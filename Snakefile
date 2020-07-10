rule all:
    input:
        auspice = "auspice/batch.json",

# Triggers the pre-analyses
rule preanalyses:
	input:
		"data/metadata.tsv",
		"data/sequences.fasta",
		"config/latlongs.tsv",
		"config/colors.tsv"


# Define file names
rule files:
	params:
		original_dataset = "pre-analyses/gisaid_hcov-19.fasta",
		new_genomes = "pre-analyses/new_genomes.fasta",
		full_metadata = "pre-analyses/metadata_nextstrain.tsv",
		metadata_lab = "pre-analyses/Sequencing-SC2.xlsx",
		cache = "config/cache_coordinates.tsv",
		keep = "config/keep.txt",
		remove = "config/remove.txt",
		reference = "config/reference.gb",
		geoscheme = "config/geoscheme.tsv",
		colour_grid = "config/colour_grid.html",
		clades = "config/clades.tsv",
		auspice_config = "config/auspice_config.json",
		dropped_strains = "config/dropped_strains.txt"



files = rules.files.params


rule add_sequences:
	message:
		"""
		Filtering sequence files to:
		- Add selected and newly sequenced genomes
		- Prevent unwanted genomes from being incorporated in the initial dataset.
		"""
	input:
		genomes = files.original_dataset,
		new_genomes = files.new_genomes,
		include = files.keep,
		exclude = files.remove
	output:
		sequences = "pre-analyses/temp_sequences.fasta"
	shell:
		"""
		python3 scripts/add_newgenomes.py \
			--genomes {input.genomes} \
			--new-genomes {input.new_genomes} \
			--keep {input.include} \
			--remove {input.exclude} \
			--output {output.sequences}
		"""


rule filter_metadata:
	message:
		"""
		Processing {input.genomes} to:
		- Filter only lines corresponding to genomes included in {input.genomes}
		- Reformat metadata by dropping or adding columns, and fixing some fields
		"""
	input:
		genomes = rules.add_sequences.output.sequences,
		metadata1 = files.full_metadata,
		metadata2 = files.metadata_lab
	output:
		filtered_metadata = "pre-analyses/metadata_filtered.tsv",
		renaming = "pre-analyses/rename.tsv",
		sequences = "data/sequences.fasta"
	shell:
		"""
		python3 scripts/filter_metadata.py \
			--genomes {input.genomes} \
			--metadata1 {input.metadata1} \
			--metadata2 {input.metadata2} \
			--output1 {output.filtered_metadata} \
			--output2 {output.renaming} \
			--output3 {output.sequences}
		"""


rule geoscheme:
	message:
		"""
		Processing {input.filtered_metadata} to:
		- Reformat columns according to geographic scheme in {input.geoscheme}
		"""
	input:
		filtered_metadata = rules.filter_metadata.output.filtered_metadata,
		geoscheme = files.geoscheme
	output:
		final_metadata = "data/metadata.tsv"
	shell:
		"""
		python3 scripts/apply_geoscheme.py \
			--metadata {input.filtered_metadata} \
			--geoscheme {input.geoscheme} \
			--output {output.final_metadata}
		"""


rule coordinates:
	message:
		"""
		Searching for coordinates (latitudes and longitudes) for samples in {input.metadata}
		"""
	input:
		metadata = rules.geoscheme.output.final_metadata,
		geoscheme = files.geoscheme,
		cache = files.cache
	params:
		columns = "region country division location"
	output:
		latlongs = "config/latlongs.tsv"
	shell:
		"""
		python3 scripts/get_coordinates.py \
			--metadata {input.metadata} \
			--geoscheme {input.geoscheme} \
			--columns {params.columns} \
			--cache {input.cache} \
			--output {output.latlongs}
		cp {output.latlongs} config/cache_coordinates.tsv
		"""


rule colours:
	message:
		"""
		Assigning colour scheme for samples in {input.metadata} based on geoscheme
		"""
	input:
		metadata = rules.geoscheme.output.final_metadata,
		latlongs = rules.coordinates.output.latlongs,
		geoscheme = files.geoscheme,
		colour_grid = files.colour_grid
	params:
		columns = "region country division location"
	output:
		colours = "config/colors.tsv"
	shell:
		"""
		python3 scripts/apply_colour_scheme.py \
		--metadata {input.metadata} \
		--coordinates {input.latlongs} \
		--geoscheme {input.geoscheme} \
		--grid {input.colour_grid} \
		--columns {params.columns} \
		--output {output.colours}
		"""






### STARTING NEXTSTRAIN PIPELINE


input_fasta = rules.filter_metadata.output.sequences,
input_metadata = rules.geoscheme.output.final_metadata,
reference = files.reference,
clades = files.clades,
lat_longs = rules.coordinates.output.latlongs,
colors = rules.colours.output.colours,
dropped_strains = files.dropped_strains,
auspice_config = "config/auspice_config.json"

### Excluding sequences included in dropped_strains
rule filter:
    message:
        """
        Filtering to
          - excluding strains in {input.exclude}
        """
    input:
        sequences = input_fasta,
        metadata = input_metadata,
        exclude = dropped_strains
    output:
        sequences = "results/filtered.fasta"
    params:
        group_by = "country year month",
        sequences_per_group = 20,
        min_date = 1900
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --output {output.sequences}
        """


### Aligning the sequences using MAFFT
rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = reference
    output:
        alignment = "results/aligned.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --remove-reference \
            --fill-gaps
        """


### Masking alignment sites
rule mask:
    message:
        """
        Mask bases in alignment
          - masking {params.mask_from_beginning} from beginning
          - masking {params.mask_from_end} from end
          - masking other sites: {params.mask_sites}
        """
    input:
        alignment = rules.align.output.alignment
    output:
        alignment = "results/masked.fasta"
    params:
        mask_from_beginning = 130,
        mask_from_end = 50,
        mask_sites = "18529 29849 29851 29853"
    shell:
        """
        python3 scripts/mask-alignment.py \
            --alignment {input.alignment} \
            --mask-from-beginning {params.mask_from_beginning} \
            --mask-from-end {params.mask_from_end} \
            --mask-sites {params.mask_sites} \
            --output {output.alignment}
        """


### Inferring Maximum Likelihood tree using the default software (IQTree)

rule tree:
    message: "Building tree"
    input:
        alignment = rules.mask.output.alignment
    output:
        tree = "results/tree_raw.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree}
        """


### Running TreeTime to estimate time for ancestral genomes

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = input_metadata
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json"
    params:
        root = "Wuhan/WH01/2019",
        coalescent = "skyline",
        clock_rate = 0.0008,
        clock_std_dev = 0.0004,
        date_inference = "marginal",
        unit = "mutations"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --root {params.root} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --clock-rate {params.clock_rate} \
            --clock-std-dev {params.clock_std_dev} \
            --divergence-units {params.unit} \
            --date-inference {params.date_inference}
        """


### Reconstructing ancestral sequences and mutations

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/nt_muts.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --inference {params.inference} \
            --output-node-data {output.node_data}
        """

## Performing amino acid translation

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = reference
    output:
        node_data = "results/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} \
        """


### Inferring ancestral locations of genomes

rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = input_metadata
    output:
        node_data = "results/traits.json",
    params:
        columns = "region"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """


### Define clades based on sets of mutations

rule clades:
    message: " Labeling clades as specified in config/clades.tsv"
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = clades
    output:
        clade_data = "results/clades.json"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output {output.clade_data}
        """


### Generating final results for visualisation with auspice

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = input_metadata,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        colors = colors,
        lat_longs = lat_longs,
        clades = rules.clades.output.clade_data,
        auspice_config = auspice_config
    output:
        auspice = rules.all.input.auspice,
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} {input.clades} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice}
        """


### Clearing the working directory (only executed when needed)

rule clean:
	message: "Removing directories: {params}"
	params:
		"results ",
		"auspice",
		"data",
		"config/colors.tsv",
		"config/latlongs.tsv",
		"pre-analyses/metadata_filtered.tsv",
		"pre-analyses/temp_sequences.fasta",
		"pre-analyses/rename.tsv"

	shell:
		"""
		rm -rfv {params}
		"""


rule delete:
	message: "Deleting directory: {params}"
	params:
		"pre-analyses"
	shell:
		"rm -rfv {params}"
