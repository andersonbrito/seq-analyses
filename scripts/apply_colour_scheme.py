# coding=utf-8
import pycountry_convert as pyCountry
import pandas as pd
import argparse
from bs4 import BeautifulSoup as BS
import pycountry
from matplotlib import cm
import numpy as np


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generate ordered colour file for nextstrain build",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", required=True, help="Reformatted nextstrain metadata file")
    parser.add_argument("--coordinates", required=True,  help="TSV coordinates file being used in the build")
    parser.add_argument("--geoscheme", required=True, help="XML file with geographic scheme")
    parser.add_argument("--grid", required=True, help="HTML file with HEX colour matrices")
    parser.add_argument("--columns", nargs='+', type=str,   help="list of columns with geographic information")
    parser.add_argument("--output", required=True, help="TSV file containing ordered HEX colours based on locations")
    args = parser.parse_args()

    metadata = args.metadata
    geoscheme = args.geoscheme
    coordinates = args.coordinates
    grid = args.grid
    columns = args.columns
    output = args.output


    # metadata = path + 'metadata_geo.tsv'
    # coordinates = path + 'lat_longs.tsv'
    # geoscheme = path + 'geoscheme.tsv'
    # grid = path + 'colour_grid.html'
    # columns = ['region', 'country', 'division', 'location']
    # output = path + 'colors.tsv'


    # pre-determined HEX colours and hues
    # force_colour = {'Family': '#0f8a3e', 'Florida': '#743d8f', 'International': '#bfbfbf', 'Staff-IN': '#0068B1', 'Staff-OUT': '#66bfff',
    #                 'Team': '#E61649', 'USA-Midwest': '#f0db75', 'USA-Northeast': '#ebcf47', 'USA-South': '#e6c419', 'USA-West': '#f5e7a3'}
    force_colour = {}
    force_hue = {'North America': 0}

    # content to be exported as final result
    latlongs = {trait: {} for trait in columns}

    # get ISO alpha3 country codes
    isos = {}
    def get_iso(country):
        global isos
        if country not in isos.keys():
            try:
                isoCode = pyCountry.country_name_to_country_alpha3(country, cn_name_format="default")
                isos[country] = isoCode
            except:
                try:
                    isoCode = pycountry.countries.search_fuzzy(country)[0].alpha_3
                    isos[country] = isoCode
                except:
                    isos[country] = ''
        return isos[country]


    # extract coordinates from latlongs file for sorting places by latitude
    for line in open(coordinates).readlines():
        if not line.startswith('\n'):
            try:
                trait, place, lat, long = line.strip().split('\t')
                if trait in latlongs.keys():
                    entry = {place: (str(lat), str(long))}
                    latlongs[trait].update(entry) # save as pre-existing result
                else:
                    print('### WARNING! ' + trait + ' is not among the pre-selected columns with geographic information!')
            except:
                pass



    ''' REORDER LOCATIONS FOR LEGEND FORMATTING '''

    # open metadata file as dataframe
    dfN = pd.read_csv(metadata, encoding='ISO-8859-1', sep='\t')

    ordered_regions = {}
    dcountries = {}
    places = []
    pinpoints = [dfN[trait].values.tolist() for trait in columns]
    for region_index in latlongs[columns[0]]:
        for address in zip(*pinpoints):
            address = list(address)
            region, country, division, location = [trait for trait in columns]
            region_name = address[0]
            country_name = address[1]
            division_name = address[2]
            location_name = address[3]
            places.append(address)

            if region_index == region_name:
                if 'subcontinent' not in ordered_regions.keys():
                    ordered_regions['subcontinent'] = {}
                if region_name not in ordered_regions['subcontinent'].keys():
                    ordered_regions['subcontinent'].update({region_name: latlongs[region][region_name]})

                ### STORE COUNTRIES WITH REGIONS AS KEYS
                if country_name in latlongs[country].keys():
                    if region_name not in dcountries.keys():
                        dcountries[region_name] = {country_name: latlongs[country][country_name]}
                    else:
                        if country_name not in dcountries[region_name].keys():
                            dcountries[region_name].update({country_name: latlongs[country][country_name]})



    # sort division entries based on sorted country entries
    ordered_countries = {}
    for region, countries in dcountries.items():
        ordered_countries[region] = {k: v for k, v in sorted(countries.items(), key=lambda item: item[1])}

    ddivisions = {}
    for country_index in [key for dict_ in ordered_countries.values() for key in dict_]:
        for address in places:
            region, country, division, location = [trait for trait in columns]
            country_name = address[1]
            division_name = address[2]
            location_name = address[3]

            if country_index == country_name:
                ### STORE DIVISIONS WITH COUNTRIES AS KEYS
                if division_name in latlongs[division].keys():
                    if country_name not in ddivisions.keys():
                        ddivisions[country_name] = {division_name: latlongs[division][division_name]}
                    else:
                        if division_name not in ddivisions[country_name].keys():
                            ddivisions[country_name].update({division_name: latlongs[division][division_name]})


    # sort division entries based on sorted country entries
    ordered_divisions = {}
    for country, divisions in ddivisions.items():
        ordered_divisions[country] = {k: v for k, v in sorted(divisions.items(), key=lambda item: item[1])}



    dlocations = {}
    for division_index in [key for dict_ in ordered_divisions.values() for key in dict_]:
        for address in places:
            address = list(address)
            region, country, division, location = [trait for trait in columns]
            division_name = address[2]
            location_name = address[3]

            if division_index == division_name:
                ### STORE LOCATIONS WITH DIVISIONS AS KEYS
                if location_name in latlongs[location].keys():
                    if division_name not in dlocations.keys():
                        dlocations[division_name] = {location_name: latlongs[location][location_name]}
                    else:
                        if location_name not in dlocations[division_name].keys():
                            dlocations[division_name].update({location_name: latlongs[location][location_name]})

    # sort locations entries based on sorted division entries
    ordered_locations = {}
    for division, locations in dlocations.items():
        ordered_locations[division] = {k: v for k, v in sorted(locations.items(), key=lambda item: item[1])}



    ''' CONVERT RGB TO HEX COLOURS '''

    # original source: https://bsou.io/posts/color-gradients-with-python

    # convert colour codes
    def RGB_to_hex(RGB):
        ''' [255,255,255] -> "#FFFFFF" '''
        # Components need to be integers for hex to make sense
        RGB = [int(x) for x in RGB]
        return "#" + "".join(["0{0:x}".format(v) if v < 16 else
                              "{0:x}".format(v) for v in RGB])

    def hex_to_RGB(hex):
        ''' "#FFFFFF" -> [255,255,255] '''
        # Pass 16 to the integer function for change of base
        return [int(hex[i:i + 2], 16) for i in range(1, 6, 2)]


    def color_dict(gradient):
        ''' Takes in a list of RGB sub-lists and returns dictionary of
          colors in RGB and hex form for use in a graphing function
          defined later on '''
        return [RGB_to_hex(RGB) for RGB in gradient]

    # create gradient
    def linear_gradient(start_hex, finish_hex, n):
        ''' returns a gradient list of (n) colors between
          two hex colors. start_hex and finish_hex
          should be the full six-digit color string,
          ("#FFFFFF") '''
        # Starting and ending colors in RGB form
        s = hex_to_RGB(start_hex)
        f = hex_to_RGB(finish_hex)
        # Initilize a list of the output colors with the starting color
        RGB_list = [s]
        # Calcuate a color at each evenly spaced value of t from 1 to n
        n = n
        for t in range(1, n):
            # Interpolate RGB vector for color at the current value of t
            curr_vector = [
                int(s[j] + (float(t) / (n - 1)) * (f[j] - s[j]))
                for j in range(3)]
            # Add it to list of output colors
            RGB_list.append(curr_vector)
        RGB_list = RGB_list
        return color_dict(RGB_list)




    ''' IMPORT GEOSCHEME '''

    # xml = BS(open(geoscheme, "r").read(), 'xml')
    # levels = xml.find('levels')

    scheme_list = open(geoscheme, "r").readlines()[1:]
    sampled_region = [key for dict_ in ordered_regions.values() for key in dict_]
    geodata = {}

    for line in scheme_list:
        if not line.startswith('\n'):
            type = line.split('\t')[0]
            if type == 'region':
                continent = line.split('\t')[1]
                region = line.split('\t')[2]
                if region in sampled_region:
                    if continent not in geodata.keys():
                        geodata[continent] = [region]
                    else:
                        geodata[continent] += [region]



    ''' IMPORT COLOUR SCHEME '''

    print('\nGenerating colour scheme...\n')
    html = BS(open(grid, "r").read(), 'html.parser')
    tables = html.find_all('table')

    # for colour_name in colour_scale.keys():
    limits = {'dark': (60, 40), 'light': (100, 80)}  # define saturation and luminance, max and min, respectively
    hue_to_hex = {}
    for table in tables:
        string_head = str(table.caption)
        if string_head == 'None':
            continue
        else:
            hue_value = int(string_head.split('"')[1])

        lux = table.tbody.find_all('tr')
        lum_value = 90  # brightest colour
        hexdark = ''
        hexligth = ''
        for row in lux:
            sat_value = 10  # unsaturated colour
            for cell in row.find_all('td'):
                # pass
                if sat_value == limits['dark'][0] and lum_value == limits['dark'][1]:
                    hexdark = cell.text.strip()
                if sat_value == limits['light'][0] and lum_value == limits['light'][1]:
                    hexligth = cell.text.strip()
                hex = (hexdark, hexligth)
                sat_value += 10
            lum_value -= 10
        hex = (hexdark, hexligth)
        hue_to_hex[hue_value] = hex


    colour_scale = {'magenta': [320, 310], 'purple': [300, 290, 280, 270, 260],
                    'blue': [250, 240, 230, 220], 'cyan': [210, 200], 'turquoise': [190, 180, 170, 160, 150],
                    'green': [140, 130, 120, 110], 'yellowgreen': [100, 90, 80, 70],
                    'yellow': [60, 50, 40], 'orange': [30, 20], 'red': [10, 0]}

    continent_hues = {'Oceania': colour_scale['magenta'], 'Asia': colour_scale['purple'],
                      'Europe': colour_scale['blue'] + colour_scale['cyan'], 'Africa': colour_scale['yellowgreen'],
                      'America': colour_scale['yellow'] + colour_scale['orange'] + colour_scale['red']}


    categories = len(sampled_region)
    colour_wheel = {}
    c = 0
    palette = []
    for continent, regions in geodata.items():
        categories = len(regions)
        hues = len(continent_hues[continent])
        for region in sampled_region:
            for position, region in zip([x for x in range(1, int(hues), int(hues/categories))], regions):
                hue = continent_hues[continent][position] # colour picker
                if hue not in palette:
                    palette.append(hue)


    for region, hue in zip(sampled_region, palette):
        if region in force_hue:
            colour_wheel[region] = hue_to_hex[force_hue[region]]
        else:
            colour_wheel[region] = hue_to_hex[hue]



    ''' SET COLOUR SCHEME FOR UPDATES '''

    # convert a hue value into an rgb colour, and then hex colour
    def hue_to_rgb(hue):
        colour = int(hue / 240 * 255)
        # print(colour)
        luminance = 0.7
        rgb = [c * 255 * luminance for c in list(cm.jet(colour))[:-1]]
        # print(rgb)
        return RGB_to_hex(rgb)

    # special colours
    dfN['batch'].fillna('X', inplace=True)
    list_updates = [up_number for up_number in sorted(set(dfN['batch'].to_list())) if up_number != 'X']
    list_hex = list([hue_to_rgb(int(x)) for x in np.linspace(30, 240, len(list_updates)*2, endpoint=True)])
    skip_hex = [h for n, h in enumerate(list_hex) if n in range(0, len(list_hex), 2)]

    results = {trait: {} for trait in columns}
    results['batch'] = {}
    for update, hex in zip(list_updates, skip_hex):
        results['batch'].update({update: hex})
        print(update, hex)


    ''' APPLY SAME HUE FOR MEMBERS OF THE SAME SUB-CONTINENT '''

    # assign countries to regions
    country_colours = {}
    reference_countries = {}
    for region, members in ordered_countries.items():
        countries = list(members.keys())
        for country in countries:
            reference_countries[country] = region
            country_colours[region] = countries


    # assign divisions to countries
    division_colours = {}
    reference_divisions = {}
    for country, members in ordered_divisions.items():
        divisions = list(members.keys())
        for division in divisions:
            reference_divisions[division] = reference_countries[country]
            if reference_countries[country] not in division_colours.keys():
                division_colours[reference_countries[country]] = divisions
            else:
                if division not in division_colours[reference_countries[country]]:
                    division_colours[reference_countries[country]].append(division)


    # assign locations to divisions
    location_colours = {}
    reference_locations = {}
    for division, members in ordered_locations.items():
        locations = list(members.keys())
        for location in locations:
            if reference_divisions[division] not in location_colours.keys():
                location_colours[reference_divisions[division]] = locations
            else:
                if location not in location_colours[reference_divisions[division]]:
                    location_colours[reference_divisions[division]].append(location)



    ''' CREATE COLOUR GRADIENT '''

    # define gradients for regions
    for continent, regions in geodata.items():
        hex_limits = []
        for region in sampled_region:
            if region in regions:
                if continent == 'America':
                    hex_limits += [list(colour_wheel[region])[0]]
                else:
                    hex_limits += list(colour_wheel[region])

        start, end = hex_limits[0], hex_limits[-1]
        if len(regions) == 1:
            gradient = linear_gradient(start, end, 4)
            gradient = [list(gradient)[2]]
        else:
            gradient = linear_gradient(start, end, len(regions))

        for region, colour in zip(regions, gradient):
            print('region', region, colour)
            results['region'].update({region: colour})

    # define gradients for country
    for hue, countries in country_colours.items():
        start, end = colour_wheel[hue]
        if len(countries) == 1:
            gradient = linear_gradient(start, end, 4)
            gradient = [list(gradient)[2]]
        else:
            gradient = linear_gradient(start, end, len(countries))
        for country, colour in zip(countries, gradient):
            print('country', country, colour)
            results['country'].update({country: colour})

    # define gradients for divisions
    for hue, divisions in division_colours.items():
        start, end = colour_wheel[hue]
        if len(divisions) == 1:
            gradient = linear_gradient(start, end, 3)
            gradient = [list(gradient)[1]]
        else:
            gradient = linear_gradient(start, end, len(divisions))
        for division, colour in zip(divisions, gradient):
            print('division', division, colour)
            results['division'].update({division: colour})

    # define gradients for locations
    for hue, locations in location_colours.items():
        start, end = colour_wheel[hue]
        if len(locations) == 1:
            gradient = linear_gradient(start, end, 3)
            gradient = [list(gradient)[1]]
        else:
            gradient = linear_gradient(start, end, len(locations))
        for location, colour in zip(locations, gradient):
            print('location', location, colour)
            results['location'].update({location: colour})

    # special colouring, hardcoded
    category = {'International': '#BFBFBF', 'USA': '#ECD26C', 'Florida': '#eca06a', 'Family': '#0F8A3E',
                'Staff-IN': '#0068B1', 'Staff-OUT': '#66BFFF', 'Team': '#E61649'}
    results['category'] = {}
    for cat, hex in category.items():
        results['category'].update({cat: hex})
        print('category', cat, hex)


    ''' EXPORT COLOUR FILE '''
    with open(output, 'w') as outfile:
        for trait, entries in results.items():
            for item, hexcolour in entries.items():
                if item in force_colour and trait not in ['location']:
                    hexcolour = force_colour[item]
                    print('* ' + item + ' is hardcode with the colour ' + hexcolour)
                line = "{}\t{}\t{}\n".format(trait, item, hexcolour.upper())
                outfile.write(line)
            outfile.write('\n')
    print('\nOrdered colour file successfully created!\n')
