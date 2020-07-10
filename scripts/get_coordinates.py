# coding=utf-8
import pandas as pd
from geopy.geocoders import Nominatim
import argparse
from bs4 import BeautifulSoup as BS
import numpy as np

geolocator = Nominatim(user_agent="email@gmail.com")  # add your email here

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generate file with latitudes and longitudes of samples listed in a metadata file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", required=True, help="Nextstrain metadata file")
    parser.add_argument("--geoscheme", required=True, help="XML file with geographic schemes")
    parser.add_argument("--columns", nargs='+', type=str, help="list of columns that need coordinates")
    parser.add_argument("--cache", required=False, help="TSV file with preexisting latitudes and longitudes")
    parser.add_argument("--output", required=True, help="TSV file containing geographic coordinates")
    args = parser.parse_args()

    metadata = args.metadata
    geoscheme = args.geoscheme
    columns = args.columns
    cache = args.cache
    output = args.output

    # metadata = path + 'metadata_geo.tsv'
    # geoscheme = path + "geoscheme.tsv"
    # columns = ['region', 'country', 'division', 'location']
    # cache = path + 'cache_coordinates.tsv'
    # output = path + 'lat_longs.tsv'

    force_coordinates = {'Puerto Rico': ('18.235853', '-66.522056'), 'Virgin Islands': ('17.727304', '-64.748327')}

    results = {trait: {} for trait in columns}  # content to be exported as final result

    # extract coordinates from cache file
    try:
        for line in open(cache).readlines():
            if not line.startswith('\n'):
                try:
                    trait, place, lat, long = line.strip().split('\t')
                    if trait in results.keys():
                        entry = {place: (str(lat), str(long))}
                        results[trait].update(entry)  # save as pre-existing result
                except:
                    pass
    except:
        pass

    # extract coordinates from TSV file
    scheme_list = open(geoscheme, "r").readlines()[1:]
    dont_search = []
    set_countries = []
    for line in scheme_list:
        if not line.startswith('\n'):
            type = line.split('\t')[0]
            if type in columns:
                coordinates = {}
                try:
                    subarea = line.split('\t')[2]  # name of the pre-defined area in the TSV file
                    lat = line.split('\t')[3]
                    long = line.split('\t')[4]
                    entry = {subarea: (lat, long)}
                    coordinates.update(entry)

                    if subarea not in results[type]:
                        results[type].update(coordinates)
                        dont_search.append(subarea)
                    country_name = subarea.split('-')[0]
                    if type == 'country' and country_name not in set_countries:
                        set_countries.append(country_name)
                except:
                    pass


    # find coordinates for locations not found in cache or XML file
    def find_coordinates(place):
        try:
            location = geolocator.geocode(place, language='en')
            lat, long = location.latitude, location.longitude
            coord = (str(lat), str(long))
            return coord
        except:
            coord = ('NA', 'NA')
            return coord


    # open metadata file as dataframe
    dfN = pd.read_csv(metadata, encoding='utf-8', sep='\t')

    queries = []
    pinpoints = [dfN[trait].values.tolist() for trait in columns if trait != 'region']
    for address in zip(*pinpoints):
        traits = [trait for trait in columns if trait != 'region']
        for position, place in enumerate(address):
            level = traits[position]
            query = list(address[0:position + 1])
            queries.append((level, query))

    not_found = []
    for unknown_place in queries:
        trait, place = unknown_place[0], unknown_place[1]
        target = place[-1]
        if target not in ['', 'NA', 'NAN', 'unknown', '-', np.nan, None]:
            try:
                if place[0].split('-')[0] in set_countries:
                    country_short = place[0].split('-')[0]  # correcting TSV pre-defined country names
                    place[0] = country_short
            except:
                pass

            if target not in results[trait]:
                new_query = []
                for name in place:
                    if name not in dont_search:
                        if place[0] == 'USA':
                            if name != 'USA':
                                name = name + ' state'
                        if name not in new_query:
                            new_query.append(name)

                item = (trait, ', '.join(new_query))
                coord = ('NA', 'NA')
                if item not in not_found:
                    coord = find_coordinates(', '.join(new_query))  # search coordinates
                if 'NA' in coord:
                    if item not in not_found:
                        not_found.append(item)
                        print('\t* WARNING! Coordinates not found for: ' + trait + ', ' + ', '.join(new_query))
                else:
                    print(trait + ', ' + target + '. Coordinates = ' + ', '.join(coord))
                    entry = {target: coord}
                    results[trait].update(entry)

    print('\n### These coordinates were found and saved in the output file:')
    with open(output, 'w') as outfile:
        for trait, lines in results.items():
            print('\n* ' + trait)
            for place, coord in lines.items():
                if place in force_coordinates:
                    lat, long = force_coordinates[place][0], force_coordinates[place][1]
                else:
                    lat, long = coord[0], coord[1]
                print(place + ': ' + lat + ', ' + long)
                line = "{}\t{}\t{}\t{}\n".format(trait, place, lat, long)
                outfile.write(line)
            outfile.write('\n')

    if len(not_found) > 1:
        print('\n### WARNING! Some coordinates were not found (see below).'
              '\nTypos or especial characters in place names my explain such errors.'
              '\nPlease fix them, and run the script again, or add coordinates manually:\n')
        for trait, address in not_found:
            print(trait + ': ' + address)

print('\nCoordinates file successfully created!\n')
