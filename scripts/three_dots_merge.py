import argparse
import json
from copy import copy, deepcopy
from os import path
from sys import stderr
import math

__author__ = 'senorrift'
# TODO: Add Kn/Ks Support

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#                                                                                                                     #
# Parse Command Line Arguments                                                                                        #
#                                                                                                                     #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

parser = argparse.ArgumentParser(description="Merge Three SynMaps")

parser.add_argument('-i1',
                    type=str,
                    required=True,
                    help='Path to pairwise SynMap (.gcoords.ks) #1 (Full Path)')

parser.add_argument('-i2',
                    type=str,
                    required=True,
                    help='Path to pairwise SynMap (.gcoords.ks) #2 (Full Path)')

parser.add_argument('-i3',
                    type=str,
                    required=True,
                    help='Path to pairwise SynMap (.gcoords.ks) #3 (Full Path)')

parser.add_argument('-o',
                    type=str,
                    required=True,
                    help="Path to output data folder")

parser.add_argument('-xid',
                    type=str,
                    required=True,
                    help="X Axis Genome ID")

parser.add_argument('-yid',
                    type=str,
                    required=True,
                    help="Y Axis Genome ID")

parser.add_argument('-zid',
                    type=str,
                    required=True,
                    help="Z Axis Genome ID")

parser.add_argument('-P',
                    action="store_true",
                    help="Parse (remove chromosomes w/o synteny)")

parser.add_argument('-M',
                    type=int,
                    help="Minimum Scaffold Length")

args = parser.parse_args()

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#                                                                                                                     #
# Functions                                                                                                           #
#                                                                                                                     #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #


def chromosome_parser(output_path, xy_input, xz_input, yz_input):
    """ Chromosome Parser
    Given three pairwise SynMap dotplot_dot.pl output JSON files,
    Generate new JSON files containing only chromosomes with synteny.

    :param xy_input: JSON file from pairwise SynMap & dotplot_dots.pl
    :param xz_input: JSON file from pairwise SynMap & dotplot_dots.pl
    :param yz_input: JSON file from pairwise SynMap & dotplot_dots.pl
    """
    from json import load, dump
    # Global Variables
    inputs = [xy_input, xz_input, yz_input]
    match_chromosomes = {}

    # Data Processing
    for input_file in inputs:
        get = load(open(input_file, 'r'))

        # Get Species IDs
        genomes = get["syntenic_points"]
        first_sp_id = genomes.keys()[0]
        second_sp_id = genomes[first_sp_id].keys()[0]

        if first_sp_id not in match_chromosomes.keys():
            match_chromosomes[first_sp_id] = []
        if second_sp_id not in match_chromosomes.keys():
            match_chromosomes[second_sp_id] = []

        # Find Chromosomes with Matches
        pairs = genomes[first_sp_id][second_sp_id]
        for first_sp_ch in pairs:
            if first_sp_ch not in match_chromosomes[first_sp_id]:
                match_chromosomes[first_sp_id].append(first_sp_ch)
            for second_sp_ch in pairs[first_sp_ch]:
                if second_sp_ch not in match_chromosomes[second_sp_id]:
                    match_chromosomes[second_sp_id].append(second_sp_ch)

    # Produce Parsed Files
    for input_file in inputs:
        output_file = "%s/parsed_%s" % (output_path, path.basename(input_file))
        parsed = load(open(input_file, 'r'))

        to_remove = {}
        for genome in parsed['genomes']:
            to_remove[genome] = []
            for chromosome in parsed['genomes'][genome]['chromosomes']:
                if chromosome['name'] not in match_chromosomes[genome]:
                    to_remove[genome].append(chromosome)

        for genome in to_remove:
            for removal in to_remove[genome]:
                parsed['genomes'][genome]['chromosomes'].remove(removal)
                # print "%s: %s removed from %s" % (input_file, removal, genome)

        dump(parsed, open(output_file, 'w'))


def get_data(first_sp_id, second_sp_id, json_file):
    """ Get Data
    Given two species IDs and a URL to the SynMap comparison JSON between them,
    Return a list of the (a, b) coordinates, the a coordinates, and the b coordinates.

    :param first_sp_id: Species ID for axis "a".
    :param second_sp_id: Species ID for axis "b".
    :param json_file: URL to SynMap comparison JSON.
    :return: lists of (a, b) coordinates, a coordinates, b coordinates.
    """

    # -----------------------------------------------------------------------------------------------------------------
    # Load JSON into Dictionary
    # -----------------------------------------------------------------------------------------------------------------
    get = json.load(open(json_file, 'r'))

    # -----------------------------------------------------------------------------------------------------------------
    # Determine Species Order
    # -----------------------------------------------------------------------------------------------------------------
    species = []
    json_sps = path.basename(json_file).lstrip("parsed_").split("_")
    if first_sp_id in json_sps and second_sp_id in json_sps:
        species = [json_sps[0], json_sps[1]]
    else:
        stderr.write("File/Species Error: Code 1")
        exit()

    # -----------------------------------------------------------------------------------------------------------------
    # Build Data Structures
    #
    #   ######## A Match Locations ########                   ######## B Match Locations ########
    #   # ab_a = {a_ch1: {b_ch1: [a, a, ...],                 # ab_b = {b_ch1: {a_ch1: [b, b, ...],
    #   #                 b_ch2: [a, a, ...],                 #                 a_ch2: [b, b, ...],
    #   #                 ...                                 #                 ...
    #   #                 b_chN: [a, a, ...]},                #                 a_chN: [b, b, ...]},
    #   #         a_ch2: {b_ch1: [a, a, ...],                 #         b_ch2: {a_ch1: [b, b, ...],
    #   #                 b_ch2: [a, a, ...],                 #                 a_ch2: [b, b, ...],
    #   #                 ...                                 #                 ...
    #   #                 b_chN: [a, a, ...]},                #                 a_chN: [b, b, ...]},
    #   #         ...                                         #         ...
    #   #         a_chN: {b_ch1: [a, a, ...],                 #         b_chN: {a_ch1: [b, b, ...],
    #   #                 b_ch2: [a, a, ...],                 #                 a_ch2: [b, b, ...],
    #   #                 ...                                 #                 ...
    #   #                 b_chN: [a, a, ...]}}                #                 a_chN: [b, b, ...]}}
    #
    #   ######## Coordinates A-Indexed ########               ######## Coordinates B-Indexed ########
    #   # ab_cords = {a_ch1: {b_ch1: [[a,b], [a,b], ...],     # ba_cords = {b_ch1: {a_ch1: [[a,b], [a,b], ...],
    #   #                     b_ch2: [[a,b], [a,b], ...],     #               a_ch2: [[a,b], [a,b], ...],
    #   #                     ...                             #               ...
    #   #                     b_chN: [[a,b], [a,b], ...]},    #               a_chN: [[a,b], [a,b], ...]},
    #   #             ...                                     #       ...
    #   #             a_chN: {b_ch1: [[a,b], [a,b], ...],     #       b_chN: {a_ch1: [[a,b], [a,b], ...]
    #   #                     ...                             #               ...
    #   #                     b_chN: [[a,b], [a,b], ...]}}    #               a_ch2: [[a,b], [a,b], ...]}}
    #
    # -----------------------------------------------------------------------------------------------------------------

    ab_a = {}
    ab_b = {}
    ab_cords = {}
    ba_cords = {}

    a_chromosomes = get["genomes"][first_sp_id]["chromosomes"]
    b_chromosomes = get["genomes"][second_sp_id]["chromosomes"]

    for a_ch in a_chromosomes:
        ab_a[a_ch["name"]] = {}
        ab_cords[a_ch["name"]] = {}
        for b_ch in b_chromosomes:
            ab_a[a_ch["name"]][b_ch["name"]] = []
            ab_cords[a_ch["name"]][b_ch["name"]] = []

    for b_ch in b_chromosomes:
        ab_b[b_ch["name"]] = {}
        ba_cords[b_ch["name"]] = {}
        for a_ch in a_chromosomes:
            ab_b[b_ch["name"]][a_ch["name"]] = []
            ba_cords[b_ch["name"]][a_ch["name"]] = []

    # -----------------------------------------------------------------------------------------------------------------
    # Get Data
    # -----------------------------------------------------------------------------------------------------------------
    data = get["syntenic_points"][species[0]][species[1]]
    # Iterate through all matches
    for sp1_chr in data:
        for sp2_chr in data[sp1_chr]:
            for match in data[sp1_chr][sp2_chr]:
                # Calculate points for first and second species in JSON
                sp1_start = int(match[0])
                sp1_end = int(match[1])
                sp1_cord = (sp1_start + sp1_end) / 2

                sp2_start = int(match[2])
                sp2_end = int(match[3])
                sp2_cord = (sp2_start + sp2_end) / 2

                kn_ks = match[4]

                # Assign correct genome/value by JSON format
                if species[0] == first_sp_id:
                    a_chr = sp1_chr
                    b_chr = sp2_chr
                    a_cord = sp1_cord
                    b_cord = sp2_cord
                elif species[0] == second_sp_id:
                    a_chr = sp2_chr
                    b_chr = sp1_chr
                    a_cord = sp2_cord
                    b_cord = sp1_cord
                else:
                    stderr.write("File/Species Error: Code 3")
                    exit()

                # Populate data sets
                coordinate = [a_cord, b_cord, kn_ks]
                ab_a[a_chr][b_chr].append(a_cord)
                ab_b[b_chr][a_chr].append(b_cord)
                ab_cords[a_chr][b_chr].append(coordinate)
                ba_cords[b_chr][a_chr].append(coordinate)

    return ab_a, ab_b, ab_cords, ba_cords


def find_matches(species_coordinate, link1, link2, link3):
    """ Find Matches
    Given a list of three species IDs (in order of X, Y, Z axis) and three pairwise comparisons covering all three,
    Returns a list of the coordinates for those matches shared between all three species.

    :param species_coordinate: list containing species ID in order [x axis species, y axis species, z axis species].
    :param link1: Link to a file with pairwise SynMap JSON.
    :param link2: Link to a file with pairwise SynMap JSON.
    :param link3: Link to a file with pairwise SynMap JSON.
    :return: Match coordinate dictionary {X Chr: {Y Chr: {ZChr : [[x_loc, y_loc, z_loc], [x_loc, y_loc, z_loc], ...]}}}
    """

    # Load genome data from each linked file
    link1_data = json.load(open(link1, 'r'))["genomes"]
    link2_data = json.load(open(link2, 'r'))["genomes"]
    link3_data = json.load(open(link3, 'r'))["genomes"]

    # Determine species contained within each file
    link1_spp = link1_data.keys()
    link2_spp = link2_data.keys()
    link3_spp = link3_data.keys()

    # Establish Correct Links/Species
    x_species_id = species_coordinate[0]
    y_species_id = species_coordinate[1]
    z_species_id = species_coordinate[2]
    xy_link = ''
    xz_link = ''
    yz_link = ''

    if x_species_id in link1_spp and y_species_id in link1_spp:
        xy_link = link1
        xy_data = link1_data
    elif x_species_id in link2_spp and y_species_id in link2_spp:
        xy_link = link2
        xy_data = link2_data
    elif x_species_id in link3_spp and y_species_id in link3_spp:
        xy_link = link3
        xy_data = link3_data

    if y_species_id in link1_spp and z_species_id in link1_spp:
        yz_link = link1
        yz_data = link1_data
    elif y_species_id in link2_spp and z_species_id in link2_spp:
        yz_link = link2
        yz_data = link2_data
    elif y_species_id in link3_spp and z_species_id in link3_spp:
        yz_link = link3
        yz_data = link3_data

    if x_species_id in link1_spp and z_species_id in link1_spp:
        xz_link = link1
    elif x_species_id in link2_spp and z_species_id in link2_spp:
        xz_link = link2
    elif x_species_id in link3_spp and z_species_id in link3_spp:
        xz_link = link3

    # Get coordinate information from each file
    xy_x, xy_y, xy_cords, yx_cords = get_data(x_species_id, y_species_id, xy_link)
    xz_x, xz_z, xz_cords, zx_cords = get_data(x_species_id, z_species_id, xz_link)
    yz_y, yz_z, yz_cords, zy_cords = get_data(y_species_id, z_species_id, yz_link)

    # Build Data Structure To Hold 3-Way Matches
    matches = {}
    for x_chr in xy_data[x_species_id]["chromosomes"]:
        matches[x_chr["name"]] = {}
        for y_chr in xy_data[y_species_id]["chromosomes"]:
            matches[x_chr["name"]][y_chr["name"]] = {}
            for z_chr in yz_data[z_species_id]["chromosomes"]:
                matches[x_chr["name"]][y_chr["name"]][z_chr["name"]] = []

    # Data Structure to Hold Histogram Information
    hist_data = {"Ks": {"xy": [], "xz": [], "yz": [], 'mean': []},
                 "Kn": {"xy": [], "xz": [], "yz": [], 'mean': []}}

    # Find 3-Way Matches
    match_count = 0
    for x_ch in xy_cords:
        for y_ch in xy_cords[x_ch]:
            for xy_match in xy_cords[x_ch][y_ch]:
                for z_ch in xz_cords[x_ch]:
                    for xz_match in xz_cords[x_ch][z_ch]:
                        if xy_match[0] == xz_match[0]:
                            for y_ch in zy_cords[z_ch]:
                                for yz_match in zy_cords[z_ch][y_ch]:
                                    if yz_match[0] == xy_match[1]:
                                        if yz_match[1] == xz_match[1]:
                                            # Build coordinate template.
                                            coordinate = [0, 0, 0, {}]
                                            # Populate XYZ coordinate fields.
                                            coordinate[0] = xy_match[0]
                                            coordinate[1] = xy_match[1]
                                            coordinate[2] = xz_match[1]
                                            # Populate Ks/Kn coordinate field.
                                            coordinate[3]['xy'] = xy_match[2]
                                            coordinate[3]['xz'] = xz_match[2]
                                            coordinate[3]['yz'] = yz_match[2]
                                            if 'NA' in xy_match[2].values() or 'NA' in xz_match[2].values() or \
                                                            'NA' in yz_match[2].values():
                                                meanKs = 'NA'
                                                meanKn = 'NA'
                                            else:
                                                meanKs = str((float(xy_match[2]['Ks']) + float(xz_match[2]['Ks']) +
                                                              float(yz_match[2]['Ks'])) / 3)
                                                meanKn = str((float(xy_match[2]['Kn']) + float(xz_match[2]['Kn']) +
                                                              float(yz_match[2]['Kn'])) / 3)
                                            coordinate[3]['mean'] = {'Ks': meanKs, 'Kn': meanKn}
                                            # Add coordinate to matches, ignoring duplicates.
                                            if coordinate not in matches[x_ch][y_ch][z_ch]:
                                                matches[x_ch][y_ch][z_ch].append(coordinate)
                                                match_count += 1
                                            # Populate Ks field in "hist_data".
                                            try:
                                                hist_data['Ks']['xy'].append(float(xy_match[2]['Ks']))
                                                hist_data['Ks']['xz'].append(float(xz_match[2]['Ks']))
                                                hist_data['Ks']['yz'].append(float(yz_match[2]['Ks']))
                                                hist_data['Ks']['mean'].append(float(meanKs))
                                            except (TypeError, ValueError):
                                                pass
                                            # Populate Kn field in "hist_data".
                                            try:
                                                hist_data['Kn']['xy'].append(float(xy_match[2]['Kn']))
                                                hist_data['Kn']['xz'].append(float(xz_match[2]['Kn']))
                                                hist_data['Kn']['yz'].append(float(yz_match[2]['Kn']))
                                                hist_data['Kn']['mean'].append(float(meanKn))
                                            except (TypeError, ValueError):
                                                pass

    # Remove any empty elements
    # matches = {xchr: { ychr: { zchr: [matches], ...}, ...} ...}
    for xchr in matches:
        for ychr in matches[xchr]:
            matches[xchr][ychr] = {k: v for k, v in matches[xchr][ychr].items() if len(v) > 0}

    for xchr in matches:
        matches[xchr] = {k: v for k, v in matches[xchr].items() if len(v) > 0}

    matches = {k: v for k, v in matches.items() if len(v) > 0}

    # Return Matches, Match Count, and Histogram Data
    return matches, match_count, hist_data


def generate_histogram(histo_data):
    # Log transform data
    log_data = {}
    for calculation in histo_data:
        log_data[calculation] = {}
        for comparison in histo_data[calculation]:
            log_data[calculation][comparison] = [math.log(n, 10) for n in histo_data[calculation][comparison] if n != 0]

    # Define 'info' sub-template.
    info_subtmpl = {'min': 0,
                    'max': 0,
                    'dif': 0,
                    'bin_size': 0,
                    'bin1': [0, 0],
                    'bin2': [0, 0],
                    'bin3': [0, 0],
                    'bin4': [0, 0],
                    'bin5': [0, 0]}
    # Define 'info' template.
    info_tmpl = {'xy': deepcopy(info_subtmpl),
                 'xz': deepcopy(info_subtmpl),
                 'yz': deepcopy(info_subtmpl),
                 'mean': deepcopy(info_subtmpl)}

    # Define 'log_info' template.
    log_info = {'Ks': deepcopy(info_tmpl),
                'Kn': deepcopy(info_tmpl)}

    # Define 'value_info' template.
    value_info = deepcopy(log_info)

    # Define a function to populate info dictionaries.
    def populate_info(info_dict, data_dict):
        for calc in data_dict:
            for comp in data_dict[calc]:
                # print calc, comp, min(log_data[calc][comp]), max(log_data[calc][comp])
                info_dict[calc][comp]['min'] = min(data_dict[calc][comp])
                info_dict[calc][comp]['max'] = max(data_dict[calc][comp])
                info_dict[calc][comp]['dif'] = abs(max(data_dict[calc][comp]) - min(data_dict[calc][comp]))
                info_dict[calc][comp]['bin_size'] = abs(max(data_dict[calc][comp]) - min(data_dict[calc][comp])) / 5
                info_dict[calc][comp]['bin1'][0] = min(data_dict[calc][comp])
                info_dict[calc][comp]['bin1'][1] = info_dict[calc][comp]['bin1'][0] + info_dict[calc][comp]['bin_size']
                info_dict[calc][comp]['bin2'][0] = copy(info_dict[calc][comp]['bin1'][1])
                info_dict[calc][comp]['bin2'][1] = info_dict[calc][comp]['bin2'][0] + info_dict[calc][comp]['bin_size']
                info_dict[calc][comp]['bin3'][0] = copy(info_dict[calc][comp]['bin2'][1])
                info_dict[calc][comp]['bin3'][1] = info_dict[calc][comp]['bin3'][0] + info_dict[calc][comp]['bin_size']
                info_dict[calc][comp]['bin4'][0] = copy(info_dict[calc][comp]['bin3'][1])
                info_dict[calc][comp]['bin4'][1] = info_dict[calc][comp]['bin4'][0] + info_dict[calc][comp]['bin_size']
                info_dict[calc][comp]['bin5'][0] = copy(info_dict[calc][comp]['bin4'][1])
                info_dict[calc][comp]['bin5'][1] = max(data_dict[calc][comp])
        return info_dict

    # Populate 'log_info'.
    log_info = populate_info(log_info, log_data)
    # Populate 'value_info'.
    value_info = populate_info(value_info, histo_data)

    # Build final histogram dictionary.
    histogram_dict = {'values': {'data': histo_data, 'info': value_info},
                      'logten': {'data': log_data, 'info': log_info}}

    return histogram_dict


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#                                                                                                                     #
# Main Script                                                                                                         #
#                                                                                                                     #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# User Inputs
x = args.xid
y = args.yid
z = args.zid
file1 = args.i1
file2 = args.i2
file3 = args.i3
ofolder = args.o.rstrip('/')
synt_out = "%s/%s_%s_%s.json" % (ofolder, x, y, z)
hist_out = "%s/%s_%s_%s_histogram.json" % (ofolder, x, y, z)

parse_file = args.P
min_length = args.M

# Option Actions
if parse_file:
    chromosome_parser(ofolder, file1, file2, file3)
    file1 = "%s/parsed_%s" % (ofolder, path.basename(file1))
    file2 = "%s/parsed_%s" % (ofolder, path.basename(file2))
    file3 = "%s/parsed_%s" % (ofolder, path.basename(file3))
    synt_out = "%s/parsed_%s_%s_%s.json" % (ofolder, x, y, z)
    hist_out = "%s/parsed_%s_%s_%s_histogram.json" % (ofolder, x, y, z)

if min_length is not None and min_length > 0:
    # TODO: Make Minimum Length Parser
    pass

# Execute match finding.
all_matches, match_number, histogram_data = find_matches([x, y, z], file1, file2, file3)

# Parse histogram_data
histogram = generate_histogram(histogram_data)

# Dump JSONs - matches .
json.dump(all_matches, open(synt_out, "wb"))
json.dump(histogram, open(hist_out, "wb"))

# Build Log & Dump to JSON
log_ks = histogram['logten']['data']['Ks']['mean']
log_kn = histogram['logten']['data']['Kn']['mean']

log = {}
log_out = "%s/log.json" % ofolder
log["status"] = "complete"
log["matches"] = str(match_number)
log["ks"] = {"min": str(min(log_ks)), "max": str(max(log_ks))}
log["kn"] = {"min": str(min(log_kn)), "max": str(max(log_kn))}
json.dump(log, open(log_out, 'wb'))
