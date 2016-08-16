#!/usr/bin/env python

__author__ = 'senorrift'

# Dependencies:
from argparse import ArgumentParser
from datetime import datetime, timedelta
from json import dump, load
from jwt import encode
from requests import get
from sys import stderr

# ---------------------------------------------------------------------------------------------------------------------
# Use Instructions
# python dotplot_dots.py --input=<input_file> --config=<path_to_config_file> --output=<output_filename_including_path> \
# --log=<log_filename_including_path> --gid1=<id1> --gid2=<id2>
#
# 3 *REQUIRED* Command Line Arguments
# --input  : SynMap Ks Blocks File (NOTE: ends with .aligncoords.gcoords.ks)
# --config : Configuration file that includes username, secret and api url
# --gid1   : Genome ID 1
# --gid2   : Genome ID 2
# --output : Output filename
# --log    : Log filename
# ---------------------------------------------------------------------------------------------------------------------

log = {}

# Parse arguments
parser = ArgumentParser()
parser.add_argument("--input", type=str, required=True, help="Input *.aligncoords.gcoords.ks file")
parser.add_argument("--config", type=str, required=True, help="Configuration file with username, API URL, and secret")
parser.add_argument("--gid1", type=str, required=True, help="CoGe Genome ID (1/2)")
parser.add_argument("--gid2", type=str, required=True, help="CoGe Genome ID (2/2)")
parser.add_argument("--output", type=str, required=True, help="Output filename")
parser.add_argument("--log", type=str, required=True, help="Log filename")
args = parser.parse_args()

# Load config file, assign user and api_base. Exit (non-zero) if unable to open or parse.
try:
    config_file = open(args.config)
    cg = load(config_file)
    config_file.close()
    user = cg["username"]
    api_base = cg["api_url"].strip().rstrip('/') + '/'
except Exception as e:
    print("dotplot_dots.py failed (Error: unable to open or parse config file) %s" % str(e))
    stderr.write("dotplot_dots.py failed (Error: unable to open or parse config file)\n" + str(e))
    exit(1)

# Load secret file, read JWT secret token. Exit (non-zero) if unable to open
try:
    sf = open(cg["secret_file"])
    secret = sf.read().strip()
    sf.close()
except Exception as e:
    print("dotplot_dots.py failed (Error: unable to open secret file) %s" % str(e))
    stderr.write("dotplot_dots.py failed (Error: unable to open secret file)\n" + str(e))
    exit(1)


# Load input into variable "get_values"
get_values = open(args.input, 'r')

# ---------------------------------------------------------------------------------------------------------------------
# Define Global Variables
# ---------------------------------------------------------------------------------------------------------------------

data = {}  # holds final data structure
# Final "data" structure
# data = { syntenic_points: { "sp1_genomeID" : {
#                                               "sp2_genomeID" : {
#                                                                 "sp1_chr" : {
#                                                                              "sp2_chr : [
#                                                                                          sp1_start,
#                                                                                          sp1_stop,
#                                                                                          sp2_start,
#                                                                                          sp2_stop,
#                                                                                          {"Kn": *,
#                                                                                           "Ks": *,
#                                                                                           "Sp1": {"name": *,
#                                                                                                   "chr" : *,
#                                                                                                   "start" : *,
#                                                                                                   "stop" : *,
#                                                                                                   "strand" : *,
#                                                                                                   "type" : *
#                                                                                                   "gene_count" : *,
#                                                                                                   "db_feature_id" : *,
#                                                                                                   "percent_id" : *
#                                                                                                  },
#                                                                                           "Sp2": {"name": *,
#                                                                                                   "chr" : *,
#                                                                                                   "start" : *,
#                                                                                                   "stop" : *,
#                                                                                                   "strand" : *,
#                                                                                                   "type" : *
#                                                                                                   "gene_count" : *,
#                                                                                                   "db_feature_id" : *,
#                                                                                                   "percent_id" : *
#                                                                                                  }
#                                                                                         ],
#                                                                                         ...
#                                                                             },
#                                                                             ...
#                                                                },
#                                                                ...
#                                              }
#                           }
#         genomes: { "sp1_genomeID" : {"chromosomes" : [{"name": <name>, "length": <length>}, ...],
#                                      <additional genome data>
#                                     },
#                    "sp2_genomeID" : {"chromosomes" : [{"name": <name>, "length": <length>}, ...],
#                                      <additional genome data>
#                                     }
#                  }
#        }
sp1_id = []  # holds all sp1 id's extracted from dataset
sp2_id = []  # holds all sp2 id's extracted from dataset
sp1_chromosomes = []  # holds all sp1 chromosomes extracted from dataset
sp2_chromosomes = []  # holds all sp2 chromosomes extracted from dataset
hits = []  # holds syntenic point data before assignment into "data"

# ---------------------------------------------------------------------------------------------------------------------
# Extract Data from Input
# ---------------------------------------------------------------------------------------------------------------------

for line in get_values:
    # Skip comment lines.
    if line[0] == '#':
        pass
    # Process non-comment lines.
    else:
        tmpl = {}  # temporary template for holding syntenic point data
        # Split contents of line into array "line_contents"
        line_contents = line.split('\t')

        # Assign template values from line contents.
        # Kn/Ks
        tmpl["Ks"] = line_contents[0]
        tmpl["Kn"] = line_contents[1]
        sp1_info = line_contents[3].split('||')
        sp2_info = line_contents[7].split('||')

        # Species 1: Genome ID, Chromosome, Start, Stop
        id_chr_1 = line_contents[2].lstrip('a').partition("_")
        tmpl["sp1_id"] = id_chr_1[0]
        if id_chr_1[0] not in sp1_id:
            sp1_id.append(id_chr_1[0])
        tmpl["sp1_ch"] = id_chr_1[2]
        if id_chr_1[2] not in sp1_chromosomes:
            sp1_chromosomes.append(id_chr_1[2])
        tmpl["sp1_start"] = line_contents[4]
        tmpl["sp1_stop"] = line_contents[5]
        tmpl["sp1_info"] = {"chr": sp1_info[0],
                            "start": sp1_info[1],
                            "stop": sp1_info[2],
                            "name": sp1_info[3],
                            "strand": sp1_info[4],
                            "type": sp1_info[5],
                            "db_feature_id": sp1_info[6],
                            "gene_count": sp1_info[7],
                            "percent_id": sp1_info[8]}

        # Species 2: Genome ID, Chromosome, Start, Stop
        id_chr_2 = line_contents[6].lstrip('b').partition("_")
        tmpl["sp2_id"] = id_chr_2[0]
        if id_chr_2[0] not in sp2_id:
            sp2_id.append(id_chr_2[0])
        tmpl["sp2_ch"] = id_chr_2[2]
        if id_chr_2[2] not in sp2_chromosomes:
            sp2_chromosomes.append(id_chr_2[2])
        tmpl["sp2_start"] = line_contents[8]
        tmpl["sp2_stop"] = line_contents[9]
        tmpl["sp2_info"] = {"chr": sp2_info[0],
                            "start": sp2_info[1],
                            "stop": sp2_info[2],
                            "name": sp2_info[3],
                            "strand": sp2_info[4],
                            "type": sp2_info[5],
                            "db_feature_id": sp2_info[6],
                            "gene_count": sp2_info[7],
                            "percent_id": sp2_info[8]}

        # Append hit (template) to "hits" list
        hits.append(tmpl)

# ---------------------------------------------------------------------------------------------------------------------
# Build 'data' structure outline & populate 'syntenic_points' component.
# ---------------------------------------------------------------------------------------------------------------------

if len(hits) < 1:  # If no hits are found, use command line gid's & create data object.
    sp1_id = args.gid1
    sp2_id = args.gid2
    data["syntenic_points"] = {sp1_id: {sp2_id: {}}}
    data["genomes"] = {}
else:  # If hits are found, determine IDs from input file & create data object.
    if len(sp1_id) != 1 or len(sp2_id) != 1:
        print("dotplot_dots.py failed (Error: wrong number of genome IDs in input file.)")
        stderr.write("dotplot_dots.py failed (Error: wrong number of genome IDs in input file.)\n")
        exit(1)
    else:
        sp1_id = sp1_id[0]
        sp2_id = sp2_id[0]

    # Build basic "data" structure outline.
    data["syntenic_points"] = {sp1_id: {sp2_id: {}}}
    for ch1 in sp1_chromosomes:
        data["syntenic_points"][sp1_id][sp2_id][ch1] = {}
        for ch2 in sp2_chromosomes:
            data["syntenic_points"][sp1_id][sp2_id][ch1][ch2] = []
    data["genomes"] = {}

    # Populate data["syntenic_points"] with hits.
    for hit in hits:
        # Entry Structure: [ch1_start, ch1_stop, ch2_start, ch2_stop, [kn, ks]]
        entry = [hit["sp1_start"],
                 hit["sp1_stop"],
                 hit["sp2_start"],
                 hit["sp2_stop"],
                 {"Kn": hit["Kn"],
                  "Ks": hit["Ks"],
                  str(sp1_id): hit["sp1_info"],
                  str(sp2_id): hit["sp2_info"]
                  }
                 ]
        data["syntenic_points"][hit["sp1_id"]][hit["sp2_id"]][hit["sp1_ch"]][hit["sp2_ch"]].append(entry)

    # Remove any empty sets.
    # Remove species 2 chromosomes with no matches.
    for chr_sp1 in data["syntenic_points"][sp1_id][sp2_id]:
        data["syntenic_points"][sp1_id][sp2_id][chr_sp1] = {chr_sp2: hits for chr_sp2, hits in data["syntenic_points"]\
            [sp1_id][sp2_id][chr_sp1].iteritems() if len(data["syntenic_points"][sp1_id][sp2_id][chr_sp1][chr_sp2]) > 0}
    # Remove species 1 chromosomes with no matches.
    data["syntenic_points"][sp1_id][sp2_id] = {sp1_ch: sp2_ch for sp1_ch, sp2_ch in data["syntenic_points"][sp1_id]\
        [sp2_id].iteritems() if len(data["syntenic_points"][sp1_id][sp2_id][sp1_ch]) > 0}

# ---------------------------------------------------------------------------------------------------------------------
# Populate 'genomes' component
# ---------------------------------------------------------------------------------------------------------------------

# Get genome information, validating if user is logged in.
if user == '':
    stderr.write("Requesting genome info: " + api_base + sp1_id + "\n")
    sp1_genome_info = get(api_base + sp1_id).json()
    stderr.write("Requesting genome info: " + api_base + sp2_id + "\n")
    sp2_genome_info = get(api_base + sp2_id).json()
else:
    # Build token & header.
    token = encode({'sub': user,
                    'exp': datetime.utcnow()+timedelta(seconds=60),
                    'iat': datetime.utcnow()},
                    secret, algorithm='HS256')
    get_headers = {'x-coge-jwt': token}

    # Request Species 1 genome info
    sp1_get = get(api_base + sp1_id, headers=get_headers)
    if sp1_get.status_code != 200:
        print(api_base + sp1_id + " returned status code " + str(sp1_get.status_code))
        stderr.write(api_base + sp1_id + " returned status code " + str(sp1_get.status_code) + "\n")
        exit(1)
    else:
        sp1_genome_info = sp1_get.json()

    # Request Species 2 Genome Info
    sp2_get = get(api_base + sp2_id, headers=get_headers)
    if sp1_get.status_code != 200:
        print(api_base + sp2_id + " returned status code " + str(sp1_get.status_code))
        stderr.write(api_base + sp2_id + " returned status code " + str(sp1_get.status_code) + "\n")
        exit(1)
    else:
        sp2_genome_info = sp2_get.json()

# Check for genome fetch errors, exit without producing outputs if detected.
try:
    sp1err = sp1_genome_info["error"]
    print("dotplot_dots.py died on genome info fetch error!")
    print("Error (gid " + sp1_id + "): " + str(sp1err))
    stderr.write("dotplot_dots.py died on genome info fetch error!" + "\n")
    stderr.write("Error (gid " + sp1_id + "): " + str(sp1err) + "\n")
    try:
        sp2err = sp2_genome_info["error"]
        print("Error (gid " + sp2_id + "): " + str(sp2err))
        stderr.write("Error (gid " + sp2_id + "): " + str(sp2err) + "\n")
    except KeyError:
        pass
    exit(1)
except KeyError:
    pass
try:
    sp2err = sp2_genome_info["error"]
    print("dotplot_dots.py died on genome info fetch error!")
    print("Error (gid " + sp2_id + "): " + str(sp2err))
    stderr.write("dotplot_dots.py died on genome info fetch error!" + "\n")
    stderr.write("Error (gid " + sp2_id + "): " + str(sp2err) + "\n")
    exit(1)
except KeyError:
    pass

# Populate data["genomes"] with genome information.
data["genomes"][sp1_id] = sp1_genome_info
data["genomes"][sp2_id] = sp2_genome_info

# ---------------------------------------------------------------------------------------------------------------------
# Write output files.
# ---------------------------------------------------------------------------------------------------------------------

# Write output.
dump(data, open(args.output, 'wb'))

# Write log
log["status"] = "complete"
log["message"] = "%s syntenic pairs identified" % str(len(hits))
dump(log, open(args.log, "wb"))
