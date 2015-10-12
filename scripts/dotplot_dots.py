__author__ = 'senorrift'
# Dependencies:

from json import dump
from os import path
from requests import get
from sys import argv, stderr

# api_base = "https://genomevolution.org/coge/api/v1/genomes/"
api_base = "https://geco.iplantc.org/asherkhb/coge/api/v1/genomes/"

# ---------------------------------------------------------------------------------------------------------------------
# Define Input/Outputs
# - 1st command line argument.
# - Synmap Output - must include ks values [this file will end with '.aligncoords.gcoords.ks').
# ---------------------------------------------------------------------------------------------------------------------

log = {}

# Assign input to variable "synmap_output", raise error and exit if not specified.
try:
    input_ksfile = argv[1]
    output_dir = path.dirname(input_ksfile).rstrip("/")
    #output_dir = argv[2].rstrip("/")
except IndexError:
    log["status"] = "failed"
    log["message"] = "Error: input/output specification"
    # TODO: Add Log Dump!
    stderr.write("dotplot_dots.py failed\n")
    stderr.write("Error: input/output specification")
    exit()

# Load input into variable "get_values"
get_values = open(input_ksfile, 'r')

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
#                                                                                          ["Kn", "Ks"]
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

        # Species 1: Genome ID, Chromosome, Start, Stop
        id_chr_1 = line_contents[2].lstrip('a').split("_")
        tmpl["sp1_id"] = id_chr_1[0]
        if id_chr_1[0] not in sp1_id:
            sp1_id.append(id_chr_1[0])
        tmpl["sp1_ch"] = id_chr_1[1]
        if id_chr_1[1] not in sp1_chromosomes:
            sp1_chromosomes.append(id_chr_1[1])
        tmpl["sp1_start"] = line_contents[4]
        tmpl["sp1_stop"] = line_contents[5]

        # Species 2: Genome ID, Chromosome, Start, Stop
        id_chr_2 = line_contents[6].lstrip('b').split("_")
        tmpl["sp2_id"] = id_chr_2[0]
        if id_chr_2[0] not in sp2_id:
            sp2_id.append(id_chr_2[0])
        tmpl["sp2_ch"] = id_chr_2[1]
        if id_chr_2[1] not in sp2_chromosomes:
            sp2_chromosomes.append(id_chr_2[1])
        tmpl["sp2_start"] = line_contents[8]
        tmpl["sp2_stop"] = line_contents[9]

        # Append hit (template) to "hits" list
        hits.append(tmpl)

# Reassign Species ID
if len(sp1_id) > 1 or len(sp2_id) > 1:
    log["status"] = "failed"
    log["message"] = "Error: Too many genome IDs"
    # TODO: Add Log Dump!
    stderr.write("dotplot_dots.py failed\n")
    stderr.write("Error: too many genome IDs")
    exit()
elif len(sp1_id) < 1 or len(sp2_id) < 1:
    log["status"] = "failed"
    log["message"] = "Error: Too few genome IDs"
    # TODO: Add Log Dump!
    stderr.write("dotplot_dots.py failed\n")
    stderr.write("Error: too few genome IDs")
    exit()
else:
    sp1_id = sp1_id[0]
    sp2_id = sp2_id[0]
    output_dir = "%s/%s/%s" % (output_dir, str(sp1_id), str(sp2_id))

# ---------------------------------------------------------------------------------------------------------------------
# Build & Dump Output JSON
# ---------------------------------------------------------------------------------------------------------------------

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
    entry = [hit["sp1_start"], hit["sp1_stop"], hit["sp2_start"], hit["sp2_stop"], {"Kn": hit["Kn"], "Ks": hit["Ks"]}]
    data["syntenic_points"][hit["sp1_id"]][hit["sp2_id"]][hit["sp1_ch"]][hit["sp2_ch"]].append(entry)

# Remove any empty sets.
# Remove species 2 chromosomes with no matches.
for chr_sp1 in data["syntenic_points"][sp1_id][sp2_id]:
    data["syntenic_points"][sp1_id][sp2_id][chr_sp1] = {chr_sp2: hits for chr_sp2, hits in data["syntenic_points"]\
        [sp1_id][sp2_id][chr_sp1].iteritems() if len(data["syntenic_points"][sp1_id][sp2_id][chr_sp1][chr_sp2]) > 0}
# Remove species 1 chromosomes with no matches.
data["syntenic_points"][sp1_id][sp2_id] = {sp1_ch: sp2_ch for sp1_ch, sp2_ch in data["syntenic_points"][sp1_id]\
    [sp2_id].iteritems() if len(data["syntenic_points"][sp1_id][sp2_id][sp1_ch]) > 0}

# Populate data["genomes"] with genome information.
sp1_genome_info = get(api_base + sp1_id).json()
sp2_genome_info = get(api_base + sp2_id).json()
data["genomes"][sp1_id] = sp1_genome_info
data["genomes"][sp2_id] = sp2_genome_info


# Dump "data" to JSON.
output_filename = "%s/%s_%s_synteny.json" % (output_dir, sp1_id, sp2_id)
dump(data, open(output_filename, 'wb'))

# Print concluding message.
log["status"] = "complete"
log["message"] = "%s syntenic pairs identified" % str(len(hits))
log_out = "%s/dotplot_dots_log.json" % output_dir
dump(log, open(log_out, "wb"))

# print "%s syntenic pairs identified!\ndotplot_dots.py Complete" % str(len(hits))