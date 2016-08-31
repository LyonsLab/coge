#!/usr/bin/env python

# FractBias, Command Line Edition.
# Contributers: Blake Joyce (bjoyce3), Asher Baltzell (asherkhb)

# ------------------------------------------------------------------------------------------------------------------- #
# Record starting time.
# ------------------------------------------------------------------------------------------------------------------- #

from datetime import datetime, timedelta
t_start = datetime.now()

# ------------------------------------------------------------------------------------------------------------------- #
# Load dependencies.
# ------------------------------------------------------------------------------------------------------------------- #

print("Loading Dependencies")

from operator import itemgetter
import csv
from itertools import islice
from natsort import natsorted
import requests
from json import dump
import jwt
from sys import stderr
import argparse

t_cp1 = datetime.now()
print('--> Complete! (' + str(t_cp1-t_start) + ')')

# ------------------------------------------------------------------------------------------------------------------- #
# Parse arguments.
# ------------------------------------------------------------------------------------------------------------------- #

parser = argparse.ArgumentParser()
parser.add_argument("--align", help="File path to SynMap syntenic alignment data output (aligncoords.gcoords file)")
parser.add_argument("--gff", help="File path to target genome GFF file")
parser.add_argument("--target",
                    help="Target genome CoGe ID number. The target genome is the one with fewer subgenomes. Used to "
                         "determine which side of the SynMap output file the target genome data exists.")
parser.add_argument("--windowsize", help="Sets the size of the sliding window for analysis")
parser.add_argument("--query",
                    help="Query genome CoGe ID number. The query genome is the one with more subgenomes or higher "
                         "syntenic depth. Used to determine which side of the SynMap output file the target genome "
                         "data exists.")
parser.add_argument("--output",
                    help="File path to output directory. Specify SynMap directory: 'storage/coge/data/diags' "
                         "file path.")
parser.add_argument("--allgenes",
                    help="Determines whether all genes in target genome gff or only conserved genes should be used in "
                         "analysis.")
parser.add_argument("--numtargetchr", help="Restricts the maximum number of target chromosomes to use in the analysis.",
                    type=int)
parser.add_argument("--numquerychr", help="Restricts the maximum number of target chromosomes to use in the analysis.",
                    type=int)
parser.add_argument("--remove_random_unknown",
                    help="Restricts the maximum number of target chromosomes to use in the analysis.", type=str)
parser.add_argument("--syndepth", help="Passes synetic depth setting from SynMap for saving a unique figure.", type=str)
parser.add_argument("--apiurl",
                    help="URL to CoGe genomes API endpoint (i.e. https://genomevolution.org/coge/api/v1/genomes)",
                    type=str)
parser.add_argument("--user", help="User requesting job (CoGe username)", type=str)
args = parser.parse_args()

# ------------------------------------------------------------------------------------------------------------------- #
# Define methods and global variables.
# ------------------------------------------------------------------------------------------------------------------- #

# Variables for importing data and parsing.
synmap_import_file = args.align
gff_import_file = args.gff

# Variables for storing data.
# for SynMap data.
d = {}
# Will contain the array of syntenic genome1_chrs, genome1_genes, genome2_chrs, and genome2_genes.
# d{target_chr
#     {target_gene
#         {query_chr
#             {query gene}
#   }}}
# for GFF data
gff_genes_target = {}  # For organization of genes on chromosomes within target genome according to start bp
gff_genes_query = {}  # For organization of genes on chromosomes withing query genome according to start bp

# Define genus/species.
genus_species = ''
with open(gff_import_file) as gff_file:
    for line in gff_file:
        if line[0:15] == '##Organism name':
            genus_species = line[17:-1]
            species_name = genus_species.replace(' ', '_')
            species_name_filter = species_name.translate(None, '(){}[]')


# Define sliding window function.
def window(seq, n):
    """
    Returns a sliding window (of width n) over data from the iterable
    s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...
    http://stackoverflow.com/questions/6822725/rolling-or-sliding-window-iterator-in-python
    """
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result

# ------------------------------------------------------------------------------------------------------------------- #
# Fetch genome info from CoGe API & store in various objects.
# ------------------------------------------------------------------------------------------------------------------- #

print("Fetching Genome Info")
t_cp2 = datetime.now()

# Build & submit requests.
if args.user == '':  # Anonymous user case.
    query_api = requests.get(args.apiurl.rstrip('/') + '/' + str(args.query))
    target_api = requests.get(args.apiurl.rstrip('/') + '/' + str(args.target))
else:  # Registered user case, attempt authentication.
    token = jwt.encode({'sub': args.user,
                        'exp': datetime.utcnow() + timedelta(seconds=60),
                        'iat': datetime.utcnow()},
                       'fracbias', algorithm='HS256')
    get_headers = {'x-coge-jwt': token}
    query_api = requests.get(args.apiurl.rstrip('/') + '/' + str(args.query), headers=get_headers)
    target_api = requests.get(args.apiurl.rstrip('/') + '/' + str(args.target), headers=get_headers)

# Store results (JSON) as python objects.
query_api = query_api.json()
target_api = target_api.json()

# Check results for errors.
try:
    err1 = query_api["error"]
    stderr.write("fractionation_bias.py died on genome info fetch error!" + "\n")
    stderr.write("Error (gid " + str(args.query) + "): " + str(err1) + "\n")
    try:
        err2 = target_api["error"]
        stderr.write("Error (gid " + str(args.target) + "): " + str(err2) + "\n")
    except KeyError:
        pass
    exit()
except KeyError:
    pass
try:
    err3 = target_api["error"]
    stderr.write("dotplot_dots.py died on genome info fetch error!" + "\n")
    stderr.write("Error (gid " + str(args.target) + "): " + str(err3) + "\n")
    exit()
except KeyError:
    pass

# Initializes chromosome lists.
target_api_chrs = []
# target_api_chrs_sorted = []
# target_api_chrs_sorted_filtered = []
# target_api_chrs_sorted_name = []
target_api_chrs_final = []
query_api_chrs = []
# query_api_chrs_sorted = []
# query_api_chrs_sorted_filtered = []
# query_api_chrs_sorted_name = []
query_api_chrs_final = []

# Parse query genome chromosomes (restricted by length largest -> smallest, and sorted according to name).
for c in query_api['chromosomes']:
    if args.remove_random_unknown == 'True':  # AKB converted to str 8/30/16
        if ("Random" in c['name']) or ("random" in c['name']) or ("Unknown" in c['name']) or ("unknown" in c['name']):
            continue
        else:
            query_api_chrs.append((c['name'], c['length']))
    if args.remove_random_unknown == 'False':  # AKB converted to str 8/30/16
        query_api_chrs.append((c['name'], c['length']))
query_api_chrs_sorted = natsorted(query_api_chrs, key=lambda ch: ch[1], reverse=True)
query_api_chrs_sorted_filtered = query_api_chrs_sorted[0:int(args.numquerychr)]
query_api_chrs_sorted_name = natsorted(query_api_chrs_sorted_filtered, key=lambda ch: ch[0])
for c in query_api_chrs_sorted_name:
    query_api_chrs_final.append(c[0])

# Parse target genome chromosomes (restricted by length largest ->smallest, and sorted according to name).
for c in target_api['chromosomes']:
    if args.remove_random_unknown == 'True':  # AKB converted to str 8/30/16
        # print c['name']  # AKB removed 8/30/16 to prevent print bloat.
        if ("Random" in c['name']) or ("random" in c['name']) or ("Unknown" in c['name']) or \
                ("unknown" in c['name']):
            continue
        else:
            target_api_chrs.append((c['name'], c['length']))
    if args.remove_random_unknown == 'False':  # AKB converted to str 8/30/16
        target_api_chrs.append((c['name'], c['length']))
target_api_chrs_sorted = natsorted(target_api_chrs, key=lambda ch: ch[1], reverse=True)
target_api_chrs_sorted_filtered = target_api_chrs_sorted[0:int(args.numtargetchr)]
target_api_chrs_sorted_name = natsorted(target_api_chrs_sorted_filtered, key=lambda ch: ch[0])
for c in target_api_chrs_sorted_name:
    target_api_chrs_final.append(c[0])

# Checkpoint and print genome fetch time.
t_cp3 = datetime.now()
print('--> Complete! (' + str(t_cp3-t_cp2) + ')')

# ------------------------------------------------------------------------------------------------------------------- #
# Read & parse SynMap and GFF files.
# ------------------------------------------------------------------------------------------------------------------- #

print "Making Synmap Data Structure"
with open(synmap_import_file, 'r') as f:  # open SynMap file containing syntenic genes
    synmap_rowcount = 0
    cols = []  # list for parsing columns from SynMap data
    decimal_strip_check_target_gene = 0
    for line in f:  # for loop to parse columns
        new_line = line.replace('||', '\t')  # converts || into tabs for universal delimination
        if line[0] != '#' and line[0] != '\n':  # sorts out columns containing syntenic block information/headings
            cols = new_line.split('\t', )  # splits all syntenic gene pair lines into parsed columns in a list
            synmap_rowcount += 1
            global target_chr
            global target_gene
            global query_chr
            global query_gene
            if synmap_rowcount == 1:
                # clean ID subgenome A from the column on left of data
                ida = cols[0]
                ida = ida[1:cols[0].index('_')]
            # Determines which are target and query genes
            if args.target == ida:
                target_chr = cols[1]
                if any(target_chr in s for s in target_api_chrs_final):
                    target_gene = str(cols[7])  # .rsplit(".", 1)[0]  #puts all genome1_genes with synteny into a list
                    # decimal_strip_check_target_gene = target_gene.find('.')
                    # if not decimal_strip_check_target_gene == (-1):
                    # target_gene = target_gene[:(decimal_strip_check_target_gene)]
                else:
                    continue
                query_chr = str(cols[13])  # puts all genome2_chrs with synteny to genes in genome1 into a list
                if any(query_chr in s for s in query_api_chrs_final):
                    # Put all genome2_genes with synteny to genes in a genome1 into a list
                    query_gene = str(cols[19])  # .rsplit(".", 1)[0]
                else:
                    continue
            else:
                target_chr = cols[13]
                if any(target_chr in s for s in target_api_chrs_final):
                    target_gene = str(cols[19])
                    # decimal_strip_check_target_gene = target_gene.find('.')
                    # if not decimal_strip_check_target_gene == (-1):
                    #    target_gene = target_gene[:(decimal_strip_check_target_gene)]
                    # puts all genome1_genes with synteny into a list
                else:
                    continue
                query_chr = str(cols[1])  # puts all genome2_chrs with synteny to genes in genome1 into a list
                if any(query_chr in s for s in query_api_chrs_final):
                    query_gene = cols[
                        7]  # .rsplit(".", 1)[0]  #puts all genome2_genes with synteny to genes in a genome1 into a list
                else:
                    continue
            # if not target_chr in d:  # AKB changed to 'not in' as per PEP 8
            if target_chr not in d:
                d[target_chr] = {}  # initializes the nested dictionary-primary level at genome1_chromosome
            # if not target_gene in d[target_chr]:  # AKB changed to 'not in' as per PEP 8
            if target_gene not in d[target_chr]:
                d[target_chr][target_gene] = {}  # initializes first nesting in dictionary-second level at genome1_genes
            # if not query_chr in d[target_chr][target_gene]:  # AKB changed to 'not in' as per PEP 8
            if query_chr not in d[target_chr][target_gene]:
                d[target_chr][target_gene][
                    query_chr] = query_gene  # initializes nested dictionary-third level at genome2_chr

# Checkpoint & print data structure construction time.
t_cp4 = datetime.now()
print('--> Complete! (' + str(t_cp4-t_cp3) + ')')

# ------------------------------------------------------------------------------------------------------------------- #
# Read target genome GFF & build gff_genes.
# ------------------------------------------------------------------------------------------------------------------- #

print("Reading and Processing GFF File")

# Create hashable data structure to look up gene/CDS name. CDS_name_dict: CDS_key_name: 1 if TRUE.
CDS_name_dict = {}  # initializes the CDS name dictionary for looking up if CDS already have been identified
# gff_genes_target = {}  # AKB removed 8/30/16 because redefined from above.

# Open & parse GFF file.
with open(gff_import_file, 'r') as g:  # opens gff file
    gffcols = []  # list of parsed gff columns
    # chr = []  # initialize list of chromosomes present in genome1 gff file  # AKB removed 8/20/16 seems unused.
    for line in g:
        new_line = line.replace(';', '\t')  # makes subdelims universal in gff file from CoGe
        new_line = new_line.replace('Name=', '')  # strips Name= off gene_name in gff file from CoGe
        if new_line[0] != '#' and new_line[0] != '\n':  # selects only lines with CDS information
            gffcols = new_line.split('\t', )  # parses all columns
            # select only 'CDS' lines for consideration
            if any(gffcols[0] in s for s in target_api_chrs_final) and (gffcols[2] == 'CDS'):
                c = gffcols[0]  # adds genome1_chrs to list
                # if not chr in gff_genes_target:  # AKB changed to 'not in' as per PEP 8
                if c not in gff_genes_target:
                    gff_genes_target[c] = []  # initializes chr list in dictionary if chr does not exist yet
                gene_name = str(gffcols[-1])
                gene_name = gene_name.rstrip("\n")
                gene_name1 = gene_name[9:]  # adds targetgenome_genes to list and removes the "ID=" added by CoGe
                start = int(gffcols[3])  # adds targetgenome_gene start bp to list for ordering as integer
                stop = int(gffcols[4])  # adds targetgenome_gene stop bp to list ?for ordering? as integer
                try:
                    CDS_name_dict[str(gene_name1)] == 1
                    continue
                except KeyError:
                    gff_genes_target[c].append(dict(gene_name=gene_name1, start=start, stop=stop))
                    CDS_name_dict[str(gene_name1)] = 1

# Sort GFF genes within target chromosomes by start position.
# gff_sort_all = {} AKB removed 8/30/16 because redefined right below
gff_sort_all = natsorted(gff_genes_target)
for c in gff_genes_target:
    # Creates dictionary for searching genes against. TODO: CONSIDER sorting on midpoint of genes
    gff_genes_sorted = sorted(gff_genes_target[c], key=itemgetter('start'))
    gff_genes_target[c] = gff_genes_sorted
    # TODO: CONSIDER WRITING A CHECK PROGRAM TO RETURN TRUE IF ALL VALUES ARE SORTED OR FALSE

# # Write SynMap dictionary and sorted GFF gene list to files. AKB removed 8/30/16 because unnecessary.
# print "Writing Out GFF and SynMap Data Structure to CSV"
# gff_sort_output_file = (args.output + "/gff_sort.txt")
# with open(str(gff_sort_output_file), 'w') as h:
#     h.write(str(gff_genes_target))
# synmap_dictionary_output_file = (args.output + "/synmap_data_structure.txt")
# with open(synmap_dictionary_output_file, 'w+') as i:
#     i.write(str(d))

# Checkpoint & print GFF processing time.
t_cp5 = datetime.now()
print('--> Complete! (' + str(t_cp5-t_cp4) + ')')

# ------------------------------------------------------------------------------------------------------------------- #
# Determine syntenic gene pairs present and output raw data CSV file from parsed data.
# ------------------------------------------------------------------------------------------------------------------- #

print "Calling Syntenic Pairs and Processing Raw SynMap Data"

windanalysis_input_dict = {}
fract_bias_raw_output_file = args.output + "/fractbias_output-" + "-TarID" + str(args.target) + "-TarChrNum" + \
                             str(args.numtargetchr) + "-SynDep" + str(args.syndepth) + "-QueryID" + str(args.query) + \
                             "-QueryChrNum" + str(args.numquerychr) + "-AllGene" + str(args.allgenes) + "-RmRnd" + \
                             str(args.remove_random_unknown) + "-WindSize" + str(args.windowsize) + ".csv"

with open(str(fract_bias_raw_output_file), 'w') as csvfile:
    csvfile.write('# Target Species: ' + target_api['organism']['name'] + ' (CoGe ID: ' + str(args.target) + ')\n')
    csvfile.write('# Query Species: ' + query_api['organism']['name'] + ' (CoGe ID: ' + str(args.query) + ')\n')
    csvfile.write('# Generated ' + str(datetime.today()) + '\n##\n')
    headers = ['Target Chromosome', 'Target Gene Name', 'Gene Order on Target Chromosome']
    headers.extend(query_api_chrs_final)
    headers.extend(query_api_chrs_final)
    writer = csv.writer(csvfile, dialect='excel', delimiter=',', lineterminator='\n')
    writer.writerow(headers)

    if args.allgenes == 'True':
        for tchr in gff_genes_target:
            col0 = chr  # writes Pineapple chr number
            count = 0
            for diction in gff_genes_target[tchr]:
                gene = diction['gene_name']
                col1 = gene  # writes gene name
                count += 1
                col2 = count  # writes pineapple gene number on pineapple chr
                # Find the query chr genes and output to columns: first gff info (col0-3), query chr (col 4-n),
                # query chr-gene (col n+1-m)
                syntenic_query_gene_presence_data = []
                syntenic_query_gene_name = []
                for qchr in query_api_chrs_final:
                    # if not tchr in windanalysis_input_dict:  # AKB changed to 'not in' as per PEP 8
                    if tchr not in windanalysis_input_dict:
                        windanalysis_input_dict[
                            tchr] = {}  # initializes the nested dictionary-primary level at genome1_chromosome
                    # if not qchr in windanalysis_input_dict[tchr]:  # AKB changed to 'not in' as per PEP 8
                    if qchr not in windanalysis_input_dict[tchr]:
                        windanalysis_input_dict[tchr][
                            qchr] = []  # initializes first nesting in dictionary-second level at genome1_genes
                    try:
                        syn_gene = d[tchr][gene][qchr]
                        syntenic_query_gene_presence_data.append(True)
                        syntenic_query_gene_name.append(syn_gene)
                        windanalysis_input_dict[tchr][qchr].append(True)
                    except KeyError:
                        syntenic_query_gene_presence_data.append(False)
                        syntenic_query_gene_name.append("")
                        windanalysis_input_dict[tchr][qchr].append(False)
                rows = [tchr, col1, col2]
                rows.extend(syntenic_query_gene_presence_data)
                rows.extend(syntenic_query_gene_name)
                writer.writerow(rows)

    elif args.allgenes == 'False':
        for tchr in gff_genes_target:
            col0 = chr  # writes Pineapple chr number
            count = 0
            for diction in gff_genes_target[tchr]:
                gene = diction['gene_name']
                col1 = gene  # writes target gene name
                count += 1
                col2 = count  # writes target gene number on target chr
                # Find the query chr genes and output to columns: first gff info (col0-3), query chr (col 4-n),
                # query chr-gene (col n+1-m)
                syntenic_query_gene_presence_data = []
                syntenic_query_gene_name = []
                for qchr in query_api_chrs_final:
                    # if not tchr in windanalysis_input_dict: # AKB changed to 'not in' as per PEP 8
                    if tchr not in windanalysis_input_dict:
                        windanalysis_input_dict[
                            tchr] = {}  # initializes the nested dictionary-primary level at genome1_chromosome
                    # if not qchr in windanalysis_input_dict[tchr]:  # AKB changed to 'not in' as per PEP 8
                    if qchr not in windanalysis_input_dict[tchr]:
                        windanalysis_input_dict[tchr][
                            qchr] = []  # initializes first nesting in dictionary-second level at genome1_genes

                    try:
                        syn_gene = d[tchr][gene][qchr]
                        syntenic_query_gene_presence_data.append(True)
                    except KeyError:
                        syntenic_query_gene_presence_data.append(False)

                if sum(syntenic_query_gene_presence_data) >= 1:
                    for qchr in query_api_chrs_final:
                        try:
                            syn_gene = d[tchr][gene][qchr]
                            syntenic_query_gene_name.append(syn_gene)
                            windanalysis_input_dict[tchr][qchr].append(True)
                        except KeyError:
                            syntenic_query_gene_name.append("")
                            windanalysis_input_dict[tchr][qchr].append(False)
                    rows = [tchr, col1, col2]
                    rows.extend(syntenic_query_gene_presence_data)
                    rows.extend(syntenic_query_gene_name)
                    writer.writerow(rows)

                elif sum(syntenic_query_gene_presence_data) < 1:
                    continue
                else:
                    print 'Target gene not classified'
                    continue

    else:
        print "Genes to be used (all genes or at least one syntenic gene) are not defined"

# Checkpoint & print synteny processing time.
t_cp6 = datetime.now()
print('--> Complete! (' + str(t_cp6-t_cp5) + ')')

# ------------------------------------------------------------------------------------------------------------------- #
# Sliding window analysis.
# ------------------------------------------------------------------------------------------------------------------- #
# Looks through each chromosome in target genome (tchr), and goes through each gene gene on that target chromosome
# (determined by the GFF uploaded above according to bp position on that target chromosome) to compare to the query
# genome (query chromsomes).
# i.e. For each chromosome in target genome, reads the genes on that chromosome and compare to all query subgenomes.

print "Running Sliding Window Analysis"

# Initialize variables.
data_output0 = []
data_output1 = []
data_output2 = []
data_output3 = []
output_dict = {}
window_size = 100
window_size_int = 0

# Process windows 100genes/sliding window and output to nested dictionary data structure:
# output_dict{target chr:}{query chr}{window count:retention%}
for tchr in windanalysis_input_dict:
    tchr_counter = tchr
    for qchr in windanalysis_input_dict[tchr]:
        counter = 0
        qchr_counter = qchr
        # if not tchr in output_dict:  # AKB changed to 'not in' as per PEP 8
        if tchr not in output_dict:
            output_dict[tchr] = {}  # initializes the nested dictionary-primary level at genome1_chromosome
        # if not qchr in output_dict[tchr]:  # AKB changed to 'not in' as per PEP 8
        if qchr not in output_dict[tchr]:
            output_dict[tchr][qchr] = {}  # initializes first nesting in dictionary-second level at genome1_genes
        # Match two arguments to int for logical comparison below, geco sees args.input as 'string'
        window_size_int = int(args.windowsize)
        try:
            if (int(len(windanalysis_input_dict[tchr][qchr]))) >= window_size_int:  # <--- changed from args.windowsize
                for each in window(windanalysis_input_dict[tchr][qchr],
                                   window_size_int):  # <--- changed from args.windowsize
                    counter += 1
                    data_output2 = float(sum(each)) / float(window_size)
                    output_dict[tchr][qchr][counter] = round(data_output2 * 100.)
        except KeyError:
            continue

# # Print some metrics.
# print "Length of d " + str(len(d))
# print "Length of gff_sort_all " + str(len(gff_sort_all))
# print "Length of gff_genes_target " + str(len(gff_genes_target))
# print "Length of windanalysis_input_dict " + str(len(windanalysis_input_dict))
# print "Length of output_dict " + str(len(output_dict))

# Checkpoint & print sliding window analysis time.
t_cp7 = datetime.now()
print('--> Complete! (' + str(t_cp7-t_cp6) + ')')

# ------------------------------------------------------------------------------------------------------------------- #
# Build & write final output files.
# ------------------------------------------------------------------------------------------------------------------- #
print("Building & Writing Final Outputs")

# Output Data for Interactive Plotting (AKB added 7/8/16)
dump_out = args.output + "/fractbias_figure-" + "-TarID" + str(args.target) + "-TarChrNum" + \
           str(args.numtargetchr) + "-SynDep" + str(args.syndepth) + "-QueryID" + str(args.query) + \
           "-QueryChrNum" + str(args.numquerychr) + "-AllGene" + str(args.allgenes) + "-RmRnd" + \
           str(args.remove_random_unknown) + "-WindSize" + str(args.windowsize) + ".json"
dump_data = {}

# Output Data for Downloading (AKB added 8/30/16)
download_out = args.output + "/fractbias_results-" + "-TarID" + str(args.target) + "-TarChrNum" + \
               str(args.numtargetchr) + "-SynDep" + str(args.syndepth) + "-QueryID" + str(args.query) + \
               "-QueryChrNum" + str(args.numquerychr) + "-AllGene" + str(args.allgenes) + "-RmRnd" + \
               str(args.remove_random_unknown) + "-WindSize" + str(args.windowsize) + ".csv"
down_data = {'targets': [], 'x_pos': []}  # [[targets...], [x_pos...], [query1_ys], ... [queryN_ys]]
for q in query_api_chrs_final:
    down_data[q] = []

# Build dump_data and initial down_data objects
ordered_targets = natsorted(output_dict.keys())
for target in ordered_targets:
    dump_data[target] = {}
    ordered_queries = natsorted(output_dict[target].keys())
    length = []
    for query in ordered_queries:
        dump_data[target][query] = []
        ordered_positions = natsorted(output_dict[target][query].keys())
        for position in ordered_positions:
            dump_data[target][query].append(output_dict[target][query][position])
            down_data['targets'].append(target)
            down_data['x_pos'].append(position)
            down_data[query].append(output_dict[target][query][position])

# Write out plotting object (.json)
dump(dump_data, open(dump_out, 'wb'))

# Zip download data (redefine down_data to zipped version)
down_header = ['Target Chr', 'Sliding Window (Count)']
down_prezip = [down_data['targets'], down_data['x_pos']]
for q in query_api_chrs_final:
    down_header.append('Query Chr ' + q + ' (%)')
    down_prezip.append(down_data[q])
down_data = zip(*down_prezip)
down_prezip = None  # Free memory.

# Write out downloadable data (.csv)
with open(download_out, 'wb') as o:
    o.write('# Target Species: ' + target_api['organism']['name'] + ' (CoGe ID: ' + str(args.target) + ')\n')
    o.write('# Query Species: ' + query_api['organism']['name'] + ' (CoGe ID: ' + str(args.query) + ')\n')
    o.write('# Generated ' + str(datetime.today()) + '\n##\n')
    w = csv.writer(o)
    w.writerow(down_header)
    for row in down_data:
        w.writerow(row)

# Checkpoint & print outputting time.
t_cp8 = datetime.now()
print('--> Complete! (' + str(t_cp8-t_cp7) + ')')

# ------------------------------------------------------------------------------------------------------------------- #
# Print concluding message.
# -------------------------------------------------------------------------------------------------------------------

print('FractBias Complete. Total time elapsed ' + str(t_cp8-t_start))

# ------------------------------------------------------------------------------------------------------------------- #
# Deprecated MatPlotLib Plotting
# ------------------------------------------------------------------------------------------------------------------- #

# # Imports for Plotting
# # Sets matplotlib to not output figures in a window (will cause an error if removed)
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
# import seaborn as sns
#
# # Check output data before graphing
# countofchrgraph = 0
# listofchrgraph = []
#
# for tchr in output_dict:
#     sumofchr = 0
#     for qchr in output_dict[tchr]:
#         try:
#             (max(output_dict[tchr][qchr].itervalues()))
#             sumofchr = sumofchr + (max(output_dict[tchr][qchr].itervalues()))
#         except ValueError:
#             continue
#     if sumofchr > 5:
#         countofchrgraph += 1
#         listofchrgraph.append(str(tchr))
#
# print listofchrgraph
# listofchrgraph = natsorted(listofchrgraph)
# print listofchrgraph
# print "Plotting FractBias Data"
#
#
# print "setting up figsize, cols, and gridspec"
# # define figure size, column layout, grid layout
# figsize = (15, (0.2 * len(query_api_chrs_final)) + 20)  # (15, (len(target_api_chrs_final)*2.4))
# cols = 2
# gs = gridspec.GridSpec(len(output_dict) // cols + 1, cols)
#
# print "setting up tableau palette"
#
# # These are the "Tableau 20" colors as RGB
# tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
#              (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
#              (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
#              (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
#              (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
#
# print "scaling palette"
# # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts
# for i in range(len(tableau20)):
#     r, g, b = tableau20[i]
#     tableau20[i] = (r / 255., g / 255., b / 255.)
#
# # Alternate color palettes Seaborn/Brewer
# print "setting palette"
# current_palette = sns.color_palette("Set2", 40)
#
# # tableau20
# print "Starting to plot figure"
#
# fig = plt.figure(figsize=figsize, frameon=False)
# print "set figure size and plt.figure"
# subplt_count = -1
# ax = []
# print "beginning subplots"
# for tchr in listofchrgraph:
#     subplt_count += 1
#     print "Plotting subplot: " + str(subplt_count)
#     count = 0
#     row = (subplt_count // cols)
#     col = subplt_count % cols
#     ax.append(fig.add_subplot(gs[row, col]))
#     for qchr in output_dict[tchr]:
#         count += 1
#         try:
#             if (max(output_dict[tchr][qchr].itervalues())) > 5:
#                 x = output_dict[tchr][qchr].keys()
#                 y = output_dict[tchr][qchr].values()
#                 # Sets up plotting conditions
#                 ax[-1].spines["top"].set_visible(False)
#                 ax[-1].spines["right"].set_visible(False)
#                 ax[-1].spines["left"].set_visible(True)
#                 ax[-1].spines["bottom"].set_visible(True)
#                 ax[-1].patch.set_visible(False)
#                 ax[-1].get_xaxis().tick_bottom()
#                 ax[-1].get_yaxis().tick_left()
#                 ax[-1].plot(x, y, color=current_palette[count], lw=2, label=str(qchr))
#                 ax[-1].set_title(label='Target Chromosome: ' + species_name_filter + " " + tchr, fontweight='bold',
#                                  fontsize=14, y=1.1, loc='left')
#                 ax[-1].set_xlabel('Window Iteration\n(Gene number)', fontsize=12, fontweight='bold')
#                 ax[-1].set_ylabel('% Retention\n(#Syntenic genes/window size)', fontsize=12, fontweight='bold')
#                 ax[-1].legend(bbox_to_anchor=(1.25, 1.1), loc=1, frameon=False, title="       Query\nChromosome",
#                               fontsize=10)
#             else:
#                 continue
#         except ValueError:
#             continue
#
# fig.subplots_adjust(wspace=0.45, hspace=0.6)
# plt.savefig(args.output + "/html/" + "fractbias_figure-" + "-TarID" + str(args.target) + "-TarChrNum" + str(
#     args.numtargetchr) + "-SynDep" + str(args.syndepth) + \
#             "-QueryID" + str(args.query) + "-QueryChrNum" + str(args.numquerychr) + "-AllGene" + str(
#     args.allgenes) + "-RmRnd" + str(args.remove_random_unknown) + "-WindSize" \
#             + str(args.windowsize) + ".png", transparent=True, bbox_inches='tight')
