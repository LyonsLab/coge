"""Geco version of FractBias SynMap Tool"""
#geco specific code highlighted by '<--------' comments

#Sets matplotlib to not output figures in a window (will cause an error if removed) <--------
import matplotlib #<--------
matplotlib.use('Agg')

#Specifically for using code on server     <--------
#Allows arguments in command line
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--align", help="File path to SynMap syntenic alignment data output (aligncoords.gcoords file)")
parser.add_argument("--gff", help="File path to target genome GFF file")
parser.add_argument("--target", help="Target genome CoGe ID number. The target genome is the one with fewer subgenomes. Used to determine which side of the SynMap output file the target genome data exists.")
parser.add_argument("--windowsize", help="Sets the size of the sliding window for analysis")
parser.add_argument("--query", help="Query genome CoGe ID number. The query genome is the one with more subgenomes or higher syntenic depth. Used to determine which side of the SynMap output file the target genome data exists.")
parser.add_argument("--output", help="File path to output directory. Specify SynMap directory: 'storage/coge/data/diags' file path.")
parser.add_argument("--allgenes", help="Determines whether all genes in target genome gff or only conserved genes should be used in analysis.")
parser.add_argument("--numtargetchr", help="Restricts the maximum number of target chromosomes to use in the analysis.", type=int)
parser.add_argument("--numquerychr", help="Restricts the maximum number of target chromosomes to use in the analysis.", type=int)
parser.add_argument("--remove_random_unknown", help="Restricts the maximum number of target chromosomes to use in the analysis.", type=bool)
parser.add_argument("--syndepth", help="Passes synetic depth setting from SynMap for saving a unique figure.", type=str)
args = parser.parse_args()

#Flag for printing debug statements



#For importing data and parsing data and benchmarking program
from operator import itemgetter
from datetime import datetime, timedelta

#Converting parsed data into raw parsed data output to csv
import csv
from itertools import islice, izip


#For analyzing raw parsed data
import collections, re
from collections import OrderedDict

    #had to uninstall python-dateutil and use old version dateutil 2.2 to avoid error
    #sudo pip uninstall python-dateutil
    #sudo pip install python-dateutil==2.2
import matplotlib.pyplot as plt
import numpy as np 
import matplotlib.gridspec as gridspec
import seaborn as sns

#had to install this using pip on local computer
from natsort import natsorted

#Library to pull query genome chromosome information from CoGe API
import requests

#PATH TO RUN IN DESKTOP
#Sorghum/Maize comparison
#python fractionation_bias_geco_from_ipython_v2.py --align /home/bjoyce3/sorghum_maize/SorghumMaize_gcoords_merge60_syntenicdepth40.txt --gff /home/bjoyce3/sorghum_maize/Sorghum_bicolor_6807_prodrun.gff --target 6807 --windowsize 100 --output /home/bjoyce3/sorghum_maize/ --allgenes True

#pina/rice comparison
#python fractionation_bias_geco.py --align /home/bjoyce3/pinarice/SynMapKsMerge120Sdepth10_Osativa-Acomosus2-1Acomv6.txt --gff /home/bjoyce3/pinarice/Ananas_comosus_pineapple_v6.gff --target 25734 --windowsize 100 --output /home/bjoyce3/pinarice


"""Methods and Global Variables"""

#For importing data and parsing <---------
synmap_import_file = args.align     #<-----------
gff_import_file = args.gff     #<-----------

# initialize dictionary to contain the array of syntenic genome1_chrs, genome1_genes, genome2_chrs, and genome2_genes
d = {}

#Parsed data and raw output to csv
gff_genes_target = {}  # initializes dictionary for organization of genes on chromosomes within target genome according to start bp
gff_genes_query = {}  #initializes dictionary for ogranization of genes on chromosomes withing query genome according to start bp


genus_species = ''
with open(gff_import_file) as gff_file:
    for line in gff_file:
        if line[0:15] == '##Organism name':
            genus_species = line[17:-1]
            species_name = genus_species.replace(' ','_')
            species_name_filter = species_name.translate(None, '(){}[]')


#Parsed data and raw output to csv <-----------
gff_sort_output_file = (args.output+"/gff_sort.txt")
synmap_dictionary_output_file = (args.output+"/synmap_data_structure.txt")
fract_bias_raw_output_file = (args.output+"/fractbias_output.csv")

#Analysis of parsed data
retention_calc_output_file = ("fractbias_sliding_window_output.csv") #<--------

#Deprecated method for determining target and query chromosomes
"""def chr_id(input_dict):
    for item in input_dict:
        if not item in target_lst:
            target_lst.append(item)
        for gene in input_dict[item]:
            for chr in input_dict[item][gene]:
                if not chr in query_lst:
                    query_lst.append(chr)"""

#http://stackoverflow.com/questions/6822725/rolling-or-sliding-window-iterator-in-python
def window(seq, n):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result
        
#DATAPaths
#synmap_import_file = '/Users/bjoyce3/Desktop/SynMapFractBiasInput/Pina-rice/SynMapKsMerge120Sdepth10_Osativa-Acomosus2-1Acomv6.txt'

#gff_import_file = '/Users/bjoyce3/Desktop/SynMapFractBiasInput/Sorghum-MaizeSynMap/Sorghum_bicolor_annos1-cds0-id_typename-nu1-upa1-add_chr0.gid6807.gff'


'''Importing Data and Making Data Structures'''
# 
# Reads SynMap and GFF files and parse data into data arrays. 
# 
# SynMap data is put into nested dictionary called 'd':
# d{target_chr
#     {target_gene
#         {query_chr
#             {query gene}
#   }}}
# 
# GFF data is put into a nested dictionary call 'gff_genes':
# 

"""Calls Genome Fetch CoGE API to get query genome and target chromosomes and orders them"""
t_api_before = datetime.now()

#retrieves api chromsome lists and length and moves json() object into Python dictionary
query_api = requests.get("https://genomevolution.org/coge/api/v1/genomes/" + str(args.query))
target_api = requests.get("https://genomevolution.org/coge/api/v1/genomes/" + str(args.target))
query_api = query_api.json()
target_api = target_api.json()

#initializes dictionaries to drop api data into
target_api_chrs = []
target_api_chrs_sorted = []
target_api_chrs_sorted_filtered = []
target_api_chrs_sorted_name = []
target_api_chrs_final = []
query_api_chrs = []
query_api_chrs_sorted = []
query_api_chrs_sorted_filtered = []
query_api_chrs_sorted_name = []
query_api_chrs_final = []

#query genome chromosomes parsed, restricted by length largest -> smallest, and sorted according to name
for chr in query_api['chromosomes']:
    if args.remove_random_unknown == True:
        if ("Random" in chr['name']) or ("random" in chr['name']) or ("Unknown" in chr['name']) or ("unknown" in chr['name']):
            continue
        else:
            query_api_chrs.append((chr['name'], chr['length']))
    if args.remove_random_unknown == False:
        query_api_chrs.append((chr['name'], chr['length']))
query_api_chrs_sorted = natsorted(query_api_chrs, key=lambda chr: chr[1], reverse=True)
query_api_chrs_sorted_filtered = query_api_chrs_sorted[0:int(args.numquerychr)]
query_api_chrs_sorted_name = natsorted(query_api_chrs_sorted_filtered, key=lambda chr: chr[0])
for chr in query_api_chrs_sorted_name:
    query_api_chrs_final.append(chr[0])

#target genome chromosomes parsed, restricted by length largest ->smallest, and sorted according to name
for chr in target_api['chromosomes']:
    if args.remove_random_unknown == True:
        print chr['name']
        if ("Random" in chr['name']) or ("random" in chr['name']) or ("Unknown" in chr['name']) or ("unknown" in chr['name']):
            continue
        else:
            target_api_chrs.append((chr['name'], chr['length']))
    if args.remove_random_unknown == False:
        target_api_chrs.append((chr['name'], chr['length']))
target_api_chrs_sorted = natsorted(target_api_chrs, key=lambda chr: chr[1], reverse=True)
target_api_chrs_sorted_filtered = target_api_chrs_sorted[0:int(args.numtargetchr)]
target_api_chrs_sorted_name = natsorted(target_api_chrs_sorted_filtered, key=lambda chr: chr[0])
for chr in target_api_chrs_sorted_name:
    target_api_chrs_final.append(chr[0])

"""print target_api_chrs 
print target_api_chrs_sorted 
print target_api_chrs_sorted_filtered 
print target_api_chrs_sorted_name
print target_api_chrs_final"""
t_api_after = datetime.now()


"""
Reads SynMap and GFF files and parse data into columns in array

"""
print "Making Synmap Data Structure"
t0 = datetime.now()

print target_api_chrs_sorted_name

with open(synmap_import_file, 'r') as f:  # open SynMap file containing syntenic genes
    synmap_rowcount = 0
    cols = []  # list for parsing columns from SynMap data
    decimal_strip_check_target_gene = 0
    for line in f:  # for loop to parse columns
        new_line = line.replace('||', '\t')  #converts || into tabs for universal delimination
        if line[0] != '#' and line[0] != '\n':  #sorts out columns containing syntenic block information/headings
            cols = new_line.split('\t', )  #splits all syntenic gene pair lines into parsed columns in a list
            synmap_rowcount += 1
            global target_chr
            global target_gene
            global query_chr
            global query_gene
            if synmap_rowcount == 1:            
                #clean ID subgenome A from the column on left of data
                ida = cols[0]
                ida = ida[1:cols[0].index('_')]
            #Determines which are target and query genes                
            if args.target == ida:
                target_chr = cols[1]
                if any(target_chr in s for s in target_api_chrs_final):
                    target_gene = str(cols[7]) #.rsplit(".", 1)[0]  #puts all genome1_genes with synteny into a list
                    #decimal_strip_check_target_gene = target_gene.find('.')
                    #if not decimal_strip_check_target_gene == (-1):
                        #target_gene = target_gene[:(decimal_strip_check_target_gene)]
                else:
                    continue
                query_chr = str(cols[13])  #puts all genome2_chrs with synteny to genes in genome1 into a list
                if any(query_chr in s for s in query_api_chrs_final):
                    query_gene = str(cols[19])  #.rsplit(".", 1)[0]  #puts all genome2_genes with synteny to genes in a genome1 into a list
                else:
                    continue
            else:
                target_chr = cols[13]
                if any(target_chr in s for s in target_api_chrs_final):
                    target_gene = str(cols[19])
                    #decimal_strip_check_target_gene = target_gene.find('.')
                    #if not decimal_strip_check_target_gene == (-1):
                    #    target_gene = target_gene[:(decimal_strip_check_target_gene)]
                    #puts all genome1_genes with synteny into a list
                else:
                    continue
                query_chr = str(cols[1])  #puts all genome2_chrs with synteny to genes in genome1 into a list
                if any(query_chr in s for s in query_api_chrs_final):
                    query_gene = cols[7]  #.rsplit(".", 1)[0]  #puts all genome2_genes with synteny to genes in a genome1 into a list
                else:
                    continue
            if not target_chr in d:
                d[target_chr] = {}  #initializes the nested dictionary-primary level at genome1_chromosome
            if not target_gene in d[target_chr]:
                d[target_chr][target_gene] = {}  #initializes first nesting in dictionary-second level at genome1_genes
            if not query_chr in d[target_chr][target_gene]:
                d[target_chr][target_gene][query_chr] = query_gene  #initializes nested dictionary-third level at genome2_chr

t1 = datetime.now()


"Create hashable data structure to look up gene/CDS name. CDS_name_dict: CDS_key_name: 1 if TRUE"
CDS_name_dict = {}  #initializes the CDS name dictionary for looking up if CDS already have been identified
gff_genes_target = {}

'''Reads GFF from target genome and puts into data structure gff_genes'''
print "Reading GFF File"

with open(gff_import_file, 'r') as g:  # opens gff file
    gffcols = []  #list of parsed gff columns
    chr = []  #initialize list of chromosomes present in genome1 gff file
    for line in g:
        new_line = line.replace(';', '\t')  #makes subdelims universal in gff file from CoGe
        new_line = new_line.replace('Name=', '')  #strips Name= off gene_name in gff file from CoGe
        if new_line[0] != '#' and new_line[0] != '\n':  #selects only lines with CDS information
            gffcols = new_line.split('\t', )  #parses all columns
            if any (gffcols[0] in s for s in target_api_chrs_final) and (gffcols[2] == 'CDS'):  #selects only 'CDS' lines for consideration
                chr = gffcols[0]  #adds genome1_chrs to list
                if not chr in gff_genes_target:
                    gff_genes_target[chr] = []  #initializes chr list in dictionary if chr does not exist yet
                gene_name = str(gffcols[-1])  
                gene_name = gene_name.rstrip("\n")
                gene_name1 = gene_name[9:] #adds targetgenome_genes to list and removes the "ID=" added by CoGe
                start = int(gffcols[3])  #adds targetgenome_gene start bp to list for ordering as integer
                stop = int(gffcols[4])  #adds targetgenome_gene stop bp to list ?for ordering? as integer
                try:
                    CDS_name_dict[str(gene_name1)] == 1
                    continue
                except KeyError:
                    gff_genes_target[chr].append(dict(gene_name=gene_name1, start=start, stop=stop))
                    CDS_name_dict[str(gene_name1)] = 1
                
'''Sorts GFF genes within target chromosomes by start position'''
gff_sort_all = {}
gff_sort_all = natsorted(gff_genes_target)
for chr in gff_genes_target:
    gff_genes_sorted = sorted(gff_genes_target[chr], key=itemgetter('start'))  #Creates dictionary for searching genes against::CONSIDER sorting on midpoint of genes rather than
    gff_genes_target[chr] = gff_genes_sorted

    #CONSIDER WRITING A CHECK PROGRAM TO RETURN TRUE IF ALL VALUES ARE SORTED OR FALSE

t2 = datetime.now()

'''Writes out SynMap dictionary and sorted GFF gene list to document for parsed output'''
print "Writing Out GFF and SynMap Data Structure to CSV"


with open(str(gff_sort_output_file), 'w') as h:
	h.write(str(gff_genes_target))
with open(synmap_dictionary_output_file, 'w+') as i:
    i.write(str(d))

t3 = datetime.now()

'''Determine syntenic gene pairs present and output Raw Data CSV file from parsed data'''
print "Calling Syntenic Pairs and Processing Raw SynMap Data"

"""chr_id(d)
target_lst = natsorted(target_lst)
query_lst = natsorted(query_lst)"""
windanalysis_input_dict = {}

with open(str(fract_bias_raw_output_file), 'w') as csvfile:
    headers = ['Target Chromosome', 'Target Gene Name', 'Gene Order on Target Chromosome']
    headers.extend(query_api_chrs_final)
    headers.extend(query_api_chrs_final)
    writer = csv.writer(csvfile, dialect='excel', delimiter=',', lineterminator='\n')
    writer.writerow(headers)

    if args.allgenes == 'True':
        for tchr in gff_genes_target:
            col0 = chr #writes Pineapple chr number
            count = 0
            for diction in gff_genes_target[tchr]:
                gene = diction['gene_name']
                col1 = gene #writes gene name
                count += 1
                col2 = count #writes pineapple gene number on pineapple chr
                #Find the query chr genes and output to columns: first gff info (col0-3), query chr (col 4-n), query chr-gene (col n+1-m)
                syntenic_query_gene_presence_data = []
                syntenic_query_gene_name = []
                for qchr in query_api_chrs_final:
                    if not tchr in windanalysis_input_dict:
                        windanalysis_input_dict[tchr] = {}  #initializes the nested dictionary-primary level at genome1_chromosome
                    if not qchr in windanalysis_input_dict[tchr]:
                        windanalysis_input_dict[tchr][qchr] = []  #initializes first nesting in dictionary-second level at genome1_genes
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
            col0 = chr #writes Pineapple chr number
            count = 0
            for diction in gff_genes_target[tchr]:
                gene = diction['gene_name']
                col1 = gene #writes target gene name
                count += 1
                col2 = count #writes target gene number on target chr
                #Find the query chr genes and output to columns: first gff info (col0-3), query chr (col 4-n), query chr-gene (col n+1-m)
                syntenic_query_gene_presence_data = []
                syntenic_query_gene_name = []
                for qchr in query_api_chrs_final:
                    if not tchr in windanalysis_input_dict:
                        windanalysis_input_dict[tchr] = {}  #initializes the nested dictionary-primary level at genome1_chromosome
                    if not qchr in windanalysis_input_dict[tchr]:
                        windanalysis_input_dict[tchr][qchr] = []  #initializes first nesting in dictionary-second level at genome1_genes
                    
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

t4 = datetime.now()

# Analysis of Data Structures
# Looks through each chromosome in target genome (tchr), and goes through each gene gene on that target chromosome
# (determined by the GFF uploaded above according to bp position on that target chromosome) to compare to the query genome (query chromsomes).



'''Analysis: for each chromosome in target genome it reads the genes on a chromosome and compares to all query subgenomes'''
print "Running Sliding Window Analysis"

data_output0 = []
data_output1 = []
data_output2 = []
data_output3 = []
output_dict = {}
window_size = 100
window_size_int = 0

#Process windows 100genes/sliding window and 
#output to nested dictionary data structure output_dict{target chr:}{query chr}{window count:retention%}
for tchr in windanalysis_input_dict:
    tchr_counter = tchr
    for qchr in windanalysis_input_dict[tchr]:
        counter = 0
        qchr_counter = qchr
        if not tchr in output_dict:
            output_dict[tchr] = {}  #initializes the nested dictionary-primary level at genome1_chromosome
        if not qchr in output_dict[tchr]:
            output_dict[tchr][qchr] = {}  #initializes first nesting in dictionary-second level at genome1_genes
        window_size_int = int(args.windowsize) #<------------ Have to match up the two arguements to int for logical comparison below, geco sees args.input as 'string'
        try:
            if (int(len(windanalysis_input_dict[tchr][qchr]))) >= window_size_int:   #<--------- changed to args.windowsize
                for each in window(windanalysis_input_dict[tchr][qchr], window_size_int):  #<--------- changed to args.windowsize
                    counter += 1
                    data_output2 = float(sum(each)) / float(window_size)
                    output_dict[tchr][qchr][counter] = round(data_output2*100.) 
        except KeyError:
            continue

#BRASSICA DATA DOES NOT PROCESS         
            
#Sort output_dict for tchr alphanumberic at top level
#if tchr are only integers (not alpha&numeric) the except statement will sort just integers
"""try:
    alphanumbsort = lambda k,v: [k, int(v)]
    output_dict = collections.OrderedDict(sorted(output_dict.items(), key=lambda t: alphanumbsort(*re.match(r'([a-zA-Z]+)(\d+)',t[0]).groups())))
except AttributeError:
    #sorted(output_dict, key=lambda x: output_dict[x])
    #natsorted(output_dict[tchr])
    print "Alphanumeric sort didn't work"
"""
t5 = datetime.now()


"""#Output processed data to a csv file for downstream analysis
with open(retention_calc_output_file, 'wb') as csvf:
    headers = ['Target Chromosome', 'Query Chromosome' 'Window Iteration (x-axis)']
    headers.extend(query_lst)
    writer = csv.writer(csvf, dialect='excel', delimiter=',', lineterminator='\n')
    writer.writerow(headers)
    for tchr in windanalysis_input_dict:
        for qchr in windanalysis_input_dict[tchr]:            
            #Prints into two columns
            writer.writerows(izip( output_dict[tchr][qchr]))"""

print "Length of d " + str(len(d))
print "Length of gff_sort_all " + str(len(gff_sort_all))
print "Length of gff_genes_target " + str(len(gff_genes_target))
print "Length of windanalysis_input_dict " + str(len(windanalysis_input_dict))
print "Length of output_dict " + str(len(output_dict))



#Check output data before graphing
countofchrgraph = 0
listofchrgraph = []

for tchr in output_dict:
    sumofchr = 0
    for qchr in output_dict[tchr]:
        try:
            (max(output_dict[tchr][qchr].itervalues()))
            sumofchr = sumofchr + (max(output_dict[tchr][qchr].itervalues()))
        except ValueError:
            continue
    if sumofchr > 5:
        countofchrgraph += 1
        listofchrgraph.append(str(tchr))

print listofchrgraph
listofchrgraph = natsorted(listofchrgraph)
print listofchrgraph
16911
##Statistics Output NEEDS FIXING

#for tchr in output_dict:
    #for qchr in output_dict[tchr]:
        #print np.mean(output_dict[tchr][qchr])
        #print np.median_grouped(output_dict[tchr][qchr])

"""Plotting"""

print "Plotting FractBias Data"

print "setting up figsize, cols, and gridspec"
#define figure size, column layout, grid layout
figsize = (15, (0.2*len(query_api_chrs_final))+20) #(15, (len(target_api_chrs_final)*2.4))
cols = 2
gs = gridspec.GridSpec(len(output_dict) // cols + 1, cols)

print "setting up tableau palette"

# These are the "Tableau 20" colors as RGB
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

print "scaling palette"
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

#Alternate color palettes Seaborn/Brewer
print "setting palette"
current_palette = sns.color_palette("Set2", 40)

#tableau20
print "Starting to plot figure"

fig = plt.figure(figsize=figsize, frameon=False)
print "set figure size and plt.figure"
subplt_count = -1
ax = []
print "beginning subplots"
for tchr in listofchrgraph:
    subplt_count += 1
    print "Plotting subplot: "+str(subplt_count)
    count = 0
    row = (subplt_count // cols)
    col = subplt_count % cols
    ax.append(fig.add_subplot(gs[row, col]))
    for qchr in output_dict[tchr]:
        count += 1
        try:
            if (max(output_dict[tchr][qchr].itervalues()))>5:
                x = output_dict[tchr][qchr].keys()
                y = output_dict[tchr][qchr].values()
                #Sets up plotting conditions
                ax[-1].spines["top"].set_visible(False)
                ax[-1].spines["right"].set_visible(False)
                ax[-1].spines["left"].set_visible(True)
                ax[-1].spines["bottom"].set_visible(True)
                ax[-1].patch.set_visible(False)
                ax[-1].get_xaxis().tick_bottom()
                ax[-1].get_yaxis().tick_left()
                ax[-1].plot(x, y, color=current_palette[count], lw=2, label=str(qchr))
                ax[-1].set_title(label='Target Chromosome: '+species_name_filter+" "+ tchr, fontweight='bold', fontsize=14, y=1.1, loc='left')
                ax[-1].set_xlabel('Window Iteration\n(Gene number)', fontsize=12, fontweight='bold')
                ax[-1].set_ylabel('% Retention\n(#Syntenic genes/window size)', fontsize=12, fontweight='bold')
                ax[-1].legend(bbox_to_anchor=(1.25, 1.1), loc=1, frameon=False, title="       Query\nChromosome", fontsize=10)
            else:
                continue
        except ValueError:
            continue


fig.subplots_adjust(wspace=0.45, hspace=0.6)
plt.savefig(args.output+"/html/"+"fractbias_figure-" + "-TarID" + str(args.target) + "-TarChrNum" + str(args.numtargetchr) + "-SynDep" + str(args.syndepth) + \
"-QueryID" + str(args.query) + "-QueryChrNum" + str(args.numquerychr) + "-AllGene" + str(args.allgenes) + "-RmRnd" + str(args.remove_random_unknown) + "-WindSize" \
+ str(args.windowsize) + ".png", transparent=True, bbox_inches='tight') 

t6 = datetime.now()


"""Printing Out Benchmarking Times"""
benchmark_api_call = t_api_after - t_api_before
benchmark_synmap_data = t1 - t0
benchmark_gff_data = t2 - t1
benchmark_output_data = t3 - t2
benchmark_processing_data = t4 - t3
benchmark_sliding_window = t5 - t4
benchmark_plotting = t6 - t5
benchmark_total = t6 - t0


"""print "API Call Time: " + str(benchmark_api_call)
print "SynMap Data Processing Time: " + str(benchmark_synmap_data)
print "GFF Data Processing Time: " + str(benchmark_gff_data)
print "Data Structure Output Time: " + str(benchmark_output_data)
print "Raw Data Processing Time: " + str(benchmark_processing_data)
print "Sliding Window Analysis Time: " + str(benchmark_sliding_window)
print "Plotting Graphs Time: " + str(benchmark_plotting)
print "Total Time: " + str(benchmark_total + benchmark_api_call)"""
