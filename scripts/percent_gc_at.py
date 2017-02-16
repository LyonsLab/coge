#!/usr/bin/python3
# coding: utf-8

import fileinput
import glob, os, sys
import collections
import itertools
from os.path import basename
from operator import itemgetter
from itertools import islice
from itertools import groupby

fasta_line_source = fileinput.FileInput(sys.argv[1])
window_size = int(sys.argv[2]);
window_step = int(sys.argv[3]);
file_name = (os.path.splitext(sys.argv[1])[0]) + '_' + sys.argv[2] + '_' + sys.argv[3] + '_out.txt'

###Read each line on the input file and splits sequence headers and sequenced (es a string) into a dictionary
def lines_to_contigs(lines):
    contig_dictionary = None
    current_contig_name = None    
    current_contig_sequence = []
    lines_with_eof_token = itertools.chain(lines, ['>EOF'])
    for line in lines_with_eof_token:
        if line.startswith('>'):
            if current_contig_name:
                contig_sequence_string = ''.join(current_contig_sequence)
                contig_sequence_string = contig_sequence_string.upper()
                contig_dictionary = {'name' : current_contig_name, 'sequence' : contig_sequence_string}
                yield contig_dictionary
                current_contig_sequence = []
            current_contig_name = line[1:].strip()
        else:  
            current_contig_sequence.append(line.strip())

###Groups A and T into AT, and C and G into GC. Other nucleotides (including ambiguous characters are added to their own category) 
def group_nucleotides(nucleotide_seqs):
    nucleotide_group_map = {'A': 'AT','T': 'AT','G': 'GC','C': 'GC','N': 'N','X': 'X', 
                            'M':'M','R':'R','Y':'Y','S':'S','W':'W','K':'K','B':'B','V':'V','D':'D','H':'H'}
    for nucleotide_seq in nucleotide_seqs:
        nucleotide_groups = (nucleotide_group_map[base] for base in nucleotide_seq['sequence'])
        yield {'name': nucleotide_seq['name'], 'sequence': nucleotide_groups}

### Creates a sliding window (size of the window is a global variable). 
### Calculates %AT, %GC, %N and %X inside the window. It also calculates %of ambiguous nucleotides if found.
### After first window, the first nucleotide inside the window is droped and the next in the sequence is added. 
### counter updates to discount the removed nucleotide and count the new nucleotide
### Percentages are recalculated.
### Each window is yielded to a dictionary.
def sliding_percentages(nucleos, window_size, window_step, unique_nucleotides):
    for nucleo in nucleos:
        counter = collections.Counter({nucleotides: 0 for nucleotides in unique_nucleotides})
        sequence = nucleo['sequence']
        counter.update(itertools.islice(sequence, window_size))
        w_start=0
        w_end = window_size
        while w_end < len(sequence):
            counter.subtract(itertools.islice(sequence, window_step))
            w_start += window_step
            w_end += window_step
            counter.update(itertools.islice(sequence, w_end - window_step, w_end))
            percentages = {nucleotides: (c / window_size)*100 for nucleotides, c in counter.items()}
            yield {'name': nucleo['name'], 'sequence': percentages, 'start': w_start, 'end': w_end}
        
### Turns all values calculated on the previous function and found on the dictionary into list
### Output list as a tab delimited format into a .txt file
def contig_list_creator(nucleos):    
    out_file = open(file_name, 'w')
    name = ''
    for nucleo in nucleos:
        if nucleo['name'] != name:
            out_file.write('#Chromosome/contig name, window_start_position, window_end_position, %AT, %GC, %N, %X\n')
            name = nucleo['name']
        out_file.write(name + '\t')
        out_file.write(str(nucleo['start']) + '\t')
        out_file.write(str(nucleo['end']) + '\t')
        out_file.write(str(nucleo['sequence']['AT']) + '\t')
        out_file.write(str(nucleo['sequence']['GC']) + '\t')
        out_file.write(str(nucleo['sequence']['N']) + '\t')
        out_file.write(str(nucleo['sequence']['X']) + '\n')
    out_file.close()

contig_list_creator(sliding_percentages(group_nucleotides(lines_to_contigs(fasta_line_source)), window_size, window_step, ('AT', 'GC', 'N', 'X')))
