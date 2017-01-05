#!/usr/bin/python
# coding: utf-8

# In[1]:

####Import library
import fileinput
#import os, sys
import glob, os, sys
import collections
import itertools
from os.path import basename
from operator import itemgetter
from itertools import islice
from itertools import groupby
import numpy as np
# from plotly import __version__
# from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
# from plotly.graph_objs import Bar, Scatter, Figure, Layout
# init_notebook_mode(connected=True)


# In[2]:

####Global variable 
window_size = 200


# In[3]:

###File input/path

#fasta_line_source = fileinput.FileInput(os.path.expanduser('Plasmodium_chabaudi_chabaudi_strain_AS.faa'))
#fasta_line_source = fileinput.FileInput(os.path.expanduser('Plasmodium_chabaudi_chabaudi_TEST2.faa'))
#fasta_line_source = fileinput.FileInput(os.path.expanduser('Plasmodium_falciparum_IT.faa'))
#fasta_line_source = fileinput.FileInput(os.path.expanduser('contigSrtPlasmoDBberghei'))
#fasta_line_source = fileinput.FileInput(os.path.expanduser('contigSrtPlasmoDBreichenowiCDC'))
#fasta_line_source = fileinput.FileInput(os.path.expanduser('coatchrom1.faa'))
#fasta_line_source = fileinput.FileInput(os.path.expanduser('4chrco.faa'))
#fasta_line_source = fileinput.FileInput(os.path.expanduser('Plasmodium_inui_San_Antonio_test3.faa'))
#fasta_line_source = fileinput.FileInput(os.path.expanduser('test_new.fasta'))
fasta_line_source = fileinput.FileInput(sys.argv[1])

# In[4]:

###Change file name and extension from path to be later used in output

file_name = (os.path.splitext(sys.argv[1])[0]) + '_out.txt'
#print(file_name)


# In[5]:

###Read each line on the input file and splits sequence headers and sequenced (es a string) into a dictionary

def lines_to_contigs(lines):
    contig_dictionary = None
    current_contig_name = None    
    current_contig_sequence = []
    lines_with_eof_token = itertools.chain(lines, ['>EOF'])
    for line in lines_with_eof_token:
        if line.startswith('>'):
            #print('line: {}'.format(line))
            if current_contig_name:
                #print(line)
                #print('current_contig_name: {}'.format(current_contig_name))
                contig_sequence_string = ''.join(current_contig_sequence)
                contig_sequence_string = contig_sequence_string.upper( )
                contig_dictionary = {'name' : current_contig_name, 'sequence' : contig_sequence_string}
                yield contig_dictionary
                current_contig_sequence = []
            current_contig_name = line[1:].strip()
        else:  
            current_contig_sequence.append(line.strip())
            
generator_lines_to_contigs = lines_to_contigs(fasta_line_source)
# for i in generator_lines_to_contigs:
#     print(i)


# In[6]:

###Splits sequence from a string to a tuple. Each nucleotide is separately listed on the tuple

def sequence_split(contigs_iterator):
    for contig in contigs_iterator:
        yield {'name': contig['name'], 'sequence': tuple(contig['sequence'])}

nucleotides2 = sequence_split(generator_lines_to_contigs)
# for i in nucleotides2:
#     print(i)


# In[7]:

###Groups A and T into AT, and C and G into GC. Other nucleotides (including ambiguous characters are added to their own category) 

def group_nucleotides(nucleotide_seqs):
    nucleotide_group_map = {'A': 'AT','T': 'AT','G': 'GC','C': 'GC','N': 'N','X': 'X', 
                            'M':'M','R':'R','Y':'Y','S':'S','W':'W','K':'K','B':'B','V':'V','D':'D','H':'H'}
    for nucleotide_seq in nucleotide_seqs:
        nucleotide_groups = (nucleotide_group_map[base] for base in nucleotide_seq['sequence'])
        yield {'name': nucleotide_seq['name'], 'sequence': nucleotide_groups}

nucleotides_seq = group_nucleotides(nucleotides2)
# for i in nucleotides_seq:
#     print(list(islice(i['sequence'], 10)))


# In[8]:

### Creates a sliding window (size of the window is a global variable). 
### Calculates %AT, %GC, %N and %X inside the window. It also calculates %of ambiguous nucleotides if found.
### After first window, the first nucleotide inside the window is droped and the next in the sequence is added. 
### counter updates to discount the removed nucleotide and count the new nucleotide
### Percentages are recalculated.
### Each window is yielded to a dictionary.

def sliding_percentages(nucleos, window_size, unique_nucleotides):
    for nucleo in nucleos:
        sliding_window = collections.deque(itertools.islice(nucleo['sequence'], window_size), maxlen=window_size)
        counter = collections.Counter({nucleotides: 0 for nucleotides in unique_nucleotides}) 
        counter.update(sliding_window)
        window_percentage = {nucleotides: (counter / window_size)*100 for nucleotides, counter in counter.items()}
        w_start=0
        for base in nucleo['sequence']:
            w_start = w_start+1
            itertools.islice(nucleo['sequence'], window_size, None)
            trailing_nucleotide = sliding_window.popleft()
            counter.subtract([trailing_nucleotide])
            w_end = w_start+window_size
            sliding_window.append(base)
            counter.update([base])
            window_percentage = {nucleotides: (counter / window_size)*100 for nucleotides, counter in counter.items()}
            yield {'name': nucleo['name'], 'sequence': window_percentage, 'start': w_start, 'end': w_end}
        
nucleotides4 = sliding_percentages(nucleotides_seq, window_size, ('AT', 'GC', 'N', 'X'))
# for i in nucleotides4:
#     print(i)


# In[9]:

### Turns all values calculated on the previous function and found on the dictionary into list
### Output list as a tab delimited format into a .txt file

def contig_list_creator(nucleos):    
    iterator_list = [] 
    for nucleo in nucleos:
        dictio = {'name': nucleo['name'],'sequence': nucleo['sequence'], 'start': nucleo['start'], 'end': nucleo['end']} 
        name = nucleo['name']
        start = nucleo['start']
        end = nucleo['end']

        iterator_list.append(dictio.copy())        

    calc_out_file = open(file_name, 'w')
    sorted_iterator_list = sorted(iterator_list, key=itemgetter('name'))
    for key, group in itertools.groupby(sorted_iterator_list, key=lambda x:x['name']):
        contig_iterator_list = list(group)
        #print(contig_iterator_list)
        
        new_AT = []
        new_GC = []
        new_N = []
        new_X = []
        list_w_start = []
        list_w_end = []
        name = key
                
        for contig in contig_iterator_list:
            AT = contig['sequence']['AT'] 
            GC = contig['sequence']['GC']
            N = contig['sequence']['N']
            X = contig['sequence']['X']
            start = contig['start']
            end = contig['end']
            #print(start)
            #print(GC)
        
            new_AT.append(AT)
            new_GC.append(GC)
            new_N.append(N)
            new_X.append(X)
            list_w_start.append(start)
            list_w_end.append(end)
                        
        #print(new_GC)
        #print(list_w_start)
            
        position_new_AT = np.asarray(new_AT)
        position_new_GC = np.asarray(new_GC)
        position_new_N = np.asarray(new_N)
        position_new_X = np.asarray(new_X)
        #print(position_new_GC)
        
        
        list_of_contig_names = []
        xle = len(position_new_X)
        for n in range(xle):
            list_of_contig_names.append(name)
        #print(list_of_contig_names)
               
        #calc_out_file = open('OUT2.txt', 'a')
        
        calc_out_file.writelines('#Chromosome/contig name, window_start_position, window_end_position, %AT, %GC, %N, %X' +'\n')
        out_list = []
        for a, b, c, d, e, f, g in zip(list_of_contig_names, list_w_start, list_w_end, position_new_AT, position_new_GC, position_new_N, position_new_X):
            all_info = (a,str(b),str(c),str(d),str(e),str(f),str(g))
            out_list.append(all_info)
        calc_out_file.writelines('\t'.join(i) + '\n' for i in out_list)        
        #calc_out_file = open('OUT.txt', 'r')
        #print (calc_out_file.read())
        
    calc_out_file.close()
                
    return

contig_list_creator(nucleotides4)

