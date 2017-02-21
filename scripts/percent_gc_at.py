#!/usr/bin/python3
# coding: utf-8

import fileinput
import glob, os, sys

in_lines = fileinput.FileInput(sys.argv[1])
window_size = int(sys.argv[2]);
window_step = int(sys.argv[3]);
file_name = (os.path.splitext(sys.argv[1])[0]) + '_' + sys.argv[2] + '_' + sys.argv[3] + '_out.txt'

### Read each line on the input file into an uppercased string
def read_fasta(lines):
    seqs = []
    lines.readline() # skip header
    for line in lines:
        seqs.append(line.strip())
    return ''.join(seqs).upper()

### Creates a sliding window (size of the window is a global variable). 
### Calculates %AT, %GC, %N and %X inside the window.
### After first window, the first window_step nucleotides inside the window are droped and the next window_step nucleotides in the sequence are added. 
### Each window is yielded as a line for the output file
def sliding_percentages(fasta, window_size, window_step):
    at = 0
    gc = 0
    n = 0
    x = 0
    for i in range(0, window_size):
        c = fasta[i]
        if c == 'A' or c == 'T' or c == 'S':
            at += 1
        elif c == 'C' or c == 'G' or c == 'W':
            gc += 1
        elif c == 'N':
            n += 1
        elif c == 'X':
            x += 1
    w_start=0
    w_end = window_size
    yield str(w_start) + '\t' + str(w_end) + '\t' + str((at/window_size)*100) + '\t' + str((gc/window_size)*100) + '\t' + str((n/window_size)*100) + '\t' + str((x/window_size)*100) + '\n'
    while w_start < len(fasta) - window_size - window_step:
        for i in range(w_start, w_start + window_step):
            c = fasta[i]
            if c == 'A' or c == 'T' or c == 'S':
                at -= 1
            elif c == 'C' or c == 'G' or c == 'W':
                gc -= 1
            elif c == 'N':
                n -= 1
            elif c == 'X':
                x -= 1
        w_start += window_step
        for i in range(w_end, w_end + window_step):
            c = fasta[i]
            if c == 'A' or c == 'T' or c == 'S':
                at += 1
            elif c == 'C' or c == 'G' or c == 'W':
                gc += 1
            elif c == 'N':
                n += 1
            elif c == 'X':
                x += 1
        w_end += window_step
        yield str(w_start) + '\t' + str(w_end) + '\t' + str((at/window_size)*100) + '\t' + str((gc/window_size)*100) + '\t' + str((n/window_size)*100) + '\t' + str((x/window_size)*100) + '\n'
        
### Output list as a tab delimited format into a .txt file
def write_windows(out_lines):    
    out_file = open(file_name, 'w')
    out_file.write('#window start\twindow end\t%AT\t%GC\t%N\t%X\n')
    for line in out_lines:
        out_file.write(line)
    out_file.close()

write_windows(sliding_percentages(read_fasta(in_lines), window_size, window_step))
