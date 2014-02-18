from __future__ import division
import argparse
import operator
import os
import math
import re
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('cufflinks', help="The input cufflinks file")
    parser.add_argument('output', help="The output filename")

    args = parser.parse_args()

    gene_exp = []
    max_fpkm = 0.0

    with open(args.cufflinks) as fh, open(args.output, "w") as out:
        fh.readline()

        for x in fh:
            y = x.strip().split('\t')
            myname = y[0]
            myloc = y[6]

            try:
                mychr, mystart, mystop = re.split(':|-', myloc);
            except:
                continue

            myfpkm = math.log(float(y[9]) + 1, 10)

            if myfpkm == 0:
                continue

            max_fpkm = max(max_fpkm, myfpkm)
            gene_exp.append([mychr, int(mystart), int(mystop), 1, myfpkm])

        # PRINT
        out.write("#CHR,START,STOP,STRAND,VALUE1(0-1),VALUE2(ANY-ANY)\n")

        for gene in sorted(gene_exp, key=operator.itemgetter(1)):
            ref, start, stop, strand, fpkm = gene

            # normalize fpkm value
            fpkm_norm  = fpkm / max_fpkm

            items = [str(item) for item in [ref, start, stop, strand, fpkm_norm, fpkm]]
            out.write(",".join(items) + "\n")


