#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging
import gzip
import os
import re
import sys

logger = logging.getLogger("synmap.gene_order")
logger.setLevel(logging.INFO)

def main(input, output, gid1, gid2, feature1, feature2):
    '''
    Returns a dag file

    If feature1 or feature2 are genomic then the genomic hits need to be
    reordered differently. Otherwise the positions will be set from the
    data within the file.
    '''
    genomic = False
    with open(output, "w") as fp:
        fp.write("#\n")

    if not(os.path.exists(input)):
        logging.error("Input file: (%s) not found.", input)
        
    if os.stat(input).st_size == 0:
        logging.error("Input file: (%s) is empty.", input)

    if not(feature1 == "genomic") and not(feature2 == "genomic"):
        logging.info("No genomic features were specified.")
    else:
        genomic = True

    logging.info("Opening %s for converting to gene order.", input)

    input_file = None
    if input.endswith("gz"):
        logging.info("Opening the input file with gzip.")
        with gzip.open(input, 'r') as fp:
            input_file = fp.readlines()
    else:
        logging.info("Opening the input file.")
        with open(input, 'r') as fp:
            input_file = fp.readlines()

    if not input_file:
        logging.error("The input file (" + input + ") was unable to be opened or is empty.")
        return 1

    if genomic:
        genomic_order = order_genes(input_file, feature1, feature2)

    logging.info("Writing %s in gene order.", args.output)

    with open(output, "a") as fp:
        for line in input_file:
            if line.startswith("#"):
                fp.write(line)
                continue

            fields = line.rstrip("\n\r").split("\t")

            if feature1 == "genomic":
                fields[2] = genomic_order[1][fields[0]][fields[1]]['order']
                fields[3] = genomic_order[1][fields[0]][fields[1]]['order']
            else:
                item = re.split("\|\|", fields[1])

                if not isinstance(item, list) or len(item) < 8:
                    msg = ("The line could not be split for feature (%s)"
                    + " in (%s).\n%s")
                    logger.critical(msg, feature1, gid1, line)
                    sys.exit(1)

                fields[2] = item[7]
                fields[3] = item[7]

            if feature2 == "genomic":
                fields[6] = genomic_order[2][fields[4]][fields[5]]['order']
                fields[7] = genomic_order[2][fields[4]][fields[5]]['order']
            else:
                item = re.split("\|\|", fields[5])

                if not isinstance(item, list) or len(item) < 8:
                    msg = ("The line could not be split for feature (%s)"
                    + " in (%s).\n%s")
                    logger.critical(msg, feature2, gid2, line)
                    sys.exit(1)

                fields[6] = item[7]
                fields[7] = item[7]

            output = reduce(lambda x, y: "\t".join([str(x),str(y)]), fields)
            fp.write("%s\n" % output)
    return 0

# TODO: Some of the work should be able to be done in parallel
# Additionally this shouldn't all be stored in memory
# @by Evan Briones
# @on 3/01/2013
def order_genes(input_file, feature1, feature2):
    '''Returns dictionary containing the genomic order

    The order of the genes in a genome are changed if the feature is genomic.
    The ordering of the genes is done by using the start position of its
    sequence.
    '''

    genomic_order = {1 : {}, 2 : {}}

    for line in input_file:
        if line.startswith("#"):
            continue

        field = line.rstrip("\n\r").split("\t")
        # FIXME: What was this intended to do? It's from the perl script
        #items = map(lambda x: re.split('\|\|', x), [fields[1], fields[5]])

        if not isinstance(field, list) or len(field) < 7:
            logger.critical("The line could not be split\n%s.", line)
            sys.exit(1)

        logger.debug("There are %d fields", len(field))

        if feature1 == "genomic":
            if not field[0] in genomic_order[1]:
                genomic_order[1][field[0]] = {field[1] : {'start' : field[2]}}
            else:
                genomic_order[1][field[0]][field[1]] = {'start' : field[2]}

        if feature2 == "genomic":
            if not field[4] in genomic_order[2]:
                genomic_order[2][field[4]] = {field[5] : {'start' : field[6]}}
            else:
                genomic_order[2][field[4]][field[5]]  = {'start' : field[6]}

    # Compares the gene by its starting position
    cmp_by_start_seq = lambda x : int(x['start'])

    if feature1 == "genomic":
        for chromosome in genomic_order[1].iterkeys():
            genes = genomic_order[1][chromosome].values();
            genes = sorted(genes, key=cmp_by_start_seq)

            for (order, gene) in enumerate(genes, start=1):
                gene['order'] = order
                logging.debug(gene)

    if feature2 == "genomic":
        for chromosome in genomic_order[2].iterkeys():
            genes = genomic_order[2][chromosome].values();
            genes = sorted(genes, key=cmp_by_start_seq)

            for (order, gene) in enumerate(genes, start=1):
                gene['order'] = order
                logging.debug(gene)

    return genomic_order

def convert_to_genomic_position(genomic_file, output):
    '''
    Return output file with genomic positions

    Convert the input file from genomic order to genomic
    position.
    '''
    with open(genomic_file, 'r') as fp, open(output, 'a') as out:
        out.write('#\n')
        for line in fp:
            if line.startswith('#'):
                data = line
            else:
                items = line.rstrip("\n\r").split("\t")
                pos1 = re.split('\|\|', items[1])
                pos2 = re.split('\|\|', items[5])

                (start, stop) = (pos1[1], pos1[2])
                if pos1[4] and re.match('-', pos1[4]):
                    (start, stop) = (stop, start)

                items[2] = start
                items[3] = stop

                (start, stop) = (pos2[1], pos2[2])

                if pos2[4] and re.match('-', pos2[4]):
                    (start, stop) = (stop, start)

                items[6] = start
                items[7] = stop
                data = "\t".join(items) + "\n"

            out.write(data)
    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='''
        Reorders the genes contained in the dag file by ascending order.
        ''')

    parser.add_argument("input",
                        help="dag file to be ordered")
    parser.add_argument("output",
                        help="dag file that is produced.")
    parser.add_argument("--gid1",
                        help="The first genome id.")
    parser.add_argument("--gid2",
                        help="The second genome id.")
    parser.add_argument("--feature1",
                        help="The feature of the first genome")
    parser.add_argument("--feature2",
                        help="The feature of the second genome")
    parser.add_argument("--positional", action="store_true",
                        help="Convert from genomic order to genomic position")

    parser.add_argument("--loglevel",
                        help="Sets the log level verboseness")
    args = parser.parse_args()

    if args.loglevel:
        log_level = getattr(logging, args.loglevel.upper(), None)

        if not isinstance(log_level, int):
            raise ValueError('Invalid log level: %s' % log_level)
        logger.setLevel(log_level)

    if args.positional:
        ret = convert_to_genomic_position(args.input, args.output)
    else:
        ret = main(args.input, args.output, args.gid1, args.gid2, args.feature1,
             args.feature2)

    sys.exit(ret)
