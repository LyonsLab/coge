#!/usr/bin/env python

import numpy as np

from argparse import ArgumentParser
from collections import defaultdict, deque
from datetime import datetime
from json import load, dump
from math import log10
from natsort import natsorted
from os import path
from sklearn.cluster import DBSCAN

#from pprint import pprint

start = datetime.now()
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                   ______                 __  _                                                      #
#                                  / ____/_  ______  _____/ /_(_)___  ____  _____                                     #
#                                 / /_  / / / / __ \/ ___/ __/ / __ \/ __ \/ ___/                                     #
#                                / __/ / /_/ / / / / /__/ /_/ / /_/ / / / (__  )                                      #
#                               /_/    \__,_/_/ /_/\___/\__/_/\____/_/ /_/____/                                       #
#                                                                                                                     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def minimizeGenomes(genomes_object, cutoff):
    """
    Remove chromosomes below cutoff size from genome_object.

    :param genomes_object:
    :param cutoff:
    :return: genomes_object with chromosomes < cutoff removed
    """
    for genome in genomes_object:
        print("Parsing %s genome (gid %s)" % (genomes_object[genome]["organism"]["name"], genome))
        print("--> Initial chromosome count: %d" % len(genomes_object[genome]["chromosomes"]))
        genomes_object[genome]["chromosomes"] = [genomes_object[genome]["chromosomes"][i] \
                                                 for i in range(len(genomes_object[genome]["chromosomes"])) \
                                                 if genomes_object[genome]['chromosomes'][i]["length"] >= cutoff]
        parsed_count = len(genomes_object[genome]["chromosomes"])
        genomes_object[genome]["chromosome_count"] = parsed_count
        print("--> Parsed chromosome count: %d" % parsed_count)

    return genomes_object

def minimizeSyntenicPoints(sp_object, genomes_object):
    """
    Remove syntenic points if containing chromosomes not in genome_object.

    :param sp_object: syntenic_points object
    :param genomes_object: genomes object (see scaleSortGenomes for full description)
    :return: syntenic points object with points from absent chromosomes removed.
    """
    # Extract GIDs from syntenic points object.
    gid1 = sp_object.keys()[0]
    gid2 = sp_object[gid1].keys()[0]

    # Extract chromosome lists from genomes object.
    sp1_chlist = [ch["name"] for ch in genomes_object[gid1]["chromosomes"]]
    sp2_chlist = [ch["name"] for ch in genomes_object[gid2]["chromosomes"]]

    # Remove any points on species 1 chromosomes if chromosome is absent in genomes object.
    sp_object[gid1][gid2] = {sp1ch: sp_object[gid1][gid2][sp1ch] \
                             for sp1ch in sp_object[gid1][gid2] \
                             if sp1ch in sp1_chlist}
    # Remove any points on species 2 chromosomes if chromosome is absent in genomes object.
    for sp1ch in sp_object[gid1][gid2]:
        sp_object[gid1][gid2][sp1ch] = {sp2ch: sp_object[gid1][gid2][sp1ch][sp2ch] \
                                        for sp2ch in sp_object[gid1][gid2][sp1ch] \
                                        if sp2ch in sp2_chlist}
    return sp_object


def syntenyParser(sp_object, genomes_object, cutoff):
    """
    Parses chromosomes under a minimum syntenic points count.

    :param sp_object: syntenic_points object
    :param genomes_object: genomes_object
    :param cutoff: minimum syntenic points
    :return:
    """
    gid1 = sp_object.keys()[0]
    gid2 = sp_object[gid1].keys()[0]

    sp1ch_hits = {}
    for ch in genomes_object[gid1]["chromosomes"]:
        sp1ch_hits[ch["name"]] = 0
    sp2ch_hits = {}
    for ch in genomes_object[gid2]["chromosomes"]:
        sp2ch_hits[ch["name"]] = 0

    for sp1ch in sp_object[gid1][gid2]:
        #if sp1ch not in sp1ch_hits:
        #    sp1ch_hits[sp1ch] = 0
        for sp2ch in sp_object[gid1][gid2][sp1ch]:
            #if sp2ch not in sp2ch_hits:
            #    sp2ch_hits[sp2ch] = 0
            if sp1ch in sp1ch_hits and sp2ch in sp2ch_hits:
                sp1ch_hits[sp1ch] += len(sp_object[gid1][gid2][sp1ch][sp2ch])
                sp2ch_hits[sp2ch] += len(sp_object[gid1][gid2][sp1ch][sp2ch])

    print("Species %s" % gid1)
    print("--> Before %d" % len(genomes_object[gid1]["chromosomes"]))
    genomes_object[gid1]["chromosomes"] = [genomes_object[gid1]["chromosomes"][i] for i in range(len(genomes_object[gid1]["chromosomes"])) if sp1ch_hits[genomes_object[gid1]["chromosomes"][i]["name"]] > cutoff]
    print("--> After %d" % len(genomes_object[gid1]["chromosomes"]))

    print("Species %s" % gid2)
    print("--> Before %d" % len(genomes_object[gid2]["chromosomes"]))
    genomes_object[gid2]["chromosomes"] = [genomes_object[gid2]["chromosomes"][i] for i in range(len(genomes_object[gid2]["chromosomes"])) if sp2ch_hits[genomes_object[gid2]["chromosomes"][i]["name"]] > cutoff]
    print("--> After %d" % len(genomes_object[gid2]["chromosomes"]))

    return genomes_object


def scaleSortGenomes(genomes_object):
    """
    Scale & Sort Genomes

    :param genomes_object: dictionary of genome objects {"<gid1>": <genome1_object>,
                                                         ... ,
                                                         "<gidN>": <genomeN_object}
    :return genome_size_absolute: dictionary of absolute (bp) genome sizes {"<gid1>": <genome1_length_bp,
                                                                            ... ,
                                                                            "<gidN>": <genomeN_length_bp}
    :return genome_size_relative: dictionary of scaled (largest = 1) genome sizes {"<gid1>": <genome1_length_scaled,
                                                                                   ... ,
                                                                                   "<gidN>": <genomeN_length_scaled}
    :return chrs_sorted: dictionary of genome keys and sorted chromosome list values.
    :return chrs_relative: dictionary of scaled chromosomes (sum = 100) {"<gid1>: {"<ch1>": <scaled_len>, ... }, ... }
    :return chrs_absolute: dictionary of chromosomes absolute (bp) length {"<gid1>: {"<ch1>": <bp_len>, ... }, ... }
    """
    ss_start = datetime.now()

    # Variables to store data.
    genome_size_absolute = {}
    genome_size_relative = {}
    chrs_sorted = {}
    chrs_relative = {}
    chrs_absolute = {}

    for gid, ginfo in genomes_object.iteritems():
        # Calculate total genome length.
        total_length = 0
        for chr in ginfo["chromosomes"]:
            total_length += chr["length"]
        # Save genome length.
        genome_size_absolute[gid] = total_length

        # Sort chromosomes.
        if sort_method == "length":
            chrs_sorted[gid] = natsorted(ginfo["chromosomes"], key=lambda chromosome: chromosome["length"], reverse=True)
        elif sort_method == "name":
            chrs_sorted[gid] = natsorted(ginfo["chromosomes"], key=lambda chromosome: chromosome["name"])
        else:
            chrs_sorted[gid] = ginfo["chromosomes"]

        # Calculate chromosome percentages.
        chrs_relative[gid] = {}
        chrs_absolute[gid] = {}
        total = 0.0
        for chr in chrs_sorted[gid]:
            length = chr["length"]
            percent = float(length)/total_length * 100
            total += percent
            chrs_relative[gid][chr["name"]] = percent
            chrs_absolute[gid][chr["name"]] = length

    ## ----- Calculate relative genome sizes. ----- ##
    genome_size_relative = genome_size_absolute.copy()
    max_genome = max(genome_size_absolute, key=lambda k: genome_size_absolute[k])
    for gid,glen in genome_size_absolute.iteritems():
        genome_size_relative[gid] = float(glen)/genome_size_absolute[max_genome]

    ss_end = datetime.now()
    print "Genome Sorting & Scaling Complete (%s)" % str(ss_end-ss_start)
    return genome_size_absolute, genome_size_relative, chrs_sorted, chrs_relative, chrs_absolute


def getPointCoords(syntenic_points, relative_chr, absolute_chr, relative_genome, sorted_chs):
    """
    Get Point Coordinates, extracts points from syntenic_points and converts to graph space (~0-100).

    :param syntenic_points: syntenic points from file.
    :param relative_chr: "chrs_relative", from scale_sort_genomes()
    :param absolute_chr: "chrs_absolute", from scale_sort_genomes()
    :param relative_genome: "genome_size_relative" from scale_sort_genomes()
    :param sorted_chs: "chrs_sorted" from scale_sort_genomes()
    :return ids: tuple of genome id's from point extraction, ordered to represent coordinates
    :return coords: coordinates np-array (n, 6) [[loc1_gid1, loc1_gid2, fid1_gid1, fid1_gid2, Kn, Ks],
                                                 ...,
                                                 [locn_gid1, locn_gid2, fidn_gid1, fidn_gid2, Kn, Ks]]
    """
    coord_start = datetime.now()

    # Assign order of gids.
    gid1 = syntenic_points.keys()[0]
    gid2 = syntenic_points[gid1].keys()[0]
    # Variable to store all coordinates.
    coords = []  # Each coord is in form [sp1_val, sp2_val, sp1_fid, sp2_fid, Kn, Ks]
    downloads = []

    # Reassign syntenic_points to just points.
    syntenic_points = syntenic_points[gid1][gid2]

    for sp1chr in syntenic_points:
        # Assign scale & length variables for current species 1 chromosome.
        sp1chr_scale = relative_chr[gid1][sp1chr]
        sp1chr_length = absolute_chr[gid1][sp1chr]
        sp1chr_relative_length = sp1chr_scale * relative_genome[gid1]
        for sp2chr in syntenic_points[sp1chr]:
            # Assign scale & length variables for current species 2 chromosome.
            sp2chr_scale = relative_chr[gid2][sp2chr]
            sp2chr_length = absolute_chr[gid2][sp2chr]
            sp2chr_relative_length = sp2chr_scale * relative_genome[gid2]

            # Calculate sum of all previous chromosomes
            sp1_sum = 0.0
            for c in sorted_chs[gid1]:
                if c["name"] == sp1chr:
                    break
                else:
                    sp1_sum += relative_chr[gid1][c["name"]] * relative_genome[gid1]
            sp2_sum = 0.0
            for c in sorted_chs[gid2]:
                if c["name"] == sp2chr:
                    break
                else:
                    sp2_sum += relative_chr[gid2][c["name"]] * relative_genome[gid2]

            # Calculate each hit coordinate, append to coordinates list.
            for hit in syntenic_points[sp1chr][sp2chr]:
                sp1_raw = (int(hit[0]) + int(hit[1])) / 2
                sp1_fid = hit[4][gid1]["db_feature_id"]
                sp2_raw = (int(hit[2]) + int(hit[3])) / 2
                sp2_fid = hit[4][gid2]["db_feature_id"]

                sp1_val = (sp1chr_relative_length * (float(sp1_raw) / sp1chr_length)) + sp1_sum
                sp2_val = (sp2chr_relative_length * (float(sp2_raw) / sp2chr_length)) + sp2_sum

                try:
                    Kn = log10(float(hit[4]["Kn"]))
                except ValueError:
                    Kn = 'NA'
                try:
                    Ks = log10(float(hit[4]["Ks"]))
                except ValueError:
                    Ks = 'NA'

                coord = [sp1_val, sp2_val, sp1_fid, sp2_fid, Kn, Ks]
                coords.append(coord)

                # Build downloadable data dictionary.
                #p1_gid||chr1||start1||stop1||name1||strand1||type1||db_feature_id1
                downloadA = "||".join([gid1, sp1chr, hit[0], hit[1], hit[4][gid1]["name"], hit[4][gid1]["strand"], hit[4][gid1]["type"], hit[4][gid1]["db_feature_id"]])
                downloadB = "||".join([gid2, sp2chr, hit[2], hit[3], hit[4][gid2]["name"], hit[4][gid2]["strand"], hit[4][gid2]["type"], hit[4][gid2]["db_feature_id"]])
                download = [downloadA, downloadB]
                downloads.append(download)

    ids = (gid1, gid2)
    coords = np.array(coords)
    downloads = np.array(downloads, dtype=str, copy=True)

    coord_end = datetime.now()
    print("Point Extraction Complete (%s)" % str(coord_end-coord_start))
    print("--> %s vs %s" % (gid1, gid2))
    return ids, coords, downloads


def reorganizeArray(id_A, id_B, gids_tuple, unordered_coordinates ):
    """
    From a tuple of genome ID's representing a coordinate array's structure (gidI, gidJ).
    and the corresponding coordinate np-array (n x 4) of design [[locI1, locJ1, fidI1, fidJ2, Kn, Ks],
                                                                  ...
                                                                 [locIn, locJn, fidIn, fidJn, Kn, Ks]]
    * Both of these are returned by getPointCoords - as "ids" and "coords", respectively.

    Restructures to represent an (A, B) order (swaps columns 0,1 and 2,3 if A=J and B=I, else unchanged)

    :param id_A: GID of genome to appear first
    :param id_B: GID of genome to appear second
    :param gids_tuple: Ordered tuple of GIDs in original coordinate array (returned by getPointCoords)
    :param unordered_coordinates: coordinate (ordered or unordered) to map order to (returned by getPointCoords)
    :return: coordinate object with [A_loc, B_loc, A_fid, B_fid] order enforced.
    """
    reorg_start = datetime.now()

    if id_A == gids_tuple[0] and id_B == gids_tuple[1]:
        ordered_coordinates = unordered_coordinates
    elif id_A == gids_tuple[1] and id_B == gids_tuple[0]:
        ordered_coordinates = unordered_coordinates[:,[1,0,3,2,4,5]]
    else:
        ordered_coordinates = None
        print "ERROR: reorganizeArray (%s, %s)" % (str(id_A), str(id_B))
        exit()

    reorg_end = datetime.now()
    print("Coordinates Reorganization Complete (%s)" % str(reorg_end-reorg_start))
    print("--> (%s,%s) to (%s,%s)" % (gids_tuple[0], gids_tuple[1], id_A, id_B))

    return ordered_coordinates


def reorganizeDownloads(id_A, id_B, downloadsList):
    """

    :param id_A:
    :param id_B:
    :param downloadsList:
    :return:
    """
    reorg_start = datetime.now()
    first = downloadsList[0][0].split("||")[0]
    second = downloadsList[0][1].split("||")[0]
    if str(id_A) == str(first):
        reorg_end = datetime.now()
        print("Downloads Reorganization Complete (%s)" % str(reorg_end-reorg_start))
        print("--> (%s,%s) to (%s,%s)" % (first, second, id_A, id_B))
        return downloadsList
    elif str(id_A) == str(second):
        reorg_end = datetime.now()
        print("Downloads Reorganization Complete (%s)" % str(reorg_end-reorg_start))
        print("--> (%s,%s) to (%s,%s)" % (first, second, id_A, id_B))
        return [[e[1], e[0]] for e in downloadsList]
    else:
        print("Reordering Error (Downloads)")
        exit()


def mergePoints(setXY, setXZ, setYZ, xid, yid, zid, downloadData):
    """
    Merge XY, XZ, YZ points.

    :param setXY: XY-points object ( (n,4) np-array [[x_loc, y_loc, x_fid, y_fid]...] ).
    :param setXZ: XZ-points object ( (n,4) np-array [[x_loc, z_loc, x_fid, z_fid]...] ).
    :param setYZ: YZ-points object ( (n,4) np-array [[y_loc, z_loc, y_fid, z_fid]...] ).
    :return: common hits ( (n,12) np-array [[x, y, z, x_fid, y_fid, z_fid, XY_Kn, XY_Ks, XZ_Kn, XZ_Ks, YZ_Kn, YZ_Ks]...] ).
    """
    merge_start = datetime.now()

    hits = []
    errors = []
    downloads = []

    # Build indexed XZ sets.
    setXZ_xindexed = defaultdict(list)
    for s in setXZ:
        setXZ_xindexed[s[2]].append(s)

    # Build indexed YZ sets.
    setYZ_yindexed = defaultdict(list)
    for s in setYZ:
        setYZ_yindexed[s[2]].append(s)

    # Find hits.
    for XY in setXY:
        XZ_X_matches = setXZ_xindexed[XY[2]]
        YZ_Y_matches = setYZ_yindexed[XY[3]]
        if len(XZ_X_matches) > 0 and len(YZ_Y_matches) > 0:
            for XZ in XZ_X_matches:
                for YZ in YZ_Y_matches:
                    if XY[2] == XZ[2] and XY[3] == YZ[2] and XZ[3] == YZ[3]:
                        if XY[0] == XZ[0] and XY[1] == YZ[0] and YZ[1] == XZ[1]:
                            #hit = [X, Y,Z, X_FID, Y_ID, Z_ID, XY_Kn, XY_Ks, XZ_Kn, XZ_Ks, YZ_Kn, YZ_Ks]
                            hit = [XY[0], YZ[0], XZ[1], XY[2], YZ[2], XZ[3],
                                   XY[4], XY[5],
                                   XZ[4], XZ[5],
                                   YZ[4], YZ[5]]
                            hits.append(hit)

                            mut = "||".join([XY[4], XY[5], XZ[4], XZ[5], YZ[4], YZ[5]])
                            down = [downloadData[xid][XY[2]], downloadData[yid][YZ[2]], downloadData[zid][XZ[3]], mut]
                            downloads.append(down)
                        else:
                            err = [XY, XZ, YZ]
                            errors.append(err)

    hits = np.array(hits)
    errors = np.array(errors)
    merge_end = datetime.now()

    print "Merge Complete (%s)" % str(merge_end-merge_start)
    print "--> Hits Found: %s" % str(hits.shape[0])
    print "--> Errors Recorded: %s" % str(errors.shape[0])
    return hits, downloads


def buildAxisTicks(sorted_chrs, relative_genome_size, relative_chromosomes):
    """
    Build Axis Ticks, for Plotly Plots.

    :param sorted_chrs: Sorted list of chromosomes.
    :param relative_genome_size: Relative total size of genome.
    :param relative_chromosomes: Relative chromosome sizes.
    :return axis: index-matched dictionary of {"positions": [], "labels": []}.
    :return length: Total relative length of axis.
    """
    axis = {"positions": [], "labels": []}
    length = 0
    for chr in sorted_chrs:
        name = chr["name"]
        axis["positions"].append(length)
        axis["labels"].append(name)
        length += relative_genome_size * relative_chromosomes[name]

    return axis, length


def buildDownloadData(dataXY, dataXZ):
    data = {}
    xgid = dataXY[0][0].split("||")[0] # For X data
    data[xgid] = {}
    ygid = dataXY[0][1].split("||")[0] # For Y data
    data[ygid] = {}
    zgid = dataXZ[0][1].split("||")[0] # For Z data
    data[zgid] = {}

    for lineitem in dataXY:
        xdata = lineitem[0]
        xfid = xdata.split("||")[7]
        data[xgid][xfid] = xdata
        ydata = lineitem[1]
        yfid = ydata.split("||")[7]
        data[ygid][yfid] = ydata
    for lineitem in dataXZ:
        zdata = lineitem[1]
        zfid = zdata.split("||")[7]
        data[zgid][zfid] = zdata

    return data


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                          ____        __  _                                                          #
#                                         / __ \____  / /_(_)___  ____  _____                                         #
#                                        / / / / __ \/ __/ / __ \/ __ \/ ___/                                         #
#                                       / /_/ / /_/ / /_/ / /_/ / / / (__  )                                          #
#                                       \____/ .___/\__/_/\____/_/ /_/____/                                           #
#                                           /_/                                                                       #
#                                                                                                                     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## ----- Argument parser. ----- ##
parser = ArgumentParser(description="Threeway SynMap Graph Builder")

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
parser.add_argument('-S',
                    type=str,
                    choices=["name", "length"],
                    default="name",
                    help="Sort method ('name' or 'length)")
parser.add_argument('-C',
                    action="store_true",
                    help="Flag to enable thinning by clustering")
parser.add_argument('-C_eps',
                    type=str,
                    default="0.5",
                    help="EPS (maximum sample distance for neighborhood) for DBSCAN clustering")
parser.add_argument('-C_ms',
                    type=int,
                    default=5,
                    help="Minimum samples in neighborhood for DBSCAN clustering")
parser.add_argument('-R',
                    type=str,
                    choices=["kn", "ks", "knks"],
                    required=False,
                    help="Enable Kn/Ks ratio cutoff ('kn', 'ks', 'knks'")
parser.add_argument('-Rby',
                    type=str,
                    choices=["xy", "xz", "yz", "mean", "median"],
                    default="xy",
                    help="Comparison to base ratio cutoff ('xy', 'xz', 'yz', 'mean', 'median')")
parser.add_argument('-Rmin',
                    type=str,
                    default="-1.0",
                    help="Minimum log10(ratio) cutoff (log10)")
parser.add_argument('-Rmax',
                    type=str,
                    default="1.0",
                    help="Maximum log10(ratio) cutoff (log10)")
parser.add_argument('-ml',
                    type=int,
                    default=0,
                    help="Minimum chromosome size (bp) for consideration")
parser.add_argument('-ms',
                    type=int,
                    default=0,
                    help="Minimum syntenic points per chromosome for consideration")
parser.add_argument('-o',
                    type=str,
                    required=True,
                    help="Path to output data folder")

args = parser.parse_args()

## -----Assign options. ----- ##
# Sorting
sort_method = args.S

# Clustering
remove_unclustered = args.C
eps = float(args.C_eps)
min_samples = args.C_ms

# Parsing
min_length = args.ml
min_synteny = args.ms
# TODO: Kn/Ks Parse

# Kn/Ks Ratio Cutoff
ratio_cutoff = args.R
ratio_by = args.Rby
rcutoff_min = float(args.Rmin)
rcutoff_max = float(args.Rmax)

#except:
#ratio_cutoff = args.R
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                 __  ___      _          _____           _       __                                  #
#                                /  |/  /___ _(_)___     / ___/__________(_)___  / /_                                 #
#                               / /|_/ / __ `/ / __ \    \__ \/ ___/ ___/ / __ \/ __/                                 #
#                              / /  / / /_/ / / / / /   ___/ / /__/ /  / / /_/ / /_                                   #
#                             /_/  /_/\__,_/_/_/ /_/   /____/\___/_/  /_/ .___/\__/                                   #
#                                                                      /_/                                            #
#                                                                                                                     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## ----- Assign GIDs. ----- ##
X_GID = args.xid
Y_GID = args.yid
Z_GID = args.zid


## ----- Load data files. ----- ##
file1 = load(open(args.i1, 'r'))
file2 = load(open(args.i2, 'r'))
file3 = load(open(args.i3, 'r'))


## ----- Consolidate genome information. ----- ##
genomes = file1["genomes"].copy()
genomes.update(file2["genomes"])
genomes.update(file3["genomes"])

## ----- Save syntenic points. ----- ##
file1_sp = file1["syntenic_points"]
file2_sp = file2["syntenic_points"]
file3_sp = file3["syntenic_points"]

## ----- Process parsing options. ----- ##
# Minimum length.
if min_length > 0:
    genomes = minimizeGenomes(genomes, min_length)
    file1_sp = minimizeSyntenicPoints(file1_sp, genomes)
    file2_sp = minimizeSyntenicPoints(file2_sp, genomes)
    file3_sp = minimizeSyntenicPoints(file3_sp, genomes)

# Blank chromosoomes.
if min_synteny > 0:
    genomes = syntenyParser(file1_sp, genomes, min_synteny)
    genomes = syntenyParser(file2_sp, genomes, min_synteny)
    genomes = syntenyParser(file3_sp, genomes, min_synteny)
    file1_sp = minimizeSyntenicPoints(file1_sp, genomes)
    file2_sp = minimizeSyntenicPoints(file2_sp, genomes)
    file3_sp = minimizeSyntenicPoints(file3_sp, genomes)

## ----- Calculate total genome sizes, relative chromosome sizes, and sort. ----- ##
g_size_abs, g_size_rel, c_sort, c_rel, c_abs = scaleSortGenomes(genomes)

## ----- Get point coordinates from each file. ----- ##
gids1, coordinates1, downloads1 = getPointCoords(file1_sp, c_rel, c_abs, g_size_rel, c_sort)
gids2, coordinates2, downloads2 = getPointCoords(file2_sp, c_rel, c_abs, g_size_rel, c_sort)
gids3, coordinates3, downloads3 = getPointCoords(file3_sp, c_rel, c_abs, g_size_rel, c_sort)

## ----- Determine order of arrays & restructure to represent correct XYZ ordering. ----- ##
# XY order.
if X_GID in gids1 and Y_GID in gids1:
    coordinatesXY = reorganizeArray(X_GID, Y_GID, gids1, coordinates1)
    downloadsXY = reorganizeDownloads(X_GID, Y_GID, downloads1)
elif X_GID in gids2 and Y_GID in gids2:
    coordinatesXY = reorganizeArray(X_GID, Y_GID, gids2, coordinates2)
    downloadsXY = reorganizeDownloads(X_GID, Y_GID, downloads2)
elif X_GID in gids3 and Y_GID in gids3:
    coordinatesXY = reorganizeArray(X_GID, Y_GID, gids3, coordinates3)
    downloadsXY = reorganizeDownloads(X_GID, Y_GID, downloads3)
else:
    coordinatesXY = None
    print "ERROR: XY coordinates unidentified."
    exit()
# XZ order.
if X_GID in gids1 and Z_GID in gids1:
    coordinatesXZ = reorganizeArray(X_GID, Z_GID, gids1, coordinates1)
    downloadsXZ = reorganizeDownloads(X_GID, Z_GID, downloads1)
elif X_GID in gids2 and Z_GID in gids2:
    coordinatesXZ = reorganizeArray(X_GID, Z_GID, gids2, coordinates2)
    downloadsXZ = reorganizeDownloads(X_GID, Z_GID, downloads2)
elif X_GID in gids3 and Z_GID in gids3:
    coordinatesXZ = reorganizeArray(X_GID, Z_GID, gids3, coordinates3)
    downloadsXZ = reorganizeDownloads(X_GID, Z_GID, downloads3)
else:
    coordinatesXZ = None
    print "ERROR: XZ coordinates unidentified."
    exit()
# YZ order.
if Y_GID in gids1 and Z_GID in gids1:
    coordinatesYZ = reorganizeArray(Y_GID, Z_GID, gids1, coordinates1)
    downloadsYZ = reorganizeDownloads(Y_GID, Z_GID, downloads1)
elif Y_GID in gids2 and Z_GID in gids2:
    coordinatesYZ = reorganizeArray(Y_GID, Z_GID, gids2, coordinates2)
    downloadsYZ = reorganizeDownloads(Y_GID, Z_GID, downloads2)
elif Y_GID in gids3 and Z_GID in gids3:
    coordinatesYZ = reorganizeArray(Y_GID, Z_GID, gids3, coordinates3)
    downloadsYZ = reorganizeDownloads(Y_GID, Z_GID, downloads3)
else:
    coordinatesYZ = None
    print "ERROR: YZ coordinates unidentified."
    exit()

## ----- Merge & calculate coordinates. ----- ##
downloadData = buildDownloadData(downloadsXY, downloadsXZ)

coordinates, downloads = mergePoints(coordinatesXY, coordinatesXZ, coordinatesYZ, X_GID, Y_GID, Z_GID, downloadData) # TODO HERE!

## ----- Clustering: DBSCAN ----- ##
# If remove_unclustered = True
if remove_unclustered:
    cluster_start = datetime.now()

    # Build vectors of only coordinate space.
    pts = coordinates[:,0:3]

    # Cluster with DBSCAN.
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(pts)
    labels = db.labels_

    # Build vectors of clustered points including all coordinate information.
    clusters = coordinates[labels != -1]
    sparse = coordinates[labels == -1]

    downloads = np.array(downloads, dtype=str)
    downloads = downloads[labels != -1]

    # Print concluding message.
    cluster_end = datetime.now()
    print("Clustering Complete (%s)" % str(cluster_end-cluster_start))
    cluster_count = len(set(labels)) - (1 if -1 in labels else 0)
    print('--> Estimated number of clusters: %d' % cluster_count)
    print('--> Total Point Count: %s') % str(coordinates.shape[0])
    print('--> Clustered Point Count: %s') % str(clusters.shape[0])
    print('--> Removed Point Count: %s') % str(sparse.shape[0])
# If remove_unclustered = False
else:
    clusters = coordinates
    sparse = np.array([[-1, -1, -1]])

## ----- Kn or Ks cutoff. ----- ##
# Value Indices
# XY: Kn=6 , Ks=7
# XZ: Kn=8 , Ks=9
# YZ: Kn=10 , Ks=11
if ratio_cutoff:
    rcut_start = datetime.now()
    start_pts = clusters.shape[0]
    clusters = clusters.tolist()
    if ratio_cutoff == 'kn':
        if ratio_by == 'xy':
            #clusters = [c for c in clusters if c[6] != 'NA' and float(c[6])>rcutoff_min and float(c[6])<rcutoff_max]
            indexes = [i for i, c in enumerate(clusters) if c[6] != 'NA' and float(c[6])>rcutoff_min and float(c[6])<rcutoff_max]
            clusters = [clusters[i] for i in indexes]
            downloads = [downloads[i] for i in indexes]
        elif ratio_by == 'xz':
            #clusters = [c for c in clusters if c[8] != 'NA' and float(c[8])>rcutoff_min and float(c[8])<rcutoff_max]
            indexes = [i for i, c in enumerate(clusters) if c[8] != 'NA' and float(c[8])>rcutoff_min and float(c[8])<rcutoff_max]
            clusters = [clusters[i] for i in indexes]
            downloads = [downloads[i] for i in indexes]
        elif ratio_by == 'yz':
            #clusters = [c for c in clusters if c[10] != 'NA' and float(c[10])>rcutoff_min and float(c[10])<rcutoff_max]
            indexes = [i for i, c in enumerate(clusters) if c[10] != 'NA' and float(c[10])>rcutoff_min and float(c[10])<rcutoff_max]
            clusters = [clusters[i] for i in indexes]
            downloads = [downloads[i] for i in indexes]

        elif ratio_by == 'mean':
            # clusters = [c for c in clusters if c[6] != 'NA' and c[8] != 'NA' and c[10] != 'NA' and \
            #             np.mean([float(c[6]), float(c[8]), float(c[10])]) > rcutoff_min and \
            #             np.mean([float(c[6]), float(c[8]), float(c[10])]) < rcutoff_max]
            indexes = [i for i, c in enumerate(clusters) if c[6] != 'NA' and c[8] != 'NA' and c[10] != 'NA' and \
                       np.mean([float(c[6]), float(c[8]), float(c[10])]) > rcutoff_min and \
                       np.mean([float(c[6]), float(c[8]), float(c[10])]) < rcutoff_max]
            clusters = [clusters[i] for i in indexes]
            downloads = [downloads[i] for i in indexes]
        elif ratio_by == 'median':
            # clusters = [c for c in clusters if c[6] != 'NA' and c[8] != 'NA' and c[10] != 'NA' and \
            #             np.median([float(c[6]), float(c[8]), float(c[10])]) > rcutoff_min and \
            #             np.median([float(c[6]), float(c[8]), float(c[10])]) < rcutoff_max]
            indexes = [i for i, c in enumerate(clusters) if c[6] != 'NA' and c[8] != 'NA' and c[10] != 'NA' and \
                       np.median([float(c[6]), float(c[8]), float(c[10])]) > rcutoff_min and \
                       np.median([float(c[6]), float(c[8]), float(c[10])]) < rcutoff_max]
            clusters = [clusters[i] for i in indexes]
            downloads = [downloads[i] for i in indexes]
        else:
            print("ERROR: Unknown ratio cutoff value. Skipping...")
    elif ratio_cutoff == 'ks':
        if ratio_by == 'xy':
            #clusters = [c for c in clusters if c[7] != 'NA' and float(c[7])>rcutoff_min and float(c[7])<rcutoff_max]
            indexes = [i for i, c in enumerate(clusters) if c[7] != 'NA' and float(c[7])>rcutoff_min and float(c[7])<rcutoff_max]
            clusters = [clusters[i] for i in indexes]
            downloads = [downloads[i] for i in indexes]
        elif ratio_by == 'xz':
            #clusters = [c for c in clusters if c[9] != 'NA' and float(c[9])>rcutoff_min and float(c[9])<rcutoff_max]
            indexes = [i for i, c in enumerate(clusters) if c[9] != 'NA' and float(c[9])>rcutoff_min and float(c[9])<rcutoff_max]
            clusters = [clusters[i] for i in indexes]
            downloads = [downloads[i] for i in indexes]
        elif ratio_by == 'yz':
            #clusters = [c for c in clusters if c[11] != 'NA' and float(c[11])>rcutoff_min and float(c[11])<rcutoff_max]
            indexes = [i for i, c in enumerate(clusters) if c[11] != 'NA' and float(c[11])>rcutoff_min and float(c[11])<rcutoff_max]
            clusters = [clusters[i] for i in indexes]
            downloads = [downloads[i] for i in indexes]
        elif ratio_by == 'mean':
            # clusters = [c for c in clusters if c[7] != 'NA' and c[9] != 'NA' and c[11] != 'NA' and \
            #             np.mean([float(c[7]), float(c[9]), float(c[11])]) > rcutoff_min and \
            #             np.mean([float(c[7]), float(c[9]), float(c[11])]) < rcutoff_max]
            indexes = [i for i, c in enumerate(clusters) if c[7] != 'NA' and c[9] != 'NA' and c[11] != 'NA' and \
                       np.mean([float(c[7]), float(c[9]), float(c[11])]) > rcutoff_min and \
                       np.mean([float(c[7]), float(c[9]), float(c[11])]) < rcutoff_max]
            clusters = [clusters[i] for i in indexes]
            downloads = [downloads[i] for i in indexes]
        elif ratio_by == 'median':
            # clusters = [c for c in clusters if c[7] != 'NA' and c[9] != 'NA' and c[11] != 'NA' and \
            #             np.median([float(c[7]), float(c[9]), float(c[11])]) > rcutoff_min and \
            #             np.median([float(c[7]), float(c[9]), float(c[11])]) < rcutoff_max]
            indexes = [i for i, c in enumerate(clusters) if c[7] != 'NA' and c[9] != 'NA' and c[11] != 'NA' and \
                       np.median([float(c[7]), float(c[9]), float(c[11])]) > rcutoff_min and \
                       np.median([float(c[7]), float(c[9]), float(c[11])]) < rcutoff_max]
            clusters = [clusters[i] for i in indexes]
            downloads = [downloads[i] for i in indexes]
        else:
            print("ERROR: Unknown ratio cutoff value. Skipping...")
    elif ratio_cutoff == 'knks':
        if ratio_by == 'xy':
            # clusters = [c for c in clusters if c[6] != 'NA' and c[7] != 'NA' and \
            #             (float(c[6]) / float(c[7])) > rcutoff_min and \
            #             (float(c[6]) / float(c[7])) < rcutoff_max]
            indexes = [c for c in clusters if c[6] != 'NA' and c[7] != 'NA' and \
                       (float(c[6]) / float(c[7])) > rcutoff_min and \
                       (float(c[6]) / float(c[7])) < rcutoff_max]
            clusters = [clusters[i] for i in indexes]
            downloads = [downloads[i] for i in indexes]
        elif ratio_by == 'xz':
            # clusters = [c for c in clusters if c[8] != 'NA' and c[9] != 'NA' and \
            #             (float(c[8]) / float(c[9])) > rcutoff_min and \
            #             (float(c[8]) / float(c[9])) < rcutoff_max]
            indexes = [c for c in clusters if c[8] != 'NA' and c[9] != 'NA' and \
                       (float(c[8]) / float(c[9])) > rcutoff_min and \
                       (float(c[8]) / float(c[9])) < rcutoff_max]
        elif ratio_by == 'yz':
            # clusters = [c for c in clusters if c[10] != 'NA' and c[11] != 'NA' and \
            #             (float(c[10]) / float(c[11])) > rcutoff_min and \
            #             (float(c[10]) / float(c[11])) < rcutoff_max]
            indexes = [c for c in clusters if c[10] != 'NA' and c[11] != 'NA' and \
                       (float(c[10]) / float(c[11])) > rcutoff_min and \
                       (float(c[10]) / float(c[11])) < rcutoff_max]
            clusters = [clusters[i] for i in indexes]
            downloads = [downloads[i] for i in indexes]
        elif ratio_by == 'mean':
            # clusters = [c for c in clusters if c[6] != 'NA' and c[7] != 'NA' and c[8] != 'NA' and \
            #             c[9] != 'NA' and c[10] != 'NA' and c[11] != 'NA' and \
            #             np.mean([float(c[6])/float(c[7]),
            #                      float(c[8])/float(c[9]),
            #                      float(c[10])/float(c[11])]) > rcutoff_min and \
            #             np.mean([float(c[6])/float(c[7]),
            #                      float(c[8])/float(c[9]),
            #                      float(c[10])/float(c[11])]) < rcutoff_max]
            indexes = [c for c in clusters if c[6] != 'NA' and c[7] != 'NA' and c[8] != 'NA' and \
                       c[9] != 'NA' and c[10] != 'NA' and c[11] != 'NA' and \
                       np.mean([float(c[6])/float(c[7]),
                                float(c[8])/float(c[9]),
                                float(c[10])/float(c[11])]) > rcutoff_min and \
                       np.mean([float(c[6])/float(c[7]),
                                float(c[8])/float(c[9]),
                                float(c[10])/float(c[11])]) < rcutoff_max]
            clusters = [clusters[i] for i in indexes]
            downloads = [downloads[i] for i in indexes]
        elif ratio_by == 'median':
            # clusters = [c for c in clusters if c[6] != 'NA' and c[7] != 'NA' and c[8] != 'NA' and \
            #             c[9] != 'NA' and c[10] != 'NA' and c[11] != 'NA' and \
            #             np.median([float(c[6])/float(c[7]),
            #                        float(c[8])/float(c[9]),
            #                        float(c[10])/float(c[11])]) > rcutoff_min and \
            #             np.median([float(c[6])/float(c[7]),
            #                        float(c[8])/float(c[9]),
            #                        float(c[10])/float(c[11])]) < rcutoff_max]
            indexes = [c for c in clusters if c[6] != 'NA' and c[7] != 'NA' and c[8] != 'NA' and \
                       c[9] != 'NA' and c[10] != 'NA' and c[11] != 'NA' and \
                       np.median([float(c[6])/float(c[7]),
                                  float(c[8])/float(c[9]),
                                  float(c[10])/float(c[11])]) > rcutoff_min and \
                       np.median([float(c[6])/float(c[7]),
                                  float(c[8])/float(c[9]),
                                  float(c[10])/float(c[11])]) < rcutoff_max]
            clusters = [clusters[i] for i in indexes]
            downloads = [downloads[i] for i in indexes]
        else:
            print("ERROR: Unknown ratio cutoff value. Skipping...")
    else:
        print("ERROR: Unknown ratio cutoff value. Skipping...")

    clusters = np.array(clusters)
    print("Thinning by %s %s Complete (%s)" % (ratio_by, ratio_cutoff, str(datetime.now()-rcut_start)))
    print("--> Removed Point Count: %s" % str(start_pts - clusters.shape[0]))
    print("--> Remaining Points: %s" % str(clusters.shape[0]))


## ----- Calculate axis data. ----- ##
x_chrs, x_len = buildAxisTicks(c_sort[X_GID], g_size_rel[X_GID], c_rel[X_GID])
y_chrs, y_len = buildAxisTicks(c_sort[Y_GID], g_size_rel[Y_GID], c_rel[Y_GID])
z_chrs, z_len = buildAxisTicks(c_sort[Z_GID], g_size_rel[Z_GID], c_rel[Z_GID])


## ----- Write out graph object(s) ----- ##
# Build graph object.
# graph_obj = {'axes': [[xsp_chname1, ... xsp_chnameN], [xsp_chstart1, ... xsp_chstartN],
#                       [ysp_chname1, ... ysp_chnameN], [ysp_chstart1, ... ysp_chstartN],
#                       [zsp_chname1, ... zsp_chnameN], [zsp_chstart1, ... zsp_chstartN]],
#              'points': [[x1, y1, z1, x1_fid, y1_fid, z1_fid, XY_Kn, XY_Ks, XZ_Kn, XZ_Ks, YZ_Kn, YZ_Ks],
#                         ...,
#                         [xN, yN, zN, xN_fid, yN_fid, zN_fid, XY_Kn, XY_Ks, XZ_Kn, XZ_Ks, YZ_Kn, YZ_Ks]],
#              'x': [gid, species_name, relative_length],
#              'y': [gid, species_name, relative_length],
#              'z': [gid, species_name, relative_length]
#             }
#
graph_obj = {}
graph_obj["points"] = clusters.tolist()
graph_obj["axes"] = [x_chrs["labels"], x_chrs["positions"],
                     y_chrs["labels"], y_chrs["positions"],
                     z_chrs["labels"], z_chrs["positions"]]
graph_obj['x'] = [X_GID, genomes[X_GID]["organism"]["name"], x_len]
graph_obj['y'] = [Y_GID, genomes[Y_GID]["organism"]["name"], y_len]
graph_obj['z'] = [Z_GID, genomes[Z_GID]["organism"]["name"], z_len]

# Build log object.
# log_obj = {'status': <status>,
#            'matches': <match_count>
#            'execution_time': <execution_time>
#           }
log_obj = {}
log_obj["status"] = "complete"
log_obj["matches"] = clusters.shape[0]
log_obj["execution_time"] = str(datetime.now() - start)

# Build output filename.
# Filename follows standard: XGID_YGID_ZGID_<options, '_' separated and in alphabetic order>_<filetype>.json
name_base = "%s_%s_%s_" % (X_GID, Y_GID, Z_GID)
name_ext = []  # deque([])
# Sort
name_ext.append('sort=%s' % sort_method)
# Parse
name_ext.append('parse.len=%s' % str(min_length))
name_ext.append('parse.syn=%s' % str(min_synteny))
# Cluster
if remove_unclustered:
    name_ext.append('cluster.eps=%s' % args.C_eps)
    name_ext.append('cluster.min=%s' % str(min_samples))
# Ratio Thinning
if ratio_cutoff:
    name_ext.append('ratio.by=%s.%s' % (ratio_by, ratio_cutoff))
    name_ext.append('ratio.min=%s' % args.Rmin)
    name_ext.append('ratio.max=%s' % args.Rmax)
# Sort name extensions, write graph & log filenames.
name_ext.sort()
graphName_out = name_base + '_'.join(name_ext) + '_graph.json'
downloadName_out = name_base + '_'.join(name_ext) + '_data.json'
logName_out = name_base + '_'.join(name_ext) + '_log.json'
# Add output path to output files.
graph_out = path.join(args.o, graphName_out)
download_out = path.join(args.o, downloadName_out)
log_out = path.join(args.o, logName_out)
# Write output files.
# graph
dump(graph_obj, open(graph_out, 'wb'))
# downloads
with open(download_out, 'wb') as dout:
    dout.write("#Xsp_gid||Xchr1||Xstart||Xstop||Xname||Xstrand||Xtype||Xdb_feature_id1\tYsp_gid||Ychr||Ystart||Ystop||Yname||Ystrand||Ytype||Ydb_feature_id\tZsp_gid||Zchr||Zstart||Zstop||Zname||Zstrand||Ztype||Zdb_feature_id\tXY_kn||XY_ks||XZ_kn||XZ_ks||YZ_kn||YZ_ks\n")
    for l in downloads:
        dout.write("\t".join(l) + '\n')
# log
dump(log_obj, open(log_out, 'wb'))

## ----- Print concluding message & total execution time. ----- ##
end = datetime.now()
print("Job Complete!")
print("Total Execution Time: %s" % str(end-start))

#Sp1_Data    Sp2_Data    Sp3_data   Mutation_Data
#sp1_gid||chr1||start1||stop1||name1||strand1||type1||db_feature_id1  sp2_gid||chr2||start2||stop2||name2||strand2||type2||db_feature_id2	sp3_gid||chr3||start3||stop3||name3||strand3||type3||db_feature_id3  1-2_kn||1-2_ks||1-3_kn||1-3_ks||2-3_kn||2-3_ks

