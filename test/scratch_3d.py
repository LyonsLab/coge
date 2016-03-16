import numpy as np

from collections import defaultdict
from datetime import datetime
from json import load
from natsort import natsorted
from plotly.plotly import plot
from plotly.graph_objs import Scatter3d, Layout, Figure, Scene
from sklearn.cluster import DBSCAN
from sys import argv

from pprint import pprint

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                          ____        __  _                                                          #
#                                         / __ \____  / /_(_)___  ____  _____                                         #
#                                        / / / / __ \/ __/ / __ \/ __ \/ ___/                                         #
#                                       / /_/ / /_/ / /_/ / /_/ / / / (__  )                                          #
#                                       \____/ .___/\__/_/\____/_/ /_/____/                                           #
#                                           /_/                                                                       #
#                                                                                                                     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## ----- Sorting ----- ##
sort_method = "name"  # "name" or "length"

## ----- Clustering ----- ##
remove_unclustered = True
eps = 0.5
min_samples = 10

## ----- Parsing ----- ##
min_length = 100000000
# no_synteny = False



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                   ______                 __  _                                                      #
#                                  / ____/_  ______  _____/ /_(_)___  ____  _____                                     #
#                                 / /_  / / / / __ \/ ___/ __/ / __ \/ __ \/ ___/                                     #
#                                / __/ / /_/ / / / / /__/ /_/ / /_/ / / / (__  )                                      #
#                               /_/    \__,_/_/ /_/\___/\__/_/\____/_/ /_/____/                                       #
#                                                                                                                     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def minimizeGenomes(genomes_object, cutoff):
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
    gid1 = sp_object.keys()[0]
    gid2 = sp_object[gid1].keys()[0]

    sp1_chlist = [ch["name"] for ch in genomes_object[gid1]["chromosomes"]]
    sp2_chlist = [ch["name"] for ch in genomes_object[gid2]["chromosomes"]]

    sp_object[gid1][gid2] = {sp1ch: sp_object[gid1][gid2][sp1ch] for sp1ch in sp_object[gid1][gid2] if sp1ch in sp1_chlist}
    for sp1ch in sp_object[gid1][gid2]:
        sp_object[gid1][gid2][sp1ch] = {sp2ch: sp_object[gid1][gid2][sp1ch][sp2ch] for sp2ch in sp_object[gid1][gid2][sp1ch] if sp2ch in sp2_chlist}

    #print sp_object[gid1][gid2].keys()

    #exit()
    #for sp1_ch in sp_object[gid1][gid2]:

    #exit()
    return sp_object


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
    :return coords: coordinates np-array (n, 4) [[loc1_gid1, loc1_gid2, fid1_gid1, fid1_gid2],
                                                 ...,
                                                 [locn_gid1, locn_gid2, fidn_gid1, fidn_gid2]]
    """
    coord_start = datetime.now()

    # Assign order of gids.
    gid1 = syntenic_points.keys()[0]
    gid2 = syntenic_points[gid1].keys()[0]
    # Variable to store all coordinates.
    coords = []  # Each coord is in form [sp1_val, sp2_val, sp1_fid, sp2_fid]

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

                coord = [sp1_val, sp2_val, sp1_fid, sp2_fid]
                coords.append(coord)

    ids = (gid1, gid2)
    coords = np.array(coords)

    coord_end = datetime.now()
    print("Point Extraction Complete (%s)" % str(coord_end-coord_start))
    print("--> %s vs %s" % (gid1, gid2))
    return ids, coords


def reorganizeArray(id_A, id_B, gids_tuple, unordered_coordinates ):
    """
    From a tuple of genome ID's representing a coordinate array's structure (gidI, gidJ).
    and the corresponding coordinate np-array (n x 4) of design [[locI1, locJ1, fidI1, fidJ2],
                                                                  ...
                                                                 [locIn, locJn, fidIn, fidJn]]
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
        ordered_coordinates = unordered_coordinates[:,[1,0,3,2]]
    else:
        ordered_coordinates = None
        print "ERROR: reorganizeArray (%s, %s)" % (str(id_A), str(id_B))
        exit()

    reorg_end = datetime.now()
    print("Reorganization Complete (%s)" % str(reorg_end-reorg_start))
    print("--> (%s,%s) to (%s,%s)" % (gids_tuple[0], gids_tuple[1], id_A, id_B))

    return ordered_coordinates


def mergePoints(setXY, setXZ, setYZ):
    """
    Merge XY, XZ, YZ points.

    :param setXY: XY-points object ( (n,4) np-array [[x_loc, y_loc, x_fid, y_fid]...] ).
    :param setXZ: XZ-points object ( (n,4) np-array [[x_loc, z_loc, x_fid, z_fid]...] ).
    :param setYZ: YZ-points object ( (n,4) np-array [[y_loc, z_loc, y_fid, z_fid]...] ).
    :return: common hits ( (n,6) np-array [[x_loc, y_loc, z_loc, x_fid, y_fid, z_fid]...] ).
    """
    merge_start = datetime.now()

    hits = []
    errors = []

    # Build indexed XY sets.
    setXY_xindexed = defaultdict(list)
    setXY_yindexed = defaultdict(list)
    for s in setXY:
        setXY_xindexed[s[2]].append(s)
        setXY_yindexed[s[3]].append(s)

    # Build indexed XZ sets.
    setXZ_xindexed = defaultdict(list)
    setXZ_zindexed = defaultdict(list)
    for s in setXZ:
        setXZ_xindexed[s[2]].append(s)
        setXZ_zindexed[s[3]].append(s)

    # Build indexed YZ sets.
    setYZ_yindexed = defaultdict(list)
    setYZ_zindexed = defaultdict(list)
    for s in setYZ:
        setYZ_yindexed[s[2]].append(s)
        setYZ_zindexed[s[3]].append(s)

    # Find hits.
    for XY in setXY:
        XZ_X_matches = setXZ_xindexed[XY[2]]
        YZ_Y_matches = setYZ_yindexed[XY[3]]
        if len(XZ_X_matches) > 0 and len(YZ_Y_matches) > 0:
            for XZ in XZ_X_matches:
                for YZ in YZ_Y_matches:
                    if XY[2] == XZ[2] and XY[3] == YZ[2] and XZ[3] == YZ[3]:
                        if XY[0] == XZ[0] and XY[1] == YZ[0] and YZ[1] == XZ[1]:
                            hit = [XY[0], YZ[0], XZ[1], XY[2], YZ[2], XZ[2]]
                            hits.append(hit)
                        else:
                            err = [XY, XZ, YZ]
                            errors.append(err)

    hits = np.array(hits)
    errors = np.array(errors)
    merge_end = datetime.now()

    print "Merge Complete (%s)" % str(merge_end-merge_start)
    print "--> Hits Found: %s" % str(hits.shape[0])
    print "--> Errors Recorded: %s" % str(errors.shape[0])
    return hits


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


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                 __  ___      _          _____           _       __                                  #
#                                /  |/  /___ _(_)___     / ___/__________(_)___  / /_                                 #
#                               / /|_/ / __ `/ / __ \    \__ \/ ___/ ___/ / __ \/ __/                                 #
#                              / /  / / /_/ / / / / /   ___/ / /__/ /  / / /_/ / /_                                   #
#                             /_/  /_/\__,_/_/_/ /_/   /____/\___/_/  /_/ .___/\__/                                   #
#                                                                      /_/                                            #
#                                                                                                                     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
start = datetime.now()

## ----- Check Command Line Arguments ----- ##
if len(argv) != 7:
    print "ERROR: Command Line Arguments"
    print("USAGE: python scratch.py <x_gid> <y_gid> <z_gid> <file1> <file2> <file3>")
    print("HINT: 'file's are the pairwise comparisons of genomes generated with dotplot_dots.py, in any order")
    exit()

## ----- Assign GIDs. ----- ##
X_GID = argv[1]
Y_GID = argv[2]
Z_GID = argv[3]

## ----- Load data files. ----- ##
file1 = load(open(argv[4], 'r'))
file2 = load(open(argv[5], 'r'))
file3 = load(open(argv[6], 'r'))


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

## ----- Calculate total genome sizes, relative chromosome sizes, and sort. ----- ##
g_size_abs, g_size_rel, c_sort, c_rel, c_abs = scaleSortGenomes(genomes)

## ----- Get point coordinates from each file. ----- ##
gids1, coordinates1 = getPointCoords(file1_sp, c_rel, c_abs, g_size_rel, c_sort)
gids2, coordinates2 = getPointCoords(file2_sp, c_rel, c_abs, g_size_rel, c_sort)
gids3, coordinates3 = getPointCoords(file3_sp, c_rel, c_abs, g_size_rel, c_sort)

## ----- Determine order of arrays & restructure to represent correct XYZ ordering. ----- ##
# XY order.
if X_GID in gids1 and Y_GID in gids1:
    coordinatesXY = reorganizeArray(X_GID, Y_GID, gids1, coordinates1)
elif X_GID in gids2 and Y_GID in gids2:
    coordinatesXY = reorganizeArray(X_GID, Y_GID, gids2, coordinates2)
elif X_GID in gids3 and Y_GID in gids3:
    coordinatesXY = reorganizeArray(X_GID, Y_GID, gids3, coordinates3)
else:
    coordinatesXY = None
    print "ERROR: XY coordinates unidentified."
    exit()
# XZ order.
if X_GID in gids1 and Z_GID in gids1:
    coordinatesXZ = reorganizeArray(X_GID, Z_GID, gids1, coordinates1)
elif X_GID in gids2 and Z_GID in gids2:
    coordinatesXZ = reorganizeArray(X_GID, Z_GID, gids2, coordinates2)
elif X_GID in gids3 and Z_GID in gids3:
    coordinatesXZ = reorganizeArray(X_GID, Z_GID, gids3, coordinates3)
else:
    coordinatesXZ = None
    print "ERROR: XZ coordinates unidentified."
    exit()
# YZ order.
if Y_GID in gids1 and Z_GID in gids1:
    coordinatesYZ = reorganizeArray(Y_GID, Z_GID, gids1, coordinates1)
elif Y_GID in gids2 and Z_GID in gids2:
    coordinatesYZ = reorganizeArray(Y_GID, Z_GID, gids2, coordinates2)
elif Y_GID in gids3 and Z_GID in gids3:
    coordinatesYZ = reorganizeArray(Y_GID, Z_GID, gids3, coordinates3)
else:
    coordinatesYZ = None
    print "ERROR: YZ coordinates unidentified."
    exit()


## ----- Merge Coordinates ----- ##
coordinates = mergePoints(coordinatesXY, coordinatesXZ, coordinatesYZ)

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

## ----- Plotting: For Confirmation ----- ##
# Build axis ticks.
x_chrs, x_len = buildAxisTicks(c_sort[X_GID], g_size_rel[X_GID], c_rel[X_GID])
y_chrs, y_len = buildAxisTicks(c_sort[Y_GID], g_size_rel[Y_GID], c_rel[Y_GID])
z_chrs, z_len = buildAxisTicks(c_sort[Z_GID], g_size_rel[Z_GID], c_rel[Z_GID])

# Define data.
kept = Scatter3d(
    x=clusters[:,0],
    y=clusters[:,1],
    z=clusters[:,2],
    name='Clustered Points',
    mode='markers'
)
removed = Scatter3d(
    x=sparse[:,0],
    y=sparse[:,1],
    z=sparse[:,2],
    name="Removed Points",
    mode='markers',
    marker=dict(
        color='red',
    )
)

# Define layout.
layout = Layout(
    title="Syntenic Dotplot: Clustering Results, EPS=%s, MinSample=%s" % (eps, min_samples),
    scene=Scene(
        xaxis=dict(
            title = genomes[X_GID]["organism"]["name"],
            range = [0, x_len],
            tickmode = 'array',
            tickvals = x_chrs["positions"],
            ticktext = x_chrs["labels"]
        ),
        yaxis=dict(
            title = genomes[Y_GID]["organism"]["name"],
            range = [0, y_len],
            tickmode = 'array',
            tickvals = y_chrs["positions"],
            ticktext = y_chrs["labels"]
        ),
        zaxis=dict(
            title = genomes[Z_GID]["organism"]["name"],
            range = [0, z_len],
            tickmode = 'array',
            tickvals = z_chrs["positions"],
            ticktext = z_chrs["labels"]
        )
    )
)

# Plot figure.
fig=Figure(data=[kept, removed], layout=layout)
plot(fig)

## ----- Print concluding message & total execution time. ----- ##
end = datetime.now()
print("Job Complete!")
print("Total Execution Time: %s" % str(end-start))