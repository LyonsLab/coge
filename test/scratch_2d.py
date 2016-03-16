import json
from sys import argv
from pprint import pprint
from natsort import natsorted

import numpy as np
from sklearn.cluster import DBSCAN

import plotly.plotly as py
from plotly.graph_objs import *

# Length Sort: natsorted(data[key]["chromosomes"], key=lambda chromosome: chromosome["length"])
# Name Sort: natsorted(data[key]["chromosomes"], key=lambda chromosome: chromosome["name"]):


## ----- OPTIONS ----- ##
# Sorting
sort_method = "name"  # "name" or "length"

# Clustering
eps = 0.5
min_samples = 10


## ----- GLOBAL VARIABLES ----- ##
genome_lengths = {}
absolute_lengths = {}
relative_lengths = {}
sorted_chs = {}
genome_scales = {}


## ----- Function Declarations ----- ##
def getPointCoords(syntenic_points):
    gid1 = syntenic_points.keys()[0]
    gid2 = syntenic_points[gid1].keys()[0]
    coords = []  # Each coord is in form [sp1_val, sp2_val, sp1_fid, sp2_fid]

    syntenic_points = syntenic_points[gid1][gid2]

    for sp1chr in syntenic_points:
        # Assign scale & length variables for current species 1 chromosome.
        sp1chr_scale = genome_scales[gid1][sp1chr]
        sp1chr_length = absolute_lengths[gid1][sp1chr]
        sp1chr_relative_length = sp1chr_scale * relative_lengths[gid1]
        for sp2chr in syntenic_points[sp1chr]:
            # Assign scale & length variables for current species 2 chromosome.
            sp2chr_scale = genome_scales[gid2][sp2chr]
            sp2chr_length = absolute_lengths[gid2][sp2chr]
            sp2chr_relative_length = sp2chr_scale * relative_lengths[gid2]

            # Calculate sum of all previous chromosomes
            sp1_sum = 0.0
            for c in sorted_chs[gid1]:
                if c["name"] == sp1chr:
                    break
                else:
                    sp1_sum += genome_scales[gid1][c["name"]] * relative_lengths[gid1]
            sp2_sum = 0.0
            for c in sorted_chs[gid2]:
                if c["name"] == sp2chr:
                    break
                else:
                    sp2_sum += genome_scales[gid2][c["name"]] * relative_lengths[gid2]

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

    coords = np.array(coords)
    return gid1, gid2, coords


## ----- Load data files. ----- ##
file1 = json.load(open(argv[1], 'r'))
file2 = json.load(open(argv[2], 'r'))
file3 = json.load(open(argv[3], 'r'))


## ----- Consolidate genome information. ----- ##
genomes = file1["genomes"].copy()
genomes.update(file2["genomes"])
genomes.update(file3["genomes"])


## ----- Calculate total genome sizes, relative chromosome sizes, and sort. ----- ##
for gid, ginfo in genomes.iteritems():
    genome_scales[gid] = {}
    absolute_lengths[gid] = {}

    # Calculate total genome length.
    total_length = 0
    for chr in ginfo["chromosomes"]:
        total_length += chr["length"]
    # Save genome length.
    genome_lengths[gid] = total_length

    # Sort chromosomes & add to sorted_chromosomes.
    if sort_method == "length":
        sorted_chromosomes = natsorted(ginfo["chromosomes"], key=lambda chromosome: chromosome["length"], reverse=True)
    elif sort_method == "name":
        sorted_chromosomes = natsorted(ginfo["chromosomes"], key=lambda chromosome: chromosome["name"])
    else:
        sorted_chromosomes = ginfo["chromosomes"]
    # Save sorted chromosomes to global variable "sorted_chs"
    sorted_chs[gid] = sorted_chromosomes

    # Calculate chromosome percentages.
    total = 0.0
    for chr in sorted_chromosomes:
        length = chr["length"]
        percent = float(length)/total_length * 100
        total += percent
        genome_scales[gid][chr["name"]] = percent
        absolute_lengths[gid][chr["name"]] = length

    # Sanity Check: Print total percentage (Should be 1, obviously.)
    #print "\t------------------\n\t%f" % total
    #print "\n"


## ----- Calculate relative genome sizes. ----- ##
relative_lengths = genome_lengths.copy()
max_genome = max(genome_lengths, key=lambda k: genome_lengths[k])
for gid,glen in genome_lengths.iteritems():
    relative_lengths[gid] = float(glen)/genome_lengths[max_genome]


## ----- Get Point Coordinates ----- ##
gid1, gid2, coordinates = getPointCoords(file1["syntenic_points"])


## ----- Clustering: DBSCAN ----- ##
pts = coordinates[:,0:2]
db = DBSCAN(eps=eps, min_samples=min_samples).fit(pts)

labels = db.labels_
clusters = coordinates[labels != -1]
sparse = coordinates[labels == -1]

cluster_count = len(set(labels)) - (1 if -1 in labels else 0)
print('Estimated number of clusters: %d' % cluster_count)
print('Total Point Count: %s') % str(coordinates.shape[0])
print('Clustered Point Count: %s') % str(clusters.shape[0])
print('Removed Point Count: %s') % str(sparse.shape[0])


## ----- Plotting: For Confirmation ----- ##
# Build x-axis ticks.
x_chrs = {"positions": [], "labels": []}
x_current = 0
for chr in sorted_chs[gid1]:
    name = chr["name"]
    x_chrs["positions"].append(x_current)
    x_chrs["labels"].append(name)
    x_current += relative_lengths[gid1] * genome_scales[gid1][name]

# Build y-axis ticks.
y_chrs = {"positions": [], "labels": []}
y_current = 0
for chr in sorted_chs[gid2]:
    name = chr["name"]
    y_chrs["positions"].append(y_current)
    y_chrs["labels"].append(name)
    y_current += relative_lengths[gid2] * genome_scales[gid2][name]

# Define data.
kept = Scatter(
    x=clusters[:,0],
    y=clusters[:,1],
    name='Clustered Points',
    mode='markers'
)

removed = Scatter(
    x=sparse[:,0],
    y=sparse[:,1],
    name="Removed Points",
    mode='markers',
    marker=dict(
        symbol='x',
        color='red',
    )
)

# Define layout.
layout = Layout(
    title="Syntenic Dotplot: Clustering Results, EPS=%s, MinSample=%s" % (eps, min_samples),
    margin=dict(
        l=125
    ),
    xaxis=dict(
        title = file1["genomes"][gid1]["organism"]["name"],
        range = [0, x_current],
        tickmode = 'array',
        tickvals = x_chrs["positions"],
        ticktext = x_chrs["labels"]
    ),
    yaxis=dict(
        title = file1["genomes"][gid2]["organism"]["name"],
        range = [0, y_current],
        tickmode = 'array',
        tickvals = y_chrs["positions"],
        ticktext = y_chrs["labels"]
    )
)

# Plot figure.
fig=Figure(data=[kept, removed], layout=layout)
py.plot(fig)
