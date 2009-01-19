"""
the diags_file contains the hits that dagchainer found, the all_file is the full list of blast hits.
this script goes through all of the hits in the dag file and adds any hit from the all_file that is within 'dist' of 
any hit in the diag.
"""


from rtree import Rtree
import sys
import os
import psyco; psyco.full()

"""
6    50
rice_1  1||36212110||36215209||OS01G61990||-1||CDS  36215209    36212110    sorghum_1   1||62493571||62496164||SB01G039040||-1||CDS 62496164    62493571    4.000000e-18    52
rice_1  1||36223907||36227093||OS01G62020||1||CDS   36223907    36227093    sorghum_1   1||62517465||62520208||SB01G039050||1||CDS  62517465    62520208    1.000000e-49    95
rice_1  1||36239128||36239914||OS01G62060||1||CDS   36239128    36239914    sorghum_1   1||62554416||62554908||SB01G039100||-1||CDS 62554908    62554416    2.000000e-10    96
rice_1  1||36293042||36293695||OS01G62130||-1||CDS  36293695    36293042    sorghum_1   1||62652716||62653309||SB01G039190||1||CDS  62652716    62653309    5.000000e-20    88
rice_1  1||36341022||36341441||OS01G62230||-1||CDS  36341441    36341022    sorghum_1   1||62699344||62699790||SB01G039260||-1||CDS 62699790    62699344    2.000000e-98    126
rice_1  1||36366694||36369819||OS01G62290||1||CDS   36366694    36369819    sorghum_1   1||62807491||62810206||SB01G039390||1||CDS  62807491    62810206    1.000000e-250   146
"""

def read_dag_to_tree(all_hits):
    """create an rtree, using query as x, subject as y
    do this for all htis, then for each diag, do an intersection
    (+bbuffer) to find nearby 
    """
    gdxs = {}
    for i, sline in enumerate(open(all_hits)):
        if sline[0] == '#': continue
        line = sline[:-1].split("\t")
        chrs = tuple(sorted([line[0], line[4]]))
        # so save the index, which will return i when queried
        # and associate i with the text line vie the dict.
        if not chrs in gdxs: gdxs[chrs] = ({}, Rtree())
        q0, q1 = sorted(map(int, line[2:4]))
        s0, s1 = sorted(map(int, line[6:8]))
        assert q0 < q1 and s0 < s1
        gdxs[chrs][1].add(i, (q0, s0, q1, s1))

        gdxs[chrs][0][i] = sline
        
    return gdxs
        



def main(dist, diags, all_hits):
    """empty docstring"""

    gdxs = read_dag_to_tree(all_hits)
    seen = {}
    for sline in open(diags):
        # reset seen for each new diagonal...
        if sline[0] == '#': seen = {}; print sline.strip(); continue
        line = sline[:-1].split("\t")
        chrs = (line[0], line[4])
        info_lines, tree = gdxs[chrs]

        q0, q1 = sorted(map(int, line[2:4]))
        s0, s1 = sorted(map(int, line[6:8]))
        assert q0 < q1 and s0 < s1

        idxs = tree.intersection((q0 - dist, s0 - dist, q1 + dist, s1 + dist))

        seen[(line[1], line[5])] = True
        print sline,
        for i in idxs:
            iline = info_lines[i]
            ikey = (iline[1], iline[5])
            if ikey in seen: continue
            seen[ikey] = True
            print iline,
        







if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("-d", "--dist", dest="dist", help="distance around each pair to look for missed pairs")
    parser.add_option("--diags", dest="diags", help="the dag aligncoords file (something like q_s.dag.nodups.filter.aligncoords")
    parser.add_option("--all", dest="all", help="the dag blast hit file containing all hits of q to s (something like q_s.dag.nodups")

    (options, _) = parser.parse_args()
    if not (options.dist and options.diags and options.all):
        sys.exit(parser.print_help())

    main(int(options.dist), options.diags, options.all)
