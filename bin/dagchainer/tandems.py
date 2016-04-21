'''
find the local duplicates given a tab-delimited blastfile.
'''
__revision__ = .15

import sys
import os
import re
import collections
try:
    import psyco
    psyco.full()
except ImportError, e:
    print >> sys.stderr, "This will run faster with psyco installed."


def tandems(locs, dupdist=4):
    """ find local dups:
    + locs: a list of integer pairs or a file of (comma separated) integer pairs. 
    + dupdist: the maxiumum distance between accn's for
            them to still be considered local duplicats.
    returns loc_dups 
       dups looks like: 234: [236, 239, 242, 243]
            where the key is the parent duplication, determined by
            order in the genome. and all other dups are within dupdist
            of the nearest dup in that dupset.

    """
    if isinstance(locs, str):
        # so assume it's a file
        locs = [map(int,xy.split(",")) for xy in open(locs)]
    assert isinstance(locs, list)
       
    all_points = [] # save a record of all dups.
    loc_dups = {} # save what we're keeping
    seendups = {} # dont use repeated pairs
    inserted = {} # map q:[s] to s:[q] and vice-versa

    for qloc, sloc in locs:
        if qloc == sloc: continue # self, self
        all_points.append((qloc, sloc))
        if abs(qloc-sloc) > dupdist: continue # too far to be local dups

        (qloc, sloc) = sorted((qloc, sloc))
        if (qloc, sloc) in seendups: continue
        seendups[(qloc, sloc)] = 1

        # qloc has a parent, append sloc to same group
        used = 0
        if qloc in inserted:
            loc_dups[inserted[qloc]].append(sloc)
            inserted[sloc] = inserted[qloc]
            used = 1
        # qloc is a parent, append sloc to it's group
        elif qloc in loc_dups:
            loc_dups[qloc].append(sloc)
            inserted[sloc] = qloc
            used = 1

        # sloc has a parent, append qloc to same group
        if sloc in inserted:
            loc_dups[inserted[sloc]].append(qloc)
            inserted[qloc] = inserted[sloc]
            used = 1
        # sloc is a parent, append qloc to it's group
        elif sloc in loc_dups:
            loc_dups[sloc].append(qloc)
            inserted[qloc] = sloc
            used = 1
        # start a new dupset
        if not used:
            loc_dups[qloc] = [sloc]

       
    # now we have the basic structure, but what if 2 fairly
    # distant dupsets later gained a common relative...
    # must collapse those to a single dupset.

    # reset 
    seendups = {}
    # get rid of repeats

    # use a set to test for intersection.
    all_points = set(all_points)
    lkeys = loc_dups.keys()
    lkeys.sort()

    kidx1, kidx2 = 0, 1
    while kidx1 + kidx2 < len(lkeys):

        k1 = lkeys[kidx1]
        vals1 = set(loc_dups[k1])

        # so any other dupset within dupdist of this one could 
        # potentially have a common dup, so look through all of them.
        while kidx2 < dupdist and kidx1 + kidx2 < len(lkeys):
            k2 = lkeys[kidx1+kidx2]

            vals2 = loc_dups[k2]
            
            # they dont have any in common.
            if not len(vals1.intersection(vals2)):
                kidx2 +=1
                continue

            # can do destructive inside loop because we delete k1 not k2
            # and we're using a separate copy of the keys.
            del loc_dups[k1]
            loc_dups[k2].extend(vals1)
            break

        kidx1 += 1
        kidx2 = 1
    
    #return loc_dups
    # remove repeats and parents and sort.
    # use new dups to make sure lowest number accn is the parent.
    new_dups = collections.defaultdict(set)
    for lkey, dups in loc_dups.iteritems():
        dups = sorted(list(set(d for d in dups if d != lkey)) + [lkey])
        new_dups[dups[0]].update(dups[1:])

    return new_dups
         

def get_self_hits(hit_file, order):

    self_hits = {}
    seen = {}
    for line in open(hit_file):
        line = line.split("\t")

        key = line[1] + line[5]
        if key in seen: continue
        seen[key] = 1

        qchr   = line[0]
        schr   = line[4]
        if qchr != schr: continue
        if not qchr in self_hits: self_hits[qchr] = []
        self_hits[qchr].append((order[line[1]], order[line[5]]))

    return self_hits


def order_from_hits(hit_file):
    
    seen = {}
    list_order = []

    bcmp = cmp # faster lookup.
    def _cmp(a, b):
        """tuples of chr, loc, name"""
        if a[0] != b[0]: return bcmp(a[0], b[0])
        return a[1] - b[1]


    for line in open(hit_file):
        line = line.split("\t")

        key = line[1] + line[5]
        if key in seen: continue
        seen[key] = 1

        qchr, qname = line[:2]
        qstart, qstop = map(int, line[2:4])
        schr, sname = line[4:6]
        sstart, sstop = map(int, line[6:8])

        qloc = (qstart + qstop)/2
        sloc = (sstart + sstop)/2
        list_order.append((qchr, qloc, qname))
        list_order.append((schr, sloc, sname))

    list_order = sorted(set(list_order), cmp=_cmp)
    return dict((name, i) for i, (qchr, loc, name) in enumerate(list_order))
        
def clean(non_parent_dups, hit_file, is_same=True):
    """print out a copy of the original hit file excluding 
    any line containing a local duplication"""
    clean_file = open(hit_file + ".nodups", "w")
#    print >>sys.stderr, "creating:" +  hit_file + ".nodups"
    dup_set = frozenset(non_parent_dups)

    seen = {}
    for sline in open(hit_file):
        line = sline.split("\t")
        if line[1] in dup_set or line[5] in dup_set: continue
        key = line[1] + line[5]
        if key in seen: continue
        seen[key] = 1
        key = line[5] + line[1]
        if key in seen: continue
        seen[key] = 1
        if is_same:
            if line[0] > line[4]:
                line[0:4], line[4:8] = line[4:8], line[0:4]
                #if int(line[1]) > int(line[2]):
                #    line[1], line[2] = line[2], line[1]
                #    line[5], line[6] = line[6], line[5]
            clean_file.write("\t".join(line))
        else:
            clean_file.write(sline)
    clean_file.close()


def main(dup_dist, hit_file, save_dups=False, clean_file=False):

    if save_dups:
        save_dups = open(hit_file + ".dups",'w')

    order = order_from_hits(hit_file)
    self_hits = get_self_hits(hit_file, order)
    invorder = dict((i, a) for a, i in order.iteritems())

    non_parent_dups = []

    for ichr, locs in self_hits.iteritems():
        locs = tandems(locs, dup_dist)
        for parent, dups in locs.iteritems():
            dup_names = [invorder[d] for d in  dups]
            parent = invorder[parent]
            if save_dups:
                print >>save_dups, "%s,%s" % (parent, ",".join(dup_names))
            non_parent_dups.extend(dup_names)
    non_parent_dups.sort()
    print "\n".join(non_parent_dups)

    # print out the original file with any hit containing a dup excluded.
    if clean_file:
        clean(non_parent_dups, hit_file)

if __name__ == "__main__":

    from optparse import OptionParser
    hit_file_help = """file with columns of:
        query_chr [tab] query_name [tab] 5_prime_query [tab] 3_prime_query [tab] subject_chr [tab] subject_name [tab] 5_prime_subject [tab] 3_prime_subject [tab] evalue (optional)"""
    dup_dist_help = "the minimum distance between 2 genes on the same chromosome which will prevent them from being called as a pair of local dups.  if they are farther apart than `dup_dist`, but share a common hit within `dup_dist` of both, then they will still be in the same set of local dups"
    save_dups_help = """if this flag is seen, a new file named HIT_FILE.dups will be created with the full dups data: the parent as the first column and dups as comma-delimited columns continuing to the left e.g.:
                        parent_dupa, dupa1, dupa2\\nparent_dupb, dupb1, dupb2, dupb3\\nparentdupc, dupc1\\n..."""
    clean_help = "if this flag is seen, then a new file will be created by the name of HIT_FILE.nodups that prints out the HITFILE contents excluding any hit that contains a non-parent local dup"

    usage = """
    python %s -d 5 -i /var/blast/rice.hits -s -r > /var/blast/rice_dups.txt

    the output to STDOUT is all dups (columns) except the parent (first). 
    this is useful to capture for later filtering of all but the parent duplication which is NOT included in the list

    the dups file will be named HIT_FILE.dups with format:
       parenta, dupa1, dupa2
       parentb, dupb1, dupb2, dupb3
       parentc, dupc1
       ...

    = Assumptions =
    the parent hit is the first gene by basepair location.
    Any e-value filtering is done before sending to this program.
    multiple hits between the same query and subject provide no additional information, only the first is used.

    = TODO =
    the order is not always with the lowest numbered gene first, if it
    gets merged to another dupset...
    """ % (sys.argv[0],)
    parser = OptionParser(usage)
    parser.add_option("-d", "--dup_dist",  dest="dup_dist",   help=dup_dist_help, type="int", default=6)
    parser.add_option("-i", "--hit_file",  dest="hit_file",   help=hit_file_help)
    parser.add_option("-s",                dest="save_dups",  help=save_dups_help, action="store_true")
    parser.add_option("-r",                dest="clean_file", help=clean_help,     action="store_true")

    

    (options, _) = parser.parse_args()
    if not options.hit_file:
        sys.exit(parser.print_help())
    main(options.dup_dist, options.hit_file, options.save_dups, options.clean_file)
