from operator import itemgetter
import numpy
import re
try: import psyco; pysco.full()
except: pass

def dag_array(dagf):
    recs = [] #collections.defaultdict(list)

    fh = open(dagf, 'r')
    qname_len = 0
    sname_len = 0
    qchr_len = 0
    schr_len = 0

    for line in fh:
        if line[0] == '#': continue
        qchr, qname, qstart, qstop, schr, sname, sstart, sstop, score  = line.rstrip("*,\n,+").split("\t")[:9]
        if len(qchr) > qchr_len: qchr_len = len(qchr)
        if len(schr) > schr_len: schr_len = len(schr)
        if len(qname) > qname_len: qname_len = len(qname)
        if len(sname) > sname_len: sname_len = len(sname)
        if not (qname, sname) in recs: recs[(qname, sname)] = []
        recs[(qname, sname)].append([qchr, qname, int(qstart), int(qstop), schr, sname, int(sstart), int(sstop), float(score)])

    fh.close()

    arr = []
    for k in sorted(recs, key=itemgetter(1)):
        arr.extend([li for li in sorted(recs[k], itemgetter(8))])

    dag_names = ('qchr', 'qname', 'qstart', 'qstop', 'schr', 'sname', 'sstart', 'sstop', 'score')
    dag_types = ['S', 'S', 'i4', 'i4', 'S', 'S', 'i4', 'i4', 'f8']
    dag_types[0] += str(qchr_len)
    dag_types[4] += str(schr_len)
    dag_types[1] += str(qname_len)
    dag_types[5] += str(sname_len)
    
    return numpy.rec.array(arr, names=dag_names, formats=dag_types)


chrre = re.compile("(\d+)")
def get_chr(line):
    try:
        return re.search(chrre, line).groups(0)[0]
    except:
        print >>sys.stderr, line
        sys.exit(2)


def blast_to_dag(blast_file, query, subject, qdups, sdups, get_chr=get_chr, condense=True):
    if qdups:
        qdups = frozenset([x.strip() for x in open(qdups)])
    if sdups:
        sdups = frozenset([x.strip() for x in open(sdups)])
    #if query == subject: subject += "2"
    qorg = query   + "_"
    sorg = subject + "_"
    seen = {}
    n_qdups = 0
    n_sdups = 0
    for line in open(blast_file):
        line = line.split("\t")

        if line[0] in qdups: n_qdups += 1; continue
        if line[1] in sdups: n_sdups += 1; continue

        if condense:
            key = line[0] + line[1]
            eval, score = map(float, line[-2:])
            if key in seen and (seen[key][0] < eval and seen[key][1] > score): continue
            seen[key] = (eval, score)

        qinfo = line[0].split("||")
        sinfo = line[1].split("||")
        # it wast just the name


        if len(qinfo) > 1:
            qchr = qinfo[0]
            qlocs = [l.lstrip('0') for l in qinfo[1:3]]
            if len(qinfo) > 4 and qinfo[4] == '-1':
                qlocs.reverse()
        else:
            # a whole chromosome, use the locs it came with.
            qlocs = line[6:8]
            qchr = line[0]
#            qchr = get_chr(line[0])

            line[0] = line[0]+"||"+qlocs[0]+"||"+qlocs[1]

        if len(sinfo) > 1:
            schr = sinfo[0]
            slocs = [l.lstrip('0') for l in sinfo[1:3]]
            if len(sinfo) > 4 and sinfo[4] == '-1':
                slocs.reverse()
        else: 
            # a whole chromosome, use the locs it came with.
            slocs = line[8:10]
            schr = line[1]
#            schr = get_chr(line[1]) 
            line[1] = line[1]+"||"+slocs[0]+"||"+slocs[1]
        print "\t".join([
             qorg + qchr, line[0] + "||" + line[2], qlocs[0], qlocs[1]
            ,sorg + schr, line[1] + "||" + line[2], slocs[0], slocs[1], line[10]])

    if qdups:
        print >>sys.stderr, "removed %i dups from query  " % n_qdups
    if sdups:
        print >>sys.stderr, "removed %i dups from subject" % n_sdups

if __name__ == "__main__":
    import sys, os
    import re
    import cPickle

    from optparse import OptionParser
    usage = """ 
    takes a tab-delimited blast file and converts it to the format used by
    dagchainer and tandems.py. output is to STDOUT.
    if (optional) files are given for query/subject_dups with format:
    dupa_name
    dupb_name
    .
    .
    dupzzz_name

    then any hits containing those are removed. from the output
    """
    parser = OptionParser(usage)

    parser.add_option("-b", "--blast_file", dest="blast_file", help="the name of the blast_file", default=False)
    parser.add_option("-q", "--query",   dest="query",   help="the name of the query organism")
    parser.add_option("-s", "--subject", dest="subject", help="the name of the subject organism")
    parser.add_option("--query_dups", dest="query_dups", help="file containing list of query dups", default=[])
    parser.add_option("--subject_dups", dest="subject_dups", help="file containing list of subject dups", default=[])
    parser.add_option("-c","--condense", dest="condense", help="condense duplicate blast hits", action="store_false")

    (options, _) = parser.parse_args()
    condense=options.condense
    if not options.blast_file:
        sys.exit(parser.print_help())
    blast_to_dag(options.blast_file, options.query, options.subject, options.query_dups, options.subject_dups, condense=condense)
