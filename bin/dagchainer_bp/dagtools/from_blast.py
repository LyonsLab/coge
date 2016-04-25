import sys
import operator
sys.path.insert(0, "/opt/src/flatfeature")
from flatfeature import Flat

def to_order(flat, outfile):
    seqids = sorted(flat.seqids)
    order = {}
    i = 0
    out = open(outfile, 'w')
    for seqid in seqids:
        feats = sorted(flat[flat['seqid'] ==seqid], key=operator.itemgetter('start'))
        for feat in feats:
            order[feat["accn"]] = i
            print >>out, feat['accn']
            i += 1
    out.close()
    return order

def blast_to_dag(blast_file, query, subject, qflat_file, sflat_file, qdups,
                 sdups, order, condense):
    """
    if order is true, we return the index of the gene, instead of the basepair position
    """
    if qdups:
        qdups = frozenset(x.strip() for x in open(qdups))
    if sdups:
        sdups = frozenset(x.strip() for x in open(sdups))

    qflat = Flat(qflat_file)
    sflat = Flat(sflat_file)

    qflat_d = {}
    for q in qflat:
        qflat_d[q["accn"]] = q
    sflat_d = {}
    for s in sflat:
        sflat_d[s["accn"]] = s

    if order:
        qorder = to_order(qflat, qflat_file + ".order")
        sorder = to_order(sflat, sflat_file + ".order")

    qorg = query   + "_"
    sorg = subject + "_"
    seen = {}
    n_qdups = 0
    n_sdups = 0
    if condense:
        for line in open(blast_file):
            line = line.split("\t")
            qname, sname = line[:2]
            key = qname + "@" + sname
            v = float(line[-2])
            if not key in seen or v < seen[key]: 
                seen[key] = v

    for line in open(blast_file):
        line = line.split("\t")
        qname, sname = line[:2]

        if qdups is not None and qname in qdups: n_qdups += 1; continue
        if sdups is not None and sname in sdups: n_sdups += 1; continue
        qfeat = qflat_d[qname]
        sfeat = sflat_d[sname]

        if condense:
            key = qname + "@" + sname
            v = float(line[-2])
            if seen[key] > v: continue
        
        if order:
            qo = qorder[qname]
            so = sorder[sname]
            print "\t".join(map(str, [
                 qorg + qfeat['seqid'], qname, qo, qo,
                 sorg + sfeat['seqid'], sname, so, so, line[-2]]))

        else:
            print "\t".join(map(str, [
                 qorg + qfeat['seqid'], qname, qfeat['start'], qfeat['end']
                ,sorg + sfeat['seqid'], sname, sfeat['start'], sfeat['end'], line[-2]]))

    if qdups:
        print >>sys.stderr, "removed %i dups from query  " % n_qdups
        if n_qdups == 0:
            raw_input("didnt remove any dups from teh query, this is Baad!!!...  press a key")
    if sdups:
        print >>sys.stderr, "removed %i dups from subject" % n_sdups
        if n_sdups == 0:
            raw_input("didnt remove any dups from the subject, this is Baad!!!...  press a key")

def genomic_blast_to_dag(blast_file, merge_at=10000):

    for line in open(blast_file):
        bline = line.split("\t")

def main(args):
    from optparse import OptionParser
    usage = """ 
    takes a tab-delimited blast file and converts it to the format used by
    dagchainer and tandems.py. output is to STDOUT.
    if (optional) files are given for query/subject_dups with format:
    dupa_name
    dupb_name

    .
    dupzzz_name

    then any hits containing those are removed. from the output
    """
    parser = OptionParser(usage)

    parser.add_option("-b", "--blast_file", dest="blast_file", 
                      help="the name of the blast_file", default=False)
    parser.add_option("-q", "--query",   dest="query",   
                      help="the name of the query organism")
    parser.add_option("-s", "--subject", dest="subject", 
                      help="the name of the subject organism")
    parser.add_option("-o", dest="order", action="store_true", default=False, 
                      help="if this is true, the relative positions are sent "
                      "to the dag file, rather than the basepair positions.")
    parser.add_option("--qflat",   dest="qflat", 
                      help="the path of the query flat")
    parser.add_option("--sflat", dest="sflat", 
                      help="the path of the subject subject")
    parser.add_option("--query_dups", dest="query_dups", 
                  help="file containing list of query dups", default=None)
    parser.add_option("--subject_dups", dest="subject_dups", 
                  help="file containing list of subject dups", default=None)
    parser.add_option("-c","--condense", dest="condense", 
                  help="condense duplicate blast hits", action="store_true")

    (options, _) = parser.parse_args()
    if not options.blast_file:
        sys.exit(parser.print_help())
    blast_to_dag(options.blast_file, options.query, options.subject,
                 options.qflat, options.sflat, 
                 options.query_dups, options.subject_dups, 
                 options.order, options.condense)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(sys.argv[1:])
    else:
        main([])

