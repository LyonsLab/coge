import collections
import sys
import os
from cdagline import DagLine


class BlastLine(object):
    __slots__ = ('query', 'subject', 'pctid', 'hitlen', 'nmismatch', 'ngaps', \
                 'qstart', 'qstop', 'sstart', 'sstop', 'evalue', 'score')

    def __init__(self, sline):
        args = sline.split("\t")
        self.query   = args[0]
        self.subject = args[1]

        self.pctid  = float(args[2])
        self.hitlen = int(args[3])

        self.nmismatch = int(args[4])
        self.ngaps     = int(args[5])

        self.qstart = int(args[6])
        self.qstop  = int(args[7])
        self.sstart = int(args[8])
        self.sstop  = int(args[9])

        self.evalue = float(args[10])
        self.score  = float(args[11])

    def __repr__(self):
        return "BlastLine('%s' to '%s', eval=%.3f, score=%.1f)" % (self.query, self.subject, self.eval, self.score)
    def __str__(self):
        return "\t".join(map(str, (getattr(self, attr) for attr in BlastLine.__slots__)))

def adjust_evalue(afile, expected_count=8, evalue_cutoff=5, oclass=DagLine):
    """
    adjust the evalues in a dagchainer/blast file by the number of times they occur.
    query/subjects that appear often will have the evalues raise (made less 
    significant).
    `afile`: path to the blast or dag file.
    `expected_count`: the number of blast hits expect per accn. a good estimate
                      is per file: [number of lines ] / [ number of unique accns ]
                      making the value higher results in less evalue adjustment.    
    `evalue_cutoff`: dont print lines with an adjusted evalue above this
    `oclass`: either DagLine or BlastLine

    yields each dag value after it has been adjusted.
    """
    assert os.path.exists(afile)
    if oclass is BlastLine:
        name1, name2 = ('query', 'subject')
    else:
        name1, name2 = ('a_accn', 'b_accn')


    expected_count = float(expected_count)
    counts = collections.defaultdict(float)
    if isinstance(afile, basestring):
        fh = sys.stdin if afile == "-" else open(afile)
    elif isinstance(afile, list):
        fh = afile

    for line in fh:
        b = oclass(line)
        counts[getattr(b, name1)] += 1.0
        counts[getattr(b, name2)] += 1.0

    for line in open(afile):
        b = oclass(line)
        count = counts[getattr(b, name1)] + counts[getattr(b, name2)]
        b.evalue = b.evalue ** (expected_count / count)

        if b.evalue <= evalue_cutoff: yield b

def main(args):
    import optparse
    p = optparse.OptionParser("""
    adjust the evalues in a dagchainer/blast file by the number of times they occur.
    query/subjects that appear often will have the evalues raise (made less 
    significant). useage:
        python %s some.dag
    """ % sys.argv[0])

    p.add_option('-t', dest='type', help="'blast' or 'dag'.\
       if not specified, the extention the extension of the file is used.",
            default=None)
    p.add_option('-c', dest='expected', type='float', help=
                "the number of blast hits expect per accn. a good estimate\n"
                "is per file: [number of lines ] / [ number of unique accns ]\n"
                "making the value higher results in less evalue adjustment.\n"
                "values between 4 and 15 seem to work well.", default=8)
    p.add_option('-o', dest='out', help="if specified, a file to write to."\
                                     " if not specified, stdout is used", 
                                    default=None)

    p.add_option('-e', dest='evalue', type='float', help="filter (dont print)"
        " lines with a value higher than this", default=5)

    opts, a_file = p.parse_args(args)
    if not a_file: sys.exit(p.print_help())
    a_file = a_file[0]
    if not opts.type:
        opts.type = "blast" if a_file.endswith(".blast") else "dag"
    oclass = DagLine if opts.type == "dag" else BlastLine

    if opts.out is not None:
        out = open(opts.out, "w")
    else: 
        out = sys.stdout 

    for d in adjust_evalue(a_file, opts.expected, opts.evalue, oclass):
        print >>out, str(d)


if __name__ == "__main__":
    main(sys.argv[1:])
