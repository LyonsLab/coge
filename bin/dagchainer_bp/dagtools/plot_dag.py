import sys
from cdagline import DagLine
import numpy as np
from itertools import cycle
import matplotlib

def main(args):

    import optparse
    p = optparse.OptionParser("plot a dagfile")
    p.add_option('-q', '--qseqid', dest='qseqid', help="seqid of the query")
    p.add_option('-s', '--sseqid', dest='sseqid', help="seqid of the subject")
    p.add_option('-p', '--png', dest='png', help="path of png to save file.")

    p.add_option('-l', '--lines', dest='lines', help="draw as lines (not dots)",
                 action='store_true')

    p.add_option('-d', '--dag', dest='dag', help='path to dag file')
    p.add_option('-b', '--background', dest='back', default=None,
                 help='path to dag file of points to plot in background')
    # TODO outfile.

    opts, _ = p.parse_args(args)
    if not (opts.qseqid and opts.sseqid and (opts.dag or opts.back) and opts.png):
        sys.exit(p.print_help())
    ax = None
    if opts.back:
        ax = plot(opts.back, opts.qseqid, opts.sseqid, opts.png, False, ax, False)

    if opts.dag:
        plot(opts.dag, opts.qseqid, opts.sseqid, opts.png, opts.lines, ax)

def plot(dagfile, qseqid, sseqid, png, lines=False, ax=None, colored=True):

    if ax is None:
        matplotlib.use('Agg')
        from matplotlib import pyplot as plt
        f = plt.figure()
        ax = f.add_subplot(1, 1, 1)
    else:
        from matplotlib import pyplot as plt

    xmax = 0
    ymax = 0

    colors = cycle('rgbcmy')
    if colored:
        alpha = 0.9
        c = 'y'
        s = 2
    else:
        alpha = 0.5
        c = 'k'
        s = 0.3

    pts = []
    for line in open(dagfile):
        if line[0] == '#': 
            if colored and qseqid in line and sseqid in line:
                c = colors.next()
            continue
        dag = DagLine(line)

        if dag.a_seqid != qseqid: continue
        if dag.b_seqid != sseqid: continue

        if lines:
            ax.plot([dag.a_start, dag.a_end], 
                    [dag.b_start, dag.b_end], c=c)
        else:
            pts.append((dag.a_start, dag.b_start, c))

        if dag.a_end > xmax: xmax = dag.a_end
        if dag.b_end > ymax: ymax = dag.b_end

    if not lines:
        pts = np.array(pts, dtype=[('x', int), ('y', int), ('c', 'S1')])
        if pts.shape[0] > 0:
            ax.scatter(pts['x'], pts['y'], edgecolor='none', c=pts['c'], s=s,
                  alpha=alpha)
    ax.axis('tight')
    plt.savefig(png)
    return ax

if __name__ == "__main__":
    if len(sys.argv) < 2: sys.argv.append('zzz')
    main(sys.argv[1:])
