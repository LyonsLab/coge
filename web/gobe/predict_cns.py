
plot = 0
if plot:
    os.environ['HOME'] = '/tmp/'
    import matplotlib
    matplotlib.use('Agg')
    import pylab

import os
from blast_misc import blast_array
import sys
import sqlite3
from tables import numexpr
from scipy.stats import linregress
import re
import commands


def find_colinear_hits(blastfile, qeval, seval, mask='query', as_str=False):
    sqlite_file = blastfile[:blastfile.rfind(".")] + ".sqlite"
    db = sqlite3.connect(sqlite_file)
    cur = db.cursor()
    
    # need these to convert the absolute coords in sqlite to match
    # the local positions in the fresh blast
    qmin, smin = [x[0] for x in cur.execute('SELECT bpmin FROM image_info ORDER BY id')]

    # so we mask the new blast with anything that's NOT an HSP
    sql = "SELECT image_id, bpmin, bpmax FROM image_data WHERE type != 'HSP'"
    if mask == 'query': sql += ' AND image_id = 1'

    b = blast_array(blastfile, dopickle=False, best_hit=0, maxkeep=99999)
    if plot:
        pylab.plot(b['qstart'], b['sstart'], "kx")
    b = b[(b['eval'] < qeval) & (b['eval'] < seval)]

    if plot:
        pylab.plot(b['qstart'], b['sstart'], "ro")

    # TODO: remove stuff that's way off the diagonal?


    for row  in cur.execute(sql).fetchall():
        (start, stop, lmin) = row[0] == 1 and ('qstart', 'qstop', qmin) or ('sstart', 'sstop', smin)
        assert (start, stop) == ('qstart', 'qstop') # for now, always only using query
        cds_start, cds_stop = row[1] - lmin + 1, row[2] - lmin + 1
        bstart, bstop  = b[start], b[stop]
        b = b[numexpr.evaluate("(((bstart < cds_start) & (bstop < cds_start)) | ((bstop  > cds_stop ) & (bstart > cds_stop)))")]
                               
    r = 0 
    if not b.shape[0]: return None
    delta = 0.2 * b['sstart'].max()

    # here, try to find a sort of line, and keep removing outliers to only get linear cnss
    for i in range(4):
        slope, intercept, r, zy, zz = linregress(b['qstart'], b['sstart'])
        #print >>sys.stderr, slope, intercept, r, zy, zz
        if r > 0.8: break
        bqstart = b['qstart']
        expected = numexpr.evaluate('intercept + slope * bqstart')
        bsstart = b['sstart']
        s = b.shape[0]
        b = b[numexpr.evaluate('bsstart - expected < delta')]
        if s == b.shape[0]: break # not removing anything.

    if plot:
        pylab.plot(b['qstart'], b['sstart'], "bo")
        pylab.savefig('/var/www/ms_tmp/d.png')
    
    cnss = []
    start_stops = [map(lambda p: int(p) + qmin - 1, pair) for pair in zip(b['qstart'], b['qstop'])]

    for qstart, qstop in start_stops:
        qres =  cur.execute('SELECT xmin, ymin, xmax, ymax, id, pair_id FROM image_data WHERE image_id = 1 AND bpmin = ? AND bpmax = ?', (qstart, qstop)).fetchone()
        this_cns = [qres[:-1]]
        if not qres: continue
        sres =  cur.execute('SELECT xmin, ymin, xmax, ymax, id FROM image_data WHERE id = ?', (qres[-1],)).fetchone()
        this_cns.append(sres)
        cnss.append(this_cns)
    return cnss

if __name__ == "__main__":
    blastfile = sys.argv[1]
    qeval = float(sys.argv[2])
    seval = float(sys.argv[3])

    colinear_starts_stops = find_colinear_hits(blastfile, qeval, seval, as_str=True)
    #print >>sys.stderr, colinear_starts_stops
    print "\n".join(colinear_starts_stops)

