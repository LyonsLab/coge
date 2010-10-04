#!/usr/bin/python

from pyamf.remoting.gateway.wsgi import WSGIGateway
import pyamf.amf3 # for some reason, have to import this...
import sys, os
import sqlite3

sys.path.extend(['/opt/apache/CoGe/gobe/'])
import predict_cns

tmpdbs = ['/opt/apache/CoGe/gobe/tmp/', '/opt/apache/CoGe/tmp/GEvo/']

dbpath = "/opt/apache/CoGe/data/sqlite/pair_tracking.db"
tracking_db = sqlite3.connect(dbpath)
tracking_db.row_factory = sqlite3.Row
tcur = tracking_db.cursor()

def bagwrap(fn):
    def newfn(*args, **kwargs):
        try:
            res = fn(*args, **kwargs)
            return res
            #if not res: return 
            #if isinstance(res, dict): return pyamf.Bag(res)
            #return res
        except Exception, e:
            #import traceback
            #rint >>sys.stderr, traceback.extract_stack()
            print >>sys.stderr, "function:" + fn.func_name
            print >>sys.stderr, "*args:" + ",".join(map(str, list(args)))
            print >>sys.stderr, sys.exc_info()

            tracking_db.rollback()
            return str(e) + str(sys.exc_info()) 
    return newfn



def get_temp_db(base_name, ext=".sqlite"):
    for t in tmpdbs:
        path = t + base_name + ext
        if os.path.exists(path):
            if ext == ".sqlite":
                print >>sys.stderr, path
                return sqlite3.connect(path)
            else:
                return path

@bagwrap
def save(*args, **kwargs):
    """called when user clicks [save] in the flash annotation swf. saves changes to db.
    args[0] looks like {'keywords':[0,2, ...], 'annos':[1,2,...], 'notes':'asdf'}
    """
    data = args[0]
    tcur.execute("UPDATE genespace SET revisit = ?, qdups = ?, sdups = ?, annotation = ?, keywords = ?, notes= ? WHERE genespace_id = ?"
                ,(  data['revisit']
                   , data['qdups']
                   , data['sdups']
                   , "|".join(map(str,data['annos']))
                   , "|".join(map(str,data['keywords']))
                   , data['notes']
                   , data['genespace_id']
                   ))
    tracking_db.commit()
    if not os.path.exists(tmp_db):
        tmp_db = tmp_db.replace('gobe','')
    print >>sys.stderr, "TEMP_DB:", tmp_db

    tmp_db = get_temp_db(data['base_name'])
    tmp_db.row_factory = sqlite3.Row
    tmp_cur = tmp_db.cursor()

    data = dict(list(data.iteritems()))
    if not 'user' in data: data['user'] = 'unknown'

    datasets = [x[0] for x in tmp_db.execute('SELECT dsid FROM image_info ORDER BY id').fetchall()]
    chroms   = [x[0] for x in tmp_db.execute('SELECT chromosome FROM image_info ORDER BY id').fetchall()]
    print >>sys.stderr, datasets

    qsql = 'SELECT * FROM image_data WHERE id IN (' + ",".join([str(p[0]) for p in data['hsp_ids']]) + ');';
    ssql = 'SELECT * FROM image_data WHERE id IN (' + ",".join([str(p[1]) for p in data['hsp_ids']]) + ');';

    # get rid of the old ones. but save history in genespace
    tcur.execute("UPDATE pair SET previous_genespace_id = genespace_id WHERE genespace_id = ?", (data['genespace_id'],))
    tcur.execute("UPDATE pair SET genespace_id = -1 WHERE genespace_id = ?", (data['genespace_id'],))

    ##################################################################################################
    # save the extents of the genespace.                                                             #
    ##################################################################################################
    tcur.execute("INSERT INTO PAIR VALUES(NULL, 'genespace', -5, ?, NULL, NULL, ?)"
            , (data['genespace_id'], data['user']))
    gspair_id = tcur.lastrowid
    tcur.execute("INSERT INTO location VALUES (NULL, ?, ?, 'NAME', NULL, 'q', ? , ?, ?, NULL)"
                , (datasets[0], gspair_id, chroms[0], data['qextents'][0], data['qextents'][1]))

    tcur.execute("INSERT INTO location VALUES (NULL, ?, ?, 'NAME', NULL, 's', ?, ?, ?, NULL)"
                , (datasets[1], gspair_id, chroms[1], data['sextents'][0], data['sextents'][1]))


    for qloc, sloc in zip(tmp_db.execute(qsql).fetchall(), tmp_db.execute(ssql).fetchall()):


        #assert sloc['bpmin'] < 454499, ('got query, need subject', sloc['bpmin'], ssql)
        # TODO: orientation (-5 ???)
        tcur.execute("INSERT INTO pair VALUES(NULL, ?, -5, ?, NULL, NULL, ?)"
                            ,('CNS', data['genespace_id'], data['user'] ))

        pair_id = tcur.lastrowid
        tcur.execute("INSERT INTO location VALUES (NULL, ?, ?, 'NAME', NULL, 'q', ?, ?, ?, NULL)"
                , (datasets[0], pair_id, chroms[0], qloc['bpmin'], qloc['bpmax']))

        tcur.execute("INSERT INTO location VALUES (NULL, ?, ?, 'NAME', NULL, 's', ?, ?, ?, NULL)"
                , (datasets[1], pair_id, chroms[1], sloc['bpmin'], sloc['bpmax']))

    tracking_db.commit()
    tmp_db.close()
    print >>sys.stderr, 'success'

    return True


def remove(genespace_id):
    #tcur.executescript("UPDATE genespace SET genespace_type = 'INVALID' where genespace_id = ?",(genespace_id,)); 
    tracking_db.commit()
    return "UPDATE genespace SET genespace_type = 'INVALID' where genespace_id = %i" % genespace_id



@bagwrap
def load(genespace_id, base_name):
    print >>sys.stderr, "IN LOAD" 
    genespace_id = int(genespace_id)
    info = tcur.execute('SELECT * FROM genespace WHERE genespace_id = ?', (genespace_id,)).fetchone()

    if info is None: 
        tracking_db.close()
        return False

    locs = tcur.execute("SELECT l.start, l.stop, l.q_or_s FROM location l, pair p WHERE p.pair_id = l.pair_id \
                                           AND p.genespace_id = ? AND p.pair_type = 'CNS' and l.q_or_s = 'q'"
                                        , (genespace_id,)).fetchall()
    bars = tcur.execute("SELECT l.start, l.stop, l.q_or_s FROM location l, pair p WHERE p.pair_id = l.pair_id \
                                           AND p.genespace_id = ? AND p.pair_type = 'genespace'"
                                        , (genespace_id,)).fetchall()

    qextents = [[b['start'], b['stop']] for b in bars if b['q_or_s'] == 'q']
    sextents = [[b['start'], b['stop']] for b in bars if b['q_or_s'] == 's']
    print >>sys.stderr, qextents

    tmp_db =  get_temp_db( base_name )
    tmp_db.row_factory = sqlite3.Row
    tmp_cur = tmp_db.cursor()
    coordslist = []
    for l in locs:
        # TODO: is order by image_id sufficient to assure we get the
        # query?
        qbps = tmp_cur.execute('SELECT xmin, ymin, xmax, ymax, id, pair_id  FROM image_data WHERE bpmin = ? AND bpmax = ? AND ABS(image_track) > 1 ORDER BY image_id', (l['start'], l['stop'])).fetchone()
        sbps = tmp_cur.execute('SELECT xmin, ymin, xmax, ymax, id  FROM image_data WHERE id = ?', (qbps['pair_id'],)).fetchone()
        print >>sys.stderr,  qbps, sbps
        coordslist.append({
                    'img1': [ qbps['xmin'], qbps['ymin'], qbps['xmax'], qbps['ymax'], qbps['id']]
                  , 'img2': [ sbps['xmin'], sbps['ymin'], sbps['xmax'], sbps['ymax'], sbps['id']]
                  });
    print >>sys.stderr, "COORDS\n"
    print >>sys.stderr, coordslist
    kwds = [int(k) for k in info['keywords'].split("|") if k ]
    anns = [int(k) for k in  info['annotation'].split("|") if k]
    print >>sys.stderr, {'notes': info['notes'], 'qdups': info['qdups'], 'sdups':info['sdups'] ,'annos':anns
            ,'keywords':kwds, 'revisit':bool(info['revisit']), 'features': coordslist
            ,'sextents': sextents, 'qextents':qextents }

    return {'notes': info['notes'], 'qdups': info['qdups'], 'sdups':info['sdups'] ,'annos':anns
            ,'keywords':kwds, 'revisit':bool(info['revisit']), 'features': coordslist
            ,'sextents': sextents, 'qextents':qextents }



def new_genespace(qanchor, sanchor, base_name):
    print >>sys.stderr, qanchor, sanchor, tmp_db
    gsid = tcur.execute('SELECT MAX(genespace_id) FROM genespace').fetchone()[0] + 1

    print >>sys.stderr, qanchor, sanchor, tmp_db
    tmp_db = get_temp_db(base_name)
    tmp_db.row_factory = sqlite3.Row
    tmp_cur = tmp_db.cursor()
    datasets = [x[0] for x in tmp_db.execute('SELECT dsid FROM image_info ORDER BY id').fetchall()]


    tcur.execute("INSERT INTO genespace VALUES (?, 'our_additional', NULL, 'user', 'f', '', '', NULL, NULL, 'Our Additional')"
            , (gsid,))
    tcur.execute("INSERT INTO pair VALUES (NULL, 'anchor', 1, ?, NULL, NULL, 'user')"
            , ( gsid,))
    pair_id = tcur.lastrowid
    tcur.execute("INSERT INTO location VALUES  (NULL, ?, ?, 'name', 0, 'q', -5, ?, ?, ?)"
            ,(datasets[0],pair_id, qanchor, sanchor, 1))

    return [qanchor, sanchor]

application = WSGIGateway({
    'save': save
    ,'remove': remove
    ,'new_genespace': new_genespace
    ,'load': load
    })

if __name__ == "__main__":
    

    #print save({'tmp_db': u'tmpdir//GEvo_Fkdb8kIf.sqlite', 'notes': u'sfasd', 'annos': [], 'sextents': [8458456, 8489597], 'revisit': False, 'qextents': [1078486, 1098951], 'sdups': 1, 'genespace_id': 5, 'qdups': 0, 'keywords': [2], 'hsp_ids': [[55, 159], [62, 161]]},)
    print load(2, 'GEvo_Fkdb8kIf')
    #print new_genespace(123, 234, 'tmpdir//GEvo_11ltNR3a.sqlite')
