#!/usr/bin/python

import web
import os
try:
    import sqlite3
except ImportError:
    from pysqlite2 import dbapi2 as sqlite3
import simplejson
import urllib

TMPDIR = "../tmp/GEvo/"
if not os.path.exists(TMPDIR):
    TMPDIR = os.path.join(os.path.dirname(__file__), TMPDIR)
DBTMPL = os.path.join(TMPDIR, "%s.sqlite")

def getdb(dbname):
    db = sqlite3.connect(DBTMPL % dbname)
    db.row_factory = sqlite3.Row
    return db


class info(object):
    def GET(self, dbname):
        web.header('Content-type', 'text/javascript')
        db = getdb(dbname)
        c = db.cursor()
        c2 = db.cursor()
        c.execute("SELECT * FROM image_info order by display_id")
        c2.execute("SELECT min(xmin) as min, max(xmax) as max, image_id FROM image_data WHERE type='anchor' GROUP BY image_id ORDER BY image_id")

        result = {}
        for i, (row, anchor) in enumerate(zip(c, c2)):
            result[row['iname']] = dict(
                title=row['title'],
                i=i,
                img_width=row['px_width'],
                bpmin=row['bpmin'],
                bpmax=row['bpmax'],
                idx=row['id'],
                xmin=anchor['min'],
                xmax=anchor['max']
            )
        return simplejson.dumps(result)


class follow(object):
    def GET(self, dbname):
        web.header('content-type', 'text/javascript')
        db = getdb(dbname)
        c = db.cursor()
        c2 = db.cursor()
        img = web.input(img=None).img

        bbox = map(float, web.input().bbox.split(","))
        ids = []
        pair_ids = []

        used_pairs = []

        def get_pair_data(pair_id):
            c.execute("""SELECT xmin, xmax, ymin, ymax, image_id, image_track FROM image_data WHERE id = ?""", (pair_id,))
            p = c.fetchone()
            return dict(
                         pair_id=pair_id,
                         pair_image_id=p['image_id'],
                          pair_track=p['image_track'],
                         pair_bbox=(p['xmin'], p['ymin'], p['xmax'], p['ymax']))

        def get_pairs(img_id, bbox):
            c.execute("""SELECT id, pair_id, image_id, xmin, xmax, ymin, ymax FROM image_data WHERE ? + 1 > xmin AND ? - 1 < xmax AND 
              ? - 1 > ymin AND ? + 1 < ymax AND image_id = ? AND pair_id != -99 AND type = 'HSP'""", \
              (bbox[2], bbox[0], bbox[3], bbox[1], img_id))
            results = c.fetchall()
            if not results: return None
            pairs = []
            for r in results:
                d = dict(id=r['id'], bbox=(r['xmin'], r['ymin'], r['xmax'], r['ymax']), image_id=r['image_id'])
                d.update(get_pair_data(r['pair_id']))
                pairs.append(d)
            
            return pairs 

        def get_pairs_for_bbox_image(xmin, xmax, img_id, exclude_track):
            c.execute("""SELECT id, pair_id, image_id, xmin, xmax, ymin, ymax 
                      FROM image_data WHERE ? + 1 > xmin AND ? - 1 < xmax AND 
                      image_id = ? AND pair_id != -99 AND image_track != ? AND type = 'HSP'""", \
              (xmax, xmin, img_id, exclude_track))
            web.debug("""SELECT id, pair_id, image_id, xmin, xmax, ymin, ymax 
                      FROM image_data WHERE ? + 1 > xmin AND ? - 1 < xmax AND 
                      image_id = ? AND pair_id != -99 AND image_track != ? AND type = 'HSP'""")
            web.debug((xmax, xmin, img_id, exclude_track))

            results = c.fetchall()
            pairs = []
            for r in results:
                d = dict(id=r['id'], bbox=(r['xmin'], r['ymin'], r['xmax'], r['ymax']), image_id=r['image_id'])
                d.update(get_pair_data(r['pair_id']))
                pairs.append(d)
            return pairs 
        
        pairs = get_pairs(img, bbox)
        i = 0
        while True:
            L = len(pairs)
            if i == L: break
            pair = pairs[i]
            new_pairs = get_pairs(pair['pair_image_id'], pair['pair_bbox'])
            for np in (_np for _np in new_pairs if not _np in pairs):
                new2 = get_pairs_for_bbox_image(np['bbox'][0], np['bbox'][2], np['image_id'], np['pair_track'])
                if new2 == []: continue
                if not new2 in pairs:
                    pairs.append(new2)
            pairs.extend([np for np in new_pairs if not np in pairs])
            i += 1
            if L == len(pairs): break
            web.debug(pairs)


class query(object):
    def GET(self, dbname):
        db = getdb(dbname)
        c = db.cursor()
        img = web.input(img=None).img

        if web.input(bbox=None).bbox:
            bbox = map(float, web.input().bbox.split(","))
            c.execute("""SELECT * FROM image_data WHERE ? + 1 > xmin AND ? - 1 < xmax AND 
                      ? - 1 > ymin AND ? + 1 < ymax AND image_id = ? AND pair_id != -99 AND type = 'HSP'""", \
                      (bbox[2], bbox[0], bbox[3], bbox[1], img))
        elif web.input(all=None).all:
            c.execute("""SELECT distinct(image_track) as image_track FROM image_data WHERE ? 
                      BETWEEN ymin AND ymax AND image_id = ? ORDER BY 
                      ABS(image_track) DESC""", (float(web.input().y), img))
            track = c.fetchone()['image_track']
            web.debug(track)
            c.execute("""SELECT id, xmin, xmax, ymin, ymax, image_id, image_track, pair_id, color, link FROM image_data 
                    WHERE ( (image_track = ?) or (image_track = (? * -1) ) ) 
                    and image_id = ? and pair_id != -99 and type = 'HSP'""", (track, track, img))

        else: # point query.
            x = float(web.input().x)
            y = float(web.input().y)
            c.execute("""SELECT * FROM image_data WHERE ? + 3 > xmin AND ? - 3
                      < xmax AND ? BETWEEN ymin and ymax and image_id = ?""", 
                      (x, x, y, img))

        c2 = db.cursor()
        # now iterate over the cursor
        results = []
        for result in c:
            c2.execute("""SELECT id, xmin, xmax, ymin, ymax, image_id, 
                       image_track, pair_id, color FROM image_data where 
                       id = ?""", (result['pair_id'], ));
            pair = c2.fetchone()
            try:
                anno = result['annotation']
                if anno.startswith('http'):
                    anno = urllib.urlopen(anno).read()
            except:
                anno = ""

            f1pts = []
            f2pts = []
            for k in ('xmin', 'ymin', 'xmax', 'ymax'):
                f1pts.append(int(round(result[k])))
                if pair:
                    f2pts.append(int(round(pair[k])))
                
            f1pts.extend([result['id'], result['image_track']])
            if pair:
                f2pts.extend([pair['id'], pair['image_track']])
            results.append(dict(
                # TODO: tell eric to add 'CoGe' to the start of his links.
                link=result['link'],
                annotation = anno,
                # TODO has_pair
                has_pair= bool(pair),
                color=(result['color'] or (pair and pair['color'])).replace('#', '0x'),
                features={
                    'key%i' % result['image_id']: f1pts,
                    'key%i' % (pair and pair['image_id'] or 999): f2pts}
            ))
	if not results:
	    results.append(dict(
		annotation = "No results found"
	    ))
        web.header('Content-type', 'text/javascript')
        return simplejson.dumps({'resultset':results})


urls = (
    # the first pattern is always the sqlite db name. e.g.: /GEVo_WxUonWBr/info
    '/([^\/]+)/info/', 'info',
    '/([^\/]+)/follow/', 'follow',
    '/([^\/]+)/query/', 'query',
)


app = web.application(urls, locals())
application = app.wsgifunc()

if __name__ == "__main__":
    app.run()
