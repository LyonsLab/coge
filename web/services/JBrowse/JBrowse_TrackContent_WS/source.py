import Cookie
import json
import math
import re
import random
import os
import sys
import string
import urllib2
from collections import defaultdict #Counter
from cgi import parse_qs, escape
from time import time

import MySQLdb as mdb

config_path = os.path.join(os.path.dirname(__file__), '../../../coge.conf')

def not_found(environ, start_response):
    """Called if no URL matches."""
    start_response('404 NOT FOUND', [('Content-Type', 'text/plain')])
    return [environ.get('PATH_INFO', '').lstrip('/')]

def db_connect():
    dir = os.path.dirname(__file__)
    path = os.path.join(dir, '../../../coge.conf')
    config = {}
    # Open config settings from config file
    f = open(path)

    # Parse config settings from opened file
    for line in f:
        if line.split('\t')[0].startswith('DB'):
            line = line.split('\t')
            config[line[0]] = line[-1].strip()
    f.close()

    # Connect to database passing config settings
    con = mdb.connect(
            host=config['DBHOST'],
            user=config['DBUSER'],
            passwd=config['DBPASS'],
            db=config['DBNAME'],
            port=int(config['DBPORT']))

    return con

def get_config():
    config = {}

    with open(config_path) as fp:
        for line in fp:
            line = line.strip()

            if line.startswith('\n') or line.startswith('#'):
                continue
            else:
                try:
                    (key, value) = line.replace('\t', ' ').split()[0:2]
                    config[key] = value
                except ValueError:
                    pass

    return config

def fetch_sequence(genome_id, chr_id, start, stop, cookie_string):
    service = '{base}{service}/sequence/{id}/{chr}?start={start};stop={stop};'

    config = get_config()
    if not config:
        return

    url = service.format(base=config['SERVER'],
        #service='services/JBrowse/service.pl',
        service='api/v1/legacy', # mdb added 2/5/15, COGE-289
        id=genome_id, chr=chr_id, start=start, stop=stop)

    sys.stderr.write(url)
    try:
        name = config['COOKIE_NAME']

        if cookie_string:
            opener = urllib2.build_opener()
            opener.addheaders.append(('Cookie', cookie_string))
            response = opener.open(url)
        else:
            response = urllib2.urlopen(url)

        sequence = response.read()
    except urllib2.URLError as e:
        return

    return sequence.lower()

# mdb added 11/4/13 issue 246 - add wobble shading
def calc_wobble(seq, start_offset, stop_offset):
    if (start_offset < 0):
        start_offset = 0
    if (stop_offset < 0 or start_offset == stop_offset):
        return 0

    # Calculate GC percent
    count = 0
    total = 0;
    for i in xrange(start_offset+2, stop_offset, 3):
        total += 1
        c = seq[i].upper()
        if (c == 'G' or c == 'C'):
            count += 1

    return (count, total)

def gc_features(environ, start_response):
    """Main feature endpoint for GC content"""
    #sys.stderr.write('gc_features\n')
    status = '200 OK'
    response_headers = [('Content-Type', 'application/json')]
    response_body = { "features" : [] }
    bucketSize = 100

    # Get the passed params of the AJAX request
    d = parse_qs(environ['QUERY_STRING'])
    start = int(d.get('start', [''])[0])
    end = int(d.get('end', [''])[0])
    scale = d.get('scale', [''])[0]
    basesPerSpan = d.get('basesPerSpan', [''])[0]
    args = environ['url_args']

    # set parsed argument variables
    genome_id = int(args['genome_id'])
    chr_id = args['chr_id']
    
    if (start < 0):
        start = 0;
    if (end <= 0):
        end = 0
        
    if (start == end):
        start_response(status, response_headers)
        return ''
        
    try:
        # Get chromosome subsequence using interbase coordinates
        cookie = get_cookie(environ)
        string = fetch_sequence(genome_id, chr_id, start+1, end+1, cookie)

        # Set bucketSize
        sizes = {'20': 1, '10': 1, '5': 2, '2': 5, '1': 25, '0.5': 75}
        try:
            bucketSize = sizes[scale]
        except KeyError:
            bucketSize = int(1 / math.pow(2, math.log10(float(scale))) * 50)

        for i in xrange(0, len(string), bucketSize):
            # Score becomes the length of the string subtracting all 'atnx'
            chunk = string[i:i+bucketSize]
            # FIX FOR NON PYTHON 2.7
            matches = defaultdict(int)
            for char in chunk:
                matches[char] += 1
            nucleotide = max(matches.iteritems(), key=lambda x: x[1])[0].lower()
            #Counter(chunk).most_common(1)[0][0].lower()
            score = len(re.sub('[atnx]', '', chunk))
            score = str(round(score / float(len(chunk)), 3))
            k = start + i
            j = start + i + bucketSize
            if (j > end):
                j = end
            response_body['features'].append({
                "start": k,
                "score": score,
                "end": j,
                "nucleotide": nucleotide,
            })

        response_body = json.dumps(response_body)

    except mdb.Error, e:
        response_body = "Error %d: %s" % (e.args[0], e.args[1])
        status = '500 Internal Server Error'

    start_response(status, response_headers)
    return response_body

def an_features(environ, start_response): # mdb rewritten 11/8/13 issue 246 - add wobble calc
    """Main feature endpoint for Annotation Feature content"""
    status = '200 OK'
    response_body = { "features" : [] }
    bucketSize = 100

    # Get the passed params of the AJAX request
    d = parse_qs(environ['QUERY_STRING'])
    start = d.get('start', [''])[0]
    end = d.get('end', [''])[0]
    try:
        show_wobble = int(d.get('showWobble', [''])[0])
    except:
        show_wobble = 0

    # set parsed argument variables
    args = environ['url_args']
    genome_id = args['genome_id']
    chr_id = args['chr_id']
    try:
        feat_type = args['feat_type']
    except KeyError:
        feat_type = ""
    # mdb added 4/21/14 issue 363
    try:
        dataset_id = args['dataset_id']
    except KeyError:
        dataset_id = 0;

    sys.stderr.write('an_features: '+str(genome_id)+' '+str(chr_id)+' '+str(start)+' '+str(end)+'\n')

    con = db_connect()
    cur = con.cursor()
    # mdb 4/24/14 - added fn.primary_name=1 constraint to keep from returning arbitrary name
    query = "SELECT l.start, l.stop, l.strand, ft.name, fn.name, \
            l.location_id, f.start, f.stop, f.feature_id, fn.primary_name, \
	    l.chromosome, f.chromosome, f.strand \
            FROM genome g \
            JOIN dataset_connector dc ON dc.genome_id = g.genome_id \
            JOIN dataset d on dc.dataset_id = d.dataset_id \
            JOIN feature f ON d.dataset_id = f.dataset_id \
            JOIN location l ON f.feature_id = l.feature_id \
            JOIN feature_name fn ON f.feature_id = fn.feature_id \
            JOIN feature_type ft ON f.feature_type_id = ft.feature_type_id \
            WHERE g.genome_id = {0} \
                AND f.chromosome = '{1}' \
                AND f.stop > {2} AND f.start <= {3} \
                AND ft.feature_type_id != 4".format(genome_id, chr_id, start, end)
    # erb 5/16/14 issue 381 - all genbank genomes have fn.primary_name=0
    #            AND fn.primary_name = 1"
    #       .format(genome_id, chr_id, start, end)

    if feat_type:
        query +=  " AND ft.name = '{0}'".format(feat_type)
    # mdb added 4/21/14 issue 363
    if dataset_id:
        query += " AND dc.dataset_id = {0}".format(dataset_id)

    try:
        # Query for features
        #sys.stderr.write(str(query))
        cur.execute(query + ";")
        results = cur.fetchall()

        if feat_type: # feature type was specified
            if results:
                response_body["features"].append({"subfeatures" : []})
            i = 0
            lastID = 0
            lastEnd = 0
            lastStart = 0
            for row in results:
                if row[8] != lastID and lastID != 0:
                    response_body["features"][i]["chr"] = chr_id
                    response_body["features"][i]["start"] = lastStart
                    response_body["features"][i]["end"] = lastEnd
                    response_body["features"][i]["uniqueID"] = lastID
                    response_body["features"][i]["name"] = lastName
                    response_body["features"][i]["strand"] = lastStrand
                    response_body["features"][i]["type"] = lastType
                    response_body["features"].append({"subfeatures" : []})
                    i += 1

                elif lastID == 0:
                    response_body["features"][i]["chr"] = chr_id
                    response_body["features"][i]["start"] = row[6]
                    response_body["features"][i]["end"] = row[7]
                    response_body["features"][i]["uniqueID"] = row[8]
                    response_body["features"][i]["name"] = row[4]
                    response_body["features"][i]["strand"] = row[2]
                    response_body["features"][i]["type"] = row[3]

		if row[2] == row[12] and row[10] == row[11] and row[0] >= row[6] and row[0] <= row[7] and row[1] >= row[6] and row[1] <= row[7]:
                    response_body["features"][i]["subfeatures"].append({
                        "chr" : chr_id, # mdb added 11/18/13 issue 240 - add chr to FeatAnno.pl onclick url
                        "start" : row[0],
                        "end" : row[1],
                        "strand" : row[2],
                        "type" : row[3],
                        "name" : row[4],
                    })

                lastStrand = row[2]
                lastType = row[3]
                lastName = row[4]
                lastStart = row[6]
                lastEnd = row[7]
                lastID = row[8]
        else: # process all feature types
            # Filter out redundant results due to multiple feature names
            feats = {}
            min_start = None
            max_stop = None
            for row in results:
		if row[2] == row[12] and row[10] == row[11] and row[0] >= row[6] and row[0] <= row[7] and row[1] >= row[6] and row[1] <= row[7]:
                    # Find bounds for sequence retrieval later (show_wobble == 1)
                    if min_start is None or min_start > row[0]:
                        min_start = row[0]
                    if max_stop is None or max_stop < row[1]:
                        max_stop = row[1]
                
                    # Hash unique features
                    location_id = row[5]
                    try:
                        feat = feats[location_id]
                        if row[9] == 1: # is feature name the primary name?
                            feats[location_id] = row
                    except KeyError:
                        feats[location_id] = row

            # Calculate wobble GC for CDS features
            wcount = {}
            wtotal = {}
            if show_wobble:
                seq = '';
                rcseq = '';
                for row in feats.values():
                    if row[3] == 'CDS' and is_overlapping(int(row[0]), int(row[1]), int(start), int(end)):
                        if not seq:
                            # Get chromosome subsequence using interbase coordinates
                            cookie = get_cookie(environ)
                            seq = fetch_sequence(genome_id, chr_id, int(min_start), int(max_stop), cookie)
                            
                        if (int(row[2]) == 1): # plus strand
                            count, total = calc_wobble(seq, int(row[0])-int(min_start), int(row[1])-int(min_start))
                        else: # minus strand
                            # Reverse complement the seq and coords
                            if not rcseq:
                                rcseq = reverse_complement(seq)
                            rcstart = int(max_stop)-int(row[1])
                            rcstop  = int(max_stop)-int(row[0])
                            count, total = calc_wobble(rcseq, rcstart, rcstop)
                            
                        id = row[8]
                        if id in wcount:
                            wcount[id] += count
                            wtotal[id] += total
                        else:
                            wcount[id] = count
                            wtotal[id] = total
            
            # Build response
            for row in feats.values():
                response_body["features"].append({
                    "chr"    : chr_id, # mdb added 11/18/13 issue 240 - add chr to FeatAnno.pl onclick url
                    "start"  : row[0],
                    "end"    : row[1],
                    "strand" : row[2],
                    "type"   : row[3],
                    "name"   : row[4] if row[3] == 'gene' else None,
                    "wobble" : round(wcount[row[8]] / float(wtotal[row[8]]), 3) if (show_wobble and row[8] in wcount) else 0
                })

    except mdb.Error, e:
        response_body = "Error %d: %s" % (e.args[0], e.args[1])
        status = '500 Internal Server Error'

    finally:
        if con:
            con.close()

    response_headers = [('Content-Type', 'application/json')]
    start_response(status, response_headers)
    response_body = json.dumps(response_body)
    return response_body

def stats(environ, start_response):
    #sys.stderr.write('global stats\n')
    start_response('200 OK', [('Content-Type', 'text/plain')])
    response_body = { "featureDensity": 0.02,
                      "scoreMin": 0,
                      "scoreMax": 1,
                      "featureDensityByRegion" : 50000,
                    }
    return json.dumps(response_body)

def region(environ, start_response):
    #sys.stderr.write('region stats\n')
    
    # Get the passed params of the AJAX request
    d = parse_qs(environ['QUERY_STRING'])
    start = int( d.get('start', [''])[0] )
    end = int( d.get('end', [''])[0] )
    bpPerBin = int( d.get('bpPerBin', [''])[0] )
    args = environ['url_args']

    response_body = { "bins" : [], "stats" : [ { "basesPerBin": bpPerBin, "max": 0, "mean": 0 } ] }
    status = '200 OK'

    # set parsed argument variables
    genome_id = args['genome_id']
    chr_id = args['chr_id']
    feat_type = args['feat_type']

    con = db_connect()
    cur = con.cursor()

    query = "SELECT l.start, l.stop \
            FROM genome g \
            JOIN dataset_connector dc ON dc.genome_id = g.genome_id \
            JOIN dataset d on dc.dataset_id = d.dataset_id \
            JOIN feature f ON d.dataset_id = f.dataset_id \
            JOIN location l ON f.feature_id = l.feature_id \
            JOIN feature_name fn ON f.feature_id = fn.feature_id \
            JOIN feature_type ft ON f.feature_type_id = ft.feature_type_id \
            WHERE g.genome_id = {0} \
            AND f.chromosome = '{1}' \
            AND f.stop > {2} AND f.start <= {3} \
            AND ft.feature_type_id != 4 \
            AND ft.name = '{4}'" \
            .format(genome_id, chr_id, start, end, feat_type)

    try:
        cur.execute(query + ";")
        results = cur.fetchall()

        if results:
            maxStop = max([int(row[1]) for row in results]) - start
            numBins = bin(maxStop, maxStop, bpPerBin) + 1
            bins = [0] * numBins
    
            for row in results:
                s = int(row[0]) - start
                e = int(row[1]) - start
                b = bin(s, e, bpPerBin)
                bins[b] += 1
                
            maxVal = max(bins)
            meanVal = sum(bins) / len(bins)
            response_body = { "bins" : bins, "stats" : [ { "basesPerBin": bpPerBin, "max": maxVal, "mean": 0 } ] }
            status = '200 OK'

    except mdb.Error, e:
        response_body = "Error %d: %s" % (e.args[0], e.args[1])
        status = '500 Internal Server Error'

    finally:
        if con:
            con.close()

    response_headers = [('Content-Type', 'application/json')]
    start_response(status, response_headers)
    response_body = json.dumps(response_body)
    return response_body

def is_overlapping(s1, e1, s2, e2):
    return (s1 <= e2 and s2 <= e1)

def reverse_complement(dna):
    complements = string.maketrans('acgtACGT', 'tgcaTGCA')
    rcseq = dna.translate(complements)[::-1]
    return rcseq

def bin(start, end, bpPerBin):
    return max(0, int((start + end) / 2 / bpPerBin))

def get_cookie(environ): # mdb added 2/5/15 for GEISHA-JBrowse, COGE-585
    try:
        return environ['HTTP_COOKIE']
    except KeyError:
        return None

urls = [
    (r'stats/global$',
        stats),
    (r'annotation/(?P<genome_id>\d+)/types/(?P<feat_type>[\w\:\-]+)/stats/region/(?P<chr_id>\S+)?$',
        region),
    (r'annotation/(?P<genome_id>\d+)/types/(?P<feat_type>[\w\:\-]+)/features/(?P<chr_id>\S+)?$',
        an_features),
    (r'annotation/(?P<genome_id>\d+)/features/(?P<chr_id>\S+)?$',
        an_features),
    (r'annotation/(?P<genome_id>\d+)/datasets/(?P<dataset_id>\S+)/features/(?P<chr_id>\S+)?$',
        an_features), # mdb added 4/21/14 issue 363
    (r'gc/(?P<genome_id>\d+)/features/(?P<chr_id>\S+)?$',
        gc_features),
]

def application(environ, start_response):
    """
    The main WSGI application. Dispatch the current request to
    the functions from above and store the regular expression
    captures in the WSGI environment as  `myapp.url_args` so that
    the functions from above can access the url placeholders.

    If nothing matches call the `not_found` function.
    """
    path = environ.get('PATH_INFO', '').lstrip('/')
    for regex, callback in urls:
        match = re.search(regex, path)
        if match is not None:
            environ['url_args'] = match.groupdict()
            return callback(environ, start_response)
    return not_found(environ, start_response)
