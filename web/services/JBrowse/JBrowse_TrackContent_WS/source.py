import MySQLdb as mdb
import json
import math
import re
import random
import os
from collections import defaultdict #Counter
from cgi import parse_qs, escape

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

def gc_features(environ, start_response):
    """Main feature endpoint for GC content"""
    status = '200 OK'
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

    con = db_connect()
    cur = con.cursor()

    try:
        cur.execute("SELECT file_path FROM genome where genome_id = %d;"
                % genome_id)

        # Open the right chromosome file derived from the pathname
        file = cur.fetchone()[0]
        file = '/'.join(file.split('/')[:-1]) + '/chr/' + chr_id
        f = open(file)
        if (start < 0):
            string = f.read(end).lower()
        else:
            f.seek(start)
            string = f.read(end - start).lower()
        f.close()

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
            if (start + i + bucketSize < start):
                k = start
            else:
                k = start + i
            if (start + i + bucketSize >= end):
                j = end
            else:
                j = start + i + bucketSize
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

    finally:
        if con:
            con.close()

    response_headers = [('Content-Type', 'application/json')]
    start_response(status, response_headers)

    return response_body

def an_features(environ, start_response):
    """Main feature endpoint for Annotation Feature content"""
    status = '200 OK'
    response_body = { "features" : [] }
    bucketSize = 100

    # Get the passed params of the AJAX request
    d = parse_qs(environ['QUERY_STRING'])
    start = d.get('start', [''])[0]
    end = d.get('end', [''])[0]
    args = environ['url_args']

    # set parsed argument variables
    genome_id = args['genome_id']
    chr_id = args['chr_id']
    try:
        feat_type = args['feat_type']
    except KeyError:
        feat_type = ""


    con = db_connect()
    cur = con.cursor()

    query = "SELECT l.start, l.stop, l.strand, ft.name, fn.name, \
            l.location_id, f.start, f.stop, f.feature_id \
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
            AND ft.feature_type_id != 4" \
            .format(genome_id, chr_id, start, end)
    if feat_type:
        query +=  " AND ft.name = '{0}'".format(feat_type)

    try:
        cur.execute(query + ";")
        results = cur.fetchall()

        if feat_type:
            if results:
                response_body["features"].append({"subfeatures" : []})
            i = 0
            lastID = 0
            lastEnd = 0
            lastStart = 0
            for row in results:
                if row[8] != lastID and lastID != 0:
                    response_body["features"][i]["start"] = lastStart
                    response_body["features"][i]["end"] = lastEnd
                    response_body["features"][i]["uniqueID"] = lastID
                    response_body["features"][i]["name"] = lastName
                    response_body["features"][i]["strand"] = lastStrand
                    response_body["features"][i]["type"] = lastType
                    response_body["features"].append({"subfeatures" : []})
                    i += 1

                elif lastID == 0:
                    response_body["features"][i]["start"] = row[6]
                    response_body["features"][i]["end"] = row[7]
                    response_body["features"][i]["uniqueID"] = row[8]
                    response_body["features"][i]["name"] = row[4]
                    response_body["features"][i]["strand"] = row[2]
                    response_body["features"][i]["type"] = row[3]

                response_body["features"][i]["subfeatures"].append({
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
        else:
            for row in results:
                response_body["features"].append({
                    "start" : row[0],
                    "end" : row[1],
                    "strand" : row[2],
                    "type" : row[3],
                    "name" : row[4],
                })

    except mdb.Error, e:
        response_body = "Error %d: %s" % (e.args[0], e.args[1])
        status = '500 Internal Server Error'

    finally:
        if con:
            con.close()

    response_headers = [('Content-Type', 'application/json')]
    start_response(status, response_headers)

    #response_body = feat_type
    response_body = json.dumps(response_body)
    return response_body

def stats(environ, start_response):
    start_response('200 OK', [('Content-Type', 'text/plain')])
    response_body = { "featureDensity": 0.02,
                      "scoreMin": 0,
                      "scoreMax": 1,
                      "featureDensityByRegion" : 50000,
                    }
    return json.dumps(response_body)

def region(environ, start_response):
    start_response('200 OK', [('Content-Type', 'text/plain')])
    return ['{}']

urls = [
    (r'stats/global$',
        stats),
    (r'stats/region/?$',
        region),
    (r'annotation/(?P<genome_id>\d+)/types/(?P<feat_type>\w+)/features/(?P<chr_id>\w+)?(.+)?$',
        an_features),
    (r'annotation/(?P<genome_id>\d+)/features/(?P<chr_id>\w+)?(.+)?$',
        an_features),
    (r'gc/(?P<genome_id>\d+)/features/(?P<chr_id>\w+)?(.+)?$',
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
