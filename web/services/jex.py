import json
import math
import re
import random
import os
import sys
from collections import defaultdict #Counter
from cgi import parse_qs, escape

import zmq

_defaults = {
    'connection' : 'tcp://localhost:5151'
}

def not_found(environ, start_response):
    """Called if no URL matches."""
    start_response('404 NOT FOUND', [('Content-Type', 'text/plain')])
    return [environ.get('PATH_INFO', '').lstrip('/')]

def stats(environ, start_response):
    sys.stderr.write('global stats\n')
    start_response('200 OK', [('Content-Type', 'text/plain')])
    response_body = { "featureDensity": 0.02,
                      "scoreMin": 0,
                      "scoreMax": 1,
                      "featureDensityByRegion" : 50000,
                    }
    return json.dumps(response_body)

def status(environ, start_response):
    start_response('200 OK', [('Content-Type', 'application/json')])

    result = {}

    try:
        args = environ['url_args']
        identity = args['id']
    except KeyError:
        return json.dumps(result)

    if not identity:
        return json.dumps(result)

    try:
        context = zmq.Context()
        socket = context.socket(zmq.REQ)
        socket.connect(_defaults['connection'])
        request = {
            'request' : 'get_status',
            'data' : {
                'id' : args['id']
            },
        }

        socket.send_json(request)
        result = socket.recv_json()
    except:
        result = {}

    return json.dumps(result)

urls = [
    (r"synmap/status/(?P<id>[a-z0-9\-]+)?", status),
    (r"synfind/status/(?P<id>[a-z0-9\-]+)?", status),
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
