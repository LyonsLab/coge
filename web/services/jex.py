import json
import itertools
import math
import re
import random
import os
import sys
import time
from collections import defaultdict #Counter
from cgi import parse_qs, escape

import zmq

_defaults = {
    'connection' : 'tcp://{}:{}',
    'max_attempts' : 30,
    'timeout' : 1000
}


def not_found(environ, start_response):
    """Called if no URL matches."""
    start_response('404 NOT FOUND', [('Content-Type', 'text/plain')])
    return [environ.get('PATH_INFO', '').lstrip('/')]

def load_config(filename):
    config = {}
    with open(filename, 'r') as fp:
        for line in fp:
            if line.startswith('#') or len(line) < 1:
                continue
            else:
                cleaned = line.strip().replace('\t', '')
                try:
                    option, value = re.split(' ', cleaned)
                    config[option] = value
                except ValueError:
                    pass

    return config

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

    result = None

    try:
        args = environ['url_args']
        identity = args['id']
    except KeyError:
        return json.dumps({})

    if not identity:
        return json.dumps({})

    context = zmq.Context()
    socket = context.socket(zmq.REQ)
    socket.setsockopt(zmq.LINGER, 0)


    root_dir = os.path.dirname(__file__)
    filepath = os.path.abspath(os.path.join(root_dir, environ['COGE_HOME'], 'coge.conf'))

    config = load_config(filepath)
    host = _defaults['connection'].format(config['JOBSERVER'],
            config['JOBPORT'])

    socket.connect(host)
    request = {
        'request' : 'get_status',
        'data' : {
            'id' : args['id']
        },
    }

    poller = zmq.Poller()
    poller.register(socket, zmq.POLLIN)
    counter = itertools.count()

    socket.send_json(request, zmq.NOBLOCK)

    while _defaults['max_attempts'] > counter.next() and not result:

        if socket in dict(poller.poll(timeout=_defaults['timeout'])):
            result = socket.recv_json(flags=zmq.NOBLOCK)
        else:
            time.sleep(1)

            try:
                socket.send_json(request, zmq.NOBLOCK)
            except zmq.ZMQError:
                pass

    socket.close()
    context.term()

    if not result:
        result = {"error" : 1}

    return json.dumps(result)

urls = [
    (r"status/(?P<id>[a-z0-9\-]+)?", status),
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
