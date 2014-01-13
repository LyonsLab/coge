#!/usr/bin/env python
import argparse
import itertools
import os
import sys
import time

import zmq
import MySQLdb as mdb

_defaults = {
    'connection' : 'tcp://localhost:5151',
    'max_attempts' : 30,
    'timeout' : 1000,
}

def db_connect(config):
    path = os.path.join(config)
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

def get_running_jobs(db):
    cursor = db.cursor()

    query = 'select job_id from job where status=1 and not page="gevo";'

    cursor.execute(query)
    return cursor.fetchall()

def fixup_jobs(db, update=False):
    cursor = db.cursor()

    print("\n" + "#" * 25)
    print("Fixing Jobs")
    print("#" * 25)

    # Fixup all jobs that are not presently in the queue
    query = "update job set end_time=NOW() where (start_time > end_time) and (status > 1)"

    # Get a list of the jobs to be fixed
    cursor.execute("select job_id from job where (start_time > end_time) and (status > 1)")

    jobs = cursor.fetchall()

    for job in jobs:
        print("Job {0} will be fixed".format(job[0]))

    if not jobs:
        print("No fixes necessary.")
        return

    if update:
        print("Executing: {0}".format(query))
        cursor.execute(query)
    else:
        print("[DRY-RUN]: {0}".format(query))


def update_jobs(db, jobs, update=False):
    query = "update job set status={status_code}, end_time=NOW() where job_id={id};"
    print("\n" + "#" * 25)
    print("Updating Jobs")
    print("#" * 25)

    if not jobs:
        print("No updates found.")

    cursor = db.cursor()

    for job in jobs:
        q1 = query.format(**job)
        if update:
            print("Executing: {0}".format(q1))
            cursor.execute(query.format(**job))
        else:
            print("[DRY-RUN]: {0}".format(q1))

def main(args={}):
    context = zmq.Context()
    socket = context.socket(zmq.REQ)
    socket.setsockopt(zmq.LINGER, 0)
    socket.connect(_defaults['connection'])

    pid = str(os.getpid())
    pidfile = "/tmp/update-jobs.pid"

    try:
        current = file(pidfile, 'r').readline()

        if os.path.exists("/proc/%s" % current):
            sys.exit()
        else:
            os.unlink(pidfile)
    except IOError:
        pass
    finally:
        file(pidfile, 'w').write(pid)

    print("#" * 25)
    print("Fetching Jobs")
    print("#" * 25)

    db = db_connect(args.config)
    jobs = get_running_jobs(db)

    poller = zmq.Poller()
    poller.register(socket, zmq.POLLIN)

    updated_jobs = []
    print("Found {0} running jobs.".format(len(jobs)))

    for job in jobs:
        request = {
            'request' : 'get_status',
            'data' : { 'id' : str(job[0]) }
        }

        socket.send_json(request, zmq.NOBLOCK)
        result = None
        counter = itertools.count()

        while _defaults['max_attempts'] > counter.next() and not result:

            if socket in dict(poller.poll(timeout=_defaults['timeout'])):
                result = socket.recv_json(flags=zmq.NOBLOCK)
            else:
                time.sleep(1)
                try:
                    socket.send_json(request, zmq.NOBLOCK)
                except zmq.ZMQError:
                    pass

        if not result:
            print("Unable to connect to the job engine.")
            sys.exit(1)

        job_status = result['status'].lower()

        if job_status  == 'running':
            print("The job {id} is still running.".format(id=job[0]))
            continue
        elif job_status == 'completed':
            code = 2
        elif job_status == 'cancelled':
            code = 3
        elif job_status == 'notfound':
            code = 4
            job_status = 'terminated'
        elif job_status == 'failed':
            code = 5
        else:
            continue

        info= { 'id' : job[0], 'status_code' : code, 'status' : job_status}
        print("Job {id} will be updated to status {status}.".format(**info))
        updated_jobs.append(info)

    update_jobs(db, updated_jobs, update=args.update)

    if args.fix:
        fixup_jobs(db, update=args.update)

    socket.close()
    context.term()
    db.close()
    os.unlink(pidfile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Update jobs in the job engine.  Run with no options for a dry-run with no database commits.')

    parser.add_argument('config', help="Database configuration file")
    parser.add_argument('--update', action='store_true', help="Commit the updates to the database.  Must be set to have commits happen.")
    parser.add_argument('--fix', action='store_true', help="Fixes a problem with time-stamps.")
    main(parser.parse_args())
