#!/bin/bash
#------------------------------------------------------------------------------
# Script to start/stop Mojolicious development web server (morbo) for this sandbox.
#------------------------------------------------------------------------------

start() {
    export PERL5LIB=./modules/perl
#    hypnotoad -f ./web/services/api.pl >> ./api.log 2>&1 &
#    ./web/services/api.pl daemon -l http://localhost:3304 >> ./api.log 2>&1 &
    morbo -l http://localhost:3304 ./web/services/api.pl >> ./api.log 2>&1 &
    echo "Started API"
}

stop() {
#    hypnotoad ./web/services/api.pl --stop
    pid=$(pgrep 'morbo')
    kill -9 $pid
    echo "Stopped API (pid " $pid ")"     
}

case "$1" in
    start)
        start
        exit 0
        ;;
    stop)
        stop
        exit 0
        ;;
    restart)
        stop
        sleep 1
        start
        exit 0
        ;;
    *)
        echo "Usage: $0 {start|stop|restart}" 1>&2
        exit 1
        ;;
esac
