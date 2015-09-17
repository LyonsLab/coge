#!/bin/bash

REMOTE_DIR="/iplant/home/coge/backup/genomic_sequence"
LOCAL_DIR="/storage2/coge/data/genomic_sequence"
/usr/local/bin/irsync -rs i:$REMOTE_DIR $LOCAL_DIR

REMOTE_DIR="/iplant/home/coge/backup/experiments"
LOCAL_DIR="/storage2/coge/data/experiments"
/usr/local/bin/irsync -rs i:$REMOTE_DIR $LOCAL_DIR

