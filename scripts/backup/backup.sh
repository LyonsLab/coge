#!/bin/bash

# Weekly cron job to backup to IRODS location
# Backed-up:
#    MySQL databases:
#       coge
#       wikidb
#       cogelinks
#    Data directories:
#       Wiki web directory
#       SVN code repo
#       CoGe data
# Created: 7/11/12 by mdb

echo `date` "Backup started"

ICMD=/usr/bin
#DAYS_UNTIL_DELETE=7
MAX_LOCAL_BACKUPS=1
MAX_REMOTE_BACKUPS=4
VERSION=$(date '+%Y%m%d')
LOCAL=/storage/coge/backup
REMOTE=backup
CONFIG=/opt/apache2/coge/coge.conf

#
# Dump databases and copy to remote IRODS location
#
DBPASS=$(grep -P "DBPASS" $CONFIG | awk '{print $2}')
LOCAL_MYSQL=$LOCAL/mysql_$VERSION
mkdir -p $LOCAL_MYSQL
echo `date` "Snapshotting MySQL databases"
#mysqldump -u root -p$DBPASS wikidb -c | gzip -9 > $LOCAL_MYSQL/wikidb.sql.gz
#mysqlhotcopy -u root -p $DBPASS cogelinks $LOCAL_MYSQL
#mysqlhotcopy --port=3307 -u root -p $DBPASS coge $LOCAL_MYSQL
innobackupex --user=root --password=$DBPASS --no-timestamp $LOCAL_MYSQL
innobackupex --apply-log $LOCAL_MYSQL
echo `date` "Pushing MySQL database backup to IRODS"
$ICMD/iput -bfr $LOCAL_MYSQL $REMOTE

#
# Remove old database backups
#
echo `date` "Deleting old LOCAL database backups"
#LOCAL_DELETIONS=`find $LOCAL/mysql_* -maxdepth 1 -type d -mtime +$DAYS_UNTIL_DELETE`
count=1
for d in `ls -td /storage/coge/backup/*/ | grep 'mysql_'` 
do
   if [ $count -gt $MAX_LOCAL_BACKUPS ];
   then
       echo deleting local $d 
       rm -rf $d
    fi
    ((count++))
done

echo `date` "Deleting old REMOTE database backups (local & remote)"
count=1
for d in `$ICMD/ils backup | grep 'mysql_' | sed 's/.*\(mysql_.*\)/\1/'`
do
    if [ $count -gt $MAX_REMOTE_BACKUPS ];
    then
        echo deleting IRODS backup/$d
        $ICMD/irm -r backup/$d
    fi
    ((count++))
done

#
# Sync data directories with remote IRODS location
#
echo `date` "Syncing data directories with IRODS"
$ICMD/icd
$ICMD/irsync -rs /opt/apache2/cogepedia i:$REMOTE/cogepedia
$ICMD/irsync -rs /opt/apache2/r i:$REMOTE/r
$ICMD/irsync -rs /storage/coge/data/genomic_sequence/ i:$REMOTE/genomic_sequence
$ICMD/irsync -rs /storage/coge/data/experiments/ i:$REMOTE/experiments

#
# Sync JEX database with remote IRODS location
#
echo `date` "Syncing JEX database with IRODS"
$ICMD/icd
$ICMD/irsync -rs /opt/Yerba/workflows.db i:$REMOTE/jex

#
# Sync /etc with remote IRODS location
#
echo `date` "Syncing /etc config files with IRODS"
$ICMD/icd
$ICMD/irsync -rs /etc i:$REMOTE/etc

echo `date` "Backup completed"
