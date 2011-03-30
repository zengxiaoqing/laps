#!/bin/sh

#Copy recent output products to remote destination

LOCAL_DATA_ROOT=$1
REMOTE_DATA_ROOT=$2

rsync -rlptgvv --rsh=ssh --delete $LOCAL_DATA_ROOT/lapsprd/www/* $REMOTE_DATA_ROOT/lapsprd/www > $LOCAL_DATA_ROOT/log/rsync_www.log.`date +\%H\%M` 2>&1

rsync -rlptgvv --exclude='log/core' --exclude='lapsprd/www' --exclude='lapsprd/lga' --rsh=ssh --delete $LOCAL_DATA_ROOT/* $REMOTE_DATA_ROOT > $LOCAL_DATA_ROOT/log/rsync.log.`date +\%H\%M` 2>&1         
