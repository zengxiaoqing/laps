#!/bin/sh

REMOTE_NODE=$1
REMOTE_DATA_ROOT=$2
DOMAIN_NAME=$3

echo "Running followup_ncarg.com on $REMOTE_NODE via do it yourself daemon...."

echo "command = /usr/nfs/common/lapb/www/followup_ncarg.com $REMOTE_DATA_ROOT $DOMAIN_NAME /usr/nfs/lapb/builds/laps"
echo           "/usr/nfs/common/lapb/www/followup_ncarg.com $REMOTE_DATA_ROOT $DOMAIN_NAME /usr/nfs/lapb/builds/laps" >> /exchange/tmp/usfsfire/queue.$REMOTE_NODE
chmod +x /exchange/tmp/usfsfire/queue.$REMOTE_NODE

echo "followup_ncarg_remote.sh finished..."

