#!/bin/sh

REMOTE_NODE=$1
REMOTE_DATA_ROOT=$2
DOMAIN_NAME=$3

echo "start followup_ncarg_remote.sh to generate web images"
date

echo "Running followup_ncarg.com on $REMOTE_NODE via ssh..."

ssh $REMOTE_NODE /usr/nfs/common/lapb/www/followup_ncarg.com $REMOTE_DATA_ROOT $DOMAIN_NAME /usr/nfs/lapb/builds/laps

date

echo "followup_ncarg_remote.sh finished..."

