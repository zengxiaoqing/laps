#!/bin/sh

LOCAL_DATA_ROOT=$1
REMOTE_DATA_ROOT=$2
DOMAIN_NAME=$3
NUM_FILES=$4

echo "start followup_ncarg_remote.sh to copy lapsprd files to brain and generate web images"
date
echo $LOCAL_DATA_ROOT

a9time=`head -2 $LOCAL_DATA_ROOT/time/systime.dat | tail -1 | cut -c 2-10`

echo "Running tar to ship the $a9time files to brain..."

cd $LOCAL_DATA_ROOT
tar cvf - time/systime.dat \
          static/pressures.nl \
          static/lapsplot.nl \
          static/nest7grid.parms \
          static/static.nest7grid \
          lapsprd/pro/$a9time.pro \
          lapsprd/lso/$a9time.lso \
          lapsprd/snd/$a9time.snd \
          lapsprd/cdw/$a9time.cdw \
          lapsprd/pin/$a9time.pin \
          lapsprd/sag/$a9time.sag \
          lapsprd/pig/$a9time.pig \
          lapsprd/prg/$a9time.prg \
          lapsprd/lwm/$a9time.lwm \
          lapsprd/lsx/$a9time.lsx \
          lapsprd/tmg/$a9time.tmg \
          lapsprd/lt1/$a9time.lt1 \
          lapsprd/pbl/$a9time.pbl \
          lapsprd/lcb/$a9time.lcb \
          lapsprd/lcv/$a9time.lcv \
          lapsprd/lh4/$a9time.lh4 \
          lapsprd/lct/$a9time.lct \
          lapsprd/lst/$a9time.lst \
          lapsprd/lmr/$a9time.lmr \
          lapsprd/lmt/$a9time.lmt \
          lapsprd/liw/$a9time.liw \
          lapsprd/lhe/$a9time.lhe \
          lapsprd/lfr/$a9time.lfr \
          lapsprd/l1s/$a9time.l1s \
          lapsprd/lm2/$a9time.lm2 \
            | ssh brain.fsl.noaa.gov "cd $REMOTE_DATA_ROOT; tar xpvf -"

date

echo "Running followup_ncarg.com on jayhawk via ssh..."

ssh brain ssh jayhawk /usr/nfs/common/lapb/www/followup_ncarg.com $REMOTE_DATA_ROOT $DOMAIN_NAME /usr/nfs/lapb/builds/laps

date

echo "Running purger.pl on brain..."
echo "perl /usr/nfs/lapb/builds/laps/etc/purger.pl -r -m $NUM_FILES $REMOTE_DATA_ROOT/lapsprd"
ssh brain "perl /usr/nfs/lapb/builds/laps/etc/purger.pl -r -m $NUM_FILES $REMOTE_DATA_ROOT/lapsprd"

