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
          static/ncarg/uscounty.dat \
          static/ncarg/state_from_counties.dat \
          static/ncarg/continent_minus_us.dat \
          lapsprd/pro/$a9time.pro \
          lapsprd/lso/$a9time.lso \
          lapsprd/lso/$a9time.lso_qc \
          lapsprd/snd/$a9time.snd \
          lapsprd/cdw/$a9time.cdw \
          lapsprd/pin/$a9time.pin \
          lapsprd/d01/$a9time.d01 \
          lapsprd/d02/$a9time.d02 \
          lapsprd/d03/$a9time.d03 \
          lapsprd/d04/$a9time.d04 \
          lapsprd/d05/$a9time.d05 \
          lapsprd/d06/$a9time.d06 \
          lapsprd/sag/$a9time.sag \
          lapsprd/pig/$a9time.pig \
          lapsprd/prg/$a9time.prg \
          lapsprd/lw3/$a9time.lw3 \
          lapsprd/lwm/$a9time.lwm \
          lapsprd/lsx/$a9time.lsx \
          lapsprd/tmg/$a9time.tmg \
          lapsprd/lt1/$a9time.lt1 \
          lapsprd/pbl/$a9time.pbl \
          lapsprd/vrz/$a9time.vrz \
          lapsprd/lps/$a9time.lps \
          lapsprd/lcb/$a9time.lcb \
          lapsprd/lcv/$a9time.lcv \
          lapsprd/lh3/$a9time.lh3 \
          lapsprd/lh4/$a9time.lh4 \
          lapsprd/lcp/$a9time.lcp \
          lapsprd/lwc/$a9time.lwc \
          lapsprd/lct/$a9time.lct \
          lapsprd/lco/$a9time.lco \
          lapsprd/lty/$a9time.lty \
          lapsprd/lst/$a9time.lst \
          lapsprd/lmr/$a9time.lmr \
          lapsprd/lmt/$a9time.lmt \
          lapsprd/liw/$a9time.liw \
          lapsprd/lhe/$a9time.lhe \
          lapsprd/lfr/$a9time.lfr \
          lapsprd/l1s/$a9time.l1s \
          lapsprd/lm2/$a9time.lm2 \
          log/sfc.wgi.$a9time \
          log/wind3d.wgi.$a9time \
          log/temp.wgi.$a9time \
          log/cloud.wgi.$a9time \
          log/hum3d.wgi.$a9time \
            | ssh brain.fsl.noaa.gov "cd $REMOTE_DATA_ROOT; tar xpvf -"

date

echo "Running followup_ncarg.com on toro via ssh..."

ssh brain ssh toro /usr/nfs/common/lapb/www/followup_ncarg.com $REMOTE_DATA_ROOT $DOMAIN_NAME /usr/nfs/lapb/builds/laps

date

echo "Running purger.pl on brain..."
echo "perl /usr/nfs/lapb/builds/laps/etc/purger.pl -r -m $NUM_FILES $REMOTE_DATA_ROOT/lapsprd"
ssh brain "perl /usr/nfs/lapb/builds/laps/etc/purger.pl -r -m $NUM_FILES $REMOTE_DATA_ROOT/lapsprd"

date

echo "followup_ncarg_remote.sh finished..."

