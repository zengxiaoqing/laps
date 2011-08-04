#!/bin/sh 

# Set up links to locally run lfmpost on an external version of WRF

LOCAL_WRF_ROOT=/lfs0/projects/hmtb/hwt_domains/hrrr_conus
REMOTE_WRF_ROOT=/whome/rtrr/hrrr
LAPS_DATA_ROOT=/pan1/projects/dlaps/analysis/wiap 

# Set model run time according to latest directory in REMOTE_WRF_ROOT
cd $REMOTE_WRF_ROOT
DATETIME=`ls -1 . | tail -2 | head -1`

# Alternatively set model run time by call to 'sched_sys.pl'

# Alternatively pass in model run time on command line

# Short term fix to force a 00Z run
DATETIME=`echo $DATETIME | cut -c1-8`00
 
LOCAL_WRF_RUN=$LOCAL_WRF_ROOT/$DATETIME
REMOTE_WRF_RUN=$REMOTE_WRF_ROOT/$DATETIME

mkdir -p $LOCAL_WRF_RUN/wrfprd
mkdir -p $LOCAL_WRF_RUN/static

cd $LOCAL_WRF_RUN/wrfprd
ln -s $REMOTE_WRF_RUN/wrfprd/wrfout* .

echo ""
pwd
ls -l  

cd $LOCAL_WRF_RUN/static
ln -s $LAPS_DATA_ROOT/static/lfmpost.nl .

echo ""
pwd
ls -l 


