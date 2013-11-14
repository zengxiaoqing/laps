#!/bin/sh 

# Set up links to locally run lfmpost on an external version of WRF

LOCAL_WRF_ROOT=$1
LAPS_DATA_ROOT=$2
REMOTE_WRF_ROOT=$3                      
DELAY=$4             # Hours

#LOCAL_WRF_ROOT=/lfs0/projects/hmtb/hwt_domains/hrrr_conus
#LAPS_DATA_ROOT=/pan1/projects/dlaps/analysis/wiap 
#REMOTE_WRF_ROOT=/home/rtrr/hrrr_dev1

# Purge local WRF_ROOT subdirectories
/usr/bin/perl /home/oplapb/builds/laps/etc/purge_w3_fcst2d.pl -f -t 72 -d $LOCAL_WRF_ROOT modelroot dummy

# Purge local WRF_ROOT logs                 
/usr/bin/perl /home/oplapb/builds/laps/etc/purger.pl -r -t 3.0 $LOCAL_WRF_ROOT/log                  

# Set model run time according to latest directory in REMOTE_WRF_ROOT
cd $REMOTE_WRF_ROOT
DATETIME=`ls -1 . | tail -2 | head -1`

# Alternatively set model run time by call to 'sched_sys.pl'
model_cycle_time=`/usr/bin/perl /home/oplapb/builds/laps/etc/read_nl.pl -d $LAPS_DATA_ROOT -n nest7grid.parms -v model_cycle_time`
DATETIME=`/usr/bin/perl /home/oplapb/builds/laps/etc/sched_sys.pl -c $model_cycle_time -d $DELAY -f yyyymmddhh`

echo "LAPS_DATA_ROOT is $LAPS_DATA_ROOT"  
echo "model_cycle_time is $model_cycle_time"
echo "DATETIME is $DATETIME"                      

# Short term fix to force a 00Z run
#DATETIME=`echo $DATETIME | cut -c1-8`00

mkdir -p $LOCAL_WRF_ROOT/log
 
LOCAL_WRF_RUN=$LOCAL_WRF_ROOT/$DATETIME
REMOTE_WRF_RUN=$REMOTE_WRF_ROOT/$DATETIME

mkdir -p $LOCAL_WRF_RUN/wrfprd
mkdir -p $LOCAL_WRF_RUN/static

cd $LOCAL_WRF_RUN/wrfprd
for file in `ls $REMOTE_WRF_RUN/wrfprd/wrfout*`; do
    ls -l $file
    ln -s $file .                              
done

echo ""
pwd
ls -l  

cd $LOCAL_WRF_RUN/static
ln -s $LAPS_DATA_ROOT/static/lfmpost.nl .

echo ""
pwd
ls -l 


