#!/bin/sh
#$ -N lfmp_5km
#$ -A hmtb
#$ -l h_rt=04:30:00
#$ -S /bin/sh
#$ -cwd
#$ -pe comp 1
#exit

export LAPSROOT=/home/oplapb/builds/laps
export LAPS_DATA_ROOT=/pan1/projects/mm5-laps/domains/LAPS5KM
export WRF_DATAROOT=/pan1/projects/mm5-laps/domains/WRFV3-5KM
export NETCDF=/opt/netcdf/3.6.2-pgi-7.1-3
export PATH=$PATH:/opt/netcdf/3.6.2-pgi-7.1-3/bin

model=wrf
numFcsts=36
fcstIncrMin=60
maxWaitSec=3600
maxHrsRun=4
project=DWR

# next 4 lines for testing
#testTime=2009-04-13-1200
#dateDir=2009041312
#maxWaitSec=0
#cd $WRF_DATAROOT/$dateDir/wrfprd
#/usr/bin/perl $LAPSROOT/etc/lfmpost.pl -m $model -f $numFcsts -i $fcstIncrMin -w $maxWaitSec -e $maxHrsRun -P $project -d $testTime -q

/usr/bin/perl $LAPSROOT/etc/lfmpost.pl -m $model -f $numFcsts -i $fcstIncrMin -w $maxWaitSec -e $maxHrsRun -P $project -q

exit 0

