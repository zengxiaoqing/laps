#!/bin/sh
#$ -N lfmp_sch
#$ -A hmtb
#$ -l h_rt=04:30:00
#$ -S /bin/sh
#$ -cwd
#$ -pe comp 1
#exit

export LAPSROOT=/home/oplapb/builds/laps
export LAPS_DATA_ROOT=/lfs0/projects/hmtb/dwr_domains/laps
export WRF_DATAROOT=/lfs0/projects/hmtb/dwr_domains/wrf-schultz
export NETCDF=/opt/netcdf/3.6.2-pgi-7.1-3
export PATH=$PATH:/opt/netcdf/3.6.2-pgi-7.1-3/bin
export phys=sch

model=wrf
numFcsts=40
fcstIncrMin=180
maxWaitSec=3600
maxHrsRun=4
project=DWR

# next 4 lines for testing
#testTime=2008-10-10-2100
#dateDir=2008101021
#cd $WRF_DATAROOT/$dateDir/wrfprd
#/usr/bin/perl $LAPSROOT/etc/lfmpost.pl -m $model -f $numFcsts -i $fcstIncrMin -w $maxWaitSec -e $maxHrsRun -d $testTime -y -sch -P $project -q

/usr/bin/perl $LAPSROOT/etc/lfmpost.pl -m $model -f $numFcsts -i $fcstIncrMin -w $maxWaitSec -e $maxHrsRun -y sch -P $project -q

exit 0

