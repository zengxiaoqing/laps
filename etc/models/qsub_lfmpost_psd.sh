#!/bin/sh
#$ -N lfmp_psd
#$ -A hmtb
#$ -l h_rt=04:30:00
#$ -S /bin/sh
#$ -cwd
#$ -pe comp 1
#exit

export LAPSROOT=/home/oplapb/builds/laps
export LAPS_DATA_ROOT=/lfs0/projects/hmtb/dwr_domains/laps_psd
export WRF_DATAROOT=/lfs0/projects/hmtb/dwr_domains/psd_gfs
export NETCDF=/opt/netcdf/3.6.2-pgi-7.1-3
export PATH=$PATH:/opt/netcdf/3.6.2-pgi-7.1-3/bin
export phys=fer

model=wrf
numFcsts=12
fcstIncrMin=60
maxWaitSec=3600
maxHrsRun=2
project=DWR

# next 4 lines for testing
#testTime=2008-10-02-1800
#dateDir=2008100218
#cd $WRF_DATAROOT/$dateDir/wrfprd
#/usr/bin/perl $LAPSROOT/etc/lfmpost.pl -m $model -f $numFcsts -i $fcstIncrMin -w $maxWaitSec -e $maxHrsRun -d $testTime -y fer -P $project -q

/usr/bin/perl $LAPSROOT/etc/lfmpost.pl -m $model -f $numFcsts -i $fcstIncrMin -w $maxWaitSec -e $maxHrsRun -y fer -P $project -q

exit 0

