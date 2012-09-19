#!/bin/csh

# Download archive data from hpss:
#
# Usage:
# Download with options: 
#   data_dir: where data will download to, 
#   yyyy, mm, dd: year, month and day
if ($1 == '' || $2 == '' || $3 == '' || $4 == '') then
  echo "Usage: hpss_arch.csh data_dir yyyy mm dd"
  echo "Example: ./hpss_model_nam.sh /scratch2/portfolios/BMC/odvars/derecho/data 2012 06 28"
  exit
endif

! current directory:
set current=`pwd`
# Switch to the data_dir directory:
cd $1
echo "Where we are now: " `pwd`

echo "Downloading: to $1 for $2 $3 $4"

# Load HPSS module:
module load hpss
echo `module list`


# Get NAM A218 background:
set dir=grib/ftp/7/0/84/0_262792_30
set modeldir=grib/nam
mkdir -p $modeldir
cd $modeldir
hsi ls -1 /BMC/fdr/$2/$3/$4/$dir >& $current/tmp.txt
grep gz $current/tmp.txt >$current/filename.txt
foreach file (`cat $current/filename.txt`)
  echo 'download' $file
  hsi get $file
  tar xvfz *.tar.gz
  rm -f *.tar.gz
end
#/scratch2/portfolios/BMC/public/data/grids/nam/A218   grib2 -> ../../../grib/ftp/7/0/84/0_262792_30/
#/BMC/fdr/2012/06/29/grib/ftp/7/0/84/0_262792_30

cd $current

echo ""
