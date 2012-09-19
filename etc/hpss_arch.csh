#!/bin/csh

# Download archive data from hpss:
#
# Usage:
# Download with options: 
#   data_dir: where data will download to, 
#   yyyy, mm, dd: year, month and day
if ($1 == '' || $2 == '' || $3 == '' || $4 == '') then
  echo "Usage: hpss_arch.csh data_dir yyyy mm dd"
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

# ACAR:
set dir=acars/qc/netcdf
mkdir -p $dir
cd $dir
hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4\0000.tar.gz
gunzip $2$3$4\0000.tar.gz
tar xvf $2$3$4\0000.tar
rm -f $2$3$4\0000.tar
cd ../../../
echo "ACAR data downloaded"
pwd

# GPSmet:
set dir=gpsmet/netcdf
mkdir -p $dir
cd $dir
hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4\0000.tar.gz 
gunzip $2$3$4\0000.tar.gz
tar xvf $2$3$4\0000.tar
rm -f $2$3$4\0000.tar
cd ../../
echo "GPSmet data downloaded"
pwd

# Metar:
set dir=metar/netcdf
mkdir -p $dir
cd $dir
hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4\0000.tar.gz
gunzip $2$3$4\0000.tar.gz
tar xvf $2$3$4\0000.tar
rm -f $2$3$4\0000.tar
cd ../../
echo "Metar data downloaded"
pwd

# Radar nowrad:
set dir=radar/wsi/nowrad/netcdf
mkdir -p $dir
cd $dir
foreach hour (00 03 06 09 12 15 18 21)
  hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4$hour\00.tar.gz
  gunzip $2$3$4$hour\00.tar.gz
  tar xvf $2$3$4$hour\00.tar
  rm -f $2$3$4$hour\00.tar
end
echo "Here are the radar file gunzipped"
cd ../../../../
echo "radar data downloaded"
pwd

# Tower:
set dir=tower/public/netcdf
mkdir -p $dir
cd $dir
hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4\0000.tar.gz
gunzip $2$3$4\0000.tar.gz
tar xvf $2$3$4\0000.tar
rm -f $2$3$4\0000.tar
cd ../../../
pwd

# RAOB:
set dir=raob/netcdf
mkdir -p $dir
cd $dir
hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4\0000.tar.gz
gunzip $2$3$4\0000.tar.gz
tar xvf $2$3$4\0000.tar
rm -f $2$3$4\0000.tar
cd ../../
pwd

# PIREP:
set dir=pirep/netcdf
mkdir -p $dir
cd $dir
hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4\0000.tar.gz
gunzip $2$3$4\0000.tar.gz
tar xvf $2$3$4\0000.tar
rm -f $2$3$4\0000.tar
cd ../../
pwd

# SAO:
set dir=sao/netcdf
mkdir -p $dir
cd $dir
hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4\0000.tar.gz
gunzip $2$3$4\0000.tar.gz
tar xvf $2$3$4\0000.tar
rm -f $2$3$4\0000.tar
cd ../../
pwd

# MADIS:
set dir=madis/LDAD/profiler/netCDF
mkdir -p $dir
cd $dir
hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4\0000.tar.gz
gunzip $2$3$4\0000.tar.gz
tar xvf $2$3$4\0000.tar
rm -f $2$3$4\0000.tar
cd ../../../../
set dir=madis/LDAD/mesonet/netCDF
mkdir -p $dir
cd $dir
foreach hour (00 06 12 18)
  hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4$hour\00.tar.gz
  gunzip $2$3$4$hour\00.tar.gz
  tar xvf $2$3$4$hour\00.tar
  rm -f $2$3$4$hour\00.tar
end
cd ../../../../
set dir=madis/LDAD/urbanet/netCDF
mkdir -p $dir
cd $dir
hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4\0000.tar.gz
gunzip $2$3$4\0000.tar.gz
tar xvf $2$3$4\0000.tar
rm -f $2$3$4\0000.tar
cd ../../../../
set dir=madis/point/HDW/netcdf
mkdir -p $dir
cd $dir
hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4\0000.tar.gz
gunzip $2$3$4\0000.tar.gz
tar xvf $2$3$4\0000.tar
rm -f $2$3$4\0000.tar
cd ../../../../
set dir=madis/point/HDW1h/netcdf
mkdir -p $dir
cd $dir
hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4\0000.tar.gz
gunzip $2$3$4\0000.tar.gz
tar xvf $2$3$4\0000.tar
rm -f $2$3$4\0000.tar
cd ../../../../
set dir=madis/point/maritime/netcdf
mkdir -p $dir
cd $dir
hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4\0000.tar.gz
gunzip $2$3$4\0000.tar.gz
tar xvf $2$3$4\0000.tar
rm -f $2$3$4\0000.tar
cd ../../../../
set dir=madis/point/POES/netcdf
mkdir -p $dir
cd $dir
hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4\0000.tar.gz
gunzip $2$3$4\0000.tar.gz
tar xvf $2$3$4\0000.tar
rm -f $2$3$4\0000.tar
cd ../../../../
set dir=madis/point/radiometer/netcdf
mkdir -p $dir
cd $dir
hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4\0000.tar.gz
gunzip $2$3$4\0000.tar.gz
tar xvf $2$3$4\0000.tar
rm -f $2$3$4\0000.tar
cd ../../../../

# sat:
set dir=sat/fsl-conus/gs/netcdf
mkdir -p $dir
cd $dir
foreach hour (00 03 06 09 12 15 18 21)
  hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4$hour\00.tar.gz
  gunzip $2$3$4$hour\00.tar.gz
  tar xvf $2$3$4$hour\00.tar
  rm -f $2$3$4$hour\00.tar 
end
cd ../../../../

# maritime:
set dir=maritime/netcdf
mkdir -p $dir
cd $dir
hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4\0000.tar.gz
gunzip $2$3$4\0000.tar..gz
tar xvf $2$3$4\0000.tar.
rm -f $2$3$4\0000.tar.
cd ../../

# profiler:
set dir=profiler/rass/external/netcdf
mkdir -p $dir
cd $dir
hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4\0000.tar.gz
gunzip $2$3$4\0000.tar..gz
tar xvf $2$3$4\0000.tar.
rm -f $2$3$4\0000.tar.
cd ../../../../
set dir=profiler/rass/noaanet/netcdf
mkdir -p $dir
cd $dir
hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4\0000.tar.gz
gunzip $2$3$4\0000.tar..gz
tar xvf $2$3$4\0000.tar.
rm -f $2$3$4\0000.tar.
cd ../../../../
set dir=profiler/wind/external/netcdf
mkdir -p $dir
cd $dir
hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4\0000.tar.gz
gunzip $2$3$4\0000.tar..gz
tar xvf $2$3$4\0000.tar.
rm -f $2$3$4\0000.tar.
cd ../../../../
set dir=profiler/wind/noaanet/netcdf
mkdir -p $dir
cd $dir
hsi get /BMC/fdr/$2/$3/$4/data/$dir/$2$3$4\0000.tar.gz
gunzip $2$3$4\0000.tar..gz
tar xvf $2$3$4\0000.tar.
rm -f $2$3$4\0000.tar.
cd ../../../../

# Get RAP background:
set dir=grib/ftp_rap_iso/7/0/105/0_67725_30
mkdir -p grib/rap
cd grib/rap
foreach hour (00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23)
  hsi get /BMC/fdr/$2/$3/$4/$dir/$2$3$4$hour\00.tar.gz
  gunzip $2$3$4$hour\00.tar.gz
  tar xvf $2$3$4$hour\00.tar
  rm -f $2$3$4$hour\00.tar 
end
#/BMC/fdr/2012/06/29/grib/ftp_rap_iso/7/0/105/0_67725_30
cd $current

echo ""
