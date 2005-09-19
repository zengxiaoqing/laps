#!/bin/csh

#Define Casedate (and hour within)
setenv YEAR $1
setenv MONTH $2
setenv DATE $3
setenv HOUR $4
setenv YYDDD $5
setenv LAPS_DATA_ROOT $6
setenv SKIP $7                 # skip processing if NetCDF radar is there?
setenv REMAP $8                # run remap_polar_netcdf.exe [yes,no]
setenv LAPSINSTALLROOT $9
setenv INPUTROOT $10           # location of Archive II data
                               # e.g. /data/ihop/$MONTH/$DATE/data/radar/wsr88d/level2

echo "Start wideband2nc.csh..."
date

#setenv OUTPUTROOT  /scratch/lapb/albers/radar
#setenv OUTPUTROOT  /data/lapb/ihop_work/raw/wsr88d/wideband

# Location of 'widebandlist.txt'
if (! -e $LAPS_DATA_ROOT/static/widebandlist.txt) then
    echo "ERROR: data file $LAPS_DATA_ROOT/static/widebandlist.txt not found..."
    exit
endif

setenv INSTALLROOT `head -1 $LAPS_DATA_ROOT/static/widebandlist.txt`
if (! -e $INSTALLROOT) then
    echo "ERROR: WIDEBAND INSTALLROOT $INSTALLROOT not found..."
    exit
endif

setenv OUTPUTROOT  $LAPS_DATA_ROOT/lapsprd/rdr/wideband

echo " "

echo "INPUTROOT=$INPUTROOT"
echo "INSTALLROOT=$INSTALLROOT"
echo "OUTPUTROOT=$OUTPUTROOT"

setenv NEXRAD_SITES $INSTALLROOT/NexradSite.cfg
if (! -e $NEXRAD_SITES) then
    echo "ERROR: $NEXRAD_SITES not found..."
    exit
endif

#Check and skip over any radars that were already processed for the hour 

foreach RADAR (`tail -1 $LAPS_DATA_ROOT/static/widebandlist.txt`)
#foreach RADAR (kama kcys kddc kfdr kftg kfws kgld kict kinx klbb ktlx kvnx)
#foreach RADAR (kcys kama)

  echo " "
  echo "Start radar $RADAR at $YEAR $MONTH $DATE $HOUR..."

# Output Directories
  mkdir -p $OUTPUTROOT
  mkdir -p $OUTPUTROOT/$RADAR
  mkdir -p $OUTPUTROOT/$RADAR/log
  mkdir -p $OUTPUTROOT/$RADAR/netcdf
  mkdir -p $OUTPUTROOT/$RADAR/netcdf/.tmp
# chmod -f u+w $OUTPUTROOT/$RADAR/netcdf/.tmp/*
# rm -f $OUTPUTROOT/$RADAR/netcdf/.tmp/*

# Log Directories
  mkdir -p $INSTALLROOT/$RADAR
  mkdir -p $INSTALLROOT/$RADAR/log

  echo " "
  echo "Look for logs in $OUTPUTROOT/$RADAR/log..."
  echo "Look for output in $OUTPUTROOT/$RADAR/netcdf..."
  echo " "

# Test for existence of output data already present?
  echo "ls -1 $OUTPUTROOT/$RADAR/netcdf/$YYDDD$HOUR* | wc -l"
  setenv COUNT `ls -1 $OUTPUTROOT/$RADAR/netcdf/$YYDDD$HOUR* | wc -l`
  echo "Number of pre-existing files is $COUNT"

  if($COUNT == "0" || $SKIP != "yes")then
      echo "Generating radar $RADAR at $YEAR $MONTH $DATE $HOUR..."

#     This works with the archived radar data for IHOP
#     find /$INPUTROOT/$RADAR -name "$YEAR$MONTH$DATE$HOUR*$RADAR*" -exec $INSTALLROOT/bin/TarNexrad2NetCDF -l $OUTPUTROOT/$RADAR/log \
#                                                 -p $RADAR -o $OUTPUTROOT/$RADAR/netcdf -c $INSTALLROOT/cdl/wsr88d_wideband.cdl -t {} \;

#     Use this filename convention for realtime data
      ls -l $INPUTROOT/$RADAR/$YEAR$MONTH$DATE$HOUR*
      find /$INPUTROOT/$RADAR -name "$YEAR$MONTH$DATE$HOUR*" -exec $INSTALLROOT/bin/TarNexrad2NetCDF -l $OUTPUTROOT/$RADAR/log \
                                                 -p $RADAR -o $OUTPUTROOT/$RADAR/netcdf -c $INSTALLROOT/cdl/wsr88d_wideband.cdl -t {} \;

      echo "ls -1 $OUTPUTROOT/$RADAR/netcdf/$YYDDD$HOUR* | wc -l"
      setenv COUNT `ls -1 $OUTPUTROOT/$RADAR/netcdf/$YYDDD$HOUR* | wc -l`
      echo "Number of files generated is $COUNT"
      echo "Finished radar $RADAR at $YEAR $MONTH $DATE $HOUR..."
  else
      echo "Pre-existing output: skipped radar $RADAR at $YEAR $MONTH $DATE $HOUR..."
  endif

  date
  
  echo " "

  if($REMAP == "yes")then
      echo "Running LAPS remapper"
      $LAPSINSTALLROOT/bin/remap_polar_netcdf.exe > $LAPS_DATA_ROOT/log/remap_polar_netcdf.log.$YYDDD$HOUR
  endif
end





