#!/bin/csh

#Define Casedate (and hour within)
setenv YEAR $1
setenv MONTH $2
setenv DATE $3
setenv HOUR $4
setenv YYDDD $5

setenv INPUTROOT   /data/ihop/$MONTH/$DATE/data/radar/wsr88d/level2

echo "Start wideband2nc.csh..."

if (-e /home/lapb/albers/bin/ihop) then
    setenv INSTALLROOT /home/lapb/albers/bin/ihop
else if (-e /home/albers/bin/ihop) then
    setenv INSTALLROOT /home/albers/bin/ihop
else
    echo "ERROR: WIDEBAND INSTALLROOT $INSTALLROOT not found..."
    exit
endif

#setenv OUTPUTROOT  /scratch/lapb/albers/radar
#setenv OUTPUTROOT  /data/lapb/ihop_work/raw/wsr88d/wideband
setenv OUTPUTROOT  $6

echo " "
date
echo "Start wideband2nc.csh..."
echo "INPUTROOT=$INPUTROOT"
echo "INSTALLROOT=$INSTALLROOT"
echo "OUTPUTROOT=$OUTPUTROOT"

setenv NEXRAD_SITES $INSTALLROOT/NexradSite.cfg
if (! -e $NEXRAD_SITES) then
    echo "ERROR: $NEXRAD_SITES not found..."
    exit
endif

#Check and skip over any radars that were already processed for the hour 
setenv SKIP $7 

foreach RADAR (kama kcys kddc kfdr kftg kfws kgld kict kinx klbb ktlx kvnx)
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
      find /$INPUTROOT/$RADAR -name "$YEAR$MONTH$DATE$HOUR*$RADAR*" -exec $INSTALLROOT/bin/TarNexrad2NetCDF -l $OUTPUTROOT/$RADAR/log \
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

end





