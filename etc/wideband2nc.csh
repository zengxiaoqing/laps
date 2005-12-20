#!/bin/csh

#Define Casedate (and hour within)
setenv LAPS_DATA_ROOT $1
setenv SKIP $2                 # skip processing if NetCDF radar is there?
setenv REMAP $3                # run remap_polar_netcdf.exe [yes,no]
setenv LAPSINSTALLROOT $4
setenv MODETIME $5 

setenv SUFFIX "_elev01"

echo "Start wideband2nc.csh..."

setenv TZ GMT
date

#Only needed for archive cases possibly
setenv YEAR `head -5 $LAPS_DATA_ROOT/time/systime.dat | tail -1 | cut -c8-11`
setenv MONTH 09
setenv DATE `head -5 $LAPS_DATA_ROOT/time/systime.dat | tail -1 | cut -c1-2`
setenv HOUR `head -3 $LAPS_DATA_ROOT/time/systime.dat | tail -1`
setenv YYDDD `head -6 $LAPS_DATA_ROOT/time/systime.dat | tail -1`

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

setenv INPUTROOT `head -2 $LAPS_DATA_ROOT/static/widebandlist.txt | tail -1`
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

# Output Directories
  mkdir -p $OUTPUTROOT
  mkdir -p $OUTPUTROOT/$RADAR
  mkdir -p $OUTPUTROOT/$RADAR/log
  mkdir -p $OUTPUTROOT/$RADAR/netcdf
  mkdir -p $OUTPUTROOT/$RADAR/netcdf/.tmp
# chmod -f u+w $OUTPUTROOT/$RADAR/netcdf/.tmp/*
# rm -f $OUTPUTROOT/$RADAR/netcdf/.tmp/*

# Log Directories
# mkdir -p $INSTALLROOT/$RADAR
# mkdir -p $INSTALLROOT/$RADAR/log

  echo " "
  echo "Look for logs in $OUTPUTROOT/$RADAR/log..."
  echo "Look for output in $OUTPUTROOT/$RADAR/netcdf..."
  echo " "

# Test for existence of output data already present?
  echo "ls -1 $OUTPUTROOT/$RADAR/netcdf/$YYDDD$HOUR* | wc -l"
  setenv COUNT `ls -1 $OUTPUTROOT/$RADAR/netcdf/$YYDDD$HOUR* | wc -l`
  echo "Number of pre-existing files is $COUNT"

  if($COUNT == "0" || $SKIP != "yes")then
#     Use this filename convention for realtime data

      if($MODETIME == "realtime")then

#         We assume that less than 24 hours of realtime Archive-II data are available on disk at any given time

          echo "Generating radar $RADAR in realtime..."

          pushd $INPUTROOT/$RADAR
              foreach file (*)
                  echo " "
                  echo "processing Archive-II file $file"

                  setenv HOUR   `echo $file | cut -c9-10`
                  setenv MINUTE `echo $file | cut -c11-12`

                 
#                 find /$OUTPUTROOT/$RADAR -name "$OUTPUTROOT/$RADAR/netcdf/*$HOUR$MINUTE$SUFFIX" -exec ls -1 "{}" \;

#                 ls -l $OUTPUTROOT/$RADAR/netcdf/*$HOUR$MINUTE$SUFFIX
#                 ls -1 $OUTPUTROOT/$RADAR/netcdf/*$HOUR$MINUTE$SUFFIX | wc -l

                  setenv OUTCOUNT `ls -1 $OUTPUTROOT/$RADAR/netcdf/*$HOUR$MINUTE$SUFFIX | wc -l`

                  echo "OUTCOUNT = "$OUTCOUNT

#                 ls -al --time-style=+%Y%m | awk '{print $6}'

#                 if (-e '$OUTPUTROOT/$RADAR/netcdf/*$HOUR$MINUTE$SUFFIX') then
#                 if ($exist == "yes") then

                  if($OUTCOUNT != "0")then
                      echo "output file already exists for     $HOUR$MINUTE$SUFFIX"
                  else
                      echo "output file does not yet exist for $HOUR$MINUTE$SUFFIX"
                      $INSTALLROOT/bin/TarNexrad2NetCDF -l $OUTPUTROOT/$RADAR/log \
                                              -p $RADAR -o $OUTPUTROOT/$RADAR/netcdf -c $INSTALLROOT/cdl/wsr88d_wideband.cdl -r $RADAR -t $INPUTROOT/$RADAR/$file
                  endif

#                 find /$INPUTROOT/$RADAR -name "$INPUTROOT/$RADAR/$file*" -exec $INSTALLROOT/bin/TarNexrad2NetCDF -l $OUTPUTROOT/$RADAR/log \
#                                         -p $RADAR -o $OUTPUTROOT/$RADAR/netcdf -c $INSTALLROOT/cdl/wsr88d_wideband.cdl -r $RADAR -t {} \;
              end
          popd
   
          echo "listing of output files for this hour "$YYDDD$HOUR
          echo "ls -1 $OUTPUTROOT/$RADAR/netcdf/$YYDDD$HOUR* | wc -l"

      else # archive case
#         This works with the archived radar data for IHOP
          echo "Generating radar $RADAR at $YEAR $MONTH $DATE $HOUR..."

          find /$INPUTROOT/$RADAR -name "$YEAR$MONTH$DATE$HOUR*$RADAR*" -exec $INSTALLROOT/bin/TarNexrad2NetCDF -l $OUTPUTROOT/$RADAR/log \
                                  -p $RADAR -o $OUTPUTROOT/$RADAR/netcdf -c $INSTALLROOT/cdl/wsr88d_wideband.cdl -r $RADAR -t {} \;

          echo "ls -1 $OUTPUTROOT/$RADAR/netcdf/$YYDDD$HOUR* | wc -l"
          setenv COUNT `ls -1 $OUTPUTROOT/$RADAR/netcdf/$YYDDD$HOUR* | wc -l`
          echo "Number of files generated is $COUNT"
          echo "Finished radar $RADAR at $YEAR $MONTH $DATE $HOUR..."

      endif

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





