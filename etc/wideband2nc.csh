#!/bin/csh

#Convert Archive-II radar files into Polar NetCDF using the GSD/ITS supplied software package.

#In realtime mode, this script can be run in cron using 'laps_driver.pl'.
#A cron line can be constructed that is analogous to this: 
#00,05,10,15,20,25,30,35,40,45,50,55   * * * *  /usr/bin/perl /usr/nfs/lapb/builds/laps/etc/laps_driver.pl wideband2nc.csh /usr/nfs/lapb/builds/laps /data/lapb/projects/smg/rt 

#In archive mode, this script can be run using 'casererun.pl' with appropriate command line arguments.  

#This script reads additional info from $LAPS_DATA_ROOT/static/widebandlist.txt

#Define Casedate (and hour within)
setenv LAPS_DATA_ROOT $1
setenv SKIP $2                 # skip processing if NetCDF radar is there?
setenv REMAP $3                # run remap_polar_netcdf.exe [yes,no]
setenv LAPSINSTALLROOT $4
setenv MODETIME $5             # run mode [realtime,archive]
setenv OUTPUTROOT_ARCHIVE $6   # location of output Polar NetCDF files (active only for archive mode)

setenv SUFFIX "_elev01"

echo "Start wideband2nc.csh..."

setenv TZ GMT
date

#setenv OUTPUTROOT  /scratch/lapb/albers/radar
#setenv OUTPUTROOT  /data/lapb/ihop_work/raw/wsr88d/wideband

# Check location of 'widebandlist.txt'
if (! -e $LAPS_DATA_ROOT/static/widebandlist.txt) then
    echo "ERROR: data file $LAPS_DATA_ROOT/static/widebandlist.txt not found..."
    exit
endif

#Obtain location of installed software that converts radar files
setenv INSTALLROOT `head -1 $LAPS_DATA_ROOT/static/widebandlist.txt`
if (! -e $INSTALLROOT) then
    echo "ERROR: WIDEBAND software INSTALLROOT $INSTALLROOT not found..."
    exit
endif

#Obtain location of Archive-II data files
setenv INPUTROOT `head -2 $LAPS_DATA_ROOT/static/widebandlist.txt | tail -1`
if (! -e $INPUTROOT) then
    echo "ERROR: WIDEBAND INPUTROOT data directory $INPUTROOT not found..."
    exit
endif

if ($MODETIME == "realtime") then
    setenv OUTPUTROOT  $LAPS_DATA_ROOT/lapsprd/rdr/wideband
    setenv HR1 `head -3 $LAPS_DATA_ROOT/time/systime.dat | tail -1`
    set HR2=$HR1
    @ HR2++
    if ($HR2 == 24) then
        set HR2=1
    endif
    echo "Processing hours $HR1 and $HR2"

#   Name of executable that converts from Nexrad format to NetCDF
    setenv NEXRAD_2_NETCDF TarNexrad2NetCDF

else # archive case

#   Access additional command line args
    setenv YEAR $7
    setenv MONTH $8
    setenv DATE $9
    setenv HOUR $10
    setenv YYDDD $11

    echo "YEAR=$YEAR"
    echo "DATE=$DATE"
    echo "HOUR=$HOUR"
    echo "YYDDD=$YYDDD"
    echo "MONTH=$MONTH"

#   This could be done if 'systime.dat' is updated actively by 'casererun.pl'
#   setenv YEAR `head -5 $LAPS_DATA_ROOT/time/systime.dat | tail -1 | cut -c8-11`
#   setenv DATE `head -5 $LAPS_DATA_ROOT/time/systime.dat | tail -1 | cut -c1-2`
#   setenv HOUR `head -3 $LAPS_DATA_ROOT/time/systime.dat | tail -1`
#   setenv YYDDD `head -6 $LAPS_DATA_ROOT/time/systime.dat | tail -1`

#   Need to run perl script to convert YEAR and DDD to MONTH
#   setenv MONTH 09 

    setenv OUTPUTROOT  $OUTPUTROOT_ARCHIVE
#   setenv OUTPUTROOT  $LAPS_DATA_ROOT/lapsprd/rdr/wideband

#   Name of executable that converts from Nexrad format to NetCDF
    setenv NEXRAD_2_NETCDF NCDC_TarNexrad2NetCDF

endif

echo " "

echo "INPUTROOT=$INPUTROOT"
echo "INSTALLROOT=$INSTALLROOT"
echo "OUTPUTROOT=$OUTPUTROOT"

setenv NEXRAD_SITES $INSTALLROOT/NexradSite.cfg
if (! -e $NEXRAD_SITES) then
    echo "ERROR: $NEXRAD_SITES not found..."
    exit
endif

if (! -e $INSTALLROOT/bin/$NEXRAD_2_NETCDF) then
    echo "ERROR: $INSTALLROOT/bin/$NEXRAD_2_NETCDF not found..."
    echo "Check to see if executable should be renamed to $NEXRAD_2_NETCDF"
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

# Use this filename convention for realtime data

  if ($MODETIME == "realtime") then

#     We assume that less than 24 hours of realtime Archive-II data are available on disk at any given time

      echo "Generating radar $RADAR in realtime..."

      pushd $INPUTROOT/$RADAR
          foreach file (*)
              echo " "
              echo "checking Archive-II file $file"

              setenv HOUR   `echo $file | cut -c9-10`
              setenv MINUTE `echo $file | cut -c11-12`
                 
#             echo "testing time $HOUR $HR1 $HR2"

              if ($HOUR == "$HR1" || $HOUR == "$HR2") then
#                 echo "matching time $HOUR $HR1 $HR2"
                  setenv OUTCOUNT `ls -1 $OUTPUTROOT/$RADAR/netcdf/*$HOUR$MINUTE$SUFFIX | wc -l`

                  echo "OUTCOUNT = "$OUTCOUNT

                  if ($OUTCOUNT != "0") then
                      echo "output file $HOUR$MINUTE$SUFFIX already exists"
                  else
                      echo "processing output file $HOUR$MINUTE$SUFFIX that does not yet exist"
                      $INSTALLROOT/bin/$NEXRAD_2_NETCDF -l $OUTPUTROOT/$RADAR/log \
                                              -p $RADAR -o $OUTPUTROOT/$RADAR/netcdf -c $INSTALLROOT/cdl/wsr88d_wideband.cdl -r $RADAR -t $INPUTROOT/$RADAR/$file
                  endif

              else
#                 echo "no matching time $HOUR $HR1 $HR2"
                  echo "outside time window $HOUR$MINUTE$SUFFIX"

              endif
          end

#         Cleanup
          rm -rf $OUTPUTROOT/$RADAR/netcdf/.tmp

      popd

  else # archive case
#     Test for existence of output data already present?
      echo "listing of output files for this hour "$YYDDD$HOUR
      echo "ls -1 $OUTPUTROOT/$RADAR/netcdf/$YYDDD$HOUR* | wc -l"

      setenv COUNT `ls -1 $OUTPUTROOT/$RADAR/netcdf/$YYDDD$HOUR* | wc -l`
      echo "Number of pre-existing files is $COUNT"

      if ($COUNT == "0" || $SKIP != "yes") then
#         This works with the archived radar data (e.g. for IHOP)
          echo "Generating radar $RADAR at $YEAR $MONTH $DATE $HOUR..."

#         Filename convention for /public
#         find /$INPUTROOT/$RADAR -name "$YEAR$MONTH$DATE$HOUR*$RADAR*" -exec $INSTALLROOT/bin/$NEXRAD_2_NETCDF -l $OUTPUTROOT/$RADAR/log \
#                                 -p $RADAR -o $OUTPUTROOT/$RADAR/netcdf -c $INSTALLROOT/cdl/wsr88d_wideband.cdl -r $RADAR -t {} \;

#         Filename convention for NCDC
          ls -l $INPUTROOT/$RADAR
          ls -l $INPUTROOT/$RADAR/*$YEAR$MONTH$DATE*
          find /$INPUTROOT/$RADAR -name "*$YEAR$MONTH$DATE*" -exec $INSTALLROOT/bin/$NEXRAD_2_NETCDF -l $OUTPUTROOT/$RADAR/log \
                                  -p $RADAR -o $OUTPUTROOT/$RADAR/netcdf -c $INSTALLROOT/cdl/wsr88d_wideband.cdl -r $RADAR -t {} \;

          echo "ls -1 $OUTPUTROOT/$RADAR/netcdf/$YYDDD$HOUR* | wc -l"
          setenv COUNT `ls -1 $OUTPUTROOT/$RADAR/netcdf/$YYDDD$HOUR* | wc -l`
          echo "Number of files generated is $COUNT"
          echo "Finished radar $RADAR at $YEAR $MONTH $DATE $HOUR..."

      else
          echo "Pre-existing output: skipped radar $RADAR at $YEAR $MONTH $DATE $HOUR..."

      endif

  endif

  date

  echo " "

end

if ($REMAP == "yes") then
    echo "Running LAPS remapper"

    if ($MODETIME == "realtime") then
        setenv TIMESTAMP `date +%y%j%H%M`
        $LAPSINSTALLROOT/bin/remap_polar_netcdf.exe > $LAPS_DATA_ROOT/log/remap_polar_netcdf.log.$TIMESTAMP
        echo "remap log file is remap_polar_netcdf.log.$TIMESTAMP"

    else # archive case
        $LAPSINSTALLROOT/bin/remap_polar_netcdf.exe > $LAPS_DATA_ROOT/log/remap_polar_netcdf.log.$YYDDD$HOUR
        echo "remap log file is remap_polar_netcdf.log.$YYDDD$HOUR"

    endif
endif




