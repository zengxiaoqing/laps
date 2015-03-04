#!/bin/csh

#Convert Archive-II radar files into Polar NetCDF using the GSD/ITS supplied software package.

#In realtime mode, this script can be run in cron using 'laps_driver.pl'.
#A cron line can be constructed that is analogous to this: 
#00,05,10,15,20,25,30,35,40,45,50,55   * * * *  /usr/bin/perl /usr/nfs/lapb/builds/laps/etc/laps_driver.pl wideband2nc.csh /usr/nfs/lapb/builds/laps /data/lapb/projects/smg/rt 

#In archive mode, this script can be run using 'casererun.pl' with appropriate command line arguments.

#NCDC files should already be untarred, though it's OK if they are compressed as .Z files.

#Bring in command line arguments
setenv LAPS_DATA_ROOT $1
setenv SKIP $2                 # skip processing if NetCDF radar is there?
setenv REMAP $3                # run remap_polar_netcdf.exe [yes,no]
setenv LAPSINSTALLROOT $4
setenv MODETIME $5             # run mode [realtime,archive]
setenv OUTPUTROOT_ARCHIVE $6   # location of output Polar NetCDF files (active only for archive mode)

#This script reads additional info from $LAPS_DATA_ROOT/static/widebandlist.txt as follows:

#First line is location of installed software that converts radar files
#Second line is location of input Archive-II data files 
#Third line is a list of radars to process. Use lower case for real-time and upper case for archive (NCDC).

setenv SUFFIX "_elev01"

echo "Start wideband2nc.csh..."

setenv TZ GMT
date

# Check location of 'widebandlist.txt'
if (! -e $LAPS_DATA_ROOT/static/widebandlist.txt) then
    echo "ERROR: data file $LAPS_DATA_ROOT/static/widebandlist.txt not found..."
    exit
endif

#Obtain location of installed software that converts radar files
setenv JAVA `head -1 $LAPS_DATA_ROOT/static/widebandlist.txt`
if (! -e $JAVA) then
    echo "ERROR: WIDEBAND software JAVA $JAVA not found..."
    exit
endif

#Obtain NXJAR
setenv NXJAR `head -2 $LAPS_DATA_ROOT/static/widebandlist.txt | tail -1`
if (! -e $NXJAR) then
    echo "ERROR: NXJAR $NXJAR not found..."
    exit
endif

#Obtain NCJAR
setenv NCJAR `head -3 $LAPS_DATA_ROOT/static/widebandlist.txt | tail -1`
if (! -e $NCJAR) then
    echo "ERROR: NCJAR $NCJAR not found..."
    exit
endif

#Obtain DASET
setenv DASET `head -4 $LAPS_DATA_ROOT/static/widebandlist.txt | tail -1`
if (! -e $NCJAR) then
    echo "ERROR: DASET $DASET not found..."
    exit
endif

#Obtain location of Archive-II data files
setenv INPUTROOT `head -5 $LAPS_DATA_ROOT/static/widebandlist.txt | tail -1`
if (! -e $INPUTROOT) then
    echo "ERROR: WIDEBAND INPUTROOT data directory $INPUTROOT not found..."
    exit
endif

#Define Casedate (and hour within)
if ($MODETIME == "realtime") then
    setenv OUTPUTROOT  $LAPS_DATA_ROOT/lapsprd/rdr/wideband
    setenv HR1 `head -3 $LAPS_DATA_ROOT/time/systime.dat | tail -1`
    set HR2=$HR1
    @ HR2++
    if ($HR2 == 24) then
        set HR2=00
    endif
    echo "Processing hours $HR1 and $HR2"

#   Name of executable that converts from Nexrad format to NetCDF
    setenv NEXRAD_2_NETCDF ArchiveNexrad2NetCDF_latest_static

    setenv A9TIME `head -2 $LAPS_DATA_ROOT/time/systime.dat | tail -1`

else # archive case

#   Access additional command line args
    setenv YEAR $7
    setenv MONTH $8
    setenv DATE $9
    setenv HOUR $10
    setenv YYDDD $11

#   Use of Java for archive:
    # Due complication with laps_driver.pl, Yuanfu decided to use environment variables to pass in $JAVA, NXJAR and NCJAR instead of through arguments:
#   setenv DASET "ucar.nc2.dataset.NetcdfDataset" # Check with JAVA script option
    echo "YEAR=$YEAR"
    echo "DATE=$DATE"
    echo "HOUR=$HOUR"
    echo "YYDDD=$YYDDD"
    echo "MONTH=$MONTH"
    echo "JAVA=$JAVA"
    echo "NXJAR=$NXJAR"
    echo "NCJAR=$NCJAR"
    echo "DASET=$DASET"

    setenv A9TIME $YYDDD$HOUR\00

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
    setenv NEXRAD_2_NETCDF ArchiveNexrad2NetCDF

endif

#Name of NetCDF CDL file
setenv CDL wsr88d_super_res_wideband.cdl

echo " "

echo "INPUTROOT=$INPUTROOT"
#echo "INSTALLROOT=$INSTALLROOT"
echo "OUTPUTROOT=$OUTPUTROOT"

#setenv NEXRAD_SITES $INSTALLROOT/NexradSite.cfg
#if (! -e $NEXRAD_SITES) then
#    echo "ERROR: $NEXRAD_SITES not found..."
#    exit
#endif

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
  echo "New radar / cycle time..."
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
                                              -p $RADAR -o $OUTPUTROOT/$RADAR/netcdf \
                                              -c $INSTALLROOT/cdl/$CDL \
                                              -t $INPUTROOT/$RADAR/$file
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
      echo "listing of pre-existing output files for this hour "$YYDDD$HOUR
      echo "ls -1 $OUTPUTROOT/$RADAR/netcdf/$YYDDD$HOUR* | wc -l"

      setenv COUNT `ls -1 $OUTPUTROOT/$RADAR/netcdf/$YYDDD$HOUR* | wc -l`
      echo "# of pre-existing output files is $COUNT"

      if ($COUNT == "0" || $SKIP != "yes") then
#         This works with the archived radar data (e.g. for IHOP)
          echo "Generating radar $RADAR at $YEAR $MONTH $DATE $HOUR..."

#         Filename convention for NCDC
#         ls -l $INPUTROOT
	  echo "$INPUTROOT $RADAR $YEAR $MONTH $DATE $HOUR"
          ls -l $INPUTROOT/$RADAR$YEAR$MONTH$DATE\_$HOUR*
          setenv COUNT_NCDC `ls -1 $INPUTROOT/*$RADAR$YEAR$MONTH$DATE\_$HOUR* | wc -l`
          if ($COUNT_NCDC != "0") then
            echo " "
            echo "running command(s) for NCDC file format(s)..."

#           Process low res data (parse compressed hhmmss files by hour)
            setenv CDL wsr88d_low_res_wideband.cdl 

            #echo find /$INPUTROOT -name "*$RADAR$YEAR$MONTH$DATE\_$HOUR????.Z" -exec $INSTALLROOT/bin/$NEXRAD_2_NETCDF -l $OUTPUTROOT/$RADAR/log 
            #echo "                        -p $RADAR -o $OUTPUTROOT/$RADAR/netcdf -c $INSTALLROOT/cdl/$CDL -t {} \;"

            # Using ITS Java script to convert NCDC data to netCDF by Yuanfu Xie

            # Check availability of environment varibles needed for java converter:
            if (! $?JAVA || ! $?NXJAR || ! $?NCJAR) then
              echo "To use Java converter, one has to set environment variables:"
              echo "JAVA to where Java executable is"
              echo "NXJAR to where nexrad-nyquist jar file is"
              echo "NCJAR to where netcdfAll*.jar file is"
              echo "These are needed for the current Java converter; otherwise"
              echo "  check the java script for the options!"
              echo "Set them and try again! Bye"
              exit
            endif
            echo find /$INPUTROOT -name "*$RADAR$YEAR$MONTH$DATE\_$HOUR????.Z" -exec $JAVA -cp $NXJAR\:$NCJAR $DASET -in inputfile -out outputfile 

            # Comment out LAPS conversion: 2 lines
            # find /$INPUTROOT -name "*$RADAR$YEAR$MONTH$DATE\_$HOUR????.Z" -exec $INSTALLROOT/bin/$NEXRAD_2_NETCDF -l $OUTPUTROOT/$RADAR/log \
            #                        -p $RADAR -o $OUTPUTROOT/$RADAR/netcdf -c $INSTALLROOT/cdl/$CDL -t {} \;

            echo ""
            echo "***Start trying files in form *$RADAR$YEAR$MONTH$DATE\_$HOUR????.Z"
            pushd $INPUTROOT
            find . -name "*$RADAR$YEAR$MONTH$DATE\_$HOUR????.Z" \
              -exec $JAVA -cp $NXJAR\:$NCJAR $DASET \
              -in {} -out $OUTPUTROOT/$RADAR/netcdf/\{}.nc \;
            popd
            echo "***Exit first try"
            echo ""

#           Process low res data (parse uncompressed hhmmss files by hour)

            # Comment out LAPS conversion: 4 lines by Yuanfu Xie
            # echo find /$INPUTROOT -name "*$RADAR$YEAR$MONTH$DATE\_$HOUR????" -exec $INSTALLROOT/bin/$NEXRAD_2_NETCDF -l $OUTPUTROOT/$RADAR/log/$YEAR$MONTH$DATE\_$HOUR.log 
            #echo "                        -p $RADAR -o $OUTPUTROOT/$RADAR/netcdf -c $INSTALLROOT/cdl/$CDL -t {} \;"

            # find /$INPUTROOT -name "*$RADAR$YEAR$MONTH$DATE\_$HOUR????" -exec $INSTALLROOT/bin/$NEXRAD_2_NETCDF -l $OUTPUTROOT/$RADAR/log \
            #                        -p $RADAR -o $OUTPUTROOT/$RADAR/netcdf -c $INSTALLROOT/cdl/$CDL -t {} \;

            # Use ITS java by Yuanfu Xie
            echo ""
            echo "***Start try files in form *$RADAR$YEAR$MONTH$DATE\_$HOUR????"
            pushd $INPUTROOT
            find . -name "*$RADAR$YEAR$MONTH$DATE\_$HOUR????" \
              -exec $JAVA -cp $NXJAR\:$NCJAR $DASET \
              -in {} -out $OUTPUTROOT/$RADAR/netcdf/\{}.nc \;
            popd
            echo "Exit second try"
            echo ""

#           Process super res data (parse hhmmss files by hour)
            setenv CDL wsr88d_super_res_wideband.cdl

            # Comment out LAPS conversion: 4 lines by Yuanfu Xie
            # echo find /$INPUTROOT -name "*$RADAR$YEAR$MONTH$DATE\_$HOUR????\_V0?*" -exec $INSTALLROOT/bin/$NEXRAD_2_NETCDF -l $OUTPUTROOT/$RADAR/log 
            # echo "                        -p $RADAR -o $OUTPUTROOT/$RADAR/netcdf -c $INSTALLROOT/cdl/$CDL -t {} \;"
            #find /$INPUTROOT -name "*$RADAR$YEAR$MONTH$DATE\_$HOUR????\_V0?*" -exec $INSTALLROOT/bin/$NEXRAD_2_NETCDF -l $OUTPUTROOT/$RADAR/log \
            #                        -p $RADAR -o $OUTPUTROOT/$RADAR/netcdf -c $INSTALLROOT/cdl/$CDL -t {} \;

            # Use ITS Java conversion: by Yuanfu Xie
            echo ""
            echo "***Start try files in form *$RADAR$YEAR$MONTH$DATE\_$HOUR????\_V0?*"
            pushd $INPUTROOT
            find . -name "*$RADAR$YEAR$MONTH$DATE\_$HOUR????\_V0?*" \
              -exec $JAVA -cp $NXJAR\:$NCJAR $DASET \
              -in {} -out $OUTPUTROOT/$RADAR/netcdf/\{}.nc \;
            popd
            echo "***Exit third try"
            echo ""

            echo "Rename and gunzip the files in $OUTPUTROOT/$RADAR/netcdf"
            pushd $OUTPUTROOT/$RADAR/netcdf
                foreach OUTFILE (*.gz.nc)
                    setenv OUTFILE2 `echo $OUTFILE  | sed 's/\.[^.]*$//'`
                    setenv OUTFILE3 `echo $OUTFILE2 | sed 's/\.[^.]*$//'`
                    setenv OUTFILE4 $OUTFILE3.nc
                    echo "mv $OUTFILE $OUTFILE4"
                          mv $OUTFILE $OUTFILE4
                end                
            popd

            echo         "ls -1 $OUTPUTROOT/$RADAR/netcdf/$RADAR$YEAR$MONTH$DATE\_$HOUR* | wc -l"
            setenv COUNT `ls -1 $OUTPUTROOT/$RADAR/netcdf/$RADAR$YEAR$MONTH$DATE\_$HOUR* | wc -l`
            echo "Updated # of output files is $COUNT"
            echo "Finished radar $RADAR at $YEAR $MONTH $DATE $HOUR..."

# If using data from the /public archive, the previous find command should return a COUNT of 0.
# Then try to use the foreach command from the realtime portion of the loop above.
#         else
#           pushd $INPUTROOT/$RADAR
#           foreach file (*)
#             echo " "
#             echo "running conversion for Archive-II (public) file $INPUTROOT/$RADAR/$file"
#             $INSTALLROOT/bin/$NEXRAD_2_NETCDF -l $OUTPUTROOT/$RADAR/log -p $RADAR -o $OUTPUTROOT/$RADAR/netcdf -c $INSTALLROOT/cdl/$CDL -t $INPUTROOT/$RADAR/$file
#           end

          endif

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




