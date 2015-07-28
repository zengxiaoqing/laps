#!/bin/bash -u

#####################################################################
#
# Wrapper script for LAPS lightning data assimilation method.  This
# script should be executed from LAPS sched.pl after mosaic_radar.exe
# and before cloud.exe
#
# In the LAPS repository this can be kept in the 'etc' directory and
# copied to the bin directory during the configure or make step
#
# Antti.Pessi@Vaisala.com 2014 & Steve Albers 2015
#####################################################################

# Set the LAPS domain. Domain is the last directory of LAPS_DATA_ROOT
# For example: LAPS_DATA_ROOT=/data/laps/domains/conus , where domain=conus

export TZ=UTC+0

export domain=`echo ${LAPS_DATA_ROOT##*/}`

echo Running runLDA.sh for domain $domain...
echo Current time:`date -u`
echo "LAPS_DATA_ROOT="$LAPS_DATA_ROOT
echo "LAPSINSTALLROOT="$LAPSINSTALLROOT

bindir=$LAPSINSTALLROOT/bin
rundir=$LAPS_DATA_ROOT/lightning/run
lapsprd=$LAPS_DATA_ROOT/lapsprd
namelists=$LAPS_DATA_ROOT/static

# Set the current cycle time from LAPS c_time.dat

epochtime=`cat $LAPS_DATA_ROOT/time/c_time.dat |tail -n1`
systime=`date -u -d "@${epochtime}" +'%Y-%m-%d %H:%M:00'`

echo "system time="$systime

# Find the lightning data path from namelist file lightning.nl
if [ -e ${namelists}/lightning.nl ]; then
    echo "${namelists}/lightning.nl was found"
else
    echo "ERROR: ${namelists}/lightning.nl does not exist - exiting"
    exit
fi

# Directory test
if [ $? -eq 0 ]; then
    echo 'Directory changed to '$LAPS_DATA_ROOT/lightning/run
else
    echo 'ERROR: failed to change directory to run - exiting'
    exit 1
fi

ligdir=`grep path_to_lightning_data ${namelists}/lightning.nl | cut -f2 -d"'"`

cd ${rundir}

# Set the current cycle

export cycle=`date -u -d "${systime} UTC" +'%y%j%H%M'`
rm -f $cycle.lig
rm -f $cycle.tmp.lig

export yyyy=`date -u -d "${systime} UTC" +'%Y'`
export   mm=`date -u -d "${systime} UTC" +'%m'`
export   dd=`date -u -d "${systime} UTC" +'%d'`
export   hh=`date -u -d "${systime} UTC" +'%H'`
export   mi=`date -u -d "${systime} UTC" +'%M'`

# Create lightning data file. This part is very specific to the
# lightning data format that the user has. The resulting file should
# contain the lightning observations for the entire assimilation
# period (e.g. 20 minutes) and the file format should be:
# date time latitude longitude
# For example: 2013-08-07 18:11:01  42.235 -124.848
# Current method for NOAA assumes 5-minute ascii files where the data timestamp
# is lagging behind the file timestamp
# data format in NOAA: yyjjjhhmm0005r.txt (e.g. 1417719450005r.txt)
# Files needed have timestamps e.g. 1855, 1900, 1905, 1910, 1915

echo "Create lightning data file"

fileprefix=`date -u -d "${systime} UTC -1 hours" +'%y%j%H'`
cat ${ligdir}/${fileprefix}550005r.txt >> ${cycle}.tmp.lig

fileprefix=`date -u -d "${systime} UTC" +'%y%j%H'`
for min in 00 05 10 15 
do
    cat $ligdir/${fileprefix}${min}0005r.txt >> $cycle.tmp.lig
done

# Choose the right lines from the .tmp file

hour=`date -u -d "${systime} UTC -1 hours" +'%H'`
grep ${hour}:5[0-9] $cycle.tmp.lig >> $cycle.lig 

hour=`date -u -d "${systime} UTC" +'%H'`
grep ${hour}:0[0-9] $cycle.tmp.lig >> $cycle.lig 


echo "Generate empty netcdf file from the template"
ncgen -o $cycle.vrz $LAPS_DATA_ROOT/cdl/vrz.cdl

echo "execute: lig2ref $yyyy-$mm-$dd $hh:$mi:00 $LAPS_DATA_ROOT"

${bindir}/lig2ref $yyyy-$mm-$dd $hh:$mi:00 $LAPS_DATA_ROOT

if [ -s $lapsprd/vrz/$cycle.vrz ] ; then
    cp -f $lapsprd/vrz/$cycle.vrz $lapsprd/vrz/$cycle.vrz.orig
fi

mv -f $cycle.vrz $lapsprd/vrz/$cycle.vrz

rm -f $cycle.lig

echo `date -u`
echo "runLDA.sh finished"

exit 0
