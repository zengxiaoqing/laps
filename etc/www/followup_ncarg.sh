#!/bin/sh

#exit

#Input LAPS_DATA_ROOT
LAPS_DATA_ROOT=$1

#Input domain name (e.g. fsld, t2)
#WEB_DATA=/w3/lapb/domains/$2
WEB_DATA=$LAPS_DATA_ROOT/lapsprd/www

#Input LAPS_ROOT ($LAPSINSTALLROOT)
LAPS_ROOT=$3

#Input LAPS A9Time (yyydddhhmm)
LAPS_A9TIME=$4

export LAPSINSTALLROOT=$LAPS_ROOT
WEB_NFS=$LAPS_ROOT/etc

SCHED="$LAPS_DATA_ROOT/time"
cd $SCHED

#utc_hour=`head -3 systime.dat | tail -1`
#utc_min=`head -4 systime.dat | tail -1`

utc_hour=`echo $LAPS_A9TIME | cut -c6-7`
utc_min=`echo $LAPS_A9TIME | cut -c8-9`

log_min=20

umask 000

echo " "                                  
echo "Running followup_ncarg.sh script" 
echo "user = "`whoami` 
echo "machine = "`uname -n` 
echo "domain name ="$2                                  
echo " "                                  

#$LAPS_ROOT/etc/wait_for_file.sh $LAPS_DATA_ROOT/log/sched.lock 50

#LAPS WWW 2-D NCAR graphics products
echo " "
date


#Determine shape (aspect ratio) of image based on domain name
#In the future this could perhaps read the nx/ny to determine this
#The threshold ratio is > 1.192 for using a "wide" window and > 1.40 for
#a "wide2" window
if test "$2" = "dwfe1-rcsv"
then
    WINDOW=0.0:0.08:1.0:0.92
    RESOLUTION=1056x885

elif test "$2" = "stmas-backup"
then
    WINDOW=0.0:0.08:1.0:0.92
    RESOLUTION=1056x885

elif test "$2" = "dwfe1-rcsv-nemo"
then
    WINDOW=0.0:0.08:1.0:0.92
    RESOLUTION=1056x885

elif test "$2" = "glaps"
then
    WINDOW=0.0:0.14:1.0:0.86
    RESOLUTION=1200x885

elif test "$2" = "glaps_sos"
then
    WINDOW=0.0:0.0796:1.0:0.9204
    RESOLUTION=3193x2681

elif test "$2" = "glaps_sos2"
then
    WINDOW=0.0:0.0796:1.0:0.9204
    RESOLUTION=3193x2681

elif test "$5" = "wide"
then
    WINDOW=0.0:0.08:1.0:0.92
    RESOLUTION=1056x885

elif test "$5" = "wide2"
then
    WINDOW=0.0:0.14:1.0:0.86
    RESOLUTION=1200x885

elif test "$2" = "f3"
then
    WINDOW=0.0:0.08:1.0:0.92
    RESOLUTION=1056x885

else
    WINDOW=0.0:0.0:1.0:1.0
    RESOLUTION=885x885

fi

echo "WWW 2-D laps_gifs started" 
echo "WEB_NFS="$WEB_NFS
echo "LAPS_DATA_ROOT="$LAPS_DATA_ROOT
echo "WEB_DATA="$WEB_DATA
echo "LAPSROOT="$LAPS_ROOT
echo "WINDOW="$WINDOW
echo "RESOLUTION="$RESOLUTION
echo "NCARG_ROOT="$NCARG_ROOT
mkdir -p $WEB_DATA 
mkdir -p $WEB_DATA/anal2d
AGE=+7 

LOG=$WEB_DATA/anal2d/laps_gifs.log.$utc_hour$utc_min

echo "Running command $WEB_NFS/www/anal2d/laps_gifs.csh $LAPS_DATA_ROOT $WEB_NFS $WINDOW $RESOLUTION $2 $AGE $LAPS_ROOT/bin $LAPS_A9TIME"
                      $WEB_NFS/www/anal2d/laps_gifs.csh $LAPS_DATA_ROOT $WEB_NFS $WINDOW $RESOLUTION $2 $AGE $LAPS_ROOT/bin $LAPS_A9TIME       1> $LOG 2>&1

echo "Additional log info in "$LOG

rm -f /tmp/GSEG0*
rm -f /tmp/GSEG1*
rm -f /tmp/GSEG2*
rm -f /tmp/GSEG3*
rm -f /tmp/GSEG4*
rm -f /tmp/GSEG5*
rm -f /tmp/GSEG6*
rm -f /tmp/GSEG7*
rm -f /tmp/GSEG8*
rm -f /tmp/GSEG9*
rm -f /tmp/GSEG*

echo " "
echo "followup_ncarg.sh finish"  
