#!/bin/sh

#exit

#Input LAPS_DATA_ROOT
LAPS_DATA_ROOT=$1

#Input domain name (e.g. fsld, t2)
#WEB_DATA=/w3/lapb/domains/$2
WEB_DATA=$LAPS_DATA_ROOT/www

#Input LAPS_ROOT ($LAPSINSTALLROOT)
LAPS_ROOT=$2

WEB_NFS=$LAPS_ROOT/etc

SCHED="$LAPS_DATA_ROOT/time"
cd $SCHED

utc_hour=`head -3 systime.dat | tail -1`
utc_min=`head -4 systime.dat | tail -1`

log_min=20

umask 000

echo " "                                  
echo "Running followup_ncarg.sh script" 
echo "user = "`whoami` 
echo "machine = "`uname -n`                                   
echo " "                                  

$LAPS_ROOT/etc/wait_for_file.sh $LAPS_DATA_ROOT/log/sched.lock 50

#LAPS WWW 2-D NCAR graphics products
echo " "
date


#Determine shape (aspect ratio) of image based on domain name
#In the future this could perhaps read the nx/ny to determine this
#The threshold ratio is > 1.192 for using a rectangular window
if test "$2" = "wiap"  
then
    WINDOW=0.0:0.08:1.0:0.92
    RESOLUTION=792x664

elif test "$2" = "nos_nb"
then
    WINDOW=0.0:0.08:1.0:0.92
    RESOLUTION=792x664

elif test "$2" = "f3"
then
    WINDOW=0.0:0.08:1.0:0.92
    RESOLUTION=792x664

else
    WINDOW=0.0:0.0:1.0:1.0
    RESOLUTION=664x664

fi

echo "WWW 2-D laps_gifs started" 
echo $WEB_NFS
echo $LAPS_DATA_ROOT
echo $WEB_DATA
echo $LAPS_ROOT
echo $WINDOW
echo $RESOLUTION
mkdir -p $WEB_DATA 
mkdir -p $WEB_DATA/anal2d
AGE=+7 
$WEB_NFS/www/anal2d/laps_gifs.csh $LAPS_DATA_ROOT $WEB_NFS $WINDOW $RESOLUTION $2 $AGE $LAPS_ROOT/bin \
       1> $WEB_DATA/anal2d/laps_gifs.log 2> $WEB_DATA/anal2d/laps_gifs.err

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
