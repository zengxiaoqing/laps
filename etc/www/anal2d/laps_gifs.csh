#!/bin/csh

echo "Starting laps_gifs.csh..."

#Set umask to 000 because oplapb is not in the same group as the rest of us
umask 000
#umask 002

#limit coredumpsize 1k

setenv LAPS_DATA_ROOT $1
setenv WEB_ROOT $2
setenv WINDOW $3
setenv RESOLUTION $4
setenv DOMAIN_SUFFIX $5
setenv AGE $6
setenv EXE_DIR $7
setenv LAPS_A9TIME $8

echo "MACHINE ="`uname -a`
echo "LAPS_DATA_ROOT ="$LAPS_DATA_ROOT
echo "LAPS_ROOT/bin ="$EXE_DIR
echo "WEB_ROOT ="$WEB_ROOT
echo "WINDOW ="$WINDOW
echo "RESOLUTION ="$RESOLUTION
echo "DOMAIN_SUFFIX ="$DOMAIN_SUFFIX
#echo "NCARG_ROOT ="$NCARG_ROOT

#export $NCARG_ROOT

setenv LAPSINSTALLROOT $7/..

if (-e /w3/lapb) then
#   In-house at NOAA/ESRL/GSD/FAB
    setenv SERVER_ROOT /w3/lapb
    setenv WWW_DOMAIN domains/$DOMAIN_SUFFIX
    setenv LAPS_GIFS     $WEB_ROOT/www/anal2d

#   Set up soft link for on-the-fly program
    echo "updating on-the-fly page softlink..."
    mkdir -p /w3/lapb/domains/$DOMAIN_SUFFIX
    cd /w3/lapb/domains/$DOMAIN_SUFFIX

    if (! -e private_data) then
        ln -s $LAPS_DATA_ROOT private_data
    endif

    ls -l /w3/lapb/domains/$DOMAIN_SUFFIX
    echo " "

#   Check whether lapsplot.exe exists for Linux runs
    if (! -e $EXE_DIR/lapsplot.exe) then
        echo "pointing to 64-bit binary for lapsplot.exe"
        setenv EXE_DIR /usr/nfs/lapb/builds64/laps/bin
    endif

else
#   External to NOAA/ESRL/GSD/FAB
    setenv SERVER_ROOT $LAPS_DATA_ROOT                
    mkdir -p $SERVER_ROOT
    setenv WWW_DOMAIN lapsprd/www              
    setenv LAPS_GIFS     $LAPSINSTALLROOT/etc/www/anal2d

endif

setenv SCRIPTDIR $LAPSINSTALLROOT/etc/www/anal2d

if (-d $LAPS_DATA_ROOT) then
    echo "$LAPS_DATA_ROOT exists"
else
    echo "$LAPS_DATA_ROOT does not exist"
    exit
endif

setenv CP cp
setenv WWW_DIR $SERVER_ROOT/$WWW_DOMAIN

echo "WWW_DIR ="$WWW_DIR

mkdir -p $WWW_DIR/anal2d
mkdir -p $WWW_DIR/anal2d/archive
mkdir -p $WWW_DIR/anal2d/recent
#mkdir -p $WWW_DIR/anal2d/loop

#Link loop directory to archive
cd    $WWW_DIR/anal2d
rm -f loop
ln -s archive loop

#Link in java looper directory
if (-e $SERVER_ROOT/looper) then
    cd    $SERVER_ROOT/looper
    rm -f pub_$DOMAIN_SUFFIX
    ln -s $WWW_DIR/anal2d/loop pub_$DOMAIN_SUFFIX
endif

#Set up scripts for pregenerated pages
setenv OPLAPB /home/fab/oplapb
if (-e $OPLAPB/www/laps_anal) then
    echo "updating pregenerated page directory..."
    echo $DOMAIN_SUFFIX > ~oplapb/www/laps_anal/domain_$DOMAIN_SUFFIX
    chmod 666 $OPLAPB/www/laps_anal/domain_$DOMAIN_SUFFIX
    ls -l $OPLAPB/www/laps_anal/domain_$DOMAIN_SUFFIX
    echo " "
endif

setenv TIME_DIR      $LAPS_DATA_ROOT/time

date -u

echo " "
echo "Obtaining time info from LAPS_A9TIME..."
#ls -l $TIME_DIR/systime.dat

#setenv utc_hour `head -3 $TIME_DIR/systime.dat | tail -1`
#setenv utc_min  `head -4 $TIME_DIR/systime.dat | tail -1`
#setenv datetime `head -2 $TIME_DIR/systime.dat | tail -1`

setenv utc_hour `echo $LAPS_A9TIME | cut -c6-7`
setenv utc_min  `echo $LAPS_A9TIME | cut -c8-9`
setenv datetime `echo $LAPS_A9TIME | cut -c1-9`
setenv latest latest

echo "utc_hour = "$utc_hour
echo "utc_min  = "$utc_min

echo " "
limit
limit cputime unlimited
echo " "
limit
echo " "

cd $LAPS_GIFS
date -u

# Wind Graphic Product
echo "Generating Graphical Wind Product"; date -u
setenv prod wd0
$SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# TT Graphic Product
echo "Generating Graphical TT Product"; date -u
setenv prod tt
$SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# TD Graphic Product
echo "Generating Graphical TD Product"; date -u
setenv prod td
$SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# Sfc RH
echo "Generating SFC RH Product"; date -u
setenv prod rh
$SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# DPD
# echo "Generating Dewpoint Depression Product"; date -u
# setenv prod pbl
# $SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# Sfc Theta(e)
echo "Generating SFC Theta(e) Product"; date -u
setenv prod fwi
$SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# Check for existence of 3D products
if (-e $LAPS_DATA_ROOT/lapsprd/lt1/$datetime.lt1 || \
    -e $LAPS_DATA_ROOT/lapsprd/balance/lt1/$datetime.lt1) then

#   PBE Graphic Product
    echo "Generating Graphical PBE Product"; date -u
    setenv prod pbe
    $SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

#   NBE Graphic Product
    echo "Generating Graphical NBE Product"; date -u
    setenv prod nbe
    $SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

#   H5 Graphic Product
    echo "Generating Graphical Wind Product"; date -u
    setenv prod h5b
    $SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

#   Storm Total Graphic Product (Precip or Snow)
    echo "Generating Graphical Storm Total Product (Precip or Snow)"; date -u
    setenv prod sto
    $SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

#   Analysis Cycle Precip (Rain or Snow)
    echo "Generating Graphical Analysis Cycle Precip Product (Rain or Snow)"; date -u
    setenv prod p01
    $SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

#   Precip Type Graphic Product
    echo "Generating Graphical Precip Type Product"; date -u
    setenv prod ptt
    $SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

#   Composite Reflectivity Graphic Product
    echo "Generating Graphical Composite Reflectivity Product"; date -u
    setenv prod crf
    $SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

#   Snow Cover Graphic Product
    echo "Generating Graphical Snow Cover Product"; date -u
    setenv prod sc
    $SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

#   Total Precipitable Water Product
    echo "Generating Graphical Total Precipitable Water Product"; date -u
    setenv prod pw
    $SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

#   LI times Omega Product
    echo "Generating Graphical LI times Omega Product"; date -u
    setenv prod liw
    $SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

#   Low Level Radar Reflectivity Product
    echo "Generating Graphical Low Level Radar Reflectivity Product"; date -u
    setenv prod llr
    $SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

#   Upslope Moisture Flux
#   echo "Generating Ventilation Index & PBL Mean Wind"; date -u
#   setenv prod umf
#   $SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

#   Solar Radiation    
    echo "Generating Solar Radiation"; date -u
    setenv prod sol
    $SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

#   11 micron brightness temperature
    echo "Generating 11 micron brightness temperature"; date -u
    setenv prod s8a
    $SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION
 
#   Ground Temperature
#   echo "Generating Ground Temperature"; date -u
#   setenv prod tgd
#   $SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION
 
#   Sfc Fireweather
    echo "Generating SFC Fireweather Product"; date -u
    setenv prod fw
    $SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

#   PBL Depth & PBL Mean Wind
    echo "Generating PBL Depth & PBL Mean Wind"; date -u
    setenv prod pbl
    $SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

#   Ventilation Index & PBL Mean Wind
    echo "Generating Ventilation Index & PBL Mean Wind"; date -u
    setenv prod vnt
    $SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

#   Haines Mid
#   echo "Mid-Level Haines Index"; date -u
#   setenv prod ham
#   $SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

#   Haines High
#   echo "High-Level Haines Index"; date -u
#   setenv prod hah
#   $SCRIPTDIR/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

#   X-sect
    echo "Generating X-Sect Product"; date -u
    setenv prod xct
#   $SCRIPTDIR/laps_gifs_sub.sh $prod 0.06:0.12:0.94:0.85 $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime 732x614
#   $SCRIPTDIR/laps_gifs_sub.sh $prod 0.06:0.12:0.94:0.85 $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime 792x664
    $SCRIPTDIR/laps_gifs_sub.sh $prod 0.06:0.12:0.94:0.85 $LAPS_GIFS $WWW_DIR $LAPS_A9TIME $LAPS_DATA_ROOT $latest $datetime 1056x885

else
    echo "$LAPS_DATA_ROOT/lapsprd/lt1/$datetime.lt1 file not present"
    echo "$LAPS_DATA_ROOT/lapsprd/balance/lt1/$datetime.lt1 file not present"
    echo "skipping non-surface based products..."

endif

# Update the date
echo "Update systime on web directory..."
head -5 $TIME_DIR/systime.dat | tail -1 > $WWW_DIR/anal2d/recent/systime_www_2d_anal.out
cat                                       $WWW_DIR/anal2d/recent/systime_www_2d_anal.out
echo `head -2 $TIME_DIR/c_time.dat | tail -1`" 7200" > $WWW_DIR/anal2d/recent/time.stamp.laps_2d_anal       

# Purge laps WWW files 
echo "Running purger of LAPS gif files, AGE = "$AGE
$SERVER_ROOT/purgers/laps_purge.com $WWW_DIR $AGE > $SERVER_ROOT/purgers/$DOMAIN_SUFFIX.out

echo "Finishing laps_gifs.csh..."
