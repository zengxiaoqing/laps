#!/bin/csh

#laps_gifs.csh

#Set umask to 000 because oplapb is not in the same group as the rest of us
umask 000
#umask 002

limit coredumpsize 1k
limit cputime 30

setenv LAPS_DATA_ROOT $1
setenv WEB_ROOT $2
setenv WINDOW $3
setenv RESOLUTION $4
setenv DOMAIN_SUFFIX $5
setenv AGE $6
setenv EXE_DIR $7

echo "MACHINE ="`uname -a`
echo "LAPS_DATA_ROOT ="$LAPS_DATA_ROOT
echo "LAPS_ROOT/bin ="$EXE_DIR
echo "WEB_ROOT ="$WEB_ROOT
echo "WINDOW ="$WINDOW
echo "RESOLUTION ="$RESOLUTION

echo " "
limit
limit cputime unlimited
echo " "
limit
echo " "

if (-e /w3/lapb) then
#   In-house at FSL
    setenv SERVER_ROOT /w3/lapb
    setenv WWW_DOMAIN domains/$DOMAIN_SUFFIX
    setenv LAPS_GIFS     $WEB_ROOT/www/anal2d

else
#   External to FSL
    setenv SERVER_ROOT $LAPS_DATA_ROOT/lapsprd/www
    mkdir -p $SERVER_ROOT
    setenv WWW_DOMAIN domains/$DOMAIN_SUFFIX
    setenv LAPS_GIFS     $EXE_DIR/..

endif

setenv CP cp
setenv WWW_DIR $SERVER_ROOT/$WWW_DOMAIN

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

if (-e /home/lapb/albers) then
    echo $DOMAIN_SUFFIX > /home/lapb/albers/www/laps_anal/domain_$DOMAIN_SUFFIX
    chmod 666 /home/lapb/albers/www/laps_anal/domain_$DOMAIN_SUFFIX
endif

setenv SCRATCH_DIR   $LAPS_GIFS
#setenv EXE_DIR       /usr/nfs/lapb/parallel/laps/bin
setenv TIME_DIR      $LAPS_DATA_ROOT/time

date

echo " "
echo "Obtaining time info from systime.dat..."
ls -l $TIME_DIR/systime.dat

setenv utc_hour `head -3 $TIME_DIR/systime.dat | tail -1`
setenv utc_min  `head -4 $TIME_DIR/systime.dat | tail -1`
setenv datetime `head -2 $TIME_DIR/systime.dat | tail -1`
setenv latest latest

echo "utc_hour = "$utc_hour
echo "utc_min  = "$utc_min

cd $LAPS_GIFS
date

# PBE Graphic Product
echo "Generating Graphical PBE Product"; date -u
setenv prod pbe
$LAPS_GIFS/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# NBE Graphic Product
echo "Generating Graphical NBE Product"; date -u
setenv prod nbe
$LAPS_GIFS/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# TT Graphic Product
echo "Generating Graphical TT Product"; date -u
setenv prod tt
$LAPS_GIFS/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# TD Graphic Product
echo "Generating Graphical TD Product"; date -u
setenv prod td
$LAPS_GIFS/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# Wind Graphic Product
echo "Generating Graphical Wind Product"; date -u
setenv prod wd0
$LAPS_GIFS/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# Storm Total Graphic Product (Precip or Snow)
echo "Generating Graphical Storm Total Product (Precip or Snow)"; date -u
setenv prod sto
$LAPS_GIFS/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# Analysis Cycle Precip (Rain or Snow)
echo "Generating Graphical Analysis Cycle Precip Product (Rain or Snow)"; date -u
setenv prod p01
$LAPS_GIFS/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# Precip Type Graphic Product
echo "Generating Graphical Precip Type Product"; date -u
setenv prod ptt
$LAPS_GIFS/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# Snow Cover Graphic Product
echo "Generating Graphical Snow Cover Product"; date -u
setenv prod sc
$LAPS_GIFS/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# Total Precipitable Water Product
echo "Generating Graphical Total Precipitable Water Product"; date -u
setenv prod pw
$LAPS_GIFS/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# LI times Omega Product
echo "Generating Graphical LI times Omega Product"; date -u
setenv prod liw
$LAPS_GIFS/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# Low Level Radar Reflectivity Product
echo "Generating Graphical Low Level Radar Reflectivity Product"; date -u
setenv prod llr
$LAPS_GIFS/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# Sfc RH
echo "Generating SFC RH Product"; date -u
setenv prod rh
$LAPS_GIFS/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# Sfc FW
echo "Generating SFC FW Product"; date -u
setenv prod fw
$LAPS_GIFS/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# Sfc FW
echo "Generating SFC Fosberg Firewx Product"; date -u
setenv prod fwi
$LAPS_GIFS/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# PBL
echo "Generating PBL Product"; date -u
setenv prod pbl
$LAPS_GIFS/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# Ventilation Index & PBL Mean Wind
echo "Generating Ventilation Index & PBL Mean Wind"; date -u
setenv prod vnt
$LAPS_GIFS/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# Haines Mid
echo "Mid-Level Haines Index"; date -u
setenv prod ham
$LAPS_GIFS/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# Haines High
echo "High-Level Haines Index"; date -u
setenv prod hah
$LAPS_GIFS/laps_gifs_sub.sh $prod $WINDOW $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime $RESOLUTION

# X-sect
echo "Generating X-Sect Product"; date -u
setenv prod xct
#$LAPS_GIFS/laps_gifs_sub.sh $prod 0.06:0.14:0.94:0.86 $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime 732x614
#$LAPS_GIFS/laps_gifs_sub.sh $prod 0.06:0.12:0.94:0.84 $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime 732x614
#$LAPS_GIFS/laps_gifs_sub.sh $prod 0.06:0.13:0.94:0.85 $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime 732x614
#$LAPS_GIFS/laps_gifs_sub.sh $prod 0.06:0.12:0.94:0.85 $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime 732x614
$LAPS_GIFS/laps_gifs_sub.sh $prod 0.06:0.12:0.94:0.85 $LAPS_GIFS $WWW_DIR $utc_hour$utc_min $LAPS_DATA_ROOT $latest $datetime 792x664

# Update the date
echo "Update systime on web directory..."
head -5 $TIME_DIR/systime.dat | tail -1 > $WWW_DIR/anal2d/recent/systime_www_2d_anal.out
cat                                       $WWW_DIR/anal2d/recent/systime_www_2d_anal.out
echo `head -2 $TIME_DIR/c_time.dat | tail -1`" 7200" > $WWW_DIR/anal2d/recent/time.stamp.laps_2d_anal       

# Purge laps WWW files 
echo "Running purger of LAPS gif files, AGE = "$AGE
$SERVER_ROOT/purgers/laps_purge.com $WWW_DIR $AGE > $SERVER_ROOT/purgers/$DOMAIN_SUFFIX.out

# Set up soft link for on-the-fly program
cd $WWW_DIR
rm -f data
ln -s $LAPS_DATA_ROOT data
