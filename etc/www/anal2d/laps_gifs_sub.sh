#!/bin/sh

#laps_gifs_sub.sh

umask 000

prod=$1
WINDOW=$2
LAPS_ETC=$3
WWW_DIR=$4
utc_hhmm=$5
LAPS_DATA_ROOT=$6
latest=$7
datetime=$8
RESOLUTION=$9

SCRATCH_DIR=$WWW_DIR/anal2d

uscore="_"
MACHINE=`uname -n`
PLATFORM=`uname -s`

#export NCARG_ROOT=/usr/local/apps/ncarg-4.0.1
#setenv NCARG_ROOT `cat /usr/nfs/lapb/bin/ncarg_root`
export NCARG_ROOT=`cat /usr/nfs/lapb/bin/ncarg_root`
SUPMAP_DATA_DIR=/home/elvis/mcdonald/data/supmap/

#alias ctrans '/usr/local/apps/ncarg-4.0.1/bin/ctrans  -verbose'

echo "prod ="$prod
echo "WINDOW ="$WINDOW
echo "NCARG_ROOT ="$NCARG_ROOT
echo "utc_hhmm ="$utc_hhmm
echo "LAPS_DATA_ROOT ="$LAPS_DATA_ROOT
echo "latest ="$latest
echo "RESOLUTION ="$RESOLUTION
echo "LAPS_ETC = "$LAPS_ETC
echo "WWW_DIR ="$WWW_DIR
echo "SCRATCH_DIR ="$SCRATCH_DIR

cd $WWW_DIR/anal2d

ulimit -t 300

ls -l $EXE_DIR/lapsplot.exe

echo "Running $EXE_DIR/lapsplot.exe"; date -u

head -2 $LAPS_DATA_ROOT/time/systime.dat | tail -1 | cut -c2-10     > $SCRATCH_DIR/lapsplot.$prod.tmp

if test -r $LAPS_DATA_ROOT/static/www/lapsplot.$prod; then
    echo "Input to lapsplot.exe (dataroot) = $LAPS_DATA_ROOT/static/www/lapsplot.$prod"
    cat     $LAPS_DATA_ROOT/static/www/lapsplot.$prod              >> $SCRATCH_DIR/lapsplot.$prod.tmp
elif test -r $SCRATCH_DIR/lapsplot.$prod; then
    echo "Input to lapsplot.exe (local) = $SCRATCH_DIR/lapsplot.$prod"
    cat     $SCRATCH_DIR/lapsplot.$prod                            >> $SCRATCH_DIR/lapsplot.$prod.tmp
else
    echo "Input to lapsplot.exe (non-local www) = $LAPS_ETC/lapsplot.$prod"
    cat     $LAPS_ETC/lapsplot.$prod                               >> $SCRATCH_DIR/lapsplot.$prod.tmp
fi

$EXE_DIR/lapsplot.exe                                           < $SCRATCH_DIR/lapsplot.$prod.tmp

ls -l gmeta

# Test whether we are at FSL on 'gizmo'
if test -r /data/lapb; then
    netpbm=no
else
    netpbm=yes
fi

# This if block is used to further assess whether we have access to netpbm
if test "$MACHINE" = "headnode.fsl.noaa.gov"; then
    ctransarg=avs
    ctransext=avs
    netpbm=no

elif test -r /usr/local/apps/ncarg-4.2.2-pgi/lib/libncarg.a; then
#   ctransarg=avs
#   ctransext=x
    ctransarg=sun
    ctransext=sun
    netpbm=yes

elif test "$PLATFORM" = "AIX"; then
    ctransarg=avs
    ctransext=x

elif test "$PLATFORM" = "HP-UX"; then
    ctransarg=sgi
    ctransext=sgi

else
    ctransarg=sgi
    ctransext=sgi
#   ctransarg=sun
#   ctransext=sun

fi

if test "$netpbm" = "yes"; then 
    date
    echo "Running $NCARG_ROOT/bin/ctrans | netpbm to make gmeta_$prod.gif file"
    $NCARG_ROOT/bin/ctrans -verbose -d sun -window $WINDOW -resolution $RESOLUTION gmeta | rasttopnm | ppmtogif > $SCRATCH_DIR/gmeta_$prod.gif
    date

else
    echo "Running $NCARG_ROOT/bin/ctrans to make gmeta_$prod.$ctransext file"

    $NCARG_ROOT/bin/ctrans -verbose -d $ctransarg -window $WINDOW -resolution $RESOLUTION gmeta > $SCRATCH_DIR/gmeta_$prod.$ctransext

    #For testing only...
    #cp gmeta gmeta_$prod

    ls -l $SCRATCH_DIR/gmeta_$prod.$ctransext
    date
    echo "Running gmeta.$ctransext to GIF converter laps_prod.sh on $MACHINE"
    $LAPS_ETC/laps_prod.sh $prod $SCRATCH_DIR $SCRATCH_DIR $ctransext

    #Cleanup
    'rm' -f $SCRATCH_DIR/gmeta_$prod.$ctransext
    echo " "
    date
    echo " "

fi

chmod 666 $SCRATCH_DIR/gmeta_$prod.gif

#Copy and link the output GIF file
mkdir -p                                                               $WWW_DIR/anal2d/archive/$prod
mkdir -p                                                               $WWW_DIR/anal2d/recent/$prod
#mkdir -p                                                              $WWW_DIR/anal2d/loop/$prod

$CP $SCRATCH_DIR/gmeta_$prod.gif                                    $WWW_DIR/anal2d/archive/$prod/$datetime.gif

rm -f                                                               $WWW_DIR/anal2d/recent/gmeta_$prod$uscore$latest.gif
ln -s $WWW_DIR/anal2d/archive/$prod/$datetime.gif                   $WWW_DIR/anal2d/recent/gmeta_$prod$uscore$latest.gif

rm -f                                                               $WWW_DIR/anal2d/recent/$prod/gmeta_$prod$uscore$latest.gif
ln -s $WWW_DIR/anal2d/archive/$prod/$datetime.gif                   $WWW_DIR/anal2d/recent/$prod/gmeta_$prod$uscore$latest.gif

rm -f                                                               $WWW_DIR/anal2d/recent/gmeta_$prod$uscore$utc_hhmm.gif
#ln -s $WWW_DIR/anal2d/archive/$prod/$datetime.gif                  $WWW_DIR/anal2d/recent/gmeta_$prod$uscore$utc_hhmm.gif

rm -f                                                               $WWW_DIR/anal2d/recent/$prod/gmeta_$prod$uscore$utc_hhmm.gif
ln -s $WWW_DIR/anal2d/archive/$prod/$datetime.gif                   $WWW_DIR/anal2d/recent/$prod/gmeta_$prod$uscore$utc_hhmm.gif

#ln -s $WWW_DIR/anal2d/archive/$prod/$datetime.gif                  $WWW_DIR/anal2d/loop/$prod/$datetime.gif

#Make the soft link for Web display of directories
cd    $WWW_DIR/anal2d/archive/$prod
rm -f                                                                  80.htaccess
ln -s $WWW_DIR/../../utilities/public-current-directory.htaccess       80.htaccess

#Create the files.txt script with the equivalent of files.sh
cd  $WWW_DIR/anal2d/loop/$prod
/bin/ls -1r *.gif > files.txt  

#Link to the files.cgi script
cd    $WWW_DIR/anal2d/loop/$prod
rm -f files.cgi
ln -s $WWW_DIR/../../looper/files.cgi files.cgi


