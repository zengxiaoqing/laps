#!/bin/ksh

#laps_gifs_sub.sh
#limit cputime 600
umask 000

prod=$1
WINDOW=$2
LAPS_ETC=$3
WWW_DIR=$4
utc_hhmm=$5
LAPS_DATA_ROOT=$6
LAPSINSTALLROOT=$7
datetime=$8
RESOLUTION=$9


echo " "
echo "start laps_gifs_sub_fcst.sh..."

uscore="_"
MACHINE=`uname -n`
PLATFORM=`uname -s`

LOOPER=/usr/nfs/common/lapb/www/fcst2d/$domain/$model/upd_files_txt.pl

#export NCARG_ROOT=/usr/local/apps/ncarg-4.0.1
#setenv NCARG_ROOT `cat /usr/nfs/lapb/bin/ncarg_root`
if test -r /usr/nfs/lapb/bin/ncarg_root; then
    export NCARG_ROOT=`cat /usr/nfs/lapb/bin/ncarg_root`
else
    export NCARG_ROOT=`ncargpath root`
fi
SUPMAP_DATA_DIR=/home/elvis/mcdonald/data/supmap/

#alias ctrans '/usr/local/apps/ncarg-4.0.1/bin/ctrans  -verbose'

LEN_DATETIME=${#datetime}
echo "LEN_DATETIME ="$LEN_DATETIME

if test $LEN_DATETIME = "13" || test $LEN_DATETIME = "14"; then
    MODEL=$utc_hhmm
    echo "MODEL ="$MODEL
    SCRATCH_DIR=$WWW_DIR/fcst2d/$MODEL/.tmp
else
    SCRATCH_DIR=$WWW_DIR/anal2d
fi
rm -f $SCRATCH_DIR/*gmeta*
rm -f $SCRATCH_DIR/*.gif
ARCHIVE=$WWW_DIR/fcst2d/$MODEL/archive
RECENT=$WWW_DIR/fcst2d/$MODEL/recent
mkdir -p $SCRATCH_DIR
mkdir -p $ARCHIVE
mkdir -p $RECENT

echo "prod ="$prod
echo "WINDOW ="$WINDOW
echo "NCARG_ROOT ="$NCARG_ROOT
echo "utc_hhmm ="$utc_hhmm
echo "LAPS_DATA_ROOT ="$LAPS_DATA_ROOT
echo "LAPSINSTALLROOT ="$LAPSINSTALLROOT
echo "datetime ="$datetime
echo "RESOLUTION ="$RESOLUTION
echo "LAPS_ETC = "$LAPS_ETC
echo "WWW_DIR ="$WWW_DIR
echo "SCRATCH_DIR ="$SCRATCH_DIR
echo "MODEL="$MODEL

echo "ARCHIVE="$ARCHIVE

cd $SCRATCH_DIR

if test $LEN_DATETIME = "13" || test $LEN_DATETIME = "14"; then
    if test $LEN_DATETIME = "13"; then
        cycle=${datetime%????}
    else # length = 14
        cycle=${datetime%?????}
    fi

    fcsthr=${datetime#?????????}
    echo $MODEL    >  $SCRATCH_DIR/datetime 
    echo $datetime >> $SCRATCH_DIR/datetime 
    echo "cycle="$cycle
    echo "fcsthr="$fcsthr
    cat /dev/null > $SCRATCH_DIR/lapsplot.$prod.tmp
#   for fname in $LAPS_ETC/lapsplot.$prod.?
    for fname in $LAPS_DATA_ROOT/static/www/fcst2d/lapsplot.$prod.?
    do
      cat $fname >> $SCRATCH_DIR/lapsplot.$prod.tmp
      cat $SCRATCH_DIR/datetime >> $SCRATCH_DIR/lapsplot.$prod.tmp
    done
    LAPSPLOT_IN=$SCRATCH_DIR/lapsplot.$prod.tmp

else
    head -2 $LAPS_DATA_ROOT/time/systime.dat | tail -1 | cut -c2-10     > $SCRATCH_DIR/lapsplot.$prod.tmp

    if test -r $SCRATCH_DIR/lapsplot.$prod; then
        echo "Input to lapsplot.exe (local) = $SCRATCH_DIR/lapsplot.$prod"
        cat     $SCRATCH_DIR/lapsplot.$prod                            >> $SCRATCH_DIR/lapsplot.$prod.tmp
    else
        echo "Input to lapsplot.exe (non-local) = $LAPS_ETC/lapsplot.$prod"
        cat     $LAPS_ETC/lapsplot.$prod                               >> $SCRATCH_DIR/lapsplot.$prod.tmp
    fi

    LAPSPLOT_IN=$SCRATCH_DIR/lapsplot.$prod.tmp

fi

export LAPS_DATA_ROOT=$LAPS_DATA_ROOT

echo "Running $LAPSINSTALLROOT/bin/lapsplot.exe < $LAPSPLOT_IN"; date -u

$LAPSINSTALLROOT/bin/lapsplot.exe < $LAPSPLOT_IN

ls -l gmeta

FILESIZE=$(stat -c%s "gmeta")
echo "Size of gmeta = $FILESIZE bytes."

netpbm=no
if test "$MACHINE" = "headnode.fsl.noaa.gov"; then
    ctransarg=avs
    ctransext=avs

elif test -r /usr/local/apps/ncarg-4.2.2-pgi/lib/libncarg.a; then
#    ctransarg=avs
#    ctransext=x
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
    echo "Running $NCARG_ROOT/bin/ctrans -verbose -d sun -window $WINDOW -resolution $RESOLUTION gmeta | rasttopnm | ppmtogif > $SCRATCH_DIR/gmeta_$prod.gif"
    $NCARG_ROOT/bin/ctrans -verbose -d sun -window $WINDOW -resolution $RESOLUTION gmeta | rasttopnm | ppmtogif > $SCRATCH_DIR/gmeta_$prod.gif
    echo "Success in ctrans to make GIF"
    ls -l $SCRATCH_DIR/gmeta_$prod.gif
    date

else # laps_prod.sh option no longer used
  echo "Running $NCARG_ROOT/bin/ctrans to make gmeta_$prod.$ctransext file"

  $NCARG_ROOT/bin/ctrans -verbose -d $ctransarg -window $WINDOW -resolution $RESOLUTION gmeta > $SCRATCH_DIR/gmeta_$prod.$ctransext

  #For testing only...
  #cp gmeta gmeta_$prod

  ls -l $SCRATCH_DIR/gmeta_$prod.$ctransext
  date
  echo "Running gmeta.$ctransext to GIF converter laps_prod.sh on $MACHINE"
  $LAPS_ETC/laps_prod.sh $prod $SCRATCH_DIR $SCRATCH_DIR $ctransext
fi
chmod 666 $SCRATCH_DIR/gmeta_$prod.gif

#Copy and link the output GIF file
mkdir -p $ARCHIVE/$cycle/$prod
mkdir -p $RECENT/$prod
if [[ $fcsthr = "0000" ]]
then
  rm -f $RECENT/$prod/*.gif
  cp $SCRATCH_DIR/datetime $RECENT/cycle.datetime
fi

touch $WWW_DIR/fcst2d/$MODEL/latest.txt

cp $SCRATCH_DIR/gmeta_$prod.gif $ARCHIVE/$cycle/$prod/$fcsthr.gif
ln -fs $ARCHIVE/$cycle/$prod/$fcsthr.gif $RECENT/$prod/$fcsthr.gif
ls -l $ARCHIVE/$cycle/$prod/$fcsthr.gif
ls -l $RECENT/$prod/$fcsthr.gif

FILENAME=$ARCHIVE/$cycle/$prod/$fcsthr.gif
FILESIZE=$(stat -c%s "$FILENAME")
echo "Size of $FILENAME = $FILESIZE bytes."

if test $FILESIZE -lt "1000"; then
    echo "ERROR: Size of $FILENAME is too small"
fi

# Set up looper
# $LOOPER $ARCHIVE/$cycle/$prod

#Cleanup

'rm' -f $SCRATCH_DIR/gmeta_$prod.$ctransext
echo " "
echo "finish laps_gifs_sub_fcst.sh..."
date
exit
