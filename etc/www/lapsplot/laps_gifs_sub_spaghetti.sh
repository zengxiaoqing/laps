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
latest=$7
datetime=$8
RESOLUTION=$9
domain=$10
model=$11

echo " "
echo "start laps_gifs_sub_spaghetti.sh..."

uscore="_"
MACHINE=`uname -n`
PLATFORM=`uname -s`

LOOPER=/usr/nfs/common/lapb/www/fcst2d/$domain/$model/upd_files_txt.pl
#export NCARG_ROOT=/usr/local/apps/ncarg-4.0.1
#setenv NCARG_ROOT `cat /usr/nfs/lapb/bin/ncarg_root`
export NCARG_ROOT=`cat /usr/nfs/lapb/bin/ncarg_root`
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
echo "latest ="$latest
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
    cat $LAPS_DATA_ROOT/static/www/fcst2d/lapsplot.$prod.1 > $SCRATCH_DIR/lapsplot.$prod.tmp

#   Append each model & datetime with product segment 2
    cd $LAPS_DATA_ROOT/lapsprd/fua 

    for modeldir in *-*
    do
      cat $LAPS_DATA_ROOT/static/www/fcst2d/lapsplot.$prod.2 >> $SCRATCH_DIR/lapsplot.$prod.tmp
#     cat $modeldir >> $SCRATCH_DIR/lapsplot.$prod.tmp
      echo "modeldir = $modeldir"
      echo $modeldir >> $SCRATCH_DIR/lapsplot.$prod.tmp
      tail -1 $SCRATCH_DIR/datetime >> $SCRATCH_DIR/lapsplot.$prod.tmp
      cat $LAPS_DATA_ROOT/static/www/fcst2d/lapsplot.$prod.3 >> $SCRATCH_DIR/lapsplot.$prod.tmp
    done

    cd $SCRATCH_DIR

    cat $LAPS_DATA_ROOT/static/www/fcst2d/lapsplot.$prod.4 >> $SCRATCH_DIR/lapsplot.$prod.tmp

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

echo "Running $EXE_DIR/lapsplot.exe < $LAPSPLOT_IN"; date -u

$EXE_DIR/lapsplot.exe < $LAPSPLOT_IN

ls -l gmeta
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
    if test -r /whome; then
        CTRANS=/opt/ncl/5.1.0_bin/bin/ctrans
    elif test -r /usr/local/apps/ncarg-4.3.1.LINUX9; then
        CTRANS=/usr/local/apps/ncarg-4.3.1.LINUX9/bin/ctrans
    else
        CTRANS=$NCARG_ROOT/bin/ctrans
    fi

    date
    echo "Running $CTRANS and netpbm programs to make gmeta_$prod.gif file"
    which rasttopnm
    which ppmtogif
    echo "$CTRANS -verbose -d sun -window $WINDOW -resolution $RESOLUTION gmeta | rasttopnm | ppmtogif > $SCRATCH_DIR/gmeta_$prod.gif"
          $CTRANS -verbose -d sun -window $WINDOW -resolution $RESOLUTION gmeta | rasttopnm | ppmtogif > $SCRATCH_DIR/gmeta_$prod.gif
    ls -l $SCRATCH_DIR/gmeta_$prod.gif
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


cp $SCRATCH_DIR/gmeta_$prod.gif $ARCHIVE/$cycle/$prod/$fcsthr.gif
ln -fs $ARCHIVE/$cycle/$prod/$fcsthr.gif $RECENT/$prod/$fcsthr.gif

# Set up looper
$LOOPER $ARCHIVE/$cycle/$prod

#Cleanup

'rm' -f $SCRATCH_DIR/gmeta_$prod.$ctransext
echo " "
date
echo " "
exit
