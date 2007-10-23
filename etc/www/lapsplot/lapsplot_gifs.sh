#!/bin/sh 

#umask 000

proc=$1
WINDOW=$2
LAPS_GIFS=$3
LAPSPLOT_IN=$4
EXE_DIR=$5
export LAPS_DATA_ROOT=$6
export NCARG_ROOT=$7
datetime=$8
RESOLUTION=$9

SCRATCH_DIR=$LAPS_GIFS/scratch

uscore="_"
MACHINE=`uname -s`
NODE=`uname -n`

#export NCARG_ROOT=/usr/local/apps/ncarg-4.2.2-pgi
SUPMAP_DATA_DIR=/home/elvis/mcdonald/data/supmap/
#alias ctrans '/usr/local/apps/ncarg-4.0.1/bin/ctrans  -verbose'

echo "proc ="$proc
echo "WINDOW ="$WINDOW
echo "LAPS_GIFS = "$LAPS_GIFS
echo "SCRATCH_DIR =  $SCRATCH_DIR"
echo "EXE_DIR =  $EXE_DIR"
echo "NCARG_ROOT ="$NCARG_ROOT
echo "LAPS_DATA_ROOT ="$LAPS_DATA_ROOT
echo "latest ="$latest
echo "RESOLUTION ="$RESOLUTION
echo "LAPSPLOT_IN ="$LAPSPLOT_IN

limit cputime 1000

mkdir -p /scratch/lapb/www
mkdir -p $SCRATCH_DIR/$proc
cd $SCRATCH_DIR/$proc

#EXE_DIR=/usr/nfs/lapb/parallel/laps/bin

date -u
echo "Running $EXE_DIR/lapsplot.exe < $LAPSPLOT_IN on $MACHINE $NODE"
#$EXE_DIR/lapsplot.exe                                          < $LAPS_GIFS/lapsplot.in
$EXE_DIR/lapsplot.exe                                           < $LAPSPLOT_IN

pwd
ls -l $SCRATCH_DIR/$proc/gmeta

if test "$MACHINE" = "AIX"; then

#   Combination for IBM
    ext1=avs
    ext2=x
    ext3=gif
    netpbm=no

#   CTRANS=/usr/local/apps/ncarg-4.0.1/bin/ctrans

else

    netpbm=yes

#   Best combination for LINUX
    ext1=avs
    ext2=avs
    ext3=gif

#   The ones below will run but produce fewer colors in color images for LINUX

#   ext1=sun
#   ext2=sun

#   ext1=xwd
#   ext2=xwd

#   Note that sgi will not work in LINUX since we are using gmeta files with WINDOW/RESOLUTION set
#   ext1=sgi
#   ext2=sgi

#   ext1=sun
#   ext2=gif

#   CTRANS=/usr/local/apps/ncarg-4.2.2-pgi/bin/ctrans

fi

CTRANS=$NCARG_ROOT/bin/ctrans


#/usr/local/apps/ncarg-4.0.1/bin/ctrans -verbose -d avs -window $WINDOW -resolution $RESOLUTION gmeta > $SCRATCH_DIR/gmeta_$proc.x
#ctrans -d avs -window 0.0:0.08:1.0:0.92 -resolution 610x512 gmeta > $SCRATCH_DIR/gmeta_$proc.x
#/usr/local/apps/ncarg-4.0.1/bin/ctrans -verbose -d $ext1 -window $WINDOW -resolution $RESOLUTION gmeta > $SCRATCH_DIR/gmeta_temp_$proc.$ext2
#$NCARG_ROOT/bin/ctrans -verbose -d $ext1 -window $WINDOW -resolution $RESOLUTION $SCRATCH_DIR/gmeta > $SCRATCH_DIR/gmeta_temp_$proc.$ext2

#$CTRANS -verbose -d $ext1 -window $WINDOW -resolution $RESOLUTION $SCRATCH_DIR/$proc/gmeta > $SCRATCH_DIR/gmeta_temp_$proc.$ext2

#ls -l $SCRATCH_DIR/gmeta_temp_$proc.$ext2

date -u

#We assume we are running this script in LINUX and convert will not properly do AVS X on LINUX
if test "$netpbm" = "yes"; then 
    date
#   echo "Running $NCARG_ROOT/bin/ctrans | netpbm to make gmeta_$proc.gif file"
    echo "Running $NCARG_ROOT/bin/ctrans -verbose -d sun -window $WINDOW -resolution $RESOLUTION gmeta | rasttopnm | ppmtogif > $SCRATCH_DIR/gmeta_$proc.gif"
    $NCARG_ROOT/bin/ctrans -verbose -d sun -window $WINDOW -resolution $RESOLUTION gmeta | rasttopnm | ppmtogif > $SCRATCH_DIR/gmeta_$proc.gif

    date -u

#   Cleanup
    echo "Cleanup"
    mv gmeta $SCRATCH_DIR/gmeta_$proc.gm;  cd ..; rmdir $SCRATCH_DIR/$proc &

elif test "$ext2" = "x"; then
    date -u
    echo "Converting $SCRATCH_DIR/$proc/gmeta file with $CTRANS to $SCRATCH_DIR/gmeta_temp_$proc.$ext2"
    $CTRANS -verbose -d $ext1 -window $WINDOW -resolution $RESOLUTION $SCRATCH_DIR/$proc/gmeta > $SCRATCH_DIR/gmeta_temp_$proc.$ext2

    ls -l $SCRATCH_DIR/gmeta_temp_$proc.$ext2

    date -u

#   echo "Running imconv.ibm.exe on ren via rsh"
#   rsh ren /usr/nfs/avs_fsl/bin/imconv.ibm.exe $SCRATCH_DIR/gmeta_temp_$proc.$ext2 $SCRATCH_DIR/gmeta_$proc.gif

#   Note that 'convert' is available on Linux but it appears to be slow for GIFs
#   echo "Running convert $SCRATCH_DIR/gmeta_temp_$proc.$ext2 $SCRATCH_DIR/gmeta_$proc.$ext3 on $MACHINE"
#   convert $SCRATCH_DIR/gmeta_temp_$proc.$ext2 $SCRATCH_DIR/gmeta_$proc.$ext3

    echo "Running imconv.ibm.exe on $MACHINE $NODE"
    /usr/nfs/avs_fsl/bin/imconv.ibm.exe $SCRATCH_DIR/gmeta_temp_$proc.$ext2 $SCRATCH_DIR/gmeta_$proc.gif

    date -u

#   Cleanup
    echo "Cleanup"
    rm -f $SCRATCH_DIR/gmeta_temp_$proc.$ext2; mv gmeta $SCRATCH_DIR/gmeta_$proc.gm;  cd ..; rmdir $SCRATCH_DIR/$proc &

fi

chmod 666 $SCRATCH_DIR/gmeta_$proc.$ext3


echo " "
date -u
echo " "

