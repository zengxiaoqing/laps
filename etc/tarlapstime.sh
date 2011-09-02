#!/bin/sh

# For an hourly cycle a good time to run this would be 110 minutes later (at 50 after the hour)
# with a MINAGE of 48 minutes and a MAXAGE of 112 minutes

# First command line argument is the LAPS_DATA_ROOT

# Second command line argument is the laps time in yydddhhmm format

# Third command line argument is year (YYYY)

# Fourth command line argument is month (MM)

# Fifth command line argument is date of the month (DD)

# Sixth argument is for expanded output [expand,noexpand]

# Seventh argument is for including static files [static,nostatic]

# Eighth argument is DOMAIN_NAME

# Ninth argument is REMOTE_DATA_ROOT                                         

LAPS_DATA_ROOT=$1
time=$2
YYYY=$3
MM=$4
DD=$5
EXPAND=$6
STATIC=$7
#MINAGE=$8
#MAXAGE=$9
DOMAIN_NAME=$8
REMOTE_DATA_ROOT=$9

echo "LAPS_DATA_ROOT = $LAPS_DATA_ROOT"
echo "time = $time"
echo "EXPAND = $EXPAND"
echo "STATIC = $STATIC"
echo "MINAGE = $MINAGE"
echo "MAXAGE = $MAXAGE"

LOGDIR=$LAPS_DATA_ROOT/log

#Create list of lapsprd output files to be potentially tarred up for Web access
echo " "
echo "Create list of lapsprd output files to be potentially tarred up for Web access"

HH=`echo $time  | cut -c6-7`
HHMM=`echo $time  | cut -c6-9`
YYDDDHH=`echo $time  | cut -c1-7`
YYDDDHHMM=$time

echo "Tarring up LAPS in $LAPS_DATA_ROOT for $time"

cd $LAPS_DATA_ROOT

rm -f lapstar.txt
touch lapstar.txt

#LAPS Data Files
ls -1 time/systime.dat                                          > lapstar.txt
ls -1 time/c_time.dat                                          >> lapstar.txt

# Verification files                           
#find ./lapsprd/verif -type f         ! -cmin +120         -print | grep -v cont | grep -v hist >> lapstar.txt
#find ./lapsprd/verif -type l         ! -cmin +120         -print | grep -v cont | grep -v hist >> lapstar.txt
#echo "log/load.png"                                                                            >> lapstar.txt
#echo "log/cloud_fcst.png"                                                                      >> lapstar.txt

# Lapsprd files (except LGA, LGB, FUA, FSF)
find ./lapsprd -type f -name "$YYDDDHHMM.*"             -print | grep -v gr2 | grep -v www | grep -v "\.v" | grep -v "\.d"  >> lapstar.txt

# LGA/LGB files (use MINAGE/MAXAGE)
#find ./lapsprd/lg?     -name "*.lg?" ! -cmin +$MAXAGE -cmin +$MINAGE -print >> lapstar.txt

# Lapsprep files (use MINAGE/MAXAGE)
#find  ./lapsprd/lapsprep    -name "LAPS:$YYYY-$MM-$DD\_$HH"      -print >> lapstar.txt

# Pregenerated web analysis files
find ./lapsprd/www/anal2d/archive/* -type f -name "$YYDDDHHMM.*"     -print   >> lapstar.txt
find ./lapsprd/www/anal2d/recent/*  -type f -name "$YYDDDHHMM.*"     -print   >> lapstar.txt
find ./lapsprd/www/anal2d/recent    -type l ! -cmin +60              -print   >> lapstar.txt

# Wgi files
find ./log     -type f -name "*.wgi.$YYDDDHHMM" -print   >> lapstar.txt

# Sfc Verification file
ls -1 log/qc/laps_sfc.ver.$HHMM                          >> lapstar.txt

# Static files
ls -1 static/*.nl                                        >> lapstar.txt
ls -1 static/nest7grid.parms                             >> lapstar.txt
ls -1 static/purger.txt                                  >> lapstar.txt
if test "$STATIC" = static; then                         
    echo "including extra static files"
    ls -1 static/static.nest7grid                        >> lapstar.txt
    ls -1 static/www/*                                   >> lapstar.txt
else
    echo "not including extra static files"
fi

ls -l $LAPS_DATA_ROOT/lapstar.txt

which tar

if test "$EXPAND" = noexpand; then
    echo "current directory is `pwd`"
    echo "making tar file $LAPS_DATA_ROOT/laps_$time.tar.gz"

    echo "tar cvfz $LAPS_DATA_ROOT/laps_$time.tar.gz -T $LAPS_DATA_ROOT/lapstar.txt"
          tar cvfz $LAPS_DATA_ROOT/laps_$time.tar.gz -T $LAPS_DATA_ROOT/lapstar.txt

    ls -l $LAPS_DATA_ROOT/laps_$time.tar.gz

elif test "$EXPAND" = expand; then
    echo "cp to $LAPS_DATA_ROOT/lapstar_$YYDDDHH expanded directory"
    rm -rf $LAPS_DATA_ROOT/lapstar_*
    mkdir -p $LAPS_DATA_ROOT/lapstar_$YYDDDHH
    tar -T lapstar.txt -cf - | (cd $LAPS_DATA_ROOT/lapstar_$YYDDDHH;  tar xfBp -)
    pwd

elif test "$EXPAND" = scp; then
    echo "scp to remote location (under construction)"  

elif test "$EXPAND" = tar; then
    echo "tar to remote location (under construction)"  
    cd $LAPS_DATA_ROOT 
    tar -T lapstar.txt -cvf - | ssh $DOMAIN_NAME "cd $REMOTE_DATA_ROOT; tar xpvf -"

else
    echo "no tar performed"                                      

fi

date -u
    
echo "finished tarlapstime.sh"
