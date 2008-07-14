#!/bin/sh
# note this is a bourne shell script

# First command line argument is the LAPS_DATA_ROOT

# Second command line argument is the laps time in yydddhhmm format

LAPS_DATA_ROOT=$1
time=$2

LOGDIR=$LAPS_DATA_ROOT/log

#Create list of lapsprd output files to be potentially tarred up for Web access
echo " "
echo "Create list of lapsprd output files to be potentially tarred up for Web access"
hour=`echo $time  | cut -c6-7`

echo "Tarring up LAPS in $LAPS_DATA_ROOT for $time"

current_dir=`pwd`

echo "Current directory is $current_dir"

pushd $LAPS_DATA_ROOT

#LAPS Data Files
ls -1 time/*.dat                      >  $current_dir/lapstar.txt
ls -1 lapsprd/lga/* | tail -1         >> $current_dir/lapstar.txt
ls -1 lapsprd/lgb/* | tail -1         >> $current_dir/lapstar.txt
ls -1 lapsprd/v01/* | tail -1         >> $current_dir/lapstar.txt
ls -1 lapsprd/rdr/001/vrc/* | tail -1 >> $current_dir/lapstar.txt
ls -1 lapsprd/*/$time*                >> $current_dir/lapstar.txt
ls -1 lapsprd/*/*/$time*              >> $current_dir/lapstar.txt
ls -1 $LOGDIR/*.log.$hour*            >> $current_dir/lapstar.txt
ls -1 $LOGDIR/*.wgi.$hour*            >> $current_dir/lapstar.txt
ls -1 static/static.nest7grid         >> $current_dir/lapstar.txt
ls -1 static/*.nl                     >> $current_dir/lapstar.txt
ls -1 static/www/*                    >> $current_dir/lapstar.txt

tar cvf $current_dir/laps.tar -T $current_dir/lapstar.txt

popd

