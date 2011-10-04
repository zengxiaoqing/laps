#!/bin/ksh



#Do for each domain/model/time
#!/bin/sh

echo " "
echo "Starting followup_fcst_xsect.sh script" 
echo "user = "`whoami`                                    
echo "machine = "`uname -n`

# Arguments
DOMAIN=$1
MODEL=$2
LAPS_DATA_ROOT=$3
DATETIME=$4
LAPSINSTALLROOT=$5

# Hard coded stuff
WINDOW=0.0:0.0:1.0:1.0
LAPS_ETC=$LAPSINSTALLROOT/etc
#LAPS_ETC=/usr/nfs/common/lapb/www/fcst2d
#LAPS_ETC=/usr/nfs/common/lapb/www/fcst2d/$domain/$MODEL

if test -r "/w3/lapb"; then
    WWW_DIR=/w3/lapb/domains/$DOMAIN
else
    WWW_DIR=$LAPS_DATA_ROOT/lapsprd/www
fi

latest=latest
RESOLUTION=730x730
export LAPS_DATA_ROOT
export LAPS_ETC

echo "LAPS_DATA_ROOT=$LAPS_DATA_ROOT"
echo "MODEL=$MODEL"

mkdir -p $WWW_DIR/fcst2d
mkdir -p $WWW_DIR/fcst2d/$MODEL

set -A products `cat $LAPS_DATA_ROOT/static/www/fcst2d/followup_xsect_prods.txt`

echo " "
echo "Start loop through products at $DATETIME"

#Call for each product
for prod in ${products[*]}
do

  startdate=`date -u +%H:%M:%S`
  echo "Running: $LAPS_ETC/www/lapsplot/laps_gifs_sub_fcst.sh $prod $WINDOW $LAPS_ETC $WWW_DIR $MODEL $LAPS_DATA_ROOT $LAPSINSTALLROOT $DATETIME $RESOLUTION"
                 $LAPS_ETC/www/lapsplot/laps_gifs_sub_fcst.sh $prod $WINDOW $LAPS_ETC $WWW_DIR $MODEL $LAPS_DATA_ROOT $LAPSINSTALLROOT $DATETIME $RESOLUTION
  enddate=`date -u +%H:%M:%S`
  echo "Timing info for $prod: $startdate $enddate"
done

echo " "
echo "Finished loop through followup_fcst_xsect.sh products at $DATETIME - end of script"

exit
