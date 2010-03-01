#!/bin/ksh



#Do for each domain/model/time
#!/bin/sh

echo "Running followup_xsect.sh script" 
echo "user = "`whoami`                                    
echo "machine = "`uname -n`

# Arguments
domain=$1
model=$2
LAPS_DATA_ROOT=$3
datetime=$4
utc_hhmm=$5
# Hard coded stuff
WINDOW=0.0:0.0:1.0:1.0
LAPS_ETC=/usr/nfs/common/lapb/www/fcst2d
#LAPS_ETC=/usr/nfs/common/lapb/www/fcst2d/$domain/$model
WWW_DIR=/w3/lapb/domains/$domain
latest=latest
RESOLUTION=730x730
export EXE_DIR=$LAPSINSTALLROOT/bin
export LAPS_DATA_ROOT
export LAPS_ETC

mkdir -p $WWW_DIR/fcst2d
mkdir -p $WWW_DIR/fcst2d/$model

set -A products `cat $LAPS_DATA_ROOT/static/www/fcst2d/followup_xsect_prods.txt`

#Call for each product
for prod in ${products[*]}
do

  startdate=`date -u +%H:%M:%S`
  $LAPS_ETC/laps_gifs_sub_fcst.sh $prod $WINDOW $LAPS_ETC $WWW_DIR $model $LAPS_DATA_ROOT $latest $datetime $RESOLUTION $domain $model
  enddate=`date -u +%H:%M:%S`
  echo "Timing info for $prod: $startdate $enddate"
done
exit
