#!/bin/sh --login

projectpath=$1 # Example: /lfs0/projects/hmtb/dwr_domains
model=$2       # Example: wrf
physics=$3     # Example: tom
lbc=$4         # Example: gfs1
nest=$5        # Example: 1 or -2

run=$model\-$physics\-$lbc
mkdir -p $projectpath/$run

log=$projectpath/$run/lfmpost_run.log

# Build qsub script
script=$projectpath/$run/qsub_lfmpost.sh
echo "#!/bin/sh"                 > $script
echo "#$ -N lfmp_$run"          >> $script
echo "#$ -A hmtb"               >> $script
echo "#$ -l h_rt=04:30:00"      >> $script
echo "#$ -S /bin/sh"            >> $script
echo "#$ -cwd"                  >> $script
echo "#$ -pe wcomp 1"           >> $script
echo "#$ -o $projectpath/$run/qlfmp.log"                          >> $script
echo "#exit"                    >> $script
echo " "                        >> $script
echo "projectpath=$projectpath" >> $script
echo "model=$model"             >> $script
echo "lbc=$lbc"                 >> $script
echo "physics=$physics"         >> $script
echo "nest=$nest"               >> $script
echo " "                        >> $script
#cat /home/oplapb/builds/laps/etc/models/qsub_lfmpost.sh.footer >> $script
echo "export LAPSROOT=/home/oplapb/builds/laps"                 >> $script
echo "export LAPS_DATA_ROOT=$projectpath/laps"                  >> $script
echo "export WRF_DATAROOT=$projectpath/$run"                    >> $script
echo "export NETCDF=/opt/netcdf/3.6.2-pgi-7.1-3"                >> $script
echo "export PATH=\$PATH:/opt/netcdf/3.6.2-pgi-7.1-3/bin"       >> $script
echo "export phys=$physics"                                     >> $script
echo " "                                                        >> $script
echo "numFcsts=40"                                              >> $script
echo "fcstIncrMin=180"                                          >> $script
echo "maxWaitSec=3600"                                          >> $script
echo "maxHrsRun=4"                                              >> $script
echo "project=DWR"                                              >> $script
echo " echo 'Running this lfmpost.pl command...'"               >> $script
echo " echo '/usr/bin/perl \$LAPSROOT/etc/lfmpost.pl -m $model -f \$numFcsts -i \$fcstIncrMin -w \$maxWaitSec -e \$maxHrsRun -y \$phys -P \$project -g \$nest -q' " >> $script
echo "       /usr/bin/perl \$LAPSROOT/etc/lfmpost.pl -m $model -f \$numFcsts -i \$fcstIncrMin -w \$maxWaitSec -e \$maxHrsRun -y \$phys -P \$project -g \$nest -q"   >> $script
echo " "                                                        >> $script
echo " "                                                        >> $script
echo " exit 0"                                                  >> $script

echo " "
echo " Running this qsub script...."
cat $script
echo " "
echo " using this command..."
echo "qsub $script > $log 2>&1"

qsub $script > $log 2>&1

