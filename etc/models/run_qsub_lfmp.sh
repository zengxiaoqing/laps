#!/bin/sh --login

modelroot=$1    # Example: /lfs0/projects/hmtb/dwr_domains/efp_conus
model=$2        # Example: wrf (arw), arw, nmm, mm5, st4
physics=$3      # Example: tom
lbc=$4          # Example: gfs1
name=$5         # Example: ewp0 or wrf-wsm6 (should match a '$LAPS_DATA_ROOT/lapsprd/fua' subdirectory name and an element of 'nest7grid.parms/fdda_model_source')
nest=$6         # Example: 1 or -2
lapsdataroot=$7 # Example: /lfs0/projects/hmtb/dwr_domains/laps_psd   
fcstIncrMin=$8  # Example: 60
numFcsts=$9     # Example: 13

mdlfilewait=200

if test "$name" != "none"; then                               # (e.g. ewp0)
    run=$name
elif test "$physics" != "none" && test "$lbc" != "none"; then # (e.g. wrf-fer-gep0)
    run=$model\-$physics\-$lbc
elif test "$physics" != "none" && test "$lbc" == "none"; then # (e.g. wrf-fer)
    run=$model\-$physics
elif test "$physics" == "none" && test "$lbc" != "none"; then # (e.g. wrf-gep0)
    run=$model\-$lbc
elif test "$physics" == "none" && test "$lbc" == "none"; then # (e.g. wrf)
    run=$model\-$lbc
fi

echo "mkdir -p $modelroot"
      mkdir -p $modelroot

log=$lapsdataroot/log/lfmpost_run.log.`date +\%H\%M`

rm -f $lapsdataroot/log/qlfmp.log.`date +\%H\%M`

# Build qsub script
script=$modelroot/qsub_lfmpost.sh
echo "#!/bin/sh --login"         > $script
echo "#$ -N lfmp_$run"          >> $script
echo "#$ -A hmtb"               >> $script
echo "#$ -l h_rt=11:00:00"      >> $script
echo "#$ -S /bin/sh"            >> $script
echo "#$ -cwd"                  >> $script
echo "#$ -pe service 1"         >> $script
echo "#$ -o $lapsdataroot/log/qlfmp.log.`date +\%H\%M`"           >> $script
echo "#$ -j y"                  >> $script
echo "#exit"                    >> $script
echo " "                        >> $script

if test "$model" = arw; then
    echo "model=wrf"                                              >> $script
elif test "$model" = nmm; then
    echo "model=$model"                                           >> $script
elif test "$model" = wrf; then
    echo "model=$model"                                           >> $script
elif test "$model" = mm5; then
    echo "model=$model"                                           >> $script
else
    echo "model=$model"                                           >> $script
fi

echo "lbc=$lbc"                 >> $script
echo "physics=$physics"         >> $script
echo "nest=$nest"               >> $script
echo " "                        >> $script

echo "export LAPSROOT=/home/oplapb/builds_lahey/laps"           >> $script
echo "export LAPS_DATA_ROOT=$lapsdataroot"                      >> $script

echo "export NETCDF=/opt/netcdf/3.6.3-pgi"                >> $script
echo "export PATH=\$PATH:/opt/netcdf/3.6.3-pgi/bin"       >> $script
echo "export phys=$physics"                                     >> $script
echo " "                                                        >> $script
echo "fcstIncrMin=$fcstIncrMin"                                 >> $script
echo "numFcsts=$numFcsts"                                       >> $script
echo "maxWaitSec=21600"                                         >> $script
echo "maxHrsRun=24"                                             >> $script
echo "name=$name"                                               >> $script
echo "project=DWR"                                              >> $script
echo " echo 'Running this lfmpost.pl command...'"               >> $script
if test "$name" != "none"; then # (e.g. ewp0)
    echo " echo '/usr/bin/perl \$LAPSROOT/etc/lfmpost.pl -m \$model -r $modelroot -f \$numFcsts -i \$fcstIncrMin -w \$maxWaitSec -e \$maxHrsRun -n \$name -P \$project -g \$nest -W $mdlfilewait -q' " >> $script
    echo "       /usr/bin/perl \$LAPSROOT/etc/lfmpost.pl -m \$model -r $modelroot -f \$numFcsts -i \$fcstIncrMin -w \$maxWaitSec -e \$maxHrsRun -n \$name -P \$project -g \$nest -W $mdlfilewait -q"   >> $script
else
    echo " echo '/usr/bin/perl \$LAPSROOT/etc/lfmpost.pl -m \$model -r $modelroot -f \$numFcsts -i \$fcstIncrMin -w \$maxWaitSec -e \$maxHrsRun -y \$phys -P \$project -g \$nest -W $mdlfilewait -q' " >> $script
    echo "       /usr/bin/perl \$LAPSROOT/etc/lfmpost.pl -m \$model -r $modelroot -f \$numFcsts -i \$fcstIncrMin -w \$maxWaitSec -e \$maxHrsRun -y \$phys -P \$project -g \$nest -W $mdlfilewait -q"   >> $script
fi
echo " "                                                        >> $script
echo " "                                                        >> $script
echo " exit 0"                                                  >> $script

echo " "
echo " Running qsub script contained in $script...."
cat $script
echo " "
echo " using this command..."
echo "qsub $script > $log 2>&1"

qsub $script > $log 2>&1

