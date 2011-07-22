#!/bin/sh --login

#Copy recent output products to remote destination

#Argument 6 is optional subdirectory (if arg 5 is specified)

LOCAL_DATA_ROOT=$1
REMOTE_DATA_ROOT=$2

if test "$6" == ""; then # copy all
  log=$LOCAL_DATA_ROOT/log/rsync.log.`date +\%H\%M`
else
  subdir=$6
  log=$LOCAL_DATA_ROOT/log/rsync.log.$subdir.`date +\%H\%M`
fi

echo " log file is $log"

echo " " > $log

#if test "$5" != "4" && test "$5" != "5"; then # copy www related files
#    echo " rsync www related files to $REMOTE_DATA_ROOT " >> $log
#    rsync -rlptgvvz --rsh=ssh --delete $LOCAL_DATA_ROOT/lapsprd/www/*   $REMOTE_DATA_ROOT/lapsprd/www   > $LOCAL_DATA_ROOT/log/rsync_www.log.`date +\%H\%M` 2>&1
#    rsync -rlptgvvz --rsh=ssh --delete $LOCAL_DATA_ROOT/lapsprd/verif/* $REMOTE_DATA_ROOT/lapsprd/verif > $LOCAL_DATA_ROOT/log/rsync_verif.log.`date +\%H\%M` 2>&1
#    rsync -rlptgvvz --rsh=ssh --delete $LOCAL_DATA_ROOT/time/*          $REMOTE_DATA_ROOT/time          > $LOCAL_DATA_ROOT/log/rsync_time.log.`date +\%H\%M` 2>&1
#fi

echo " " >> $log
echo " rsync non-www files to $REMOTE_DATA_ROOT " >> $log

if test "$3" = qsub; then

    HHMM=00:15
    if test "$4" != ""; then
        HHMM=$4
    fi

#   Set up qsub script
    if test "$subdir" == ""; then # copy all
        script=$LOCAL_DATA_ROOT/log/qsub_rsync_out.sh
        echo "#!/bin/sh"                 > $script
        echo "#$ -N qsub_rsync_out"     >> $script
    else
        script=$LOCAL_DATA_ROOT/log/qsub_rsync_out_$subdir.sh
        echo "#!/bin/sh"                 > $script
        echo "#$ -N qsub_rsync_out_$subdir"  >> $script
    fi
    echo "#$ -A dlaps"              >> $script
    echo "#$ -l h_rt=$HHMM:00"      >> $script
    echo "#$ -S /bin/sh"            >> $script
    echo "#$ -cwd"                  >> $script
    echo "#$ -pe service 1"         >> $script
    if test "$subdir" == ""; then # copy all
        echo "#$ -o $LOCAL_DATA_ROOT/log/qsub_rsync_out.log.`date +\%H\%M`"      >> $script
    else
        echo "#$ -o $LOCAL_DATA_ROOT/log/qsub_rsync_out_$subdir.log.`date +\%H\%M`"   >> $script
    fi
    echo "#$ -j y"                  >> $script
    echo "#exit"                    >> $script
    echo " "                        >> $script
    echo "LOCAL_DATA_ROOT=$LOCAL_DATA_ROOT" >> $script
    echo "REMOTE_DATA_ROOT=$REMOTE_DATA_ROOT" >> $script

    echo "cd $LOCAL_DATA_ROOT"      >> $script

    echo "date -u"                  >> $script

    if test "$5" == ""; then # copy all
        echo " "                        >> $script
        echo "rsync -rlptgvvz --exclude='log/core' --exclude='time' --exclude='lapsprd/www' --exclude='lapsprd/verif' --exclude='lapsprd/lga' --exclude='lapsprd/bigfile' --exclude='lapsprd/lapsprep' --rsh=ssh --delete \$LOCAL_DATA_ROOT/* \$REMOTE_DATA_ROOT >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
    fi

    if test "$5" == "1"; then # copy all except fua
        echo " "                        >> $script
        echo "rsync -rlptgvvz --exclude='log/core' --exclude='time' --exclude='lapsprd/www' --exclude='lapsprd/verif' --exclude='lapsprd/lga' --exclude='lapsprd/bigfile' --exclude='lapsprd/lapsprep' --exclude='lapsprd/fua' --rsh=ssh --delete \$LOCAL_DATA_ROOT/* \$REMOTE_DATA_ROOT >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
    fi

    if test "$5" == "2"; then # copy all except fua/fsf
        echo " "                        >> $script

        echo "rsync -rlptgvvz --rsh=ssh --delete $LOCAL_DATA_ROOT/lapsprd/www/*   $REMOTE_DATA_ROOT/lapsprd/www   >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
        echo "rsync -rlptgvvz --rsh=ssh --delete $LOCAL_DATA_ROOT/lapsprd/verif/* $REMOTE_DATA_ROOT/lapsprd/verif >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
        echo "rsync -rlptgvvz --rsh=ssh --delete $LOCAL_DATA_ROOT/time/*          $REMOTE_DATA_ROOT/time          >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script

        echo "rsync -rlptgvvz --exclude='log/core' --exclude='time' --exclude='lapsprd/www' --exclude='lapsprd/verif' --exclude='lapsprd/lga' --exclude='lapsprd/bigfile' --exclude='lapsprd/lapsprep' --exclude='lapsprd/fua' --exclude='lapsprd/fsf' --rsh=ssh --delete \$LOCAL_DATA_ROOT/* \$REMOTE_DATA_ROOT >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
#       echo "rsync -rlptgvvz --exclude='log/core'                                                                    --exclude='lapsprd/lga' --exclude='lapsprd/bigfile' --exclude='lapsprd/lapsprep' --exclude='lapsprd/fua' --exclude='lapsprd/fsf' --rsh=ssh --delete \$LOCAL_DATA_ROOT/* \$REMOTE_DATA_ROOT >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script

    fi

    if test "$5" == "3"; then # copy just fua
        echo " "                        >> $script
        echo "rsync -rlptgvvz --rsh=ssh --delete \$LOCAL_DATA_ROOT/lapsprd/fua/* \$REMOTE_DATA_ROOT/lapsprd/fua >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
    fi

    if test "$5" == "4"; then # loop through fua/fsf
      cd $LOCAL_DATA_ROOT/lapsprd/fua
      pwd >> $log                                              
      for subdir in `ls`; do
        echo "build command for $subdir subdirectory" >> $log
        echo " "                        >> $script
        echo "rsync -rlptgvvz --rsh=ssh --delete \$LOCAL_DATA_ROOT/lapsprd/fua/$subdir \$REMOTE_DATA_ROOT/lapsprd/fua >> \$LOCAL_DATA_ROOT/log/rsync_qsub_fuafsf.log.`date +\%H\%M` 2>&1" >> $script
        echo "rsync -rlptgvvz --rsh=ssh --delete \$LOCAL_DATA_ROOT/lapsprd/fsf/$subdir \$REMOTE_DATA_ROOT/lapsprd/fsf >> \$LOCAL_DATA_ROOT/log/rsync_qsub_fuafsf.log.`date +\%H\%M` 2>&1" >> $script
      done
    fi

    if test "$5" == "5"; then # copy individual fua/fsf subdirectory
      cd $LOCAL_DATA_ROOT/lapsprd/fua
      pwd >> $log                                              
      echo "build command for $subdir subdirectory" >> $log
      echo " "                        >> $script
      echo "rsync -rlptgvvz --rsh=ssh --delete \$LOCAL_DATA_ROOT/lapsprd/fsf/$subdir \$REMOTE_DATA_ROOT/lapsprd/fsf >> \$LOCAL_DATA_ROOT/log/rsync_qsub_fuafsf_$subdir.log.`date +\%H\%M` 2>&1" >> $script
      echo "rsync -rlptgvvz --rsh=ssh --delete \$LOCAL_DATA_ROOT/lapsprd/fua/$subdir \$REMOTE_DATA_ROOT/lapsprd/fua >> \$LOCAL_DATA_ROOT/log/rsync_qsub_fuafsf_$subdir.log.`date +\%H\%M` 2>&1" >> $script
    fi

    echo " "                                 >> $log
    echo " Running this qsub script...."     >> $log
    cat $script                              >> $log
    echo " "                                 >> $log
    echo " using this command..."            >> $log
    echo "/usr/local/fsl/bin/qsub $script >> $log 2>&1"     >> $log

          /usr/local/fsl/bin/qsub $script >> $log 2>&1      >> $log

else

    echo " non-qsub case with direct rsync... "         >> $log

    echo " rsync -rlptgvvz --exclude='log/core' --exclude='time' --exclude='lapsprd/www' --exclude='lapsprd/verif' --exclude='lapsprd/lga' --exclude='lapsprd/bigfile' --exclude='lapsprd/lapsprep' --rsh=ssh --delete $LOCAL_DATA_ROOT/* $REMOTE_DATA_ROOT >> $log 2>&1" >> $log

           rsync -rlptgvvz --exclude='log/core' --exclude='time' --exclude='lapsprd/www' --exclude='lapsprd/verif' --exclude='lapsprd/lga' --exclude='lapsprd/bigfile' --exclude='lapsprd/lapsprep' --rsh=ssh --delete $LOCAL_DATA_ROOT/* $REMOTE_DATA_ROOT >> $log 2>&1         

fi

