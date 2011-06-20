#!/bin/sh --login

#Copy recent output products to remote destination

LOCAL_DATA_ROOT=$1
REMOTE_DATA_ROOT=$2

rsync -rlptgvvz --rsh=ssh --delete $LOCAL_DATA_ROOT/lapsprd/www/*   $REMOTE_DATA_ROOT/lapsprd/www   > $LOCAL_DATA_ROOT/log/rsync_www.log.`date +\%H\%M` 2>&1
rsync -rlptgvvz --rsh=ssh --delete $LOCAL_DATA_ROOT/lapsprd/verif/* $REMOTE_DATA_ROOT/lapsprd/verif > $LOCAL_DATA_ROOT/log/rsync_verif.log.`date +\%H\%M` 2>&1

log=$LOCAL_DATA_ROOT/log/rsync.log.`date +\%H\%M`

echo " rsync non-www files to $REMOTE_DATA_ROOT " > $log

if test "$3" = qsub; then

    HHMM=00:15
    if test "$4" != ""; then
        HHMM=$4
    fi

#   Set up qsub script
    script=$LOCAL_DATA_ROOT/log/qsub_rsync_out.sh
    echo "#!/bin/sh"                 > $script
    echo "#$ -N qsub_rsync_out"     >> $script
    echo "#$ -A dlaps"              >> $script
    echo "#$ -l h_rt=$HHMM:00"      >> $script
    echo "#$ -S /bin/sh"            >> $script
    echo "#$ -cwd"                  >> $script
    echo "#$ -pe service 1"         >> $script
    echo "#$ -o $LOCAL_DATA_ROOT/log/qsub_rsync_out.log"              >> $script
    echo "#exit"                    >> $script
    echo " "                        >> $script
    echo "LOCAL_DATA_ROOT=$LOCAL_DATA_ROOT" >> $script
    echo "REMOTE_DATA_ROOT=$REMOTE_DATA_ROOT" >> $script

    echo "rsync -rlptgvvz --exclude='log/core' --exclude='lapsprd/www' --exclude='lapsprd/verif' --exclude='lapsprd/lga' --exclude='lapsprd/bigfile' --exclude='lapsprd/lapsprep' --rsh=ssh --delete \$LOCAL_DATA_ROOT/* \$REMOTE_DATA_ROOT >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script

    echo " "                                 >> $log
    echo " Running this qsub script...."     >> $log
    cat $script                              >> $log
    echo " "                                 >> $log
    echo " using this command..."            >> $log
    echo "/usr/local/fsl/bin/qsub_wait $script >> $log 2>&1"     >> $log

          /usr/local/fsl/bin/qsub_wait $script >> $log 2>&1      >> $log

else

    echo " non-qsub case with direct rsync... "         >> $log

    echo " rsync -rlptgvvz --exclude='log/core' --exclude='lapsprd/www' --exclude='lapsprd/verif' --exclude='lapsprd/lga' --exclude='lapsprd/bigfile' --exclude='lapsprd/lapsprep' --rsh=ssh --delete $LOCAL_DATA_ROOT/* $REMOTE_DATA_ROOT >> $log 2>&1" >> $log

           rsync -rlptgvvz --exclude='log/core' --exclude='lapsprd/www' --exclude='lapsprd/verif' --exclude='lapsprd/lga' --exclude='lapsprd/bigfile' --exclude='lapsprd/lapsprep' --rsh=ssh --delete $LOCAL_DATA_ROOT/* $REMOTE_DATA_ROOT >> $log 2>&1         

fi

