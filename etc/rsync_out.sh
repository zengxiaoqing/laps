#!/bin/sh

#Copy recent output products to remote destination

LOCAL_DATA_ROOT=$1
REMOTE_DATA_ROOT=$2

rsync -rlptgvv --rsh=ssh --delete $LOCAL_DATA_ROOT/lapsprd/www/* $REMOTE_DATA_ROOT/lapsprd/www > $LOCAL_DATA_ROOT/log/rsync_www.log.`date +\%H\%M` 2>&1

log=$LOCAL_DATA_ROOT/log/rsync.log.`date +\%H\%M`

if test "$3" = qsub; then

#   Set up qsub script
    script=$LOCAL_DATA_ROOT/log/qsub_rsync_out.sh
    echo "#!/bin/sh"                 > $script
    echo "#$ -N qsub_rsync_out"     >> $script
    echo "#$ -A dlaps               >> $script
    echo "#$ -l h_rt=00:30:00"      >> $script
    echo "#$ -S /bin/sh"            >> $script
    echo "#$ -cwd"                  >> $script
    echo "#$ -pe service 1"         >> $script
    echo "#$ -o $LOCAL_DATA_ROOT/log/qsub_rsync_out.log"              >> $script
    echo "#exit"                    >> $script
    echo " "                        >> $script
    echo "LOCAL_DATA_ROOT=$LOCAL_DATA_ROOT" >> $script
    echo "REMOTE_DATA_ROOT=$REMOTE_DATA_ROOT" >> $script

    echo "rsync -rlptgvv --exclude='log/core' --exclude='lapsprd/www' --exclude='lapsprd/lga' --rsh=ssh --delete \$LOCAL_DATA_ROOT/* \$REMOTE_DATA_ROOT > \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script

    echo " "
    echo " Running this qsub script...."
    cat $script
    echo " "
    echo " using this command..."
    echo "qsub $script > $log 2>&1"

    qsub $script > $log 2>&1

else

    rsync -rlptgvv --exclude='log/core' --exclude='lapsprd/www' --exclude='lapsprd/lga' --rsh=ssh --delete $LOCAL_DATA_ROOT/* $REMOTE_DATA_ROOT > $log 2>&1         

fi

