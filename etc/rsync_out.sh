#!/bin/sh --login

#Copy recent output products to remote destination

#Argument 1 is local DATA_ROOT

#Argument 2 is remote DATA_ROOT (include the node unless arg 5 is set to 5r,6,7r,8r)

#Argument 3 tells whether to use 'qsub'

#Argument 4 is optional wall clock run time (e.g. 02:00)

#Argument 5 is optional number that controls the qsub action
#           2:  copy all except fua/fsa via rsync in prioritized sequence
#           5:  copy individual fua/fsf subdirectory via rsync
#           5r: reverse copy individual fua/fsf subdirectory via rsync
#           6:  copy individual fua/fsf subdirectory via scp (and remote purge)
#           7:  copy all except fua/fsa and verif via rsync in prioritized sequence
#           7r: reverse copy all except fua/fsa and verif via rsync in prioritized sequence
#           8:  copy just verif via rsync

#Argument 6 is optional subdirectory (if arg 5 is set to "5", "5r" or "6")

#Argument 7 is optional modelroot subdirectory (should be used if arg 5 is set to "6")

#Argument 8 is optional purge time (should be used if arg 5 is set to "6")

#Argument 9 is optional remote node (can be used if arg 5 is set to "6", default is 'clank')

umask 002

LOCAL_DATA_ROOT=$1
REMOTE_DATA_ROOT=$2

if test "$REMOTE_DATA_ROOT" == "`echo $REMOTE_DATA_ROOT | grep -v :`"; then # does not contain a semicolon
  RSH=""
else
  RSH=--rsh=ssh
fi

DELETE=--delete

if test "$5" == "8"; then # copy verif
  log=$LOCAL_DATA_ROOT/log/rsync.log.verif.`date +\%H\%M`
  DELETE=""
elif test "$5" == "7"; then 
  DELETE=""
  log=$LOCAL_DATA_ROOT/log/rsync.log.`date +\%H\%M`
elif test "$5" == "7r"; then 
  DELETE=""
  log=$LOCAL_DATA_ROOT/log/rsync.log.`date +\%H\%M`
elif test "$6" == ""; then # copy all
  log=$LOCAL_DATA_ROOT/log/rsync.log.`date +\%H\%M`
else
  subdir=$6
  log=$LOCAL_DATA_ROOT/log/rsync.log.$subdir.`date +\%H\%M`
fi

echo " log file is $log"

echo " " > $log

#if test "$5" != "4" && test "$5" != "5"; then # copy www related files
#    echo " rsync www related files to $REMOTE_DATA_ROOT " >> $log
#    rsync -rlptgvvz $RSH $DELETE $LOCAL_DATA_ROOT/lapsprd/www/*   $REMOTE_DATA_ROOT/lapsprd/www   > $LOCAL_DATA_ROOT/log/rsync_www.log.`date +\%H\%M` 2>&1
#    rsync -rlptgvvz $RSH $DELETE $LOCAL_DATA_ROOT/lapsprd/verif/* $REMOTE_DATA_ROOT/lapsprd/verif > $LOCAL_DATA_ROOT/log/rsync_verif.log.`date +\%H\%M` 2>&1
#    rsync -rlptgvvz $RSH $DELETE $LOCAL_DATA_ROOT/time/*          $REMOTE_DATA_ROOT/time          > $LOCAL_DATA_ROOT/log/rsync_time.log.`date +\%H\%M` 2>&1
#fi

echo " " >> $log
echo " rsync files to $REMOTE_DATA_ROOT " >> $log

echo " " >> $log
echo " RSH variable is: $RSH " >> $log

echo " " >> $log
echo " DELETE variable is: $DELETE " >> $log

echo " " >> $log
echo " Option 5 is: $5 " >> $log

date -u > $LOCAL_DATA_ROOT/time/rsynctime.dat

REMOTE_NODE=oplapb@clank   

if test "$3" = qsub; then

    HHMM=00:15
    if test "$4" != ""; then
        HHMM=$4
    fi

#   Set up qsub script
    if test "$subdir" == ""; then # copy all
        script=$LOCAL_DATA_ROOT/log/qsub_rsync_out.sh
        echo "#!/bin/sh"                 > $script
        echo "#PBS -N qsub_rsync_out"     >> $script
    else
        script=$LOCAL_DATA_ROOT/log/qsub_rsync_out_$subdir.sh
        echo "#!/bin/sh"                 > $script
        echo "#PBS -N rsync_$subdir"  >> $script
    fi
    echo "#PBS -A dlaps"           >> $script
    echo "#PBS -l procs=1,walltime=$HHMM:00"      >> $script
    echo "#PBS -m n"               >> $script
    echo "#PBS -q service"         >> $script
    if test "$subdir" == ""; then # copy all
        echo "#PBS -o $LOCAL_DATA_ROOT/log/qsub_rsync_out.log.`date +\%H\%M`"      >> $script
    else
        echo "#PBS -o $LOCAL_DATA_ROOT/log/rsync_qsub_fuafsf_$subdir.log.`date +\%H\%M`"   >> $script
    fi
    echo "#exit"                    >> $script
    echo " "                        >> $script
    echo "LOCAL_DATA_ROOT=$LOCAL_DATA_ROOT" >> $script
    echo "REMOTE_DATA_ROOT=$REMOTE_DATA_ROOT" >> $script

    echo "cd $LOCAL_DATA_ROOT"      >> $script

    echo "date -u"                  >> $script

    if test "$5" == ""; then # copy all (except www verif)
        echo " "                        >> $script
        echo "rsync -rlptgvvz --exclude='log/core' --exclude='time' --exclude='lapsprd/www' --exclude='lapsprd/verif' --exclude='lapsprd/lga' --exclude='lapsprd/bigfile' --exclude='lapsprd/lapsprep' $RSH $DELETE \$LOCAL_DATA_ROOT/* \$REMOTE_DATA_ROOT  > \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
    fi

    if test "$5" == "1"; then # copy all except fua
        echo " "                        >> $script
        echo "rsync -rlptgvvz --exclude='log/core' --exclude='time' --exclude='lapsprd/www' --exclude='lapsprd/verif' --exclude='lapsprd/lga' --exclude='lapsprd/bigfile' --exclude='lapsprd/lapsprep' --exclude='lapsprd/fua' $RSH $DELETE \$LOCAL_DATA_ROOT/* \$REMOTE_DATA_ROOT >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
    fi

    if test "$5" == "2"; then # copy all except fua/fsf
        echo " "                                                                   >> $script
        echo "date -u  > \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script

        echo " "                                                                   >> $script
        echo "Start copy of www directories"                                       >> $script
        echo "echo rsync -rlptgvvz $RSH $DELETE $LOCAL_DATA_ROOT/lapsprd/www/*                        $REMOTE_DATA_ROOT/lapsprd/www"                                                               >> $script
        echo "     rsync -rlptgvvz $RSH $DELETE $LOCAL_DATA_ROOT/lapsprd/www/*                        $REMOTE_DATA_ROOT/lapsprd/www   >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
        echo " "                                                                   >> $script
        echo "date -u >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script

        echo " "                                                                   >> $script
        echo "Start copy of verif directories"                                     >> $script
        echo "rsync -rlptgvvz $RSH         $LOCAL_DATA_ROOT/lapsprd/verif/* --exclude='REF/cont' $REMOTE_DATA_ROOT/lapsprd/verif >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
        echo " "                                                                   >> $script
        echo "date -u >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script

        echo " "                                                                   >> $script
        echo "Start copy of time directories"                                      >> $script
        echo "rsync -rlptgvvz $RSH $DELETE $LOCAL_DATA_ROOT/time/*                               $REMOTE_DATA_ROOT/time          >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
        echo " "                                                                   >> $script
        echo "date -u >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script

        echo " "                                                                   >> $script
        echo "Start copy of most overall dataroot directories"                     >> $script

#       --exclude-from=$LOCAL_DATA_ROOT/static/exclude.txt

        echo "echo rsync -rlptgvvz --exclude-from=$LOCAL_DATA_ROOT/static/exclude.txt $RSH $DELETE \$LOCAL_DATA_ROOT/* \$REMOTE_DATA_ROOT"                                                             >> $script
        echo "     rsync -rlptgvvz --exclude-from=$LOCAL_DATA_ROOT/static/exclude.txt $RSH $DELETE \$LOCAL_DATA_ROOT/* \$REMOTE_DATA_ROOT >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
        echo " "                                                                   >> $script
        echo "date -u >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script

    fi

    if test "$5" == "3"; then # copy just fua
        echo " "                        >> $script
        echo "rsync -rlptgvvz $RSH $DELETE \$LOCAL_DATA_ROOT/lapsprd/fua/* \$REMOTE_DATA_ROOT/lapsprd/fua >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
    fi

    if test "$5" == "4"; then # loop through fua/fsf
      cd $LOCAL_DATA_ROOT/lapsprd/fua
      pwd >> $log                                              
      for subdir in `ls`; do
        echo "build command for $subdir subdirectory" >> $log
        echo " "                        >> $script
        echo "rsync -rlptgvvz $RSH $DELETE \$LOCAL_DATA_ROOT/lapsprd/fua/$subdir \$REMOTE_DATA_ROOT/lapsprd/fua >> \$LOCAL_DATA_ROOT/log/rsync_qsub_fuafsf.log.`date +\%H\%M` 2>&1" >> $script
        echo "rsync -rlptgvvz $RSH $DELETE \$LOCAL_DATA_ROOT/lapsprd/fsf/$subdir \$REMOTE_DATA_ROOT/lapsprd/fsf >> \$LOCAL_DATA_ROOT/log/rsync_qsub_fuafsf.log.`date +\%H\%M` 2>&1" >> $script
      done
    fi

    if test "$5" == "5"; then # copy individual fua/fsf subdirectory via rsync
      cd $LOCAL_DATA_ROOT/lapsprd/fua
      pwd >> $log                                              
      echo "build command for $subdir subdirectory" >> $log
      echo " "                        >> $script
      echo "rsync -rlptgvvz $RSH \$LOCAL_DATA_ROOT/lapsprd/fsf/$subdir \$REMOTE_DATA_ROOT/lapsprd/fsf  > \$LOCAL_DATA_ROOT/log/rsync_qsub_fuafsf_$subdir.log.`date +\%H\%M` 2>&1" >> $script
      echo "rsync -rlptgvvz $RSH \$LOCAL_DATA_ROOT/lapsprd/fua/$subdir \$REMOTE_DATA_ROOT/lapsprd/fua >> \$LOCAL_DATA_ROOT/log/rsync_qsub_fuafsf_$subdir.log.`date +\%H\%M` 2>&1" >> $script
    fi

    if test "$5" == "5r"; then # reverse copy individual fua/fsf subdirectory via rsync
      cd $LOCAL_DATA_ROOT/lapsprd/fua
      pwd >> $log                                              
      echo "build command for $subdir subdirectory" >> $log
      echo " "                        >> $script
      echo "ssh clank.fsl.noaa.gov 'ssh pinky rsync -rlptgvvz $RSH $DELETE jetscp.rdhpcs.noaa.gov:$LOCAL_DATA_ROOT/lapsprd/fsf/$subdir $REMOTE_DATA_ROOT/lapsprd/fsf'  > \$LOCAL_DATA_ROOT/log/rsync_qsub_fuafsf_$subdir.log.`date +\%H\%M` 2>&1" >> $script
      echo "ssh clank.fsl.noaa.gov 'ssh pinky rsync -rlptgvvz $RSH $DELETE jetscp.rdhpcs.noaa.gov:$LOCAL_DATA_ROOT/lapsprd/fua/$subdir $REMOTE_DATA_ROOT/lapsprd/fua'  > \$LOCAL_DATA_ROOT/log/rsync_qsub_fuafsf_$subdir.log.`date +\%H\%M` 2>&1" >> $script
    fi

    if test "$5" == "6"; then # copy individual fua/fsf subdirectory via scp (and remote purge)                          
      MODEL_DATA_ROOT=$7
      REMOTE_PURGE_TIME=$8
      MODELTYPE=`echo $subdir | cut -c1-3`
      MODELCONFIG=`echo $subdir | cut -c5-10`

      if test "$9" != ""; then
        REMOTE_NODE=$9
      fi

      MODEL_CYCLE_TIME=`/usr/bin/perl /home/oplapb/builds/laps/etc/read_nl.pl -d $LOCAL_DATA_ROOT -n nest7grid.parms -v model_cycle_time`
      MODEL_INIT_TIME=`/usr/bin/perl /home/oplapb/builds/laps/etc/sched_sys.pl -c $MODEL_CYCLE_TIME -f yyyymmddhh`
      MODEL_FCST_INTVL=`/usr/bin/perl /home/oplapb/builds/laps/etc/read_nl.pl -d $LOCAL_DATA_ROOT -n nest7grid.parms -v model_fcst_intvl`
      N_FCST_STEPS_P1=`/usr/bin/perl /home/oplapb/builds/laps/etc/read_nl.pl -d $LOCAL_DATA_ROOT -n nest7grid.parms -v n_fcst_steps_p1`

      echo "MODEL_CYCLE_TIME = $MODEL_CYCLE_TIME" >> $log
      echo "MODEL_INIT_TIME  = $MODEL_INIT_TIME"  >> $log

      if test "$MODEL_CYCLE_TIME" > 21600; then
          SCP_WAIT_TIME=21600
      else
          SCP_WAIT_TIME=$MODEL_CYCLE_TIME
      fi

      echo "SCP_WAIT_TIME  = $SCP_WAIT_TIME"  >> $log

      MODEL_FCST_INTVL_MIN=`echo $MODEL_FCST_INTVL / 60 | bc`
      echo "MODEL_FCST_INTVL_MIN  = $MODEL_FCST_INTVL_MIN"  >> $log

      cd $LOCAL_DATA_ROOT/lapsprd/fua
      pwd >> $log                                              
      echo "build commands for $subdir subdirectory" >> $log
      echo " "                        >> $script

      echo "ssh $REMOTE_NODE /usr/nfs/lapb/builds/laps/etc/purger.pl -t $REMOTE_PURGE_TIME $REMOTE_DATA_ROOT/lapsprd/fua/$subdir           > \$LOCAL_DATA_ROOT/log/rsync_qsub_fuafsf_$subdir.log.`date +\%H\%M` 2>&1" >> $script
      echo "ssh $REMOTE_NODE /usr/nfs/lapb/builds/laps/etc/purger.pl -t $REMOTE_PURGE_TIME $REMOTE_DATA_ROOT/lapsprd/fsf/$subdir          >> \$LOCAL_DATA_ROOT/log/rsync_qsub_fuafsf_$subdir.log.`date +\%H\%M` 2>&1" >> $script
      echo " "                        >> $script

      echo "date -u"                  >> $script
      echo " "                        >> $script
#     lfmpost_scp_fuafsf.pl may need an option allowing the input of LAPS_DATA_ROOT or this script would need MODEL_DATA_ROOT
#                    perl /home/oplapb/builds/laps/etc/models/lfmpost_scp_fuafsf.pl -l oplapb@clank:/w3/jet/fab/wrf5km -m wrf -r /pan1/projects/mm5-laps/domains/WRFV3-5KM -f 37 -i 15 -w 10800               -a wsm6  
      echo "/usr/bin/perl /home/oplapb/builds/laps/etc/models/lfmpost_scp_fuafsf.pl -l $REMOTE_NODE:$REMOTE_DATA_ROOT -L      -r $MODEL_DATA_ROOT -d $MODEL_INIT_TIME      -f $N_FCST_STEPS_P1 -i $MODEL_FCST_INTVL_MIN -w $SCP_WAIT_TIME -m $MODELTYPE -a $MODELCONFIG >> \$LOCAL_DATA_ROOT/log/rsync_qsub_fuafsf_$subdir.log.`date +\%H\%M` 2>&1" >> $script
      echo " "                        >> $script

      echo "date -u"                  >> $script

#     exit
    fi

    if test "$5" == "7"; then # copy all in prioritized list except fua/fsf and verif
        echo " "                                                                   >> $script
        echo "date -u  > \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script

        echo " "                                                                   >> $script
        echo "Start copy of www directories"                                       >> $script
        echo "echo rsync -rlptgvvz $RSH $DELETE $LOCAL_DATA_ROOT/lapsprd/www/*                        $REMOTE_DATA_ROOT/lapsprd/www"                                                               >> $script
        echo "     rsync -rlptgvvz $RSH $DELETE $LOCAL_DATA_ROOT/lapsprd/www/*                        $REMOTE_DATA_ROOT/lapsprd/www   >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
        echo " "                                                                   >> $script
        echo "date -u >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script

        echo " "                                                                   >> $script
        echo "Start copy of time directory"                                        >> $script
        echo "rsync -rlptgvvz $RSH $DELETE $LOCAL_DATA_ROOT/time/systime.dat                          $REMOTE_DATA_ROOT/time          >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
        echo " "                                                                   >> $script
        echo "date -u >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script

        echo " "                                                                   >> $script
        echo "Start copy of most overall dataroot directories"                     >> $script

#       --exclude-from=$LOCAL_DATA_ROOT/static/exclude.txt

        echo "echo rsync -rlptgvvz --exclude-from=$LOCAL_DATA_ROOT/static/exclude.txt $RSH $DELETE \$LOCAL_DATA_ROOT/* \$REMOTE_DATA_ROOT"                                                             >> $script
        echo "     rsync -rlptgvvz --exclude-from=$LOCAL_DATA_ROOT/static/exclude.txt $RSH $DELETE \$LOCAL_DATA_ROOT/* \$REMOTE_DATA_ROOT >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
        echo " "                                                                   >> $script
        echo "date -u >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
    fi

    if test "$5" == "7r"; then # reverse copy all in prioritized list except fua/fsf and verif
        echo "under construction"            >> $log
        echo " "                                                                   >> $script
        echo "date -u  > \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script

        echo " "                                                                   >> $script
        echo "Start copy of www directories"                                       >> $script
        echo "echo rsync -rlptgvvz $RSH $DELETE $LOCAL_DATA_ROOT/lapsprd/www/*                $REMOTE_NODE:$REMOTE_DATA_ROOT/lapsprd/www"                                                               >> $script
        echo "     rsync -rlptgvvz $RSH $DELETE $LOCAL_DATA_ROOT/lapsprd/www/*                $REMOTE_NODE:$REMOTE_DATA_ROOT/lapsprd/www   >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
        echo " "                                                                   >> $script
        echo "date -u >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script

        echo " "                                                                   >> $script
        echo "Start copy of time directory"                                        >> $script
        echo "rsync -rlptgvvz $RSH $LOCAL_DATA_ROOT/time/systime.dat                          $REMOTE_NODE:$REMOTE_DATA_ROOT/time          >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
        echo " "                                                                   >> $script
        echo "date -u >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script

        echo " "                                                                   >> $script
        echo "Start reverse copy of most overall dataroot directories"             >> $script
        echo "rsync -rlptgvvz $RSH $LOCAL_DATA_ROOT/static/exclude.txt                        $REMOTE_NODE:$REMOTE_DATA_ROOT/static        >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
        echo "echo ssh clank.fsl.noaa.gov 'ssh pinky rsync -rlptgvvz --exclude-from=$REMOTE_DATA_ROOT/static/exclude.txt $RSH --timeout=500 jetscp.rdhpcs.noaa.gov:$LOCAL_DATA_ROOT/ $REMOTE_DATA_ROOT'                                                            " >> $script
        echo "     ssh clank.fsl.noaa.gov 'ssh pinky rsync -rlptgvvz --exclude-from=$REMOTE_DATA_ROOT/static/exclude.txt $RSH --timeout=500 jetscp.rdhpcs.noaa.gov:$LOCAL_DATA_ROOT/ $REMOTE_DATA_ROOT' >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
#       echo "echo rsync -rlptgvvz --exclude-from=$LOCAL_DATA_ROOT/static/exclude.txt $RSH $DELETE \$LOCAL_DATA_ROOT/* \$REMOTE_DATA_ROOT"                                                             >> $script
#       echo "     rsync -rlptgvvz --exclude-from=$LOCAL_DATA_ROOT/static/exclude.txt $RSH $DELETE \$LOCAL_DATA_ROOT/* \$REMOTE_DATA_ROOT >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
        echo " "                                                                   >> $script
        echo "date -u >> \$LOCAL_DATA_ROOT/log/rsync_qsub.log.`date +\%H\%M` 2>&1" >> $script
    fi

    if test "$5" == "8"; then # copy just verif
        echo " "                                                                   >> $script
        echo "date -u  > \$LOCAL_DATA_ROOT/log/rsync_qsub_verif.log.`date +\%H\%M` 2>&1" >> $script

        echo " "                                                                   >> $script
        echo "Start copy of verif directories"                                     >> $script
        echo "rsync -rlptgvvz $RSH $DELETE $LOCAL_DATA_ROOT/lapsprd/verif/* --exclude='REF/cont' --exclude='REF/hist' $REMOTE_DATA_ROOT/lapsprd/verif >> \$LOCAL_DATA_ROOT/log/rsync_qsub_verif.log.`date +\%H\%M` 2>&1" >> $script
        echo "rsync -rlptgvvz $RSH $DELETE $LOCAL_DATA_ROOT/log/load.png                             $REMOTE_DATA_ROOT/log           >> \$LOCAL_DATA_ROOT/log/rsync_qsub_verif.log.`date +\%H\%M` 2>&1" >> $script
        echo "rsync -rlptgvvz $RSH $DELETE $LOCAL_DATA_ROOT/log/cloud_fcst.png                       $REMOTE_DATA_ROOT/log           >> \$LOCAL_DATA_ROOT/log/rsync_qsub_verif.log.`date +\%H\%M` 2>&1" >> $script
        echo "rsync -rlptgvvz $RSH $DELETE $LOCAL_DATA_ROOT/time/modelvtime.dat                      $REMOTE_DATA_ROOT/time          >> \$LOCAL_DATA_ROOT/log/rsync_qsub_verif.log.`date +\%H\%M` 2>&1" >> $script
        echo " "                                                                   >> $script
        echo "date -u >> \$LOCAL_DATA_ROOT/log/rsync_qsub_verif.log.`date +\%H\%M` 2>&1" >> $script
    fi

    if test "$5" == "8r"; then # reverse copy just verif
        echo "under construction"            >> $log
    fi

    echo " "                                 >> $log
    echo " Running this qsub script...."     >> $log
    cat $script                              >> $log
    echo " "                                 >> $log
    echo " using this command..."            >> $log
    echo "qsub $script >> $log 2>&1"     >> $log

          qsub $script >> $log 2>&1      >> $log

else

    echo " non-qsub case with direct rsync... "         >> $log

    echo " rsync www related files to $REMOTE_DATA_ROOT " >> $log

    echo " rsync -rlptgvvz $RSH $DELETE $LOCAL_DATA_ROOT/lapsprd/www/*   $REMOTE_DATA_ROOT/lapsprd/www >> $log 2>&1" >> $log
           rsync -rlptgvvz $RSH $DELETE $LOCAL_DATA_ROOT/lapsprd/www/*   $REMOTE_DATA_ROOT/lapsprd/www >> $log 2>&1 

    echo " rsync -rlptgvvz $RSH $DELETE $LOCAL_DATA_ROOT/lapsprd/verif/* --exclude='REF/cont'  $REMOTE_DATA_ROOT/lapsprd/verif >> $log 2>&1" >> $log
           rsync -rlptgvvz $RSH $DELETE $LOCAL_DATA_ROOT/lapsprd/verif/* --exclude='REF/cont'  $REMOTE_DATA_ROOT/lapsprd/verif >> $log 2>&1 

    echo " rsync -rlptgvvz $RSH $DELETE $LOCAL_DATA_ROOT/time/*   $REMOTE_DATA_ROOT/time >> $log 2>&1" >> $log
           rsync -rlptgvvz $RSH $DELETE $LOCAL_DATA_ROOT/time/*   $REMOTE_DATA_ROOT/time >> $log 2>&1 

    echo " rsync non-www related files to $REMOTE_DATA_ROOT " >> $log

    echo " rsync -rlptgvvz --exclude-from=$LOCAL_DATA_ROOT/static/exclude.txt $RSH $DELETE $LOCAL_DATA_ROOT/* $REMOTE_DATA_ROOT >> $log 2>&1" >> $log
           rsync -rlptgvvz --exclude-from=$LOCAL_DATA_ROOT/static/exclude.txt $RSH $DELETE $LOCAL_DATA_ROOT/* $REMOTE_DATA_ROOT >> $log 2>&1         

fi

