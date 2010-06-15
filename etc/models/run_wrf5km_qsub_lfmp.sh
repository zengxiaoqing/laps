#!/bin/sh --login

qsub /home/oplapb/builds/laps/etc/models/qsub_lfmpost_wrf5km.sh > /pan1/projects/mm5-laps/domains/WRFV3-5KM/log/lfmpost_run.log 2>&1
