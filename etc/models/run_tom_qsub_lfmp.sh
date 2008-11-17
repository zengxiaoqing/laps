#!/bin/sh --login

qsub /home/oplapb/builds/laps/etc/models/qsub_lfmpost_tom.sh > /lfs0/projects/hmtb/dwr_domains/wrf-tom/lfmpost_run.log 2>&1
