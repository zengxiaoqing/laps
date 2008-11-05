#!/bin/sh --login

qsub /home/oplapb/builds/laps/etc/qsub_lfmpost_nmm.sh > /lfs0/projects/hmtb/dwr_domains/wrf-nmm/lfmpost_run.log 2>&1
