#!/bin/sh --login

qsub /home/oplapb/builds/laps/etc/models/qsub_lfmpost_psd.sh > /lfs0/projects/hmtb/dwr_domains/psd_gfs/lfmpost_run.log 2>&1
