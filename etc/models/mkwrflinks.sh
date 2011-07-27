#!/bin/sh

# Set up links to locally run lfmpost on an external version of WRF

setenv DATETIME $1

setenv LOCAL_WRF_ROOT /lfs0/projects/hmtb/hwt_domains/hrrr_conus
setenv REMOTE_WRF_ROOT /whome/rtrr/hrrr
 
setenv LOCAL_WRF_RUN  $LOCAL_WRF_ROOT/$DATETIME
setenv REMOTE_WRF_RUN $REMOTE_WRF_ROOT/$DATETIME

cd $LOCAL_WRF_RUN/wrfprd
ln -s $REMOTE_WRF_RUN/wrfprd/wrfout* 

.

