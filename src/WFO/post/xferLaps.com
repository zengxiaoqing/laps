#!/bin/perl

# Linda Wharton and Paul Schultz -- Feb 97
# Note that in all the system calls we invoke the bourne shell (sh)
# explicitly, because this perl script is being invoked from a csh.

# Import names of LAPS bigfile depositories.
$FXA_DATA = $ENV{"FXA_DATA"};
$FXA_DATA_BACKUP = $ENV{"FXA_DATA_BACKUP"};

# For local use.
$XFR_HOME="/usr/local/laps/makeWFObigfile";
$FCSTPRD="/usr/local/laps/laps/nest7grid/lapsprd";
$BIGFILE_PATH="Grid/FSL/netCDF/LAPS_Grid/LAPS";
$WFO_FCSTPRD="$FXA_DATA/$BIGFILE_PATH";
$WFO_FCSTPRD_BACKUP="$FXA_DATA_BACKUP/$BIGFILE_PATH";
$TRANS_TBL="$XFR_HOME/public_laps2wfo.tbl";
$LAPSTIME="/usr/local/laps/laps/nest7grid/sched/systime.dat";
$CDL_PATH="$FXA_DATA/Grid/FSL/CDL";

#Export environment variables for use in xfer_laps.
%ENV = (
       "XFR_HOME",    $XFR_HOME,
       "FCSTPRD",     $FCSTPRD,
       "WFO_FCSTPRD", $WFO_FCSTPRD,
       "TRANS_TBL",   $TRANS_TBL,
       "LAPSTIME",    $LAPSTIME,
       "CDL_PATH",    $CDL_PATH
       );

open(LAPSTIME,$LAPSTIME) or die "Can't open $LAPSTIME";
$filetime = <LAPSTIME>;
$dummy = <LAPSTIME>; 
$utc_hour = <LAPSTIME>; 
$utc_min = <LAPSTIME>; 
close(LAPSTIME);
chomp($filetime); chomp($utc_hour); chomp($utc_min);
$filetime =~ s/ *//;  #fortran write may have added leading spaces
$filetime= $filetime - 315619200;

chdir($XFR_HOME) or die "Cannot cd to $XFR_HOME";
system "sh './xfer_laps 1> $XFR_HOME/logs/xfer.log.$utc_hour$utc_min 2> $XFR_HOME/logs/xfer.err.$utc_hour$utc_min'";

open(LOGFILE, ">>$XFR_HOME/logs/xfer.log.$utc_hour$utc_min");
$oldhandle = select(LOGFILE);
print " \n";
select($oldhandle);
close(LOGFILE);

system "sh './create_filename $filetime >filename'";
open(FILENAME,"filename") or die "Can't open temporary output file for create_filename";
$filename = <FILENAME>;
chomp($filename);
close(FILENAME);
unlink("filename");

system "sh 'cp $WFO_FCSTPRD/$filename $WFO_FCSTPRD_BACKUP/$filename'";
system "sh '~fxa/bin/GridNotify.script Laps $filetime 0'";
