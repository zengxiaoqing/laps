#!/usr/local/perl5/bin/perl

package laps_tools;
use strict;
umask 002;

sub laps_domain_name{
   my $LAPS_DATA_ROOT = shift(@_);
   my @components = split("/",$LAPS_DATA_ROOT); 
   my $i = 0; my $isave = 0;
   foreach (@components){
   if($_ eq "data"){
#print "location in list = $i\n";
      if($i != 1){$isave=$i-1;} }
      $i++;}
      $isave=$i-1 if($isave == 0);
   return my $DOMAIN_NAME = @components[$isave];
}
1;

sub laps_data_root{
   my $LAPS_DATA_ROOT = shift(@_);
   my $DATAROOT;
   my $DOMAIN_NAME = &laps_domain_name($LAPS_DATA_ROOT);
   my @components = split("/",$LAPS_DATA_ROOT);
   my $i = 0; my $isave=0;
   foreach (@components){
   if($_ eq "$DOMAIN_NAME"){
      if($i != 1){$isave=$i-1;} }
      $i++;}
   $isave=$i-1 if($isave == 0);
   $i=0;
   while ($i <= $isave) {
      $DATAROOT="$DATAROOT"."@components[$i]"."/";
      $i++;}
   return $DATAROOT;
}
1;

#===========================================================================
#
# laps_data_root is location of static, lapsprd, cdl, time,
# and log, and is where the data dirs are being created if
# they don't already exist.
# J. Smart 1-10-00:
#    "     1-18-00: product subdirectories added. Removed
#                   reference to lapssrcroot.

sub mkdatadirs{

  my $LAPS_DATA_ROOT = shift;
  $LAPS_DATA_ROOT = $ENV{LAPS_DATA_ROOT} if ! $LAPS_DATA_ROOT;

  my (@datadirs) = ('cdl','lapsprd','log','log/qc','static','time');
  my (@lapsprddirs) = ('l1s','lc3','lcb','lco','lcp','lct','lcv',
'lf1','lga','lh3','lh4','lhe','lil','liw','lm1','lm2','lmd','lmr','lmt',
'lpbl','lps','lq3','lrp','lrs','lso','lsx','lt1','lty',
'lvd','lvd/goes08','lvd/goes10','lvd/meteos','lvd/gmssat','lw3','lwc','lwm',
'msg','pig','pin','prg','pro','sag','vrc','snd',
'v01','v02','v03','v04','v05','v06','v07','v08','v09','v10','v11','v12',
'v13','v14','v15','v16','v17','v18','v19','vdr',
'd01','d02','d03','d04','d05','d06','d07','d08','d09','d10','d11','d12',
'd13','d14','d15','d16','d17','d18','d19','d20',
'ln3','lsr','lsr/dmsp01','lsr/dmsp02','lsr/goes08','lsr/goes10','lsr/tros12','lsr/tros14','cdw','rdr',
'rdr/001','rdr/002','rdr/003','rdr/004','rdr/005','rdr/006','rdr/007','rdr/008','rdr/009',
'rdr/001/raw','rdr/001/vrc','rdr/002/vrc','rdr/003/vrc','rdr/004/vrc',
'rdr/005/vrc','rdr/006/vrc','rdr/007/vrc','rdr/008/vrc','rdr/009/vrc',
'lgb','ls2','dprep','fsf','stats','balance','balance/lt1','balance/lw3','fua',
'grid','ram','rsf','lsq','tmg','lst','pbl','model','model/varfiles','model/output','model/sfc');

  my $datadirs; my $lapsprddirs;

  foreach (@datadirs){
     mkdir "$LAPS_DATA_ROOT/$_",0777 if(! -e "$LAPS_DATA_ROOT/$_");}

# this perhaps can be used once (if) the lapsprd subdirectories are checked-in to CVS.
#    if( -e "$LAPSSRCROOT/data/lapsprd")    {
#        opendir(DATADIRS,"$LAPSSRCROOT/data/lapsprd");
#        @lapsprddirs = readdir DATADIRS;   }

  foreach (@lapsprddirs) {
     mkdir "$LAPS_DATA_ROOT/lapsprd/$_",0777;}

  return;
}
1;
