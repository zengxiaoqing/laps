#!/usr/local/apps/bin/perl

package laps_tools;
use strict;
umask 002;

#
# routine get_nl_values returns value of namelist variable given
# LAPS_DATA_ROOT -> location of namelists (-d opt or environ var);
# namelist filename -> commaind line (-n) input;
# namelist variable name -> command line (-v) input;
# namelist to open is a template or not (template is defined or not)
#          in such case LAPS_DATA_ROOT = path to template files.
#
# returned namelist value is in an array.
#
# J. Smart 4-28-00

sub get_nl_value{

   my ($namelist_file, $namelist_var, $LAPS_DATA_ROOT, $template) = @_;

   if(!defined $LAPS_DATA_ROOT) {$LAPS_DATA_ROOT = $ENV{LAPS_DATA_ROOT};}
   if(!defined $LAPS_DATA_ROOT) {print "LAPS_DATA_ROOT not defined "; exit;}

   defined $namelist_file || die "namelist filename required input\n";
   defined $namelist_var || die "namelist variable required input\n";

   if(!defined $template){
      open(NLF,"$LAPS_DATA_ROOT/static/$namelist_file") or die "Can't open $LAPS_DATA_ROOT/$namelist_file";
   }else{
      open(NLF,"$LAPS_DATA_ROOT/$namelist_file") or die "Can't open $LAPS_DATA_ROOT/$namelist_file";
   }

   my @nlf=<NLF>;
   my $nlf=@nlf;
   close (NLF);

   my @nlcomps;
   my @nl_values;
   my $nlcomps;
   my $search_for_equal = 0;
   my $continue_adding_lines = 1; 
   my $pattern = 0;

   foreach (@nlf){

      if($continue_adding_lines == 1){

         if($search_for_equal == 1){
            $pattern = grep /=/, $_;
            if($pattern != 1){
               $pattern = grep /^\s*\/\s*$/, $_;
            }
         }
                        
         chomp; s/\s+$//; s/\,$//;
         if(/^\s+$namelist_var\s*/i) {
           @nlcomps = split('=',$_);
           @nl_values = split /\,/,$nlcomps[1];
           $search_for_equal = 1;
         }elsif($search_for_equal == 1 && $pattern == 0){

                @nlcomps = split /\,/,$_;
                $nlcomps = @nlcomps;
                push @nl_values, split(',',$nlcomps[0]);

         }elsif($search_for_equal == 1 && $pattern == 1){
              $continue_adding_lines = 0;
         }
      }
   }
# strip off quotes and leading white space.
   foreach (@nl_values){
      s/\'//g; s/\s*//;}

   return @nl_values;
}
1;
#
#================================================================================
#
sub update_nl{
    my($LAPS_DATA_ROOT,$nl_file,$nl_var,@new_values) = @_;

    print "LAPS_DATA_ROOT = $LAPS_DATA_ROOT\n";
    print "nl file = $nl_file\n";
    print "nl var = $nl_var \n";
    foreach (@new_values){
       print "new values = $_\n";}
    print "\n";
#
# open and compare nest7grid.parms and *.nl in $srcroot/static and $dataroot/static
# Add variables and files found in $srcroot/static but not $dataroot/static 
#
# ***   Correction for this code ... do not add files (like from srcroot/static)
#
# Retain the values of variables found in $dataroot/static for variables in both files
#
# some of this software was from "laps" etc/laps_localization.pl and is used
# here to replace namelist variables with new values.
#
#
# save the original copy of $nl_file unless it already exists
#
    unless (-e "$LAPS_DATA_ROOT/static/$nl_file.bak"){
      system("cp $LAPS_DATA_ROOT/static/$nl_file $LAPS_DATA_ROOT/static/$nl_file.bak");}

# save the existing namelist file for later merger operation with template 
# ------------------------------------------------------------------------
    if(!-e "$LAPS_DATA_ROOT/static/tmp"){
       print "make tmp directory and save $nl_file \n";
       mkdir "$LAPS_DATA_ROOT/static/tmp", 0777 or die "Can't make directory $LAPS_DATA_ROOT/static/tmp";
       system("mv $LAPS_DATA_ROOT/static/$nl_file  $LAPS_DATA_ROOT/static/tmp/."); 
    }else{
       system("mv $LAPS_DATA_ROOT/static/$nl_file  $LAPS_DATA_ROOT/static/tmp/.");
    }

# write a "template" namelist file in static with new variable info
# -----------------------------------------------------------------
    my $filename = "$LAPS_DATA_ROOT/static/$nl_file";
    my $nl_line;

    if($nl_file eq "nest7grid.parms"){
      $nl_line = "lapsparms_nl";
    }else{
       my @fname_part = split /\./, $nl_file;
       $nl_line = $fname_part[0]."_".$fname_part[1];
    }

    &write_namelist($filename,$nl_var,$nl_line,@new_values);

    my($var, $val, $line, $eon);
    my %nl_vals;
    my %comments;

# open the namelist file for which the variable is to be replaced with a new value. 

        print "First pass of namelist parser $nl_file\n";
        open(FILE,"$LAPS_DATA_ROOT/static/$nl_file");
        my @template = <FILE>;
        close(FILE);

        $var='';
        my $mark=0;
        foreach $line (@template){
            if($line =~ /^\s*\&/){
                $mark=1;
                next;
            }elsif($line =~ /^\s*\//){
                $mark=2;
                next;
            }elsif($line =~ /^[!cC]/){
                $comments{$nl_file} .= $line;
                $mark=3;
                next;
            }elsif($line =~ /^\s*(\S+)\s*=\s*(.*)$/){
                $var = $1;
                $var =~ tr/a-z/A-Z/;
                $nl_vals{$var} = $2;
                next;
            }elsif($line =~ /^(.*)$/){
                $nl_vals{$var} .= "\n$1";
                next;
            }
            if($mark>0){
                $var = '';
                $mark=0;
            }
        }

        open(INFILE,"$LAPS_DATA_ROOT/static/tmp/$nl_file");
        my @infile = <INFILE>;
        close(INFILE);

        
#--- here only open the saved (in tmp) file that represent the original and
#    merge the new "template" just written with the saved namelist file in tmp

        print "merging $LAPS_DATA_ROOT/static/tmp/$nl_file into $LAPS_DATA_ROOT/static/$nl_file\n";
        open(OUTFILE,">$LAPS_DATA_ROOT/static/$nl_file") or die "Could not open $LAPS_DATA_ROOT/static/$nl_file to write";

        my @comments = split("\n",$comments{$nl_file});

        foreach $line (@infile){
            next if($line eq "\n");
#           print ">$line< ".length($line)."\n";
            if($line =~ /^\s*\//){
                print "End of namelist found\n";
                $eon = 1;
                next;
            }
            
            if($line =~ /^\s*(\S+)\s*=\s*(.*)$/){
                $var = $1;
                
                $var =~ tr/a-z/A-Z/;
                $val = $2;
                if(exists $nl_vals{$var}){
#                   print "Found $var = $val\n";
                    $val = $nl_vals{$var};
                }
                $val =~ s/\n$//;
                print OUTFILE " $var = $val\n";
                next;
            }elsif($line =~ /^[!cC]/){
                chomp($line);
                my $tmpline = $line;
                $tmpline =~ s/[(\[\]\\\/\(\)\!\$\^)]/\$1/g;
                next if(grep(/$tmpline/,@comments)>0);
                push(@comments,$line);
                next;
#               print OUTFILE $line;            ;
            }elsif($line =~ /^\s*&/){
                print OUTFILE $line;
                next;
            }elsif(($line =~ /^(\s*[^&\/].*)$/) && exists $nl_vals{$var}){
                next;
            }
            print OUTFILE $line;

        }
        print OUTFILE " \/\n";
        foreach(@comments){
            print OUTFILE "$_\n";
        }
        close(OUTFILE);
#   }

}

#=========================================================================
#
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
#
# =========================================================================
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
#    "     6-28-01: added wrfsi functionality using $domain_type argument.

sub mkdatadirs{

  my ($LAPS_DATA_ROOT,$LAPS_SRC_ROOT,$domain_type);

  $LAPS_DATA_ROOT = shift or die "no $LAPS_DATA_ROOT specified to mkdatadirs\n";
  $LAPS_SRC_ROOT  = shift or die "no $LAPS_SRC_ROOT input to mkdatadirs\n";
  $domain_type    = shift or die "no domain_type specified to mkdatadirs\n";

  if(! -e $LAPS_DATA_ROOT){die "$LAPS_DATA_ROOT does not exist\n";}
  if(! -e $LAPS_SRC_ROOT) {die "$LAPS_SRC_ROOT  does not exist\n";}

  if($domain_type ne "laps" && $domain_type ne "wrfsi"){
     die "Stop: Unknown domain type input to mkdatadirs = $domain_type\n";
  }

  my ($datadirs, $lapsprddirs);
  my (@datadirs, @lapsprddirs);
  my (@fua_dirs, @fsf_dirs);
  my (@fdda_dirs);
  my (@lga_dirs, @lgb_dirs);
  my (@bkgd_dirs);

  if($domain_type eq "laps"){

     (@datadirs) = qw (cdl lapsprd log log/qc static time);
     (@lapsprddirs) = qw (l1s lc3 lcb lco lcp lct lcv 
lf1  lh3 lh4 lhe lil liw lm1 lm2 lmd lmr lmt 
lpbl lps lq3 lrp lrs lso lsx lt1 lty lfr
lvd lvd/goes08 lvd/goes09 lvd/goes10 lvd/goes12 lvd/meteos lvd/gmssat
lw3 lwc lwm ctp msg pig pin prg pro sag vrc vrz snd 
v01 v02 v03 v04 v05 v06 v07 v08 v09 v10 v11 v12 
v13 v14 v15 v16 v17 v18 v19 vdr 
d01 d02 d03 d04 d05 d06 d07 d08 d09 d10 d11 d12 
d13 d14 d15 d16 d17 d18 d19 d20 
ln3 
lsr lsr/dmsp01 lsr/dmsp02 lsr/goes08 lsr/goes09 lsr/goes10 lsr/goes12
lsr/tros12 lsr/tros14 cdw rdr 
rdr/001 rdr/002 rdr/003 rdr/004 rdr/005 rdr/006 rdr/007 rdr/008 rdr/009 
rdr/001/vrc rdr/001/raw rdr/002/vrc rdr/002/raw rdr/003/vrc rdr/003/raw
rdr/004/vrc rdr/004/raw rdr/005/vrc rdr/005/raw rdr/006/vrc rdr/006/raw
rdr/007/vrc rdr/007/raw rdr/008/vrc rdr/008/raw rdr/009/vrc rdr/009/raw 
ls2 lapsprep lapsprep/mm5 lapsprep/rams lapsprep/wrf lapsprep/cdf 
dprep stats balance balance/lt1 balance/lw3 balance/lh3 balance/lq3 balance/air
grid ram rsf lsq tmg lst pbl model model/varfiles model/output model/sfc
verif verif/noBal verif/Bal verif/Bkgd);

     if(-e "$LAPS_DATA_ROOT/static/nest7grid.parms"){
        print "using LAPS_DATA_ROOT nest7grid.parms for fdda dirs\n";
        @fdda_dirs = &laps_tools::get_nl_value('nest7grid.parms','fdda_model_source_cmn',$LAPS_DATA_ROOT);
     }elsif(-e "$LAPS_SRC_ROOT/data"){
            print "using LAPS_SRC_ROOT nest7grid.parms for fdda dirs\n";
            @fdda_dirs = &laps_tools::get_nl_value('nest7grid.parms','fdda_model_source_cmn',"$LAPS_SRC_ROOT/data");
     }else{
           print "file nest7grid.parms not available in dataroot or source root\n";
           print "for acquiring fdda directories.  Terminating.\n";
           exit;
     }

     @bkgd_dirs = &get_bkgd_models($LAPS_SRC_ROOT);

     print "adding fdda_model_source subdirectories to lapsprd dirs\n";
     my $ii = 0;
     @fua_dirs[$ii] = 'fua';
     @fsf_dirs[$ii] = 'fsf';
     foreach (@fdda_dirs){
        if($_ ne "lga"){
              $ii++;
              @fua_dirs[$ii]=@fua_dirs[0]."/".$_;
              @fsf_dirs[$ii]=@fsf_dirs[0]."/".$_;
        }
     }
     print "fua dirs: @fua_dirs\n";
     print "fsf dirs: @fsf_dirs\n";

     print "adding background model subdirectories to lapsprd dirs\n";
     $lga_dirs[0]='lga';
     $lgb_dirs[0]='lgb';

     $ii=0;
     foreach (@bkgd_dirs){
              $ii++;
              $lga_dirs[$ii]=$lga_dirs[0]."/".$_;
              $lgb_dirs[$ii]=$lgb_dirs[0]."/".$_;
     }
     print "lga dirs: @lga_dirs\n";
     print "lgb dirs: @lgb_dirs\n";

  }else{
     (@datadirs) = qw (cdl siprd log static)
  }

  foreach (@datadirs){
     mkdir "$LAPS_DATA_ROOT/$_",0777 if(! -e "$LAPS_DATA_ROOT/$_");
  }

# this perhaps can be used once (if) the lapsprd subdirectories are checked-in to CVS.
#    if( -e "$LAPSSRCROOT/data/lapsprd")    {
#        opendir(DATADIRS,"$LAPSSRCROOT/data/lapsprd");
#        @lapsprddirs = readdir DATADIRS;   }
#        close(DATADIRS);

  foreach (@lapsprddirs) {
     mkdir "$LAPS_DATA_ROOT/lapsprd/$_",0777 if(! -e "$LAPS_DATA_ROOT/lapsprd/$_");
  }
  foreach (@lga_dirs) {
     mkdir "$LAPS_DATA_ROOT/lapsprd/$_",0777 if(! -e "$LAPS_DATA_ROOT/lapsprd/$_");
  }
  foreach (@lgb_dirs) {
     mkdir "$LAPS_DATA_ROOT/lapsprd/$_",0777 if(! -e "$LAPS_DATA_ROOT/lapsprd/$_");
  }
  foreach (@fua_dirs) {
     mkdir "$LAPS_DATA_ROOT/lapsprd/$_",0777 if(! -e "$LAPS_DATA_ROOT/lapsprd/$_");
  }
  foreach (@fsf_dirs) {
     mkdir "$LAPS_DATA_ROOT/lapsprd/$_",0777 if(! -e "$LAPS_DATA_ROOT/lapsprd/$_");
  }

  return;
}
1;
#
# ---------------------------------------------------------------------
sub get_bkgd_models{

  my ($SRCROOT) = @_;
  my $bgdata_inc = "$SRCROOT/src/include/bgdata.inc";
  open(BGD,$bgdata_inc) or die "Can't open $bgdata_inc file";
  my @bgdata_inc = <BGD>;
  my ($models,$i,@models);
  foreach (@bgdata_inc){
   if(/data/){
      $models = substr($_,26,length($_)); 
      chomp $models;
      @models = split(",",$models);
      for ($i=0; $i<=$#models; $i++){
         @models[$i]=~s/'//g;
      }
      last; 
   }
  }
# print "\n # of models = $#models + 1: @models\n";
  return @models;
}
1;
#
# ---------------------------------------------------------------------
sub write_namelist {

    my ($filename,$nl_var,$nl_line,@nl_values)=@_;
    open(NLF,">$filename");
    my $first_time = 1;
    print NLF " &".$nl_line."\n";
    foreach (@nl_values){
       print "value = $_\n";
       if($first_time == 1){
          if(/^\d+/ || /^\s*\.\D+\./ || /^\s*\-\d+/ || /^\'/){
             print NLF " ".$nl_var."=".$_.",";          #this for digits, .true./.false., neg #'s
             $first_time = 0;
          }elsif($first_time eq 1){
             print NLF " ".$nl_var."='".$_."',";        #this for character strings
             print NLF "\n" if(length($_)>25);          #separate long strings with line feed.
             $first_time = 0;
          }
       }elsif(/^\d+/ || /^\s*\.\D+\./ || /^\s*\-\d+/){  #this for namelist arrays
          print NLF $_.",";
       }else{
          print NLF "'".$_."',";
          print NLF "\n" if(length($_)>25);             #separate long strings with line feed.
       }
    }

    print NLF "\n /\n";
    close NLF;
}
#
#---------------------------------------------------------------------
#
sub get_pressures {

    my $LAPS_DATA_ROOT = shift(@_);
    my @pressures;
    if(-e "$LAPS_DATA_ROOT/static/pressures.nl")
    {
     open(PRES,"$LAPS_DATA_ROOT/static/pressures.nl") or die "Can't open $LAPS_DATA_ROOT/static/pressures.nl";
     my @plines = <PRES>;
     close PRES;
     my $i=0;
     foreach (@plines){
       if(/\s*(\d+)/){
          @pressures[$i]=$1;
          $i++
       }
     }
    }else{
     print "Warning: pressures.nl does not exist in $LAPS_DATA_ROOT/static\n";
    }
    return @pressures;
}
1;
#$i=0;
#foreach (@pressures){
#   print "$_\n";
#}
#
#
#---------------------------------------------------------------------
#
sub julian {
    my($yr,$mo,$dy) = @_;

    my($b,$g,$d,$e,$f,$today,$first_of_year);
    $yr = ($yr < 70) ? ($yr + 2000) : ($yr + 1900);

    # Use temporary vars to compute num of days since Oct 1, 1582 to today
    $b = int ( ($mo - 14) / 12 );
    $g = $yr + 4900 + $b;
    $b = $mo - 2 - 12*$b;

    $d = int( (1461*($g-100))/4);
    $e = int( (367*$b)/12);
    $f = int( (3*int($g/100))/4);

    $today = $d + $e - $f + $dy - 2432076;

    # Now compute number of days from Oct 1, 1582 to Jan 1, $yr
    $mo = 1;
    $dy = 1;
    $b = int( ($mo - 14) / 12);
    $g = $yr + 4900 + $b;
    $b = $mo - 2 - 12*$b;

    $d = int( (1461*($g-100))/4);
    $e = int( (367*$b)/12);
    $f = int( (3*int($g/100))/4);

    $first_of_year = $d + $e - $f + $dy - 2432076;

    # Julian day from 1st of year is $today-$first_of_year+1

    $today - $first_of_year + 1;
}
#
#---------------------------------------------------------------------
#
sub get_system_type {

    my $dataroot = shift(@_);

    if(!defined $dataroot){
       $dataroot = $ENV{LAPS_DATA_ROOT}  if( $ENV{LAPS_DATA_ROOT} );
       $dataroot = $ENV{LAPSINSTALLROOT} if( $ENV{LAPSINSTALLROOT} && !defined $dataroot);
       $dataroot = $ENV{LAPS_SRC_ROOT}   if( $ENV{LAPS_SRC_ROOT}   && !defined $dataroot);
       $dataroot = $ENV{MOAD_DATAROOT}   if( $ENV{MOAD_DATAROOT}   && !defined $dataroot);
       $dataroot = $ENV{INSTALLROOT}     if( $ENV{INSTALLROOT}     && !defined $dataroot);
       $dataroot = $ENV{SOURCE_ROOT}     if( $ENV{SOURCE_ROOT}     && !defined $dataroot);
    }
    print "The dataroot = $dataroot \n";

    my @dirs;
    opendir(DIR,$dataroot) or die "Can't open $dataroot";
    @dirs = readdir DIR;
    close(DIR);

    my $wrfsystem=0;
    my $lapssystem=0;

    foreach (@dirs) {
#      print "'$_',\n";
       if($_ eq "siprd"   || $_ eq "wrfprd"){$wrfsystem=1;}
       if($_ eq "lapsprd" || $_ eq "log")  {$lapssystem=1;}
       }
    if($wrfsystem==1 && $lapssystem==1){
       print "Found ambiguous system \n";
       print "Both siprd and lapsprd in dataroot \n";
    }elsif($wrfsystem==1){
       print "Welcome to WRF\n";
    }else{
       print "Welcome to LAPS\n";
    }
}
1;
#
#---------------------------------------------------------------------
# given coordinates of two places in radians, compute distance in meters
#
sub great_circle_distance {
    my ($lat1,$long1,$lat2,$long2) = @_;

    # approx radius of Earth in meters.  True radius varies from
    # 6357km (polar) to 6378km (equatorial).
    #my $earth_radius = 6367000;
#JS: modified to be in synch with the LAPS definition
    my $earth_radius = 6371200;

    my $dlon = $long2 - $long1;
    my $dlat = $lat2 - $lat1;
    my $a = (sin($dlat / 2)) ** 2 
            + cos($lat1) * cos($lat2) * (sin($dlon / 2)) ** 2;
    my $d = 2 * atan2(sqrt($a), sqrt(1 - $a));

    # This is a simpler formula, but it's subject to rounding errors
    # for small distances.  See http://www.census.gov/cgi-bin/geo/gisfaq?Q5.1
    # my $d = &acos(sin($lat1) * sin($lat2)
    #               + cos($lat1) * cos($lat2) * cos($long1-$long2));

    return $earth_radius * $d;
}
1;
#
#---------------------------------------------------------------------
#
sub date_to_i4time {
 
# This subroutine will accept as input the elements of the date (i.e.
# day, month, year, hours, minutes, seconds)
# and return the corresponding integer time in seconds since
# January 1, 1970. (Brian Jamison)
#
# example call:
#
#   $yr = 2000;
#   $mo = 6;
#   $dy = 3;
#   $hr = 12;
#   $mn = 0;
#   $sc = 0;
#   ($i4time) = &date_to_i4time($yr,$mo,$dy,$hr,$mn,$sc);
#   print "after calling sub, i4time is $i4time\n";
#
# Notes on input:
#
#   year ($yr)  -  must be 4 digit integer (i.e. 1992 instead of 92)
#
 
    my ($year,$month,$day,$hours,$minutes,$seconds)=@_;
    my ($i4time);
    my ($start,$jul,$jul_minus_1);
    my ($secs_per_day,$secs_per_hour,$secs_per_minute);
    my ($ly);
 
    $start = 1970;
    $i4time = 0;
 
    while ($start != $year) {
     
          ($ly) = &leapyear_tf($start);
 
          if ($ly) {
              $i4time = $i4time + 31622400;
          } else {
              $i4time = $i4time + 31536000;
          }
 
          $start++;
       
    }
 
#
# Get the julian day to add the seconds per day
#
    ($jul) = &get_julian_day($day,$month,$year);
 
    $jul_minus_1 = $jul - 1;
    $secs_per_day = $jul_minus_1 * 86400;
    $i4time = $i4time + $secs_per_day;
#
# Add in the seconds per hour, seconds per minute, and the seconds
#
    $secs_per_hour = $hours * 3600;
    $secs_per_minute = $minutes * 60;
    $i4time = $i4time + $secs_per_hour + $secs_per_minute + $seconds;

    return ($i4time);
 
}
#
#-------------------------------------------------------------------
#
sub i4time_to_date {
   
# This subroutine will accept as input the integer time in seconds since
# January 1, 1970 and return the elements of the date (i.e.
# day, month, year, hours, minutes, seconds) (Brian Jamison)
#
# example call:
#
#   $i4time = 1055894400;
#   ($yr,$mo,$dy,$hr,$mn,$sc) = &i4time_to_date($i4time);
#   print "after calling sub,\n
#            yr is $yr\n
#            mo is $mo\n
#            dy is $dy\n
#            hr is $hr\n
#            mn is $mn\n
#            sc is $sc\n";
#
# Notes on output:
#
#   year ($yr)  -  will be output as a 4 digit integer
#
                                                                                
  my ($i4time)=@_;
  my ($year,$month,$day,$hours,$minutes,$seconds);
  my ($count,$lastcount,$jday);
  my ($start,$jul,$jul_minus_1);
  my ($secs_per_day,$secs_per_hour,$secs_per_minute);
  my ($ly);
                                                                                
  my @nmonth = qw(31 28 31 30 31 30 31 31 30 31 30 31);
  my @nmonth_ly = qw(31 29 31 30 31 30 31 31 30 31 30 31);
   
  $count = 0;
  $year = 1970;
   
  while ($i4time >= $count) {
       
    ($ly) = &leapyear_tf($year);
   
    if ($ly) {
      $count += 31622400;
      $lastcount = 31622400;
    } else {
      $count += 31536000;
      $lastcount = 31536000;
    }
   
    $year++;
         
  }
                                                                                
  $count -= $lastcount;
  $year--;
   
  my $daycount = 0;
   
  while ($i4time >= $count) {
   
    $daycount++;
    $count += 86400;
   
  }
   
  $jday = $daycount;
  $daycount--;
  $count -= 86400;
   
  $month = 0;
  ($ly) = &leapyear_tf($year);
   
  while ($daycount >= 0) {
   
    $month++;
   
    if ($ly) {
      $daycount -= $nmonth_ly[$month - 1];
      if ($daycount <= 0) {
        $day = $nmonth_ly[$month - 1] + $daycount + 1;
      }
    } else {
      $daycount -= $nmonth[$month - 1];
      if ($daycount <= 0) {
        $day = $nmonth[$month - 1] + $daycount + 1;
      }
    }
  }
                                                                                
  $hours = 0;
   
  while ($i4time >= $count) {
    $count += 3600;
    $hours++;
  }
     
  $count -= 3600;
  $hours--;
   
  $minutes = 0;
   
  while ($i4time >= $count) {
    $count += 60;
    $minutes++;
  }
   
  $count -= 60;
  $minutes--;
   
  $seconds = $i4time - $count;
                                                                                
  return ($year,$month,$day,$hours,$minutes,$seconds);
   
}
#
#-------------------------------------------------------------------
#
sub leapyear_tf {
 
# This subroutine will accept as input the 4 digit year and
# return the variable "ly" which will have a value of 1 for leap years
# and 0 for other years. (Brian Jamison)
 
  my ($year)=@_;
  my ($yeardiv4,$yeardiv100,$yeardiv400);
 
# Initialize some logical variables to be false (0=false,1=true)
 
  $yeardiv4 = 0;
  $yeardiv100 = 0;
  $yeardiv400 = 0;
 
# Test to see if the year is a leap year
# Leap year definition: If the year is evenly divisible by 4, it is a
# leap year unless it is a centenary year (i.e. 1800,1900, etc.).  However
# centenary years evenly divisible by 400 are leap years (e.g. 2000).
 
  if ($year % 4 == 0) {$yeardiv4 = 1;}
  if ($year % 100 == 0) {$yeardiv100 = 1;}
  if ($year % 400 == 0) {$yeardiv400 = 1;}
 
  my $ly = 0;
  if ($yeardiv4) {
    $ly = 1;
    if ($yeardiv100 && !$yeardiv400) {$ly = 0;}
  }
 
  return ($ly);
 
}
# 
#-------------------------------------------------------------------
#
sub get_julian_day {
 
# This subroutine will accept as input the day, month, and year and
# return the corresponding julian day. (Brian Jamison)
#
# Notes on input:
#
#   iyear  -  must be 4 digit integer (i.e. 1992 instead of 92)
 
  my ($iday,$imonth,$iyear)=@_;
 
  my ($ly) = &leapyear_tf($iyear);
 
  my @noleap = (0,31,59,90,120,151,181,212,243,273,304,334);
  my @leap = (0,31,60,91,121,152,182,213,244,274,305,335);
 
  my $ijul = $iday + $noleap[$imonth-1];
  if ($ly) {$ijul = $iday + $leap[$imonth-1]};
 
  return ($ijul);
 
}
1;
#
# ------------------------------------------------------------------------
#
sub systime{
    use Time::Local;
    my($DATAROOT,$delay,$cycle_time,$archive_time,$write_systime_dat) = @_;

#   See if -t option was used to pass in archive time in calling routine
    my $ctime = 0;    # Initial declaration outside the scope of the if test
    if($archive_time > 0) {
#       print "Archive data time is being set to a ctime of $archive_time...\n";
        $ctime = $archive_time;
    }else{
#       print "Setting ctime based on clock time and delay of $delay in hours...\n";
        $ctime = time - $delay*3600;
    }

#   print "ctime = $ctime\n";

    my @MON = qw(JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC);
    my($sec,$min,$hour,$mday,$mon,$year,$yyyy); 
    if($archive_time > 0){
      ($year,$mon,$mday,$hour,$min,$sec) = &laps_tools::i4time_to_date($archive_time);
       $yyyy=$year;
       $year=substr($year,2,2);
       $mon=$mon-1;
    }else{
       ($sec,$min,$hour,$mday,$mon,$year) = gmtime($ctime);
        $yyyy = 1900+$year;
        $year= $year-100 if($year>99);
        $year='0'.$year if(length($year)<2);
        $mday = '0'.$mday if(length($mday)<2);
    }
#
# This resets to the top of the cycle
# ---------------------------------------
    if($archive_time==0){
       my $minute=$min;
       my $i=0;
       while($i<60){
       my $mod = $i%int($cycle_time/60);
        $min = $i if($i<$minute && $mod==0);
        $i++;
        #     print "$i $mod\n";
       }
    }

    $ctime = timegm(0,$min,$hour,$mday,$mon,$yyyy-1900);

    my $ftime = $ctime +  315619200;
    $min = '0'.$min if(length($min)<2);
    $hour = '0'.$hour if(length($hour)<2);

    my $jjj = &laps_tools::julian($year,$mon+1,$mday);

    $jjj="0".$jjj while(length($jjj)< 3);

    my $yyjjjhhmm = "$year$jjj$hour$min";

    if(defined $write_systime_dat){
       open(TFILE,">$DATAROOT/time/c_time.dat");
       print TFILE " $year$jjj$hour$min\n";
       print TFILE "   $ctime\n";
       close(TFILE);

       open(TFILE,">$DATAROOT/time/systime.dat");
       print TFILE "  $ftime\n";
       print TFILE " $year$jjj$hour$min\n";
       print TFILE "$hour\n";
       print TFILE "$min\n";
       print TFILE "$mday-$MON[$mon]-$yyyy $hour$min\n";
       print TFILE "$year$jjj\n";
       close(TFILE);
    }

    return ($yyjjjhhmm);
}
