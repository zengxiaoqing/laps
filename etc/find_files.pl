#!/usr/bin/perl

use strict;
use English;
use Time::Local;
use Time::JulianDay;
use Date::Manip;
use vars qw($opt_d $opt_t $opt_x $opt_p $opt_i);

use Getopt::Std;

getopts('d:c:a:rm:i:t:e:q:Q:KV:M:n:A:f:F:wslu:vp:TE:P:S:XR:x:y:C:NY:');

my $numfiles = 0;
my $filenamelength = 0;
my @files_found = ();
my $dt = 0;
my ($year_t,$month_t,$day_t,$hour_t,$minute_t) = (0,0,0,0,0);
my $pattern = "";
my $maxfilelength = 0;

#Time pattern:
if ($opt_p =~ "yyyymmddhhmm") {       # standard date
  $pattern = "\\d\\d\\d\\d\\d\\d\\d\\d\\d\\d\\d\\d";
  if ($opt_t =~ m/(\d\d\d\d)(\d\d)(\d\d)(\d\d)(\d\d)/) {
    ($year_t,$month_t,$day_t,$hour_t,$minute_t) = ($1,$2,$3,$4,$5);
  }
} elsif ($opt_p =~ "yyyydddhhmm")  { # Days in the year
  $pattern = "\\d\\d\\d\\d\\d\\d\\d\\d\\d\\d\\d";
  if ($opt_t =~ m/(\d\d\d\d)(\d\d\d)(\d\d)(\d\d)/) {
    ($year_t,$day_t,$hour_t,$minute_t) = ($1,$2,$3,$4);

    # Re-calculate standard date:
    ($year_t,$month_t,$day_t) = Date_NthDayOfYear($year_t,$day_t);
  }
} elsif ($opt_p =~ "yydddhhmm")  { # 2 digits for year and Days in the year
  $pattern = "\\d\\d\\d\\d\\d\\d\\d\\d\\d";
  if ($opt_t =~ m/(\d\d)(\d\d\d)(\d\d)(\d\d)/) {
    ($year_t,$day_t,$hour_t,$minute_t) = ($1,$2,$3,$4);
    # Assuming this code is used between year 2000 to 2999:
    $year_t = 2000+$year_t;

    # Re-calculate standard date:
    ($year_t,$month_t,$day_t) = Date_NthDayOfYear($year_t,$day_t);
  }
} else {
  print "Pattern is not implemented!\n";
  exit;
}

# Epoach time of opt_t:
# timelocal uses month of (0-11) and years since 1900:
my $epoch_t = timelocal(0,$minute_t,$hour_t,$day_t,$month_t-1,$year_t-1900);

# Open file for outputing filenames within the given time window:
open(OUT,">/tmp/temp.files");

# Check director existing:
if (defined $opt_d && -d $opt_d) {
  opendir (DIR, $opt_d) or die $!;
  while (my $file = readdir(DIR)) {
    if ($file =~ m/$opt_x/) {
      if ($file =~ m/($pattern)/) {       # () matching allows saving the matching string in $1
        my ($year,$month,$day,$hour,$minute) = get_file_time($1,$opt_p);

        if ($year =~ m/\d\d\d\d/) {
          # Convert to Epoch time: No seconds counted: 0 for the first argument:
          # timelocal uses month of (0-11) and years since 1900:
          my $epoch_f = timelocal(0,$minute,$hour,$day,$month-1,$year-1900);

          $dt = $epoch_f-$epoch_t;
          if ($dt < 0) {
            $dt = -$dt;
          }

          if ($dt <= $opt_i) {
            $numfiles = $numfiles+1;
            push(@files_found,$file);
            if (length($file) > $maxfilelength) {
              $maxfilelength = length($file);
            }
            "";
          }
        }
      }
    } 
  }
} else {
  print "Option is not defined or search dir does not exist\n";
  "";
}

# Output:
print OUT "$numfiles $maxfilelength\n";
my $i=0;
while ($i < $numfiles) {
  print OUT "@files_found[$i]\n";
  $i++;
}

close(OUT);

# This function is to parse $1 for year, month, day, hour, minute.
# Usage: get_file_time($1, $2) : $1 is the string to parse; $2 is the pattern to decode time
sub get_file_time {

  if ($_[1] eq "yyyymmddhhmm") { #yyyymmddhhmm
    if ($_[0] =~ m/(\d\d\d\d)(\d\d)(\d\d)(\d\d)(\d\d)/) {
      ($1,$2,$3,$4,$5);
    }
  } elsif ($_[1] eq "yyyydddhhmm") { #yyyydddhhmm
    if ($_[0] =~ m/(\d\d\d\d)(\d\d\d)(\d\d)(\d\d)/) {
      my ($y, $m, $d) = Date_NthDayOfYear($1,$2);
      ($y, $m, $d, $3, $4);
    }
  } elsif ($_[1] eq "yydddhhmm") { #yyyydddhhmm
    if ($_[0] =~ m/(\d\d)(\d\d\d)(\d\d)(\d\d)/) {
      my $y = $_[0]+2000;
      my ($y, $m, $d) = Date_NthDayOfYear($y,$2);
      ($y, $m, $d, $3, $4);
    }
  } else {
    print "Pattern is not yet implemented: $_[1]\n";
    "";
  }
}
