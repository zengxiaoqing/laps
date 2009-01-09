package namelist;
use strict;

#===============================================================================

sub read {

   my ($file) = @_;
   my ($line,$id,$value,@value,$nvalue,%var);

   open(NL,"$file") or die "Namelist not found:  $file\n";
   while(<NL>) {
      $line = substr($_,0,-1);       # Remove end-of-line c/r
      $line =~ s/\s+//g;             # Remove white space
      ($line) = split(/!|#/,$line);  # Remove comments
      $line = substr($line,0,-1) if (substr($line,-1) eq ',');
      ($id,$value) = split(/=/,$line);
      if ($id ne "") {
         $id =~ tr/a-z/A-Z/;         # Make sure $id is upper case
         $var{"$id"} = "$value";
      }
   }
   close(NL);

   return %var;
}
1;

#===============================================================================

package time;
use strict;

#===============================================================================

# &unix_to_time: Calculate month, day, year, hour, min, sec from unix time
# Arguments: unix time
# Returns: month, day, year, hour, min, sec

sub unix_to_time {
    my($utm) = @_;
    my($i, $j, $n, $l, $d, $m, $y);

    $n = int($utm/86400);

    $utm -= 86400*$n;
    $l    = $n + 2509157;
    $n    = int((4*$l)/146097);
    $l   -= int( (146097*$n + 3)/4);
    $i    = int( (4000*($l+1))/1461001);
    $l   += 31 - int((1461*$i)/4);
    $j    = int((80*$l)/2447);
    $d    = $l - int( (2447*$j)/80);
    $l    = int($j/11);
    $m    = $j + 2 - 12*$l;
    $y    = 100*($n-49) + $i + $l;

    my @answer = ($y, $m, $d, int($utm/3600), int(($utm%3600)/60), int($utm%60));

    @answer
}
1;

#===============================================================================

# &time_to_unix: Calculate unix time from month, day, year, hour, min, sec
# Arguments: month, day, year (4 digit), hour, min, sec
# Returns: unix time

sub time_to_unix {
   my($month, $day, $year, $hour, $minute, $second) = @_;
   my($b, $g, $d, $e, $f, $today);

   $b = int ( ($month - 14) / 12 );
   $g = $year + 4900 + $b;
   $b = $month - 2 - 12*$b;

   $d = int( (1461*($g-100))/4);
   $e = int( (367*$b)/12);
   $f = int( (3*int($g/100))/4);

   $today = $d + $e - $f + $day - 2432076;

   86400*($today-40587) + 3600*$hour + 60*$minute + $second;
}
1;

#===============================================================================

package oputil;
use strict;

#===============================================================================

sub purge {

   my ($runs,$purge,$time) = @_;
   my (@rundir,$nrundir,$n,$py,$pm,$pd,$ph,$runtime);

   opendir(RUNS,$runs);
   @rundir = sort(grep(/^\d{4}-\d{2}-\d{2}-\d{4}/, readdir(RUNS)));
   closedir(RUNS);
   $nrundir = @rundir;
   if ($purge < 0) {
      $purge = int(-$purge);
      print "    > Keep $purge model runs.\n";
      for ($n=0;$n<$nrundir-$purge;$n++) {
         print "      Removing $runs/$rundir[$n]\n";
         system("rm -rf $runs/$rundir[$n]");
      }
      return;
   } else {
      print "    > Keep model runs that are less than $purge hours old.\n";
      foreach (@rundir) {
         if (/(\d{4})-(\d{2})-(\d{2})-(\d{2})/) {
            $py = $1;
            $pm = $2;
            $pd = $3;
            $ph = $4;
            $runtime = &time::time_to_unix($pm,$pd,$py,$ph,0,0);
            if ($time-$runtime > $purge*3600) {
               print "      Removing $runs/$_\n";
               system("rm -rf $runs/$_");
            } else { return; }
         }
      }
   }
}
1;
