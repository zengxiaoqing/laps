# *****************  fxa.pm ****************************************
# Copyright (C) 1998  James P. Edwards
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
#
package fxa;
sub Get_env'fxa{ #'
    local($fxa_env_file) = "/awips/fxa/.environs";
    open(FXA,"$fxa_env_file") || return 0;
    @fxaenv = <FXA>;
    close(FXA);
    my $host = `hostname`;
    $host =~ s/(^[^\.]+)\./$1/;
    $fxa_env_file.=".".$host;
    
    chomp($fxa_env_name);
    open(FXA,"$fxa_env_file");
    push(@fxaenv,<FXA>);
    close(FXA);
    my $i;
    foreach(@fxaenv){
        $i++;

        next if /^\#/;
        next unless(/\s*([^\s]+)\s+([^\s]+)\s*$/);
        $evar = $1;
        $eval = $2;

	$eval =~ s/\$$evar/$ENV{$evar}/g;
	$eval =~ s/\$\{$evar\}/$ENV{$evar}/g;

       $eval =~ s/\$\{/\$/g;
        $eval =~ s/\}//g;

        $ENV{$evar} = $eval;
    }
    close(FXA);
    foreach(keys %ENV){
	while($ENV{$_} =~ /\$([^\:\s\/]+)/){
	  my $sub = $1;
#          print "$_ $sub $ENV{$sub}\n";
	  $ENV{$_} =~ s/\$$sub/$ENV{$sub}/;
        }
    }
    return 1;
}

1;
sub Set_logdir'fxa{ #'
    
    if( -d $ENV{LOG_DIR}){
#       local($yymmdd) = `date -u +%y%m%d`;
#       chomp($yymmdd);
        local(@gmtime) = gmtime;
        $gmtime[4]++;
#  Modified to make logfile directory to  $LOG_DIR/yyyymmdd/laps LW 1/7/00
#       $gmtime[5] -= 100 if($gmtime[5]>99);
        $gmtime[5] += 1900;
#       for($i=3;$i<=5;$i++){
        for($i=3;$i<=4;$i++){
          $gmtime[$i]="0".$gmtime[$i] if(length($gmtime[$i])<2);
        }
        local($yymmdd) = $gmtime[5].$gmtime[4].$gmtime[3];
#        print "hera $yymmdd\n"; chomp($yymmdd); print $yymmdd;


        $LAPS_LOG_PATH= "$ENV{LOG_DIR}/$yymmdd";

        if(! -d $LAPS_LOG_PATH){
            mkdir $LAPS_LOG_PATH, 0777 || 
                die "Could not create Log directory $LAPS_LOG_PATH";
        }
        $LAPS_LOG_PATH .= "/laps";
        if(! -d $LAPS_LOG_PATH){
            mkdir $LAPS_LOG_PATH, 0777 || 
                die "Could not create Log directory $LAPS_LOG_PATH";
        }
    }else{
        die "Could not find log dir $ENV{LOG_DIR}";
    }
    return $LAPS_LOG_PATH;
}

1;

