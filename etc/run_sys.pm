# *****************  run_sys.pm ****************************************
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
package run_sys;
require 5.002;
require Exporter;

@run_sys::ISA = qw(Exporter);

@EXPORT = qw(run_sys);
use strict;

sub run_sys{
    my($sys,$noexit) = @_;
    my $rc = 0xffff & system($sys);
    if($rc == 0){
	print "$sys completed \n";
    }elsif($rc == 0xff00){
	print "ERROR: $rc Command $sys failed: $! ";
	exit unless($noexit);
    }elsif($rc > 0x80){
	$rc >>= 8; 
	print "STATUS: $sys returned non-zero exit status $rc\n";
        exit unless($noexit);
    }else{
	print "ERROR: $sys ran with ";
	if($rc & 0x80){
	    $rc &= ~0x80;
	    print "coredump from ";
	}
	print "signal $rc\n";
	exit -1 unless($noexit);
    }
    return 1;
}
1;
