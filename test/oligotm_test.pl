# Copyright (c) 2007
# Triinu Koressaar and Maido Remm
# All rights reserved.
# 
#     This file is part of the oligotm library.
# 
#     The oligotm library is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.
# 
#     The oligotm library is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with the oligtm library (file gpl.txt in the source
#     distribution); if not, write to the Free Software
#     Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

use constant EPSILON => 1e-5;
$nr=1;
$failure=0;

die "Cannot execute ../src/oligotm" unless -x '../src/oligotm';

print STDERR "Tests nr $nr-";

open F, "./oligotm.txt" or die "Cannot open oligotm.txt\n";
while(<F>){
    chomp;
    next if $. == 1;
    @tmp=split(/\t/);
    $cmd="../src/oligotm -tp $tmp[1] -sc $tmp[2] -mv $tmp[3] -dv $tmp[4] -n $tmp[5] $tmp[0] |";
    open(CMD, $cmd); # execute the command and take the output
    $tm=<CMD>;
    chomp $tm;
    close(CMD);
    if($tm){
	if(($tm-$tmp[6])>EPSILON) {
	    $failure++;
	    print STDERR "$cmd FAILED (expected $tm, got $tmp[6])\n";
	}
    } else {
	$failure++;
	print STDERR  "$cmd FAILED (no output)\n";
    }
    $nr++;
}
close F;

$nr--;
if(!$failure){
    print STDERR "$nr: OK\n";
}
else{
    print STDERR "$nr: $failure FAILURES\n";
}
