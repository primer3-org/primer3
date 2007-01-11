#!/usr/bin/perl -w

use constant EPSILON => 1e-5;
$nr=1;
$failure=0;

die "Cannot execute ../src/oligotm" unless -x '../src/oligotm';

print STDERR "Tests nr $nr-";

open F, "./oligotm.txt" or die "Can't open oligotm.txt\n";
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
	$failure++ if(($tm-$tmp[6])>EPSILON);
    }
    else{
	$failure++;
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
