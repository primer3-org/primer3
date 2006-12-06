#!/usr/local/bin/perl

# Script that re-implements the calculation for long sequences.
# Compare output to that of long_seq_tm_test.pl

my ($seq, $start, $len) = @ARGV;
my $s = substr($seq, $start, $len);
my $GC_count = $s =~ tr/gGcC/gGcC/;

$r = log(50/1000.0)/log(10);
$tm = 81.5 + 16.6*$r + 41*$GC_count/$len - 600./$len;
print "tm = $tm\n";
