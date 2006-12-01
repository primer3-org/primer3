#!/usr/local/bin/perl

# Script for checking output of executable src/long_seq_tm_test
# To make long_seq_tm_test go to ../src, and do
# make long_seq_tm_test

my ($seq, $start, $len) = @ARGV;
my $s = substr($seq, $start, $len);
my $GC_count = $s =~ tr/gGcC/gGcC/;

$r = log(50/1000.0)/log(10);
$tm = 81.5 + 16.6*$r + 41*$GC_count/$len - 600./$len;
print "tm = $tm\n";
