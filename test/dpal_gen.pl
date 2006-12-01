#!/usr/local/bin/perl5 -w

$exe = shift @ARGV;

$ENV{TC_SILENT} = '1';

while (<>) {
    chop;
    $ENV{TC_COMMENT} = "$exe $_";
    system "$exe $_";
}
