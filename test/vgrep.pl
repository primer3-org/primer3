#!/usr/local/bin/perl

system "grep ERROR             *.vg */*.vg | grep -v 'ERROR SUMMARY: 0 errors'";
system "grep 'definitely lost' *.vg */*.vg | grep -v ' 0 bytes'";
system "grep 'possibly lost'   *.vg */*.vg | grep -v ' 0 bytes'";

system "grep ERROR             *.valg      | grep -v 'ERROR SUMMARY: 0 errors'";
system "grep 'definitely lost' *.valg      | grep -v ' 0 bytes'";
system "grep 'possibly lost'   *.valg      | grep -v ' 0 bytes'";
