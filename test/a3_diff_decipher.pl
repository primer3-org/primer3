#!/usr/bin/perl

# Regression test driver for the amplicon3_core executable.
#
# For usage, see the usage statement in the code, below.
#
# ======================================================================
# (c) Copyright 2021
# Whitehead Institute for Biomedical Research, Steve Rozen, 
# Andreas Untergasser and Helen Skaletsky
# All rights reserved.
# 
#   This file is part of the primer3 suite.
#
#   The primer3 suite is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License as
#   published by the Free Software Foundation; either version 2 of the
#   License, or (at your option) any later version.
#
#   This software is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this file (file gpl-2.0.txt in the source distribution); if
#   not, write to the Free Software Foundation, Inc., 51 Franklin St,
#   Fifth Floor, Boston, MA 02110-1301 USA
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# ======================================================================

use warnings 'all';
use strict;
use Cwd;
use Getopt::Long;
use POSIX; # For testing wait / system return value, e.g. WIFSIGNALED, WTERMSIG.....
use Config; # To get the number for SIGINT

sub perldiff($$);
sub perldifffloat($$);
sub main();

# Call system() with warnings turned off; needed for ActiveState / MS Windows.
sub _nowarn_system($); 

our $set_files = 'amplicon3/';

our %signo;


# Global variable for Errors
my $all_ok;

main();

sub main() {
    my %args;

    # select STDERR;
    $| = 1;
    my $exit_stat = 0;

    print "$0 compare with DECIPHER output\n";
    
    $all_ok = 1;
    my $overall_start_time = time();

    if (defined $Config{sig_name}) {
          my $i = 0;
          for my $name (split(' ', $Config{sig_name})) {
            $signo{$name} = $i;
            $i++;
          }
    }

    my $amplicons = $set_files . 'amplicons.csv';

    open(AMPLISEQ, '<', $amplicons) or die "Cannot read $amplicons";
    while(<AMPLISEQ>) {
        my @cur_amplicon = split(";", $_);
        my $r;  # Return value for tests

        my $amplicon3_output = $set_files . 'heli_ts_' . $cur_amplicon[0] . '_output';
        my $decipher_output = $set_files . 'decipher_' . $cur_amplicon[0] . '_heli.csv';

        $r = perldifffloat($amplicon3_output, $decipher_output);

        print ($cur_amplicon[0] . ": Helix values ");
        if ($r == 0) {
            print "[OK]\n";
        } else {
            $all_ok = 0;
            print "[FAILED]\n";
            $exit_stat = -1;
        }
        $r = 0;

        $amplicon3_output = $set_files . 'only_mv_ts_' . $cur_amplicon[0] . '_output';
        open(AMPOUT, '<', $amplicon3_output) or die "open $amplicon3_output: $!";
        my @out_arr = <AMPOUT>;
        my @temps;
        my @melt;
        my @deriv;
        for my $i (0 .. $#out_arr) {
            my @line = split("=", $out_arr[$i]);
            $line[1] =~ s/\n//g;
            if ($line[0] eq "AMPLICON_TEMPERATURES") {
                @temps = split(",", $line[1]);
            }
            if ($line[0] eq "AMPLICON_MELT_CURVE") {
                @melt = split(",", $line[1]);
            }
            if ($line[0] eq "AMPLICON_DERIVATIVE_CURVE") {
                @deriv = split(",", $line[1]);
            }
        }
        my $sum_melt = "";
        my $sum_deriv = "";
        for my $k (0 .. $#temps) {
            $sum_melt .= $temps[$k];
            $sum_deriv .= $temps[$k];
            $sum_melt .= "\t";
            $sum_deriv .= "\t";
            $sum_melt .= $melt[$k];
            $sum_deriv .= $deriv[$k];
            $sum_melt .= "\n";
            $sum_deriv .= "\n";
        }
        close(AMPOUT);

        $decipher_output = $set_files . 'decipher_' . $cur_amplicon[0] . '_melt.csv';
        $r = perldiff $sum_melt, $decipher_output;
        print ($cur_amplicon[0] . ": Melting values ");
        if ($r == 0) {
            print "[OK]\n";
        } else {
            $all_ok = 0;
            print "[FAILED]\n";
            $exit_stat = -1;
        }
        $r = 0;
        $decipher_output = $set_files . 'decipher_' . $cur_amplicon[0] . '_deriv.csv';
        $r = perldiff $sum_deriv, $decipher_output;
        print ($cur_amplicon[0] . ": Derivate values ");
        if ($r == 0) {
            print "[OK]\n";
        } else {
            $all_ok = 0;
            print "[FAILED]\n";
            $exit_stat = -1;
        }
    }
    close(AMPLISEQ);

    print "\nDONE ";

    print $all_ok ? "Passed all tests - [OK]\n\n\n" : "At least one test failed - [FAILED]\n\n\n";

    if ($all_ok != 1 ) {
        $exit_stat = -1;
    }

    exit 0;  #  Generally we want the testing to continue.
}

# Usage: perldiff("string", "filename2")
# Return 0 if no differences are found,
# othewise return 1;
sub perldiff($$) {
    my ($f1, $f2) = @_;

    open(F2, '<', $f2) or die "open $f2: $!";

    my @f1 = split("\n", $f1);
    my @f2 = <F2>;

    my $test_count = 0;
    my $all_count = 0;

    # If different number of lines, return FAIL.
    if (@f1 != @f2) {
        print "\n   Different number of lines;\n";
        print "   Check file $f2\n\n";
        return 1;
    }
    # check for differences on the lines, themselves
    my $linenumber = 0;
    my $line_end_diff = 0;
    my $linErr = 0;
    # iterate using lines in file1
    while (@f1) {
        $linenumber++;

        # get the lines from each respective file
        my $l1 = shift @f1;
        my $l2 = shift @f2;
        my @arr1 = split("\t", $l1);
        my @arr2 = split("\t", $l2);
        if ($#arr1 != $#arr2) {
            $linErr = 1;
        }
        for my $i (0 .. $#arr1) {
            if ($arr1[$i] - $arr2[$i] > 0.000003) {
                $linErr = 1;
            }
            if (($arr1[$i] < 0.000003) || ($arr2[$i] < 0.000003)) {
               $test_count++;
            }
            $all_count++;
        }

        $linenumber++;
        # Check for difference between two edited lines (line by line)
        if ($linErr != 0) {
            print "\n   Difference found at line $linenumber!\n";
            print "   There may be additional differences.\n";
            print "   Check file $f2\n\n";
            return 1;
        }
    }
    close(F1);
    close(F2);
    my $percent = ($test_count / $all_count) * 100;
    print "\n   $test_count of $all_count (" . sprintf("%.1f", $percent) . "%) below test limit of 0.000003.\n";

    return 0;
}


sub perldifffloat($$) {
    my ($f1, $f2) = @_;

    open(F1, '<', $f1) or die "open $f1: $!";
    open(F2, '<', $f2) or die "open $f2: $!";

    my @f1 = <F1>;
    my @f2 = <F2>;

    shift @f1;
    shift @f1;
    shift @f2;

    my $test_count = 0;
    my $all_count = 0;

    # If different number of lines, return FAIL.
    if (@f1 != @f2) {
        print "\n   Different number of lines;\n";
        print "   on command line in test directory run:\n";
        print "   diff $f1 $f2\n\n";
        return 1;
    }
    # check for differences on the lines, themselves
    my $linenumber = 0;
    my $line_end_diff = 0;
    my $linErr = 0;
    # iterate using lines in file1
    while (@f1) {
        $linenumber++;

        # get the lines from each respective file
        my $l1 = shift @f1;
        my $l2 = shift @f2;
        my @arr1 = split("\t", $l1);
        my @arr2 = split("\t", $l2);
        if ($#arr1 != $#arr2) {
            $linErr = 1;
        }
        for my $i (0 .. $#arr1) {
            if ($arr1[$i] - $arr2[$i] > 0.000003) {
                $linErr = 1;
            }
            if (($arr1[$i] < 0.000003) || ($arr2[$i] < 0.000003)) {
               $test_count++;
            }
            $all_count++;
        }

        $linenumber++;
        # Check for difference between two edited lines (line by line)
        if ($linErr != 0) {
            print "\n   Difference found at line $linenumber!\n";
            print "   There may be additional differences.\n";
            print "   On command line in test directory run\n\n";
            print "   diff $f1 $f2\n\n";
            return 1;
        }
    }
    close(F1);
    close(F2);
    my $percent = ($test_count / $all_count) * 100;
    print "\n   $test_count of $all_count (" . sprintf("%.1f", $percent) . "%) below test limit of 0.000003.\n";

    return 0;
}
