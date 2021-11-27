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
sub test_fatal_errors();
sub main();

# Call system() with warnings turned off; needed for ActiveState / MS Windows.
sub _nowarn_system($); 

our $def_executable = "../src/amplicon3_core";
our $exe = '../src/amplicon3_core';
our $set_files = '../test/amplicon3/';
our ($verbose, $do_valgrind, $do_valgrinda, $do_valgrindb, $return_action,
     $winFlag, $onetest);

our %signo;

our $valgrind_format;
                   
# Global variable for Errors
my $all_ok;

main();

sub main() {
    my %args;

    # select STDERR;
    $| = 1;

    print "$0 for amplicon3_core version 2.5.0\n";
    
    $all_ok = 1;
    my $overall_start_time = time();

    if (defined $Config{sig_name}) {
          my $i = 0;
          for my $name (split(' ', $Config{sig_name})) {
            $signo{$name} = $i;
            $i++;
          }
    }

    # GetOptions handles various flag abbreviations and formats,
    # such as  -e ../src/amplicon3_core, --exe ../src/amplicon3_core,
    # --exe=.../src/amplicon3_core
    if (!GetOptions(\%args,
                    'valgrind',
                    'valgrinda',
                    'valgrindb',
                    'action',
                    'verbose',
                    'windows',
                    'executable=s',
                    )) {
        print "Usage: perl p3test.pl \\\n",
        "    [ --executable <amplicon3 executable> ] [ --valgrind ] [ --windows ]\n",
        "\n",
        "    where <amplicon3 executable> defaults to ../src/amplicon3_core\n";
        exit -1;
    }

    $exe = $args{'executable'} if defined$ args{'executable'};
    $winFlag = defined $args{'windows'};
    $verbose = defined $args{'verbose'};
    $do_valgrind = $args{'valgrind'};
    $do_valgrinda = $args{'valgrinda'};
    $do_valgrindb = $args{'valgrindb'};
    $return_action = $args{'action'};
    if ($do_valgrinda) {
        $do_valgrind = 1;
    }

    if ($do_valgrindb) {
        exit 0;
    }

    if ($winFlag && $do_valgrind) {
        print "$0: Cannot specify both --valgrind and --windows\n";
        exit -1;
    }

    my $valgrind_exe = "/usr/bin/valgrind";
    my $log_file_arg_for_valgrind = "--log-file-exactly";
    if ($do_valgrind) {
        if (!-x $valgrind_exe) {
            print "Cannot find $valgrind_exe; will try `/usr/local/bin/valgrind`\n";
            $valgrind_exe = "/usr/local/bin/valgrind";
        }
    if (!-x $valgrind_exe) {
        print "Cannot find $valgrind_exe; will try `which valgrind`\n";
        $valgrind_exe= `which valgrind`;
        chomp($valgrind_exe);
        if (!$valgrind_exe || ! -x $valgrind_exe) {
        die "Cannot execute $valgrind_exe";
        }
    }
    # Need to deal with different arguments in
    # different version of valgrind
    my $valgrind_version = `$valgrind_exe --version`;
    if ($valgrind_version =~ /([0-9]+)\.([0-9]+)\./) {
        if (($1 > 3) || (($1 == 3) && ($2 > 2))) {
        $log_file_arg_for_valgrind = "--log-file";
        }
    }
    }

    $valgrind_format  = $valgrind_exe
        . " --leak-check=yes --show-reachable=yes "  # --track-origins=yes
        . "$log_file_arg_for_valgrind=%s.vg ";

    if ($winFlag) {
        $exe = '..\\src\\amplicon3_core.exe';
        $set_files = '..\\test\\amplicon3\\';
    }

    my $exit_stat = 0;

    die "Cannot execute $exe" unless -x $exe;

    print 
        "\n\n$0: testing $exe\n\n",
        "START, ", scalar(localtime), "\n";
    print "valgrind mode\n" if $do_valgrind;

    my $amplicons = $set_files . 'amplicons.csv';

    my @TESTS = ( 
                  'plain',
                  'plain_on_command_line',
                  'classic',
                  'heli',
                  'long',
                  'optimal',
                  'only_mv',
                  'no_dmso'
                );

    # The range of this for loop is a set of test names
    # that get translated into file names inside the loop.
    for my $test (@TESTS) {

        open(AMPLISEQ, '<', $amplicons) or die "Cannot read $amplicons";
        while(<AMPLISEQ>) {
            my @cur_amplicon = split(";", $_);
            if ($#cur_amplicon < 1) {
                print(" No Sequence found!\n");
                next;
            }
            print ($test . "_" . $cur_amplicon[0] . "...");
            my $test_start_time = time();

            my $valgrind_prefix
                = $do_valgrind ? sprintf $valgrind_format, $test . '_ts_'. $cur_amplicon[0] : '';

            # Figure out what the files are called for a particular test
            my $testx = $test;
            $testx =~ s/_formatted$//;
            my $input = $set_files . $testx . '_input';
            open(READCOM, '<', $input) or die "Cannot read $input";
            my $com_line = <READCOM>;
            close(READCOM);
            chop($com_line);
            my $output = $set_files . $test  . '_ts_'. $cur_amplicon[0] . '_output';
            my $tmp = $set_files . $test . '_ts_'. $cur_amplicon[0] . '.tmp';

            # Make sure that needed files are present and readable
            # die "Cannot read $output"  unless -r $output;

            my $r;  # Return value for tests

            chop($cur_amplicon[1]);
            my $cmd = "$valgrind_prefix$exe $com_line $cur_amplicon[1] > $tmp";
            $r = _nowarn_system($cmd);

            unless ($r == 0) {
                $all_ok = 0;
                print "NON-0 EXIT: $r [FAILED]\n";
                $exit_stat = -1;
                next;
            }

            $r = perldiff $output, $tmp;
            print ((time() - $test_start_time), "s ");

            if ($r == 0) {
                print "[OK]\n";
            } else {
                $all_ok = 0;
                print "[FAILED]\n";
                $exit_stat = -1;
            }
        }
        close(AMPLISEQ);
    }                 # End of long for loop, for my $test in (.....) 

    # ================================================== 
    # If we were running under valgrind to look for memory-related
    # errors (reading uninitialized memory, writing off the end of
    # an array, etc) or leaks, then we look through the valgrind
    # logs to summarize errors and leaks.
    if ($do_valgrind) {
        # Assume this is Unix/Linux envrionment, so
        # we have grep.
        my $r = system "grep ERROR *.vg */*.vg | grep -v 'ERROR SUMMARY: 0 errors'";
        if (!$r) { # !$r because grep returns 0 if something is found,
            # and if something is found, we have a problem.
            $all_ok = 0;
            print "\nValgrind ERRORs found [FAILED]";
            $exit_stat = -1;
        }
        $r = system "grep 'definitely lost' *.vg */*.vg | grep -v ' 0 bytes'";
        if (!$r) {
            $all_ok = 0;
            print "\nValgrind LEAKSs found [WARNING]";
            $exit_stat = -1;
        }
        $r = system "grep 'possibly lost' *.vg */*.vg   | grep -v ' 0 bytes'";
        if (!$r) {
            $all_ok = 0;
            print "\nValgrind LEAKSs found [WARNING]";
            $exit_stat = -1;
        }
    }
    print "\nTests in $0 ran for ", (time() - $overall_start_time), " s \(", ((time() - $overall_start_time)/60), " min\)\n";
    print "\nDONE ", scalar(localtime), " ";

    print $all_ok ? "Passed all tests - [OK]\n\n\n" : "At least one test failed - [FAILED]\n\n\n";

    if ($all_ok != 1 ) {
        $exit_stat = -1;
    }

    if ($return_action || $do_valgrinda) {
        exit $exit_stat;
    } else {
        exit 0;  #  Generally we want the testing to continue.
    }
}

# Usage: perldiff("filename1", "filename2")
# Return 0 if no differences are found,
# othewise return 1;
sub perldiff($$) {
    my ($f1, $f2) = @_;
    open F1, $f1 or die "open $f1: $!"; 
    open F2, $f2 or die "open $f2: $!";
    my @f1 = <F1>;
    my @f2 = <F2>;
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
    # iterate using lines in file1
    while (@f1) {
        $linenumber++;

        # get the lines from each respective file
        my $l1 = shift @f1;
        my $l2 = shift @f2;
        my $l1_orig = $l1;
        my $l2_orig = $l2;

        # Handle the diff in empty_1.{out2,tmp2} due to
        # different executable names.
        if ($exe ne $def_executable && $l1 =~/^USAGE:\s+\S+/) {
            $l1 =~ s/^USAGE:\s+\S+/USAGE: ... /i;
            $l2 =~ s/^USAGE:\s+\S+/USAGE: ... /i;
            if ($verbose) {
                print "removing executable name from\n",
                "$l1_orig\n$l2_orig\n";
            }
        }

        # Remove everything up to and including the first
        # colon, which removes the executable name.
        # This substitution must follow the USAGE
        # substitution, above
        my $regex = "[^:]+:";

        my $quoteexe = quotemeta($def_executable);

        # Edit executable name
        if ($exe ne $def_executable && ($l1 =~ /$quoteexe/ || $l2  =~ /$quoteexe/)) {
            $l1 =~ s/$regex//g;
            $l2 =~ s/$regex//g;
            if ($verbose) {
                print "removing <executable>: from\n",
                "$l1_orig\n$l2_orig\n";
            }
        }

        # Edit release number
        if ($l1 ne $l2) {
            if ($l1 =~ /(libamplicon3|amplicon3) release \d+\.\d+\.\d+/
                && $l2 =~ /(libamplicon3|amplicon3) release \d+\.\d+\.\d+/) {
                $l1 =~ s/(libamplicon3|amplicon3) release \d+\.\d+\.\d+//;
                $l2 =~ s/(libamplicon3|amplicon3) release \d+\.\d+\.\d+//;
            }
        }

        # If this is the tag with the settings file path, replace \ by / to make it
        # the same on both Linux and Windows
        if ($l1 =~ /^P3_SETTINGS_FILE_USED/ && $l2 =~ /^P3_SETTINGS_FILE_USED/) {
            $l1 =~ s/\\/\//g;
            $l2 =~ s/\\/\//g;
        }

        # Fix exponent on windows
        $l2 =~ s/e-0(\d\d)/e-$1/g;
        $l2 =~ s/e0(\d\d)/e$1/g;

        $linenumber++;
        # Check for difference between two edited lines (line by line)
        if ($l1 ne $l2) {
            print 
                "\n   Difference found at line $linenumber:\n",
                "   <  $l1_orig\n",
                "   >  $l2_orig\n";
            print "   There may be additional differences.\n";
            print "   On command line in test directory run\n\n";
            print "   diff $f1 $f2\n\n";
            return 1;
        }
    }
    return 0;
}

sub _nowarn_system($) {
    my $cmd = shift;
    if ($verbose) { print "\n$cmd\n" }
    no warnings 'all';
    my $r = system $cmd;
    my $r2 = $?;
    if (!$winFlag && WIFEXITED($r2)) {
        $r = WEXITSTATUS($r2);
        if (defined $signo{'INT'}) {
            if ($r == $signo{'INT'}) {
                print "\nCommand: $cmd\n";
                print "generated exit value $r\n";
                print "Presumably caused by catching SIGINT\n";
                print "Tests halted\n\n";
                print "\nWARNING: not all tests executed ... [FAILED]\n";
                if (!$all_ok) {
                    print "\n\nIn addition, at least 1 test FAILED ... [FAILED]\n\n";
                }
                exit;
            }
        }
    }
    if (!$winFlag && WIFSIGNALED($r2)) {
        my $r2 = WTERMSIG($r2);
        print "Exited with signal $r2\n";
    }
    $r;
}
