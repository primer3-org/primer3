
# Command line test driver for the primer3_core executable.
#
# For usage, see the usage statement in the code, below.
#
# ======================================================================
# (c) Copyright 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008,
#  2010,2011,2012
# Whitehead Institute for Biomedical Research, Steve Rozen, 
# Andreas Untergasser and Helen Skaletsky
# All rights reserved.
# 
#   This file is part of the primer3 suite.
#
#   The primer3 suite is free software; you can
#   redistribute them and/or modify them under the terms of the GNU
#   General Public License as published by the Free Software Foundation;
#   either version 2 of the License, or (at your option) any later
#   version.
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
sub runtest($$$$);
sub main();

# Call system() with warnings turned off; needed for ActiveState / MS Windows.
sub _nowarn_system($); 

our $def_executable = "../src/primer3_core";
our $exe = '../src/primer3_core';
our $set_files = '../test/';
our ($verbose, $do_valgrind, $winFlag, $fastFlag);

our %signo;

our $valgrind_format;
                   
# Global variable for Errors
my $all_ok;

main();

sub main() {
    my %args;

    # select STDERR;
    $| = 1;

    $all_ok = 1;

    if (defined $Config{sig_name}) {
          my $i = 0;
          for my $name (split(' ', $Config{sig_name})) {
            $signo{$name} = $i;
            $i++;
          }
    }

    print "\n\nTESTING command line arguments\n\n";

    # GetOptions handles various flag abbreviations and formats,
    # such as  -e ../src/primer3_core, --exe ../src/primer3_core, 
    # --exe=.../src/primer3_core
    if (!GetOptions(\%args,
                    'valgrind',
                    'windows',
                    'fast',
                    'verbose',
                    'executable=s',
                    )) {
        print "Usage: perl cmdline_test.pl \\\n",
        "    [--executable <primer3 executable>] [ --valgrind ] [  --verbose ] [--windows] [--fast]\n",
        "\n",
        "    where <primer3 executable> defaults to ../src/primer3_core\n";
        exit -1;
    }

    $exe = $args{'executable'} if defined$ args{'executable'};
    $winFlag = defined $args{'windows'};
    $verbose = defined $args{'verbose'};
    $fastFlag = defined $args{'fast'};
    $do_valgrind = $args{'valgrind'};
    if ($winFlag && $do_valgrind) {
        print "$0: Cannot specify both --valgrind and --windows\n";
        exit -1;
    }

    my $valgrind_exe = "/usr/local/bin/valgrind";
    my $log_file_arg_for_valgrind = "--log-file-exactly";
    if ($do_valgrind) {
	if (!-x $valgrind_exe) { 
	    warn "Cannot find $valgrind_exe; will try `which valgrind`\n";
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
        . " --leak-check=yes --show-reachable=yes "
        . "$log_file_arg_for_valgrind=%s.vg ";

    if ($winFlag) {
        $exe = '..\\src\\primer3_core.exe';
        $set_files = '..\\test\\';
    }

    my $exit_stat = 0;

    die "Cannot execute $exe" unless -x $exe;
    
    my $valgrind_prefix;
    my $cmd;
    
    # Use primer_check_input and primer_check_output files
    my $input = 'primer_check' . '_input';
    my $output = 'primer_check' . '_output';
    
    # Make sure that needed files are present and readable
    die "Cannot read $input"  unless -r $input;
    die "Cannot read $output" unless -r $output;
    
    # Test 1 : check '-output' and 'input_file'
    $valgrind_prefix
       = $do_valgrind ? sprintf $valgrind_format, 'cmd_test1' : '';
    $cmd = "$valgrind_prefix$exe -default_version=1 --strict -o cmd_test1.tmp $input";
    $exit_stat = runtest(1, $cmd, $output, 0);
    
    # Test 2 : check '> output' and '< input_file'
    $valgrind_prefix
       = $do_valgrind ? sprintf $valgrind_format, 'cmd_test2' : '';
    $cmd = "$valgrind_prefix$exe -default_version=1 -st > cmd_test2.tmp < $input";
    $exit_stat = runtest(2, $cmd, $output, 0);
    
    # Test 3 : check incorrect flag and '-error_file'
    $valgrind_prefix
       = $do_valgrind ? sprintf $valgrind_format, 'cmd_test3' : '';
    $cmd = "$valgrind_prefix$exe -flag -err=cmd_test3.tmp < $input";
    $exit_stat = runtest(3, $cmd, 'cmd_test3_output', 255);
    
    # Test 4 : check nonexistent input file and '2> error_file'
    $valgrind_prefix
       = $do_valgrind ? sprintf $valgrind_format, 'cmd_test4' : '';
    $cmd = "$valgrind_prefix$exe invalid_input 2> cmd_test4.tmp";
    $exit_stat = runtest(4, $cmd, 'cmd_test4_output', 255);

    # Test 5 : check that io_version=3 fails
    $valgrind_prefix
       = $do_valgrind ? sprintf $valgrind_format, 'cmd_test5' : '';
    $cmd = "$valgrind_prefix$exe -io_version=3 2> cmd_test5.tmp";
    $exit_stat = runtest(5, $cmd, 'cmd_test5_output', 255);

    # ================================================== 
    # If we were running under valgrind to look for memory-related
    # errors (reading uninitialized memory, writing off the end of
    # an array, etc) or leaks, then we look through the valgrind
    # logs to summarize errors and leaks.
    if ($do_valgrind) {
        # Assume this is Unix/Linux environment, so
        # we have grep.
        my $r = system "grep ERROR *.vg | grep -v 'ERROR SUMMARY: 0 errors'";
        if (!$r) { # !$r because grep returns 0 if something is found,
            # and if something is found, we have a problem.
            $exit_stat = -1;
        }
        $r = system "grep 'definitely lost' *.vg | grep -v ' 0 bytes'";
        if (!$r) {
            $exit_stat = -1;
        }
        $r = system "grep 'possibly lost' *.vg   | grep -v ' 0 bytes'";
        if (!$r) {
            $exit_stat = -1;
        }
    }

    print $all_ok ? "Passed all tests - [OK]\n\n\n" : "At least one test failed - [FAILED]\n\n";

    exit $exit_stat;
}

# Usage: runtest(test_number, command, output_file, expected_exit_code)
sub runtest($$$$) {
    my $r;
    my $exit_stat = 0;
    my ($i, $cmd, $output, $exit) = @_;
    
    print "Test $i... ";
    
    $r = _nowarn_system($cmd);
    if ($winFlag) {
	$r = $r >> 8;
    }
        
    unless ($r == $exit) {
        $all_ok = 0;
        print "incorrect exit code [FAILED]\n";
        $exit_stat = -1;
    } else {
        $r = perldiff $output, "cmd_test$i.tmp";
    	if ($r == 0) {
    	   print "[OK]\n";
    	} else {
    	   $all_ok = 0;
    	   print "[FAILED]\n";
    	   $exit_stat = -1;
        }
    }
    return $exit_stat;
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
        print "Different number of lines ($f1 versus $f2)\n";
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
            if ($l1 =~ /(libprimer3|primer3) release \d+\.\d+\.\d+/
                && $l2 =~ /(libprimer3|primer3) release \d+\.\d+\.\d+/) {
                $l1 =~ s/(libprimer3|primer3) release \d+\.\d+\.\d+//;
                $l2 =~ s/(libprimer3|primer3) release \d+\.\d+\.\d+//;
            }
        }

	# If this is the tag with the settings file path, replace \ by / to make it
	# the same on both Linux and Windows
	if ($l1 =~ /^P3_SETTINGS_FILE_USED/ && $l2 =~ /^P3_SETTINGS_FILE_USED/) {
	    $l1 =~ s/\\/\//g;
	    $l2 =~ s/\\/\//g;
	}

        $linenumber++;
        # Check for difference between two edited lines (line by line)
        if ($l1 ne $l2) {
            print 
                "Difference found at line $linenumber:\n<  $l1_orig\n>  $l2_orig\n";
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
