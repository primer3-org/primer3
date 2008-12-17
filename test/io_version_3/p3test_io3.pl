
# Regression test driver for the primer3_core executable
# FOR -io_verions=3 THIS IS A TEMPORARY TEST DRIVER
#
# For usage, see the usage statement in the code, below.
#
# ======================================================================
# (c) Copyright 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008
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
sub test_fatal_errors();
sub main();

# Call system() with warnings turned off; needed for ActiveState / MS Windows.
sub _nowarn_system($); 

our $def_executable = "../src/primer3_core";
our $exe = "../../src/primer3_core";
our ($verbose, $do_valgrind, $winFlag, $fastFlag);

our %signo;

our $valgrind_format;
                   
# Global variable for Errors
my $all_ok;

main();

sub main() {
    my %args;

    select STDERR; $| = 1;

    $all_ok = 1;
    my $start_time;
    $start_time = time();

    if (defined $Config{sig_name}) {
          my $i = 0;
          for my $name (split(' ', $Config{sig_name})) {
            $signo{$name} = $i;
            $i++;
          }
    }

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
        print "Usage: perl p3test.pl \\\n",
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
	if ($valgrind_version =~ /3\.3\./) {
	    $log_file_arg_for_valgrind = "--log-file";
	}
    }

    $valgrind_format  = $valgrind_exe
        . " --leak-check=yes --show-reachable=yes "
        . "$log_file_arg_for_valgrind=%s.vg ";

    if ($winFlag) {
        $exe = '..\\..\\src\\primer3_core.exe';
    }

    my $exit_stat = 0;

    die "Cannot execute $exe" unless -x $exe;

    print 
        "\n\n$0: testing $exe\n\n",
        "START, ", scalar(localtime), "\n";
    print "verbose mode\n" if $verbose;
    print "valgrind mode\n" if $do_valgrind;

    test_fatal_errors();

    # The range of this for loop is a set of test names
    # that get translated into file names inside the loop.
    for my $test (
                  'primer_boundary', # Put the quickest tests first.
                  'primer_boundary_formatted',
                  'primer_boundary1',
                  'primer_boundary1_formatted',

                  'primer_internal',
                  'primer_internal1',
                  'primer_internal_formatted',
                  'primer_internal1_formatted',

                  'primer_tm_lc_masking',
                  'primer_tm_lc_masking_formatted',

                  'primer_start_codon',

                  'primer_task',
                  'primer_task_formatted',

                  'primer_check',
                  'primer_must_use',
                  'primer_must_use_formatted',

                  'primer_syntax',

                  'primer_end_pathology',
                  'primer_num_best',
                  'primer_quality_boundary',

                  'primer',
                  'primer1',

                  'primer_mispriming',
                  'primer_mispriming_formatted',
                  'primer_mispriming_boundary1',
                  'primer_mispriming_boundary1_formatted',
                  'primer_mispriming_boundary2',
                  'primer_mispriming_boundary2_formatted',
                  'primer_mispriming_long_lib',

                  'primer_rat',
                  'primer_human',
                  'long_seq',
                  'primer_position_penalty',
                  'primer_position_penalty_formatted',
                  'p3-tmpl-mispriming',
                  # Put slow tests last
                  'primer_obj_fn',
                  'primer_lib_amb_codes',
                  ) {

        # We are inside the for loop here....
        print "$test...";

        if ($fastFlag && (($test eq 'p3_3_prime_n')
                || ($test eq 'primer_obj_fn'))) {
            print "[skiped in fast mode]\n";
            next;
        }

        if ($test eq 'primer_lib_amb_codes') {
            if ($fastFlag) {
                print "[skiped in fast mode]\n";       
                next;
            }
            print 
                "\nNOTE: this test takes _much_ longer than the others ",
                "(10 to 20 minutes or more).\n",
                "starting $test at ", scalar(localtime), "...";
        }
        my $valgrind_prefix
            = $do_valgrind ? sprintf $valgrind_format, $test : '';

        # Figure out what the files are called for a particular test
        my $testx = $test;
        $testx =~ s/_formatted$//;
        my $input = $testx . '_input';
        my $output = $test . '_output';
        my $tmp = $test . '.tmp';

        # Make sure that needed files are present and readable
        die "Cannot read $input"  unless -r $input;
        die "Cannot read $output"  unless -r $output;

        my $r;                  # Return value for tests

        if (0 && ($test eq 'primer' || $test eq 'primer1')) { # Skipping these tests
            # These tests generate primer lists, which
            # need to be checked separately.

            my $list_tmp = $test.'_list_tmp';
            # We need to chdir below because primer3 puts the 'list' files
            # in the current working directory. 

            # get a list of the files to remove (if any) in this directory 
            my @tempList = glob('./' . $test.'_list_tmp/*');
            # delete them
            for my $t (@tempList) {
                unlink $t;
            }
            # go down into primer|primer1 list_tmp directory
            chdir $list_tmp;

            my $tmpCmd;

	    $tmpCmd = "$valgrind_prefix ../$exe -strict_tags -io_version=3 <../$input >../$tmp";
            $r = _nowarn_system($tmpCmd);
            # back to main directory
            chdir "../";
        } elsif ($test =~ /formatted$/) {
            my $cmd = "$valgrind_prefix$exe -strict_tags -io_version=3 -format_output <$input >$tmp";
            $r = _nowarn_system($cmd);
        } else {
            my $cmd = "$valgrind_prefix$exe -strict_tags -io_version=3 <$input >$tmp";
            $r = _nowarn_system($cmd);
        }

        unless ($r == 0) {
            $all_ok = 0;
            print "NON-0 EXIT: $r [FAILED]\n";
            $exit_stat = -1;
            next;
        }

        $r = perldiff $output, $tmp;
        if ($r == 0) {
            print "[OK]\n";
        } else {
            $all_ok = 0;
            print "[FAILED]\n";
            $exit_stat = -1;
        }

        if ($test eq 'primer' || $test eq 'primer1') {
            my $list_tmp = $test.'_list_tmp';
            my $list_last = $test.'_list_last';
            my @tempList = glob($list_last.'/*');
            # do _file by file_ comparison within primer_list_last to primer_list_tmp
            # since we are not using diff - diff used to do directory comparison
            for my $t (@tempList) {
                # sneakily get correct paths since glob of primer_list_last returns
                # primer_list_last/filename_within_primer_list_last
                my $regex = "[^/]+/";
                $t=~ s/$regex//g;
                if (perldiff $list_tmp."/".$t, $list_last."/".$t) {
                    $r = 1;
                }
            }
            print $test. "_list_files...";
            if ($r == 0) {
                print "[OK]\n";
            } 
            else {
                $all_ok = 0;
                print "[FAILED]\n";
                $exit_stat = -1;
            }
        }
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
            $exit_stat = -1;
        }
        $r = system "grep 'definitely lost' *.vg */*.vg | grep -v ' 0 bytes'";
        if (!$r) {
            $exit_stat = -1;
        }
        $r = system "grep 'possibly lost' *.vg */*.vg   | grep -v ' 0 bytes'";
        if (!$r) {
            $exit_stat = -1;
        }
    }
    print "Tests ran for ", (time() - $start_time), " seconds.\n";
    print "\n\nDONE ", scalar(localtime), " ";

    print $all_ok ? "Passed all tests - [OK]\n\n\n" : "At least one test failed - [FAILED]\n\n\n";

    exit $exit_stat;
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
        print "Different number of lines\n";
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
            if ($l1 =~ /primer3 release \d+\.\d+\.\d+/
                && $l2 =~ /primer3 release \d+\.\d+\.\d+/) {
                $l1 =~ s/primer3 release \d+\.\d+\.\d+//;
                $l2 =~ s/primer3 release \d+\.\d+\.\d+//;
            }
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

sub test_fatal_errors() {

    # Get all the files that match primer_global_err/*.in:
    my @inputs = glob("./primer_global_err/*.in");
    my $r;
    my $problem = 0;
    print "\ntesting fatal errors...\n";
    for (@inputs) {
        my ($root) = /(.*)\.in$/;  # Hint, the parens around $root give
                                   # the result of the match in
                                   # an array context.
        print "  $root\n";
        my $valgrind_prefix
            = $do_valgrind ? sprintf $valgrind_format, $root : '';

        my $cmd = "$valgrind_prefix$exe -io_version=3 <$_ > $root.tmp 2> $root.tmp2";

	$r = _nowarn_system($cmd);

        if ($? == 0) {
            my $r = $? >> 8;
            print
                "\nErroneous 0 exit status ($?) from command $cmd\n";
            $problem = 1;
        }
        if (perldiff "$root.tmp", "$root.out") {
            print
                "Difference found between $root.out and $root.tmp\nfrom $cmd\n\n";
            $problem = 1;
        }
        if (perldiff "$root.tmp2", "$root.out2") {
            print 
                "\nDifference found between $root.out2 and $root.tmp2\nfrom $cmd\n\n";
            $problem = 1;
        }
    }
    if ($problem == 1){
        $all_ok = 0;
    }
    print $problem ? "[FAILED]" : "[OK]" ,"\n";
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
