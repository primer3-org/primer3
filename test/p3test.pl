
# Regression test driver for the primer3_core executable.
#
# For usage, see the usage statement in the code, below.
#
# ======================================================================
# (c) Copyright 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008,2010,
#  2011,2012
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

our $def_executable = "../src/primer3_core";
our $exe = '../src/primer3_core';
our $set_files = '../test/';
our ($verbose, $do_valgrind, $winFlag, $fastFlag, $onetest);

our %signo;

our $valgrind_format;
                   
# Global variable for Errors
my $all_ok;

main();

sub main() {
    my %args;

    # select STDERR;
    $| = 1;

    print "$0 for primer3_core version 2.3.6\n";
    
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
    # such as  -e ../src/primer3_core, --exe ../src/primer3_core, 
    # --exe=.../src/primer3_core
    if (!GetOptions(\%args,
                    'valgrind',
                    'windows',
                    'fast',
                    'verbose',
                    'onetest=s',
                    'executable=s',
                    )) {
        print "Usage: perl p3test.pl \\\n",
        "    [--executable <primer3 executable>] [ --onetest <test_name> ] ",
	"[ --valgrind ] [  --verbose ] [--windows] [--fast]\n",
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
        . " --leak-check=yes --show-reachable=yes "
        . "$log_file_arg_for_valgrind=%s.vg ";

    if ($winFlag) {
        $exe = '..\\src\\primer3_core.exe';
        $set_files = '..\\test\\';
    }

    my $exit_stat = 0;

    die "Cannot execute $exe" unless -x $exe;

    print 
        "\n\n$0: testing $exe\n\n",
        "START, ", scalar(localtime), "\n";
    print "verbose mode\n" if $verbose;
    print "valgrind mode\n" if $do_valgrind;

    # Tests in %default_version2 use --default_version=2.
    # Others use -default_version=1.
    my %default_version2 = ('th-w-other-tasks' => 1,
			    'primer_thermod_align' => 1,
			    'primer_thermod_align_formatted' => 1,
			    'primer1_th' => 1,
			    'primer_mispriming_th' => 1,
			    'primer_must_use_th' => 1,
			    'primer_new_tasks_th' => 1,
			    'primer_task_th' => 1,
	                    'primer_thal_args' => 1,
	
			    # 'test_compl_error' => 1,
			    # 'test_left_to_right_of_right' => 1,

			    'primer_boundary' => 1, # Put the quickest tests first.

			    # 'primer_boundary1' => 1,
			    # 'primer_boundary_formatted' => 1,
			    # 'primer_boundary1_formatted' => 1,
 
                            'primer3_v1_1_4_default_settings',
                            'primer3web_v0_4_0_default_settings',
   	                    'primer3web_v3_0_0_default_settings',
                            'primer3web_v4_0_0_default_settings',
			    
			    # 'primer_internal' => 1,
			    # 'primer_internal1' => 1,
			    # 'primer_internal_formatted' => 1,
			    # 'primer_internal1_formatted' => 1,
			    
			    # 'primer_ok_regions' => 1,
			    # 'primer_ok_regions_formatted' => 1,
			    # 'primer_ok_regions2' => 1,

			    # 'primer_tm_lc_masking' => 1,
			    # 'primer_tm_lc_masking_formatted' => 1,

			    # 'primer_start_codon' => 1,
			    
			    # 'primer_task' => 1,
			    # 'primer_task_formatted' => 1,
			    # 'primer_renewed_tasks' => 1,
			    # 'primer_new_tasks' => 1,
			    # 'primer_new_tasks_formatted' => 1,
			    
			    # 'primer_must_overlap_point' => 1,
			    # 'primer_overlap_junction' => 1,
			    
			    # 'primer_all_settingsfiles' => 1,
			    # 'primer_high_tm_load_set' => 1,
			    # 'primer_high_gc_load_set' => 1,

			    # 'primer_gc_end' => 1,
			    # 'primer_check' => 1,
			    # 'primer_must_use' => 1,
			    # 'primer_must_use_formatted' => 1,
			    # 'primer_syntax' => 1,
			    # 'primer_end_pathology' => 1,
			    # 'primer_num_best' => 1,
			    # 'primer_quality_boundary' => 1,
			    # 'primer' => 1,
			    # 'primer1' => 1,
			    # 'primer_mispriming' => 1,
			    # 'primer_mispriming_formatted' => 1,
			    # 'primer_mispriming_boundary1' => 1,
			    # 'primer_mispriming_boundary1_formatted' => 1,
			    # 'primer_mispriming_boundary2' => 1,
			    # 'primer_mispriming_boundary2_formatted' => 1,
			    # 'primer_mispriming_long_lib' => 1,
			    # 'primer_rat' => 1,
			    # 'primer_human' => 1,
			    # 'long_seq' => 1,
			    # 'primer_position_penalty' => 1,
			    # 'primer_position_penalty_formatted' => 1,
			    # 'p3-tmpl-mispriming' => 1,

			    # Put slow tests last

			    # 'p3_3_prime_n' => 1,
			    # 'primer_three_prime_distance' => 1,
			    
			    # 'primer_obj_fn' => 1,
			    
			    # 'primer_lib_amb_codes' => 1,
			    );

    my @TESTS = ( 
			# New tests that use new melting temperature
			# and thermodynamic alignments
			# Put the quickest tests first.
			'primer_must_use_th',
			'primer_task_th',
			'primer_thal_args',
	                'primer_thal_max_seq_error',
	                'primer_first_base_index',
                        'primer_must_match',

			'test_compl_error',
			'test_left_to_right_of_right',
                        'dv_conc_vs_dntp_conc',
			'primer_boundary',
			'primer_boundary1',
			'primer_boundary_formatted',
			'primer_boundary1_formatted',

                        'primer3_v1_1_4_default_settings',
                        'primer3web_v0_4_0_default_settings',
	                'primer3web_v3_0_0_default_settings',
                        'primer3web_v4_0_0_default_settings',

			'primer_internal',
			'primer_internal1',
			'primer_internal_formatted',
			'primer_internal1_formatted',

			'primer_ok_regions',
			'primer_ok_regions_formatted',
			'primer_ok_regions2',

			'primer_tm_lc_masking',
			'primer_tm_lc_masking_formatted',

			'primer_start_codon',

			'primer_task',
			'primer_task_formatted',
			'primer_renewed_tasks',
			'primer_new_tasks',
			'primer_new_tasks_formatted',

			'primer_must_overlap_point',
			'primer_overlap_junction',

			'primer_all_settingsfiles',
			'primer_high_tm_load_set',
			'primer_high_gc_load_set',

			'primer_gc_end',
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

			'primer_three_prime_distance',

			'primer_obj_fn',
			'p3_3_prime_n',

			# Put slow tests last
			'primer_mispriming_th',
			'th-w-other-tasks',
			'primer_new_tasks_th',
			'primer_thermod_align',             
			'primer_thermod_align_formatted',
			'primer1_th',
			'primer_lib_amb_codes',

		  );


    if (!$args{'onetest'}) {
	test_fatal_errors();
    } else {
	@TESTS = $args{'onetest'};
    }

    # The range of this for loop is a set of test names
    # that get translated into file names inside the loop.
    for my $test (@TESTS) {

        # We are inside the for loop here....
        print "$test...";
	my $test_start_time = time();

        if ($fastFlag && (($test eq 'th-w-other-tasks')
            || ($test eq 'primer_obj_fn')
            || ($test eq 'primer1_th')
            || ($test eq 'primer_mispriming_th')
            || ($test eq 'primer_new_tasks_th')
	        || ($test eq 'primer_thermod_align')
	        || ($test eq 'primer_thermod_align_formatted'))) {
            print "[skiped in fast mode]\n";
            next;
        }

	my $default_version;
	if ($default_version2{$test}) {
	    $default_version = '-default_version=2';
	} else {
	    $default_version = '-default_version=1';
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

	my %list_test 
	    = ('primer' => 1,
	       'primer1' => 1,
	       'primer1_th' => 1,
	       'th-w-other-tasks' => 1);

        if (exists($list_test{$test})) {
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
            if (!chdir $list_tmp) { die "chdir $list_tmp: $!\n" }

            my $tmpCmd;
            # generate the necc. files; If $winFlag is 
            # set, run command with Windows backslashes
            # in path.
	    my $settings = '';
	
            if ($winFlag) {
		if ($test eq 'th-w-other-tasks') {
		    $settings = '-p3_settings ..\\th-w-other-tasks-settings.txt -echo_settings';
		}
		$tmpCmd = "..\\$exe $default_version -strict_tags $settings ../$input >../$tmp";
            }  else {
		if ($test eq 'th-w-other-tasks') {
		    $settings = '-p3_settings ../th-w-other-tasks-settings.txt -echo_settings';
		}
                $tmpCmd = "$valgrind_prefix ../$exe $default_version -strict_tags $settings ../$input >../$tmp";
            }
            $r = _nowarn_system($tmpCmd);
            # back to main directory
            if (!chdir "../") { die "chdir \"..\": $!\n" }

        } elsif ($test =~ /settings$/) {
	    my $cmd = "$valgrind_prefix$exe $default_version -strict_tags -p3_settings_file=../$test.txt -echo <$input >$tmp";
            $r = _nowarn_system($cmd);

	} elsif ($test =~ /formatted$/) {
            my $cmd = "$valgrind_prefix$exe $default_version -strict_tags -format_output <$input >$tmp";
            $r = _nowarn_system($cmd);

        } elsif ($test =~ /_load_set/) {
            my $cmd = "$valgrind_prefix$exe $default_version -strict_tags -p3_settings_file=$set_files$test.set -echo <$input >$tmp";
            $r = _nowarn_system($cmd);

        } else {
            my $cmd = "$valgrind_prefix$exe $default_version -strict_tags <$input >$tmp";
            $r = _nowarn_system($cmd);
        }

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

        if (exists($list_test{$test})) {
	    # Special processing for tests that generate files containing
	    # primer lists.
            my $list_tmp = $test.'_list_tmp';
            my $list_last = $test.'_list_last';
	    if (!-r $list_last) {
		print "Configuration error, $list_last does not exist or is not readable ";
		$r = 1;
	    }
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
            } else {
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
	    $all_ok = 0;
	    print "\nValgrind ERRORs found [FAILED]";
            $exit_stat = -1;
        }
        $r = system "grep 'definitely lost' *.vg */*.vg | grep -v ' 0 bytes'";
        if (!$r) {
	    print "\nValgrind LEAKSs found [WARNING]";
            $exit_stat = -1;
        }
        $r = system "grep 'possibly lost' *.vg */*.vg   | grep -v ' 0 bytes'";
        if (!$r) {
	    print "\nValgrind LEAKSs found [WARNING]";
            $exit_stat = -1;
        }
    }
    print "\nTests in $0 ran for ", (time() - $overall_start_time), " s \(", ((time() - $overall_start_time)/60), " min\)\n";
    print "\nDONE ", scalar(localtime), " ";

    print $all_ok ? "Passed all tests - [OK]\n\n\n" : "At least one test failed - [FAILED]\n\n\n";

    exit 0 # $exit_stat; Change here, generally we want the testing to continue.
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

sub test_fatal_errors() {

    # Get all the files that match primer_global_err/*.in:
    my @inputs = glob("./primer_global_err/*.in");
    my $r;
    my $fatal_error_problem = 0;
    print "\ntesting fatal errors...\n";
    my $output_and_err_tested = 0;
    for (@inputs) {
        my ($root) = /(.*)\.in$/;  # Hint, the parens around $root give
                                   # the result of the match in
                                   # an array context.
        print "  $root\n";
        my $valgrind_prefix
            = $do_valgrind ? sprintf $valgrind_format, $root : '';

        my $cmd;
	if ($_ =~ /bad_settings\d\.in/) {
	    # For testing the settings files we need the
	    # names of tests and the settings to be parallel
	    $cmd = "$valgrind_prefix$exe -strict_tags -p3_settings_file $_ primer_global_err/input_for_settings_tests.txt  > $root.tmp 2> $root.tmp2";
	    # print STDERR $cmd, "\n";
	} else {
	    if ($output_and_err_tested) {
		$cmd = "$valgrind_prefix$exe -strict_tags <$_ > $root.tmp 2> $root.tmp2";
	    } else {
		# We need to test the -output and -error command line arguments plus
                # simply taking the file name as an argument (no "<")
		$cmd = "$valgrind_prefix$exe -strict_tags $_ -output $root.tmp -error $root.tmp2";
		$output_and_err_tested = 1;
		print "Testing -output and -error flags on\n$cmd\n";
	    }
	}

	$r = _nowarn_system($cmd);

        if ($? == 0) {
            my $r = $? >> 8;
            print
                "\nErroneous 0 exit status ($?) from command $cmd\n";
            $fatal_error_problem = 1;
        }
        if (perldiff "$root.tmp", "$root.out") {
            print
                "Difference found between $root.out and $root.tmp\nfrom $cmd\n\n";
            $fatal_error_problem = 1;
        }
        if (perldiff "$root.tmp2", "$root.out2") {
            print 
                "\nDifference found between $root.out2 and $root.tmp2\nfrom $cmd\n\n";
            $fatal_error_problem = 1;
        }
    }
    if ($fatal_error_problem == 1){
        $all_ok = 0;
    }
    print(($fatal_error_problem ? "[FAILED]" : "[OK]") ,"\n");
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
