# Regression test driver for the primer3_core executable.
#
# For usage, see the usage statement in the code, below.
#
use warnings 'all';
use strict;
use Cwd;
use Getopt::Long;

sub perldiff($$);
sub test_fatal_errors();
sub main();

# Call system() with warnings turned off; needed for ActiveState / MS Windows.
sub _nowarn_system($); 

our $def_executable = "../src/primer3_core";
our $exe = '../src/primer3_core';
our ($verbose, $do_valgrind, $winFlag);

# The following format works with valgrind-3.2.3:
our $valgrind_format  = "/usr/local/bin/valgrind "
		   . " --leak-check=yes --show-reachable=yes "
		   . " --log-file-exactly=%s.vg ";

main();

sub main() {
    my %args;

    # GetOptions handles various flag abbreviations and formats,
    # such as  -e ../src/primer3_core, --exe ../src/primer3_core, 
    # --exe=.../src/primer3_core
    if (!GetOptions(\%args,
		    'valgrind',
		    'windows',
		    'verbose',
		    'executable=s',
		    )) {
	print STDERR "Usage: $0 [--executable <primer3 executable>] [ --valgrind ] [  --verbose ] [--windows]\n";
	exit -1;
    }

    $exe = $args{'executable'} if defined$ args{'executable'};
    $winFlag = defined $args{'windows'};
    $verbose = defined $args{'verbose'};
    $do_valgrind = $args{'valgrind'};
    if ($winFlag && $do_valgrind) {
	print STDERR "$0: Cannot specify both --valgrind and --windows\n";
	exit -1;
    }

    if ($winFlag) {
	$exe = '..\\src\\primer3_core.exe';
	# $def_executable = $exe; # keep things happy @ line 237 Brant probably not necessary any more
                                  # Also, line numbers of particular statements are not very stable.
    }

    my $exit_stat = 0;

    die "Cannot execute $exe" unless -x $exe;

    print STDERR 
	"\n\n$0: testing $exe\n\n",
	"START, ", scalar(localtime), "\n";
    print STDERR "verbose mode\n" if $verbose;
    print STDERR "valgrind mode\n" if $do_valgrind;

    test_fatal_errors;

    for my $test (
		  'primer_boundary', # Put the quickest tests first.
		  'primer_internal',
		  'primer_boundary_formatted',
		  'primer_internal_formatted',
		  'primer_start_codon',
		  'primer_boundary1',
		  'primer_internal1',
		  'primer_task',
		  'primer_task_formatted',
		  'primer_boundary1_formatted',
		  'primer_internal1_formatted',
		  'primer_check',
		  'primer_must_use',
		  'primer_must_use_formatted',
		  'primer_syntax',
		  'primer_end_pathology',
		  'primer_num_best',
		  'primer_quality_boundary',
		  'primer_obj_fn',
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
		  'primer_tm_lc_masking',
		  'primer_tm_lc_masking_formatted',
		  # Put primer_lib_amb_codes last because it is slow
		  'primer_lib_amb_codes',
		  ) {
	print STDERR "$test...";
	if ($test eq 'primer_lib_amb_codes') {
	    print STDERR 
		"\nNOTE: this test takes _much_ longer than the others ",
		"(10 to 20 minutes or more).\n",
		"starting $test at ", scalar(localtime), "...";
	}
	my $valgrind_prefix
	    = $do_valgrind ? sprintf $valgrind_format, $test : '';

	my $testx = $test;
	$testx =~ s/_formatted$//;
	my $input = $testx . '_input';
	my $output = $test . '_output';
	my $tmp = $test . '.tmp';
	die "Cannot read $input"  unless -r $input;
	die "Cannot read $output"  unless -r $output;

	my $r; # Return value for tests
	if ($test eq 'primer' || $test eq 'primer1') {
	    my $list_tmp = $test.'_list_tmp';
	    # We need to chdir below because primer3 puts the 'list' files
	    # in the current working directory.  Therefore we adjust
	    # the TestCenter result directory.

	    # get a list of the files to remove (if any) in this directory 
	    my @tempList = glob('./' . $test.'_list_tmp/*');
	    # delete them
	    for my $t (@tempList) {
        	unlink $t;
	    }
	    # go down into primer|primer1 list_tmp directory
	    chdir $list_tmp;

	    my $tmpCmd;
	    # generate the necc. files; If $winFlag is 
	    # set, run command with Windows backslashes
	    # in path.
	    if ($winFlag) {
	        $tmpCmd = "..\\$exe -strict_tags <../$input >../$tmp";
	    }  else {
		$tmpCmd = "$valgrind_prefix ../$exe -strict_tags <../$input >../$tmp";
	    }

	    $r = _nowarn_system($tmpCmd);
	    # back to main directory
	    chdir "../";
	} elsif ($test =~ /formatted$/) {
	    my $cmd = "$valgrind_prefix$exe -strict_tags -format_output <$input >$tmp";
	    $r = _nowarn_system($cmd);
	} else {
	    my $cmd = "$valgrind_prefix$exe -strict_tags <$input >$tmp";
	    $r = _nowarn_system($cmd);
	}

	unless ($r == 0) {
	    print STDERR "NON-0 EXIT: $r\n";
	    $exit_stat = -1;
	    next;
	}

	$r = perldiff $output, $tmp;

	if ($r == 0) {
	    print STDERR "[OK]\n";
	} else {
	    print STDERR "[FAILED]\n";
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
		$r = perldiff $list_tmp."/".$t, $list_last."/".$t;
	    }
	    print STDERR $test. "_list_files ";
	    if ($r == 0) {
		print STDERR "[OK]\n";
	    } 
	    else { 
		print STDERR "[FAILED]\n";
		$exit_stat = -1;
	    }
	}
    }
    unlink("./core") if -e "./core";
    if ($do_valgrind) {
	# Assume this is Unix/Linux envrionment, so
	# we have grep.
	my $r = system "grep ERROR *.vg */*.vg | grep -v 'ERROR SUMMARY: 0 errors'";
	if (!$r) { # !$r because grep returns 0 if something is found,
	           # and if something is found, we have a problem.
	    $exit_stat = -1;
	}
	$r = system "grep 'definitely lost' *.vg */*.vg | grep -v '0 bytes'";
	if (!$r) {
	    $exit_stat = -1;
	}
	$r = system "grep 'possibly lost' *.vg */*.vg   | grep -v '0 bytes'";
	if (!$r) {
	    $exit_stat = -1;
	}
    }
    print STDERR "DONE ", scalar(localtime), "\n";
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
	print STDERR "Different number of lines\n";
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
		print STDERR "removing executable name from\n",
		"$l1_orig\n$l2_orig\n";
	    }
	}

	# Remove everything up to and including the first
	# colon, which removes the executable name.
        # This substitution must follow the USAGE
	# substitution, above
	my $regex = "[^:]+:";

	my $quoteexe = quotemeta($def_executable);

	if ($exe ne $def_executable && ($l1 =~ /$quoteexe/ || $l2  =~ /$quoteexe/)) {
	    $l1 =~ s/$regex//g;
	    $l2 =~ s/$regex//g;
	    if ($verbose) {
		print STDERR "removing <executable>: from\n",
		"$l1_orig\n$l2_orig\n";
	    }
	}

        $linenumber++;
	# check for differences in edited lines (line by line)
	if ($l1 ne $l2) {
	    print STDERR 
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
    print STDERR "\ntesting fatal errors...\n";
    for (@inputs) {
        my ($root) = /(.*)\.in$/;  # Hint, the parens around $root give
                                   # the result of the match in
                                   # an array context.
	print STDERR "  $root\n";
	my $valgrind_prefix
	    = $do_valgrind ? sprintf $valgrind_format, $root : '';

	my $cmd = "$valgrind_prefix$exe <$_ > $root.tmp 2> $root.tmp2";
	if ($winFlag) {
	    $r = _nowarn_system($cmd);  # FIX ME --- both branches are the same
	}
	else {
	    $r = _nowarn_system($cmd);
	}
	if ($? == 0) {
	    my $r = $? >> 8;
	    print STDERR
		"\nErroneous 0 exit status ($?) from command $cmd\n";
	    $problem = 1;
	}
	if (perldiff "$root.tmp", "$root.out") {
	    print STDERR
		"Difference found between $root.out and $root.tmp\nfrom $cmd\n\n";
	    $problem = 1;
	}
	if (perldiff "$root.tmp2", "$root.out2") {
	    print STDERR 
		"\nDifference found between $root.out2 and $root.tmp2\nfrom $cmd\n\n";
	    $problem = 1;
	}
    }
    print STDERR $problem ? "[FAILED]" : "[OK]" ,"\n";
}

sub _nowarn_system($) {
    my $cmd = shift;
    no warnings 'all';
    system $cmd;
}
