# Regression test driver for the primer3_core executable.
#
# Usage: perl p3test.pl [<primer3>] [-w|--windows] [-v|--valgrind]
#
# <primer3> defaults to ../src/primer3_core
#
# If <primer3> is specified, the executable run is <primer3>.
#
# Stderr difference tests for fatal errors are performed only if
# <primer3> is '<any dir>/primer3_core' (because the executable
# name is part of the text written to stderr).

#use warnings 'all';
use strict;
use Cwd;

sub perldiff($$);
sub test_fatal_errors($$);
sub main();

our $def_executable = "../src/primer3_core";
our $exe = '../src/primer3_core';

main();

sub main() {

    $exe = $ARGV[0] if defined $ARGV[0];

    my ($winFlag, $valgrind_prefix) = (0, '');

    # catch ARGV for windows/valgrind commands...
    for my $arg (@ARGV) {
	if ($arg eq "--windows" || $arg eq "-w") {
	    # need to automatically get absolute path to windows executable
	    #my $winPath = getcwd();
	    # find where path ends with "something/something/test[/]"
	    #my $regex = "test[/]*\$";
	    # replace this with "/src/primer3_core.exe"
	    #$winPath =~ s/$regex/src\/primer3_core\.exe/g;
	    my $winPath = "..\\src\\primer3_core\.exe";
	    # shove this into the $exe variable
	    $exe = $winPath;
	    $def_executable = $winPath; # keep things happy @ line 237
	    # set a binary windows flag to 1 (used later - keeps from calling @ARGV again)
	    $winFlag = 1;
	}
	elsif ($arg eq "--valgrind" || $arg eq "-v") {
	    $valgrind_prefix = 
		"valgrind --leak-check=yes --show-reachable=yes --logfile=p3vg ";
	}
    }

    my $exit_stat = 0;

    die "Cannot execute $exe" unless -x $exe;

    print STDERR 
	"\n\n$0: testing $valgrind_prefix$exe\n\nSTART, ", scalar(localtime), "\n";

    test_fatal_errors($winFlag, $valgrind_prefix);

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
	    print STDERR "\nNOTE: this test takes _much_ longer than the others ";
	    print STDERR "(10 to 20 minutes or more).\n";
	    print STDERR "starting $test at ", scalar(localtime), "...";
	}
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
	    # generate the necc. files; If winFlag is 
	    # set, run command with absolute
	    # path to primer3_core.exe
	    if ($winFlag) {
	        $tmpCmd = "..\\$exe -strict_tags <../$input >../$tmp";
	    }  else {
		$tmpCmd = "../$exe -strict_tags <../$input >../$tmp";
	    }
	    $r = system $tmpCmd;
	    # back to main directory
	    chdir "../";
	} elsif ($test =~ /formatted$/) {
	    my $cmd = "$valgrind_prefix$exe -strict_tags -format_output <$input >$tmp";
	    $r = system $cmd;
	} else {
	    my $cmd = "$valgrind_prefix$exe -strict_tags <$input >$tmp";
	    $r = system $cmd;
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
    # check for line number differences - if so, return FAIL
    if (@f1 != @f2) {
        return 1;
    }
    # check for differences on the lines, themselves
    my $linenumber = 0;
    my $line_end_diff = 0;
    # iterate using lines in file1
    while (@f1) {
	# get the lines from each respective file
	my $l1 = shift @f1;
	my $l2 = shift @f2;

	# strip out directory goofiness containing executable name
	# (remove up to and including the colon from each file being compared)
	my $regex = "[^:]+:";
        $l1 =~ s/$regex//g;
        $l2 =~ s/$regex//g;
        $linenumber++;
	# check for differences in edited lines (line by line)
	if ($l1 ne $l2) {
	    return 1;
	}
    }
    return 0;
}

sub test_fatal_errors($$) {
    my ($winFlag, $valgrind_prefix) = @_;

    # Get all the files that match primer_global_err/*.in:
    my @inputs = glob("./primer_global_err/*.in");
    my $r;
    my $problem = 0;
    print STDERR "\ntesting fatal errors...";
    for (@inputs) {
        my ($root) = /(.*)\.in$/;  # Hint, the parens around $root give
                                   # the result of the match in
                                   # an array context.
	my $cmd = "$valgrind_prefix$exe <$_ > $root.tmp 2> $root.tmp2";
	if ($winFlag) {
	    system $cmd;
	}
	else {
	    system $cmd;
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
