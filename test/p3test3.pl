sub perldiff($$);
sub runtst($);
sub main();

#our $exe = '../src/primer3_core';
our $exe = './p3core.pl';
	   our $input ;
	   our $tmp ;
	   our $output ;
main();

sub main() {

    $cmd = 'rm -f *.tmp';
    my $r = system $cmd;
    

		  runtst('primer_boundary'); # Put the quickest tests first.
#		  runtst('p3_3_prime_n');
		  runtst('primer_tm_lc_masking');
#		  runtst('primer_tm_lc_masking_formatted');
		  runtst('primer_internal');
#		  runtst('primer_boundary_formatted');
#		  runtst('primer_internal_formatted');
		  runtst('primer_start_codon');
		  runtst('primer_boundary1');
		  runtst('primer_internal1');
		  runtst('primer_task');
#		  runtst('primer_task_formatted');
#		  runtst('primer_boundary1_formatted');
#		  runtst('primer_internal1_formatted');
		  runtst('primer_check');
		  runtst('primer_must_use');
#		  runtst('primer_must_use_formatted');
#		  runtst('primer_syntax');
		  runtst('primer_end_pathology');
		  runtst('primer_num_best');
#		  runtst('primer_quality_boundary');
		  runtst('primer_obj_fn');
#		  runtst('primer');
		  runtst('primer1');
		  runtst('primer_mispriming');
#		  runtst('primer_mispriming_formatted');
		  runtst('primer_mispriming_boundary1');
#		  runtst('primer_mispriming_boundary1_formatted');
		  runtst('primer_mispriming_boundary2');
#		  runtst('primer_mispriming_boundary2_formatted');
		  runtst('primer_mispriming_long_lib');
		  runtst('primer_rat');
		  runtst('primer_human');
#		  runtst('long_seq');
		  runtst('primer_position_penalty');
#		  runtst('primer_position_penalty_formatted');
#		  runtst('p3-tmpl-mispriming');
#		  runtst('v4_old_tasks');
#		  runtst('v4_renewed_tasks');
#		  runtst('v4_new_tasks');
		  # Put primer_lib_amb_codes last because it is slow
#		  runtst('primer_lib_amb_codes');

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


sub runtst($) {
    $test = shift;
    print "$test";
    $testx = $test ;
    $testx =~s/_formatted//;
    $input = $testx . '_input';
    $output = $test . '_output';
    $tmp = $test . '.tmp';
    $cmd = "$exe < $input > $tmp";
    system  "$cmd";
	$r = perldiff $output, $tmp;
	if ($r == 0) {
	    print "[OK]\n";
	} else {
	    print "[FAILED]\n";
	}

}
 

