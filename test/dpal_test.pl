
# Usage: perl dpal_test.pl

use strict;
use warnings 'all';
use Getopt::Long;
use Carp;

our $test_count = 0;
our $exit_status = 0;
our $do_valgrind;
our $valgrind_format
    = "/usr/local/bin/valgrind --leak-check=yes "
    . " --show-reachable=yes --log-file-exactly=ntdpal.%0.4d.valg ";

sub runtest($$$$$) {
    my ($desc, $ntdpal_args, $infile, $benchfile, $outfile) = @_;
    print $desc, '...';
    open Q, $infile or confess "open $infile: $!";

    # Reopen STDOUT to get all ntdpal output in one file
    open OLDOUT, ">&STDOUT" or confess "Cannot dup STDOUT: $!";
    open STDOUT, '>', $outfile or confess "open STDOUT '>' $outfile: $!";
    while (my $in = <Q>) {
	$test_count++;
	my $valgrind_prefix 
	    = $do_valgrind ? sprintf($valgrind_format, $test_count) : '';
	system "$valgrind_prefix ../src/ntdpal $ntdpal_args $in";
    }
    close Q;
    open STDOUT, ">&OLDOUT" or confess "Cannot dup OLDOUT: $!";
    close OLDOUT;
    my $r = system "diff $benchfile $outfile";
    print $r == 0 ? "OK\n" : "FAILED\n";
    if ($r) { $exit_status = -1; }
}

select STDERR;

if (!GetOptions('valgrind', \$do_valgrind)) {
    print "Usage: $0 [ --valgrind ]\n";
    exit -1;
}

die "Cannot execute ../src/ntdpal" unless -x '../src/ntdpal';

# Test error handling on over-long input sequence:
my $valgrind_prefix = $do_valgrind ? sprintf($valgrind_format, $test_count) : '';
print 'Error handling of too-long sequence...';
my $r = system '$valgrind_prefeix ../src/ntdpal ACGTGTTCGTCGTAAATAACATGCTATATT GACGTAGACAACCCTGTGTTTAGCCTGCGTTTTGTGCCATCCTAATGCTTTACTAGATCACTGAGCCACCTCCCAAGGACTACACCTAGCGGTATTTCGTACATTAACTAGGATCCTTTTCCACATGGACTACAATGTCTGCCGAGCATGCGATGGGGTACCGCGCCCGCGCACATACGCGCGCAGAGCTTTTGGAGGCATACCTACACCGGCGAGGGGCTGCGGTTTATTGACACTGAAACGGGATAACGAGTCGCTGAATTGAGCCAAAAATATGCAAGCGTCACAAATTGTGACAAAAATTTTAAAGGAAAAATTAGACCATTGATTCTGAAGTGGTGCGTATAGGACCAGTCGTGGCAATGAGACCGATTTGAGTAGCACTAGCTCAAACACTGTCTGGGTCGCCATCAAGGCCACAAGAACTTAAGCAGCCGTCACCCTATAGAAGGTTAAGCGACGGTTAGGGCTTCTGGCAACGAAAGTTGTCGGTTCGTCCTGTGCCAACGTGTGGCAAAGTCTACTATGATTCGATTGTTGACGTGTCGACAGGCTGTTTCGCTGGATACCCCACCTTGATAATTTTTCTCGTCGAACGCTAGCAGTTTTTTTTTCAACGGCCCGGAATCTGTAAGAGGCCGTTGCAGGAACGCGTGTGTATGTAAATGCCCACTACTTCTGTTATGTACCCAAATGGCGTGCGGCGTGGATGTATAGTGTCGACCCTCCATAATCGGGCGGACGGTCGTGGGGTATGTATGATCTTCGGCACTGATTCGCCTCGAGTCTATATGTTCTTAATCCAGACCTTCGGGGAAAGCCTACTTTCCATCCGTTGTCTAGCGTCATGCCAGTGACTACTGTTGTATTGTCTGGTTCCTAAGATAGCCATGGATTCCGGACATCGACGATGCACAAGAGCGTTAGCGCTGGTGTGCAACGCAACGTCGCGAAGGCTGGGTTACAGCGTGATCTCCTGGCTGCACCCAGATGCAGAGGGACATACCTACGATGAATAGGTGCGTCTGTTTATAAACGCCCAATCCTAGCAAAAATCACAACTAAGACAGTGTATGGAAGACCCACCAGTTGTGGGCGAATGGTCAGGTATACAAGATCGTGTCAAGACGGAACTTAAGCTTCTGTGCGCTCTCCATGCGAGCTGGTACGTCTGGACGGCGAGGTATGAGTGAATGACCATCCATGGCAACTTTCGTGTTCTACGACAGATACGAGCTCGACGGACGACCTGGTGACCAGTAGTATATGCGCGTCCGTCGGCCAGACTTTCCAAACGCCCTTTCAACGAGATACATGCGAACACGCTACAATTTCTCGTTCCGTCTAAAGTCGATACTCGCAAGCCCAGGCCCGTTACTACAACGCTGTTAATAGGATCAGAAGGGCCATAAGACTTTGGCAGCGGTAGCTAGGAAAGTGATGGTTGTGATGGCCCTAGTAAGGAGTCAGCCATCTACCCAACTATTTGAATGGGACCATAGCCAAGGGACCCAGCTGTTCCTTAGAAACCTGGTGACTCCCTTAGCCAATTGTGTAACTTCGTGCGTGCCAGTATTACACCTATAATCACAAGACCCCTTCAATACGAGTCCTGTGGCGTAGTGTTCCATCAAAACAATCAAGAACAGATTTCCGGTCCCCGTTGTGTTGGGATCTAGCGGACGTTGTCGGTAGATCAATAACGTAAATGCGAATCGAAGTTCTCTGGCCTAAAACAACTGCGCGCAGGGCCTCCGGTCATTGCATCTTTCTTGTCTCTCGTGAGGGCGTGATTCGTTTACCTGGAGCGAGCCGGGCACAAGAGCTATGGATTATTGGCTGGTGCAAAAACCATTCTAGCTACAATTATACTCGCGTGTCGACGATAAGAGTGAAATCACTGCGTAGGCAAACTGCCGGGTCACCAAGAGAGGCTGATACCGCGGTTCACCC l > dpal.tmp 2>&1';
open X, 'dpal.tmp';
my @foo = <X>;
close X;
if ($foo[0] eq "Error: Sequence 2 longer than DPAL_MAX_ALIGN and alignment is requested\n") {
    print "OK\n"
} else {
    print "FAILED\n";
    $exit_status = -1;
}

runtest('Default implementations + alignment',
	"",  "dpal_input", 'dpal_output',
	'dpal_output1.tmp');

runtest('Default implementations + NO alignment 1',
	"-s",  'dpal_input', 'dpal_score_output',
	'dpal_score_output1.tmp');

runtest('Default implementations + NO alignment 2',
	"-s", "dpal_long_input", 'dpal_long_score_output',
	"dpal_long_score_output1.tmp");

runtest('Force _dpal_generic',
	'-s -f1', "dpal_input", 'dpal_score_output',
	"dpal_score_output4.tmp");

runtest('Force _dpal_long_nopath_generic 1',
	'-s -f2', "dpal_input", 'dpal_score_output',
	"dpal_score_output2.tmp");

runtest('Force _dpal_long_nopath_generic 2',
	'-s -f2', "dpal_long_input", 'dpal_long_score_output',
	"dpal_long_score_output2.tmp");

runtest('Force long maxgap1 functions 1',
	'-s -f3',  "dpal_input", 'dpal_score_output',
	"dpal_score_output3.tmp");

runtest('Force long maxgap1 functions 2',
	'-s -f3',  "dpal_long_input", 'dpal_long_score_output', 
	"dpal_long_score_output3.tmp");

if ($do_valgrind) {
    # Assume this is Unix/Linux envrionment, so
    # we have grep.
    my $r = system "grep ERROR ntdpal.*.*valg | grep -v 'ERROR SUMMARY: 0 errors'";
    if (!$r) {	   # !$r because grep returns 0 if something is found,
	# and if something is found, we have a problem.
	$exit_status = -1;
    }
    $r = system "grep 'definitely lost' ntdpal.*.*valg | grep -v '0 bytes'";
    if (!$r) {
	$exit_status = -1;
    }
    $r = system "grep 'possibly lost' ntdpal.*.*valg  | grep -v '0 bytes'";
    if (!$r) {
	$exit_status = -1;
    }
}

exit $exit_status;

# sub runtest($$$$$) {
#     my ($desc, $ntdpal_args, $infile, $benchfile, $outfile) = @_;
#     print $desc, '...';
#     open Q, $infile or confess "open $infile: $!";
# 
#     # Reopen STDOUT to get all ntdpal output in one file
#     open OLDOUT, ">&STDOUT" or confess "Cannot dup STDOUT: $!";
#     open STDOUT, '>', $outfile or confess "open STDOUT '>' $outfile: $!";
#     while (my $in = <Q>) {
# 	$test_count++;
# 	my $valgrind_prefix 
# 	    = $do_valgrind ? sprintf($valgrind_format, $test_count) : '';
# 	system "$valgrind_prefix ../src/ntdpal $ntdpal_args $in";
#     }
#     close Q;
#     open STDOUT, ">&OLDOUT" or confess "Cannot dup OLDOUT: $!";
#     close OLDOUT;
#     $r = system "diff $benchfile $outfile";
#     print $r == 0 ? "OK\n" : "FAILED\n";
# 
# }
