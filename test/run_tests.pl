#!/usr/bin/perl

# ======================================================================
# Runs all test suites one after the other, also when one of them
# fails, collects the failures and prints a final summary with the
# time each suite took.
#
# Usage (from the test directory):
#   perl run_tests.pl [ --windows ] [ --valgrind | --valgrinda | --valgrindb ]
#                     [ --fast ] [ --only=SUITE[,SUITE...] ] [ --help ]
#
# --valgrinda and --valgrindb split the valgrind run in two parts, as
# the CI jobs do (make valgrnda / make valgrndb):
#   --valgrinda: the small suites with --valgrind, a3test.pl and
#                p3test.pl with --valgrinda (first part of the p3 tests)
#   --valgrindb: only p3test.pl with --valgrindb (second part)
#
# Exit status: 0 if every suite passed, 1 otherwise.
#
# On GitHub Actions (GITHUB_ACTIONS=true) the output of each suite is
# folded into a group, every failure becomes an error annotation and
# the summary is added to the job summary ($GITHUB_STEP_SUMMARY).
#
# This file is part of primer3 and is distributed under the terms of
# the GNU General Public License, version 2 or (at your option) any
# later version.  See the file gpl-2.0.txt in the source distribution.
# ======================================================================

use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(time);

$| = 1;

# Name, script and the options each script understands.  Suites with
# split_valgrind understand --valgrinda (and p3test.pl --valgrindb).
my @SUITES = (
    { name => 'cmdline', script => 'cmdline_test.pl', opts => [qw(windows valgrind)] },
    { name => 'oligotm', script => 'oligotm_test.pl', opts => [qw(windows valgrind)] },
    { name => 'dpal',    script => 'dpal_test.pl',    opts => [qw(windows valgrind)] },
    { name => 'thal',    script => 'thal_test.pl',    opts => [qw(windows valgrind)] },
    { name => 'masker',  script => 'masker_test.pl',  opts => [qw(windows valgrind)],
      no_windows => 1 },
    { name => 'a3',      script => 'a3test.pl',       opts => [qw(windows valgrind)],
      extra => ['--action'], split_valgrind => 1 },
    { name => 'p3',      script => 'p3test.pl',       opts => [qw(windows valgrind fast)],
      extra => ['--action'], split_valgrind => 1 },
);

my %args;
if (!GetOptions(\%args, 'windows', 'valgrind', 'valgrinda', 'valgrindb',
                'fast', 'only=s', 'help')
    || $args{'help'}) {
    print "Usage: perl run_tests.pl [ --windows ] [ --valgrind | --valgrinda | --valgrindb ]\n",
          "                         [ --fast ] [ --only=SUITE[,SUITE...] ]\n",
          "Suites: ", join(' ', map { $_->{name} } @SUITES), "\n",
          "--fast is passed to p3test.pl only.\n",
          "--valgrinda: small suites with --valgrind, a3 and p3 with --valgrinda.\n",
          "--valgrindb: only p3 with --valgrindb.\n";
    exit($args{'help'} ? 0 : 1);
}
my $n_valgrind = grep { $args{$_} } qw(valgrind valgrinda valgrindb);
if ($n_valgrind > 1) {
    print "run_tests.pl: Specify only one of --valgrind, --valgrinda and --valgrindb\n";
    exit 1;
}
if ($args{'windows'} && $n_valgrind) {
    print "run_tests.pl: Cannot specify both valgrind and --windows\n";
    exit 1;
}
if ($args{'valgrindb'} && defined $args{'only'}) {
    print "run_tests.pl: --valgrindb runs only p3; do not combine it with --only\n";
    exit 1;
}

my %only;
if (defined $args{'only'}) {
    my %known = map { $_->{name} => 1 } @SUITES;
    for my $s (split /,/, $args{'only'}) {
        if (!$known{$s}) {
            print "run_tests.pl: Unknown suite '$s'. Suites: ",
                  join(' ', map { $_->{name} } @SUITES), "\n";
            exit 1;
        }
        $only{$s} = 1;
    }
}
if ($args{'valgrindb'}) {
    %only = (p3 => 1);
}

my $gha = defined $ENV{'GITHUB_ACTIONS'} && $ENV{'GITHUB_ACTIONS'} eq 'true';

my $interrupted = 0;
$SIG{'INT'} = sub { $interrupted = 1; };

my @results;
my $total_start = time();

for my $suite (@SUITES) {
    last if $interrupted;
    next if %only && !$only{$suite->{name}};

    my %res = (name => $suite->{name}, failures => [], time => 0);
    if ($suite->{no_windows} && $args{'windows'}) {
        $res{status} = 'SKIPPED';
        push @results, \%res;
        next;
    }

    my @cmd_args = map { "--$_" } grep { $args{$_} } @{$suite->{opts}};
    if ($args{'valgrinda'}) {
        push @cmd_args, $suite->{split_valgrind} ? '--valgrinda' : '--valgrind';
    } elsif ($args{'valgrindb'}) {
        push @cmd_args, '--valgrindb';
    }
    push @cmd_args, @{$suite->{extra}} if $suite->{extra};
    my $cmd = join(' ', qq{"$^X"}, $suite->{script}, @cmd_args) . ' 2>&1';

    if ($gha) {
        print "::group::$suite->{name} ($suite->{script})\n";
    } else {
        print "\n", '=' x 70, "\n";
        print "TESTING $suite->{name}: $cmd\n";
        print '=' x 70, "\n";
    }

    my $start = time();
    my $exit;
    my $current = '';   # name of the test that is running
    my $last = '';      # last non-empty output line
    if (open(my $fh, '-|', $cmd)) {
        while (my $line = <$fh>) {
            print $line;
            $line =~ s/[\r\n]+$//;
            my $name = test_name($line);
            $current = $name if defined $name;
            if ($line =~ /\[FAILED\]/) {
                my $d = failure_description($line, $current);
                push @{$res{failures}}, $d if defined $d;
            }
            $last = $line if $line =~ /\S/;
        }
        close($fh);
        $exit = $?;
    } else {
        print "Cannot run $cmd: $!\n";
        $exit = -1;
    }
    $res{time} = time() - $start;

    print "::endgroup::\n" if $gha;

    if ($exit == -1) {
        $res{exit} = 'not started';
    } elsif ($exit & 127) {
        $res{exit} = 'signal ' . ($exit & 127);
        $interrupted = 1 if ($exit & 127) == 2;   # SIGINT
    } else {
        $res{exit} = $exit >> 8;
    }

    if ($res{exit} ne '0' || @{$res{failures}}) {
        $res{status} = 'FAILED';
        if (!@{$res{failures}}) {
            my $d = "exit status $res{exit}";
            $last =~ s/^\s+|\s+$//g;
            $d .= ", last output: $last" if $last ne '';
            push @{$res{failures}}, $d;
        }
    } else {
        $res{status} = 'OK';
    }
    push @results, \%res;
}

my $total_time = time() - $total_start;
my @failed = grep { $_->{status} eq 'FAILED' } @results;
my $n_fail = 0;
$n_fail += scalar @{$_->{failures}} for @failed;

# ---------------------------------------------------------------- summary
print "\n", '=' x 70, "\nSUMMARY\n", '=' x 70, "\n";
for my $r (@results) {
    my $extra = '';
    if ($r->{status} eq 'FAILED') {
        my $n = scalar @{$r->{failures}};
        $extra = "  ($n failure" . ($n == 1 ? '' : 's') . ")";
    }
    printf "%-10s %-8s %10s%s\n", $r->{name}, $r->{status},
           fmt_time($r->{time}), $extra;
}
printf "%-10s %-8s %10s\n", 'total', '', fmt_time($total_time);

if (@failed) {
    print "\nFailures:\n";
    for my $r (@failed) {
        for my $f (@{$r->{failures}}) {
            print "  $r->{name}: $f\n";
        }
    }
}

my $final;
if ($interrupted) {
    $final = 'INTERRUPTED - not all suites were run';
} elsif (@failed) {
    $final = scalar(@failed) . ' of ' . scalar(@results)
           . " suites FAILED ($n_fail failure" . ($n_fail == 1 ? "" : "s") . ")";
} else {
    $final = 'ALL TESTS PASSED';
}
print "\nRESULT: $final in ", fmt_time($total_time), "\n\n";

if ($gha) {
    for my $r (@failed) {
        for my $f (@{$r->{failures}}) {
            print "::error title=$r->{name}::", gha_escape($f), "\n";
        }
    }
    write_step_summary();
}

exit(($interrupted || @failed) ? 1 : 0);

# ======================================================================

# The scripts print "NAME..." before running a test; diff output and
# the result follow on the same or on later lines.  Returns NAME or
# undef.  Names need a letter, so "1......" (a sub-step) is no name.
sub test_name {
    my ($line) = @_;
    return undef if $line =~ /^\s*[<>]/;   # diff output
    if ($line =~ /^\s*(.*?[A-Za-z].*?)\s*\.{3,}/) {
        return $1;
    }
    return undef;
}

# Turns a line containing [FAILED] into a description of the failed
# test: the test name plus what the line says besides [FAILED].
# Returns undef for the closing lines of a suite, which only repeat
# that something failed.
sub failure_description {
    my ($line, $current) = @_;
    return undef if $line =~ /At least one test failed/;
    return undef if $line =~ /In addition, at least 1 test FAILED/;

    my $name = test_name($line);
    my $rest = $line;
    if (defined $name) {
        $rest =~ s/^\s*.*?[A-Za-z].*?\s*\.{3,}//;
    } elsif ($line =~ /valgrind/i) {
        # The valgrind check at the end of a script covers all its tests,
        # not the test that ran last.
        $name = '';
    } else {
        $name = $current;
    }
    $rest =~ s/\[FAILED\]//;
    $rest =~ s/(^|\s)\d+(\.\d+)?s(?=\s|$)/ /g;   # run time like "3s"
    $rest =~ s/^[\s.\d]+$//;                      # "1......" sub-steps
    $rest =~ s/^\s+|\s+$//g;

    if ($name eq '') {
        return $rest eq '' ? 'failed' : $rest;
    }
    return $rest eq '' ? $name : "$name: $rest";
}

sub fmt_time {
    my ($s) = @_;
    my $m = int($s / 60);
    return sprintf("%d:%04.1f", $m, $s - 60 * $m);
}

sub gha_escape {
    my ($s) = @_;
    $s =~ s/%/%25/g;
    $s =~ s/\r/%0D/g;
    $s =~ s/\n/%0A/g;
    return $s;
}

sub write_step_summary {
    my $file = $ENV{'GITHUB_STEP_SUMMARY'};
    return if !defined $file || $file eq '';
    open(my $fh, '>>', $file) or return;
    print $fh "### Primer3 tests: $final in ", fmt_time($total_time), "\n\n";
    print $fh "| Suite | Result | Time | Failures |\n|---|---|---|---|\n";
    for my $r (@results) {
        my $n = $r->{status} eq 'FAILED' ? scalar @{$r->{failures}} : '';
        print $fh "| $r->{name} | $r->{status} | ", fmt_time($r->{time}),
                  " | $n |\n";
    }
    if (@failed) {
        print $fh "\n**Failures**\n\n";
        for my $r (@failed) {
            for my $f (@{$r->{failures}}) {
                my $e = $f;
                $e =~ s/([\\`*_|<>])/\\$1/g;
                print $fh "- $r->{name}: $e\n";
            }
        }
    }
    print $fh "\n";
    close($fh);
}
