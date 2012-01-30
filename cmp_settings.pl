#!/usr/bin/perl

use strict;
use warnings "all";
use Getopt::Long;

my $file1;
my $file2;

my %tags1;
my %tags2;

sub read_file($$)
{
    my ($file, $tags) = @_;
    open IN, $file or die "open $file: $!\n";
    # skip first line
    my $line = <IN>;
    while ($line = <IN>) {
	chomp($line);
	# skip empty lines and comments
	if ($line =~ /^\s*$/) { next; }
	if ($line =~ "^#") { next; }
	unless ($line =~ /(\S*)=(\S*)/) { print STDERR "wrong line format: $line\n"; next; }
	my $tag = $1;
	my $value = $2;
	$tags->{$tag} = $value;
    }
    close IN;
}

if (!GetOptions("file1=s" => \$file1, "file2=s" => \$file2) || !defined($file1) || !defined($file2)) {
    die "\nUsage:\n"
        . "$0 -file1 <filename1> -file2 <filename2>\n";
}

read_file($file1, \%tags1);
read_file($file2, \%tags2);

my @only1;

# cmp tags1 with tags2 print anything common different
print "Common tags with different values:\n";
foreach my $tag (keys %tags1) {
    if (defined($tags2{$tag})) {
	if ($tags1{$tag} ne $tags2{$tag}) {
	    print "\t$tag:\n\t\t$file1: $tag=$tags1{$tag}\n\t\t$file2: $tag=$tags2{$tag}\n"
	}
    } else {
	push(@only1, $tag);
    }
}

print "Tags that exist only in $file1:\n";
foreach my $tag (@only1) {
    print "\t$tag=$tags1{$tag}\n";
}

print "Tags that exist only in $file2:\n";
foreach my $tag (keys %tags2) {
    if (!defined($tags1{$tag})) {
	print "\t$tag=$tags2{$tag}\n";
    }
}






