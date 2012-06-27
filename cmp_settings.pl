#!/usr/bin/perl
#
# Logical compare (i.e. diff) of two settings files.
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
	unless ($line =~ /(\S*)=(.*)/) { print STDERR "wrong line format: $line\n"; next; }
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
foreach my $tag (sort (keys %tags1)) {
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






