#!/usr/bin/perl -w

# (c) Copyright 1996,1997,1998,1999,2000,2001,2004,2006,2007, 2008
# Whitehead Institute for Biomedical Research, Steve Rozen, 
# Andreas Untergasser and Helen Skaletsky
# All rights reserved.
# 
#     This file is part of the primer3 suite and libraries.
# 
#     The primer3 suite and libraries are free software;
#     you can redistribute them and/or modify them under the terms
#     of the GNU General Public License as published by the Free
#     Software Foundation; either version 2 of the License, or (at
#     your option) any later version.
# 
#     This software is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this software (file gpl-2.0.txt in the source
#     distribution); if not, write to the Free Software
#     Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
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

use strict;
use Cwd;
use File::Copy;

# This is what wil be changed:
my %docTags = (PRIMER_NUM_NS_ACCEPTED => "PRIMER_MAX_NS_ACCEPTED",
PRIMER_INTERNAL_NUM_NS_ACCEPTED => "PRIMER_INTERNAL_MAX_NS_ACCEPTED",
PRIMER_MAX_DIFF_TM => "PRIMER_PAIR_MAX_DIFF_TM",
PRIMER_PAIR_COMPL_ANY => "PRIMER_PAIR_MAX_COMPL_ANY",
PRIMER_PAIR_COMPL_END => "PRIMER_PAIR_MAX_COMPL_END",
PRIMER_SELF_ANY => "PRIMER_MAX_SELF_ANY",
PRIMER_INTERNAL_SELF_ANY => "PRIMER_INTERNAL_MAX_SELF_ANY",
PRIMER_SELF_END => "PRIMER_MAX_SELF_END",
PRIMER_INTERNAL_SELF_END => "PRIMER_INTERNAL_MAX_SELF_END");

# Directory Name
my $directory = cwd;

my @directoryFiles;

if ( opendir(DIR, $directory)) { 
	push @directoryFiles, readdir(DIR);
	closedir(DIR);
} else {
	print "Couldn't open $directory for reading!\n";
}


print ("Start Translation...\n");
my $input;
my $output;

foreach my $si_file (@directoryFiles) {
	if (($si_file eq ".") or ($si_file eq "..")
		or ($si_file eq ".svn")
		or ($si_file =~ /_formatted_output$/)
		or ($si_file =~ /^dpal_/)) {
		next;
	}
#	print $si_file."\n";
	if ($si_file =~ /.in$/) {
#		print $si_file."\n";
		$input = $si_file;
		$output = $si_file;
		translate_file($input,$output);
	}
	if ($si_file =~ /.out$/) {
#		print $si_file."\n";
		$input = $si_file;
		$output = $si_file;
		translate_file($input,$output);
	}
    if ($si_file =~ /_input$/) {
#       print $si_file."\n";
        $input = $si_file;
        $output = $si_file;
        translate_file($input,$output);
    }
    if ($si_file =~ /_output$/) {
#       print $si_file."\n";
        $input = $si_file;
        $output = $si_file;
        translate_file($input,$output);
    }
}
#translate_file("primer_task_input","primer_task_input_tr");


print ("\nTranslation finished\n");


##########################################
# Reads a File into a string and returns #
#      0 for success                     #
#      or error message                  #
##########################################
sub translate_file {
	my $file = shift;	# File to read
	my $target = shift; # File to save in
	my $fin_target = $target;
	my $error = 0;

	my $text_input = "";
	my $text_output = "";
	my @string_array;
	my @lineA;
	
	# Read the file in a string
	$error = read2String(\$file, \$text_input);
	if ($error ne "0"){
		print $error;
		return $error;
	}
	my $backup = $file."_bakup";
#	$error = string2file($backup, $text_input);
	if ($error ne "0"){
		print $error;
	}
	
	# Split the string into single lines
	@string_array = split '\n', $text_input;
	
	# Handle each line separate
	foreach my $line (@string_array) {
		# First handle the = only lines
		if ($line eq "="){
			$text_output .= "=\n";
		} 
		# Now handle lines ending with =
		elsif ($line =~ /=$/) {
			$line =~ s/=$//;
			$text_output .= $line."=\n";
		}
		# Catch lines without a = (should not happen)
		elsif (!($line =~ /=/)) {
			print "\nError in $file at Line: ". $line."\n\n";
			$fin_target = $target."_error";
		}
		# Translate the old tags into new tags
		else {
			@lineA = split '=', $line;
			if (defined $docTags{$lineA[0]}){
				$text_output .= $docTags{$lineA[0]}."=".$lineA[1]."\n";
			} else {
				$text_output .= $lineA[0]."=".$lineA[1]."\n";
			}
		}
	}
	
	# $text_output
	
	# Save the file
	$error = string2file($fin_target, $text_output);
	if ($error ne "0"){
		print $error;
	}

    return $error;
}


##########################################
# Reads a File into a string and returns #
#      0 for success                     #
#      or error message                  #
##########################################
sub read2String {
	my $file = shift;	# File Name
	my $target = shift; # String to save in
	my $error = "0";
	
	if (open (FILE, ("<".${$file}))) { 
		while (<FILE>) {
		  	${$target} .= $_;
		}
		close(FILE);
		
		# solve the newline problem with other platforms
		if ( ${$target} =~ /\r\n/ ) {
			${$target} =~ s/\r\n/\n/g;
		}
		if ( ${$target} =~ /\r/ ) {
			${$target} =~ s/\r/\n/g;
		}
	} 
	else {
		$error = "Couldn't open ${$file} for reading!\n";
	}

    return $error;
}


##########################################
# Writes a string int a File and returns #
#      0 for success                     #
#      or error message                  #
##########################################
sub string2file {
	my $file = shift;	# File Name
	my $string = shift; # String to save in
	my $error = 0;
	
	if (open FILE, ">", $file) { 
		print FILE $string;
		close(FILE);
		
	} 
	else {
		$error = "Couldn't open $file for writing!\n";
	}

    return $error;
}
