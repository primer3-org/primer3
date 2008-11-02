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

# This is output of create_tag_files.pl
my %docTags = (PRIMER_SEQUENCE_ID => "SEQUENCE_ID",
MARKER_NAME => "SEQUENCE_ID",
SEQUENCE_ID => "SEQUENCE_ID",
PRIMER_ERROR => "PRIMER_ERROR",
SEQUENCE => "SEQUENCE_TEMPLATE",
SEQUENCE_TEMPLATE => "SEQUENCE_TEMPLATE",
INCLUDED_REGION => "SEQUENCE_INCLUDED_REGION",
SEQUENCE_INCLUDED_REGION => "SEQUENCE_INCLUDED_REGION",
TARGET => "SEQUENCE_TARGET",
SEQUENCE_TARGET => "SEQUENCE_TARGET",
EXCLUDED_REGION => "SEQUENCE_EXCLUDED_REGION",
SEQUENCE_EXCLUDED_REGION => "SEQUENCE_EXCLUDED_REGION",
PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION => "SEQUENCE_INTERNAL_EXCLUDED_REGION",
SEQUENCE_INTERNAL_EXCLUDED_REGION => "SEQUENCE_INTERNAL_EXCLUDED_REGION",
PRIMER_LEFT_INPUT => "SEQUENCE_PRIMER",
SEQUENCE_PRIMER => "SEQUENCE_PRIMER",
PRIMER_INTERNAL_OLIGO_INPUT => "SEQUENCE_INTERNAL_OLIGO",
SEQUENCE_INTERNAL_OLIGO => "SEQUENCE_INTERNAL_OLIGO",
PRIMER_RIGHT_INPUT => "SEQUENCE_PRIMER_REVCOMP",
SEQUENCE_PRIMER_REVCOMP => "SEQUENCE_PRIMER_REVCOMP",
PRIMER_START_CODON_POSITION => "SEQUENCE_START_CODON_POSITION",
SEQUENCE_START_CODON_POSITION => "SEQUENCE_START_CODON_POSITION",
PRIMER_SEQUENCE_QUALITY => "SEQUENCE_QUALITY",
SEQUENCE_QUALITY => "SEQUENCE_QUALITY",
SEQUENCE_FORCE_LEFT_START => "SEQUENCE_FORCE_LEFT_START",
SEQUENCE_FORCE_LEFT_END => "SEQUENCE_FORCE_LEFT_END",
SEQUENCE_FORCE_RIGHT_START => "SEQUENCE_FORCE_RIGHT_START",
SEQUENCE_FORCE_RIGHT_END => "SEQUENCE_FORCE_RIGHT_END",
PRIMER_TASK => "PRIMER_TASK",
PRIMER_PICK_LEFT_PRIMER => "PRIMER_PICK_LEFT_PRIMER",
PRIMER_PICK_INTERNAL_OLIGO => "PRIMER_PICK_INTERNAL_OLIGO",
PRIMER_PICK_RIGHT_PRIMER => "PRIMER_PICK_RIGHT_PRIMER",
PRIMER_NUM_RETURN => "PRIMER_NUM_RETURN",
PRIMER_DEFAULT_PRODUCT => "PRIMER_PRODUCT_SIZE_RANGE",
PRIMER_PRODUCT_SIZE_RANGE => "PRIMER_PRODUCT_SIZE_RANGE",
PRIMER_PRODUCT_OPT_SIZE => "PRIMER_PRODUCT_OPT_SIZE",
PRIMER_PAIR_WT_PRODUCT_SIZE_LT => "PRIMER_PAIR_WT_PRODUCT_SIZE_LT",
PRIMER_PAIR_WT_PRODUCT_SIZE_GT => "PRIMER_PAIR_WT_PRODUCT_SIZE_GT",
PRIMER_MIN_SIZE => "PRIMER_MIN_SIZE",
PRIMER_INTERNAL_OLIGO_MIN_SIZE => "PRIMER_INTERNAL_MIN_SIZE",
PRIMER_INTERNAL_MIN_SIZE => "PRIMER_INTERNAL_MIN_SIZE",
PRIMER_DEFAULT_SIZE => "PRIMER_OPT_SIZE",
PRIMER_OPT_SIZE => "PRIMER_OPT_SIZE",
PRIMER_INTERNAL_OLIGO_OPT_SIZE => "PRIMER_INTERNAL_OPT_SIZE",
PRIMER_INTERNAL_OPT_SIZE => "PRIMER_INTERNAL_OPT_SIZE",
PRIMER_MAX_SIZE => "PRIMER_MAX_SIZE",
PRIMER_INTERNAL_OLIGO_MAX_SIZE => "PRIMER_INTERNAL_MAX_SIZE",
PRIMER_INTERNAL_MAX_SIZE => "PRIMER_INTERNAL_MAX_SIZE",
PRIMER_WT_SIZE_LT => "PRIMER_WT_SIZE_LT",
PRIMER_IO_WT_SIZE_LT => "PRIMER_INTERNAL_WT_SIZE_LT",
PRIMER_INTERNAL_WT_SIZE_LT => "PRIMER_INTERNAL_WT_SIZE_LT",
PRIMER_WT_SIZE_GT => "PRIMER_WT_SIZE_GT",
PRIMER_IO_WT_SIZE_GT => "PRIMER_INTERNAL_WT_SIZE_GT",
PRIMER_INTERNAL_WT_SIZE_GT => "PRIMER_INTERNAL_WT_SIZE_GT",
PRIMER_MIN_GC => "PRIMER_MIN_GC",
PRIMER_INTERNAL_OLIGO_MIN_GC => "PRIMER_INTERNAL_MIN_GC",
PRIMER_INTERNAL_MIN_GC => "PRIMER_INTERNAL_MIN_GC",
PRIMER_OPT_GC_PERCENT => "PRIMER_OPT_GC_PERCENT",
PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT => "PRIMER_INTERNAL_OPT_GC_PERCENT",
PRIMER_INTERNAL_OPT_GC_PERCENT => "PRIMER_INTERNAL_OPT_GC_PERCENT",
PRIMER_MAX_GC => "PRIMER_MAX_GC",
PRIMER_INTERNAL_OLIGO_MAX_GC => "PRIMER_INTERNAL_MAX_GC",
PRIMER_INTERNAL_MAX_GC => "PRIMER_INTERNAL_MAX_GC",
PRIMER_WT_GC_PERCENT_LT => "PRIMER_WT_GC_PERCENT_LT",
PRIMER_IO_WT_GC_PERCENT_LT => "PRIMER_INTERNAL_WT_GC_PERCENT_LT",
PRIMER_INTERNAL_WT_GC_PERCENT_LT => "PRIMER_INTERNAL_WT_GC_PERCENT_LT",
PRIMER_WT_GC_PERCENT_GT => "PRIMER_WT_GC_PERCENT_GT",
PRIMER_IO_WT_GC_PERCENT_GT => "PRIMER_INTERNAL_WT_GC_PERCENT_GT",
PRIMER_INTERNAL_WT_GC_PERCENT_GT => "PRIMER_INTERNAL_WT_GC_PERCENT_GT",
PRIMER_GC_CLAMP => "PRIMER_GC_CLAMP",
PRIMER_MIN_TM => "PRIMER_MIN_TM",
PRIMER_INTERNAL_OLIGO_MIN_TM => "PRIMER_INTERNAL_MIN_TM",
PRIMER_INTERNAL_MIN_TM => "PRIMER_INTERNAL_MIN_TM",
PRIMER_OPT_TM => "PRIMER_OPT_TM",
PRIMER_INTERNAL_OLIGO_OPT_TM => "PRIMER_INTERNAL_OPT_TM",
PRIMER_INTERNAL_OPT_TM => "PRIMER_INTERNAL_OPT_TM",
PRIMER_MAX_TM => "PRIMER_MAX_TM",
PRIMER_INTERNAL_OLIGO_MAX_TM => "PRIMER_INTERNAL_MAX_TM",
PRIMER_INTERNAL_MAX_TM => "PRIMER_INTERNAL_MAX_TM",
PRIMER_MAX_DIFF_TM => "PRIMER_MAX_DIFF_TM",
PRIMER_WT_TM_LT => "PRIMER_WT_TM_LT",
PRIMER_IO_WT_TM_LT => "PRIMER_INTERNAL_WT_TM_LT",
PRIMER_INTERNAL_WT_TM_LT => "PRIMER_INTERNAL_WT_TM_LT",
PRIMER_WT_TM_GT => "PRIMER_WT_TM_GT",
PRIMER_IO_WT_TM_GT => "PRIMER_INTERNAL_WT_TM_GT",
PRIMER_INTERNAL_WT_TM_GT => "PRIMER_INTERNAL_WT_TM_GT",
PRIMER_PRODUCT_MIN_TM => "PRIMER_PRODUCT_MIN_TM",
PRIMER_PRODUCT_OPT_TM => "PRIMER_PRODUCT_OPT_TM",
PRIMER_PRODUCT_MAX_TM => "PRIMER_PRODUCT_MAX_TM",
PRIMER_PAIR_WT_PRODUCT_TM_LT => "PRIMER_PAIR_WT_PRODUCT_TM_LT",
PRIMER_PAIR_WT_PRODUCT_TM_GT => "PRIMER_PAIR_WT_PRODUCT_TM_GT",
PRIMER_TM_SANTALUCIA => "PRIMER_TM_FORMULA",
PRIMER_TM_FORMULA => "PRIMER_TM_FORMULA",
PRIMER_SALT_CONC => "PRIMER_SALT_MONOVALENT",
PRIMER_SALT_MONOVALENT => "PRIMER_SALT_MONOVALENT",
PRIMER_INTERNAL_OLIGO_SALT_CONC => "PRIMER_INTERNAL_SALT_MONOVALENT",
PRIMER_INTERNAL_SALT_MONOVALENT => "PRIMER_INTERNAL_SALT_MONOVALENT",
PRIMER_DIVALENT_CONC => "PRIMER_SALT_DIVALENT",
PRIMER_SALT_DIVALENT => "PRIMER_SALT_DIVALENT",
PRIMER_INTERNAL_OLIGO_DIVALENT_CONC => "PRIMER_INTERNAL_SALT_DIVALENT",
PRIMER_INTERNAL_SALT_DIVALENT => "PRIMER_INTERNAL_SALT_DIVALENT",
PRIMER_DNTP_CONC => "PRIMER_DNTP_CONC",
PRIMER_INTERNAL_OLIGO_DNTP_CONC => "PRIMER_INTERNAL_DNTP_CONC",
PRIMER_INTERNAL_DNTP_CONC => "PRIMER_INTERNAL_DNTP_CONC",
PRIMER_SALT_CORRECTIONS => "PRIMER_SALT_CORRECTIONS",
PRIMER_DNA_CONC => "PRIMER_DNA_CONC",
PRIMER_INTERNAL_OLIGO_DNA_CONC => "PRIMER_INTERNAL_DNA_CONC",
PRIMER_INTERNAL_DNA_CONC => "PRIMER_INTERNAL_DNA_CONC",
PRIMER_SELF_ANY => "PRIMER_SELF_ANY",
PRIMER_INTERNAL_OLIGO_SELF_ANY => "PRIMER_INTERNAL_SELF_ANY",
PRIMER_INTERNAL_SELF_ANY => "PRIMER_INTERNAL_SELF_ANY",
PRIMER_SELF_END => "PRIMER_SELF_END",
PRIMER_INTERNAL_OLIGO_SELF_END => "PRIMER_INTERNAL_SELF_END",
PRIMER_INTERNAL_SELF_END => "PRIMER_INTERNAL_SELF_END",
PRIMER_MAX_END_STABILITY => "PRIMER_MAX_END_STABILITY",
PRIMER_NUM_NS_ACCEPTED => "PRIMER_NUM_NS_ACCEPTED",
PRIMER_INTERNAL_OLIGO_NUM_NS => "PRIMER_INTERNAL_NUM_NS_ACCEPTED",
PRIMER_INTERNAL_NUM_NS_ACCEPTED => "PRIMER_INTERNAL_NUM_NS_ACCEPTED",
PRIMER_WT_NUM_NS => "PRIMER_WT_NUM_NS",
PRIMER_IO_WT_NUM_NS => "PRIMER_INTERNAL_WT_NUM_NS",
PRIMER_INTERNAL_WT_NUM_NS => "PRIMER_INTERNAL_WT_NUM_NS",
PRIMER_MAX_POLY_X => "PRIMER_MAX_POLY_X",
PRIMER_INTERNAL_OLIGO_MAX_POLY_X => "PRIMER_INTERNAL_MAX_POLY_X",
PRIMER_INTERNAL_MAX_POLY_X => "PRIMER_INTERNAL_MAX_POLY_X",
PRIMER_MIN_THREE_PRIME_DISTANCE => "PRIMER_MIN_THREE_PRIME_DISTANCE",
PRIMER_PICK_ANYWAY => "PRIMER_PICK_ANYWAY",
PRIMER_SHOW_OLIGO_PROBLEMS => "PRIMER_SHOW_OLIGO_PROBLEMS",
PRIMER_LOWERCASE_MASKING => "PRIMER_LOWERCASE_MASKING",
PRIMER_EXPLAIN_FLAG => "PRIMER_EXPLAIN_FLAG",
PRIMER_LIBERAL_BASE => "PRIMER_LIBERAL_BASE",
PRIMER_FIRST_BASE_INDEX => "PRIMER_FIRST_BASE_INDEX",
PRIMER_MISPRIMING_LIBRARY => "PRIMER_MISPRIMING_LIBRARY",
PRIMER_INTERNAL_OLIGO_MISHYB_LIBRARY => "PRIMER_INTERNAL_MISHYB_LIBRARY",
PRIMER_INTERNAL_MISHYB_LIBRARY => "PRIMER_INTERNAL_MISHYB_LIBRARY",
PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS => "PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS",
PRIMER_MAX_MISPRIMING => "PRIMER_MAX_MISPRIMING",
PRIMER_INTERNAL_OLIGO_MAX_MISHYB => "PRIMER_INTERNAL_MAX_MISHYB",
PRIMER_INTERNAL_MAX_MISHYB => "PRIMER_INTERNAL_MAX_MISHYB",
PRIMER_MAX_TEMPLATE_MISPRIMING => "PRIMER_MAX_TEMPLATE_MISPRIMING",
PRIMER_INTERNAL_OLIGO_MAX_TEMPLATE_MISHYB => "PRIMER_INTERNAL_MAX_TEMPLATE_MISHYB",
PRIMER_INTERNAL_MAX_TEMPLATE_MISHYB => "PRIMER_INTERNAL_MAX_TEMPLATE_MISHYB",
PRIMER_PAIR_MAX_MISPRIMING => "PRIMER_PAIR_MAX_MISPRIMING",
PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING => "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING",
PRIMER_MIN_QUALITY => "PRIMER_MIN_QUALITY",
PRIMER_INTERNAL_OLIGO_MIN_QUALITY => "PRIMER_INTERNAL_MIN_QUALITY",
PRIMER_INTERNAL_MIN_QUALITY => "PRIMER_INTERNAL_MIN_QUALITY",
PRIMER_MIN_END_QUALITY => "PRIMER_MIN_END_QUALITY",
PRIMER_QUALITY_RANGE_MIN => "PRIMER_QUALITY_RANGE_MIN",
PRIMER_QUALITY_RANGE_MAX => "PRIMER_QUALITY_RANGE_MAX",
PRIMER_INSIDE_PENALTY => "PRIMER_INSIDE_PENALTY",
PRIMER_OUTSIDE_PENALTY => "PRIMER_OUTSIDE_PENALTY",
PRIMER_SEQUENCING_LEAD => "PRIMER_SEQUENCING_LEAD",
PRIMER_SEQUENCING_SPACING => "PRIMER_SEQUENCING_SPACING",
PRIMER_SEQUENCING_INTERVAL => "PRIMER_SEQUENCING_INTERVAL",
PRIMER_SEQUENCING_ACCURACY => "PRIMER_SEQUENCING_ACCURACY",
PRIMER_WT_COMPL_ANY => "PRIMER_WT_SELF_ANY",
PRIMER_WT_COMPL_END => "PRIMER_WT_SELF_END",
PRIMER_WT_SELF_ANY => "PRIMER_WT_SELF_ANY",
PRIMER_WT_SELF_END => "PRIMER_WT_SELF_END",
PRIMER_WT_REP_SIM => "PRIMER_WT_REP_SIM",
PRIMER_WT_SEQ_QUAL => "PRIMER_WT_SEQ_QUAL",
PRIMER_WT_END_QUAL => "PRIMER_WT_END_QUAL",
PRIMER_WT_POS_PENALTY => "PRIMER_WT_POS_PENALTY",
PRIMER_WT_END_STABILITY => "PRIMER_WT_END_STABILITY",
PRIMER_WT_TEMPLATE_MISPRIMING => "PRIMER_WT_TEMPLATE_MISPRIMING",
PRIMER_INTERNAL_WT_TEMPLATE_MISHYB => "PRIMER_INTERNAL_WT_TEMPLATE_MISHYB",
PRIMER_PAIR_WT_PR_PENALTY => "PRIMER_PAIR_WT_PR_PENALTY",
PRIMER_PAIR_WT_IO_PENALTY => "PRIMER_PAIR_WT_IO_PENALTY",
PRIMER_PAIR_WT_DIFF_TM => "PRIMER_PAIR_WT_DIFF_TM",
PRIMER_PAIR_WT_COMPL_ANY => "PRIMER_PAIR_WT_COMPL_ANY",
PRIMER_PAIR_WT_COMPL_END => "PRIMER_PAIR_WT_COMPL_END",
PRIMER_PAIR_WT_REP_SIM => "PRIMER_PAIR_WT_REP_SIM",
PRIMER_PAIR_WT_TEMPLATE_MISPRIMING => "PRIMER_PAIR_WT_TEMPLATE_MISPRIMING",
PRIMER_IO_WT_COMPL_ANY => "PRIMER_INTERNAL_WT_SELF_ANY",
PRIMER_INTERNAL_WT_SELF_ANY => "PRIMER_INTERNAL_WT_SELF_ANY",
PRIMER_IO_WT_COMPL_END => "PRIMER_INTERNAL_WT_SELF_END",
PRIMER_INTERNAL_WT_SELF_END => "PRIMER_INTERNAL_WT_SELF_END",
PRIMER_IO_WT_REP_SIM => "PRIMER_INTERNAL_WT_REP_SIM",
PRIMER_INTERNAL_WT_REP_SIM => "PRIMER_INTERNAL_WT_REP_SIM",
PRIMER_IO_WT_SEQ_QUAL => "PRIMER_INTERNAL_WT_SEQ_QUAL",
PRIMER_INTERNAL_WT_SEQ_QUAL => "PRIMER_INTERNAL_WT_SEQ_QUAL",
PRIMER_IO_WT_END_QUAL => "PRIMER_INTERNAL_WT_END_QUAL",
PRIMER_INTERNAL_WT_END_QUAL => "PRIMER_INTERNAL_WT_END_QUAL",
PRIMER_PAIR_ANY => "PRIMER_PAIR_ANY",
PRIMER_PAIR_END => "PRIMER_PAIR_END",
PRIMER_PRIMERS_IN_PAIRS_MUST_BE_UNIQUE => "PRIMER_PRIMERS_IN_PAIRS_MUST_BE_UNIQUE",
P3_FILE_ID => "P3_FILE_ID",
PRIMER_FILE_FLAG => "P3_FILE_FLAG",
P3_FILE_FLAG => "P3_FILE_FLAG",
PRIMER_COMMENT => "P3_COMMENT",
COMMENT => "P3_COMMENT",
P3_COMMENT => "P3_COMMENT");

my @files = ( 'primer_boundary', # Put the quickest tests first.
		  'p3_3_prime_n',
		  'primer_tm_lc_masking',
		  'primer_tm_lc_masking_formatted',
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
		  'v4_old_tasks',
		  'v4_renewed_tasks',
		  'v4_new_tasks',
		  # Put primer_lib_amb_codes last because it is slow
		  'primer_lib_amb_codes');


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
		or ($si_file eq ".svn")) {
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
			$text_output .= $docTags{$line}."=\n";
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
				print "Error - in $file unknown Tag ".$lineA[0]." is used.\n";
				$fin_target = $target."_error";
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
