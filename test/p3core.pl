#!/usr/bin/perl

# Approximate perl re-implementation of the primer3_core executable
# for testing.
#
# For usage, see the usage statement in the code, below.
#
# ======================================================================
# (c) Copyright 1996,1997,1998,1999,2000,2001,2004,2006,2007, 2008
# Whitehead Institute for Biomedical Research, Steve Rozen, 
# Andreas Untergasser and Helen Skaletsky
# All rights reserved.
# 
#   This file is part of the primer3 suite.
#
#   The primer3 suite is free software; you can
#   redistribute them and/or modify them under the terms of the GNU
#   General Public License as published by the Free Software Foundation;
#   either version 2 of the License, or (at your option) any later
#   version.
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

use Carp;
use blib '../../../p3perl/trunk';
use Primer3 ':all';

sub set_setters();
sub main();

our %dispatch;
our $sa ;
our $gs;
our $tag;
our $value;
our $err = "0" ;
our $retval = 0 ;
our $explain_flag = 0;
our $show_oligo_problems = 0;
our $printargs = 0 ; # set to 1 if dumping arguments to choose_primers

main();

sub main() {

    select STDOUT; $| = 1;  # Not sure we need this ....
    set_setters();
    $gs = pl_create_global_settings();
    if ($printargs) {
     	pl_set_dump($gs) ;
        }

    $/ = "\n=\n";
    while (1) {
	my $rec = <>;
	last if !defined $rec;
	chomp $rec;
	if (!$rec) {
	    confess "Record $. is empty\n";
	}
	$sa = pl_create_seq_arg();
	
	my %rec;
	my @rec = split /\n/, $rec;
	my $tag_found = 0;
	for my $line (@rec) {
	    print "$line\n";
	    if ($line !~ /(^\w+\s*)=(.*)/) {
		confess "Bad line $line\n";
	    }
	    ($tag, $value) = ($1, $2);
	    $tag_found = 1;
	    if ($dispatch{$tag})  {
		&{$dispatch{$tag}}($value);
	    } else { confess "no call for $tag" }
	}
	if (!$tag_found) {
	    confess "Record $. is empty\n";
	}
	if (!$err) {
	    $retval = pl_choose_primers($gs, $sa); 
	    p3_print_boulder($gs, $sa, $retval, $explain_flag, $show_oligo_problems);
           # boulder_print generates the final '='
	} else 	{
	    print "$err=\n" ;
	}
	pl_destroy_seq_args($sa);
	$err = "0" ;
    }
}

sub set_setters() {

# PER SEQUENCE INPUTS ==================================================
     $dispatch{'SEQUENCE'} = 
sub ($) { my $v = shift; pl_set_sa_sequence($sa, $v) }; 
     $dispatch{'SEQUENCE_DNA'} = 
sub ($) { my $v = shift; pl_set_sa_sequence($sa, $v) }; 
     $dispatch{'PRIMER_SEQUENCE_QUALITY'} = 
sub ($) { my $v = shift;
	  my @nums = split( /\s+/, $v) ;
	  my $n = shift @nums ;
	  pl_set_sa_empty_quality($sa);     
	  while(defined $n && $n ne "") {
	      pl_sa_add_to_quality_array($sa, $n) ;
	      $n = shift @nums ;
	    } 
          } ;
    $dispatch{'PRIMER_SEQUENCE_ID'} = 
sub ($) { my $v = shift; my $i = pl_set_sa_sequence_name($sa, $v) };     
    $dispatch{'SEQUENCE_ID'} = 
sub ($) { my $v = shift; my $i = pl_set_sa_sequence_name($sa, $v) };     
     $dispatch{'MARKER_NAME'} = 
sub ($) { my $v = shift;  my $i = pl_set_sa_sequence_name($sa, $v) };     
    $dispatch{'PRIMER_LEFT_INPUT'} = 
sub ($) { my $v = shift;  my $i = pl_set_sa_left_input($sa, $v) };     
    $dispatch{'PRIMER_RIGHT_INPUT'} = 
sub ($) { my $v = shift;  my $i = pl_set_sa_right_input($sa, $v) };     
    $dispatch{'PRIMER_INTERNAL_OLIGO_INPUT'} = 
sub ($) { my $v = shift;  my $i = pl_set_sa_internal_input($sa, $v) };     
    $dispatch{'TARGET'} = 
sub ($) { my $v = shift;
	  my @nums = split /[, ]/, $v ;
	  my $f = shift @nums ;
	  my $s = shift @nums ;
	  if ($f =~m/a-zA-Z/) { undef $f ;}
	  while (defined $f && defined $s) {
	      if (pl_add_to_sa_tar2($sa, $f, $s)) {$err = "PRIMER_ERROR=Too many elements for tag $tag\n"; return ;}     
	      $f = shift @nums ;
	      $s = shift @nums ;
	    }};
    $dispatch{'EXCLUDED_REGION'} = 
sub ($) { my $v = shift;
	  my @nums = split /[, ]/, $v ;
	  my $f = shift @nums ;
	  my $s = shift @nums ;
	  while (defined $f) {
	      if ($s eq "" || $f eq "") { return ; }
	      if (pl_add_to_sa_excl2($sa, $f, $s)) {$err = "PRIMER_ERROR=Too many elements for tag $tag\n"; return ;}     
	      $f = shift @nums ;
	      $s = shift @nums ;
	    }};
    $dispatch{'PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION'} = 
sub ($) { my $v = shift;
	  my @nums = split /[, ]/, $v ;
	  my $f = shift @nums ;
	  my $s = shift @nums ;
	  while (defined $f && defined $s) {
	      if (pl_add_to_sa_excl_internal2($sa, $f, $s)) {$err =  "PRIMER_ERROR=Too many elements for tag $tag\n"; return ;}     
	      $f = shift @nums ;
	      $s = shift @nums ;
	    }};
   $dispatch{'INCLUDED_REGION'} = 
sub ($) { my $v = shift;
	  my @nums = split /[, -]/, $v ;
	  my $f = shift @nums ;
	  my $s = shift @nums ;
	  if ($f eq "" || $s eq "") { return ; }
	  pl_set_sa_incl_s($sa, $f);     
	  pl_set_sa_incl_l($sa, $s)} ;     
   $dispatch{'PRIMER_START_CODON_POSITION'} = 
sub ($) { my $v = shift;  my $i = pl_set_sa_start_codon_pos($sa, $v) };     

# GLOBAL SETTINGS

   $dispatch{'PRIMER_PRODUCT_SIZE_RANGE'} = 
sub ($) { my $v = shift;
          pl_empty_gs_product_size_range($gs);     
	  my @nums = split /[-\" ]/, $v ;
	  my $f = "";
	  my $s = "";
	  while (1) {
	  $f = shift @nums ;
	  if (defined $f && $f eq "") {$f = shift @nums};
	  $s = shift @nums ;
	  if (defined $s && $s eq "") {$s = shift @nums};
	  if (!defined $f && !defined $s) { return ; }
	  if ($f ne "" && $s ne "")
	      {if (pl_add_to_gs_product_size_range($gs, $f, $s)) {$err = "PRIMER_ERROR=Too many elements for tag $tag\n"; return;}} ;
	  }} ;

   $dispatch{'PRIMER_DEFAULT_PRODUCT'} = 
sub ($) { my $v = shift;
	  my @nums = split /[-\"]/, $v ;
	  my $f = shift @nums ;
          if ($f eq "") { $f = shift @nums; }

	  my $s = shift @nums ;
          if ($s eq "") { $s = shift @nums; }
	  pl_set_gs_prmin($gs, $f, 0);     
	  pl_set_gs_prmax($gs, $s, 0)} ;     

   $dispatch{'PRIMER_DEFAULT_SIZE'} = 
sub ($) { my $v = shift;  my $i = pl_set_sa_primer_opt_size($gs, $v) };
   $dispatch{'PRIMER_OPT_SIZE'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_opt_size($gs, $v) };          
   $dispatch{'PRIMER_MIN_SIZE'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_min_size($gs, $v) };          
   $dispatch{'PRIMER_MAX_SIZE'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_max_size($gs, $v) };          
   $dispatch{'PRIMER_MAX_POLY_X'} =
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_max_poly_x($gs, $v) };          
   $dispatch{'PRIMER_OPT_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_opt_tm($gs, $v) };          
   $dispatch{'PRIMER_OPT_GC_PERCENT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_opt_gc_percent($gs, $v) };          
   $dispatch{'PRIMER_MIN_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_min_tm($gs, $v) };          
   $dispatch{'PRIMER_MAX_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_max_tm($gs, $v) };          
   $dispatch{'PRIMER_MAX_DIFF_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_max_diff_tm($gs, $v) };          
   $dispatch{'PRIMER_TM_SANTALUCIA'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_tm_santalucia($gs, $v) };          
   $dispatch{'PRIMER_SALT_CORRECTIONS'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_salt_corrections($gs, $v) };          
   $dispatch{'PRIMER_MIN_GC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_min_gc($gs, $v) };          
   $dispatch{'PRIMER_MAX_GC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_max_gc($gs, $v) };          
   $dispatch{'PRIMER_SALT_CONC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_salt_conc($gs, $v) };          
   $dispatch{'PRIMER_DIVALENT_CONC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_divalent_conc($gs, $v) };          
   $dispatch{'PRIMER_DNTP_CONC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_dntp_conc($gs, $v) };          
   $dispatch{'PRIMER_DNA_CONC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_dna_conc($gs, $v) };          
   $dispatch{'PRIMER_NUM_NS_ACCEPTED'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_num_ns_accepted($gs, $v) };          
   $dispatch{'PRIMER_PRODUCT_OPT_SIZE'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_product_opt_size($gs, $v) };          
   $dispatch{'PRIMER_SELF_ANY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_self_any($gs, $v) };          
   $dispatch{'PRIMER_SELF_END'} = 

# THIS NEEDS TO BE REMOVED:
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_self_end($gs, $v) };          
   $dispatch{'PRIMER_FILE_FLAG'} = 

sub ($) { my $v = shift;  my $i = pl_set_gs_primer_file_flag($gs, $v) };          
   $dispatch{'PRIMER_PICK_ANYWAY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pick_anyway($gs, $v) };     
   $dispatch{'PRIMER_GC_CLAMP'} =
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_gc_clamp($gs, $v) };          
   $dispatch{'PRIMER_EXPLAIN_FLAG'} = 
sub ($) { my $v = shift;  $explain_flag = 1 };          
   $dispatch{'PRIMER_LIBERAL_BASE'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_liberal_base($gs, $v) };          
   $dispatch{'PRIMER_FIRST_BASE_INDEX'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_first_base_index($gs, $v) };          
   $dispatch{'PRIMER_NUM_RETURN'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_num_return($gs, $v) };          
   $dispatch{'PRIMER_MIN_QUALITY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_min_quality($gs, $v) };          
   $dispatch{'PRIMER_MIN_END_QUALITY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_min_end_quality($gs, $v) };          
   $dispatch{'PRIMER_QUALITY_RANGE_MIN'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_quality_range_min($gs, $v) };          
   $dispatch{'PRIMER_QUALITY_RANGE_MAX'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_quality_range_max($gs, $v) };          
   $dispatch{'PRIMER_PRODUCT_MAX_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_product_max_tm($gs, $v) };          
   $dispatch{'PRIMER_PRODUCT_MIN_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_product_min_tm($gs, $v) };          
   $dispatch{'PRIMER_PRODUCT_OPT_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_product_opt_tm($gs, $v) };          
   $dispatch{'PRIMER_TASK'} = 
sub ($) { my $v = shift; my $i = pl_set_gs_primer_task($gs, $v); };
   $dispatch{'PRIMER_PICK_RIGHT_PRIMER'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pick_right_primer($gs, $v) };          
   $dispatch{'PRIMER_PICK_INTERNAL_OLIGO'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pick_internal_oligo($gs, $v) };          
   $dispatch{'PRIMER_PICK_LEFT_PRIMER'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pick_left_primer($gs, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_OPT_SIZE'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_opt_size($gs, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_MAX_SIZE'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_max_size($gs, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_MIN_SIZE'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_min_size($gs, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_MAX_POLY_X'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_max_poly_x($gs, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_OPT_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_opt_tm($gs, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_opt_gc_percent($gs, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_MAX_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_max_tm($gs, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_MIN_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_min_tm($gs, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_MIN_GC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_min_gc($gs, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_MAX_GC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_max_gc($gs, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_SALT_CONC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_salt_conc($gs, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_DIVALENT_CONC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_divalent_conc($gs, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_DNTP_CONC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_dntp_conc($gs, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_DNA_CONC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_dna_conc($gs, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_NUM_NS'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_num_ns($gs, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_MIN_QUALITY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_min_quality($gs, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_SELF_ANY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_self_any($gs, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_SELF_END'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_self_end($gs, $v) };          
   $dispatch{'PRIMER_MAX_LIBRARY_MISPRIMING'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_max_mispriming($gs, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_MAX_MISHYB'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_max_mishyb($gs, $v) };          
   $dispatch{'PRIMER_PAIR_MAX_LIBRARY_MISPRIMING'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_max_mispriming($gs, $v) };          
   $dispatch{'PRIMER_MAX_TEMPLATE_MISPRIMING'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_max_template_mispriming($gs, $v) };          
   $dispatch{'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_max_template_mispriming($gs, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_MAX_TEMPLATE_MISHYB'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_max_template_mishyb($gs, $v) };          
   $dispatch{'PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_ambiguity_codes_consensus($gs, $v) };          
   $dispatch{'PRIMER_INSIDE_PENALTY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_inside_penalty($gs, $v) };          
   $dispatch{'PRIMER_OUTSIDE_PENALTY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_outside_penalty($gs, $v) };          
   $dispatch{'PRIMER_MISPRIMING_LIBRARY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_mispriming_library($gs, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_MISHYB_LIBRARY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_mishyb_library($gs, $v) };          
   $dispatch{'PRIMER_MAX_END_STABILITY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_max_end_stability($gs, $v) };          
   $dispatch{'PRIMER_LOWERCASE_MASKING'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_lowercase_masking($gs, $v) };          
   $dispatch{'PRIMER_WT_TM_GT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_tm_gt($gs, $v) };          
   $dispatch{'PRIMER_WT_TM_LT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_tm_lt($gs, $v) };          
   $dispatch{'PRIMER_WT_GC_PERCENT_GT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_gc_percent_gt($gs, $v) };          
   $dispatch{'PRIMER_WT_GC_PERCENT_LT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_gc_percent_lt($gs, $v) };          
   $dispatch{'PRIMER_WT_SIZE_LT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_size_lt($gs, $v) };          
   $dispatch{'PRIMER_WT_SIZE_GT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_size_gt($gs, $v) };          
   $dispatch{'PRIMER_WT_COMPL_ANY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_compl_any($gs, $v) };          
   $dispatch{'PRIMER_WT_COMPL_END'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_compl_end($gs, $v) };          
   $dispatch{'PRIMER_WT_NUM_NS'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_num_ns($gs, $v) };          
   $dispatch{'PRIMER_WT_LIBRARY_MISPRIMING'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_rep_sim($gs, $v) };          
   $dispatch{'PRIMER_WT_SEQ_QUAL'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_seq_qual($gs, $v) };          
   $dispatch{'PRIMER_WT_END_QUAL'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_end_qual($gs, $v) };          
   $dispatch{'PRIMER_WT_POS_PENALTY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_pos_penalty($gs, $v) };          
   $dispatch{'PRIMER_WT_END_STABILITY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_end_stability($gs, $v) };          
   $dispatch{'PRIMER_WT_TEMPLATE_MISPRIMING'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_template_mispriming($gs, $v) };          
   $dispatch{'PRIMER_IO_WT_TM_GT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_tm_gt($gs, $v) };          
   $dispatch{'PRIMER_IO_WT_TM_LT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_tm_lt($gs, $v) };          
   $dispatch{'PRIMER_IO_WT_GC_PERCENT_GT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_gc_percent_gt($gs, $v) };          
   $dispatch{'PRIMER_IO_WT_GC_PERCENT_LT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_gc_percent_lt($gs, $v) };          
   $dispatch{'PRIMER_IO_WT_SIZE_LT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_size_lt($gs, $v) };          
   $dispatch{'PRIMER_IO_WT_SIZE_GT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_size_gt($gs, $v) };          
   $dispatch{'PRIMER_IO_WT_COMPL_ANY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_compl_any($gs, $v) };          
   $dispatch{'PRIMER_IO_WT_COMPL_END'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_compl_end($gs, $v) };          
   $dispatch{'PRIMER_IO_WT_NUM_NS'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_num_ns($gs, $v) };          
   $dispatch{'PRIMER_IO_WT_LIBRARY_MISPRIMING'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_rep_sim($gs, $v) };          
   $dispatch{'PRIMER_IO_WT_SEQ_QUAL'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_seq_qual($gs, $v) };          
   $dispatch{'PRIMER_IO_WT_END_QUAL'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_end_qual($gs, $v) };          
   $dispatch{'PRIMER_IO_WT_TEMPLATE_MISHYB'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_template_mishyb($gs, $v) };          
   $dispatch{'PRIMER_PAIR_WT_PR_PENALTY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_pr_penalty($gs, $v) };          
   $dispatch{'PRIMER_PAIR_WT_IO_PENALTY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_io_penalty($gs, $v) };          
   $dispatch{'PRIMER_PAIR_WT_DIFF_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_diff_tm($gs, $v) };          
   $dispatch{'PRIMER_PAIR_WT_COMPL_ANY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_compl_any($gs, $v) };          
   $dispatch{'PRIMER_PAIR_WT_COMPL_END'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_compl_end($gs, $v) };          
   $dispatch{'PRIMER_PAIR_WT_PRODUCT_TM_LT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_product_tm_lt($gs, $v) };          
   $dispatch{'PRIMER_PAIR_WT_PRODUCT_TM_GT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_product_tm_gt($gs, $v) };          
   $dispatch{'PRIMER_PAIR_WT_PRODUCT_SIZE_GT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_product_size_gt($gs, $v) };          
   $dispatch{'PRIMER_PAIR_WT_PRODUCT_SIZE_LT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_product_size_lt($gs, $v) };          
   $dispatch{'PRIMER_PAIR_WT_LIBRARY_MISPRIMING'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_rep_sim($gs, $v) };          
   $dispatch{'PRIMER_PAIR_WT_TEMPLATE_MISPRIMING'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_template_mispriming($gs, $v) };
$dispatch{'COMMENT'} = sub ($) {};
$dispatch{'PRIMER_COMMENT'} = sub ($) {};
$dispatch{'PRIMER_SHOW_OLIGO_PROBLEMS'} = sub { $show_oligo_problems = 1; };
}
