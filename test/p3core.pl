#!/usr/bin/perl

use warnings 'all';
use strict;

use Carp;
use blib '../../../p3perl/trunk';
use Primer3 ':all';

sub set_setters();
sub main();

    our %dispatch;
    our $sa ;
    our $tag;
    our $gs;

main();

sub main() {

    select STDOUT; $| = 1;  # Not sure we need this ....
    set_setters();
    $gs = pl_create_global_settings();

    $/ = "\n=\n";
    while (1) {
	my $rec = <>;
	last if !defined $rec;
	chomp $rec;
	$sa = pl_create_seq_arg();
	my %rec;
	my @rec = split /\n/, $rec;
	for my $line (@rec) {
	    print "$line\n";
	    if ($line !~ /(^\w+\s*)=(.*)/) {
		confess "Bad line $line\n";
	    }
	    my ($tag, $value) = ($1, $2);
	    if ($dispatch{$tag})  {
		&{$dispatch{$tag}}($value);
	    } else { confess "no call for $tag" }
	}
	my $retval = pl_choose_primers($gs, $sa);
	pl_boulder_print($gs, $sa, $retval) ;  # boulder_print generates the final '='
	pl_destroy_seq_args($sa);
    }
}

sub set_setters() {
     $dispatch{'SEQUENCE'} = 
sub ($) { my $v = shift; pl_set_sa_sequence($sa, $v) }; 
     $dispatch{'PRIMER_SEQUENCE_QUALITY'} = 
sub ($) { my $v = shift;
	  my @nums = split /[\b]/, $v ;
	  my $n = shift @nums ;
	  my $nq = pl_sa_get_n_quality($sa);     
	  while(defined $n) {
	      Mytest2:pl_set_sa_quality($sa, $nq, $n) ;
	      $nq++;
	      $n = shift @nums ;
	    } 
          pl_set_sa_n_quality($sa, $nq)} ;
    $dispatch{'PRIMER_SEQUENCE_ID'} = 
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
	      pl_add_to_sa_tar2($sa, $f, $s);     
	      $f = shift @nums ;
	      $s = shift @nums ;
	    }};
    $dispatch{'EXCLUDED_REGION'} = 
sub ($) { my $v = shift;
	  my @nums = split /[, ]/, $v ;
	  my $f = shift @nums ;
	  my $s = shift @nums ;
	  while (defined $f) {
	      pl_add_to_sa_excl2($sa, $f, $s);     
	      $f = shift @nums ;
	      $s = shift @nums ;
	    }};
    $dispatch{'PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION'} = 
sub ($) { my $v = shift;
	  my @nums = split /[, ]/, $v ;
	  my $f = shift @nums ;
	  my $s = shift @nums ;
	  while (defined $f) {
	      pl_add_to_sa_excl_internal2($sa, $f, $s);     
	      $f = shift @nums ;
	      $s = shift @nums ;
	    }};
   $dispatch{'INCLUDED_REGION'} = 
sub ($) { my $v = shift;
	  my @nums = split /[, ]/, $v ;
	  my $f = shift @nums ;
	  my $s = shift @nums ;
	  pl_set_sa_incl_s($sa, $f);     
	  pl_set_sa_incl_l($sa, $s)} ;     
   $dispatch{'PRIMER_START_CODON_POSITION'} = 
sub ($) { my $v = shift;  my $i = pl_set_sa_start_codon_pos($sa, $v) };     
   $dispatch{'PRIMER_PRODUCT_SIZE_RANGE'} = 
sub ($) { my $v = shift;
	  my @nums = split /[-" ]/, $v ;
	  my $f = shift @nums ;
          if ($f eq "") { $f = shift @nums; }

	  my $s = shift @nums ;
          if ($s eq "") { $s = shift @nums; }
          
          pl_empty_gs_product_size_range($gs);     
	  pl_add_to_gs_product_size_range($gs, $f, $s) };     

   $dispatch{'PRIMER_DEFAULT_PRODUCT'} = 
sub ($) { my $v = shift;
	  my @nums = split /[-"]/, $v ;
	  my $f = shift @nums ;
          if ($f eq "") { $f = shift @nums; }

	  my $s = shift @nums ;
          if ($s eq "") { $s = shift @nums; }

	  pl_set_gs_prmin($sa, $f, 1);     
	  pl_set_gs_prmax($sa, $s, 1)} ;     
   $dispatch{'PRIMER_DEFAULT_SIZE'} = 
sub ($) { my $v = shift;  my $i = pl_set_sa_primer_opt_size($sa, $v) };          
   $dispatch{'PRIMER_OPT_SIZE'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_opt_size($sa, $v) };          
   $dispatch{'PRIMER_MIN_SIZE'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_min_size($sa, $v) };          
   $dispatch{'PRIMER_MAX_SIZE'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_max_size($sa, $v) };          
   $dispatch{'PRIMER_MAX_POLY_X'} =
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_max_poly_x($sa, $v) };          
   $dispatch{'PRIMER_OPT_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_opt_tm($sa, $v) };          
   $dispatch{'PRIMER_OPT_GC_PERCENT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_opt_gc_percent($sa, $v) };          
   $dispatch{'PRIMER_MIN_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_min_tm($sa, $v) };          
   $dispatch{'PRIMER_MAX_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_max_tm($sa, $v) };          
   $dispatch{'PRIMER_MAX_TM_DIFF_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_max_tm_diff_tm($sa, $v) };          
   $dispatch{'PRIMER_TM_SANTALUCIA'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_tm_santalucia($sa, $v) };          
   $dispatch{'PRIMER_SALT_CORRECTIONS'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_salt_corrections($sa, $v) };          
   $dispatch{'PRIMER_MIN_GC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_min_gc($sa, $v) };          
   $dispatch{'PRIMER_MAX_GC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_max_gc($sa, $v) };          
   $dispatch{'PRIMER_SALT_CONC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_salt_conc($sa, $v) };          
   $dispatch{'PRIMER_DIVALENT_CONC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_divalent_conc($sa, $v) };          
   $dispatch{'PRIMER_DNTP_CONC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_dntp_conc($sa, $v) };          
   $dispatch{'PRIMER_DNA_CONC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_dna_conc($sa, $v) };          
   $dispatch{'PRIMER_NUM_NS_ACCEPTED'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_num_ns_accepted($sa, $v) };          
   $dispatch{'PRIMER_PRODUCT_OPT_SIZE'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_product_opt_size($sa, $v) };          
   $dispatch{'PRIMER_SELF_ANY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_self_any($sa, $v) };          
   $dispatch{'PRIMER_SELF_END'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_self_end($sa, $v) };          
   $dispatch{'PRIMER_FILE_FLAG'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_file_flag($sa, $v) };          
   $dispatch{'PRIMER_PICK_ANYWAY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pick_anyway($sa, $v) };     
   $dispatch{'PRIMER_GC_CLAMP'} =
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_gc_clamp($sa, $v) };          
   $dispatch{'PRIMER_EXPLAIN_FLAG'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_explain_flag($sa, $v) };          
   $dispatch{'PRIMER_LIBERAL_BASE'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_liberal_base($sa, $v) };          
   $dispatch{'PRIMER_FIRST_BASE_INDEX'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_first_base_index($sa, $v) };          
   $dispatch{'PRIMER_NUM_RETURN'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_num_return($sa, $v) };          
   $dispatch{'PRIMER_MIN_QUALITY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_min_quality($sa, $v) };          
   $dispatch{'PRIMER_MIN_END_QUALITY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_min_end_quality($sa, $v) };          
   $dispatch{'PRIMER_QUALITY_RANGE_MIN'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_quality_range_min($sa, $v) };          
   $dispatch{'PRIMER_QUALITY_RANGE_MAX'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_quality_range_max($sa, $v) };          
   $dispatch{'PRIMER_PRODUCT_MAX_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_product_max_tm($sa, $v) };          
   $dispatch{'PRIMER_PRODUCT_MIN_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_product_min_tm($sa, $v) };          
   $dispatch{'PRIMER_PRODUCT_OPT_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_product_opt_tm($sa, $v) };          
   $dispatch{'PRIMER_TASK'} = 
# WARNING this needs to be done for the other possible tasks
sub ($) { my $v = shift;  
	  my $i = pl_set_gs_primer_task($sa, $v);
	      if ($v =~ /^\s*PICK_PCR_PRIMERS\s*$/) {
		   pl_set_gs_primer_pick_left($sa, $v);
		   pl_set_gs_primer_pick_right($sa, $v);
	      }
 };
   $dispatch{'PRIMER_PICK_RIGHT_PRIMER'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pick_right($sa, $v) };          
   $dispatch{'PRIMER_PICK_INTERNAL_OLIGO'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pick_internal_oligo($sa, $v) };          
   $dispatch{'PRIMER_PICK_LEFT_PRIMER'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pick_left_primer($sa, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_OPT_SIZE'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_opt_size($sa, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_MAX_SIZE'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_max_size($sa, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_MIN_SIZE'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_min_size($sa, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_MAX_POLY_X'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_max_poly_x($sa, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_OPT_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_opt_tm($sa, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_opt_gc_percent($sa, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_MAX_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_max_tm($sa, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_MIN_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_min_tm($sa, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_MIN_GC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_min_gc($sa, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_MAX_GC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_max_gc($sa, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_SALT_CONC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_salt_conc($sa, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_DIVALENT_CONC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_divalent_conc($sa, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_DNTP_CONC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_dntp_conc($sa, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_DNA_CONC'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_dna_conc($sa, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_NUM_NS'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_num_ns($sa, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_MIN_QUALITY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_min_quality($sa, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_SELF_ANY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_self_any($sa, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_SELF_END'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_self_end($sa, $v) };          
   $dispatch{'PRIMER_MAX_MISPRIMING'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_max_mispriming($sa, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_MAX_MISHYB'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_max_mishyb($sa, $v) };          
   $dispatch{'PRIMER_PAIR_MAX_MISPRIMING'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_max_mispriming($sa, $v) };          
   $dispatch{'PRIMER_MAX_TEMPLATE_MISPRIMING'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_max_template_mispriming($sa, $v) };          
   $dispatch{'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_max_template_mispriming($sa, $v) };          
   $dispatch{'PRIMER_INTERNAL_OLIGO_MAX_TEMPLATE_MISHYB'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_internal_oligo_max_template_mishyb($sa, $v) };          
   $dispatch{'PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_ambiguity_codes_consensus($sa, $v) };          
   $dispatch{'PRIMER_INSIDE_PENALTY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_inside_penalty($sa, $v) };          
   $dispatch{'PRIMER_OUTSIDE_PENALTY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_outside_penalty($sa, $v) };          
   $dispatch{'PRIMER_MISPRIMING_LIBRARY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_mispriming_library($sa, $v) };          
   $dispatch{'PRIMER_MISHYB_LIBRARY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_mishyb_library($sa, $v) };          
   $dispatch{'PRIMER_MAX_END_STABILITY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_max_end_stability($sa, $v) };          
   $dispatch{'PRIMER_LOWERCASE_MASKING'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_lowercase_masking($sa, $v) };          
   $dispatch{'PRIMER_WT_TM_GT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_tm_gt($sa, $v) };          
   $dispatch{'PRIMER_WT_TM_LT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_tm_lt($sa, $v) };          
   $dispatch{'PRIMER_WT_GC_PERCENT_GT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_gc_percent_gt($sa, $v) };          
   $dispatch{'PRIMER_WT_GC_PERCENT_LT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_gc_percent_lt($sa, $v) };          
   $dispatch{'PRIMER_WT_SIZE_LT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_size_lt($sa, $v) };          
   $dispatch{'PRIMER_WT_SIZE_GT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_size_gt($sa, $v) };          
   $dispatch{'PRIMER_WT_COMPL_ANY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_compl_any($sa, $v) };          
   $dispatch{'PRIMER_WT_COMPL_END'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_compl_end($sa, $v) };          
   $dispatch{'PRIMER_WT_NUM_NS'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_num_ns($sa, $v) };          
   $dispatch{'PRIMER_WT_REP_SIM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_rep_sim($sa, $v) };          
   $dispatch{'PRIMER_WT_SEQ_QUAL'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_seq_qual($sa, $v) };          
   $dispatch{'PRIMER_WT_END_QUAL'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_end_qual($sa, $v) };          
   $dispatch{'PRIMER_WT_POS_PENALTY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_pos_penalty($sa, $v) };          
   $dispatch{'PRIMER_WT_END_STABILITY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_end_stability($sa, $v) };          
   $dispatch{'PRIMER_WT_TEMPLATE_MISPRIMING'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_wt_template_mispriming($sa, $v) };          
   $dispatch{'PRIMER_IO_WT_TM_GT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_tm_gt($sa, $v) };          
   $dispatch{'PRIMER_IO_WT_TM_LT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_tm_lt($sa, $v) };          
   $dispatch{'PRIMER_IO_WT_GC_PERCENT_GT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_gc_percent_gt($sa, $v) };          
   $dispatch{'PRIMER_IO_WT_GC_PERCENT_LT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_gc_percent_lt($sa, $v) };          
   $dispatch{'PRIMER_IO_WT_SIZE_LT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_size_lt($sa, $v) };          
   $dispatch{'PRIMER_IO_WT_SIZE_GT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_size_gt($sa, $v) };          
   $dispatch{'PRIMER_IO_WT_COMPL_ANY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_compl_any($sa, $v) };          
   $dispatch{'PRIMER_IO_WT_COMPL_END'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_compl_end($sa, $v) };          
   $dispatch{'PRIMER_IO_WT_NUM_NS'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_num_ns($sa, $v) };          
   $dispatch{'PRIMER_IO_WT_REP_SIM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_rep_sim($sa, $v) };          
   $dispatch{'PRIMER_IO_WT_SEQ_QUAL'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_seq_qual($sa, $v) };          
   $dispatch{'PRIMER_IO_WT_END_QUAL'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_end_qual($sa, $v) };          
   $dispatch{'PRIMER_IO_WT_TEMPLATE_MISHYB'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_io_wt_template_mishyb($sa, $v) };          
   $dispatch{'PRIMER_PAIR_WT_PR_PENALTY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_pr_penalty($sa, $v) };          
   $dispatch{'PRIMER_PAIR_WT_IO_PENALTY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_io_penalty($sa, $v) };          
   $dispatch{'PRIMER_PAIR_WT_DIFF_TM'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_diff_tm($sa, $v) };          
   $dispatch{'PRIMER_PAIR_WT_COMPL_ANY'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_compl_any($sa, $v) };          
   $dispatch{'PRIMER_PAIR_WT_COMPL_END'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_compl_end($sa, $v) };          
   $dispatch{'PRIMER_PAIR_WT_PRODUCT_TM_LT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_product_tm_lt($sa, $v) };          
   $dispatch{'PRIMER_PAIR_WT_PRODUCT_TM_GT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_product_tm_gt($sa, $v) };          
   $dispatch{'PRIMER_PAIR_WT_PRODUCT_SIZE_GT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_product_size_gt($sa, $v) };          
   $dispatch{'PRIMER_PAIR_WT_PRODUCT_SIZE_LT'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_product_size_lt($sa, $v) };          
   $dispatch{'PRIMER_PAIR_WT_REP_SIM,'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_rep_sim($sa, $v) };          
   $dispatch{'PRIMER_PAIR_WT_TEMPLATE_MISPRIMING'} = 
sub ($) { my $v = shift;  my $i = pl_set_gs_primer_pair_wt_template_mispriming($sa, $v) };
}
