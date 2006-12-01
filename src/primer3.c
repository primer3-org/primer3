/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006
Whitehead Institute for Biomedical Research, Steve Rozen
(http://jura.wi.mit.edu/rozen), and Helen Skaletsky
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

   * Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the
distribution.
   * Neither the names of the copyright holders nor contributors may
be used to endorse or promote products derived from this software
without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <string.h>
#include <float.h>
#include "primer3.h"
#include "primer3_release.h"

#define MACRO_VALUE_AS_STRING(A) MACRO_STRING(A)

static int    check_intervals(const char *, const int,
			      interval_array_t, const int, seq_args *);
static char   dna_to_upper(char *, int);
static char   *strstr_nocase(char *, char *);


/* Default parameter values.  */
#define OPT_SIZE            20
#define MIN_SIZE             18
#define MAX_SIZE             27

#define OPT_TM             60.0
#define MIN_TM             57.0
#define MAX_TM             63.0
#define MAX_DIFF_TM       100.0
#define DEFAULT_OPT_GC_PERCENT PR_UNDEFINED_INT_OPT
#define MIN_GC             20.0
#define MAX_GC             80.0
#define SALT_CONC          50.0
#define DNA_CONC           50.0
#define NUM_NS_ACCEPTED       0
#define MAX_POLY_X            5
#define SELF_ANY            800
#define SELF_END            300
#define PAIR_COMPL_ANY      800
#define PAIR_COMPL_END      300
#define FILE_FLAG             0
#define EXPLAIN_FLAG          0
#define GC_CLAMP              0
#define LIBERAL_BASE          0
#define PICK_INTERNAL_OLIGO   0
#define PRIMER_TASK           0
#define INTERNAL_OLIGO_OPT_SIZE   20
#define INTERNAL_OLIGO_MIN_SIZE   18
#define INTERNAL_OLIGO_MAX_SIZE   27
#define INTERNAL_OLIGO_OPT_TM     60.0
#define INTERNAL_OLIGO_MIN_TM     57.0
#define INTERNAL_OLIGO_MAX_TM     63.0
#define INTERNAL_OLIGO_MIN_GC     20.0
#define INTERNAL_OLIGO_MAX_GC     80.0
#define INTERNAL_OLIGO_SALT_CONC         50.0
#define INTERNAL_OLIGO_DNA_CONC          50.0
#define INTERNAL_OLIGO_NUM_NS               0
#define INTERNAL_OLIGO_MAX_POLY_X           5 
#define INTERNAL_OLIGO_SELF_ANY          1200
#define INTERNAL_OLIGO_SELF_END          1200
#define INTERNAL_OLIGO_REPEAT_SIMILARITY 1200
#define REPEAT_SIMILARITY                1200
#define PAIR_REPEAT_SIMILARITY           2400
#define FIRST_BASE_INDEX            0
#define NUM_RETURN                  5
#define MIN_QUALITY                 0
#define QUALITY_RANGE_MIN           0
#define QUALITY_RANGE_MAX         100
#define DEFAULT_MAX_END_STABILITY    100.0
#define PRIMER_PRODUCT_OPT_SIZE      PR_UNDEFINED_INT_OPT
#define PRIMER_PRODUCT_OPT_TM        PR_UNDEFINED_DBL_OPT
#define MAX_TEMPLATE_MISPRIMING      PR_UNDEFINED_ALIGN_OPT
#define PAIR_MAX_TEMPLATE_MISPRIMING PR_UNDEFINED_ALIGN_OPT
#define IO_MAX_TEMPLATE_MISHYB       PR_UNDEFINED_ALIGN_OPT

#define LIB_AMBIGUITY_CODES_CONSENSUS 1
/*  For backward compatibility. It turns out that
    this _not_ what one normally wants, since many
    libraries contain strings of N, which then match
    every oligo (very bad).
*/

/* Weights for objective functions for oligos and pairs. */
#define PRIMER_WT_TM_GT          1
#define PRIMER_WT_TM_LT          1
#define PRIMER_WT_SIZE_LT        1
#define PRIMER_WT_SIZE_GT        1
#define PRIMER_WT_GC_PERCENT_LT  0
#define PRIMER_WT_GC_PERCENT_GT  0
#define PRIMER_WT_COMPL_ANY      0
#define PRIMER_WT_COMPL_END      0
#define PRIMER_WT_NUM_NS         0
#define PRIMER_WT_REP_SIM        0
#define PRIMER_WT_SEQ_QUAL       0
#define PRIMER_WT_END_QUAL       0
#define PRIMER_WT_POS_PENALTY    1
#define PRIMER_WT_END_STABILITY  0

#define IO_WT_TM_GT          1
#define IO_WT_TM_LT          1
#define IO_WT_SIZE_LT        1
#define IO_WT_SIZE_GT        1
#define IO_WT_GC_PERCENT_LT  0
#define IO_WT_GC_PERCENT_GT  0
#define IO_WT_COMPL_ANY      0
#define IO_WT_COMPL_END      0
#define IO_WT_NUM_NS         0
#define IO_WT_REP_SIM        0
#define IO_WT_SEQ_QUAL       0
#define IO_WT_END_QUAL       0

#define PAIR_WT_PRIMER_PENALTY      1
#define PAIR_WT_IO_PENALTY          0
#define PAIR_WT_DIFF_TM             0
#define PAIR_WT_COMPL_ANY           0
#define PAIR_WT_COMPL_END           0
#define PAIR_WT_REP_SIM             0
#define PAIR_WT_PRODUCT_TM_LT       0
#define PAIR_WT_PRODUCT_TM_GT       0
#define PAIR_WT_PRODUCT_SIZE_LT     0
#define PAIR_WT_PRODUCT_SIZE_GT     0

void
pr_set_default_global_args(a)
    primer_args *a;
{
    memset(a, 0, sizeof(*a));  
    a->primer_opt_size  = OPT_SIZE;
    a->primer_min_size  = MIN_SIZE;
    a->primer_max_size  = MAX_SIZE;
    a->opt_tm           = OPT_TM;
    a->min_tm           = MIN_TM;
    a->max_tm           = MAX_TM;
    a->max_diff_tm      = MAX_DIFF_TM;
    a->min_gc           = MIN_GC;
    a->opt_gc_content   = DEFAULT_OPT_GC_PERCENT;
    a->max_gc           = MAX_GC;
    a->salt_conc        = SALT_CONC;
    a->dna_conc         = DNA_CONC;
    a->num_ns_accepted  = NUM_NS_ACCEPTED;
    a->self_any         = SELF_ANY;
    a->self_end         = SELF_END;
    a->pair_compl_any   = PAIR_COMPL_ANY;
    a->pair_compl_end   = PAIR_COMPL_END;
    a->file_flag        = FILE_FLAG;
    a->explain_flag     = EXPLAIN_FLAG;
    a->gc_clamp         = GC_CLAMP;
    a->max_poly_x       = MAX_POLY_X;
    a->liberal_base      = LIBERAL_BASE;
    a->primer_task       = PRIMER_TASK;
    a->first_base_index  = FIRST_BASE_INDEX;
    a->num_return        = NUM_RETURN;
    a->pr_min[0]         = 100;
    a->pr_max[0]         = 300;
    a->num_intervals     = 1;
    a->repeat_compl      = REPEAT_SIMILARITY;
    a->pair_repeat_compl = PAIR_REPEAT_SIMILARITY;
    a->min_quality       = MIN_QUALITY;
    a->min_end_quality   = MIN_QUALITY;
    a->quality_range_min = QUALITY_RANGE_MIN;
    a->quality_range_max = QUALITY_RANGE_MAX;
    a->outside_penalty   = PR_DEFAULT_OUTSIDE_PENALTY;
    a->inside_penalty    = PR_DEFAULT_INSIDE_PENALTY;
    a->max_end_stability = DEFAULT_MAX_END_STABILITY;
    a->product_max_tm    = PR_DEFAULT_PRODUCT_MAX_TM;
    a->product_min_tm    = PR_DEFAULT_PRODUCT_MIN_TM;
    a->product_opt_tm    = PRIMER_PRODUCT_OPT_TM;
    a->product_opt_size  = PRIMER_PRODUCT_OPT_SIZE;
    a->max_template_mispriming
                          = MAX_TEMPLATE_MISPRIMING;
    a->pair_max_template_mispriming
                          = PAIR_MAX_TEMPLATE_MISPRIMING;

    a->io_primer_opt_size = INTERNAL_OLIGO_OPT_SIZE;
    a->io_primer_min_size = INTERNAL_OLIGO_MIN_SIZE;
    a->io_primer_max_size = INTERNAL_OLIGO_MAX_SIZE;
    a->io_opt_tm          = INTERNAL_OLIGO_OPT_TM;
    a->io_min_tm          = INTERNAL_OLIGO_MIN_TM;
    a->io_max_tm          = INTERNAL_OLIGO_MAX_TM;
    a->io_min_gc          = INTERNAL_OLIGO_MIN_GC;
    a->io_opt_gc_content  = DEFAULT_OPT_GC_PERCENT;
    a->io_max_gc          = INTERNAL_OLIGO_MAX_GC;
    a->io_max_poly_x      = INTERNAL_OLIGO_MAX_POLY_X;
    a->io_salt_conc       = INTERNAL_OLIGO_SALT_CONC;
    a->io_dna_conc        = INTERNAL_OLIGO_DNA_CONC;
    a->io_num_ns_accepted = INTERNAL_OLIGO_NUM_NS;
    a->io_self_any        = INTERNAL_OLIGO_SELF_ANY;
    a->io_self_end        = INTERNAL_OLIGO_SELF_END;
    a->io_repeat_compl    = INTERNAL_OLIGO_REPEAT_SIMILARITY;
    a->io_min_quality     = MIN_QUALITY;
    a->io_min_end_quality = MIN_QUALITY;
    a->io_max_template_mishyb
                          = IO_MAX_TEMPLATE_MISHYB;

    a->primer_weights.temp_gt       = PRIMER_WT_TM_GT;
    a->primer_weights.temp_lt       = PRIMER_WT_TM_LT;
    a->primer_weights.length_gt     = PRIMER_WT_SIZE_GT;
    a->primer_weights.length_lt     = PRIMER_WT_SIZE_LT;
    a->primer_weights.gc_content_gt = PRIMER_WT_GC_PERCENT_GT;
    a->primer_weights.gc_content_lt = PRIMER_WT_GC_PERCENT_LT;
    a->primer_weights.compl_any     = PRIMER_WT_COMPL_ANY;
    a->primer_weights.compl_end     = PRIMER_WT_COMPL_END;
    a->primer_weights.num_ns        = PRIMER_WT_NUM_NS;
    a->primer_weights.repeat_sim    = PRIMER_WT_REP_SIM;
    a->primer_weights.seq_quality   = PRIMER_WT_SEQ_QUAL;
    a->primer_weights.end_quality   = PRIMER_WT_END_QUAL;
    a->primer_weights.pos_penalty   = PRIMER_WT_POS_PENALTY;
    a->primer_weights.end_stability = PRIMER_WT_END_STABILITY;

    a->io_weights.temp_gt     = IO_WT_TM_GT;
    a->io_weights.temp_lt     = IO_WT_TM_LT;
    a->io_weights.length_gt   = IO_WT_SIZE_GT;
    a->io_weights.length_lt   = IO_WT_SIZE_LT;
    a->io_weights.gc_content_gt = IO_WT_GC_PERCENT_GT;
    a->io_weights.gc_content_lt = IO_WT_GC_PERCENT_LT;
    a->io_weights.compl_any   = IO_WT_COMPL_ANY;
    a->io_weights.compl_end   = IO_WT_COMPL_END;
    a->io_weights.num_ns      = IO_WT_NUM_NS;
    a->io_weights.repeat_sim  = IO_WT_REP_SIM;
    a->io_weights.seq_quality = IO_WT_SEQ_QUAL;
    a->io_weights.end_quality = IO_WT_END_QUAL;

    a->pr_pair_weights.primer_quality  = PAIR_WT_PRIMER_PENALTY;
    a->pr_pair_weights.io_quality      = PAIR_WT_IO_PENALTY;
    a->pr_pair_weights.diff_tm         = PAIR_WT_DIFF_TM;
    a->pr_pair_weights.compl_any       = PAIR_WT_COMPL_ANY;
    a->pr_pair_weights.compl_end       = PAIR_WT_COMPL_END;
    a->pr_pair_weights.repeat_sim      = PAIR_WT_REP_SIM;
    a->pr_pair_weights.product_tm_lt   = PAIR_WT_PRODUCT_TM_LT;
    a->pr_pair_weights.product_tm_gt   = PAIR_WT_PRODUCT_TM_GT;
    a->pr_pair_weights.product_size_lt = PAIR_WT_PRODUCT_SIZE_LT;
    a->pr_pair_weights.product_size_gt = PAIR_WT_PRODUCT_SIZE_GT;

    /* a->short_match = 0; not used */

    a->lib_ambiguity_codes_consensus   = LIB_AMBIGUITY_CODES_CONSENSUS;
}

/*
 * Return 1 on error, 0 on success.  Set sa->trimmed_seq and possibly modify
 * sa->tar.  Upcase and check all bases in sa->trimmed_seq.
 */
int
_pr_data_control(pa, sa)
    primer_args *pa;
    seq_args *sa;
{
    static char s1[MAX_PRIMER_LENGTH+1];
    int i, pr_min;
    int seq_len = strlen(sa->sequence);
    char offending_char = '\0';

    if (pa->io_max_template_mishyb >= 0)
      pr_append_new_chunk(&pa->glob_err,
			  "PRIMER_INTERNAL_OLIGO_MAX_TEMPLATE_MISHYB is not supported");

    if (pa->primer_min_size < 1)
      pr_append_new_chunk(&pa->glob_err, "PRIMER_MIN_SIZE must be >= 1");

    if (pa->primer_max_size > MAX_PRIMER_LENGTH) {
      pr_append_new_chunk(&pa->glob_err,
			  "PRIMER_MAX_SIZE exceeds built-in maximum of ");
      pr_append(&pa->glob_err, MACRO_VALUE_AS_STRING(MAX_PRIMER_LENGTH));
      return 1;
    }

    if (pa->primer_opt_size > pa->primer_max_size) {
	pr_append_new_chunk(&pa->glob_err,
			    "PRIMER_{OPT,DEFAULT}_SIZE > PRIMER_MAX_SIZE");
	return 1;
    }

    if (pa->primer_opt_size < pa->primer_min_size) {
	pr_append_new_chunk(&pa->glob_err,
			    "PRIMER_{OPT,DEFAULT}_SIZE < PRIMER_MIN_SIZE");
	return 1;
    }

    if (pa->io_primer_max_size > MAX_PRIMER_LENGTH) {
	pr_append_new_chunk(&pa->glob_err,
		  "PRIMER_INTERNAL_OLIGO_MAX_SIZE exceeds built-in maximum");
        return 1;
    }

    if (pa->io_primer_opt_size > pa->io_primer_max_size) {
	pr_append_new_chunk(&pa->glob_err,
		  "PRIMER_INTERNAL_OLIGO_{OPT,DEFAULT}_SIZE > MAX_SIZE");
        return 1;
    }

    if (pa->io_primer_opt_size < pa->io_primer_min_size) {
	pr_append_new_chunk(&pa->glob_err,
		  "PRIMER_INTERNAL_OLIGO_{OPT,DEFAULT}_SIZE < MIN_SIZE");
        return 1;
    }

    if (pa->gc_clamp > pa->primer_min_size) {
	pr_append_new_chunk(&pa->glob_err,
			    "PRIMER_GC_CLAMP > PRIMER_MIN_SIZE");
	return 1;
    }

    if (NULL == sa->sequence_name && pa->file_flag) {
	pr_append_new_chunk(&sa->error,
			    "Need PRIMER_SEQUENCE_ID if PRIMER_FILE_FLAG != 0");
	return 1;
    }

    if (0 == pa->num_intervals) {
	pr_append_new_chunk(&pa->glob_err,
			    "Empty value for PRIMER_PRODUCT_SIZE_RANGE");
	return 1;
    }
    for (i = 0; i < pa->num_intervals; i++) {
	if (pa->pr_min[i] > pa->pr_max[i] || pa->pr_min[i] < 0) {
	    pr_append_new_chunk(&pa->glob_err,
				"Illegal element in PRIMER_PRODUCT_SIZE_RANGE");
	    return 1;
	}
    }

    pr_min = INT_MAX;
    for(i=0;i<pa->num_intervals;i++)
	if(pa->pr_min[i]<pr_min) pr_min=pa->pr_min[i];

    if (pa->primer_max_size > pr_min) {
	pr_append_new_chunk(&pa->glob_err,
			    "PRIMER_MAX_SIZE > min PRIMER_PRODUCT_SIZE_RANGE");
	return 1;
    }

    if ((pick_pcr_primers_and_hyb_probe == pa->primer_task 
	 || pick_hyb_probe_only == pa->primer_task)
	&& pa->io_primer_max_size > pr_min) {
	pr_append_new_chunk(&pa->glob_err,
		 "PRIMER_INTERNAL_OLIGO_MAX_SIZE > min PRIMER_PRODUCT_SIZE_RANGE");
        return 1;
    }

    if (pa->num_return < 1) {
	pr_append_new_chunk(&pa->glob_err,
			    "PRIMER_NUM_RETURN < 1");
        return 1;
    }

    if (sa->incl_l >= INT_MAX) {
	pr_append_new_chunk(&sa->error, "Value for INCLUDED_REGION too large");
	return 1;
    }

    if (sa->incl_s < 0 || sa->incl_l < 0 
	|| sa->incl_s + sa->incl_l > seq_len) {
	pr_append_new_chunk(&sa->error, "Illegal value for INCLUDED_REGION");
	return 1;
    }
    
    if (sa->incl_l < pr_min && pa->primer_task != pick_hyb_probe_only
	&& pa->primer_task != pick_left_only
	&& pa->primer_task != pick_right_only) {
	pr_append_new_chunk(&sa->error,
	   "INCLUDED_REGION length < min PRIMER_PRODUCT_SIZE_RANGE");
	return 1;
    }

    if (pa->max_end_stability < 0) {
        pr_append_new_chunk(&sa->error,
			    "PRIMER_MAX_END_STABILITY must be non-negative");
	return 1;
    }

    if (!PR_START_CODON_POS_IS_NULL(sa)) {
      if (!PR_POSITION_PENALTY_IS_NULL(pa)) {
	pr_append_new_chunk(&sa->error,
	   "Cannot accept both PRIMER_START_CODON_POSITION and non-default ");
	pr_append(&sa->error,
	   "arguments for PRIMER_INSIDE_PENALTY or PRIMER_OUTSIDE_PENALTY");
      }
      if (sa->start_codon_pos  > (sa->incl_s + sa->incl_l - 3))
	pr_append_new_chunk(&sa->error,
	   "Start codon position not contained in INCLUDED_REGION");
      else {
	if (sa->start_codon_pos >= 0
	    && ((sa->sequence[sa->start_codon_pos] != 'A'
		 && sa->sequence[sa->start_codon_pos] != 'a')
		|| (sa->sequence[sa->start_codon_pos + 1] != 'T'
		    && sa->sequence[sa->start_codon_pos + 1] != 't')
		|| (sa->sequence[sa->start_codon_pos + 2] != 'G'
		    && sa->sequence[sa->start_codon_pos + 2] != 'g')))
	  pr_append_new_chunk(&sa->error,
			      "No start codon at PRIMER_START_CODON_POSITION");
      }
    }

    sa->trimmed_seq = pr_safe_malloc(sa->incl_l + 1);
    _pr_substr(sa->sequence, sa->incl_s, sa->incl_l, sa->trimmed_seq);

    sa->upcased_seq = pr_safe_malloc(strlen(sa->sequence) + 1);
    strcpy(sa->upcased_seq, sa->sequence);
    if ((offending_char = dna_to_upper(sa->upcased_seq, 1))) {
      offending_char = '\0';
      /* NEW */
      /* TODO add warning or error (depending on liberal base)
         here. */
    }
    sa->upcased_seq_r = pr_safe_malloc(strlen(sa->sequence) + 1);
    _pr_reverse_complement(sa->upcased_seq, sa->upcased_seq_r);

    if (check_intervals("TARGET", sa->num_targets, sa->tar, seq_len, sa)
	== 1) return 1;
    sa->start_codon_pos -= sa->incl_s;

    if (check_intervals("EXCLUDED_REGION", sa->num_excl, sa->excl,
			seq_len, sa)
	== 1) return 1;

    if (check_intervals("PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION",
			sa->num_internal_excl, sa->excl_internal,
			seq_len, sa)
	== 1) return 1;

    if (NULL != sa->quality){
	if(pa->min_quality != 0 && pa->min_quality < pa->quality_range_min) {
	   pr_append_new_chunk(&pa->glob_err,
	       "PRIMER_MIN_QUALITY < PRIMER_QUALITY_RANGE_MIN");
           return 1;
        }
	if(pa->min_quality != 0 && pa->min_quality > pa->quality_range_max) {
	   pr_append_new_chunk(&pa->glob_err,
	       "PRIMER_MIN_QUALITY > PRIMER_QUALITY_RANGE_MAX");
           return 1;
        }
	if(pa->io_min_quality != 0 && pa->io_min_quality <pa->quality_range_min) {
	   pr_append_new_chunk(&pa->glob_err,
	    "PRIMER_INTERNAL_OLIGO_MIN_QUALITY < PRIMER_QUALITY_RANGE_MIN");
           return 1;
        }
	if(pa->io_min_quality != 0 && pa->io_min_quality > pa->quality_range_max) {
	   pr_append_new_chunk(&pa->glob_err,
	     "PRIMER_INTERNAL_OLIGO_MIN_QUALITY > PRIMER_QUALITY_RANGE_MAX");
           return 1;
        }
	for(i=0; i< seq_len; i++) {
	   if(sa->quality[i] < pa->quality_range_min ||
		 sa->quality[i] > pa->quality_range_max) {
             pr_append_new_chunk(&sa->error,
		"Sequence quality score out of range");
             return 1;
           }
        }
    }
    else if (pa->primer_weights.seq_quality || pa->io_weights.seq_quality) {
	 pr_append_new_chunk(&sa->error,
	      "Sequence quality is part of objective function but sequence quality is not defined");
         return 1;
    }

    if ((offending_char = dna_to_upper(sa->trimmed_seq, 0))) {
      if (pa->liberal_base) {
	pr_append_new_chunk(&sa->warning,
			    "Unrecognized base in input sequence");
      }
      else {
	pr_append_new_chunk(&sa->error, "Unrecognized base in input sequence");
	return 1;
      }
    }

    if(pa->opt_tm < pa->min_tm || pa->opt_tm > pa->max_tm) {
	 pr_append_new_chunk(&pa->glob_err,
			     "Optimum primer Tm lower than minimum or higher than maximum");
	 return 1;
    }
    if(pa->io_opt_tm < pa->io_min_tm || pa->io_opt_tm > pa->io_max_tm) {
	 pr_append_new_chunk(&pa->glob_err,
			     "Illegal values for PRIMER_INTERNAL_OLIGO_TM");
	 return 1;
    }
    if(pa->min_gc>pa->max_gc||pa->min_gc>100||pa->max_gc<0){
	 pr_append_new_chunk(&pa->glob_err,
			     "Illegal value for PRIMER_MAX_GC and PRIMER_MIN_GC");
	 return 1;
    }
    if(pa->io_min_gc>pa->io_max_gc||pa->io_min_gc>100||pa->io_max_gc<0){
	 pr_append_new_chunk(&pa->glob_err,
			     "Illegal value for PRIMER_INTERNAL_OLIGO_GC");
	 return 1;
    }
    if(pa->num_ns_accepted<0){
	 pr_append_new_chunk(&pa->glob_err,
			     "Illegal value for PRIMER_NUM_NS_ACCEPTED");
	 return 1;
    }
    if(pa->io_num_ns_accepted<0){ 
	 pr_append_new_chunk(&pa->glob_err,
			     "Illegal value for PRIMER_INTERNAL_OLIGO_NUM_NS");
	 return 1;
    }
    if(pa->self_any<0||pa->self_end<0
		       ||pa->pair_compl_any<0||pa->pair_compl_end<0){
         pr_append_new_chunk(&pa->glob_err,
	     "Illegal value for primer complementarity restrictions");
	 return 1;
    }
    if(pa->io_self_any<0||pa->io_self_end<0){
	  pr_append_new_chunk(&pa->glob_err,
	      "Illegal value for internal oligo complementarity restrictions");
	  return 1;
    }
    if(pa->salt_conc<=0||pa->dna_conc<=0){
	  pr_append_new_chunk(&pa->glob_err,
	      "Illegal value for primer salt or dna concentration");
	  return 1;
    }
    if(pa->io_salt_conc<=0||pa->io_dna_conc<=0){
	  pr_append_new_chunk(&pa->glob_err,
	      "Illegal value for internal oligo salt or dna concentration");
	  return 1;
    }
    if (!_PR_DEFAULT_POSITION_PENALTIES(pa) && sa->num_targets > 1) {
      pr_append_new_chunk(&sa->error,
			  "Non-default inside penalty or outside penalty ");
      pr_append(&sa->error,
		"is valid only when number of targets <= 1");
    }
    if (!_PR_DEFAULT_POSITION_PENALTIES(pa) && 0 == sa->num_targets) {
      pr_append_new_chunk(&sa->warning,
			  "Non-default inside penalty or outside penalty ");
      pr_append(&sa->warning,
		"has no effect when number of targets is 0");
    }
    if (pa->primer_task != pick_pcr_primers_and_hyb_probe 
	&& pa->primer_task != pick_hyb_probe_only
	&& sa->internal_input) {
      pr_append_new_chunk(&sa->error,
			  "Not specified to pick internal oligos");
      pr_append(&sa->error,
		" but a specific internal oligo is provided");
    }
    if (sa->internal_input) {
      if (strlen(sa->internal_input) > pa->io_primer_max_size)
	pr_append_new_chunk(&sa->error, "Specified internal oligo too long");

      if (strlen(sa->internal_input) < pa->io_primer_min_size)
	pr_append_new_chunk(&sa->error, "Specified internal oligo too short");

      if (!strstr_nocase(sa->sequence, sa->internal_input))
	pr_append_new_chunk(&sa->error,
			    "Specified internal oligo not in sequence");
      else if (!strstr_nocase(sa->trimmed_seq, sa->internal_input))
	pr_append_new_chunk(&sa->error,
			    "Specified internal oligo not in Included Region");
    }
    if (sa->left_input) {
      if (strlen(sa->left_input) > pa->primer_max_size)
	pr_append_new_chunk(&sa->error, "Specified left primer too long");
      if (strlen(sa->left_input) < pa->primer_min_size)
	pr_append_new_chunk(&sa->error, "Specified left primer too short");
      if (!strstr_nocase(sa->sequence, sa->left_input))
	pr_append_new_chunk(&sa->error,
			    "Specified left primer not in sequence");
      else if (!strstr_nocase(sa->trimmed_seq, sa->left_input))
	pr_append_new_chunk(&sa->error,
			    "Specified left primer not in Included Region");
    }
    if (sa->right_input) {
      if (strlen(sa->right_input) < pa->primer_min_size)
	pr_append_new_chunk(&sa->error, "Specified right primer too short");
      if (strlen(sa->right_input) > pa->primer_max_size) {
	pr_append_new_chunk(&sa->error, "Specified right primer too long");
       } else { /* We do not want to overflow s1. */
	_pr_reverse_complement(sa->right_input,s1);
	if (!strstr_nocase(sa->sequence, s1))
	  pr_append_new_chunk(&sa->error,
			      "Specified right primer not in sequence");
	else if (!strstr_nocase(sa->trimmed_seq, s1))
	  pr_append_new_chunk(&sa->error,
			      "Specified right primer not in Included Region");
      }
    }

    if ((pa->pr_pair_weights.product_tm_lt || 
	 pa->pr_pair_weights.product_tm_gt)
	&& pa->product_opt_tm == PR_UNDEFINED_DBL_OPT) {
        pr_append_new_chunk(&pa->glob_err, 
	   "Product temperature is part of objective function while optimum temperature is not defined");
        return 1;
     }
	
    if((pa->pr_pair_weights.product_size_lt ||
	pa->pr_pair_weights.product_size_gt) 
       && pa->product_opt_size == PR_UNDEFINED_INT_OPT){
       pr_append_new_chunk(&pa->glob_err,
	  "Product size is part of objective function while optimum size is not defined");
       return 1;
    }

    if ((pa->primer_weights.gc_content_lt || 
	 pa->primer_weights.gc_content_gt)
	&& pa->opt_gc_content == DEFAULT_OPT_GC_PERCENT) {
        pr_append_new_chunk(&pa->glob_err, 
	   "Primer GC content is part of objective function while optimum gc_content is not defined");
        return 1;
     }
	
    if ((pa->io_weights.gc_content_lt || 
	 pa->io_weights.gc_content_gt)
	&& pa->io_opt_gc_content == DEFAULT_OPT_GC_PERCENT) {
        pr_append_new_chunk(&pa->glob_err, 
	   "Hyb probe GC content is part of objective function while optimum gc_content is not defined");
        return 1;
     }
	
    if ((pa->primer_task != pick_pcr_primers_and_hyb_probe 
	 && pa->primer_task != pick_hyb_probe_only ) &&
			(pa->pr_pair_weights.io_quality)) {
       pr_append_new_chunk(&pa->glob_err,
	  "Internal oligo quality is part of objective function while internal oligo choice is not required");
       return 1;
    }

    if (pa->primer_weights.repeat_sim && (!pa->repeat_lib.seq_num)) {
       pr_append_new_chunk(&pa->glob_err,
	  "Mispriming score is part of objective function, but mispriming library is not defined");
       return 1;
    }

    if(pa->io_weights.repeat_sim && (!pa->io_mishyb_library.seq_num)){
      pr_append_new_chunk(&pa->glob_err,
      "Internal oligo mispriming score is part of objective function while mishyb library is not defined");
      return 1;
    }

    if(pa->pr_pair_weights.repeat_sim && (!pa->repeat_lib.seq_num)){
      pr_append_new_chunk(&pa->glob_err,
	"Mispriming score is part of objective function, but mispriming library is not defined");
      return 1;
    }

    if(pa->pr_pair_weights.io_quality 
	&& pa->primer_task != pick_pcr_primers_and_hyb_probe ) {
	  pr_append_new_chunk(&pa->glob_err,
	   "Internal oligo quality is part of objective function while internal oligo choice is not required");
        return 1;
    }

    return (NULL == sa->error.data && NULL == pa->glob_err.data) ? 0 : 1;
}

/* 
 * Check intervals, and add to sa->error.  Update the start of each interval to
 * be relative to the start of the included region.
 */ 
static int
check_intervals(tag_name, num_intervals, intervals, seq_len, sa)
    const char *tag_name;
    const int num_intervals;
    interval_array_t intervals;
    const int seq_len;
    seq_args *sa;
{
    int i;
    int outside_warning_issued = 0;
    for (i=0; i < num_intervals; i++) {
	if (intervals[i][0] + intervals[i][1] > seq_len) {
	    pr_append_new_chunk(&sa->error, tag_name);
	    pr_append(&sa->error, " beyond end of sequence");
	    return 1;
	}
	/* Cause the interval start to be relative to the included region. */
	intervals[i][0] -= sa->incl_s;
	/* Check that intervals are within the included region. */
	if (intervals[i][0] < 0
	    || intervals[i][0] + intervals[i][1] > sa->incl_l) {
	    if (!outside_warning_issued) {
		pr_append_new_chunk(&sa->warning, tag_name);
		pr_append(&sa->warning,
			  " outside of INCLUDED_REGION");
		outside_warning_issued = 1;
	    }
	}
	if (intervals[i][1] < 0) {
	    pr_append_new_chunk(&sa->error, "Negative ");
	    pr_append(&sa->error, tag_name);
	    pr_append(&sa->error, " length");
	    return 1;
	}
    }
    return 0;
}

static char * strstr_nocase(s1, s2)
char *s1, *s2;
{
   int  n1, n2;
   char *p, q, *tmp;

   if(s1 == NULL || s2 == NULL) return NULL;
   n1 = strlen(s1); n2 = strlen(s2);
   if(n1 < n2) return NULL;

   tmp = pr_safe_malloc(n1 + 1);
   strcpy(tmp, s1);

   q = *tmp; p = tmp;
   while(q != '\0' && q != '\n'){
      q = *(p + n2);
      *(p + n2) = '\0';
      if(strcmp_nocase(p, s2)){
	 *(p + n2) = q; p++; continue;
      }
      else {free(tmp); return p;}
   }
   free(tmp); return NULL;
}

/* Upcase a DNA string, s, in place.  If amibiguity_code_ok is false the
   string can consist of acgtnACGTN.  If it is true then the IUB/IUPAC
   ambiguity codes are are allowed.  Return the first unrecognized letter if
   any is seen (and turn it to 'N' in s).  Otherwise return '\0'.  */
static char
dna_to_upper(s, ambiguity_code_ok)
    char * s;
    int ambiguity_code_ok;
{
  char *p = s;
  int unrecognized_base = '\0';
  while (*p) {
    switch (*p) {
      case 'a': case 'A': *p='A'; break;
      case 'c': case 'C': *p='C'; break;
      case 'g': case 'G': *p='G'; break;
      case 't': case 'T': *p='T'; break;
      case 'n': case 'N': *p='N'; break;
      default: 
	if (ambiguity_code_ok) {
	  switch (*p) {
	  case 'r': case 'R': *p = 'R'; break;
	  case 'y': case 'Y': *p = 'Y'; break;
	  case 'm': case 'M': *p = 'M'; break;
	  case 'w': case 'W': *p = 'W'; break;
	  case 's': case 'S': *p = 'S'; break;
	  case 'k': case 'K': *p = 'K'; break;
	  case 'd': case 'D': *p = 'D'; break;
	  case 'h': case 'H': *p = 'H'; break;
	  case 'v': case 'V': *p = 'V'; break;
	  case 'b': case 'B': *p = 'B'; break;
	  }
	} else {
	  if (!unrecognized_base) unrecognized_base = *p;
	  *p = 'N';
	}
	break;
      }
    p++;
  }
  return unrecognized_base;
}

int
_pr_need_template_mispriming(pa)
  const primer_args *pa;
{
  return 
    pa->max_template_mispriming >= 0
    || pa->primer_weights.template_mispriming > 0.0
    || _pr_need_pair_template_mispriming(pa);
}

int
_pr_need_pair_template_mispriming(pa)
  const primer_args *pa;
{
  return 
    pa->pair_max_template_mispriming >= 0
    || pa->pr_pair_weights.template_mispriming > 0.0;
}

/* Takes substring of seq starting from n with length m and puts it to s.    */
void
_pr_substr(const char *seq, int n, int m, char *s)
{
	int i;
	for(i=n;i<n+m;i++)s[i-n]=seq[i];
	s[m]='\0';
}

/* Reverse and complement the sequence seq and put the result in s. */ 
void
_pr_reverse_complement(const char *seq, char *s)
{
    const char *p = seq;
    char *q = s;

    while (*p != '\0') p++;
    p--;
    while (p >= seq) {
	switch (*p)
	{
	case 'A': *q='T'; break;
	case 'C': *q='G'; break;
	case 'G': *q='C'; break;
	case 'T': *q='A'; break;
	case 'U': *q='A'; break;

	case 'B': *q='V'; break;
	case 'D': *q='H'; break;
        case 'H': *q='D'; break;
        case 'V': *q='B'; break;
        case 'R': *q='Y'; break;
        case 'Y': *q='R'; break;
	case 'K': *q='M'; break;
        case 'M': *q='K'; break;
        case 'S': *q='S'; break;
        case 'W': *q='W'; break;

	case 'N': *q='N'; break;

	case 'a': *q='t'; break;
	case 'c': *q='g'; break;
	case 'g': *q='c'; break;
	case 't': *q='a'; break;
	case 'u': *q='a'; break;

	case 'b': *q='v'; break;
	case 'd': *q='h'; break;
        case 'h': *q='d'; break;
        case 'v': *q='b'; break;
        case 'r': *q='y'; break;
        case 'y': *q='r'; break;
	case 'k': *q='m'; break;
        case 'm': *q='k'; break;
        case 's': *q='s'; break;
        case 'w': *q='w'; break;

	case 'n': *q='n'; break;
	}
	p--; q++;
    }
    *q = '\0';
}
