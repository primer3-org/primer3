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

#include <limits.h>
#include <stdlib.h>  /* strtod, strtol,... */
#include <ctype.h> /* toupper */
#include "boulder_input.h"
#include "primer3_release.h"

#define INIT_BUF_SIZE 1024
#define INIT_LIB_SIZE  500

/* Static functions. */
static void   adjust_base_index_interval_list(interval_array_t, int, int);
static void   parse_align_score(const char *, const char *, short *,
				pr_append_str *);
static void   parse_double(const char *, const char *, double *,
			   pr_append_str *);
static void   parse_int(const char *, const char *, int *, pr_append_str *);
static const char *parse_int_pair(const char *, const char *, char, int *, int *,
			    pr_append_str *);
static void   parse_interval_list(const char *, const char *, int*,
				  interval_array_t, pr_append_str *);
static void   parse_product_size(const char *, char *, primer_args *,
				 pr_append_str *);
static void   tag_syntax_error(const char *, const char *,  pr_append_str *);
static void   read_seq_lib(seq_lib *, const char *, const char *);
static char   upcase_and_check_char(char *);
static char*  read_line(FILE *);
static double parse_seq_name(char *);
static void   reverse_complement_seq_lib(seq_lib *);
static int    parse_seq_quality(char *, int **);

/* 
 * Hack to support old SunOS headers.  (We do not try to declare _all_
 * undeclared functions; only those with non-int return types.)
 */
#ifndef __cplusplus
extern double strtod();
#endif

/* 
 * Read data from input stream until a "=" line occurs.  Assign parameter
 * values for primer picking and perform primary data control. Return 0 for
 * end of data and 1 otherwise.  If sa->error is not NULL the data is
 * erroneous and should not be processed. Echo the input lines to stdout.
 */
#define COMPARE(TAG) (!strncmp(s, TAG, tag_len) \
                      && ('=' == s[tag_len] || ' ' == s[tag_len]) \
                      && '\0' == TAG[tag_len])

#define COMPARE_AND_MALLOC(TAG,T)                  \
   if (COMPARE(TAG)) {                             \
       if (T) {                                    \
           pr_append_new_chunk(parse_err,          \
                               "Duplicate tag: "); \
           pr_append(parse_err, TAG);              \
       } else {                                    \
           T = pr_safe_malloc(datum_len + 1);      \
           strcpy(T, datum);                       \
       }                                           \
       continue;                                   \
   }

#define COMPARE_ALIGN_SCORE(TAG,T)                     \
   if (COMPARE(TAG)) {                                 \
       parse_align_score(TAG, datum, &(T), parse_err); \
       continue;                                       \
   }

#define COMPARE_FLOAT(TAG,T)                      \
   if (COMPARE(TAG)) {                            \
       parse_double(TAG, datum, &(T), parse_err); \
       continue;                                  \
   }

#define COMPARE_INT(TAG,T)                     \
   if (COMPARE(TAG)) {                         \
       parse_int(TAG, datum, &(T), parse_err); \
       continue;                               \
   }

#define COMPARE_INTERVAL_LIST(TAG, SIZE, LIST)                   \
   if (COMPARE(TAG)) {                                           \
       parse_interval_list(TAG, datum, &SIZE, LIST, parse_err);  \
       continue;                                                 \
   }

int
read_record(prog_args, pa, sa) 
    const program_args *prog_args;
    primer_args *pa;
    seq_args *sa; 
{ 
    int line_len, seq_len, n_quality;
    int tag_len, datum_len;
    int data_found = 0;
    int pick_internal_oligo = 2;
    char *s, *n, *datum, *task_tmp = NULL;
    const char *p;
    pr_append_str *parse_err;
    char *repeat_file = NULL, *int_repeat_file = NULL;

    memset(&sa->error, 0, sizeof(sa->error));
    memset(&pa->glob_err, 0, sizeof(pa->glob_err));
    memset(sa, 0, sizeof(*sa));
    sa->start_codon_pos = PR_DEFAULT_START_CODON_POS;
    sa->incl_l = -1; /* Indicates logical NULL. */
    n_quality = 0;

    while ((s = read_line(stdin)) != NULL && strcmp(s,"=")) {
	data_found = 1;
	if (0 == prog_args->format_output) printf("%s\n", s);
	line_len = strlen(s);
	if ((n=strchr(s,'=')) == NULL) {
	    /* 
	     * The input line is illegal, but we still have to read to the end
	     * of the record.
	     */
	    pr_append_new_chunk(&pa->glob_err, "Input line with no '=': ");
	    pr_append(&pa->glob_err, s);
	} else {	
	    tag_len = n - s;
	    datum = n + 1;
	    datum_len = line_len - tag_len - 1;
	    
	    /* 
	     * "Sequence" (i.e. Per-Record) Arguments".
	     */
	    parse_err = &sa->error;
	    COMPARE_AND_MALLOC("SEQUENCE", sa->sequence);
	    if (COMPARE("PRIMER_SEQUENCE_QUALITY")) {
	       if((n_quality = parse_seq_quality(datum, &sa->quality)) == 0){
		   pr_append_new_chunk(&sa->error, 
				     "Error in sequence quality data");
		   continue;
               }
	       continue;
            }

	    COMPARE_AND_MALLOC("PRIMER_SEQUENCE_ID", sa->sequence_name);
	    COMPARE_AND_MALLOC("MARKER_NAME", sa->sequence_name);
            COMPARE_AND_MALLOC("PRIMER_LEFT_INPUT", sa->left_input);
            COMPARE_AND_MALLOC("PRIMER_RIGHT_INPUT", sa->right_input);
            COMPARE_AND_MALLOC("PRIMER_INTERNAL_OLIGO_INPUT", sa->internal_input);

	    COMPARE_INTERVAL_LIST("TARGET", sa->num_targets, sa->tar) ;
	    COMPARE_INTERVAL_LIST("EXCLUDED_REGION", sa->num_excl,
				  sa->excl);
	    COMPARE_INTERVAL_LIST("PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION",
				   sa->num_internal_excl,  
				   sa->excl_internal);
	    if (COMPARE("INCLUDED_REGION")) {
		p = parse_int_pair("INCLUDED_REGION", datum, ',',
				   &sa->incl_s, &sa->incl_l, parse_err);
		if (NULL == p) /* 
                                * An error; the message is already
                                * in parse_err.
                                */
		  continue;

		while (' ' == *p || '\t' == *p) p++;
		if (*p != '\n' && *p != '\0')
		    tag_syntax_error("INCLUDED_REGION", datum,
				     parse_err);
		continue;
	    }

	    COMPARE_INT("PRIMER_START_CODON_POSITION", sa->start_codon_pos);

	    /* 
	     * "Global" Arguments (those that persist between boulder
	     * records).
	     */
	    parse_err = &pa->glob_err;
	    if (COMPARE("PRIMER_PRODUCT_SIZE_RANGE")
		|| COMPARE("PRIMER_DEFAULT_PRODUCT")) {
		parse_product_size("PRIMER_PRODUCT_SIZE_RANGE", datum, pa,
				   parse_err);
		continue;
	    }
	    COMPARE_INT("PRIMER_DEFAULT_SIZE", pa->primer_opt_size);
	    COMPARE_INT("PRIMER_OPT_SIZE", pa->primer_opt_size);
	    COMPARE_INT("PRIMER_MIN_SIZE", pa->primer_min_size);
	    COMPARE_INT("PRIMER_MAX_SIZE", pa->primer_max_size);
	    COMPARE_INT("PRIMER_MAX_POLY_X", pa->max_poly_x);
	    COMPARE_FLOAT("PRIMER_OPT_TM", pa->opt_tm);
	    COMPARE_FLOAT("PRIMER_OPT_GC_PERCENT", pa->opt_gc_content);
	    COMPARE_FLOAT("PRIMER_MIN_TM", pa->min_tm);
	    COMPARE_FLOAT("PRIMER_MAX_TM", pa->max_tm);
	    COMPARE_FLOAT("PRIMER_MAX_DIFF_TM", pa->max_diff_tm);

	    COMPARE_INT("PRIMER_TM_SANTALUCIA",
			pa->tm_santalucia);    /* added by T.Koressaar */
	    COMPARE_INT("PRIMER_SALT_CORRECTIONS",
			pa->salt_corrections); /* added by T.Koressaar */
	    COMPARE_FLOAT("PRIMER_MIN_GC", pa->min_gc);
	    COMPARE_FLOAT("PRIMER_MAX_GC", pa->max_gc);
	    COMPARE_FLOAT("PRIMER_SALT_CONC", pa->salt_conc);
	    COMPARE_FLOAT("PRIMER_DIVALENT_CONC", pa->divalent_conc); /* added by T.Koressaar */
	    COMPARE_FLOAT("PRIMER_DNTP_CONC", pa->dntp_conc); /* added by T.Koressaar */
	    COMPARE_FLOAT("PRIMER_DNA_CONC", pa->dna_conc);
	    COMPARE_INT("PRIMER_NUM_NS_ACCEPTED", pa->num_ns_accepted);
	    COMPARE_INT("PRIMER_PRODUCT_OPT_SIZE", pa->product_opt_size);
	    COMPARE_ALIGN_SCORE("PRIMER_SELF_ANY", pa->self_any);
	    COMPARE_ALIGN_SCORE("PRIMER_SELF_END", pa->self_end);
	    COMPARE_INT("PRIMER_FILE_FLAG", pa->file_flag);
	    COMPARE_INT("PRIMER_PICK_ANYWAY", pa->pick_anyway);
	    COMPARE_INT("PRIMER_GC_CLAMP", pa->gc_clamp);
	    COMPARE_INT("PRIMER_EXPLAIN_FLAG", pa->explain_flag);
	    COMPARE_INT("PRIMER_LIBERAL_BASE", pa->liberal_base);
	    COMPARE_INT("PRIMER_FIRST_BASE_INDEX", pa->first_base_index);
	    COMPARE_INT("PRIMER_NUM_RETURN", pa->num_return);
	    COMPARE_INT("PRIMER_MIN_QUALITY", pa->min_quality);
	    COMPARE_INT("PRIMER_MIN_END_QUALITY", pa->min_end_quality);
	    COMPARE_INT("PRIMER_QUALITY_RANGE_MIN",
				   pa->quality_range_min);
            COMPARE_INT("PRIMER_QUALITY_RANGE_MAX",
				   pa->quality_range_max);

	    COMPARE_FLOAT("PRIMER_PRODUCT_MAX_TM", pa->product_max_tm);
	    COMPARE_FLOAT("PRIMER_PRODUCT_MIN_TM", pa->product_min_tm);
	    COMPARE_FLOAT("PRIMER_PRODUCT_OPT_TM", pa->product_opt_tm);

	    COMPARE_INT("PRIMER_PICK_INTERNAL_OLIGO", pick_internal_oligo);

	    COMPARE_AND_MALLOC("PRIMER_TASK", task_tmp);

	    COMPARE_INT("PRIMER_INTERNAL_OLIGO_OPT_SIZE", pa->io_primer_opt_size);
	    COMPARE_INT("PRIMER_INTERNAL_OLIGO_MAX_SIZE", pa->io_primer_max_size);
	    COMPARE_INT("PRIMER_INTERNAL_OLIGO_MIN_SIZE", pa->io_primer_min_size);
	    COMPARE_INT("PRIMER_INTERNAL_OLIGO_MAX_POLY_X", pa->io_max_poly_x);
	    COMPARE_FLOAT("PRIMER_INTERNAL_OLIGO_OPT_TM", pa->io_opt_tm);
	    COMPARE_FLOAT("PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT", pa->io_opt_gc_content);
	    COMPARE_FLOAT("PRIMER_INTERNAL_OLIGO_MAX_TM", pa->io_max_tm);
	    COMPARE_FLOAT("PRIMER_INTERNAL_OLIGO_MIN_TM", pa->io_min_tm);
	    COMPARE_FLOAT("PRIMER_INTERNAL_OLIGO_MIN_GC", pa->io_min_gc);
            COMPARE_FLOAT("PRIMER_INTERNAL_OLIGO_MAX_GC", pa->io_max_gc);
	    COMPARE_FLOAT("PRIMER_INTERNAL_OLIGO_SALT_CONC",pa->io_salt_conc);
           COMPARE_FLOAT("PRIMER_INTERNAL_OLIGO_DIVALENT_CONC",pa->io_divalent_conc);
	   COMPARE_FLOAT("PRIMER_INTERNAL_OLIGO_DNTP_CONC",pa->io_dntp_conc); /* added by T.Koressaar */
	   COMPARE_FLOAT("PRIMER_INTERNAL_OLIGO_DNA_CONC", pa->io_dna_conc); /* added by T.Koressaar */
	    COMPARE_INT("PRIMER_INTERNAL_OLIGO_NUM_NS", pa->io_num_ns_accepted);
	    COMPARE_INT("PRIMER_INTERNAL_OLIGO_MIN_QUALITY", pa->io_min_quality);

	    COMPARE_ALIGN_SCORE("PRIMER_INTERNAL_OLIGO_SELF_ANY",
				pa->io_self_any);
	    COMPARE_ALIGN_SCORE("PRIMER_INTERNAL_OLIGO_SELF_END", 
				pa->io_self_end);
	    COMPARE_ALIGN_SCORE("PRIMER_MAX_MISPRIMING",
				pa->repeat_compl);
	    COMPARE_ALIGN_SCORE("PRIMER_INTERNAL_OLIGO_MAX_MISHYB",
				pa->io_repeat_compl);
	    COMPARE_ALIGN_SCORE("PRIMER_PAIR_MAX_MISPRIMING",
				pa->pair_repeat_compl);

	    /* Mispriming / mishybing in the template. */
	    COMPARE_ALIGN_SCORE("PRIMER_MAX_TEMPLATE_MISPRIMING",
				pa->max_template_mispriming);
	    COMPARE_ALIGN_SCORE("PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING",
				pa->pair_max_template_mispriming);
	    COMPARE_ALIGN_SCORE("PRIMER_INTERNAL_OLIGO_MAX_TEMPLATE_MISHYB",
				pa->io_max_template_mishyb);

            /* Control interpretation of ambiguity codes in mispriming
               and mishyb libraries. */
	    COMPARE_INT("PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS",
			pa->lib_ambiguity_codes_consensus);


	    COMPARE_FLOAT("PRIMER_INSIDE_PENALTY", pa->inside_penalty);
	    COMPARE_FLOAT("PRIMER_OUTSIDE_PENALTY", pa->outside_penalty);
            if (COMPARE("PRIMER_MISPRIMING_LIBRARY")) {
		if (repeat_file != NULL) {
		    pr_append_new_chunk(&pa->glob_err,
					"Duplicate PRIMER_MISPRIMING_LIBRARY tag");
		    free(repeat_file);
		    repeat_file = NULL;
		} else {
		    repeat_file = pr_safe_malloc(strlen(datum) + 1);
		    strcpy(repeat_file, datum);
		}
		continue;
	    }
            if (COMPARE("PRIMER_INTERNAL_OLIGO_MISHYB_LIBRARY")) {
		if (int_repeat_file != NULL) {
		    pr_append_new_chunk(&pa->glob_err,
					"Duplicate PRIMER_INTERNAL_OLIGO_MISHYB_LIBRARY tag");
		    free(int_repeat_file);
		    int_repeat_file = NULL;
		} else {
		    int_repeat_file = pr_safe_malloc(strlen(datum) + 1);
		    strcpy(int_repeat_file, datum);
		}
		continue;
	    }
	    if (COMPARE("PRIMER_COMMENT") || COMPARE("COMMENT")) continue;
	    COMPARE_FLOAT("PRIMER_MAX_END_STABILITY", pa->max_end_stability);

	    COMPARE_INT("PRIMER_LOWERCASE_MASKING",
			pa->lowercase_masking); /* added by T. Koressaar */

	    /* weights for objective functions  */
            /* CHANGE TEMP/temp -> TM/tm */
	    COMPARE_FLOAT("PRIMER_WT_TM_GT", pa->primer_weights.temp_gt);
	    COMPARE_FLOAT("PRIMER_WT_TM_LT", pa->primer_weights.temp_lt);
	    COMPARE_FLOAT("PRIMER_WT_GC_PERCENT_GT", pa->primer_weights.gc_content_gt);
	    COMPARE_FLOAT("PRIMER_WT_GC_PERCENT_LT", pa->primer_weights.gc_content_lt);
	    COMPARE_FLOAT("PRIMER_WT_SIZE_LT", pa->primer_weights.length_lt);
	    COMPARE_FLOAT("PRIMER_WT_SIZE_GT", pa->primer_weights.length_gt);
	    COMPARE_FLOAT("PRIMER_WT_COMPL_ANY", pa->primer_weights.compl_any);
	    COMPARE_FLOAT("PRIMER_WT_COMPL_END", pa->primer_weights.compl_end);
	    COMPARE_FLOAT("PRIMER_WT_NUM_NS", pa->primer_weights.num_ns);
	    COMPARE_FLOAT("PRIMER_WT_REP_SIM", pa->primer_weights.repeat_sim);
	    COMPARE_FLOAT("PRIMER_WT_SEQ_QUAL", pa->primer_weights.seq_quality);
	    COMPARE_FLOAT("PRIMER_WT_END_QUAL", pa->primer_weights.end_quality);
	    COMPARE_FLOAT("PRIMER_WT_POS_PENALTY", pa->primer_weights.pos_penalty);
	    COMPARE_FLOAT("PRIMER_WT_END_STABILITY",
			  pa->primer_weights.end_stability);
	    COMPARE_FLOAT("PRIMER_WT_TEMPLATE_MISPRIMING",
			  pa->primer_weights.template_mispriming);

	    COMPARE_FLOAT("PRIMER_IO_WT_TM_GT", pa->io_weights.temp_gt);
	    COMPARE_FLOAT("PRIMER_IO_WT_TM_LT", pa->io_weights.temp_lt);
	    COMPARE_FLOAT("PRIMER_IO_WT_GC_PERCENT_GT", pa->io_weights.gc_content_gt);
	    COMPARE_FLOAT("PRIMER_IO_WT_GC_PERCENT_LT", pa->io_weights.gc_content_lt);
	    COMPARE_FLOAT("PRIMER_IO_WT_SIZE_LT", pa->io_weights.length_lt);
	    COMPARE_FLOAT("PRIMER_IO_WT_SIZE_GT", pa->io_weights.length_gt);
	    COMPARE_FLOAT("PRIMER_IO_WT_COMPL_ANY", pa->io_weights.compl_any);
	    COMPARE_FLOAT("PRIMER_IO_WT_COMPL_END", pa->io_weights.compl_end);
	    COMPARE_FLOAT("PRIMER_IO_WT_NUM_NS", pa->io_weights.num_ns);
	    COMPARE_FLOAT("PRIMER_IO_WT_REP_SIM", pa->io_weights.repeat_sim);
	    COMPARE_FLOAT("PRIMER_IO_WT_SEQ_QUAL", pa->io_weights.seq_quality);
	    COMPARE_FLOAT("PRIMER_IO_WT_END_QUAL", pa->io_weights.end_quality);
	    COMPARE_FLOAT("PRIMER_IO_WT_TEMPLATE_MISHYB",
			  pa->io_weights.template_mispriming);

	    COMPARE_FLOAT("PRIMER_PAIR_WT_PR_PENALTY", 
					      pa->pr_pair_weights.primer_quality);
            COMPARE_FLOAT("PRIMER_PAIR_WT_IO_PENALTY",
					      pa->pr_pair_weights.io_quality);
            COMPARE_FLOAT("PRIMER_PAIR_WT_DIFF_TM",
					      pa->pr_pair_weights.diff_tm);
            COMPARE_FLOAT("PRIMER_PAIR_WT_COMPL_ANY",
					      pa->pr_pair_weights.compl_any);
            COMPARE_FLOAT("PRIMER_PAIR_WT_COMPL_END",
					      pa->pr_pair_weights.compl_end);

	    COMPARE_FLOAT("PRIMER_PAIR_WT_PRODUCT_TM_LT",
					      pa->pr_pair_weights.product_tm_lt);
	    COMPARE_FLOAT("PRIMER_PAIR_WT_PRODUCT_TM_GT",
					      pa->pr_pair_weights.product_tm_gt);
	    COMPARE_FLOAT("PRIMER_PAIR_WT_PRODUCT_SIZE_GT",
					   pa->pr_pair_weights.product_size_gt);
            COMPARE_FLOAT("PRIMER_PAIR_WT_PRODUCT_SIZE_LT",
					   pa->pr_pair_weights.product_size_lt);

	    COMPARE_FLOAT("PRIMER_PAIR_WT_REP_SIM",
					   pa->pr_pair_weights.repeat_sim);

	    COMPARE_FLOAT("PRIMER_PAIR_WT_TEMPLATE_MISPRIMING",
			  pa->pr_pair_weights.template_mispriming);
	}
	if (1 == prog_args->strict_tags) {
	    pr_append_new_chunk(&pa->glob_err, "Unrecognized tag: ");
	    pr_append(&pa->glob_err, s);
	    fprintf(stderr, "Unrecognized tag: %s\n", s);
	}
    }

    if (NULL == s) { /* End of file. */
	if (data_found) {
	    pr_append_new_chunk(&pa->glob_err, 
				"Final record not terminated by '='");
	    return 1;
	} else return 0;
    }

    if(task_tmp != NULL) {
          if (!strcmp_nocase(task_tmp, "pick_pcr_primers"))
	    pa->primer_task = pick_pcr_primers;
          else if (!strcmp_nocase(task_tmp, "pick_pcr_primers_and_hyb_probe"))
	    pa->primer_task = pick_pcr_primers_and_hyb_probe;
	  else if (!strcmp_nocase(task_tmp, "pick_left_only"))
	    pa->primer_task = pick_left_only;
          else if (!strcmp_nocase(task_tmp, "pick_right_only"))
	    pa->primer_task = pick_right_only;
          else if (!strcmp_nocase(task_tmp, "pick_hyb_probe_only"))
	    pa->primer_task = pick_hyb_probe_only;
          else   pr_append_new_chunk(&pa->glob_err,
				     "Unrecognized PRIMER_TASK");
	  free(task_tmp);
    }

    if (NULL == sa->sequence)
	pr_append_new_chunk(&sa->error, "Missing SEQUENCE tag");
    else {
	seq_len = strlen(sa->sequence);
	if (sa->incl_l == -1) {
	     sa->incl_l = seq_len;
	     sa->incl_s = pa->first_base_index;
	}
	if(n_quality !=0 && n_quality != seq_len)
	     pr_append_new_chunk(&sa->error, "Error in sequence quality data");
        if((pa->min_quality != 0 || pa->io_min_quality != 0) && n_quality == 0) 
	     pr_append_new_chunk(&sa->error, "Sequence quality data missing");
	if(pa->min_quality != 0 && pa->min_end_quality < pa->min_quality)
	     pa->min_end_quality = pa->min_quality;
    }

    /* 
     * WARNING: read_seq_lib uses read_line, so repeat files cannot be read
     * inside the while ((s = read_line(stdin))...)  loop above.
     */

    if (NULL != repeat_file) {
      if ('\0' == *repeat_file) {
	free_seq_lib(&pa->repeat_lib);
      }
      else {
	read_seq_lib(&pa->repeat_lib, repeat_file, "mispriming library");
	if(pa->repeat_lib.error.data != NULL) {
	  pr_append_new_chunk(&pa->glob_err, pa->repeat_lib.error.data);
	}
      }
      free(repeat_file);
      repeat_file = NULL;
    }

    if (NULL != int_repeat_file) {
      if ('\0' == *int_repeat_file) {
	free_seq_lib(&pa->io_mishyb_library);
      }
      else {
	read_seq_lib(&pa->io_mishyb_library, int_repeat_file,
		     "internal oligo mishyb library");
	if(pa->io_mishyb_library.error.data != NULL) {
	  pr_append_new_chunk(&pa->glob_err, pa->io_mishyb_library.error.data);
	}
      }
      free(int_repeat_file);
      int_repeat_file = NULL;
    }
    
    if((pick_internal_oligo == 1 || pick_internal_oligo == 0) &&
       (pa->primer_task == pick_left_only || 
	pa->primer_task == pick_right_only ||
	pa->primer_task == pick_hyb_probe_only)) 
	  pr_append_new_chunk(&pa->glob_err, 
	    "Contradiction in primer_task definition");
    else if(pick_internal_oligo == 1) 
	  pa->primer_task = 1;
    else if(pick_internal_oligo == 0)pa->primer_task = 0;


    /* Adjust base indexes in sa. */
    sa->incl_s -= pa->first_base_index;
    sa->start_codon_pos -= pa->first_base_index;
    adjust_base_index_interval_list(sa->tar, sa->num_targets,
				    pa->first_base_index);
    adjust_base_index_interval_list(sa->excl, sa->num_excl,
				    pa->first_base_index);
    adjust_base_index_interval_list(sa->excl_internal,
				    sa->num_internal_excl,
				    pa->first_base_index);
    return 1;
}
#undef COMPARE
#undef COMPARE_AND_MALLOC
#undef COMPARE_INT
#undef COMPARE_FLOAT
#undef COMPARE_INTERVAL_LIST

static void
adjust_base_index_interval_list(intervals, num, first_index)
    interval_array_t intervals;
    int num, first_index;
{
    int i;
    for (i = 0; i < num; i++) intervals[i][0] -= first_index;
}

/* 
 * Read a line of any length from stdin.  Return NULL on end of file,
 * otherwise return a pointer to static storage containing the line.  Any
 * trailing newline is stripped off.
 */
static char*
read_line(file)
FILE *file;
{
    static size_t ssz;
    static char *s = NULL;

    size_t remaining_size;
    char *p, *n;

    if (NULL == s) {
	ssz = INIT_BUF_SIZE;
	s = pr_safe_malloc(ssz);
    }
    p = s;
    remaining_size = ssz;
    while (1) {
	if (fgets(p, remaining_size, file) == NULL) /* End of file. */
	    return p == s ? NULL : s;

	if ((n = strchr(p, '\n')) != NULL) {
	    *n = '\0';
	    return s;
	}

	/* We did not get the whole line. */
	
	/* 
         * The following assertion is a bit of hack, a at least for 32-bit
         * machines, because we will usually run out of address space first.
         * Really we should treat an over-long line as an input error, but
         * since an over-long line is unlikely and we do want to provide some
         * protection....
	 */
	PR_ASSERT(ssz <= INT_MAX);
	if (ssz >= INT_MAX / 2)
	    ssz = INT_MAX;
	else {
	    ssz *= 2;
	}
	s = pr_safe_realloc(s, ssz);
	p = strchr(s, '\0');
	remaining_size = ssz - (p - s);
    }
}

static void
tag_syntax_error(tag_name, datum, err)
    const char *tag_name, *datum;
    pr_append_str *err;
{
    pr_append_new_chunk(err, "Illegal ");
    pr_append(err, tag_name);
    pr_append(err, " value: ");
    pr_append(err, datum);
}

static void
parse_align_score(tag_name, datum, out, err)
    const char *datum, *tag_name;
    short *out;
    pr_append_str *err;
{
    double d;

    parse_double(tag_name, datum, &d, err);
    d *= PR_ALIGN_SCORE_PRECISION;
    if (d > SHRT_MAX) {
	pr_append_new_chunk(err, "Value too large at tag ");
	pr_append(err, tag_name);
    } else {
	/* Should we be rounding here? */
	*out = (short)d;
    }
}    

static void
parse_double(tag_name, datum, out, err)
    const char *datum, *tag_name;
    double *out;
    pr_append_str *err;
{
    char *nptr;
    *out = strtod(datum, &nptr);
    if (nptr == datum) {
	/* Empty string or complete junk. */
	tag_syntax_error(tag_name, datum, err);
	*out = 0.0;
	return;
    }
    /* Look for trailing junk. */
    while (*nptr != '\n' && *nptr != '\0') {
	if (*nptr != ' ' && *nptr != '\t') {
	    tag_syntax_error(tag_name, datum, err);
	    break;
	}
	nptr++;
    }
}

static void
parse_int(tag_name, datum, out, err)
    const char *datum, *tag_name;
    int *out;
    pr_append_str *err;
{
    char *nptr;
    long tlong;
    tlong = strtol(datum, &nptr, 10);
    if (tlong > INT_MAX || tlong < INT_MIN) {
	tag_syntax_error(tag_name, datum, err);
	pr_append(err, " (value too large or too small)");
	return;
    }
    *out = tlong;
    if (nptr == datum) {
	/* Empty string or complete junk. */
	tag_syntax_error(tag_name, datum, err);
	return;
    }
    /* Look for trailing junk. */
    while (*nptr != '\n' && *nptr != '\0') {
	if (*nptr != ' ' && *nptr != '\t') {
	    tag_syntax_error(tag_name, datum, err);
	    break;
	}
	nptr++;
    }
}

/* 
 * For correct input, return a pointer to the first non-tab, non-space
 * character after the second integer, and place the integers in out1 and
 * out2.  On incorrect input, return NULL;
 */
static const char *
parse_int_pair(tag_name, datum, sep, out1, out2, err)
    const char    *tag_name, *datum;
    char          sep;          /* The separator, e.g. ',' or '-'. */
    int           *out1, *out2; /* The 2 integers. */
    pr_append_str *err;         /* Error messages. */
{
    char *nptr, *tmp;
    long tlong;
    tlong = strtol(datum, &nptr, 10);
    if (tlong > INT_MAX || tlong < INT_MIN) {
	tag_syntax_error(tag_name, datum, err);
	pr_append(err, " (value too large or too small)");
	return NULL;
    }
    *out1 = tlong;
    if (nptr == datum) {
	tag_syntax_error(tag_name, datum, err);
	return NULL;
    }
    while (' ' == *nptr || '\t' == *nptr) nptr++;
    if (sep != *nptr) {
	tag_syntax_error(tag_name, datum, err);
	return NULL;
    }
    nptr++; /* Advance past separator. */
    while (' ' == *nptr || '\t' == *nptr) nptr++;
    tmp = nptr;
    tlong = strtol(tmp, &nptr, 10);
    if (tlong > INT_MAX || tlong < INT_MIN) {
	tag_syntax_error(tag_name, datum, err);
	pr_append(err, " (value too large or too small)");
	return NULL;
    }
    *out2 = tlong;
    if (nptr == tmp) {
	tag_syntax_error(tag_name, datum, err);
	return NULL;
    }
    while (' ' == *nptr || '\t' == *nptr) nptr++;

    /* A hack to live with the old TARGET syntax. */
    if (',' == *nptr && !strcmp(tag_name, "TARGET")) {
	/* Skip the old-fashioned "description". */
	while(' ' != *nptr && '\t' != *nptr 
	      && '\0' != *nptr && '\n' != *nptr) nptr++;
	/* Advance to non-space, non-tab. */
	while (' ' == *nptr || '\t' == *nptr) nptr++;
    }
    return nptr;
}

static void
parse_interval_list(tag_name, datum, count, interval_array, err)
    const char *tag_name;
    const char *datum;
    int *count;
    interval_array_t interval_array;
    pr_append_str *err;
{
    const char *p = datum;
    while (' ' == *p || '\t' == *p) p++;
    while (*p != '\0' && *p != '\n') {
	if (*count >= PR_MAX_INTERVAL_ARRAY) {
	    pr_append_new_chunk(err, "Too many elements for tag ");
	    pr_append(err, tag_name);
	    return;
	}
	p = parse_int_pair(tag_name, p, ',', 
			   &interval_array[*count][0],
			   &interval_array[*count][1],
			   err);
	if (NULL == p) return;
	(*count)++;
    }
}

static void
parse_product_size(tag_name, in, pa, err)
    const char *tag_name;
    char *in;
    primer_args *pa;
    pr_append_str *err;
{
    char *q, *s = in;
    const char *p;
    int i;
    /* 
     * Handle possible double quotes around the value.
     * (This handling is needed for backward compatibility with v2.)
     */
    if ('"' == *s)  {
      s++;
      in++;
      q = strchr(s, '"');
      if (NULL == q) {
	pr_append_new_chunk(err, tag_name);
	pr_append(err, " begins but does not end with a quote");
	return;
      }
      /* Ignore the " and everything after it. */
      *q = '\0';
    }
    p = in;
    while (' ' == *p || '\t' == *p) p++;
    i = 0;
    while (*p != '\0' && *p != '\n') {
	if (i >= PR_MAX_INTERVAL_ARRAY) {
	    pr_append_new_chunk(err, "Too many values for ");
	    pr_append(err, tag_name);
	    return;
	}
	p = parse_int_pair(tag_name, p, '-',
			   &pa->pr_min[i], &pa->pr_max[i], err);
	if (NULL == p) return;
	i++;
    }
    pa->num_intervals = i;
}

/* 
 * Reads any file in fasta format and fills in *lib.  Sets lib->error to a
 * non-empty string on error.
 */
static void
read_seq_lib(lib, filename, errfrag)
    seq_lib *lib;
    const char *filename;
    const char *errfrag;
{
    char  *p;
    FILE *file;
    int i, m, k;
    size_t j, n;
    char buf[2];
    char offender = '\0', tmp;

    PR_ASSERT(NULL != lib);
    PR_ASSERT(NULL != filename);

    free_seq_lib(lib);

    lib->repeat_file = pr_safe_malloc(strlen(filename) + 1);
    strcpy(lib->repeat_file, filename);

    if((file = fopen(lib->repeat_file,"r")) == NULL) {
	pr_append_new_chunk(&lib->error,
			    "Cannot open ");
	goto ERROR;
    }

    j = INIT_BUF_SIZE;
    n = INIT_LIB_SIZE;
    lib->names = pr_safe_malloc(INIT_LIB_SIZE*sizeof(*lib->names));
    lib->seqs  = pr_safe_malloc(INIT_LIB_SIZE*sizeof(*lib->seqs));
    lib->weight= pr_safe_malloc(INIT_LIB_SIZE*sizeof(*lib->weight));
    lib->seq_num = 0;

    i = -1;  m = 0; k = 0;
    while((p = read_line(file))) {
	if(*p == '>'){
	    i++;
	    if(i >= n) {
		n += INIT_LIB_SIZE;
		lib->names = pr_safe_realloc(lib->names,n*sizeof(*lib->names));
		lib->seqs  = pr_safe_realloc(lib->seqs ,n*sizeof(*lib->seqs));
		lib->weight= pr_safe_realloc(lib->weight,
					     n*sizeof(*lib->weight));
	    }
	    p++;
	    lib->names[i] = pr_safe_malloc(strlen(p) + 1);
	    strcpy(lib->names[i],p);
	    lib->weight[i] = parse_seq_name(lib->names[i]);
	    lib->seqs[i] = pr_safe_malloc(INIT_BUF_SIZE);
	    lib->seqs[i][0] = '\0';
	    lib->seq_num = i+1;
	    if(lib->weight[i] < 0) {
		pr_append_new_chunk(&lib->error, "Illegal weight in ");
		goto ERROR;
	    }
	    j = INIT_BUF_SIZE;
	    k = 0;
	    if(i > 0) {
		/* We are actually testing the previous sequence. */
		if(strlen(lib->seqs[i-1]) == 0) {
		    pr_append_new_chunk(&lib->error, "Empty sequence in ");
		    goto ERROR;
		}
		tmp = upcase_and_check_char(lib->seqs[i-1]);
		m += tmp;
		if (tmp && '\0' == offender) offender = tmp;
	    }
	    p--;
	}
	else {
	    if(i < 0){ 
		pr_append_new_chunk(&lib->error,
				    "Missing id line (expected '>') in ");
		goto ERROR;
	    } else {
		if(k+strlen(p) > j-2){
		    while(j-2 < k+ strlen(p))j += INIT_BUF_SIZE;
		    lib->seqs[i] = pr_safe_realloc(lib->seqs[i], j);

		}
		strcat(lib->seqs[i], p);
		k += strlen(p);
	    }
	}
    }
    if(i < 0) {
	pr_append_new_chunk(&lib->error, "Empty ");
	goto ERROR;
    }
    else if(strlen(lib->seqs[i]) == 0) {
	pr_append_new_chunk(&lib->error, "Empty sequence in ");
	goto ERROR;
    }
    tmp = upcase_and_check_char(lib->seqs[i]);
    m += tmp;
    if (tmp && '\0' == offender) offender = tmp;
    if (offender) {
	pr_append_new_chunk(&lib->warning,
			    "Unrecognized character (");
	buf[0] = offender;
	buf[1] = '\0';
	pr_append(&lib->warning, buf);
	pr_append(&lib->warning, ") in ");
	pr_append(&lib->warning, errfrag);
	pr_append(&lib->warning, " ");
	pr_append(&lib->warning, lib->repeat_file);
    }
    fclose(file);
    reverse_complement_seq_lib(lib);
    return;

 ERROR:
    pr_append(&lib->error, errfrag);
    pr_append(&lib->error, " ");
    pr_append(&lib->error, lib->repeat_file);
    if (file) fclose(file);
}

/* 
 * Free exogenous storage associated with a seq_lib (but not the seq_lib
 * itself).  Silently ignore NULL p.  Set *p to 0 bytes.
 */
void
free_seq_lib(p)
    seq_lib *p;
{
    int i;
    if (NULL == p) return;

    if ( NULL != p->repeat_file) free(p->repeat_file);
    if (NULL != p->seqs) { 
	for(i = 0; i < p->seq_num; i++)
	    if (NULL != p->seqs[i]) free(p->seqs[i]);
	free(p->seqs);
    }
    if (NULL != p->names) {
	for(i = 0; i < p->seq_num; i++)
	    if (NULL != p->names[i]) free(p->names[i]);
	free(p->names);
    }
    if (NULL != p->weight) free(p->weight);
    if (NULL != p->error.data) free(p->error.data);
    if (NULL != p->warning.data) free(p->warning.data);
    if (NULL != p->rev_compl_seqs) free(p->rev_compl_seqs);
    memset(p, 0, sizeof(*p));
}

/* 
 * Removes spaces and "end-of-line" characters
 * from the sequence, replaces all other
 * characters except A, T, G, C and IUB/IUPAC
 * codes with N.  Returns 0 if there were no such
 * replacements and the first non-ACGT IUB
 * character otherwise. 
 */
static char
upcase_and_check_char(s)
    char *s;
{
    int i, j, n, m;

    j = 0; m = 0;
    n = strlen(s);
    for(i=0; i<n; i++){
      
	switch(s[i])
	{
	case 'a' : s[i-j] = 'A'; break;
	case 'g' : s[i-j] = 'G'; break;
	case 'c' : s[i-j] = 'C'; break;
	case 't' : s[i-j] = 'T'; break;
	case 'n' : s[i-j] = 'N'; break;
	case 'A' : s[i-j] = 'A'; break;
	case 'G' : s[i-j] = 'G'; break;
	case 'C' : s[i-j] = 'C'; break;
	case 'T' : s[i-j] = 'T'; break;
	case 'N' : s[i-j] = 'N'; break;

        case 'b' : case 'B': 
        case 'd' : case 'D':
        case 'h' : case 'H':
        case 'v' : case 'V':
        case 'r' : case 'R':
        case 'y' : case 'Y':
        case 'k' : case 'K':
        case 'm' : case 'M':
	case 's' : case 'S':
	case 'w' : case 'W':
	  s[i-j] = toupper(s[i]); break;

	case '\n': j++;          break;
	case ' ' : j++;          break;
	case '\t': j++;          break;
	case '\r': j++;          break;
	default  : if (!m) m = s[i]; s[i-j] = 'N'; 
	}
    }
    s[n-j] = '\0';
    return m;
}

static double
parse_seq_name(s)
char *s;
{
    char *p, *q;
    double n;

    p = s;
    while( *p != '*' && *p != '\0' ) p++;
    if (*p == '\0' ) return 1;
    else {
	 p++;
	 n = strtod( p, &q );
	 if( q == p ) return -1;
    }
    if(n > PR_MAX_LIBRARY_WT) return -1;

    return n;
}

static void
reverse_complement_seq_lib(lib)
seq_lib  *lib;
{
    int i, n, k;
    if((n = lib->seq_num) == 0) return;
    else {
	lib->names = pr_safe_realloc(lib->names, 2*n*sizeof(*lib->names));
	lib->seqs = pr_safe_realloc(lib->seqs, 2*n*sizeof(*lib->seqs));
	lib->weight = pr_safe_realloc(lib->weight, 2*n*sizeof(*lib->weight));
	lib->rev_compl_seqs = pr_safe_malloc(2*n*sizeof(*lib->seqs));

	lib->seq_num *= 2;
	for(i=n; i<lib->seq_num; i++){
	    k = strlen(lib->names[i-n]);
	    lib->names[i] = pr_safe_malloc(k + 9);
	    strcpy(lib->names[i], "reverse ");
	    strcat(lib->names[i], lib->names[i-n]);
	    lib->seqs[i] = pr_safe_malloc(strlen(lib->seqs[i-n]) + 1);
	    _pr_reverse_complement(lib->seqs[i-n], lib->seqs[i]);
	    lib->weight[i] = lib->weight[i-n];
	    lib->rev_compl_seqs[i-n] = lib->seqs[i];
	    lib->rev_compl_seqs[i] = lib->seqs[i-n];
       }
    }
    return;
}

static int
parse_seq_quality(s, num)
   char *s;
   int **num;
{
   int k, i=0, *g;
   long t;
   char *p, *q;

   p = q = s;
   k = strlen(s);
   g = *num = malloc(sizeof(int)*k);
   while(*p == ' ' || *p == '\t'){
      p++;
      if(*p == '\0' || *p == '\n') return 0;
   }
   while(*q != '\0' && *q != '\n'){
      t = strtol(p, &q, 10);
      if(q == p) return i;
      p = q;
      *g = t;
      g++;
      i++;
   }
   return i;
}

