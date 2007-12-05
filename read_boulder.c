/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007
Whitehead Institute for Biomedical Research, Steve Rozen
(http://jura.wi.mit.edu/rozen), and Helen Skaletsky
All rights reserved.

    This file is part of primer3 and the primer3 suite.

    Primer3 and the primer3 suite are free software;
    you can redistribute them and/or modify them under the terms
    of the GNU General Public License as published by the Free
    Software Foundation; either version 2 of the License, or (at
    your option) any later version.

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this file (file gpl-2.0.txt in the source
    distribution); if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

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
#include <string.h> /* memset, strlen,  strcmp, ... */
#include <unistd.h> /* write */
#include "read_boulder.h"

#define INIT_BUF_SIZE 1024
#define INIT_LIB_SIZE  500
#define PR_MAX_LIBRARY_WT 100.0

/* Static functions. */
static void *_rb_safe_malloc(size_t x);

static void   parse_align_score(const char *, const char *, short *,
				pr_append_str *);
static void   parse_double(const char *, const char *, double *,
			   pr_append_str *);
static void   parse_int(const char *, const char *, int *, pr_append_str *);
static const char *parse_int_pair(const char *, const char *, char, int *, int *,
			    pr_append_str *);

/* static void   parse_interval_list_old(const char *, const char *, int*,
   interval_array_t, pr_append_str *); */

static void   parse_interval_list(const char *tag_name,
				   const char *datum,
				   interval_array_t2 *interval_arr,
				   pr_append_str *err);


static void   parse_product_size(const char *, char *, p3_global_settings *,
				 pr_append_str *);
static void   tag_syntax_error(const char *, const char *,  pr_append_str *);
static int    parse_seq_quality(char *, int **);

/* 
 * Hack to support old SunOS headers.  (We do not try to declare _all_
 * undeclared functions; only those with non-int return types.)
 */
#ifndef __cplusplus
extern double strtod();
#endif

/* 
 * See read_boulder.h for description.
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
           T = _rb_safe_malloc(datum_len + 1);      \
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

#if 0
#define COMPARE_INTERVAL_LIST(TAG, SIZE, LIST)                   \
   if (COMPARE(TAG)) {                                           \
       parse_interval_list(TAG, datum, &SIZE, LIST, parse_err);  \
       continue;                                                 \
   }
#endif

#define COMPARE_INTERVAL_LIST(TAG, PLACE)                   \
   if (COMPARE(TAG)) {                                           \
       parse_interval_list(TAG, datum, PLACE, parse_err);  \
       continue;                                                 \
   }


/* 
 * See read_boulder.h for description.
 */
int
read_record(const program_args *prog_args, 
	    int   echo_output,
	    p3_global_settings *pa, 
	    seq_args *sa, 
	    pr_append_str *glob_err,  /* Really should be called fatal_parse_err */
	    pr_append_str *nonfatal_parse_err) { 
  int line_len; /* seq_len; n_quality; */
    int tag_len, datum_len;
    int data_found = 0;
    int pick_internal_oligo = 2;
    char *s, *n, *datum, *task_tmp = NULL;
    const char *p;
    pr_append_str *parse_err;
    pr_append_str *non_fatal_err;
    char *repeat_file_path = NULL, *int_repeat_file_path = NULL;

    /* FIX ME call p3_create_seq_args inside read_boulder? */

    non_fatal_err = nonfatal_parse_err;

    while ((s = p3_read_line(stdin)) != NULL && strcmp(s,"=")) {
	data_found = 1;
	if (echo_output) printf("%s\n", s);
	line_len = strlen(s);
	if ((n=strchr(s,'=')) == NULL) {
	    /* 
	     * The input line is illegal because it has no
	     * "=" in it, but we still will read to the end
	     * of the record.
	     */
	    pr_append_new_chunk(glob_err, "Input line with no '=': ");
	    pr_append(glob_err, s);
	} else {	
	    tag_len = n - s;
	    datum = n + 1;
	    datum_len = line_len - tag_len - 1;
	    
	    /* 
	     * Process "Sequence" (i.e. Per-Record) Arguments".
	     */
	    parse_err = non_fatal_err;

	    /* COMPARE_AND_MALLOC("SEQUENCE", sa->sequence); */
	    if (COMPARE("SEQUENCE")) {   /* NEW WAY */
	      if (/* p3_get_seq_arg_sequence(sa) */ sa->sequence) {
		pr_append_new_chunk(parse_err,
				    "Duplicate tag: ");
		pr_append(parse_err, "SEQUENCE"); 
	      } else {
		/* p3_set_seq_arg_sequence */
		if (p3_set_seq_args_sequence(sa, datum)) exit(-2);
	      }
	      continue;
	    }

	    if (COMPARE("PRIMER_SEQUENCE_QUALITY")) {
	       if ((sa->n_quality = parse_seq_quality(datum, &sa->quality)) == 0) {
		 pr_append_new_chunk(parse_err, /*&sa->error, */ 
				     "Error in sequence quality data");
		 continue;  /* FIX ME superfluous ? */
               }
	       continue;
            }

	    COMPARE_AND_MALLOC("PRIMER_SEQUENCE_ID", sa->sequence_name);
	    COMPARE_AND_MALLOC("MARKER_NAME", sa->sequence_name);
            COMPARE_AND_MALLOC("PRIMER_LEFT_INPUT", sa->left_input);
            COMPARE_AND_MALLOC("PRIMER_RIGHT_INPUT", sa->right_input);
            COMPARE_AND_MALLOC("PRIMER_INTERNAL_OLIGO_INPUT", sa->internal_input);

	    COMPARE_INTERVAL_LIST("TARGET", &sa->tar2);
	    COMPARE_INTERVAL_LIST("EXCLUDED_REGION", &sa->excl2);
	    COMPARE_INTERVAL_LIST("PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION",
				  &sa->excl_internal2);

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
	     * Process "Global" Arguments (those that persist between boulder
	     * records).
	     */
	    parse_err = glob_err;  /* These errors are considered fatal. */
	    if (COMPARE("PRIMER_PRODUCT_SIZE_RANGE")
		|| COMPARE("PRIMER_DEFAULT_PRODUCT")) {
		parse_product_size("PRIMER_PRODUCT_SIZE_RANGE", datum, pa,
				   parse_err);
		continue;
	    }
	    COMPARE_INT("PRIMER_DEFAULT_SIZE", pa->p_args.opt_size);
	    COMPARE_INT("PRIMER_OPT_SIZE", pa->p_args.opt_size);
	    COMPARE_INT("PRIMER_MIN_SIZE", pa->p_args.min_size);
	    COMPARE_INT("PRIMER_MAX_SIZE", pa->p_args.max_size);
	    COMPARE_INT("PRIMER_MAX_POLY_X", pa->p_args.max_poly_x);
	    COMPARE_FLOAT("PRIMER_OPT_TM", pa->p_args.opt_tm);
	    COMPARE_FLOAT("PRIMER_OPT_GC_PERCENT", pa->p_args.opt_gc_content);
	    COMPARE_FLOAT("PRIMER_MIN_TM", pa->p_args.min_tm);
	    COMPARE_FLOAT("PRIMER_MAX_TM", pa->p_args.max_tm);
	    COMPARE_FLOAT("PRIMER_MAX_DIFF_TM", pa->max_diff_tm);

	    COMPARE_INT("PRIMER_TM_SANTALUCIA",
			pa->tm_santalucia);    /* added by T.Koressaar */
	    COMPARE_INT("PRIMER_SALT_CORRECTIONS",
			pa->salt_corrections); /* added by T.Koressaar */
	    COMPARE_FLOAT("PRIMER_MIN_GC", pa->p_args.min_gc);
	    COMPARE_FLOAT("PRIMER_MAX_GC", pa->p_args.max_gc);
	    COMPARE_FLOAT("PRIMER_SALT_CONC", pa->p_args.salt_conc);

	    /* added by T.Koressaar: */
	    COMPARE_FLOAT("PRIMER_DIVALENT_CONC", pa->p_args.divalent_conc);

	    /* added by T.Koressaar: */
	    COMPARE_FLOAT("PRIMER_DNTP_CONC", pa->p_args.dntp_conc);

	    COMPARE_FLOAT("PRIMER_DNA_CONC", pa->p_args.dna_conc);
	    COMPARE_INT("PRIMER_NUM_NS_ACCEPTED", pa->p_args.num_ns_accepted);
	    COMPARE_INT("PRIMER_PRODUCT_OPT_SIZE", pa->product_opt_size);
	    COMPARE_ALIGN_SCORE("PRIMER_SELF_ANY", pa->p_args.max_self_any);
	    COMPARE_ALIGN_SCORE("PRIMER_SELF_END", pa->p_args.max_self_end);
	    COMPARE_INT("PRIMER_FILE_FLAG", pa->file_flag);
	    COMPARE_INT("PRIMER_PICK_ANYWAY", pa->pick_anyway);
	    COMPARE_INT("PRIMER_GC_CLAMP", pa->gc_clamp);
	    COMPARE_INT("PRIMER_EXPLAIN_FLAG", pa->explain_flag);
	    COMPARE_INT("PRIMER_LIBERAL_BASE", pa->liberal_base);
	    COMPARE_INT("PRIMER_FIRST_BASE_INDEX", pa->first_base_index);
	    COMPARE_INT("PRIMER_NUM_RETURN", pa->num_return);
	    COMPARE_INT("PRIMER_MIN_QUALITY", pa->p_args.min_quality);
	    COMPARE_INT("PRIMER_MIN_END_QUALITY", pa->p_args.min_end_quality);
	    COMPARE_INT("PRIMER_QUALITY_RANGE_MIN",
				   pa->quality_range_min);
            COMPARE_INT("PRIMER_QUALITY_RANGE_MAX",
				   pa->quality_range_max);

	    COMPARE_FLOAT("PRIMER_PRODUCT_MAX_TM", pa->product_max_tm);
	    COMPARE_FLOAT("PRIMER_PRODUCT_MIN_TM", pa->product_min_tm);
	    COMPARE_FLOAT("PRIMER_PRODUCT_OPT_TM", pa->product_opt_tm);

	    COMPARE_INT("PRIMER_PICK_INTERNAL_OLIGO", pick_internal_oligo);

	    COMPARE_AND_MALLOC("PRIMER_TASK", task_tmp);

	    COMPARE_INT("PRIMER_INTERNAL_OLIGO_OPT_SIZE", pa->o_args.opt_size);
	    COMPARE_INT("PRIMER_INTERNAL_OLIGO_MAX_SIZE", pa->o_args.max_size);
	    COMPARE_INT("PRIMER_INTERNAL_OLIGO_MIN_SIZE", pa->o_args.min_size);
	    COMPARE_INT("PRIMER_INTERNAL_OLIGO_MAX_POLY_X", pa->o_args.max_poly_x);
	    COMPARE_FLOAT("PRIMER_INTERNAL_OLIGO_OPT_TM", pa->o_args.opt_tm);
	    COMPARE_FLOAT("PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT",
			  pa->o_args.opt_gc_content);
	    COMPARE_FLOAT("PRIMER_INTERNAL_OLIGO_MAX_TM", pa->o_args.max_tm);
	    COMPARE_FLOAT("PRIMER_INTERNAL_OLIGO_MIN_TM", pa->o_args.min_tm);
	    COMPARE_FLOAT("PRIMER_INTERNAL_OLIGO_MIN_GC", pa->o_args.min_gc);
            COMPARE_FLOAT("PRIMER_INTERNAL_OLIGO_MAX_GC", pa->o_args.max_gc);
	    COMPARE_FLOAT("PRIMER_INTERNAL_OLIGO_SALT_CONC",  pa->o_args.salt_conc);

	    COMPARE_FLOAT("PRIMER_INTERNAL_OLIGO_DIVALENT_CONC",
			  pa->o_args.divalent_conc); /* added by T.Koressaar */
	    COMPARE_FLOAT("PRIMER_INTERNAL_OLIGO_DNTP_CONC",
			  pa->o_args.dntp_conc); /* added by T.Koressaar */

	    COMPARE_FLOAT("PRIMER_INTERNAL_OLIGO_DNA_CONC", pa->o_args.dna_conc);
	    COMPARE_INT("PRIMER_INTERNAL_OLIGO_NUM_NS", pa->o_args.num_ns_accepted);
	    COMPARE_INT("PRIMER_INTERNAL_OLIGO_MIN_QUALITY", pa->o_args.min_quality);

	    COMPARE_ALIGN_SCORE("PRIMER_INTERNAL_OLIGO_SELF_ANY",
				pa->o_args.max_self_any);
	    COMPARE_ALIGN_SCORE("PRIMER_INTERNAL_OLIGO_SELF_END", 
				pa->o_args.max_self_end);
	    COMPARE_ALIGN_SCORE("PRIMER_MAX_MISPRIMING",
				pa->p_args.max_repeat_compl);
	    COMPARE_ALIGN_SCORE("PRIMER_INTERNAL_OLIGO_MAX_MISHYB",
				pa->o_args.max_repeat_compl);
	    COMPARE_ALIGN_SCORE("PRIMER_PAIR_MAX_MISPRIMING",
				pa->pair_repeat_compl);

	    /* Mispriming / mishybing in the template. */
	    COMPARE_ALIGN_SCORE("PRIMER_MAX_TEMPLATE_MISPRIMING",
				pa->p_args.max_template_mispriming);
	    COMPARE_ALIGN_SCORE("PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING",
				pa->pair_max_template_mispriming);
	    COMPARE_ALIGN_SCORE("PRIMER_INTERNAL_OLIGO_MAX_TEMPLATE_MISHYB",
				pa->o_args.max_template_mispriming);

            /* Control interpretation of ambiguity codes in mispriming
               and mishyb libraries. */
	    COMPARE_INT("PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS",
			pa->lib_ambiguity_codes_consensus);


	    COMPARE_FLOAT("PRIMER_INSIDE_PENALTY", pa->inside_penalty);
	    COMPARE_FLOAT("PRIMER_OUTSIDE_PENALTY", pa->outside_penalty);
            if (COMPARE("PRIMER_MISPRIMING_LIBRARY")) {
		if (repeat_file_path != NULL) {
		    pr_append_new_chunk(glob_err,
					"Duplicate PRIMER_MISPRIMING_LIBRARY tag");
		    free(repeat_file_path);
		    repeat_file_path = NULL;
		} else {
		    repeat_file_path = _rb_safe_malloc(strlen(datum) + 1);
		    strcpy(repeat_file_path, datum);
		}
		continue;
	    }
            if (COMPARE("PRIMER_INTERNAL_OLIGO_MISHYB_LIBRARY")) {
		if (int_repeat_file_path != NULL) {
		    pr_append_new_chunk(glob_err,
					"Duplicate PRIMER_INTERNAL_OLIGO_MISHYB_LIBRARY tag");
		    free(int_repeat_file_path);
		    int_repeat_file_path = NULL;
		} else {
		    int_repeat_file_path = _rb_safe_malloc(strlen(datum) + 1);
		    strcpy(int_repeat_file_path, datum);
		}
		continue;
	    }
	    if (COMPARE("PRIMER_COMMENT") || COMPARE("COMMENT")) continue;
	    COMPARE_FLOAT("PRIMER_MAX_END_STABILITY", pa->max_end_stability);

	    COMPARE_INT("PRIMER_LOWERCASE_MASKING",
			pa->lowercase_masking); /* added by T. Koressaar */

	    /* weights for objective functions  */
            /* CHANGE TEMP/temp -> TM/tm */
	    COMPARE_FLOAT("PRIMER_WT_TM_GT", pa->p_args.weights.temp_gt);
	    COMPARE_FLOAT("PRIMER_WT_TM_LT", pa->p_args.weights.temp_lt);
	    COMPARE_FLOAT("PRIMER_WT_GC_PERCENT_GT", pa->p_args.weights.gc_content_gt);
	    COMPARE_FLOAT("PRIMER_WT_GC_PERCENT_LT", pa->p_args.weights.gc_content_lt);
	    COMPARE_FLOAT("PRIMER_WT_SIZE_LT", pa->p_args.weights.length_lt);
	    COMPARE_FLOAT("PRIMER_WT_SIZE_GT", pa->p_args.weights.length_gt);
	    COMPARE_FLOAT("PRIMER_WT_COMPL_ANY", pa->p_args.weights.compl_any);
	    COMPARE_FLOAT("PRIMER_WT_COMPL_END", pa->p_args.weights.compl_end);
	    COMPARE_FLOAT("PRIMER_WT_NUM_NS", pa->p_args.weights.num_ns);
	    COMPARE_FLOAT("PRIMER_WT_REP_SIM", pa->p_args.weights.repeat_sim);
	    COMPARE_FLOAT("PRIMER_WT_SEQ_QUAL", pa->p_args.weights.seq_quality);
	    COMPARE_FLOAT("PRIMER_WT_END_QUAL", pa->p_args.weights.end_quality);
	    COMPARE_FLOAT("PRIMER_WT_POS_PENALTY", pa->p_args.weights.pos_penalty);
	    COMPARE_FLOAT("PRIMER_WT_END_STABILITY",
			  pa->p_args.weights.end_stability);
	    COMPARE_FLOAT("PRIMER_WT_TEMPLATE_MISPRIMING",
			  pa->p_args.weights.template_mispriming);

	    COMPARE_FLOAT("PRIMER_IO_WT_TM_GT", pa->o_args.weights.temp_gt);
	    COMPARE_FLOAT("PRIMER_IO_WT_TM_LT", pa->o_args.weights.temp_lt);
	    COMPARE_FLOAT("PRIMER_IO_WT_GC_PERCENT_GT", pa->o_args.weights.gc_content_gt);
	    COMPARE_FLOAT("PRIMER_IO_WT_GC_PERCENT_LT", pa->o_args.weights.gc_content_lt);
	    COMPARE_FLOAT("PRIMER_IO_WT_SIZE_LT", pa->o_args.weights.length_lt);
	    COMPARE_FLOAT("PRIMER_IO_WT_SIZE_GT", pa->o_args.weights.length_gt);
	    COMPARE_FLOAT("PRIMER_IO_WT_COMPL_ANY", pa->o_args.weights.compl_any);
	    COMPARE_FLOAT("PRIMER_IO_WT_COMPL_END", pa->o_args.weights.compl_end);
	    COMPARE_FLOAT("PRIMER_IO_WT_NUM_NS", pa->o_args.weights.num_ns);
	    COMPARE_FLOAT("PRIMER_IO_WT_REP_SIM", pa->o_args.weights.repeat_sim);
	    COMPARE_FLOAT("PRIMER_IO_WT_SEQ_QUAL", pa->o_args.weights.seq_quality);
	    COMPARE_FLOAT("PRIMER_IO_WT_END_QUAL", pa->o_args.weights.end_quality);
	    COMPARE_FLOAT("PRIMER_IO_WT_TEMPLATE_MISHYB",
			  pa->o_args.weights.template_mispriming);

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
	    pr_append_new_chunk(glob_err, "Unrecognized tag: ");
	    pr_append(glob_err, s);
	    fprintf(stderr, "Unrecognized tag: %s\n", s);
	}
    }  /* while ((s = p3_read_line(stdin)) != NULL && strcmp(s,"=")) { */

    if (NULL == s) { /* End of file. */
	if (data_found) {
	    pr_append_new_chunk(glob_err, 
				"Final record not terminated by '='");
	    return 1;
	} else return 0;
    }

    if (task_tmp != NULL) {
      pa->pick_left_primer = 0;
      pa->pick_right_primer = 0;
      pa->pick_internal_oligo = 0;
      if (!strcmp_nocase(task_tmp, "pick_pcr_primers")) {
	pa->primer_task = pick_pcr_primers;
	pa->pick_left_primer = 1;
	pa->pick_right_primer = 1;
      } else if (!strcmp_nocase(task_tmp, "pick_pcr_primers_and_hyb_probe")) {
	pa->primer_task = pick_pcr_primers_and_hyb_probe;
	pa->pick_left_primer = 1;
	pa->pick_right_primer = 1;
	pa->pick_internal_oligo = 1;
      } else if (!strcmp_nocase(task_tmp, "pick_left_only")) {
	pa->primer_task = pick_left_only;
	pa->pick_left_primer = 1;
      } else if (!strcmp_nocase(task_tmp, "pick_right_only")) {
	pa->primer_task = pick_right_only;
	pa->pick_right_primer = 1;
      } else if (!strcmp_nocase(task_tmp, "pick_hyb_probe_only")) {
	pa->primer_task = pick_hyb_probe_only;
	pa->pick_internal_oligo = 1;
      } else 
	pr_append_new_chunk(glob_err,
			    "Unrecognized PRIMER_TASK");
	  free(task_tmp);
    }

   /* 
     * WARNING: read_seq_lib uses p3_read_line, so repeat files cannot be read
     * inside the while ((s = p3_read_line(stdin))...)  loop above.
     * FIX ME, in fact the reading of the library contents probably
     * belongs inside primer3_boulder_main.c or libprimer3.c.
     */

    if (NULL != repeat_file_path) {
      destroy_seq_lib(pa->p_args.repeat_lib);
      if ('\0' == *repeat_file_path) {
	/* Input now specifies no repeat library. */
	pa->p_args.repeat_lib = NULL;
      }
      else {
	pa->p_args.repeat_lib
	  = read_and_create_seq_lib(repeat_file_path, 
				    "mispriming library");
	if(pa->p_args.repeat_lib->error.data != NULL) {
	  pr_append_new_chunk(glob_err, pa->p_args.repeat_lib->error.data);
	}
      }
      free(repeat_file_path);
      repeat_file_path = NULL;
    }

    if (NULL != int_repeat_file_path) {
      destroy_seq_lib(pa->o_args.repeat_lib);
      if ('\0' == *int_repeat_file_path) {
	/* Input now specifies no mishybridization library. */
	pa->o_args.repeat_lib = NULL;
      }
      else {
	pa->o_args.repeat_lib = 
	  read_and_create_seq_lib(int_repeat_file_path,
				  "internal oligo mishyb library");
	if(pa->o_args.repeat_lib->error.data != NULL) {
	  pr_append_new_chunk(glob_err, pa->o_args.repeat_lib->error.data);
	}
      }
      free(int_repeat_file_path);
      int_repeat_file_path = NULL;
    }
    
    /* This next belongs here rather than libprimer3, because it deals
       with potential incompatibility with old tags (kept for backward
       compatibility, and new tags.  */
    if((pick_internal_oligo == 1 || pick_internal_oligo == 0) &&
       (pa->primer_task == pick_left_only || 
	pa->primer_task == pick_right_only ||
	pa->primer_task == pick_hyb_probe_only)) 
	  pr_append_new_chunk(glob_err, 
	    "Contradiction in primer_task definition");
    else if (pick_internal_oligo == 1) 
	  pa->primer_task = 1;
    else if (pick_internal_oligo == 0) pa->primer_task = 0;

    return 1;
}
#undef COMPARE
#undef COMPARE_AND_MALLOC
#undef COMPARE_INT
#undef COMPARE_FLOAT
#undef COMPARE_INTERVAL_LIST

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
      /* FIX ME --> libprimer3.c?  */
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
parse_interval_list(const char *tag_name,
		    const char *datum,
		    interval_array_t2 *interval_arr,
		    pr_append_str *err)
{
  const char *p = datum;
  int i1, i2;
  int ret = 0;
  while (' ' == *p || '\t' == *p) p++;
  while (*p != '\0' && *p != '\n') {
    p = parse_int_pair(tag_name, p, ',', &i1, &i2, err);
    if (NULL == p) return;
    ret = p3_add_to_interval_array(interval_arr, i1, i2);
    if (ret) {
      pr_append_new_chunk(err, "Too many elements for tag ");
      pr_append(err, tag_name);
    }
  }
}

#if 0
static void
parse_interval_list_old(tag_name, datum, count, interval_array, err)
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
#endif

static void
parse_product_size(tag_name, in, pa, err)
    const char *tag_name;
    char *in;
    p3_global_settings *pa;
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

/* FIX ME -- TEST THIS A BIT MORE WITH
   off-by-one errors on length of the 
   quality score */
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
   g = *num = _rb_safe_malloc(sizeof(int)*k);
   memset(g, 0, k * sizeof(int));

   /* Skip leading blanks or tabs; FIX ME, TEST, probably not needed */
   while(*p == ' ' || *p == '\t'){
      p++;
      if (*p == '\0' || *p == '\n') return 0;
   }

   while (*q != '\0' && *q != '\n') {
      t = strtol(p, &q, 10);
      if (q == p) return i;
      p = q;
      *g = t;
      g++;
      i++;
   }
   return i;  /* FIX ME, is this branch ever taken? */
}

/* =========================================================== */
/* Fail-stop wrapper for memory allocation.                   */
/* =========================================================== */
/* 
 * Panic messages for when the program runs out of memory.
 */
#define OOM_MESSAGE      "out of memory in read_boulder\n"
#define OOM_MESSAGE_LEN  30
#define OOM_ERROR write(2, OOM_MESSAGE, OOM_MESSAGE_LEN), exit(-2)

static void *
_rb_safe_malloc(size_t x)
{
    void *r = malloc(x);
    if (NULL == r) OOM_ERROR;
    return r;
}

/* End of fail-stop wrappers. */
/* =========================================================== */
