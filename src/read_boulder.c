/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008,2009,
              2010,2011,2012
Whitehead Institute for Biomedical Research, Steve Rozen
(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky.
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
#include <ctype.h> /* toupper, isspace */
#include <string.h> /* memset, strlen,  strcmp, ... */
#include <float.h>
#include "read_boulder.h"

#define INIT_BUF_SIZE 1024
#define INIT_LIB_SIZE  500
#define PR_MAX_LIBRARY_WT 100.0

char *thermodynamic_params_path = NULL;
int   thermodynamic_path_changed = 1;

/* Static functions. */
static void *_rb_safe_malloc(size_t x);
static void  out_of_memory_error();

static void   parse_double(const char *, const char *, double *,
                           pr_append_str *);

static void   parse_int(const char *, const char *, int *, pr_append_str *);

static const char *parse_int_pair(const char *, const char *,
                                  char, int *, int *,
                                  pr_append_str *);

static char  *parse_2_int_pair(const char*, char*, char,
			       char, int*, int*, int*, int*,
			       pr_append_str*); 

static void   parse_interval_list(const char *tag_name,
                                  const char *datum,
                                  interval_array_t2 *interval_arr,
                                  pr_append_str *err);

static void   parse_2_interval_list(const char *tag_name,
				    char *datum,
				    interval_array_t4 *interval_arr,
				    pr_append_str *err);   

static int    parse_intron_list(char *, int *, int *);

static void   parse_product_size(const char *, char *, p3_global_settings *,
                                 pr_append_str *);

static void   pr_append(pr_append_str *, const char *);
static void   pr_append_new_chunk(pr_append_str *x, const char *s);

static void   tag_syntax_error(const char *, const char *,  pr_append_str *);
static int    parse_seq_quality(char *, seq_args *);

static const char *pr_program_name = "TMP";

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
           T = (char*) _rb_safe_malloc(datum_len + 1); \
           strcpy(T, datum);                       \
       }                                           \
       continue;                                   \
   }

#define OVERWRITE_COMPARE_AND_MALLOC(TAG,T)        \
   if (COMPARE(TAG)) {                             \
       if (T) {                                    \
	 free(T);				   \
       }                                           \
       T = (char*) _rb_safe_malloc(datum_len + 1); \
       strcpy(T, datum);                           \
       continue;                                   \
   }

#define COMPARE_FLOAT(TAG,T)                      \
   if (COMPARE(TAG)) {                            \
       parse_double(TAG, datum, &(T), parse_err); \
       continue;                                  \
   }

/* This macro added 2011 08 28 to allow use of 
   of the 'setting' functions, such as 
   p3_set_gs_primer_self_end, in setting
   arguments. At this point, this macro is
   only used twice, and there are no analogous
   macros for integers,etc.
*/
#define COMPARE_FLOAT_USE_FN(TAG, set_function)        \
   if (COMPARE(TAG)) {                                 \
     parse_double(TAG, datum, &tmp_double, parse_err); \
     set_function(pa, tmp_double);                     \
     continue;                                         \
   }

#define COMPARE_INT(TAG,T)                     \
   if (COMPARE(TAG)) {                         \
       parse_int(TAG, datum, &(T), parse_err); \
       continue;                               \
   }

#define COMPARE_INTERVAL_LIST(TAG, PLACE)                  \
   if (COMPARE(TAG)) {                                     \
       parse_interval_list(TAG, datum, PLACE, parse_err);  \
       continue;                                           \
   }

#define COMPARE_2_INTERVAL_LIST(TAG, PLACE)                \
   if (COMPARE(TAG)) {                                     \
       parse_2_interval_list(TAG, datum, PLACE, parse_err);\
       continue;                                           \
   }

/* c-basic-offset in emacs is set to 2 */
/* 
 * See read_boulder.h for description.
 */
int
read_boulder_record(FILE *file_input,
                    const int *strict_tags,
                    const int *io_version,
                    int   echo_output, /* should be echo_input */
                    const p3_file_type file_type,
                    p3_global_settings *pa, 
                    seq_args *sa, 
                    pr_append_str *glob_err,  /* Really should be called fatal_parse_err */
                    pr_append_str *nonfatal_parse_err,
		    pr_append_str *warnings,
                    read_boulder_record_results *res) 
{
  int line_len;
  int tag_len, datum_len;
  int data_found = 0;
  char *s, *n, *datum, *task_tmp = NULL;
  const char *p;
  pr_append_str *parse_err;
  pr_append_str *non_fatal_err;
  char *repeat_file_path = NULL, *int_repeat_file_path = NULL;
  int tmp_int;
  /* int min_3_prime = 0, min_5_prime = 0; Removed 10/20/2010 */
  int min_3_prime_distance_global = 0;    /* needed to check if both global and specific*/
  int min_3_prime_distance_specific = 0;  /* are given in same boulder record */
  int min_three_prime_distance;           /* holder for the value of this tag */
  double tmp_double;

  non_fatal_err = nonfatal_parse_err;

  while (((s = p3_read_line(file_input)) != NULL) && (strcmp(s,"="))) {
    /* If we are reading from a settings file, then skip every 
       line except those begining "PRIMER_" or "P3_FILE_ID".
       Hint: strncomp returns 0 if both strings are equal */
     if (file_type == settings 
	 && strncmp(s, "PRIMER_", 7) /* "s does not begin with 'PRIMER_'" */
	 && strncmp(s, "P3_FILE_ID", 10) /* "s does not begin with
					    'P3_FILE_ID'" */
	 ) {
      continue;
    }
    /* Silently ignore all primer3plus tags */
    /* Removed 10/20/10 because P3P_ tags
       are already ignored in settings files. */
    /* if (!(strncmp(s, "P3P_", 4))) {
      continue;
      } */

    data_found = 1;
    /* Print out the input */
    if (echo_output) printf("%s\n", s);
    line_len = strlen(s);
    /* If the line has an "=" read the tag in the right place */
    if ((n=strchr(s,'=')) == NULL) {
      /* The input line is illegal because it has no
       * "=" in it, but we still will read to the end
       * of the record. */
      pr_append_new_chunk(glob_err, "Input line with no '=': ");
      pr_append(glob_err, s);
      continue;
    } 
    /* Read in the new tags used from primer3 version 2.0 on */    
    else {
      /* Get the tag and the value pointers */
      tag_len = n - s;
      datum = n + 1;
      datum_len = line_len - tag_len - 1;
            
      /* Process "Sequence" (i.e. Per-Record) Arguments". */
      parse_err = non_fatal_err;
              
      /* COMPARE_AND_MALLOC("SEQUENCE", sa->sequence); */
      if (COMPARE("SEQUENCE_TEMPLATE")) {   /* NEW WAY */
        if (/* p3_get_seq_arg_sequence(sa) */ sa->sequence) {
          pr_append_new_chunk(parse_err, "Duplicate tag: ");
          pr_append(parse_err, "SEQUENCE_TEMPLATE"); 
        } else {
          if (p3_set_sa_sequence(sa, datum)) exit(-2);
        }
        continue;
      }
      if (COMPARE("SEQUENCE_QUALITY")) {
        if ((sa->n_quality = parse_seq_quality(datum, sa)) == 0) {
          pr_append_new_chunk(parse_err,
                              "Error in sequence quality data");
        }
        continue;
      }
      COMPARE_AND_MALLOC("SEQUENCE_ID", sa->sequence_name);
      COMPARE_AND_MALLOC("SEQUENCE_PRIMER", sa->left_input);
      COMPARE_AND_MALLOC("SEQUENCE_PRIMER_REVCOMP", sa->right_input);
      COMPARE_AND_MALLOC("SEQUENCE_INTERNAL_OLIGO", sa->internal_input);
      COMPARE_2_INTERVAL_LIST("SEQUENCE_PRIMER_PAIR_OK_REGION_LIST", &sa->ok_regions);
      COMPARE_INTERVAL_LIST("SEQUENCE_TARGET", &sa->tar2);
      COMPARE_INTERVAL_LIST("SEQUENCE_EXCLUDED_REGION", &sa->excl2);
      COMPARE_INTERVAL_LIST("SEQUENCE_INTERNAL_EXCLUDED_REGION",
                              &sa->excl_internal2);
      if (COMPARE("SEQUENCE_OVERLAP_JUNCTION_LIST")) {
	if (parse_intron_list(datum, sa->primer_overlap_junctions, 
			      &sa->primer_overlap_junctions_count) == 0) {
          pr_append_new_chunk(parse_err,
			      "Error in SEQUENCE_PRIMER_OVERLAP_JUNCTION_LIST");
        }
        continue;
      }
      if (COMPARE("SEQUENCE_INCLUDED_REGION")) {
        p = parse_int_pair("SEQUENCE_INCLUDED_REGION", datum, ',',
                           &sa->incl_s, &sa->incl_l, parse_err);
        if (NULL == p) /* An error; the message is already
                        * in parse_err.  */
          continue;
        while (' ' == *p || '\t' == *p) p++;
        if (*p != '\n' && *p != '\0')
            tag_syntax_error("SEQUENCE_INCLUDED_REGION", datum,
                             parse_err);
          continue;
      }
      COMPARE_INT("SEQUENCE_START_CODON_POSITION", sa->start_codon_pos);
      COMPARE_INT("SEQUENCE_FORCE_LEFT_START", sa->force_left_start);
      COMPARE_INT("SEQUENCE_FORCE_LEFT_END", sa->force_left_end);
      COMPARE_INT("SEQUENCE_FORCE_RIGHT_START", sa->force_right_start);
      COMPARE_INT("SEQUENCE_FORCE_RIGHT_END", sa->force_right_end);
      /* Process "Global" Arguments (those that persist between boulder
       * records). */
      parse_err = glob_err;  /* These errors are considered fatal. */
      if (COMPARE("PRIMER_PRODUCT_SIZE_RANGE")) {
        parse_product_size("PRIMER_PRODUCT_SIZE_RANGE", datum, pa,
                           parse_err);
        continue;
      }
      COMPARE_INT("PRIMER_OPT_SIZE", pa->p_args.opt_size);
      COMPARE_INT("PRIMER_MIN_SIZE", pa->p_args.min_size);
      COMPARE_INT("PRIMER_MAX_SIZE", pa->p_args.max_size);
      COMPARE_INT("PRIMER_MAX_POLY_X", pa->p_args.max_poly_x);
      COMPARE_FLOAT("PRIMER_OPT_TM", pa->p_args.opt_tm);
      COMPARE_FLOAT("PRIMER_OPT_GC_PERCENT", pa->p_args.opt_gc_content);
      COMPARE_FLOAT("PRIMER_MIN_TM", pa->p_args.min_tm);
      COMPARE_FLOAT("PRIMER_MAX_TM", pa->p_args.max_tm);
      COMPARE_FLOAT("PRIMER_PAIR_MAX_DIFF_TM", pa->max_diff_tm);
      if (COMPARE("PRIMER_TM_FORMULA")) {
          parse_int("PRIMER_TM_FORMULA", datum, &tmp_int, parse_err);
          pa->tm_santalucia = (tm_method_type) tmp_int;    /* added by T.Koressaar */
          continue;
      }
      if (COMPARE("PRIMER_SALT_CORRECTIONS")) {
        parse_int("PRIMER_SALT_CORRECTIONS", datum, &tmp_int, parse_err);
        pa->salt_corrections = (salt_correction_type) tmp_int; /* added by T.Koressaar */
        continue;
      }
      COMPARE_FLOAT("PRIMER_MIN_GC", pa->p_args.min_gc);
      COMPARE_FLOAT("PRIMER_MAX_GC", pa->p_args.max_gc);
      COMPARE_FLOAT("PRIMER_SALT_MONOVALENT", pa->p_args.salt_conc);
      COMPARE_FLOAT("PRIMER_SALT_DIVALENT", pa->p_args.divalent_conc);
      COMPARE_FLOAT("PRIMER_DNTP_CONC", pa->p_args.dntp_conc);
      COMPARE_FLOAT("PRIMER_DNA_CONC", pa->p_args.dna_conc);
      COMPARE_INT("PRIMER_MAX_NS_ACCEPTED", pa->p_args.num_ns_accepted);
      COMPARE_INT("PRIMER_PRODUCT_OPT_SIZE", pa->product_opt_size);
      COMPARE_FLOAT("PRIMER_MAX_SELF_ANY", pa->p_args.max_self_any);

      /* COMPARE_FLOAT("PRIMER_MAX_SELF_END", pa->p_args.max_self_end); */
      /* NEW */ COMPARE_FLOAT_USE_FN("PRIMER_MAX_SELF_END", p3_set_gs_primer_self_end)
	/* if (COMPARE("PRIMER_MAX_SELF_END")) {
	  parse_double("PRIMER_MAX_SELF_END", datum, &tmp_double, parse_err);
	  p3_set_gs_primer_self_end(pa, tmp_double);
	  continue;
	  } */

      COMPARE_FLOAT("PRIMER_MAX_SELF_ANY_TH", pa->p_args.max_self_any_th);
      COMPARE_FLOAT("PRIMER_MAX_SELF_END_TH", pa->p_args.max_self_end_th);
      COMPARE_FLOAT("PRIMER_MAX_HAIRPIN_TH", pa->p_args.max_hairpin_th);   
      COMPARE_FLOAT("PRIMER_PAIR_MAX_COMPL_ANY", pa->pair_compl_any);
      COMPARE_FLOAT("PRIMER_PAIR_MAX_COMPL_END", pa->pair_compl_end);
      COMPARE_FLOAT("PRIMER_PAIR_MAX_COMPL_ANY_TH", pa->pair_compl_any_th);
      COMPARE_FLOAT("PRIMER_PAIR_MAX_COMPL_END_TH", pa->pair_compl_end_th);
      COMPARE_INT("P3_FILE_FLAG", res->file_flag);
      COMPARE_INT("PRIMER_PICK_ANYWAY", pa->pick_anyway);
      COMPARE_INT("PRIMER_GC_CLAMP", pa->gc_clamp);
      COMPARE_INT("PRIMER_MAX_END_GC", pa->max_end_gc);
      COMPARE_INT("PRIMER_EXPLAIN_FLAG", res->explain_flag);
      COMPARE_INT("PRIMER_LIBERAL_BASE", pa->liberal_base);
      COMPARE_INT("PRIMER_FIRST_BASE_INDEX", pa->first_base_index);
      COMPARE_INT("PRIMER_NUM_RETURN", pa->num_return);
      COMPARE_INT("PRIMER_MIN_QUALITY", pa->p_args.min_quality);
      COMPARE_INT("PRIMER_MIN_END_QUALITY", pa->p_args.min_end_quality);
      if (COMPARE("PRIMER_MIN_THREE_PRIME_DISTANCE")) {
	parse_int("PRIMER_MIN_THREE_PRIME_DISTANCE", datum, &(min_three_prime_distance), parse_err);
	/* check if specific tag also specified - error in this case */
	if (min_3_prime_distance_specific == 1) {
	  pr_append_new_chunk(glob_err,
                              "Both PRIMER_MIN_THREE_PRIME_DISTANCE and PRIMER_{LEFT/RIGHT}_MIN_THREE_PRIME_DISTANCE specified");
	} else {
	  min_3_prime_distance_global = 1;
	  /* set up individual flags */
	  pa->min_left_three_prime_distance = min_three_prime_distance;
	  pa->min_right_three_prime_distance = min_three_prime_distance;
	}
	continue;
      }
      if (COMPARE("PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE")) {
	parse_int("PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE", datum, &(pa->min_left_three_prime_distance), parse_err);
	/* check if global tag also specified - error in this case */
	if (min_3_prime_distance_global == 1) {
	  pr_append_new_chunk(glob_err,
                              "Both PRIMER_MIN_THREE_PRIME_DISTANCE and PRIMER_{LEFT/RIGHT}_MIN_THREE_PRIME_DISTANCE specified");
	} else {
	  min_3_prime_distance_specific = 1;
	}
	continue;
      }
      if (COMPARE("PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE")) {
	parse_int("PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE", datum, &(pa->min_right_three_prime_distance), parse_err);
	/* check if global tag also specified - error in this case */
	if (min_3_prime_distance_global == 1) {
	  pr_append_new_chunk(glob_err,
                              "Both PRIMER_MIN_THREE_PRIME_DISTANCE and PRIMER_{LEFT/RIGHT}_MIN_THREE_PRIME_DISTANCE specified");
	} else {
	  min_3_prime_distance_specific = 1;
	}
	continue;
      }
      if (file_type == settings) {
        if (COMPARE("P3_FILE_ID")) continue;
      }
      COMPARE_INT("PRIMER_QUALITY_RANGE_MIN", pa->quality_range_min);
      COMPARE_INT("PRIMER_QUALITY_RANGE_MAX", pa->quality_range_max);
      COMPARE_FLOAT("PRIMER_PRODUCT_MAX_TM", pa->product_max_tm);
      COMPARE_FLOAT("PRIMER_PRODUCT_MIN_TM", pa->product_min_tm);
      COMPARE_FLOAT("PRIMER_PRODUCT_OPT_TM", pa->product_opt_tm);
      COMPARE_INT("PRIMER_SEQUENCING_LEAD", pa->sequencing.lead);
      COMPARE_INT("PRIMER_SEQUENCING_SPACING", pa->sequencing.spacing);
      COMPARE_INT("PRIMER_SEQUENCING_INTERVAL", pa->sequencing.interval);
      COMPARE_INT("PRIMER_SEQUENCING_ACCURACY", pa->sequencing.accuracy);
      if (COMPARE("PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION")) {
	parse_int("PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION", datum, &pa->min_5_prime_overlap_of_junction, parse_err);
	/* min_5_prime = 1; Removed 10/20/2010 */
	continue;
      }
      if (COMPARE("PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION")) {
	parse_int("PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION", datum, &pa->min_3_prime_overlap_of_junction, parse_err);
	/* min_3_prime = 1; Removed 10/20/2010 */
	continue;
      }
      COMPARE_AND_MALLOC("PRIMER_TASK", task_tmp);
      COMPARE_INT("PRIMER_PICK_RIGHT_PRIMER", pa->pick_right_primer);
      COMPARE_INT("PRIMER_PICK_INTERNAL_OLIGO", pa->pick_internal_oligo);
      COMPARE_INT("PRIMER_PICK_LEFT_PRIMER", pa->pick_left_primer);
      COMPARE_INT("PRIMER_INTERNAL_OPT_SIZE", pa->o_args.opt_size);
      COMPARE_INT("PRIMER_INTERNAL_MAX_SIZE", pa->o_args.max_size);
      COMPARE_INT("PRIMER_INTERNAL_MIN_SIZE", pa->o_args.min_size);
      COMPARE_INT("PRIMER_INTERNAL_MAX_POLY_X", pa->o_args.max_poly_x);
      COMPARE_FLOAT("PRIMER_INTERNAL_OPT_TM", pa->o_args.opt_tm);
      COMPARE_FLOAT("PRIMER_INTERNAL_OPT_GC_PERCENT",
                    pa->o_args.opt_gc_content);
      COMPARE_FLOAT("PRIMER_INTERNAL_MAX_TM", pa->o_args.max_tm);
      COMPARE_FLOAT("PRIMER_INTERNAL_MIN_TM", pa->o_args.min_tm);
      COMPARE_FLOAT("PRIMER_INTERNAL_MIN_GC", pa->o_args.min_gc);
      COMPARE_FLOAT("PRIMER_INTERNAL_MAX_GC", pa->o_args.max_gc);
      COMPARE_FLOAT("PRIMER_INTERNAL_SALT_MONOVALENT",
                        pa->o_args.salt_conc);
      COMPARE_FLOAT("PRIMER_INTERNAL_SALT_DIVALENT",
                    pa->o_args.divalent_conc);
      COMPARE_FLOAT("PRIMER_INTERNAL_DNTP_CONC",
                    pa->o_args.dntp_conc);
      COMPARE_FLOAT("PRIMER_INTERNAL_DNA_CONC", pa->o_args.dna_conc);
      COMPARE_INT("PRIMER_INTERNAL_MAX_NS_ACCEPTED", pa->o_args.num_ns_accepted);
      COMPARE_INT("PRIMER_INTERNAL_MIN_QUALITY", pa->o_args.min_quality);
      COMPARE_FLOAT("PRIMER_INTERNAL_MAX_SELF_ANY",
		    pa->o_args.max_self_any);

      /* COMPARE_FLOAT("PRIMER_INTERNAL_MAX_SELF_END", 
	 pa->o_args.max_self_end); */
      /* NEW  */ COMPARE_FLOAT_USE_FN("PRIMER_INTERNAL_MAX_SELF_END", 
				      p3_set_gs_primer_internal_oligo_self_end);
      /* if (COMPARE("PRIMER_INTERNAL_MAX_SELF_END")) { 
	 parse_double("PRIMER_INTERNAL_MAX_SELF_END", datum, &tmp_double, parse_err);
	 p3_set_gs_primer_internal_oligo_self_end(pa, tmp_double);
	 continue;
	 } */

      COMPARE_FLOAT("PRIMER_INTERNAL_MAX_SELF_ANY_TH",
			   pa->o_args.max_self_any_th);
      COMPARE_FLOAT("PRIMER_INTERNAL_MAX_SELF_END_TH",
				 pa->o_args.max_self_end_th);
      COMPARE_FLOAT("PRIMER_INTERNAL_MAX_HAIRPIN_TH",
			       pa->o_args.max_hairpin_th);
      COMPARE_FLOAT("PRIMER_MAX_LIBRARY_MISPRIMING",
                    pa->p_args.max_repeat_compl);
      COMPARE_FLOAT("PRIMER_INTERNAL_MAX_LIBRARY_MISHYB",
                    pa->o_args.max_repeat_compl);
      COMPARE_FLOAT("PRIMER_PAIR_MAX_LIBRARY_MISPRIMING",
                    pa->pair_repeat_compl);
      /* Mispriming / mishybing in the template. */
      COMPARE_FLOAT("PRIMER_MAX_TEMPLATE_MISPRIMING",
                    pa->p_args.max_template_mispriming);
      COMPARE_FLOAT("PRIMER_MAX_TEMPLATE_MISPRIMING_TH",
			  pa->p_args.max_template_mispriming_th);
      COMPARE_FLOAT("PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING",
                    pa->pair_max_template_mispriming);
      COMPARE_FLOAT("PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH",
		    pa->pair_max_template_mispriming_th);
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
          repeat_file_path = (char*) _rb_safe_malloc(strlen(datum) + 1);
          strcpy(repeat_file_path, datum);
        }
        continue;
      }
      if (COMPARE("PRIMER_INTERNAL_MISHYB_LIBRARY")) {
        if (int_repeat_file_path != NULL) {
          pr_append_new_chunk(glob_err,
                              "Duplicate PRIMER_INTERNAL_MISHYB_LIBRARY tag");
          free(int_repeat_file_path);
          int_repeat_file_path = NULL;
        } else {
          int_repeat_file_path = (char*) _rb_safe_malloc(strlen(datum) + 1);
          strcpy(int_repeat_file_path, datum);
        }
        continue;
      }
      if (COMPARE("P3_COMMENT")) continue;
      COMPARE_FLOAT("PRIMER_MAX_END_STABILITY", pa->max_end_stability);

      COMPARE_INT("PRIMER_LOWERCASE_MASKING",
                  pa->lowercase_masking); 
      /* added by T. Koressaar */
      COMPARE_INT("PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT", pa->thermodynamic_oligo_alignment);
      COMPARE_INT("PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT", pa->thermodynamic_template_alignment);
      if (COMPARE("PRIMER_THERMODYNAMIC_PARAMETERS_PATH")) {
        if (thermodynamic_params_path == NULL) {
          thermodynamic_params_path = (char*) _rb_safe_malloc(datum_len + 1);
          strcpy(thermodynamic_params_path, datum);
          thermodynamic_path_changed = 1;
        }
        /* check if path changes */
	else if (strcmp(thermodynamic_params_path, datum)) {
	  free(thermodynamic_params_path);
	  thermodynamic_params_path = (char*) _rb_safe_malloc(datum_len + 1); 
	  strcpy(thermodynamic_params_path, datum);
	  thermodynamic_path_changed = 1;
        }
        continue;
      }
      /* added by A. Untergasser */
      OVERWRITE_COMPARE_AND_MALLOC("PRIMER_MUST_MATCH_FIVE_PRIME", pa->p_args.must_match_five_prime);
      OVERWRITE_COMPARE_AND_MALLOC("PRIMER_MUST_MATCH_THREE_PRIME", pa->p_args.must_match_three_prime);
      OVERWRITE_COMPARE_AND_MALLOC("PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME", pa->o_args.must_match_five_prime);
      OVERWRITE_COMPARE_AND_MALLOC("PRIMER_INTERNAL_MUST_MATCH_THREE_PRIME", pa->o_args.must_match_three_prime);
      /* weights for objective functions  */
      /* CHANGE TEMP/temp -> TM/tm */
      COMPARE_FLOAT("PRIMER_WT_TM_GT", pa->p_args.weights.temp_gt);
      COMPARE_FLOAT("PRIMER_WT_TM_LT", pa->p_args.weights.temp_lt);
      COMPARE_FLOAT("PRIMER_WT_GC_PERCENT_GT", pa->p_args.weights.gc_content_gt);
      COMPARE_FLOAT("PRIMER_WT_GC_PERCENT_LT", pa->p_args.weights.gc_content_lt);
      COMPARE_FLOAT("PRIMER_WT_SIZE_LT", pa->p_args.weights.length_lt);
      COMPARE_FLOAT("PRIMER_WT_SIZE_GT", pa->p_args.weights.length_gt);
      COMPARE_FLOAT("PRIMER_WT_SELF_ANY", pa->p_args.weights.compl_any);
      COMPARE_FLOAT("PRIMER_WT_SELF_END", pa->p_args.weights.compl_end);
      COMPARE_FLOAT("PRIMER_WT_SELF_ANY_TH", pa->p_args.weights.compl_any_th);
      COMPARE_FLOAT("PRIMER_WT_SELF_END_TH", pa->p_args.weights.compl_end_th);
      COMPARE_FLOAT("PRIMER_WT_HAIRPIN_TH", pa->p_args.weights.hairpin_th);
      COMPARE_FLOAT("PRIMER_WT_NUM_NS", pa->p_args.weights.num_ns);
      COMPARE_FLOAT("PRIMER_WT_LIBRARY_MISPRIMING", pa->p_args.weights.repeat_sim);
      COMPARE_FLOAT("PRIMER_WT_SEQ_QUAL", pa->p_args.weights.seq_quality);
      COMPARE_FLOAT("PRIMER_WT_END_QUAL", pa->p_args.weights.end_quality);
      COMPARE_FLOAT("PRIMER_WT_POS_PENALTY", pa->p_args.weights.pos_penalty);
      COMPARE_FLOAT("PRIMER_WT_END_STABILITY",
                    pa->p_args.weights.end_stability);
      COMPARE_FLOAT("PRIMER_WT_TEMPLATE_MISPRIMING",
                    pa->p_args.weights.template_mispriming);
      COMPARE_FLOAT("PRIMER_WT_TEMPLATE_MISPRIMING_TH",
		                        pa->p_args.weights.template_mispriming_th);
      COMPARE_FLOAT("PRIMER_INTERNAL_WT_TM_GT", pa->o_args.weights.temp_gt);
      COMPARE_FLOAT("PRIMER_INTERNAL_WT_TM_LT", pa->o_args.weights.temp_lt);
      COMPARE_FLOAT("PRIMER_INTERNAL_WT_GC_PERCENT_GT", pa->o_args.weights.gc_content_gt);
      COMPARE_FLOAT("PRIMER_INTERNAL_WT_GC_PERCENT_LT", pa->o_args.weights.gc_content_lt);
      COMPARE_FLOAT("PRIMER_INTERNAL_WT_SIZE_LT", pa->o_args.weights.length_lt);
      COMPARE_FLOAT("PRIMER_INTERNAL_WT_SIZE_GT", pa->o_args.weights.length_gt);
      COMPARE_FLOAT("PRIMER_INTERNAL_WT_SELF_ANY", pa->o_args.weights.compl_any);
      COMPARE_FLOAT("PRIMER_INTERNAL_WT_SELF_END", pa->o_args.weights.compl_end);
      COMPARE_FLOAT("PRIMER_INTERNAL_WT_SELF_ANY_TH", pa->o_args.weights.compl_any_th);
      COMPARE_FLOAT("PRIMER_INTERNAL_WT_SELF_END_TH", pa->o_args.weights.compl_end_th);
      COMPARE_FLOAT("PRIMER_INTERNAL_WT_HAIRPIN_TH", pa->o_args.weights.hairpin_th);
      COMPARE_FLOAT("PRIMER_INTERNAL_WT_NUM_NS", pa->o_args.weights.num_ns);
      COMPARE_FLOAT("PRIMER_INTERNAL_WT_LIBRARY_MISHYB", pa->o_args.weights.repeat_sim);
      COMPARE_FLOAT("PRIMER_INTERNAL_WT_SEQ_QUAL", pa->o_args.weights.seq_quality);
      COMPARE_FLOAT("PRIMER_INTERNAL_WT_END_QUAL", pa->o_args.weights.end_quality);
      COMPARE_FLOAT("PRIMER_WT_TEMPLATE_MISPRIMING_TH",
		    pa->o_args.weights.template_mispriming_th);
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
      COMPARE_FLOAT("PRIMER_PAIR_WT_COMPL_ANY_TH",
		    pa->pr_pair_weights.compl_any_th);
      COMPARE_FLOAT("PRIMER_PAIR_WT_COMPL_END_TH",
		    pa->pr_pair_weights.compl_end_th);
      COMPARE_FLOAT("PRIMER_PAIR_WT_PRODUCT_TM_LT",
                    pa->pr_pair_weights.product_tm_lt);
      COMPARE_FLOAT("PRIMER_PAIR_WT_PRODUCT_TM_GT",
                    pa->pr_pair_weights.product_tm_gt);
      COMPARE_FLOAT("PRIMER_PAIR_WT_PRODUCT_SIZE_GT",
                    pa->pr_pair_weights.product_size_gt);
      COMPARE_FLOAT("PRIMER_PAIR_WT_PRODUCT_SIZE_LT",
                    pa->pr_pair_weights.product_size_lt);
      COMPARE_FLOAT("PRIMER_PAIR_WT_LIBRARY_MISPRIMING",
                    pa->pr_pair_weights.repeat_sim);
      COMPARE_FLOAT("PRIMER_PAIR_WT_TEMPLATE_MISPRIMING",
                    pa->pr_pair_weights.template_mispriming);
      COMPARE_FLOAT("PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH",
		    pa->pr_pair_weights.template_mispriming_th);
    }
    /* End of reading the tags in the right place */
        
    /*  Complain about unrecognized tags */
     if (*strict_tags == 1) {
	pr_append_new_chunk(glob_err, "Unrecognized tag: ");
	pr_append(glob_err, s);
	fprintf(stderr, "Unrecognized tag: %s\n", s);
     }
  }  /* while ((s = p3_read_line(stdin)) != NULL && strcmp(s,"=")) { */

  /* Check if the record was terminated by "=" */
  if (NULL == s) { /* End of file. */
    if (data_found) {
       pr_append_new_chunk(glob_err, 
                          "Final record not terminated by '='");
      return 1;
    } else  return 0;
  }
    
  /* Figure out the right settings for the tasks*/
  if (task_tmp != NULL) {
    if (!strcmp_nocase(task_tmp, "pick_pcr_primers")) {
      pa->primer_task = generic;
      pa->pick_left_primer = 1;
      pa->pick_right_primer = 1;
      pa->pick_internal_oligo = 0;
    } else if (!strcmp_nocase(task_tmp, "pick_pcr_primers_and_hyb_probe")) {
      pa->primer_task = generic; 
      pa->pick_left_primer = 1;
      pa->pick_right_primer = 1;
      pa->pick_internal_oligo = 1;
    } else if (!strcmp_nocase(task_tmp, "pick_left_only")) {
      pa->primer_task = generic;
      pa->pick_left_primer = 1;
      pa->pick_right_primer = 0;
      pa->pick_internal_oligo = 0;
    } else if (!strcmp_nocase(task_tmp, "pick_right_only")) {
      pa->primer_task = generic;
      pa->pick_left_primer = 0;
      pa->pick_right_primer = 1;
      pa->pick_internal_oligo = 0;
    } else if (!strcmp_nocase(task_tmp, "pick_hyb_probe_only")) {
      pa->primer_task = generic;
      pa->pick_left_primer = 0;
      pa->pick_right_primer = 0;
      pa->pick_internal_oligo = 1;
    } else if (!strcmp_nocase(task_tmp, "generic")) {
      pa->primer_task = generic;
    } else if (!strcmp_nocase(task_tmp, "pick_detection_primers")) {
      pa->primer_task = generic; /* Deliberate duplication for
				    backward compatibility. */
    } else if (!strcmp_nocase(task_tmp, "pick_cloning_primers")) {
      pa->primer_task = pick_cloning_primers;
    } else if (!strcmp_nocase(task_tmp, "pick_discriminative_primers")) {
      pa->primer_task = pick_discriminative_primers;
    } else if (!strcmp_nocase(task_tmp, "pick_sequencing_primers")) {
      pa->primer_task = pick_sequencing_primers;
    } else if (!strcmp_nocase(task_tmp, "pick_primer_list")) {
      pa->primer_task = pick_primer_list;
    } else if (!strcmp_nocase(task_tmp, "check_primers")) {
      pa->primer_task = check_primers;
      /* check_primers sets the picking flags itself */
      pa->pick_left_primer = 0;
      pa->pick_right_primer = 0;
      pa->pick_internal_oligo = 0;
      if (sa->left_input){
        pa->pick_left_primer = 1;
      }
      if (sa->right_input){
        pa->pick_right_primer = 1;
      }
      if (sa->internal_input){
        pa->pick_internal_oligo = 1;
      }
    } else 
      pr_append_new_chunk(glob_err,
                          "Unrecognized PRIMER_TASK");
    free(task_tmp);
  }

  /* 
   * WARNING: read_and_create_seq_lib uses p3_read_line, so repeat
   * library files cannot be read inside the 
   * while ((s = p3_read_line(stdin))...) loop above.
   *
   * FIX ME, in fact the reading
   * of the library contents probably belongs inside
   * primer3_boulder_main.c or libprimer3.c.
   */

  /* Read in the repeat libraries */
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

  /* Read in the repeat libraries for internal oligo */
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

  if (pa->primer_task == pick_pcr_primers_and_hyb_probe) {
    PR_ASSERT(pa->pick_internal_oligo);
  }

  /* Removed 10/20/2010 as excessively compulsive, especially in the context
     of taking input from web pages.
     if ((min_3_prime || min_5_prime) && (sa->primer_overlap_junctions_count == 0)) {
     pr_append_new_chunk(warnings,
                         "SEQUENCE_OVERLAP_JUNCTION_LIST not given, but "
			 "PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION or "
			 "PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION specified");
			 } */

  return 1;
}
#undef COMPARE
#undef COMPARE_AND_MALLOC
#undef COMPARE_INT
#undef COMPARE_FLOAT
#undef COMPARE_INTERVAL_LIST

static void
pr_append2(pr_append_str *err,
	   const char* s1, 
	   const char* s2) {
    pr_append_new_chunk(err, s1);
    pr_append(err, s2);
}  

int 
read_p3_file(const char *file_name,
	     const p3_file_type expected_file_type,
	     int echo_output,
	     int strict_tags,
	     p3_global_settings *pa, 
	     seq_args *sa,
	     pr_append_str *fatal_err,
	     pr_append_str *nonfatal_err,
	     pr_append_str *warnings,
	     read_boulder_record_results *read_boulder_record_res) 
{
  FILE *file;
  int ret_par = 0;
  int io_version = 4;
  char *line1;
  char *line2;
  char *line3;
  p3_file_type file_type = all_parameters;
    
  /* Check if a file name was provided */
  PR_ASSERT(NULL != file_name);

  /* Open the file */
  if (!(file = fopen(file_name,"r"))) {
    pr_append2(fatal_err, "Cannot open ", file_name);
    return ret_par;
  }

  /* Parse and interpret the first 3 lines */

  /* Line 1 */
  line1 = p3_read_line(file);
  if (!line1) {
    pr_append2(fatal_err, 
	       "Settings file is empty: ", 
	       file_name);
    return ret_par;
  }
  if (strcmp(line1,
	     "Primer3 File - http://primer3.sourceforge.net")) {
      pr_append2(fatal_err,
		 "First line must be \"Primer3 File - http://primer3.sourceforge.net\" in ",
		 file_name);
      return ret_par;
  }

  /* Line 2 */
  line2 = p3_read_line(file);
  if (!line2) {
    pr_append2(fatal_err, 
	       "Incorrect file format (too few lines) in ", 
	       file_name);
    return ret_par;
  }
  if (!strcmp(line2,"P3_FILE_TYPE=all_parameters")) {
    file_type = all_parameters;
  } else if (!strcmp(line2,"P3_FILE_TYPE=sequence")) {
    file_type = sequence;
  } else if (!strcmp(line2,"P3_FILE_TYPE=settings")) {
    file_type = settings;
  } else {
    pr_append2(fatal_err, 
	       "Unknown file type in at line 2 (line2='", line2);
    pr_append(fatal_err, "') in ");
    pr_append(fatal_err, file_name);
    return ret_par;
  }
  if (echo_output) {
    printf("P3_SETTINGS_FILE_USED=%s\n", file_name);
    printf("%s\n", line2);
  }

  /* Line 3 */
  line3 = p3_read_line(file);
  if (!line3) {
    pr_append2(fatal_err, 
	       "Incorrect file format (too few lines) in ", 
	       file_name);
      return ret_par;
  }
  if (strcmp(line3, "")) {
    pr_append2(fatal_err, "Line 3 must be empty in ",
	       file_name);
    return ret_par;
  }

  /* Check if the file type matches the expected type */
  if (file_type != expected_file_type){
    pr_append_new_chunk(nonfatal_err, 
			"Unexpected P3 file type parsed");
  }

  /* read the file */
  ret_par = read_boulder_record(file, &strict_tags, &io_version, 
				echo_output, expected_file_type,
				pa, sa, fatal_err, 
				nonfatal_err, warnings, 
				read_boulder_record_res);

  if (echo_output) printf("P3_SETTINGS_FILE_END=\n");
  if (file) fclose(file);
         
  return ret_par;
}

static void
tag_syntax_error(const char *tag_name, const char *datum,
                 pr_append_str *err)
{
    pr_append_new_chunk(err, "Illegal ");
    pr_append(err, tag_name);
    pr_append(err, " value: ");
    pr_append(err, datum);
}

static void
parse_double(const char *tag_name, const char *datum,
             double *out, pr_append_str *err)
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
parse_int(const char *tag_name, const char *datum,
          int *out, pr_append_str *err)
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
parse_int_pair(const char    *tag_name, const char *datum,
               char          sep,              /* The separator, e.g. ',' or '-'. */
               int           *out1, int *out2, /* The 2 integers. */
               pr_append_str *err)             /* Error messages. */
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
      pr_append(err, tag_name); return;
    }
  }
}

/*
 * For correct input, return a pointer to the first non-tab, non-space
 * character after the forth integer and after the separator sep2, and
 * place the integers in out1, out2, out3 and out4. On incorrect input,
 * return NULL; If any of the 4 integers is not specified, the
 * corresponding output value will be -1.
 */
static char *
parse_2_int_pair(const char    *tag_name, char *datum,
		 char          sep,              /* The separator between 2 numbers, e.g. ',' or '-'. */
		 char          sep2,             /* Separator between 2 intervals */
		 int           *out1, int *out2, 
		 int           *out3, int *out4, /* The 4 integers. */
		 pr_append_str *err)             /* Error messages. */
{
    char *nptr, *tmp;
    long tlong;

    nptr = datum;

    if (*nptr == sep) {
      *out1 = -1;
    } else {
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
    }
    while (' ' == *nptr || '\t' == *nptr) nptr++;
    if (sep != *nptr) {
        tag_syntax_error(tag_name, datum, err);
        return NULL;
    }
    nptr++; /* Advance past separator. */
    while (' ' == *nptr || '\t' == *nptr) nptr++;
    tmp = nptr;
    if (*nptr == sep) {
      *out2 = -1;
    } else {
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
    }
    while (' ' == *nptr || '\t' == *nptr) nptr++;
    if (sep != *nptr) {
        tag_syntax_error(tag_name, datum, err);
        return NULL;
    }
    nptr++; /* Advance past separator. */
    while (' ' == *nptr || '\t' == *nptr) nptr++;
    tmp = nptr;
    if (*nptr == sep) {
      *out3 = -1;
    } else {
      tlong = strtol(tmp, &nptr, 10);
      if (tlong > INT_MAX || tlong < INT_MIN) {
        tag_syntax_error(tag_name, datum, err);
        pr_append(err, " (value too large or too small)");
        return NULL;
      }
      *out3 = tlong;
      if (nptr == tmp) {
	tag_syntax_error(tag_name, datum, err);
        return NULL;
      }
    }
    while (' ' == *nptr || '\t' == *nptr) nptr++;
    if (sep != *nptr) {
        tag_syntax_error(tag_name, datum, err);
        return NULL;
    }
    nptr++; /* Advance past separator. */
    while (' ' == *nptr || '\t' == *nptr) nptr++;
    tmp = nptr;
    if ((*nptr == '\n') || (*nptr == '\0') || (*nptr == sep2)) {
      *out4 = -1;
    } else {
      tlong = strtol(tmp, &nptr, 10);
      if (tlong > INT_MAX || tlong < INT_MIN) {
        tag_syntax_error(tag_name, datum, err);
        pr_append(err, " (value too large or too small)");
        return NULL;
      }
      *out4 = tlong;
      if (nptr == tmp) {
	tag_syntax_error(tag_name, datum, err);
        return NULL;
      }
    }
    while (' ' == *nptr || '\t' == *nptr) nptr++;
    /* Should have reached an interval separator or end of line/string */
    if ((*nptr != sep2) && (*nptr != '\0') && (*nptr != '\n')) {
      tag_syntax_error(tag_name, datum, err);
      return NULL;
    }
    if (*nptr == sep2) {
      /* Advance past the interval separator */
      nptr++;
    }
    while (' ' == *nptr || '\t' == *nptr) nptr++;

    return nptr;
}

static void
parse_2_interval_list(const char *tag_name,
		      char *datum,
		      interval_array_t4 *interval_arr,
		      pr_append_str *err)
{
  char *p = datum;
  int i1, i2, i3, i4;
  int ret = 0;

  interval_arr->count = 0;
  interval_arr->any_pair = 0;
  interval_arr->any_left = interval_arr->any_right = 0;

  while (' ' == *p || '\t' == *p) p++;
  while (*p != '\0' && *p != '\n') {
    p = parse_2_int_pair(tag_name, p, ',', ';', &i1, &i2, &i3, &i4, err);
    if (NULL == p) return;
    ret = p3_add_to_2_interval_array(interval_arr, i1, i2, i3, i4);
    if (ret == 1) {
      pr_append_new_chunk(err, "Too many elements for tag ");
      pr_append(err, tag_name); return;
    } else if (ret == 2) {
      pr_append_new_chunk(err, "Invalid range at tag  ");
      pr_append(err, tag_name); return;
    }
  }
}

static int
parse_intron_list(char *s,
		  int *list,
		  int *count) 
{
  long t;
  char *p, *q;

  *count = 0;

  p = q = s;

  while (*q != '\0' && *q != '\n') {
    t = strtol(p, &q, 10);
    if (q == p) {
      while (*q != '\0') {
        if (!isspace(*q)) {
	  *count = 0;
          return 0; 
        }
        q++;
      }
      return *count;
    }
    if (t > INT_MAX || t < INT_MIN) {
      return 0;
    }
    list[*count] = t;
    (*count)++;

    p = q;
  }
  return *count;
}

static void
parse_product_size(const char *tag_name, char *in,
                   p3_global_settings *pa,
                   pr_append_str *err)
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

/* Would be worthwhile to test this a bit more off-by-one errors on
   length of the quality score */

/* This function returns the number of elements in
   the sequence quality vector, and updates
   sargs->quality,  sargs->n_quality,
   and sargs->quality_storage_size */
static int
parse_seq_quality(char *s,
                  seq_args *sargs) 
{
  long t;
  char *p, *q;

  p3_set_sa_empty_quality(sargs);

  p = q = s;

  while (*q != '\0' && *q != '\n') {
    t = strtol(p, &q, 10);
    if (q == p) {
      while (*q != '\0') {
        if (!isspace(*q)) {
          p3_set_sa_empty_quality(sargs);
          return 0; 
        }
        q++;
      }
      return sargs->n_quality;
    }
    p3_sa_add_to_quality_array(sargs, t);

    p = q;
  }
  return sargs->n_quality;
}

/* =========================================================== */
/* Fail-stop wrapper for memory allocation.                    */
/* =========================================================== */
/* 
 * Panic messages for when the program runs out of memory.
 */

static void *
_rb_safe_malloc(size_t x)
{
    void *r = malloc(x);
    if (NULL == r)
            out_of_memory_error();
    return r;
}

static void
pr_append(pr_append_str *x,
          const char *s)
{
  if (pr_append_external(x, s))
          out_of_memory_error();
}

static void
pr_append_new_chunk(pr_append_str *x,
                    const char *s)
{
  if (pr_append_new_chunk_external(x, s))
    out_of_memory_error();
}

static void
out_of_memory_error() 
{
  fprintf(stderr, "out of memory in read_boulder\n");
  exit(-2);
}

/* End of fail-stop wrappers. */
/* =========================================================== */
