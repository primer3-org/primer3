/*
 * Copyright (c) 1996, Whitehead Institute for Biomedical Research. All rights
 * reserved.  Please see full software use agreement in primer3_main.c or by
 * executing primer3 with -h.
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
static void   parse_align_score(primer_error *, const char *, const char *,
				short *, pr_append_str *);
static void   parse_double(primer_error *, const char *, const char *,
			   double *, pr_append_str *);
static void   parse_int(primer_error *, const char *, const char *, int *,
			pr_append_str *);
static const char *parse_int_pair(primer_error *,
				  const char *, const char *, char,
				  int *, int *, pr_append_str *);
static void   parse_interval_list(primer_error *, const char *, const char *,
				  int*, interval_array_t, pr_append_str *);
static void   parse_product_size(primer_error *,
				 const char *, char *, primer_args *,
				 pr_append_str *);
static void   tag_syntax_error(primer_error *, const char *, const char *,  pr_append_str *);
static void   read_seq_lib(primer_error *, seq_lib *, const char *, const char *);
static short  upcase_and_check_char(char *);
static char*  read_line(primer_error *, FILE *);
static double parse_seq_name(char *);
static void   reverse_complement_seq_lib(primer_error *, seq_lib *);
static int    parse_seq_quality(char *, int **);


/* 
 * Hack to support old SunOS headers.  (We do not try to declare _all_
 * undeclared functions; only those with non-int return types.)
 */
#ifndef __cplusplus
extern double strtod();
#endif

/*
 * ==========================================================================
 * External APIs
 * ==========================================================================
 */

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
           jump_append_new_chunk(perr, parse_err,  \
                               "Duplicate tag: "); \
           jump_append(perr, parse_err, TAG);      \
       } else {                                    \
           T = pr_jump_malloc(perr, datum_len + 1);\
           strcpy(T, datum);                       \
       }                                           \
       continue;                                   \
   }

#define COMPARE_ALIGN_SCORE(TAG,T)                           \
   if (COMPARE(TAG)) {                                       \
       parse_align_score(perr, TAG, datum, &(T), parse_err); \
       continue;                                             \
   }

#define COMPARE_FLOAT(TAG,T)                            \
   if (COMPARE(TAG)) {                                  \
       parse_double(perr, TAG, datum, &(T), parse_err); \
       continue;                                        \
   }

#define COMPARE_INT(TAG,T)                           \
   if (COMPARE(TAG)) {                               \
       parse_int(perr, TAG, datum, &(T), parse_err); \
       continue;                                     \
   }

#define COMPARE_INTERVAL_LIST(TAG, SIZE, LIST      )                   \
   if (COMPARE(TAG)) {                                                 \
       parse_interval_list(perr, TAG, datum, &SIZE, LIST, parse_err);  \
       continue;                                                       \
   }

int
read_record(const program_args *prog_args,
	    primer_args *pa,
	    seq_args *sa)
{ 
    int line_len, seq_len, n_quality;
    int tag_len, datum_len;
    int data_found = 0;
    int pick_internal_oligo = 2;
    char *s, *n, *datum, *task_tmp = NULL;
    const char *p;
    pr_append_str *parse_err;
    char *repeat_file = NULL, *int_repeat_file = NULL;
    primer_error errstruct;
    primer_error *perr = &errstruct;

    memset(&sa->error, 0, sizeof(sa->error));
    memset(&pa->glob_err, 0, sizeof(pa->glob_err));
    memset(sa, 0, sizeof(*sa));
    sa->start_codon_pos = PR_DEFAULT_START_CODON_POS;
    sa->incl_l = -1; /* Indicates logical NULL. */
    n_quality = 0;

    /* non-local error recovery */
    perr->system_errno = 0;
    perr->local_errno = 0;
    perr->error_msg = NULL;
    if (setjmp(perr->jmpenv) != 0)
	return 0;

    while ((s = read_line(perr, stdin)) != NULL && strcmp(s,"=")) {
	data_found = 1;
	if (0 == prog_args->format_output) printf("%s\n", s);
	line_len = strlen(s);
	if ((n=strchr(s,'=')) == NULL) {
	    /* 
	     * The input line is illegal, but we still have to read to the end
	     * of the record.
	     */
	    jump_append_new_chunk(perr, &pa->glob_err, "Input line with no '=': ");
	    jump_append(perr, &pa->glob_err, s);
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
		   jump_append_new_chunk(perr, &sa->error, 
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
		p = parse_int_pair(perr, "INCLUDED_REGION", datum, ',',
				   &sa->incl_s, &sa->incl_l, parse_err);
		if (NULL == p) /* 
                                * An error; the message is already
                                * in parse_err.
                                */
		  continue;

		while (' ' == *p || '\t' == *p) p++;
		if (*p != '\n' && *p != '\0')
		    tag_syntax_error(perr, "INCLUDED_REGION", datum,
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
		parse_product_size(perr, "PRIMER_PRODUCT_SIZE_RANGE", datum, pa,
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
	    COMPARE_FLOAT("PRIMER_MIN_GC", pa->min_gc);
	    COMPARE_FLOAT("PRIMER_MAX_GC", pa->max_gc);
	    COMPARE_FLOAT("PRIMER_SALT_CONC", pa->salt_conc);
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
            COMPARE_FLOAT("PRIMER_INTERNAL_OLIGO_DNA_CONC", pa->io_dna_conc);
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
	    COMPARE_FLOAT("PRIMER_INSIDE_PENALTY", pa->inside_penalty);
	    COMPARE_FLOAT("PRIMER_OUTSIDE_PENALTY", pa->outside_penalty);
            if (COMPARE("PRIMER_MISPRIMING_LIBRARY")) {
		if (repeat_file != NULL) {
		    jump_append_new_chunk(perr, &pa->glob_err,
					  "Duplicate PRIMER_MISPRIMING_LIBRARY"
					  " tag");
		    free(repeat_file);
		    repeat_file = NULL;
		} else {
		    repeat_file = pr_jump_malloc(perr, strlen(datum) + 1);
		    strcpy(repeat_file, datum);
		}
		continue;
	    }
            if (COMPARE("PRIMER_INTERNAL_OLIGO_MISHYB_LIBRARY")) {
		if (int_repeat_file != NULL) {
		    jump_append_new_chunk(perr, &pa->glob_err,
					  "Duplicate PRIMER_INTERNAL_OLIGO_"
					  "MISHYB_LIBRARY tag");
		    free(int_repeat_file);
		    int_repeat_file = NULL;
		} else {
		    int_repeat_file = pr_jump_malloc(perr, strlen(datum) + 1);
		    strcpy(int_repeat_file, datum);
		}
		continue;
	    }
	    if (COMPARE("PRIMER_COMMENT") || COMPARE("COMMENT")) continue;
	    COMPARE_FLOAT("PRIMER_MAX_END_STABILITY", pa->max_end_stability);

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

            /* CHANGE Add PRIMER_PAIR_WT_PRODUCT_TM. */

	}
	if (1 == prog_args->strict_tags) {
	    jump_append_new_chunk(perr, &pa->glob_err, "Unrecognized tag: ");
	    jump_append(perr, &pa->glob_err, s);
	    fprintf(stderr, "Unrecognized tag: %s\n", s);
	}
    }

    if (NULL == s) { /* End of file. */
	if (data_found) {
	    jump_append_new_chunk(perr, &pa->glob_err, 
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
          else   jump_append_new_chunk(perr, &pa->glob_err,
				       "Unrecognized PRIMER_TASK");
	  free(task_tmp);
    }

    if (NULL == sa->sequence)
	jump_append_new_chunk(perr, &sa->error, "Missing SEQUENCE tag");
    else {
	seq_len = strlen(sa->sequence);
	if (sa->incl_l == -1) {
	     sa->incl_l = seq_len;
	     sa->incl_s = pa->first_base_index;
	}
	if(n_quality !=0 && n_quality != seq_len)
	    jump_append_new_chunk(perr, &sa->error,
				  "Error in sequence quality data");
        if((pa->min_quality != 0 || pa->io_min_quality != 0) && n_quality == 0) 
	    jump_append_new_chunk(perr, &sa->error,
				  "Sequence quality data missing");
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
	read_seq_lib(perr, &pa->repeat_lib, repeat_file, "mispriming library");
	if(pa->repeat_lib.error.data != NULL) {
	    jump_append_new_chunk(perr, &pa->glob_err,
				  pa->repeat_lib.error.data);
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
	read_seq_lib(perr, &pa->io_mishyb_library, int_repeat_file,
		     "internal oligo mishyb library");
	if(pa->io_mishyb_library.error.data != NULL) {
	    jump_append_new_chunk(perr, &pa->glob_err,
				  pa->io_mishyb_library.error.data);
	}
      }
      free(int_repeat_file);
      int_repeat_file = NULL;
    }
    
    if((pick_internal_oligo == 1 || pick_internal_oligo == 0) &&
       (pa->primer_task == pick_left_only || 
	pa->primer_task == pick_right_only ||
	pa->primer_task == pick_hyb_probe_only)) 
	jump_append_new_chunk(perr, &pa->glob_err, 
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

void
free_record(seq_args *sa) 
{ 
    if (NULL != sa->internal_input) free(sa->internal_input);
    if (NULL != sa->left_input) free(sa->left_input);
    if (NULL != sa->right_input) free(sa->right_input);
    if (NULL != sa->sequence) free(sa->sequence);
    if (NULL != sa->quality)  free(sa->quality);
    if (NULL != sa->trimmed_seq) free(sa->trimmed_seq);
    if (NULL != sa->sequence_name) free(sa->sequence_name);
    if (NULL != sa->error.data) free(sa->error.data);
    if (NULL != sa->warning.data) free(sa->warning.data);
}

/* 
 * Free exogenous storage associated with a seq_lib (but not the seq_lib
 * itself).  Silently ignore NULL p.  Set *p to 0 bytes.
 */
void
free_seq_lib(seq_lib *p)
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
 * ==========================================================================
 * Internal functions - mostly children of the read_record function.
 * These should all be static.
 * ==========================================================================
 */

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
read_line(perr, file)
     primer_error *perr;
     FILE *file;
{
    static size_t ssz;
    static char *s = NULL;

    size_t remaining_size;
    char *p, *n;

    if (NULL == s) {
	ssz = INIT_BUF_SIZE;
	s = pr_jump_malloc(perr, ssz);
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
	s = pr_jump_realloc(perr, s, ssz);
	p = strchr(s, '\0');
	remaining_size = ssz - (p - s);
    }
}

static void
tag_syntax_error(perr, tag_name, datum, err)
    primer_error *perr;
    const char *tag_name, *datum;
    pr_append_str *err;
{
    jump_append_new_chunk(perr, err, "Illegal ");
    jump_append(perr, err, tag_name);
    jump_append(perr, err, " value: ");
    jump_append(perr, err, datum);
}

static void
parse_align_score(perr, tag_name, datum, out, err)
    primer_error *perr;
    const char *datum, *tag_name;
    short *out;
    pr_append_str *err;
{
    double d;

    parse_double(perr, tag_name, datum, &d, err);
    d *= PR_ALIGN_SCORE_PRECISION;
    if (d > SHRT_MAX) {
	jump_append_new_chunk(perr, err, "Value too large at tag ");
	jump_append(perr, err, tag_name);
    } else {
	/* Should we be rounding here? */
	*out = (short)d;
    }
}    

static void
parse_double(perr, tag_name, datum, out, err)
    primer_error *perr;
    const char *datum, *tag_name;
    double *out;
    pr_append_str *err;
{
    char *nptr;
    *out = strtod(datum, &nptr);
    if (nptr == datum) {
	/* Empty string or complete junk. */
	tag_syntax_error(perr, tag_name, datum, err);
	*out = 0.0;
	return;
    }
    /* Look for trailing junk. */
    while (*nptr != '\n' && *nptr != '\0') {
	if (*nptr != ' ' && *nptr != '\t') {
	    tag_syntax_error(perr, tag_name, datum, err);
	    break;
	}
	nptr++;
    }
}

static void
parse_int(perr, tag_name, datum, out, err)
    primer_error *perr;
    const char *datum, *tag_name;
    int *out;
    pr_append_str *err;
{
    char *nptr;
    long tlong;
    tlong = strtol(datum, &nptr, 10);
    if (tlong > INT_MAX || tlong < INT_MIN) {
	tag_syntax_error(perr, tag_name, datum, err);
	jump_append(perr, err, " (value too large or too small)");
	return;
    }
    *out = tlong;
    if (nptr == datum) {
	/* Empty string or complete junk. */
	tag_syntax_error(perr, tag_name, datum, err);
	return;
    }
    /* Look for trailing junk. */
    while (*nptr != '\n' && *nptr != '\0') {
	if (*nptr != ' ' && *nptr != '\t') {
	    tag_syntax_error(perr, tag_name, datum, err);
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
parse_int_pair(primer_error  *perr,
	       const char    *tag_name,
	       const char    *datum,
	       char           sep,	/* The separator, e.g. ',' or '-'. */
	       int           *out1,	/* The 2 integers. */
	       int           *out2,
	       pr_append_str *err)	/* Error messages. */
{
    char *nptr, *tmp;
    long tlong;
    tlong = strtol(datum, &nptr, 10);
    if (tlong > INT_MAX || tlong < INT_MIN) {
	tag_syntax_error(perr, tag_name, datum, err);
	jump_append(perr, err, " (value too large or too small)");
	return NULL;
    }
    *out1 = tlong;
    if (nptr == datum) {
	tag_syntax_error(perr, tag_name, datum, err);
	return NULL;
    }
    while (' ' == *nptr || '\t' == *nptr) nptr++;
    if (sep != *nptr) {
	tag_syntax_error(perr, tag_name, datum, err);
	return NULL;
    }
    nptr++; /* Advance past separator. */
    while (' ' == *nptr || '\t' == *nptr) nptr++;
    tmp = nptr;
    tlong = strtol(tmp, &nptr, 10);
    if (tlong > INT_MAX || tlong < INT_MIN) {
	tag_syntax_error(perr, tag_name, datum, err);
	jump_append(perr, err, " (value too large or too small)");
	return NULL;
    }
    *out2 = tlong;
    if (nptr == tmp) {
	tag_syntax_error(perr, tag_name, datum, err);
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
parse_interval_list(perr, tag_name, datum, count, interval_array, err)
    primer_error *perr;
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
	    jump_append_new_chunk(perr, err, "Too many elements for tag ");
	    jump_append(perr, err, tag_name);
	    return;
	}
	p = parse_int_pair(perr, tag_name, p, ',', 
			   &interval_array[*count][0],
			   &interval_array[*count][1],
			   err);
	if (NULL == p) return;
	(*count)++;
    }
}

static void
parse_product_size(perr, tag_name, in, pa, err)
    primer_error *perr;
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
	jump_append_new_chunk(perr, err, tag_name);
	jump_append(perr, err, " begins but does not end with a quote");
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
	    jump_append_new_chunk(perr, err, "Too many values for ");
	    jump_append(perr, err, tag_name);
	    return;
	}
	p = parse_int_pair(perr, tag_name, p, '-',
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
read_seq_lib(perr, lib, filename, errfrag)
    primer_error *perr;
    seq_lib *lib;
    const char *filename;
    const char *errfrag;
{
    char  *p;
    FILE *file;
    int i, m, k;
    size_t j, n;

    PR_ASSERT(NULL != lib);
    PR_ASSERT(NULL != filename);

    free_seq_lib(lib);

    lib->repeat_file = pr_jump_malloc(perr, strlen(filename) + 1);
    strcpy(lib->repeat_file, filename);

    if((file = fopen(lib->repeat_file,"r")) == NULL) {
	jump_append_new_chunk(perr, &lib->error,
			    "Cannot open ");
	goto ERROR;
    }

    j = INIT_BUF_SIZE;
    n = INIT_LIB_SIZE;
    lib->names = pr_jump_malloc(perr, INIT_LIB_SIZE*sizeof(*lib->names));
    lib->seqs  = pr_jump_malloc(perr, INIT_LIB_SIZE*sizeof(*lib->seqs));
    lib->weight= pr_jump_malloc(perr, INIT_LIB_SIZE*sizeof(*lib->weight));
    lib->seq_num = 0;

    i = -1;  m = 0; k = 0;
    while((p = read_line(perr, file))) {
	if(*p == '>'){
	    i++;
	    if(i >= n) {
		n += INIT_LIB_SIZE;
		lib->names = pr_jump_realloc(perr,
					     lib->names,n*sizeof(*lib->names));
		lib->seqs  = pr_jump_realloc(perr,
					     lib->seqs ,n*sizeof(*lib->seqs));
		lib->weight= pr_jump_realloc(perr,
					     lib->weight,
					     n*sizeof(*lib->weight));
	    }
	    p++;
	    lib->names[i] = pr_jump_malloc(perr, strlen(p) + 1);
	    strcpy(lib->names[i],p);
	    lib->weight[i] = parse_seq_name(lib->names[i]);
	    lib->seqs[i] = pr_jump_malloc(perr, INIT_BUF_SIZE);
	    lib->seqs[i][0] = '\0';
	    lib->seq_num = i+1;
	    if(lib->weight[i] < 0) {
		jump_append_new_chunk(perr, &lib->error, "Illegal weight in ");
		goto ERROR;
	    }
	    j = INIT_BUF_SIZE;
	    k = 0;
	    if(i > 0) {
		/* We are actually testing the previous sequence. */
		if(strlen(lib->seqs[i-1]) == 0) {
		    jump_append_new_chunk(perr, &lib->error, "Empty sequence in ");
		    goto ERROR;
		}
		m += upcase_and_check_char(lib->seqs[i-1]);
	    }
	    p--;
	}
	else {
	    if(i < 0){ 
		jump_append_new_chunk(perr, &lib->error,
				    "Missing id line (expected '>') in ");
		goto ERROR;
	    } else {
		if(k+strlen(p) > j-2){
		    while(j-2 < k+ strlen(p))j += INIT_BUF_SIZE;
		    lib->seqs[i] = pr_jump_realloc(perr, lib->seqs[i], j);

		}
		strcat(lib->seqs[i], p);
		k += strlen(p);
	    }
	}
    }
    if(i < 0) {
	jump_append_new_chunk(perr, &lib->error, "Empty ");
	goto ERROR;
    }
    else if(strlen(lib->seqs[i]) == 0) {
	jump_append_new_chunk(perr, &lib->error, "Empty sequence in ");
	goto ERROR;
    }
    m += upcase_and_check_char(lib->seqs[i]);
    if (m > 0) {
	jump_append_new_chunk(perr, &lib->warning,
			    "Unrecognized character in ");
	jump_append(perr, &lib->warning, errfrag);
	jump_append(perr, &lib->warning, " ");
	jump_append(perr, &lib->warning, lib->repeat_file);
    }
    fclose(file);
    reverse_complement_seq_lib(perr, lib);
    return;

 ERROR:
    jump_append(perr, &lib->error, errfrag);
    jump_append(perr, &lib->error, " ");
    jump_append(perr, &lib->error, lib->repeat_file);
    fclose(file);
}

/* 
 * Removes spaces and "end-of-line" characters from the sequence, replaces all
 * other characters except A, T, G and C with N. Returns 0 if there were not
 * such replacements and 1 otherwise.
 */
static short upcase_and_check_char(s)
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

        case 'b': case 'B': 
        case 'd': case 'D':
        case 'h': case 'H':
        case 'v': case 'V':
        case 'r': case 'R':
        case 'y': case 'Y':
        case 'k': case 'K':
        case 'm': case 'M':
	case 's': case 'S':
	case 'w': case 'W':
	  s[i-j] = toupper(s[i]); break;

	case '\n': j++;          break;
	case ' ' : j++;          break;
	case '\t': j++;          break;
	default  : s[i-j] = 'N';
	    m = 1;
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
reverse_complement_seq_lib(perr, lib)
     primer_error *perr;
     seq_lib  *lib;
{
    int i, n, k;
    if((n = lib->seq_num) == 0) return;
    else {
	lib->names = pr_jump_realloc(perr,
				     lib->names,
				     2*n*sizeof(*lib->names));
	lib->seqs = pr_jump_realloc(perr,
				    lib->seqs,
				    2*n*sizeof(*lib->seqs));
	lib->weight = pr_jump_realloc(perr,
				      lib->weight,
				      2*n*sizeof(*lib->weight));
	lib->rev_compl_seqs = pr_jump_malloc(perr,
					     2*n*sizeof(*lib->seqs));

	lib->seq_num *= 2;
	for(i=n; i<lib->seq_num; i++){
	    k = strlen(lib->names[i-n]);
	    lib->names[i] = pr_jump_malloc(perr, k + 9);
	    strcpy(lib->names[i], "reverse ");
	    strcat(lib->names[i], lib->names[i-n]);
	    lib->seqs[i] = pr_jump_malloc(perr, strlen(lib->seqs[i-n]) + 1);
	    reverse_complement(lib->seqs[i-n], lib->seqs[i]);
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

