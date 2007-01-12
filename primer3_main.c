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

#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <float.h>
#include "primer3_release.h"
#include "format_output.h"
#include "dpal.h"
#include "oligotm.h"
#include "primer3.h"
#include "boulder_input.h"

/* #define's */

/* 
 * Panic messages for when the program runs out of memory.  pr_program_name and
 * pr_program_name_len must be set at the beginning of main.
 */
#define OOM_MESSAGE      ": out of memory\n"
#define OOM_MESSAGE_LEN  16
#define OOM_STMT1 write(2, pr_program_name, pr_program_name_len)
#define OOM_STMT2 write(2, OOM_MESSAGE, OOM_MESSAGE_LEN), exit(-2)
#define OOM_ERROR OOM_STMT1, OOM_STMT2

#ifndef MAX_PRIMER_LENGTH
#error "Define MAX_PRIMER_LENGTH in Makefile..."
  /* to ensure that MAX_PRIMER_LENGTH <= DPAL_MAX_ALIGN. */
#endif
#if (MAX_PRIMER_LENGTH > DPAL_MAX_ALIGN) 
#error "MAX_PRIMER_LENGTH must be <= DPAL_MAX_ALIGN"
#endif

#define MAX_NN_TM_LENGTH 36 /* The maxium length for which to use the
			       nearest neighbor model when calculating
			       oligo Tms. */

#define ALIGN_SCORE_UNDEF     SHRT_MIN

#define MACRO_CAT_2(A,B) A##B

#define INITIAL_LIST_LEN     2000 /* Initial size of oligo lists. */
#define INITIAL_NUM_RETURN   5    /* Initial space to allocate for pairs to
				     return. */

#define PAIR_OK 1
#define PAIR_FAILED 0

#define OK_OR_MUST_USE(H) ((H)->ok == OV_OK || (H)->must_use)

/* Function declarations. */
static void   add_must_use_warnings(seq_args *, const char *,
				    const oligo_stats *);
static void   add_pair(const primer_pair *, pair_array_t *);
static short  align(const char *, const char*, const dpal_args *a);
static int    choose_pair(const primer_args *,
			  seq_args *,
			  int, int, int, int,
			  const dpal_args*, const dpal_args*,
			  const dpal_args*,
			  pair_array_t *);
static void   check_sequence_quality(const primer_args *, primer_rec *,
				     oligo_type, const seq_args *, int, int,
				     int *, int *);
static int    choose_internal_oligo(const primer_rec *, const primer_rec *,
				    int, int *, seq_args *,
				    const primer_args *,
				    const dpal_args*,
				    const dpal_args *, const dpal_args *);
void          compute_position_penalty(const primer_args *, const seq_args *, 
				       primer_rec *, oligo_type);

static void   create_and_print_file(const seq_args *, int, const primer_rec[],
				    const oligo_type, const int, const int,
				    const char *);
static int    find_stop_codon(const char *, int, int);
static void   gc_and_n_content(const int, const int, const char *, primer_rec *);
static int    make_internal_oligos_list(const primer_args *,
					seq_args *,
					int *, const dpal_args *,
					const dpal_args*,
					const dpal_args *);
static int    make_lists(const primer_args *,
			 seq_args *,
			 int *, int *, const dpal_args *, 
			 const dpal_args*, const dpal_args *);
static double obj_fn(const primer_args *, primer_pair *);
static short  oligo_max_template_mispriming(const primer_rec *);
static int    oligo_overlaps_interval(const int, const int,
				      interval_array_t, const int);
static int    oligo_pair_seen(const primer_pair *, const pair_array_t *);
static void   oligo_param(const primer_args *pa,
			  primer_rec *, oligo_type,
			  const dpal_args*,
			  const dpal_args *, const dpal_args *,
			  seq_args *, oligo_stats *);
static int    pair_param(const primer_args *, seq_args *,
			 int, int, int, primer_pair *,
			 const dpal_args*,
			 const dpal_args *, const dpal_args*);
static int    pair_spans_target(const primer_pair *, const seq_args *);
static int    primer_pair_comp(const void *, const void*);
static int    primer_rec_comp(const void *, const void *);
static void   print_list(const seq_args *, const primer_args *,
			 int, int, int);
static void   print_list_header(FILE *, oligo_type, int, int);
static void   print_oligo(FILE *, const seq_args *, int, const primer_rec *,
			  oligo_type, int, int);
static void   print_usage();
static void   boulder_print_pairs(const program_args *, const primer_args *, 
				  const seq_args *, const pair_array_t *);
static void   print_explain(const oligo_stats *, oligo_type);
static void   print_all_explain(const primer_args *, const seq_args *);
static FILE   *safe_fopen(const char*, const char*);
static void   set_dpal_args(dpal_args *);
static void   sig_handler(int);
static double p_obj_fn(const primer_args *, primer_rec *, int );
static void   oligo_compl(primer_rec *, const primer_args *, seq_args *,
			  oligo_type, const dpal_args *, 
			  const dpal_args*, const dpal_args *);
static void   oligo_mispriming(primer_rec *,
			       const primer_args *,
			       seq_args *,
			       oligo_type, 
			       const dpal_args *);
static int    pair_repeat_sim(primer_pair *, const primer_args *);
static void   boulder_print_oligos(const primer_args *, 
				   const seq_args *, int, oligo_type);
static void   free_repeat_sim_score(int, int, int);

/* edited by T. Koressaar for lowercase masking:  */
static void   check_if_lowercase_masked(const int position,
					const char *sequence,
					primer_rec *h);

/* Global static variables. */
static const char * copyright[] = {
"",
"Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006",
"Whitehead Institute for Biomedical Research, Steve Rozen",
"(http://jura.wi.mit.edu/rozen), and Helen Skaletsky",
"All rights reserved.",
"",
"Redistribution and use in source and binary forms, with or without",
"modification, are permitted provided that the following conditions are",
"met:",
"",
"   * Redistributions of source code must retain the above copyright",
"notice, this list of conditions and the following disclaimer.",
"   * Redistributions in binary form must reproduce the above",
"copyright notice, this list of conditions and the following disclaimer",
"in the documentation and/or other materials provided with the",
"distribution.",
"   * Neither the names of the copyright holders nor contributors may",
"be used to endorse or promote products derived from this software",
"without specific prior written permission.",
"",
"THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS",
"\"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT",
"LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR",
"A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT",
"OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,",
"SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT",
"LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,",
"DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY",
"THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT",
"(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE",
"OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.",
NULL
};

/* Primer lists and their lengths. */
static primer_rec *f   = NULL;
static primer_rec *r   = NULL;
static primer_rec *mid = NULL;
static int f_len=0, r_len=0, mid_len=0;

/* Other global variables. */
const char *pr_program_name;
int pr_program_name_len;

/* Argument structures for dpal (too cumbersome to pass around
   individually.) */
static dpal_args *lib_local_dpal_args, *lib_local_end_dpal_args;

typedef struct oligo_array {
  int len;
  primer_rec *data;
} oligo_array;

#define FREE_STUFF { \
    free(sa); free(local_args); free(end_args);               \
    free(local_end_args);                                     \
    free(lib_local_dpal_args); free(lib_local_end_dpal_args); \
    if (0 != best_pairs.storage_size) free(best_pairs.pairs); \
    free_seq_lib(&pa->repeat_lib);                            \
    free_seq_lib(&pa->io_mishyb_library); free(pa);           \
}

int
main(argc,argv)
    int argc;
    char *argv[]; 
{ 
    program_args prog_args;
    primer_args *pa;
    seq_args *sa;
    dpal_args *local_args;
    dpal_args *end_args;
    dpal_args *local_end_args;
    int input_found=0;
    pair_array_t best_pairs;
    int n_f, n_r, n_m;  /* Will be initialized in pr_choice. */

    pr_program_name = argv[0];
    pr_program_name_len = strlen(argv[0]);
    /* 
     * We set up some signal handlers in case someone starts up the program
     * from the command line, wonders why nothing is happening, and then kills
     * the program.
     */
    signal(SIGINT, sig_handler);
    signal(SIGTERM, sig_handler);

    /* 
     * We allocate the following structures on the heap rather than on the
     * stack in order to take advantage of testcenter memory access checking
     * during testing.
     */
    pa = pr_safe_malloc(sizeof(*pa));
    sa = pr_safe_malloc(sizeof(*sa));
    local_args = pr_safe_malloc(sizeof(*local_args));
    end_args = pr_safe_malloc(sizeof(*end_args));
    local_end_args = pr_safe_malloc(sizeof(*local_end_args));
    lib_local_end_dpal_args = pr_safe_malloc(sizeof(*lib_local_end_dpal_args));
    lib_local_dpal_args = pr_safe_malloc(sizeof(*lib_local_dpal_args));

    best_pairs.storage_size = best_pairs.num_pairs = 0;
    pr_set_default_global_args(pa);

    memset(&prog_args, 0, sizeof(prog_args));
    while (--argc > 0) {
	argv++;
	if (!strcmp(*argv, "-format_output"))
	    prog_args.format_output = 1;
	else if (!strcmp(*argv, "-2x_compat"))
	    prog_args.twox_compat = 1;
	else if (!strcmp(*argv, "-strict_tags"))
	    prog_args.strict_tags = 1;
	else  {
	    print_usage();
	    FREE_STUFF;
	    exit(-1);
	}
    }

    set_dpal_args(local_args);
    local_args->flag = DPAL_LOCAL;

    set_dpal_args(end_args);
    end_args->flag = DPAL_GLOBAL_END;

    set_dpal_args(local_end_args);
    local_end_args->flag = DPAL_LOCAL_END;

    /* 
     * Read the data from input stream record by record and process it if
     * there are no errors.
     */
    while ((read_record(&prog_args, pa, sa)) > 0) {
	input_found = 1;

	*lib_local_dpal_args = *local_args;
	if (pa->lib_ambiguity_codes_consensus) {
	  PR_ASSERT(dpal_set_ambiguity_code_matrix(lib_local_dpal_args));
	}

	*lib_local_end_dpal_args = *local_end_args;
	if (pa->lib_ambiguity_codes_consensus) {
	  PR_ASSERT(dpal_set_ambiguity_code_matrix(lib_local_end_dpal_args));
	}

	n_f = n_m = n_r = 0;
	if (NULL == sa->error.data && NULL == pa->glob_err.data) {
	    pr_choice(pa, sa, local_args, end_args, local_end_args, 
			 &best_pairs, &n_f, &n_r, &n_m);
        }
	if (NULL != pa->glob_err.data) 
	    pr_append_new_chunk(&sa->error, pa->glob_err.data);

	if (pick_pcr_primers == pa->primer_task
	    || pick_pcr_primers_and_hyb_probe == pa->primer_task) {
	   if (prog_args.format_output) {
	       format_pairs(stdout, pa, sa, &best_pairs);
	   }
	   else {
	       boulder_print_pairs(&prog_args, pa, sa, &best_pairs);
	   }
	}
	else if(pa->primer_task == pick_left_only) {
	   if (prog_args.format_output) 
		     format_oligos(stdout, pa, sa, f, n_f, OT_LEFT);
	   else boulder_print_oligos(pa, sa, n_f, OT_LEFT);
        }
	else if(pa->primer_task == pick_right_only) {
	   if (prog_args.format_output) 
		     format_oligos(stdout, pa, sa, r, n_r, OT_RIGHT);
	   else boulder_print_oligos(pa, sa, n_r, OT_RIGHT);
        }
	else if(pa->primer_task == pick_hyb_probe_only) {
	   if(prog_args.format_output) 
		     format_oligos(stdout, pa, sa, mid, n_m, OT_INTL);
	   else boulder_print_oligos(pa, sa, n_m, OT_INTL);
        }

	best_pairs.num_pairs = 0;

	if(pa->repeat_lib.seq_num > 0 || pa->io_mishyb_library.seq_num > 0) 
	  free_repeat_sim_score(n_f, n_r, n_m);
	if (NULL != sa->internal_input) free(sa->internal_input);
	if (NULL != sa->left_input) free(sa->left_input);
	if (NULL != sa->right_input) free(sa->right_input);
	if (NULL != sa->sequence) free(sa->sequence);
	if (NULL != sa->quality)  free(sa->quality);
	if (NULL != sa->trimmed_seq) free(sa->trimmed_seq);

	/* edited by T. Koressaar for lowercase masking */
	if (NULL != sa->trimmed_orig_seq) free(sa->trimmed_orig_seq);

	if (NULL != sa->upcased_seq) free(sa->upcased_seq);
	if (NULL != sa->upcased_seq_r) free(sa->upcased_seq_r);
	if (NULL != sa->sequence_name) free(sa->sequence_name);
	if (NULL != sa->error.data) free(sa->error.data);
	if (NULL != sa->warning.data) free(sa->warning.data);
	
	if (NULL != pa->glob_err.data) {
	    fprintf(stderr, "%s: %s\n", pr_program_name, pa->glob_err.data);
	    free(pa->glob_err.data);
	    FREE_STUFF;
	    exit(-4);
	}
    }
    FREE_STUFF;   /* Free memory */
    if (0 == input_found) {
	print_usage();
	exit(-3);
    }
    return 0;
}
#undef FREE_STUFF

/* 
 * Find up to pa->num_return primer pairs for the sequence seq with t targets.
 * Set sa->error and return 1 on error; otherwise return 0.
 */
void
pr_choice(pa, sa, local_args, end_args, local_end_args, best_pairs,
	  n_f, n_r, n_m)
  primer_args *pa;
  seq_args *sa;
  const dpal_args *local_args, *end_args, *local_end_args;
  pair_array_t *best_pairs;
  int *n_f;     /* Number of acceptable left primers found. */
  int *n_r;     /* Number of acceptable right primers found. */
  int *n_m;     /* Number of acceptable internal oligos found. */
  /* *n_f, *n_r, *n_m must initialized before calling this function. */
{
    int i;    /* Loop index. */
    int int_num; /* Product size range counter. */
    pair_array_t p;

    PR_ASSERT(NULL != pa);
    PR_ASSERT(NULL != sa);
    PR_ASSERT(NULL != local_args);
    PR_ASSERT(NULL != end_args);
    PR_ASSERT(NULL != local_end_args);
    PR_ASSERT(NULL != n_f);
    PR_ASSERT(NULL != n_r);
    PR_ASSERT(NULL != n_m);
    PR_ASSERT(0 == *n_f);
    PR_ASSERT(0 == *n_r);
    PR_ASSERT(0 == *n_m);
    
    if (_pr_data_control(pa, sa) !=0 ) return;

    if (NULL == f) {
	f = pr_safe_malloc(sizeof(*f) * INITIAL_LIST_LEN);
	r = pr_safe_malloc(sizeof(*r) * INITIAL_LIST_LEN);
	f_len = r_len = INITIAL_LIST_LEN;
    }

    if (make_lists(pa, sa, n_f, n_r,local_args, 
		   end_args, local_end_args)!=0 ) {
	return;
    }

    if ((pa->primer_task == pick_hyb_probe_only 
	 || pa->primer_task == pick_pcr_primers_and_hyb_probe)
        && make_internal_oligos_list(pa, sa, n_m, local_args, 
				     end_args, local_end_args) != 0)
	return;

    /* Creates files with left, right, and internal oligos. */
    if (pa->file_flag) print_list(sa, pa, *n_f, *n_r, *n_m);

    /* We sort _after_ printing lists to maintain the order of test output. */
    if(pa->primer_task != pick_left_only 
       && pa->primer_task != pick_hyb_probe_only) 
      qsort(&r[0], *n_r, sizeof(*r), primer_rec_comp);
    if(pa->primer_task != pick_right_only
       && pa->primer_task != pick_hyb_probe_only) 
      qsort(&f[0], *n_f, sizeof(*f), primer_rec_comp);
    if(pa->primer_task == pick_hyb_probe_only)
      qsort(&mid[0], *n_m, sizeof(*mid), primer_rec_comp);

    p.storage_size = p.num_pairs = 0;
    if (pa->primer_task == pick_pcr_primers 
	|| pa->primer_task == pick_pcr_primers_and_hyb_probe){

      /* Look for pa->num_return best primer pairs. */
      for(int_num=0; int_num < pa->num_intervals; int_num++) {
	if(choose_pair(pa, sa, *n_f, *n_r, *n_m, int_num,
		       local_args, end_args, local_end_args, &p)!=0)
	  continue;

	for (i = 0; i < p.num_pairs && best_pairs->num_pairs < pa->num_return;
	     i++)
	  if (!oligo_pair_seen(&p.pairs[i], best_pairs))
	    add_pair(&p.pairs[i], best_pairs);

	if (pa->num_return == best_pairs->num_pairs) break;
	p.num_pairs = 0;
      }
    }

    /* If it was necessary to use a left_input, right_input,
       or internal_oligo_input primer that was
       unacceptable, then add warnings. */
    if (pa->pick_anyway) {
      if (sa->left_input) {
	add_must_use_warnings(sa, "Left primer", &sa->left_expl);
      }
      if (sa->right_input) {
	add_must_use_warnings(sa, "Right primer", &sa->right_expl);
      }
      if (sa->internal_input) {
	add_must_use_warnings(sa, "Hybridization probe", &sa->intl_expl);
      }
    }

    if (0 != p.storage_size) free(p.pairs);
    return;
}

/* Call this function only if the 'stat's contains
   the _errors_ associated with a given primer
   i.e. that primer was supplied by the caller
   and pick_anyway is set. */
static void
add_must_use_warnings(sa, text, stats)
  seq_args *sa;
  const char* text;
  const oligo_stats *stats;
{
  const char *sep = "/";
  pr_append_str s;

  s.data = NULL;
  s.storage_size = 0;

  if (stats->ns) pr_append_w_sep(&s, sep, "Too many Ns");
  if (stats->target) pr_append_w_sep(&s, sep, "Overlaps Target");
  if (stats->excluded) pr_append_w_sep(&s, sep, "Overlaps Excluded Region");
  if (stats->gc) pr_append_w_sep(&s, sep, "Unacceptable GC content");
  if (stats->gc_clamp) pr_append_w_sep(&s, sep, "No GC clamp");
  if (stats->temp_min) pr_append_w_sep(&s, sep, "Tm too low");
  if (stats->temp_max) pr_append_w_sep(&s, sep, "Tm too high");
  if (stats->compl_any) pr_append_w_sep(&s, sep, "High self complementarity");
  if (stats->compl_end) 
    pr_append_w_sep(&s, sep, "High end self complementarity");
  if (stats->repeat_score)
    pr_append_w_sep(&s, sep, "High similarity to mispriming or mishyb library");
  if (stats->poly_x) pr_append_w_sep(&s, sep, "Long poly-X");
  if (stats->seq_quality) pr_append_w_sep(&s, sep, "Low sequence quality");
  if (stats->stability) pr_append_w_sep(&s, sep, "High 3' stability");
  if (stats->no_orf) pr_append_w_sep(&s, sep, "Would not amplify any ORF");

  /* edited by T. Koressaar for lowercase masking: */
  if (stats->gmasked)
    pr_append_w_sep(&s, sep, "Masked with lowercase letter");

  if (s.data) {
    pr_append_new_chunk(&sa->warning, text);
    pr_append(&sa->warning, " is unacceptable: ");
    pr_append(&sa->warning, s.data);
    free(s.data);
  }

}

/* Return 1 iff pair is already in the first num_pairs elements of retpair. */
static int
oligo_pair_seen(pair, retpair)
    const primer_pair *pair;
    const pair_array_t *retpair;
{
    const primer_pair *q = &retpair->pairs[0], 
		      *stop = &retpair->pairs[retpair->num_pairs];
    for (; q < stop; q++) {
	if (q->left->start == pair->left->start
	    && q->left->length == pair->left->length
	    && q->right->start == pair->right->start
	    && q->right->length == pair->right->length) return 1;
    }
    return 0;
}

/* Add 'pair' to 'retpair'. */
static void
add_pair(pair, retpair)
    const primer_pair *pair;
    pair_array_t *retpair;
{
    if (0 == retpair->storage_size) {
	retpair->storage_size = INITIAL_NUM_RETURN;
	retpair->pairs 
	    = pr_safe_malloc(retpair->storage_size * sizeof(*retpair->pairs));
    } else if (retpair->storage_size == retpair->num_pairs) {
	retpair->storage_size *= 2;
	retpair->pairs
	    = pr_safe_realloc(retpair->pairs,
			      retpair->storage_size * sizeof(*retpair->pairs));
    }
    retpair->pairs[retpair->num_pairs] = *pair;
    retpair->num_pairs++;
}

/* Make lists of acceptable left and right primers. */
static int
make_lists(pa, sa, n_f, n_r, local_args, end_args, local_end_args)
    const primer_args *pa;
    seq_args *sa;
    int *n_f, *n_r;
    const dpal_args *local_args, *end_args, *local_end_args;
{
    int left, right;
    int i,j,n,k,pr_min;
    int tar_l, tar_r, f_b, r_b;

    /* 
     * The position of the intial base of the rightmost stop codon that is
     * to the left of sa->start_codon_pos; valid only if sa->start_codon_pos
     * is "not null".  We will not want to include a stop codon to the right
     * of the the start codon in the amplicon.
     */
    int stop_codon1 = -1; 
       
    char s[MAX_PRIMER_LENGTH+1],s1[MAX_PRIMER_LENGTH+1];
    primer_rec h;

    left = right = 0;
    if (!PR_START_CODON_POS_IS_NULL(sa)) {
      stop_codon1 = find_stop_codon(sa->trimmed_seq, 
				    sa->start_codon_pos, -1);

      sa->stop_codon_pos = find_stop_codon(sa->trimmed_seq, 
					   sa->start_codon_pos,  1);
      sa->stop_codon_pos += sa->incl_s;
    }

    pr_min = INT_MAX;
    for(i=0;i<pa->num_intervals;i++)if(pa->pr_min[i]<pr_min)
       pr_min= pa->pr_min[i];

    PR_ASSERT(INT_MAX > (n=strlen(sa->trimmed_seq)));
    tar_r = 0;
    tar_l = n;
    for(i=0;i<sa->num_targets;i++) {
	if(sa->tar[i][0]>tar_r)tar_r = sa->tar[i][0];
	if(sa->tar[i][0]+sa->tar[i][1]-1<tar_l)tar_l=
	    sa->tar[i][0]+sa->tar[i][1]-1;
    }

    if (_PR_DEFAULT_POSITION_PENALTIES(pa)) {
      if (0 == tar_r) tar_r = n;
      if (tar_l == n) tar_l = 0;
    } else {
      tar_r = n;
      tar_l = 0;
    }

    if (pa->primer_task == pick_left_only)
      f_b = n - 1;
    else if (tar_r - 1 < n - pr_min + pa->primer_max_size - 1 
	&& !(pa->pick_anyway && sa->left_input))
      f_b=tar_r - 1;
    else 
      f_b = n - pr_min + pa->primer_max_size-1;
    k = 0;
    if(pa->primer_task != pick_right_only && pa->primer_task != pick_hyb_probe_only){
    left=n; right=0;
    for (i = f_b; i >= pa->primer_min_size - 1; i--) {
	s[0]='\0';
	for (j = pa->primer_min_size; j <= pa->primer_max_size; j++) {
	    if (i-j > n-pr_min-1 && pick_left_only != pa->primer_task) continue;
	    if (i-j+1>=0) {
		if (k >= f_len) {
		    f_len += (f_len >> 1);
		    f = pr_safe_realloc(f, f_len * sizeof(*f));
		}
		h.start=i-j+1;
		h.length=j;
		h.repeat_sim.score = NULL;
		_pr_substr(sa->trimmed_seq,h.start,h.length,s);

		if (sa->left_input && strcmp_nocase(sa->left_input, s))
		  continue;

		h.must_use = (sa->left_input && pa->pick_anyway);

		if (pa->explain_flag) sa->left_expl.considered++;

		if (!PR_START_CODON_POS_IS_NULL(sa)
		/* Make sure the primer would amplify at least part of
		   the ORF. */
		    && (0 != (h.start - sa->start_codon_pos) % 3
			|| h.start <= stop_codon1
			|| (sa->stop_codon_pos != -1 
			    && h.start >= sa->stop_codon_pos))) {
		  if (pa->explain_flag) sa->left_expl.no_orf++;
		  if (!pa->pick_anyway) continue;
		}

		h.repeat_sim.score = NULL;
		oligo_param(pa, &h, OT_LEFT, local_args, end_args, 
			    local_end_args, sa, &sa->left_expl);
		if (OK_OR_MUST_USE(&h)) {
		  h.quality = p_obj_fn(pa, &h, 0);
		  f[k] = h;
		  if(f[k].start < left) left=f[k].start;
		  k++;
		} else if (h.ok==OV_TOO_MANY_NS || h.ok==OV_INTERSECT_TARGET
			   || h.ok==OV_SELF_ANY || h.ok==OV_END_STAB
			   || h.ok==OV_POLY_X || h.ok==OV_EXCL_REGION || h.ok==OV_GC_CLAMP
			   || h.ok == OV_SEQ_QUALITY || h.ok == OV_LIB_SIM ) {
		  /* Break from the inner for loop, because there is no
		     legal longer oligo with the same 3' sequence. */
		  break;
		}
	    }
	    else break;
	}
    }
    }
    *n_f = k;


    if (pa->primer_task == pick_right_only)
      r_b = 0;
    else if (tar_l+1>pr_min - pa->primer_max_size
	&& !(pa->pick_anyway && sa->right_input))
      r_b = tar_l+1;
    else 
      r_b = pr_min - pa->primer_max_size;
    k = 0;
    if(pa->primer_task != pick_left_only && pa->primer_task != pick_hyb_probe_only) {
    for(i=r_b; i<=n-pa->primer_min_size; i++) {
	s[0]='\0';
	for(j = pa->primer_min_size; j <= pa->primer_max_size; j++) {
	    if (i+j<pr_min && pa->primer_task != pick_right_only) continue;
	    if(i+j-1<n) {
		if (k >= r_len) {
		    r_len += (r_len >> 1);
		    r = pr_safe_realloc(r, r_len * sizeof(*r));
		}
		h.start=i+j-1;
		h.length=j;
		h.repeat_sim.score = NULL;
		_pr_substr(sa->trimmed_seq,i,j,s);
		_pr_reverse_complement(s,s1);

		if (sa->right_input && strcmp_nocase(sa->right_input, s1))
		  continue;
		h.must_use = (sa->right_input && pa->pick_anyway);

		h.repeat_sim.score = NULL;
		oligo_param(pa, &h, OT_RIGHT, local_args, end_args, 
					  local_end_args, sa, &sa->right_expl);
		sa->right_expl.considered++;
		if (OK_OR_MUST_USE(&h)) {
		  h.quality = p_obj_fn(pa, &h, 0);
		  r[k] = h;
		  if (r[k].start > right) right=r[k].start;
		  k++;
		} else if (h.ok==OV_TOO_MANY_NS || h.ok==OV_INTERSECT_TARGET
		    || h.ok==OV_SELF_ANY || h.ok==OV_END_STAB
		    || h.ok==OV_POLY_X || h.ok==OV_EXCL_REGION || h.ok==OV_GC_CLAMP
		    || h.ok == OV_SEQ_QUALITY || h.ok == OV_LIB_SIM ) {
		  /* Break from the inner for loop, because there is no
		     legal longer oligo with the same 3' sequence. */
		  break;
		}
	    }
	    else break;
	}
    }
    }
    *n_r=k;

    /* 
     * Return 1 if one of lists is empty or if leftmost left primer and
     * rightmost right primer do not provide sufficient product size.
     */
    sa->left_expl.ok = *n_f;
    sa->right_expl.ok = *n_r;
    if ((pa->primer_task != pick_right_only &&
		  pa->primer_task != pick_hyb_probe_only && 0 == *n_f) || 
       ((pa->primer_task != pick_left_only && 
		  pa->primer_task != pick_hyb_probe_only) && 0 == *n_r)) return 1;
    else if ((pa->primer_task == pick_pcr_primers || 
		  pa->primer_task == pick_pcr_primers_and_hyb_probe) 
		  && right - left < pr_min - 1) {
	sa->pair_expl.product    = 1;
	sa->pair_expl.considered = 1;
	return 1;
    } else return 0;
}

/* 
 * Make complete list of acceptable internal oligos in mid.  Place the number
 * of valid elements in mid in *n_m.  Return 1 if there are no acceptable
 * internal oligos; otherwise return 0.
 */
static int
make_internal_oligos_list(pa, sa, n_m, local_args, end_args, local_end_args)
    const primer_args *pa;
    seq_args *sa;
    int *n_m;
    const dpal_args *local_args, *end_args, *local_end_args;
{
  int i, j, n, k;

  char s[MAX_PRIMER_LENGTH+1];
  primer_rec h;

  if (NULL == mid) {
    mid_len = INITIAL_LIST_LEN;
    mid = pr_safe_malloc(sizeof(*mid) * mid_len);
  }

  n = strlen(sa->trimmed_seq);
  k = 0;
  for(i = n - 1; i >= pa->io_primer_min_size-1; i--) {
    s[0] = '\0';
    for(j = pa->io_primer_min_size; j <=pa->io_primer_max_size; j++) {
      if(i-j < -1) break;
      if (k >= mid_len) {
	mid_len += (mid_len >> 1);
	mid = pr_safe_realloc(mid, mid_len * sizeof(*mid));
      }
      h.start = i - j +1;
      h.length = j;
      h.repeat_sim.score = NULL;
      _pr_substr(sa->trimmed_seq, h.start, h.length,s);

      if (sa->internal_input && strcmp_nocase(sa->internal_input, s))
	continue;
      h.must_use = (sa->internal_input && pa->pick_anyway);

      h.repeat_sim.score = NULL;
      oligo_param(pa, &h, OT_INTL, local_args, end_args, 
		  local_end_args, sa, &sa->intl_expl);
      sa->intl_expl.considered++;
      if (OK_OR_MUST_USE(&h)) {
	h.quality = p_obj_fn(pa, &h, 2);
	mid[k] = h;
	k++;
      } else if (h.ok==OV_TOO_MANY_NS || h.ok==OV_INTERSECT_TARGET
		 || h.ok==OV_SELF_ANY || h.ok==OV_POLY_X
		 || h.ok==OV_EXCL_REGION || h.ok==OV_GC_CLAMP
		 || h.ok==OV_SEQ_QUALITY || h.ok==OV_LIB_SIM ) {
	/* Break from the inner for loop, because there is no
           legal longer oligo with the same 3' sequence. */
	break;
      }
    }
  }
  *n_m = k;
  sa->intl_expl.ok = *n_m;
  if(*n_m==0)return 1;
  else return 0;
}

/*
 * Compute various characteristics of the oligo, and determine
 * if it is acceptable.
 */
#define OUTSIDE_START_WT  30.0
#define INSIDE_START_WT   20.0
#define INSIDE_STOP_WT   100.0
#define OUTSIDE_STOP_WT    0.5
static void
oligo_param(pa, h, l, local_args, end_args, local_end_args, sa, stats)
    const primer_args *pa;
    primer_rec *h;
    oligo_type l;
    const dpal_args *local_args, *end_args, *local_end_args;
    seq_args *sa;
    oligo_stats *stats;
{
    int i,j,k, min_q, min_end_q;
    int poly_x, max_poly_x;
    int must_use = h->must_use;
    const char *seq = sa->trimmed_seq;
    char s1[MAX_PRIMER_LENGTH+1], s1_rev[MAX_PRIMER_LENGTH+1];



    h->ok = OV_UNINITIALIZED;
    h->target = h->gc_content = h->num_ns=h->excl=0;

    h->template_mispriming = h->template_mispriming_r = ALIGN_SCORE_UNDEF;

    PR_ASSERT(OT_LEFT == l || OT_RIGHT == l || OT_INTL == l);
    
    if (OT_LEFT == l || OT_INTL == l) {j = h->start; k=j+h->length-1;}
    else {j = h->start-h->length+1; k=h->start;}

    PR_ASSERT(k >= 0);
    PR_ASSERT(k < TRIMMED_SEQ_LEN(sa));
   
    /* edited by T. Koressaar for lowercase masking */
    if(pa->lowercase_masking==1) {
      if(l==OT_LEFT) {
	 check_if_lowercase_masked(k, sa->trimmed_orig_seq,h);
      }
      if(l==OT_RIGHT) {
	 check_if_lowercase_masked(j, sa->trimmed_orig_seq,h);
      }
      if(h->ok==OV_GMASKED) {
	 stats->gmasked++;
	 if (!must_use) return;
      }
    }
    /* end T. Koressar's changes */

    gc_and_n_content(j, k-j+1, sa->trimmed_seq, h);

    if (((OT_LEFT == l || OT_RIGHT == l) && 
	 h->num_ns > pa->num_ns_accepted) || 
	(OT_INTL == l && h->num_ns > pa->io_num_ns_accepted) ) {
      h->ok = OV_TOO_MANY_NS;
      stats->ns++;
      if (!must_use) return;
    }

    /* Upstream error checking has ensured that we use non-default position
       penalties only when there is 0 or 1 target. */
    PR_ASSERT(sa->num_targets <= 1 || _PR_DEFAULT_POSITION_PENALTIES(pa));
    if (l < 2 
	&& _PR_DEFAULT_POSITION_PENALTIES(pa)
	&& oligo_overlaps_interval(j, k-j+1, sa->tar, sa->num_targets)) {
      h->position_penalty = 0.0;
      h->position_penalty_infinite = '\1';
      h->target = 1;
    } else if (l < 2 && !_PR_DEFAULT_POSITION_PENALTIES(pa)
	     && 1 == sa->num_targets) {
      compute_position_penalty(pa, sa, h, l);
      if (h->position_penalty_infinite) h->target = 1;
    } else {
      h->position_penalty = 0.0;
      h->position_penalty_infinite = '\0';
    }

    if (!PR_START_CODON_POS_IS_NULL(sa)) {
      if (OT_LEFT == l) {
	if (sa->start_codon_pos > h->start)
	  h->position_penalty
	    = (sa->start_codon_pos - h->start) * OUTSIDE_START_WT;
	else
	  h->position_penalty 
	    = (h->start - sa->start_codon_pos) * INSIDE_START_WT;
      }
      else if (OT_RIGHT == l) {
	if (-1 == sa->stop_codon_pos) {
	  h->position_penalty = (TRIMMED_SEQ_LEN(sa) - h->start - 1) * INSIDE_STOP_WT;
	} else if (sa->stop_codon_pos < h->start)
	  h->position_penalty
	    = (h->start - sa->stop_codon_pos) * OUTSIDE_STOP_WT;
	else
	  h->position_penalty
	    = (sa->stop_codon_pos - h->start) * INSIDE_STOP_WT;
      }
    }

    if (l < 2 && oligo_overlaps_interval(j, k-j+1, sa->excl, sa->num_excl))
	h->excl = 1;

    if (l == 2 && oligo_overlaps_interval(j, k-j+1, sa->excl_internal,
					  sa->num_internal_excl))
	h->excl = 1;

    if(l < 2 && h->target==1) {
	h->ok = OV_INTERSECT_TARGET;
	stats->target++;
	if (!must_use) return;
    }

    if(h->excl==1){
	h->ok = OV_EXCL_REGION;
	stats->excluded++;
	if (!must_use) return;
    }

    if((l<2 && (h->gc_content< pa->min_gc || h->gc_content > pa->max_gc)) ||
       (l==2 && (h->gc_content< pa->io_min_gc || h->gc_content > pa->io_max_gc))){
	h->ok = OV_GC_CONTENT;
	stats->gc++;
	if (!must_use) return;
    }
    if(pa->gc_clamp != 0){
       if(OT_LEFT == l){
	 for(i=k-pa->gc_clamp+1; i<= k; i++)if(seq[i] !='G'&&seq[i] !='C'){
	   h->ok = OV_GC_CLAMP;
	   stats->gc_clamp++;
	   if (!must_use) return; else break;
	 }
       }
       if(OT_RIGHT == l){
	 for(i=j; i<j+pa->gc_clamp; i++)if(seq[i] != 'G' && seq[i] != 'C'){
	   h->ok = OV_GC_CLAMP;
	   stats->gc_clamp++;
	   if (!must_use) return; else break;
	 }
       }
    }
            
    check_sequence_quality(pa, h, l, sa, j, k, &min_q, &min_end_q);
    if(OT_LEFT == l || OT_RIGHT == l) {
      if (min_q < pa->min_quality) {
	h->ok = OV_SEQ_QUALITY;
	stats->seq_quality++;
	if (!must_use) return;
      } else if (min_end_q < pa->min_end_quality) {
	h->ok = OV_SEQ_QUALITY;
	stats->seq_quality++;
	if (!must_use) return;
      }
    } else if (OT_INTL == l) {
      if (min_q < pa->io_min_quality) {
	h->ok = OV_SEQ_QUALITY;
	stats->seq_quality++;
	if (!must_use) return;
      }
    } else {
      PR_ASSERT(0); /* Programming error. */
    }

    if(OT_LEFT == l || OT_RIGHT == l) max_poly_x = pa->max_poly_x;     
    else max_poly_x = pa->io_max_poly_x;
    if(max_poly_x > 0) {
          poly_x = 1;
	  for(i=j+1;i<=k;i++){
                if(seq[i] == seq[i-1]||seq[i] == 'N'){
		       poly_x++;
		       if(poly_x > max_poly_x){
			     h->ok = OV_POLY_X;
			     stats->poly_x++;
			     if (!must_use) return; else break;
                       }
                }
		else poly_x = 1;
          }
     }

    _pr_substr(seq,j,k-j+1,s1);
                   
   if(OT_LEFT == l || OT_RIGHT == l)
     h->temp 
     = seqtm(s1, pa->dna_conc, pa->salt_conc, pa->divalent_conc, pa->dntp_conc, 
	     MAX_NN_TM_LENGTH,
	     pa->tm_santalucia,
	     pa->salt_corrections); 
   else
     h->temp
     = seqtm(s1, pa->io_dna_conc, pa->io_salt_conc, pa->io_divalent_conc, pa->io_dntp_conc,
	     MAX_NN_TM_LENGTH,
	     pa->tm_santalucia,
	     pa->salt_corrections);
         
    if (((l == OT_LEFT || l == OT_RIGHT) && h->temp < pa->min_tm)
	|| (l==OT_INTL && h->temp<pa->io_min_tm)) {
	h->ok = OV_TM_LOW;
	stats->temp_min++;
	if (!must_use) return;
    }
    if (((OT_LEFT == l || OT_RIGHT == l) && h->temp>pa->max_tm) 
	|| (OT_INTL == l && h->temp>pa->io_max_tm)) {
	h->ok = OV_TM_HIGH;
	stats->temp_max++;
	if (!must_use) return;
    }
    if (OT_LEFT == l) {
      if ((h->end_stability = end_oligodg(s1, 5,
					  pa->tm_santalucia))
	  > pa->max_end_stability) {
	h->ok = OV_END_STAB;
	stats->stability++;
	if (!must_use) return;
      }
    } else if (OT_RIGHT == l) {
      _pr_reverse_complement(s1, s1_rev);
      if ((h->end_stability = end_oligodg(s1_rev, 5,
					  pa->tm_santalucia))
	  > pa->max_end_stability) {
	  h->ok = OV_END_STAB;
	  stats->stability++;
	  if (!must_use) return;
      }
    }

    if (must_use
	|| pa->file_flag 
	|| (pa->primer_task != pick_pcr_primers && 
	    pa->primer_task != pick_pcr_primers_and_hyb_probe)
	|| ((OT_RIGHT == l || OT_LEFT == l) &&
	    (pa->primer_weights.compl_any || pa->primer_weights.compl_end))
	|| (OT_INTL == l 
	    && (pa->io_weights.compl_any || pa->io_weights.compl_end))) {

      oligo_compl(h, pa, sa, l, local_args, end_args, local_end_args);

      if (h->ok != OV_UNINITIALIZED && !must_use) {
	PR_ASSERT(h->ok != OV_OK);
	return; 
      }

    } else
      h->self_any = h->self_end  = ALIGN_SCORE_UNDEF;

    if (must_use
	|| pa->file_flag
	||(pa->primer_task != pick_pcr_primers && 
	   pa->primer_task != pick_pcr_primers_and_hyb_probe)
	|| ((OT_RIGHT == l || OT_LEFT == l) && pa->primer_weights.repeat_sim)
	|| ((OT_RIGHT == l || OT_LEFT == l) 
	    && pa->primer_weights.template_mispriming)
	|| (OT_INTL == l && pa->io_weights.repeat_sim)) {

      oligo_mispriming(h, pa, sa, l, local_end_args);

    }

    if (OV_UNINITIALIZED == h->ok) h->ok = OV_OK;
}
#undef OUTSIDE_START_WT
#undef INSIDE_START_WT
#undef INSIDE_STOP_WT
#undef OUTSIDE_STOP_WT

static void
check_sequence_quality(pa, h, l, sa, j, k, r_min_q, r_min_q_end)
    const primer_args *pa;
    primer_rec *h;
    oligo_type l;
    const seq_args *sa;
    int j, k;
    int *r_min_q, *r_min_q_end;
{
  int i,min_q, min_q_end, m, q;

  q = min_q = min_q_end = pa->quality_range_max;

  if (NULL != sa->quality) {

    if(OT_LEFT == l || OT_RIGHT == l){
      min_q = pa->min_quality;
      min_q_end = pa->min_end_quality;
    }
    else {min_q = min_q_end = pa->io_min_quality;} 

    if (OT_LEFT == l || OT_INTL == l) {

      for(i = k-4; i <= k; i++) {
	if(i < j) continue;
	m = sa->quality[i + sa->incl_s];
	if (m < q) q = m;
      }
      min_q_end = q;

      for(i = j; i<=k-5; i++) {
	m = sa->quality[i + sa->incl_s];
	if (m < q) q = m;
      }
      min_q = q;

    } else if (OT_RIGHT == l) {
      for(i = j; i < j+5; i++) {
	if(i > k) break;
	m = sa->quality[i + sa->incl_s];
	if (m < q) q = m;
      }
      min_q_end = q;
       
      for(i = j+5; i <= k; i++) {
	m = sa->quality[i + sa->incl_s];
	if (m < q) q = m;
      }
      min_q = q;
    } else {
      PR_ASSERT(0); /* Programming error. */
    }
  }
  h->seq_quality = *r_min_q = min_q;
  *r_min_q_end = min_q_end;
}

/* 
 * Set h->gc_content to the GC % of the 'len' bases starting at 'start' in
 * 'sequence' (excluding 'N's).  Set h->num_ns to the number of 'N's in that
 * subsequence.
 */
static void
gc_and_n_content(start, len, sequence, h)
    const int start, len;
    const char *sequence;
    primer_rec *h;
{
    const char* p = &sequence[start];
    const char* stop = p + len;
    int num_gc = 0, num_gcat = 0, num_n = 0;
    while (p < stop) {
	if ('N' == *p) 
	    num_n++;
	else {
	    num_gcat++;
	    if ('C' == *p || 'G' == *p) num_gc++;
	}
	p++;
    }
    h->num_ns = num_n;
    if (0 == num_gcat) h->gc_content= 0.0;
    else h->gc_content = 100.0 * ((double)num_gc)/num_gcat;
}

static int
oligo_overlaps_interval(start, len, intervals, num_intervals)
    const int start, len;
    interval_array_t intervals;
    const int num_intervals;
{
    int i;
    int last = start + len - 1;
    for (i = 0; i < num_intervals; i++)
	if (!(last < intervals[i][0] 
	      || start > (intervals[i][0] + intervals[i][1] - 1)))
	    return 1;
    return 0;
}

/* Calculate the part of the objective function due to one primer. */
static double
p_obj_fn(pa, h, j)
    const primer_args *pa;
    primer_rec *h;
    int  j;
{
  double sum;

  sum = 0;
  if (j == OT_LEFT || j == OT_RIGHT) {
      if(pa->primer_weights.temp_gt && h->temp > pa->opt_tm)
	   sum += pa->primer_weights.temp_gt * (h->temp - pa->opt_tm);
      if(pa->primer_weights.temp_lt && h->temp < pa->opt_tm)
	   sum += pa->primer_weights.temp_lt * (pa->opt_tm - h->temp);

      if(pa->primer_weights.gc_content_gt && h->gc_content > pa->opt_gc_content)
	   sum += pa->primer_weights.gc_content_gt 
	     * (h->gc_content - pa->opt_gc_content);
      if(pa->primer_weights.gc_content_lt && h->gc_content < pa->opt_gc_content)
	   sum += pa->primer_weights.gc_content_lt 
	     * (pa->opt_gc_content - h->gc_content);

      if(pa->primer_weights.length_lt && h->length < pa->primer_opt_size)
	   sum += pa->primer_weights.length_lt * (pa->primer_opt_size - h->length);
      if(pa->primer_weights.length_gt && h->length > pa->primer_opt_size)
	   sum += pa->primer_weights.length_gt * (h->length - pa->primer_opt_size);
      if(pa->primer_weights.compl_any)
	   sum += pa->primer_weights.compl_any * h->self_any
	     / PR_ALIGN_SCORE_PRECISION;
      if(pa->primer_weights.compl_end)
	   sum += pa->primer_weights.compl_end * h->self_end
	     / PR_ALIGN_SCORE_PRECISION;
      if(pa->primer_weights.num_ns)
	   sum += pa->primer_weights.num_ns * h->num_ns;
      if(pa->primer_weights.repeat_sim)
	   sum += pa->primer_weights.repeat_sim 
	     * h->repeat_sim.score[h->repeat_sim.max]
	     / PR_ALIGN_SCORE_PRECISION;
      if (!h->target) {
	/* We might be evaluating p_obj_fn with h->target if
	   the client supplied 'pick_anyway' and specified
	   a primer or oligo. */
	PR_ASSERT(!h->position_penalty_infinite);
	if(pa->primer_weights.pos_penalty)
	  sum += pa->primer_weights.pos_penalty * h->position_penalty;
      }
      if(pa->primer_weights.end_stability)
	   sum += pa->primer_weights.end_stability * h->end_stability;
      if(pa->primer_weights.seq_quality)
	   sum += pa->primer_weights.seq_quality * 
			     (pa->quality_range_max - h->seq_quality);

      if (pa->primer_weights.template_mispriming) {
	PR_ASSERT(oligo_max_template_mispriming(h) != ALIGN_SCORE_UNDEF);
	sum += pa->primer_weights.template_mispriming * 
	  oligo_max_template_mispriming(h);
      }

      return sum;
  } else if (j == OT_INTL) {
      if(pa->io_weights.temp_gt && h->temp > pa->io_opt_tm)
	 sum += pa->io_weights.temp_gt * (h->temp - pa->io_opt_tm);
      if(pa->io_weights.temp_lt && h->temp < pa->io_opt_tm)
	 sum += pa->io_weights.temp_lt * (pa->io_opt_tm - h->temp);

      if(pa->io_weights.gc_content_gt && h->gc_content > pa->io_opt_gc_content)
	   sum += pa->io_weights.gc_content_gt
	     * (h->gc_content - pa->io_opt_gc_content);
      if(pa->io_weights.gc_content_lt && h->gc_content < pa->io_opt_gc_content)
	   sum += pa->io_weights.gc_content_lt
	     * (pa->io_opt_gc_content - h->gc_content);

      if(pa->io_weights.length_lt && h->length < pa->io_primer_opt_size)
	 sum += pa->io_weights.length_lt * (pa->io_primer_opt_size - h->length);
      if(pa->io_weights.length_gt && h->length  > pa->io_primer_opt_size)
	 sum += pa->io_weights.length_gt * (h->length - pa->io_primer_opt_size);
      if(pa->io_weights.compl_any)
	 sum += pa->io_weights.compl_any * h->self_any / PR_ALIGN_SCORE_PRECISION;
      if(pa->io_weights.compl_end)
	 sum += pa->io_weights.compl_end * h->self_end / PR_ALIGN_SCORE_PRECISION;
      if(pa->io_weights.num_ns)
	 sum += pa->io_weights.num_ns * h->num_ns;
      if(pa->io_weights.repeat_sim)
	 sum += pa->io_weights.repeat_sim
	   * h->repeat_sim.score[h->repeat_sim.max]
	   / PR_ALIGN_SCORE_PRECISION;
      if(pa->io_weights.seq_quality)
	 sum += pa->io_weights.seq_quality * 
		     (pa->quality_range_max - h->seq_quality);
      return sum;
  } else {
    PR_ASSERT(0); /* Programmig error. */
  }
}

/* Return max of h->template_mispriming and h->template_mispriming_r (max
   template mispriming on either strand). */
static short
oligo_max_template_mispriming(h)
  const primer_rec *h;
{
  return h->template_mispriming > h->template_mispriming_r ?
    h->template_mispriming : h->template_mispriming_r;
}

static void
print_list(sa, pa, n_f, n_r, n_m)
    const seq_args *sa;
    const primer_args *pa;
    int n_f, n_r, n_m;
{
    int first_base_index = pa->first_base_index;

    if(pa->primer_task != pick_right_only && pa->primer_task != pick_hyb_probe_only)
       create_and_print_file(sa, n_f, f, OT_LEFT, first_base_index, 
			  NULL != pa->repeat_lib.repeat_file, ".for");
    if(pa->primer_task != pick_left_only && pa->primer_task != pick_hyb_probe_only)
       create_and_print_file(sa, n_r, r, OT_RIGHT, first_base_index,
			  NULL != pa->repeat_lib.repeat_file, ".rev");
    if ( pa->primer_task == pick_pcr_primers_and_hyb_probe 
				|| pa->primer_task == pick_hyb_probe_only)
      create_and_print_file(sa, n_m, mid, OT_INTL,
			    first_base_index,
			    NULL != pa->io_mishyb_library.repeat_file,
			    ".int");
}

static void
create_and_print_file(sa, n, oligo_arr, o_type,
		      first_base_index, print_lib_sim, ext)
    const seq_args *sa;
    int n;
    const primer_rec oligo_arr[];
    const oligo_type o_type;
    const int first_base_index, print_lib_sim;
    const char *ext;
{
    int i;
    char *file = pr_safe_malloc(strlen(sa->sequence_name) + strlen(ext) + 1);
    FILE *fh;

    strcpy(file, sa->sequence_name);
    strcat(file,ext);
    fh = safe_fopen(file,"w");
    print_list_header(fh, o_type, first_base_index, print_lib_sim);
    for(i=0; i<n; i++) 
	print_oligo(fh, sa, i, &oligo_arr[i], o_type,
		    first_base_index, print_lib_sim);
    fclose(fh);
    free(file);
}

static void
print_list_header(fh, type, first_base_index, print_lib_sim)
    FILE *fh;
    oligo_type type;
    int first_base_index, print_lib_sim;
{
    fprintf(fh, "ACCEPTABLE %s\n",
            OT_LEFT == type ? "LEFT PRIMERS"
            : OT_RIGHT == type ? "RIGHT PRIMERS" : "INTERNAL OLIGOS");
    fprintf(fh, "                               %4d-based     ", 
	    first_base_index);
    if (print_lib_sim)
	fprintf(fh, "#               self  self   lib  qual-\n");
    else
	fprintf(fh, "#               self  self  qual-\n");
    fprintf(fh, "   # sequence                       start ln  "); 
    if (print_lib_sim)
	fprintf(fh, "N   GC%%     Tm   any   end   sim   lity\n");
    else 
	fprintf(fh, "N   GC%%     Tm   any   end   lity\n");
}

static void
print_oligo(fh, sa, index, h, type, first_base_index, print_lib_sim)
    FILE *fh;
    const seq_args *sa;
    int index;
    const primer_rec *h;
    oligo_type type;
    int first_base_index, print_lib_sim;
{
    char *p = (OT_RIGHT != type) 
	? pr_oligo_sequence(sa, h) : pr_oligo_rev_c_sequence(sa, h);
    if (print_lib_sim)
	fprintf(fh,
		"%4d %-30s %5d %2d %2d %5.2f %5.3f %5.2f %5.2f %5.2f %6.3f\n",
		index, p, h->start+sa->incl_s + first_base_index, 
		h->length,
		h->num_ns, h->gc_content, h->temp,
		h->self_any / PR_ALIGN_SCORE_PRECISION,
		h->self_end / PR_ALIGN_SCORE_PRECISION,
		h->repeat_sim.score[h->repeat_sim.max] /PR_ALIGN_SCORE_PRECISION,
		h->quality);
    else
	fprintf(fh,
		"%4d %-30s %5d %2d %2d %5.2f %5.3f %5.2f %5.2f %6.3f\n",
		index, p, h->start+sa->incl_s + first_base_index, 
		h->length,
		h->num_ns, h->gc_content, h->temp,
		h->self_any / PR_ALIGN_SCORE_PRECISION,
		h->self_end / PR_ALIGN_SCORE_PRECISION,
		h->quality);
}

static int
choose_pair(pa, sa, n_f, n_r, n_m, int_num, local_args, end_args, 
					      local_end_args, p)
    const primer_args *pa;
    seq_args *sa;
    int n_f;
    int n_r;
    int n_m;
    int int_num;
    const dpal_args *local_args, *end_args, *local_end_args;
    pair_array_t *p;
{
  int i,j;
  int k; /* 
	  * The number of acceptable primer pairs saved in p.
	  * (k <= pa->num_return.)
	  */
  int i0, n_last, n_int;
  primer_pair worst_pair; /* The worst pair among those being "remembered". */
  int i_worst;             /* The index within p of worst_pair. */

  primer_pair h;

  k=0; 

  i_worst = 0;
  n_last = n_f;
  for(i=0;i<n_r;i++) {
    /* 
     * Make a quick cut based on the the quality of the best left
     * primer.
     *
     * worst_pair must be defined if k >= pa->num_return.
     */
    if (!OK_OR_MUST_USE(&r[i])) continue;
    if (k >= pa->num_return 
	&& (r[i].quality+f[0].quality > worst_pair.pair_quality 
	    || worst_pair.pair_quality == 0))
      break;

    for(j=0;j<n_last;j++) {
      /* 
       * Invariant: if 2 pairs in p have the same pair_quality, then the
       * pair discovered second will have a higher index within p.
       *
       * This invariant is needed to produce consistent results when 2
       * pairs have the same pair_quality.
       */

      if (!OK_OR_MUST_USE(&r[i])) break;
      if (!OK_OR_MUST_USE(&f[j])) continue;
      if(k>= pa->num_return 
	 && (f[j].quality+r[i].quality > worst_pair.pair_quality
	     || worst_pair.pair_quality == 0)) {
	/* worst_pair must be defined if k >= pa->num_return. */
	n_last=j;
	break;
      }

      if (PAIR_OK ==
	  pair_param(pa, sa, j, i, int_num, &h, local_args, end_args,
		     local_end_args)) {

	if (!pa->pr_pair_weights.io_quality) {
	  h.pair_quality = obj_fn(pa, &h);
	  PR_ASSERT(h.pair_quality >= 0.0);
	}


	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe
	     && choose_internal_oligo(h.left, h.right,
				      n_m, &n_int, sa, pa, local_args,
				      end_args, local_end_args)!=0) {
	  sa->pair_expl.internal++;
	  continue;
	}
	sa->pair_expl.ok++;
	if (k < pa->num_return) {
	  if ( pa->primer_task == pick_pcr_primers_and_hyb_probe) 
	    h.intl = &mid[n_int];
	  if(pa->pr_pair_weights.io_quality) {
	    h.pair_quality = obj_fn(pa, &h);
	    PR_ASSERT(h.pair_quality >= 0.0);
	  }

	  add_pair(&h, p);
	  /* p->pairs[k] = h; */
	  if (k == 0 || primer_pair_comp(&h, &worst_pair) > 0){
	    worst_pair = h;
	    i_worst = k;
	  }
	  k++;
	} else {
	  if ( pa->primer_task == pick_pcr_primers_and_hyb_probe)
	    h.intl = &mid[n_int];
	  if(pa->pr_pair_weights.io_quality) {
	    h.pair_quality = obj_fn(pa, &h);
	    PR_ASSERT(h.pair_quality >= 0.0);
	  }

	  if (primer_pair_comp(&h, &worst_pair) < 0) {
	    /* 
	     * There are already pa->num_return results, and vl is better than
	     * the pa->num_return the quality found so far.
	     */
	    p->pairs[i_worst] = h;
	    worst_pair = h; /* h is a lower bound on the worst pair. */
	    for (i0 = 0; i0<pa->num_return; i0++)
	      if(primer_pair_comp(&p->pairs[i0], &worst_pair) > 0) {
		i_worst = i0;
		worst_pair = p->pairs[i0];

	      }
	  }
	}
      }
    }	
  }
  if(k!=0) qsort(p->pairs, k, sizeof(primer_pair), primer_pair_comp);
  p->num_pairs = k;
  if (k==0) return 1;
  else return 0;
}

/* Choose best internal oligo for given pair of left and right primers. */
static int
choose_internal_oligo(left, right, n_m, nm, sa, pa, local_args, 
					end_args,  local_end_args)
     const primer_rec *left, *right;
     int n_m;
     int *nm;
     seq_args *sa;
     const primer_args *pa;
     const dpal_args *local_args, *end_args, *local_end_args;
{
   int i,k;
   double min;

   min = 1000000.;
   i = -1;
   for(k=0; k < n_m; k++){
     if ((mid[k].start > (left->start + (left->length-1))) 
	&& ((mid[k].start+(mid[k].length-1)) < (right->start-right->length+1)) 
	&& (mid[k].quality < min) 
	&& (OK_OR_MUST_USE(&mid[k]))) {
       if(mid[k].self_any == ALIGN_SCORE_UNDEF){
	 oligo_compl(&mid[k], pa, sa, OT_INTL, local_args,
		     end_args, local_end_args);
	 if (!OK_OR_MUST_USE(&mid[k])) continue;
       }
       if(mid[k].repeat_sim.score == NULL) {
         oligo_mispriming(&mid[k], pa, sa, OT_INTL, local_args);
	 if (!OK_OR_MUST_USE(&mid[k])) continue;
       }
       min = mid[k].quality;
       i=k;
     }    
   }
   *nm = i;
   if(*nm < 0) return 1;
   return 0;
}

/* Compare function for sorting primer records. */
static int
primer_rec_comp(x1, x2)
    const void *x1, *x2;
{
    const primer_rec *a1 = x1, *a2 = x2;

    if(a1->quality < a2->quality) return -1;
    if (a1->quality > a2->quality) return 1;

    /* 
     * We want primer_rec_comp to always return a non-0 result, because
     * different implementations of qsort differ in how they treat "equal"
     * elements, making it difficult to compare test output on different
     * systems.
     */
     if(a1->start > a2->start) return -1;
     if(a1->start < a2->start) return 1;
     if(a1->length < a2->length) return -1;
     return 1;
}

/* Compare function for sorting primer records. */
static int
primer_pair_comp(x1, x2)
    const void *x1, *x2;
{
    const primer_pair *a1 = x1, *a2 = x2;
    int y1, y2;

    if(a1->pair_quality < a2->pair_quality) return -1;
    if (a1->pair_quality > a2->pair_quality) return 1;

    if (a1->compl_measure < a2->compl_measure) return -1;
    if (a1->compl_measure > a2->compl_measure) return 1;

    /* 
     * The following statements ensure that sorting
     * produces the same order on all systems regardless
     * of whether the sorting function is stable.
     */

    y1 = a1->left->start;
    y2 = a2->left->start;
    if (y1 > y2) return -1;  /* prefer left primers to the right. */
    if (y1 < y2) return 1;
    
    y1 = a1->right->start;
    y2 = a2->right->start;
    if (y1 < y2) return -1; /* prefer right primers to the left. */
    if (y1 > y2) return 1; 

    y1 = a1->left->length;
    y2 = a2->left->length;
    if (y1 < y2) return -1;  /* prefer shorter primers. */
    if (y1 > y2) return 1;
    
    y1 = a1->right->length;
    y2 = a2->right->length;
    if (y1 < y2) return -1; /* prefer shorter primers. */
    if (y1 > y2) return 1; 
    
    return 0;
}

/* 
 * Defines parameter values for given primer pair. Returns PAIR_OK if the pair is
 * acceptable; PAIR_FAILED otherwise.
 */
static int
pair_param(pa, sa, m, n, int_num, h, local_args, end_args, local_end_args)
    const primer_args *pa;
    seq_args *sa;
    int m, n, int_num;
    primer_pair *h;
    const dpal_args *local_args, *end_args, *local_end_args;
{
    char s1[MAX_PRIMER_LENGTH+1], s2[MAX_PRIMER_LENGTH+1], 
    s1_rev[MAX_PRIMER_LENGTH+1], s2_rev[MAX_PRIMER_LENGTH+1];
    short compl_end;

    /* FUTURE CODE: we must use the pair if the caller specifed
       both the left and the right primer. */
    int must_use = 0;

    int pair_failed_flag = 0;
    double min_oligo_tm;

    h->left = &f[m];
    h->right = &r[n];
    h->product_size = r[n].start - f[m].start+1;
    h->target = 0;
    h->compl_any = h->compl_end = 0;

    sa->pair_expl.considered++;

    if(h->product_size < pa->pr_min[int_num] || 
		h->product_size > pa->pr_max[int_num]) {
	sa->pair_expl.product++;
	h->product_size = -1;
	if (!must_use) return PAIR_FAILED;
	else pair_failed_flag = 1;
    }

    if (sa->num_targets > 0) {
	if (pair_spans_target(h, sa))
	    h->target = 1;
	else {
	    h->target = -1;
	    sa->pair_expl.target++;
	    if (!must_use) return PAIR_FAILED;
	    else pair_failed_flag = 1;
	}
    }

    /* Compute product Tm and related parameters; check constraints. */
   h->product_tm 
     = long_seq_tm(sa->trimmed_seq, h->left->start,
		   h->right->start - h->left->start + 1, pa->salt_conc, pa->divalent_conc, pa->dntp_conc);
      
    PR_ASSERT(h->product_tm != OLIGOTM_ERROR);

    min_oligo_tm 
      = h->left->temp > h->right->temp ? h->right->temp : h->left->temp;
    h->product_tm_oligo_tm_diff = h->product_tm - min_oligo_tm;
    h->t_opt_a  = 0.3 * min_oligo_tm + 0.7 * h->product_tm - 14.9;

    if (pa->product_min_tm != PR_DEFAULT_PRODUCT_MIN_TM
	&& h->product_tm < pa->product_min_tm) {
      sa->pair_expl.low_tm++;
      if (!must_use) return PAIR_FAILED;
      else pair_failed_flag = 1;
    }

    if (pa->product_max_tm != PR_DEFAULT_PRODUCT_MAX_TM
	&& h->product_tm > pa->product_max_tm) {
      sa->pair_expl.high_tm++;
      if (!must_use) return PAIR_FAILED;
      else pair_failed_flag = 1;
    }
      
    h->diff_tm = fabs(f[m].temp - r[n].temp);
    if (h->diff_tm > pa->max_diff_tm) {
	sa->pair_expl.temp_diff++;
	if (!must_use) return PAIR_FAILED;
	else pair_failed_flag = 1;
    }

    _pr_substr(sa->trimmed_seq,f[m].start,f[m].length,s1);
    _pr_substr(sa->trimmed_seq,r[n].start-r[n].length+1,r[n].length,s2);
    if(f[m].self_any == ALIGN_SCORE_UNDEF){
       oligo_compl(&f[m], pa, sa, OT_LEFT, local_args, end_args, local_end_args);
       if (!OK_OR_MUST_USE(&f[m])) {
	  sa->pair_expl.considered--;
	  return PAIR_FAILED;
       }
    }
    if(f[m].repeat_sim.score == NULL){
       oligo_mispriming(&f[m], pa, sa, OT_LEFT, local_end_args);
       if (!OK_OR_MUST_USE(&f[m])) {
	   sa->pair_expl.considered--;
	   return PAIR_FAILED;
       }
    }
    if(r[n].self_any == ALIGN_SCORE_UNDEF){
       oligo_compl(&r[n], pa, sa, OT_RIGHT, local_args, end_args, local_end_args);
       if (!OK_OR_MUST_USE(&r[n])) {
	  sa->pair_expl.considered--;
	  return PAIR_FAILED;
       }
    }
    if(r[n].repeat_sim.score == NULL){
       oligo_mispriming(&r[n], pa, sa, OT_RIGHT, local_end_args);
       if (!OK_OR_MUST_USE(&r[n])) {
	  sa->pair_expl.considered--;
	  return PAIR_FAILED;
       }
    }
	
    /* 
     * Similarity between s1 and s2 is equivalent to complementarity between
     * s2's complement and s1.  (Both s1 and s2 are taken from the same strand.)
     */
    h->compl_any = align(s1,s2,local_args);
    if (h->compl_any > pa->self_any) {
	sa->pair_expl.compl_any++;
	return PAIR_FAILED;
    }

    if ((h->compl_end = align(s1,s2,end_args)) > pa->self_end) {
	    sa->pair_expl.compl_end++;
	    return PAIR_FAILED;
    }
    /*
     * It is conceivable (though very unlikely in practice) that
     * align(s2_rev, s1_rev, end_args) > align(s1,s2,end_args).
     */
    _pr_reverse_complement(s1, s1_rev);
    _pr_reverse_complement(s2, s2_rev);
    if((compl_end = align(s2_rev, s1_rev, end_args)) > h->compl_end)  {
	if (compl_end > pa->self_end) {
	    sa->pair_expl.compl_end++;
	    return PAIR_FAILED;
	}
	h->compl_end = compl_end;
    }
    h->compl_measure = 	
	(h->right->self_end  + h->left->self_end + h->compl_end) * 1.1
	    + h->right->self_any + h->left->self_any + h->compl_any;
    if((h->repeat_sim = pair_repeat_sim(h, pa)) > pa->pair_repeat_compl){
	 sa->pair_expl.repeat_sim++;
	 return PAIR_FAILED;
    }

    if (!_pr_need_pair_template_mispriming(pa))
      h->template_mispriming = ALIGN_SCORE_UNDEF;
    else {
      PR_ASSERT(h->left->template_mispriming != ALIGN_SCORE_UNDEF);
      PR_ASSERT(h->left->template_mispriming_r != ALIGN_SCORE_UNDEF);
      PR_ASSERT(h->right->template_mispriming != ALIGN_SCORE_UNDEF);
      PR_ASSERT(h->right->template_mispriming_r != ALIGN_SCORE_UNDEF);
      h->template_mispriming =
	h->left->template_mispriming + h->right->template_mispriming_r;
      if ((h->left->template_mispriming_r + h->right->template_mispriming)
	  > h->template_mispriming)
      h->template_mispriming 
	= h->left->template_mispriming_r + h->right->template_mispriming;

      if (pa->pair_max_template_mispriming >= 0.0
	  && h->template_mispriming > pa->pair_max_template_mispriming) {
	sa->pair_expl.template_mispriming++;
	return PAIR_FAILED;
      }

    }
    return PAIR_OK;
}

/* compute_position_penalty is experimental code. */
void
compute_position_penalty(pa, sa, h, o_type)
  const primer_args *pa;
  const seq_args *sa;
  primer_rec *h;
  oligo_type o_type;
{
  int three_prime_base;
  int inside_flag = 0;
  int target_begin, target_end;

  PR_ASSERT(OT_LEFT == o_type || OT_RIGHT == o_type);
  PR_ASSERT(1 == sa->num_targets);
  target_begin = sa->tar[0][0];
  target_end = target_begin + sa->tar[0][1] - 1;

  three_prime_base = OT_LEFT == o_type
    ? h->start + h->length - 1 : h->start - h->length + 1;
  h->position_penalty_infinite = '\1';
  h->position_penalty = 0.0;
  
  if (OT_LEFT == o_type) {
    if (three_prime_base <= target_end) {
      h->position_penalty_infinite = '\0';
      if (three_prime_base < target_begin)
	h->position_penalty = target_begin - three_prime_base - 1;
      else {
	h->position_penalty = three_prime_base - target_begin + 1;
	inside_flag = 1;
      }
    }
  } else { /* OT_RIGHT == o_type */
    if (three_prime_base >= target_begin) {
      h->position_penalty_infinite = '\0';
      if (three_prime_base > target_end) {
	h->position_penalty = three_prime_base - target_end - 1;
      } else {
	h->position_penalty = target_end - three_prime_base + 1;
	inside_flag = 1;
      }
    }
  }
  if (!inside_flag)
    h->position_penalty *= pa->outside_penalty;
  else {
    h->position_penalty *= pa->inside_penalty;
  }
}

/* 
 * Return 1 if 'pair' spans any target (in sa->tar); otherwise return 0; An
 * oligo pair spans a target, t, if the last base of the left primer is to
 * left of the last base of t and the first base of the right primer is to
 * the right of the first base of t.  Of course the primers must
 * still be in a legal position with respect to each other.
 */
static int
pair_spans_target(pair, sa)
    const primer_pair *pair;
    const seq_args *sa;
{
    int i;
    int last_of_left = pair->left->start + pair->left->length - 1;
    int first_of_right = pair->right->start - pair->right->length + 1;
    int target_first, target_last;
    for (i = 0; i < sa->num_targets; i++) {
      target_first = sa->tar[i][0];
      target_last = target_first + sa->tar[i][1] - 1;
	if (last_of_left <= target_last
	    && first_of_right >= target_first
	    && last_of_left < first_of_right)
	    return 1;
    }
    return 0;
}

/* 
 * Return the value of the objective function value for given primer pair.
 * We must know that the pair is acceptable before calling obj_fn.
 */
static double
obj_fn(pa, h)
    const primer_args *pa;
    primer_pair *h;
{
    double sum;

    sum = 0.0;

    if(pa->pr_pair_weights.primer_quality)
       sum += pa->pr_pair_weights.primer_quality * (h->left->quality + h->right->quality);

    if(pa->pr_pair_weights.io_quality && 
        pa->primer_task == pick_pcr_primers_and_hyb_probe)
       sum += pa->pr_pair_weights.io_quality * h->intl->quality;

    if(pa->pr_pair_weights.diff_tm)
       sum += pa->pr_pair_weights.diff_tm * h->diff_tm;

    if(pa->pr_pair_weights.compl_any)
       sum += pa->pr_pair_weights.compl_any * h->compl_any / PR_ALIGN_SCORE_PRECISION;

    if(pa->pr_pair_weights.compl_end)
       sum += pa->pr_pair_weights.compl_end * h->compl_end / PR_ALIGN_SCORE_PRECISION;

    if(pa->pr_pair_weights.product_tm_lt && h->product_tm < pa->product_opt_tm)
	sum += pa->pr_pair_weights.product_tm_lt * 
			      (pa->product_opt_tm - h->product_tm);

    if(pa->pr_pair_weights.product_tm_gt && h->product_tm > pa->product_opt_tm)
	sum += pa->pr_pair_weights.product_tm_gt *
			       (h->product_tm - pa->product_opt_tm);

    if(pa->pr_pair_weights.product_size_lt &&
	    h->product_size < pa->product_opt_size) 
	    sum += pa->pr_pair_weights.product_size_lt * 
		 (pa->product_opt_size - h->product_size);

    if(pa->pr_pair_weights.product_size_gt &&
	    h->product_size > pa->product_opt_size)
	    sum += pa->pr_pair_weights.product_size_gt *
		  (h->product_size - pa->product_opt_size);

    if(pa->pr_pair_weights.repeat_sim)
      sum += pa->pr_pair_weights.repeat_sim * h->repeat_sim;

    if (pa->pr_pair_weights.template_mispriming) {
      PR_ASSERT(pa->pr_pair_weights.template_mispriming >= 0.0);
      PR_ASSERT(h->template_mispriming >= 0);
      sum += pa->pr_pair_weights.template_mispriming * h->template_mispriming;
    }

    PR_ASSERT(sum >= 0.0);

    return sum;
}

char *
pr_gather_warnings(sa, pa)
    const seq_args *sa;
    const primer_args *pa;
{
    pr_append_str warning;

    PR_ASSERT(NULL != sa);
    PR_ASSERT(NULL != pa);

    warning.data = NULL;
    warning.storage_size = 0;

    if (pa->repeat_lib.warning.data)
	pr_append_new_chunk(&warning, pa->repeat_lib.warning.data);

    if(pa->io_mishyb_library.warning.data != NULL) {
	pr_append_new_chunk(&warning, pa->io_mishyb_library.warning.data); 
	pr_append(&warning, " (for internal oligo)");
    }

    if (sa->warning.data) pr_append_new_chunk(&warning, sa->warning.data);
    return pr_is_empty(&warning) ? NULL : warning.data;
}

/* Print the data for chosen primer pairs to stdout in "boulderio" format. */
static void
boulder_print_pairs(prog_args, pa, sa, best_pairs)
    const program_args *prog_args;
    const primer_args *pa;
    const seq_args *sa;
    const pair_array_t *best_pairs;
{
    const char *left_tag, *right_tag, *intl_tag, *prod_size_tag;
    char *warning;
    char suffix [3];
    primer_rec *fwd, *rev, *intl;
    int i, incl_s = sa->incl_s;

    PR_ASSERT(NULL != pa);
    PR_ASSERT(NULL != sa);
    PR_ASSERT(NULL != prog_args);

    if (prog_args->twox_compat) {
	left_tag = "FORWARD_PRIMER";
	right_tag = "REVERSE_PRIMER";
	intl_tag = "MIDDLE_OLIGO";
	prod_size_tag = "PRODUCT_SIZE";
    } else {
	left_tag = "PRIMER_LEFT";
	right_tag = "PRIMER_RIGHT";
	intl_tag = "PRIMER_INTERNAL_OLIGO";
	prod_size_tag = "PRIMER_PRODUCT_SIZE";
    }

    if ((warning = pr_gather_warnings(sa, pa)) != NULL) {
      
	printf("PRIMER_WARNING=%s\n", warning);
	free(warning);
    }
    
    if (sa->error.data != NULL) {
	printf("PRIMER_ERROR=%s\n=\n", sa->error.data);
	if (fflush(stdout) == EOF) {
	    perror("fflush(stdout) failed");
	    exit(-1);
	}
	return;
    }

    if (pa->explain_flag) print_all_explain(pa, sa);

    if (!PR_START_CODON_POS_IS_NULL(sa))
      printf("PRIMER_STOP_CODON_POSITION=%d\n", sa->stop_codon_pos);

    for(i=0; i<best_pairs->num_pairs; i++) {
	fwd = best_pairs->pairs[i].left;
	rev = best_pairs->pairs[i].right;
	intl = best_pairs->pairs[i].intl;

	if (i == 0) suffix[0] = '\0';
	else sprintf(suffix, "_%d", i);

	printf("PRIMER_PAIR_PENALTY%s=%.4f\n", suffix,
	       best_pairs->pairs[i].pair_quality);
	printf("PRIMER_LEFT%s_PENALTY=%f\n", suffix,
	       fwd->quality);
	printf("PRIMER_RIGHT%s_PENALTY=%f\n", suffix,
		   rev->quality);
	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe)
	  printf("PRIMER_INTERNAL_OLIGO%s_PENALTY=%f\n", suffix,
		 intl->quality);

        /* Print sequences. */
	printf("PRIMER_LEFT%s_SEQUENCE=%s\n", suffix,
	       pr_oligo_sequence(sa, fwd));
	printf("PRIMER_RIGHT%s_SEQUENCE=%s\n", suffix,
	       pr_oligo_rev_c_sequence(sa, rev));
	if( pa->primer_task == pick_pcr_primers_and_hyb_probe)
#if USE_OLD_FORMAT_MISTAKE
	    printf("PRIMER_INTERNAL%s_OLIGO_SEQUENCE=%s\n", suffix,
		   pr_oligo_sequence(sa,intl));
#else
	    printf("PRIMER_INTERNAL_OLIGO%s_SEQUENCE=%s\n", suffix,
		   pr_oligo_sequence(sa,intl));
#endif

	printf("%s%s=%d,%d\n", left_tag, suffix,
	       fwd->start + incl_s + pa->first_base_index,
	       fwd->length);
	printf("%s%s=%d,%d\n", right_tag, suffix,
	       rev->start + incl_s + pa->first_base_index,
	       rev->length);
	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe)
	    printf("%s%s=%d,%d\n", intl_tag, suffix,
		   intl->start + incl_s + pa->first_base_index,
		   intl->length);

	printf("PRIMER_LEFT%s_TM=%.3f\n", suffix, fwd->temp);
	printf("PRIMER_RIGHT%s_TM=%.3f\n", suffix, rev->temp);
	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe)
	    printf("PRIMER_INTERNAL_OLIGO%s_TM=%.3f\n",suffix, intl->temp);


	printf("PRIMER_LEFT%s_GC_PERCENT=%.3f\n", suffix, fwd->gc_content);
	printf("PRIMER_RIGHT%s_GC_PERCENT=%.3f\n", suffix, rev->gc_content);
	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe)
	    printf("PRIMER_INTERNAL_OLIGO%s_GC_PERCENT=%.3f\n",suffix,
		   intl->gc_content);

	printf("PRIMER_LEFT%s_SELF_ANY=%.2f\n", suffix,
	       fwd->self_any / PR_ALIGN_SCORE_PRECISION);
	printf("PRIMER_RIGHT%s_SELF_ANY=%.2f\n", suffix,
	       rev->self_any / PR_ALIGN_SCORE_PRECISION);
	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe)
	    printf("PRIMER_INTERNAL_OLIGO%s_SELF_ANY=%.2f\n", suffix,
		   intl->self_any / PR_ALIGN_SCORE_PRECISION);
	printf("PRIMER_LEFT%s_SELF_END=%.2f\n", suffix,
	       fwd->self_end / PR_ALIGN_SCORE_PRECISION);
	printf("PRIMER_RIGHT%s_SELF_END=%.2f\n",
	       suffix,rev->self_end / PR_ALIGN_SCORE_PRECISION);
	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe)
	    printf("PRIMER_INTERNAL_OLIGO%s_SELF_END=%.2f\n", suffix,
		   intl->self_end / PR_ALIGN_SCORE_PRECISION);
        if (pa->repeat_lib.seq_num > 0) {
	    printf("PRIMER_LEFT%s_MISPRIMING_SCORE=%.2f, %s\n", suffix,
		   fwd->repeat_sim.score[fwd->repeat_sim.max] / PR_ALIGN_SCORE_PRECISION,
		   fwd->repeat_sim.name);
            printf("PRIMER_RIGHT%s_MISPRIMING_SCORE=%.2f, %s\n", suffix,
		   rev->repeat_sim.score[rev->repeat_sim.max] / PR_ALIGN_SCORE_PRECISION,
                   rev->repeat_sim.name);
            printf("PRIMER_PAIR%s_MISPRIMING_SCORE=%.2f, %s\n", suffix,
		   best_pairs->pairs[i].repeat_sim / PR_ALIGN_SCORE_PRECISION,
		   best_pairs->pairs[i].rep_name);
        }
	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe
	    && pa->io_mishyb_library.seq_num > 0)
	    printf("PRIMER_INTERNAL_OLIGO%s_MISHYB_SCORE=%.2f, %s\n", suffix,
		   intl->repeat_sim.score[intl->repeat_sim.max]
		   / PR_ALIGN_SCORE_PRECISION,
		   intl->repeat_sim.name);
	if (NULL != sa->quality){
	   printf("PRIMER_LEFT%s_MIN_SEQ_QUALITY=%d\n", suffix,
		   fwd->seq_quality);
           printf("PRIMER_RIGHT%s_MIN_SEQ_QUALITY=%d\n", suffix,
		   rev->seq_quality);
        }

	if (!_PR_DEFAULT_POSITION_PENALTIES(pa)
            || !PR_START_CODON_POS_IS_NULL(sa)) {
	  printf("PRIMER_LEFT%s_POSITION_PENALTY=%f\n", suffix,
		 fwd->position_penalty);
	  printf("PRIMER_RIGHT%s_POSITION_PENALTY=%f\n", suffix,
		 rev->position_penalty);
	}

	printf("PRIMER_LEFT%s_END_STABILITY=%.4f\n",
	       suffix, fwd->end_stability);
	printf("PRIMER_RIGHT%s_END_STABILITY=%.4f\n",
	       suffix, rev->end_stability);

	if (oligo_max_template_mispriming(fwd) != ALIGN_SCORE_UNDEF)
	  printf("PRIMER_LEFT%s_TEMPLATE_MISPRIMING=%.4f\n", suffix,
		 oligo_max_template_mispriming(fwd)
		 / PR_ALIGN_SCORE_PRECISION);

	if (oligo_max_template_mispriming(rev) != ALIGN_SCORE_UNDEF)
	  printf("PRIMER_RIGHT%s_TEMPLATE_MISPRIMING=%.4f\n", suffix,
		 oligo_max_template_mispriming(rev)
		 / PR_ALIGN_SCORE_PRECISION);


	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe
	    && NULL != sa->quality)
	   printf("PRIMER_INTERNAL_OLIGO%s_MIN_SEQ_QUALITY=%d\n",
		   suffix, intl->seq_quality);
	printf("PRIMER_PAIR%s_COMPL_ANY=%.2f\n", suffix,
	       best_pairs->pairs[i].compl_any / PR_ALIGN_SCORE_PRECISION);
	printf("PRIMER_PAIR%s_COMPL_END=%.2f\n", suffix,
	       best_pairs->pairs[i].compl_end  / PR_ALIGN_SCORE_PRECISION);

	printf("%s%s=%d\n", prod_size_tag, suffix,
				     best_pairs->pairs[i].product_size);

	if (pa->product_max_tm != PR_DEFAULT_PRODUCT_MAX_TM
	    || pa->product_min_tm != PR_DEFAULT_PRODUCT_MIN_TM) {
	  printf("PRIMER_PRODUCT_TM%s=%.4f\n", suffix,
		 best_pairs->pairs[i].product_tm);

	  printf("PRIMER_PRODUCT_TM_OLIGO_TM_DIFF%s=%.4f\n", suffix,
	      best_pairs->pairs[i].product_tm_oligo_tm_diff);

	   printf("PRIMER_PAIR%s_T_OPT_A=%.4f\n", suffix,
	      best_pairs->pairs[i].t_opt_a);
	}

	if (best_pairs->pairs[i].template_mispriming != ALIGN_SCORE_UNDEF)
	  printf("PRIMER_PAIR%s_TEMPLATE_MISPRIMING=%.2f\n", suffix,
		 best_pairs->pairs[i].template_mispriming 
		 / PR_ALIGN_SCORE_PRECISION);

    }
    printf("=\n");
    if (fflush(stdout) == EOF) {
	perror("fflush(stdout) failed");
	exit(-1);
    }
}

static void 
boulder_print_oligos(pa, sa, n, l)
    const primer_args *pa;
    const seq_args *sa;
    int n;
    oligo_type l;
{
    char *warning;
    int i, j;
    char suffix [3], type[256];
    /* type must be larger than the length of "PRIMER_INTERNAL_OLIGO". */

    primer_rec *oligo;
    int incl_s = sa->incl_s;


    if ((warning = pr_gather_warnings(sa, pa)) != NULL) {
	printf("PRIMER_WARNING=%s\n", warning);
	free(warning);
    }
    if (sa->error.data != NULL) {
	printf("PRIMER_ERROR=%s\n=\n", sa->error.data);
	if (fflush(stdout) == EOF) {
	    perror("fflush(stdout) failed");
	    exit(-1);
        }
	return;
    }

    if(l == OT_LEFT) strcpy(type, "PRIMER_LEFT");
    else if(l == OT_RIGHT) strcpy(type, "PRIMER_RIGHT");
    else strcpy(type, "PRIMER_INTERNAL_OLIGO");

    if (pa->explain_flag) print_all_explain(pa, sa);

    i = 0;
    j = (pa->num_return < n) ? pa->num_return : n;
    if(l == OT_LEFT) oligo = f;
    else if (l == OT_RIGHT) oligo = r;
    else  oligo = mid;

    while(i < j) {
	if (i == 0) suffix[0] = '\0';
	else sprintf(suffix, "_%d", i);

	printf("%s%s_PENALTY=%.4f\n", type, suffix, oligo[i].quality);
	if(l == OT_RIGHT)
	   printf("%s%s_SEQUENCE=%s\n", type, suffix,
		  pr_oligo_rev_c_sequence(sa, &oligo[i]));
        else printf("%s%s_SEQUENCE=%s\n", type, suffix,
		    pr_oligo_sequence(sa, &oligo[i])); 
	printf("%s%s=%d,%d\n", type, suffix, 
		       oligo[i].start + incl_s + pa->first_base_index,
		       oligo[i].length);
	printf("%s%s_TM=%.3f\n", type, suffix, oligo[i].temp);
	printf("%s%s_GC_PERCENT=%.3f\n", type, suffix, oligo[i].gc_content);
        printf("%s%s_SELF_ANY=%.2f\n", type, suffix,
		       oligo[i].self_any / PR_ALIGN_SCORE_PRECISION);
        printf("%s%s_SELF_END=%.2f\n", type, suffix,
		       oligo[i].self_end / PR_ALIGN_SCORE_PRECISION);
        if ((l == OT_LEFT || l == OT_RIGHT) && pa->repeat_lib.seq_num > 0 ) 
	    printf("%s%s_MISPRIMING_SCORE=%.2f, %s\n", type, suffix,
		    oligo[i].repeat_sim.score[oligo[i].repeat_sim.max] /PR_ALIGN_SCORE_PRECISION,
		    oligo[i].repeat_sim.name);
        if (l == OT_INTL && pa->io_mishyb_library.seq_num > 0)
	    printf("%s%s_MISHYB_SCORE=%.2f,%s\n", type, suffix,
	           oligo[i].repeat_sim.score[oligo[i].repeat_sim.max]/ PR_ALIGN_SCORE_PRECISION,
                   oligo[i].repeat_sim.name);
	if (NULL != sa->quality)printf("%s%s_MIN_SEQ_QUALITY=%d\n",
		   type, suffix, oligo[i].seq_quality);
        if (PR_DEFAULT_INSIDE_PENALTY != pa->inside_penalty
	    || PR_DEFAULT_OUTSIDE_PENALTY != pa->outside_penalty != 0.0)
	    printf("%s%s_POSITION_PENALTY=%f\n", type, suffix,
		    oligo[i].position_penalty);
	if(l == OT_LEFT || l == OT_RIGHT)
	    printf("%s%s_END_STABILITY=%.4f\n", type, suffix,
		    oligo[i].end_stability);

	if (oligo_max_template_mispriming(&oligo[i]) != ALIGN_SCORE_UNDEF)
	  printf("%s%s_TEMPLATE_MISPRIMING=%.2f\n", type, suffix,
		 oligo_max_template_mispriming(&oligo[i])
		 / PR_ALIGN_SCORE_PRECISION);

        i++;
    }
    printf("=\n");
}

static void
print_all_explain(pa, sa)
    const primer_args *pa;
    const seq_args *sa;
{
  if (pa->explain_flag) {
    if (pa->primer_task != pick_right_only
	&& pa->primer_task != pick_hyb_probe_only
	&& !(pa->pick_anyway && sa->left_input))
      print_explain(&sa->left_expl,OT_LEFT);

    if (pa->primer_task != pick_left_only 
	&& pa->primer_task != pick_hyb_probe_only
	&& !(pa->pick_anyway && sa->right_input))
      print_explain(&sa->right_expl,OT_RIGHT);

    if ((pa->primer_task == pick_hyb_probe_only
	 || pa->primer_task == pick_pcr_primers_and_hyb_probe)
	&& !(pa->pick_anyway && sa->internal_input)) 
      print_explain(&sa->intl_expl, OT_INTL);

    if (pa->primer_task  == pick_pcr_primers
	|| pa->primer_task == pick_pcr_primers_and_hyb_probe) {
      printf("PRIMER_PAIR_EXPLAIN=");
      pr_print_pair_explain(stdout, sa);
    }
  }
}

void
pr_print_pair_explain(f, sa)
  FILE *f;
  const seq_args *sa;
{
    fprintf(f, "considered %d",sa->pair_expl.considered);
    if (sa->pair_expl.target)
      fprintf(f, ", no target %d", sa->pair_expl.target);
    if (sa->pair_expl.product)
      fprintf(f, ", unacceptable product size %d", sa->pair_expl.product);
    if (sa->pair_expl.low_tm)
      fprintf(f, ", low product Tm %d", sa->pair_expl.low_tm);
    if (sa->pair_expl.high_tm)
      fprintf(f, ", high product Tm %d", sa->pair_expl.high_tm);
    if (sa->pair_expl.temp_diff) 
      fprintf(f, ", tm diff too large %d",sa->pair_expl.temp_diff);
    if (sa->pair_expl.compl_any) 
      fprintf(f, ", high any compl %d", sa->pair_expl.compl_any);
    if (sa->pair_expl.compl_end) 
      fprintf(f, ", high end compl %d", sa->pair_expl.compl_end);
    if (sa->pair_expl.internal) 
      fprintf(f, ", no internal oligo %d", sa->pair_expl.internal);
    if (sa->pair_expl.repeat_sim)
      fprintf(f, ", high mispriming library similarity %d",
	      sa->pair_expl.repeat_sim);
    if (sa->pair_expl.template_mispriming)
      fprintf(f, ", high template mispriming score %d",
	      sa->pair_expl.template_mispriming);
    fprintf(f, ", ok %d\n", sa->pair_expl.ok);
}

static void
print_explain(stat, l)
    const oligo_stats *stat;
    oligo_type l;
{
    if(OT_LEFT == l)printf("PRIMER_LEFT_EXPLAIN=");
    else if(OT_RIGHT == l)printf("PRIMER_RIGHT_EXPLAIN=");
	 else printf("PRIMER_INTERNAL_OLIGO_EXPLAIN=");

    printf("considered %d", stat->considered);
    if (stat->no_orf) printf(", would not amplify any of the ORF %d", stat->no_orf);
    if(stat->ns)printf(", too many Ns %d", stat->ns);
    if(stat->target)printf(", overlap target %d", stat->target);
    if(stat->excluded)printf(", overlap excluded region %d", stat->excluded);
    if(stat->gc)printf(", GC content failed %d", stat->gc);
    if(stat->gc_clamp)printf(", GC clamp failed %d", stat->gc_clamp);
    if(stat->temp_min)printf(", low tm %d", stat->temp_min);
    if(stat->temp_max)printf(", high tm %d", stat->temp_max);
    if(stat->compl_any)printf(", high any compl %d", stat->compl_any);
    if(stat->compl_end)printf(", high end compl %d", stat->compl_end);
    if(stat->repeat_score) printf(", high repeat similarity %d", stat->repeat_score);
    if(stat->poly_x)printf(", long poly-x seq %d", stat->poly_x);
    if(stat->seq_quality)printf(",low sequence quality %d", stat->seq_quality);
    if (stat->stability) printf(",high 3' stability %d", stat->stability);
    if (stat->template_mispriming) printf(",high template mispriming score %d",
					  stat->template_mispriming);
    /* edited by T. Koressaar for lowercase masking */
    if(stat->gmasked) printf(",lowercase masking of 3' end %d",stat->gmasked);
   
    printf(", ok %d\n", stat->ok);
}

void
pr_append(x, s)
    pr_append_str *x;
    const char *s;
{
    int xlen, slen;
    if (NULL == x->data) {
	x->storage_size = 24;
	x->data = pr_safe_malloc(x->storage_size);
	*x->data = '\0';
    }
    xlen = strlen(x->data);
    slen = strlen(s);
    if (xlen + slen + 1 > x->storage_size) {
	x->storage_size += 2 * (slen + 1);
	x->data = pr_safe_realloc(x->data, x->storage_size);
    }
    strcpy(x->data + xlen, s);
}

void
pr_append_new_chunk(x, s)
    pr_append_str *x;
    const char *s;
{
  pr_append_w_sep(x, "; ", s);
}

void
pr_append_w_sep(x, sep, s)
    pr_append_str *x;
    const char *sep;
    const char *s;
{
    if (pr_is_empty(x))
	pr_append(x, s);
    else {
	pr_append(x, sep);
	pr_append(x, s);
    }
}

void
pr_set_empty(x)
    pr_append_str *x;
{
    PR_ASSERT(NULL != x);
    if (NULL != x->data) *x->data = '\0';
}

int
pr_is_empty(x)
    const pr_append_str *x;
{
    PR_ASSERT(NULL != x);
    return  NULL == x->data || '\0' == *x->data;
}

static short
align(s1, s2, a)
    const char *s1, *s2;
    const dpal_args *a;
{
    dpal_results r;
    dpal(s1, s2, a, &r);
    PR_ASSERT(r.score <= SHRT_MAX);
    return ((r.score<0) ? 0 : (short)r.score);
}

/* Set dpal args to appropriate values for primer picking. */
static void
set_dpal_args(a)
    dpal_args *a;
{
    unsigned int i, j;

    memset(a, 0, sizeof(*a));
    for (i = 0; i <= UCHAR_MAX; i++)
	for (j = 0; j <= UCHAR_MAX; j++)
	  if (('A' == i || 'C' == i || 'G' == i || 'T' == i || 'N' == i)
	      && ('A' == j || 'C' == j || 'G' == j || 'T' == j 
		  || 'N' == j)) {
	    if (i == 'N' || j == 'N') 
	      a->ssm[i][j] = -25;
	    else if (i == j)
	      a->ssm[i][j] = 100;
	    else 
	      a->ssm[i][j] = -100;
	  } else
	    a->ssm[i][j] = INT_MIN;

    a->gap                = -200;
    a->gapl               = -200;
    a->flag               = DPAL_LOCAL;
    a->max_gap            = 1;
    a->fail_stop          = 1;
    a->check_chars        = 0;
    a->debug              = 0;
    a->score_only         = 1;
    a->force_generic      = 0;
    a->force_long_generic = 0;
    a->force_long_maxgap1 = 0;
}

char *
pr_oligo_sequence(sa, o)
    const seq_args *sa;
    const primer_rec *o;
{
    static char s[MAX_PRIMER_LENGTH+1];
    int seq_len;
    PR_ASSERT(NULL != sa);
    PR_ASSERT(NULL != o);
    seq_len = strlen(sa->sequence);
    PR_ASSERT(o->start + sa->incl_s >= 0);
    PR_ASSERT(o->start + sa->incl_s + o->length <= seq_len);
    _pr_substr(sa->sequence, sa->incl_s + o->start, o->length, s);
    return &s[0];
}

char *
pr_oligo_rev_c_sequence(sa, o)
    const seq_args *sa;
    const primer_rec *o;
{
    static char s[MAX_PRIMER_LENGTH+1], s1[MAX_PRIMER_LENGTH+1];
    int seq_len, start;
    PR_ASSERT(NULL != sa);
    PR_ASSERT(NULL != o);
    seq_len = strlen(sa->sequence);
    start = sa->incl_s + o->start - o->length + 1;
    PR_ASSERT(start >= 0);
    PR_ASSERT(start + o->length <= seq_len);
    _pr_substr(sa->sequence, start, o->length, s);
    _pr_reverse_complement(s,s1);
    return &s1[0];
}

static void
print_usage()
{
    const char **p;
    p = &copyright[0];
    while (NULL != *p) fprintf(stderr, "%s\n", *p++);
    fprintf(stderr, 
	    "\n\nUSAGE: %s %s %s\n", pr_program_name,
	    "[-format_output] [-2x_compat]",
	    "[-strict_tags]");
    fprintf(stderr, "This is primer3 (%s)\n", pr_release);
    fprintf(stderr, "Input must be provided on standard input.\n");
    fprintf(stderr, "For example:\n");
    fprintf(stderr, "$ primer3_core < my_input_file\n");
}

static void
sig_handler(signal)
    int signal;
{
    print_usage();
    fprintf(stderr, "%s: received signal %d\n", pr_program_name, signal);
    exit(signal);
}

/* Calculate self complementarity. */
static void
oligo_compl(h, ha, sa, l, local_args, end_args, local_end_args)
    primer_rec *h;
    const primer_args *ha;
    seq_args *sa;
    oligo_type l;
    const dpal_args *local_args, *end_args, *local_end_args;
{
    char s[MAX_PRIMER_LENGTH+1], s1[MAX_PRIMER_LENGTH+1];
    int j;
    short self_any, self_end;

    if (OT_INTL == l) {
      self_any = ha->io_self_any;
      self_end = ha->io_self_end;
    } else {
      self_any = ha->self_any;
      self_end = ha->self_end;
    }

    j =  (OT_LEFT == l || OT_INTL == l) ? h->start : h->start-h->length+1;

    _pr_substr(sa->trimmed_seq, j, h->length, s1);
    _pr_reverse_complement(s1, s);

    h->self_any = align(s1, s, local_args);
    if(h->self_any > self_any){
	h->ok = OV_SELF_ANY;
	if      (OT_LEFT  == l) {
	     sa->left_expl.compl_any++;
	     sa->left_expl.ok--;
        }
	else if (OT_RIGHT == l) {
	     sa->right_expl.compl_any++;
	     sa->right_expl.ok--;
        }
	else {
	     sa->intl_expl.compl_any++;
	     sa->intl_expl.ok--;
        }
	if (!h->must_use) return;
    }

    h->self_end = (l != OT_RIGHT) ? align(s1,s,end_args) : align(s,s1,end_args);
    if(h->self_end > self_end) {
	h->ok = OV_SELF_END;
	if      (OT_LEFT  == l) {
	    sa->left_expl.compl_end++;
	    sa->left_expl.ok--;
        }
	else if (OT_RIGHT == l) {
	    sa->right_expl.compl_end++;
	    sa->right_expl.ok--;
        }
	else {
	    sa->intl_expl.compl_end++;
	    sa->intl_expl.ok--;
        }
	return;
    }
}

static void 
oligo_mispriming(h, pa, sa, l, align_args)
   primer_rec *h;
   const primer_args *pa;
   seq_args *sa;
   oligo_type l;
   const dpal_args *align_args;
{
  char 
    s[MAX_PRIMER_LENGTH+1],     /* Will contain the oligo sequence. */
    s_tmp[MAX_PRIMER_LENGTH+1], /* Scratch buffer. */
    s_r[MAX_PRIMER_LENGTH+1];   /* Will contain s reverse complemented. */

  char *oseq, *target, *target_r;
  double w;
  const seq_lib *lib;
  int i;
  int first, last; /* Indexes of first and last bases of the oligo in sa->trimmed_seq,
		     that is, WITHIN THE INCLUDED REGION. */
  int first_untrimmed, last_untrimmed;  
                  /* Indexes of first and last bases of the oligo in sa->seq,
		     that is, WITHIN THE TOTAL SEQUENCE INPUT. */
  int min, max, tmp;
  int seqlen;
  int debug = 0;
  int match_length;
  short  lib_compl;
  short  tmp_score;
  char   tmp_char;

  if (OT_INTL == l) {
    lib = &(pa->io_mishyb_library);
    lib_compl = pa->io_repeat_compl;
  } else {
    lib = &(pa->repeat_lib);
    lib_compl = pa->repeat_compl;
  }

  first =  (OT_LEFT == l || OT_INTL == l)
    ? h->start 
    : h->start - h->length + 1;
  last  =  (OT_LEFT == l || OT_INTL == l)
    ? h->start + h->length - 1
    : h->start;

  match_length = h->length;

  _pr_substr(sa->trimmed_seq, first, h->length, s_tmp);
  _pr_substr(s_tmp, 0, match_length, s);
  _pr_reverse_complement(s, s_tmp);  /* FIX ME -- is s_tmp needed? */
  _pr_substr(s_tmp, 0, match_length, s_r);

  /*
   * Calculate maximum similarity to sequences from user defined repeat
   * library. Compare it with maximum allowed repeat similarity.
   */
  if(lib->seq_num > 0) {
    h->repeat_sim.score = 
      pr_safe_malloc(lib->seq_num * sizeof(short));
    h->repeat_sim.max = h->repeat_sim.min = 0;
    max = min = 0;
    h->repeat_sim.name = lib->names[0];
    for(i = 0; i < lib->seq_num; i++){
      if (OT_LEFT == l) w = lib->weight[i] *
			  align(s, lib->seqs[i], lib_local_end_dpal_args);
      else if (OT_INTL == l) w = lib->weight[i] *
			       align(s, lib->seqs[i], lib_local_dpal_args);
      else w = lib->weight[i] *
	     align(s_r, lib->rev_compl_seqs[i], lib_local_end_dpal_args);

      h->repeat_sim.score[i] = w;
      if(w > max){
	max = w;
	h->repeat_sim.max = i;
	h->repeat_sim.name = lib->names[i];
      }
      if(w < min){
	min = w;
	h->repeat_sim.min = i;
      }
      if (w > lib_compl) {
	h->ok = OV_LIB_SIM;
	if (OT_LEFT  == l) {
	  sa->left_expl.repeat_score++;
	  sa->left_expl.ok--;
	}
	else if (OT_RIGHT == l) {
	  sa->right_expl.repeat_score++;
	  sa->right_expl.ok--;
	}
	else {
	  sa->intl_expl.repeat_score++;
	  sa->intl_expl.ok--;
	}
	if (!h->must_use) return;
      }
    }
  }

  if (_pr_need_template_mispriming(pa) && (l == OT_RIGHT || l == OT_LEFT)) {
    /* Calculate maximum similarity to ectopic sites in the template. */

    seqlen = strlen(sa->upcased_seq);
    first_untrimmed = sa->incl_s + first;
    last_untrimmed = sa->incl_s + last;

    if (l == OT_LEFT) {
      oseq = &s[0];
      target = &sa->upcased_seq[0];
      target_r = &sa->upcased_seq_r[0];
    } else {  /* l == OT_RIGHT */
      if (debug) 
	fprintf(stderr, "first_untrimmed = %d, last_untrimmed = %d\n",
		first_untrimmed, last_untrimmed);
      oseq = &s_r[0];
      target = &sa->upcased_seq_r[0];
      target_r = &sa->upcased_seq[0];
      /* We need to adjust first_untrimmed and last_untrimmed so that
	 they are correct in the reverse-complemented
	 sequence.
      */
      tmp = (seqlen - last_untrimmed) - 1;
      last_untrimmed  = (seqlen - first_untrimmed) - 1;
      first_untrimmed = tmp;
    }

    /* 1. Align to the template 5' of the oligo. */
    tmp_char = target[first_untrimmed];
    target[first_untrimmed] = '\0';

    tmp_score = align(oseq, target, align_args);

    if (debug) {
      if (l == OT_LEFT) fprintf(stderr, "\n************ OLIGO = LEFT\n");
      else fprintf(stderr,              "\n************ OLIGO = RIGHT\n");
      fprintf(stderr, "first_untrimmed = %d, last_untrimmed = %d\n",
	      first_untrimmed, last_untrimmed);
			
      fprintf(stderr, "5' of oligo: Score %d aligning %s against %s\n\n", tmp_score,
	      oseq, target);
    }

    target[first_untrimmed] = tmp_char;

    /* 2. Align to the template 3' of the oligo. */
    h->template_mispriming
      = align(oseq, &target[0] + last_untrimmed + 1, align_args);

    if (debug)
      fprintf(stderr, "3' of oligo Score %d aligning %s against %s\n\n",
	      h->template_mispriming, oseq, &target[0] + last_untrimmed + 1);

    /* 3. Take the max of 1. and 2. */
    if (tmp_score > h->template_mispriming)
      h->template_mispriming = tmp_score;

    /* 4. Align to the reverse strand of the template. */
    h->template_mispriming_r
      = align(oseq, target_r, align_args);

    if (debug)
      fprintf(stderr, "other strand Score %d aligning %s against %s\n\n", 
	      h->template_mispriming_r, oseq, target_r);


    if (pa->max_template_mispriming >= 0 
	 && oligo_max_template_mispriming(h) > pa->max_template_mispriming) {
      h->ok = OV_TEMPLATE_MISPRIMING;
      if (OT_LEFT == l) {
	sa->left_expl.template_mispriming++;
	sa->left_expl.ok--;
      } else if (OT_RIGHT == l) {
	sa->right_expl.template_mispriming++;
	sa->right_expl.ok--;
      } else PR_ASSERT(0); /* Should not get here. */
      if (!h->must_use) return;
    }
  }
}

static int
pair_repeat_sim(h, pa)
  primer_pair *h;
  const primer_args *pa;
{
  int i, n, max, w;
  primer_rec *fw, *rev;
     
  fw = h->left;
  rev = h->right;

  max = 0;
  n = pa->repeat_lib.seq_num;
  if(n == 0) return 0;
  h->rep_name =  pa->repeat_lib.names[0] ;
  for(i = 0; i < n; i++) {
    if((w=(fw->repeat_sim.score[i] +
	   rev->repeat_sim.score[i])) > max) {
      max = w;
      h->rep_name =  pa->repeat_lib.names[i] ;
    }
  }
  return max;
}

/* 
 * 's' is the sequence in which to find the stop codon.
 * 'start' is the position of a start codon.
 * 
 * There are two modes depending on 'direction'.
 *
 * If direction is 1:
 *
 * Return the index of the first stop codon to the right of
 * 'start' in 's' (in same frame).
 * If there is no such stop codon return -1.
 *
 * If direction is -1:
 *
 * If 'start' is negative then return -1.
 * Otherwise return the index the first stop codon to left of 'start'.
 * If there is no such stop codon return -1.
 *  
 * Note: we don't insist that the start codon be ATG, since in some
 * cases the caller will not have the full sequence in 's', nor even
 * know the postion of the start codon relative to s.
 */
int
find_stop_codon(s, start, direction)
  const char* s;
{
  const char *p, *q;
  int increment = 3 * direction;
  int len = strlen(s);

  PR_ASSERT(s != NULL);
  PR_ASSERT(direction == 1 || direction == -1);
  PR_ASSERT(len >= 3);
  PR_ASSERT(start <= (len - 3));

  if (start < 0) {
    if (direction == 1)
      while (start < 0) start += increment;
    else
      return -1;
  }

  for (p = &s[start];
       p >= &s[0]
	 && *p
	 && *(p + 1)
	 && *(p + 2);
       p += increment) {
    if ('T' != *p && 't' != *p) continue;
    q = p + 1;
    if ('A' == *q || 'a' == *q) {
      q++;
      if  ('G' == *q || 'g' == *q || 'A' == *q || 'a' == *q)
	return p - &s[0];
    } else if ('G' == *q || 'g' == *q) {
      q++;
      if ('A' == *q || 'a' == *q)
	return p - &s[0];
    }
  }
  return -1;
}

int
strcmp_nocase(s1, s2)
char *s1, *s2;
{
   static char M[UCHAR_MAX];
   static int f = 0;
   int i;
   char *p, *q;

   if(f != 1){
      for(i = 0; i < UCHAR_MAX; i++) M[i] = i;
      i = 'a'; M[i] = 'A'; i = 'b'; M[i] = 'B'; i = 'c'; M[i] = 'C';
      i = 'A'; M[i] = 'a'; i = 'B'; M[i] = 'b'; i = 'C'; M[i] = 'c';
      i = 'd'; M[i] = 'D'; i = 'e'; M[i] = 'E'; i = 'f'; M[i] = 'F';
      i = 'D'; M[i] = 'd'; i = 'E'; M[i] = 'e'; i = 'F'; M[i] = 'f';
      i = 'g'; M[i] = 'G'; i = 'h'; M[i] = 'H'; i = 'i'; M[i] = 'I';
      i = 'G'; M[i] = 'g'; i = 'H'; M[i] = 'h'; i = 'I'; M[i] = 'i';
      i = 'k'; M[i] = 'K'; i = 'l'; M[i] = 'L'; i = 'm'; M[i] = 'M';
      i = 'K'; M[i] = 'k'; i = 'L'; M[i] = 'l'; i = 'M'; M[i] = 'm';
      i = 'n'; M[i] = 'N'; i = 'o'; M[i] = 'O'; i = 'p'; M[i] = 'P';
      i = 'N'; M[i] = 'n'; i = 'O'; M[i] = 'o'; i = 'P'; M[i] = 'p';
      i = 'q'; M[i] = 'Q'; i = 'r'; M[i] = 'R'; i = 's'; M[i] = 'S';
      i = 'Q'; M[i] = 'q'; i = 'R'; M[i] = 'r'; i = 'S'; M[i] = 's';
      i = 't'; M[i] = 'T'; i = 'u'; M[i] = 'U'; i = 'v'; M[i] = 'V';
      i = 'T'; M[i] = 't'; i = 'U'; M[i] = 'u'; i = 'V'; M[i] = 'v';
      i = 'w'; M[i] = 'W'; i = 'x'; M[i] = 'X'; i = 'y'; M[i] = 'Y';
      i = 'W'; M[i] = 'w'; i = 'X'; M[i] = 'x'; i = 'Y'; M[i] = 'y';
      i = 'z'; M[i] = 'Z'; i = 'Z'; M[i] = 'z'; i = 'j'; M[i] = 'J';
      i = 'J'; M[i] = 'j';
      f = 1;
   }

   if(s1 == NULL || s2 == NULL) return 1;
   if(strlen(s1) != strlen(s2)) return 1;
   p = s1; q = s2;
   while(*p != '\0' && *p != '\n' && *q != '\0' && *q != '\n'){
      i = *p;
      if(*p == *q || M[i] == *q ) {p++; q++; continue;}

      return 1;
   }
   return 0;
}

static void
free_repeat_sim_score(n_f, n_r, n_m)
int n_f, n_r, n_m;
{
   int i;

   for(i = 0; i < n_f; i++) if(f[i].repeat_sim.score != NULL)
     { free(f[i].repeat_sim.score); f[i].repeat_sim.score = NULL; }

   for(i = 0; i < n_r; i++) if(r[i].repeat_sim.score != NULL)
     { free(r[i].repeat_sim.score); r[i].repeat_sim.score = NULL; }

   for(i = 0; i < n_m; i++) if(mid[i].repeat_sim.score != NULL)
     { free(mid[i].repeat_sim.score); mid[i].repeat_sim.score = NULL; }
}

/*  Edited by T. Koressaar for lowercase masking. This function checks
 if the 3' end of the primer has been masked by lowercase letter.
 Function created/Added by Eric Reppo, July 9, 2002
 */
static void
check_if_lowercase_masked(position, sequence, h)
     const int position;
     const char *sequence;
     primer_rec *h;
{   
   const char* p = &sequence[position];
   if ('a' == *p || 'c' == *p ||'g' == *p || 't' == *p) {
      h->ok=OV_GMASKED;
   }
}

/* =========================================================== */
/* Various fail-stop wrappers for standard library functions.  */

void *
pr_safe_malloc(x)
    size_t x;
{
    void *r = malloc(x);
    if (NULL == r) OOM_ERROR;
    return r;
}

void *
pr_safe_realloc(p, x)
    void *p;
    size_t x;
{
    void *r = realloc(p, x);
    if (NULL == r) OOM_ERROR;
    return r;
}

static FILE *
safe_fopen(path, mode)
    const char *path, *mode;
{
    FILE *r = fopen(path, mode);
    if (NULL == r) {
	fprintf(stderr, "%s: unable to open file %s:",
		pr_program_name, path);
	perror("");
	exit (-1);
    }
    return r;
}

/* End of fail-stop wrappers. */
/* =========================================================== */
