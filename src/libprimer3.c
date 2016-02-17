/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008,2009,
              2010,2011,2012,2016
Whitehead Institute for Biomedical Research, Steve Rozen
(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky.
All rights reserved.

    This file is part of the primer3 suite and libraries.

    The primer3 suite and libraries are free software;
    you can redistribute them and/or modify them under the terms
    of the GNU General Public License as published by the Free
    Software Foundation; either version 2 of the License, or (at
    your option) any later version.

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this software (file gpl-2.0.txt in the source
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

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <signal.h>
#include <float.h>
#include <string.h>
#include <ctype.h> /* toupper */

#ifdef __GNUC__
#include <ext/hash_map>
#else
#include <hash_map>
#endif

namespace std
{
  using namespace __gnu_cxx;
}

#include "dpal.h"
#include "thal.h"
#include "oligotm.h"
#include "libprimer3.h"

/* #define's */

/*
 * OPTIMIZE_OK_REGIONS 1 allows _optimize_ok_regions_list() to use the
 * max/min product size info and the max/min oligo length to reduce
 * the sizes of the ok regions (while still generating the same primer
 * pairs in the same order).  Set OPTIMIZE_OK_REGIONS to 0 confirm
 * that the results do not change.  (The output tags
 * PRIMER_{LEFT,RIGHT,PAIR}_EXPLAIN _will_ likely change.)
 */
#define OPTIMIZE_OK_REGIONS 1

#ifndef MAX_PRIMER_LENGTH
#define MAX_PRIMER_LENGTH 36
#endif
#if (MAX_PRIMER_LENGTH > DPAL_MAX_ALIGN)
#error "MAX_PRIMER_LENGTH must be <= DPAL_MAX_ALIGN"
#endif
#if (MAX_PRIMER_LENGTH > THAL_MAX_ALIGN)
# error "MAX_PRIMER_LENGTH must be <= THAL_MAX_ALIGN"
#endif
#define MAX_NN_TM_LENGTH 36 /* The maxium length for which to use the
                               nearest neighbor model when calculating
                               oligo Tms. */

#define MACRO_CAT_2(A,B) A##B
#define MACRO_VALUE_AS_STRING(A) MACRO_STRING(A)

#define PR_POSITION_PENALTY_IS_NULL(PA) \
(PR_DEFAULT_INSIDE_PENALTY == (PA)->inside_penalty \
 && PR_DEFAULT_OUTSIDE_PENALTY == (PA)->outside_penalty)

#define INITIAL_LIST_LEN     2000 /* Initial size of oligo lists. */
#define INITIAL_NUM_RETURN   5    /* Initial space to allocate for pairs to
                                     return. */

#define PAIR_OK 1
#define PAIR_FAILED 0

#define OK_OR_MUST_USE(H) (!p3_ol_has_any_problem(H) || (H)->must_use)

#define PR_UNDEFINED_INT_OPT          INT_MIN
#define PR_UNDEFINED_DBL_OPT          DBL_MIN

/* Undefined value for alignment score (meaning do not check) used for
   maximum template mispriming or mishyb. */
#define PR_UNDEFINED_ALIGN_OPT        -100.0

#define TRIMMED_SEQ_LEN(X) ((X)->incl_l)

typedef struct dpal_arg_holder {
  dpal_args *local;
  dpal_args *end;
  dpal_args *local_end;
  dpal_args *local_ambig;
  dpal_args *local_end_ambig;
} dpal_arg_holder;
typedef struct thal_arg_holder {
  thal_args *any;
  thal_args *end1;
  thal_args *end2;
  thal_args *hairpin_th;
} thal_arg_holder;

static jmp_buf _jmp_buf;

static int thermodynamic_alignment_length_error = 0;
static char *thermodynamic_alignment_length_error_msg = NULL;

/* Variables needed in choose_pair_or_triple, need to be global
   for memory de-allocation reasons. Full description given in the function. */
static int *max_j_seen;
static std::hash_map<int, primer_pair*> **pairs;

/* Function declarations. */

static void  pr_set_default_global_args_1(p3_global_settings *);
static void  pr_set_default_global_args_2(p3_global_settings *);
static void _adjust_seq_args(const p3_global_settings *pa,
                             seq_args *sa,
                             pr_append_str *nonfatal_err,
                             pr_append_str *warning);

static void _optimize_ok_regions_list(const p3_global_settings *pa,
				      seq_args *sa);

static int any_5_prime_ol_extension_has_problem(const primer_rec *);

static int p3_ol_is_uninitialized(const primer_rec *);

static int fake_a_sequence(seq_args *sa, const p3_global_settings *pa);

static int    _pr_data_control(const p3_global_settings *,
                               const seq_args *,
                               pr_append_str *glob_err,
                               pr_append_str *nonfatal_err,
                               pr_append_str *warning);

static int    _pr_need_pair_template_mispriming(const p3_global_settings *pa);
static int    _pr_need_pair_template_mispriming_thermod(const p3_global_settings *pa);

static int    _pr_need_template_mispriming(const p3_global_settings *);
static int    _pr_need_template_mispriming_thermod(const p3_global_settings *);

static void   _pr_substr(const char *, int, int, char *);


static int    _check_and_adjust_intervals(seq_args *sa,
                                          int seq_len,
                                          int first_index,
                                          pr_append_str * nonfatal_err,
                                          pr_append_str *warning);

static int    _check_and_adjust_overlap_pos(seq_args *sa,
                                            int *list,
                                            int *count,
                                            const char *tag,
                                            int seq_len,
                                            int first_index,
                                            pr_append_str * nonfatal_err,
                                            pr_append_str *warning);

static int    _check_and_adjust_1_interval(const char *,
                                           int num,
                                           interval_array_t,
                                           int,
                                           int first_index,
                                           pr_append_str *err,
                                           seq_args *,
                                           pr_append_str
                                           *warning, int);

static void   sort_primer_array(oligo_array *);

static void   add_pair(const primer_pair *, pair_array_t *);

static double align(const char *, const char*, const dpal_args *a);

static double align_thermod(const char *, const char *, const thal_args *a);

static int    characterize_pair(p3retval *p,
                                const p3_global_settings *,
                                const seq_args *,
                                int, int, int,
                                primer_pair *,
                                const dpal_arg_holder*,
                                const thal_arg_holder*,
                                int update_stats);

static void    choose_pair_or_triple(p3retval *,
                                    const p3_global_settings *,
                                    const seq_args *,
                                    const dpal_arg_holder *,
                                    const thal_arg_holder *,
				    const thal_arg_holder *,
                                    pair_array_t *);

static int    sequence_quality_is_ok(const p3_global_settings *, primer_rec *,
                                     oligo_type,
                                     const seq_args *, int, int,
                                     oligo_stats *global_oligo_stats,
                                     const args_for_one_oligo_or_primer *);

static int    choose_internal_oligo(p3retval *,
                                    const primer_rec *, const primer_rec *,
                                    int *,
                                    const seq_args *,
                                    const p3_global_settings *,
                                    const dpal_arg_holder *,
                                    const thal_arg_holder *);

void          compute_position_penalty(const p3_global_settings *,
                                       const seq_args *,
                                       primer_rec *, oligo_type);

static p3retval *create_p3retval(void);

static char   dna_to_upper(char *, int);

static int    find_stop_codon(const char *, int, int);

static void   gc_and_n_content(int, int, const char *, primer_rec *);

static int    make_detection_primer_lists(p3retval *,
                                const p3_global_settings *,
                                const seq_args *,
                                const dpal_arg_holder *,
                                const thal_arg_holder *);

static int    make_complete_primer_lists(p3retval *retval,
                                const p3_global_settings *pa,
                                const seq_args *sa,
                                const dpal_arg_holder *dpal_arg_to_use,
                                const thal_arg_holder *thal_arg_to_use,
				const thal_arg_holder *thal_oligo_arg_to_use);

static int    add_primers_to_check(p3retval *retval,
                                const p3_global_settings *pa,
                                const seq_args *sa,
                                const dpal_arg_holder *dpal_arg_to_use,
                                const thal_arg_holder *thal_arg_to_use,
				const thal_arg_holder *thal_oligo_arg_to_use);

static int    pick_sequencing_primer_list(p3retval *retval,
                                const p3_global_settings *pa,
                                const seq_args *sa,
                                const dpal_arg_holder *dpal_arg_to_use,
                                const thal_arg_holder *thal_arg_to_use);

static int    make_internal_oligo_list(p3retval *,
                                       const p3_global_settings *,
                                       const seq_args *,
                                       const dpal_arg_holder *,
                                       const thal_arg_holder *);

static int    pick_only_best_primer(const int, const int,
                                    oligo_array *,
                                    const p3_global_settings *,
                                    const seq_args *,
                                    const dpal_arg_holder *,
                                    const thal_arg_holder *,
                                    p3retval *);

static int    pick_primer_range(const int, const int, int *,
                                oligo_array *,
                                const p3_global_settings *,
                                const seq_args *,
                                const dpal_arg_holder *,
                                const thal_arg_holder *,
                                p3retval *retval);

static int    add_one_primer(const char *, int *, oligo_array *,
                             const p3_global_settings *,
                             const seq_args *,
                             const dpal_arg_holder *,
                             const thal_arg_holder *,
                             p3retval *);

static int   add_one_primer_by_position(int, int, int *,
                                        oligo_array *,
                                        const p3_global_settings *,
                                        const seq_args *,
                                        const dpal_arg_holder *,
                                        const thal_arg_holder *,
                                        p3retval *);

static int   pick_primers_by_position(const int, const int,
                                      int *,
                                      oligo_array *,
                                      const p3_global_settings *,
                                      const seq_args *,
                                      const dpal_arg_holder *,
                                      const thal_arg_holder *,
                                      p3retval *);

static double obj_fn(const p3_global_settings *, primer_pair *);

static int    left_oligo_in_pair_overlaps_used_oligo(const primer_rec *left,
                                                     const primer_pair *best_pair,
                                                     int min_dist);

static int    right_oligo_in_pair_overlaps_used_oligo(const primer_rec *right,
                                                      const primer_pair *best_pair,
                                                      int min_dist);

static int    oligo_overlaps_interval(int, int,
                                      const interval_array_t, int);

static void   calc_and_check_oligo_features(const p3_global_settings *pa,
                                       primer_rec *,
                                       oligo_type,
                                       const dpal_arg_holder*,
                                       const thal_arg_holder*,
                                       const seq_args *, oligo_stats *,
                                       p3retval *,
                                       const char *);

static void   pr_append(pr_append_str *, const char *);

static void   pr_append_new_chunk(pr_append_str *x, const char *s);

static int    pair_spans_target(const primer_pair *, const seq_args *);
static void   pr_append_w_sep(pr_append_str *, const char *, const char *);
static void*  pr_safe_malloc(size_t x);
static void*  pr_safe_realloc(void *p, size_t x);

static int    compare_primer_pair(const void *, const void*);

static int    primer_rec_comp(const void *, const void *);
static int    print_list_header(FILE *, oligo_type, int, int, int);
static int    print_oligo(FILE *, const seq_args *, int, const primer_rec *,
                          oligo_type, int, int,int);
static char   *strstr_nocase(char *, char *);

static double p_obj_fn(const p3_global_settings *, primer_rec *, int );

static void   oligo_compl(primer_rec *,
                          const args_for_one_oligo_or_primer *po_args,
                          oligo_stats *,
                          const dpal_arg_holder *,
                          const char *oligo_seq,
                          const char *revc_oligo_seq
                          );
static void   oligo_compl_thermod(primer_rec *, /* TO DO 2012-06-01 -- update by removing last argument. */
                          const args_for_one_oligo_or_primer *po_args,
                          oligo_stats *,
                          const thal_arg_holder *,
                          const char *oligo_seq,
                          const char *revc_oligo_seq
                          );
static void   oligo_hairpin(primer_rec *,
                            const args_for_one_oligo_or_primer *po_args,
                            oligo_stats *,
                            const thal_arg_holder *,
                            const char *oligo_seq
                            );

static void   oligo_compute_sequence_and_reverse(primer_rec *,
                                                 const seq_args *,
                                                 oligo_type,  
                                                 int*, int*,
                                                 char*, char*);

static void   oligo_repeat_library_mispriming(primer_rec *,
                                              const p3_global_settings *,
                                              const seq_args *,
                                              oligo_type,
                                              oligo_stats *,
                                              const dpal_arg_holder *);

static void   oligo_template_mispriming(primer_rec *,
                                        const p3_global_settings *,
                                        const seq_args *,
                                        oligo_type,
                                        oligo_stats *,
                                        const dpal_args *,
                                        const thal_args *);

static int    pair_repeat_sim(primer_pair *, const p3_global_settings *);

static void   free_repeat_sim_score(p3retval *);
static void   free_primer_repeat_sim_score(primer_rec *); 

/* edited by T. Koressaar for lowercase masking:  */
static int    is_lowercase_masked(int position,
                                   const char *sequence,
                                   primer_rec *h,
                                   oligo_stats *);

static int    primer_must_match(const p3_global_settings *pa,
				primer_rec *h,
				oligo_stats *stats,
				const char *input_oligo_seq,
                                char *match_three_prime,
                                char *match_five_prime);
static int    compare_nucleotides(const char a, const char b);
static int    test_must_match_parameters(char *test);

static void set_retval_both_stop_codons(const seq_args *sa, p3retval *retval);

/* Functions to set bitfield parameters for oligos (or primers) */
static void bf_set_overlaps_target(primer_rec *, int);
static int bf_get_overlaps_target(const primer_rec *);
static void bf_set_overlaps_excl_region(primer_rec *, int);
static int bf_get_overlaps_excl_region(const primer_rec *);
static void bf_set_infinite_pos_penalty(primer_rec *, int);
static int bf_get_infinite_pos_penalty(const primer_rec *);

/* Functions to record problems with oligos (or primers) */
static void initialize_op(primer_rec *);
static void op_set_completely_written(primer_rec *);
static void op_set_must_match_err(primer_rec *);
static void op_set_too_many_ns(primer_rec *);
static void op_set_overlaps_target(primer_rec *);
static void op_set_high_gc_content(primer_rec *);
static void op_set_low_gc_content(primer_rec *);
static void op_set_high_tm(primer_rec *);
static void op_set_low_tm(primer_rec *);
static void op_set_overlaps_excluded_region(primer_rec *);
static void op_set_not_in_any_ok_region(primer_rec *);
static void op_set_high_self_any(primer_rec *oligo);
static void op_set_high_self_end(primer_rec *oligo);
static void op_set_high_hairpin_th(primer_rec *oligo);
static void op_set_no_gc_glamp(primer_rec *);
static void op_set_too_many_gc_at_end(primer_rec *);
static void op_set_high_end_stability(primer_rec *);
static void op_set_high_poly_x(primer_rec *);
static void op_set_low_sequence_quality(primer_rec *);
static void op_set_low_end_sequence_quality(primer_rec *oligo);
static void op_set_high_similarity_to_non_template_seq(primer_rec *);
static void op_set_high_similarity_to_multiple_template_sites(primer_rec *);
static void op_set_overlaps_masked_sequence(primer_rec *);
static void op_set_too_long(primer_rec *);
static void op_set_too_short(primer_rec *);
static void op_set_does_not_amplify_orf(primer_rec *);
/* End functions to set problems in oligos */

/* Global static variables. */
static const char *primer3_copyright_char_star = "\n"
"Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008\n"
"  2009,2010,2011,2012\n"
"Whitehead Institute for Biomedical Research, Steve Rozen\n"
"(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky\n"
"All rights reserved.\n"
"\n"
"    This file is part of the primer3 suite and libraries.\n"
"\n"
"    The primer3 suite and libraries are free software;\n"
"    you can redistribute them and/or modify them under the terms\n"
"    of the GNU General Public License as published by the Free\n"
"    Software Foundation; either version 2 of the License, or (at\n"
"    your option) any later version.\n"
"\n"
"    This software is distributed in the hope that it will be useful,\n"
"    but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"    GNU General Public License for more details.\n"
"\n"
"    You should have received a copy of the GNU General Public License\n"
"    along with this software (file gpl-2.0.txt in the source\n"
"    distribution); if not, write to the Free Software\n"
"    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA\n"
"\n"
"THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n"
"\"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n"
"LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR\n"
"A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT\n"
"OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,\n"
"SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT\n"
"LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,\n"
"DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY\n"
"THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT\n"
"(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE\n"
"OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n"
"\n";

static const char *pr_program_name = "probably primer3_core";

static const int use_end_for_th_template_mispriming = 1;

#define DEFAULT_OPT_GC_PERCENT PR_UNDEFINED_INT_OPT

/* Set the program name */
void
p3_set_program_name(const char *pname) 
{
  pr_program_name = pname;
}

/* ============================================================ */
/* BEGIN functions for global settings                          */
/* ============================================================ */

/* Allocate space for global settings and fill in defaults */
p3_global_settings *
p3_create_global_settings() 
{
  p3_global_settings *r;

  if (!(r = (p3_global_settings *) malloc(sizeof(*r)))) {
    return NULL;
  }

  pr_set_default_global_args_2(r);

  return r;
}

/* Allocate space for global settings and fill in defaults */
p3_global_settings *
p3_create_global_settings_default_version_1() 
{
  p3_global_settings *r;

  if (!(r = (p3_global_settings *) malloc(sizeof(*r)))) {
    return NULL;
  }

  pr_set_default_global_args_1(r);

  return r;
}

/* Free the space of global settings */
void
p3_destroy_global_settings(p3_global_settings *a) 
{
  if (NULL != a) {
    if (NULL != a->p_args.must_match_five_prime) {
    	free(a->p_args.must_match_five_prime);
    }
    if (NULL != a->p_args.must_match_three_prime) {
    	free(a->p_args.must_match_three_prime);
    }
    if (NULL != a->o_args.must_match_five_prime) {
    	free(a->o_args.must_match_five_prime);
    }
    if (NULL != a->o_args.must_match_three_prime) {
    	free(a->o_args.must_match_three_prime);
    }
    destroy_seq_lib(a->p_args.repeat_lib);
    destroy_seq_lib(a->o_args.repeat_lib);
    free(a);
  }
}

static void
pr_set_default_global_args_2(p3_global_settings *a) 
/* Write the default values for default_values=2 into
   the p3_global_settings struct */
{
  pr_set_default_global_args_1(a);
  a->tm_santalucia                    = santalucia_auto;
  a->salt_corrections                 = santalucia;
  a->thermodynamic_oligo_alignment    = 1;
  a->thermodynamic_template_alignment = 0;
  a->p_args.divalent_conc             = 1.5;
  a->p_args.dntp_conc                 = 0.6;
  a->lib_ambiguity_codes_consensus    = 0;
}

/* Write the default values for default_values=1 into
   the p3_global_settings struct */
static void
pr_set_default_global_args_1(p3_global_settings *a) 
{
  memset(a, 0, sizeof(*a));

  /* Arguments for primers ================================= */
  a->p_args.opt_size          = 20;
  a->p_args.min_size          = 18;
  a->p_args.max_size          = 27;

  a->p_args.opt_tm            = 60;
  a->p_args.min_tm            = 57;
  a->p_args.max_tm            = 63;

  a->p_args.min_gc            = 20.0;
  a->p_args.opt_gc_content    = DEFAULT_OPT_GC_PERCENT;
  a->p_args.max_gc            = 80.0;
  a->p_args.salt_conc         = 50.0;
  a->p_args.divalent_conc     = 0.0;
  a->p_args.dntp_conc         = 0.0;

  a->p_args.dna_conc          = 50.0;
  a->p_args.num_ns_accepted   = 0;
  a->p_args.max_self_any      = 8.0;
  a->p_args.max_self_end      = 3.0;
  a->p_args.max_self_any_th   = 47.0;
  a->p_args.max_self_end_th   = 47.0;
  a->p_args.max_hairpin_th    = 47.0;
  a->p_args.max_poly_x        = 5;
  a->p_args.max_repeat_compl  = 12.0;
  a->p_args.min_quality       = 0;
  a->p_args.min_end_quality   = 0;
  a->p_args.max_template_mispriming = PR_UNDEFINED_ALIGN_OPT;
  a->p_args.max_template_mispriming_th = PR_UNDEFINED_ALIGN_OPT;
  /* The following apply only to primers (and not to internal
     oligos). */
  a->gc_clamp                 = 0;
  a->max_end_gc               = 5;

/* Weights for objective functions for oligos and pairs. */
  a->p_args.weights.compl_any     = 0;
  a->p_args.weights.compl_any_th  = 0;
  a->p_args.weights.compl_end     = 0;
  a->p_args.weights.compl_end_th  = 0;
  a->p_args.weights.end_quality   = 0;
  a->p_args.weights.end_stability = 0;
  a->p_args.weights.gc_content_gt = 0;
  a->p_args.weights.gc_content_lt = 0;
  a->p_args.weights.hairpin_th    = 0;
  a->p_args.weights.length_gt     = 1;
  a->p_args.weights.length_lt     = 1;
  a->p_args.weights.num_ns        = 0;
  a->p_args.weights.pos_penalty   = 1;
  a->p_args.weights.repeat_sim    = 0;
  a->p_args.weights.seq_quality   = 0;
  a->p_args.weights.temp_cutoff   = 5;
  a->p_args.weights.temp_gt       = 1;
  a->p_args.weights.temp_lt       = 1;
  a->p_args.weights.template_mispriming = 0.0;
  a->p_args.weights.template_mispriming_th = 0.0;
  a->p_args.must_match_five_prime  = NULL;
  a->p_args.must_match_three_prime = NULL;
  /* End of weights for objective functions for oligos and pairs. */

  /* End of arguments for primers =========================== */

  a->max_diff_tm         = 100.0;
  a->tm_santalucia       = breslauer_auto;
  a->salt_corrections    = schildkraut;
  a->pair_compl_any      = 8.0;
  a->pair_compl_end      = 3.0;
  a->pair_compl_any_th   = 47.0;
  a->pair_compl_end_th   = 47.0;
  a->thermodynamic_oligo_alignment = 0;
  a->thermodynamic_template_alignment = 0;
  a->liberal_base        = 0;
  a->primer_task         = generic;
  a->pick_left_primer    = 1;
  a->pick_right_primer   = 1;
  a->pick_internal_oligo = 0;
  a->first_base_index    = 0;
  a->num_return          = 5;
  a->pr_min[0]           = 100;
  a->pr_max[0]           = 300;
  a->num_intervals       = 1;
  a->pair_repeat_compl   = 24.0;
  a->quality_range_min   = 0;
  a->quality_range_max   = 100;
  a->outside_penalty     = PR_DEFAULT_OUTSIDE_PENALTY;
  a->inside_penalty      = PR_DEFAULT_INSIDE_PENALTY;
  a->max_end_stability   = 100.0;
  a->lowercase_masking   = 0;
  a->product_max_tm      = PR_DEFAULT_PRODUCT_MAX_TM;
  a->product_min_tm      = PR_DEFAULT_PRODUCT_MIN_TM;
  a->product_opt_tm      = PR_UNDEFINED_DBL_OPT;
  a->product_opt_size    = PR_UNDEFINED_INT_OPT;
  a->pair_max_template_mispriming = PR_UNDEFINED_ALIGN_OPT;
  a->pair_max_template_mispriming_th = PR_UNDEFINED_ALIGN_OPT;
  a->o_args.opt_size        = 20;
  a->o_args.min_size        = 18;
  a->o_args.max_size        = 27;
  a->o_args.opt_tm          = 60.0;
  a->o_args.min_tm          = 57.0;
  a->o_args.max_tm          = 63.0;
  a->o_args.min_gc          = 20.0;
  a->o_args.max_gc          = 80.0;
  a->o_args.opt_gc_content  = DEFAULT_OPT_GC_PERCENT;
  a->o_args.max_poly_x      = 5;
  a->o_args.salt_conc       = 50.0;
  a->o_args.divalent_conc   = 0.0;
  a->o_args.dntp_conc       = 0.0;
  a->o_args.dna_conc        = 50.0;
  a->o_args.num_ns_accepted = 0;
  a->o_args.max_self_any    = 12.0;
  a->o_args.max_self_end    = 12.0;
  a->o_args.max_self_any_th = 47.0;
  a->o_args.max_self_end_th = 47.0;
  a->o_args.max_hairpin_th  = 47.0;
  a->o_args.max_repeat_compl= 12.0;

  a->o_args.min_quality           = 0;
  a->o_args.min_end_quality       = 0;
  a->o_args.max_template_mispriming = PR_UNDEFINED_ALIGN_OPT;
  a->o_args.max_template_mispriming_th = PR_UNDEFINED_ALIGN_OPT;
  a->o_args.weights.temp_gt       = 1;
  a->o_args.weights.temp_lt       = 1;
  a->o_args.weights.length_gt     = 1;
  a->o_args.weights.length_lt     = 1;
  a->o_args.weights.gc_content_gt = 0;
  a->o_args.weights.gc_content_lt = 0;
  a->o_args.weights.compl_any     = 0;
  a->o_args.weights.compl_end     = 0;
  a->o_args.weights.compl_any_th  = 0;
  a->o_args.weights.compl_end_th  = 0;
  a->o_args.weights.hairpin_th    = 0;
  a->o_args.weights.num_ns        = 0;
  a->o_args.weights.repeat_sim    = 0;
  a->o_args.weights.seq_quality   = 0;
  a->o_args.weights.end_quality   = 0;
  a->o_args.must_match_five_prime  = NULL;
  a->o_args.must_match_three_prime = NULL;

  a->pr_pair_weights.primer_quality  = 1;
  a->pr_pair_weights.io_quality      = 0;
  a->pr_pair_weights.diff_tm         = 0;
  a->pr_pair_weights.compl_any       = 0;
  a->pr_pair_weights.compl_end       = 0;
  a->pr_pair_weights.compl_any_th    = 0;
  a->pr_pair_weights.compl_end_th    = 0;
  a->pr_pair_weights.temp_cutoff     = 5;
  a->pr_pair_weights.repeat_sim      = 0;

  a->pr_pair_weights.product_tm_lt   = 0;
  a->pr_pair_weights.product_tm_gt   = 0;
  a->pr_pair_weights.product_size_lt = 0;
  a->pr_pair_weights.product_size_gt = 0;

  a->lib_ambiguity_codes_consensus   = 1;
  /*  Set to 1 for backward compatibility. This _NOT_ what
      one normally wants, since many libraries contain
      strings of N, which then match every oligo (very bad).
  */

  a->min_left_three_prime_distance   = -1;
  a->min_right_three_prime_distance  = -1;

  a->sequencing.lead                 = 50;
  a->sequencing.spacing              = 500;
  a->sequencing.interval             = 250;
  a->sequencing.accuracy             = 20;
  
  a->min_5_prime_overlap_of_junction = 7;
  a->min_3_prime_overlap_of_junction = 4;
    
}

/* Add a pair of integers to an  array of intervals */
int
p3_add_to_interval_array(interval_array_t2 *interval_arr, int i1, int i2)
{
  int c = interval_arr->count;
  if (c >= PR_MAX_INTERVAL_ARRAY) return 1;
  interval_arr->pairs[c][0] = i1;
  interval_arr->pairs[c][1] = i2;
  interval_arr->count++;
  return  0;
}

/* Add 2 pairs of integers to an array of 2 intervals */
int
p3_add_to_2_interval_array(interval_array_t4 *interval_arr, int i1, int i2, int i3, int i4)
{
  int c = interval_arr->count;
  if (c >= PR_MAX_INTERVAL_ARRAY) return 1;
  /* for a region either both values are given, or none is given */
  if (((i1 == -1) && (i2 != -1)) || ((i1 != -1) && (i2 == -1)))
    return 2;
  if (((i3 == -1) && (i4 != -1)) || ((i3 != -1) && (i4 == -1)))
    return 2;
  interval_arr->left_pairs[c][0] = i1;
  interval_arr->left_pairs[c][1] = i2;
  interval_arr->right_pairs[c][0] = i3;
  interval_arr->right_pairs[c][1] = i4;
  if ((i1 == -1) && (i2 == -1))
    interval_arr->any_left = 1;
  if ((i3 == -1) && (i4 == -1))
    interval_arr->any_right = 1;
  if ((i1 == -1) && (i2 == -1) && (i3 == -1) && (i4 == -1))
    interval_arr->any_pair = 1;
  interval_arr->count++;
  return  0;
}

/* ============================================================ */
/* END functions for global settings                            */
/* ============================================================ */

int
interval_array_t2_count(const interval_array_t2 *array) 
{
  return array->count;
}

const int *
interval_array_t2_get_pair(const interval_array_t2 *array, int i) 
{
  if (i > array->count) abort();
  if (i < 0) abort();
  return array->pairs[i];
}

/* ============================================================ */
/* BEGIN functions for p3retval                                 */
/* ============================================================ */

/* Allocate a new primer3 state. Return NULL if out of memory. Assuming
   malloc sets errno to ENOMEM according to Unix98, set errno to ENOMEM
   on out-of-memory error. */
static p3retval *
create_p3retval(void) 
{
  p3retval *state = (p3retval *)malloc(sizeof(*state));
  if (!state)
    return NULL;

  state->fwd.oligo
    = (primer_rec *) malloc(sizeof(*state->fwd.oligo)  * INITIAL_LIST_LEN);

  state->rev.oligo
    = (primer_rec *) malloc(sizeof(*state->rev.oligo)  * INITIAL_LIST_LEN);

  state->intl.oligo
    = (primer_rec *) malloc(sizeof(*state->intl.oligo) * INITIAL_LIST_LEN);

  if (state->fwd.oligo == NULL
      || state->rev.oligo == NULL
      || state->intl.oligo == NULL) {
    free(state);
    return NULL;
  }

  state->fwd.storage_size  = INITIAL_LIST_LEN;
  state->rev.storage_size  = INITIAL_LIST_LEN;
  state->intl.storage_size = INITIAL_LIST_LEN;

  state->fwd.num_elem  = 0;
  state->rev.num_elem  = 0;
  state->intl.num_elem = 0;

  state->fwd.type  = OT_LEFT;
  state->intl.type = OT_INTL;
  state->rev.type  = OT_RIGHT;

  state->best_pairs.storage_size = 0;
  state->best_pairs.pairs = NULL;
  state->best_pairs.num_pairs = 0;

  init_pr_append_str(&state->glob_err);
  init_pr_append_str(&state->per_sequence_err);
  init_pr_append_str(&state->warnings);

  memset(&state->fwd.expl,  0, sizeof(state->fwd.expl));
  memset(&state->rev.expl,  0, sizeof(state->rev.expl));
  memset(&state->intl.expl, 0, sizeof(state->intl.expl));

  memset(&state->best_pairs.expl, 0, sizeof(state->best_pairs.expl));

  return state;
}
/* ============================================================ */
/* BEGIN functions for dpal_arg_holder                          */
/* ============================================================ */

/* Create the dpal arg holder */
dpal_arg_holder *
create_dpal_arg_holder () 
{

  dpal_arg_holder *h
    = (dpal_arg_holder *) pr_safe_malloc(sizeof(dpal_arg_holder));

  h->local = (dpal_args *) pr_safe_malloc(sizeof(*h->local));
  set_dpal_args(h->local);
  h->local->flag = DPAL_LOCAL;

  h->end = (dpal_args *) pr_safe_malloc(sizeof(*h->end));
  set_dpal_args(h->end);
  h->end->flag = DPAL_GLOBAL_END;

  h->local_end = (dpal_args *) pr_safe_malloc(sizeof(*h->local_end));
  set_dpal_args(h->local_end);
  h->local_end->flag = DPAL_LOCAL_END;

  h->local_ambig  = (dpal_args *) pr_safe_malloc(sizeof(*h->local_ambig));
  *h->local_ambig = *h->local;
  PR_ASSERT(dpal_set_ambiguity_code_matrix(h->local_ambig));

  h->local_end_ambig = (dpal_args *) pr_safe_malloc(sizeof(*h->local_end_ambig));
  *h->local_end_ambig = *h->local_end;
  PR_ASSERT(dpal_set_ambiguity_code_matrix(h->local_end_ambig));

  return h;
}

/* Free the dpal arg holder */
void
destroy_dpal_arg_holder(dpal_arg_holder *h) 
{
  if (NULL != h) {
    free(h->local);
    free(h->end);
    free(h->local_end);
    free(h->local_ambig);
    free(h->local_end_ambig);
    free(h);
  }
}

/* This is a static variable that is initialized once
   in choose_primers().  We make this variable have
   'file' scope, we do not have to remember to free the associated
   storage after each call to choose_primers(). */
static dpal_arg_holder *dpal_arg_to_use = NULL;

/* ============================================================ */
/* END functions for dpal_arg_holder                            */
/* ============================================================ */

/* ============================================================ */
/* BEGIN functions for thal_arg_holder                           */
/* ============================================================ */
/* Create the thal arg holder */
thal_arg_holder *
create_thal_arg_holder (const args_for_one_oligo_or_primer *po_args)
{
   
   thal_arg_holder *h
     = (thal_arg_holder *) pr_safe_malloc(sizeof(thal_arg_holder));
   
   h->any = (thal_args *) pr_safe_malloc(sizeof(*h->any));
   set_thal_default_args(h->any);
   h->any->type = thal_any;
   h->any->mv = po_args->salt_conc;
   h->any->dv = po_args->divalent_conc;
   h->any->dntp = po_args->dntp_conc;
   h->any->dna_conc = po_args->dna_conc;
   
   h->end1 = (thal_args *) pr_safe_malloc(sizeof(*h->end1));
   set_thal_default_args(h->end1);
   h->end1->type = thal_end1;
   h->end1->mv = po_args->salt_conc;
   h->end1->dv = po_args->divalent_conc;
   h->end1->dntp = po_args->dntp_conc;
   h->end1->dna_conc = po_args->dna_conc;
   
   h->end2 = (thal_args *) pr_safe_malloc(sizeof(*h->end2));
   set_thal_default_args(h->end2);
   h->end2->type = thal_end2;
   h->end2->mv = po_args->salt_conc;
   h->end2->dv = po_args->divalent_conc;
   h->end2->dntp = po_args->dntp_conc;
   h->end2->dna_conc = po_args->dna_conc;
   
   h->hairpin_th  = (thal_args *) pr_safe_malloc(sizeof(*h->hairpin_th));
   set_thal_default_args(h->hairpin_th); 
   h->hairpin_th->type = thal_hairpin;
   h->hairpin_th->mv = po_args->salt_conc;
   h->hairpin_th->dv = po_args->divalent_conc;
   h->hairpin_th->dntp = po_args->dntp_conc;
   h->hairpin_th->dna_conc = po_args->dna_conc;
   h->hairpin_th->dimer = 0;
   
   return h;
}

/* Free the thal arg holder */
void
destroy_thal_arg_holder(thal_arg_holder *h) 
{  
  if (NULL != h) {
    free(h->any);
    free(h->end1);
    free(h->end2);
    free(h->hairpin_th);
    free(h);
  }
}

static thal_arg_holder *thal_arg_to_use = NULL;
static thal_arg_holder *thal_oligo_arg_to_use = NULL;
/* ============================================================ */
/* END functions for thal_arg_holder                            */
/* ============================================================ */

/* Deallocate a primer3 state */
void
destroy_p3retval(p3retval *state)
{
  if (!state)
    return;

  free_repeat_sim_score(state);
  free(state->fwd.oligo);
  free(state->rev.oligo);
  free(state->intl.oligo);
  if (state->best_pairs.storage_size != 0 && state->best_pairs.pairs)
    free(state->best_pairs.pairs);

  destroy_pr_append_str_data(&state->glob_err);
  destroy_pr_append_str_data(&state->per_sequence_err);
  destroy_pr_append_str_data(&state->warnings);

  free(state);
}

void 
destroy_dpal_thal_arg_holder() {
   if(dpal_arg_to_use)
     destroy_dpal_arg_holder(dpal_arg_to_use);
   if(thal_arg_to_use)
     destroy_thal_arg_holder(thal_arg_to_use);
   if(thal_oligo_arg_to_use)
     destroy_thal_arg_holder(thal_oligo_arg_to_use);
}

const oligo_array *
p3_get_rv_fwd(const p3retval *r) {
  return &r->fwd;
}

const oligo_array *
p3_get_rv_intl(const p3retval *r) {
  return &r->intl;
}

const oligo_array *
p3_get_rv_rev(const p3retval *r) {
  return &r->rev;
}

const pair_array_t *
p3_get_rv_best_pairs(const p3retval *r) {
  return &r->best_pairs;
}

/* ============================================================ */
/* END functions for p3retval                                   */
/* ============================================================ */

/* ============================================================ */
/* BEGIN functions for seq_arg_holder                           */
/* ============================================================ */

/* Create and initialize a seq_args data structure */
seq_args *
create_seq_arg() 
{
  seq_args *r = (seq_args *) malloc(sizeof(*r));
  if (NULL == r) return NULL; /* Out of memory */
  memset(r, 0, sizeof(*r));
  r->start_codon_pos = PR_DEFAULT_START_CODON_POS;
  r->incl_l = -1; /* Indicates logical NULL. */

  r->force_left_start = PR_NULL_FORCE_POSITION; /* Indicates logical NULL. */
  r->force_left_end = PR_NULL_FORCE_POSITION; /* Indicates logical NULL. */
  r->force_right_start = PR_NULL_FORCE_POSITION; /* Indicates logical NULL. */
  r->force_right_end = PR_NULL_FORCE_POSITION; /* Indicates logical NULL. */
  r->primer_overlap_junctions_count = 0;

  r->n_quality = 0;
  r->quality = NULL;

  return r;
}

/* Free a seq_arg data structure */
void
destroy_seq_args(seq_args *sa) 
{
  if (NULL == sa) return;
  free(sa->internal_input);
  free(sa->left_input);
  free(sa->right_input);
  free(sa->sequence);
  free(sa->quality);
  free(sa->trimmed_seq);

  /* edited by T. Koressaar for lowercase masking */
  free(sa->trimmed_orig_seq);

  free(sa->upcased_seq);
  free(sa->upcased_seq_r);
  free(sa->sequence_name);
  free(sa);
}

/* ============================================================ */
/* END functions for seq_arg_holder                             */
/* ============================================================ */


/* Function used to de-allocate the memory used by the hash maps
   allocated in choose_pair_or_triple. */
static void free_pair_memory(int rev_num_elem)
{
  std::hash_map<int, primer_pair*> *hmap;
  std::hash_map<int, primer_pair*>::iterator it;
  primer_pair *pp;
  int i;

  free(max_j_seen);
  for (i=0; i<rev_num_elem; i++) {
    hmap = pairs[i];
    if (hmap) {
      for (it=hmap->begin(); it!=hmap->end(); it++) {
        pp = it->second;
        if (pp != NULL)
	  delete pp;
      }
      delete hmap;
    }
  }
  free(pairs);
}


/* ============================================================ */
/* BEGIN choose_primers()                                       */
/* The main primer3 interface                                   */
/* See libprimer3.h for documentation.                          */
/* ============================================================ */
p3retval *
choose_primers(const p3_global_settings *pa,
               seq_args *sa)
{
  /* Create retval and set were to find the results */
  p3retval *retval = create_p3retval();
  if (retval == NULL)  return NULL;

  PR_ASSERT(NULL != pa);
  PR_ASSERT(NULL != sa);

  if (pa->dump) {
    printf("Start of choose_primers:\n");
    p3_print_args(pa, sa) ;
  }

  /* Set the general output type */
  if (pa->pick_left_primer && pa->pick_right_primer) {
    retval->output_type = primer_pairs;
  } else {
    retval->output_type = primer_list;
  }
  if (pa->primer_task == pick_primer_list ||
      pa->primer_task == pick_sequencing_primers) {
    retval->output_type = primer_list;
  }

  /*
   * For catching ENOMEM.  WARNING: We can only use longjmp to escape
   * from errors that have been called through choose_primers().
   * Therefore, if we subsequently update other static functions in
   * this file to have external linkage then we need to check whether
   * they call (or use functions that in turn call) longjmp to handle
   * ENOMEM.
   */
  if (setjmp(_jmp_buf) != 0) {
    /* Check if this was a thermodynamic alignment length error. */
    if (thermodynamic_alignment_length_error == 1) {
      thermodynamic_alignment_length_error = 0;
      /* Set per sequence error */
      pr_append_new_chunk(&retval->per_sequence_err, 
			  thermodynamic_alignment_length_error_msg);
      free(thermodynamic_alignment_length_error_msg);
      thermodynamic_alignment_length_error_msg = NULL;
      /* Other necessary cleanup. */
      free_pair_memory(retval->rev.num_elem);
      return retval;
    }
    /* This was a memory error. */
    destroy_p3retval(retval);
    return NULL;  /* If we get here, that means errno should be ENOMEM. */
  }

  /* Change some parameters to fit the task */
  _adjust_seq_args(pa, sa, &retval->per_sequence_err,&retval->warnings);

  if (pa->dump) {
    printf("After _adjust_seq_args\n");
    p3_print_args(pa, sa) ;
  }

  if (!pr_is_empty(&retval->per_sequence_err)) return retval;

  /* TO DO -- move the check below to pr_data_control and issue a
     warning if it fails. */
  /* if (pa->p_args.min_quality != 0
     && pa->p_args.min_end_quality < pa->p_args.min_quality)
       ... issue warning ...
     pa->p_args.min_end_quality = pa->p_args.min_quality; */

  /* Check if the input in sa and pa makes sense */
  if (_pr_data_control(pa, sa,
                       &retval->glob_err,  /* Fatal errors */
                       &retval->per_sequence_err, /* Non-fatal errors */
                       &retval->warnings
                       ) !=0 ) {
    return retval;
  }


  set_retval_both_stop_codons(sa, retval);

  /* Set the parameters for alignment functions
     if dpal_arg_to_use, a static variable that has 'file'
     scope, has not yet been initialized. */
  if (dpal_arg_to_use == NULL)
    dpal_arg_to_use = create_dpal_arg_holder();
   if(thal_arg_to_use == NULL) {
      thal_arg_to_use = create_thal_arg_holder(&pa->p_args);
   } else {
      destroy_thal_arg_holder(thal_arg_to_use);
      thal_arg_to_use = create_thal_arg_holder(&pa->p_args);
   }
   if(thal_oligo_arg_to_use == NULL) {
      thal_oligo_arg_to_use = create_thal_arg_holder(&pa->o_args);
   } else {
      destroy_thal_arg_holder(thal_oligo_arg_to_use);
      thal_oligo_arg_to_use = create_thal_arg_holder(&pa->o_args);
   } 
  if (pa->primer_task == pick_primer_list) {
    make_complete_primer_lists(retval, pa, sa,
                               dpal_arg_to_use,thal_arg_to_use,thal_oligo_arg_to_use);
  } else if (pa->primer_task == pick_sequencing_primers) {
    pick_sequencing_primer_list(retval, pa, sa,
                                dpal_arg_to_use,thal_arg_to_use);
  } else if (pa->primer_task == check_primers) {
    add_primers_to_check(retval, pa, sa,
                         dpal_arg_to_use, thal_arg_to_use, thal_oligo_arg_to_use);
  } else { /* The general way to pick primers */
    /* Populate the forward and reverse primer lists */
    if (make_detection_primer_lists(retval, pa, sa,
                                    dpal_arg_to_use,thal_arg_to_use) != 0) {
      /* There was an error */
      return retval;
    }
    /* Populate the internal oligo lists */
    if ( pa->pick_internal_oligo) {
      if (make_internal_oligo_list(retval, pa, sa,
                                   dpal_arg_to_use,thal_oligo_arg_to_use) != 0) {
        /* There was an error*/
        return retval;
      }
    }
  }

  if (pa->pick_right_primer &&
      (pa->primer_task != pick_sequencing_primers))
    sort_primer_array(&retval->rev);

  if (pa->pick_left_primer &&
      (pa->primer_task != pick_sequencing_primers))
    sort_primer_array(&retval->fwd);

  /* If we are returning a list of internal oligos, sort them by their
     'goodness'. We do not care if these are sorted if we end up in
     choose_pair_or_triple(), since this calls
     choose_internal_oligo(), which selects the best internal oligo
     for a given primer pair. */
  if (retval->output_type == primer_list && pa->pick_internal_oligo == 1)
    sort_primer_array(&retval->intl);

  /* Select primer pairs if needed */
  if (retval->output_type == primer_pairs) {
    choose_pair_or_triple(retval, pa, sa, dpal_arg_to_use, thal_arg_to_use,
       thal_oligo_arg_to_use, &retval->best_pairs);
  }

  if (pa->dump) {
    printf("End of choose_primers:\n");
    p3_print_args(pa, sa) ;
  }

  return retval;
}
/* ============================================================ */
/* END choose_primers()                                         */
/* ============================================================ */

/* ============================================================ */
/* BEGIN choose_pair_or_triple

   This function uses retval->fwd and retval->rev and
   updates the oligos in these array.

   This function posibly uses retval->intl and updates its
   elements via choose_internal_oligo().

   This function examines primer pairs or triples to find
   pa->num_return pairs or triples to return.

   Results are returned in best_pairs and in
   retval->best_pairs.expl */
/* ============================================================ */
static void
choose_pair_or_triple(p3retval *retval,
                      const p3_global_settings *pa,
                      const seq_args *sa,
                      const dpal_arg_holder *dpal_arg_to_use,
                      const thal_arg_holder *thal_arg_to_use,
		      const thal_arg_holder *thal_oligo_arg_to_use,
                      pair_array_t *best_pairs) {
  int i,j;   /* Loop index. */
  int n_int; /* Index of the internal oligo */
  /*int *max_j_seen; */  /* The maximum value of j (loop index for forward primers)
                            that has been examined for every reverse primer
                            index (i) -- global variable now */
  int update_stats = 1;  /* Flag to indicate whether pair_stats
                            should be updated. */
  primer_pair h;             /* The current pair which is being evaluated. */
  primer_pair the_best_pair; /* The best pair is being "remembered". */
  pair_stats *pair_expl = &retval->best_pairs.expl; /* For statistics */

  int product_size_range_index = 0;
  int trace_me = 0;
  int the_best_i, the_best_j;

  /* Hash maps used to store pairs that were computed */

  /* std::hash_map<int, primer_pair*> **pairs; */
  /* pairs is an array of pointers to hash maps.  It will be indexed
     by the indices of the reverse primers in retval->rev. -- global var now */

  std::hash_map<int, primer_pair*> *hmap, *best_hmap = NULL;
  /* hmap and best_hmap will be pointers to hash maps also pointed to
     by elements of pairs. */

  std::hash_map<int, primer_pair*>::iterator it;
  primer_pair *pp, *best_pp = NULL;
  int pair_found = 0;

  pairs =
    (std::hash_map<int, primer_pair*>**) 
    calloc (retval->rev.num_elem, 
	    sizeof(std::hash_map<int, primer_pair*>*));
  if (!pairs) longjmp(_jmp_buf, 1);
  
  memset(&the_best_pair, 0, sizeof(the_best_pair));
  max_j_seen = (int *) malloc(sizeof(int) * retval->rev.num_elem);
  for (i = 0; i < retval->rev.num_elem; i++) max_j_seen[i] = -1;

  /* Pick pairs till we have enough. */     
  while(1) {

    /* FIX ME, is this memset really needed?  Apparently slow */
    memset(&the_best_pair, 0, sizeof(the_best_pair));

    the_best_i = -1;
    the_best_j = -1;
    /* To start put penalty to the maximum */
    the_best_pair.pair_quality = DBL_MAX;

    /* Iterate over the reverse primers. */
    for (i = 0; i < retval->rev.num_elem; i++) {

      hmap = pairs[i];  
      /* Pairs[i] is NULL if there has never been an assignment to
	 pairs[i] because pairs was allocated by calloc, which
	 sets the allocated memory to 0. */

      /* Only use a primer that *might be* legal or that the caller
         has provided and specified as "must use".  Primers are *NOT*
         FULLY ASSESSED until the call to characterize_pair(), in
         order to avoid expensive computations (mostly alignments)
         unless necessary. */
      if (!OK_OR_MUST_USE(&retval->rev.oligo[i])) {
        /* Can free the memory used by the hmap associated to this
	   reverse primer */
        if (hmap) {
          for (it=hmap->begin(); it!=hmap->end(); it++) {
	    /* it->second is the second element (i.e. the 'value', as
	       opposed to the 'key'). */
            pp = it->second;
            delete pp;
          }
	  if (hmap == best_hmap) best_hmap = NULL;
          delete hmap;
	  hmap = NULL;
          pairs[i] = NULL;
        }
        continue;
      }

      /* If the pair cannot be better than the one already 
       * selected, then we can skip remaining reverse primers */
      if (pa->pr_pair_weights.primer_quality *
           (retval->rev.oligo[i].quality + retval->fwd.oligo[0].quality)
              > the_best_pair.pair_quality) {
        break;
      }

      if (retval->rev.oligo[i].overlaps) {
        /* The stats will not keep track of the pair correctly
           after the first pass, because an oligo might
           have been legal on one pass but become illegal on
           a subsequent pass. */
        if (update_stats) {
          if (trace_me)
            fprintf(stderr,
                    "i=%d, j=%d, overlaps_oligo_in_better_pair++\n",
                    i, j);
          pair_expl->overlaps_oligo_in_better_pair++;
        }
        /* Can free the memory used by the hmap associated to this reverse primer */
        if (hmap) {
          for (it=hmap->begin(); it!=hmap->end(); it++) {
	    /* it->second is the second element (i.e. the 'value', as
	       opposed to the 'key'). */
            pp = it->second;
            delete pp;
          }
	  if (hmap == best_hmap) best_hmap = NULL;
          delete hmap;
	  hmap = NULL;
          pairs[i] = NULL;
        }
        continue;
      }

      /* Loop over forward primers */
      for (j=0; j<retval->fwd.num_elem; j++) {

        /* We check the reverse oligo again, because we may
           have determined that it is "not ok", even though
           (as a far as we knew), it was ok above. */
        if (!OK_OR_MUST_USE(&retval->rev.oligo[i])) {
          /* Can free the memory used by the hmap associated to this reverse primer */
          if (hmap) {
            for (it=hmap->begin(); it!=hmap->end(); it++) {
	      /* it->second is the second element (i.e. the 'value', as
		 opposed to the 'key'). */
              pp = it->second;
              delete pp;
            }
	    if (hmap == best_hmap) best_hmap = NULL;
	    delete hmap;
	    hmap = NULL;
            pairs[i] = NULL;
          }
          break;
        }

        /* Only use a primer that is legal, or that the caller
           has provided and specified as "must use". */
        if (!OK_OR_MUST_USE(&retval->fwd.oligo[j])) continue;

        /* If the pair cannot be better than the one already 
         * selected, then we can skip remaining forward primers 
         * for this reverse primer */
        if (pa->pr_pair_weights.primer_quality *
            (retval->fwd.oligo[j].quality + retval->rev.oligo[i].quality) 
            > the_best_pair.pair_quality) {
          break;
        }

        /* Need to have this here because if we break just above, then,
           at a later iteration, we may need to examine the oligo
           pair with reverse oligo at i and forward oligo at j. */
        update_stats = 0;
        if (j > max_j_seen[i]) {
          if (trace_me)
            fprintf(stderr, 
		    "updates ON: i=%d, j=%d, max_j_seen[%d]=%d\n", 
		    i, j, i, max_j_seen[i]);
          max_j_seen[i] = j;
          if (trace_me)
            fprintf(stderr, "max_j_seen[%d] --> %d\n", i, max_j_seen[i]);
          if (trace_me) fprintf(stderr, "updates on\n");
          update_stats = 1;
        }

        if (retval->fwd.oligo[j].overlaps) {
          /* The stats will not keep track of the pair correctly
             after the first pass, because an oligo might
             have been legal on one pass but become illegal on
             a subsequent pass. */
          if (update_stats) {
            if (trace_me)
              fprintf(stderr,
                      "i=%d, j=%d, overlaps_oligo_in_better_pair++\n",
                      i, j);
            pair_expl->overlaps_oligo_in_better_pair++;
          }
          continue;
        }

	/* Some simple checks first, before searching the hashmap */
	int must_use = 0;
	if ((pa->primer_task == check_primers) || 
	    ((retval->fwd.oligo[j].must_use != 0) &&
	     (retval->rev.oligo[i].must_use != 0))) {
	  must_use = 1;
	}
	
	/* Determine if overlap with an overlap point is required, and
	   if so, whether one of the primers in the pairs overlaps
	   that point. */
	if ((sa->primer_overlap_junctions_count > 0)
	    && !(retval->rev.oligo[i].overlaps_overlap_position
		 || retval->fwd.oligo[j].overlaps_overlap_position)
	    ) {
	  if (update_stats) { 
	    pair_expl->considered++;
	    pair_expl->does_not_overlap_a_required_point++; 
	  }
	  if (!must_use) continue;
	}

	/* Check product size now */
	double product_size 
	  = retval->rev.oligo[i].start - retval->fwd.oligo[j].start+1;

	if (product_size < pa->pr_min[product_size_range_index] ||
	    product_size > pa->pr_max[product_size_range_index]) {
	  if (update_stats) {
	    /* This line NEW */ if (!must_use)
	                           pair_expl->considered++;
	    pair_expl->product++; 
	  }
	  if (!must_use) continue;
	}

        /* Check if pair was already computed */
        pair_found = 0;
        if (hmap) {
          it = hmap->find(j);
          if (it != hmap->end()) {
            pair_found = 1;

	    /* it->second is the second element (i.e. the 'value', as
	       opposed to the 'key'). */
            pp = it->second; 
            if (pp) { 
	      /* The pair was computed, it isn't illegal and it wasn't
		 selected yet */
	      if (update_stats) {
		pair_expl->considered++;
		if (trace_me)
		  fprintf(stderr, "ok++\n");
		pair_expl->ok++;
	      }
	      /* Check if this is a better pair */
	      if (compare_primer_pair(pp, &the_best_pair) < 0) {
		the_best_pair = *pp;
		the_best_i = i;
		the_best_j = j;
		best_hmap = hmap;
		best_pp = pp;
	      }
	      
	      /* There cannot be a better pair */
	      if (the_best_pair.pair_quality == 0) {
		break;
	      } 
            } /* else - pp is NULL - it's illegal or already selected */
          }
        } else {
          /* Create this hashmap */
          hmap = new std::hash_map<int, primer_pair*>;
          if (!hmap)
            longjmp(_jmp_buf, 1);
          pairs[i] = hmap;
        }

        if (!pair_found) {
          /* Characterize the pair. h is initialized by this
             call. */
          int tmp = 
	    characterize_pair(retval, pa, sa, j, i,
			      product_size_range_index,
			      &h, dpal_arg_to_use,
			      thal_arg_to_use,
			      update_stats);
          if (tmp == PAIR_OK) {

            /* Choose internal oligo if needed */
            if (pa->pick_right_primer && pa->pick_left_primer
                && pa->pick_internal_oligo) {
              if (choose_internal_oligo(retval, h.left, h.right,
                                        &n_int, sa, pa,
                                        dpal_arg_to_use, thal_oligo_arg_to_use)!=0) {

                /* We were UNable to choose an internal oligo. */
                if (update_stats) { 
                  pair_expl->internal++;
                }

                /* Mark the pair as not good - the entry in the hash
		   map will be a NULL */
                (*hmap)[j] = NULL;
                continue;
              } else {
                /* We DID choose an internal oligo, and we
                   set h.intl to point to it. */
                h.intl = &retval->intl.oligo[n_int];
              }
            }

            if (update_stats) { 
              if (trace_me)
                fprintf(stderr, "ok++\n");
              pair_expl->ok++;
            }

            /* Calculate the pair penalty */
            h.pair_quality = obj_fn(pa, &h);
            PR_ASSERT(h.pair_quality >= 0.0);
 
            /* Save the pair */
            pp = new primer_pair;
            if (!pp)
              longjmp(_jmp_buf, 1);
            *pp = h;
            (*hmap)[j] = pp;
              
            /* The current pair (h) is the new best pair if it is
	       better than the best pair so far. */
            if (compare_primer_pair(&h, &the_best_pair) < 0) {
              the_best_pair = h;
              the_best_i = i;
              the_best_j = j;
              best_hmap = hmap;
              best_pp = pp;
            }

            /* There cannot be a better pair */
            if (the_best_pair.pair_quality == 0) {
              break;
            }
          } else if (tmp == PAIR_FAILED) {
            /* Illegal pair */
            (*hmap)[j] = NULL;
          } 
        }
      }  /* for (j=0; j<retval->fwd.num_elem; j++) -- inner loop */

      /* Check if there cannot be a better pair than best found */
      if (the_best_pair.pair_quality == 0) {
        break;
      }

    } /* for (i = 0; i < retval->rev.num_elem; i++) --- outer loop */

    if (the_best_pair.pair_quality == DBL_MAX){
      /* No pair was found. Try another product-size-range,
       if one exists. */

      product_size_range_index++;
      /* Re-set the high-water marks for the indices
         for reverse and forward primers:  */
      for (i = 0; i < retval->rev.num_elem; i++) max_j_seen[i] = -1;

      if (!(product_size_range_index < pa->num_intervals)) {
        /* We ran out of product-size-ranges. Exit the while loop. */
        break;

        /* Our bookkeeping was incorrect unless the assertion below is
	   true. If num_intervals > 1 or min_three_*_prime_distance >
	   -1 the assertion below might not be true. */
        PR_ASSERT(!((pa->num_intervals == 1) && 
		    ((pa->min_left_three_prime_distance == -1) ||
		     (pa->min_right_three_prime_distance == -1)))
                  || (best_pairs->num_pairs == pair_expl->ok));
      }

    } else {
      /* Store the best primer for output */

      if (trace_me)
        fprintf(stderr, "ADD pair i=%d, j=%d\n", the_best_i, the_best_j);

      add_pair(&the_best_pair, best_pairs);

      /* Mark the pair as already selected */
      delete best_pp;
      (*best_hmap)[the_best_j] = NULL;

      /* Update the overlaps flags */
      for (i = 0; i < retval->rev.num_elem; i++) {
	if (right_oligo_in_pair_overlaps_used_oligo(&retval->rev.oligo[i],
						    &the_best_pair,
						    pa->min_right_three_prime_distance)) {
	  retval->rev.oligo[i].overlaps = 1;
	}
      }
      for (j = 0; j < retval->fwd.num_elem; j++) {
	if (left_oligo_in_pair_overlaps_used_oligo(&retval->fwd.oligo[j],
						   &the_best_pair,
						   pa->min_left_three_prime_distance)) {
	  retval->fwd.oligo[j].overlaps = 1;
	}
      }

      /* If we have enough then stop the while loop */
      if (pa->num_return == best_pairs->num_pairs) {
        break;
      }
    }
  } /* end of while(continue_trying == 1) */

  /* Final cleanup of dynamically allocated storage for this
     function. */
  free_pair_memory(retval->rev.num_elem);
}
/* ============================================================ */
/* END choose_pair_or_triple                                    */
/* ============================================================ */

/* ============================================================ */
/* BEGIN choose_internal_oligo */
/* Choose best internal oligo for given pair of left and right
   primers. */
/* ============================================================ */
static int
choose_internal_oligo(p3retval *retval,
                      const primer_rec *left,
                      const primer_rec *right,
                      int *nm,
                      const seq_args *sa,
                      const p3_global_settings *pa,
                      const dpal_arg_holder *dpal_arg_to_use,
                      const thal_arg_holder *thal_arg_to_use
                      )
{
  int i,k;
  double min;
  char oligo_seq[MAX_PRIMER_LENGTH+1], revc_oligo_seq[MAX_PRIMER_LENGTH+1];
  primer_rec *h;
  min = 1000000.;
  i = -1;

  for (k=0; k < retval->intl.num_elem; k++) {
    h = &retval->intl.oligo[k];  /* h is the record for the oligo currently
                                    under consideration */

    if ((h->start > (left->start + (left->length-1)))
        && ((h->start + (h->length-1))
            < (right->start-right->length+1))
        && (h->quality < min)
        && (OK_OR_MUST_USE(h))) {

      if (h->self_any == ALIGN_SCORE_UNDEF && pa->thermodynamic_oligo_alignment==0) {

        _pr_substr(sa->trimmed_seq, h->start, h->length, oligo_seq);
        p3_reverse_complement(oligo_seq, revc_oligo_seq);

        oligo_compl(h, &pa->o_args, &retval->intl.expl,
                    dpal_arg_to_use, oligo_seq, revc_oligo_seq);
        if (!OK_OR_MUST_USE(h)) continue;
      }
       
      if (h->self_any == ALIGN_SCORE_UNDEF && pa->thermodynamic_oligo_alignment==1) {
        _pr_substr(sa->trimmed_seq, h->start, h->length, oligo_seq);
	p3_reverse_complement(oligo_seq, revc_oligo_seq);
          
        oligo_compl_thermod(h, &pa->o_args, &retval->intl.expl,
                            thal_oligo_arg_to_use, oligo_seq, oligo_seq);
        if (!OK_OR_MUST_USE(h)) continue;
      }
      if(h->hairpin_th == ALIGN_SCORE_UNDEF && pa->thermodynamic_oligo_alignment==1) {
	_pr_substr(sa->trimmed_seq, h->start, h->length, oligo_seq);
	oligo_hairpin(h, &pa->o_args,
                      &retval->intl.expl, thal_oligo_arg_to_use,
                      oligo_seq);
        if (!OK_OR_MUST_USE(h)) continue;
      }
       
      if (h->repeat_sim.score == NULL) {
        oligo_repeat_library_mispriming(h, pa, sa, OT_INTL, &retval->intl.expl,
                                        dpal_arg_to_use);
        if (!OK_OR_MUST_USE(h)) continue;
      }

      min = h->quality;
      i=k;

    } /* if ((h->start.... */

  }  /* for (k=0;..... */

  *nm = i;
  if(*nm < 0) return 1;
  return 0;
}
/* ============================================================ */
/* END choose_internal_oligo */
/* ============================================================ */



/* 
   Translate the values in the stats struct into an warning string.
   Call this function only if the 'stat's contains the _errors_
   associated with a given primer i.e. that primer was supplied by the
   caller and pick_anyway is set.
*/
void
add_must_use_warnings(pr_append_str *warning,
                      const char* text,
                      const oligo_stats *stats)
{
  const char *sep = "/";
  pr_append_str s;

  s.data = NULL;
  s.storage_size = 0;

  if (stats->size_min) pr_append_w_sep(&s, sep, "Too short");
  if (stats->size_max) pr_append_w_sep(&s, sep, "Too long");
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
  if (stats->hairpin_th) pr_append_w_sep(&s, sep, "High hairpin stability (thermod. approach)");
  if (stats->repeat_score)
    pr_append_w_sep(&s, sep, "High similarity to mispriming or mishyb library");
  if (stats->poly_x) pr_append_w_sep(&s, sep, "Long poly-X");
  if (stats->seq_quality) pr_append_w_sep(&s, sep, "Low sequence quality");
  if (stats->stability) pr_append_w_sep(&s, sep, "High 3' stability");
  if (stats->no_orf) pr_append_w_sep(&s, sep, "Would not amplify any ORF");
  if (stats->not_in_any_left_ok_region) pr_append_w_sep(&s, sep, "Not in any ok left region");
  if (stats->not_in_any_right_ok_region) pr_append_w_sep(&s, sep, "Not in any ok right region");

  /* edited by T. Koressaar for lowercase masking: */
  if (stats->gmasked)
    pr_append_w_sep(&s, sep, "Masked with lowercase letter");

  if (stats->must_match_fail)
    pr_append_w_sep(&s, sep, "Failed must_match requirements");

  if (s.data) {
    pr_append_new_chunk(warning, text);
    pr_append(warning, " is unacceptable: ");
    pr_append(warning, s.data);
    free(s.data);
  }
}

/*
   Also see the documentation in librimer3.h,
   p3_global_settings.min_{left/right_}three_prime_distance.

   If min_dist == -1, return 0

   If min_dist ==  0:
     Return 1 if the left or the right
     primer is identical to the left or right
     primer (respectively) in best_pair.

     Otherwise return 0.

   If min_dist > 0:
     Return 1 if EITHER:
     (1) The 3' end of the left primer in pair is
         < min_dist from the 3' end of the left primer
         in best_pair
     OR
     (2) The 3' end of the right primer in pair is
         < min_dist from the 3' end of the right primer
         in best_pair

*/

static int
right_oligo_in_pair_overlaps_used_oligo(const primer_rec *right,
                                        const primer_pair *best_pair,
                                        int min_dist)
{
  int best_pos, pair_pos;

  if (min_dist == -1)
    return 0;

  best_pos = best_pair->right->start - best_pair->right->length + 1;

  pair_pos = right->start - right->length + 1;

  if ((abs(best_pos - pair_pos) < min_dist)
      && (min_dist != 0)) { return 1; }
  
  if ((best_pair->right->length == right->length)
      && (best_pair->right->start == right->start)
      && (min_dist == 0)) {
    return 1;
  }

  return 0;
}

static int
left_oligo_in_pair_overlaps_used_oligo(const primer_rec *left,
                                       const primer_pair *best_pair,
                                       int min_dist)
{
  int best_pos, pair_pos;

  if (min_dist == -1)
    return 0;

  best_pos =  best_pair->left->start + best_pair->left->length - 1;

  pair_pos = left->start + left->length -  1;

  if ((abs(best_pos - pair_pos) < min_dist)
      && (min_dist != 0)) { return 1; }
  
  if ((best_pair->left->length == left->length)
      && (best_pair->left->start == left->start)
      && (min_dist == 0)) {
    return 1;
  }

  return 0;
}

/* Add 'pair' to 'retpair'; always appends. */
static void
add_pair(const primer_pair *pair,
         pair_array_t *retpair)
{
  /* If array is not initialised, initialise it with the size of pairs
     to return */
  if (0 == retpair->storage_size) {
    retpair->storage_size = INITIAL_NUM_RETURN;
    retpair->pairs
      = (primer_pair *)
      pr_safe_malloc(retpair->storage_size * sizeof(*retpair->pairs));
  } else {
    if (retpair->storage_size == retpair->num_pairs) {
      /* We need more space, so realloc double the space*/
      retpair->storage_size *= 2;
      retpair->pairs
        = (primer_pair *)
        pr_safe_realloc(retpair->pairs,
                        retpair->storage_size * sizeof(*retpair->pairs));
    }
  }
  /* Copy the pair into the storage place */
  retpair->pairs[retpair->num_pairs] = *pair;
  retpair->num_pairs++;
}

/*
 * Make lists of acceptable left and right primers.  After return, the
 * lists are stored in retval->fwd.oligo and retval->rev.oligo and the
 * coresponding list sizes are stored in retval->fwd.num_elem and
 * retval->rev.num_elem.  Return 1 if one of lists is empty or if
 * leftmost left primer and rightmost right primer do not provide
 * sufficient product size.
 */
static int
make_detection_primer_lists(p3retval *retval,
                            const p3_global_settings *pa,
                            const seq_args *sa,
                            const dpal_arg_holder *dpal_arg_to_use,
                            const thal_arg_holder *thal_arg_to_use)
{
  int left, right;
  int length, start;
  int i,n,/*k,*/ pr_min;
  int tar_l, tar_r, f_b, r_b;
  pair_stats *pair_expl = &retval->best_pairs.expl; /* To store the statistics for pairs */

  /* Var to save the very left and very right primer */
  left = right = 0;

  /* Set pr_min to the very smallest
     allowable product size. */
  pr_min = INT_MAX;
  for (i=0; i < pa->num_intervals; i++)
    if(pa->pr_min[i] < pr_min)
      pr_min = pa->pr_min[i];

  /* Get the length of the sequence */
  PR_ASSERT(INT_MAX > (n=strlen(sa->trimmed_seq)));

  tar_r = 0; /* Target start position */
  tar_l = n; /* Target length */

  /* Iterate over target array */
  for (i=0; i < sa->tar2.count; i++) {

    /* Select the rightmost target start */
    if (sa->tar2.pairs[i][0] > tar_r)
      tar_r = sa->tar2.pairs[i][0];

    /* Select the rightmost target end */
    if (sa->tar2.pairs[i][0] + sa->tar2.pairs[i][1] - 1 < tar_l)
      tar_l = sa->tar2.pairs[i][0] + sa->tar2.pairs[i][1] - 1;
  }

  if (_PR_DEFAULT_POSITION_PENALTIES(pa)) {
    if (0 == tar_r) tar_r = n;
    if (tar_l == n) tar_l = 0;
  } else {
    tar_r = n;
    tar_l = 0;
  }

  /* We use some global information to restrict the region
     of the input sequence in which we generate candidate
     oligos. */
  if (retval->output_type == primer_list && pa->pick_left_primer == 1)
    f_b = n - 1;
  else if (tar_r - 1 < n - pr_min + pa->p_args.max_size - 1
           && !(pa->pick_anyway && sa->left_input))
    f_b=tar_r - 1;
  else
    f_b = n - pr_min + pa->p_args.max_size-1;

  if (pa->pick_left_primer) {
    /* We will need a left primer. */
    left=n; right=0;
    length = f_b - pa->p_args.min_size + 1;
    start = pa->p_args.min_size - 1;

    /* Use the primer provided */
     if (sa->left_input) {
      add_one_primer(sa->left_input, &left, &retval->fwd,
                     pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
         
    }
    /* Pick primers at one position */
    else if(sa->force_left_start > -1 ||
                sa->force_left_end > -1) {
      pick_primers_by_position(sa->force_left_start, sa->force_left_end,
                               &left, &retval->fwd, pa, sa,
                               dpal_arg_to_use, thal_arg_to_use, retval);
    }
    /* Or pick all good in the given range */
    else {
       pick_primer_range(start, length, &left, &retval->fwd,
                        pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
    }

  }  /* if (pa->pick_left_primer) */
   if (retval->output_type == primer_list && pa->pick_right_primer == 1)
    r_b = 0;
  else if (tar_l+1>pr_min - pa->p_args.max_size
           && !(pa->pick_anyway && sa->right_input))
    r_b = tar_l+1;
  else
    r_b = pr_min - pa->p_args.max_size;

  if ( pa->pick_right_primer ) {

    /* We will need a right primer */
    length = n-pa->p_args.min_size - r_b + 1;
    start = r_b;

    /* Use the primer provided */
    if (sa->right_input) {
      add_one_primer(sa->right_input, &right, &retval->rev,
                     pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
      /*  pick_right_primers(start, length, &right, &retval->rev,
          pa, sa, dpal_arg_to_use, retval);*/

    }
    /* Pick primers at one position */
    else if(sa->force_right_start > -1 ||
                sa->force_right_end > -1) {
      pick_primers_by_position(sa->force_right_start, sa->force_right_end,
                               &right, &retval->rev, pa, sa,
                               dpal_arg_to_use, thal_arg_to_use, retval);
    }
    /* Or pick all good in the given range */
    else {
      pick_primer_range(start, length, &right, &retval->rev,
                        pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
    }
  }

  /*
   * Return 1 if either the left primer list or the right primer
   * list is empty or if leftmost left primer and
   * rightmost right primer do not provide sufficient product size.
   */

  if ((pa->pick_left_primer && 0 == retval->fwd.num_elem)
      || ((pa->pick_right_primer)  && 0 == retval->rev.num_elem)) {
    return 1;
  } else if (!((sa->right_input)
                && (sa->left_input))
             && pa->pick_left_primer
             && pa->pick_right_primer
             && (right - left) < (pr_min - 1)) {
    pair_expl->product    = 1;
    pair_expl->considered = 1;
    return 1;
  } else return 0;
} /* make_detection_primer_lists */

/*
 * Make complete list of acceptable internal oligos in retval->intl.oligo.
 * and place the number of valid elements in mid in retval->intl.num_elem.  Return 1 if
 * there are no acceptable internal oligos; otherwise return 0.
 */
static int
make_internal_oligo_list(p3retval *retval,
                         const p3_global_settings *pa,
                         const seq_args *sa,
                         const dpal_arg_holder *dpal_arg_to_use,
                         const thal_arg_holder *thal_arg_to_use)
{
  int ret;
  int left = 0;

  /* Use the primer provided */
  if ((sa->internal_input) || (pa->primer_task == check_primers)){
          ret = add_one_primer(sa->internal_input, &left, &retval->intl,
                               pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
  }
  else {
    /* Pick all good in the given range */

    /* Use the settings to select a proper range */
    int length = strlen(sa->trimmed_seq) - pa->o_args.min_size;
    int start = pa->o_args.min_size - 1;
    int left = 0;

    ret = pick_primer_range(start, length, &left, &retval->intl,
                            pa, sa, dpal_arg_to_use, thal_arg_to_use,
                            retval);
  }
  return ret;
} /* make_internal_oligo_list */

/*
 * Make lists of acceptable left and right primers.  After return, the
 * lists are stored in retval->fwd.oligo and retval->rev.oligo and the
 * coresponding list sizes are stored in retval->fwd.num_elem and
 * retval->rev.num_elem.  Return 1 if one of lists is empty or if
 * leftmost left primer and rightmost right primer do not provide
 * sufficient product size.
 */
static int
make_complete_primer_lists(p3retval *retval,
                  const p3_global_settings *pa,
                  const seq_args *sa,
                  const dpal_arg_holder *dpal_arg_to_use,
                  const thal_arg_holder *thal_arg_to_use,
		  const thal_arg_holder *thal_oligo_arg_to_use)
{
  int extreme;
  int length, start;
  int n;

  /* Get the length of the sequence */
   PR_ASSERT(INT_MAX > (n=strlen(sa->trimmed_seq)));
   if (pa->pick_left_primer) {
    /* We will need a left primer. */
    extreme = 0;
    length = n - pa->p_args.min_size;
    start = pa->p_args.min_size - 1;

    /* Pick all good in the given range */
    pick_primer_range(start, length, &extreme, &retval->fwd,
                      pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);

  }  /* if (pa->pick_left_primer) */
  if ( pa->pick_right_primer ) {
    /* We will need a right primer */
    extreme = n;
    length = n - pa->p_args.min_size + 1;
    start = 0;

    /* Pick all good in the given range */
    pick_primer_range(start, length, &extreme, &retval->rev,
                      pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
  }

  if ( pa->pick_internal_oligo ) {
    /* We will need a internal oligo */
    length = n - pa->o_args.min_size;
    start = pa->o_args.min_size - 1;
    extreme = 0;

    /* Pick all good in the given range */
    pick_primer_range(start, length, &extreme, &retval->intl,
                      pa, sa, dpal_arg_to_use, thal_oligo_arg_to_use, retval);
  }

  return 0;
} /* make_complete_primer_lists */

/*
 * Add caller-specified primers to retval->fwd, retval->rev, and/or
 * retval->intl.  The primers do not have to be "legal". It is possible
 * to add more than one primer to retval->fwd, etc, if the input
 * primer is found in multiple locations in the template. Similar to
 * make_complete_primer_lists.
 */
static int
add_primers_to_check(p3retval *retval,
		     const p3_global_settings *pa,
		     const seq_args *sa,
		     const dpal_arg_holder *dpal_arg_to_use,
		     const thal_arg_holder *thal_arg_to_use,
		     const thal_arg_holder *thal_oligo_arg_to_use)
{
  int extreme = 0; /* Required when calling add_one_primer
		      but not used in the current function. */

  if (sa->left_input) {
    add_one_primer(sa->left_input, &extreme, &retval->fwd,
                   pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
  }

  if (sa->right_input) {
    add_one_primer(sa->right_input, &extreme, &retval->rev,
                   pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
  }

  if (sa->internal_input) {
    add_one_primer(sa->internal_input, &extreme, &retval->intl,
                   pa, sa, dpal_arg_to_use, thal_oligo_arg_to_use, retval);
  }

  return 0;
} /* add_primers_to_check */

/*
 * Make lists of acceptable left and right primers.  After return, the
 * lists are stored in retval->fwd.oligo and retval->rev.oligo and the
 * coresponding list sizes are stored in retval->fwd.num_elem and
 * retval->rev.num_elem.  Return 1 if one of lists is empty or if
 * leftmost left primer and rightmost right primer do not provide
 * sufficient product size.
 */
static int
pick_sequencing_primer_list(p3retval *retval,
                            const p3_global_settings *pa,
                            const seq_args *sa,
                            const dpal_arg_holder *dpal_arg_to_use,
                            const thal_arg_holder *thal_arg_to_use)
{
  int length, start;
  int n, rest_accuracy;
  int primer_nr; /* number of primers we need to pick */
  int tar_n; /* counter for the targets */
  int step_nr; /* counter to step through the targets */
  int sequenced_len; /* bp sequenced in good quality */
  int extra_seq; /* bp sequenced additionally on both sides */

  /* best location for the 3' end of the fwd primer: */
  int pr_position_f;

  /* best location for the 3' end of the rev primer: */
  int pr_position_r;

  /* Get the length of the sequence */
  PR_ASSERT(INT_MAX > (n=strlen(sa->trimmed_seq)));

  /* For each target needed loop*/
  for (tar_n=0; tar_n < sa->tar2.count; tar_n++) {

    /* Calculate the amount of primers needed */
    primer_nr = 1;
    if ((pa->pick_left_primer) && (pa->pick_right_primer)){
      sequenced_len = pa->sequencing.interval;
      while(sequenced_len < sa->tar2.pairs[tar_n][1]) {
        primer_nr++;
        sequenced_len = pa->sequencing.spacing * (primer_nr - 1)
          + pa->sequencing.interval;
      }
    } else {
      sequenced_len = pa->sequencing.spacing;
      while(sequenced_len < sa->tar2.pairs[tar_n][1]) {
        primer_nr++;
        sequenced_len = pa->sequencing.spacing * primer_nr;
      }
    }
    /* Calculate the overlap on the sides */
    extra_seq = (sequenced_len - sa->tar2.pairs[tar_n][1]) / 2;

    /* Pick primers for each position */
    for ( step_nr = 0 ; step_nr < primer_nr ; step_nr++ ) {
      pr_position_f = sa->tar2.pairs[tar_n][0] - extra_seq
        + ( pa->sequencing.spacing * step_nr )
        - pa->sequencing.lead;
      if ((pa->pick_left_primer) && (pa->pick_right_primer)) {
        pr_position_r = sa->tar2.pairs[tar_n][0] - extra_seq
          + ( pa->sequencing.spacing * step_nr )
          + pa->sequencing.interval
          + pa->sequencing.lead;
      } else {
          pr_position_r = sa->tar2.pairs[tar_n][0] - extra_seq
            + ( pa->sequencing.spacing * (step_nr+1))
            + pa->sequencing.lead;        
      }
      /* Check if calculated positions make sense */
      /* position_f cannot be outside included region */
      if (pr_position_f < (pa->p_args.min_size -1)) {
        pr_position_f = pa->p_args.min_size - 1;
      }
      if (pr_position_f > (n - pa->p_args.min_size - 1)) {
        pr_position_f = n - pa->p_args.min_size - 1;
        /* Actually this should never happen */
        pr_append_new_chunk(&retval->warnings,
                            "Calculation error in forward "
                            "sequencing position calculation");
      }
      /* position_r cannot be outside included region */
      if (pr_position_r < (pa->p_args.min_size - 1)) {
        pr_position_r = pa->p_args.min_size - 1;
        /* Actually this should never happen */
        pr_append_new_chunk(&retval->warnings,
                            "Calculation error in reverse "
                            "sequencing position calculation");
      }
      if (pr_position_r > (n - pa->p_args.min_size - 1)) {
        pr_position_r = n - pa->p_args.min_size - 1;
      }
      /* Now all pr_positions are within the sequence */
      if (pa->pick_left_primer) {
        /* Set the start and length for the regions */
        start = pr_position_f - pa->sequencing.accuracy;
        if (start < 0) {
          rest_accuracy = pr_position_f + 1;
          start = 0;
        } else {
          rest_accuracy = pa->sequencing.accuracy;
        }
        length = rest_accuracy + pa->sequencing.accuracy ;
        if ((start + length) > n) {
          length = n - start;
        }
        /* Pick all good in the given range */
        pick_only_best_primer(start, length, &retval->fwd,
                              pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
      }
      if (pa->pick_right_primer) {
        start = pr_position_r - pa->sequencing.accuracy;
        if (start < 0) {
          rest_accuracy = pr_position_r + 1;
          start = 0;
        } else {
          rest_accuracy = pa->sequencing.accuracy;
        }
        length = rest_accuracy + pa->sequencing.accuracy ;
        if ((start + length) > n) {
          length = n - start;
        }
        /* Pick all good in the given range */
        pick_only_best_primer(start, length, &retval->rev,
                              pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
      }
    }

  } /* End of Target Loop */

  /* Print an error if not all primers will be printed */
  if (retval->fwd.num_elem > pa->num_return 
                  || retval->rev.num_elem > pa->num_return) {
    pr_append_new_chunk( &retval->warnings,
     "Increase PRIMER_NUM_RETURN to obtain all sequencing primers");
  }
    
  return 0;
} /* pick_sequencing_primer_list */

static void
add_oligo_to_oligo_array(oligo_array *oarray, primer_rec orec) {
  /* Allocate some space for primers if needed */
  if (NULL == oarray->oligo) {
    oarray->storage_size = INITIAL_LIST_LEN;
    oarray->oligo
      = (primer_rec *)
      pr_safe_malloc(sizeof(*oarray->oligo) * oarray->storage_size);
  }
  /* If there is no space on the array, allocate new space */
  if ((oarray->num_elem + 1) >= oarray->storage_size) { /* TO DO, is +1 really needed? */
    oarray->storage_size += (oarray->storage_size >> 1);
    oarray->oligo
      = (primer_rec *)
      pr_safe_realloc(oarray->oligo,
                      oarray->storage_size * sizeof(*oarray->oligo));
  }
  oarray->oligo[oarray->num_elem] = orec;
  oarray->num_elem++;
}

/* pick_primer_range picks all primers which have their 3' end in the range
 * from start to start+length and store only the best *oligo  */
static int
pick_only_best_primer(const int start,
                      const int length,
                      oligo_array *oligo,
                      const p3_global_settings *pa,
                      const seq_args *sa,
                      const dpal_arg_holder *dpal_arg_to_use,
                      const thal_arg_holder *thal_arg_to_use,
                      p3retval *retval)
{
  /* Variables for the loop */
  int i, j, primer_size_small, primer_size_large;
  int n, found_primer;
  char number[20];
  char *p_number = &number[0];
  int temp_value;

  /* Array to store one primer sequences in */
  char oligo_seq[MAX_PRIMER_LENGTH+1];

  /* Struct to store the primer parameters in */
  primer_rec h;
  primer_rec best;
  memset(&h, 0, sizeof(primer_rec));
  memset(&best, 0, sizeof(primer_rec));
  best.quality = 1000.00;
  found_primer = 0;

  /* Set n to the length of included region */
  PR_ASSERT(INT_MAX > (n=strlen(sa->trimmed_seq)));

  /* Conditions for primer length */
  if (oligo->type == OT_INTL) {
    primer_size_small=pa->o_args.min_size;
    primer_size_large=pa->o_args.max_size;
  }
  else {
    primer_size_small=pa->p_args.min_size;
    primer_size_large=pa->p_args.max_size;
  }

  /* Loop over locations in the sequence */
  for(i = start + length - 1; i >= start; i--) {
    oligo_seq[0] = '\0';

    /* Loop over possible primer lengths, from min to max */
    for (j = primer_size_small; j <= primer_size_large; j++) {
      /* Set the length of the primer */
      h.length = j;

      /* Set repeat_sim to NULL as indicator that the repeat_sim
         struct is not initialized. */
      h.repeat_sim.score = NULL;

      /* Figure out positions for left primers and internal oligos */
      if (oligo->type != OT_RIGHT) {
        /* Break if the primer is bigger than the sequence left */
        if(i-j < -1) continue;

        /* Set the start of the primer */
        h.start = i - j + 1;

        /* Put the real primer sequence in oligo_seq */
        _pr_substr(sa->trimmed_seq, h.start, j, oligo_seq);
      } else {
        /* Figure out positions for reverse primers */
        /* Break if the primer is bigger than the sequence left*/
        if(i+j > n) continue;

        /* Set the start of the primer */
        h.start = i+j-1;

        /* Put the real primer sequence in s */
        _pr_substr(sa->trimmed_seq,  i, j, oligo_seq);

      }

      /* Do not force primer3 to use this oligo */
      h.must_use = 0;

      h.overlaps = 0;

      /* Add it to the considered statistics */
      oligo->expl.considered++;
      /* Calculate all the primer parameters */
      calc_and_check_oligo_features(pa, &h, oligo->type, dpal_arg_to_use, thal_arg_to_use,
                                    sa, &oligo->expl, retval, oligo_seq);

      /* If primer has to be used or is OK */
      if (OK_OR_MUST_USE(&h)) {
        /* Calculate the penalty */
        h.quality = p_obj_fn(pa, &h, oligo->type);
        /* Save the primer in the array */
        if (h.quality < best.quality) {
	  /* Free memory used by previous best primer. */
	  free_primer_repeat_sim_score(&best);
          best = h;
          found_primer = 1;
        } else {
	  /* Free memory used by this primer. */
	  free_primer_repeat_sim_score(&h);
	}
      }
      else {
	/* Free memory used by this primer. */
	free_primer_repeat_sim_score(&h);
	if (any_5_prime_ol_extension_has_problem(&h)) {
	  /* Break from the inner for loop, because there is no
	     legal longer oligo with the same 3' sequence. */
	  break;
	}
      }
    } /* j: Loop over possible primer length from min to max */
  } /* i: Loop over the sequence */

  if (found_primer == 1) {
    /* Add the best to the array */
    add_oligo_to_oligo_array(oligo, best);
    /* Update statistics with how many primers are good */
    oligo->expl.ok = oligo->expl.ok + 1;
  } else {
    if (oligo->type == OT_RIGHT) {
      pr_append_new_chunk(&retval->warnings, "No right primer found in range ");
    } else {
      pr_append_new_chunk(&retval->warnings, "No left primer found in range ");
    }
    temp_value = start + pa->first_base_index;
    sprintf(p_number, "%d", temp_value);
    pr_append(&retval->warnings, p_number);
    pr_append(&retval->warnings, " - ");
    temp_value = start + length + pa->first_base_index;
    sprintf(p_number, "%d", temp_value);
    pr_append(&retval->warnings, p_number);
  }
  if (oligo->num_elem == 0) return 1;
  else return 0;
} /* pick_only_best_primer */

/* pick_primer_range picks all legal primers in the range from start
   to start+length and stores them in *oligo  */
static int
pick_primer_range(const int start, const int length, int *extreme,
                  oligo_array *oligo, const p3_global_settings *pa,
                  const seq_args *sa,
                  const dpal_arg_holder *dpal_arg_to_use,
                  const thal_arg_holder *thal_arg_to_use,
                  p3retval *retval)
{
  /* Variables for the loop */
  int i, j;
  int primer_size_small, primer_size_large;
  int pr_min, n;

  /* Array to store one primer sequences in */
  char oligo_seq[MAX_PRIMER_LENGTH+1];

  /* Struct to store the primer parameters in */
  primer_rec h;
  memset(&h, 0, sizeof(primer_rec));

  /* Set pr_min to the very smallest
     allowable product size. */
  pr_min = INT_MAX;
  for (i=0; i < pa->num_intervals; i++)
    if(pa->pr_min[i] < pr_min)
      pr_min = pa->pr_min[i];

  /* Set n to the length of included region */
  PR_ASSERT(INT_MAX > (n=strlen(sa->trimmed_seq)));

  if (oligo->type == OT_INTL) {
    primer_size_small=pa->o_args.min_size; 
    primer_size_large=pa->o_args.max_size;
  }
  else {
    primer_size_small=pa->p_args.min_size;
    primer_size_large=pa->p_args.max_size;
  }

  /* Loop over locations in the sequence */
  for(i = start + length; i >= start; i--) {
    oligo_seq[0] = '\0';

    /* Loop over possible primer lengths, from min to max */
    for (j = primer_size_small; j <= primer_size_large; j++) {

      /* Set the length of the primer */
      h.length = j;

      /* Figure out positions for left primers and internal oligos */
      if (oligo->type != OT_RIGHT) {
        /* Check if the product is of sufficient size */
        if (i-j > n-pr_min-1 && retval->output_type == primer_pairs
            && oligo->type == OT_LEFT) continue;

        /* Break if the primer is bigger than the sequence left */
        if (i-j < -1) break;

        /* Set the start of the primer */
        h.start = i - j + 1;

        /* Put the real primer sequence in oligo_seq */
        _pr_substr(sa->trimmed_seq, h.start, j, oligo_seq);
      } else {
        /* Figure out positions for reverse primers */
        /* Check if the product is of sufficient size */
        if (i+j < pr_min && retval->output_type == primer_pairs) continue;

        /* Break if the primer is bigger than the sequence left*/
        if (i+j > n) break;

        /* Set the start of the primer */
        h.start = i+j-1;

        /* Put the real primer sequence in s */
        _pr_substr(sa->trimmed_seq,  i, j, oligo_seq);

      }

      /* Do not force primer3 to use this oligo */
      h.must_use = 0;
 
      h.overlaps = 0;

      /* Add it to the considered statistics */
      oligo->expl.considered++;
      /* Calculate all the primer parameters */
      calc_and_check_oligo_features(pa, &h, oligo->type, dpal_arg_to_use, thal_arg_to_use,
                                    sa, &oligo->expl, retval, oligo_seq);
       /* If primer has to be used or is OK */
      if (OK_OR_MUST_USE(&h)) {
        /* Calculate the penalty */
        h.quality = p_obj_fn(pa, &h, oligo->type);
        /* Save the primer in the array */
        add_oligo_to_oligo_array(oligo, h);
        /* Update the most extreme primer variable */
        if (( h.start < *extreme) && (oligo->type != OT_RIGHT))
          *extreme = h.start;
        /* Update the most extreme primer variable */
        if ((h.start > *extreme) && (oligo->type == OT_RIGHT))
          *extreme = h.start;
      } else {
	/* Free memory used by this primer. */
	free_primer_repeat_sim_score(&h);
	if (any_5_prime_ol_extension_has_problem(&h)) {
	  /* Break from the inner for loop, because there is no
	     legal longer oligo with the same 3' sequence. */
	  break;
	}
      }
    } /* j: Loop over possible primer length from min to max */
  } /* i: Loop over the sequence */

  /* Update statistics with how many primers are good */
  oligo->expl.ok = oligo->num_elem;

  if (oligo->num_elem == 0) return 1;
  else return 0;
} /* pick_primer_range */

/* add_one_primer finds one primer in the trimmed sequence and stores
 * it in *oligo The main difference to the general fuction is that it
 * calculates its length and it will add a primer of any length to the
 * list */
static int
add_one_primer(const char *primer, int *extreme, oligo_array *oligo,
               const p3_global_settings *pa,
               const seq_args *sa,
               const dpal_arg_holder *dpal_arg_to_use,
               const thal_arg_holder *thal_arg_to_use,
               p3retval *retval) {
  /* Variables for the loop */
  int i, j;
  int n;

  /* Array to store one primer sequences in */
  char oligo_seq[MAX_PRIMER_LENGTH+1] , test_oligo[MAX_PRIMER_LENGTH+1];

  /* Struct to store the primer parameters in */
  primer_rec h;
  memset(&h, 0, sizeof(primer_rec));

  /* Copy *primer into test_oligo */
  test_oligo[0] = '\0';
  if (oligo->type != OT_RIGHT) {
    strncat(test_oligo, primer, MAX_PRIMER_LENGTH);
  } else {
    p3_reverse_complement(primer, test_oligo);
  }

  PR_ASSERT(INT_MAX > (n=strlen(sa->trimmed_seq)));

  /* This time we already know the size of the primer */
  j = strlen(primer);

  /* Loop over the whole sequence */
  for(i = strlen(sa->trimmed_seq); i >= 0; i--) {
    oligo_seq[0] = '\0';

    /* Set the length of the primer */
    h.length = j;

    /* Figure out positions for forward primers */
    if (oligo->type != OT_RIGHT) {
      /* Break if the primer is bigger than the sequence left*/
      if(i-j < -1) continue;

      /* Set the start of the primer */
      h.start = i - j +1;

      /* Put the real primer sequence in s */
      _pr_substr(sa->trimmed_seq, h.start, j, oligo_seq);
    }
    /* Figure out positions for reverse primers */
    else {
      /* Break if the primer is bigger than the sequence left*/
      if(i+j>n) continue;

      /* Set the start of the primer */
      h.start=i+j-1;

      /* Put the real primer sequence in s */
      _pr_substr(sa->trimmed_seq,  i, j, oligo_seq);
    }

    /* Compare the primer with the sequence */
    if (strcmp_nocase(test_oligo, oligo_seq))
      continue;

    /* Force primer3 to use this oligo */
    h.must_use = (1 && pa->pick_anyway);

    h.overlaps = 0;

    /* Add it to the considered statistics */
    oligo->expl.considered++;

    /* Calculate all the primer parameters */
    calc_and_check_oligo_features(pa, &h, oligo->type, dpal_arg_to_use, thal_arg_to_use,
                                  sa, &oligo->expl, retval, oligo_seq);

    /* If primer has to be used or is OK */
    if (OK_OR_MUST_USE(&h)) {
      /* Calculate the penalty */
      h.quality = p_obj_fn(pa, &h, oligo->type);
      /* Save the primer in the array */
      add_oligo_to_oligo_array(oligo,  h);
      /* Update the most extreme primer variable */
      if ((h.start < *extreme) && (oligo->type != OT_RIGHT))
        *extreme = h.start;
      /* Update the most extreme primer variable */
      if ((h.start > *extreme) && (oligo->type == OT_RIGHT))
        *extreme = h.start;
    } else {
      /* Free memory used by this primer. */
      free_primer_repeat_sim_score(&h);
    }
  } /* i: Loop over the sequence */
    /* Update array with how many primers are good */
  /* Update statistics with how many primers are good */
  oligo->expl.ok = oligo->num_elem;

  if (oligo->num_elem == 0) return 1;
  else {
    if (oligo->num_elem > 1) {
      pr_append_new_chunk(&retval->warnings, 
			  "More than one position in template for input oligo ");
      pr_append(&retval->warnings, primer);
    }
    return 0; /* Success */
  }
}

/* add_one_primer finds one primer in the trimmed sequence and stores
 * it in *oligo The main difference to the general fuction is that it
 * calculates its length and it will add a primer of any length to the
 * list */
static int
add_one_primer_by_position(int start, int length, int *extreme, oligo_array *oligo,
                           const p3_global_settings *pa,
                           const seq_args *sa,
                           const dpal_arg_holder *dpal_arg_to_use,
                           const thal_arg_holder *thal_arg_to_use,
                           p3retval *retval) {
  /* Variables for the loop */
  int i;
  int n, found_primer;

  /* Array to store one primer sequences in */
  char oligo_seq[MAX_PRIMER_LENGTH+1];

  /* Struct to store the primer parameters in */
  primer_rec h;
  memset(&h, 0, sizeof(primer_rec));

  /* Retun 1 for no primer found */
  found_primer = 1;

  PR_ASSERT(INT_MAX > (n=strlen(sa->trimmed_seq)));

  /* Just to be sure */
  if (start < 0) {
    return 1;
  }
  if (start >= n) {
    return 1;
  }
  if (oligo->type != OT_RIGHT) {
    if ((start + length) > n) {
      return 1;
    }
  } else {
    if ((start - length + 1) < 0) {
      return 1;
    }
  }

  oligo_seq[0] = '\0';

  /* Set the length of the primer */
  h.length = length;

  /* Figure out positions for forward primers */
  if (oligo->type != OT_RIGHT) {
    /* Set the start of the primer */
    h.start = start;

    /* Put the real primer sequence in s */
    _pr_substr(sa->trimmed_seq, h.start, length, oligo_seq);
  }
  /* Figure out positions for reverse primers */
  else {
    i = start - length + 1;
    /* Set the start of the primer */
    h.start = start;

    /* Put the real primer sequence in s */
    _pr_substr(sa->trimmed_seq, i, length, oligo_seq);
  }

  /* Force primer3 to use this oligo */
  h.must_use = (1 && pa->pick_anyway);

  h.overlaps = 0;

  /* Add it to the considered statistics */
  oligo->expl.considered++;

  /* Calculate all the primer parameters */
  calc_and_check_oligo_features(pa, &h, oligo->type, dpal_arg_to_use, thal_arg_to_use,
                                sa, &oligo->expl, retval, oligo_seq);

  /* If primer has to be used or is OK */
  if (OK_OR_MUST_USE(&h)) {
    /* Calculate the penalty */
      h.quality = p_obj_fn(pa, &h, oligo->type);
      /* Save the primer in the array */
      add_oligo_to_oligo_array(oligo, h);
      found_primer = 0;
      /* Update the most extreme primer variable */
      if (( h.start < *extreme) &&
          (oligo->type != OT_RIGHT))
        *extreme =  h.start;
      /* Update the most extreme primer variable */
      if (( h.start > *extreme) &&
          (oligo->type == OT_RIGHT))
        *extreme =  h.start;
      /* Update the number of primers */
  } else {
    /* Free memory used by this primer. */
    free_primer_repeat_sim_score(&h);
  }
  /* Update array with how many primers are good */
  /* Update statistics with how many primers are good */
  oligo->expl.ok = oligo->num_elem;

  /* return 0 for success */
  return found_primer;
}

static int
pick_primers_by_position(const int start, const int end, int *extreme,
                         oligo_array *oligo, const p3_global_settings *pa,
                         const seq_args *sa,
                         const dpal_arg_holder *dpal_arg_to_use,
                         const thal_arg_holder *thal_arg_to_use,
                         p3retval *retval)
{
  int found_primer, length, j, ret, new_start;
  found_primer = 1;
  ret = 1;

  if(start > -1 && end > -1) {
    if (oligo->type != OT_RIGHT) {
      length = end - start + 1;
    } else {
      length = start - end + 1;
    }

    found_primer = add_one_primer_by_position(start, length, extreme, oligo,
                                              pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
    return found_primer;
  } else if (start > -1) {
    /* Loop over possible primer lengths, from min to max */
    ret = 0;
    for (j = pa->p_args.min_size; j <= pa->p_args.max_size; j++) {
      ret += add_one_primer_by_position(start, j, extreme, oligo,
                                        pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
      if (ret == 0) {
        found_primer = 0;
      }
    }
    return found_primer;
  } else if (end > -1) {
    /* Loop over possible primer lengths, from min to max */
    ret = 0;
    for (j = pa->p_args.min_size; j <= pa->p_args.max_size; j++) {
      if (oligo->type != OT_RIGHT) {
        new_start = end - j + 1;
      } else {
        new_start = end + j - 1;
      }
      ret += add_one_primer_by_position(new_start, j, extreme, oligo,
                                        pa, sa, dpal_arg_to_use, thal_arg_to_use, retval);
      if (ret == 0) {
        found_primer = 0;
      }
    }
    return found_primer;
  } else {
    /* Actually this should never happen */
    pr_append_new_chunk(&retval->warnings,
           "Calculation error in forced primer position calculation");
    return 1;
  }
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
calc_and_check_oligo_features(const p3_global_settings *pa,
                              primer_rec *h,
                              oligo_type otype,
                              const dpal_arg_holder *dpal_arg_to_use,
                              const thal_arg_holder *thal_arg_to_use,
                              const seq_args *sa,
                              oligo_stats *stats,
                              p3retval *retval,

                              /* This is 5'->3' on the template sequence: */
                              const char *input_oligo_seq
                              )
{
  int i, j, k, for_i, gc_count;
  int three_prime_pos; /* position of 3' base of oligo */
  oligo_type l = otype;
  int poly_x, max_poly_x;
  int must_use = h->must_use;
  int three_conditions 
    = (must_use || pa->file_flag || retval->output_type == primer_list);
  const char *seq = sa->trimmed_seq;
  const thal_args *thal_args_for_template_mispriming 
    = use_end_for_th_template_mispriming 
    ? thal_arg_to_use->end1
    : thal_arg_to_use->any;

  char s1_rev[MAX_PRIMER_LENGTH+1];
  const char *oligo_seq;
  const char *revc_oligo_seq;

  const args_for_one_oligo_or_primer *po_args;
   
  /* Initialize slots in h */
  initialize_op(h);
  h->overlaps = 0;

  /* Set repeat_sim to NULL as indicator that the repeat_sim
     struct is not initialized. */
  h->repeat_sim.score = NULL;

  h->gc_content = h->num_ns = 0;
  h->overlaps_overlap_position = 0;
  h->template_mispriming = h->template_mispriming_r = ALIGN_SCORE_UNDEF;
  h->template_mispriming_ok = 0;
   
  PR_ASSERT(OT_LEFT == l || OT_RIGHT == l || OT_INTL == l);
   
  p3_reverse_complement(input_oligo_seq, s1_rev);
  if (OT_RIGHT == l) {
    oligo_seq = s1_rev;
    revc_oligo_seq = input_oligo_seq;
 } else {
    oligo_seq = input_oligo_seq;
    revc_oligo_seq = s1_rev;
  }

  if (OT_INTL == l) {
    po_args = &pa->o_args;
  } else {
    po_args = &pa->p_args;
  }

  /* Set j and k, and sanity check */
  if (OT_LEFT == otype || OT_INTL == otype) {
    j = h->start;
    three_prime_pos = k = j+h->length-1;
  }  else {
    three_prime_pos = j = h->start-h->length+1;
    k = h->start;
  }
  PR_ASSERT(k >= 0);
  PR_ASSERT(k < TRIMMED_SEQ_LEN(sa));

  if ((otype == OT_LEFT) && !PR_START_CODON_POS_IS_NULL(sa)
      /* Make sure the primer would amplify at least part of
         the ORF. */
      && (0 != (h->start - sa->start_codon_pos) % 3
          || h->start <= retval->upstream_stop_codon
          || (retval->stop_codon_pos != -1
              && h->start >= retval->stop_codon_pos))) {
    stats->no_orf++;
    op_set_does_not_amplify_orf(h);
    if (!pa->pick_anyway) return;
  }

  /* edited by T. Koressaar for lowercase masking */
  if(pa->lowercase_masking == 1) {
    if (is_lowercase_masked(three_prime_pos,
                            sa->trimmed_orig_seq,
                            h, stats)) {
      if (!must_use) return;
    }
  }
  /* end T. Koressar's changes */

  /* edited by A. Untergasser for forcing sequence use */
  if ((po_args->must_match_five_prime != NULL) or
      (po_args->must_match_three_prime != NULL)) {
    if (primer_must_match(pa, h, stats, oligo_seq,
                          po_args->must_match_three_prime,
                          po_args->must_match_five_prime)) {
      if (!must_use) {
	op_set_must_match_err(h);
	return;
      }
    }
  }
  /* end A. Untergasser's changes */

  gc_and_n_content(j, k-j+1, sa->trimmed_seq, h);

  if (h->num_ns > po_args->num_ns_accepted) {
    op_set_too_many_ns(h);
    stats->ns++;
    if (!must_use) return;
  }

  /* Upstream error checking has ensured that we use non-default position
     penalties only when there is 0 or 1 target. */
  PR_ASSERT(sa->tar2.count <= 1 || _PR_DEFAULT_POSITION_PENALTIES(pa));

  if (pa->primer_task == pick_sequencing_primers) {
    h->position_penalty = 0.0;
  } else if (l != OT_INTL
             && _PR_DEFAULT_POSITION_PENALTIES(pa)
             && oligo_overlaps_interval(j, k-j+1, sa->tar2.pairs, sa->tar2.count)) {
    h->position_penalty = 0.0;
    bf_set_infinite_pos_penalty(h,1);
    bf_set_overlaps_target(h,1);
  } else if (l != OT_INTL && !_PR_DEFAULT_POSITION_PENALTIES(pa)
             && 1 == sa->tar2.count) {
    compute_position_penalty(pa, sa, h, l);
    if (bf_get_infinite_pos_penalty(h)) {
        bf_set_overlaps_target(h,1);
    }
  } else {
    h->position_penalty = 0.0;
    bf_set_infinite_pos_penalty(h,0);
  }

  if (!PR_START_CODON_POS_IS_NULL(sa)) {
    if (OT_LEFT == l) {
      if (sa->start_codon_pos > h->start)
        h->position_penalty
          = (sa->start_codon_pos - h->start) * OUTSIDE_START_WT;
      else
        h->position_penalty
          = (h->start - sa->start_codon_pos) * INSIDE_START_WT;
    } else if (OT_RIGHT == l) {

      if (-1 == retval->stop_codon_pos) {
        h->position_penalty = (TRIMMED_SEQ_LEN(sa) - h->start - 1) * INSIDE_STOP_WT;
      } else if (retval->stop_codon_pos < h->start) {
        h->position_penalty
          = (h->start - retval->stop_codon_pos) * OUTSIDE_STOP_WT;
      } else {
        h->position_penalty
          = (retval->stop_codon_pos - h->start) * INSIDE_STOP_WT;
      }
    }
  }

  /* TO DO Simplify logic here */
  if (l != OT_INTL
      && oligo_overlaps_interval(j, k-j+1, sa->excl2.pairs, sa->excl2.count))
    bf_set_overlaps_excl_region(h,1);

  if (l == OT_INTL
      && oligo_overlaps_interval(j, k-j+1, sa->excl_internal2.pairs,
                                              sa->excl_internal2.count))
    bf_set_overlaps_excl_region(h,1);

  if(l != OT_INTL && bf_get_overlaps_target(h)) {
    op_set_overlaps_target(h);
    stats->target++;
    if (!must_use) return;
  }

  if(bf_get_overlaps_excl_region(h)){
    op_set_overlaps_excluded_region(h);
    stats->excluded++;
    if (!must_use) return;
  }

  /* Check if the oligo is included in any ok region */
  if ((l == OT_LEFT) && (sa->ok_regions.count > 0) && (!sa->ok_regions.any_left)) {
    int included = 0;
    for (i=0; i<sa->ok_regions.count; i++) {
      if ((j >= sa->ok_regions.left_pairs[i][0]) && 
          (k <= sa->ok_regions.left_pairs[i][0] + sa->ok_regions.left_pairs[i][1] - 1)) {
        included = 1;
        break;
      }
    }
    if (!included) {
      op_set_not_in_any_ok_region(h);
      stats->not_in_any_left_ok_region++;
      if (!must_use) return;
    }
  }
  if ((l == OT_RIGHT) && (sa->ok_regions.count > 0) && (!sa->ok_regions.any_right)) {
    int included = 0;
    for (i=0; i<sa->ok_regions.count; i++) {
      if ((j >= sa->ok_regions.right_pairs[i][0]) && 
          (k <= sa->ok_regions.right_pairs[i][0] + sa->ok_regions.right_pairs[i][1] - 1)) {
        included = 1;
        break;
      }
    }
    if (!included) {
      op_set_not_in_any_ok_region(h);
      stats->not_in_any_right_ok_region++;
      if (!must_use) return;
    }
  }

  if (h->gc_content < po_args->min_gc) {
    op_set_low_gc_content(h);
    stats->gc++;
    if (!must_use) return;
  } else if (h->gc_content > po_args->max_gc) {
    op_set_high_gc_content(h);
    stats->gc++;
    if (!must_use) return;
  }

  if (OT_LEFT == l || OT_RIGHT == l) {
    /* gc_clamp is applicable only to primers (as opposed
       to primers and hybridzations oligos. */
    for (i = 0; i < pa->gc_clamp; i++) {
      /* We want to look at the 3' end of the oligo being
         assessed, so we look in the 5' end of its reverse-complement */
      if (revc_oligo_seq[i] != 'G' && revc_oligo_seq[i] != 'C') {
        op_set_no_gc_glamp(h);
        stats->gc_clamp++;
        if (!must_use) return; else break;
      }
    }
  }

#if 0
  /* Old code, works */
  if(pa->gc_clamp != 0){
    if(OT_LEFT == l) {
      for(i=k-pa->gc_clamp+1; i<= k; i++)if(seq[i] !='G'&&seq[i] !='C'){
        op_set_no_gc_glamp(h);
        stats->gc_clamp++;
        if (!must_use) return; else break;
      }
    }
    if(OT_RIGHT == l){
      for(i=j; i<j+pa->gc_clamp; i++)if(seq[i] != 'G' && seq[i] != 'C'){
        op_set_no_gc_glamp(h);
        stats->gc_clamp++;
        if (!must_use) return; else break;
      }
    }
  }
#endif

  if (pa->max_end_gc < 5 
      && (OT_LEFT == l || OT_RIGHT == l)) {
    /* The CGs are only counted in the END of primers (as opposed
       to primers and hybridzations oligos. */
    gc_count=0;
    for (i = 0; i < 5; i++) {
      /* We want to look at the 3' end of the oligo being
         assessed, so we look in the 5' end of its reverse-complement */
      if (revc_oligo_seq[i] == 'G' || revc_oligo_seq[i] == 'C') {
        gc_count++;
      }
    }
    if (gc_count > pa->max_end_gc) {
      op_set_too_many_gc_at_end(h);
      stats->gc_end_high++;
      if (!must_use) return;
    }
  }

  /* Tentative interface for generic oligo  testing
     function */
  if (!sequence_quality_is_ok(pa, h, l, sa, j, k,
                              stats, po_args)
      && !must_use) return;

  max_poly_x = po_args->max_poly_x;
  if (max_poly_x > 0) {
    poly_x = 1;
    for(i=j+1;i<=k;i++){
      if(seq[i] == seq[i-1]||seq[i] == 'N'){
        poly_x++;
        if(poly_x > max_poly_x){
          op_set_high_poly_x(h);
          stats->poly_x++;
          if (!must_use) return; else break;
        }
      }
      else poly_x = 1;
    }
  }

  h->temp  /* Oligo/primer melting temperature */
    = seqtm(oligo_seq, po_args->dna_conc,
            po_args->salt_conc,
            po_args->divalent_conc,
            po_args->dntp_conc,
            MAX_NN_TM_LENGTH,
            pa->tm_santalucia,
            pa->salt_corrections);

  if (h->temp < po_args->min_tm) {
    op_set_low_tm(h);
    stats->temp_min++;
    if (!must_use) return;
  }

  if (h->temp > po_args->max_tm) {
    op_set_high_tm(h);
    stats->temp_max++;
    if (!must_use) return;
  }

  /* End stability is applicable only to primers
     (not to oligos) */
  if (OT_LEFT == l || OT_RIGHT == l) {
    if ((h->end_stability = end_oligodg(oligo_seq, 5,
                                        pa->tm_santalucia))
        > pa->max_end_stability) {
      /* Must not be called on a hybridization probe / internal oligo: */
      op_set_high_end_stability(h);
      stats->stability++;
      if (!must_use) return;
    }
  }

  if ((must_use
      || pa->file_flag
      || retval->output_type == primer_list
      || po_args->weights.compl_any
      || po_args->weights.compl_end)
      && pa->thermodynamic_oligo_alignment==0) {

    oligo_compl(h, po_args,
                stats, dpal_arg_to_use,
                oligo_seq, revc_oligo_seq);

    if ((!(p3_ol_is_uninitialized(h))) && !must_use) {
      PR_ASSERT(!p3_ol_is_ok(h));
      return;
    }
  } else {
    /* Thermodynamical approach: for primers only  */
    if ((must_use
       || pa->file_flag
       || retval->output_type == primer_list
       || po_args->weights.compl_any_th
       || po_args->weights.compl_end_th)
       && pa->thermodynamic_oligo_alignment==1
       ) 
     {
        oligo_compl_thermod(h, po_args,
                            stats, thal_arg_to_use,
                            oligo_seq, oligo_seq);
                            /* input_oligo_seq, input_oligo_seq);*/
                            /* oligo_seq, revc_oligo_seq); */
        if ((!(p3_ol_is_uninitialized(h))) && !must_use) {
           PR_ASSERT(!p3_ol_is_ok(h));
           return;
        }
     } else  {
        h->self_any = h->self_end  = ALIGN_SCORE_UNDEF;
     }
  }

  if ((three_conditions
       || po_args->weights.hairpin_th
#if 0
        || po_args->weights.compl_any_th /* Triinu, is this needed? */
        || po_args->weights.compl_end_th /* Triinu, is this needed? */
#endif
       )
       && pa->thermodynamic_oligo_alignment==1
       ) {
      oligo_hairpin(h, po_args,
                    stats, thal_arg_to_use,
                    /* input_oligo_seq);*/
                    oligo_seq);
      if ((!(p3_ol_is_uninitialized(h))) && !must_use) {
         PR_ASSERT(!p3_ol_is_ok(h));
         return;
      }
   } else   {
    /* This will get calculated later if necessary, in characterize_pair. */
    h->hairpin_th = ALIGN_SCORE_UNDEF;
   }
   /* end of thermod. approach */


#if 0
  /* FIX ME FIXME */
  if (((must_use
        || pa->file_flag
        || retval->output_type == primer_list
        || po_args->weights.repeat_sim
        || ((OT_RIGHT == l || OT_LEFT == l)
            && pa->p_args.weights.template_mispriming))

       && pa->thermodynamic_oligo_alignment==0) 
      || (pa->thermodynamic_oligo_alignment==1 && po_args->weights.repeat_sim)) {

    oligo_repeat_library_mispriming(h, pa, sa, l, stats,
                                    dpal_arg_to_use);
    if (pa->thermodynamic_template_alignment==0 && OK_OR_MUST_USE(h)) {
      oligo_template_mispriming(h, pa, sa, l, stats,
				dpal_arg_to_use->local_end,
				thal_args_for_template_mispriming);
    }
  }
  
  if ((must_use
       || pa->file_flag
       || retval->output_type == primer_list
       || ((OT_RIGHT == l || OT_LEFT == l)
	   && pa->p_args.weights.template_mispriming_th))

      && pa->thermodynamic_template_alignment==1)
 {

    oligo_template_mispriming(h, pa, sa, l, stats,
                              dpal_arg_to_use->local_end,
                              thal_args_for_template_mispriming);
  }
  /* END FIX ME FIXME */   
#else

  if (three_conditions || po_args->weights.repeat_sim) {
    oligo_repeat_library_mispriming(h, pa, sa, l, stats,
                                    dpal_arg_to_use);
  }


  if (!OK_OR_MUST_USE(h)) return;

  if (three_conditions
      || 
      ( /* Do we need template mispriming for the penalty function? */
       (OT_RIGHT == l || OT_LEFT == l)
       && 
       (
	(pa->p_args.weights.template_mispriming && !pa->thermodynamic_template_alignment)
	||
	(pa->p_args.weights.template_mispriming_th && pa->thermodynamic_template_alignment)
	)
       )
      ) {
    if (OK_OR_MUST_USE(h)) {
      oligo_template_mispriming(h, pa, sa, l, stats,
				dpal_arg_to_use->local_end,
				thal_args_for_template_mispriming);
    }
  }

#endif

  if (h->length > po_args->max_size ) {
    op_set_too_long(h);
    stats->size_max ++;
    if (!must_use) return;
  }

  if (h->length < po_args->min_size ) {
    op_set_too_short(h);
    stats->size_min ++;
    if (!must_use) return;
  }
  
  for (for_i=0; for_i < sa->primer_overlap_junctions_count; for_i++) {
    if (OT_LEFT == l 
        && ((h->start + pa->min_5_prime_overlap_of_junction - 1) 
            <= sa->primer_overlap_junctions[for_i])
        && ((h->start + h->length - pa->min_3_prime_overlap_of_junction))
        > sa->primer_overlap_junctions[for_i]) {
      h->overlaps_overlap_position = 1;
      /* no need to continue checking */
      break;
    }
    if (OT_RIGHT == l
        && ((h->start - h->length + pa->min_3_prime_overlap_of_junction) 
            <= sa->primer_overlap_junctions[for_i])
        && ((h->start - pa->min_5_prime_overlap_of_junction + 1))
        > sa->primer_overlap_junctions[for_i]) {
      h->overlaps_overlap_position = 1;
      /* no need to continue checking */
      break;
    }
  }

  /* FIX ME FIXME Steve, is this really needed? */
  op_set_completely_written(h);

} /* calc_and_check_oligo_features */
#undef OUTSIDE_START_WT
#undef INSIDE_START_WT
#undef INSIDE_STOP_WT
#undef OUTSIDE_STOP_WT

/* Calculate the minimum sequence quality and the minimum sequence
   quality of the 3' end of a primer or oligo.  Set h->seq_quality
   and h->seq_end_quality with these values.  Return 1 (== ok)
   if sa->quality is undefined.  Otherwise, return 1 (ok)
   if h->seq_quality and h->seq_end_quality are within
   range, or else return 0.
*/
static int
sequence_quality_is_ok(const p3_global_settings *pa,
                       primer_rec *h,
                       oligo_type l,
                       const seq_args *sa,
                       int j, int k,
                       oligo_stats *global_oligo_stats,
                       const args_for_one_oligo_or_primer *po_args) {
  int i, min_q, min_q_end, m, q;
  int  retval = 1;

  if (NULL == sa->quality) {
    h->seq_end_quality = h->seq_quality = pa->quality_range_max;
    return 1;
  }

  q = pa->quality_range_max;

  min_q = po_args->min_quality;
  if (OT_LEFT == l || OT_RIGHT == l) {
    min_q_end = po_args->min_end_quality;
  }  else {
    min_q_end = min_q;
  }

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

  h->seq_quality = min_q;
  h->seq_end_quality = min_q_end;

  if (h->seq_quality < po_args->min_quality) {
    op_set_low_sequence_quality(h);
    global_oligo_stats->seq_quality++;
    retval = 0;
    return retval;
  }

  if (OT_LEFT == l || OT_RIGHT == l) {
    if (h->seq_end_quality < po_args->min_end_quality) {
      op_set_low_end_sequence_quality(h);
      global_oligo_stats->seq_quality++;
      retval = 0;
    }
  }
  return retval;
}

/*
 * Set h->gc_content to the GC % of the 'len' bases starting at 'start' in
 * 'sequence' (excluding 'N's).  Set h->num_ns to the number of 'N's in that
 * subsequence.
 */
static void
gc_and_n_content(int start, int len,
                 const char *sequence,
                 primer_rec *h)
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
oligo_overlaps_interval(int start,
                        int len,
                        const interval_array_t intervals,
                        int num_intervals)
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
p_obj_fn(const p3_global_settings *pa,
         primer_rec *h,
         int  j) {
  double sum;

  sum = 0;
  if (j == OT_LEFT || j == OT_RIGHT) {
      if(pa->p_args.weights.temp_gt && h->temp > pa->p_args.opt_tm)
           sum += pa->p_args.weights.temp_gt * (h->temp - pa->p_args.opt_tm);
      if (pa->p_args.weights.temp_lt && h->temp < pa->p_args.opt_tm)
           sum += pa->p_args.weights.temp_lt * (pa->p_args.opt_tm - h->temp);

      if (pa->p_args.weights.gc_content_gt && h->gc_content > pa->p_args.opt_gc_content)
           sum += pa->p_args.weights.gc_content_gt
             * (h->gc_content - pa->p_args.opt_gc_content);
      if (pa->p_args.weights.gc_content_lt && h->gc_content < pa->p_args.opt_gc_content)
           sum += pa->p_args.weights.gc_content_lt
             * (pa->p_args.opt_gc_content - h->gc_content);

      if (pa->p_args.weights.length_lt && h->length < pa->p_args.opt_size)
           sum += pa->p_args.weights.length_lt * (pa->p_args.opt_size - h->length);
      if (pa->p_args.weights.length_gt && h->length > pa->p_args.opt_size)
           sum += pa->p_args.weights.length_gt * (h->length - pa->p_args.opt_size);

      /* BEGIN: secondary structures */
      if (pa->thermodynamic_oligo_alignment==0) {

	if (pa->p_args.weights.compl_any)
           sum += pa->p_args.weights.compl_any * h->self_any;
	if (pa->p_args.weights.compl_end)
           sum += pa->p_args.weights.compl_end * h->self_end;

      } else if (pa->thermodynamic_oligo_alignment==1) {

	if (pa->p_args.weights.compl_any_th) {
	  if ((h->temp - pa->p_args.weights.temp_cutoff) <= h->self_any)
	    sum += pa->p_args.weights.compl_any_th 
	      * (h->self_any - (h->temp - pa->p_args.weights.temp_cutoff - 1.0)); /* -1.0 is added for the case where == */
	  else
	    sum += pa->p_args.weights.compl_any_th 
	      * (1/(h->temp - pa->p_args.weights.temp_cutoff + 1.0 - h->self_any));
	}
     
	if (pa->p_args.weights.compl_end_th) {
	  if ((h->temp - pa->p_args.weights.temp_cutoff) <= h->self_end)
	    sum += pa->p_args.weights.compl_end_th 
	      * (h->self_end - (h->temp - pa->p_args.weights.temp_cutoff - 1.0));
	  else
	    sum += pa->p_args.weights.compl_end_th 
	      * (1/(h->temp - pa->p_args.weights.temp_cutoff + 1.0 - h->self_end));
	}
     
	if (pa->p_args.weights.hairpin_th) {
	  if ((h->temp - pa->p_args.weights.temp_cutoff) <= h->hairpin_th)
	    sum += pa->p_args.weights.hairpin_th 
	      * (h->hairpin_th - (h->temp - pa->p_args.weights.temp_cutoff - 1.0));
	  else
	    sum += pa->p_args.weights.hairpin_th 
	      * (1/(h->temp - pa->p_args.weights.temp_cutoff + 1.0 - h->hairpin_th));
	}

      } else {
	PR_ASSERT(0) 
	  /* Programming error, 
	     illegal value for pa->thermodynamic_oligo_alignment */
      }
       /* END: secondary structures */
     
      if (pa->p_args.weights.num_ns)
           sum += pa->p_args.weights.num_ns * h->num_ns;
      if(pa->p_args.weights.repeat_sim)
           sum += pa->p_args.weights.repeat_sim
             * h->repeat_sim.score[h->repeat_sim.max];
      if (!bf_get_overlaps_target(h)) {
        /* We might be evaluating p_obj_fn with h->target if
           the client supplied 'pick_anyway' and specified
           a primer or oligo. */
        PR_ASSERT(!(bf_get_infinite_pos_penalty(h)));
        if(pa->p_args.weights.pos_penalty)
          sum += pa->p_args.weights.pos_penalty * h->position_penalty;
      }
      if(pa->p_args.weights.end_stability)
           sum += pa->p_args.weights.end_stability * h->end_stability;

      /* FIX ME QUALITY WT Change here */
      if(pa->p_args.weights.seq_quality)
           sum += pa->p_args.weights.seq_quality *
	     (pa->quality_range_max - h->seq_quality);  /* Look for end seq quality */

      if (pa->p_args.weights.template_mispriming && pa->thermodynamic_template_alignment==0) {
        PR_ASSERT(oligo_max_template_mispriming(h) != ALIGN_SCORE_UNDEF);
        sum += pa->p_args.weights.template_mispriming *
          oligo_max_template_mispriming(h);
      }
                                            
      if (pa->p_args.weights.template_mispriming_th && pa->thermodynamic_template_alignment==1) {

	PR_ASSERT(oligo_max_template_mispriming_thermod(h) != ALIGN_SCORE_UNDEF);

	if((h->temp - pa->p_args.weights.temp_cutoff) <= oligo_max_template_mispriming_thermod(h))
	  sum += pa->p_args.weights.template_mispriming_th * 
	    (oligo_max_template_mispriming_thermod(h) - (h->temp - pa->p_args.weights.temp_cutoff - 1.0));
	if((h->temp - pa->p_args.weights.temp_cutoff) > oligo_max_template_mispriming_thermod(h))
	  sum += pa->p_args.weights.template_mispriming_th *
	    (1/(h->temp - pa->p_args.weights.temp_cutoff + 1.0 - oligo_max_template_mispriming_thermod(h)));
      }
      return sum;
  } else if (j == OT_INTL) {
      if(pa->o_args.weights.temp_gt && h->temp > pa->o_args.opt_tm)
         sum += pa->o_args.weights.temp_gt * (h->temp - pa->o_args.opt_tm);
      if(pa->o_args.weights.temp_lt && h->temp < pa->o_args.opt_tm)
         sum += pa->o_args.weights.temp_lt * (pa->o_args.opt_tm - h->temp);

      if(pa->o_args.weights.gc_content_gt && h->gc_content > pa->o_args.opt_gc_content)
           sum += pa->o_args.weights.gc_content_gt
             * (h->gc_content - pa->o_args.opt_gc_content);
      if(pa->o_args.weights.gc_content_lt && h->gc_content < pa->o_args.opt_gc_content)
           sum += pa->o_args.weights.gc_content_lt
             * (pa->o_args.opt_gc_content - h->gc_content);

      if(pa->o_args.weights.length_lt && h->length < pa->o_args.opt_size)
         sum += pa->o_args.weights.length_lt * (pa->o_args.opt_size - h->length);
      if(pa->o_args.weights.length_gt && h->length  > pa->o_args.opt_size)
         sum += pa->o_args.weights.length_gt * (h->length - pa->o_args.opt_size);
      if(pa->o_args.weights.compl_any && pa->thermodynamic_oligo_alignment==0)
         sum += pa->o_args.weights.compl_any * h->self_any;
      if(pa->o_args.weights.compl_end && pa->thermodynamic_oligo_alignment==0)
         sum += pa->o_args.weights.compl_end * h->self_end;

      /* begin thermodynamical approach */
      if (pa->thermodynamic_oligo_alignment==1) {

	if (pa->o_args.weights.compl_any_th) {
	  if ((h->temp - pa->o_args.weights.temp_cutoff) <= h->self_any)
	    sum += 
	      pa->o_args.weights.compl_any_th 
	      * (h->self_any - (h->temp - pa->o_args.weights.temp_cutoff - 1.0)); /* -1.0 is added for the case where == */
	  else 
	    sum += 
	      pa->o_args.weights.compl_any_th 
	      * (1/(h->temp - pa->o_args.weights.temp_cutoff + 1.0 - h->self_any));
	}

	if (pa->o_args.weights.compl_end_th) {
	  if ((h->temp - pa->o_args.weights.temp_cutoff) <= h->self_end) 
	    sum += 
	      pa->o_args.weights.compl_end_th 
	      * (h->self_end - (h->temp - pa->o_args.weights.temp_cutoff - 1.0));
	  else
	    sum += 
	      pa->o_args.weights.compl_end_th 
	      * (1/(h->temp - pa->o_args.weights.temp_cutoff + 1.0 - h->self_end));
	}
     
	if (pa->o_args.weights.hairpin_th) {
	  if ((h->temp - pa->o_args.weights.temp_cutoff) <= h->hairpin_th) 
	    sum += pa->o_args.weights.hairpin_th * (h->hairpin_th - (h->temp - pa->o_args.weights.temp_cutoff - 1.0));
	  else
	    sum += pa->o_args.weights.hairpin_th * (1/(h->temp - pa->o_args.weights.temp_cutoff + 1.0 - h->hairpin_th));
	}
      }
      /* end thermodynamical approach */

      if(pa->o_args.weights.num_ns)
         sum += pa->o_args.weights.num_ns * h->num_ns;
      if(pa->o_args.weights.repeat_sim)
         sum += pa->o_args.weights.repeat_sim
           * h->repeat_sim.score[h->repeat_sim.max];

      /* FIX ME QUALITY WT */
      if(pa->o_args.weights.seq_quality)
         sum += pa->o_args.weights.seq_quality *
                     (pa->quality_range_max - h->seq_quality);

      return sum;
  } else {
    PR_ASSERT(0); /* Programmig error. */
  }
}

/* Return max of h->template_mispriming and h->template_mispriming_r (max
   template mispriming on either strand). */
double 
oligo_max_template_mispriming(const primer_rec *h) {
   return h->template_mispriming > h->template_mispriming_r ?
     h->template_mispriming : h->template_mispriming_r;
}

double 
oligo_max_template_mispriming_thermod(const primer_rec *h) {
   return h->template_mispriming > h->template_mispriming_r ?
     h->template_mispriming : h->template_mispriming_r;
}
   
/* Sort a given primer array by penalty */
static void
sort_primer_array(oligo_array *oligo)
{
  qsort(&oligo->oligo[0], oligo->num_elem, sizeof(*oligo->oligo),
        primer_rec_comp);
}

/* Compare function for sorting primer records. */
static int
primer_rec_comp(const void *x1, const void *x2)
{
  const primer_rec *a1 = (const primer_rec *) x1;
  const primer_rec *a2 = (const primer_rec *) x2;

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
compare_primer_pair(const void *x1, const void *x2)
{
  const primer_pair *a1 = (const primer_pair *) x1;
  const primer_pair *a2 = (const primer_pair *) x2;
  const double epsilon = 1e-6;
  int y1, y2;

  if ((a1->pair_quality + epsilon) < a2->pair_quality) return -1;
  if (a1->pair_quality > (a2->pair_quality + epsilon)) return 1;

  /*
   * The following statements ensure that we get a stable order that
   * is the same on all systems.
   */

  y1 = a1->left->start;
  y2 = a2->left->start;
  if (y1 > y2) return -1; /* prefer left primers to the right. */
  if (y1 < y2) return 1;

  y1 = a1->right->start;
  y2 = a2->right->start;
  if (y1 < y2) return -1; /* prefer right primers to the left. */
  if (y1 > y2) return 1;

  y1 = a1->left->length;
  y2 = a2->left->length;
  if (y1 < y2) return -1; /* prefer shorter primers. */
  if (y1 > y2) return 1;

  y1 = a1->right->length;
  y2 = a2->right->length;
  if (y1 < y2) return -1; /* prefer shorter primers. */
  if (y1 > y2) return 1;

  return 0;
}

/*
 * Defines parameter values for given primer pair. Returns PAIR_OK if the pair is
 * acceptable; PAIR_FAILED otherwise.  This function sets the
 * various elements of *ppair, and also calculates some key EXPENSIVE
 * parameters for individual primers if they have not been set yet.
 */
static int
characterize_pair(p3retval *retval,
                  const p3_global_settings *pa,
                  const seq_args *sa,
                  int m, int n, int int_num,
                  primer_pair *ppair,
                  const dpal_arg_holder *dpal_arg_to_use,
                  const thal_arg_holder *thal_arg_to_use,
                  int update_stats) {
  char s1[MAX_PRIMER_LENGTH+1], s2[MAX_PRIMER_LENGTH+1],
    s1_rev[MAX_PRIMER_LENGTH+1], s2_rev[MAX_PRIMER_LENGTH+1];
  double compl_end;
  pair_stats *pair_expl = &retval->best_pairs.expl;
  int must_use = 0;
  double min_oligo_tm;
  int i;
  const thal_args *thal_args_for_template_mispriming 
    = use_end_for_th_template_mispriming 
    ? thal_arg_to_use->end1
    : thal_arg_to_use->any;

  memset(ppair, 0, sizeof(*ppair));

  ppair->left = &retval->fwd.oligo[m];
  ppair->right = &retval->rev.oligo[n];
  ppair->product_size = retval->rev.oligo[n].start - retval->fwd.oligo[m].start+1;

  ppair->target = 0;
  ppair->compl_any = ppair->compl_end = 0;
  if (update_stats) { 
    pair_expl->considered++;
  }

  if (pa->primer_task == check_primers) {
    must_use = 1;
  }

  /* We must use the pair if the caller specifed
     both the left and the right primer. */
  if (ppair->left->must_use != 0 &&
      ppair->right->must_use != 0) {

    /* But if caller specified primers are in
       reversed order we cannot analyze them
       as a pair. */
    if (ppair->product_size < 1) {
      pair_expl->reversed++;
      return PAIR_FAILED;
    }

    must_use = 1;
  }

  ppair->must_use = must_use;

  if (sa->tar2.count > 0) {
    if (pair_spans_target(ppair, sa)) {
      ppair->target = 1;
    } else {
      ppair->target = -1;
      if (update_stats) { pair_expl->target++; }
      if (!must_use) return PAIR_FAILED;
    }
  }

  /* ============================================================= */

  /* Check that the pair is included in one of the specified ok regions */
  if ((sa->ok_regions.count > 0) && (!sa->ok_regions.any_pair)) {
    int included = 0;
    int l_start = ppair->left->start, l_end = ppair->left->start + ppair->left->length - 1;
    int r_start = ppair->right->start - ppair->right->length + 1, r_end = ppair->right->start;
    for (i=0; i<sa->ok_regions.count; i++) {
      if (sa->ok_regions.left_pairs[i][0] == -1) {
        /* any left primer, check right primer */
        if ((r_start >= sa->ok_regions.right_pairs[i][0]) &&
            (r_end <= sa->ok_regions.right_pairs[i][0] + sa->ok_regions.right_pairs[i][1] - 1)) {
          included = 1;
          break;
        }
      } else if (sa->ok_regions.right_pairs[i][0] == -1) {
        /* any right primer, check the left primer */
        if ((l_start >= sa->ok_regions.left_pairs[i][0]) && 
            (l_end <= sa->ok_regions.left_pairs[i][0] + sa->ok_regions.left_pairs[i][1] - 1)) {
          included = 1;
          break;
        }
      } else {
        /* check both primers */
        if ((l_start >= sa->ok_regions.left_pairs[i][0]) && 
            (l_end <= sa->ok_regions.left_pairs[i][0] + sa->ok_regions.left_pairs[i][1] - 1) &&
            (r_start >= sa->ok_regions.right_pairs[i][0]) &&
            (r_end <= sa->ok_regions.right_pairs[i][0] + sa->ok_regions.right_pairs[i][1] - 1)) {
          included = 1;
          break;
        }  
      }
    }
    if (!included) {
      if (update_stats) { pair_expl->not_in_any_ok_region++; }
      if (!must_use) return PAIR_FAILED;
    }
  }

  /* ============================================================= */
  /* Compute product Tm and related parameters; check constraints. */

  if (ppair->right->start - ppair->left->start + 1 <= 0) {
    fprintf(stderr, "temporary");
  }
  PR_ASSERT(ppair->right->start - ppair->left->start + 1 > 0)
  ppair->product_tm
    = long_seq_tm(sa->trimmed_seq, ppair->left->start,
                  ppair->right->start - ppair->left->start + 1,
                  /* TO DO -- skewed, it would be better to not use p_args elements here */
                  pa->p_args.salt_conc,
                  pa->p_args.divalent_conc,
                  pa->p_args.dntp_conc);

  PR_ASSERT(ppair->product_tm != OLIGOTM_ERROR);

  min_oligo_tm
    = ppair->left->temp > ppair->right->temp ? ppair->right->temp : ppair->left->temp;
  ppair->product_tm_oligo_tm_diff = ppair->product_tm - min_oligo_tm;
  ppair->t_opt_a  = 0.3 * min_oligo_tm + 0.7 * ppair->product_tm - 14.9;

  if (pa->product_min_tm != PR_DEFAULT_PRODUCT_MIN_TM
      && ppair->product_tm < pa->product_min_tm) {
    if (update_stats) { pair_expl->low_tm++; }
    if (!must_use) return PAIR_FAILED;
  }

  if (pa->product_max_tm != PR_DEFAULT_PRODUCT_MAX_TM
      && ppair->product_tm > pa->product_max_tm) {
    if (update_stats) { pair_expl->high_tm++; }
    if (!must_use) return PAIR_FAILED;
  }

  ppair->diff_tm = fabs(retval->fwd.oligo[m].temp - retval->rev.oligo[n].temp);
  if (ppair->diff_tm > pa->max_diff_tm) {
    if (update_stats) { pair_expl->temp_diff++; }
    if (!must_use) return PAIR_FAILED;
  }

  /* End of product-temperature related computations. */
  /* ============================================================= */


  /* ============================================================= */
  /* BEGIN "EXPENSIVE" computations on _individual_ primers
     in the pair.  These have been postponed until
     this point in the interests of efficiency.
  */

  /* ============================================================= */
  /* Estimate secondary structure and primer-dimer of _individual_
     primers if not already calculated. */

  /* s1 is the forward oligo. */
  _pr_substr(sa->trimmed_seq,
             retval->fwd.oligo[m].start,
             retval->fwd.oligo[m].length,
             s1);

  /* s2 is the reverse oligo. */
  _pr_substr(sa->trimmed_seq,
             retval->rev.oligo[n].start - retval->rev.oligo[n].length + 1,
             retval->rev.oligo[n].length,
             s2);

  p3_reverse_complement(s1, s1_rev);
  p3_reverse_complement(s2, s2_rev);


  if (retval->fwd.oligo[m].self_any == ALIGN_SCORE_UNDEF 
      && pa->thermodynamic_oligo_alignment==0) {
    /* We have not yet computed the 'self_any' paramter,
       which is an estimate of self primer-dimer and secondary
       structure propensity. */
    oligo_compl(&retval->fwd.oligo[m], &pa->p_args,
                &retval->fwd.expl, dpal_arg_to_use, s1, s1_rev);

    if (!OK_OR_MUST_USE(&retval->fwd.oligo[m])) {
      pair_expl->considered--;
      if (!must_use) return PAIR_FAILED;
    }
  }
   /* Thermodynamic approach, fwd-primer */
   if (retval->fwd.oligo[m].self_any == ALIGN_SCORE_UNDEF 
       && pa->thermodynamic_oligo_alignment==1) {
      oligo_compl_thermod(&retval->fwd.oligo[m], &pa->p_args,
                          &retval->fwd.expl, thal_arg_to_use, s1, s1); /* ! s1, s1_rev */
      
      if (!OK_OR_MUST_USE(&retval->fwd.oligo[m])) {
         pair_expl->considered--;
         if (!must_use) return PAIR_FAILED;
      }
   }   
   if (retval->fwd.oligo[m].hairpin_th == ALIGN_SCORE_UNDEF 
       && pa->thermodynamic_oligo_alignment==1) {
      oligo_hairpin(&retval->fwd.oligo[m], &pa->p_args,
                          &retval->fwd.expl, thal_arg_to_use, s1);
      if (!OK_OR_MUST_USE(&retval->fwd.oligo[m])) {
         pair_expl->considered--;
         if (!must_use) return PAIR_FAILED;
      }
   }
   
   /* End of thermodynamic approach */
   
  if (retval->rev.oligo[n].self_any == ALIGN_SCORE_UNDEF && pa->thermodynamic_oligo_alignment==0) {
    oligo_compl(&retval->rev.oligo[n], &pa->p_args,
                &retval->rev.expl, dpal_arg_to_use, s2_rev, s2);

    if (!OK_OR_MUST_USE(&retval->rev.oligo[n])) {
      pair_expl->considered--;
      if (!must_use) return PAIR_FAILED;
    }
  }
   /* Thermodynamic approach */
   if (retval->rev.oligo[n].self_any == ALIGN_SCORE_UNDEF && pa->thermodynamic_oligo_alignment==1) {
      oligo_compl_thermod(&retval->rev.oligo[n], &pa->p_args,
                          &retval->rev.expl, thal_arg_to_use, s2_rev, s2_rev); /* s2_rev, s2 */
      
      if (!OK_OR_MUST_USE(&retval->rev.oligo[n])) {
         pair_expl->considered--;
         if (!must_use) return PAIR_FAILED;
      }  
   }
   if (retval->rev.oligo[n].hairpin_th == ALIGN_SCORE_UNDEF && pa->thermodynamic_oligo_alignment==1) {
      oligo_hairpin(&retval->rev.oligo[n], &pa->p_args,
                    &retval->rev.expl, thal_arg_to_use, s2_rev);
      if (!OK_OR_MUST_USE(&retval->rev.oligo[n])) {
         pair_expl->considered--;
         if (!must_use) return PAIR_FAILED;
      }
   }
   
   /* End of thermodynamic approach */
  /* End of secondary structure and primer-dimer of _individual_ 
     primers. */
  /* ============================================================= */


  /* ============================================================= */
  /* Mispriming of _individual_ primers to template and 
     to repeat libraries. */

  if (retval->fwd.oligo[m].repeat_sim.score == NULL) {
    /* We have not yet checked the oligo against the repeat library. */
    oligo_repeat_library_mispriming(&retval->fwd.oligo[m], pa, sa, OT_LEFT,
                                    &retval->fwd.expl,dpal_arg_to_use);
    if (OK_OR_MUST_USE(&retval->fwd.oligo[m])) {
      oligo_template_mispriming(&retval->fwd.oligo[m], pa, sa, OT_LEFT,
				&retval->fwd.expl,
				dpal_arg_to_use->local_end,
				thal_args_for_template_mispriming);
    }
    if (!OK_OR_MUST_USE(&retval->fwd.oligo[m])) {
      pair_expl->considered--;
      if (!must_use) return PAIR_FAILED;
    }
  }
   
  if (retval->rev.oligo[n].repeat_sim.score == NULL) {
    oligo_repeat_library_mispriming(&retval->rev.oligo[n], pa, sa, OT_RIGHT,
                                    &retval->rev.expl, dpal_arg_to_use);
    if (OK_OR_MUST_USE(&retval->rev.oligo[n])) {
      oligo_template_mispriming(&retval->rev.oligo[n], pa, sa, OT_RIGHT,
				&retval->rev.expl, 
				dpal_arg_to_use->local_end,
				thal_args_for_template_mispriming);
    }
    if (!OK_OR_MUST_USE(&retval->rev.oligo[n])) {
      pair_expl->considered--;
      if (!must_use) return PAIR_FAILED;
    }
  }
   
  /* End of mispriming of _individual_ primers to template and
     to repeat libraries. */
  /* ============================================================= */


  /* ============================================================= */
  /*
   * Similarity between s1 and s2 is equivalent to complementarity between
   * s2's complement and s1.  (Both s1 and s2 are taken from the same strand.)
   */
   if(pa->thermodynamic_oligo_alignment==0) {
      ppair->compl_any = align(s1, s2, dpal_arg_to_use->local);
      if (ppair->compl_any > pa->pair_compl_any) {
         if (update_stats) { pair_expl->compl_any++; }
         if (!must_use) return PAIR_FAILED;
      }
      
      ppair->compl_end = align(s1, s2, dpal_arg_to_use->end);
      if (ppair->compl_end > pa->pair_compl_end) {
         if (update_stats) { pair_expl->compl_end++; }
         if (!must_use) return PAIR_FAILED;
      }
   } else {
      /* thermodynamical approach */
      ppair->compl_any = align_thermod(s1, s2_rev, thal_arg_to_use->any);
      if (ppair->compl_any > pa->pair_compl_any_th) {
         if (update_stats) {
            pair_expl->compl_any++; 
         }
         if (!must_use) return PAIR_FAILED;
      }
      ppair->compl_end = align_thermod(s1, s2_rev, thal_arg_to_use->end1);
      compl_end        = align_thermod(s1, s2_rev, thal_arg_to_use->end2); /* Triinu Please check */
      if (ppair->compl_end < compl_end) {
         ppair->compl_end = compl_end;
      }
      if (ppair->compl_end > pa->pair_compl_end_th) {
         if (update_stats) {
            pair_expl->compl_end++; 
         }
         if (!must_use) return PAIR_FAILED;
      }
   }
   
  /*
   * It is conceivable (though unlikely) that
   * align(s2_rev, s1_rev, end_args) > align(s1,s2,end_args).
   */
  if (pa->thermodynamic_oligo_alignment==0 && (compl_end = align(s2_rev, s1_rev, dpal_arg_to_use->end))
      > ppair->compl_end) {
    if (compl_end > pa->p_args.max_self_end) {
      if (update_stats) { pair_expl->compl_end++; }
      if (!must_use) return PAIR_FAILED;
    }
    ppair->compl_end = compl_end;
  }

  if ((ppair->repeat_sim = pair_repeat_sim(ppair, pa))
      > pa->pair_repeat_compl) {
    if (update_stats) { pair_expl->repeat_sim++; }
    if (!must_use) return PAIR_FAILED;
  }
   /* thermodynamic approach */
   if (pa->thermodynamic_oligo_alignment==1 && 
       ((compl_end = align_thermod(s2, s1_rev, thal_arg_to_use->end1)) > ppair->compl_end || 
	(compl_end = align_thermod(s2, s1_rev, thal_arg_to_use->end2)) > ppair->compl_end)) {
      if (compl_end > pa->p_args.max_self_end_th) {
         if (update_stats) {
            pair_expl->compl_end++; 
         }
         if (!must_use) return PAIR_FAILED;
      }
      ppair->compl_end = compl_end;
   }
   
  /* ============================================================= */


  /* ============================================================= */
  /* Calculate _pair_ mispriming, if necessary. */
 
   if (pa->thermodynamic_template_alignment == 0) {
     if (!_pr_need_pair_template_mispriming(pa))
        ppair->template_mispriming = ALIGN_SCORE_UNDEF;
     else {
        PR_ASSERT(ppair->left->template_mispriming != ALIGN_SCORE_UNDEF);
        PR_ASSERT(ppair->left->template_mispriming_r != ALIGN_SCORE_UNDEF);
        PR_ASSERT(ppair->right->template_mispriming != ALIGN_SCORE_UNDEF);
        PR_ASSERT(ppair->right->template_mispriming_r != ALIGN_SCORE_UNDEF);
        ppair->template_mispriming =
           ppair->left->template_mispriming + ppair->right->template_mispriming_r;
        if ((ppair->left->template_mispriming_r + ppair->right->template_mispriming)
             > ppair->template_mispriming)
           ppair->template_mispriming
           = ppair->left->template_mispriming_r + ppair->right->template_mispriming;

        if (pa->pair_max_template_mispriming >= 0.0
             && ppair->template_mispriming > pa->pair_max_template_mispriming) {
            if (update_stats) { pair_expl->template_mispriming++; }
            if (!must_use) return PAIR_FAILED;
        }
     }
   } else { /* thermodynamic approach */
     if (!_pr_need_pair_template_mispriming_thermod(pa))
       ppair->template_mispriming = ALIGN_SCORE_UNDEF;
     else {
       PR_ASSERT(ppair->left->template_mispriming != ALIGN_SCORE_UNDEF);
       PR_ASSERT(ppair->left->template_mispriming_r != ALIGN_SCORE_UNDEF);
       PR_ASSERT(ppair->right->template_mispriming != ALIGN_SCORE_UNDEF);
       PR_ASSERT(ppair->right->template_mispriming_r != ALIGN_SCORE_UNDEF);
       ppair->template_mispriming =
           ppair->left->template_mispriming + ppair->right->template_mispriming_r;
       if ((ppair->left->template_mispriming_r + ppair->right->template_mispriming)
             > ppair->template_mispriming)
           ppair->template_mispriming
           = ppair->left->template_mispriming_r + ppair->right->template_mispriming;

       if (pa->pair_max_template_mispriming_th
           && ppair->template_mispriming > pa->pair_max_template_mispriming_th) {
            if (update_stats) {
               pair_expl->template_mispriming++;
            }
            if (!must_use) return PAIR_FAILED;
       }
     }
   }
      
   /* End of calculating _pair_ mispriming if necessary. */
   /* ============================================================= */

  return PAIR_OK;
} /* characterize_pair */


/*
 */
void
compute_position_penalty(const p3_global_settings *pa,
                         const seq_args *sa,
                         primer_rec *h,
                         oligo_type o_type)
{
  int three_prime_base;
  int inside_flag = 0;
  int target_begin, target_end;

  PR_ASSERT(OT_LEFT == o_type || OT_RIGHT == o_type);
  PR_ASSERT(1 == sa->tar2.count);
  target_begin = sa->tar2.pairs[0][0];
  target_end = target_begin + sa->tar2.pairs[0][1] - 1;

  three_prime_base = OT_LEFT == o_type
    ? h->start + h->length - 1 : h->start - h->length + 1;
  bf_set_infinite_pos_penalty(h,1);
  h->position_penalty = 0.0;

  if (OT_LEFT == o_type) {
    if (three_prime_base <= target_end) {
      bf_set_infinite_pos_penalty(h,0);
      if (three_prime_base < target_begin)
        h->position_penalty = target_begin - three_prime_base - 1;
      else {
        h->position_penalty = three_prime_base - target_begin + 1;
        inside_flag = 1;
      }
    }
  } else { /* OT_RIGHT == o_type */
    if (three_prime_base >= target_begin) {
      bf_set_infinite_pos_penalty(h,0);
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
pair_spans_target(const primer_pair *pair, const seq_args *sa)
{
  int i;
  int last_of_left = pair->left->start + pair->left->length - 1;
  int first_of_right = pair->right->start - pair->right->length + 1;
  int target_first, target_last;
  for (i = 0; i < sa->tar2.count; i++) {
    target_first = sa->tar2.pairs[i][0];
    target_last = target_first + sa->tar2.pairs[i][1] - 1;
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
obj_fn(const p3_global_settings *pa, primer_pair *h)
{
  double sum;
  double lower_tm;

  sum = 0.0;
  lower_tm = h->right->temp;
  if(h->left->temp < h->right->temp) lower_tm = h->left->temp;
     
  if(pa->pr_pair_weights.primer_quality)   /*  HERE 1 */
    sum += pa->pr_pair_weights.primer_quality * (h->left->quality + h->right->quality);

  if(pa->pr_pair_weights.io_quality && pa->pick_right_primer
     && pa->pick_left_primer && pa->pick_internal_oligo)
    sum += pa->pr_pair_weights.io_quality * h->intl->quality;

  if(pa->pr_pair_weights.diff_tm)
    sum += pa->pr_pair_weights.diff_tm * h->diff_tm;

  if (pa->thermodynamic_oligo_alignment==0) {
    if (pa->pr_pair_weights.compl_any)
      sum += pa->pr_pair_weights.compl_any * h->compl_any;
    if (pa->pr_pair_weights.compl_end)
      sum += pa->pr_pair_weights.compl_end * h->compl_end;
  } else if (pa->thermodynamic_oligo_alignment==1) {

    if (pa->pr_pair_weights.compl_any_th && ((lower_tm - pa->pr_pair_weights.temp_cutoff) <= h->compl_any))
      sum += pa->pr_pair_weights.compl_any_th * (h->compl_any - (lower_tm - pa->pr_pair_weights.temp_cutoff - 1.0));
    if (pa->pr_pair_weights.compl_any_th && ((lower_tm - pa->pr_pair_weights.temp_cutoff) > h->compl_any))
      sum += pa->pr_pair_weights.compl_any_th * (1/(lower_tm - pa->pr_pair_weights.temp_cutoff + 1.0 - h->compl_any));
   
    if (pa->pr_pair_weights.compl_end_th && ((lower_tm - pa->pr_pair_weights.temp_cutoff) <= h->compl_end))
      sum += pa->pr_pair_weights.compl_end_th * (h->compl_end - (lower_tm - pa->pr_pair_weights.temp_cutoff - 1.0));
    if (pa->pr_pair_weights.compl_end_th && ((lower_tm - pa->pr_pair_weights.temp_cutoff) > h->compl_end))
      sum += pa->pr_pair_weights.compl_end_th * (1/(lower_tm - pa->pr_pair_weights.temp_cutoff + 1.0 - h->compl_end));
  } else {
    PR_ASSERT(0);
    /* Programming error */
  }

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

  if (pa->pr_pair_weights.template_mispriming && pa->thermodynamic_template_alignment==0) {
    PR_ASSERT(pa->pr_pair_weights.template_mispriming >= 0.0);
    PR_ASSERT(h->template_mispriming >= 0.0);
    sum += pa->pr_pair_weights.template_mispriming * h->template_mispriming;
  }
  if (pa->pr_pair_weights.template_mispriming_th && pa->thermodynamic_template_alignment==1) {
      PR_ASSERT(pa->pr_pair_weights.template_mispriming_th >= 0.0);
      PR_ASSERT(h->template_mispriming >= 0.0);
      if((lower_tm - pa->pr_pair_weights.temp_cutoff) <= h->template_mispriming)
        sum += pa->pr_pair_weights.template_mispriming_th * 
        (h->template_mispriming - (lower_tm - pa->pr_pair_weights.temp_cutoff - 1.0));
      if((lower_tm - pa->pr_pair_weights.temp_cutoff) > h->template_mispriming)
        sum += pa->pr_pair_weights.template_mispriming_th *
        (1/(lower_tm - pa->pr_pair_weights.temp_cutoff + 1.0 - h->template_mispriming));
   }
   PR_ASSERT(sum >= 0.0);
   return sum;
}

static double
align(const char *s1,
      const char *s2,
      const dpal_args *a) {
   dpal_results r;
   if(a->flag == DPAL_LOCAL || a->flag == DPAL_LOCAL_END) {
    if (strlen(s2) < 3) {
      /* For extremely short alignments we simply
         max out the score, because the dpal subroutines
         for these cannot handle this case.
         TO DO: this can probably be corrected in dpal. */
      return strlen(s2);
    }
  }
  dpal((const unsigned char *) s1, (const unsigned char *) s2, a, &r);
  PR_ASSERT(r.score <= SHRT_MAX);
  if (r.score == DPAL_ERROR_SCORE) {
    /* There was an error. */
    if (errno == ENOMEM) {
      longjmp(_jmp_buf, 1);
    } else {
      /* This branch is taken only if there is a programming error, in
         that s1 or s2 were NULL or contained an illegal character. We
         try to print some debugging information before aborting. */
      fprintf(stderr, "%s", r.msg);
      /* Always false, causes an abort: */
      PR_ASSERT(r.score != DPAL_ERROR_SCORE);
    }
  }
  return ((r.score < 0.0) ? 0.0 : r.score / PR_ALIGN_SCORE_PRECISION);
}

/* Helper function */
static int
_set_string(char **loc, const char *new_string) {
  if (*loc) {
    free(*loc);
  }
  if (!(*loc = (char *) malloc(strlen(new_string) + 1)))
    return 1; /* ENOMEM */
  strcpy(*loc, new_string);
  return 0;
}

static double
align_thermod(const char *s1,
              const char *s2,
              const thal_args *a) 
{  
   int thal_trace=0;
   thal_results r;
   thal((const unsigned char *) s1, (const unsigned char *) s2, a, &r);
   if (thal_trace) {
     fprintf(stdout, 
	     "thal, thal_args, type=%d maxLoop=%d mv=%f dv=%f "
	     "dntp=%f dna_conc=%f, temp=%f, temponly=%d dimer=%d\n",
	     a->type, a->maxLoop, a->mv, a->dv, a->dntp, a->dna_conc, 
	     a->temp, a->temponly, a->dimer);
     fprintf(stdout, 
	     "thal: s1=%s s2=%s temp=%f msg=%s end1=%d end2=%d\n", 
	     s1, s2, r.temp, r.msg, r.align_end_1, r.align_end_2);
   }
   PR_ASSERT(r.temp <= DBL_MAX);
   if (r.temp == THAL_ERROR_SCORE) 
     {
        /* There was an error. */
        if (errno == ENOMEM) {
	  longjmp(_jmp_buf, 1);
	} else {
	  /* This branch is taken if there is an error other than
	     ENOMEM, in which case we treat it as a per sequence error.
	  */
	  /* Save the error string */
	  if (!_set_string(&thermodynamic_alignment_length_error_msg, r.msg)) {
	    thermodynamic_alignment_length_error = 1;
	  }
	  longjmp(_jmp_buf, 1);
	}
     }
   return ((r.temp < 0.0) ? 0.0 : (double) (r.temp));
}

/* Return the sequence of oligo in
   a ****static***** buffer.  The
   sequence returned is changed at the
   next call to pr_oligo_sequence. */
char *
pr_oligo_sequence(const seq_args *sa,
    const primer_rec *oligo)
{
  static char s[MAX_PRIMER_LENGTH+1];
  int seq_len;
  PR_ASSERT(NULL != sa);
  PR_ASSERT(NULL != oligo);
  seq_len = strlen(sa->sequence);
  PR_ASSERT(oligo->start + sa->incl_s >= 0);
  PR_ASSERT(oligo->start + sa->incl_s + oligo->length <= seq_len);
  _pr_substr(sa->sequence, sa->incl_s + oligo->start, oligo->length, s);
  return &s[0];
}

char *
pr_oligo_rev_c_sequence(const seq_args *sa,
    const primer_rec *o)
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
  p3_reverse_complement(s,s1);
  return &s1[0];
}

/* Calculate self complementarity (both h->self_any and h->self_end)
   which we use as an approximation for both secondary structure and
   self primer-dimer. */
static void
oligo_compl(primer_rec *h,
              const args_for_one_oligo_or_primer *po_args,
              oligo_stats *ostats,
              const dpal_arg_holder *dpal_arg_to_use,
              const char *oligo_seq,
              const char *revc_oligo_seq)
{
   PR_ASSERT(h != NULL);
   
   h->self_any = align(oligo_seq, revc_oligo_seq, dpal_arg_to_use->local);
   if (h->self_any > po_args->max_self_any) {
      op_set_high_self_any(h);
      ostats->compl_any++;
      ostats->ok--;     
      if (!h->must_use) return;
   }
   
   h->self_end = align(oligo_seq, revc_oligo_seq, dpal_arg_to_use->end);
   if (h->self_end > po_args->max_self_end) {
      op_set_high_self_end(h);
      ostats->compl_end++;
      ostats->ok--;
      if (!h->must_use) return;
   }
}

static void
oligo_compl_thermod(primer_rec *h,
                    const args_for_one_oligo_or_primer *po_args,
                    oligo_stats *ostats,
                    const thal_arg_holder *thal_arg_to_use,
                    const char *oligo_seq,
                    const char *revc_oligo_seq)
{
   
   PR_ASSERT(h != NULL);
   h->self_any = align_thermod(oligo_seq, revc_oligo_seq, thal_arg_to_use->any);
   if(h->self_any > po_args->max_self_any_th) { 
      op_set_high_self_any(h);
      ostats->compl_any++;
      ostats->ok--;
      if (!h->must_use) return;
   }  
   h->self_end = align_thermod(oligo_seq, revc_oligo_seq, thal_arg_to_use->end1);  
   if(h->self_end > po_args->max_self_end_th){  
      op_set_high_self_end(h);
      ostats->compl_end++;
      ostats->ok--;
      if (!h->must_use) return;
   }   
}

static void 
oligo_hairpin(primer_rec *h,
              const args_for_one_oligo_or_primer *po_args,
              oligo_stats *ostats,
              const thal_arg_holder *thal_arg_to_use,
              const char *oligo_seq)
{
   PR_ASSERT(h != NULL);
   h->hairpin_th = align_thermod(oligo_seq, oligo_seq, thal_arg_to_use->hairpin_th);
   if(h->hairpin_th > po_args->max_hairpin_th) {
      op_set_high_hairpin_th(h);
      ostats->hairpin_th++;
      ostats->ok--;
      return;
   }
}

static void
primer_mispriming_to_template(primer_rec *h,
                              const p3_global_settings *pa,
                              const seq_args *sa,
                              oligo_type l,
                              oligo_stats *ostats,
                              int first,
                              int last,
                              /* The oligo sequence: */
                              const char *s,
                              /* s reverse complemented: */
                              const char *s_r,
                              const dpal_args *align_args
                              )
{
  const char *oseq;
  char *target, *target_r;
  int tmp, seqlen;
  int debug = 0;
  int first_untrimmed, last_untrimmed;
                  /* Indexes of first and last bases of the oligo in sa->seq,
                     that is, WITHIN THE TOTAL SEQUENCE INPUT. */

  /* first, last are indexes of first and last bases of the oligo in
     sa->trimmed_seq, that is, WITHIN THE INCLUDED REGION. */

  char   tmp_char;
  double tmp_score;

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
    fprintf(stderr, "first_untrimmed = %d, last_untrimmed = %d,first = %d, last = %d\n",
            first_untrimmed, last_untrimmed, first, last);

    fprintf(stderr, "5' of oligo: Score %f aligning %s against %s\n\n", tmp_score,
            oseq, target);
  }

  target[first_untrimmed] = tmp_char;

  /* 2. Align to the template 3' of the oligo. */
  h->template_mispriming
    = align(oseq, &target[0] + last_untrimmed + 1, align_args);

  if (debug)
    fprintf(stderr, "3' of oligo Score %f aligning %s against %s\n\n",
            h->template_mispriming, oseq, &target[0] + last_untrimmed + 1);

  /* 3. Take the max of 1. and 2. */
  if (tmp_score > h->template_mispriming)
    h->template_mispriming = tmp_score;

  /* 4. Align to the reverse strand of the template. */
  h->template_mispriming_r
    = align(oseq, target_r, align_args);

  if (debug)
    fprintf(stderr, "other strand Score %f aligning %s against %s\n\n",
            h->template_mispriming_r, oseq, target_r);

  if (pa->p_args.max_template_mispriming >= 0) {
    if (oligo_max_template_mispriming(h)
	> pa->p_args.max_template_mispriming) {
      op_set_high_similarity_to_multiple_template_sites(h);
      if (OT_LEFT == l || OT_RIGHT == l ) {
        ostats->template_mispriming++;
        ostats->ok--;
      } else PR_ASSERT(0); /* Should not get here. */
    } else {
      /* This oligo is ok, mark it so we do not do this check again. */
      h->template_mispriming_ok = 1;
    }
  }
}

static void
primer_mispriming_to_template_thermod(primer_rec *h,
                                      const p3_global_settings *pa,
                                      const seq_args *sa,
                                      oligo_type l,
                                      oligo_stats *ostats,
                                      int first,
                                      int last,
                                      /* The oligo sequence: */
                                      const char *s,
                                      /* s reverse complemented: */
                                      const char *s_r,
                                      const thal_args *align_args)
{   
   const char *oseq;
   char *target, *target_r;
   int tmp, seqlen;
   int debug = 0;
   int first_untrimmed, last_untrimmed;
   /* Indexes of first and last bases of the oligo in sa->seq,
    *                      that is, WITHIN THE TOTAL SEQUENCE INPUT. */
   
   /* first, last are indexes of first and last bases of the oligo in
    *      sa->trimmed_seq, that is, WITHIN THE INCLUDED REGION. */
   
   char   tmp_char;
   double  tmp_score;
   
   seqlen = strlen(sa->upcased_seq);
   first_untrimmed = sa->incl_s + first;
   last_untrimmed = sa->incl_s + last;
   
   if (l == OT_RIGHT) {
      oseq = &s_r[0];
      target = &sa->upcased_seq[0];
      target_r = &sa->upcased_seq_r[0];
   } else { /* l == OT_RIGHT */
      if (debug)
        fprintf(stdout, "first_untrimmed = %d, last_untrimmed = %d\n",
                first_untrimmed, last_untrimmed);
      oseq = &s[0];
      target = &sa->upcased_seq_r[0];
      target_r = &sa->upcased_seq[0];
      /* We need to adjust first_untrimmed and last_untrimmed so that
       * they are correct in the reverse-complemented
       * sequence.
             *     */
            tmp = (seqlen - last_untrimmed) - 1;
            last_untrimmed  = (seqlen - first_untrimmed) - 1;
            first_untrimmed = tmp;
   }
   
   /* 1. Align to the template 5' of the oligo. */
   /* tmp_char = target_r[first_untrimmed]; Corrected 2012-05-30 */
   tmp_char = target[first_untrimmed];
   target[first_untrimmed] = '\0';
   
   tmp_score = align_thermod(oseq, target, align_args);
     
   if (debug) {
        if (l == OT_LEFT) fprintf(stdout, "\n************ OLIGO = LEFT\n");
        else fprintf(stdout,              "\n************ OLIGO = RIGHT\n");
        fprintf(stdout, "first_untrimmed = %d, last_untrimmed = %d, first = %d, last = %d\n",
                first_untrimmed, last_untrimmed, first, last);
        
        fprintf(stdout, "5' of oligo: Score %f aligning %s against %s\n\n", tmp_score,
                oseq, target);
   }
   target[first_untrimmed] = tmp_char;
   
   /* 2. Align to the template 3' of the oligo. */
   h->template_mispriming
     = align_thermod(oseq, &target[0] + last_untrimmed + 1, align_args);
   
   if (debug)
     fprintf(stdout, "3' of oligo Score %f aligning %s against %s\n\n",
	     h->template_mispriming, oseq, &target[0] + last_untrimmed + 1);
   
   /* 3. Take the max of 1. and 2. */
   if (tmp_score > h->template_mispriming)
     h->template_mispriming = tmp_score;
   
   /* 4. Align to the reverse strand of the template. */
   h->template_mispriming_r
     = align_thermod(oseq, target_r, align_args);
   
   if (debug)
     fprintf(stdout,
	     "In primer_mispriming_to_template_thermod\n"
	     " other strand Score %f aligning %s against %s\n\n",
	     h->template_mispriming_r, oseq, target_r);

   if (pa->p_args.max_template_mispriming_th >= 0) {
     if (oligo_max_template_mispriming_thermod(h)
	 > pa->p_args.max_template_mispriming_th) {
       op_set_high_similarity_to_multiple_template_sites(h);
       if (OT_LEFT == l || OT_RIGHT == l ) {
	 ostats->template_mispriming++;
	 ostats->ok--;
       }
       else PR_ASSERT(0); /* Should not get here. */
     } else {
       /* This oligo is ok, mark it so we do not do this check again. */
       h->template_mispriming_ok = 1;
     }
   }
}

static void
oligo_compute_sequence_and_reverse(primer_rec *h,
                                   const seq_args *sa,
                                   oligo_type l,  
                                   int *first, int *last,
                                   char *s, char *s_r)
{
  *first =  (OT_LEFT == l || OT_INTL == l)
    ? h->start
    : h->start - h->length + 1;
  *last  =  (OT_LEFT == l || OT_INTL == l)
    ? h->start + h->length - 1
    : h->start;

  _pr_substr(sa->trimmed_seq, *first, h->length, s);
  p3_reverse_complement(s, s_r);
}

/* Possible improvement -- pass in the oligo sequences */
static void
oligo_repeat_library_mispriming(primer_rec *h,
                                const p3_global_settings *pa,
                                const seq_args *sa,
                                oligo_type l,
                                oligo_stats *ostats,
                                const dpal_arg_holder *dpal_arg_to_use)
{
  char
    s[MAX_PRIMER_LENGTH+1],     /* Will contain the oligo sequence. */
    s_r[MAX_PRIMER_LENGTH+1];   /* Will contain s reverse complemented. */

  double w;
  const seq_lib *lib;
  int i;
  int first, last; /* Indexes of first and last bases of the oligo in
                     sa->trimmed_seq, that is, WITHIN THE INCLUDED
                     REGION. */
  int   min, max;
  short max_lib_compl;

  /* First, check the oligo against the repeat library. */
  if (OT_INTL == l) {
    lib = pa->o_args.repeat_lib;
    max_lib_compl = (short) pa->o_args.max_repeat_compl;
  } else {
    lib = pa->p_args.repeat_lib;
    max_lib_compl = (short) pa->p_args.max_repeat_compl;
  }

  oligo_compute_sequence_and_reverse(h, sa, l, &first, &last, s, s_r);

  /*
   * Calculate maximum similarity to sequences from user defined repeat
   * library. Compare it with maximum allowed repeat similarity.
   */

  if (seq_lib_num_seq(lib) > 0) {
    /* Library exists and is non-empty. */

    h->repeat_sim.score =
      (double *) pr_safe_malloc(lib->seq_num * sizeof(double));
    h->repeat_sim.max = h->repeat_sim.min = 0;
    max = min = 0;
    h->repeat_sim.name = lib->names[0];

    for (i = 0; i < lib->seq_num; i++) {
      if (OT_LEFT == l)
        w = lib->weight[i] *
          align(s, lib->seqs[i],
                (pa->lib_ambiguity_codes_consensus
                 ? dpal_arg_to_use->local_end_ambig
                 : dpal_arg_to_use->local_end));

      else if (OT_INTL == l)
        w = lib->weight[i] *
          align(s, lib->seqs[i],
                (pa->lib_ambiguity_codes_consensus
                 ? dpal_arg_to_use->local_ambig
                 : dpal_arg_to_use->local));

      else
        w = lib->weight[i] *
          align(s_r, lib->rev_compl_seqs[i],
                (pa->lib_ambiguity_codes_consensus
                 ? dpal_arg_to_use->local_end_ambig
                 : dpal_arg_to_use->local));


      if (w > SHRT_MAX || w < SHRT_MIN) {
        abort(); /* TO DO, propagate error */
        /* This check is necessary for the next 9 lines */
      }
      h->repeat_sim.score[i] = w;
      if(w > max){
        max = (int) w;
        h->repeat_sim.max = i;
        h->repeat_sim.name = lib->names[i];
      }
      if(w < min){
        min = (int) w;
        h->repeat_sim.min = i;
      }

      if (w > max_lib_compl) {
        op_set_high_similarity_to_non_template_seq(h);
        ostats->repeat_score++;
        ostats->ok--;
        if (!h->must_use) return;
      } /* if w > max_lib_compl */
    } /* for */
  } /* if library exists and is non-empty */
  /* End of checking against the repeat library */  
}

/* This function carries out either the old "dpal" alignment or the
   thermodynamic alignment, depending on the value of
   pa->thermodynamic_template_alignment. */
static void
oligo_template_mispriming(primer_rec *h,
                          const p3_global_settings *pa,
                          const seq_args *sa,
                          oligo_type l,
                          oligo_stats *ostats,
                          const dpal_args *d_align_args,
                          const thal_args *t_align_args)
{
  char
    s[MAX_PRIMER_LENGTH+1],     /* Will contain the oligo sequence. */
    s_r[MAX_PRIMER_LENGTH+1];   /* Will contain s reverse complemented. */

  int first, last; /* Indexes of first and last bases of the oligo in
                     sa->trimmed_seq, that is, WITHIN THE INCLUDED
                     REGION. */
  
  /* Check if we already did this and the oligo was ok. */
  if (h->template_mispriming_ok) {
    return;
  }

  oligo_compute_sequence_and_reverse(h, sa, l, &first, &last, s, s_r);

  /* Calculate maximum similarity to ectopic sites in the template. */
  if (l == OT_RIGHT || l == OT_LEFT) {
    if (pa->thermodynamic_template_alignment == 0 && _pr_need_template_mispriming(pa))
      primer_mispriming_to_template(h, pa, sa, l,
                                    ostats, first,
                                    last, s, s_r, d_align_args);
    if (pa->thermodynamic_template_alignment == 1 && _pr_need_template_mispriming_thermod(pa))
      primer_mispriming_to_template_thermod(h, pa, sa, l,
                                            ostats, first,
                                            last, s, s_r, t_align_args);
  }
}

static int
pair_repeat_sim(primer_pair *h,
                const p3_global_settings *pa) {
  int i, n, max, w;
  primer_rec *fw, *rev;

  fw = h->left;
  rev = h->right;

  max = 0;
  n = seq_lib_num_seq(pa->p_args.repeat_lib);
  if(n == 0) return 0;
  h->rep_name =  pa->p_args.repeat_lib->names[0] ;
  for (i = 0; i < n; i++) {
    if ((w = (int) (fw->repeat_sim.score[i] +
                    rev->repeat_sim.score[i])) > max) {
      max = w;
      h->rep_name =  pa->p_args.repeat_lib->names[i] ;
    }
  }
  return max;
}

static void
set_retval_both_stop_codons(const seq_args *sa, p3retval *retval) {
  /* The position of the intial base of the rightmost stop codon that
   * is to the left of sa->start_codon_pos; valid only if
   * sa->start_codon_pos is "not null".  We will not want to include
   * a stop codon to the right of the the start codon in the
   * amplicon. */
  retval->upstream_stop_codon = find_stop_codon(sa->trimmed_seq,
                                                sa->start_codon_pos, -1);
  retval->upstream_stop_codon += sa->incl_s;
  retval->stop_codon_pos = find_stop_codon(sa->trimmed_seq,
                                             sa->start_codon_pos,  1);
  retval->stop_codon_pos += sa->incl_s;
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
static int
find_stop_codon(const char* s,
                int start,
                int direction) {
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
strcmp_nocase(const char *s1, const char *s2)
{
  static char M[UCHAR_MAX];
  static int f = 0;
  int i;
  const char *p, *q;

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
free_repeat_sim_score(p3retval *state)
{
  int i;

  for (i = 0; i < state->fwd.num_elem; i++) {
    if (state->fwd.oligo[i].repeat_sim.score != NULL) {
      free(state->fwd.oligo[i].repeat_sim.score);
      state->fwd.oligo[i].repeat_sim.score = NULL;
    }
  }

  for (i = 0; i < state->rev.num_elem; i++) {
    if (state->rev.oligo[i].repeat_sim.score != NULL) {
      free(state->rev.oligo[i].repeat_sim.score);
      state->rev.oligo[i].repeat_sim.score = NULL;
    }
  }

  for (i = 0; i < state->intl.num_elem; i++) {
    if (state->intl.oligo[i].repeat_sim.score != NULL) {
      free(state->intl.oligo[i].repeat_sim.score);
      state->intl.oligo[i].repeat_sim.score = NULL;
    }
  }
}

static void
free_primer_repeat_sim_score(primer_rec *h) 
{ 
  if (h->repeat_sim.score != NULL) { 
    free(h->repeat_sim.score); 
    h->repeat_sim.score = NULL;
  } 
}

/*
 Edited by T. Koressaar for lowercase masking. This function checks
 if the 3' end of the primer has been masked by lowercase letter.
 Function created/Added by Eric Reppo, July 9, 2002
 */
static int
is_lowercase_masked(int position,
                    const char *sequence,
                    primer_rec *h,
                    oligo_stats *global_oligo_stats)
{
  const char* p = &sequence[position];
  if ('a' == *p || 'c' == *p ||'g' == *p || 't' == *p) {
    op_set_overlaps_masked_sequence(h);
    global_oligo_stats->gmasked++;
    return 1;
  }
  return 0;
}

/*
 Edited by A. Untergasser
 */
static int
primer_must_match(const p3_global_settings *pa, primer_rec *h, oligo_stats *global_oligo_stats,
		  /* This is 5'->3' on the template sequence: */
		  const char *input_oligo_seq, char *match_three_prime, char *match_five_prime)
{
  const char *seq;
  char *test;
  int length = h->length - 5;
  if (match_five_prime != NULL) {
    seq = input_oligo_seq;
    test = match_five_prime;
    for (int i = 0; i < 5; i++) {
      if (!compare_nucleotides(*seq, *test)) {
	global_oligo_stats->must_match_fail++;
	return 1;
      }
      seq++;
      test++;
    }
  }
  if (match_three_prime != NULL) {
    seq = input_oligo_seq;
    test = match_three_prime;
    seq = seq + length;
    for (int i = 0; i < 5; i++) {
      if (!compare_nucleotides(*seq, *test)) {
	global_oligo_stats->must_match_fail++;
	return 1;
      }
      seq++;
      test++;
		}
  }
  return 0;
}

/* For a [NACTG] is allowed, for b [NACTGRYWSMKBHDV].*/
static int compare_nucleotides(const char a, const char b) 
{
  char x = a;
  char y = b;
  /* Convert to uppercase */
  if(a >= 'a' && a <= 'z'){
    x = ('A' + a - 'a');
  }
  if(b >= 'a' && b <= 'z'){
    y = ('A' + b - 'a');
  }
  
  if ( x == y ) {
    return 1;
  }
  if (( x == 'N') or (y == 'N')) {
    return 1;
  }
  if (x == 'A') {
    if ((y == 'R') or (y == 'W') or (y == 'M') or
	(y == 'H') or (y == 'D') or (y == 'V')){
      return 1;
    }
  }
  if (x == 'G') {
    if ((y == 'R') or (y == 'S') or (y == 'K') or
	(y == 'B') or (y == 'D') or (y == 'V')){
      return 1;
    }
  }
  if (x == 'C') {
    if ((y == 'Y') or (y == 'S') or (y == 'M') or
	(y == 'B') or (y == 'H') or (y == 'V')){
      return 1;
    }
  }
  if (x == 'T') {
    if ((y == 'Y') or (y == 'W') or (y == 'K') or
	(y == 'B') or (y == 'H') or (y == 'D')){
      return 1;
    }
  }
  return 0;
}

static int test_must_match_parameters(char *test) 
{
  int i = 0;
  char x;
  while (*test != '\0') {
    /* First to UPPER letters */
    if ((*test >= 'a') && (*test <= 'z')){
      x = ('A' + (*test) - 'a');
    } else {
      x = *test;
    }
    /* Second check that it is in Range A-Z */
    if ((x < 'A') || (x > 'Z')){
      return 1;
    }
    /* Check it is NACTGRYWSMKBHDV */
    if ((x == 'N') ||
	(x == 'A') || (x == 'C') ||
	(x == 'T') || (x == 'G') ||
	(x == 'R') || (x == 'Y') ||
	(x == 'W') || (x == 'S') ||
	(x == 'M') || (x == 'K') ||
	(x == 'B') || (x == 'H') ||
	(x == 'D') || (x == 'V')){
      test++;
      i++;
    } else {
      return 1;
    }
  }
  /* Check it is 5 letters */
  if (i != 5) {
    return 1;
  }
  return 0;
}

/* Put substring of seq starting at n with length m into s. */
void
_pr_substr(const char *seq, int n, int m, char *s)
{
  int i;
  for(i=n;i<n+m;i++)s[i-n]=seq[i];
  s[m]='\0';
}

/* Reverse and complement the sequence seq and put the result in s.
   WARNING: It is up the caller to ensure that s points to enough
   space. */
void
p3_reverse_complement(const char *seq, char *s)
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

int
_pr_need_template_mispriming(const p3_global_settings *pa) {
  return
    pa->p_args.max_template_mispriming >= 0
    || pa->p_args.weights.template_mispriming > 0.0
    || _pr_need_pair_template_mispriming(pa);
}
int
_pr_need_template_mispriming_thermod(const p3_global_settings *pa) {
     return
       pa->p_args.max_template_mispriming_th >= 0
       || pa->p_args.weights.template_mispriming_th > 0.0
       || _pr_need_pair_template_mispriming_thermod(pa);
  }

int
_pr_need_pair_template_mispriming(const p3_global_settings *pa)
{
  return
    pa->pair_max_template_mispriming >= 0
    || pa->pr_pair_weights.template_mispriming > 0.0;
}

int
_pr_need_pair_template_mispriming_thermod(const p3_global_settings *pa)
{   
   return
     pa->pair_max_template_mispriming_th >= 0
     || pa->pr_pair_weights.template_mispriming_th > 0.0;
}

/* Upcase a DNA string, s, in place.  If amibiguity_code_ok is false,
   then the characters acgtnACGTN are 'recognized' and upcased.  If
   ambiguity_code_ok is true, then the IUB/IUPAC ambiguity codes are
   also 'recognized' and upcased.

   Change any unrecognized letters in *s to 'N"

   Return the first unrecognized letter, if
   any, that is seen. Otherwise return '\0'. */
static char
dna_to_upper(char * s, int ambiguity_code_ok)
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

static char *
strstr_nocase(char *s1, char *s2) {
  int  n1, n2;
  char *p, q, *tmp;

  if(s1 == NULL || s2 == NULL) return NULL;
  n1 = strlen(s1); n2 = strlen(s2);
  if(n1 < n2) return NULL;

  tmp = (char *) pr_safe_malloc(n1 + 1);
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

#define CHECK  if (r > bsize || r < 0) return "Internal error, not enough space for \"explain\" string";  bufp += r; bsize -= r
#define SP_AND_CHECK(FMT, VAL) { r = sprintf(bufp, FMT, VAL); CHECK; }
#define IF_SP_AND_CHECK(FMT, VAL) { if (VAL) { SP_AND_CHECK(FMT, VAL) } }
const char *
p3_pair_explain_string(const pair_stats *pair_expl)
{
  /* WARNING WARNING WARNING

  Static, fixed-size buffer.  If you add more calls to SP_AND_CHECK or
  IF_SP_AND_CHECK, you need to make sure the buffer cannot overflow.
  */

  static char buf[10000];
  char *bufp = buf;
  size_t bsize = 10000;
  size_t r;

  SP_AND_CHECK("considered %d", pair_expl->considered)
  IF_SP_AND_CHECK(", no target %d", pair_expl->target)
  IF_SP_AND_CHECK(", unacceptable product size %d", pair_expl->product)
  IF_SP_AND_CHECK(", low product Tm %d", pair_expl->low_tm)
  IF_SP_AND_CHECK(", high product Tm %d", pair_expl->high_tm)
  IF_SP_AND_CHECK(", tm diff too large %d",pair_expl->temp_diff)
  IF_SP_AND_CHECK(", high any compl %d", pair_expl->compl_any)
  IF_SP_AND_CHECK(", high end compl %d", pair_expl->compl_end)
  IF_SP_AND_CHECK(", no internal oligo %d", pair_expl->internal)
  IF_SP_AND_CHECK(", high mispriming library similarity %d",
                  pair_expl->repeat_sim)
  IF_SP_AND_CHECK(", no overlap of required point %d",
                  pair_expl->does_not_overlap_a_required_point)
  IF_SP_AND_CHECK(", primer in pair overlaps a primer in a better pair %d",
                  pair_expl->overlaps_oligo_in_better_pair)
  IF_SP_AND_CHECK(", high template mispriming score %d",
                  pair_expl->template_mispriming)
  IF_SP_AND_CHECK(", not in any ok region %d", 
                  pair_expl->not_in_any_ok_region);
  IF_SP_AND_CHECK(", left primer to right of right primer %d",
		  pair_expl->reversed);

  SP_AND_CHECK(", ok %d", pair_expl->ok)

  return buf;
}

const char *
p3_oligo_explain_string(const oligo_stats *stat)
{
  /* WARNING WARNING WARNING

  Static, fixed-size buffer.  If you add more calls to SP_AND_CHECK or
  IF_SP_AND_CHECK, you need to make sure the buffer cannot overflow.
  */

  static char buf[10000];
  char *bufp = buf;
  size_t bsize = 10000;
  size_t r;

  SP_AND_CHECK("considered %d", stat->considered)
  IF_SP_AND_CHECK(", would not amplify any of the ORF %d", stat->no_orf)
  IF_SP_AND_CHECK(", too many Ns %d", stat->ns)
  IF_SP_AND_CHECK(", overlap target %d", stat->target)
  IF_SP_AND_CHECK(", overlap excluded region %d", stat->excluded)
  IF_SP_AND_CHECK(", GC content failed %d", stat->gc)
  IF_SP_AND_CHECK(", GC clamp failed %d", stat->gc_clamp)
  IF_SP_AND_CHECK(", low tm %d", stat->temp_min)
  IF_SP_AND_CHECK(", high tm %d", stat->temp_max)
  IF_SP_AND_CHECK(", high any compl %d", stat->compl_any)
  IF_SP_AND_CHECK(", high end compl %d", stat->compl_end)
  IF_SP_AND_CHECK(", high hairpin stability %d", stat->hairpin_th)
  IF_SP_AND_CHECK(", high repeat similarity %d", stat->repeat_score)
  IF_SP_AND_CHECK(", long poly-x seq %d", stat->poly_x)
  IF_SP_AND_CHECK(", low sequence quality %d", stat->seq_quality)
  IF_SP_AND_CHECK(", high 3' stability %d", stat->stability)
  IF_SP_AND_CHECK(", high template mispriming score %d",
                  stat->template_mispriming)
  IF_SP_AND_CHECK(", lowercase masking of 3' end %d",stat->gmasked)
  IF_SP_AND_CHECK(", failed must_match requirements %d",stat->must_match_fail)
  IF_SP_AND_CHECK(", not in any ok left region %d", 
                  stat->not_in_any_left_ok_region)
  IF_SP_AND_CHECK(", not in any ok right region %d", 
                  stat->not_in_any_right_ok_region)
  SP_AND_CHECK(", ok %d", stat->ok)
  return buf;
}
#undef CHECK
#undef SP_AND_CHECK
#undef IF_SP_AND_CHEC

const char *
p3_get_oligo_array_explain_string(const oligo_array *oligo_array)
{
  return p3_oligo_explain_string(&oligo_array->expl);
}

const char *
p3_get_pair_array_explain_string(const pair_array_t *pair_array)
{
  return p3_pair_explain_string(&pair_array->expl);
}

const char *
libprimer3_release(void) 
{
  return "libprimer3 release 2.3.7";
}

const char *
primer3_copyright(void) 
{
  return primer3_copyright_char_star;
}

/* ============================================================ */
/* BEGIN Internal and external functions for pr_append_str      */
/* ============================================================ */

void
init_pr_append_str(pr_append_str *s) 
{
  s->data = NULL;
  s->storage_size = 0;
}

const char *
pr_append_str_chars(const pr_append_str *x) 
{
  return x->data;
}

pr_append_str *
create_pr_append_str() 
{
  /* We cannot use pr_safe_malloc here
     because this function will be called outside
     of the setjmp(....) */
  pr_append_str *ret;

  ret = (pr_append_str *) malloc(sizeof(pr_append_str));
  if (NULL == ret) return NULL;
  init_pr_append_str(ret);
  return ret;
}

void
destroy_pr_append_str_data(pr_append_str *str) 
{
  if (NULL == str) return;
  free(str->data);
  str->data = NULL;
}

void
destroy_pr_append_str(pr_append_str *str) 
{
  if (str == NULL) return;
  destroy_pr_append_str_data(str);
  free(str);
}

int
pr_append_external(pr_append_str *x,  const char *s) 
{
  int xlen, slen;

  PR_ASSERT(NULL != s);
  PR_ASSERT(NULL != x);

  if (NULL == x->data) {
    x->storage_size = 24;
    x->data = (char *) malloc(x->storage_size);
    if (NULL == x->data) return 1; /* out of memory */
    *x->data = '\0';
  }
  xlen = strlen(x->data);
  slen = strlen(s);
  if (xlen + slen + 1 > x->storage_size) {
    x->storage_size += 2 * (slen + 1);
    x->data = (char *) realloc(x->data, x->storage_size);
    if (NULL == x->data) return 1; /* out of memory */
  }
  strcpy(x->data + xlen, s);
  return 0;
}

static void
pr_append_new_chunk(pr_append_str *x, const char *s)
{
  PR_ASSERT(NULL != x)
  if (NULL == s) return;
  pr_append_w_sep(x, "; ", s);
}

int
pr_append_new_chunk_external(pr_append_str *x, const char *s)
{
  PR_ASSERT(NULL != x)
    if (NULL == s) return 0;
  return(pr_append_w_sep_external(x, "; ", s));
}

int
pr_append_w_sep_external(pr_append_str *x,
                         const char *sep,
                         const char *s)
{
  PR_ASSERT(NULL != x)
  PR_ASSERT(NULL != s)
  PR_ASSERT(NULL != sep)
    if (pr_is_empty(x)) {
      return(pr_append_external(x, s));
    } else {
        return(pr_append_external(x, sep)
               || pr_append_external(x, s));
    }
}

void
pr_set_empty(pr_append_str *x)
{
  PR_ASSERT(NULL != x);
  if (NULL != x->data) *x->data = '\0';
}

int
pr_is_empty( const pr_append_str *x)
{
  PR_ASSERT(NULL != x);
  return  NULL == x->data || '\0' == *x->data;
}

static void
pr_append_w_sep(pr_append_str *x,
                const char *sep,
                const char *s)
{
  if (pr_append_w_sep_external(x, sep, s)) longjmp(_jmp_buf, 1);
}

static void
pr_append(pr_append_str *x,
                 const char *s)
{
  if (pr_append_external(x, s)) longjmp(_jmp_buf, 1);
}

/* ============================================================ */
/* END internal and external functions for pr_append_str        */
/* ============================================================ */


/* =========================================================== */
/* Malloc and realloc wrappers that longjmp() on failure       */
/* =========================================================== */
static void *
pr_safe_malloc(size_t x)
{
  void *r = malloc(x);
  if (NULL == r) longjmp(_jmp_buf, 1);
  return r;
}

static void *
pr_safe_realloc(void *p, size_t x)
{
  void *r = realloc(p, x);
  if (NULL == r) longjmp(_jmp_buf, 1);
  return r;
}

/* =========================================================== */
/* End of malloc/realloc wrappers.                             */
/* =========================================================== */


/* ============================================================ */
/* START functions which check and modify the input             */
/* ============================================================ */

/* Fuction to set the included region and fix the start positions */
static void
_adjust_seq_args(const p3_global_settings *pa,
                 seq_args *sa,
                 pr_append_str *nonfatal_err,
                 pr_append_str *warning)
{
  int seq_len, inc_len;

  /* Create a seq for check primers if needed */
  if (pa->primer_task == check_primers) {
    if (NULL == sa->sequence) {
      fake_a_sequence(sa, pa);
    }
  }
  if (pa->primer_task == pick_sequencing_primers && sa->incl_l != -1) {
    pr_append_new_chunk(nonfatal_err,
			"Task pick_sequencing_primers cannot be combined with included region");
    return;
  }

  /*
     Complain if there is no sequence; We need to check this
     error here, because this function cannot do its work if
     sa->sequence == NULL
  */
  if (NULL == sa->sequence) {
    if (pa->primer_task == check_primers) {
      pr_append_new_chunk(nonfatal_err, "No primers provided");
    } else {
      pr_append_new_chunk(nonfatal_err, "Missing SEQUENCE tag");
    }
    return;
  }

  seq_len = strlen(sa->sequence);

  /* For pick_cloning_primers set the forced positions */
  if (pa->primer_task == pick_cloning_primers) {
    if(sa->incl_l == -1) {
      sa->force_left_start = pa->first_base_index;
      sa->force_right_start = seq_len + pa->first_base_index - 1;
    } else {
      sa->force_left_start = sa->incl_s;
      sa->force_right_start = sa->incl_s + sa->incl_l - 1;
    }
  }

  /* For pick_discriminative_primers set the forced positions */
  if (pa->primer_task == pick_discriminative_primers) {
    /* Changed here from incl_s and incl_l to sa->tar2->pairs[0][0/1] */
    if (sa->tar2.count != 1) {
      pr_append_new_chunk(nonfatal_err,
			  "Task pick_discriminative_primers requires exactly one SEQUENCE_TARGET");
    }
    sa->force_left_end = sa->tar2.pairs[0][0];
    sa->force_right_end = sa->tar2.pairs[0][0] + sa->tar2.pairs[0][1] - 1;
  }

  /* If no included region is specified,
   * use the whole sequence as included region */
  if (sa->incl_l == -1) {
    sa->incl_l = seq_len;
    sa->incl_s = pa->first_base_index;
  }

  /* Generate at least one target */
  if (pa->primer_task == pick_sequencing_primers && sa->tar2.count == 0) {
    sa->tar2.pairs[0][0] = pa->first_base_index;
    sa->tar2.pairs[0][1] = seq_len;
    sa->tar2.count = 1;
  }

  /* Fix the start of the included region and start codon */
  sa->incl_s -= pa->first_base_index;
  sa->start_codon_pos -= pa->first_base_index;

  /* Fix the start */
  sa->force_left_start -= pa->first_base_index;
  sa->force_left_end -= pa->first_base_index;
  sa->force_right_start -= pa->first_base_index;
  sa->force_right_end -= pa->first_base_index;

  /* Make it relative to included region */
  sa->force_left_start -= sa->incl_s;
  sa->force_left_end -= sa->incl_s;
  sa->force_right_start -= sa->incl_s;
  sa->force_right_end -= sa->incl_s;

  inc_len = sa->incl_s + sa->incl_l - 1;

  if ((sa->incl_l < INT_MAX) && (sa->incl_s > -1)
      && (sa->incl_l > -1) && (inc_len < seq_len) ) {
    /* Copies inluded region into trimmed_seq */
    sa->trimmed_seq = (char *) pr_safe_malloc(sa->incl_l + 1);
    _pr_substr(sa->sequence, sa->incl_s, sa->incl_l, sa->trimmed_seq);

    /* Copies inluded region into trimmed_orig_seq */
    /* edited by T. Koressaar for lowercase masking */
    sa->trimmed_orig_seq = (char *) pr_safe_malloc(sa->incl_l + 1);
    _pr_substr(sa->sequence, sa->incl_s, sa->incl_l, sa->trimmed_orig_seq);

    /* Copies the whole sequence into upcased_seq */
    sa->upcased_seq = (char *) pr_safe_malloc(strlen(sa->sequence) + 1);
    strcpy(sa->upcased_seq, sa->sequence);
    dna_to_upper(sa->upcased_seq, 1);
    /* We do not need to check for illegal characters in the return
       from dna_to_upper(), because errors are checked in
       _pr_data_control sa->trimmed_seq. */

    /* Copies the reverse complement of the whole sequence into upcased_seq_r */
    sa->upcased_seq_r = (char *) pr_safe_malloc(strlen(sa->sequence) + 1);
    p3_reverse_complement(sa->upcased_seq, sa->upcased_seq_r);
  }

  if (_check_and_adjust_intervals(sa,
                                  seq_len,
                                  pa->first_base_index,
                                  nonfatal_err, warning)) {
    return;
  }

  if (_check_and_adjust_overlap_pos(sa,
                                    sa->primer_overlap_junctions,
                                    &sa->primer_overlap_junctions_count,
                                    "SEQUENCE_OVERLAP_JUNCTION_LIST",
                                    seq_len,
                                    pa->first_base_index,
                                    nonfatal_err, warning)) {
    return;
  }

  /* Update ok regions, if non empty */
  if (sa->ok_regions.count > 0) {
    _optimize_ok_regions_list(pa, sa);
  }

}

/*
 * This function uses the max/min product size info and the max/min
 * oligo length in order to reduce the ranges of the ok regions. On
 * some imputs this improves speed dramatically.
 */
static void
_optimize_ok_regions_list(const p3_global_settings *pa,
			  seq_args *sa)
{
  /* We do this only if we enabled the optimization and
   * the primers were NOT specified. */
  if (!OPTIMIZE_OK_REGIONS || (sa->left_input) || (sa->right_input)) {
    return;
  }

  /* If any pair is allowed, no point in doing this */
  if (sa->ok_regions.any_pair) {
    return;
  }

  int pmin = INT_MAX;
  int pmax = 0;
  int omin = pa->p_args.min_size;
  int omax = pa->p_args.max_size;

  /* Determine min/max product size */
  for (int i=0; i<pa->num_intervals; i++) {
    if (pa->pr_min[i] < pmin) { pmin = pa->pr_min[i]; }
    if (pa->pr_max[i] > pmax) { pmax = pa->pr_max[i]; }
  }

  /* Update each region */
  for (int i=0; i<sa->ok_regions.count; i++) {
    int ls = -1, le = -1, rs = -1, re = -1;
    int new_ls = -1, new_le = -1, new_rs = -1, new_re = -1;
    if (sa->ok_regions.left_pairs[i][0] != -1) {
      ls = sa->ok_regions.left_pairs[i][0];
      le = sa->ok_regions.left_pairs[i][0]
	+ sa->ok_regions.left_pairs[i][1] - 1;
    }
    if (sa->ok_regions.right_pairs[i][0] != -1) {
      rs = sa->ok_regions.right_pairs[i][0];
      re = sa->ok_regions.right_pairs[i][0]
	+ sa->ok_regions.right_pairs[i][1] - 1;
    }
    /* Compute new right region based on left range and min/max values
       of product size and oligo length */
    if (ls != -1) {
      new_rs = ls + pmin - omax - 1; /* -1 just to be safe */
      new_re = le - omin + pmax + 1; /* +1 just to be safe */
      /* Adjust the ranges */
      if ((rs == -1) || (new_rs > rs)) { rs = new_rs; }
      if ((re == -1) || (new_re < re)) { re = new_re; }
      if (rs < 0) { rs = 0; }
      if (re > (signed) strlen(sa->sequence)) { re = strlen(sa->sequence); }
    }
    /* Compute new left region based on right range and min/max values
       of product size and oligo length */
    if (rs != -1) {
      new_ls = rs + omin - pmax - 1; /* -1 just to be safe */
      new_le = re - pmin + omax + 1; /* +1 just to be safe */
      /* Adjust the ranges */
      if ((ls == -1) || (new_ls > ls)) { ls = new_ls; }
      if ((le == -1) || (new_le < le)) { le = new_le; }
      if (ls < 0) { ls = 0; }
      if (le > (signed) strlen(sa->sequence)) { le = strlen(sa->sequence); }
    }
    /* Temporary testing fprintf: */
    /* fprintf(stderr, "Adjusted range [%d,%d,%d,%d] to [%d,%d,%d,%d],
	    pmin is %d, pmax is %d, omin is %d, omax is %d\n",
	    sa->ok_regions.left_pairs[i][0],
	    sa->ok_regions.left_pairs[i][0] +
	    sa->ok_regions.left_pairs[i][1] - 1,
	    sa->ok_regions.right_pairs[i][0],
	    sa->ok_regions.right_pairs[i][0] +
	    sa->ok_regions.right_pairs[i][1] - 1, ls, le, rs, re,
	    pmin, pmax, omin, omax);
    */
    sa->ok_regions.left_pairs[i][0] = ls;
    sa->ok_regions.left_pairs[i][1] = le - ls + 1;
    sa->ok_regions.right_pairs[i][0] = rs;
    sa->ok_regions.right_pairs[i][1] = re - rs + 1;
  }
  /* any_left and any_right not true anymore */
  sa->ok_regions.any_left = 0;
  sa->ok_regions.any_right = 0;
}

/*
 * Return 1 on error, 0 on success.  fake_a_sequence creates a sequence
 * out of the provided primers and puts them in sa.
 */

static int
fake_a_sequence(seq_args *sa, const p3_global_settings *pa)
{
  int i, product_size, space, ns_to_fill;
  char *rev = NULL;
  int ns_to_fill_first, ns_to_fill_second;

  /* Determine the product size */
  if ( pa->product_opt_size == PR_UNDEFINED_INT_OPT){
    product_size = pa->pr_max[0] - pa->pr_min[0];
  } else {
    product_size = pa->product_opt_size;
  } 

  space = product_size + 1;
  ns_to_fill = product_size;


  /* Calculate how many Ns have to be added */
  if (sa->left_input){
    ns_to_fill = ns_to_fill - strlen(sa->left_input);
  }
  if (sa->right_input){
    ns_to_fill = ns_to_fill - strlen(sa->right_input);
    rev = (char *) pr_safe_malloc(strlen(sa->right_input) + 1);
    p3_reverse_complement(sa->right_input, rev);
  }
  if (sa->internal_input){
          ns_to_fill = ns_to_fill - strlen(sa->internal_input);
  }
  /* Return if there are no primers provided */
  if (ns_to_fill == product_size + 1){
    return 0;
  }
  ns_to_fill_first = ns_to_fill / 2;
  ns_to_fill_second = ns_to_fill - ns_to_fill_first;
  /* Allocate the space for the sequence */  
  sa->sequence = (char *) pr_safe_malloc(space);
  *sa->sequence = '\0';
  /* Copy over the primers */
  if (sa->left_input){
    strcat(sa->sequence, sa->left_input);
  }

  /* Add the Ns*/
  for (i = 0; i < ns_to_fill_first; i++ ) {
    strcat(sa->sequence, "N\0");
  }
  if (sa->internal_input){
    strcat(sa->sequence, sa->internal_input);
  }
  for (i = 0; i < ns_to_fill_second; i++ ) {
    strcat(sa->sequence, "N\0");
  }
  if (sa->right_input){
    strcat(sa->sequence, rev);
  }
  free(rev);
  return 0;
}

/*
 * Return 1 on error, 0 on success.  Set sa->trimmed_seq and possibly modify
 * sa->tar.  Upcase and check all bases in sa->trimmed_seq.
 * TO DO -- this would probably be cleaner if it only
 * checked, rather than updated, sa.
 */
/* Check if the input in sa and pa makes sense */
int
_pr_data_control(const p3_global_settings *pa,
                 const seq_args *sa,
                 pr_append_str *glob_err,
                 pr_append_str *nonfatal_err,
                 pr_append_str *warning)
{
  static char s1[MAX_PRIMER_LENGTH+1];
  int i, pr_min, seq_len;
  char offending_char = '\0';

  seq_len = strlen(sa->sequence);

  /* If sequence quality is provided, is it as long as the sequence? */
  if (sa->n_quality !=0 && sa->n_quality != seq_len)
    pr_append_new_chunk(nonfatal_err,
                        "Error in sequence quality data");

  if ((pa->min_left_three_prime_distance < -1) ||
      (pa->min_right_three_prime_distance < -1))
    pr_append_new_chunk(nonfatal_err,
                        "Minimum 3' distance must be >= -1 "
                        "(min_*_three_prime_distance)");

  if ((pa->p_args.min_quality != 0 || pa->o_args.min_quality != 0)
      && sa->n_quality == 0)
    pr_append_new_chunk(nonfatal_err, "Sequence quality data missing");

  if (pa->first_base_index < PR_NULL_FORCE_POSITION) {
    pr_append_new_chunk(glob_err, "Value too small at tag PRIMER_FIRST_BASE_INDEX");
    return 1;
  }

  if (pa->p_args.max_template_mispriming > SHRT_MAX && pa->thermodynamic_template_alignment == 0) {
    pr_append_new_chunk(glob_err, "Value too large at tag PRIMER_MAX_TEMPLATE_MISPRIMING");
    return 1;
  }

  if (pa->pair_max_template_mispriming > SHRT_MAX && pa->thermodynamic_template_alignment == 0) {
    pr_append_new_chunk(glob_err, "Value too large at tag PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING");
    return 1;
  }

  if (pa->p_args.max_repeat_compl > SHRT_MAX && pa->thermodynamic_oligo_alignment == 0) {
    pr_append_new_chunk(glob_err, "Value too large at tag PRIMER_MAX_LIBRARY_MISPRIMING");
    return 1;
  }

  if (pa->o_args.max_repeat_compl > SHRT_MAX && pa->thermodynamic_oligo_alignment == 0) {
    pr_append_new_chunk(glob_err, "Value too large at tag PRIMER_INTERNAL_MAX_LIBRARY_MISHYB");
    return 1;
  }

  if (pa->pair_repeat_compl > SHRT_MAX && pa->thermodynamic_oligo_alignment == 0) {
    pr_append_new_chunk(glob_err, "Value too large at tag PRIMER_PAIR_MAX_LIBRARY_MISPRIMING");
    return 1;
  }

  if (pa->o_args.max_template_mispriming >= 0 && pa->thermodynamic_template_alignment==0)
    pr_append_new_chunk(glob_err,
                        "PRIMER_INTERNAL_MAX_TEMPLATE_MISHYB is not supported");
  if (pa->o_args.max_template_mispriming_th >= 0 && pa->thermodynamic_template_alignment==1)
    pr_append_new_chunk(glob_err,
                        "PRIMER_INTERNAL_MAX_TEMPLATE_MISHYB_TH is not supported");
  if (pa->p_args.min_size < 1)
    pr_append_new_chunk(glob_err, "PRIMER_MIN_SIZE must be >= 1");

  if (pa->p_args.max_size > MAX_PRIMER_LENGTH) {
    pr_append_new_chunk(glob_err,
                        "PRIMER_MAX_SIZE exceeds built-in maximum of ");
    pr_append(glob_err, MACRO_VALUE_AS_STRING(MAX_PRIMER_LENGTH));
    return 1;
  }

  if (pa->p_args.opt_size > pa->p_args.max_size) {
    pr_append_new_chunk(glob_err,
                        "PRIMER_{OPT,DEFAULT}_SIZE > PRIMER_MAX_SIZE");
    return 1;
  }

  if (pa->p_args.opt_size < pa->p_args.min_size) {
    pr_append_new_chunk(glob_err,
                        "PRIMER_{OPT,DEFAULT}_SIZE < PRIMER_MIN_SIZE");
    return 1;
  }

  if (pa->o_args.max_size > MAX_PRIMER_LENGTH) {
    pr_append_new_chunk(glob_err,
                        "PRIMER_INTERNAL_MAX_SIZE exceeds built-in maximum");
    return 1;
  }

  if (pa->o_args.opt_size > pa->o_args.max_size) {
    pr_append_new_chunk(glob_err,
                        "PRIMER_INTERNAL_{OPT,DEFAULT}_SIZE > MAX_SIZE");
    return 1;
  }

  if (pa->o_args.opt_size < pa->o_args.min_size) {
    pr_append_new_chunk(glob_err,
                        "PRIMER_INTERNAL_{OPT,DEFAULT}_SIZE < MIN_SIZE");
    return 1;
  }

  /* A GC clamp can not be bigger then the primer */
  if (pa->gc_clamp > pa->p_args.min_size) {
    pr_append_new_chunk(glob_err,
                        "PRIMER_GC_CLAMP > PRIMER_MIN_SIZE");
    return 1;
  }

  /* Must be >= 0 and <= 5 */
  if ((pa->max_end_gc < 0) 
      || (pa->max_end_gc > 5)) {
    pr_append_new_chunk(glob_err,
                        "PRIMER_MAX_END_GC must be between 0 to 5");
    return 1;
  }

  /* Product size must be provided */
  if (0 == pa->num_intervals) {
    pr_append_new_chunk( glob_err,
                         "Empty value for PRIMER_PRODUCT_SIZE_RANGE");
    return 1;
  }
  for (i = 0; i < pa->num_intervals; i++) {
    if (pa->pr_min[i] > pa->pr_max[i] || pa->pr_min[i] < 0) {
      pr_append_new_chunk(glob_err,
                          "Illegal element in PRIMER_PRODUCT_SIZE_RANGE");
      return 1;
    }
  }

  pr_min = INT_MAX;
  /* Check if the primer is bigger then the product */
  for(i=0;i<pa->num_intervals;i++)
    if(pa->pr_min[i]<pr_min) pr_min=pa->pr_min[i];

  if (pa->p_args.max_size > pr_min) {
    pr_append_new_chunk(glob_err,
                        "PRIMER_MAX_SIZE > min PRIMER_PRODUCT_SIZE_RANGE");
    return 1;
  }

  if ((pa->pick_internal_oligo == 1 )
      && pa->o_args.max_size > pr_min) {
    pr_append_new_chunk(glob_err,
                        "PRIMER_INTERNAL_MAX_SIZE > min PRIMER_PRODUCT_SIZE_RANGE");
    return 1;
  }

  /* There must be at least one primer returned */
  if (pa->num_return < 1) {
    pr_append_new_chunk(glob_err,
                        "PRIMER_NUM_RETURN < 1");
    return 1;
  }

  if ((pa->p_args.must_match_five_prime != NULL) && (strlen(pa->p_args.must_match_five_prime) != 5)) {
    pr_append_new_chunk(glob_err,
                        "PRIMER_MUST_MATCH_FIVE_PRIME must have 5 characters");
    return 1;
  }
  if ((pa->p_args.must_match_three_prime != NULL) && (strlen(pa->p_args.must_match_three_prime) != 5)) {
    pr_append_new_chunk(glob_err,
                        "PRIMER_MUST_MATCH_THREE_PRIME must have 5 characters");
    return 1;
  }
  if ((pa->o_args.must_match_five_prime != NULL) && (strlen(pa->o_args.must_match_five_prime) != 5)) {
    pr_append_new_chunk(glob_err,
                        "PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME must have 5 characters");
    return 1;
  }
  if ((pa->o_args.must_match_three_prime != NULL) && (strlen(pa->o_args.must_match_three_prime) != 5)) {
    pr_append_new_chunk(glob_err,
                        "PRIMER_INTERNAL_MUST_MATCH_THREE_PRIME must have 5 characters");
    return 1;
  }

  if (sa->incl_l >= INT_MAX) {
    pr_append_new_chunk(nonfatal_err, "Value for SEQUENCE_INCLUDED_REGION too large");
    return 1;
  }

  if (sa->incl_s < 0 || sa->incl_l < 0
      || sa->incl_s + sa->incl_l > seq_len) {
    pr_append_new_chunk(nonfatal_err, "Illegal value for SEQUENCE_INCLUDED_REGION");
    return 1;
  }

  /* The product must fit in the included region */
  if (sa->incl_l < pr_min && pa->pick_left_primer == 1
      && pa->pick_right_primer == 1) {
    if (pa->primer_task == check_primers) {
      pr_append_new_chunk(warning,
                          "SEQUENCE_INCLUDED_REGION length < min PRIMER_PRODUCT_SIZE_RANGE");
    } else if (pa->primer_task != pick_primer_list) {
      pr_append_new_chunk(nonfatal_err,
                          "SEQUENCE_INCLUDED_REGION length < min PRIMER_PRODUCT_SIZE_RANGE");
    }

    if (pa->primer_task == generic) {
      return 1;
    }
  }

  if (pa->max_end_stability < 0) {
    pr_append_new_chunk(nonfatal_err,
                        "PRIMER_MAX_END_STABILITY must be non-negative");
    return 1;
  }

  /* Is the startodon ATG and in the incl. region */
  if (!PR_START_CODON_POS_IS_NULL(sa)) {
    if (!PR_POSITION_PENALTY_IS_NULL(pa)) {
      pr_append_new_chunk(nonfatal_err,
                          "Cannot accept both SEQUENCE_START_CODON_POSITION and non-default ");
      pr_append(nonfatal_err,
                "arguments for PRIMER_INSIDE_PENALTY or PRIMER_OUTSIDE_PENALTY");
    }
    if (sa->start_codon_pos  > (sa->incl_s + sa->incl_l - 3)) {
      pr_append_new_chunk(nonfatal_err,
                          "Start codon position not contained in SEQUENCE_INCLUDED_REGION");
      return 1;
    } else {
      if (sa->start_codon_pos >= 0
          && ((sa->sequence[sa->start_codon_pos] != 'A'
               && sa->sequence[sa->start_codon_pos] != 'a')
              || (sa->sequence[sa->start_codon_pos + 1] != 'T'
                  && sa->sequence[sa->start_codon_pos + 1] != 't')
              || (sa->sequence[sa->start_codon_pos + 2] != 'G'
                  && sa->sequence[sa->start_codon_pos + 2] != 'g'))) {
        pr_append_new_chunk(nonfatal_err,
                            "No start codon at SEQUENCE_START_CODON_POSITION");
        return 1;
      }
    }
  }

  if (NULL != sa->quality) {
    if(pa->p_args.min_quality != 0 && pa->p_args.min_quality < pa->quality_range_min) {
      pr_append_new_chunk(glob_err,
                          "PRIMER_MIN_QUALITY < PRIMER_QUALITY_RANGE_MIN");
      return 1;
    }
    if (pa->p_args.min_quality != 0 && pa->p_args.min_quality > pa->quality_range_max) {
      pr_append_new_chunk(glob_err,
                          "PRIMER_MIN_QUALITY > PRIMER_QUALITY_RANGE_MAX");
      return 1;
    }
    if (pa->o_args.min_quality != 0 && pa->o_args.min_quality <pa->quality_range_min) {
      pr_append_new_chunk(glob_err,
                          "PRIMER_INTERNAL_MIN_QUALITY < PRIMER_QUALITY_RANGE_MIN");
      return 1;
    }
    if (pa->o_args.min_quality != 0 && pa->o_args.min_quality > pa->quality_range_max) {
      pr_append_new_chunk(glob_err,
                          "PRIMER_INTERNAL_MIN_QUALITY > PRIMER_QUALITY_RANGE_MAX");
      return 1;
    }
    for(i=0; i < sa->n_quality; i++) {
      if(sa->quality[i] < pa->quality_range_min ||
         sa->quality[i] > pa->quality_range_max) {
        pr_append_new_chunk(nonfatal_err,
                            "Sequence quality score out of range");
        return 1;
      }
    }
  }
  else if (pa->p_args.weights.seq_quality || pa->o_args.weights.seq_quality) {
    pr_append_new_chunk(nonfatal_err,
                        "Sequence quality is part of objective function but sequence quality is not defined");
    return 1;
  }

  if ((offending_char = dna_to_upper(sa->trimmed_seq, 0))) {
    if (pa->liberal_base) {
      pr_append_new_chunk(warning,
                          "Unrecognized base in input sequence");
    }
    else {
      pr_append_new_chunk(nonfatal_err,
                          "Unrecognized base in input sequence");
      return 1;
    }
  }

  if (pa->p_args.opt_tm < pa->p_args.min_tm
      || pa->p_args.opt_tm > pa->p_args.max_tm) {
    pr_append_new_chunk(glob_err,
                        "Optimum primer Tm lower than minimum or higher than maximum");
    return 1;
  }

  if (pa->o_args.opt_tm < pa->o_args.min_tm
      || pa->o_args.opt_tm > pa->o_args.max_tm) {
    pr_append_new_chunk(glob_err,
                        "Optimum internal oligo Tm lower than minimum or higher than maximum");
    return 1;
  }

  if (pa->p_args.min_gc > pa->p_args.max_gc
      || pa->p_args.min_gc > 100
      || pa->p_args.max_gc < 0){
    pr_append_new_chunk(glob_err,
                        "Illegal value for PRIMER_MAX_GC and PRIMER_MIN_GC");
    return 1;
  }

  if (pa->o_args.min_gc > pa->o_args.max_gc
      || pa->o_args.min_gc > 100
      || pa->o_args.max_gc < 0) {
    pr_append_new_chunk(glob_err,
                        "Illegal value for PRIMER_INTERNAL_OLIGO_GC");
    return 1;
  }
  if (pa->p_args.num_ns_accepted < 0) {
    pr_append_new_chunk(glob_err,
                        "Illegal value for PRIMER_MAX_NS_ACCEPTED");
    return 1;
  }
  if (pa->o_args.num_ns_accepted < 0){
    pr_append_new_chunk(glob_err,
                        "Illegal value for PRIMER_INTERNAL_MAX_NS_ACCEPTED");
    return 1;
  }
  if (pa->p_args.max_self_any < 0 || pa->p_args.max_self_any > SHRT_MAX
      || pa->p_args.max_self_end < 0 || pa->p_args.max_self_end > SHRT_MAX
      || pa->pair_compl_any < 0 || pa->pair_compl_any > SHRT_MAX 
      || pa->pair_compl_end < 0 || pa->pair_compl_end > SHRT_MAX) {
    pr_append_new_chunk(glob_err,
                        "Illegal value for primer complementarity restrictions");
    return 1;
  }
  
  if (pa->p_args.max_self_any_th < 0
      || pa->p_args.max_self_end_th < 0 || pa->p_args.max_hairpin_th < 0
      || pa->pair_compl_any_th < 0 || pa->pair_compl_end_th < 0) {
     pr_append_new_chunk(glob_err,
                         "Illegal value for primer complementarity restrictions (thermod. approach)");
     return 1;
  }
   
  if( pa->o_args.max_self_any < 0 || pa->o_args.max_self_any > SHRT_MAX
      || pa->o_args.max_self_end < 0 || pa->o_args.max_self_end > SHRT_MAX) {
    pr_append_new_chunk(glob_err,
                        "Illegal value for internal oligo complementarity restrictions");
    return 1;
  }
  
  if( pa->o_args.max_self_any_th < 0
      || pa->o_args.max_self_end_th < 0 || pa->o_args.max_hairpin_th < 0) {
     pr_append_new_chunk(glob_err,
                         "Illegal value for internal oligo complementarity restrictions");
     return 1;
  }
   
  if (pa->p_args.salt_conc <= 0 || pa->p_args.dna_conc<=0){
    pr_append_new_chunk(glob_err,
                        "Illegal value for primer salt or dna concentration");
    return 1;
  }

  if((pa->p_args.dntp_conc < 0 && pa->p_args.divalent_conc !=0 )
     || pa->p_args.divalent_conc<0){ /* added by T.Koressaar */
    pr_append_new_chunk(glob_err, "Illegal value for primer divalent salt or dNTP concentration");
    return 1;
  }

  if(pa->o_args.salt_conc<=0||pa->o_args.dna_conc<=0){
    pr_append_new_chunk(glob_err,
                        "Illegal value for internal oligo salt or dna concentration");
    return 1;
  }

  if((pa->o_args.dntp_conc<0 && pa->o_args.divalent_conc!=0)
     || pa->o_args.divalent_conc < 0) { /* added by T.Koressaar */
    pr_append_new_chunk(glob_err,
                        "Illegal value for internal oligo divalent salt or dNTP concentration");
    return 1;
  }

  if (!_PR_DEFAULT_POSITION_PENALTIES(pa) && sa->tar2.count > 1) {
    pr_append_new_chunk(nonfatal_err,
                        "Non-default inside penalty or outside penalty ");
    pr_append(nonfatal_err,
              "is valid only when number of targets <= 1");
  }
  if (!_PR_DEFAULT_POSITION_PENALTIES(pa) && 0 == sa->tar2.count) {
    pr_append_new_chunk(warning,
                        "Non-default inside penalty or outside penalty ");
    pr_append(warning,
              "has no effect when number of targets is 0");
  }
  if (pa->pick_internal_oligo != 1 && sa->internal_input) {
    pr_append_new_chunk(nonfatal_err,
                        "Not specified to pick internal oligos");
    pr_append(nonfatal_err,
              " but a specific internal oligo is provided");
  }
  if (sa->internal_input) {
    if (strlen(sa->internal_input) > MAX_PRIMER_LENGTH) {
      pr_append_new_chunk(glob_err,
                          "Specified internal oligo exceeds built-in maximum of ");
      pr_append(glob_err, MACRO_VALUE_AS_STRING(MAX_PRIMER_LENGTH));
      return 1;
    }
    if (strlen(sa->internal_input) > (unsigned) pa->o_args.max_size)
      pr_append_new_chunk(warning,
                          "Specified internal oligo > PRIMER_INTERNAL_MAX_SIZE");

    if (strlen(sa->internal_input) < (unsigned) pa->o_args.min_size)
      pr_append_new_chunk(warning,
                          "Specified internal oligo < PRIMER_INTERNAL_MIN_SIZE");

    if (!strstr_nocase(sa->sequence, sa->internal_input))
      pr_append_new_chunk(nonfatal_err,
                          "Specified internal oligo not in sequence");
    else if (!strstr_nocase(sa->trimmed_seq, sa->internal_input))
      pr_append_new_chunk(nonfatal_err,
                          "Specified internal oligo not in Included Region");
  }
  if (sa->left_input) {
    if (strlen(sa->left_input) > MAX_PRIMER_LENGTH) {
      pr_append_new_chunk(glob_err,
                          "Specified left primer exceeds built-in maximum of ");
      pr_append(glob_err, MACRO_VALUE_AS_STRING(MAX_PRIMER_LENGTH));
      return 1;
    }
    if (strlen(sa->left_input) > (unsigned) pa->p_args.max_size)
      pr_append_new_chunk(warning,
                          "Specified left primer > PRIMER_MAX_SIZE");
    if (strlen(sa->left_input) < (unsigned) pa->p_args.min_size)
      pr_append_new_chunk(warning,
                          "Specified left primer < PRIMER_MIN_SIZE");
    if (!strstr_nocase(sa->sequence, sa->left_input))
      pr_append_new_chunk(nonfatal_err,
                          "Specified left primer not in sequence");
    else if (!strstr_nocase(sa->trimmed_seq, sa->left_input))
      pr_append_new_chunk(nonfatal_err,
                          "Specified left primer not in Included Region");
  }
  if (sa->right_input) {
    if (strlen(sa->right_input) > MAX_PRIMER_LENGTH) {
      pr_append_new_chunk(glob_err,
                          "Specified right primer exceeds built-in maximum of ");
      pr_append(glob_err, MACRO_VALUE_AS_STRING(MAX_PRIMER_LENGTH));
      return 1;
    }
    if (strlen(sa->right_input) < (unsigned) pa->p_args.min_size)
      pr_append_new_chunk(warning,
                          "Specified right primer < PRIMER_MIN_SIZE");
    if (strlen(sa->right_input) > (unsigned) pa->p_args.max_size) {
      pr_append_new_chunk(warning,
                          "Specified right primer > PRIMER_MAX_SIZE");
    } else { /* We do not want to overflow s1. */
      p3_reverse_complement(sa->right_input,s1);
      if (!strstr_nocase(sa->sequence, s1))
        pr_append_new_chunk(nonfatal_err,
                            "Specified right primer not in sequence");
      else if (!strstr_nocase(sa->trimmed_seq, s1))
        pr_append_new_chunk(nonfatal_err,
                            "Specified right primer not in Included Region");
    }
  }

  if ((pa->pr_pair_weights.product_tm_lt ||
       pa->pr_pair_weights.product_tm_gt)
      && pa->product_opt_tm == PR_UNDEFINED_DBL_OPT) {
    pr_append_new_chunk(glob_err,
                        "Product temperature is part of objective function while optimum temperature is not defined");
    return 1;
  }

  if((pa->pr_pair_weights.product_size_lt ||
      pa->pr_pair_weights.product_size_gt)
     && pa->product_opt_size == PR_UNDEFINED_INT_OPT){
    pr_append_new_chunk(glob_err,
                        "Product size is part of objective function while optimum size is not defined");
    return 1;
  }

  if ((pa->p_args.weights.gc_content_lt ||
       pa->p_args.weights.gc_content_gt)
      && pa->p_args.opt_gc_content == DEFAULT_OPT_GC_PERCENT) {
    pr_append_new_chunk(glob_err,
                        "Primer GC content is part of objective function while optimum gc_content is not defined");
    return 1;
  }

  if ((pa->o_args.weights.gc_content_lt ||
       pa->o_args.weights.gc_content_gt)
      && pa->o_args.opt_gc_content == DEFAULT_OPT_GC_PERCENT) {
    pr_append_new_chunk(glob_err,
                        "Hyb probe GC content is part of objective function "
                        "while optimum gc_content is not defined");
    return 1;
  }

  if ((pa->pick_internal_oligo != 1) &&
      (pa->pr_pair_weights.io_quality)) {
    pr_append_new_chunk(glob_err,
                        "Internal oligo quality is part of objective function "
                        "while internal oligo choice is not required");
    return 1;
  }

  if (pa->p_args.weights.repeat_sim && (!seq_lib_num_seq(pa->p_args.repeat_lib))) {
    pr_append_new_chunk(glob_err,
                        "Mispriming score is part of objective function, but mispriming library is not defined");
    return 1;
  }

  if (pa->o_args.weights.repeat_sim && (!seq_lib_num_seq(pa->o_args.repeat_lib))) {
    pr_append_new_chunk(glob_err,
                        "Internal oligo mispriming score is part of objective function while mishyb library is not defined");
    return 1;
  }

  if (pa->pr_pair_weights.repeat_sim && (!(seq_lib_num_seq(pa->p_args.repeat_lib)))) {
    pr_append_new_chunk(glob_err,
                        "Mispriming score is part of objective function, "
                        "but mispriming library is not defined");
    return 1;
  }

  if(pa->pr_pair_weights.io_quality
     && pa->pick_internal_oligo == 0 ) {
    pr_append_new_chunk(glob_err,
                        "Internal oligo quality is part of objective function"
                        " while internal oligo choice is not required");
    return 1;
  }

  if (pa->sequencing.lead < 0) {
    pr_append_new_chunk(glob_err,
                        "Illegal value for PRIMER_SEQUENCING_LEAD");
    return 1;
  }

  if (pa->sequencing.interval < 0) {
    pr_append_new_chunk(glob_err,
                        "Illegal value for PRIMER_SEQUENCING_INTERVAL");
    return 1;
  }

  if (pa->sequencing.accuracy < 0) {
    pr_append_new_chunk(glob_err,
                        "Illegal value for PRIMER_SEQUENCING_ACCURACY");
    return 1;
  }

  if (pa->sequencing.spacing < 0) {
    pr_append_new_chunk(glob_err,
                        "Illegal value for PRIMER_SEQUENCING_SPACING");
    return 1;
  }

  if(pa->sequencing.interval > pa->sequencing.spacing) {
    pr_append_new_chunk(glob_err,
                        "PRIMER_SEQUENCING_INTERVAL > PRIMER_SEQUENCING_SPACING");
    return 1;
  }

  if(pa->sequencing.accuracy > pa->sequencing.spacing) {
    pr_append_new_chunk(glob_err,
                        "PRIMER_SEQUENCING_ACCURACY > PRIMER_SEQUENCING_SPACING");
    return 1;
  }

  if(pa->sequencing.lead > pa->sequencing.spacing) {
    pr_append_new_chunk(glob_err,
                        "PRIMER_SEQUENCING_LEAD > PRIMER_SEQUENCING_SPACING");
    return 1;
  }

  if((sa->force_left_start > -1) && (sa->force_left_end > -1)
     && (sa->force_left_start > sa->force_left_end)) {
    pr_append_new_chunk(glob_err,
                        "SEQUENCE_FORCE_LEFT_START > SEQUENCE_FORCE_LEFT_END");
    return 1;
  }

  if((sa->force_right_end > -1) && (sa->force_right_start > -1)
     && (sa->force_right_end > sa->force_right_start)) {
    pr_append_new_chunk(glob_err,
                        "SEQUENCE_FORCE_RIGHT_END > SEQUENCE_FORCE_RIGHT_START");
    return 1;
  }

  if (pa->min_5_prime_overlap_of_junction < 1) {
    pr_append_new_chunk(glob_err,
                        "Illegal value for PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION");
    return 1;
  }

  if (pa->min_3_prime_overlap_of_junction < 1) {
    pr_append_new_chunk(glob_err,
                        "Illegal value for PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION");
    return 1;
  }

  if ((sa->primer_overlap_junctions_count > 0) &&
      (pa->min_5_prime_overlap_of_junction > (pa->p_args.max_size / 2))) {
    pr_append_new_chunk(glob_err,
                        "PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION > PRIMER_MAX_SIZE / 2");
    return 1;
  }

  if ((sa->primer_overlap_junctions_count > 0) && 
      (pa->min_3_prime_overlap_of_junction > (pa->p_args.max_size / 2))) {
    pr_append_new_chunk(glob_err,
                        "PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION > PRIMER_MAX_SIZE / 2");
    return 1;
  }
  
  if (pa->p_args.divalent_conc > 0.0 && pa->p_args.dntp_conc <= 0.0) {
    pr_append_new_chunk(warning,
			"PRIMER_SALT_DIVALENT > 0.0 "
			"but PRIMER_DNTP_CONC <= 0.0; "
			"use reasonable value for PRIMER_DNTP_CONC");
  }

  if ((pa->p_args.must_match_five_prime != NULL) &&
      (test_must_match_parameters(pa->p_args.must_match_five_prime))) {
    pr_append_new_chunk(glob_err,
                        "Illegal values for PRIMER_MUST_MATCH_FIVE_PRIME");
    return 1;
  }

  if ((pa->p_args.must_match_three_prime != NULL) &&
      (test_must_match_parameters(pa->p_args.must_match_three_prime))) {
    pr_append_new_chunk(glob_err,
                        "Illegal values for PRIMER_MUST_MATCH_THREE_PRIME");
    return 1;
  }
  
  if ((pa->o_args.must_match_five_prime != NULL) &&
      (test_must_match_parameters(pa->o_args.must_match_five_prime))) {
    pr_append_new_chunk(glob_err,
                        "Illegal values for PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME");
    return 1;
  }
  
  if ((pa->o_args.must_match_three_prime != NULL) &&
      (test_must_match_parameters(pa->o_args.must_match_three_prime))) {
    pr_append_new_chunk(glob_err,
                        "Illegal values for PRIMER_INTERNAL_MUST_MATCH_THREE_PRIME");
    return 1;
  }
  
  return (NULL == nonfatal_err->data && NULL == glob_err->data) ? 0 : 1;
} /* _pr_data_control */

static int
_check_and_adjust_overlap_pos(seq_args *sa,
                              int *list,
                              int *count, 
                              const char *tag,
                              int seq_len,
                              int first_index,
                              pr_append_str *nonfatal_err,
                              pr_append_str *warning) {
  int i;
  int outside_warning_issued = 0;
  char buffer[255];

  if (*count == 0) {
    return 0;
  }

  for (i = 0; i < *count; i++) {
    /* Subtract first_index from the "must overlap" positions */
    list[i] -= first_index;

    if (list[i] >= seq_len) {
      sprintf(buffer, "%s beyond end of sequence", tag);
      pr_append_new_chunk(nonfatal_err, buffer);
      return 1;
    }

    if (list[i] < 0) {
      sprintf(buffer, "Negative %s length", tag);
      pr_append_new_chunk(nonfatal_err, buffer);
      return 1;
    }

    /* Cause the intron positions to be relative to the included region. */
    list[i] -= sa->incl_s;

    /* Check that intron positions are within the included region. */
    if (list[i] < 0 
        || list[i] > sa->incl_l) {
      if (!outside_warning_issued) {
        sprintf(buffer, "%s outside of INCLUDED_REGION", tag);
        pr_append_new_chunk(warning, buffer);
        outside_warning_issued = 1;
      }
    }
  }

  return 0;
}

static int
_check_and_adjust_intervals(seq_args *sa,
                            int seq_len,
                            int first_index,
                            pr_append_str * nonfatal_err,
                            pr_append_str *warning) {

  if (_check_and_adjust_1_interval("TARGET",
                                   sa->tar2.count, 
                                   sa->tar2.pairs,
                                   seq_len,
                                   first_index,
                                   nonfatal_err, sa, warning, 0)
      == 1) return 1;
  sa->start_codon_pos -= sa->incl_s;

  if (_check_and_adjust_1_interval("EXCLUDED_REGION",
                                   sa->excl2.count, sa->excl2.pairs,
                                   seq_len, first_index, 
                                   nonfatal_err, sa, warning, 0)
      == 1) return 1;

  if (_check_and_adjust_1_interval("PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION",
                                   sa->excl_internal2.count,
                                   sa->excl_internal2.pairs,
                                   seq_len,
                                   first_index,
                                   nonfatal_err, sa, warning, 0)
      == 1) return 1;
  if (_check_and_adjust_1_interval("PRIMER_PAIR_OK_REGION_LIST",
                                   sa->ok_regions.count, 
                                   sa->ok_regions.left_pairs,
                                   seq_len,
                                   first_index,
                                   nonfatal_err, sa, warning, 1)
      == 1) return 1;
  if (_check_and_adjust_1_interval("PRIMER_PAIR_OK_REGION_LIST",
                                   sa->ok_regions.count, 
                                   sa->ok_regions.right_pairs,
                                   seq_len,
                                   first_index,
                                   nonfatal_err, sa, warning, 1)
      == 1) return 1;
  return 0;
}

/*
 * Check intervals, and add any errors to err.
 * Update the start of each interval to
 * be relative to the start of the included region.
 */
static int
_check_and_adjust_1_interval(const char *tag_name,
                             int num_intervals,
                             interval_array_t intervals,
                             int seq_len,
                             int first_index,
                             pr_append_str *err,
                             seq_args *sa,
                             pr_append_str *warning,
                             int empty_allowed)
{
  int i;
  int outside_warning_issued = 0;

  /* Subtract the first_index from the start positions in the interval
     array */
  for (i = 0; i < num_intervals; i++) {
    if (empty_allowed && (intervals[i][0] == -1) && (intervals[i][1] == -1))
      continue;
    if (empty_allowed && (((intervals[i][0] == -1) && (intervals[i][1] != -1)) 
                          || ((intervals[i][0] != -1) && (intervals[i][1] == -1)))) {
      pr_append_new_chunk(err, tag_name);
      pr_append(err, " illegal interval");
      return 1;
    }
    intervals[i][0] -= first_index;
  }

  for (i=0; i < num_intervals; i++) {
    if (empty_allowed && (intervals[i][0] == -1) && (intervals[i][1] == -1))
      continue;
    if (intervals[i][0] + intervals[i][1] > seq_len) {
      pr_append_new_chunk(err, tag_name);
      pr_append(err, " beyond end of sequence");
      return 1;
    }
    /* Cause the interval start to be relative to the included region. */
    intervals[i][0] -= sa->incl_s;
    /* Check that intervals are within the included region. */
    if (intervals[i][0] < 0
        || intervals[i][0] + intervals[i][1] > sa->incl_l) {
      if (!outside_warning_issued) {
        pr_append_new_chunk(warning, tag_name);
        pr_append(warning,
                  " outside of INCLUDED_REGION");
        outside_warning_issued = 1;
      }
    }
    if (intervals[i][1] < 0) {
      pr_append_new_chunk(err, "Negative ");
      pr_append(err, tag_name);
      pr_append(err, " length");
      return 1;
    }
  }
  return 0;
} /* _check_and_adjust_intervals  */

/* ============================================================ */
/* END functions that check and modify the input                */
/* ============================================================ */

/* ============================================================ */
/* BEGIN create primer files functions                          */
/* ============================================================ */

/* START OF WHAT SHOULD GO TO IO_PRIMER_FILES.C */

/* 
   Creates up to three files, the names of which are based on the
   argument 'file_name'.  One file is a table of forward primers, one
   a table of reverse primers, and one a table of internal
   hybridization oligos, depending on what the caller to
   choose_primers() requested (in the 'pa' argument).  Returns 0 on
   success, 1 on error.  On error, check errno for ENOMEM. Used to
   implement P3_FILE_FLAG=1.
 */
int
p3_print_oligo_lists(const p3retval *retval,
                     const seq_args *sa,
                     const p3_global_settings *pa,
                     pr_append_str *err,
                     const char *file_name)
{
  /* Figure out if the sequence starts at position 1 or 0 */
  int   first_base_index = pa->first_base_index;
  int   ret;
  /* Start building up a filename */
  char *file;
  FILE *fh;

  if (setjmp(_jmp_buf) != 0) {
    return 1;  /* If we get here, that means we returned via a longjmp.
                  In this case errno should be ENOMEM. */

  }

  file = (char *) malloc(strlen(sa->sequence_name) + 5);
  if (NULL == file) return 1; /* ENOMEM */

  /* Check if the left primers have to be printed */
  if( pa->pick_left_primer ) {
    /* Create the file name and open file*/
    strcpy(file, sa->sequence_name);
    strcat(file, ".for");
    if (!(fh = fopen(file,"w"))) {
      if (pr_append_new_chunk_external(err, "Unable to open file "))
        return 1;
      if (pr_append_external(err, file)) return 1;
      if (pr_append_external(err, " for writing")) return 1;
      free(file);
      return 1;
    }

    /* Print the content to the file */
    ret = 
      p3_print_one_oligo_list(sa, retval->fwd.num_elem,
			      retval->fwd.oligo, OT_LEFT,
			      first_base_index, NULL != pa->p_args.repeat_lib, 
			      fh,pa->thermodynamic_oligo_alignment);
    fclose(fh);
    if (ret) return 1;
  }

  /* Check if the right primers have to be printed */
  if (pa->pick_right_primer) {
    strcpy(file, sa->sequence_name);
    strcat(file, ".rev");
    if (!(fh = fopen(file,"w"))) {
      pr_append_new_chunk(err, "Unable to open file ");
      pr_append(err, file);
      pr_append(err, " for writing");
      free(file);
      return 1;
    }
    /* Print the content to the file */
    ret = p3_print_one_oligo_list(sa, retval->rev.num_elem,
				  retval->rev.oligo, OT_RIGHT,
				  first_base_index, NULL != pa->p_args.repeat_lib, 
				  fh, pa->thermodynamic_oligo_alignment);

    fclose(fh);
    if (ret) return 1;
  }

  /* Check if the internal oligos have to be printed */
  if (pa->pick_internal_oligo) {
    /* Create the file name and open file*/
    strcpy(file, sa->sequence_name);
    strcat(file, ".int");
    if (!(fh = fopen(file,"w"))) {
      if (pr_append_new_chunk_external(err, "Unable to open file "))
        return 1;
      if (pr_append_external(err, file)) return 1;
      if (pr_append_external(err, " for writing")) return 1;
      free(file);
      return 1;
    }
    /* Print the content to the file */
    ret = p3_print_one_oligo_list(sa, retval->intl.num_elem,
                                  retval->intl.oligo, OT_INTL,
                                  first_base_index, 
                                  NULL != pa->o_args.repeat_lib, 
				  fh, pa->thermodynamic_oligo_alignment);
    fclose(fh);
    if (ret) return 1;
  }
  free(file);
  return 0;
}

/* Print out the content of one primer array */
/* Return 1 on error, otherwise 0. */
int
p3_print_one_oligo_list(const seq_args *sa,
                        int n,
                        const primer_rec *oligo_arr,
                        const oligo_type o_type,
                        const int first_base_index,
                        const int print_lib_sim,
                        FILE  *fh,
                        const int thermodynamic_oligo_alignment)
{
  int i;

  /* Print out the header for the table */
  if (print_list_header(fh, o_type, first_base_index, print_lib_sim, thermodynamic_oligo_alignment))
    return 1; /* error */
  /* Iterate over the array */
  for (i = 0; i < n; i++) {
    /* Print each single oligo */
    if (print_oligo(fh, sa, i, &oligo_arr[i], o_type,
                    first_base_index, print_lib_sim, thermodynamic_oligo_alignment))
      return 1; /* error */
  }
  return 0; /* success */
}

static int
print_list_header(FILE *fh,
                  oligo_type type,
                  int first_base_index,
                  int print_lib_sim,
                  int thermodynamic_oligo_alignment)
{
  int ret;
  ret = fprintf(fh, "ACCEPTABLE %s\n",
                OT_LEFT == type ? "LEFT PRIMERS"
                : OT_RIGHT == type ? "RIGHT PRIMERS" : "INTERNAL OLIGOS");
  if (ret < 0) return 1;

  ret = fprintf(fh, "                               %4d-based     ",
                first_base_index);
  if (ret < 0) return 1;

   if(thermodynamic_oligo_alignment)
     ret = fprintf(fh, "#                self   self hair-");
   else 
     ret = fprintf(fh, "#               self  self");
   if (ret < 0) return 1;
   if (print_lib_sim)
     /* ret = fprintf(fh, "#               self  self   lib  qual-\n"); */
     ret = fprintf(fh, "   lib");
   if (ret < 0) return 1;
   ret = fprintf(fh, "  qual-\n");
   if (ret < 0) return 1;
   
   ret = fprintf(fh, "   # sequence                       start ln  ");
   if (ret < 0) return 1;
   
   ret = fprintf(fh, "N   GC%%     Tm");
   if (ret < 0) return 1;
   if(thermodynamic_oligo_alignment)
     ret = fprintf(fh, " any_th end_th   pin");
   else 
     ret = fprintf(fh, "   any   end");
   if (ret < 0) return 1;
   if (print_lib_sim)
     ret = fprintf(fh, "   sim   lity\n");
   else
     ret = fprintf(fh, "   lity\n");
   
  if (ret < 0) return 1;
  return 0;
}

static int
print_oligo(FILE *fh,
            const seq_args *sa,
            int index,
            const primer_rec *h,
            oligo_type type,
            int first_base_index,
            int print_lib_sim,
            int thermodynamical_approach)
{
  int ret;

  /* WARNING, *p points to static storage that is overwritten on next
     call to pr_oligo_sequence or pr_oligo_rev_c_sequence. */
  char *p = 
    (OT_RIGHT != type)
    ? pr_oligo_sequence(sa, h)
    : pr_oligo_rev_c_sequence(sa, h);

  ret = fprintf(fh,
		"%4d %-30s %5d %2d %2d %5.2f %5.3f %5.2f %5.2f", 
		index, p, h->start+sa->incl_s + first_base_index,
		h->length,
		h->num_ns, h->gc_content, h->temp,
		h->self_any,
		h->self_end);
  if (ret < 0) return 1;

  if (1==thermodynamical_approach) {
    ret = fprintf(fh, " %5.2f", h->hairpin_th);
    if (ret < 0) return 1;
  }

  if (ret < 0) return 1;
  if (print_lib_sim) {
    PR_ASSERT(h->repeat_sim.score != NULL);
    ret = fprintf(fh, " %5.2f",
		  h->repeat_sim.score[h->repeat_sim.max]);
    if (ret < 0) return 1;
  }
  ret = fprintf(fh, " %6.3f\n", h->quality);
  if (ret < 0) return 1;
  else return 0;
}

/* END OF WHAT SHOULD GO TO IO_PRIMER_FILES.C */

/* ============================================================ */
/* END create primer files functions                            */
/* ============================================================ */

/* ============================================================
  * BEGIN p3_read_line, a utility used
  * both for reading boulder input and for reading sequence
  * libraries.
  * ============================================================ */

/* Read a line of any length from file.  Return NULL on end of file,
 * otherwise return a pointer to static storage containing the line.  Any
 * trailing newline is stripped off.  Be careful -- there is only
 * one static buffer, so subsequent calls will replace the string
 * pointed to by the return value.
 * Used in seq_lib functions and also exported
 * used in p3_read_line and seq_lib functions
 * */
#define INIT_BUF_SIZE 1024

char*
p3_read_line(FILE *file)
{
  static size_t ssz;
  static char *s = NULL;

  size_t remaining_size;
  char *p, *n;

  if (NULL == s) {
    ssz = INIT_BUF_SIZE;
    s = (char *) pr_safe_malloc(ssz);
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
    s = (char *) pr_safe_realloc(s, ssz);
    p = strchr(s, '\0');
    remaining_size = ssz - (p - s);
  }
}

/* ============================================================ */
/* END p3_read_line                                             */
/* ============================================================ */

/* ============================================================ */
/* BEGIN 'get' functions for seq_args                           */
/* ============================================================ */

const interval_array_t2 *
p3_get_sa_tar2(const seq_args *sargs) {
  return &sargs->tar2 ;
}

const interval_array_t2 *
p3_get_sa_excl2(const seq_args *sargs) {
  return &sargs->excl2 ;
}

const interval_array_t2 *
p3_get_sa_excl_internal2(const seq_args *sargs) {
  return &sargs->excl_internal2 ;
}

const interval_array_t4 *
p3_get_sa_ok_regions(const seq_args *sargs)
{
  return &sargs->ok_regions;
}

const int* 
p3_get_sa_overlap_junctions(const seq_args *sargs)
{
  return sargs->primer_overlap_junctions;
}

/* ============================================================ */
/* END 'get' functions for seq_args                             */
/* ============================================================ */

int 
p3_get_sa_n_quality(seq_args *sargs) {
  return sargs->n_quality ;
}

/* ============================================================ */
/* BEGIN 'set' functions for seq_args                           */
/* ============================================================ */

int
p3_set_sa_sequence(seq_args *sargs, const char *sequence) {
  return _set_string(&sargs->sequence, sequence) ;
}

void
p3_set_sa_primer_sequence_quality(seq_args *sargs, int quality) {
  sargs->quality[sargs->n_quality++] = quality ;
}

int
p3_set_sa_sequence_name(seq_args *sargs, const char* sequence_name) {
 return _set_string(&sargs->sequence_name, sequence_name);
}

int 
p3_set_sa_left_input(seq_args *sargs, const char *s) {
 return _set_string(&sargs->left_input, s);
}

int 
p3_set_sa_right_input(seq_args *sargs, const char *s) {
 return _set_string(&sargs->right_input, s);
}

int 
p3_set_sa_internal_input(seq_args *sargs, const char *s) {
 return _set_string(&sargs->internal_input, s);
}

void
p3_set_sa_incl_s(seq_args *sargs, int incl_s) {
  sargs->incl_s = incl_s;
}

void
p3_set_sa_incl_l(seq_args *sargs, int incl_l) {
  sargs->incl_l = incl_l;
}

void 
p3_set_sa_empty_quality(seq_args *sargs) {
  sargs->n_quality = 0;
}

void 
p3_sa_add_to_quality_array(seq_args *sargs, int quality) {
  int n = sargs->n_quality;
  if (sargs->quality_storage_size == 0) {
    sargs->quality_storage_size = 3000;
    sargs->quality
      = (int *)
      pr_safe_malloc(sizeof(*sargs->quality)
                     * sargs->quality_storage_size);
  }
  if (n > sargs->quality_storage_size) {
    sargs->quality_storage_size *= 2;
    sargs->quality
      = (int *)
      pr_safe_realloc(sargs->quality,
                      sizeof(*sargs->quality)
                      * sargs->quality_storage_size);
  }
  sargs->quality[n] = quality;
  sargs->n_quality++;
}

int
p3_sa_add_to_overlap_junctions_array(seq_args *sargs, int overlap)
{
  int c = sargs->primer_overlap_junctions_count;
  if (c >= PR_MAX_INTERVAL_ARRAY) return 1;
  sargs->primer_overlap_junctions[sargs->primer_overlap_junctions_count++] = overlap;
  return 0;
}

void
p3_set_sa_start_codon_pos(seq_args *sargs, int start_codon_pos) {
  sargs->start_codon_pos = start_codon_pos;
}

int
p3_add_to_sa_tar2(seq_args *sargs, int n1, int n2) {
  return p3_add_to_interval_array(&sargs->tar2, n1, n2);
}

int
p3_add_to_sa_excl2(seq_args *sargs, int n1, int n2) {
  return p3_add_to_interval_array(&sargs->excl2, n1, n2);
}

int
p3_add_to_sa_excl_internal2(seq_args *sargs, int n1, int n2) {
  return p3_add_to_interval_array(&sargs->excl_internal2, n1, n2);
}

int
p3_add_to_sa_ok_regions(seq_args *sargs, int l1, int l2, int r1, int r2)
{
  return p3_add_to_2_interval_array(&sargs->ok_regions, l1, l2, r1, r2);
}

/* ============================================================ */
/* END 'set' functions for seq_args                             */
/* ============================================================ */

/* ============================================================ */
/* BEGIN 'get' functions for p3_global_settings                 */
/* ============================================================ */

args_for_one_oligo_or_primer *
p3_get_gs_p_args(p3_global_settings * p) {
  return &p->p_args;
}

args_for_one_oligo_or_primer *
p3_get_gs_o_args(p3_global_settings * p) {
  return &p->o_args;
}

/* ============================================================ */
/* END 'get' functions for p3_global_settings                   */
/* ============================================================ */


/* ============================================================ */
/* BEGIN 'set' functions for p3_global_settings                 */
/* ============================================================ */

void
p3_set_gs_prmin (p3_global_settings * p , int val, int i) {
  p->pr_min[i] = val ;
}

void
p3_set_gs_prmax (p3_global_settings * p , int val, int i) {
  p->pr_max[i] = val ;
}

void
p3_set_gs_primer_opt_size(p3_global_settings * p , int val) {
  p->p_args.opt_size = val ;
}

void
p3_set_gs_primer_min_size(p3_global_settings * p , int val) {
  p->p_args.min_size = val ;
}

void
p3_set_gs_primer_max_size(p3_global_settings * p , int val) {
  p->p_args.max_size = val ;
}

void
p3_set_gs_primer_max_poly_x(p3_global_settings * p , int val) {
  p->p_args.max_poly_x = val;

}

void
p3_set_gs_primer_opt_tm(p3_global_settings * p , double d) {
  p->p_args.opt_tm = d ;
}

void
p3_set_gs_primer_opt_gc_percent(p3_global_settings * p , double d) {
  p->p_args.opt_gc_content = d ;
}

void
p3_set_gs_primer_min_tm(p3_global_settings * p , double d) {
  p->p_args.min_tm = d ;
}
void
p3_set_gs_primer_max_tm(p3_global_settings * p , double d) {
  p->p_args.max_tm = d;
}

void
p3_set_gs_primer_max_diff_tm(p3_global_settings * p , double val) {
  p->max_diff_tm = val;
}

void
p3_set_gs_primer_tm_santalucia(p3_global_settings * p,
                               tm_method_type val) {
  p->tm_santalucia = val ;
}

void
p3_set_gs_primer_salt_corrections(p3_global_settings * p,
                                  salt_correction_type val) {
  p->salt_corrections = val;
}

void
p3_set_gs_primer_min_gc(p3_global_settings * p , double val) {
  p->p_args.min_gc = val ;
}

void
p3_set_gs_primer_max_gc(p3_global_settings * p , double val) {
   p->p_args.max_gc = val ;
}

void
p3_set_gs_primer_salt_conc(p3_global_settings * p , double val) {
  p->p_args.salt_conc = val ;
}

void
p3_set_gs_primer_divalent_conc(p3_global_settings * p , double val)
{
 p->p_args.divalent_conc = val;
}

void
p3_set_gs_primer_dntp_conc(p3_global_settings * p , double val)
{
  p->p_args.dntp_conc = val ;
}

void
p3_set_gs_primer_dna_conc(p3_global_settings * p , double val)
{
  p->p_args.dna_conc = val ;
}

void
p3_set_gs_primer_num_ns_accepted(p3_global_settings * p , int val)
{
  p->p_args.num_ns_accepted = val ;
}

void
p3_set_gs_primer_product_opt_size(p3_global_settings * p , int val)
{
  p->product_opt_size = val ;
}

void
p3_set_gs_primer_self_any(p3_global_settings * p , double val)
{
  p->p_args.max_self_any = val;
}

void
p3_set_gs_primer_self_any_th(p3_global_settings * p , double val)
{
  p->p_args.max_self_any_th = val;
}

void
p3_set_gs_primer_self_end(p3_global_settings * p , double val)
{
  p->p_args.max_self_end = val; /* (short) (val * 100); */
}

void
p3_set_gs_primer_self_end_th(p3_global_settings * p , double val)
{
   p->p_args.max_self_end_th = val;
}

void
p3_set_gs_primer_hairpin_th(p3_global_settings * p , double val) {
   p->p_args.max_hairpin_th = val;
}

void   /* Called in primer3_boulder_main.c. */
p3_set_gs_primer_file_flag(p3_global_settings * p , int file_flag)
{
  p->file_flag = file_flag;
}

void
p3_set_gs_pick_anyway(p3_global_settings * p , int pick_anyway)
{
  p->pick_anyway = pick_anyway ;
}

void
p3_set_gs_gc_clamp(p3_global_settings * p , int gc_clamp)
{
  p->gc_clamp = gc_clamp;
}

void
p3_set_gs_primer_liberal_base(p3_global_settings * p , int val) {
  p->liberal_base = val;
}

void
p3_set_gs_primer_first_base_index(p3_global_settings * p , int first_base_index){
  p->first_base_index = first_base_index;
}

void
p3_set_gs_primer_num_return(p3_global_settings * p , int val) {
  p->num_return = val ;
}

void
p3_set_gs_primer_min_quality(p3_global_settings * p , int val) {
  p->p_args.min_quality = val ;
}

void
p3_set_gs_primer_min_end_quality(p3_global_settings * p , int val) {
  p->p_args.min_end_quality = val ;
}

void
p3_set_gs_primer_quality_range_min(p3_global_settings * p , int val) {
  p->quality_range_min = val ;
}

void
p3_set_gs_primer_quality_range_max(p3_global_settings * p , int val) {
  p->quality_range_max = val ;
}

void
p3_set_gs_primer_product_max_tm(p3_global_settings * p , double val) {
  p->product_max_tm = val ;
}

void
p3_set_gs_primer_product_min_tm(p3_global_settings * p , double val) {
  p->product_min_tm = val ;
}

void
p3_set_gs_primer_product_opt_tm(p3_global_settings * p , double val) {
  p->product_opt_tm = val ;
}

void
p3_set_gs_primer_task(p3_global_settings * pa , char * task_tmp)
{
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
  }
}

void
p3_set_gs_primer_pick_right_primer(p3_global_settings * p , int pick_right_primer)
{
  p->pick_right_primer = pick_right_primer;
}

void
p3_set_gs_primer_pick_internal_oligo(p3_global_settings * p , int pick_internal_oligo)
{
  p->pick_internal_oligo = pick_internal_oligo;
}

void
p3_set_gs_primer_pick_left_primer(p3_global_settings * p , int pick_left_primer) {
  p->pick_left_primer = pick_left_primer;
}

void
p3_set_gs_primer_internal_oligo_opt_size(p3_global_settings * p , int val) {
  p->o_args.opt_size = val ;
}

void
p3_set_gs_primer_internal_oligo_max_size(p3_global_settings * p , int val) {
  p->o_args.max_size = val ;
}

void
p3_set_gs_primer_internal_oligo_min_size(p3_global_settings * p , int val) {
  p->o_args.min_size = val ;
}

void
p3_set_gs_primer_internal_oligo_max_poly_x(p3_global_settings * p , int val) {
  p->o_args.max_poly_x = val ;
}

void
p3_set_gs_primer_internal_oligo_opt_tm(p3_global_settings * p , double val) {
  p->o_args.opt_tm = val ;
}

void
p3_set_gs_primer_internal_oligo_opt_gc_percent(p3_global_settings * p , double val) {
  p->o_args.opt_gc_content = val ;
}

void
p3_set_gs_primer_internal_oligo_max_tm(p3_global_settings * p , double val) {
  p->o_args.max_tm = val ;
}

void
p3_set_gs_primer_internal_oligo_min_tm(p3_global_settings * p , double val) {
  p->o_args.min_tm = val ;
}

void
p3_set_gs_primer_internal_oligo_min_gc(p3_global_settings * p , double val) {
  p->o_args.min_gc = val ;
}

void
p3_set_gs_primer_internal_oligo_max_gc(p3_global_settings * p , double val) {
  p->o_args.max_gc = val ;
}

void
p3_set_gs_primer_internal_oligo_salt_conc(p3_global_settings * p , double val) {
  p->o_args.salt_conc = val ;
}

void
p3_set_gs_primer_internal_oligo_divalent_conc(p3_global_settings * p , double val) {
  p->o_args.divalent_conc = val ;
}

void
p3_set_gs_primer_internal_oligo_dntp_conc(p3_global_settings * p , double val) {
  p->o_args.dntp_conc = val ;
}

void
p3_set_gs_primer_internal_oligo_dna_conc(p3_global_settings * p , double val) {
  p->o_args.dna_conc = val ;
}

void
p3_set_gs_primer_internal_oligo_num_ns(p3_global_settings * p , int val) {
  p->o_args.num_ns_accepted = val ;
}

void
p3_set_gs_primer_internal_oligo_min_quality(p3_global_settings * p , int val) {
  p->o_args.min_quality = val ;
}

void
p3_set_gs_primer_internal_oligo_self_any(p3_global_settings * p , double val) {
  p->o_args.max_self_any = val;
}

void
p3_set_gs_primer_internal_oligo_self_any_th(p3_global_settings * p , double val) {   
  p->o_args.max_self_any_th = val;
}

void
p3_set_gs_primer_internal_oligo_self_end(p3_global_settings * p , double val) {
  p->o_args.max_self_end = val; /* (short) (val * 100); */
}

void
p3_set_gs_primer_internal_oligo_self_end_th(p3_global_settings * p , double val) {  
   p->o_args.max_self_end_th = val;
}

void
p3_set_gs_primer_internal_oligo_hairpin_th(p3_global_settings * p , double val) {
      p->o_args.max_hairpin_th = val;
}


void
p3_set_gs_primer_max_mispriming(p3_global_settings * p , double val) {
  p->p_args.max_repeat_compl = val;
}

void
p3_set_gs_primer_internal_oligo_max_mishyb(p3_global_settings * p , double val) {
  p->o_args.max_repeat_compl = val;
}

void
p3_set_gs_primer_pair_max_mispriming(p3_global_settings * p , double val) {
  p->pair_repeat_compl = val;
}

void
p3_set_gs_primer_max_template_mispriming(p3_global_settings * p , double val) {
  p->p_args.max_template_mispriming = val;
}

void
p3_set_gs_primer_max_template_mispriming_th(p3_global_settings * p , double val) {  
  p->p_args.max_template_mispriming_th = val;
}

void
p3_set_gs_primer_lib_ambiguity_codes_consensus(p3_global_settings * p , int val) {
  p->lib_ambiguity_codes_consensus = val ;
}

void
p3_set_gs_primer_inside_penalty(p3_global_settings * p , double val) {
  p->inside_penalty = val ;
}

void
p3_set_gs_primer_outside_penalty(p3_global_settings * p , double val) {
  p->outside_penalty = val ;
}

void
p3_set_gs_primer_mispriming_library(p3_global_settings * p , char * path) {
  /* TO DO: Re-factor this? TO DO: Check for errors */
  p->p_args.repeat_lib = read_and_create_seq_lib(path, "mispriming library") ;
}

void
p3_set_gs_primer_internal_oligo_mishyb_library(p3_global_settings * p , char * path) {
  /* TO DO: Re-factor this? TO DO: Check for errors */
  p->o_args.repeat_lib = read_and_create_seq_lib(path, "internal oligio mishyb library") ;
}

void
p3_set_gs_primer_max_end_stability(p3_global_settings * p , double val) {
  p->max_end_stability = val ;
}

void
p3_set_gs_primer_lowercase_masking(p3_global_settings * p , int val) {
  p->lowercase_masking = val ;
}

void
p3_set_gs_primer_thermodynamic_oligo_alignment(p3_global_settings * p , int val) {   
   p->thermodynamic_oligo_alignment = val ;
}

void
p3_set_gs_primer_thermodynamic_template_alignment(p3_global_settings * p , int val) {   
   p->thermodynamic_template_alignment = val ;
}

void
p3_set_gs_primer_wt_tm_gt(p3_global_settings * p , double val) {
  p->p_args.weights.temp_gt = val ;
}

void
p3_set_gs_primer_wt_tm_lt(p3_global_settings * p , double val) {
  p->p_args.weights.temp_lt = val ;
}

void
p3_set_gs_primer_wt_gc_percent_gt(p3_global_settings * p , double val) {
  p->p_args.weights.gc_content_gt = val ;
}

void
p3_set_gs_primer_wt_gc_percent_lt(p3_global_settings * p , double val) {
  p->p_args.weights.gc_content_lt = val ;
}

void
p3_set_gs_primer_wt_size_lt(p3_global_settings * p , double val) {
  p->p_args.weights.length_lt = val ;
}

void
p3_set_gs_primer_wt_size_gt(p3_global_settings * p , double val) {
  p->p_args.weights.length_gt = val ;
}

void
p3_set_gs_primer_wt_compl_any(p3_global_settings * p , double val) {
  p->p_args.weights.compl_any = val ;
}

void
p3_set_gs_primer_wt_compl_any_th(p3_global_settings * p , double val) {   
   p->p_args.weights.compl_any_th = val ;
}

void
p3_set_gs_primer_wt_compl_end(p3_global_settings * p , double val) {
  p->p_args.weights.compl_end = val ;
}

void
p3_set_gs_primer_wt_compl_end_th(p3_global_settings * p , double val) {   
     p->p_args.weights.compl_end_th = val ;
}

void
p3_set_gs_primer_wt_hairpin_th(p3_global_settings * p , double val) {
   p->p_args.weights.hairpin_th = val ;
}

void
p3_set_gs_primer_wt_num_ns(p3_global_settings * p , double val) {
  p->p_args.weights.num_ns = val ;
}

void
p3_set_gs_primer_wt_rep_sim(p3_global_settings * p , double val) {
  p->p_args.weights.repeat_sim = val ;
}

void
p3_set_gs_primer_wt_seq_qual(p3_global_settings * p , double val) {
  p->p_args.weights.seq_quality = val ;
}

void
p3_set_gs_primer_wt_end_qual(p3_global_settings * p , double val) {
  p->p_args.weights.end_quality = val ;
}

void
p3_set_gs_primer_wt_pos_penalty(p3_global_settings * p , double val){
  p->p_args.weights.pos_penalty = val ;
}

void
p3_set_gs_primer_wt_end_stability(p3_global_settings * p , double val) {
  p->p_args.weights.end_stability = val ;
}

void
p3_set_gs_primer_wt_template_mispriming(p3_global_settings * p , double val) {
  p->p_args.weights.template_mispriming = val ;
}

void
p3_set_gs_primer_wt_template_mispriming_th(p3_global_settings * p , double val) {   
     p->p_args.weights.template_mispriming_th = val ;
}

void
p3_set_gs_primer_io_wt_tm_gt(p3_global_settings * p , double val) {
  p->o_args.weights.temp_gt = val ;
}

void
p3_set_gs_primer_io_wt_tm_lt(p3_global_settings * p , double val){
  p->o_args.weights.temp_lt = val ;
}

void
p3_set_gs_primer_io_wt_gc_percent_gt(p3_global_settings * p , double val) {
  p->o_args.weights.gc_content_gt = val ;
}

void
p3_set_gs_primer_io_wt_gc_percent_lt(p3_global_settings * p , double val) {
  p->o_args.weights.gc_content_lt = val ;

}

void
p3_set_gs_primer_io_wt_size_lt(p3_global_settings * p , double val) {
  p->o_args.weights.length_lt = val ;
}

void
p3_set_gs_primer_io_wt_size_gt(p3_global_settings * p , double val) {
  p->o_args.weights.length_gt = val ;
}

void
p3_set_gs_primer_io_wt_wt_compl_any(p3_global_settings * p , double val) {
  p->o_args.weights.compl_any = val ;
}

void
p3_set_gs_primer_io_wt_wt_compl_any_th(p3_global_settings * p , double val) {
     p->o_args.weights.compl_any_th = val ;
}

void
p3_set_gs_primer_io_wt_compl_end(p3_global_settings * p , double val) {
  p->o_args.weights.compl_end = val ;
}

void
p3_set_gs_primer_io_wt_compl_end_th(p3_global_settings * p , double val) {
     p->o_args.weights.compl_end_th = val ;
}

void
p3_set_gs_primer_io_wt_hairpin_th(p3_global_settings * p , double val) {   
   p->o_args.weights.hairpin_th = val ;
}

void
p3_set_gs_primer_io_wt_num_ns(p3_global_settings * p , double val) {
  p->o_args.weights.num_ns = val ;
}

void
p3_set_gs_primer_io_wt_rep_sim(p3_global_settings * p , double val) {
  p->o_args.weights.repeat_sim = val ;
}

void
p3_set_gs_primer_io_wt_seq_qual(p3_global_settings * p , double val) {
  p->o_args.weights.seq_quality = val ;
}

void
p3_set_gs_primer_io_wt_end_qual(p3_global_settings * p , double val) {
  p->o_args.weights.end_quality = val ;
}

void
p3_set_gs_primer_pair_wt_pr_penalty(p3_global_settings * p , double val) {
  p->pr_pair_weights.primer_quality = val ;
}

void
p3_set_gs_primer_pair_wt_io_penalty(p3_global_settings * p , double val) {
  p->pr_pair_weights.io_quality = val ;
}

void
p3_set_gs_primer_pair_wt_diff_tm(p3_global_settings * p , double val) {
  p->pr_pair_weights.diff_tm = val ;
}

void
p3_set_gs_primer_pair_wt_compl_any(p3_global_settings * p , double val) {
  p->pr_pair_weights.compl_any = val ;
}

void
p3_set_gs_primer_pair_wt_compl_any_th(p3_global_settings * p , double val) {   
   p->pr_pair_weights.compl_any_th = val ;
}

void
p3_set_gs_primer_pair_wt_compl_end(p3_global_settings * p , double val) {
  p->pr_pair_weights.compl_end = val ;
}

void
p3_set_gs_primer_pair_wt_compl_end_th(p3_global_settings * p , double val) {   
   p->pr_pair_weights.compl_end_th = val ;
}

void
p3_set_gs_primer_pair_wt_product_tm_lt(p3_global_settings * p , double val) {
  p->pr_pair_weights.product_tm_lt = val ;
}

void
p3_set_gs_primer_pair_wt_product_tm_gt(p3_global_settings * p , double val) {
  p->pr_pair_weights.product_tm_gt = val ;
}

void
p3_set_gs_primer_pair_wt_product_size_gt(p3_global_settings * p , double val) {
  p->pr_pair_weights.product_size_gt = val ;
}

void
p3_set_gs_primer_pair_wt_product_size_lt(p3_global_settings * p , double val) {
  p->pr_pair_weights.product_size_lt = val ;
}

void
p3_set_gs_primer_pair_wt_rep_sim(p3_global_settings * p , double val) {
  p->pr_pair_weights.repeat_sim = val ;
}

void
p3_set_gs_primer_pair_wt_template_mispriming(p3_global_settings * p , double val) {
  p->pr_pair_weights.template_mispriming = val ;
}

void
p3_set_gs_primer_pair_wt_template_mispriming_th(p3_global_settings * p , double val){   
     p->pr_pair_weights.template_mispriming_th = val ;
}

void
p3_set_gs_lib_ambiguity_codes_consensus(p3_global_settings * p,
                                        int lib_ambiguity_codes_consensus)
{
  p->lib_ambiguity_codes_consensus = lib_ambiguity_codes_consensus;
}

void
p3_set_gs_quality_range_min(p3_global_settings * p , int quality_range_min){
  p->quality_range_min = quality_range_min;
}

void
p3_set_gs_quality_range_max(p3_global_settings * p , int quality_range_max){
  p->quality_range_max = quality_range_max ;
}

void 
p3_set_gs_max_end_gc(p3_global_settings *p, int max_end_gc) {
  p->max_end_gc = max_end_gc;
}

void
p3_set_gs_max_end_stability(p3_global_settings * p , int max_end_stability){
  p->max_end_stability = max_end_stability;
}

void
p3_set_gs_lowercase_masking(p3_global_settings * p , int lowercase_masking){
  p->lowercase_masking = lowercase_masking;
}

void
p3_set_gs_outside_penalty(p3_global_settings * p , double outside_penalty){
  p->outside_penalty = outside_penalty;
}

void
p3_set_gs_inside_penalty(p3_global_settings * p , double inside_penalty){
  p->inside_penalty = inside_penalty;
}

void
p3_set_gs_pair_max_template_mispriming(p3_global_settings * p,
                                       double  pair_max_template_mispriming)
{
  p->pair_max_template_mispriming = pair_max_template_mispriming;
}

void
p3_set_gs_pair_max_template_mispriming_th(p3_global_settings * p,
                                          double  pair_max_template_mispriming_th)
{
  p->pair_max_template_mispriming_th = pair_max_template_mispriming_th;
}


void
p3_set_gs_pair_repeat_compl(p3_global_settings * p, double pair_repeat_compl)
{
  p->pair_repeat_compl = pair_repeat_compl;
}

void
p3_set_gs_pair_compl_any(p3_global_settings * p , double pair_compl_any)
{
  p->pair_compl_any = pair_compl_any;
}

void
p3_set_gs_pair_compl_any_th(p3_global_settings * p , double pair_compl_any_th) 
{
  p->pair_compl_any_th = pair_compl_any_th;
}

void
p3_set_gs_pair_compl_end(p3_global_settings * p , double  pair_compl_end)
{
  p->pair_compl_end = pair_compl_end;
}

void
p3_set_gs_pair_compl_end_th(p3_global_settings * p , double  pair_compl_end_th) 
{   
  p->pair_compl_end = pair_compl_end_th;
}

void
p3_set_gs_min_left_three_prime_distance(p3_global_settings *p, int min_distance) 
{
  p->min_left_three_prime_distance = min_distance;
}

void
p3_set_gs_min_right_three_prime_distance(p3_global_settings *p, int min_distance) 
{
  p->min_right_three_prime_distance = min_distance;
}

void 
p3_set_gs_min_5_prime_overlap_of_junction(p3_global_settings *p, int min_5_prime)
{
  p->min_5_prime_overlap_of_junction = min_5_prime;
}

void 
p3_set_gs_min_3_prime_overlap_of_junction(p3_global_settings *p, int min_3_prime)
{
  p->min_3_prime_overlap_of_junction = min_3_prime;
}

void
p3_set_gs_max_diff_tm(p3_global_settings * p , double max_diff_tm) {
  p->max_diff_tm = max_diff_tm;
}

void
p3_set_gs_primer_pick_anyway(p3_global_settings * p , int val) {
  p->pick_anyway = val ;
}

void
p3_set_gs_primer_gc_clamp(p3_global_settings * p , int val) {
  p->gc_clamp = val;
}

void
p3_empty_gs_product_size_range(p3_global_settings *pgs) {
  pgs->num_intervals = 0;
}

int
p3_add_to_gs_product_size_range(p3_global_settings *pgs,
                                int n1, int n2) {
  if (pgs->num_intervals >= PR_MAX_INTERVAL_ARRAY)
    return 1;
  pgs->pr_min[pgs->num_intervals]  = n1;
  pgs->pr_max[pgs->num_intervals]  = n2;
  pgs->num_intervals++;
  return 0;
}

/* ============================================================ */
/* END 'set' functions for p3_global_settings                   */
/* ============================================================ */


/* ============================================================ */
/* START functions for setting and getting oligo problems       */
/* ============================================================ */

static void
initialize_op(primer_rec *oligo) {
  oligo->problems.prob = 0UL;  /* 0UL is the unsigned long zero */
}

/* We use bitfields to store all primer data. The idea is that 
 * a primer is uninitialized if all bits are 0. It was evaluated 
 * if OP_PARTIALLY_WRITTEN is true. If OP_COMPLETELY_WRITTEN is
 * also true the primer was evaluated till the end - meaning if 
 * all OP_... are false, the primer is OK, if some are true 
 * primer3 was forced to use it (must_use). */


#define OP_PARTIALLY_WRITTEN                   (1UL <<  0)
#define OP_COMPLETELY_WRITTEN                  (1UL <<  1)
#define BF_OVERLAPS_TARGET                     (1UL <<  2)
#define BF_OVERLAPS_EXCL_REGION                (1UL <<  3)
#define BF_INFINITE_POSITION_PENALTY           (1UL <<  4)
/* Space for more bitfields */

#define OP_NOT_IN_ANY_OK_REGION                (1UL <<  7)
#define OP_TOO_MANY_NS                         (1UL <<  8) /* 3prime problem*/
#define OP_OVERLAPS_TARGET                     (1UL <<  9) /* 3prime problem*/
#define OP_HIGH_GC_CONTENT                     (1UL << 10)
#define OP_LOW_GC_CONTENT                      (1UL << 11)
#define OP_HIGH_TM                             (1UL << 12)
#define OP_LOW_TM                              (1UL << 13)
#define OP_OVERLAPS_EXCL_REGION                (1UL << 14) /* 3prime problem*/
#define OP_HIGH_SELF_ANY                       (1UL << 15) /* 3prime problem*/
#define OP_HIGH_SELF_END                       (1UL << 16)
#define OP_NO_GC_CLAMP                         (1UL << 17) /* 3prime problem*/
#define OP_HIGH_END_STABILITY                  (1UL << 18) /* 3prime problem*/
#define OP_HIGH_POLY_X                         (1UL << 19) /* 3prime problem*/
#define OP_LOW_SEQUENCE_QUALITY                (1UL << 20) /* 3prime problem*/
#define OP_LOW_END_SEQUENCE_QUALITY            (1UL << 21) /* 3prime problem*/
#define OP_HIGH_SIM_TO_NON_TEMPLATE_SEQ        (1UL << 22) /* 3prime problem*/
#define OP_HIGH_SIM_TO_MULTI_TEMPLATE_SITES    (1UL << 23)
#define OP_OVERLAPS_MASKED_SEQ                 (1UL << 24)
#define OP_TOO_LONG                            (1UL << 25) /* 3prime problem*/
#define OP_TOO_SHORT                           (1UL << 26)
#define OP_DOES_NOT_AMPLIFY_ORF                (1UL << 27)
#define OP_TOO_MANY_GC_AT_END                  (1UL << 28) /* 3prime problem*/
#define OP_HIGH_HAIRPIN                        (1UL << 29) /* 3prime problem*/
#define OP_MUST_MATCH_ERR                      (1UL << 30)

/* Space for more Errors */

/* Tip: the calculator of windows (in scientific mode) can easily convert 
   between decimal and binary numbers */

/* all bits 1 except bits 0 to 6 */
/* (~0UL) ^ 127UL = 1111 1111  1111 1111  1111 1111  1000 0000 */
static unsigned long int any_problem = (~0UL) ^ 127UL;
/* 310297344UL = 0001 0010  0111 1110  1100 0011  0000 0000 */
static unsigned long int five_prime_problem = 310297344UL;

int
p3_ol_is_uninitialized(const primer_rec *oligo) {
  return (oligo->problems.prob == 0UL);
}

int
p3_ol_is_ok(const primer_rec *oligo) {
  return (oligo->problems.prob & OP_COMPLETELY_WRITTEN) != 0;
}
 
/* 
   Return 1 iff the argument 'oligo' has any problems -- i.e.
   violations of design constraints.
 */
int
p3_ol_has_any_problem(const primer_rec *oligo) {
  return (oligo->problems.prob & any_problem) != 0;
}

static int
any_5_prime_ol_extension_has_problem(const primer_rec *oligo) {
  return (oligo->problems.prob & five_prime_problem) != 0;
}

/*  
    Return a string details the the problems in 'oligo', i.e. the
    constraints that 'oligo' violates.  WARNING: Returns a pointer to
    static storage, which is over-written on next call.
 */
#define ADD_OP_STR(COND, STR) if (prob & COND) { tmp = STR; strcpy(next, tmp); next += strlen(tmp); }
const char *
p3_get_ol_problem_string(const primer_rec *oligo) {
  static char output[4096*2];  /* WARNING, make sure all error strings
                                  concatenated will not overflow this
                                  buffer. */
  char *next = output;
  const char *tmp;
  unsigned long int prob = oligo->problems.prob;
  output[0] = '\0';

  if (
      (prob & OP_PARTIALLY_WRITTEN)
      &&
      !(prob & OP_COMPLETELY_WRITTEN)
      ) {
    tmp = " Not completely checked;";
    strcpy(next, tmp);
    next += strlen(tmp);
  }

  {
    ADD_OP_STR(OP_TOO_MANY_NS,
               " Too many Ns;")
    ADD_OP_STR(OP_OVERLAPS_TARGET,
               " Overlaps target;")
    ADD_OP_STR(OP_HIGH_GC_CONTENT,
               " GC content too high;")
    ADD_OP_STR(OP_LOW_GC_CONTENT,
               " GC content too low;")
    ADD_OP_STR(OP_HIGH_TM,
               " Temperature too high;")
    ADD_OP_STR(OP_LOW_TM,
               " Temperature too low;")
    ADD_OP_STR(OP_OVERLAPS_EXCL_REGION,
               " Overlaps an excluded region;")
    ADD_OP_STR(OP_NOT_IN_ANY_OK_REGION,
               " Not in any ok region;")
    ADD_OP_STR(OP_HIGH_SELF_ANY,
               " Similarity to self too high;")
    ADD_OP_STR(OP_HIGH_SELF_END,
               " Similary to 3' end of self too high;")
    ADD_OP_STR(OP_HIGH_HAIRPIN,
               " Hairpin stability too high;")
    ADD_OP_STR(OP_NO_GC_CLAMP,
               " No 3' GC clamp;")
    ADD_OP_STR(OP_TOO_MANY_GC_AT_END,
               " Too many GCs at 3' end;")
    ADD_OP_STR(OP_HIGH_END_STABILITY,
               " 3' end too stable (delta-G too high);")
    ADD_OP_STR(OP_HIGH_POLY_X,
               " Contains too-long poly nucleotide tract;")
    ADD_OP_STR(OP_LOW_SEQUENCE_QUALITY,
               " Template sequence quality too low;")
    ADD_OP_STR(OP_LOW_END_SEQUENCE_QUALITY,
               " Template sequence quality at 3' end too low;")
    ADD_OP_STR(OP_HIGH_SIM_TO_NON_TEMPLATE_SEQ,
               " Similarity to non-template sequence too high;")
    ADD_OP_STR(OP_HIGH_SIM_TO_MULTI_TEMPLATE_SITES,
               " Similarity to multiple sites in template;")
    ADD_OP_STR(OP_OVERLAPS_MASKED_SEQ,
               " 3' base overlaps masked sequence;")
    ADD_OP_STR(OP_TOO_LONG,
               " Too long;")
    ADD_OP_STR(OP_TOO_SHORT,
               " Too short;")
    ADD_OP_STR(OP_DOES_NOT_AMPLIFY_ORF,
               " Would not amplify an open reading frame;")
    ADD_OP_STR(OP_MUST_MATCH_ERR,
               " Failed must_match requirements;")
  }
  return output;
}
#undef ADD_OP_STR

static void 
bf_set_overlaps_target(primer_rec *oligo, int val){
  if (val == 0) {
    oligo->problems.prob |= BF_OVERLAPS_TARGET;
    oligo->problems.prob ^= BF_OVERLAPS_TARGET;  
  } else {
    oligo->problems.prob |= BF_OVERLAPS_TARGET;
  }
}

static int
bf_get_overlaps_target(const primer_rec *oligo) {
  return (oligo->problems.prob & BF_OVERLAPS_TARGET) != 0;
}

static void 
bf_set_overlaps_excl_region(primer_rec *oligo, int val){
  if (val == 0) {
    oligo->problems.prob |= BF_OVERLAPS_EXCL_REGION;
    oligo->problems.prob ^= BF_OVERLAPS_EXCL_REGION;  
  } else {
    oligo->problems.prob |= BF_OVERLAPS_EXCL_REGION;
  }
}

static int
bf_get_overlaps_excl_region(const primer_rec *oligo) {
  return (oligo->problems.prob & BF_OVERLAPS_EXCL_REGION) != 0;
}

static void 
bf_set_infinite_pos_penalty(primer_rec *oligo, int val){
  if (val == 0) {
    oligo->problems.prob |= BF_INFINITE_POSITION_PENALTY;
    oligo->problems.prob ^= BF_INFINITE_POSITION_PENALTY;  
  } else {
    oligo->problems.prob |= BF_INFINITE_POSITION_PENALTY;
  }
}

static int
bf_get_infinite_pos_penalty(const primer_rec *oligo) {
  return (oligo->problems.prob & BF_INFINITE_POSITION_PENALTY) != 0;
}

static void
op_set_does_not_amplify_orf(primer_rec *oligo) {
  oligo->problems.prob |= OP_DOES_NOT_AMPLIFY_ORF;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

static void
op_set_completely_written(primer_rec *oligo) {
  oligo->problems.prob |= OP_COMPLETELY_WRITTEN;
}

static void
op_set_must_match_err(primer_rec *oligo) {
  oligo->problems.prob |= OP_MUST_MATCH_ERR;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

static void
op_set_too_many_ns(primer_rec *oligo) {
  oligo->problems.prob |= OP_TOO_MANY_NS;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

static void
op_set_overlaps_target(primer_rec *oligo) {
  oligo->problems.prob |= OP_OVERLAPS_TARGET;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

static void
op_set_high_gc_content(primer_rec *oligo) {
  oligo->problems.prob |= OP_HIGH_GC_CONTENT;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

static void
op_set_low_gc_content(primer_rec *oligo) {
  oligo->problems.prob |= OP_LOW_GC_CONTENT;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

static void
op_set_high_tm(primer_rec *oligo) {
  oligo->problems.prob |= OP_HIGH_TM;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

static void
op_set_low_tm(primer_rec *oligo) {
  oligo->problems.prob |= OP_LOW_TM;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

static void
op_set_overlaps_excluded_region(primer_rec *oligo) {
  oligo->problems.prob |= OP_OVERLAPS_EXCL_REGION;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

static void
op_set_not_in_any_ok_region(primer_rec *oligo) {
  oligo->problems.prob |= OP_NOT_IN_ANY_OK_REGION;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

static void
op_set_high_self_any(primer_rec *oligo) {
  oligo->problems.prob |= OP_HIGH_SELF_ANY;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

static void
op_set_high_self_end(primer_rec *oligo) {
  oligo->problems.prob |= OP_HIGH_SELF_END;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

static void
op_set_high_hairpin_th(primer_rec *oligo) {
   oligo->problems.prob |= OP_HIGH_HAIRPIN;
   oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

static void
op_set_no_gc_glamp(primer_rec *oligo) {
  oligo->problems.prob |= OP_NO_GC_CLAMP;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

static void
op_set_too_many_gc_at_end(primer_rec *oligo) {
  oligo->problems.prob |= OP_TOO_MANY_GC_AT_END;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

/* Must not be called on a hybridization probe / internal oligo */
static void
op_set_high_end_stability(primer_rec *oligo) {
  oligo->problems.prob |= OP_HIGH_END_STABILITY;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

static void
op_set_high_poly_x(primer_rec *oligo) {
  oligo->problems.prob |= OP_HIGH_POLY_X;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

static void
op_set_low_sequence_quality(primer_rec *oligo) {
  oligo->problems.prob |= OP_LOW_SEQUENCE_QUALITY;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

static void
op_set_low_end_sequence_quality(primer_rec *oligo) {
  oligo->problems.prob |= OP_LOW_END_SEQUENCE_QUALITY;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

static void
op_set_high_similarity_to_non_template_seq(primer_rec *oligo) {
  oligo->problems.prob |= OP_HIGH_SIM_TO_NON_TEMPLATE_SEQ;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

static void
op_set_high_similarity_to_multiple_template_sites(primer_rec *oligo) {
  oligo->problems.prob |= OP_HIGH_SIM_TO_MULTI_TEMPLATE_SITES;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}
   
static void
op_set_overlaps_masked_sequence(primer_rec *oligo) {
  oligo->problems.prob |= OP_OVERLAPS_MASKED_SEQ;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

static void
op_set_too_long(primer_rec *oligo) {
  oligo->problems.prob |= OP_TOO_LONG;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

static void
op_set_too_short(primer_rec *oligo) {
  oligo->problems.prob |= OP_TOO_SHORT;
  oligo->problems.prob |= OP_PARTIALLY_WRITTEN;
}

/* ============================================================ */
/* END functions for setting and getting oligo problems         */
/* ============================================================ */


/* ============================================================ */
/* START functions for getting values from p3retvals            */
/* ============================================================ */

char *
p3_get_rv_and_gs_warnings(const p3retval *retval,
                          const p3_global_settings *pa) {

  pr_append_str warning;

  PR_ASSERT(NULL != pa);

  init_pr_append_str(&warning);

  if (seq_lib_warning_data(pa->p_args.repeat_lib))
    pr_append_new_chunk(&warning, seq_lib_warning_data(pa->p_args.repeat_lib));

  if(seq_lib_warning_data(pa->o_args.repeat_lib)) {
    pr_append_new_chunk(&warning, seq_lib_warning_data(pa->o_args.repeat_lib));
    pr_append(&warning, " (for internal oligo)");
  }

  if (!pr_is_empty(&retval->warnings))
    pr_append_new_chunk(&warning,  retval->warnings.data);

  return pr_is_empty(&warning) ? NULL : warning.data;
}

const char *
p3_get_rv_global_errors(const p3retval *retval) {
  return retval->glob_err.data;
}

const char *
p3_get_rv_per_sequence_errors(const p3retval *retval) {
  return retval->per_sequence_err.data;
}

p3_output_type
p3_get_rv_output_type(const p3retval *r) {
  return r->output_type;
}

const char *
p3_get_rv_warnings(const p3retval *r) {
  return pr_append_str_chars(&r->warnings);
}

int
p3_get_rv_stop_codon_pos(p3retval *r) {
  return r->stop_codon_pos;
}

/* ============================================================ */
/* END functions for getting values from p3retvals            */
/* ============================================================ */

void 
p3_print_args(const p3_global_settings *p, seq_args *s)
{
  int i;

  if (p != NULL) {
    printf("=============\n");
    printf("BEGIN GLOBAL ARGS\n") ;
    printf("  primer_task %i\n", p->primer_task);
    printf("  pick_left_primer %i\n", p->pick_left_primer);
    printf("  pick_right_primer %i\n", p->pick_right_primer);
    printf("  pick_internal_oligo %i\n", p->pick_internal_oligo);
    printf("  file_flag %i\n", p->file_flag) ;
    printf("  first_base_index %i\n", p->first_base_index);
    printf("  liberal_base %i\n", p->liberal_base );
    printf("  num_return %i\n", p->num_return) ;
    printf("  pick_anyway %i\n", p->pick_anyway);
    printf("  lib_ambiguity_codes_consensus %i\n",
           p->lib_ambiguity_codes_consensus) ;
    printf("  quality_range_min %i\n", p->quality_range_min) ;
    printf("  quality_range_max %i\n", p->quality_range_max) ;

    printf("  tm_santalucia %i\n", p->tm_santalucia) ;
    printf("  salt_corrections %i\n", p->salt_corrections) ;
    printf("  max_end_stability %f\n", p->max_end_stability) ;
    printf("  gc_clamp %i\n", p->gc_clamp) ;
    printf("  max_end_gc %i\n", p->max_end_gc);
    printf("  lowercase_masking %i\n", p->lowercase_masking) ;
    printf("  thermodynamic_oligo_alignment %i\n", p->thermodynamic_oligo_alignment);
    printf("  thermodynamic_template_alignment %i\n", p->thermodynamic_template_alignment);
    printf("  outside_penalty %f\n", p->outside_penalty) ;
    printf("  inside_penalty %f\n", p->inside_penalty) ;
    printf("  number of product size ranges: %d\n", p->num_intervals);
    printf("  product size ranges:\n");
    for (i = 0; i < p->num_intervals; i++) {
      printf("  %d - %d \n", p->pr_min[i], p->pr_max[i]);
    }
    printf("  product_opt_size %i\n", p->product_opt_size) ;
    printf("  product_max_tm %f\n", p->product_max_tm) ;
    printf("  product_min_tm %f\n", p->product_min_tm) ;
    printf("  product_opt_tm %f\n", p->product_opt_tm) ;
    printf("  pair_max_template_mispriming %f\n", p->pair_max_template_mispriming) ;
    printf("  pair_max_template_mispriming_th %f\n", p->pair_max_template_mispriming_th) ;
    printf("  pair_repeat_compl %f\n", p->pair_repeat_compl) ;
    printf("  pair_compl_any %f\n", p->pair_compl_any) ;
    printf("  pair_compl_end %f\n", p->pair_compl_end) ;
    printf("  pair_compl_any_th %f\n", p->pair_compl_any_th) ;
    printf("  pair_compl_end_th %f\n", p->pair_compl_end_th) ;
     
    printf("  min_left_three_prime_distance %i\n", p->min_left_three_prime_distance) ;
    printf("  min_right_three_prime_distance %i\n", p->min_right_three_prime_distance) ;
    printf("  min_5_prime_overlap_of_junction %i\n", p->min_5_prime_overlap_of_junction);
    printf("  min_3_prime_overlap_of_junction %i\n", p->min_3_prime_overlap_of_junction);
    printf("  dump %i\n", p->dump);

    printf("  begin pr_pair_weights\n") ;
    printf("    primer_quality %f\n", p->pr_pair_weights.primer_quality) ;
    printf("    io_quality %f\n", p->pr_pair_weights.io_quality) ;
    printf("    diff_tm %f\n", p->pr_pair_weights.diff_tm) ;
    printf("    compl_any %f\n", p->pr_pair_weights.compl_any) ;
    printf("    compl_end %f\n", p->pr_pair_weights.compl_end) ;
    printf("    compl_any_th %f\n", p->pr_pair_weights.compl_any_th) ;
    printf("    compl_end_th %f\n", p->pr_pair_weights.compl_end_th) ;
    printf("    product_tm_lt %f\n", p->pr_pair_weights.product_tm_lt) ;
    printf("    product_tm_gt %f\n", p->pr_pair_weights.product_tm_gt) ;
    printf("    product_size_lt %f\n", p->pr_pair_weights.product_size_lt) ;
    printf("    product_size_gt %f\n", p->pr_pair_weights.product_size_gt) ;
    printf("    repeat_sim %f\n", p->pr_pair_weights.repeat_sim) ;
    printf("    template_mispriming %f\n", p->pr_pair_weights.template_mispriming) ;
    printf("    template_mispriming_th %f\n", p->pr_pair_weights.template_mispriming_th) ;
    printf("  end pair_weights\n") ;


    printf("\n\n");
    printf("=============\n");
    printf("BEGIN primer_args\n");
    printf("begin oligo_weights\n");
    printf("temp_gt %f\n", p->p_args.weights.temp_gt ) ;
    printf("temp_gt %f\n", p->p_args.weights.temp_gt) ;
    printf("temp_lt %f\n", p->p_args.weights.temp_lt) ;
    printf("gc_content_gt %f\n", p->p_args.weights.gc_content_gt) ;
    printf("gc_content_lt %f\n", p->p_args.weights.gc_content_lt) ;
    printf("compl_any %f\n", p->p_args.weights.compl_any) ;
    printf("compl_end %f\n", p->p_args.weights.compl_end) ;
    printf("compl_any_th %f\n", p->p_args.weights.compl_any_th) ;
    printf("compl_end_th %f\n", p->p_args.weights.compl_end_th) ;
    printf("hairpin %f\n", p->p_args.weights.hairpin_th) ;
    printf("repeat_sim %f\n", p->p_args.weights.repeat_sim) ;
    printf("length_lt %f\n", p->p_args.weights.length_lt) ;
    printf("length_gt %f\n", p->p_args.weights.length_gt) ;
    printf("seq_quality %f\n", p->p_args.weights.seq_quality) ;
    printf("end_quality %f\n", p->p_args.weights.end_quality) ;
    printf("pos_penalty %f\n", p->p_args.weights.pos_penalty) ;
    printf("end_stability %f\n", p->p_args.weights.end_stability) ;
    printf("num_ns %f\n", p->p_args.weights.num_ns) ;
    printf("template_mispriming %f\n", p->p_args.weights.template_mispriming) ;
    printf("template_mispriming_th %f\n", p->p_args.weights.template_mispriming_th) ;
    printf("end oligo_weights\n") ;

    printf("opt_tm %f\n", p->p_args.opt_tm) ;
    printf("min_tm %f\n", p->p_args.min_tm) ;
    printf("max_tm %f\n", p->p_args.max_tm) ;
    printf("opt_gc_content %f\n", p->p_args.opt_gc_content) ;
    printf("max_gc %f\n", p->p_args.max_gc) ;
    printf("min_gc %f\n", p->p_args.min_gc) ;
    printf("divalent_conc %f\n", p->p_args.divalent_conc) ;
    printf("dntp_conc %f\n", p->p_args.dntp_conc) ;
    printf("dna_conc %f\n", p->p_args.dna_conc) ;
    printf("num_ns_accepted %i\n", p->p_args.num_ns_accepted) ;
    printf("opt_size %i\n", p->p_args.opt_size) ;
    printf("min_size %i\n", p->p_args.min_size) ;
    printf("max_size %i\n", p->p_args.max_size) ;
    printf("max_poly_x %i\n", p->p_args.max_poly_x) ;
    printf("min_end_quality %i\n", p->p_args.min_end_quality) ;
    printf("min_quality %i\n", p->p_args.min_quality) ;
    printf("max_self_any %f\n", p->p_args.max_self_any) ;
    printf("max_self_end %f\n", p->p_args.max_self_end) ;
    printf("max_self_any_th %f\n", p->p_args.max_self_any_th) ;
    printf("max_self_end_th %f\n", p->p_args.max_self_end_th) ;
    printf("max_hairpin %f\n", p->p_args.max_hairpin_th) ;
    printf("max_repeat_compl %f\n", p->p_args.max_repeat_compl) ;
    printf("max_template_mispriming %f\n", p->p_args.max_template_mispriming) ;
    printf("max_template_mispriming_th %f\n", p->p_args.max_template_mispriming_th) ;
    printf("end primer args\n") ;

    printf("begin internal oligo args (p->o_args.)\n") ;

    printf("  begin internal oligo_weights (p->o_args.weights.)\n") ;
    printf("    temp_gt %f\n", p->o_args.weights.temp_gt) ;
    printf("    temp_lt %f\n", p->o_args.weights.temp_lt) ;
    printf("    gc_content_gt %f\n", p->o_args.weights.gc_content_gt) ;
    printf("    gc_content_lt %f\n", p->o_args.weights.gc_content_lt) ;
    printf("    compl_any %f\n", p->o_args.weights.compl_any) ;
    printf("    compl_end %f\n", p->o_args.weights.compl_end) ;
    printf("    compl_any_th %f\n", p->o_args.weights.compl_any_th) ;
    printf("    compl_end_th %f\n", p->o_args.weights.compl_end_th) ;
    printf("    hairpin %f\n", p->o_args.weights.hairpin_th) ;
    printf("    repeat_sim %f\n", p->o_args.weights.repeat_sim) ;
    printf("    length_lt %f\n", p->o_args.weights.length_lt) ;
    printf("    length_gt %f\n", p->o_args.weights.length_gt) ;
    printf("    seq_quality %f\n", p->o_args.weights.seq_quality) ;
    printf("    end_quality %f\n", p->o_args.weights.end_quality) ;
    printf("    pos_penalty %f\n", p->o_args.weights.pos_penalty) ;
    printf("    end_stability %f\n", p->o_args.weights.end_stability) ;
    printf("    num_ns %f\n", p->o_args.weights.num_ns) ;
    printf("  end internal oligo_weights\n") ;

    printf("  opt_tm %f\n", p->o_args.opt_tm) ;
    printf("  min_tm %f\n", p->o_args.min_tm) ;
    printf("  max_tm %f\n", p->o_args.max_tm) ;
    printf("  opt_gc_content %f\n", p->o_args.opt_gc_content) ;
    printf("  max_gc %f\n", p->o_args.max_gc) ;
    printf("  min_gc %f\n", p->o_args.min_gc) ;
    printf("  divalent_conc %f\n", p->o_args.divalent_conc) ;
    printf("  dntp_conc %f\n", p->o_args.dntp_conc) ;
    printf("  dna_conc %f\n", p->o_args.dna_conc) ;
    printf("  num_ns_accepted %i\n", p->o_args.num_ns_accepted) ;
    printf("  opt_size %i\n", p->o_args.opt_size) ;
    printf("  min_size %i\n", p->o_args.min_size) ;
    printf("  max_size %i\n", p->o_args.max_size) ;
    printf("  max_poly_x %i\n", p->o_args.max_poly_x) ;
    printf("  min_end_quality %i\n", p->o_args.min_end_quality) ;
    printf("  min_quality %i\n", p->o_args.min_quality) ;
    printf("  max_self_any %f\n", p->o_args.max_self_any) ;
    printf("  max_self_end %f\n", p->o_args.max_self_end) ;
    printf("  max_repeat_compl %f\n", p->o_args.max_repeat_compl) ;
    printf("  end internal oligo args\n");
    printf("\n");
    printf("END GLOBAL ARGS\n");
    printf("=============\n");
    printf("\n");
  }

  if (s != NULL) {
    printf("=============\n");
    printf("BEGIN SEQUENCE ARGS\n") ;
    /* TO DO: complete the statments for dumping this data
       printf("interval_array_t2 tar2 %i\n",
       int pairs[PR_MAX_INTERVAL_ARRAY][2]) ;
       int count) ;
       printf("interval_array_t2 excl2 %i\n",
       int pairs[PR_MAX_INTERVAL_ARRAY][2]) ;
       int count) ;
       printf("interval_array_t2 excl_internal2 %i\n",
       int pairs[PR_MAX_INTERVAL_ARRAY][2]) ;
       int count) ;
       printf("ok_regions \n");
    */

    if (s->primer_overlap_junctions_count > 0) {
      printf("primer_overlap_junctions_count %i\n",
             s->primer_overlap_junctions_count);
      printf("primer_overlap_junctions_list [\n");
      for (i = 0; i < s->primer_overlap_junctions_count; i++) {
        printf("   %i\n", s->primer_overlap_junctions[i]);
      }
      printf("]\n");
    }

    printf("incl_s %i\n", s->incl_s) ;
    printf("incl_l %i\n", s->incl_l) ;
    printf("start_codon_pos %i\n", s->start_codon_pos) ;
    printf("n_quality %i\n", s->n_quality) ;
    /* TO DO printf("quality%i\", s->quality) ; */
    printf("quality_storage_size %i\n", s->quality_storage_size) ;
    printf("*sequence %s\n", s->sequence) ;
    printf("*sequence_name %s\n", s->sequence_name) ;
    printf("*sequence_file %s\n", s->sequence_file) ;
    printf("*trimmed_seq %s\n", s->trimmed_seq) ;
    printf("*trimmed_orig_seq %s\n", s->trimmed_orig_seq) ;
    printf("*upcased_seq %s\n", s->upcased_seq) ;
    printf("*upcased_seq_r %s\n", s->upcased_seq_r) ;
    printf("*left_input %s\n", s->left_input) ;
    printf("*right_input %s\n", s->right_input) ;
    printf("*internal_input %s\n", s->internal_input) ;
    printf("force_left_start %i\n", s->force_left_start) ;
    printf("force_left_end %i\n", s->force_left_end) ;
    printf("force_right_start %i\n", s->force_right_start) ;
    printf("force_right_end %i\n", s->force_right_end) ;
    printf("END SEQUENCE ARGS\n") ;
    printf("=============\n");
    printf("\n");
  }
}
