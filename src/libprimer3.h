/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007
Whitehead Institute for Biomedical Research, Steve Rozen
(http://jura.wi.mit.edu/rozen), and Helen Skaletsky
All rights reserved.

    This file is part of primer3 and the libprimer3 library.

    Primer3 and the libprimer3 library are free software;
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

#ifndef LIBPRIMER3_H
#define LIBPRIMER3_H

#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h> /* SHRT_MIN */

/* ALIGN_SCORE_UNDEF is used only libprimer3 and clients, not in dpal */
#define ALIGN_SCORE_UNDEF             SHRT_MIN

/* These next 5 are exposed for format_output.c -- probabaly should be reviewed. */
#define PR_INFINITE_POSITION_PENALTY -1.0
#define PR_DEFAULT_INSIDE_PENALTY     PR_INFINITE_POSITION_PENALTY
#define PR_DEFAULT_OUTSIDE_PENALTY    0.0
#define PR_DEFAULT_PRODUCT_MAX_TM     1000000.0
#define PR_DEFAULT_PRODUCT_MIN_TM     -1000000.0

/*  Exposed in the read_boulder input routine.... */
#define PR_NULL_START_CODON_POS       -1000000
#define PR_DEFAULT_START_CODON_POS    PR_NULL_START_CODON_POS
#define PR_START_CODON_POS_IS_NULL(SA) ((SA)->start_codon_pos <= PR_NULL_START_CODON_POS)

#define _PR_DEFAULT_POSITION_PENALTIES(PA) \
    (PR_DEFAULT_INSIDE_PENALTY == pa->inside_penalty \
     && PR_DEFAULT_OUTSIDE_PENALTY == pa->outside_penalty)

#define PR_ALIGN_SCORE_PRECISION 100.0

#define MACRO_STRING(X) #X
/* pr_progam_name must be set in main(). */
#define PR_ASSERT(COND)                                  \
if (!(COND)) {                                           \
    fprintf(stderr, "%s:%s:%d, assertion (%s) failed\n", \
	   pr_program_name, __FILE__, __LINE__,          \
	   MACRO_STRING(COND));                          \
    abort();                                             \
}

typedef enum task { pick_pcr_primers               = 0,
		    pick_pcr_primers_and_hyb_probe = 1,
		    pick_left_only                 = 2,
                    pick_right_only                = 3,
                    pick_hyb_probe_only            = 4,
                } task;

/* pr_append_str is an append-only string ADT. */
typedef struct pr_append_str {
    int storage_size;
    char *data;
} pr_append_str;

/* The seq_lib struct represents a library of sequences. */
/* Clients do not need to know the details of this structure. */
typedef struct seq_lib {
    char **names;         /* An array of sequence names. */
    char **seqs;          /* An array of sequences. */
    char **rev_compl_seqs;/* An array of reversed-complemented sequences.
                             x->rev_compl_seqs[i] is the reverse complement
                             of x->seqs[i], which lets us keep track of pairwise
                             mispriming.  See reverse_complement_seq_lib(). */
    double *weight;       /* An array of weights. */
    char   *repeat_file;  /* The path of the file containing the library. */
    pr_append_str error;  /* Global error message if any.  */
    pr_append_str warning;/* Warning message. */
    int seq_num;          /* The number of names, sequences, and weights. */
} seq_lib;

/* Maxima needed for interface data structures. */
#define PR_MAX_INTERVAL_ARRAY 200 /* 
				   * Maximum number of input intervals
				   * supported; used for targets, excluded
				   * regions, product-size intervals, etc.
				   */

/* 
 * Arguments to the primer program as a whole.  Values for these arguments are
 * retained _across_ different input records.  (These are the so-called
 * "Global" arguments in the documentation.)
 */
typedef struct oligo_weights {
    double temp_gt;
    double temp_lt;
    double gc_content_gt;
    double gc_content_lt;
    double compl_any;
    double compl_end;
    double repeat_sim;
    double length_lt;
    double length_gt;
    double seq_quality;
    double end_quality;
    double pos_penalty;
    double end_stability;
    double num_ns;
    double template_mispriming;
} oligo_weights;

typedef struct pair_weights {
    double primer_quality;
    double io_quality;
    double diff_tm;
    double compl_any;
    double compl_end;
    double product_tm_lt;
    double product_tm_gt;
    double product_size_lt;
    double product_size_gt;
    double repeat_sim;
    double template_mispriming;
} pair_weights;

typedef struct args_for_one_oligo_or_primer {
  seq_lib       *repeat_lib;
  oligo_weights weights;
  double opt_tm;
  double min_tm;
  double max_tm;
  double opt_gc_content;
  double max_gc;
  double min_gc;

  /* FIX ME -- skewed -- also used for product Tm */
  double salt_conc;
  double divalent_conc; /* added by T.Koressaar, divalent salt concentration mmol/l */
  double dntp_conc; /* added by T.Koressaar, for considering divalent salt concentration */

  double dna_conc;
  int    num_ns_accepted;
  int    opt_size;
  int    min_size;
  int    max_size;
  int    max_poly_x;      /* 
			   * Maximum length of mononucleotide sequence in an
			   * oligo.
			   */

  int    min_end_quality;
  int    min_quality;       /* Minimum quality permitted for oligo sequence.*/

  short  max_self_any;  
  short  max_self_end;
  short  max_repeat_compl;   /* 
			      * Acceptable complementarity with repeat
			      * sequences.
			      */

  short  max_template_mispriming;
} args_for_one_oligo_or_primer;


typedef struct primargs {
  /* ================================================== */
  /* Arguments that control behavior of choose_primers() */
  task   primer_task;          /* 2 if left primer only, 3 if right primer only,
				* 4 if internal oligo only.    */

  int    file_flag;
  int    explain_flag;

  int    first_base_index;  /* 
			     * The index of the first base in the input
			     * sequence.  This parameter is ignored within
			     * pr_choice; pr_choice's caller must assure that
			     * all indexes are 0-based.  However, this
			     * parameter should used by output routines to
			     * adjust base indexes.
			     */

  int    liberal_base;   /* 
			  * If non-0 then turn characters other than
			  * [ATGCNatgcn] into N.
			  */

  int    num_return; /* The number of best primer pairs to return. */

  int    pick_anyway;    /* Pick even if input primer or oligos
			    violate constraints. */

  int    lib_ambiguity_codes_consensus;
  /* If non-0, treat ambiguity codes in a mispriming/mishyb
     library as representing a consensus.  So, for example,
     S would match C or G.  N would match any nucleotide.
     It turns out that this _not_ what one normally wants,
     since many libraries contain strings of N, which then
     match every oligo (very bad).
  */

  int    quality_range_min;
  int    quality_range_max;

  /* ================================================== */
  /* Writable return argument for errors. */
  pr_append_str glob_err;
  /* FIX ME (maybe) make this a return argument. */

  /* ================================================== */
  /* Arguments for individual oligos and/or primers */
  args_for_one_oligo_or_primer p_args;
  args_for_one_oligo_or_primer o_args;

  /* Added by T.Koressaar. Specify the table of thermodynamic
     parameters to use (Options are ... SantaLucia 1998) */
  int tm_santalucia;  

  /* Added by T.Koressaar. Specify the salt correction formula for Tm
     calculations. (Options are ... */
  int salt_corrections; 

  /* Arguments applicable to primers but not oligos */
  double max_end_stability;
  /* The maximum value allowed for the delta
   * G of disruption for the 5 3' bases of
   * a primer.
   */
  int    gc_clamp;              /* Required number of GCs at *3' end. */

  /* ================================================== */
  /* Arguments related to primer and/or oligo
     location in the template. */

  int lowercase_masking; 
  /* added by T.Koressaar for primer design from lowercase masked
     template */

  double outside_penalty; /* Multiply this value times the number of NTs
			   * from the 3' end to the the (unique) target to
			   * get the 'position penalty'.
			   * Meaningless if there are multiple targets
			   * or if the primer cannot be part of a pair
			   * that spans the target.
			   */

  double inside_penalty;  /* Multiply this value times the number of NT
			   * positions by which the primer overlaps
			   * the (unique) target to the 'position penalty'.
			   * Meaningless if there are multiple targets
			   * or if the primer cannot be part of a pair
			   * that spans the target.
			   */


  /* ================================================== */
  /* Arguments for primer pairs and products. */

  int    pr_min[PR_MAX_INTERVAL_ARRAY]; /* Minimum product sizes. */
  int    pr_max[PR_MAX_INTERVAL_ARRAY]; /* Maximum product sizes. */
  int    num_intervals;         /* 
				 * Number of product size intervals
				 * (i.e. number of elements in pr_min and
				 * pr_max)
				 */

  int    product_opt_size;
  double product_max_tm;
  double product_min_tm;
  double product_opt_tm;
  short  pair_max_template_mispriming;
  short  pair_repeat_compl;
  short  pair_compl_any;
  short  pair_compl_end;

  /* Max diff between tm of primer and tm of product.
     Cannot be calculated until product is known. */
  double max_diff_tm; 

  pair_weights  pr_pair_weights;

} primer_args;

typedef enum oligo_type { OT_LEFT = 0, OT_RIGHT = 1, OT_INTL = 2 }
  oligo_type;

typedef enum oligo_violation { OV_UNINITIALIZED = -1,
			       OV_OK=0, 
			       OV_TOO_MANY_NS=1, 
			       OV_INTERSECT_TARGET=2,
                               OV_GC_CONTENT=3, 
			       OV_TM_LOW=4, 
			       OV_TM_HIGH=5, 
			       OV_SELF_ANY=6,
                               OV_SELF_END=7,
			       OV_EXCL_REGION=8,
                               OV_GC_CLAMP=9,
			       OV_END_STAB=10, 
			       OV_POLY_X=11,
			       OV_SEQ_QUALITY=12,
                               OV_LIB_SIM=13,
			       OV_TEMPLATE_MISPRIMING=14,
                               OV_GMASKED=14 /* edited by T. Koressaar for lowercase masking */ 
} oligo_violation;

typedef struct rep_sim {
  char *name;      /* Name of the sequence from given file in fasta
		    * format with maximum similarity to the oligo.
		    */
  short min;       /* 
                    * The minimum score in slot 'score' (below).
                    * (Used when the objective function involves
                    * minimization of mispriming possibilities.)
                    */
  short max;       /* The maximum score in slot 'score' (below). */
  short *score;    /* 
                    * Array of similarity (i.e. false-priming) scores,
                    * one for each entry in the 'repeat_lib' slot
                    * of the primargs struct. 
                    */
} rep_sim;

typedef struct primrec {

  rep_sim repeat_sim;
                   /* Name of the sequence from given file in fasta
		    * format with maximum similarity to the oligo
		    * and corresponding alignment score.
		    */

  double temp;     /* 
		    * The oligo melting temperature calculated for the
		    * primer.
		    */

  double gc_content;

  double position_penalty; 
                  /*
                   * Penalty for distance from "ideal" position as specified
	           * by inside_penalty and outside_penalty.
                   */

  double quality;  /* Part of objective function due to this primer. */

  double end_stability;
                   /* Delta G of disription of 5 3' bases. */
  int    start;    /* The 0-based index of the leftmost base of the primer
                      WITH RESPECT TO THE seq_args FIELD trimmed_seq. */
  int    seq_quality; /* Minimum quality score of bases included. */   
  short  self_any; /* Self complementarity as local alignment * 100. */
  short  self_end; /* Self complementarity at 3' end * 100. */
  short  template_mispriming;
                   /* Max 3' complementarity to any ectopic site in template
		      on the given template strand. */
  short  template_mispriming_r;
                   /* Max 3' complementarity to any ectopic site in the
		      template on the reverse complement of the given template
		      strand. */
  char   target;   /* 
		    * 0 if this primer does not overlap any target, 1 if it
		    * does.
		    */
  char   excl;     /* 
		    * 0 if does not overlap any excluded region, 1 if it
		    * does.
		    */
  oligo_violation ok;
  char   length;   /* Length of the oligo. */
  char   num_ns;   /* Number of Ns in the oligo. */
  char   position_penalty_infinite; 
                   /* Non-0 if the position penalty is infinite. */
  char   must_use; /* Non-0 if the oligo must be used even if it is illegal. */
} primer_rec;

/* 
 * The structure for a pair of primers. (So that we can have a function 
 * that returns a pair of primers.)
 */
typedef struct primpair {
  double pair_quality;
  double compl_measure; /* 
			 * A measure of self-complementarity of left and right
			 * primers in the pair, as well as complementarity
			 * between left and right primers.  The function
			 * choice returns pairs with the minimal value for
			 * this field when 2 pairs have the same
			 * pair_quality.
			 */
  double diff_tm;       /* Absolute value of the difference between melting
			 * temperatures for left and right primers. 
			 */
   
  double product_tm;    /* Estimated melting temperature of the product. */

  double product_tm_oligo_tm_diff;
                        /* Difference in Tm between the primer with lowest Tm
			   the product Tm. */

  double t_opt_a;

  int    compl_any;     /* 
			 * Local complementarity score between left and right
			 * primers (* 100).
			 */

  int    compl_end;     /* 
		         * 3'-anchored global complementatory score between *
		         * left and right primers (* 100).
			 */

  int    template_mispriming;
                        /* Maximum total mispriming score of both primers
			   to ectopic sites in the template, on "same"
			   strand (* 100). */

  short  repeat_sim;    /* Maximum total similarity of both primers to the
			 * sequence from given file in fasta format.
			 */
  primer_rec *left;     /* Left primer. */
  primer_rec *right;    /* Right primer. */
  primer_rec *intl;     /* Internal oligo. */

  int    product_size;  /* Product size. */
  int    target;        /* 
			 * 1 if there is a target between the right and left
			 * primers.
			 */
  char   *rep_name;
} primer_pair;

typedef struct pair_array_t {
    int         storage_size;
    int         num_pairs;
    primer_pair *pairs;
} pair_array_t;

typedef int interval_array_t[PR_MAX_INTERVAL_ARRAY][2];

typedef struct oligo_stats {
  int considered;          /* Total number of tested oligos of given type   */
  int ns;                  /* Number of oligos rejected because of Ns       */
  int target;              /* Overlapping targets.                          */
  int excluded;            /* Overlapping excluded regions.                 */
  int gc;                  /* Unacceptable GC content.                      */
  int gc_clamp;            /* Don't have required number of GCs at 3' end.  */
  int temp_min;            /* Melting temperature below t_min.              */
  int temp_max;            /* Melting temperature more than t_max.          */
  int compl_any;           /* Self-complementarity too high.                */
  int compl_end;           /* Self-complementarity at 3' end too high.      */
  int repeat_score;        /* Complementarity with repeat sequence too high.*/
  int poly_x;              /* Long mononucleotide sequence inside.          */
  int seq_quality;         /* Low quality of bases included.                */
  int stability;           /* Stability of 5 3' bases too high.             */
  int no_orf;              /* Would not amplify any of the specified ORF
			     (valid for left primers only).                 */
  int template_mispriming; /* Template mispriming score too high.           */
  int ok;                  /* Number of acceptable oligos.                  */
   int gmasked;            /* edited by T. Koressaar, number of gmasked oligo*/
} oligo_stats;

typedef struct pair_stats {
  int considered;          /* Total number of pairs or triples tested.      */
  int product;             /* Pairs providing incorrect product size.       */
  int target;              /* Pairs without any target between primers.     */
  int temp_diff;           /* Melting temperature difference too high.      */
  int compl_any;           /* Pairwise complementarity larger than allowed. */
  int compl_end;           /* The same for 3' end complementarity.          */
  int internal;            /* Internal oligo was not found.                 */
  int repeat_sim;          /* Complementarity with repeat sequence too high.*/
  int high_tm;             /* Product Tm too high.                          */
  int low_tm;              /* Product Tm too low.                           */
  int template_mispriming; /* Sum of template mispriming scores too hihg.   */
  int ok;                  /* Number that were ok.                          */
} pair_stats;

/*
 * Arguments relating to a single particular source sequence (for which
 * we will pick primer(s), etc.
 */
typedef struct seq_args {
  int num_targets;        /* The number of targets. */
  interval_array_t tar;   /*
			   * The targets themselves; tar[i][0] is the start
			   * of the ith target, tar[i][1] its length.  These
			   * are presented as indexes within the sequence
			   * slot, but during the execution of choice() they
			   * are recalculated to be indexes within
			   * trimmed_seq.
			   */
  int num_excl;           /* The number of excluded regions.  */
  interval_array_t excl;  /* The same as for targets.
			   * These are presented as indexes within
			   * the sequence slot, but during the
			   * execution of choice() they are recalculated
			   * to be indexes within trimmed_seq.
			   */
  int num_internal_excl;  /* Number of excluded regions for internal oligo.*/
  interval_array_t excl_internal;
  /* Similar to excl. */
  int incl_s;             /* The 0-based start of included region. */
  int incl_l;             /* 
			   * The length of the included region, which is
			   * also the length of the trimmed_seq field.
			   */
  int  start_codon_pos;   /* Index of first base of the start codon. */
  int  stop_codon_pos;    /* 
			   * An optional _output_, meaninful if a
			   * start_codon_pos is "not nul".  The position of
			   * the intial base of the leftmost stop codon that
			   * is to the right of sa->start_codon_pos.
			   */
  int  *quality;          /* Vector of quality scores. */
  int n_quality;          /* Length 'quality' */
  char *sequence;         /* The template sequence itself as input, 
			     not trimmed, not up-cased. */
  char *sequence_name;    /* An identifier for the sequence. */
  char *sequence_file;    /* Another identifer for the sequence. */
  char *trimmed_seq;      /* The included region only, _UPCASED_. */

  /* Element add by T. Koressaar support lowercase masking: */
  char *trimmed_orig_seq; /* Trimmed version of the original,
			     mixed-case sequence. */

  char *upcased_seq;      /* Upper case version of sequence
			     (_not_ trimmed). */
  char *upcased_seq_r;    /* Upper case version of sequence, 
			     other strand (_not_ trimmed). */
  char *left_input;       /* A left primer to check or design around. */
  char *right_input;      /* A right primer to check or design around. */
  char *internal_input;   /* An internal oligo to check or design around. */

  /* ================================================== */
  /*  Output (writable) arguments. */  
  pr_append_str error;    /* Error messages. */
  pr_append_str warning;  /* Warning messages. */
  oligo_stats left_expl;  /* Left primers statistics. */
  oligo_stats right_expl; /* Right primers statistics. */
  oligo_stats intl_expl;  /* Internal oligos statistics. */
  pair_stats  pair_expl;  /* Pair statistics. */
  /*  END output arguments. */
  /* ================================================== */

} seq_args;

/*
 * A global 'state' for primer3. All formerly global variables are held here.
 * Free the memory after usage with destroy_p3retval
 */
typedef struct p3retval {
  /* Arrays of oligo (primer) records. */
  primer_rec *f, *r, *mid;

  /* Number of initialized elements in f, r and mid */
  int n_f, n_r, n_m;		

  /* Storage lengths of f, r, mid */
  int f_len, r_len, mid_len;

  /* Array of best primer pairs */
  pair_array_t best_pairs;

  /* Temporary */
  int other_error;
} p3retval;

/* Allocate and initialize a new primer3 state.  Return NULL on
   ENOMEM. */
p3retval *create_p3retval(void);

/* Deallocate a primer3 state */
void destroy_p3retval(p3retval *);

void destroy_seq_args(seq_args *);

/* 
 * Choose individual primers or oligos, or primer pairs, or primer
 * pairs with internal oligos. On ENOMEM return -1 and set errno. On
 * other error return 1 and set sa->error and/or pa->glob_err. On no
 * error return 0.
 */
p3retval *choose_primers(p3retval *retval,  primer_args *pa, seq_args *sa);

char  *pr_oligo_sequence(const seq_args *, const primer_rec *);

char  *pr_oligo_rev_c_sequence(const seq_args *, const primer_rec *);

void  pr_set_default_global_args(primer_args *);

char  *pr_gather_warnings(const seq_args *, const primer_args *);

void  pr_append(pr_append_str *, const char *);
void  pr_append_new_chunk(pr_append_str *, const char *);

void  pr_print_pair_explain(FILE *, const seq_args *);

const char  *libprimer3_release(void);
const char  **libprimer3_copyright(void);

/* An accessor function for a primer_rec *. */
short oligo_max_template_mispriming(const primer_rec *);

int   strcmp_nocase(char *, char *);

/* ======================================================= */
/* Functions for creating and destroying a seq_lib object. */
/* ======================================================= */

/*  
 * Reads any file in fasta format and returns a newly allocated
 * seq_lib, lib.  Sets lib.error to a non-empty string on any error
 * other than ENOMEM.  Returns NULL on ENOMEM.
 */
seq_lib *
read_and_create_seq_lib(const char *filename, const char* errfrag);

void
destroy_seq_lib(seq_lib *lib);

/* number of sequences in a seq_lib* */
int
seq_lib_num_seq(const seq_lib* lib);

char *
seq_lib_warning_data(const seq_lib *lib);

#endif
