#ifndef PR_PRIMER_H
#define PR_PRIMER_H

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <setjmp.h>
#include "dpal.h"

#define PR_INFINITE_POSITION_PENALTY -1.0
#define PR_DEFAULT_OUTSIDE_PENALTY    0.0
#define PR_DEFAULT_INSIDE_PENALTY     PR_INFINITE_POSITION_PENALTY
#define PR_DEFAULT_PRODUCT_MAX_TM     1000000.0
#define PR_DEFAULT_PRODUCT_MIN_TM     -1000000.0

#define PR_POSITION_PENALTY_IS_NULL(PA) \
(PR_DEFAULT_INSIDE_PENALTY == (PA)->inside_penalty \
 && PR_DEFAULT_OUTSIDE_PENALTY == (PA)->outside_penalty)

#define PR_NULL_START_CODON_POS       -1000000
#define PR_DEFAULT_START_CODON_POS    PR_NULL_START_CODON_POS

#define PR_START_CODON_POS_IS_NULL(SA) \
((SA)->start_codon_pos <= PR_NULL_START_CODON_POS)

/* Maxima needed for interface data structures. */
#define PR_MAX_INTERVAL_ARRAY 200 /* 
				   * Maximum number of input intervals
				   * supported; used for targets, excluded
				   * regions, product-size intervals, etc.
				   */

#define PR_ALIGN_SCORE_PRECISION 100.0
#define PR_MAX_LIBRARY_WT 100.0

#define TRIMMED_SEQ_LEN(X) ((X)->incl_l)

#define MACRO_STRING(X) #X
#if defined(NDEBUG)
#define PR_ASSERT(COND)                                          \
if (!(COND)) {                                                   \
    fprintf(stderr, "libprimer3:%s:%d, assertion (%s) failed\n", \
	   __FILE__, __LINE__,                                   \
	   MACRO_STRING(COND));                                  \
}
#else
#define PR_ASSERT(COND)                                          \
if (!(COND)) {                                                   \
    fprintf(stderr, "libprimer3:%s:%d, assertion (%s) failed\n", \
	   __FILE__, __LINE__,                                   \
	   MACRO_STRING(COND));                                  \
    abort();                                                     \
}
#endif

#define DEFAULT_POSITION_PENALTIES(PA) \
    (PR_DEFAULT_INSIDE_PENALTY == pa->inside_penalty \
     && PR_DEFAULT_OUTSIDE_PENALTY == pa->outside_penalty)

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
                               OV_LIB_SIM=13
} oligo_violation;

typedef enum primer_errno {
    PR_ERR_NONE = 0,
    PR_ERR_OUT_OF_MEMORY = 1,
    PR_ERR_CANNOT_OPEN_FILE = 2,
    PR_ERR_ALIGNMENT_FAILED = 3
} primer_errno;

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
  rep_sim repeat_sim; /* Name of the sequence from given file in fasta
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
  double end_stability; /* Delta G of disription of 5 3' bases. */
  int    start;    /* The 0-based index of the leftmost base of the primer. */
  int    seq_quality; /* Minimum quality score of bases included. */   
  short  self_any; /* Self complementarity as local alignment * 100. */
  short  self_end; /* Self complementarity at 3' end * 100. */
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
 * Error handling information.
 * This contains error codes and string along with a jmp_buf for error
 * recovery.
 */
typedef struct primer_error {
    int system_errno;			/* A copy of the system 'errno' */
    primer_errno local_errno;		/* Primer3 error code */
    char *error_msg;			/* text version of local_errno */
    jmp_buf jmpenv;			/* errors caught in API funcs */
} primer_error;


/* 
 * The structure for a pair of primers. (So that we can have a function which
 * returns a pair of primers.)
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

/*
 * A global 'state' for primer3. All global variables are held within here.
 * Before using primer3 you need to create a new state using primer3_create.
 * Free the memory after usage with primer3_destroy.
 */
typedef struct primer_state {
    /* Alignment arguments */
    dpal_args local_args;		/* local, noends, noambig */
    dpal_args local_end_args;		/* local,   ends, noambig */
    dpal_args end_args;			/* global,  ends, noambig */
    dpal_args local_args_ambig;		/* local, noends,   ambig */
    dpal_args local_end_args_ambig;	/* local,   ends,   ambig */

    /* Primer info */
    primer_rec *f, *r, *mid;	/* The primers */
    int n_f, n_r, n_m;		/* Number of elements in f, r and mid */
    int f_len, r_len, mid_len;	/* and their lengths */
    pair_array_t best_pairs;	/* The best primer pairs */

    primer_error err;		/* Error handling */
} primer_state;

/*
 * This typedef poses a problem for const. "interval_array_t" cannot be
 * passed into a function requiring "const interval_array_t" as the const
 * binds to the wrong point. An alternative is to specifying the function
 * prototype as accepting "int (* const tar)[2]", but that's somewhat ugly.
 * This could be simulated with a typedef:
 *     typedef int (*const const_interval_array_t)[2];
 * That's also ugly, so for now, const has been removed from
 * interval_array_t arguments.
 */
typedef int interval_array_t[PR_MAX_INTERVAL_ARRAY][2];

/* pr_append_str is an append-only string ADT. */
typedef struct pr_append_str {
    int storage_size;
    char *data;
} pr_append_str;

/* The seq_lib struct represents a library of sequences. */
typedef struct seq_lib {
    char **names;         /* An array of sequence names. */
    char **seqs;          /* An array of sequences. */
    char **rev_compl_seqs;/* An array of reversed-complemented sequences. */
    double *weight;       /* An array of weights. */
    char   *repeat_file;  /* The path of the file containing the library. */
    pr_append_str error;  /* Global error message if any.  */
    pr_append_str warning;/* Warning message. */
    int seq_num;          /* The number of names, sequences, and weights. */
} seq_lib;

/* 
 * Arguments to the primer program as a whole.  Values for these arguments are
 * retained _across_ different input records.  (These are the so-called
 * "Global" arguments in the documentation.)
 */
typedef struct oligo_weigts {
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
} pair_weights;

typedef enum task { pick_pcr_primers = 0,
		    pick_pcr_primers_and_hyb_probe = 1,
		    pick_left_only = 2,
                    pick_right_only = 3,
                    pick_hyb_probe_only =4,
                } task;

typedef struct primargs {
    int    pr_min[PR_MAX_INTERVAL_ARRAY]; /* Minimum product sizes. */
    int    pr_max[PR_MAX_INTERVAL_ARRAY]; /* Maximum product sizes. */
    seq_lib repeat_lib;  /* Library of sequences to avoid. */

    seq_lib io_mishyb_library;

    oligo_weights primer_weights;
    oligo_weights io_weights;
    pair_weights  pr_pair_weights;

    pr_append_str glob_err;

    double opt_tm;
    double min_tm;
    double max_tm;
    double max_diff_tm;
    double opt_gc_content;
    double max_gc;
    double min_gc;
    double salt_conc;
    double dna_conc;

    double io_opt_tm;
    double io_min_tm;
    double io_max_tm;
    double io_opt_gc_content;
    double io_max_gc;
    double io_min_gc;
    double io_salt_conc;
    double io_dna_conc;

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

    double product_max_tm;
    double product_min_tm;
    double product_opt_tm;
    double max_end_stability;
                            /* The maximum value allowed for the delta
                             * G of disruption for the 5 3' bases of
                             * a primer.
                             */
    int    num_intervals;         /* 
				   * Number of product size intervals
				   * (i.e. number of elements in pr_min and
				   * pr_max)
                                   */
    int    num_ns_accepted;
    task   primer_task;          /* 2 if left primer only, 3 if right primer only,
				   * 4 if internal oligo only.    */

    int    file_flag;
    int    explain_flag;
    int    primer_opt_size;
    int    primer_min_size;
    int    primer_max_size;
    int    product_opt_size;

  /*internal oligo*/
    int    io_num_ns_accepted;
    int    io_primer_opt_size;
    int    io_primer_min_size;
    int    io_primer_max_size;

    int    gc_clamp;              /* Required number of GCs at *3' end. */
   

    int    liberal_base;   /* 
			    * If non-0 then turn characters other than
			    * [ATGCNatgcn] into N.
			    */

    int    max_poly_x;      /* 
			     * Maximum length of mononucleotide sequence in an
			     * oligo.
			     */
    int    io_max_poly_x;


    int    first_base_index;  /* 
			       * The index of the first base in the input
			       * sequence.  This parameter is ignored within
			       * primer3_choose; caller must assure that
			       * all indexes are 0-based.  However, this
			       * parameter should used by output routines to
			       * adjust base indexes.
			       */
    int    num_return; /* The number of best primer pairs to return. */
    int    min_quality;       /* Minimum quality permitted for oligo sequence.*/
    int    min_end_quality;   /* Minimum quality permitted at 3' end. */
    int    quality_range_min;
    int    quality_range_max;

    int    io_min_quality;
    int    io_min_end_quality;

    int    pick_anyway;    /* Pick even if input primer or oligos
                              violate constraints. */
                                

    short  repeat_compl;   /* 
			    * Acceptable complementarity with repeat
			    * sequences.
			    */
    short  io_repeat_compl;

    short  pair_repeat_compl;

    short  self_any;  
    short  self_end;

    short  io_self_any;  
    short  io_self_end;

    short  pair_compl_any;
    short  pair_compl_end;
} primer_args;

typedef struct pair_stats {
    int considered;         /* Total number of pairs or triples tested.      */
    int product;            /* Pairs providing incorrect product size.       */
    int target;             /* Pairs without any target between primers.     */
    int temp_diff;          /* Melting temperature difference too high.      */
    int compl_any;          /* Pairwise complementarity larger than allowed. */
    int compl_end;          /* The same for 3' end complementarity.          */
    int internal;           /* Internal oligo was not found.                 */
    int repeat_sim;         /* Complementarity with repeat sequence too high.*/
    int high_tm;            /* Product Tm too high.                          */
    int low_tm;             /* Product Tm too low.                           */
    int ok;                 /* Number that were ok.                          */
} pair_stats;

typedef struct oligo_stats {
    int considered;         /* Total number of tested oligos of given type   */
    int ns;                 /* Number of oligos rejected because of Ns       */
    int target;             /* Overlapping targets.                          */
    int excluded;           /* Overlapping excluded regions.                 */
    int gc;                 /* Unappropriate GC content.                     */
    int gc_clamp;           /* Don't have required number of GCs at 3' end.  */
    int temp_min;           /* Melting temperature below t_min.              */
    int temp_max;           /* Melting temperature more than t_max.          */
    int compl_any;          /* Self-complementarity too high.                */
    int compl_end;          /* Self-complementarity at 3' end too high.      */
    int repeat;             /* Complementarity with repeat sequence too hign.*/
    int poly_x;             /* Long mononucleotide sequence inside.          */
    int seq_quality;        /* Low quality of bases included.                */
    int stability;          /* Stability of 5 3' bases too high.             */
    int no_orf;             /* Would not amplify any of the specified orf
                               (valid for left primers only).                 */
    int ok;                 /* Number of acceptable oligos in this category. */
} oligo_stats;

/*
 * Arguments relating to a particular source sequence for which
 * we will pick primers.  Values for seq_args are _not_ retained
 * across different input records.
 */
typedef struct seqargs {
    pr_append_str error;    /* Error messages. */
    pr_append_str warning;  /* Warning messages. */
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
    char *sequence;         /* The sequence itself. */
    char *sequence_name;    /* An identifier for the sequence. */
    char *sequence_file;    /* Another identifer for the sequence. */

    char *trimmed_seq;      /* The included region only. */
    char *left_input;       /* A left primer to check or design around. */
    char *right_input;      /* A right primer to check or design around. */
    char *internal_input;   /* An internal oligo to check or design around. */

    oligo_stats left_expl;  /* Left primers statistics. */
    oligo_stats right_expl; /* Right primers statistics. */
    oligo_stats intl_expl;  /* Internal oligos statistics. */
    pair_stats  pair_expl;  /* Pair statistics. */
} seq_args;

/*
 * Return the ('\0'-terminated) sequence of oligo o in static storage that is
 * overwritten on each call.
 */
char *pr_oligo_sequence(const seq_args *, const primer_rec *o);

/*
 * Return the ('\0'-terminated) reverse-complemented sequence of oligo o in
 * static storage that is overwritten on each call.
 */
char *pr_oligo_rev_c_sequence(const seq_args *, const primer_rec *o);

/* Append s to x. Returns 0 for success, -1 for failure */
int pr_append(pr_append_str *, const char *);

/* Append s to x if x is empty.  Otherwise append concat("; ",s) to x. */
int pr_append_new_chunk(pr_append_str *x, const char *s);

/* Append s to x if x is empty. Otherwise append concat(sep, s) to x. */
int pr_append_w_sep(pr_append_str *x, const char *sep, const char *);

/* Set the argument to the empty string. */
void pr_set_empty(pr_append_str *);

/* Return true iff the argument is empty. */
int pr_is_empty(const pr_append_str *);

/* 
 * Gather all warnings and return a char * containing them.  The return value
 * points to storage that the caller must free.  Return NULL if there are no
 * errors.
 */
char *pr_gather_warnings(const seq_args *, const primer_args *);

/* Choose primer pairs or triples */
int prime3_choose(primer_state *, primer_args *, seq_args *, pair_array_t *);

/* Print an "explanation string" for primer pair search. */
void pr_print_pair_explain(FILE *, const seq_args *);

void reverse_complement(const char *, char *);

int    strcmp_nocase(char *, char *);

void   set_default_global_primer_args(primer_args *);

/* Allocate and deallocate a new primer3 state */
primer_state *primer3_create(void);
void primer3_destroy(primer_state *);

/* 
 * Finds up to pa->num_return primer pairs for the sequence seq with t targets.
 * Set sa->error and return 1 on error; otherwise return 0.
 */
int
primer3_choose(primer_state *state,
	       primer_args *pa,
	       seq_args *sa);

/*
 * -----------------------------------------------------------------
 * longjmp utilising versions of the pr_append* functions and malloc
 * functions. These are only used from within the primer3_choose (and
 * children) and read_record functions as primer3_choose and read_record
 * set up the jmpbuf handler required for use in these jump_* functions.
 */
FILE   *jump_fopen(primer_error *, const char*, const char*);
void *pr_jump_malloc(primer_error *, size_t);
void *pr_jump_realloc(primer_error *, void *, size_t);
void jump_error(primer_error *, primer_errno);
void jump_append(primer_error *, pr_append_str *, const char *);
void jump_append_new_chunk(primer_error *, pr_append_str *x, const char *s);
void jump_append_w_sep(primer_error *, pr_append_str *x, const char *sep, const char *);

#endif
