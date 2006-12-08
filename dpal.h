#ifndef _DPAL_H
#define _DPAL_H
#include <limits.h>

#ifndef DPAL_MAX_ALIGN
#define DPAL_MAX_ALIGN   1600 /* 
			       * The maximum size of a string that can be
			       * aligned with with generic dpal and for which
			       * we can return a "path".  Several arrays of
			       * size DPAL_MAX_ALIGN X DPAL_MAX_ALIGN are
			       * statically allocated in dpal.o
			       */
#endif

#define DPAL_LOCAL        0  /* Return a local alignment. */
#define DPAL_GLOBAL_END   1  /* 
			      * Return a global alignment _anchored at the end
			      * of the first sequence_.
			      */
#define DPAL_GLOBAL       2  /* 
			      * Return an arbitrary global alignment, that is
			      * one anchored at the end of either the first or
			      * the second sequence.
			      */
#define DPAL_LOCAL_END    3   /* 
                               * Return a local alignment that includes the
			       * end (but not necessarily the beginning) of
			       * the first sequence.
			       */

/*
 * It is not possible to specify end-gap penalties for the DPAL_GLOBAL_END
 * and DPLAL_GLOBAL flags.
 */

/* 
 * The data structure that stores the "scoring system matrix". (The socring
 * system matrix data structure is of size UCHAR_MAX + 1 by UCHAR_MAX + 1.)
 */
typedef int dpal_ssm[UCHAR_MAX + 1][UCHAR_MAX + 1];

/* Structure for passing in arguments to the main function, dpal. */
typedef struct {
    int check_chars;        /* 
		             * If non-0, check for and raise an error on an
		             * illegal character in the input strings.
			     */
    int debug;              /* 
			     * If non-0, print debugging information to
			     * stderr.
			     */
    int fail_stop;           /* Exit with -1 on error. */
    int flag;                /* 
			      * One of DPAL_GLOBAL, DPAL_LOCAL,
			      * DPAL_GLOBAL_END, DPAL_LOCAL_END
			      */
    int force_generic;      /* Force the use of the generic function. */
    int force_long_generic; /* 
			     * Force the use of the long generic no-path
			     * function.
			     */
    int force_long_maxgap1; /* Force the use of the long maxgap 1 functions. */
    int gap;                 /* The "gap opening" penalty. */
    int gapl;                /* The "gap extension" penalty. */
    int max_gap;             /* 
		              * The maximum allowable size for a gap. -1
		              * indicates that the gap can be of any size.
			      */
    int score_max;           /* If greater than 0 stop search as soon as
			      * score > score_max.
			      */
    int score_only;          /* 
			      * If non-0, only print the score on
			      * stdout. (Incompatible with debug.)
			      */
    dpal_ssm ssm;            /* The scoring system matrix. */
} dpal_args;

/* Structure for receiving results from the main function, dpal. */
typedef struct {
    const char *msg;
    int   path[DPAL_MAX_ALIGN][2];
    int   path_length;
    int   align_end_1; /* Last alignment position in the 1st sequence. */
    int   align_end_2; /* Last alignment position in the 2nd sequence. */
    int   score;
} dpal_results;

/* Initialize the argument to the default matrix for nucleotide matches. */
void dpal_set_default_nt_args(dpal_args *);
/* Routine primarily for testing: sets CC & GG matches to 3, AA & TT 
   matches to 2. */
void dpal_set_h_nt_matrix(dpal_args *);

/* The argument a must be a DNA scoring matrix.  Modifies a so that it for a
   match between any two ambiguity codes (or between ambiguity code and base),
   e.g. B and S, the score will be the maximum of score between any base in B
   and any base in S, in the example between any pair in {C, G, T} X {C, G}.
   This function overwrites any scores already associated with pairs of
   ambiguity codes.  Return 0 on error, 1 on success.
*/
int dpal_set_ambiguity_code_matrix(dpal_args *);

/* 
 * Align the first 2 arguments, using the scoring matix mat and the additional
 * arguments supplied in the align_args.  Return the results in the
 * align_results struct.
 *
 * Returns 0 for success
 *         1 for failure
 */
int dpal(const unsigned char *, const unsigned char*,
	 const dpal_args *, dpal_results *);
#endif
