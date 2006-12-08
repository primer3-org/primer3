/*
 * Copyright (c) 1996, Whitehead Institute for Biomedical Research. All rights
 * reserved.  Please see full software use agreement in primer3_main.c or by
 * executing primer3 with -h.
 */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "dpal.h"
#include "primer3_release.h"

/*
 * We should probably remove the DPAL_FORGET_PATH compile-time option.
 * Efficiency now derives primarily for specialized versions of _dpal* for
 * particular parameter values.
 */
#ifndef DPAL_FORGET_PATH
/* 
 * Print an alignment on stderr, given the 2 aligned sequences, the "trace"
 * matrix, and the coordinate of the end point of the alignment to print.
 */
static void print_align(const unsigned char *, const unsigned char *,
			int[DPAL_MAX_ALIGN][DPAL_MAX_ALIGN][3], int, int,
			const dpal_args*);
#endif

/* 
 * Return 1 if there is an illegal character in the first argument, and place
 * the illegal character in the address contained in the last argument.
 */
static int  illegal_char(const unsigned char *, const dpal_ssm, char *);

static const char *xlate_ambiguity_code(char);

/*
 * The next function headers are for various versions of the
 * dynamic-programming alignment code optimized for particular input argument
 * values. 
 */
static int _dpal_generic(const unsigned char *,
			 const unsigned char *,
			 const int,
			 const int,
			 const dpal_args *,
			 dpal_results *);

static void _dpal_long_nopath_maxgap1_local_end(const unsigned char *,
						const unsigned char *,
						const int,
						const int,
						const dpal_args *,
						dpal_results *);

static void _dpal_long_nopath_generic(const unsigned char *,
			  const unsigned char *,
                          const int,
			  const int,
                          const dpal_args *,
			  dpal_results *);

static void _dpal_long_nopath_maxgap1_local(const unsigned char *,
					    const unsigned char *,
					    const int,
					    const int,
					    const dpal_args *,
					    dpal_results *);

static void _dpal_long_nopath_maxgap1_global_end(const unsigned char *,
						 const unsigned char *,
						 const int,
						 const int,
						 const dpal_args *,
						 dpal_results *);

void
dpal_set_default_nt_args(a)
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

    a->check_chars        = 1;
    a->debug              = 0;
    a->fail_stop          = 1;
    a->flag               = DPAL_LOCAL;
    a->force_generic      = 0;
    a->force_long_generic = 0;
    a->force_long_maxgap1 = 0;
    a->gap                = -100;
    a->gapl               = -100;
    a->max_gap            = 3;
    a->score_only         = 0;
}

void
dpal_set_h_nt_matrix(a)
    dpal_args *a;
{
    unsigned int i, j;

    for (i = 0; i <= UCHAR_MAX; i++)
	for (j = 0; j <= UCHAR_MAX; j++)
	    if (('A' == i || 'C' == i || 'G' == i || 'T' == i || 'N' == i)
		&& ('A' == j || 'C' == j || 'G' == j || 'T' == j 
		    || 'N' == j)) {
		    if (i == 'N' || j == 'N') 
			a->ssm[i][j] = -50;
		    else if (i == j) {
		      if ('C' == i || 'G' == i)
			a->ssm[i][j] = 300;
		      else
			a->ssm[i][j] = 200;
		    }
		    else 
			a->ssm[i][j] = -50;
		} else
		    a->ssm[i][j] = INT_MIN;
}

/* The argument a must be a DNA scoring matrix.
   Modify a so that it for a match between
   any two ambiguity codes (or between ambiguity code and base),
   e.g. B and S, the score will be the maximum of
   score between any base in B and any base in S,
   in the example between any pair in {C, G, T} X {C, G}.
*/
int
dpal_set_ambiguity_code_matrix(a)
  dpal_args *a;
{
  const char *c1, *c2;
  const char *amb_codes = "BDHVRYKMSWN";
  const char *all_bases = "ACGT";
  const char *bases1, *bases2, *b1, *b2;
  int extreme;
  for (c1 = amb_codes; *c1; c1++) {
    bases1 = xlate_ambiguity_code(*c1);
    if (!bases1) return 0;

    /* Do matches between c1 and all other
       ambiguity codes. */
    for (c2 = amb_codes; *c2; c2++) {
      bases2 = xlate_ambiguity_code(*c2);
      if (!bases2) return 0;
      extreme = INT_MIN;
      for (b1 = bases1; *b1; b1++) {
	for (b2 = bases2; *b2; b2++) {
	  if (a->ssm[*b1][*b2] > extreme) {
	    extreme = a->ssm[*b1][*b2];
	  }
	}
      }
      /* extreme is now the maximum score
	 for a match between any 2 bases
         represented *c1, *c2. */
      a->ssm[*c1][*c2] = extreme;
    }

    /* Do matches between c1 and all bases. */
    for (b2 = all_bases; *b2; b2++) {
      extreme = INT_MIN;
      for (b1 = bases1; *b1; b1++) {
	if (a->ssm[*b1][*b2] > extreme) {
	  extreme = a->ssm[*b1][*b2];
	}
      }
      a->ssm[*c1][*b2] = extreme;
      a->ssm[*b2][*c1] = extreme;
    }
  }
  return 1;
}

static const char *
xlate_ambiguity_code(char c)
{
  if ('N' == c)      return "ACGT";
  else if ('B' == c) return "CGT";
  else if ('D' == c) return "AGT";
  else if ('H' == c) return "ACT";
  else if ('V' == c) return "ACG";
  else if ('R' == c) return "AG";
  else if ('Y' == c) return "CT";
  else if ('K' == c) return "GT";
  else if ('M' == c) return "AC";
  else if ('S' == c) return "CG";
  else if ('W' == c) return "AT";
  else return NULL; /* Error condition */
}

#define CHECK_ERROR(COND,MSG) if (COND) { out->msg = MSG; goto FAIL; }

int
dpal(X, Y, in, out)
    const unsigned char *X,*Y;
    const dpal_args *in;
    dpal_results *out;
{
    int xlen, ylen;
    char msg[] = "Illegal character in input: ?";

    /*
     * strlen's prototype is strlen(const char *s) which means compiler
     * warnings when passing "const unsigned char *". This means we need lots
     * of casts throughout this code, but we hide this ugliness behind a
     * macro.
     */
#define ustrlen(X) (strlen((const char *)(X)))

    out->score = INT_MIN;
    out->path_length = 0;
    out->msg = NULL;

    CHECK_ERROR(NULL == X,   "NULL first sequence");
    CHECK_ERROR(NULL == Y,   "NULL second sequence");
    CHECK_ERROR(NULL == in,  "NULL 'in' pointer");
    CHECK_ERROR(NULL == out, "NULL 'out' pointer");
    CHECK_ERROR(in->flag != DPAL_GLOBAL
		&& in->flag != DPAL_GLOBAL_END
		&& in->flag != DPAL_LOCAL_END
		&& in->flag != DPAL_LOCAL,
		"Illegal flag");
    if (in->check_chars) {
	CHECK_ERROR(illegal_char(X, in->ssm, &msg[28]), msg);
	CHECK_ERROR(illegal_char(Y, in->ssm, &msg[28]), msg);
    }

    xlen = ustrlen(X);
    ylen = ustrlen(Y);

    out->align_end_1 = -1;
    out->align_end_2 = -1;

    if ('\0' == *X) {
	out->msg = "Empty first sequence";
	out->score = 0;
	return 1;
    }
    if ('\0' == *Y) {
	out->msg = "Empty second sequence";
	out->score = 0;
	return 1;
    }
    CHECK_ERROR(in->debug != 0 && in->score_only != 0,
		"score_only must be 0 if debug is non-0");
    if (1 == in->force_generic || in->debug == 1 || 0 == in->score_only) {
	/* 
	 * A true value of in->debug really means "print alignment on stderr"
	 * and implies 0 == a.score_only.
	 */
	CHECK_ERROR(xlen > DPAL_MAX_ALIGN,
		    "Sequence 1 longer than DPAL_MAX_ALIGN and alignment is requested");
	CHECK_ERROR(ylen > DPAL_MAX_ALIGN,
		    "Sequence 2 longer than DPAL_MAX_ALIGN and alignment is requested");
	if (_dpal_generic(X, Y, xlen, ylen, in, out))
	    goto FAIL;
    } else if (1 == in->force_long_generic) {
	_dpal_long_nopath_generic(X, Y, xlen, ylen, in, out);
    } else if (1 == in->max_gap ) {
	if (DPAL_LOCAL == in->flag)
	    _dpal_long_nopath_maxgap1_local(X, Y, xlen, ylen, in, out);
        else if (DPAL_GLOBAL_END == in->flag)
	    _dpal_long_nopath_maxgap1_global_end(X, Y, xlen, ylen, in, out);
	else if (DPAL_LOCAL_END == in->flag)
	    _dpal_long_nopath_maxgap1_local_end(X, Y, xlen, ylen, in, out);
        else if (xlen <= DPAL_MAX_ALIGN && ylen <= DPAL_MAX_ALIGN) {
	    if (_dpal_generic(X, Y, xlen, ylen, in, out))
		goto FAIL;
        } else _dpal_long_nopath_generic(X, Y, xlen, ylen, in, out);
    }
    else if (xlen < DPAL_MAX_ALIGN && ylen < DPAL_MAX_ALIGN) {
	if (_dpal_generic(X, Y, xlen, ylen, in, out))
	    goto FAIL;
    } else
      _dpal_long_nopath_generic(X, Y, xlen, ylen, in, out);

    return 0;
 FAIL:
    if (in->fail_stop) {
	fprintf(stderr, "\n%s\n", out->msg);
	exit(-1);
    }
    return 1;
}

static int
_dpal_generic(X, Y, xlen, ylen, in, out)
    const unsigned char *X,*Y;
    const int xlen, ylen;
    const dpal_args *in;
    dpal_results *out;
{

    /* The "score matrix" (matrix of best scores). */
    static int S[DPAL_MAX_ALIGN][DPAL_MAX_ALIGN];

#if NEW_TEST
    static int * S[DPAL_MAX_ALIGN] = NULL;
    if (NULL == S) 
      S = safe_malloc(ylen * sizeof(int) * DPAL_MAX_ALIGN);
    else
      S = safe_realloc(ylen * sizeof(int) * DPAL_MAX_ALIGN);
#endif

#if NEW_TEST2
    static int ** S = NULL;
    if (NULL == S) 
      S = safe_malloc(xlen * sizeof(int *));
    else {
      S = safe_realloc(xlen * sizeof(int *));
    }
    for (i = 0; i < xlen; i++)
      S[i] = safe_malloc(ylen * sizeof(int));
    /* Add free at end of function. */
#endif

#ifndef DPAL_FORGET_PATH
    /* The matrix of "trace" pointers */
    static int P[DPAL_MAX_ALIGN][DPAL_MAX_ALIGN][3];
#endif

    register int i, j, k, mg, c;
    register int gap = in->gap, gapl = in->gapl, max_gap = in->max_gap;

#ifndef DPAL_FORGET_PATH
    int i0,j0;
    int saved_k;
#endif 

    int I,J;		/* Coordinates of the maximum score. */
    int smax;           /* The optimum score. */
    int score;          /* Current score. */

    int a,b,max;

#ifdef DPAL_PRINT_COVERAGE
    fprintf(stderr, "_dpal_generic called\n");
#endif

    CHECK_ERROR(xlen > DPAL_MAX_ALIGN,
		"First sequence too long for _dpal_generic");
    CHECK_ERROR(ylen > DPAL_MAX_ALIGN,
		"Second sequence too long for _dpal_generic");

    /* Initialize the 0th column of the score matrix. */
    smax = INT_MIN;
    for(i=0; i < xlen; i++) {
	score = in->ssm[X[i]][Y[0]]; 
	if (DPAL_LOCAL == in->flag) {
	    if (score < 0) score = 0;
	    if(score > smax) {
		smax = score;
		I=i; J=0;
	    }
	}
	else if (DPAL_LOCAL_END == in->flag) {if (score < 0) score = 0;}
	S[i][0] = score;
    }	
    /* Move code for find global-alignment and end-anchored
       alignment below? */
    if (DPAL_LOCAL != in->flag) {
	/* 
	 * For a non-local alignment we restrict our search for the maximum
	 * score to the last row.
	 */
	smax = S[xlen-1][0]; I=xlen-1; J=0;
    }
	   
    /* Initialize the 0th row of the score matrix. */
    for(j=0; j<ylen; j++) { 
	score = in->ssm[X[0]][Y[j]]; 

	if(DPAL_LOCAL == in->flag){
	    if (score < 0) score = 0;
	    if(score > smax){
		smax = score;
		I=0; J=j;
	    }
	}
	else if (DPAL_LOCAL_END == in->flag) {if (score < 0) score = 0;}
	S[0][j] = score;
    }	
    if(DPAL_GLOBAL == in->flag&&S[0][ylen-1]>smax){
		smax = S[0][ylen-1];
		I=0; J=ylen-1;
    }

    /* Further is the solution for dynamic programming problem. */
    for(i=1; i<xlen; i++) {
	for(j=1; j<ylen; j++) {

	    a=S[i-1][j-1];

	    b = c = INT_MIN;
	    if (1 == max_gap) {
		if (i > 1) {
		    b = S[i-2][j-1] + gap;
#ifndef DPAL_FORGET_PATH
		    i0 = i - 2;
#endif
		}
		if (j > 1) {
		    c = S[i-1][j-2] + gap;
#ifndef DPAL_FORGET_PATH
		    j0 = j - 2;
#endif
		}
	    } else if (max_gap > 1) {
		max=INT_MIN;
		mg=(max_gap+1>i||max_gap<0)?i:max_gap+1;
		for(k=2; k<=mg; k++) {
		    c = S[i-k][j-1] + gap + gapl*(k-2);
		    if(c>max){
			max=c;
#ifndef DPAL_FORGET_PATH
			i0 = i-k;
#endif
		    }
		}
		b=max;

		max=INT_MIN;
		mg=(max_gap+1>j||max_gap<0)?j:max_gap+1;
		for(k=2;k<=mg;k++) {
		    c = S[i-1][j-k] + gap + gapl*(k-2);
		    if(c>max){
			max=c;
#ifndef DPAL_FORGET_PATH
			j0 = j-k;
#endif
		    }
		}
		c=max;
	    }

	    if(a>=b && a>=c) {
		score = a + in->ssm[X[i]][Y[j]];
#ifndef DPAL_FORGET_PATH
		P[i][j][1] = i-1;
		P[i][j][2] = j-1;
#endif
	    } else if (b > a && b >= c) {
		score = b + in->ssm[X[i]][Y[j]];
#ifndef DPAL_FORGET_PATH
		P[i][j][1] = i0;
		P[i][j][2] = j-1;
#endif
	    } else if (c > a && c > b) {
		score = c + in->ssm[X[i]][Y[j]];
#ifndef DPAL_FORGET_PATH
		P[i][j][1] = i-1;
		P[i][j][2] = j0;
#endif
	    }

	    if (score >= smax)
		/* 
		 * Because of comparison '>=' immediately above, dpal reports
		 * ungapped (i.e. diagonal) alignments if there is a choice
		 * of more than one optimum alignment.
		 */
		/* Move code to get 'g' and 'e' maxima to a separate loop ? */
		if (DPAL_LOCAL == in->flag 
		    || (DPAL_GLOBAL_END == in->flag && i == xlen-1)
		    || (DPAL_LOCAL_END   == in->flag && i == xlen-1)
		    || (DPAL_GLOBAL == in->flag&& (i==xlen-1||j==ylen-1))) {
		    /*  
		     * If in->flag is DPAL_LOCAL, then a cell anywhere within
		     * S may be the endpoint of the alignment.  If in->flag is
		     * DPAL_GLOBAL_END, then only cells in the last row may be
		     * the endpoint.  If in->flag is DPAL_GLOBAL cells in the
		     * last row _or_ the last column may be the endpoint.
		     */
		    smax = score;
		    I = i;
		    J = j;
	    } /*  put else here ? */
	    if (score < 0 && (DPAL_LOCAL == in->flag 
			      || DPAL_LOCAL_END == in->flag))
		/* 
		 * For a local alignment, 0 is the lowest score that we record
		 * in S.
		 */
		score = 0;

	    S[i][j]=score;
	}
    }
    /* I and J now specify the last pair of an optimum alignment. */

#ifndef DPAL_FORGET_PATH    
    k = (I > J) ? I+1 : J+1;
    saved_k=k;

    out->path[k][0]=I; out->path[k][1]=J;
    while(out->path[k][0]!=0&&out->path[k][1]!=0) {
	if ((in->flag== DPAL_LOCAL || in->flag == DPAL_LOCAL_END)
		 &&S[out->path[k][0]][out->path[k][1]]==0) {
	  k++; break;
	}
	out->path[k-1][0] = P[out->path[k][0]][out->path[k][1]][1];
	out->path[k-1][1] = P[out->path[k][0]][out->path[k][1]][2];
	k--;
    }
    if (k>0) {
	for (i=0;i<=saved_k-k;i++) {
	    out->path[i][0] = out->path[i+k][0];
	    out->path[i][1] = out->path[i+k][1];
	}
    }
#endif

    if ((DPAL_LOCAL == in->flag 
	 || DPAL_LOCAL_END == in->flag)&& S[I][J] <= 0) {
	/* There is no alignment at all. */
	out->score = 0;
	out->path_length = 0;
    } else {
	out->score = smax;
	out->align_end_1 = I;
	out->align_end_2 = J;
#ifndef DPAL_FORGET_PATH	
	out->path_length = saved_k - k + 1;
#else
	out->path_length = 0;
#endif
    }
#ifndef DPAL_FORGET_PATH
    if (in->debug) print_align(X,Y,P,I,J, in);
#endif
    return 0;
 FAIL:
    if (in->fail_stop) {
	fprintf(stderr, "\n%s\n", out->msg);
	exit(-1);
    }
    return 1;
} /* _dpal_generic */

/* Linear space, no path, for any value of maxgap and for any alignment. */
static void
_dpal_long_nopath_generic(X, Y, xlen, ylen, in, out)
    const unsigned char *X,*Y;
    const int xlen, ylen;
    const dpal_args *in;
    dpal_results *out;
{
    /* The "score matrix" (matrix of best scores). */
    int **S, *SI, **P; 

    register int i, j, k, mg, mgy, c;
    register int gap = in->gap, gapl = in->gapl, max_gap = in->max_gap;


    int I,J;		/* Coordinates of the maximum score. */
    int smax;           /* The optimum score. */
    int score;          /* Current score. */

#ifdef DPAL_PRINT_COVERAGE
    fprintf(stderr, "_dpal_long_nopath_generic called\n");
#endif

    out->score = INT_MIN;
    out->path_length = 0;
    out->msg = NULL;

    P = malloc(sizeof(*S)*(max_gap+2));
    S = malloc(sizeof(*S)*(max_gap+2));
    for(i=0; i<max_gap+2; i++){
	    P[i] = malloc(sizeof(SI)*xlen);
	    S[i] = P[i];
    }

    /* Initialize the 0th column of the score matrix. */
    smax = INT_MIN;
    for(i=0; i < xlen; i++) {
	score = in->ssm[X[i]][Y[0]]; 
	if (DPAL_LOCAL == in->flag) {
	    if (score < 0) score = 0;
	    if(score > smax) {
		smax = score;
		I=i; J=0;
	    }
	}
	else if (DPAL_LOCAL_END == in->flag) {if (score < 0) score = 0;}
	S[0][i] = score;
    }	
    /* Move code for find global-alignment and end-anchored
       alignment below? */
    if (DPAL_LOCAL != in->flag) {
	/* 
	 * For a non-local alignment we restrict our search for the maximum
	 * score to the last row.
	 */
	smax = S[0][xlen-1]; I=xlen-1; J=0;
    }
	   
    /* Initialize the 0th row of the score matrix. 
    for(j=0; j<ylen; j++) { 
	score = in->ssm[X[0]][Y[j]]; 

	if(DPAL_LOCAL == in->flag){
	    if (score < 0) score = 0;
	    if(score > smax){
		smax = score;
		I=0; J=j;
	    }
	}
	S[0][j] = score;
    }	
    
    if(DPAL_GLOBAL == in->flag&&S[0][ylen-1]>smax){
		smax = S[0][ylen-1];
		I=0; J=ylen-1;
    }
    */

    /* Further is the solution for dynamic programming problem. */
    for(j=1; j<ylen; j++) {
	mgy=(max_gap+1>j||max_gap<0)?j:max_gap+1;
	score = in->ssm[X[0]][Y[j]];
	 if (DPAL_LOCAL == in->flag) {
	     if (score < 0) score = 0;
	     if(score > smax) smax = score;
         }    
	 else if (DPAL_LOCAL_END == in->flag) { if (score < 0) score = 0;}
	 else if (DPAL_GLOBAL == in->flag && j == ylen-1 && score > smax)
			smax = score;
	S[mgy][0] = score;
	for(i=1; i<xlen; i++) {

	    score=S[mgy-1][i-1];

		mg=(max_gap+1>i||max_gap<0)?i:max_gap+1;
		for(k=2; k<=mg; k++) 
		    if((c = S[mgy-1][i-k] + gap + gapl*(k-2)) > score)score = c;

		for(k=2;k<=mgy;k++) 
		    if((c = S[mgy-k][i-1] + gap + gapl*(k-2)) > score)score=c;

		score += in->ssm[X[i]][Y[j]];

	    if (score >= smax)
		/* 
		 * Because of comparison '>=' immediately above, dpal reports
		 * ungapped (i.e. diagonal) alignments if there is a choice
		 * of more than one optimum alignment.
		 */
		/* Move code to get 'g' and 'e' maxima to a separate loop ? */
		if (DPAL_LOCAL == in->flag 
		    || ((DPAL_GLOBAL_END == in->flag
			 || DPAL_LOCAL_END == in->flag) 
			&& i == xlen-1)
		    || (DPAL_GLOBAL == in->flag&& (i==xlen-1||j==ylen-1))) {
		    /*  
		     * If in->flag is DPAL_LOCAL, then a cell anywhere within
		     * S may be the endpoint of the alignment.  If in->flag is
		     * DPAL_GLOBAL_END, then only cells in the last row may be
		     * the endpoint.  If in->flag is DPAL_GLOBAL cells in the
		     * last row _or_ the last column may be the endpoint.
		     */
		    smax = score;
		    I = i;
		    J = j;
	    } /*  put else here ? */
	    if (score < 0 && (DPAL_LOCAL == in->flag
			      || DPAL_LOCAL_END == in->flag))
		/* 
		 * For a local alignment, 0 is the lowest score that we record
		 * in S.
		 */
		score = 0;

	    S[mgy][i]=score;
	}
	if(mgy == max_gap + 1){
 	    SI = S[0];
	    for(i=0; i<mgy; i++)S[i] = S[i+1];
	    S[mgy] = SI;
        }
    }
    /* I and J now specify the last pair of an optimum alignment. */

    if (DPAL_LOCAL == in->flag && smax <= 0) {
	/* There is no alignment at all. */
	out->score = 0;
	out->path_length = 0;
    } else {
	out->score = smax;
	out->align_end_1 = I;
	out->align_end_2 = J;
    }
    for(i=0; i< max_gap + 2; i++) free(P[i]);
    free(S);
    free(P);
} /* _dpal_long_nopath_generic */

static void
_dpal_long_nopath_maxgap1_local(X, Y, xlen, ylen, in, out)
    const unsigned char *X,*Y;
    const int xlen, ylen;
    const dpal_args *in;
    dpal_results *out;
{
    /* The "score matrix" (matrix of best scores). */
    int *S0, *S1, *S2; 
    int *P0, *P1, *P2;
    int *S;

    register int i, j;
    register int gap = in->gap;
    register int smax;           /* The optimum score. */
    register int score;          /* Current score. */
    register int a;

#ifdef DPAL_PRINT_COVERAGE
    fprintf(stderr, "_dpal_long_nopath_maxgap1_local called\n");
#endif

    P0 = malloc(sizeof(int)*ylen);
    P1 = malloc(sizeof(int)*ylen);
    P2 = malloc(sizeof(int)*ylen);

    S0 = P0; S1 = P1; S2 = P2;

    smax = 0; /* For local alignment score can never be less than 0. */

    /* Initialize the 0th row of the score matrix. */
    for(j=0; j < ylen; j++) { 
	score = in->ssm[X[0]][Y[j]]; 
	if (score < 0) score = 0;
	else if (score > smax) smax = score;
	/*S[0][j] = score;*/
	S0[j] = score;
    }	

    /* Set the 1st row of the score matrix. */
    score = in->ssm[X[1]][Y[0]];
    if(score < 0) score = 0;
    else if (score > smax) smax = score;
    S1[0] = score;
    for(j=1; j < ylen; j++) {
	score = S0[j-1];
	if(j>1 && (a=S0[j-2] + gap) > score)score = a;
	score += in->ssm[X[1]][Y[j]];
	if (score < 0) score = 0;
	else if(score > smax) smax = score;
	S1[j] = score;
    }

    for(i=2; i < xlen; i++) {
	score = in->ssm[X[i]][Y[0]];
	if (score < 0) score = 0;
	else if (score > smax) smax = score;
	S2[0] = score;
	score = S1[0];
	if((a=S0[0] + gap) > score) score = a;
	score += in->ssm[X[i]][Y[1]];
	if(score < 0) score = 0;
	else if (score > smax) smax = score;
	S2[1] = score;
	for(j=2; j < ylen; j++) {
	    score = S0[j-1];
	    if((a=S1[j-2])>score) score = a;
	    score +=gap;
	    if((a=S1[j-1]) >score) score = a;

	    score += in->ssm[X[i]][Y[j]];	
	    if (score < 0 ) score = 0;
	    else if (score > smax) smax = score;
	    S2[j]=score;
	}
	S = S0; S0 = S1; S1 = S2; S2 = S;
    }
    out->score = smax;
    out->path_length=0;
    free(P0); free(P1); free(P2);
} /* _dpal_long_nopath_maxgap1_local */

static void
_dpal_long_nopath_maxgap1_global_end(X, Y, xlen, ylen, in, out)
    const unsigned char *X,*Y;
    const int xlen, ylen;
    const dpal_args *in;
    dpal_results *out;
{
    /* The "score matrix" (matrix of best scores). */
    int *S0, *S1, *S2, *S; 
    int *P0, *P1, *P2;

    register int i, j, k;
    register int gap = in->gap;
    register int smax;           /* The optimum score. */
    register int score;          /* Current score. */
    register int a, t;

#ifdef DPAL_PRINT_COVERAGE
    fprintf(stderr, "_dpal_long_nopath_maxgap1_global_end called\n");
#endif

    P0 = malloc(sizeof(int)*xlen);
    P1 = malloc(sizeof(int)*xlen);
    P2 = malloc(sizeof(int)*xlen);

    S0 = P0; S1 = P1; S2 = P2;

    smax = in->ssm[X[xlen-1]][Y[0]];
	   
    /* Set the 0th row of the score matrix. */
    for(j=0; j<xlen; j++) S0[j] = in->ssm[X[j]][Y[0]]; 

    /* Set the 1st row of the score matrix. */
    S1[0] = in->ssm[X[0]][Y[1]];
    for(j=1; j < xlen; j++){
	score = S0[j-1];
	if(j>1 && (a=S0[j-2] + gap)> score)score = a;
	score += in->ssm[X[j]][Y[1]];
	if(score > smax && j == xlen-1) smax = score;
        S1[j] = score;
    }

    k = ylen - (int)(xlen / 2) + 1;
    if (k<1) k = 1;

    /* Set the rectangular part of almost the remainder of the matrix. */
    for(j=2; j<k+1; j++) {
	S2[0] = in->ssm[X[0]][Y[j]];
	score = S1[0];
	if((a=S0[0]+gap) > score) score = a;
	score += in->ssm[X[1]][Y[j]];
	S2[1] = score;
       for(i=2; i<xlen-1; i++) {
	   score = S1[i-2];
	   if((a=S0[i-1]) > score)score = a;
	   score += gap;
	   if((a=S1[i-1]) > score)score = a;
	   score += in->ssm[X[i]][Y[j]];
	   S2[i] = score;
        }
	score = S1[xlen-3];
	if((a=S0[xlen-2]) > score)score = a;
	score += gap;
	if((a=S1[xlen-2]) > score)score = a;
	score += in->ssm[X[xlen-1]][Y[j]];
	S2[xlen-1] = score;
	if(score > smax) smax = score;
	S = S0; S0 = S1; S1 = S2; S2 = S;
    }

    /* Set the triangular part of almost the remainder of the matrix. */
    t = 2;
    for(j=k+1; j<ylen; j++) {
       for(i=t; i<xlen-1; i++) {
	    score = S1[i-2];
	    if((a=S0[i-1]) > score) score = a;
	    score += gap;
	    if((a=S1[i-1]) > score) score = a;
	    score += in->ssm[X[i]][Y[j]];
	    S2[i] = score;
       }
       t += 2;
       score = S1[xlen-3];
       if((a=S0[xlen-2]) > score)score = a;
       score += gap;
       if((a=S1[xlen-2]) > score)score = a;
       score += in->ssm[X[xlen-1]][Y[j]];
       S2[xlen-1] = score;
       if(score > smax) smax = score;
       S = S0; S0 = S1; S1 = S2; S2 = S;
    }

    free(P0); free(P1); free(P2);
    out->score = smax;
    out->path_length=0;
} /* _dpal_long_nopath_maxgap_global_end */


#ifndef DPAL_FORGET_PATH
/* Reconstruct the best path and print the alignment. */
void
print_align(X,Y,P,I,J, dargs)
    const unsigned char *X,*Y;
    int P[DPAL_MAX_ALIGN][DPAL_MAX_ALIGN][3], I, J;
    const dpal_args *dargs;
{
	int JX[DPAL_MAX_ALIGN],JY[DPAL_MAX_ALIGN];
	int k,i,j,n,m;
	char sx[3*DPAL_MAX_ALIGN],sy[3*DPAL_MAX_ALIGN],sxy[3*DPAL_MAX_ALIGN];

	for(i=0;i<3*DPAL_MAX_ALIGN;i++){
		sx[i] = ' ';sy[i] = ' '; sxy[i] = ' ';
	}
	if(I>J)k=I+1;
	else k=J+1;

	n=k;
	JX[k] = I;
	JY[k] = J;
	while(JX[k]!=0&&JY[k]!=0){
		JX[k-1] = P[JX[k]][JY[k]][1];
		JY[k-1] = P[JX[k]][JY[k]][2];
		k--;
	}
	if(JX[k]>JY[k]){
		for(i=0;i<JX[k];i++)sx[i] = X[i];
		for(i=0;i<JX[k]-JY[k];i++)sy[i] = ' ';
		j = JX[k]-JY[k];
		for(i=JX[k]-JY[k];i<JX[k];i++)sy[i] = Y[i-j];
		m = JX[k];
	}
	else{
		for(i=0;i<JY[k];i++)sy[i] = Y[i];
		for(i=0;i<JY[k]-JX[k];i++)sx[i] = ' ';
		j= JY[k]-JX[k];
		for(i=j;i<JY[k];i++)sx[i] = X[i-j];
		m = JY[k];
	}
	for(i=0;i<m;i++)sxy[i] = ' ';
	for(i=k;i<n;i++){
		sx[m] = X[JX[i]];
		sy[m] = Y[JY[i]];
		/* if(sx[m]==sy[m]&&sx[m]!='N') sxy[m] = '|'; */
		if (dargs->ssm[(unsigned char)sx[m]][(unsigned char)sy[m]] > 0)
		  sxy[m] = '|';
		else sxy[m]=' ';
		if(JX[i+1]-JX[i]>JY[i+1]-JY[i]){
			for(j=1;j<JX[i+1]-JX[i];j++){
				sy[m+j] = '-';
				sx[m+j] = X[JX[i]+j];
				sxy[m+j] = ' ';
			}
			m += JX[i+1]-JX[i]-1;
		}
		if(JY[i+1]-JY[i]>JX[i+1]-JX[i]){
			for(j=1;j<JY[i+1]-JY[i];j++){
				sx[m+j] = '-';
				sy[m+j] = Y[JY[i]+j];
				sxy[m+j] = ' ';
			}
			m += JY[i+1]-JY[i]-1;
		}
		m++;
	}
	sx[m] = X[I];
	sy[m] = Y[J];
	for(i=m+1;i<m+ustrlen(X)-I;i++)sx[i]=X[i-m+I];
	for(i=m+1;i<m+ustrlen(Y)-J;i++)sy[i]=Y[i-m+J];
	/* if (sx[m]==sy[m] && sx[m]!='N') sxy[m]='|'; */
	if (dargs->ssm[(unsigned char)sx[m]][(unsigned char)sy[m]] > 0)
	  sxy[m] = '|';
	else sxy[m]=' ';
	m++;
	if(ustrlen(X)-I>ustrlen(Y)-J)k=m+ustrlen(X)-I;
	else k=m+ustrlen(Y)-J;
	
	j=0;
	while(j<k){
		for(i=j;i<j+70;i++) fprintf(stderr, "%c",sx[i]);
		fprintf(stderr, "\n");
		for(i=j;i<j+70;i++) fprintf(stderr, "%c",sxy[i]);
		fprintf(stderr, "\n");
		for(i=j;i<j+70;i++) fprintf(stderr, "%c",sy[i]); fprintf(stderr,"\n");
		for(i=0;i<70;i++)   fprintf(stderr, "_");
		fprintf(stderr, "\n");
		j +=70;
	}
}
#endif

static int
illegal_char(X, ssm, out)
    const unsigned char *X;
    const dpal_ssm ssm;
    char *out;
{
    register const unsigned char *p;
    for (p = X; *p != '\0' && ssm[*p][*p] != INT_MIN; p++);
    if (*p == '\0')
	return 0;
    else {
	*out = *p;
	return 1;
    }
}
static void
_dpal_long_nopath_maxgap1_local_end(X, Y, xlen, ylen, in, out)
    const unsigned char *X,*Y;
    const int xlen, ylen;
    const dpal_args *in;
    dpal_results *out;
{
    /* The "score matrix" (matrix of best scores). */
    int *S0, *S1, *S2; 
    int *P0, *P1, *P2;
    int *S;

    register int i, j;
    register int gap = in->gap;
    register int smax;           /* The optimum score. */
    register int score;          /* Current score. */
    register int a;

#ifdef DPAL_PRINT_COVERAGE
    fprintf(stderr, "_dpal_long_nopath_maxgap1_local_end called\n");
#endif

    P0 = malloc(sizeof(int)*ylen);
    P1 = malloc(sizeof(int)*ylen);
    P2 = malloc(sizeof(int)*ylen);

    S0 = P0; S1 = P1; S2 = P2;

    smax = 0; /* For local alignment score can never be less than 0. */

    /* Initialize the 0th row of the score matrix. */
    for(j=0; j < ylen; j++) { 
	score = in->ssm[X[0]][Y[j]]; 
	if (score < 0) score = 0;
	/*S[0][j] = score;*/
	S0[j] = score;
    }	

    /* Set the 1st row of the score matrix. */
    score = in->ssm[X[1]][Y[0]];
    if(score < 0) score = 0;
    S1[0] = score;
    for(j=1; j < ylen; j++) {
	score = S0[j-1];
	if(j>1 && (a=S0[j-2] + gap) > score)score = a;
	score += in->ssm[X[1]][Y[j]];
	if (score < 0) score = 0;
	S1[j] = score;
    }

    for(i=2; i < xlen - 1; i++) {
	score = in->ssm[X[i]][Y[0]];
	if (score < 0) score = 0;
	S2[0] = score;
	score = S1[0];
	if((a=S0[0] + gap) > score) score = a;
	score += in->ssm[X[i]][Y[1]];
	if(score < 0) score = 0;
	S2[1] = score;
	for(j=2; j < ylen; j++) {
	    score = S0[j-1];
	    if((a=S1[j-2])>score) score = a;
	    score +=gap;
	    if((a=S1[j-1]) >score) score = a;

	    score += in->ssm[X[i]][Y[j]];	
	    if (score < 0 ) score = 0;
	    S2[j]=score;
	}
	S = S0; S0 = S1; S1 = S2; S2 = S;
    }
    /* Calculate scores for last row (i = xlen-1) and find smax */
    i = xlen - 1;
    score = in->ssm[X[i]][Y[0]];
    if (score < 0) score = 0;
    else if (score > smax) smax = score;
    S2[0] = score;
    score = S1[0];
    if((a=S0[0] + gap) > score) score = a;
    score += in->ssm[X[i]][Y[1]];
    if(score < 0) score = 0;
    else if (score > smax) smax = score;
    S2[1] = score;
    for(j=2; j < ylen; j++) {
	score = S0[j-1];
	if((a=S1[j-2])>score) score = a;
	score +=gap;
	if((a=S1[j-1]) >score) score = a;
	score += in->ssm[X[i]][Y[j]];
	if (score < 0 ) score = 0;
	else if (score > smax) smax = score;
	S2[j]=score;
    }
    out->score = smax;
    out->path_length=0;
    free(P0); free(P1); free(P2);
} /* _dpal_long_nopath_maxgap1_local */
