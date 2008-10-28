/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008
Whitehead Institute for Biomedical Research, Steve Rozen
(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky
All rights reserved.

    This file is part the primer3 software suite.

    This software suite is free software;
    you can redistribute is and/or modify it under the terms
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

#ifndef _DPAL_H
#define _DPAL_H
#include <limits.h> /* INT_MIN */

#ifdef __cplusplus
    extern "C" {
#endif

#define DPAL_ERROR_SCORE INT_MIN

#define DPAL_EXIT_ON_ERROR 0
  /* 0 means do not exit on error. */

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
void dpal_set_h_nt_matrix(dpal_args *a);

/* The argument, a, must be a DNA scoring matrix.  Modifies a so that it for a
   match between any two ambiguity codes (or between ambiguity code and base),
   e.g. B and S, the score will be the maximum of score between any base in B
   and any base in S, in the example between any pair in {C, G, T} X {C, G}.
   This function overwrites any scores already associated with pairs of
   ambiguity codes.  Return 0 on error, 1 on success.
*/
int dpal_set_ambiguity_code_matrix(dpal_args *);

/* 
 * Align the first 2 arguments, using the scoring
 * matix and other arguments supplied in the dpal_args
 * argument.  Return results in argument 'out'.
 * On error, sets out->score to DPAL_ERROR_SCORE
 * and puts additional information in out->msg.
 */
void dpal(const unsigned char *, const unsigned char*,
	  const dpal_args *, dpal_results *out);

void set_dpal_args(dpal_args *);

#ifdef __cplusplus
    }
#endif

#endif
