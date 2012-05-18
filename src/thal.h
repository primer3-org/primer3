/* Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008,2009,2010
                 2011,2012
 Whitehead Institute for Biomedical Research, Steve Rozen
 (http://purl.com/STEVEROZEN/), and Helen Skaletsky
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
 THEORY OF LIABILITY, WHETHER IN CONTRACT, ST
 RICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */

#ifndef _THAL_H
#define _THAL_H

#include <float.h> /* ! mul ei ole float.h-d includes DBL_MAX */
#include <math.h>
#include <limits.h>


#ifndef THAL_ERROR_SCORE
# define THAL_ERROR_SCORE -_INFINITY
#endif

/* The maximum length of _one_ of the two sequences being aligned in a
   thermodynamic alignment. In other words, the length of one sequence
   must be <= THAL_MAX_ALIGN, but the other sequence can be longer.
   The rationale behind this value (60) is that this is the maxium
   reasonable length for nearest neighbor models. It is the maxium
   length at which we can restrict our model to only two states of
   melting: fully intact duplex or completely dissociated single
   strands. */
#ifndef THAL_MAX_ALIGN
#define THAL_MAX_ALIGN 60
#endif

/* The maxium length of the other sequence in a thermodynamic
   alignment. This value can be increased, though alignments against
   very long sequences will be quite slow. As of 2012-05-18, we only
   potentially see sequences longer this when checking for mispriming
   in the template ('max_template_mispriming') in libprimer3.c, which
   is really designed to find sites of ectopic primer very close (a
   few kilobases) from the location of the cadidate primer. */
#ifndef THAL_MAX_SEQ
#define THAL_MAX_SEQ   10000
#endif

/*** BEGIN CONSTANTS ***/

extern const double _INFINITY;
extern const double ABSOLUTE_ZERO;
extern const int MAX_LOOP; /* the maximum size of loop that can be calculated;
			    for larger loops formula must be implemented */
extern const int MIN_LOOP;

/*** END CONSTANTS ***/

/* BEGIN TYPEDEFs */

typedef enum thal_alignment_type {
  thal_any = 1,
  thal_end1 = 2,
  thal_end2 = 3,
  thal_hairpin = 4,
} thal_alignment_type;

/* Structure for passing arguments to THermodynamic ALignment calculation */
typedef struct {
   int debug; /* if non zero, print debugging info to stderr */
   thal_alignment_type type; /* one of the
	      1 THAL_ANY, (by default)
	      2 THAL_END1,
	      3 THAL_END2,
	      4 THAL_HAIRPIN */
   int maxLoop;  /* maximum size of loop to consider; longer than 30 bp are not allowed */
   double mv; /* concentration of monovalent cations */
   double dv; /* concentration of divalent cations */
   double dntp; /* concentration of dNTP-s */
   double dna_conc; /* concentration of oligonucleotides */
   double temp; /* temperature from which hairpin structures will be calculated */
   int temponly; /* if non zero, print only temperature to stderr */
   int dimer; /* if non zero, dimer structure is calculated */
} thal_args;

/* Structure for receiving results from the thermodynamic alignment calculation */
typedef struct {
   char msg[255];
   double temp;
   int align_end_1;
   int align_end_2;
} thal_results;

/*** END OF TYPEDEFS ***/

void set_thal_default_args(thal_args *a);
void set_thal_oligo_default_args(thal_args *a);

/* Read the thermodynamic values (parameters) from the parameter files
   in the directory specified by 'path'.  Return 0 on success and -1
   on error. The thermodynamic values are stored in multiple static
   variables. */
/* Here is an example of how this function is used in 
   primer3_boulder_main.c: */
#if 0
  if (get_thermodynamic_values(thermodynamic_params_path, &o)) {
    fprintf(stderr, "%s\n", o.msg);
    exit(-1);
  }
#endif
int  get_thermodynamic_values(const char* path, thal_results *o);

void destroy_thal_structures();

/* Central method for finding the best alignment.  On error, o->temp
   is set to THAL_ERROR_SCORE and a message is put in o->msg.  The
   error might be caused by ENOMEM. To determine this it is necessary
   to check errno.
*/

void thal(const unsigned char *oligo1, 
	  const unsigned char *oligo2, 
	  const thal_args* a, 
	  thal_results* o);

#endif
