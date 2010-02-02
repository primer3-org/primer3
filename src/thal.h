/* Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008,2009
 Whitehead Institute for Biomedical Research, Steve Rozen
 (http://jura.wi.mit.edu/rozen), and Helen Skaletsky
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

#ifndef THAL_MAX_ALIGN
#define THAL_MAX_ALIGN 60
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

int  get_thermodynamic_values(const char* path, thal_results *o);

void destroy_thal_structures();

/* Central method for finding the best alignment */
void thal(const unsigned char *oligo1, const unsigned char *oligo2, const thal_args* a, thal_results* o);

#endif
