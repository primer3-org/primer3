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

#ifndef PRECISION
# define PRECISION 1
#endif

#ifndef THAL_ERROR_SCORE
# define THAL_ERROR_SCORE -INFINITY
#endif

#ifndef THAL_EXIT_ON_ERROR
#define THAL_EXIT_ON_ERROR 0
#endif

#ifndef THAL_MAX_ALIGN
#define THAL_MAX_ALIGN 60
#endif

#ifndef THAL_MAX_SEQ
#define THAL_MAX_SEQ   1600
#endif

/*** BEGIN CONSTANTS ***/

extern const double INFINITY;
extern const double R; /* cal/Kmol */
extern const double ILAS; /* Internal Loop Entropy ASymmetry correction -0.3kcal/mol*/
extern const double ILAH; /* Internal Loop EntHalpy Asymmetry correction */
extern const double AT_H; /* Enthalpy of AT penalty */
extern const double AT_S; /* Entropy of AT penalty */
extern const double MinEntropyCutoff; /* to filter out unreal entropies - positive entropies 
				       - it can affect calculating Tm to unreal values */
extern const double MinEntropy; /* initiation */
extern const double G2; /* structures w higher G are considered to be more unstabile */
extern const double ABSOLUTE_ZERO;
extern const int MAX_LOOP; /* the maximum size of loop that can be calculated; 
			    for larger loops formula must be implemented */
extern const int MIN_LOOP;
extern const char BASES[5]; /* bases to be considered - N is every symbol that is not A, G, C, T */
extern const char BASE_PAIRS[4][4]; /* allowed basepairs */
extern const int BPI[5][5]; /* matrix for allowed; bp 0 - no bp, watson crick bp - 1 */
/* extern const char *msg; I think this is not needed */

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
   int fail_stop; /* exit with -1 on error */
} thal_args;

/* Structure for receiving results from the thermodynamic alignment calculation */
typedef struct {
   const char *msg;
   double temp;
   int align_end_1;
   int align_end_2;
} thal_results;

/*** END OF TYPEDEFS ***/

/*** BEGIN STRUCTs ***/

struct triloop { 
   char loop[5]; 
   double value; };

struct tetraloop { 
   char loop[6]; 
   double value; };


struct tracer /* structure for tracebacku - unimolecular str */ {
   int i;
   int j;
   int mtrx; /* [0 1] EntropyDPT/EnthalpyDPT*/
   struct tracer* next;
};

/*** END STRUCTs ***/

void set_thal_default_args(thal_args *a);

int length_unsig_char(const unsigned char * str); /* returns length of unsigned char; to avoid warnings while compiling */

/* Central method for finding the best alignment */
void thal(const unsigned char *oligo1, const unsigned char *oligo2, const thal_args* a, thal_results* o);

unsigned char str2int(char c); /* converts DNA sequence to int; 0-A, 1-C, 2-G, 3-T, 4-whatever */

int seqcmp(unsigned char* seq1, unsigned char* seq2, int length);

double saltCorrectS (double mv, double dv, double dntp); /* part of calculating salt correction 
							  for Tm by SantaLucia et al */
int smallest(int a, int b, int c);

int identic(unsigned char* a, unsigned char* b, int len); /* checks if two seq-s are identical */

FILE* openParamFile(char* name,thal_results* o,int a); /* file of thermodynamic params */

/* get thermodynamic tables */
void getStack(double stackEntropies[5][5][5][5], double stackEnthalpies[5][5][5][5],thal_results* o,int a);

void verifyStackTable(double stack[5][5][5][5], char* type); /* just for debugging; the method is turned off by default */

void getStackint2(double stackEntropiesint2[5][5][5][5], double stackint2Enthalpies[5][5][5][5],thal_results* o,int a);

void getDangle(double dangleEntropies3[5][5][5], double dangleEnthalpies3[5][5][5], double dangleEntropies5[5][5][5], 
		double dangleEnthalpies5[5][5][5],thal_results* o,int a);

void getLoop(double hairpinLoopEnntropies[30], double interiorLoopEntropies[30], double bulgeLoopEntropiess[30],
	      double hairpinLoopEnthalpies[30], double interiorLoopEnthalpies[30], double bulgeLoopEnthalpies[30],thal_results* o,int a);

void getTstack(double tstackEntropies[5][5][5][5], double tstackEnthalpies[5][5][5][5],thal_results* o,int a);

void getTstack2(double tstack2Entropies[5][5][5][5], double tstack2Enthalpies[5][5][5][5],thal_results* o,int a);

void getTriloop(struct triloop**, struct triloop**, int* num,thal_results* o,int a);

void getTetraloop(struct tetraloop**, struct tetraloop**, int* num,thal_results* o,int a);

void tableStartATS(double atp_value, double atp[5][5]); /* creates table of entropy values for nucleotides 
							 to which AT-penlty must be applied */

void tableStartATH(double atp_value, double atp[5][5]);

int comp3loop(const void*, const void*); /* checks if sequnece consists of specific triloop */

int comp4loop(const void*, const void*); /* checks if sequnece consists of specific tetraloop */

void initMatrix(); /* initiates thermodynamic parameter tables of entropy and enthalpy for dimer */

void initMatrix2(); /* initiates thermodynamic parameter tables of entropy and enthalpy for monomer */

void fillMatrix(int maxLoop,thal_results* o, int a); /* calc-s thermod values into dynamic progr table (dimer) */

void fillMatrix2(int maxLoop, thal_results* o, int a); /* calc-s thermod values into dynamic progr table (monomer) */

void maxTM(int i, int j); /* finds max Tm while filling the dyn progr table using stacking S and stacking H (dimer) */

void maxTM2(int i, int j); /* finds max Tm while filling the dyn progr table using stacking S and stacking H (monomer) */

/* calculates bulges and internal loops for dimer structures */
void calc_bulge_internal(int ii, int jj, int i, int j, double* EntropyEnthalpy, int traceback, int maxLoop);

/* calculates bulges and internal loops for monomer structures */
void calc_bulge_internal2(int ii, int jj, int i, int j, double* EntropyEnthalpy, int traceback, int maxLoop);

/* carries out Bulge and Internal loop and stack calculations to hairpin */
void CBI(int i, int j, double* EntropyEnthalpy, int traceback, int maxLoop);

/* finds monomer structure that has maximum Tm */
void calc_hairpin(int i, int j, double* EntropyEnthalpy, int traceback);

double Ss(int i, int j, int k); /* returns stack entropy */
double Hs(int i, int j, int k); /* returns stack enthalpy */

/* calculate terminal entropy S and terminal enthalpy H starting reading from 5'end (Left hand/3' end - Right end) */
void LSH(int i, int j, double* EntropyEnthalpy);
void RSH(int i, int j, double* EntropyEnthalpy);

void reverse(unsigned char *s);

int max5(double, double, double, double, double);

/* Is sequence symmetrical */
int symmetry_thermo(const unsigned char* seq);

/* traceback for dimers */
void traceback(int i, int j, double RT, int* ps1, int* ps2, int maxLoop,thal_results* o, int a);

/* traceback for hairpins */
void tracebacku(int*, int,thal_results*, int);

/* prints ascii output of dimer structure */
void drawDimer(int*, int*, double, double, double, int, double, thal_results *, int);

/* prints ascii output of hairpin structure */
void drawHairpin(int*, double, double, int, double, thal_results *, int);

int equal(double a, double b);

void strcatc(char*, char);

void push(struct tracer**, int, int, int,thal_results*,int); /* to add elements to struct */

/* terminal bp for monomer structure */
void calc_terminal_bp(double temp);

/* executed in calc_terminal_bp; to find structure that corresponds to max Tm for terminal bp */
double END5_1(int,int); /* END5_1(X,1/2) - 1=Enthalpy, 2=Entropy*/
double END5_2(int,int);
double END5_3(int,int);
double END5_4(int,int);

double Hd5(int,int); /* returns thermodynamic value (H) for 5' dangling end */
double Hd3(int,int); /* returns thermodynamic value (H) for 3' dangling end */
double Sd5(int,int); /* returns thermodynamic value (S) for 5' dangling end */
double Sd3(int,int); /* returns thermodynamic value (S) for 3' dangling end */
double Ststack(int,int); /* returns entropy value for terminal stack */
double Htstack(int,int); /* returns enthalpy value for terminal stack */

/* memory stuff */
#ifndef XMALLOC_H
# define XMALLOC_H
void* safe_calloc(size_t, size_t,thal_results* o, int a);
void* safe_malloc(size_t,thal_results* o, int a);
void* safe_realloc(void*, size_t,thal_results* o, int a);
double* safe_recalloc(double* ptr, int m, int n,thal_results* o, int a);
#endif

#endif
