/*
 Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2009
 Whitehead Institute for Biomedical Research, Steve Rozen
 (http://jura.wi.mit.edu/rozen), and Helen Skaletsky
 All rights reserved.
 
       This file is part of primer3 software suite.
   
       This software suite is is free software;
       you can redistribute it and/or modify it under the terms
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
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON A THEORY 
 OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <limits.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <setjmp.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>

#include "thal.h"
#define DEBUG
#ifndef MIN_HRPN_LOOP
#define MIN_HRPN_LOOP 3 /*  minimum size of hairpin loop */
#endif

/* table where bp-s enthalpies, that retrieve to the most stable Tm, are saved */
#ifdef EnthalpyDPT
# undef EnthalpyDPT
#endif
#define EnthalpyDPT(i, j) enthalpyDPT[(j) + ((i-1)*len3) - (1)]

/* table where bp-s entropies, that retrieve to the most stable Tm, are saved */
#ifdef EntropyDPT
# undef EntropyDPT
#endif
#define EntropyDPT(i, j) entropyDPT[(j) + ((i-1)*len3) - (1)]

/* entropies of most stable hairpin terminal bp */
#ifndef SEND5
# define SEND5(i) send5[i]
#endif

/* enthalpies of most stable hairpin terminal bp */
#ifndef HEND5
# define HEND5(i) hend5[i]
#endif

#define CHECK_ERROR(COND,MSG) if (COND){o->msg = MSG; goto FAIL;}
#define THAL_OOM_ERROR { o->msg = "Out of memory", errno=ENOMEM; goto FAIL; }

#define bpIndx(a, b) BPI[a][b] /* for traceing matrix BPI */
#define min2(a, b) ((a) < (b) ? (a) : (b))
#define max2(a, b) ((a) < (b) ? (2) : (1))
#define atPenaltyS(a, b) atpS[a][b]
#define atPenaltyH(a, b) atpH[a][b]

#define SMALL_NON_ZERO 0.000001
#define DBL_EQ(X,Y) (((X) - (Y)) < (SMALL_NON_ZERO) ? (1) : (2)) /* 1 when numbers are equal */

#ifdef INTEGER
# define isFinite(x) (x < INFINITY / 2)
#else
# define isFinite(x) finite(x)
#endif

/*** BEGIN CONSTANTS ***/
# ifdef INTEGER
const double INFINITY = 999999.0;
# else
const double INFINITY = 1.0 / 0.0;
# endif

const double R = 1.9872; /* cal/Kmol */
const double ILAS = (-300 / 310.15); /* Internal Loop Entropy ASymmetry correction -0.3kcal/mol*/
const double ILAH = 0.0; /* Internal Loop EntHalpy Asymmetry correction */
const double AT_H = 2200.0; /* AT penalty */
const double AT_S = 6.9; /* AT penalty */
const double MinEntropyCutoff = -2500.0; /* to filter out non-existing entropies */
const double MinEntropy = -3224.0; /* initiation */
const double G2 = 0.0; /* structures w higher G are considered to be unstabile */
const double ABSOLUTE_ZERO = 273.15;
const int MAX_LOOP = 30; /* the maximum size of loop that can be calculated; for larger loops formula must be implemented */
const int MIN_LOOP = 0;
const char BASES[5] = {'A', 'C', 'G', 'T', 'N'}; /* bases to be considered - N is every symbol that is not A, G, C,$
						  */ 
const char BASE_PAIRS[4][4] = {"A-T", "C-G", "G-C", "T-A" }; /* allowed basepairs */
/* matrix for allowed; bp 0 - no bp, watson crick bp - 1 */
const int BPI[5][5] =  {
     {0, 0, 0, 1, 0}, /* A, C, G, T, N; */
     {0, 0, 1, 0, 0},
     {0, 1, 0, 0, 0},
     {1, 0, 0, 0, 0},
     {0, 0, 0, 0, 0}};

/*** END OF CONSTANTS ***/

int numTriloops; /* hairpin triloop penalties */
int numTetraloops; /* hairpin tetraloop penalties */
double atpS[5][5]; /* AT penalty */
double atpH[5][5]; /* AT penalty */
double *send5, *hend5; /* calc 5'  */
/* w/o init not constant anymore, cause for unimolecular and bimolecular foldings there are different values */
double dplx_init_H; /* initiation enthalpy; for duplex 200, for unimolecular structure 0 */
double dplx_init_S; /* initiation entropy; for duplex -5.7, for unimoleculat structure 0 */
double saltCorrection; /* value calculated by saltCorrectS, includes correction for monovalent and divalent cations */
double RC; /* universal gas constant multiplied w DNA conc - for melting temperature */
double SHleft; /* var that helps to find str w highest melting temperature */
int bestI, bestJ; /* starting position of most stable str */
double* enthalpyDPT; /* matrix for values of enthalpy */
double* entropyDPT; /* matrix for values of entropy */
int i;
unsigned char *oligo1, *oligo2; /* inserted oligo sequenced */
unsigned char *numSeq1, *numSeq2; /* same as oligo1 and oligo2 but converted to numbers */
int len1, len2, len3; /* length of sequense 1 and 2 *//* 17.02.2009 int temponly;*/ /* print only temperature of the predicted structure */
double dangleEntropies3[5][5][5]; /* thermodynamic paramteres for 3' dangling ends */
double dangleEnthalpies3[5][5][5]; /* ther params for 3' dangling ends */
double dangleEntropies5[5][5][5];  /* ther params for 5' dangling ends */
double dangleEnthalpies5[5][5][5]; /* ther params for 5' dangling ends */
double stackEntropies[5][5][5][5]; /* ther params for perfect match pairs */
double stackEnthalpies[5][5][5][5]; /* ther params for perfect match pairs */
double stackint2Entropies[5][5][5][5]; /*ther params for perfect match and internal mm */
double stackint2Enthalpies[5][5][5][5]; /* ther params for perfect match and internal mm*/
double interiorLoopEntropies[30]; /* interior loop params according to length of the loop */
double bulgeLoopEntropies[30]; /* bulge loop params according to length of the loop */
double hairpinLoopEntropies[30]; /* hairpin loop params accordint to length of the loop */
double interiorLoopEnthalpies[30]; /* same as interiorLoopEntropies but values of entropy */
double bulgeLoopEnthalpies[30]; /* same as bulgeLoopEntropies but values of entropy */
double hairpinLoopEnthalpies[30]; /* same as hairpinLoopEntropies but values of entropy */
double tstackEntropies[5][5][5][5]; /* ther params for terminal mismatches */
double tstackEnthalpies[5][5][5][5]; /* ther params for terminal mismatches */
double tstack2Entropies[5][5][5][5]; /* ther params for internal terminal mismatches */
double tstack2Enthalpies[5][5][5][5]; /* ther params for internal terminal mismatches */
struct triloop* triloopEntropies; /* ther penalties for given triloop seq-s */
struct triloop* triloopEnthalpies; /* ther penalties for given triloop seq-s */
struct tetraloop* tetraloopEntropies; /* ther penalties for given tetraloop seq-s */
struct tetraloop* tetraloopEnthalpies; /* ther penalties for given tetraloop seq-s */

static void fail_action(int a, thal_results *o)  {
   if (a) {
      fprintf(stderr, "\n%s\n", o->msg);
      exit(-1);
   }
   o->temp=THAL_ERROR_SCORE;
}

/* central method: execute all sub-methods for calculating secondary structure for dimer or for monomer */
void thal(const unsigned char *oligo_f, const unsigned char *oligo_r, const thal_args *a, thal_results *o) {
   double* SH;
   int i, j;
   int len_f, len_r;
   double T1;
   int k;
   int *bp;
   unsigned char *oligo2_rev;
   double mh, ms;
   if(a->debug == 0) {
#ifdef DEBUG
#undef DEBUG
#endif
   } else {
#define DEBUG
   }
   send5 = hend5 = NULL;
   enthalpyDPT = entropyDPT = NULL;
   numSeq1 = numSeq2 = NULL;
   oligo1 = oligo2 = NULL;
   o->msg = NULL;
   o->temp = THAL_ERROR_SCORE;
   CHECK_ERROR(NULL == oligo_f, "NULL first sequence");
   CHECK_ERROR(NULL == oligo_r, "NULL second sequence");
   len_f = length_unsig_char(oligo_f);
   len_r = length_unsig_char(oligo_r);
   CHECK_ERROR((len_f > THAL_MAX_ALIGN) && (len_r > THAL_MAX_ALIGN), 
	       "Sequences longer than THAL_MAX_ALIGN for thermodynamical alignment (nearest-neighbor approach)");
   CHECK_ERROR((len_f > THAL_MAX_SEQ), "Sequence 1 longer than THAL_MAX_SEQ and alignment is requested");
   CHECK_ERROR((len_r > THAL_MAX_SEQ), "Sequence 2 longer than THAL_MAX_SEQ and alignment is requested");

   CHECK_ERROR(NULL == a,  "NULL 'in' pointer");
   if (NULL == o) return; /* Leave it to the caller to crash */
   CHECK_ERROR(a->type != thal_any
	       && a->type != thal_end1
	       && a->type != thal_end2
	       && a->type != thal_hairpin,
	       "Illegal type");
   o->align_end_1 = -1;
   o->align_end_2 = -1;
   if ('\0' == oligo_f) { 
      o->msg = "Empty first sequence";
      o->temp = 0.0;
      return;
   }
   if ('\0' == oligo_r) {
      o->msg = "Empty second sequence";
      o->temp = 0.0;
      return;
   }
   if (0 == len_f) {
      o->temp = 0.0;
      return;
   }
   if (0 == len_r) {
      o->temp = 0.0;
      return;
   }
   if(a->type!=3) {	
      oligo1 = safe_malloc((len_f + 1) * sizeof(unsigned char),o,a->fail_stop);
      oligo2 = safe_malloc((len_r + 1) * sizeof(unsigned char),o,a->fail_stop);
      strcpy((char*)oligo1,(const char*)oligo_f);
      strcpy((char*)oligo2,(const char*)oligo_r);
   } else  {
      oligo1 = safe_malloc((len_r + 1) * sizeof(unsigned char),o,a->fail_stop);
      oligo2 = safe_malloc((len_f + 1) * sizeof(unsigned char),o,a->fail_stop);
      strcpy((char*)oligo1,(const char*)oligo_r);
      strcpy((char*)oligo2,(const char*)oligo_f);
   }
   /*** INIT values for unimolecular and bimolecular structures ***/
   if (a->type==4) { /* unimolecular folding */
      len2 = length_unsig_char(oligo2);
      len3 = len2 -1;
      dplx_init_H = 0.0;
      dplx_init_S = -0.00000000001;
      RC=0;
   } else if(a->type!=4) {
      /* hybridization of two oligos */
      dplx_init_H = 200;
      dplx_init_S = -5.7;
      if(symmetry_thermo(oligo1) && symmetry_thermo(oligo2)) {
	 RC = R  * log(a->dna_conc/1000000000.0);
      } else {
	 RC = R  * log(a->dna_conc/4000000000.0);
      }
      if(a->type!=3) {
	 oligo2_rev=safe_malloc((length_unsig_char(oligo_r) + 1) * sizeof(unsigned char),o,a->fail_stop);
	 strcpy((char*)oligo2_rev,(const char*)oligo_r);
      } else {
	 oligo2_rev=safe_malloc((length_unsig_char(oligo_f) + 1) * sizeof(unsigned char),o,a->fail_stop);
	 strcpy((char*)oligo2_rev,(const char*)oligo_f);
      }   
      reverse(oligo2_rev); /* REVERSE oligo2, so it goes to dpt 3'->5' direction */
      free(oligo2);
      oligo2=NULL;
      oligo2=&oligo2_rev[0];
   } else {
      o->msg = "Wrong alignment type!";
      o->temp = THAL_ERROR_SCORE;
#ifdef DEBUG
      fprintf(stderr, o->msg);
#endif
      exit(-1);
   }
   len1 = length_unsig_char(oligo1);
   len2 = length_unsig_char(oligo2);
   /* convert nucleotides to numbers */
   numSeq1 = safe_realloc(numSeq1, len1 + 2,o,a->fail_stop);
   numSeq2 = safe_realloc(numSeq2, len2 + 2,o,a->fail_stop); 
   
   /*** BEGIN: get thermodynamic values ***/
   
   getStack(stackEntropies, stackEnthalpies,o,a->fail_stop);
   /* verifyStackTable(stackEntropies, "entropy");
    verifyStackTable(stackEnthalpies, "enthalpy"); */ /* this is for code debugging */
   getStackint2(stackint2Entropies, stackint2Enthalpies,o,a->fail_stop);
   getDangle(dangleEntropies3, dangleEnthalpies3, dangleEntropies5, dangleEnthalpies5,o,a->fail_stop);
   getLoop(hairpinLoopEntropies, interiorLoopEntropies, bulgeLoopEntropies, hairpinLoopEnthalpies,
	    interiorLoopEnthalpies, bulgeLoopEnthalpies,o,a->fail_stop);
   getTstack(tstackEntropies, tstackEnthalpies,o,a->fail_stop);
   getTstack2(tstack2Entropies, tstack2Enthalpies,o,a->fail_stop);
   getTriloop(&triloopEntropies, &triloopEnthalpies, &numTriloops,o,a->fail_stop);
   getTetraloop(&tetraloopEntropies, &tetraloopEnthalpies, &numTetraloops,o,a->fail_stop);
   /* getting the AT-penalties */
   tableStartATS(AT_S,atpS);
   tableStartATH(AT_H,atpH);
   
   /*** END: get thermodynamic values ***/
   
   /*** Calc part of the salt correction ***/
   saltCorrection=saltCorrectS(a->mv,a->dv,a->dntp); /* salt correction for entropy, must be multiplied with N, which is
						   the total number of phosphates in the duplex divided by 2; 8bp dplx N=7 */
   
   if(a->type == 4){ /* monomer */
      /* terminal basepairs */
      send5 = safe_realloc(send5, (len1 + 1) * sizeof(double),o,a->fail_stop);
      hend5 = safe_realloc(hend5, (len1 + 1) * sizeof(double),o,a->fail_stop);
   }
   for(i = 0; i < len1; i++) oligo1[i] = toupper(oligo1[i]);
   for(i = 0; i < len2; i++) oligo2[i] = toupper(oligo2[i]);
   for(i = 1; i <= len1; ++i) numSeq1[i] = str2int(oligo1[i - 1]);
   for(i = 1; i <= len2; ++i) numSeq2[i] = str2int(oligo2[i - 1]);
   numSeq1[0] = numSeq1[len1 + 1] = numSeq2[0] = numSeq2[len2 + 1] = 4; /* mark as N-s */
   if (a->type==4) { /* calculate structure of monomer */
      enthalpyDPT = safe_recalloc(enthalpyDPT, len1, len2,o,a->fail_stop);
      entropyDPT = safe_recalloc(entropyDPT, len1, len2,o,a->fail_stop);
      initMatrix2();
      fillMatrix2(a->maxLoop,o,a->fail_stop);
      calc_terminal_bp(a->temp);
      mh = HEND5(len1);
      ms = SEND5(len1);
      o->align_end_1=mh;
      o->align_end_2=ms;
      bp = safe_calloc(len1, sizeof(int),o,a->fail_stop);
      for (k = 0; k < len1; ++k) bp[k] = 0;
      if(isFinite(mh)) {
	 tracebacku(bp, a->maxLoop,o,a->fail_stop);
	 /* traceback for unimolecular structure */
	 drawHairpin(bp, mh, ms, a->temponly,a->temp,o,a->fail_stop); /* if temponly=1 then return after printing basic therm data */
      } else if(a->temponly==0) {
	 fputs("No secondary structure could be calculated\n",stderr);
      }
      
      if(o->temp==-INFINITY && o->msg==NULL) o->temp=0.0; 
      free(bp);
      free(triloopEntropies);
      free(triloopEnthalpies);
      free(tetraloopEntropies);
      free(tetraloopEnthalpies);
      free(enthalpyDPT);
      free(entropyDPT);
      free(numSeq1);
      free(numSeq2);
      free(send5);
      free(hend5);
      free(oligo1);
      free(oligo2);
      return;
   } else if(a->type!=4) { /* Hybridization of two moleculs */
      len3 = len2;
      enthalpyDPT = safe_recalloc(enthalpyDPT, len1, len2,o,a->fail_stop); /* dyn. programming table for dS and dH */
      entropyDPT = safe_recalloc(entropyDPT, len1, len2,o,a->fail_stop); /* enthalpyDPT is 3D array represented as 1D array */
      initMatrix();
      fillMatrix(a->maxLoop,o,a->fail_stop);
      SHleft = -INFINITY;
      SH = safe_malloc(2 * sizeof(double),o,a->fail_stop);
      /* calculate terminal basepairs */
      bestI = bestJ = 0;
      if(a->type==1)
	for (i = 1; i <= len1; i++) {
	   for (j = 1; j <= len2; j++) {
	      RSH(i, j, SH);
	      SH[0] = SH[0]+SMALL_NON_ZERO; /* this adding is done for compiler, optimization -O2 vs -O0 */
	      SH[1] = SH[1]+SMALL_NON_ZERO;
	      T1 = ((EnthalpyDPT(i, j)+ SH[1] + dplx_init_H) / ((EntropyDPT(i, j)) + SH[0] +
								dplx_init_S + RC)) - ABSOLUTE_ZERO;
	      if (T1 > SHleft  && ((EntropyDPT(i, j) + SH[0])<0 && (SH[1] + EnthalpyDPT(i, j))<0)) {
		 SHleft = T1;
		 bestI = i;
		 bestJ = j;
	      }
	   }
	}
      int *ps1, *ps2;
      ps1 = safe_calloc(len1, sizeof(int),o,a->fail_stop);
      ps2 = safe_calloc(len2, sizeof(int),o,a->fail_stop);
      for (i = 0; i < len1; ++i)
	ps1[i] = 0;
      for (j = 0; j < len2; ++j)
	ps2[j] = 0;
      if(a->type == 2 || a->type == 3)	{
	 /* THAL_END1 */
	 bestI = bestJ = 0;
	 bestI = len1;
	 i = len1;
	 SHleft = -INFINITY;
	 for (j = 1; j <= len2; ++j) {
	    RSH(i, j, SH);
	    SH[0] = SH[0]+SMALL_NON_ZERO; /* this adding is done for compiler, optimization -O2 vs -O0, 
					   that compiler could understand that SH is changed in this cycle */
	    SH[1] = SH[1]+SMALL_NON_ZERO;
	    T1 = ((EnthalpyDPT(i, j)+ SH[1] + dplx_init_H) / ((EntropyDPT(i, j)) + SH[0] +
							      dplx_init_S + RC)) - ABSOLUTE_ZERO;
	    if (T1 > SHleft && ((SH[0] + EntropyDPT(i, j))<0 && (SH[1] + EnthalpyDPT(i, j))<0)) {
	       SHleft = T1;
	       bestJ = j;
	    }
	 }
      }
      if (!isFinite(SHleft)) bestI = bestJ = 1;
      double dH, dS;
      RSH(bestI, bestJ, SH);
      dH = EnthalpyDPT(bestI, bestJ)+ SH[1] + dplx_init_H;
      dS = (EntropyDPT(bestI, bestJ) + SH[0] + dplx_init_S);
      /* tracebacking */
      for (i = 0; i < len1; ++i)
	ps1[i] = 0;
      for (j = 0; j < len2; ++j)
	ps2[j] = 0;
      if(isFinite(EnthalpyDPT(bestI, bestJ))){
	 traceback(bestI, bestJ, RC, ps1, ps2, a->maxLoop,o,a->fail_stop);
	 drawDimer(ps1, ps2, SHleft, dH, dS, a->temponly,a->temp,o,a->fail_stop);
	 o->align_end_1=bestI;
	 o->align_end_2=bestJ;
      } else  {
	 o->temp = 0.0;
      }
      free(ps1);
      free(ps2);
      free(SH);
      free(oligo2_rev);
      free(triloopEntropies);
      free(triloopEnthalpies);
      free(tetraloopEntropies);
      free(tetraloopEnthalpies);
      free(enthalpyDPT);
      free(entropyDPT);
      free(numSeq1);
      free(numSeq2);
      free(oligo1);
      return;
   }
   return;
FAIL:
   fail_action(a->fail_stop, o);
} 
/*** END thal() ***/

/* Set default args */
void set_thal_default_args(thal_args *a) {
   memset(a, 0, sizeof(*a));
   a->debug = 0;
   a->type = 1; /* thal_alignment_type THAL_ANY */
   a->maxLoop = MAX_LOOP;
   a->mv = 50; /* mM */
   a->dv = 0.0; /* mM */
   a->dntp = 0.8; /* mM */
   a->dna_conc = 50; /* nM */
   a->temp = 310.15; /* Kelvin */
   a->temponly = 1; /* return only melting temperature of predicted structure */
   a->dimer = 1; /* by default dimer structure is calculated */
   a->fail_stop = THAL_EXIT_ON_ERROR;
}

unsigned char str2int(char c) {
   switch (c) {
    case 'A': case '0':
      return 0;
    case 'C': case '1':
      return 1;
    case 'G': case '2':
      return 2;
    case 'T': case '3':
      return 3;
   }
   return 4;
}

int seqcmp(unsigned char* seq1, unsigned char* seq2, int length) {
   int i;
   for (i = 0; i < length; ++i)
     if (seq1[i] < seq2[i])
       return -1;
   else if (seq1[i] > seq2[i])
     return 1;
   return 0;
}

/* memory stuff */

double* safe_recalloc(double* ptr, int m, int n,thal_results* o, int a) {
   return safe_realloc(ptr, m * n * sizeof(double),o,a);
}

void* safe_calloc(size_t m, size_t n,thal_results *o, int a) {
   void* ptr;
   if (!(ptr = calloc(m, n))) {
#ifdef DEBUG
      fputs("Error in calloc()\n", stderr);
#endif
      THAL_OOM_ERROR;
      exit(EXIT_FAILURE);
   }
   return ptr;
   
FAIL:
   fail_action(a, o);
   exit(EXIT_FAILURE);
}

void* safe_malloc(size_t n,thal_results *o, int a) {
   void* ptr;
   if (!(ptr = malloc(n))) {
#ifdef DEBUG
      fputs("Error in malloc()\n", stderr);
#endif
      THAL_OOM_ERROR;
      exit(EXIT_FAILURE);
   }
   return ptr;
   
FAIL:
   fail_action(a, o);
   exit(EXIT_FAILURE);
}

void* safe_realloc(void* ptr, size_t n,thal_results *o, int a) {
   ptr = realloc(ptr, n);
   if (ptr==NULL) {
#ifdef DEBUG
      fputs("Error in realloc()\n", stderr);
#endif
      THAL_OOM_ERROR;
      exit(EXIT_FAILURE);
   }
   return ptr;

FAIL:
   fail_action(a, o);
   exit(EXIT_FAILURE);
}

int smallest(int a, int b, int c) {
   if (a <= b && a <= c)
     return a;
   if (b <= c)
     return b;
   return c;
}

int identic(unsigned char* a, unsigned char* b, int len) {
   int i;
   for (i = 1; i <= len; ++i)
     if (a[i] != b[i])
       return 0;
   return 1;
}

int max5(double a, double b, double c, double d, double e){
   if(a > b && a > c && a > d && a > e) return 1;
   else if(b > c && b > d && b > e) return 2;
   else if(c > d && c > e) return 3;
   else if(d > e) return 4;
   else return 5;
}

void push(struct tracer** stack, int i, int j, int mtrx,thal_results* o,int a)
{
   struct tracer* new_top;  
   new_top = safe_malloc(sizeof(struct tracer),o,a);
   new_top->i = i;
   new_top->j = j;
   new_top->mtrx = mtrx;
   new_top->next = *stack;
   *stack = new_top;
}

void reverse(unsigned char *s) {
   int i,j;
   char c;
   for (i = 0, j = length_unsig_char(s)-1; i < j; i++, j--) {
      c = s[i];
      s[i] = s[j];
      s[j] = c;
   }
}

FILE* openParamFile(char* fname,thal_results* o,int a) {
   FILE* file;
   char* paramdir;
   file = fopen(fname, "rt");
   if(!file) {
      paramdir = safe_malloc(strlen(PARAMFILES) + strlen(fname) + 1,o,a);
      strcpy(paramdir,PARAMFILES);
      strcat(paramdir,fname);
      if(!(file = fopen(paramdir, "rt"))) {
#ifdef DEBUG
	 perror(paramdir);
#endif
	 free(paramdir);
	 exit(EXIT_FAILURE);
      }
      free(paramdir);
   }
   return file;
}

double saltCorrectS (double mv, double dv, double dntp) {
   if(dv<=0) dntp=dv;
   return 0.368*((log((mv+120*(sqrt(dv-dntp)))/1000)));
}

	
void getStack(double stackEntropies[5][5][5][5], double stackEnthalpies[5][5][5][5],thal_results* o,int a) {
   int i, j, ii, jj;
   FILE *sFile, *hFile;
   sFile = openParamFile("stack.ds",o,a);
   hFile = openParamFile("stack.dh",o,a);
   for (i = 0; i < 5; ++i) {
      for (ii = 0; ii < 5; ++ii) {
	 for (j = 0; j < 5; ++j) {
	    for (jj = 0; jj < 5; ++jj) {
	       if (i == 4 || j == 4 || ii == 4 || jj == 4) { 
		  stackEntropies[i][ii][j][jj] = -1.0;
		  stackEnthalpies[i][ii][j][jj] = INFINITY;
	       } else {
		  fscanf(sFile, "%lf", &stackEntropies[i][ii][j][jj]);
		  fscanf(hFile, "%lf", &stackEnthalpies[i][ii][j][jj]);
		  if (!isFinite(stackEntropies[i][ii][j][jj]) || !isFinite(stackEnthalpies[i][ii][j][jj])) {
		     stackEntropies[i][ii][j][jj] = -1.0;
		     stackEnthalpies[i][ii][j][jj] = INFINITY;
		  }
	       }
	    }
	 }
      }
   }
   fclose(sFile);
   fclose(hFile);
}

void getStackint2(double stackint2Entropies[5][5][5][5], double stackint2Enthalpies[5][5][5][5],thal_results* o,int a) {
   int i, j, ii, jj;
   FILE *sFile, *hFile;
   sFile = openParamFile("stackmm.ds",o,a);
   hFile = openParamFile("stackmm.dh",o,a);
   for (i = 0; i < 5; ++i) {
      for (ii = 0; ii < 5; ++ii) {
	 for (j = 0; j < 5; ++j) {
	    for (jj = 0; jj < 5; ++jj) {
	       if (i == 4 || j == 4 || ii == 4 || jj == 4) {
		  stackint2Entropies[i][ii][j][jj] = -1.0;
		  stackint2Enthalpies[i][ii][j][jj] = INFINITY;
	       } else {
		  fscanf(sFile, "%lf", &stackint2Entropies[i][ii][j][jj]);
		  fscanf(hFile, "%lf", &stackint2Enthalpies[i][ii][j][jj]);
		  if (!isFinite(stackint2Entropies[i][ii][j][jj]) || !isFinite(stackint2Enthalpies[i][ii][j][jj])) {
		     stackint2Entropies[i][ii][j][jj] = -1.0;
		     stackint2Enthalpies[i][ii][j][jj] = INFINITY;
		  }
	       }
	    }
	 }	 
      }	
   }  
   fclose(sFile);
   fclose(hFile);
}


void verifyStackTable(double stack[5][5][5][5], char* type) {
   
   int i, j, ii, jj;
   for (i = 0; i < 4; ++i)
     for (j = 0; j < 4; ++j)
       for (ii = 0; ii < 4; ++ii)
	 for (jj = 0; jj < 4; ++jj)
	   if (stack[i][j][ii][jj] != stack[jj][ii][j][i])
#ifdef DEBUG
	     fprintf(stderr, "Warning: symmetrical stacks _are_ _not_ equal: %c-%c/%c-%c stack %s is %g; %c-%c/%c-%c stack %s is %g\n", 
#endif
		     BASES[i], BASES[j], BASES[ii], BASES[jj], type, stack[i][j][ii][jj], BASES[jj], 
		     BASES[ii], BASES[j], BASES[i], type, stack[jj][ii][j][i]);
}


void getDangle(double dangleEntropies3[5][5][5], double dangleEnthalpies3[5][5][5], double dangleEntropies5[5][5][5], 
		double dangleEnthalpies5[5][5][5],thal_results* o,int a) {
   int i, j, k;
   FILE *sFile, *hFile;
   sFile = openParamFile("dangle.ds",o,a);
   hFile = openParamFile("dangle.dh",o,a);
   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       for (k = 0; k < 5; ++k) {	    
	  if (i == 4 || j == 4) {
	     dangleEntropies3[i][k][j] = -1.0;
	     dangleEnthalpies3[i][k][j] = INFINITY;
	  } else if (k == 4) {
	     dangleEntropies3[i][k][j] = -1.0;
	     dangleEnthalpies3[i][k][j] = INFINITY;
	  } else {
	     fscanf(sFile, "%lf", &dangleEntropies3[i][k][j]);
	     fscanf(hFile, "%lf", &dangleEnthalpies3[i][k][j]);
	     if(!isFinite(dangleEntropies3[i][k][j]) || !isFinite(dangleEnthalpies3[i][k][j])) {
		dangleEntropies3[i][k][j] = -1.0;
		dangleEnthalpies3[i][k][j] = INFINITY;	     }

	  }
       }
 
   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       for (k = 0; k < 5; ++k) {
	  if (i == 4 || j == 4) {
	     dangleEntropies5[i][j][k] = -1.0; 
	     dangleEnthalpies5[i][j][k] = INFINITY;
	       } else if (k == 4) {
		  dangleEntropies5[i][j][k] = -1.0;
		  dangleEnthalpies5[i][j][k] = INFINITY;
	       } else {
		  fscanf(sFile, "%lf", &dangleEntropies5[i][j][k]);
		  fscanf(hFile, "%lf", &dangleEnthalpies5[i][j][k]);
	     if(!isFinite(dangleEntropies5[i][j][k]) || !isFinite(dangleEnthalpies5[i][j][k])) {
		dangleEntropies5[i][j][k] = -1.0;
		dangleEnthalpies5[i][j][k] = INFINITY;
	     }
	  }
       }
   fclose(sFile);
   fclose(hFile);
}

void getLoop(double hairpinLoopEntropies[30], double interiorLoopEntropies[30], double bulgeLoopEntropies[30], 
	      double hairpinLoopEnthalpies[30], double interiorLoopEnthalpies[30], double bulgeLoopEnthalpies[30],thal_results* o,int a) {
   int k;
   FILE *sFile, *hFile;
   sFile = openParamFile("loops.ds",o,a);
   hFile = openParamFile("loops.dh",o,a);
   for (k = 0; k < 30; ++k) {
      fscanf(sFile, "%*f%lf%lf%lf", &interiorLoopEntropies[k], &bulgeLoopEntropies[k], &hairpinLoopEntropies[k]);
      fscanf(hFile, "%*f%lf%lf%lf", &interiorLoopEnthalpies[k], &bulgeLoopEnthalpies[k], &hairpinLoopEnthalpies[k]);
   }
   fclose(sFile);
   fclose(hFile);
}

void getTstack(double tstackEntropies[5][5][5][5], double tstackEnthalpies[5][5][5][5],thal_results* o,int a) {
   int i1, j1, i2, j2;
   FILE *sFile, *hFile;
   sFile = openParamFile("tstack_tm_inf.ds",o,a);
   hFile = openParamFile("tstack.dh",o,a);
   for (i1 = 0; i1 < 5; ++i1)
     for (i2 = 0; i2 < 5; ++i2)
       for (j1 = 0; j1 < 5; ++j1)
	 for (j2 = 0; j2 < 5; ++j2)
	   if (i1 == 4 || j1 == 4) {
	      tstackEnthalpies[i1][i2][j1][j2] = INFINITY;
	      tstackEntropies[i1][i2][j1][j2] = -1.0;
	   } else if (i2 == 4 || j2 == 4) {
	      tstackEntropies[i1][i2][j1][j2] = 0.00000000001;
	      tstackEnthalpies[i1][i2][j1][j2] = 0.0;
	   } else {
	      fscanf(sFile, "%lf", &tstackEntropies[i1][i2][j1][j2]);
	      fscanf(hFile, "%lf", &tstackEnthalpies[i1][i2][j1][j2]);
	      if (!isFinite(tstackEntropies[i1][i2][j1][j2]) || !isFinite(tstackEnthalpies[i1][i2][j1][j2])) {	   
		 tstackEntropies[i1][i2][j1][j2] = -1.0;
		 tstackEnthalpies[i1][i2][j1][j2] = INFINITY;
	      }
	   }
   fclose(sFile);
   fclose(hFile);
}

void getTstack2(double tstack2Entropies[5][5][5][5], double tstack2Enthalpies[5][5][5][5],thal_results* o,int a){
   
   int i1, j1, i2, j2;
   FILE *sFile, *hFile;
   sFile = openParamFile("tstack2.ds",o,a);
   hFile = openParamFile("tstack2.dh",o,a);
   for (i1 = 0; i1 < 5; ++i1)
     for (i2 = 0; i2 < 5; ++i2)
       for (j1 = 0; j1 < 5; ++j1)
	 for (j2 = 0; j2 < 5; ++j2)
	   if (i1 == 4 || j1 == 4)  {
	      tstack2Enthalpies[i1][i2][j1][j2] = INFINITY; 
	      tstack2Entropies[i1][i2][j1][j2] = -1.0;
	   } else if (i2 == 4 || j2 == 4) {
	      tstack2Entropies[i1][i2][j1][j2] = 0.00000000001;
	      tstack2Enthalpies[i1][i2][j1][j2] = 0.0;
	   } else {
	      
	      fscanf(sFile, "%lf", &tstack2Entropies[i1][i2][j1][j2]);
	      fscanf(hFile, "%lf", &tstack2Enthalpies[i1][i2][j1][j2]);
	      if (!isFinite(tstack2Entropies[i1][i2][j1][j2]) || !isFinite(tstack2Enthalpies[i1][i2][j1][j2])) {
		 tstack2Entropies[i1][i2][j1][j2] = -1.0;
		 tstack2Enthalpies[i1][i2][j1][j2] = INFINITY;
	      }
	   }   
   fclose(sFile);
   fclose(hFile);
}

void getTriloop(struct triloop** triloopEntropies, struct triloop** triloopEnthalpies, int* num,thal_results* o,int a) {
   FILE *sFile, *hFile;
   int i, size;
   double value;
   sFile = openParamFile("triloop.ds",o,a);
   *num = 0;
   size = 16;
   *triloopEntropies = calloc(16, sizeof(struct triloop));
   while (fscanf(sFile, "%5s %lg", (*triloopEntropies)[*num].loop, &value) == 2) {
      for (i = 0; i < 5; ++i)
	(*triloopEntropies)[*num].loop[i] = str2int((*triloopEntropies)[*num].loop[i]);
      (*triloopEntropies)[*num].value = value;
      ++*num;
      if (*num == size)	{
	 size *= 2;
	 *triloopEntropies = realloc(*triloopEntropies, size * sizeof(struct triloop));
      }
   }
   *triloopEntropies = realloc(*triloopEntropies, *num * sizeof(struct triloop));
      
   fclose(sFile);
   
   hFile = openParamFile("triloop.dh",o,a);
   *num = 0;
   size = 16;
   *triloopEnthalpies = calloc(16, sizeof(struct triloop));
   	
   while (fscanf(hFile, "%5s %lg", (*triloopEnthalpies)[*num].loop, &value) == 2) {
      
      for (i = 0; i < 5; ++i)
	(*triloopEnthalpies)[*num].loop[i] = str2int((*triloopEnthalpies)[*num].loop[i]);
      (*triloopEnthalpies)[*num].value = value;
      ++*num;
      if (*num == size) {
	 size *= 2;
	 *triloopEnthalpies = realloc(*triloopEnthalpies, size * sizeof(struct triloop));
      }      
   }  
   *triloopEnthalpies = realloc(*triloopEnthalpies, *num * sizeof(struct triloop));
   
   fclose(hFile);
}

void getTetraloop(struct tetraloop** tetraloopEntropies, struct tetraloop** tetraloopEnthalpies, int* num,thal_results* o,int a) {
   
   FILE *sFile, *hFile;
   int i, size;
   double value;
   sFile = openParamFile("tetraloop.ds",o,a);
   *num = 0;
   size = 16;
   *tetraloopEntropies = calloc(16, sizeof(struct tetraloop));
   while (fscanf(sFile, "%6s %lg", (*tetraloopEntropies)[*num].loop, &value) == 2) {
      for (i = 0; i < 6; ++i)
	(*tetraloopEntropies)[*num].loop[i] = str2int((*tetraloopEntropies)[*num].loop[i]);
      (*tetraloopEntropies)[*num].value = value;
      ++*num;
      if (*num == size) {
	 size *= 2;
	 *tetraloopEntropies = realloc(*tetraloopEntropies, size * sizeof(struct tetraloop));
      }
   }
   *tetraloopEntropies = realloc(*tetraloopEntropies, *num * sizeof(struct tetraloop));
   fclose(sFile);
   
   hFile = openParamFile("tetraloop.dh",o,a);
   *num = 0;
   size = 16;
   *tetraloopEnthalpies = calloc(16, sizeof(struct tetraloop));
   while (fscanf(hFile, "%6s %lg", (*tetraloopEnthalpies)[*num].loop, &value) == 2) {
      for (i = 0; i < 6; ++i)
	(*tetraloopEnthalpies)[*num].loop[i] = str2int((*tetraloopEnthalpies)[*num].loop[i]);
      (*tetraloopEnthalpies)[*num].value = value;
      ++*num;
      if (*num == size) {
	 size *= 2;
	 *tetraloopEnthalpies = realloc(*tetraloopEnthalpies, size * sizeof(struct tetraloop));
      }
   }
   *tetraloopEnthalpies = realloc(*tetraloopEnthalpies, *num * sizeof(struct tetraloop));
   	
   fclose(hFile);
}

void tableStartATS(double atp_value, double atpS[5][5]){
   
   int i, j;
   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       atpS[i][j] = 0.00000000001;
     atpS[0][3] = atpS[3][0] = atp_value;
}


void tableStartATH(double atp_value, double atpH[5][5]) {
   
   int i, j;
   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       atpH[i][j] = 0.0;
   
     atpH[0][3] = atpH[3][0] = atp_value;
}

int comp3loop(const void* loop1, const void* loop2)
{
   
     int i;
     const unsigned char* h1 = loop1;
     const struct triloop *h2 = loop2;
   
     for (i = 0; i < 5; ++i)
         if (h1[i] < h2->loop[i])
             return -1;
       else if (h1[i] > h2->loop[i])
           return 1;
   
     return 0;
}

int comp4loop(const void* loop1, const void* loop2) {
   int i;
   const unsigned char* h1 = loop1;
   const struct tetraloop *h2 = loop2;
   
   for (i = 0; i < 6; ++i)
     if (h1[i] < h2->loop[i])
       return -1;
   else if (h1[i] > h2->loop[i])
     return 1;
   
   return 0;
}


void initMatrix() {   
   int i, j;
   for (i = 1; i <= len1; ++i) {
      for (j = 1; j <= len2; ++j) {
	 if (bpIndx(numSeq1[i], numSeq2[j]) == 0)  {
	    EnthalpyDPT(i, j) = INFINITY;
	    EntropyDPT(i, j) = -1.0;
	 } else {
	    EnthalpyDPT(i, j) = 0.0;
	    EntropyDPT(i, j) = MinEntropy;
	 }
      }
   }
}

void initMatrix2() {
   int i, j;   
   for (i = 1; i <= len1; ++i)
     for (j = i; j <= len2; ++j)
       if (j - i < MIN_HRPN_LOOP + 1 || (bpIndx(numSeq1[i], numSeq1[j]) == 0)) {
	  EnthalpyDPT(i, j) = INFINITY;
	  EntropyDPT(i, j) = -1.0;
       } else {
	  EnthalpyDPT(i, j) = 0.0;
	  EntropyDPT(i, j) = MinEntropy;
       }  
   
}

void fillMatrix(int maxLoop,thal_results *o, int a) {
   int d, i, j, ii, jj;
   double* SH;
   SH = safe_malloc(2 * sizeof(double),o,a);
   for (i = 1; i <= len1; ++i) {
      for (j = 1; j <= len2; ++j) {
	 if(isFinite(EnthalpyDPT(i, j))) { /* if finite */
	    SH[0] = -1.0;
	    SH[1] = INFINITY;
	    LSH(i,j,SH);
	    if(isFinite(SH[1])) {
	       EntropyDPT(i,j) = SH[0];
	       EnthalpyDPT(i,j) = SH[1]; 
	    }
	    if (i > 1 && j > 1) { 
	       maxTM(i, j); /* stack: sets EntropyDPT(i, j) and EnthalpyDPT(i, j) */	
	       for(d = 3; d <= maxLoop + 2; d++) { /* max=30, length over 30 is not allowed */
		  ii = i - 1; 
		  jj = - ii - d + (j + i); 
		  if (jj < 1) {
		     ii -= abs(jj-1);
		     jj = 1; 
		  }
		  for (; ii > 0 && jj < j; --ii, ++jj) {
		     if (isFinite(EnthalpyDPT(ii, jj))) {
			SH[0] = -1.0;
			SH[1] = INFINITY;
			calc_bulge_internal(ii, jj, i, j, SH,0,maxLoop); 
			if(SH[0] < MinEntropyCutoff) {
			   /* to not give dH any value if dS is unreasonable */
			   SH[0] = MinEntropy;
			   SH[1] = 0.0;
			}			
			if(isFinite(SH[1])) {
			   EnthalpyDPT(i, j) = SH[1];
			   EntropyDPT(i, j) = SH[0]; 
			}
		     }
		  }
	       }
	    } /* if */
	 }
      } /* for */  
   } /* for */
   free(SH);
}

void fillMatrix2(int maxLoop,thal_results* o, int a) {
   int i, j;
   double* SH;
   SH = safe_malloc(2 * sizeof(double),o,a);
   for (j = 2; j <= len2; ++j)
     for (i = j - MIN_HRPN_LOOP - 1; i >= 1; --i) {
	if (isFinite(EnthalpyDPT(i, j))) {
	   SH[0] = -1.0;
	   SH[1] = INFINITY;
	   maxTM2(i,j); /* calculate stack */
	   CBI(i, j, SH, 0,maxLoop); /* calculate Bulge and Internal loop and stack */
	   SH[0] = -1.0; 
	   SH[1] = INFINITY;
	   calc_hairpin(i, j, SH, 0);
	   if(isFinite(SH[1])) {
	      if(SH[0] < MinEntropyCutoff){ /* to not give dH any value if dS is unreasonable */
		 SH[0] = MinEntropy;
		 SH[1] = 0.0;
	      }	      
	      EntropyDPT(i,j) = SH[0];
	      EnthalpyDPT(i,j) = SH[1];
	   }
	}
     }
   free(SH);
}


void maxTM(int i, int j) {
   double T0, T1;
   double S0, S1;
   double H0, H1;
   T0 = T1 = -INFINITY;
   S0 = EntropyDPT(i, j);
   H0 = EnthalpyDPT(i, j);
   T0 = (H0 + dplx_init_H) /(S0 + dplx_init_S + RC); /* at current position */
   if(isFinite(EnthalpyDPT(i - 1, j - 1)) && isFinite(Hs(i - 1, j - 1, 1))) {
      S1 = (EntropyDPT(i - 1, j - 1) + Ss(i - 1, j - 1, 1));
      H1 = (EnthalpyDPT(i - 1, j - 1) + Hs(i - 1, j - 1, 1));
   } else {
      S1 = -1.0;
      H1 = INFINITY;
   }
   T1 = (H1 + dplx_init_H) /(S1 + dplx_init_S + RC);
   if(S1 < MinEntropyCutoff) {
      /* to not give dH any value if dS is unreasonable */
      S1 = MinEntropy;
      H1 = 0.0;
   }   
   if(S0 < MinEntropyCutoff) {
      /* to not give dH any value if dS is unreasonable */
      S0 = MinEntropy;
      H0 = 0.0;
   }
   if(max2(T0, T1)==2 || (S0>0 && H0>0)){ /* T1 on suurem */
      EntropyDPT(i, j) = S1;
      EnthalpyDPT(i, j) = H1;
   } else if(max2(T0, T1)==1) {
      EntropyDPT(i, j) = S0;
      EnthalpyDPT(i, j) = H0;
   }
}

void maxTM2(int i, int j) {
   double T0, T1;
   double S0, S1;
   double H0, H1;
   T0 = T1 = -INFINITY;
   S0 = EntropyDPT(i, j);
   H0 = EnthalpyDPT(i, j);
   T0 = (H0 + dplx_init_H) /(S0 + dplx_init_S + RC);
   if(isFinite(EnthalpyDPT(i, j))) {	
      S1 = (EntropyDPT(i + 1, j - 1) + Ss(i, j, 2));
      H1 = (EnthalpyDPT(i + 1, j - 1) + Hs(i, j, 2));
   } else {
      S1 = -1.0;
      H1 = INFINITY;
   }  
   T1 = (H1 + dplx_init_H) /(S1 + dplx_init_S + RC);
   if(S1 < MinEntropyCutoff) {
      S1 = MinEntropy;
      H1 = 0.0;
   }
   if(S0 < MinEntropyCutoff) {
      S0 = MinEntropy;
      H0 = 0.0;
   }
   
   if(max2(T0, T1)==2) { 
      EntropyDPT(i, j) = S1;
      EnthalpyDPT(i, j) = H1;
   } else if(max2(T0, T1)==1) {
      EntropyDPT(i, j) = S0;
      EnthalpyDPT(i, j) = H0;
   }  
}


void LSH(int i, int j, double* EntropyEnthalpy) { 
   double S1, H1, T1;
   double S2, H2, T2;
   S1 = S2 = -1.0;
   H1 = H2 = -INFINITY;
   T1 = T2 = -INFINITY;
   if (bpIndx(numSeq1[i], numSeq2[j]) == 0) {
      EntropyDPT(i, j) = -1.0;
      EnthalpyDPT(i, j) = INFINITY;
      return;
   }
   S1 = atPenaltyS(numSeq1[i], numSeq2[j]) + tstack2Entropies[numSeq2[j]][numSeq2[j-1]][numSeq1[i]][numSeq1[i-1]];
   H1 = atPenaltyH(numSeq1[i], numSeq2[j]) + tstack2Enthalpies[numSeq2[j]][numSeq2[j-1]][numSeq1[i]][numSeq1[i-1]];
   if(!isFinite(H1)) {
      H1 = INFINITY;
      S1 = -1.0;
   }   
   /** If there is two dangling ends at the same end of duplex **/
   if(isFinite(dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]]) && isFinite(dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]] +
	dangleEntropies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]] +
	dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      if(!isFinite(H2)) {
	 H2 = INFINITY;
	 S2 = -1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1)) { 
	 T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
	 if(T1<T2) {
	    S1 = S2;
	    H1 = H2;
	    T1 = T2;
	 }
      } else {
	 S1 = S2;
	 H1 = H2;
	 T1 = T2;
      }     
   } else if (isFinite(dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]];
      if(!isFinite(H2)) {
	 H2 = INFINITY;
	 S2 = -1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1)) {
	 T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
	 if(T1<T2) {
	    S1 = S2;
	    H1 = H2;
	    T1 = T2;
	 }
      } else {
	 S1 = S2;
	 H1 = H2;
	 T1 = T2; 
      }
   } else if (isFinite(dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      if(!isFinite(H2)) {
	 H2 = INFINITY;
	 S2 = -1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1)) {
	 T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
	 if(T1 < T2) {
	    S1 = S2;
	    H1 = H2;
	    T1 = T2; 
	 }
      } else {
	 S1 = S2;
	 H1 = H2;
	 T1 = T2; 
      }
   }
   S2 = atPenaltyS(numSeq1[i], numSeq2[j]);
   H2 = atPenaltyH(numSeq1[i], numSeq2[j]);   
   T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
   if(isFinite(H1)) {  
      if(T1 < T2) {
	 EntropyEnthalpy[0] = S2;
	 EntropyEnthalpy[1] = H2;
      } else {
	 EntropyEnthalpy[0] = S1;
	 EntropyEnthalpy[1] = H1;
      }
   } else {
      EntropyEnthalpy[0] = S2;
      EntropyEnthalpy[1] = H2; 
   }    
   return;
}

void RSH(int i, int j, double* EntropyEnthalpy) { 
   double S1, S2;
   double H1, H2;
   double T1, T2;
   S1 = S2 = -1.0;
   H1 = H2 = INFINITY;
   T1 = T2 = -INFINITY;
   if (bpIndx(numSeq1[i], numSeq2[j]) == 0) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = INFINITY;
      return;
   }
   S1 = atPenaltyS(numSeq1[i], numSeq2[j]) + tstack2Entropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   H1 = atPenaltyH(numSeq1[i], numSeq2[j]) + tstack2Enthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   if(!isFinite(H1)) {
      H1 = INFINITY;
      S1 = -1.0;
   }
   if(isFinite(dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]]) && isFinite(dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +
	dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +
	dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
      if(!isFinite(H2)) {
	 H2 = INFINITY;
	 S2 = -1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1)) {
	 T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
	 if(T1 < T2) {
	    S1 = S2;
	    H1 = H2;
	    T1 = T2;
	 }
      } else {
	 S1 = S2;
	 H1 = H2;
	 T1 = T2;
      }
   } 
   
   if(isFinite(dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]];
      if(!isFinite(H2)) {
	 H2 = INFINITY;
	 S2 = -1.0;
      }      
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1)) {
	 T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);	
	 if(T1 < T2) {
	    S1 = S2;
	    H1 = H2;
	    T1 = T2;
	 }
      } else {
	 S1 = S2;
	 H1 = H2;
	 T1 = T2;
      }
   }
   
   if(isFinite(dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
      if(!isFinite(H2)) {
	 H2 = INFINITY;
	 S2 = -1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1)) {
	 T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
	 if(T1 < T2) {
	    S1 = S2;
	    H1 = H2;
	    T1 = T2;
	 }
      } else { 
	 S1 = S2;
	 H1 = H2;
	 T1 = T2;
      }
   }
   S2 = atPenaltyS(numSeq1[i], numSeq2[j]);
   H2 = atPenaltyH(numSeq1[i], numSeq2[j]);
   T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
   if(isFinite(H1)) {
      if(T1 < T2) {
	 EntropyEnthalpy[0] = S2;
	 EntropyEnthalpy[1] = H2;
      } else {
	 EntropyEnthalpy[0] = S1;
	 EntropyEnthalpy[1] = H1;
      }
   } else {
      EntropyEnthalpy[0] = S2;
      EntropyEnthalpy[1] = H2;
   }
   return;
}

double Ss(int i, int j, int k) {   
   if(k==2) {
      if (i >= j)
	return -1.0;
      if (i == len1 || j == len2 + 1)
	return -1.0;
      
      if (i > len1)
	i -= len1;
      if (j > len2)
	j -= len2;
      return stackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]];
   } else {
      return stackEntropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   }
}


double Hs(int i, int j, int k) {
   if(k==2) {
      if (i >= j)
	return INFINITY;
      if (i == len1 || j == len2 + 1)
	return INFINITY;
      
      if (i > len1)
	i -= len1;
      if (j > len2)
	j -= len2;
      if(isFinite(stackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]])) {
	 return stackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]];
      } else {
	 return INFINITY;
      }
   } else {
      return stackEnthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   }
}

void CBI(int i, int j, double* EntropyEnthalpy, int traceback, int maxLoop) {
   int d, ii, jj;
   for (d = j - i - 3; d >= MIN_HRPN_LOOP + 1 && d >= j - i - 2 - maxLoop; --d)
     for (ii = i + 1; ii < j - d && ii <= len1; ++ii) {
	jj = d + ii;
	if(traceback==0) {
	   EntropyEnthalpy[0] = -1.0;
	   EntropyEnthalpy[1] = INFINITY;
	}
	if (isFinite(EnthalpyDPT(ii, jj))) {
	   calc_bulge_internal2(i, j, ii, jj, EntropyEnthalpy, traceback,maxLoop);
	   if(isFinite(EntropyEnthalpy[1])) {
	      if(EntropyEnthalpy[0] < MinEntropyCutoff) {
		 EntropyEnthalpy[0] = MinEntropy;
		 EntropyEnthalpy[1] = 0.0;
	      }
	      if(traceback==0) {
		 EnthalpyDPT(i, j) = EntropyEnthalpy[1];
		 EntropyDPT(i, j) = EntropyEnthalpy[0];
	      }
	   }
	}
     }
   return;
}

void calc_hairpin(int i, int j, double* EntropyEnthalpy, int traceback) {
   int loopSize = j - i - 1;
   double T1, T2;
   T1 = T2 = -INFINITY;
   if(loopSize < MIN_HRPN_LOOP) { 
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = INFINITY;
      return;
   }
   if (i <= len1 && len2 < j) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = INFINITY;
      return;
   } else if (i > len2) {
      i -= len1;
      j -= len2;
   }
   if(loopSize <= 30) {	
      EntropyEnthalpy[1] = hairpinLoopEnthalpies[loopSize - 1];
      EntropyEnthalpy[0] = hairpinLoopEntropies[loopSize - 1];
   } else {
      EntropyEnthalpy[1] = hairpinLoopEnthalpies[29];
      EntropyEnthalpy[0] = hairpinLoopEntropies[29];
   }
      
   if (loopSize > 3) { /* for loops 4 bp and more in length, terminal mm are accounted */
      EntropyEnthalpy[1] += tstack2Enthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq1[j]][numSeq1[j - 1]];
      EntropyEnthalpy[0] += tstack2Entropies[numSeq1[i]][numSeq1[i + 1]][numSeq1[j]][numSeq1[j - 1]];
   } else if(loopSize == 3){ /* for loops 3 bp in length at-penalty is considered */
      EntropyEnthalpy[1] += atPenaltyH(numSeq1[i], numSeq1[j]); 
      EntropyEnthalpy[0] += atPenaltyS(numSeq1[i], numSeq1[j]);
   }
   
   if (loopSize == 3) {	 /* closing AT-penalty (+), triloop bonus, hairpin of 3 (+) */
      struct triloop* loop;
      if (numTriloops) {
	 if ((loop = bsearch(numSeq1 + i, triloopEnthalpies, numTriloops, sizeof(struct triloop), comp3loop)))
	   EntropyEnthalpy[1] += loop->value;
	 if ((loop = bsearch(numSeq1 + i, triloopEntropies, numTriloops, sizeof(struct triloop), comp3loop)))
	   EntropyEnthalpy[0] += loop->value;
      }
   } else if (loopSize == 4) { /* terminal mismatch, tetraloop bonus, hairpin of 4 */
      struct tetraloop* loop;
      if (numTetraloops) {
	 if ((loop = bsearch(numSeq1 + i, tetraloopEnthalpies, numTetraloops, sizeof(struct tetraloop), comp4loop))) {
	    EntropyEnthalpy[1] += loop->value;
	 }
	 if ((loop = bsearch(numSeq1 + i, tetraloopEntropies, numTetraloops, sizeof(struct tetraloop), comp4loop))) {
	    EntropyEnthalpy[0] += loop->value;
	 }
      }
   }
   if(!isFinite(EntropyEnthalpy[1])) {
      EntropyEnthalpy[1] = INFINITY;
      EntropyEnthalpy[0] = -1.0;
   }
   T1 = (EntropyEnthalpy[1] + dplx_init_H) / ((EntropyEnthalpy[0] + dplx_init_S + RC));
   T2 = (EnthalpyDPT(i, j) + dplx_init_H) / ((EntropyDPT(i, j)) + dplx_init_S + RC);
   if(T1 < T2 && traceback == 0) {
      EntropyEnthalpy[0] = EntropyDPT(i, j);
      EntropyEnthalpy[1] = EnthalpyDPT(i, j);
   }   
   return;
}


void calc_bulge_internal(int i, int j, int ii, int jj, double* EntropyEnthalpy, int traceback, int maxLoop) {
   int loopSize1, loopSize2, loopSize;
   double T1, T2;
   double S,H;
   int N, N_loop;
   T1 = T2 = -INFINITY;
   S = MinEntropy;
   H = 0;
   loopSize1 = ii - i - 1;
   loopSize2 = jj - j - 1;
   if(ii < jj) {
      N = ((2 * i)/2);
      N_loop = N;
      if(loopSize1 > 2) N_loop -= (loopSize1 - 2);
      if(loopSize2 > 2) N_loop -= (loopSize2 - 2);
   } else {
      N = ((2 * j)/2);
      N_loop = 2 * jj;
      if(loopSize1 > 2) N_loop -= (loopSize1 - 2);
      if(loopSize2 > 2) N_loop -= (loopSize2 - 2);
      N_loop = (N_loop/2) - 1;
   }
#ifdef DEBUG
   if (ii <= i){
      fputs("Error in calc_bulge_internal(): ii is not greater than i\n", stderr);
   }
   if (jj <= j)
     fputs("Error in calc_bulge_internal(): jj is not greater than j\n", stderr);
#endif
   
#ifdef DEBUG
   if (loopSize1 + loopSize2 > maxLoop) {
      fputs("Error: calc_bulge_internal() called with loopSize1 + loopSize2 > maxLoop\n", stderr);
      return;
   }
#endif
#ifdef DEBUG
   if (loopSize1 == 0 && loopSize2 == 0) {
      fputs("Error: calc_bulge_internal() called with nonsense\n", stderr);
      return;
   }
#endif
   loopSize = loopSize1 + loopSize2 -1;
   if((loopSize1 == 0 && loopSize2 > 0) || (loopSize2 == 0 && loopSize1 > 0)) { /* only bulges have to be considered */
      if(loopSize2 == 1 || loopSize1 == 1) { /* bulge loop of size one is treated differently
					      the intervening nn-pair must be added */
	 
	 if((loopSize2 == 1 && loopSize1 == 0) || (loopSize2 == 0 && loopSize1 == 1)) {
	    H = bulgeLoopEnthalpies[loopSize] +
	      stackEnthalpies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
	    S = bulgeLoopEntropies[loopSize] +
	      stackEntropies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
	 }  
	 H += EnthalpyDPT(i, j);
	 S += EntropyDPT(i, j);
	 if(!isFinite(H)) {
	    H = INFINITY;
	    S = -1.0;
	 }
	 
	 T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
	 T2 = (EnthalpyDPT(ii, jj) + dplx_init_H) / ((EntropyDPT(ii, jj)) + dplx_init_S + RC);
	 
	 if((T1 > T2) || ((traceback && T1 >= T2) || (traceback==1))) {
	    EntropyEnthalpy[0] = S;
	    EntropyEnthalpy[1] = H;
	 }
      } else { /* we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30 */
	 
	 H = bulgeLoopEnthalpies[loopSize] + atPenaltyH(numSeq1[i], numSeq2[j]) + atPenaltyH(numSeq1[ii], numSeq2[jj]);
	 H += EnthalpyDPT(i, j);
	 
	 S = bulgeLoopEntropies[loopSize] + atPenaltyS(numSeq1[i], numSeq2[j]) + atPenaltyS(numSeq1[ii], numSeq2[jj]);
	 S += EntropyDPT(i, j);
	 if(!isFinite(H)) {
	    H = INFINITY;
	    S = -1.0;
	 }
	 
	 T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
	 T2 = (EnthalpyDPT(ii, jj) + dplx_init_H) / (EntropyDPT(ii, jj) + dplx_init_S + RC);
	 if((T1 > T2) || ((traceback && T1 >= T2) || (traceback==1))) {
	    EntropyEnthalpy[0] = S;
	    EntropyEnthalpy[1] = H;
	 }
      }
   } else if (loopSize1 == 1 && loopSize2 == 1) { 
      S = stackint2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] + 
	stackint2Entropies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]]; 
      S += EntropyDPT(i, j);
      
      H = stackint2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
	stackint2Enthalpies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]]; 
      H += EnthalpyDPT(i, j);
      if(!isFinite(H)) {
	 H = INFINITY;
	 S = -1.0;
      }
      
      T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
      T2 = (EnthalpyDPT(ii, jj) + dplx_init_H) / (EntropyDPT(ii, jj) + dplx_init_S + RC);
      
      if((DBL_EQ(T1,T2) == 2) || traceback==1) {
	 if((T1 > T2) || (traceback && T1 >= T2)) {
	    EntropyEnthalpy[0] = S;
	    EntropyEnthalpy[1] = H;
	 }
      }      
      return;
   } else { /* only internal loops */
      H = interiorLoopEnthalpies[loopSize] + tstackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
	tstackEnthalpies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]]
	+ (ILAH * abs(loopSize1 - loopSize2));
      H += EnthalpyDPT(i, j);
      
      S = interiorLoopEntropies[loopSize] + tstackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
	tstackEntropies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]] + (ILAS * abs(loopSize1 - loopSize2));
      S += EntropyDPT(i, j);
      if(!isFinite(H)) {
	 H = INFINITY;
	 S = -1.0;
      }
      T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
      T2 = (EnthalpyDPT(ii, jj) + dplx_init_H) / ((EntropyDPT(ii, jj)) + dplx_init_S + RC);
      if((T1 > T2) || ((traceback && T1 >= T2) || (traceback==1))) {
	 EntropyEnthalpy[0] = S;
	 EntropyEnthalpy[1] = H;
      }
   }
   return;
}

void calc_bulge_internal2(int i, int j, int ii, int jj, double* EntropyEnthalpy, int traceback, int maxLoop) {
   int loopSize1, loopSize2, loopSize;
   double T1, T2;
   double S,H;
   int N, N_loop;
   T1 = T2 = -INFINITY;
   S = MinEntropy;
   H = 0.0;
   loopSize1 = ii - i - 1;
   loopSize2 = j - jj - 1;
   if (loopSize1 + loopSize2 > maxLoop) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = INFINITY;
      return;
   }
   if(i < (len1 -j)) {
      N  = i;
      N_loop = (i - 1);
   } else {
      N = len1-j;
      N_loop = len1 - j - 1;
   }
#ifdef DEBUG
   if (ii <= i)
     fputs("Error in calc_bulge_internal(): ii isn't greater than i\n", stderr);
   if (jj >= j)
     fputs("Error in calc_bulge_internal(): jj isn't less than j\n", stderr);
   if (ii >= jj)
     fputs("Error in calc_bulge_internal(): jj isn't greater than ii\n", stderr);
   
   if ((i <= len1 && len1 < ii) || (jj <= len2 && len2 < j))  {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = INFINITY;
      return;
   }
#endif
   
#ifdef DEBUG
   if (loopSize1 + loopSize2 > maxLoop) {
      fputs("Error: calc_bulge_internal() called with loopSize1 + loopSize2 > maxLoop\n", stderr);
      return;
   }
#endif
#ifdef DEBUG
   if (loopSize1 == 0 && loopSize2 == 0) {
      fputs("Error: calc_bulge_internal() called with nonsense\n", stderr);
      return;
   }
#endif
   
#ifdef DEBUG
   if (i > len1)
     i -= len1;
   if (ii > len1)
     ii -= len1;
   if (j > len2)
     j -= len2;
   if (jj > len2)
     jj -= len2;
#endif
   loopSize = loopSize1 + loopSize2 -1; /* for indx only */
   if((loopSize1 == 0 && loopSize2 > 0) || (loopSize2 == 0 && loopSize1 > 0)) { /* only bulges have to be considered */
      if(loopSize2 == 1 || loopSize1 == 1) { /* bulge loop of size one is treated differently
					      the intervening nn-pair must be added */
	 if((loopSize2 == 1 && loopSize1 == 0) || (loopSize2 == 0 && loopSize1 == 1)) {
	    H = bulgeLoopEnthalpies[loopSize] + 
	      stackEnthalpies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]]; 		 
	    S = bulgeLoopEntropies[loopSize] +
	      stackEntropies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
	 }
	 if(traceback!=1) {
	    H += EnthalpyDPT(ii, jj); /* bulge koos otsaga, st bulge i,j-ni */
	    S += EntropyDPT(ii, jj);
	 }
	 
	 if(!isFinite(H)) {
	    H = INFINITY;
	    S = -1.0;
	 }
	 
	 T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
	 T2 = (EnthalpyDPT(i, j) + dplx_init_H) / ((EntropyDPT(i, j)) + dplx_init_S + RC);
	 
	 if((T1 > T2) || ((traceback && T1 >= T2) || traceback==1)) {
	    EntropyEnthalpy[0] = S;
	    EntropyEnthalpy[1] = H;
	 }
	 
      } else { /* we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30 */
	 
	 H = bulgeLoopEnthalpies[loopSize] + atPenaltyH(numSeq1[i], numSeq2[j]) + atPenaltyH(numSeq1[ii], numSeq2[jj]);
	 if(traceback!=1)
	   H += EnthalpyDPT(ii, jj);
	 
	 S = bulgeLoopEntropies[loopSize] + atPenaltyS(numSeq1[i], numSeq2[j]) + atPenaltyS(numSeq1[ii], numSeq2[jj]);
	 if(traceback!=1)
	   S += EntropyDPT(ii, jj);
	 if(!isFinite(H)) {
	    H = INFINITY;
	    S = -1.0;
	 }
	 T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
	 T2 = (EnthalpyDPT(i, j) + dplx_init_H) / (EntropyDPT(i, j) + dplx_init_S + RC);
	 
	 if((T1 > T2) || ((traceback && T1 >= T2) || (traceback==1))) {
	    EntropyEnthalpy[0] = S;
	    EntropyEnthalpy[1] = H;
	 }
      }      
   } /* end of calculating bulges */
   else if (loopSize1 == 1 && loopSize2 == 1) {
      /* mismatch nearest neighbor parameters */
         
      S = stackint2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] + 
	stackint2Entropies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]];     
      if(traceback!=1)
	S += EntropyDPT(ii, jj);
      
      H = stackint2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
	stackint2Enthalpies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]];
      if(traceback!=1)
	H += EnthalpyDPT(ii, jj);
      if(!isFinite(H)) {
	 H = INFINITY;
	 S = -1.0;
      }      
      T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
      T2 = (EnthalpyDPT(i, j) + dplx_init_H) / (EntropyDPT(i, j) + dplx_init_S + RC);
      if((DBL_EQ(T1,T2) == 2) || traceback) {
	 if((T1 > T2) || ((traceback && T1 >= T2) || traceback==1)) {
	    EntropyEnthalpy[0] = S;
	    EntropyEnthalpy[1] = H;
	 }
      }      
      return;
   } else { /* only internal loops */ 
      
      H = interiorLoopEnthalpies[loopSize] + tstackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
	tstackEnthalpies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]]
	+ (ILAH * abs(loopSize1 - loopSize2));
      if(traceback!=1)
	H += EnthalpyDPT(ii, jj);
      
      S = interiorLoopEntropies[loopSize] + tstackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
	tstackEntropies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]] + (ILAS * abs(loopSize1 - loopSize2));
      if(traceback!=1)
	S += EntropyDPT(ii, jj);
      if(!isFinite(H)) {
	 H = INFINITY;
	 S = -1.0;
      }
      
      T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
      T2 = (EnthalpyDPT(i, j) + dplx_init_H) / ((EntropyDPT(i, j)) + dplx_init_S + RC);
      if((T1 > T2) || ((traceback && T1 >= T2) || (traceback==1))) {
	 EntropyEnthalpy[0] = S;
	 EntropyEnthalpy[1] = H;
      }
   }
   return;
}

void calc_terminal_bp(double temp) { /* compute exterior loop */
   int i;
   int max;
   SEND5(0) = SEND5(1) = -1.0;
   HEND5(0) = HEND5(1) = INFINITY;
   for(i = 2; i<=(len1); i++) {
      SEND5(i) = MinEntropy;
      HEND5(i) = 0;
   }
   
   double T1, T2, T3, T4, T5;
   T1 = T2 = T3 = T4 = T5 = -INFINITY;
   double G;
   /* adding terminal penalties to 3' end and to 5' end */
   for(i = 2; i <= len1; ++i) {
      max = 0;
      T1 = T2 = T3 = T4 = T5 = -INFINITY;      
      T1 = (HEND5(i - 1) + dplx_init_H) / (SEND5(i - 1) + dplx_init_S + RC);
      T2 = (END5_1(i,1) + dplx_init_H) / (END5_1(i,2) + dplx_init_S + RC);
      T3 = (END5_2(i,1) + dplx_init_H) / (END5_2(i,2) + dplx_init_S + RC);
      T4 = (END5_3(i,1) + dplx_init_H) / (END5_3(i,2) + dplx_init_S + RC);
      T5 = (END5_4(i,1) + dplx_init_H) / (END5_4(i,2) + dplx_init_S + RC); 
      max = max5(T1,T2,T3,T4,T5);
      switch (max) {
       case 1: 
	 SEND5(i) = SEND5(i - 1);
	 HEND5(i) = HEND5(i - 1);
	 break;
       case 2:
	 G = END5_1(i,1) - (temp * (END5_1(i,2)));
	 if(G < G2) {
	    SEND5(i) = END5_1(i,2);
	    HEND5(i) = END5_1(i,1);
	 } else {
	    SEND5(i) = SEND5(i - 1);
	    HEND5(i) = HEND5(i - 1);
	 }
	 break;
       case 3:
	 G = END5_2(i,1) - (temp * (END5_2(i,2)));
	 if(G < G2) {
	    SEND5(i) = END5_2(i,2);
	    HEND5(i) = END5_2(i,1);
	 } else {
	    SEND5(i) = SEND5(i - 1);
	    HEND5(i) = HEND5(i - 1);
	 }
	 break;
       case 4:
	 G = END5_3(i,1) - (temp * (END5_3(i,2)));
	 if(G < G2) {
	    SEND5(i) = END5_3(i,2);
	    HEND5(i) = END5_3(i,1);
	 } else {
	    SEND5(i) = SEND5(i - 1);
	    HEND5(i) = HEND5(i - 1);
	 }
	 break;
       case 5:
	 G = END5_4(i,1) - (temp * (END5_4(i,2)));
	 if(G < G2) {
	    SEND5(i) = END5_4(i,2);
	    HEND5(i) = END5_4(i,1);
	 } else {
	    SEND5(i) = SEND5(i - 1);
	    HEND5(i) = HEND5(i - 1);
	 }
	 break;
       default:
#ifdef DEBUG
	 printf ("WARNING: max5 returned character code %d ??\n", max);
#endif
      }
   }
}

double END5_1(int i,int hs) {
   int k, max_tm_flag;
   double max_tm; /* energy min */
   double T1, T2;
   double H, S;
   double H_max, S_max;
   H_max = H = INFINITY;
   S_max = S = -1.0;
   T1 = T2 = -INFINITY;
   max_tm = -INFINITY;
   for(k = 0; k <= i - MIN_HRPN_LOOP - 2; ++k) {
      T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
      T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC); 
      max_tm_flag= max2(T1,T2);
      if(max_tm_flag==1) {
	 H = HEND5(k) + atPenaltyH(numSeq1[k + 1], numSeq1[i]) + EnthalpyDPT(k + 1, i);
	 S = SEND5(k) + atPenaltyS(numSeq1[k + 1], numSeq1[i]) + EntropyDPT(k + 1, i);
	 if(!isFinite(H) || H > 0 || S > 0) { /* H and S must be greater than 0 to avoid BS */
	    H = INFINITY;
	    S = -1.0;
	 } 
	 T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);	
      } else {
	 H = 0 + atPenaltyH(numSeq1[k + 1], numSeq1[i]) + EnthalpyDPT(k + 1, i);
	 S = 0 + atPenaltyS(numSeq1[k + 1], numSeq1[i]) + EntropyDPT(k + 1, i);
	 if(!isFinite(H) || H > 0 || S > 0) {
	    H = INFINITY;
	    S = -1.0;
	 }
	 T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
      }
      max_tm_flag = max2(max_tm, T1);
      if(max_tm_flag == 2) {
	 if(S > MinEntropyCutoff) {
	    H_max = H;
	    S_max = S;
	    max_tm = T1;
	 }
      }
   }
   if (hs == 1) return H_max;
   return S_max;
}

double END5_2(int i,int hs) {
   int k,max_tm_flag;
   double max_tm;
   double T1, T2;
   double H, S;
   double H_max, S_max;
   H_max = H = INFINITY;
   T1 = T2 = max_tm = -INFINITY;
   S_max = S = -1.0;
   for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k) {
      T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
      T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC); 
      max_tm_flag = max2(T1,T2);
      if(max_tm_flag == 1) {
	 H = HEND5(k) + atPenaltyH(numSeq1[k + 2], numSeq1[i]) + Hd5(i, k + 2) + EnthalpyDPT(k + 2, i);
	 S = SEND5(k) + atPenaltyS(numSeq1[k + 2], numSeq1[i]) + Sd5(i, k + 2) + EntropyDPT(k + 2, i);
	 if(!isFinite(H) || H > 0 || S > 0) {
	    H = INFINITY;
	    S = -1.0;
	 }
	 T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);	 
      } else {
	 H = 0 + atPenaltyH(numSeq1[k + 2], numSeq1[i]) + Hd5(i, k + 2) + EnthalpyDPT(k + 2, i);
	 S = 0 + atPenaltyS(numSeq1[k + 2], numSeq1[i]) + Sd5(i, k + 2) + EntropyDPT(k + 2, i);
	 if(!isFinite(H) || H > 0 || S > 0) {
	    H = INFINITY;
	    S = -1.0;
	 }
	 T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
      }
      max_tm_flag = max2(max_tm, T1);
      if(max_tm_flag == 2) {
	 if(S > MinEntropyCutoff) {
	    H_max = H;
	    S_max = S;
	    max_tm = T1;
	 }
      }
   }
   if (hs == 1) return H_max;
   return S_max;
}

double END5_3(int i,int hs) {
   int k, max_tm_flag;
   double max_tm;
   double T1, T2;
   double H, S;
   double H_max, S_max;
   H_max = H = INFINITY;;
   T1 = T2 = max_tm = -INFINITY;
   S_max = S = -1.0;
   for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k) {
      T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
      T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC);       
      max_tm_flag = max2(T1,T2);
      if(max_tm_flag==1) {
	 H = HEND5(k) + atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) + Hd3(i - 1, k + 1) + EnthalpyDPT(k + 1, i - 1);
	 S = SEND5(k) + atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) + Sd3(i - 1, k + 1) + EntropyDPT(k + 1, i - 1);
	 if(!isFinite(H) || H > 0 || S > 0) {
	    H = INFINITY;
	    S = -1.0;
	 }	 
	 T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
      } else {
	 H = 0 + atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) + Hd3(i - 1, k + 1) + EnthalpyDPT(k + 1, i - 1);
	 S = 0 + atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) + Sd3(i - 1, k + 1) + EntropyDPT(k + 1, i - 1);
	 if(!isFinite(H) || H > 0 || S > 0) {
	    H = INFINITY;
	    S = -1.0;
	 }
	 T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
      }
      max_tm_flag = max2(max_tm, T1);
      if(max_tm_flag == 2) {
	 if(S > MinEntropyCutoff) {
	    H_max = H;
	    S_max = S;
	    max_tm = T1;
	 }
      }
   }
   if (hs == 1) return H_max;
   return S_max;
}

double END5_4(int i,int hs){
   int k, max_tm_flag;
   double max_tm;
   double T1, T2;
   double H, S;
   double H_max, S_max;
   H_max = H = INFINITY;
   T1 = T2 = max_tm = -INFINITY;
   S_max = S = -1.0;
   for(k = 0; k <= i - MIN_HRPN_LOOP - 4; ++k) {
      T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
      T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC); 
      max_tm_flag = max2(T1,T2);
      if(max_tm_flag == 1) {
	 H = HEND5(k) + atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) + Htstack(i - 1, k + 2) + EnthalpyDPT(k + 2, i - 1);
	 S = SEND5(k) + atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) + Ststack(i - 1, k + 2) + EntropyDPT(k + 2, i - 1);
	 if(!isFinite(H) || H > 0 || S > 0) {
	    H = INFINITY;
	    S = -1.0;
	 }
	 T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
      } else {
	 H = 0 + atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) + Htstack(i - 1, k + 2) + EnthalpyDPT(k + 2, i - 1);
	 S = 0 + atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) + Ststack(i - 1, k + 2) + EntropyDPT(k + 2, i - 1);
	 if(!isFinite(H) || H > 0 || S > 0) {
	    H = INFINITY;
	    S = -1.0;
	 }
	 T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
      }
      max_tm_flag = max2(max_tm, T1);
      if(max_tm_flag == 2) {
	 if(S > MinEntropyCutoff) {
	    H_max = H;
	    S_max = S;
	    max_tm = T1;
	 }
      }
   }
   if (hs == 1) return H_max;
   return S_max;
}


double Sd5(int i, int j) {
   return dangleEntropies5[numSeq1[i]][numSeq1[j]][numSeq1[j - 1]];
}

double Hd5(int i, int j) {
   return dangleEnthalpies5[numSeq1[i]][numSeq1[j]][numSeq1[j - 1]];
}

double Sd3(int i, int j) {
   return dangleEntropies3[numSeq1[i]][numSeq1[i+1]][numSeq1[j]];
}

double Hd3(int i, int j) {
   return dangleEnthalpies3[numSeq1[i]][numSeq1[i+1]][numSeq1[j]];
}

double Ststack(int i, int j) {
   return tstack2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq1[j]][numSeq1[j-1]];
}

double Htstack(int i, int j) { /* e.g AG_TC 210 */
   return tstack2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq1[j]][numSeq1[j-1]];
}

/* Return 1 if string is symmetrical, 0 otherwise. */
int symmetry_thermo(const unsigned char* seq) {
   register char s;
   register char e;
   const unsigned char *seq_end=seq;
   int i = 0;
   int seq_len=length_unsig_char(seq);
   int mp = seq_len/2;
   if(seq_len%2==1) {
      return 0;
   }
   seq_end+=seq_len;
   seq_end--;
   while(i<mp) {
      i++;
      s=*seq;
      e=*seq_end;
      if ((s=='A' && e!='T')
	  || (s=='T' && e!='A')
	  || (e=='A' && s!='T')
	  || (e=='T' && s!='A')) {
	 return 0;
      }
      if ((s=='C' && e!='G')
	  || (s=='G' && e!='C')
	  || (e=='C' && s!='G')
	  || (e=='G' && s!='C')) {
	 return 0;
      }
      seq++;
      seq_end--;
   }   
   return 1;
}

int length_unsig_char(const unsigned char * str) {
   int i = 0;
   while(*(str++)) {
      i++;
      if(i == INT_MAX)
	return -1;
   }
   return i;
}

void tracebacku(int* bp, int maxLoop,thal_results* o, int a) /* traceback for unimolecular structure */
{
   int i, j;
   i = j = 0;
   int ii, jj, k;
   struct tracer *top, *stack = NULL;
   double* SH1;
   double* SH2;
   double* EntropyEnthalpy;
   SH1 = safe_malloc(2 * sizeof(double),o,a);
   SH2 = safe_malloc(2 * sizeof(double),o,a);
   EntropyEnthalpy = safe_malloc(2 * sizeof(double),o,a);
   push(&stack,len1, 0, 1,o,a);
   while(stack) {
      top = stack;
      stack = stack->next;
      i = top->i;
      j = top->j;
      if(top->mtrx==1) {
	 while (equal(SEND5(i), SEND5(i - 1)) && equal(HEND5(i), HEND5(i - 1))) /* if previous structure is the same as this one */
	   --i;
	 if (i == 0)
	   continue;
	 if (equal(SEND5(i), END5_1(i,2)) && equal(HEND5(i), END5_1(i,1))) {
	    for (k = 0; k <= i - MIN_HRPN_LOOP - 2; ++k)
	      if (equal(SEND5(i), atPenaltyS(numSeq1[k + 1], numSeq1[i]) + EntropyDPT(k + 1, i)) &&
		  equal(HEND5(i), atPenaltyH(numSeq1[k + 1], numSeq1[i]) + EnthalpyDPT(k + 1, i))) {
		 push(&stack, k + 1, i,0,o,a);
		 break;
	      }
	    else if (equal(SEND5(i), SEND5(k) + atPenaltyS(numSeq1[k + 1], numSeq1[i]) + EntropyDPT(k + 1, i)) &&
		     equal(HEND5(i), HEND5(k) + atPenaltyH(numSeq1[k + 1], numSeq1[i]) + EnthalpyDPT(k + 1, i))) {
	       push(&stack, k + 1, i, 0,o,a); 
	       push(&stack, k, 0, 1,o,a);
	       break;
	    }
	 }
	 else if (equal(SEND5(i), END5_2(i,2)) && equal(HEND5(i), END5_2(i,1))) {
	    for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k)
	      if (equal(SEND5(i), atPenaltyS(numSeq1[k + 2], numSeq1[i]) + Sd5(i, k + 2) + EntropyDPT(k + 2, i)) &&
		  equal(HEND5(i), atPenaltyH(numSeq1[k + 2], numSeq1[i]) + Hd5(i, k + 2) + EnthalpyDPT(k + 2, i))) {
		 push(&stack, k + 2, i, 0,o,a);
		 break;
	      }
	    else if (equal(SEND5(i), SEND5(k) + atPenaltyS(numSeq1[k + 2], numSeq1[i]) + Sd5(i, k + 2) + EntropyDPT(k + 2, i)) &&
		     equal(HEND5(i), HEND5(k) + atPenaltyH(numSeq1[k + 2], numSeq1[i]) + Hd5(i, k + 2) + EnthalpyDPT(k + 2, i))) {
	       push(&stack, k + 2, i, 0,o,a);
	       push(&stack, k, 0, 1,o,a);
	       break;
	    }
	 }
	 else if (equal(SEND5(i), END5_3(i,2)) && equal(HEND5(i), END5_3(i,1))) {
	    for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k)
	      if (equal(SEND5(i), atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) + Sd3(i - 1, k + 1) + EntropyDPT(k + 1, i - 1))
		  && equal(HEND5(i), atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) + Hd3(i - 1, k + 1) + EnthalpyDPT(k + 1, i - 1))) {
		 push(&stack, k + 1, i - 1, 0,o,a);
		 break;
	      }
	    else if (equal(SEND5(i), SEND5(k) + atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) + Sd3(i - 1, k + 1) + EntropyDPT(k + 1, i - 1)) &&
		     equal(HEND5(i), HEND5(k) + atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) + Hd3(i - 1, k + 1) + EnthalpyDPT(k + 1, i - 1))) {
	       push(&stack, k + 1, i - 1, 0,o,a); /* matrix 0  */
	       push(&stack, k, 0, 1,o,a); /* matrix 3 */
	       break;
	    }
	 }
	 else if(equal(SEND5(i), END5_4(i,2)) && equal(HEND5(i), END5_4(i,1))) {
	    for (k = 0; k <= i - MIN_HRPN_LOOP - 4; ++k)
	      if (equal(SEND5(i), atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) + Ststack(i - 1, k + 2) + EntropyDPT(k + 2, i - 1)) &&
		  equal(HEND5(i), atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) + Htstack(i - 1, k + 2) + EnthalpyDPT(k + 2, i - 1))) {
		 push(&stack, k + 2, i - 1, 0,o,a);
		 break;
	      }
	    else if (equal(SEND5(i), SEND5(k) + atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) + Ststack(i - 1, k + 2) + EntropyDPT(k + 2, i - 1)) &&
		     equal(HEND5(i), HEND5(k) + atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) + Htstack(i - 1, k + 2) + EnthalpyDPT(k + 2, i - 1)) ) {
	       push(&stack, k + 2, i - 1, 0,o,a);
	       push(&stack, k, 0, 1,o,a);
	       break;
	    }
	 }	      
      }
      else if(top->mtrx==0) {
	 bp[i - 1] = j;
	 bp[j - 1] = i;
	 SH1[0] = -1.0;
	 SH1[1] = INFINITY;
	 calc_hairpin(i, j, SH1, 1); /* 1 means that we use this method in traceback */
	 SH2[0] = -1.0;
	 SH2[1] = INFINITY;
	 CBI(i,j,SH2,2,maxLoop);
	 if (equal(EntropyDPT(i, j), Ss(i, j, 2) + EntropyDPT(i + 1, j - 1)) &&
	     equal(EnthalpyDPT(i, j), Hs(i, j, 2) + EnthalpyDPT(i + 1, j - 1))) {
	    push(&stack, i + 1, j - 1, 0,o,a);
	 }
	 else if (equal(EntropyDPT(i, j), SH1[0]) && equal(EnthalpyDPT(i,j), SH1[1]));
	 else if (equal(EntropyDPT(i, j), SH2[0]) && equal(EnthalpyDPT(i, j), SH2[1])) {
	    int d, done;
	    for (done = 0, d = j - i - 3; d >= MIN_HRPN_LOOP + 1 && d >= j - i - 2 - maxLoop && !done; --d)
	      for (ii = i + 1; ii < j - d; ++ii) {
		 jj = d + ii;
		 EntropyEnthalpy[0] = -1.0;
		 EntropyEnthalpy[1] = INFINITY;
		 calc_bulge_internal2(i, j, ii, jj,EntropyEnthalpy,1,maxLoop);
		 if (equal(EntropyDPT(i, j), EntropyEnthalpy[0] + EntropyDPT(ii, jj)) && 
		     equal(EnthalpyDPT(i, j), EntropyEnthalpy[1] + EnthalpyDPT(ii, jj))) {
		    push(&stack, ii, jj, 0,o,a);
		    ++done;
		    break;
		 }
	      }
	 } else { 
	 }	 
      }
      free(top);
   }   
   free(SH1);
   free(SH2);
   free(EntropyEnthalpy);
}


void traceback(int i, int j, double RT, int* ps1, int* ps2, int maxLoop,thal_results* o, int a)
{
   int d, ii, jj, done;
   double* SH;
   SH = safe_malloc(2 * sizeof(double),o,a);
   ps1[i - 1] = j;
   ps2[j - 1] = i;
   while(1) {
      SH[0] = -1.0;
      SH[1] = INFINITY;
      LSH(i,j,SH);
      if(equal(EntropyDPT(i,j),SH[0]) && equal(EnthalpyDPT(i,j),SH[1])) {
	 break;
      }
      done = 0;
      if (i > 1 && j > 1 && equal(EntropyDPT(i,j), Ss(i - 1, j - 1, 1) + EntropyDPT(i - 1, j - 1))) {
	 i = i - 1;
	 j = j - 1;
	 ps1[i - 1] = j;
	 ps2[j - 1] = i;
	 done = 1;
      }
      for (d = 3; !done && d <= maxLoop + 2; ++d) {
	 ii = i - 1;
	 jj = -ii - d + (j + i);
	 if (jj < 1) {
	    ii -= abs(jj-1);
	    jj = 1;
	 }
	 for (; !done && ii > 0 && jj < j; --ii, ++jj) {
	    SH[0] = -1.0;
	    SH[1] = INFINITY;
	    calc_bulge_internal(ii, jj, i, j, SH,1,maxLoop);
	    if (equal(EntropyDPT(i, j), SH[0]) && equal(EnthalpyDPT(i, j), SH[1])) {
	       i = ii;
	       j = jj;
	       ps1[i - 1] = j;
	       ps2[j - 1] = i;
	       done = 1;
	       break;
	    }
	 }
      }
   }
   free(SH);
}

void drawDimer(int* ps1, int* ps2, double temp, double H, double S, int temponly, double t37, thal_results *o, int a) {
   int i, j, k, numSS1, numSS2, N;
   char* duplex[4];
   double G, t;
   t = G = 0;
   if (!isFinite(temp)){
      if(temponly==0) {
	 printf("No predicted secondary structures for given sequences\n");	
      }
      o->temp = 0.0; /* lets use generalization here; this should rather be very negative value */
      o->msg = "No predicted sec struc for given seq";
      return;
   } else {
      N=0;
      for(i=0;i<len1;i++){
	 if(ps1[i]>0) ++N;      
      }
      for(i=0;i<len2;i++) {
	 if(ps2[i]>0) ++N;
      }
      N = (N/2) -1;
      t = ((H) / (S + (N * saltCorrection) + RC)) - ABSOLUTE_ZERO;
      if(temponly==0) {
	 G = (H) - (t37 * (S + (N * saltCorrection)));
	 S = S + (N * saltCorrection);
	 o->temp = (double) t / PRECISION;
	 /* maybe user does not need as precise as that */
	 /* printf("Thermodynamical values:\t%d\tdS = %g\tdH = %g\tdG = %g\tt = %g\tN = %d, SaltC=%f, RC=%f\n", 
		len1, (double) S / PRECISION, (double) H / PRECISION, (double) G / PRECISION, (double) t / PRECISION, (int) N, saltCorrection, RC); */
	 printf("Calculated thermodynamical parameters for dimer:\tdS = %g\tdH = %g\tdG = %g\tt = %g\n", 
		(double) S / PRECISION, (double) H / PRECISION, (double) G / PRECISION, (double) t / PRECISION);
      } else {
	 o->temp = (double) t / PRECISION;
	 return;
      }
   }

   duplex[0] = safe_malloc(len1 + len2 + 1,o,a);
   duplex[1] = safe_malloc(len1 + len2 + 1,o,a);
   duplex[2] = safe_malloc(len1 + len2 + 1,o,a);
   duplex[3] = safe_malloc(len1 + len2 + 1,o,a);
   duplex[0][0] = duplex[1][0] = duplex[2][0] = duplex[3][0] = 0;
   
   i = 0;
   numSS1 = 0;
   while (ps1[i++] == 0) ++numSS1;
   j = 0;
   numSS2 = 0;
   while (ps2[j++] == 0) ++numSS2;
   
   if (numSS1 >= numSS2){
      for (i = 0; i < numSS1; ++i) {
	 strcatc(duplex[0], oligo1[i]);
	 strcatc(duplex[1], ' ');
	 strcatc(duplex[2], ' ');
      }
      for (j = 0; j < numSS1 - numSS2; ++j) strcatc(duplex[3], ' ');
      for (j = 0; j < numSS2; ++j) strcatc(duplex[3], oligo2[j]);
   } else {
      for (j = 0; j < numSS2; ++j) {
	 strcatc(duplex[3], oligo2[j]);
	 strcatc(duplex[1], ' ');
	 strcatc(duplex[2], ' ');
      }	
      for (i = 0; i < numSS2 - numSS1; ++i)
	strcatc(duplex[0], ' ');
      for (i = 0; i < numSS1; ++i)
	strcatc(duplex[0], oligo1[i]);
   }
   i = numSS1 + 1;
   j = numSS2 + 1;
   
   while (i <= len1) {
      while (i <= len1 && ps1[i - 1] != 0 && j <= len2 && ps2[j - 1] != 0) {
	 strcatc(duplex[0], ' ');
	 strcatc(duplex[1], oligo1[i - 1]);
	 strcatc(duplex[2], oligo2[j - 1]);
	 strcatc(duplex[3], ' ');
	 ++i;
	 ++j;
      }
      numSS1 = 0;
      while (i <= len1 && ps1[i - 1] == 0) {
	 strcatc(duplex[0], oligo1[i - 1]);
	 strcatc(duplex[1], ' ');
	 ++numSS1;
	 ++i;
      }
      numSS2 = 0;
      while (j <= len2 && ps2[j - 1] == 0) {
	 strcatc(duplex[2], ' ');
	 strcatc(duplex[3], oligo2[j - 1]);
	 ++numSS2;
	 ++j;
      }
      if (numSS1 < numSS2)
	for (k = 0; k < numSS2 - numSS1; ++k) {
	   strcatc(duplex[0], '-');
	   strcatc(duplex[1], ' ');
	}
      else if (numSS1 > numSS2)
	for (k = 0; k < numSS1 - numSS2; ++k) {
	   strcatc(duplex[2], ' ');
	   strcatc(duplex[3], '-');
	}
   }
   printf("SEQ\t");
   printf("%s\n", duplex[0]);
   printf("SEQ\t");
   printf("%s\n", duplex[1]);
   printf("STR\t");
   printf("%s\n", duplex[2]);
   printf("STR\t");
   printf("%s\n", duplex[3]);
   
   free(duplex[0]);
   free(duplex[1]);
   free(duplex[2]);
   free(duplex[3]);
   
   return;
}

void drawHairpin(int* bp, double mh, double ms, int temponly, double temp, thal_results *o, int a)
{
   /* Plain text */
   int i, N;
   N = 0;
   double mg, t;
   if (!isFinite(ms) || !isFinite(mh)) {
      if(temponly == 0) {
	 printf("0\tdS = %g\tdH = %g\tinf\tinf\n", (double) ms / PRECISION,(double) mh / PRECISION);
#ifdef DEBUG
	 fputs("No temperature could be calculated\n",stderr);
#endif
      } else {
	 o->temp = 0.0; /* lets use generalization here */
	 o->msg = "No predicted sec struc for given seq\n";
      }
   } else {
      if(temponly == 0) {
	 for (i = 1; i < len1; ++i) {
	    if(bp[i-1] > 0) N++;
	 }
      } else {
	 for (i = 1; i < len1; ++i) {
	    if(bp[i-1] > 0) N++;
	 }	 
      }      
      t = (mh / (ms + (((N/2)-1) * saltCorrection))) - ABSOLUTE_ZERO;
      if(temponly == 0) {
	 mg = mh - (temp * (ms + (((N/2)-1) * saltCorrection)));      
	 ms = ms + (((N/2)-1) * saltCorrection);
	 o->temp = (double) t / PRECISION;
	 printf("Calculated thermodynamical parameters for dimer:\t%d\tdS = %g\tdH = %g\tdG = %g\tt = %g\n", 
		len1, (double) ms / PRECISION, (double) mh / PRECISION, (double) mg / PRECISION, (double) t / PRECISION);
      } else {
	 o->temp = (double) t / PRECISION;
	 return;
      }      
   } 
   /* plain-text output */
   char* asciiRow;
   asciiRow = safe_malloc(len1,o,a);
   for(i = 0; i < len1; ++i) asciiRow[i] = '0'; 
   for(i = 1; i < len1+1; ++i) {
      if(bp[i-1] == 0) {
	 asciiRow[(i-1)] = '-';
      } else {
	 if(bp[i-1] > (i-1)) {
	    asciiRow[(bp[i-1]-1)]='\\';
	 } else  {
	    asciiRow[(bp[i-1]-1)]='/';
	 }
      }  
   }
   printf("SEQ\t");
   for(i = 0; i < len1; ++i) printf("%c",asciiRow[i]);
   printf("\nSTR\t%s\n", oligo1);
   free(asciiRow);
   return;
}


int equal(double a, double b) {
#ifdef INTEGER
   return a == b;
#endif
   
   if (!finite(a) || !finite(b))
     return 0;
   return fabs(a - b) < 1e-5;
   
   if (a == 0 && b == 0)
     return 1;
}

void strcatc(char* str, char c) {
   str[strlen(str) + 1] = 0;
   str[strlen(str)] = c;
}
