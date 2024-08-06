/*
 Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2009,2010,
               2011,2012
 Whitehead Institute for Biomedical Research, Steve Rozen
 (http://purl.com/STEVEROZEN/), and Helen Skaletsky
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

#if defined(__sun)
#include <ieeefp.h>
#endif

#include "thal.h"
#include "thal_default_params.h"

#define STR(X) #X
#define LONG_SEQ_ERR_STR(MAX_LEN) "Target sequence length > maximum allowed (" STR(MAX_LEN) ") in thermodynamic alignment"
#define XSTR(X) STR(X)

#define INIT_BUF_SIZE 1024

#ifdef INTEGER
# define isFinite(x) (x < _INFINITY / 2)
#else
# define isFinite(x) isfinite(x)
#endif

/*** BEGIN CONSTANTS ***/
// static const double _INFINITY is defined in thal_default_params.h
static const int min_hrpn_loop = 3;
static const double R = 1.9872; /* cal/Kmol */
static const double ILAS = (-300 / 310.15); /* Internal Loop Entropy ASymmetry correction -0.3kcal/mol*/
static const double ILAH = 0.0; /* Internal Loop EntHalpy Asymmetry correction */
static const double AT_H = 2200.0; /* AT penalty */
static const double AT_S = 6.9; /* AT penalty */
static const double MinEntropyCutoff = -2500.0; /* to filter out non-existing entropies */
static const double MinEntropy = -3224.0; /* initiation */
const double ABSOLUTE_ZERO = 273.15;
const double TEMP_KELVIN = 310.15;
const int MAX_LOOP = 30; /* the maximum size of loop that can be calculated; for larger loops formula must be implemented */
const int MIN_LOOP = 0;
//static const char BASES[5] = {'A', 'C', 'G', 'T', 'N'}; /* bases to be considered - N is every symbol that is not A, G, C,$
//                                                  */
//static const char BASE_PAIRS[4][4] = {"A-T", "C-G", "G-C", "T-A" }; /* allowed basepairs */
/* matrix for allowed; bp 0 - no bp, watson crick bp - 1 */
static const int bpIndx[5][5] =  {
     {0, 0, 0, 1, 0}, /* A, C, G, T, N; */
     {0, 0, 1, 0, 0},
     {0, 1, 0, 0, 0},
     {1, 0, 0, 0, 0},
     {0, 0, 0, 0, 0}};

/*** END OF CONSTANTS ***/

/*** BEGIN STRUCTs ***/
/*
Defined in thal_default_params.h:
struct triloop
struct tetraloop
*/

struct tracer /* structure for traceback_monomer - unimolecular str */ {
  int i;
  int j;
  int mtrx; /* [0 1] EntropyDPT/EnthalpyDPT*/
  struct tracer* next;
};

/*** END STRUCTs ***/


//=====================================================================================
//Functions for dimer calculation
//=====================================================================================
 /* initiates thermodynamic parameter tables of entropy and enthalpy for dimer */
static void initMatrix_dimer(double **entropyDPT, double **enthalpyDPT, const unsigned char *numSeq1, 
                                 const unsigned char *numSeq2, int oligo1_len, int oligo2_len);
 /* calc-s thermod values into dynamic progr table (dimer) */
static void fillMatrix_dimer(int maxLoop, double **entropyDPT, double **enthalpyDPT, double RC,
                                 double dplx_init_S, double dplx_init_H, const unsigned char *numSeq1,
                                 const unsigned char *numSeq2, int oligo1_len, int oligo2_len, thal_results* o);
static void calc_bulge_internal_dimer(int ii, int jj, int i, int j, double* EntropyEnthalpy, const double *saved_RSH,
                                 int traceback, int maxLoop, const double *const *entropyDPT, const double *const *enthalpyDPT,
                                 double RC, double dplx_init_S, double dplx_init_H,
                                 const unsigned char *numSeq1, const unsigned char *numSeq2);
static void traceback_dimer(int i, int j, double RC, int* ps1, int* ps2, int maxLoop, const double *const *entropyDPT,
                                 const double *const *enthalpyDPT, double dplx_init_S, double dplx_init_H, const unsigned char *numSeq1,
                                 const unsigned char *numSeq2, int oligo1_len, int oligo2_len, thal_results* o);
char *drawDimer(int*, int*, double, double, const thal_mode mode, double, const unsigned char *oligo1, const unsigned char *oligo2,
               double saltCorrection, double RC, int oligo1_len, int oligo2_len, jmp_buf, thal_results *);
/* calculate terminal entropy S and terminal enthalpy H starting reading from 5'end */
static void LSH(int i, int j, double* EntropyEnthalpy, double RC, double dplx_init_S,
               double dplx_init_H, const unsigned char *numSeq1, const unsigned char *numSeq2);

//=====================================================================================
//Functions for dimer and hairpin calculation
//=====================================================================================
/* calculate terminal entropy S and terminal enthalpy H starting reading from 3'end */
static void RSH(int i, int j, double* EntropyEnthalpy, double RC, double dplx_init_S,
               double dplx_init_H, const unsigned char *numSeq1, const unsigned char *numSeq2);

//=====================================================================================
//Functions for hairpin calculation
//=====================================================================================
static void initMatrix_monomer(double **entropyDPT, double **enthalpyDPT, const unsigned char *numSeq1, const unsigned char *numSeq2,
                              int oligo1_len, int oligo2_len); /* initiates thermodynamic parameter tables of entropy and enthalpy for monomer */
 /* calcs thermod values into dynamic progr table (monomer) */
static void fillMatrix_monomer(int maxLoop, double **entropyDPT, double **enthalpyDPT, double RC, const unsigned char *numSeq1,
                              const unsigned char *numSeq2, int oligo1_len, int oligo2_len, thal_results* o);
static void calc_bulge_internal_monomer(int ii, int jj, int i, int j, double* EntropyEnthalpy, int traceback, int maxLoop,
                                    const  double *const *entropyDPT, const double *const *enthalpyDPT, double RC,
                                    const unsigned char *numSeq1, const unsigned char *numSeq2);
static void calc_terminal_bp(double temp, const double *const *entropyDPT, const double *const *enthalpyDPT, double *send5, double *hend5,
                              double RC, const unsigned char *numSeq1, const unsigned char *numSeq2, int oligo1_len, int oligo2_len);
/* finds monomer structure that has maximum Tm */
static void calc_hairpin(int i, int j, double* EntropyEnthalpy, int traceback, const double *const *entropyDPT, const double *const *enthalpyDPT,
                        double RC, const unsigned char *numSeq1, const unsigned char *numSeq2, int oligo1_len, int oligo2_len);
static void push(struct tracer**, int, int, int, jmp_buf, thal_results*); /* to add elements to struct */
static void traceback_monomer(int*, int, const double *const *entropyDPT, const double *const *enthalpyDPT, double *send5, double *hend5, double RC,
                              double dplx_init_S, double dplx_init_H,  const unsigned char *numSeq1, const unsigned char *numSeq2,
                              int oligo1_len, int oligo2_len, jmp_buf, thal_results*);
char *drawHairpin(int*, double, double, const thal_mode mode, double, const unsigned char *oligo1, const unsigned char *oligo2, double saltCorrection,
                  int oligo1_len, int oligo2_len, jmp_buf, thal_results *);

//=====================================================================================
//Misc helper functions
//=====================================================================================
static int comp3loop(const void*, const void*); /* checks if sequnece consists of specific triloop */
static int comp4loop(const void*, const void*); /* checks if sequnece consists of specific tetraloop */
static int equal(double a, double b);

//=====================================================================================
//Initializing functions
//=====================================================================================
static int symmetry_thermo(const unsigned char* seq);
static double saltCorrectS (double mv, double dv, double dntp); /* part of calculating salt correction
                                                                   for Tm by SantaLucia et al */
int thal_check_errors(const unsigned char *oligo_f, const unsigned char *oligo_r, int *len_f, int *len_r, const thal_args *a, thal_results *o);

//=====================================================================================
//Functions for allocating memory
//=====================================================================================
static void* safe_calloc(size_t, size_t, jmp_buf _jmp_buf, thal_results* o);
static void* safe_malloc(size_t, jmp_buf, thal_results* o);
static void* safe_realloc(void*, size_t, jmp_buf, thal_results* o);
double **allocate_DPT(int oligo1_len, int oligo2_len, jmp_buf _jmp_buf, thal_results *o);
void free_DPT(double **dpt);

//=====================================================================================
//Functions for string manipulation
//=====================================================================================
static int length_unsig_char(const unsigned char * str); /* returns length of unsigned char; to avoid warnings while compiling */
static unsigned char str2int(char c); /* converts DNA sequence to int; 0-A, 1-C, 2-G, 3-T, 4-whatever */
static void reverse(unsigned char *s);
/* Is sequence symmetrical */
static void save_append_string(char** ret, int *space, thal_results *o, const char *str, jmp_buf);
static void save_append_char(char** ret, int *space, thal_results *o, const char str, jmp_buf);
static void strcatc(char*, char);
static double readDouble(char **str, jmp_buf, thal_results* o);

//=====================================================================================
//Functions for changing thermodynamic parameters
//=====================================================================================
static char* readParamFile(const char* dirname, const char* fname, jmp_buf, thal_results* o); /* file of thermodynamic params */
static void readLoop(char **str, double *v1, double *v2, double *v3, jmp_buf, thal_results *o);
static int readTLoop(char **str, char *s, double *v, int triloop, jmp_buf, thal_results *o);
static void getStack(double stackEntropies[5][5][5][5], double stackEnthalpies[5][5][5][5], const thal_parameters *tp, jmp_buf, thal_results* o);
static void getStackint2(double stackEntropiesint2[5][5][5][5], double stackint2Enthalpies[5][5][5][5], const thal_parameters *tp, jmp_buf, thal_results* o);
static void getDangle(double dangleEntropies3[5][5][5], double dangleEnthalpies3[5][5][5], double dangleEntropies5[5][5][5],
                      double dangleEnthalpies5[5][5][5], const thal_parameters *tp, jmp_buf, thal_results* o);
static void getTstack(double tstackEntropies[5][5][5][5], double tstackEnthalpies[5][5][5][5], const thal_parameters *tp, jmp_buf, thal_results* o);
static void getTstack2(double tstack2Entropies[5][5][5][5], double tstack2Enthalpies[5][5][5][5], const thal_parameters *tp, jmp_buf, thal_results* o);
static void getTriloop(struct triloop**, struct triloop**, int* num, const thal_parameters *tp, jmp_buf _jmp_buf, thal_results* o);
static void getTetraloop(struct tetraloop**, struct tetraloop**, int* num, const thal_parameters *tp, jmp_buf _jmp_buf, thal_results* o);
static void getLoop(double hairpinLoopEnntropies[30], double interiorLoopEntropies[30], double bulgeLoopEntropiess[30],
             double hairpinLoopEnthalpies[30], double interiorLoopEnthalpies[30], double bulgeLoopEnthalpies[30], const thal_parameters *tp, jmp_buf, thal_results* o);
static void tableStartATS(double atp_value, double atp[5][5]); /* creates table of entropy values for nucleotides
                                                                  to which AT-penlty must be applied */
static void tableStartATH(double atp_value, double atp[5][5]);

/*
Thermodynamic parameters from thal_default_params.h:

static double atpS[5][5];  AT penalty 
static double atpH[5][5];  AT penalty 
static int numTriloops;  hairpin triloop penalties 
static int numTetraloops;  hairpin tetraloop penalties 
static double dangleEntropies3[5][5][5]; thermodynamic paramteres for 3' dangling ends 
static double dangleEnthalpies3[5][5][5];  ther params for 3' dangling ends 
static double dangleEntropies5[5][5][5];   ther params for 5' dangling ends 
static double dangleEnthalpies5[5][5][5];  ther params for 5' dangling ends 
static double stackEntropies[5][5][5][5];  ther params for perfect match pairs 
static double stackEnthalpies[5][5][5][5];  ther params for perfect match pairs 
static double stackint2Entropies[5][5][5][5]; ther params for perfect match and internal mm 
static double stackint2Enthalpies[5][5][5][5];  ther params for perfect match and internal mm
static double interiorLoopEntropies[30];  interior loop params according to length of the loop 
static double bulgeLoopEntropies[30];  bulge loop params according to length of the loop 
static double hairpinLoopEntropies[30];  hairpin loop params accordint to length of the loop 
static double interiorLoopEnthalpies[30];  same as interiorLoopEntropies but values of entropy 
static double bulgeLoopEnthalpies[30];  same as bulgeLoopEntropies but values of entropy 
static double hairpinLoopEnthalpies[30];  same as hairpinLoopEntropies but values of entropy 
static double tstackEntropies[5][5][5][5];  ther params for terminal T0mismatches 
static double tstackEnthalpies[5][5][5][5];  ther params for terminal mismatches 
static double tstack2Entropies[5][5][5][5];  ther params for internal terminal mismatches 
static double tstack2Enthalpies[5][5][5][5];  ther params for internal terminal mismatches 
static struct triloop* triloopEntropies;  ther penalties for given triloop seq-s 
static struct triloop* triloopEnthalpies;  ther penalties for given triloop seq-s 
static struct tetraloop* tetraloopEntropies;  ther penalties for given tetraloop seq-s 
static struct tetraloop* tetraloopEnthalpies;  ther penalties for given tetraloop seq-s 
*/
//static double *send5, *hend5; /* calc 5'  */
/* w/o init not constant anymore, cause for unimolecular and bimolecular foldings there are different values */
//static double dplx_init_H; /* initiation enthalpy; for duplex 200, for unimolecular structure 0 */
//static double dplx_init_S; /* initiation entropy; for duplex -5.7, for unimoleculat structure 0 */
//static double saltCorrection; /* value calculated by saltCorrectS, includes correction for monovalent and divalent cations */
//static double RC; /* universal gas constant multiplied w DNA conc - for melting temperature */
//static int bestI, bestJ; /* starting position of most stable str */
//static double** enthalpyDPT; /* matrix for values of enthalpy */
//static double** entropyDPT; /* matrix for values of entropy */
//static unsigned char *oligo1, *oligo2; /* inserted oligo sequenced */
//static unsigned char *numSeq1, *numSeq2; /* same as oligo1 and oligo2 but converted to numbers */
//static int oligo1_len, oligo2_len; /* length of sequense 1 and 2 *//* 17.02.2009 int temponly;*/ /* print only temperature of the predicted structure */

/* central method: execute all sub-methods for calculating secondary
   structure for dimer or for monomer */
void 
thal(const unsigned char *oligo_f, 
     const unsigned char *oligo_r, 
     const thal_args *a,
     const thal_mode mode,
     thal_results *o)
{
   double SH[2];
   int i, j;
   int len_f, len_r;
   int k;
   int *bp;
   unsigned char *oligo2_rev = NULL;
   double mh, ms;
   double G1, bestG;
   jmp_buf _jmp_buf;
   double saltCorrection;
   double RC;
   double dplx_init_S;
   double dplx_init_H;
   int bestI;
   int bestJ;
   int oligo1_len;
   int oligo2_len;

   unsigned char *numSeq1 = NULL;
   unsigned char *numSeq2 = NULL;
   unsigned char *oligo1 = NULL;
   unsigned char *oligo2 = NULL;
   strcpy(o->msg, "");
   o->temp = THAL_ERROR_SCORE;
   errno = 0; 

   if (setjmp(_jmp_buf) != 0) {
     o->temp = THAL_ERROR_SCORE;
     return;  /* If we get here, that means we returned via a
                 longjmp.  In this case errno might be ENOMEM,
                 but not necessarily. */
   }

   if(thal_check_errors(oligo_f, oligo_r, &len_f, &len_r, a, o))
      return;

   if(a->type!=3) {
      oligo1 = (unsigned char*) safe_malloc((len_f + 1) * sizeof(unsigned char), _jmp_buf, o);
      oligo2 = (unsigned char*) safe_malloc((len_r + 1) * sizeof(unsigned char), _jmp_buf, o);
      strcpy((char*)oligo1,(const char*)oligo_f);
      strcpy((char*)oligo2,(const char*)oligo_r);
   } else  {
      oligo1 = (unsigned char*) safe_malloc((len_r + 1) * sizeof(unsigned char), _jmp_buf, o);
      oligo2 = (unsigned char*) safe_malloc((len_f + 1) * sizeof(unsigned char), _jmp_buf, o);
      strcpy((char*)oligo1,(const char*)oligo_r);
      strcpy((char*)oligo2,(const char*)oligo_f);
   }
   oligo1_len = length_unsig_char(oligo1);
   oligo2_len = length_unsig_char(oligo2);
   double **enthalpyDPT = allocate_DPT(oligo1_len, oligo2_len, _jmp_buf, o);
   double **entropyDPT = allocate_DPT(oligo1_len, oligo2_len, _jmp_buf, o);
   /*** INIT values for unimolecular and bimolecular structures ***/
   if (a->type==4) { /* unimolecular folding */
      dplx_init_H = 0.0;
      dplx_init_S = -0.00000000001;
      RC=0;
   } else  {
      /* hybridization of two oligos */
      dplx_init_H = 200;
      dplx_init_S = -5.7;
      if(symmetry_thermo(oligo1) && symmetry_thermo(oligo2)) {
         RC = R  * log(a->dna_conc/1000000000.0);
      } else {
         RC = R  * log(a->dna_conc/4000000000.0);
      }
      if(a->type!=3) {
         oligo2_rev = (unsigned char*) safe_malloc((length_unsig_char(oligo_r) + 1) * sizeof(unsigned char), _jmp_buf, o);
         strcpy((char*)oligo2_rev,(const char*)oligo_r);
      } else {
         oligo2_rev = (unsigned char*) safe_malloc((length_unsig_char(oligo_f) + 1) * sizeof(unsigned char), _jmp_buf, o);
         strcpy((char*)oligo2_rev,(const char*)oligo_f);
      }
      reverse(oligo2_rev); /* REVERSE oligo2, so it goes to dpt 3'->5' direction */
      free(oligo2);
      oligo2=NULL;
      oligo2=&oligo2_rev[0];
   }
   /* convert nucleotides to numbers */
   numSeq1 = (unsigned char*) safe_realloc(numSeq1, oligo1_len + 2, _jmp_buf, o);
   numSeq2 = (unsigned char*) safe_realloc(numSeq2, oligo2_len + 2, _jmp_buf, o);

   /*** Calc part of the salt correction ***/
   saltCorrection=saltCorrectS(a->mv,a->dv,a->dntp); /* salt correction for entropy, must be multiplied with N, which is
                                                   the total number of phosphates in the duplex divided by 2; 8bp dplx N=7 */

   for(i = 0; i < oligo1_len; i++) oligo1[i] = toupper(oligo1[i]);
   for(i = 0; i < oligo2_len; i++) oligo2[i] = toupper(oligo2[i]);
   for(i = 1; i <= oligo1_len; ++i) numSeq1[i] = str2int(oligo1[i - 1]);
   for(i = 1; i <= oligo2_len; ++i) numSeq2[i] = str2int(oligo2[i - 1]);
   numSeq1[0] = numSeq1[oligo1_len + 1] = numSeq2[0] = numSeq2[oligo2_len + 1] = 4; /* mark as N-s */

   if (a->type==4) { /* calculate structure of monomer */
      double *hend5 = NULL;
      double *send5 = NULL;
      send5 = (double*) safe_realloc(send5, (oligo1_len + 1) * sizeof(double), _jmp_buf, o);
      hend5 = (double*) safe_realloc(hend5, (oligo1_len + 1) * sizeof(double), _jmp_buf, o);
      initMatrix_monomer(entropyDPT, enthalpyDPT, numSeq1, numSeq2, oligo1_len, oligo2_len);
      fillMatrix_monomer(a->maxLoop, entropyDPT, enthalpyDPT, RC, numSeq1, numSeq2, oligo1_len, oligo2_len, o);
      calc_terminal_bp(a->temp, (const double **)entropyDPT, (const double **)enthalpyDPT, send5, hend5, RC, numSeq1, numSeq2, oligo1_len, oligo2_len);
      mh = hend5[oligo1_len];
      ms = send5[oligo1_len];
      o->align_end_1 = (int) mh;
      o->align_end_2 = (int) ms;
      bp = (int*) safe_calloc(oligo1_len, sizeof(int), _jmp_buf, o);
      for (k = 0; k < oligo1_len; ++k) bp[k] = 0;
      if(isFinite(mh)) {
        traceback_monomer(bp, a->maxLoop, (const double **)entropyDPT, (const double **)enthalpyDPT, send5, hend5, RC, dplx_init_S, dplx_init_H, numSeq1, numSeq2, oligo1_len, oligo2_len, _jmp_buf, o);
        /* traceback for unimolecular structure */
        o->sec_struct=drawHairpin(bp, mh, ms, mode,a->temp, oligo1, oligo2, saltCorrection, oligo1_len, oligo2_len, _jmp_buf, o); /* if mode=THL_FAST or THL_DEBUG_F then return after printing basic therm data */
      } else if((mode != THL_FAST) && (mode != THL_DEBUG_F) && (mode != THL_STRUCT)) {
        fputs("No secondary structure could be calculated\n",stderr);
      }

      if(o->temp==-_INFINITY && (!strcmp(o->msg, ""))) o->temp=0.0;
      free(bp);
      free_DPT(enthalpyDPT);
      free_DPT(entropyDPT);
      free(numSeq1);
      free(numSeq2);
      free(send5);
      free(hend5);
      free(oligo1);
      free(oligo2);
      return;
   } else if(a->type!=4) { /* Hybridization of two moleculs */
      initMatrix_dimer(entropyDPT, enthalpyDPT, numSeq1, numSeq2, oligo1_len, oligo2_len);
      fillMatrix_dimer(a->maxLoop, entropyDPT, enthalpyDPT, RC, dplx_init_S, dplx_init_H, numSeq1, numSeq2, oligo1_len, oligo2_len, o);
      /* calculate terminal basepairs */
      bestI = bestJ = 0; 
      G1 = bestG = _INFINITY;
      if(a->type==1)
        for (i = 1; i <= oligo1_len; i++) {
           for (j = 1; j <= oligo2_len; j++) {
              RSH(i, j, SH, RC, dplx_init_S, dplx_init_H, numSeq1, numSeq2);
              G1 = (enthalpyDPT[i][j]+ SH[1] + dplx_init_H) - TEMP_KELVIN*(entropyDPT[i][j] + SH[0] + dplx_init_S);  
              if(G1<bestG){
                 bestG = G1;
                 bestI = i;
                 bestJ = j;
              }
           }
        }
      int *ps1, *ps2;
      ps1 = (int*) safe_calloc(oligo1_len, sizeof(int), _jmp_buf, o);
      ps2 = (int*) safe_calloc(oligo2_len, sizeof(int), _jmp_buf, o);
      for (i = 0; i < oligo1_len; ++i)
        ps1[i] = 0;
      for (j = 0; j < oligo2_len; ++j)
        ps2[j] = 0;
      if(a->type == 2 || a->type == 3)        {
         /* THAL_END1 */
         bestI = bestJ = 0;
         bestI = oligo1_len;
         i = oligo1_len;
         G1 = bestG = _INFINITY;
         for (j = 1; j <= oligo2_len; ++j) {
            RSH(i, j, SH, RC, dplx_init_S, dplx_init_H, numSeq1, numSeq2);
            G1 = (enthalpyDPT[i][j]+ SH[1] + dplx_init_H) - TEMP_KELVIN*(entropyDPT[i][j] + SH[0] + dplx_init_S);  
                if(G1<bestG){
                   bestG = G1;
                         bestJ = j;
            }
         }
      }
      if (!isFinite(bestG)) bestI = bestJ = 1;
      double dH, dS;
      RSH(bestI, bestJ, SH, RC, dplx_init_S, dplx_init_H, numSeq1, numSeq2);
      dH = enthalpyDPT[bestI][bestJ]+ SH[1] + dplx_init_H;
      dS = (entropyDPT[bestI][bestJ] + SH[0] + dplx_init_S);
      /* tracebacking */
      for (i = 0; i < oligo1_len; ++i)
        ps1[i] = 0;
      for (j = 0; j < oligo2_len; ++j)
        ps2[j] = 0;
      if(isFinite(enthalpyDPT[bestI][bestJ])){
         traceback_dimer(bestI, bestJ, RC, ps1, ps2, a->maxLoop, (const double **)entropyDPT, (const double **)enthalpyDPT, dplx_init_S, dplx_init_H, numSeq1, numSeq2, oligo1_len, oligo2_len, o);
         o->sec_struct=drawDimer(ps1, ps2, dH, dS, mode, a->temp, oligo1, oligo2, saltCorrection, RC, oligo1_len, oligo2_len, _jmp_buf, o);
         o->align_end_1=bestI;
         o->align_end_2=bestJ;
      } else  {
         o->temp = 0.0;
         /* fputs("No secondary structure could be calculated\n",stderr); */
      }
      free(ps1);
      free(ps2);
      free(oligo2_rev);
      free_DPT(enthalpyDPT);
      free_DPT(entropyDPT);
      free(numSeq1);
      free(numSeq2);
      free(oligo1);
      return;
   }
   return;
}
/*** END thal() ***/

//=====================================================================================
//Functions for dimer calculation
//=====================================================================================

static void 
initMatrix_dimer(double **entropyDPT, double** enthalpyDPT, const unsigned char *numSeq1, const unsigned char *numSeq2, int oligo1_len, int oligo2_len)
{
   int i, j;
   for (i = 1; i <= oligo1_len; ++i) {
      for (j = 1; j <= oligo2_len; ++j) {
         if (bpIndx[numSeq1[i]][numSeq2[j]] == 0)  {
            enthalpyDPT[i][j] = _INFINITY;
            entropyDPT[i][j] = -1.0;
         } else {
            enthalpyDPT[i][j] = 0.0;
            entropyDPT[i][j] = MinEntropy;
         }
      }
   }
}

static void 
fillMatrix_dimer(int maxLoop, double **entropyDPT, double **enthalpyDPT, double RC, double dplx_init_S, double dplx_init_H, const unsigned char *numSeq1, const unsigned char *numSeq2, int oligo1_len, int oligo2_len, thal_results *o)
{
   int d, i, j, ii, jj;
   double SH[2];
   double saved_RSH[2];
   double T0, T1;
   double S0, S1;
   double H0, H1;

   for (i = 1; i <= oligo1_len; ++i) {
      for (j = 1; j <= oligo2_len; ++j) {
         if(isFinite(enthalpyDPT[i][j])) { /* if finite */
            SH[0] = -1.0;
            SH[1] = _INFINITY;
            LSH(i, j, SH, RC, dplx_init_S, dplx_init_H, numSeq1, numSeq2);
            if(isFinite(SH[1])) {
               entropyDPT[i][j] = SH[0];
               enthalpyDPT[i][j] = SH[1];
            }
            if (i > 1 && j > 1) {
               T0 = T1 = -_INFINITY;
               S0 = entropyDPT[i][j];
               H0 = enthalpyDPT[i][j];
               RSH(i, j, SH, RC, dplx_init_S, dplx_init_H, numSeq1, numSeq2);
               saved_RSH[0] = SH[0];
               saved_RSH[1] = SH[1];
               T0 = (H0 + dplx_init_H + SH[1]) /(S0 + dplx_init_S + SH[0] + RC); /* at current position */
               if(isFinite(enthalpyDPT[i-1][j-1]) && isFinite(stackEnthalpies[numSeq1[i-1]][numSeq1[i]][numSeq2[j-1]][numSeq2[j]])) {
                  S1 = entropyDPT[i-1][j-1] + stackEntropies[numSeq1[i-1]][numSeq1[i]][numSeq2[j-1]][numSeq2[j]];
                  H1 = enthalpyDPT[i-1][j-1] + stackEnthalpies[numSeq1[i-1]][numSeq1[i]][numSeq2[j-1]][numSeq2[j]];
                  T1 = (H1 + dplx_init_H + SH[1]) /(S1 + dplx_init_S + SH[0] + RC);
               } else {
                  S1 = -1.0;
                  H1 = _INFINITY;
                  T1 = (H1 + dplx_init_H) /(S1 + dplx_init_S + RC);
               }
               
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
               if(T1 > T0) { 
                  entropyDPT[i][j] = S1;
                  enthalpyDPT[i][j] = H1;
               } else if(T0 >= T1) {
                  entropyDPT[i][j] = S0;
                  enthalpyDPT[i][j] = H0;
               }

               for(d = 3; d <= maxLoop + 2; d++) { /* max=30, length over 30 is not allowed */
                  ii = i - 1;
                  jj = - ii - d + (j + i);
                  if (jj < 1) {
                     ii -= abs(jj-1);
                     jj = 1;
                  }
                  for (; ii > 0 && jj < j; --ii, ++jj) {
                     SH[0] = -1.0;
                     SH[1] = _INFINITY;
                     if (isFinite(enthalpyDPT[ii][jj])) {
                        calc_bulge_internal_dimer(ii, jj, i, j, SH, saved_RSH, 0, maxLoop, (const double **)entropyDPT, (const double **)enthalpyDPT, RC, dplx_init_S, dplx_init_H, numSeq1, numSeq2);
                        if(SH[0] < MinEntropyCutoff) {
                           /* to not give dH any value if dS is unreasonable */
                           SH[0] = MinEntropy;
                           SH[1] = 0.0;
                        }
                        if(isFinite(SH[1])) {
                           enthalpyDPT[i][j] = SH[1];
                           entropyDPT[i][j] = SH[0];
                        }
                     }
                  }
               }
            } /* if */
         }
      } /* for */
   } /* for */
}

static void 
calc_bulge_internal_dimer(int i, int j, int ii, int jj, double* EntropyEnthalpy, const double *saved_RSH,
                        int traceback, int maxLoop, const double *const *entropyDPT,  const double *const *enthalpyDPT, double RC,
                        double dplx_init_S, double dplx_init_H, const unsigned char *numSeq1, const unsigned char *numSeq2)
{
   int loopSize1, loopSize2, loopSize;
   double S,H,G1,G2;
   double SH[2];
   SH[0] = -1.0;
   SH[1] = _INFINITY;
   S = -1.0;
   H = _INFINITY;
   loopSize1 = ii - i - 1;
   loopSize2 = jj - j - 1;
   loopSize = loopSize1 + loopSize2-1;
   if((loopSize1 == 0 && loopSize2 > 0) || (loopSize2 == 0 && loopSize1 > 0)) { /* only bulges have to be considered */
      if(loopSize2 == 1 || loopSize1 == 1) { /* bulge loop of size one is treated differently
                                              the intervening nn-pair must be added */

         if((loopSize2 == 1 && loopSize1 == 0) || (loopSize2 == 0 && loopSize1 == 1)) {
            H = bulgeLoopEnthalpies[loopSize] +
              stackEnthalpies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
            S = bulgeLoopEntropies[loopSize] +
              stackEntropies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
         }
         if((H > 0) || (S > 0)){
            H = _INFINITY;
            S = -1.0;
         }
         H += enthalpyDPT[i][j];
         S += entropyDPT[i][j];
         if(!isFinite(H)) {
            H = _INFINITY;
            S = -1.0;
         } 
         if((isFinite(H)) || (traceback==1)) {
            SH[0] = saved_RSH[0];
            SH[1] = saved_RSH[1];
            G1 = H+SH[1] -TEMP_KELVIN*(S+SH[0]);
            G2 = enthalpyDPT[ii][jj]+SH[1]-TEMP_KELVIN*((entropyDPT[ii][jj]+SH[0]));
            if((G1< G2) || (traceback==1)) {
               EntropyEnthalpy[0] = S;
               EntropyEnthalpy[1] = H;
            }
         }
      } else { /* we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30 */

         H = bulgeLoopEnthalpies[loopSize] + atpH[numSeq1[i]][numSeq2[j]] + atpH[numSeq1[ii]][numSeq2[jj]];
         H += enthalpyDPT[i][j];

         S = bulgeLoopEntropies[loopSize] + atpS[numSeq1[i]][numSeq2[j]] + atpS[numSeq1[ii]][numSeq2[jj]];
         S += entropyDPT[i][j];
         if(!isFinite(H)) {
            H = _INFINITY;
            S = -1.0;
         }
         if((H > 0) && (S > 0)){ 
            H = _INFINITY;
            S = -1.0;
         }
         if((isFinite(H)) || (traceback==1)){ 
            SH[0] = saved_RSH[0];
            SH[1] = saved_RSH[1];
            G1 = H+SH[1] -TEMP_KELVIN*(S+SH[0]);
            G2 = enthalpyDPT[ii][jj]+SH[1]-TEMP_KELVIN*(entropyDPT[ii][jj]+SH[0]);
            if(G1< G2 || (traceback==1)){
               EntropyEnthalpy[0] = S;
               EntropyEnthalpy[1] = H;
            }
         }
      }
   } else if (loopSize1 == 1 && loopSize2 == 1) {
      S = stackint2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
        stackint2Entropies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]];
      S += entropyDPT[i][j];

      H = stackint2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
        stackint2Enthalpies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]];
      H += enthalpyDPT[i][j];
      if(!isFinite(H)) {
         H = _INFINITY;
         S = -1.0;
      }
      if((H > 0) && (S > 0)){
         H = _INFINITY;
         S = -1.0;
      }    
      if((isFinite(H)) || (traceback==1)){
         SH[0] = saved_RSH[0];
         SH[1] = saved_RSH[1];
         G1 = H+SH[1] -TEMP_KELVIN*(S+SH[0]);
         G2 = enthalpyDPT[ii][jj]+SH[1]-TEMP_KELVIN*(entropyDPT[ii][jj]+SH[0]);
         if((G1< G2) || traceback==1) {
            EntropyEnthalpy[0] = S;
            EntropyEnthalpy[1] = H;
         }
      }
      return;
   } else { /* only internal loops */
      H = interiorLoopEnthalpies[loopSize] + tstackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
        tstackEnthalpies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]]
        + (ILAH * abs(loopSize1 - loopSize2));
      H += enthalpyDPT[i][j];

      S = interiorLoopEntropies[loopSize] + tstackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
        tstackEntropies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]] + (ILAS * abs(loopSize1 - loopSize2));
      S += entropyDPT[i][j];
      if(!isFinite(H)) {
         H = _INFINITY;
         S = -1.0;
      }
   if((H > 0) && (S > 0)){ 
         H = _INFINITY;
         S = -1.0;
      }
      if((isFinite(H)) || (traceback==1)){
         SH[0] = saved_RSH[0];
         SH[1] = saved_RSH[1];
         G1 = H+SH[1] -TEMP_KELVIN*(S+SH[0]);
         G2 = enthalpyDPT[ii][jj]+SH[1]-TEMP_KELVIN*(entropyDPT[ii][jj]+SH[0]);
         if((G1< G2) || (traceback==1)){
                  EntropyEnthalpy[0] = S;
                  EntropyEnthalpy[1] = H;
            }
         }
   }
   return;
}

static void 
traceback_dimer(int i, int j, double RC, int* ps1, int* ps2, int maxLoop, const double *const *entropyDPT, const double *const *enthalpyDPT,
               double dplx_init_S, double dplx_init_H, const unsigned char *numSeq1, const unsigned char *numSeq2, int oligo1_len,
               int oligo2_len, thal_results* o)
{
   int d, ii, jj, done;
   double SH[2];
   double saved_RSH[2];
   ps1[i - 1] = j;
   ps2[j - 1] = i;
   while(1) {
      SH[0] = -1.0;
      SH[1] = _INFINITY;
      LSH(i,j,SH, RC, dplx_init_S, dplx_init_H, numSeq1, numSeq2);
      if(equal(entropyDPT[i][j],SH[0]) && equal(enthalpyDPT[i][j],SH[1])) {
         break;
      }
      done = 0;
      if (i > 1 && j > 1 && equal(entropyDPT[i][j], stackEntropies[numSeq1[i-1]][numSeq1[i]][numSeq2[j-1]][numSeq2[j]] + entropyDPT[i-1][j-1]) 
         && equal(enthalpyDPT[i][j], stackEnthalpies[numSeq1[i-1]][numSeq1[i]][numSeq2[j-1]][numSeq2[j]] + enthalpyDPT[i-1][j-1])) {
         i = i - 1;
         j = j - 1;
         ps1[i - 1] = j;
         ps2[j - 1] = i;
         done = 1;
      }
      for (d = 3; !done && d <= maxLoop + 2; ++d) {
         RSH(i, j, saved_RSH, RC, dplx_init_S, dplx_init_H, numSeq1, numSeq2);
         ii = i - 1;
         jj = -ii - d + (j + i);
         if (jj < 1) {
            ii -= abs(jj-1);
            jj = 1;
         }
         for (; !done && ii > 0 && jj < j; --ii, ++jj) {
            SH[0] = -1.0;
            SH[1] = _INFINITY;
            calc_bulge_internal_dimer(ii, jj, i, j, SH, saved_RSH, 1, maxLoop, entropyDPT, enthalpyDPT, RC, dplx_init_S, dplx_init_H, numSeq1, numSeq2);
            if (equal(entropyDPT[i][j], SH[0]) && equal(enthalpyDPT[i][j], SH[1])) {
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
}

char * 
drawDimer(int* ps1, int* ps2, double H, double S, const thal_mode mode, double t37, const unsigned char *oligo1, const unsigned char *oligo2,
         double saltCorrection, double RC, int oligo1_len, int oligo2_len, jmp_buf _jmp_buf, thal_results *o)
{
   int  ret_space = 0;
   char *ret_ptr = NULL;
   int ret_nr, ret_pr_once;
   char ret_para[400];
   char* ret_str[4];
   int i, j, k, numSS1, numSS2, N;
   char* duplex[4];
   double G, t;
   t = G = 0;
   N=0;
   for(i=0;i<oligo1_len;i++){
      if(ps1[i]>0) ++N;
   }
   for(i=0;i<oligo2_len;i++) {
      if(ps2[i]>0) ++N;
   }
   N = (N/2) -1;
   t = ((H) / (S + (N * saltCorrection) + RC)) - ABSOLUTE_ZERO;
   G = (H) - (t37 * (S + (N * saltCorrection)));
   S = S + (N * saltCorrection);
   o->dg = G;
   o->ds = S;
   o->dh = H;
   o->temp = (double) t;
   if((mode != THL_FAST) && (mode != THL_DEBUG_F)) {
      /* maybe user does not need as precise as that */
      /* printf("Thermodynamical values:\t%d\tdS = %g\tdH = %g\tdG = %g\tt = %g\tN = %d, SaltC=%f, RC=%f\n",
               oligo1_len, (double) S, (double) H, (double) G, (double) t, (int) N, saltCorrection, RC); */
      if (mode != THL_STRUCT) {
         printf("Calculated thermodynamical parameters for dimer:\tdS = %g\tdH = %g\tdG = %g\tt = %g\n",
               (double) S, (double) H, (double) G, (double) t);
      } else {
         snprintf(ret_para, 400, "Tm: %.1f&deg;C  dG: %.0f cal/mol  dH: %.0f cal/mol  dS: %.0f cal/mol*K\\n",
                  (double) t, (double) G, (double) H, (double) S);
      }
   } else {
      return NULL;
   }

   duplex[0] = (char*) safe_malloc(oligo1_len + oligo2_len + 1, _jmp_buf, o);
   duplex[1] = (char*) safe_malloc(oligo1_len + oligo2_len + 1, _jmp_buf, o);
   duplex[2] = (char*) safe_malloc(oligo1_len + oligo2_len + 1, _jmp_buf, o);
   duplex[3] = (char*) safe_malloc(oligo1_len + oligo2_len + 1, _jmp_buf, o);
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

   while (i <= oligo1_len) {
      while (i <= oligo1_len && ps1[i - 1] != 0 && j <= oligo2_len && ps2[j - 1] != 0) {
         strcatc(duplex[0], ' ');
         strcatc(duplex[1], oligo1[i - 1]);
         strcatc(duplex[2], oligo2[j - 1]);
         strcatc(duplex[3], ' ');
         ++i;
         ++j;
      }
      numSS1 = 0;
      while (i <= oligo1_len && ps1[i - 1] == 0) {
         strcatc(duplex[0], oligo1[i - 1]);
         strcatc(duplex[1], ' ');
         ++numSS1;
         ++i;
      }
      numSS2 = 0;
      while (j <= oligo2_len && ps2[j - 1] == 0) {
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
   if ((mode == THL_GENERAL) || (mode == THL_DEBUG)) {
     printf("SEQ\t");
     printf("%s\n", duplex[0]);
     printf("SEQ\t");
     printf("%s\n", duplex[1]);
     printf("STR\t");
     printf("%s\n", duplex[2]);
     printf("STR\t");
     printf("%s\n", duplex[3]);
   }
   if (mode == THL_STRUCT) {
     ret_str[3] = NULL;
     ret_str[0] = (char*) safe_malloc(oligo1_len + oligo2_len + 10, _jmp_buf, o);
     ret_str[1] = (char*) safe_malloc(oligo1_len + oligo2_len + 10, _jmp_buf, o);
     ret_str[2] = (char*) safe_malloc(oligo1_len + oligo2_len + 10, _jmp_buf, o);
     ret_str[0][0] = ret_str[1][0] = ret_str[2][0] = '\0';

     /* Join top primer */
     strcpy(ret_str[0], "   ");
     strcat(ret_str[0], duplex[0]);
     ret_nr = 0;
     while (duplex[1][ret_nr] != '\0') {
       if (duplex[1][ret_nr] == 'A' || duplex[1][ret_nr] == 'T' || 
           duplex[1][ret_nr] == 'C' || duplex[1][ret_nr] == 'G' || 
           duplex[1][ret_nr] == '-') {
         ret_str[0][ret_nr + 3] = duplex[1][ret_nr];
       }
       ret_nr++;
     }
     if (strlen(duplex[1]) > strlen(duplex[0])) {
       ret_str[0][strlen(duplex[1]) + 3] = '\0';
     }
     /* Clean Ends */
     ret_nr = strlen(ret_str[0]) - 1;
     while (ret_nr > 0 && (ret_str[0][ret_nr] == ' ' || ret_str[0][ret_nr] == '-')) {
       ret_str[0][ret_nr--] = '\0';
     }
     /* Write the 5' */
     ret_nr = 3;
     ret_pr_once = 1;
     while (ret_str[0][ret_nr] != '\0' && ret_pr_once == 1) {
       if (ret_str[0][ret_nr] == 'A' || ret_str[0][ret_nr] == 'T' ||
           ret_str[0][ret_nr] == 'C' || ret_str[0][ret_nr] == 'G' ||
           ret_str[0][ret_nr] == '-') {
         ret_str[0][ret_nr - 3] = '5';
         ret_str[0][ret_nr - 2] = '\'';
         ret_pr_once = 0;
       }
       ret_nr++;
     }

     /* Create the align tics */
     strcpy(ret_str[1], "     ");
     for (i = 0 ; i < strlen(duplex[1]) ; i++) {
       if (duplex[1][i] == 'A' || duplex[1][i] == 'T' || 
           duplex[1][i] == 'C' || duplex[1][i] == 'G' ) {
         ret_str[1][i + 3] = '|';
       } else {
         ret_str[1][i + 3] = ' ';
       }
       ret_str[1][i + 4] = '\0';
     }
     /* Clean Ends */
     ret_nr = strlen(ret_str[1]) - 1;
     while (ret_nr > 0 && ret_str[1][ret_nr] == ' ') {
       ret_str[1][ret_nr--] = '\0';
     }
     /* Join bottom primer */
     strcpy(ret_str[2], "   ");
     strcat(ret_str[2], duplex[2]);
     ret_nr = 0;
     while (duplex[3][ret_nr] != '\0') {
       if (duplex[3][ret_nr] == 'A' || duplex[3][ret_nr] == 'T' ||
           duplex[3][ret_nr] == 'C' || duplex[3][ret_nr] == 'G' ||
           duplex[3][ret_nr] == '-') {
         ret_str[2][ret_nr + 3] = duplex[3][ret_nr];
       }
       ret_nr++;
     }
     if (strlen(duplex[3]) > strlen(duplex[2])) {
       ret_str[2][strlen(duplex[3]) + 3] = '\0';
     }
     /* Clean Ends */
     ret_nr = strlen(ret_str[2]) - 1;
     while (ret_nr > 0 && (ret_str[2][ret_nr] == ' ' || ret_str[2][ret_nr] == '-')) {
       ret_str[2][ret_nr--] = '\0';
     }
     /* Write the 5' */
     ret_nr = 3;
     ret_pr_once = 1;
     while (ret_str[2][ret_nr] != '\0' && ret_pr_once == 1) {
       if (ret_str[2][ret_nr] == 'A' || ret_str[2][ret_nr] == 'T' ||
           ret_str[2][ret_nr] == 'C' || ret_str[2][ret_nr] == 'G' ||
           ret_str[2][ret_nr] == '-') {
         ret_str[2][ret_nr - 3] = '3';
         ret_str[2][ret_nr - 2] = '\'';
         ret_pr_once = 0;
       }
       ret_nr++;
     }

     save_append_string(&ret_str[3], &ret_space, o, ret_para, _jmp_buf);
     save_append_string(&ret_str[3], &ret_space, o, ret_str[0], _jmp_buf);
     save_append_string(&ret_str[3], &ret_space, o, " 3\'\\n", _jmp_buf);
     save_append_string(&ret_str[3], &ret_space, o, ret_str[1], _jmp_buf);
     save_append_string(&ret_str[3], &ret_space, o, "\\n", _jmp_buf);
     save_append_string(&ret_str[3], &ret_space, o, ret_str[2], _jmp_buf);
     save_append_string(&ret_str[3], &ret_space, o, " 5\'\\n", _jmp_buf);

     ret_ptr = (char *) safe_malloc(strlen(ret_str[3]) + 1, _jmp_buf, o);
     strcpy(ret_ptr, ret_str[3]);
     if (ret_str[3]) {
       free(ret_str[3]);
     }
     free(ret_str[0]);
     free(ret_str[1]);
     free(ret_str[2]);
   }
   free(duplex[0]);
   free(duplex[1]);
   free(duplex[2]);
   free(duplex[3]);

   return ret_ptr;
}

static void 
LSH(int i, int j, double* EntropyEnthalpy, double RC, double dplx_init_S, double dplx_init_H, const unsigned char *numSeq1, const unsigned char *numSeq2)
{
   double S1, H1, T1, G1;
   double S2, H2, T2, G2;
   S1 = S2 = -1.0;
   H1 = H2 = -_INFINITY;
   T1 = T2 = -_INFINITY;
   if (bpIndx[numSeq1[i]][numSeq2[j]] == 0) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      return;
   }
   S1 = atpS[numSeq1[i]][numSeq2[j]] + tstack2Entropies[numSeq2[j]][numSeq2[j-1]][numSeq1[i]][numSeq1[i-1]];
   H1 = atpH[numSeq1[i]][numSeq2[j]] + tstack2Enthalpies[numSeq2[j]][numSeq2[j-1]][numSeq1[i]][numSeq1[i-1]];
   G1 = H1 - TEMP_KELVIN*S1;
   if(G1>0) {
      H1 = _INFINITY;
      S1 = -1.0;
      G1 = 1.0;
   }
   /** If there is two dangling ends at the same end of duplex **/
   if((bpIndx[numSeq1[i-1]][numSeq2[j-1]] != 1 ) && isFinite(dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]]) && isFinite(dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]])) {
      S2 = atpS[numSeq1[i]][numSeq2[j]] + dangleEntropies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]] +
        dangleEntropies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      H2 = atpH[numSeq1[i]][numSeq2[j]] + dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]] +
        dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(G2>0) {
         H2 = _INFINITY;
         S2 = -1.0;
     G2 = 1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
         T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
         if(T1<T2 && G2<0) {
            S1 = S2;
            H1 = H2;
            T1 = T2;
         }
      } else if(G2<0) {
         S1 = S2;
         H1 = H2;
         T1 = T2;
      }
   } else if ((bpIndx[numSeq1[i-1]][numSeq2[j-1]] != 1) && isFinite(dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]])) {
      S2 = atpS[numSeq1[i]][numSeq2[j]] + dangleEntropies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]];
      H2 = atpH[numSeq1[i]][numSeq2[j]] + dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(G2>0) {
         H2 = _INFINITY;
         S2 = -1.0;
     G2 = 1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(G1<0) {
         T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
         if(T1<T2 && G2<0) {
            S1 = S2;
            H1 = H2;
            T1 = T2;
         }
      } else if (G2<0){
         S1 = S2;
         H1 = H2;
         T1 = T2;
      }
   } else if ((bpIndx[numSeq1[i-1]][numSeq2[j-1]] != 1) && isFinite(dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]])) {
      S2 = atpS[numSeq1[i]][numSeq2[j]] + dangleEntropies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      H2 = atpH[numSeq1[i]][numSeq2[j]] + dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(G2>0) {
         H2 = _INFINITY;
         S2 = -1.0;
     G2 = 1.0;     
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(G1<0) {
         T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
         if(T1 < T2  && G2<0) {
            S1 = S2;
            H1 = H2;
            T1 = T2;
         }
      } else if(G2<0) {
         S1 = S2;
         H1 = H2;
         T1 = T2;
      }
   }
   S2 = atpS[numSeq1[i]][numSeq2[j]];
   H2 = atpH[numSeq1[i]][numSeq2[j]];
   T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
   G1 = H1 -TEMP_KELVIN*S1;   
   G2 = H2 -TEMP_KELVIN*S2;
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

//=====================================================================================
//Functions used in both dimer and hairpin calculations
//=====================================================================================

static void 
RSH(int i, int j, double* EntropyEnthalpy, double RC, double dplx_init_S, double dplx_init_H, const unsigned char *numSeq1, const unsigned char *numSeq2)
{
   double G1, G2;
   double S1, S2;
   double H1, H2;
   double T1, T2;
   S1 = S2 = -1.0;
   H1 = H2 = _INFINITY;
   T1 = T2 = -_INFINITY;
   if (bpIndx[numSeq1[i]][numSeq2[j]] == 0) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      return;
   }
   S1 = atpS[numSeq1[i]][numSeq2[j]] + tstack2Entropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   H1 = atpH[numSeq1[i]][numSeq2[j]] + tstack2Enthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   G1 = H1 - TEMP_KELVIN*S1;
   if(G1>0) {
      H1 = _INFINITY;
      S1 = -1.0;
      G1 = 1.0;
   }
   
   if(bpIndx[numSeq1[i+1]][numSeq2[j+1]] == 0 && isFinite(dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]]) && isFinite(dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]])) {
      S2 = atpS[numSeq1[i]][numSeq2[j]] + dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +
        dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
      H2 = atpH[numSeq1[i]][numSeq2[j]] + dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +
        dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
    G2 = H2 - TEMP_KELVIN*S2;
      if(G2>0) {
         H2 = _INFINITY;
         S2 = -1.0;
     G2 = 1.0;
      }
      
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(G1<0) {
         T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
         if(T1 < T2 && G2<0) {
            S1 = S2;
            H1 = H2;
            T1 = T2;
         }
      } else if(G2<0){
         S1 = S2;
         H1 = H2;
         T1 = T2;
      }
   }

   else if(bpIndx[numSeq1[i+1]][numSeq2[j+1]] == 0 && isFinite(dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]])) {
      S2 = atpS[numSeq1[i]][numSeq2[j]] + dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]];
      H2 = atpH[numSeq1[i]][numSeq2[j]] + dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(G2 >0) {
         H2 = _INFINITY;
         S2 = -1.0;
         G2 = 1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(G1<0) {
         T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
         if(T1 < T2 && G2<0) {
            S1 = S2;
            H1 = H2;
            T1 = T2;
         }
      } else if(G2<0){
         S1 = S2;
         H1 = H2;
         T1 = T2;
      }
   }

   else if(bpIndx[numSeq1[i+1]][numSeq2[j+1]] == 0 && isFinite(dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]])) {
      S2 = atpS[numSeq1[i]][numSeq2[j]] + dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
      H2 = atpH[numSeq1[i]][numSeq2[j]] + dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(G2>0) {
         H2 = _INFINITY;
         S2 = -1.0;
     G2 = 1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(G1<0) {
         T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
         if(T1 < T2 && G2<0) {
            S1 = S2;
            H1 = H2;
            T1 = T2;
         }
      } else if (G2<0){
         S1 = S2;
         H1 = H2;
         T1 = T2;
      }
   }
   S2 = atpS[numSeq1[i]][numSeq2[j]];
   H2 = atpH[numSeq1[i]][numSeq2[j]];
   T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
   G1 = H1 -TEMP_KELVIN*S1;
   G2 =  H2 -TEMP_KELVIN*S2;
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

//=====================================================================================
//Functions used in hairpin calculations
//=====================================================================================

static void 
initMatrix_monomer(double **entropyDPT, double **enthalpyDPT, const unsigned char *numSeq1, const unsigned char *numSeq2, int oligo1_len, int oligo2_len)
{
   int i, j;
   for (i = 1; i <= oligo1_len; ++i)
     for (j = i; j <= oligo2_len; ++j)
       if (j - i < min_hrpn_loop + 1 || (bpIndx[numSeq1[i]][numSeq1[j]] == 0)) {
          enthalpyDPT[i][j] = _INFINITY;
          entropyDPT[i][j] = -1.0;
       } else {
          enthalpyDPT[i][j] = 0.0;
          entropyDPT[i][j] = MinEntropy;
       }

}
static void 
fillMatrix_monomer(int maxLoop, double **entropyDPT, double **enthalpyDPT, double RC, const unsigned char *numSeq1, const unsigned char *numSeq2, int oligo1_len, int oligo2_len, thal_results* o)
{
   int i, j;
   double SH[2];
   double T0, T1;
   double S0, S1;
   double H0, H1;

   for (j = 2; j <= oligo2_len; ++j)
      for (i = j - min_hrpn_loop - 1; i >= 1; --i) {
         if (isFinite(enthalpyDPT[i][j])) {
            SH[0] = -1.0;
            SH[1] = _INFINITY;
            T0 = T1 = -_INFINITY;
            S0 = entropyDPT[i][j];
            H0 = enthalpyDPT[i][j];
            T0 = (H0) /(S0 + RC);
            if(isFinite(enthalpyDPT[i][j])) {
               S1 = (entropyDPT[i + 1][j - 1] + stackEntropies[numSeq2[i]][numSeq2[i+1]][numSeq2[j]][numSeq2[j-1]]);
               H1 = (enthalpyDPT[i + 1][j - 1] + stackEnthalpies[numSeq2[i]][numSeq2[i+1]][numSeq2[j]][numSeq2[j-1]]);
            } else {
               S1 = -1.0;
               H1 = _INFINITY;
            }
            T1 = (H1) /(S1 + RC);
            if(S1 < MinEntropyCutoff) {
               S1 = MinEntropy;
               H1 = 0.0;
            }
            if(S0 < MinEntropyCutoff) {
               S0 = MinEntropy;
               H0 = 0.0;
            }

            if(T1 > T0) {
               entropyDPT[i][j] = S1;
               enthalpyDPT[i][j] = H1;
            } else {
               entropyDPT[i][j] = S0;
               enthalpyDPT[i][j] = H0;
            }

            int d, ii, jj;
            for (d = j - i - 3; d >= min_hrpn_loop + 1 && d >= j - i - 2 - maxLoop; --d){
               for (ii = i + 1; ii < j - d && ii <= oligo1_len; ++ii) {
                  jj = d + ii;
                  SH[0] = -1.0;
                  SH[1] = _INFINITY;
                  if (isFinite(enthalpyDPT[ii][jj]) && isFinite(enthalpyDPT[i][j])) {
                     calc_bulge_internal_monomer(i, j, ii, jj, SH, 0, maxLoop, (const double **)entropyDPT, (const double **)enthalpyDPT, RC, numSeq1, numSeq2);
                     if(isFinite(SH[1])) {
                        if(SH[0] < MinEntropyCutoff) {
                           SH[0] = MinEntropy;
                           SH[1] = 0.0;
                        }
                        enthalpyDPT[i][j] = SH[1];
                        entropyDPT[i][j] = SH[0];
                     }
                  }
               }
            }
           SH[0] = -1.0;
           SH[1] = _INFINITY;
           calc_hairpin(i, j, SH, 0, (const double **)entropyDPT, (const double **)enthalpyDPT, RC, numSeq1, numSeq2, oligo1_len, oligo2_len);
           if(isFinite(SH[1])) {
              if(SH[0] < MinEntropyCutoff){ /* to not give dH any value if dS is unreasonable */
                 SH[0] = MinEntropy;
                 SH[1] = 0.0;
              }
              entropyDPT[i][j] = SH[0];
              enthalpyDPT[i][j] = SH[1];
           }
        }
     }
}

static void 
calc_hairpin(int i, int j, double* EntropyEnthalpy, int traceback, const double *const *entropyDPT, const double *const *enthalpyDPT, double RC, const unsigned char *numSeq1, const unsigned char *numSeq2, int oligo1_len, int oligo2_len)
{
   int loopSize = j - i - 1;
   double G1, G2;
   G1 = G2 = -_INFINITY;
   double SH[2];
   SH[0] = -1.0;
   SH[1] = _INFINITY;
   if(loopSize < min_hrpn_loop) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      return;
   }
   if (i <= oligo1_len && oligo2_len < j) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      return;
   } else if (i > oligo2_len) {
      i -= oligo1_len;
      j -= oligo2_len;
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
      EntropyEnthalpy[1] += atpH[numSeq1[i]][numSeq1[j]];
      EntropyEnthalpy[0] += atpS[numSeq1[i]][numSeq1[j]];
   }

   if (loopSize == 3) {         /* closing AT-penalty (+), triloop bonus, hairpin of 3 (+) */
      struct triloop* loop;
      if (numTriloops) {
         if ((loop = (struct triloop*) bsearch(numSeq1 + i, triloopEnthalpies, numTriloops, sizeof(struct triloop), comp3loop)))
           EntropyEnthalpy[1] += loop->value;
         if ((loop = (struct triloop*) bsearch(numSeq1 + i, triloopEntropies, numTriloops, sizeof(struct triloop), comp3loop)))
           EntropyEnthalpy[0] += loop->value;
      }
   } else if (loopSize == 4) { /* terminal mismatch, tetraloop bonus, hairpin of 4 */
      struct tetraloop* loop;
      if (numTetraloops) {
         if ((loop = (struct tetraloop*) bsearch(numSeq1 + i, tetraloopEnthalpies, numTetraloops, sizeof(struct tetraloop), comp4loop))) {
            EntropyEnthalpy[1] += loop->value;
         }
         if ((loop = (struct tetraloop*) bsearch(numSeq1 + i, tetraloopEntropies, numTetraloops, sizeof(struct tetraloop), comp4loop))) {
            EntropyEnthalpy[0] += loop->value;
         }
      }
   }
   if(!isFinite(EntropyEnthalpy[1])) {
      EntropyEnthalpy[1] = _INFINITY;
      EntropyEnthalpy[0] = -1.0;
   }
   if((EntropyEnthalpy[1] > 0) && (EntropyEnthalpy[0] > 0) && ((enthalpyDPT[i][j] <= 0) || (entropyDPT[i][j] <= 0))) { /* if both, S and H are positive */
      EntropyEnthalpy[1] = _INFINITY;
      EntropyEnthalpy[0] = -1.0;
   }
   RSH(i, j, SH, RC, 0.0, 0.0, numSeq1, numSeq2);
   G1 = EntropyEnthalpy[1]+SH[1] -TEMP_KELVIN*(EntropyEnthalpy[0]+SH[0]);
   G2 = enthalpyDPT[i][j]+SH[1] -TEMP_KELVIN*(entropyDPT[i][j]+SH[0]);
     if(G2 < G1 && traceback == 0) {
      EntropyEnthalpy[0] = entropyDPT[i][j];
      EntropyEnthalpy[1] = enthalpyDPT[i][j];
   }
   return;
}

static void 
calc_bulge_internal_monomer(int i, int j, int ii, int jj, double* EntropyEnthalpy, int traceback, int maxLoop,
                        const double *const *entropyDPT, const double *const *enthalpyDPT, double RC, 
                        const unsigned char *numSeq1, const unsigned char *numSeq2)
{
   int loopSize1, loopSize2, loopSize;
   double T1, T2;
   double S,H;
   /* int N, N_loop; Triinu, please review */
   T1 = T2 = -_INFINITY;
   S = MinEntropy;
   H = 0.0;
   loopSize1 = ii - i - 1;
   loopSize2 = j - jj - 1;
   if (loopSize1 + loopSize2 > maxLoop) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      return;
   }
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
            H += enthalpyDPT[ii][jj]; /* bulge koos otsaga, st bulge i,j-ni */
            S += entropyDPT[ii][jj];
         }

         if(isFinite(H)) {
            T1 = (H) / (S + RC);
            T2 = (enthalpyDPT[i][j]) / ((entropyDPT[i][j]) + RC);

            if((T1 > T2) || ((traceback && T1 >= T2) || traceback==1)) {
               EntropyEnthalpy[0] = S;
               EntropyEnthalpy[1] = H;
            }
         }

      } else { /* we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30 */

         H = bulgeLoopEnthalpies[loopSize] + atpH[numSeq1[i]][numSeq2[j]] + atpH[numSeq1[ii]][numSeq2[jj]];
         if(traceback!=1)
           H += enthalpyDPT[ii][jj];

         S = bulgeLoopEntropies[loopSize] + atpS[numSeq1[i]][numSeq2[j]] + atpS[numSeq1[ii]][numSeq2[jj]];
         if(traceback!=1)
           S += entropyDPT[ii][jj];
         if(isFinite(H)) {
            T1 = (H) / (S + RC);
            T2 = (enthalpyDPT[i][j]) / (entropyDPT[i][j] + RC);
            if((T1 > T2) || ((traceback && T1 >= T2) || (traceback==1))) {
               EntropyEnthalpy[0] = S;
               EntropyEnthalpy[1] = H;
            }
         }
      }
   } /* end of calculating bulges */
   else if (loopSize1 == 1 && loopSize2 == 1) {
      /* mismatch nearest neighbor parameters */

      S = stackint2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
        stackint2Entropies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]];
      if(traceback!=1)
        S += entropyDPT[ii][jj];

      H = stackint2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
        stackint2Enthalpies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]];
      if(traceback!=1)
        H += enthalpyDPT[ii][jj];
      if(isFinite(H)) {
         T1 = (H) / (S + RC);
         T2 = (enthalpyDPT[i][j]) / (entropyDPT[i][j] + RC);
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
        H += enthalpyDPT[ii][jj];

      S = interiorLoopEntropies[loopSize] + tstackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
        tstackEntropies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]] + (ILAS * abs(loopSize1 - loopSize2));
      if(traceback!=1)
        S += entropyDPT[ii][jj];
      if(isFinite(H)) {
         T1 = (H) / (S + RC);
         T2 = (enthalpyDPT[i][j]) / ((entropyDPT[i][j]) + RC);
         if((T1 > T2) || ((traceback && T1 >= T2) || (traceback==1))) {
            EntropyEnthalpy[0] = S;
            EntropyEnthalpy[1] = H;
         }
      }
   }
   return;
}

static void 
calc_terminal_bp(double temp, const double *const *entropyDPT, const double *const *enthalpyDPT, double *send5, double *hend5, double RC,
                  const unsigned char *numSeq1, const unsigned char *numSeq2, int oligo1_len, int oligo2_len) { /* compute exterior loop */
   int i, k;
   send5[0] = send5[1] = -1.0;
   hend5[0] = hend5[1] = _INFINITY;
   for(i = 2; i<=(oligo1_len); i++) {
      send5[i] = MinEntropy;
      hend5[i] = 0;
   }

   double max_tm, S_max, H_max, H, S, T0, T1, G;
   for (i = 2; i <= oligo1_len; i++){
      max_tm = (hend5[i - 1]) / (send5[i - 1] + RC);
      S_max = send5[i-1];
      H_max = hend5[i-1];
      for(k = 0; k <= i - min_hrpn_loop - 2; ++k) {
         //END5_1
         T0 = (hend5[k]) /(send5[k] + RC);
         H = atpH[numSeq1[k + 1]][numSeq1[i]] + enthalpyDPT[k + 1][i];
         S = atpS[numSeq1[k + 1]][numSeq1[i]] + entropyDPT[k + 1][i];
         if(T0 >= 0.0) {
            H += hend5[k];
            S += send5[k];
         }
         if(H <= 0 && S <= 0) {
            T1 = (H) / (S + RC);
            G = H - temp*S;
            if((max_tm < T1) && (G<0.0)) {
               if(S > MinEntropyCutoff) {
                  H_max = H;
                  S_max = S;
                  max_tm = T1;
               }
            }
         }

         //END5_2
         H = atpH[numSeq1[k + 2]][numSeq1[i]] + dangleEnthalpies5[numSeq1[i]][numSeq1[k+2]][numSeq1[k+1]] + enthalpyDPT[k + 2][i];
         S = atpS[numSeq1[k + 2]][numSeq1[i]] + dangleEntropies5[numSeq1[i]][numSeq1[k+2]][numSeq1[k+1]] + entropyDPT[k + 2][i];
         if(T0 >= 0.0) {
            H += hend5[k];
            S += send5[k];
         }
         if(H <= 0 && S <= 0) {
            T1 = (H) / (S + RC);
            G = H - temp*S;
            if((max_tm < T1) && (G<0.0)) {
               if(S > MinEntropyCutoff) {
                  H_max = H;
                  S_max = S;
                  max_tm = T1;
               }
            }
         }

         //END5_3
         H = atpH[numSeq1[k + 1]][numSeq1[i - 1]] + dangleEnthalpies3[numSeq1[i-1]][numSeq1[i]][numSeq1[k+1]] + enthalpyDPT[k + 1][i - 1];
         S = atpS[numSeq1[k + 1]][numSeq1[i - 1]] + dangleEntropies3[numSeq1[i-1]][numSeq1[i]][numSeq1[k+1]] + entropyDPT[k + 1][i - 1];
         if(T0 >= 0.0) {
            H += hend5[k];
            S += send5[k];
         }
         if(H <= 0 && S <= 0) {
            T1 = (H) / (S + RC);
            G = H - temp*S;
            if((max_tm < T1) && (G<0.0)) {
               if(S > MinEntropyCutoff) {
                  H_max = H;
                  S_max = S;
                  max_tm = T1;
               }
            }
         }

         //END5_4
         H =atpH[numSeq1[k + 2]][numSeq1[i - 1]] + tstack2Enthalpies[numSeq1[i-1]][numSeq1[i]][numSeq1[k+2]][numSeq1[k+1]] + enthalpyDPT[k + 2][i - 1];
         S =atpS[numSeq1[k + 2]][numSeq1[i - 1]] + tstack2Entropies[numSeq1[i-1]][numSeq1[i]][numSeq1[k+2]][numSeq1[k+1]] + entropyDPT[k + 2][i - 1];
         if(T0 >= 0.0) {
            H += hend5[k];
            S += send5[k];
         }
         if(H <= 0 && S <= 0) {
            T1 = (H) / (S + RC);
            G = H - temp*S;
            if((max_tm < T1) && (G<0.0)) {
               if(S > MinEntropyCutoff) {
                  H_max = H;
                  S_max = S;
                  max_tm = T1;
               }
            }
         }
      }
      hend5[i] = H_max;
      send5[i] = S_max;
   }
}

static void 
push(struct tracer** stack, int i, int j, int mtrx, jmp_buf _jmp_buf, thal_results* o)
{
   struct tracer* new_top;
   new_top = (struct tracer*) safe_malloc(sizeof(struct tracer), _jmp_buf, o);
   new_top->i = i;
   new_top->j = j;
   new_top->mtrx = mtrx;
   new_top->next = *stack;
   *stack = new_top;
}

static void 
traceback_monomer(int* bp, int maxLoop, const double *const *entropyDPT,  const double *const *enthalpyDPT,
               double *send5, double *hend5, double RC, double dplx_init_S, double dplx_init_H,
               const unsigned char *numSeq1, const unsigned char *numSeq2, int oligo1_len,
               int oligo2_len, jmp_buf _jmp_buf, thal_results* o) /* traceback for unimolecular structure */
{
   int i, j;
   i = j = 0;
   int ii, jj, k;
   struct tracer *top, *stack = NULL;
   double SH1[2];
   double EntropyEnthalpy[2];
   double T1, H, S, max_tm;

   double T0 = -_INFINITY;
   push(&stack,oligo1_len, 0, 1, _jmp_buf, o);
   while(stack) {
      top = stack;
      stack = stack->next;
      i = top->i;
      j = top->j;
      if(top->mtrx==1) {
         while (equal(send5[i], send5[i - 1]) && equal(hend5[i], hend5[i - 1])) /* if previous structure is the same as this one */
           --i;
         if (i == 0)
           continue;
         //END5_1
         max_tm = (hend5[i - 1]) / (send5[i - 1] + RC);
         for (k = 0; k <= i - min_hrpn_loop - 2; ++k){
            T0 = (hend5[k]) /(send5[k] + RC);
            H = atpH[numSeq1[k + 1]][numSeq1[i]] + enthalpyDPT[k + 1][i];
            S = atpS[numSeq1[k + 1]][numSeq1[i]] + entropyDPT[k + 1][i];
            if(T0 >= 0.0) {
               H += hend5[k];
               S += send5[k];
            }
            if(H <= 0 && S <= 0) {
               T1 = (H) / (S + RC);
               if(max_tm < T1) {
                  if(S > MinEntropyCutoff) {
                     if (equal(send5[i], S) && equal(hend5[i], H)){
                        if (T0 >= 0.0){
                           push(&stack, k + 1, i, 0, _jmp_buf, o);
                           push(&stack, k, 0, 1, _jmp_buf, o);
                        }
                        else {
                           push(&stack, k + 1, i, 0, _jmp_buf, o);
                        }
                        break;
                     }
                  }
               }
            }
            //END5_2
            H = atpH[numSeq1[k + 2]][numSeq1[i]] + dangleEnthalpies5[numSeq1[i]][numSeq1[k+2]][numSeq1[k+1]] + enthalpyDPT[k + 2][i];
            S = atpS[numSeq1[k + 2]][numSeq1[i]] + dangleEntropies5[numSeq1[i]][numSeq1[k+2]][numSeq1[k+1]] + entropyDPT[k + 2][i];
            if(T0 >= 0.0) {
               H += hend5[k];
               S += send5[k];
            }
            if(H <= 0 && S <= 0) {
               T1 = (H) / (S + RC);
               if(max_tm < T1) {
                  if(S > MinEntropyCutoff) {
                     if (equal(send5[i], S) && equal(hend5[i], H)){
                        if (T0 >= 0.0){
                           push(&stack, k + 2, i, 0, _jmp_buf, o);
                           push(&stack, k, 0, 1, _jmp_buf, o);
                        }
                        else {
                           push(&stack, k + 2, i, 0, _jmp_buf, o);
                        }
                        break;
                     }
                  }
               }
            }
            //END5_3
            H = atpH[numSeq1[k + 1]][numSeq1[i - 1]] + dangleEnthalpies3[numSeq1[i-1]][numSeq1[i]][numSeq1[k+1]] + enthalpyDPT[k + 1][i - 1];
            S = atpS[numSeq1[k + 1]][numSeq1[i - 1]] + dangleEntropies3[numSeq1[i-1]][numSeq1[i]][numSeq1[k+1]] + entropyDPT[k + 1][i - 1];
            if(T0 >= 0.0) {
               H += hend5[k];
               S += send5[k];
            }
            if(H <= 0 && S <= 0) {
               T1 = (H) / (S + RC);
               if(max_tm < T1) {
                  if(S > MinEntropyCutoff) {
                     if (equal(send5[i], S) && equal(hend5[i], H)){
                        if (T0 >= 0.0){
                           push(&stack, k + 1, i - 1, 0, _jmp_buf, o);
                           push(&stack, k, 0, 1, _jmp_buf, o);
                        }
                        else {
                           push(&stack, k + 1, i - 1, 0, _jmp_buf, o);
                        }
                        break;
                     }
                  }
               }
            }
            //END5_4
            H =atpH[numSeq1[k + 2]][numSeq1[i - 1]] + tstack2Enthalpies[numSeq1[i-1]][numSeq1[i]][numSeq1[k+2]][numSeq1[k+1]] + enthalpyDPT[k + 2][i - 1];
            S =atpS[numSeq1[k + 2]][numSeq1[i - 1]] + tstack2Entropies[numSeq1[i-1]][numSeq1[i]][numSeq1[k+2]][numSeq1[k+1]] + entropyDPT[k + 2][i - 1];
            if(T0 >= 0.0) {
               H += hend5[k];
               S += send5[k];
            }
            if(H <= 0 && S <= 0) {
               T1 = (H) / (S + RC);
               if(max_tm < T1) {
                  if(S > MinEntropyCutoff) {
                     if (equal(send5[i], S) && equal(hend5[i], H)){
                        if (T0 >= 0.0){
                           push(&stack, k + 2, i - 1, 0, _jmp_buf, o);
                           push(&stack, k, 0, 1, _jmp_buf, o);
                        }
                        else {
                           push(&stack, k + 2, i - 1, 0, _jmp_buf, o);
                        }
                        break;
                     }
                  }
               }
            }
         }
      }
      if(top->mtrx==0) {
         bp[i - 1] = j;
         bp[j - 1] = i;
         SH1[0] = -1.0;
         SH1[1] = _INFINITY;
         calc_hairpin(i, j, SH1, 1, entropyDPT, enthalpyDPT, RC, numSeq1, numSeq2, oligo1_len, oligo2_len); /* 1 means that we use this method in traceback */
         if (equal(entropyDPT[i][j], stackEntropies[numSeq2[i]][numSeq2[i+1]][numSeq2[j]][numSeq2[j-1]] + entropyDPT[i + 1][j - 1]) &&
             equal(enthalpyDPT[i][j], stackEnthalpies[numSeq2[i]][numSeq2[i+1]][numSeq2[j]][numSeq2[j-1]] + enthalpyDPT[i + 1][j - 1])) {
            push(&stack, i + 1, j - 1, 0, _jmp_buf, o);
         }
         else if (!equal(entropyDPT[i][j], SH1[0]) || !equal(enthalpyDPT[i][j], SH1[1])) {
            int d, done;
            for (done = 0, d = j - i - 3; d >= min_hrpn_loop + 1 && d >= j - i - 2 - maxLoop && !done; --d)
              for (ii = i + 1; ii < j - d; ++ii) {
                 jj = d + ii;
                 EntropyEnthalpy[0] = -1.0;
                 EntropyEnthalpy[1] = _INFINITY;
                 calc_bulge_internal_monomer(i, j, ii, jj, EntropyEnthalpy, 1, maxLoop, entropyDPT, enthalpyDPT, RC, numSeq1, numSeq2);
                 if (equal(entropyDPT[i][j], EntropyEnthalpy[0] + entropyDPT[ii][jj]) &&
                     equal(enthalpyDPT[i][j], EntropyEnthalpy[1] + enthalpyDPT[ii][jj])) {
                    push(&stack, ii, jj, 0, _jmp_buf, o);
                    ++done;
                    break;
                 }
              }
         }
      }
      free(top);
   }
}

char * 
drawHairpin(int* bp, double mh, double ms, const thal_mode mode, double temp, const unsigned char *oligo1, const unsigned char *oligo2,
            double saltCorrection, int oligo1_len, int oligo2_len, jmp_buf _jmp_buf, thal_results *o)
{
   int  ret_space = 0;
   char *ret_ptr;
   int ret_last_l, ret_first_r, ret_center, ret_left_end, ret_right_start, ret_left_len, ret_right_len;
   int ret_add_sp_l, ret_add_sp_r;
   char ret_center_char;
   char ret_para[400];
   char* ret_str;
   /* Plain text */
   int i, N;
   ret_ptr = NULL;
   N = 0;
   double mg, t;
   if (!isFinite(ms) || !isFinite(mh)) {
      if((mode != THL_FAST) && (mode != THL_DEBUG_F)) {
        if (mode != THL_STRUCT) {
          printf("0\tdS = %g\tdH = %g\tinf\tinf\n", (double) ms,(double) mh);
#ifdef DEBUG
          fputs("No temperature could be calculated\n",stderr);
#endif
        }
      } else {
         o->temp = 0.0; /* lets use generalization here */
         strcpy(o->msg, "No predicted sec struc for given seq\n");
      }
   } else {
      if((mode != THL_FAST) && (mode != THL_DEBUG_F)) {
         for (i = 1; i < oligo1_len; ++i) {
            if(bp[i-1] > 0) N++;
         }
      } else {
         for (i = 1; i < oligo1_len; ++i) {
            if(bp[i-1] > 0) N++;
         }
      }
      t = (mh / (ms + (((N/2)-1) * saltCorrection))) - ABSOLUTE_ZERO;
      mg = mh - (temp * (ms + (((N/2)-1) * saltCorrection)));
      ms = ms + (((N/2)-1) * saltCorrection);
      o->dg = mg;
      o->ds = ms;
      o->dh = mh;
      o->temp = (double) t;
      if((mode != THL_FAST) && (mode != THL_DEBUG_F)) {
         if (mode != THL_STRUCT) {
           printf("Calculated thermodynamical parameters for dimer:\t%d\tdS = %g\tdH = %g\tdG = %g\tt = %g\n",
                  oligo1_len, (double) ms, (double) mh, (double) mg, (double) t);
         } else {
           snprintf(ret_para, 400, "Tm: %.1f&deg;C  dG: %.0f cal/mol  dH: %.0f cal/mol  dS: %.0f cal/mol*K\\n",
                   (double) t, (double) mg, (double) mh, (double) ms);
         }
      } else {
         return NULL;
      }
   }
   /* plain-text output */
   char* asciiRow;
   asciiRow = (char*) safe_malloc(oligo1_len, _jmp_buf, o);
   for(i = 0; i < oligo1_len; ++i) asciiRow[i] = '0';
   for(i = 1; i < oligo1_len+1; ++i) {
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
   if ((mode == THL_GENERAL) || (mode == THL_DEBUG)) {
     printf("SEQ\t");
     for(i = 0; i < oligo1_len; ++i) printf("%c",asciiRow[i]);
     printf("\nSTR\t%s\n", oligo1);
   }
   if (mode == THL_STRUCT) {
     ret_str = NULL;

     save_append_string(&ret_str, &ret_space, o, ret_para, _jmp_buf);

     ret_last_l = -1;
     ret_first_r = -1;
     ret_center_char = '|';
     for(i = 0; i < oligo1_len; ++i) {
       if (asciiRow[i] == '/') {
         ret_last_l = i;
       }
       if ((ret_first_r == -1) && (asciiRow[i] == '\\')) {
         ret_first_r = i;
       }
     }
     ret_center = ret_first_r - ret_last_l;
     if (ret_center % 2 == 0) { 
       /* ret_center is odd */
       ret_left_end = ret_last_l + (ret_first_r - ret_last_l) / 2 - 1;
       ret_center_char = (char) oligo1[ret_left_end + 1]; 
       ret_right_start = ret_left_end + 2;
     } else {
       /* ret_center is even */
       ret_left_end = ret_last_l + (ret_first_r - ret_last_l - 1) / 2;
       ret_right_start = ret_left_end + 1;
     }
     ret_left_len = ret_left_end + 1;
     ret_right_len = oligo1_len - ret_right_start;
     ret_add_sp_l = 0;
     ret_add_sp_r = 0;
     if (ret_left_len > ret_right_len) {
       ret_add_sp_r = ret_left_len - ret_right_len + 1;
     }
     if (ret_right_len > ret_left_len) {
       ret_add_sp_l = ret_right_len - ret_left_len;
     }
     for (i = 0 ; i < ret_add_sp_l ; i++) {
       save_append_char(&ret_str, &ret_space, o, ' ', _jmp_buf);
     }
     save_append_string(&ret_str, &ret_space, o, "5' ", _jmp_buf);
     for (i = 0 ; i < ret_left_len ; i++) {
       save_append_char(&ret_str, &ret_space, o, (char) oligo1[i], _jmp_buf);
     }
     save_append_string(&ret_str, &ret_space, o, "U+2510\\n   ", _jmp_buf);
     for (i = 0 ; i < ret_add_sp_l ; i++) {
       save_append_char(&ret_str, &ret_space, o, ' ', _jmp_buf);
     }
     for (i = 0 ; i < ret_left_len ; i++) {
       if (asciiRow[i] == '/') {
         save_append_char(&ret_str, &ret_space, o, '|', _jmp_buf);
       } else {
         save_append_char(&ret_str, &ret_space, o, ' ', _jmp_buf);
       }
     }
     if (ret_center_char == '|' ) {
       save_append_string(&ret_str, &ret_space, o, "U+2502", _jmp_buf);
     } else {
       save_append_char(&ret_str, &ret_space, o, ret_center_char, _jmp_buf);
     }
     save_append_string(&ret_str, &ret_space, o, "\\n", _jmp_buf);
     for (i = 0 ; i < ret_add_sp_r - 1 ; i++) {
       save_append_char(&ret_str, &ret_space, o, ' ', _jmp_buf);
     }
     save_append_string(&ret_str, &ret_space, o, "3' ", _jmp_buf);
     for (i = oligo1_len ; i > ret_right_start - 1; i--) {
       save_append_char(&ret_str, &ret_space, o, (char) oligo1[i], _jmp_buf);
     }
     save_append_string(&ret_str, &ret_space, o, "U+2518\\n", _jmp_buf);

     ret_ptr = (char *) safe_malloc(strlen(ret_str) + 1, _jmp_buf, o);
     strcpy(ret_ptr, ret_str);
     if (ret_str != NULL) {
       free(ret_str);
     }
   }
   free(asciiRow);
   return ret_ptr;
}

//=====================================================================================
//Misc helper functions
//=====================================================================================

static int 
comp3loop(const void* loop1, const void* loop2)
{

     int i;
     const unsigned char* h1 = (const unsigned char*) loop1;
     const struct triloop *h2 = (const struct triloop*) loop2;

     for (i = 0; i < 5; ++i)
         if (h1[i] < h2->loop[i])
             return -1;
       else if (h1[i] > h2->loop[i])
           return 1;

     return 0;
}

static int 
comp4loop(const void* loop1, const void* loop2)
{
   int i;
   const unsigned char* h1 = (const unsigned char*) loop1;
   const struct tetraloop *h2 = (const struct tetraloop*) loop2;

   for (i = 0; i < 6; ++i)
     if (h1[i] < h2->loop[i])
       return -1;
   else if (h1[i] > h2->loop[i])
     return 1;

   return 0;
}

static int 
equal(double a, double b)
{
#ifdef INTEGER
   return a == b;
#endif

   if (!isfinite(a) || !isfinite(b))
     return 0;
   return fabs(a - b) < 1e-5;

   if (a == 0 && b == 0)
     return 1;
}

//=====================================================================================
//Initializing functions
//=====================================================================================

static double 
saltCorrectS (double mv, double dv, double dntp)
{
   if(dv<=0) dntp=dv;
   return 0.368*((log((mv+120*(sqrt(fmax(0.0, dv-dntp))))/1000)));
}

int thal_check_errors(const unsigned char *oligo_f, const unsigned char *oligo_r, int *len_f, int *len_r, const thal_args *a, thal_results *o){
   if (oligo_f == NULL){
      strcpy(o->msg, "NULL first sequence");
      return 1;
   }
   if (oligo_r == NULL){
      strcpy(o->msg, "NULL second sequence");
      return 1;
   }
   *len_f = length_unsig_char(oligo_f);
   *len_r = length_unsig_char(oligo_r);

   /* The following error messages will be seen by end users and will
      not be easy to understand. */
   if((*len_f > THAL_MAX_ALIGN) && (*len_r > THAL_MAX_ALIGN)){
      strcpy(o->msg, "Both sequences longer than " XSTR(THAL_MAX_ALIGN)
         " for thermodynamic alignment");
      return 1;
   }
   if((*len_f > THAL_MAX_SEQ)){ 
      strcpy(o->msg, LONG_SEQ_ERR_STR(THAL_MAX_SEQ) " (1)");
      return 1;
   }
   if((*len_r > THAL_MAX_SEQ)){ 
      strcpy(o->msg, LONG_SEQ_ERR_STR(THAL_MAX_SEQ) " (2)");
      return 1;
   }

   if(NULL == a){
      strcpy(o->msg, "NULL 'in' pointer");
      return 1;
   }
   if (NULL == o) return 1; /* Leave it to the caller to crash */
   if((a->type != thal_any) && (a->type != thal_end1) && (a->type != thal_end2) && (a->type != thal_hairpin)){
      strcpy(o->msg, "Illegal type");
      return 1;
   }
   o->align_end_1 = -1;
   o->align_end_2 = -1;
   if (oligo_f && '\0' == *oligo_f) {
      strcpy(o->msg, "Empty first sequence");
      o->temp = 0.0;
      return 1;
   }
   if (oligo_r && '\0' == *oligo_r) {
      strcpy(o->msg, "Empty second sequence");
      o->temp = 0.0;
      return 1;
   }
   if (0 == *len_f) {
      o->temp = 0.0;
      return 1;
   }
   if (0 == *len_r) {
      o->temp = 0.0;
      return 1;
   }
   return 0;
}

/* Set default args */
void 
set_thal_default_args(thal_args *a)
{
   memset(a, 0, sizeof(*a));
   a->type = thal_any; /* thal_alignment_type THAL_ANY */
   a->maxLoop = MAX_LOOP;
   a->mv = 50; /* mM */
   a->dv = 0.0; /* mM */
   a->dntp = 0.8; /* mM */
   a->dna_conc = 50; /* nM */
   a->temp = TEMP_KELVIN; /* Kelvin */
   a->dimer = 1; /* by default dimer structure is calculated */
}

/* Set default args for oligo */
void
set_thal_oligo_default_args(thal_args *a)    
{
   memset(a, 0, sizeof(*a));
   a->type = thal_any; /* thal_alignment_type THAL_ANY */
   a->maxLoop = MAX_LOOP;
   a->mv = 50; /* mM */
   a->dv = 0.0; /* mM */
   a->dntp = 0.0; /* mM */
   a->dna_conc = 50; /* nM */
   a->temp = TEMP_KELVIN; /* Kelvin */
   a->dimer = 1; /* by default dimer structure is calculated */
}

/* Return 1 if string is symmetrical, 0 otherwise. */
static int 
symmetry_thermo(const unsigned char* seq)
{
   char s;
   char e;
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
      s=toupper(*seq);
      e=toupper(*seq_end);
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

//=====================================================================================
//Functions for allocating memory
//=====================================================================================

static void* 
safe_calloc(size_t m, size_t n, jmp_buf _jmp_buf, thal_results *o)
{
   void* ptr;
   if (!(ptr = calloc(m, n))) {
#ifdef DEBUG
      fputs("Error in calloc()\n", stderr);
#endif
   strcpy(o->msg, "Out of memory");
   errno = ENOMEM;
   longjmp(_jmp_buf, 1);
   }
   return ptr;
}

static void* 
safe_malloc(size_t n, jmp_buf _jmp_buf, thal_results *o)
{
   void* ptr;
   if (!(ptr = malloc(n))) {
#ifdef DEBUG
      fputs("Error in malloc()\n", stderr);
#endif
   strcpy(o->msg, "Out of memory");
   errno = ENOMEM;
   longjmp(_jmp_buf, 1);
   }
   return ptr;
}

static void* 
safe_realloc(void* ptr, size_t n, jmp_buf _jmp_buf, thal_results *o)
{
   ptr = realloc(ptr, n);
   if (ptr == NULL) {
#ifdef DEBUG
      fputs("Error in realloc()\n", stderr);
#endif
   strcpy(o->msg, "Out of memory");
   errno = ENOMEM;
   longjmp(_jmp_buf, 1);
   }
   return ptr;
}

double **allocate_DPT(int oligo1_len, int oligo2_len, jmp_buf _jmp_buf, thal_results *o){
   //Add one to each dimension due to the way the DPT is indexed
   //Row i=0 and column j=0 are never used, but the wasted memory is negligible
   double *dpt = (double *) safe_malloc(sizeof(double) * (oligo1_len + 1) * (oligo2_len +1), _jmp_buf, o);
   double **rows = (double **) safe_malloc(sizeof(double *) * (oligo1_len + 1), _jmp_buf, o);
   for(int i = 0; i < oligo1_len+1; i++){
      rows[i] = &dpt[i * (oligo2_len +1)];
   }
   return rows;
}

void free_DPT(double **dpt){
   free(dpt[0]);
   free(dpt);
}

//=====================================================================================
//Functions for string manipulation
//=====================================================================================

static unsigned char 
str2int(char c)
{
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

static char*
th_read_str_line(char **str, jmp_buf _jmp_buf, thal_results* o)
{
  if (*str == NULL) {
    return NULL;
  }
  char *ptr = *str;
  char *ini = *str;
  while(1) {
    if ((*ptr == '\n') || (*ptr == '\0')) {
      char *ret = NULL;
      if (!(ret = (char *) malloc(sizeof(char) * (ptr - ini + 1)))) {
#ifdef DEBUG
        fputs("Error in malloc()\n", stderr);
#endif
         strcpy(o->msg, "Out of memory");
         errno = ENOMEM;
         longjmp(_jmp_buf, 1);
      }
      /* copy line */
      strncpy(ret, ini, (ptr - ini + 1));
      ret[ptr - ini] = '\0';

      if (*ptr == '\0') { /* End of String */
        *str = NULL;
      } else {
        ptr++;
        if (*ptr == '\0') { /* End of String */
          *str = NULL;
        } else {
          *str = ptr;
        }
      }
      if (ptr == ini) {
        if (ret != NULL) {
          free(ret);
        }
        return NULL;
      } else {  
        return ret;
      }
    }
    ptr++;
  }
}

static void 
reverse(unsigned char *s)
{
   int i,j;
   char c;
   for (i = 0, j = length_unsig_char(s)-1; i < j; i++, j--) {
      c = s[i];
      s[i] = s[j];
      s[j] = c;
   }
}

/* These functions are needed as "inf" cannot be read on Windows directly */
static double 
readDouble(char **str, jmp_buf _jmp_buf, thal_results* o)
{
  double result;
  char *line = th_read_str_line(str, _jmp_buf, o);
  /* skip any spaces at beginning of the line */
  while (isspace(*line)) line++;
  if (!strncmp(line, "inf", 3)) {
    free(line);
    return _INFINITY;
  }
  sscanf(line, "%lf", &result);
  if (line != NULL) {
    free(line);
  }
  return result;
}
static void 
strcatc(char* str, char c)
{
   str[strlen(str) + 1] = 0;
   str[strlen(str)] = c;
}

static void
save_append_string(char** ret, int *space, thal_results *o, const char *str, jmp_buf _jmp_buf) {
  int xlen, slen;
  if (str == NULL) {
    return;
  }
  if (*ret == NULL) {
    *ret = (char *) safe_malloc(sizeof(char)*500, _jmp_buf, o);
    *ret[0] = '\0';
    *space = 500;
  }
  xlen = strlen(*ret);
  slen = strlen(str);
  if (xlen + slen + 1 > *space) {
    *space += 4 * (slen + 1);
    *ret = (char *) safe_realloc(*ret, *space, _jmp_buf, o);
  }
  strcpy(*ret + xlen, str);
  return;
}

static void
save_append_char(char** ret, int *space, thal_results *o, const char str, jmp_buf _jmp_buf) {
  char fix[3];
  fix[0] = str;
  fix[1] = '\0';
  save_append_string(ret, space, o, fix, _jmp_buf);
}

static int 
length_unsig_char(const unsigned char * str)
{
   int i = 0;
   while(*(str++)) {
      i++;
      if(i == INT_MAX)
        return -1;
   }
   return i;
}


//=====================================================================================
//Functions for reading and parsing thermodynamic parameters
//=====================================================================================

/* Initialize the thermodynamic values (parameters) */
int  thal_set_null_parameters(thal_parameters *a) {
  a->dangle_dh = NULL;
  a->dangle_ds = NULL;
  a->loops_dh = NULL;
  a->loops_ds = NULL;
  a->stack_dh = NULL;
  a->stack_ds = NULL;
  a->stackmm_dh = NULL;
  a->stackmm_ds = NULL;
  a->tetraloop_dh = NULL;
  a->tetraloop_ds = NULL;
  a->triloop_dh = NULL;
  a->triloop_ds = NULL;
  a->tstack_tm_inf_ds = NULL;
  a->tstack_dh = NULL;
  a->tstack2_dh = NULL;
  a->tstack2_ds = NULL;
  return 0;
}

/* Free the thermodynamic values (parameters) */
int  thal_free_parameters(thal_parameters *a) {
  if (NULL != a->dangle_dh) {
    free(a->dangle_dh);
    a->dangle_dh = NULL;
  }
  if (NULL != a->dangle_ds) {
    free(a->dangle_ds);
    a->dangle_ds = NULL;
  }
  if (NULL != a->loops_dh) {
    free(a->loops_dh);
    a->loops_dh = NULL;
  }
  if (NULL != a->loops_ds) {
    free(a->loops_ds);
    a->loops_ds = NULL;
  }
  if (NULL != a->stack_dh) {
    free(a->stack_dh);
    a->stack_dh = NULL;
  }
  if (NULL != a->stack_ds) {
    free(a->stack_ds);
    a->stack_ds = NULL;
  }
  if (NULL != a->stackmm_dh) {
    free(a->stackmm_dh);
    a->stackmm_dh = NULL;
  }
  if (NULL != a->stackmm_ds) {
    free(a->stackmm_ds);
    a->stackmm_ds = NULL;
  }
  if (NULL != a->tetraloop_dh) {
    free(a->tetraloop_dh);
    a->tetraloop_dh = NULL;
  }
  if (NULL != a->tetraloop_ds) {
    free(a->tetraloop_ds);
    a->tetraloop_ds = NULL;
  }
  if (NULL != a->triloop_dh) {
    free(a->triloop_dh);
    a->triloop_dh = NULL;
  }
  if (NULL != a->triloop_ds) {
    free(a->triloop_ds);
    a->triloop_ds = NULL;
  }
  if (NULL != a->tstack_tm_inf_ds) {
    free(a->tstack_tm_inf_ds);
    a->tstack_tm_inf_ds = NULL;
  }
  if (NULL != a->tstack_dh) {
    free(a->tstack_dh);
    a->tstack_dh = NULL;
  }
  if (NULL != a->tstack2_dh) {
    free(a->tstack2_dh);
    a->tstack2_dh = NULL;
  }
  if (NULL != a->tstack2_ds) {
    free(a->tstack2_ds);
    a->tstack2_ds = NULL;
  }
  return 0;
}

/* Read the thermodynamic values (parameters) from the parameter files
   in the directory specified by 'path'.  Return 0 on success and -1
   on error. The thermodynamic values are stored in multiple static
   variables. */
int 
get_thermodynamic_values(const thal_parameters *tp, thal_results *o)
{
   jmp_buf _jmp_buf;
  if (setjmp(_jmp_buf) != 0) {
     return -1;
  }
  getStack(stackEntropies, stackEnthalpies, tp, _jmp_buf, o);
  /* verifyStackTable(stackEntropies, "entropy");
     verifyStackTable(stackEnthalpies, "enthalpy"); */ /* this is for code debugging */
  getStackint2(stackint2Entropies, stackint2Enthalpies, tp, _jmp_buf, o);
  getDangle(dangleEntropies3, dangleEnthalpies3, dangleEntropies5, dangleEnthalpies5, tp, _jmp_buf, o);
  getLoop(hairpinLoopEntropies, interiorLoopEntropies, bulgeLoopEntropies, hairpinLoopEnthalpies,
          interiorLoopEnthalpies, bulgeLoopEnthalpies, tp, _jmp_buf, o);
  getTstack(tstackEntropies, tstackEnthalpies, tp, _jmp_buf, o);
  getTstack2(tstack2Entropies, tstack2Enthalpies, tp, _jmp_buf, o);
  getTriloop(&triloopEntropies, &triloopEnthalpies, &numTriloops, tp, _jmp_buf, o);
  getTetraloop(&tetraloopEntropies, &tetraloopEnthalpies, &numTetraloops, tp, _jmp_buf, o);
  /* getting the AT-penalties */
  tableStartATS(AT_S, atpS);
  tableStartATH(AT_H, atpH);

  return 0;
}

void 
destroy_thal_structures()
{
  if ((triloopEntropies != NULL) && (triloopEntropies != defaultTriloopEntropies)){
    free(triloopEntropies);
    triloopEntropies = NULL;
  }
  if ((triloopEnthalpies != NULL) && (triloopEnthalpies != defaultTriloopEnthalpies)){
    free(triloopEnthalpies);
    triloopEnthalpies = NULL;
  }
  if ((tetraloopEntropies != NULL) && (tetraloopEntropies != defaultTetraloopEntropies)){
    free(tetraloopEntropies);
    tetraloopEntropies = NULL;
  }
  if ((tetraloopEnthalpies != NULL) && (tetraloopEnthalpies != defaultTetraloopEnthalpies)){
    free(tetraloopEnthalpies);
    tetraloopEnthalpies = NULL;
  }
}


static char* 
readParamFile(const char* dirname, const char* fname, jmp_buf _jmp_buf, thal_results* o)
{
  FILE* file;
  char* ret = NULL;
  char* paramdir = NULL;
  paramdir = (char*) safe_malloc(strlen(dirname) + strlen(fname) + 2, _jmp_buf, o);
  strcpy(paramdir, dirname);
#ifdef OS_WIN
  if (paramdir[strlen(paramdir) - 1] != '\\') {
    strcat(paramdir, "\\\0");
  }
#else
  if (paramdir[strlen(paramdir) - 1] != '/') {
    strcat(paramdir, "/\0");
  }
#endif
  strcat(paramdir, fname);
  if (!(file = fopen(paramdir, "r"))) {
    snprintf(o->msg, 255, "Unable to open file %s", paramdir);
    if (paramdir != NULL) {
      free(paramdir);
      paramdir = NULL;
    }
    longjmp(_jmp_buf, 1);
    return NULL;
  }
  if (paramdir != NULL) {
    free(paramdir);
    paramdir = NULL;
  }
  char c;
  int i = 0;
  size_t ssz = INIT_BUF_SIZE;
  size_t remaining_size;
  remaining_size = ssz;
  ret = (char*) safe_malloc(ssz, _jmp_buf, o);
  while (1) {
    if (feof(file)) {
      ret[i] = '\0';
      fclose(file);
      return ret;
    }
    c = fgetc(file);
    remaining_size -= sizeof(char);
    if (remaining_size <= 0) {
      if (ssz >= INT_MAX / 2) {
        strcpy(o->msg, "Out of memory");
        free(ret);
        longjmp(_jmp_buf, 1);
        return NULL;
      } else {
        ssz += INIT_BUF_SIZE;
        remaining_size += INIT_BUF_SIZE;
      }
      ret = (char *) safe_realloc(ret, ssz, _jmp_buf, o);
    }
    ret[i] = c;
    i++;
  }
}

int
thal_load_parameters(const char *path, thal_parameters *a, thal_results* o)
{
   jmp_buf _jmp_buf;
  thal_free_parameters(a);
  if (setjmp(_jmp_buf) != 0) {
    printf("longjump\n");
    return -1;
  }
  a->dangle_dh = readParamFile(path, "dangle.dh", _jmp_buf, o);
  a->dangle_ds = readParamFile(path, "dangle.ds", _jmp_buf, o);
  a->loops_dh = readParamFile(path, "loops.dh", _jmp_buf, o);
  a->loops_ds = readParamFile(path, "loops.ds", _jmp_buf, o);
  a->stack_dh = readParamFile(path, "stack.dh", _jmp_buf, o);
  a->stack_ds = readParamFile(path, "stack.ds", _jmp_buf, o);
  a->stackmm_dh = readParamFile(path, "stackmm.dh", _jmp_buf, o);
  a->stackmm_ds = readParamFile(path, "stackmm.ds", _jmp_buf, o);
  a->tetraloop_dh = readParamFile(path, "tetraloop.dh", _jmp_buf, o);
  a->tetraloop_ds = readParamFile(path, "tetraloop.ds", _jmp_buf, o);
  a->triloop_dh = readParamFile(path, "triloop.dh", _jmp_buf, o);
  a->triloop_ds = readParamFile(path, "triloop.ds", _jmp_buf, o);
  a->tstack_tm_inf_ds = readParamFile(path, "tstack_tm_inf.ds", _jmp_buf, o);
  a->tstack_dh = readParamFile(path, "tstack.dh", _jmp_buf, o);
  a->tstack2_dh = readParamFile(path, "tstack2.dh", _jmp_buf, o);
  a->tstack2_ds = readParamFile(path, "tstack2.ds", _jmp_buf, o);
  return 0;
}

/* Reads a line containing 4 doubles, which can be specified as "inf". */
static void
readLoop(char **str, double *v1, double *v2, double *v3, jmp_buf _jmp_buf, thal_results *o)
{
  char *line = th_read_str_line(str, _jmp_buf, o);
  char *p = line, *q;
  /* skip first number on the line */
  while (isspace(*p)) p++;
  while (isdigit(*p)) p++;
  while (isspace(*p)) p++;
  /* read second number */
  q = p;
  while (!isspace(*q)) q++;
  *q = '\0'; q++;
  if (!strcmp(p, "inf")) *v1 = _INFINITY;
  else sscanf(p, "%lf", v1);
  while (isspace(*q)) q++;
  /* read third number */
  p = q;
  while (!isspace(*p)) p++;
  *p = '\0'; p++;
  if (!strcmp(q, "inf")) *v2 = _INFINITY;
  else sscanf(q, "%lf", v2);
  while (isspace(*p)) p++;
  /* read last number */
  q = p;
  while (!isspace(*q) && (*q != '\0')) q++;
  *q = '\0';
  if (!strcmp(p, "inf")) *v3 = _INFINITY;
  else sscanf(p, "%lf", v3);
  if (line != NULL) {
    free(line);
  }
}

/* Reads a line containing a short string and a double, used for reading a triloop or tetraloop. */
static int
readTLoop(char **str, char *s, double *v, int triloop, jmp_buf _jmp_buf, thal_results *o)
{
  char *line = th_read_str_line(str, _jmp_buf, o);
  if (!line) return -1;
  char *p = line, *q;
  /* skip first spaces */
  while (isspace(*p)) p++;
  /* read the string */
  q = p;
  while (isalpha(*q)) q++;
  *q = '\0'; q++;
  if (triloop) {
    strncpy(s, p, 5);   /*triloop string has 5 characters*/
  } else {
    strncpy(s, p, 6);   /*tetraloop string has 6 characters*/
  }
  /* skip all spaces */
  while (isspace(*q)) q++;
  p = q;
  while (!isspace(*p) && (*p != '\0')) p++;
  *p = '\0';
  if (!strcmp(q, "inf")) *v = _INFINITY;
  else sscanf(q, "%lg", v);
  if (line != NULL) {
    free(line);
  }
  return 0;
} 

static void 
getStack(double stackEntropies[5][5][5][5], double stackEnthalpies[5][5][5][5], const thal_parameters *tp, jmp_buf _jmp_buf, thal_results* o)
{
   int i, j, ii, jj;
   char *pt_ds = tp->stack_ds;
   char *pt_dh = tp->stack_dh;
   for (i = 0; i < 5; ++i) {
      for (ii = 0; ii < 5; ++ii) {
         for (j = 0; j < 5; ++j) {
            for (jj = 0; jj < 5; ++jj) {
               if (i == 4 || j == 4 || ii == 4 || jj == 4) {
                  stackEntropies[i][ii][j][jj] = -1.0;
                  stackEnthalpies[i][ii][j][jj] = _INFINITY;
               } else {
                  stackEntropies[i][ii][j][jj] = readDouble(&pt_ds, _jmp_buf, o);
                  stackEnthalpies[i][ii][j][jj] = readDouble(&pt_dh, _jmp_buf, o);
                  if (!isFinite(stackEntropies[i][ii][j][jj]) || !isFinite(stackEnthalpies[i][ii][j][jj])) {
                     stackEntropies[i][ii][j][jj] = -1.0;
                     stackEnthalpies[i][ii][j][jj] = _INFINITY;
                  }
               }
            }
         }
      }
   }
}

static void 
getStackint2(double stackint2Entropies[5][5][5][5], double stackint2Enthalpies[5][5][5][5], const thal_parameters *tp, jmp_buf _jmp_buf, thal_results* o)
{
   int i, j, ii, jj;
   char *pt_ds = tp->stackmm_ds;
   char *pt_dh = tp->stackmm_dh;
   for (i = 0; i < 5; ++i) {
      for (ii = 0; ii < 5; ++ii) {
         for (j = 0; j < 5; ++j) {
            for (jj = 0; jj < 5; ++jj) {
               if (i == 4 || j == 4 || ii == 4 || jj == 4) {
                  stackint2Entropies[i][ii][j][jj] = -1.0;
                  stackint2Enthalpies[i][ii][j][jj] = _INFINITY;
               } else {
                  stackint2Entropies[i][ii][j][jj] = readDouble(&pt_ds, _jmp_buf, o);
                  stackint2Enthalpies[i][ii][j][jj] = readDouble(&pt_dh, _jmp_buf, o);
                  if (!isFinite(stackint2Entropies[i][ii][j][jj]) || !isFinite(stackint2Enthalpies[i][ii][j][jj])) {
                     stackint2Entropies[i][ii][j][jj] = -1.0;
                     stackint2Enthalpies[i][ii][j][jj] = _INFINITY;
                  }
               }
            }
         }
      }
   }
}

static void 
getDangle(double dangleEntropies3[5][5][5], double dangleEnthalpies3[5][5][5], double dangleEntropies5[5][5][5],
          double dangleEnthalpies5[5][5][5], const thal_parameters *tp, jmp_buf _jmp_buf, thal_results* o)
{
   int i, j, k;
   char *pt_ds = tp->dangle_ds;
   char *pt_dh = tp->dangle_dh;
   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       for (k = 0; k < 5; ++k) {
          if (i == 4 || j == 4) {
             dangleEntropies3[i][k][j] = -1.0;
             dangleEnthalpies3[i][k][j] = _INFINITY;
          } else if (k == 4) {
             dangleEntropies3[i][k][j] = -1.0;
             dangleEnthalpies3[i][k][j] = _INFINITY;
          } else {
             dangleEntropies3[i][k][j] = readDouble(&pt_ds, _jmp_buf, o);
             dangleEnthalpies3[i][k][j] = readDouble(&pt_dh, _jmp_buf, o);
             if(!isFinite(dangleEntropies3[i][k][j]) || !isFinite(dangleEnthalpies3[i][k][j])) {
                dangleEntropies3[i][k][j] = -1.0;
                dangleEnthalpies3[i][k][j] = _INFINITY;             
             }
          }
       }

   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       for (k = 0; k < 5; ++k) {
          if (i == 4 || j == 4) {
             dangleEntropies5[i][j][k] = -1.0;
             dangleEnthalpies5[i][j][k] = _INFINITY;
          } else if (k == 4) {
             dangleEntropies5[i][j][k] = -1.0;
             dangleEnthalpies5[i][j][k] = _INFINITY;
          } else {
             dangleEntropies5[i][j][k] = readDouble(&pt_ds, _jmp_buf, o);
             dangleEnthalpies5[i][j][k] = readDouble(&pt_dh, _jmp_buf, o);
             if(!isFinite(dangleEntropies5[i][j][k]) || !isFinite(dangleEnthalpies5[i][j][k])) {
                dangleEntropies5[i][j][k] = -1.0;
                dangleEnthalpies5[i][j][k] = _INFINITY;
             }
          }
       }
}

static void 
getLoop(double hairpinLoopEntropies[30], double interiorLoopEntropies[30], double bulgeLoopEntropies[30],
        double hairpinLoopEnthalpies[30], double interiorLoopEnthalpies[30], double bulgeLoopEnthalpies[30], 
        const thal_parameters *tp, jmp_buf _jmp_buf, thal_results* o)
{
   int k;
   char *pt_ds = tp->loops_ds;
   char *pt_dh = tp->loops_dh;
   for (k = 0; k < 30; ++k) {
      readLoop(&pt_ds, &interiorLoopEntropies[k], &bulgeLoopEntropies[k], &hairpinLoopEntropies[k], _jmp_buf, o);
      readLoop(&pt_dh, &interiorLoopEnthalpies[k], &bulgeLoopEnthalpies[k], &hairpinLoopEnthalpies[k], _jmp_buf, o);
   }
}

static void 
getTstack(double tstackEntropies[5][5][5][5], double tstackEnthalpies[5][5][5][5], const thal_parameters *tp, jmp_buf _jmp_buf, thal_results* o)
{
   int i1, j1, i2, j2;
   char *pt_ds = tp->tstack_tm_inf_ds;
   char *pt_dh = tp->tstack_dh;
   for (i1 = 0; i1 < 5; ++i1)
     for (i2 = 0; i2 < 5; ++i2)
       for (j1 = 0; j1 < 5; ++j1)
         for (j2 = 0; j2 < 5; ++j2)
           if (i1 == 4 || j1 == 4) {
              tstackEnthalpies[i1][i2][j1][j2] = _INFINITY;
              tstackEntropies[i1][i2][j1][j2] = -1.0;
           } else if (i2 == 4 || j2 == 4) {
              tstackEntropies[i1][i2][j1][j2] = 0.00000000001;
              tstackEnthalpies[i1][i2][j1][j2] = 0.0;
           } else {
              tstackEntropies[i1][i2][j1][j2] = readDouble(&pt_ds, _jmp_buf, o);
              tstackEnthalpies[i1][i2][j1][j2] = readDouble(&pt_dh, _jmp_buf, o);
              if (!isFinite(tstackEntropies[i1][i2][j1][j2]) || !isFinite(tstackEnthalpies[i1][i2][j1][j2])) {
                 tstackEntropies[i1][i2][j1][j2] = -1.0;
                 tstackEnthalpies[i1][i2][j1][j2] = _INFINITY;
              }
           }
}

static void 
getTstack2(double tstack2Entropies[5][5][5][5], double tstack2Enthalpies[5][5][5][5], const thal_parameters *tp, jmp_buf _jmp_buf, thal_results* o)
{

   int i1, j1, i2, j2;
   char *pt_ds = tp->tstack2_ds;
   char *pt_dh = tp->tstack2_dh;
   for (i1 = 0; i1 < 5; ++i1)
     for (i2 = 0; i2 < 5; ++i2)
       for (j1 = 0; j1 < 5; ++j1)
         for (j2 = 0; j2 < 5; ++j2)
           if (i1 == 4 || j1 == 4)  {
              tstack2Enthalpies[i1][i2][j1][j2] = _INFINITY;
              tstack2Entropies[i1][i2][j1][j2] = -1.0;
           } else if (i2 == 4 || j2 == 4) {
              tstack2Entropies[i1][i2][j1][j2] = 0.00000000001;
              tstack2Enthalpies[i1][i2][j1][j2] = 0.0;
           } else {
              tstack2Entropies[i1][i2][j1][j2] = readDouble(&pt_ds, _jmp_buf, o);
              tstack2Enthalpies[i1][i2][j1][j2] = readDouble(&pt_dh, _jmp_buf, o);
              if (!isFinite(tstack2Entropies[i1][i2][j1][j2]) || !isFinite(tstack2Enthalpies[i1][i2][j1][j2])) {
                 tstack2Entropies[i1][i2][j1][j2] = -1.0;
                 tstack2Enthalpies[i1][i2][j1][j2] = _INFINITY;
              }
           }
}

static void 
getTriloop(struct triloop** triloopEntropies, struct triloop** triloopEnthalpies, int* num, const thal_parameters *tp, jmp_buf _jmp_buf, thal_results* o)
{
   int i, size;
   double value;
   char *pt_ds = tp->triloop_ds;
   *num = 0;
   size = 16;
   if ((*triloopEntropies != NULL) && (*triloopEntropies != defaultTriloopEntropies)) {
     free(*triloopEntropies);
     *triloopEntropies = NULL;
   }
   *triloopEntropies = (struct triloop*) safe_calloc(16, sizeof(struct triloop), _jmp_buf, o);
   while (readTLoop(&pt_ds, (*triloopEntropies)[*num].loop, &value, 1, _jmp_buf, o) != -1) {
      for (i = 0; i < 5; ++i)
        (*triloopEntropies)[*num].loop[i] = str2int((*triloopEntropies)[*num].loop[i]);
      (*triloopEntropies)[*num].value = value;
      ++*num;
      if (*num == size)        {
         size *= 2;
         *triloopEntropies = (struct triloop*) safe_realloc(*triloopEntropies, size * sizeof(struct triloop), _jmp_buf, o);
      }
   }
   *triloopEntropies = (struct triloop*) safe_realloc(*triloopEntropies, *num * sizeof(struct triloop), _jmp_buf, o);

   char *pt_dh = tp->triloop_dh;
   *num = 0;
   size = 16;

   if ((*triloopEnthalpies != NULL) && (*triloopEnthalpies != defaultTriloopEnthalpies)) {
     free(*triloopEnthalpies);
     *triloopEnthalpies = NULL;
   }
   *triloopEnthalpies = (struct triloop*) safe_calloc(16, sizeof(struct triloop), _jmp_buf, o);
   while (readTLoop(&pt_dh, (*triloopEnthalpies)[*num].loop, &value, 1, _jmp_buf, o) != -1) {
      for (i = 0; i < 5; ++i)
        (*triloopEnthalpies)[*num].loop[i] = str2int((*triloopEnthalpies)[*num].loop[i]);
      (*triloopEnthalpies)[*num].value = value;
      ++*num;
      if (*num == size) {
         size *= 2;
         *triloopEnthalpies = (struct triloop*) safe_realloc(*triloopEnthalpies, size * sizeof(struct triloop), _jmp_buf, o);
      }
   }
   *triloopEnthalpies = (struct triloop*) safe_realloc(*triloopEnthalpies, *num * sizeof(struct triloop), _jmp_buf, o);
}

static void 
getTetraloop(struct tetraloop** tetraloopEntropies, struct tetraloop** tetraloopEnthalpies, int* num, const thal_parameters *tp, jmp_buf _jmp_buf, thal_results* o)
{
   int i, size;
   double value;
   char *pt_ds = tp->tetraloop_ds;
   *num = 0;
   size = 16;
   if ((*tetraloopEntropies != NULL) && (*tetraloopEntropies != defaultTetraloopEntropies)) {
     free(*tetraloopEntropies);
     *tetraloopEntropies = NULL;
   }
   *tetraloopEntropies = (struct tetraloop*) safe_calloc(16, sizeof(struct tetraloop), _jmp_buf, o);
   while (readTLoop(&pt_ds, (*tetraloopEntropies)[*num].loop, &value, 0, _jmp_buf, o) != -1) {
      for (i = 0; i < 6; ++i)
        (*tetraloopEntropies)[*num].loop[i] = str2int((*tetraloopEntropies)[*num].loop[i]);
      (*tetraloopEntropies)[*num].value = value;
      ++*num;
      if (*num == size) {
         size *= 2;
         *tetraloopEntropies = (struct tetraloop*) safe_realloc(*tetraloopEntropies, size * sizeof(struct tetraloop), _jmp_buf, o);
      }
   }
   *tetraloopEntropies = (struct tetraloop*) safe_realloc(*tetraloopEntropies, *num * sizeof(struct tetraloop), _jmp_buf, o);

   char *pt_dh = tp->tetraloop_dh;
   *num = 0;
   size = 16;
   if ((*tetraloopEnthalpies != NULL) && (*tetraloopEnthalpies != defaultTetraloopEnthalpies)) {
     free(*tetraloopEnthalpies);
     *tetraloopEnthalpies = NULL;
   }
   *tetraloopEnthalpies = (struct tetraloop*) safe_calloc(16, sizeof(struct tetraloop), _jmp_buf, o);
   while (readTLoop(&pt_dh, (*tetraloopEnthalpies)[*num].loop, &value, 0, _jmp_buf, o) != -1) {
      for (i = 0; i < 6; ++i)
        (*tetraloopEnthalpies)[*num].loop[i] = str2int((*tetraloopEnthalpies)[*num].loop[i]);
      (*tetraloopEnthalpies)[*num].value = value;
      ++*num;
      if (*num == size) {
         size *= 2;
         *tetraloopEnthalpies = (struct tetraloop*) safe_realloc(*tetraloopEnthalpies, size * sizeof(struct tetraloop), _jmp_buf, o);
      }
   }
   *tetraloopEnthalpies = (struct tetraloop*) safe_realloc(*tetraloopEnthalpies, *num * sizeof(struct tetraloop), _jmp_buf, o);
}

static void 
tableStartATS(double atp_value, double atpS[5][5])
{

   int i, j;
   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       atpS[i][j] = 0.00000000001;
   atpS[0][3] = atpS[3][0] = atp_value;
}


static void 
tableStartATH(double atp_value, double atpH[5][5])
{

   int i, j;
   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       atpH[i][j] = 0.0;

   atpH[0][3] = atpH[3][0] = atp_value;
}
