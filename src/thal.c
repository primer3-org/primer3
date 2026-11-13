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
#include "thal_internal.h"
/* The thermodynamic parameter tables live in thal_params.c (the sole TU
   that includes thal_default_params.h).  We see them via the extern
   declarations in thal_internal.h. */

#define STR(X) #X
#define LONG_SEQ_ERR_STR(MAX_LEN) "Target sequence length > maximum allowed (" STR(MAX_LEN) ") in thermodynamic alignment"
#define XSTR(X) STR(X)

/* file-private constants */
static const int min_hrpn_loop = 3;
static const double R = 1.9872; /* cal/Kmol */
static const double ILAS = (-300 / 310.15); /* Internal Loop Entropy ASymmetry correction -0.3kcal/mol*/
static const double ILAH = 0.0; /* Internal Loop EntHalpy Asymmetry correction */
static const double dplx_init_H_dimer = 200;
static const double dplx_init_S_dimer = -5.7;
static const int is_complement[5][5] = { /* watson-crick bp -> 1; else 0 */
     {0, 0, 0, 1, 0}, /* A, C, G, T, N */
     {0, 0, 1, 0, 0},
     {0, 1, 0, 0, 0},
     {1, 0, 0, 0, 0},
     {0, 0, 0, 0, 0}};

/* constants exported through thal_internal.h */
const double MinEntropy = -3224.0;
const double ABSOLUTE_ZERO = 273.15;
const double TEMP_KELVIN = 310.15;
/* MAX_LOOP / MIN_LOOP are exported through thal.h */
const int MAX_LOOP = 30;
const int MIN_LOOP = 0;

/* tracer-stack node for traceback_monomer */
struct tracer {
   int i;
   int j;
   int mtrx; /* 0 = stem (dpt), 1 = prefix (send5/hend5) */
   struct tracer *next;
};


//=====================================================================================
//Functions for dimer calculation
//=====================================================================================
 /* calc-s thermod values into dynamic progr table (dimer) */
static void fillMatrix_dimer(int maxLoop, struct dpt_entry **dpt, double RC,
                                 const unsigned char *numSeq1, const unsigned char *numSeq2, int oligo1_len,
                                 int oligo2_len, thal_results* o);
static void calc_bulge_internal_dimer(int ii, int jj, int i, int j, double* EntropyEnthalpy,
                                 const struct dpt_entry* const *dpt, const unsigned char *numSeq1, const unsigned char *numSeq2);
static void traceback_dimer(int i, int j, int* ps1, int* ps2, const struct dpt_entry* const *dpt);
/* calculate terminal entropy S and terminal enthalpy H starting reading from 5'end */
static void LSH(int i, int j, double* EntropyEnthalpy, double RC,
               const unsigned char *numSeq1, const unsigned char *numSeq2);

//=====================================================================================
//Functions for dimer and hairpin calculation
//=====================================================================================
/* calculate terminal entropy S and terminal enthalpy H starting reading from 3'end */
static void RSH(int i, int j, double* EntropyEnthalpy, double RC, double dplx_init_S,
               double dplx_init_H, const unsigned char *numSeq1, const unsigned char *numSeq2);

//=====================================================================================
//Functions for hairpin calculation
//=====================================================================================
static void initMatrix_monomer(struct dpt_entry **dpt, const unsigned char *numSeq1,
                              int oligo1_len); /* initiates thermodynamic parameter tables of entropy and enthalpy for monomer */
 /* calcs thermod values into dynamic progr table (monomer) */
static void fillMatrix_monomer(int maxLoop, struct dpt_entry **dpt, double RC, const unsigned char *numSeq1,
                              int oligo1_len, thal_results* o);
static void calc_bulge_internal_monomer(int ii, int jj, int i, int j, double* EntropyEnthalpy, int traceback, int maxLoop,
                                    const struct dpt_entry* const *dpt, double RC,
                                    const unsigned char *numSeq1);
static void calc_terminal_bp(double temp, const struct dpt_entry* const *dpt, double *send5, double *hend5,
                              double RC, const unsigned char *numSeq1, int oligo1_len);
/* (S, H) for one of the four exterior-loop attachment patterns
   END5_1..END5_4; shared between calc_terminal_bp and traceback_monomer */
static double end5_candidate(int variant, int i, int k, double RC,
                             const struct dpt_entry* const *dpt,
                             const unsigned char *numSeq1,
                             const double *hend5, const double *send5,
                             double *EntropyEnthalpy);
/* finds monomer structure that has maximum Tm */
static void calc_hairpin(int i, int j, double* EntropyEnthalpy, const struct dpt_entry* const *dpt,
                        double RC, const unsigned char *numSeq1, int oligo1_len);
static void push(struct tracer**, int, int, int, jmp_buf, thal_results*); /* to add elements to struct */
static void traceback_monomer(int*, int, const struct dpt_entry* const *dpt, double *send5, double *hend5, double RC,
                              const unsigned char *numSeq1, int oligo1_len, jmp_buf, thal_results*);

/* drawDimer / drawHairpin are declared in thal_internal.h (defined in thal_draw.c) */

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
static double saltCorrectS(double mv, double dv, double dntp); /* SantaLucia-style salt correction */
static int thal_check_errors(const unsigned char *oligo_f, const unsigned char *oligo_r, int *len_f, int *len_r, const thal_args *a, thal_results *o);

//=====================================================================================
//Functions for allocating memory
//=====================================================================================
/* safe_calloc/malloc/realloc live in thal_util.c (declared in thal_internal.h) */
static struct dpt_entry **allocate_DPT(int oligo1_len, int oligo2_len, jmp_buf _jmp_buf, thal_results *o);
static void free_DPT(struct dpt_entry **dpt);

//=====================================================================================
//Functions for string manipulation
//=====================================================================================
static int length_unsig_char(const unsigned char * str); /* returns length of unsigned char; to avoid warnings while compiling */
/* str2int — declared in thal_internal.h (shared with thal_params.c) */
static void reverse(unsigned char *s);

/* Parameter file loaders (readDouble, readParamFile, readLoop, readTLoop,
   get*, tableStart*) live in thal_params.c.  Output-formatting helpers
   (strcatc, save_append_*) live in thal_draw.c.  Anything called from
   here is declared in thal_internal.h. */

/* central method: execute all sub-methods for calculating secondary
   structure for dimer or for monomer */
void 
thal(const unsigned char *oligo_f, 
     const unsigned char *oligo_r, 
     const thal_args *a,
     const thal_mode mode,
     thal_results *o)
{
   int len_f, len_r;
   jmp_buf _jmp_buf;
   double RC;
   int oligo1_len;
   int oligo2_len;

   unsigned char *numSeq1 = NULL;
   unsigned char *numSeq2 = NULL;
   unsigned char *oligo1 = NULL;
   unsigned char *oligo2 = NULL;
   struct dpt_entry **dpt = NULL;
   if (NULL == o) return; /* Leave it to the caller to crash */
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
   /*** INIT values for unimolecular and bimolecular structures ***/
   if (a->type==4) { /* unimolecular folding */
      RC=0;
   } else  {
      /* hybridization of two oligos */
      if(symmetry_thermo(oligo1) && symmetry_thermo(oligo2)) {
         RC = R  * log(a->dna_conc/1000000000.0);
      } else {
         RC = R  * log(a->dna_conc/4000000000.0);
      }
      reverse(oligo2); /* REVERSE oligo2, so it goes to dpt 3'->5' direction */
   }
   /* convert nucleotides to numbers */
   numSeq1 = (unsigned char*) safe_malloc(oligo1_len + 2, _jmp_buf, o);
   numSeq2 = (unsigned char*) safe_malloc(oligo2_len + 2, _jmp_buf, o);

   //Allocate Dynamic Programming Table
   dpt = allocate_DPT(oligo1_len, oligo2_len, _jmp_buf, o);

   for(int i = 0; i < oligo1_len; i++) oligo1[i] = toupper(oligo1[i]);
   for(int i = 0; i < oligo2_len; i++) oligo2[i] = toupper(oligo2[i]);
   for(int i = 1; i <= oligo1_len; ++i) numSeq1[i] = str2int(oligo1[i - 1]);
   for(int i = 1; i <= oligo2_len; ++i) numSeq2[i] = str2int(oligo2[i - 1]);
   numSeq1[0] = numSeq1[oligo1_len + 1] = numSeq2[0] = numSeq2[oligo2_len + 1] = 4; /* mark as N-s */

   if (a->type==4) { /* calculate structure of monomer */
      double mh, ms;
      double *hend5 = (double*) safe_malloc((oligo1_len + 1) * sizeof(double), _jmp_buf, o);
      double *send5 = (double*) safe_malloc((oligo1_len + 1) * sizeof(double), _jmp_buf, o);
      int *bp = (int*) safe_calloc(oligo1_len, sizeof(int), _jmp_buf, o);
      initMatrix_monomer(dpt, numSeq1, oligo1_len);
      fillMatrix_monomer(a->maxLoop, dpt, RC, numSeq1, oligo1_len, o);
      calc_terminal_bp(a->temp, (const struct dpt_entry **)dpt, send5, hend5, RC, numSeq1, oligo1_len);
      send5_dump(send5, hend5, oligo1_len, "after calc_terminal_bp");
      mh = hend5[oligo1_len];
      ms = send5[oligo1_len];
      o->align_end_1 = (int) mh;
      o->align_end_2 = (int) ms;
      for (int i = 0; i < oligo1_len; ++i) bp[i] = 0;
      if(isFinite(mh)) {
        traceback_monomer(bp, a->maxLoop, (const struct dpt_entry **) dpt, send5, hend5, RC, numSeq1, oligo1_len, _jmp_buf, o);
        /* traceback for unimolecular structure */
        o->sec_struct=drawHairpin(bp, mh, ms, mode,a->temp, oligo1, oligo2, saltCorrectS(a->mv,a->dv,a->dntp), oligo1_len, oligo2_len, _jmp_buf, o); /* if mode=THL_FAST or THL_DEBUG_F then return after printing basic therm data */
      } else if((mode != THL_FAST) && (mode != THL_DEBUG_F) && (mode != THL_STRUCT)) {
        fputs("No secondary structure could be calculated\n",stderr);
      }

      if(o->temp==-_INFINITY && (!strcmp(o->msg, ""))) o->temp=0.0;
      free(bp);
      free_DPT(dpt);
      free(numSeq1);
      free(numSeq2);
      free(send5);
      free(hend5);
      free(oligo1);
      free(oligo2);
      dpt_dump_close();
      return;
   } else if(a->type!=4) { /* Hybridization of two moleculs */
      int bestI;
      int bestJ;
      double G1, bestG;
      double SH[2];
      int i, j;
      int *ps1, *ps2;
      ps1 = (int*) safe_calloc(oligo1_len, sizeof(int), _jmp_buf, o);
      ps2 = (int*) safe_calloc(oligo2_len, sizeof(int), _jmp_buf, o);
      fillMatrix_dimer(a->maxLoop, dpt, RC, numSeq1, numSeq2, oligo1_len, oligo2_len, o);
      /* calculate terminal basepairs */
      bestI = bestJ = 1; 
      G1 = bestG = _INFINITY;

      if(a->type==1)
         i = 1;
      else
         i = oligo1_len; //THAL_END1: oligo1 3' end must be part of terminal pair

      for (; i <= oligo1_len; i++) {
         for (j = 1; j <= oligo2_len; j++) {
            if(is_complement[numSeq1[i]][numSeq2[j]]){
               RSH(i, j, SH, RC, dplx_init_S_dimer, dplx_init_H_dimer, numSeq1, numSeq2);
               G1 = (dpt[i][j].h+ SH[1] + dplx_init_H_dimer) - TEMP_KELVIN*(dpt[i][j].s + SH[0] + dplx_init_S_dimer);  
               if(G1<bestG){
                  bestG = G1;
                  bestI = i;
                  bestJ = j;
               }
            }
         }
      }

      /* tracebacking */
      if(is_complement[numSeq1[bestI]][numSeq2[bestJ]]){
         RSH(bestI, bestJ, SH, RC, dplx_init_S_dimer, dplx_init_H_dimer, numSeq1, numSeq2);
         traceback_dimer(bestI, bestJ, ps1, ps2, (const struct dpt_entry **)dpt);
         int N=0;
         for(i=0;i<oligo1_len;i++)
            if(ps1[i]>0) ++N;
         N--;
         o->dh = dpt[bestI][bestJ].h+ SH[1] + dplx_init_H_dimer;
         o->ds = (dpt[bestI][bestJ].s + SH[0] + dplx_init_S_dimer);
         o->ds = o->ds + (N * saltCorrectS(a->mv,a->dv,a->dntp));
         o->dg = (o->dh) - (a->temp * o->ds);
         o->temp = ((o->dh) / (o->ds + RC)) - ABSOLUTE_ZERO;
         o->align_end_1=bestI;
         o->align_end_2=bestJ;
         if(mode != THL_FAST)
            o->sec_struct=drawDimer(ps1, ps2, mode, a->temp, oligo1, oligo2, oligo1_len, oligo2_len, _jmp_buf, o);
      } else  {
         o->temp = 0.0;
         /* fputs("No secondary structure could be calculated\n",stderr); */
      }
      free(ps1);
      free(ps2);
      free(oligo2);
      free_DPT(dpt);
      free(numSeq1);
      free(numSeq2);
      free(oligo1);
      dpt_dump_close();
      return;
   }
   dpt_dump_close();
   return;
}
/*** END thal() ***/

//=====================================================================================
//Functions for dimer calculation
//=====================================================================================

static void fillMatrix_dimer(int maxLoop, struct dpt_entry **dpt, double RC,
                                 const unsigned char *numSeq1, const unsigned char *numSeq2, int oligo1_len,
                                 int oligo2_len, thal_results* o)
{
   int d, i, j, ii, jj;
   double SH[2];
   double saved_RSH[2];
   double bestG, newG;

   for (i = 1; i <= oligo1_len; ++i) {
      for (j = 1; j <= oligo2_len; ++j) {
         dpt[i][j].tb_i = -1;
         dpt[i][j].tb_j = -1;
         if(!is_complement[numSeq1[i]][numSeq2[j]]){
            dpt[i][j].h = _INFINITY;
            dpt[i][j].s = -1.0;
         } else {
            LSH(i, j, SH, RC, numSeq1, numSeq2);
            dpt[i][j].s = SH[0];
            dpt[i][j].h = SH[1];
            if (i > 1 && j > 1){
               RSH(i, j, SH, RC, dplx_init_S_dimer, dplx_init_H_dimer, numSeq1, numSeq2);
               saved_RSH[0] = SH[0];
               saved_RSH[1] = SH[1];

               if(is_complement[numSeq1[i-1]][numSeq2[j-1]]){
                  dpt[i][j].s = dpt[i-1][j-1].s + stackEntropies[numSeq1[i-1]][numSeq1[i]][numSeq2[j-1]][numSeq2[j]];
                  dpt[i][j].h = dpt[i-1][j-1].h + stackEnthalpies[numSeq1[i-1]][numSeq1[i]][numSeq2[j-1]][numSeq2[j]];
                  dpt[i][j].tb_i = i-1;
                  dpt[i][j].tb_j = j-1;
               }
               bestG = dpt[i][j].h + saved_RSH[1] - TEMP_KELVIN * (dpt[i][j].s + saved_RSH[0]);
               for(d = 3; d <= maxLoop + 2; d++) { /* max=30, length over 30 is not allowed */
                  ii = i - 1;
                  jj = - ii - d + (j + i);
                  if (jj < 1) {
                     ii += jj-1;
                     jj = 1;
                  }
                  for (; ii > 0 && jj < j; --ii, ++jj) {
                     if(is_complement[numSeq1[ii]][numSeq2[jj]]){
                        calc_bulge_internal_dimer(ii, jj, i, j, SH, (const struct dpt_entry **)dpt, numSeq1, numSeq2);
                        newG = SH[1]+saved_RSH[1] -TEMP_KELVIN*(SH[0]+saved_RSH[0]);
                        if(newG < bestG ) {
                           bestG = newG;
                           //This condition is only here to pick the same alignment as the old traceback
                           //when there are multiple, equal alignments.
                           if(!(equal(dpt[i][j].h, SH[1]) && equal(dpt[i][j].s, SH[0]))){
                              dpt[i][j].tb_i = ii;
                              dpt[i][j].tb_j = jj;
                           }
                           dpt[i][j].h = SH[1];
                           dpt[i][j].s = SH[0];
                        }
                     }
                  }
               }
            } /* if */
         }
      } /* for */
      {
         char _lbl[64];
         snprintf(_lbl, sizeof _lbl, "fillMatrix_dimer after row i=%d", i);
         dpt_dump(dpt, oligo1_len, oligo2_len, 0, _lbl);
      }
   } /* for */
}

static void 
calc_bulge_internal_dimer(int i, int j, int ii, int jj, double* EntropyEnthalpy,
                        const struct dpt_entry* const *dpt,
                        const unsigned char *numSeq1, const unsigned char *numSeq2)
{
   int loopSize1, loopSize2, loopSize;
   double S,H;
   S = -1.0;
   H = _INFINITY;
   loopSize1 = ii - i - 1;
   loopSize2 = jj - j - 1;
   loopSize = loopSize1 + loopSize2-1;
   if(loopSize1 == 0 || loopSize2 == 0) { /* only bulges have to be considered */
      //bulge loop of size one is treated differently. the intervening nn-pair must be added
      if(loopSize2 == 1 || loopSize1 == 1) {
         H = bulgeLoopEnthalpies[loopSize] +
             stackEnthalpies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
         S = bulgeLoopEntropies[loopSize] +
             stackEntropies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
      } else { /* we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30 */
         H = bulgeLoopEnthalpies[loopSize] + atpH[numSeq1[i]][numSeq2[j]] + atpH[numSeq1[ii]][numSeq2[jj]];
         S = bulgeLoopEntropies[loopSize] + atpS[numSeq1[i]][numSeq2[j]] + atpS[numSeq1[ii]][numSeq2[jj]];
      }
   } else if (loopSize1 == 1 && loopSize2 == 1) {
      S = stackint2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
          stackint2Entropies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]];
      H = stackint2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
          stackint2Enthalpies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]];
   } else { /* only internal loops */
      //Only calculate H and S if there is a mismatch in both nearest neighbor stacks. 
      //This removes the need for the tstack table and improves performance.
      //At this point, i,j and ii,jj are known to be complementary, so only need to check neighbors.
      //The only difference between tstack and tstack2 tables is when both pairs are complementary
      if((!is_complement[numSeq1[ii-1]][numSeq2[jj-1]]) && (!is_complement[numSeq1[i+1]][numSeq2[j+1]])){
         H = interiorLoopEnthalpies[loopSize] + tstack2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
             tstack2Enthalpies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]] +
             (ILAH * abs(loopSize1 - loopSize2));
         S = interiorLoopEntropies[loopSize] + tstack2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
             tstack2Entropies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]] + 
             (ILAS * abs(loopSize1 - loopSize2));
      }
   }
   EntropyEnthalpy[0] = S + dpt[i][j].s;
   EntropyEnthalpy[1] = H + dpt[i][j].h;
   return;
}

/*Each time the DPTs are modified, the traceback matrix stores the i and j values of the last
complementary pair used in the S and H calculations. To find the alignment, just start at the
i and j with the lowest dG and follow the tracback matrix.
Example: if the terminal i and j are 10,10 you go to traceback_matrix[10][10] to see
what the previous pairing is. If there is not an internal loop it would be 9,9. If it
was 7,9 that would mean that there was a bulge loop, etc.
After following the chain of indices, you know that i,j is the 5' terminal pair when
traceback_matrix[i][j] == -1,-1*/
static void 
traceback_dimer(int i, int j, int* ps1, int* ps2, const struct dpt_entry* const *dpt)
{
   int tmp_i;
   ps1[i-1] = j;
   ps2[j-1] = i;
   while(dpt[i][j].tb_i >= 0){
      tmp_i = dpt[i][j].tb_i;
      j = dpt[i][j].tb_j;
      i = tmp_i;
      ps1[i-1] = j;
      ps2[j-1] = i;
   }
}

static void 
LSH(int i, int j, double* EntropyEnthalpy, double RC, const unsigned char *numSeq1, const unsigned char *numSeq2)
{
   double S1, H1, T1, G1;
   double S2, H2, T2, G2;
   S1 = S2 = -1.0;
   H1 = H2 = -_INFINITY;
   T1 = T2 = -_INFINITY;

   S1 = atpS[numSeq1[i]][numSeq2[j]] + tstack2Entropies[numSeq2[j]][numSeq2[j-1]][numSeq1[i]][numSeq1[i-1]];
   H1 = atpH[numSeq1[i]][numSeq2[j]] + tstack2Enthalpies[numSeq2[j]][numSeq2[j-1]][numSeq1[i]][numSeq1[i-1]];
   G1 = H1 - TEMP_KELVIN*S1;
   
   if(!is_complement[numSeq1[i-1]][numSeq2[j-1]]){ 
      S2 = atpS[numSeq1[i]][numSeq2[j]];
      H2 = atpH[numSeq1[i]][numSeq2[j]];
      if ((numSeq1[i-1] == 4)){
         S2 += dangleEntropies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]];
         H2 += dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]];
      } else if ((numSeq2[j-1] == 4)){
         S2 += dangleEntropies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
         H2 += dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      } else {
         S2 += dangleEntropies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]] +
               dangleEntropies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
         H2 += dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]] +
               dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      }
      G2 = H2 - TEMP_KELVIN*S2;
      T2 = (H2 + dplx_init_H_dimer) / (S2 + dplx_init_S_dimer + RC);
      if(G1<0) {
         T1 = (H1 + dplx_init_H_dimer) / (S1 + dplx_init_S_dimer + RC);
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
   T2 = (H2 + dplx_init_H_dimer) / (S2 + dplx_init_S_dimer + RC);
   G1 = H1 -TEMP_KELVIN*S1;   
   G2 = H2 -TEMP_KELVIN*S2;
   if(T1 < T2) {
      EntropyEnthalpy[0] = S2;
      EntropyEnthalpy[1] = H2;
   } else {
      EntropyEnthalpy[0] = S1;
      EntropyEnthalpy[1] = H1;
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
   if (is_complement[numSeq1[i]][numSeq2[j]] == 0) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      return;
   }
   S1 = atpS[numSeq1[i]][numSeq2[j]] + tstack2Entropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   H1 = atpH[numSeq1[i]][numSeq2[j]] + tstack2Enthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   G1 = H1 - TEMP_KELVIN*S1;
   
   if(!is_complement[numSeq1[i+1]][numSeq2[j+1]]){
      S2 = atpS[numSeq1[i]][numSeq2[j]];
      H2 = atpH[numSeq1[i]][numSeq2[j]];
      if((numSeq2[j+1] == 4)){
         S2 += dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]];
         H2 += dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]];
         G2 = H2 - TEMP_KELVIN*S2;
      } else if((numSeq1[i+1] == 4)){
         S2 += dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
         H2 += dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
         G2 = H2 - TEMP_KELVIN*S2;
      } else {
         S2 += dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +
               dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
         H2 += dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +
               dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
         G2 = H2 - TEMP_KELVIN*S2;
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

   S2 = atpS[numSeq1[i]][numSeq2[j]];
   H2 = atpH[numSeq1[i]][numSeq2[j]];
   T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
   G1 = H1 -TEMP_KELVIN*S1;
   G2 =  H2 -TEMP_KELVIN*S2;
   if(T1 < T2) {
      EntropyEnthalpy[0] = S2;
      EntropyEnthalpy[1] = H2;
   } else {
      EntropyEnthalpy[0] = S1;
      EntropyEnthalpy[1] = H1;
   }
   return;
}

//=====================================================================================
//Functions used in hairpin calculations
//=====================================================================================

static void 
initMatrix_monomer(struct dpt_entry **dpt, const unsigned char *numSeq1, int oligo1_len)
{
   int i, j;
   for (i = 1; i <= oligo1_len; ++i){
     for (j = i; j <= oligo1_len; ++j){
       if (j - i < min_hrpn_loop + 1 || (is_complement[numSeq1[i]][numSeq1[j]] == 0)) {
          dpt[i][j].h = _INFINITY;
          dpt[i][j].s = -1.0;
       } else {
          dpt[i][j].h = 0.0;
          dpt[i][j].s = MinEntropy;
       }
       dpt[i][j].tb_i = -1;
       dpt[i][j].tb_j = -1;
      }
   }
   dpt_dump(dpt, oligo1_len, oligo1_len, 1, "initMatrix_monomer");
}
static void
fillMatrix_monomer(int maxLoop, struct dpt_entry **dpt, double RC, const unsigned char *numSeq1, int oligo1_len, thal_results* o)
{
   int i, j;
   double SH[2];
   double T0, T1;
   double S1;
   double H1;

   for (j = 2; j <= oligo1_len; ++j) {
      for (i = j - min_hrpn_loop - 1; i >= 1; --i) {
         if (is_complement[numSeq1[i]][numSeq1[j]]) {
            T0 = (dpt[i][j].h) /(dpt[i][j].s + RC);
            S1 = (dpt[i+1][j-1].s + stackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq1[j]][numSeq1[j-1]]);
            H1 = (dpt[i+1][j-1].h + stackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq1[j]][numSeq1[j-1]]);
            T1 = (H1) /(S1 + RC);

            if(T1 > T0) {
               dpt[i][j].s = S1;
               dpt[i][j].h = H1;
            }

            int d, ii, jj;
            for (d = j - i - 3; d >= min_hrpn_loop + 1 && d >= j - i - 2 - maxLoop; --d){
               for (ii = i + 1; ii < j - d && ii <= oligo1_len; ++ii) {
                  jj = d + ii;
                  SH[0] = -1.0;
                  SH[1] = _INFINITY;
                  if (is_complement[numSeq1[ii]][numSeq1[jj]]) {
                     calc_bulge_internal_monomer(i, j, ii, jj, SH, 0, maxLoop, (const struct dpt_entry **) dpt, RC, numSeq1);
                     if(isFinite(SH[1])) {
                        dpt[i][j].h = SH[1];
                        dpt[i][j].s = SH[0];
                     }
                  }
               }
            }
            calc_hairpin(i, j, SH, (const struct dpt_entry **)dpt, RC, numSeq1, oligo1_len);
            dpt[i][j].s = SH[0];
            dpt[i][j].h = SH[1];
        }
      }
      {
         char _lbl[64];
         snprintf(_lbl, sizeof _lbl, "fillMatrix_monomer after outer j=%d", j);
         dpt_dump(dpt, oligo1_len, oligo1_len, 1, _lbl);
      }
   }
}

static void 
calc_hairpin(int i, int j, double* EntropyEnthalpy, const struct dpt_entry* const *dpt, double RC, const unsigned char *numSeq1, int oligo1_len)
{
   int loopSize = j - i - 1;
   double G1, G2;
   G1 = G2 = -_INFINITY;
   double SH[2];
   SH[0] = -1.0;
   SH[1] = _INFINITY;

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
      if (loopSize == 4) { /* terminal mismatch, tetraloop bonus, hairpin of 4 */
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
   } else if(loopSize == 3){ /* for loops 3 bp in length at-penalty is considered */
      EntropyEnthalpy[1] += atpH[numSeq1[i]][numSeq1[j]];
      EntropyEnthalpy[0] += atpS[numSeq1[i]][numSeq1[j]];
      struct triloop* loop;
      if (numTriloops) {
         if ((loop = (struct triloop*) bsearch(numSeq1 + i, triloopEnthalpies, numTriloops, sizeof(struct triloop), comp3loop)))
           EntropyEnthalpy[1] += loop->value;
         if ((loop = (struct triloop*) bsearch(numSeq1 + i, triloopEntropies, numTriloops, sizeof(struct triloop), comp3loop)))
           EntropyEnthalpy[0] += loop->value;
      }
   }

   RSH(i, j, SH, RC, 0.0, 0.0, numSeq1, numSeq1);
   G1 = EntropyEnthalpy[1]+SH[1] -TEMP_KELVIN*(EntropyEnthalpy[0]+SH[0]);
   G2 = dpt[i][j].h+SH[1] -TEMP_KELVIN*(dpt[i][j].s+SH[0]);
     if(G2 < G1) {
      EntropyEnthalpy[0] = dpt[i][j].s;
      EntropyEnthalpy[1] = dpt[i][j].h;
   }
   return;
}

static void 
calc_bulge_internal_monomer(int i, int j, int ii, int jj, double* EntropyEnthalpy, int traceback, int maxLoop,
                        const struct dpt_entry* const *dpt, double RC, 
                        const unsigned char *numSeq1)
{
   int loopSize1, loopSize2, loopSize;
   double T1, T2;
   double S,H;
   /* int N, N_loop; Triinu, please review */
   T1 = T2 = -_INFINITY;
   S = -1.0;
   H = _INFINITY;
   loopSize1 = ii - i - 1;
   loopSize2 = j - jj - 1;
   loopSize = loopSize1 + loopSize2 -1; /* for indx only */
   if(loopSize1 == 0 || loopSize2 == 0) { /* only bulges have to be considered */
      if(loopSize2 == 1 || loopSize1 == 1) { /* bulge loop of size one is treated differently
                                              the intervening nn-pair must be added */
         H = bulgeLoopEnthalpies[loopSize] +
            stackEnthalpies[numSeq1[i]][numSeq1[ii]][numSeq1[j]][numSeq1[jj]];
         S = bulgeLoopEntropies[loopSize] +
            stackEntropies[numSeq1[i]][numSeq1[ii]][numSeq1[j]][numSeq1[jj]];
      } else { /* we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30 */
         H = bulgeLoopEnthalpies[loopSize] + atpH[numSeq1[i]][numSeq1[j]] + atpH[numSeq1[ii]][numSeq1[jj]];
         S = bulgeLoopEntropies[loopSize] + atpS[numSeq1[i]][numSeq1[j]] + atpS[numSeq1[ii]][numSeq1[jj]];
      }
   } /* end of calculating bulges */
   else if (loopSize1 == 1 && loopSize2 == 1) {
      /* mismatch nearest neighbor parameters */
      S = stackint2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq1[j]][numSeq1[j-1]] +
        stackint2Entropies[numSeq1[jj]][numSeq1[jj+1]][numSeq1[ii]][numSeq1[ii-1]];
      H = stackint2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq1[j]][numSeq1[j-1]] +
        stackint2Enthalpies[numSeq1[jj]][numSeq1[jj+1]][numSeq1[ii]][numSeq1[ii-1]];
   } else { /* only internal loops */
      if((!is_complement[numSeq1[ii-1]][numSeq1[jj+1]]) && (!is_complement[numSeq1[i+1]][numSeq1[j-1]])){
         H = interiorLoopEnthalpies[loopSize] + tstack2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq1[j]][numSeq1[j-1]] +
         tstack2Enthalpies[numSeq1[jj]][numSeq1[jj+1]][numSeq1[ii]][numSeq1[ii-1]]
         + (ILAH * abs(loopSize1 - loopSize2));
         S = interiorLoopEntropies[loopSize] + tstack2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq1[j]][numSeq1[j-1]] +
         tstack2Entropies[numSeq1[jj]][numSeq1[jj+1]][numSeq1[ii]][numSeq1[ii-1]] + (ILAS * abs(loopSize1 - loopSize2));
      }
   }

   if(traceback!=1) {
      H += dpt[ii][jj].h; /* bulge koos otsaga, st bulge i,j-ni */
      S += dpt[ii][jj].s;
   }

   if(isFinite(H)) {
      T1 = (H) / (S + RC);
      T2 = (dpt[i][j].h) / ((dpt[i][j].s) + RC);

      if((T1 > T2) || traceback==1) {
         EntropyEnthalpy[0] = S;
         EntropyEnthalpy[1] = H;
      }
   }
   return;
}

/* (S, H) for one of the four exterior-loop attachment patterns used in
   calc_terminal_bp / traceback_monomer.

     variant 1 = END5_1: stem (k+1..i), no flanks
         5' k+1  k+2  3'
         3'  i   i-1  5'

     variant 2 = END5_2: stem (k+2..i) with 5' overhang at k+1
         5' k+1  k+2       3'
         3'       i   i-1  5'

     variant 3 = END5_3: stem (k+1..i-1) with 3' overhang at i
         5'      k+1  k+2  3'
         3'  i   i-1       5'

     variant 4 = END5_4: stem (k+2..i-1) with both flanks (terminal mismatch)
         5' k+1 k+2 3'
         3'  i  i-1 5'

   The (S, H) value includes the prefix structure (hend5[k], send5[k])
   iff its melting point T0 = H/(S+RC) is non-negative — otherwise the
   prefix is assumed to be the empty structure. T0 is returned because
   traceback_monomer uses it to decide whether to push a prefix-DP
   traceback frame.

   No is_complement validity check is performed. calc_terminal_bp
   applies one per variant; traceback_monomer relies on _INFINITY in
   invalid dpt cells to reject candidates via the downstream
   H<=0, S<=0 test. */
static double
end5_candidate(int variant, int i, int k, double RC,
               const struct dpt_entry* const *dpt,
               const unsigned char *numSeq1,
               const double *hend5, const double *send5,
               double *EntropyEnthalpy)
{
   double T0 = hend5[k] / (send5[k] + RC);
   double H = _INFINITY;
   double S = -1.0;
   switch (variant) {
   case 1: /* END5_1 */
      H = atpH[numSeq1[k + 1]][numSeq1[i]] + dpt[k + 1][i].h;
      S = atpS[numSeq1[k + 1]][numSeq1[i]] + dpt[k + 1][i].s;
      break;
   case 2: /* END5_2 */
      H = atpH[numSeq1[k + 2]][numSeq1[i]]
        + dangleEnthalpies5[numSeq1[i]][numSeq1[k + 2]][numSeq1[k + 1]]
        + dpt[k + 2][i].h;
      S = atpS[numSeq1[k + 2]][numSeq1[i]]
        + dangleEntropies5[numSeq1[i]][numSeq1[k + 2]][numSeq1[k + 1]]
        + dpt[k + 2][i].s;
      break;
   case 3: /* END5_3 */
      H = atpH[numSeq1[k + 1]][numSeq1[i - 1]]
        + dangleEnthalpies3[numSeq1[i - 1]][numSeq1[i]][numSeq1[k + 1]]
        + dpt[k + 1][i - 1].h;
      S = atpS[numSeq1[k + 1]][numSeq1[i - 1]]
        + dangleEntropies3[numSeq1[i - 1]][numSeq1[i]][numSeq1[k + 1]]
        + dpt[k + 1][i - 1].s;
      break;
   case 4: /* END5_4 — Probably should not consider the case where k+1 and i
              are complementary; in that case k+1 and i should be the
              terminal bp instead. */
      H = atpH[numSeq1[k + 2]][numSeq1[i - 1]]
        + tstack2Enthalpies[numSeq1[i - 1]][numSeq1[i]][numSeq1[k + 2]][numSeq1[k + 1]]
        + dpt[k + 2][i - 1].h;
      S = atpS[numSeq1[k + 2]][numSeq1[i - 1]]
        + tstack2Entropies[numSeq1[i - 1]][numSeq1[i]][numSeq1[k + 2]][numSeq1[k + 1]]
        + dpt[k + 2][i - 1].s;
      break;
   }
   if (T0 >= 0.0) {
      H += hend5[k];
      S += send5[k];
   }
   EntropyEnthalpy[0] = S;
   EntropyEnthalpy[1] = H;
   return T0;
}

static void
calc_terminal_bp(double temp, const struct dpt_entry* const *dpt, double *send5, double *hend5, double RC,
                  const unsigned char *numSeq1, int oligo1_len) { /* compute exterior loop */
   int i, k;
   send5[0] = send5[1] = -1.0;
   hend5[0] = hend5[1] = _INFINITY;

   double max_tm, S_max, H_max, H, S, T1, G;
   double SH[2];
   for (i = 2; i <= oligo1_len; i++){
      max_tm = (hend5[i - 1]) / (send5[i - 1] + RC);
      S_max = send5[i-1];
      H_max = hend5[i-1];
      for(k = 0; k <= i - min_hrpn_loop - 2; ++k) {
         /* END5_1: stem (k+1..i), no flanks */
         if (is_complement[numSeq1[i-1]][numSeq1[k+2]] && is_complement[numSeq1[i]][numSeq1[k+1]]) {
            (void) end5_candidate(1, i, k, RC, dpt, numSeq1, hend5, send5, SH);
            S = SH[0]; H = SH[1];
            if (H <= 0 && S <= 0) {
               T1 = H / (S + RC);
               G = H - temp * S;
               if ((max_tm < T1) && (G < 0.0)) {
                  H_max = H; S_max = S; max_tm = T1;
               }
            }
         }

         /* END5_2: stem (k+2..i) with 5' overhang at k+1 */
         if (is_complement[numSeq1[i]][numSeq1[k+2]]) {
            (void) end5_candidate(2, i, k, RC, dpt, numSeq1, hend5, send5, SH);
            S = SH[0]; H = SH[1];
            if (H <= 0 && S <= 0) {
               T1 = H / (S + RC);
               G = H - temp * S;
               if ((max_tm < T1) && (G < 0.0)) {
                  H_max = H; S_max = S; max_tm = T1;
               }
            }
         }

         /* END5_3: stem (k+1..i-1) with 3' overhang at i */
         if (is_complement[numSeq1[i-1]][numSeq1[k+1]]) {
            (void) end5_candidate(3, i, k, RC, dpt, numSeq1, hend5, send5, SH);
            S = SH[0]; H = SH[1];
            if (H <= 0 && S <= 0) {
               T1 = H / (S + RC);
               G = H - temp * S;
               if ((max_tm < T1) && (G < 0.0)) {
                  H_max = H; S_max = S; max_tm = T1;
               }
            }
         }

         /* END5_4: stem (k+2..i-1) with both flanks (terminal mismatch) */
         if (is_complement[numSeq1[i-1]][numSeq1[k+2]]) {
            (void) end5_candidate(4, i, k, RC, dpt, numSeq1, hend5, send5, SH);
            S = SH[0]; H = SH[1];
            if (H <= 0 && S <= 0) {
               T1 = H / (S + RC);
               G = H - temp * S;
               if ((max_tm < T1) && (G < 0.0)) {
                  H_max = H; S_max = S; max_tm = T1;
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
traceback_monomer(int* bp, int maxLoop, const struct dpt_entry* const *dpt,
               double *send5, double *hend5, double RC,
               const unsigned char *numSeq1, int oligo1_len,
               jmp_buf _jmp_buf, thal_results* o) /* traceback for unimolecular structure */
{
   int i, j;
   i = j = 0;
   int ii, jj, k;
   struct tracer *top, *stack = NULL;
   double SH[2];
   double T1, H, S, max_tm;

   double T0 = -_INFINITY;
   push(&stack,oligo1_len, 0, 1, _jmp_buf, o);
   while(stack) {
      top = stack;
      stack = stack->next;
      i = top->i;
      j = top->j;
      if(top->mtrx==1) {
         while (equal(send5[i], send5[i - 1]) && equal(hend5[i], hend5[i - 1])) //if previous structure is the same as this one 
           --i;
         if (i == 0)
           continue;
         max_tm = (hend5[i - 1]) / (send5[i - 1] + RC);
         for (k = 0; k <= i - min_hrpn_loop - 2; ++k){
            /* END5_1: stem (k+1..i), no flanks */
            T0 = end5_candidate(1, i, k, RC, dpt, numSeq1, hend5, send5, SH);
            S = SH[0]; H = SH[1];
            if (H <= 0 && S <= 0) {
               T1 = H / (S + RC);
               if (max_tm < T1) {
                  if (equal(send5[i], S) && equal(hend5[i], H)) {
                     if (T0 >= 0.0) {
                        push(&stack, k + 1, i, 0, _jmp_buf, o);
                        push(&stack, k, 0, 1, _jmp_buf, o);
                     } else {
                        push(&stack, k + 1, i, 0, _jmp_buf, o);
                     }
                     break;
                  }
               }
            }

            /* END5_2: stem (k+2..i) with 5' overhang at k+1 */
            (void) end5_candidate(2, i, k, RC, dpt, numSeq1, hend5, send5, SH);
            S = SH[0]; H = SH[1];
            if (H <= 0 && S <= 0) {
               T1 = H / (S + RC);
               if (max_tm < T1) {
                  if (equal(send5[i], S) && equal(hend5[i], H)) {
                     if (T0 >= 0.0) {
                        push(&stack, k + 2, i, 0, _jmp_buf, o);
                        push(&stack, k, 0, 1, _jmp_buf, o);
                     } else {
                        push(&stack, k + 2, i, 0, _jmp_buf, o);
                     }
                     break;
                  }
               }
            }

            /* END5_3: stem (k+1..i-1) with 3' overhang at i */
            (void) end5_candidate(3, i, k, RC, dpt, numSeq1, hend5, send5, SH);
            S = SH[0]; H = SH[1];
            if (H <= 0 && S <= 0) {
               T1 = H / (S + RC);
               if (max_tm < T1) {
                  if (equal(send5[i], S) && equal(hend5[i], H)) {
                     if (T0 >= 0.0) {
                        push(&stack, k + 1, i - 1, 0, _jmp_buf, o);
                        push(&stack, k, 0, 1, _jmp_buf, o);
                     } else {
                        push(&stack, k + 1, i - 1, 0, _jmp_buf, o);
                     }
                     break;
                  }
               }
            }

            /* END5_4: stem (k+2..i-1) with both flanks (terminal mismatch) */
            (void) end5_candidate(4, i, k, RC, dpt, numSeq1, hend5, send5, SH);
            S = SH[0]; H = SH[1];
            if (H <= 0 && S <= 0) {
               T1 = H / (S + RC);
               if (max_tm < T1) {
                  if (equal(send5[i], S) && equal(hend5[i], H)) {
                     if (T0 >= 0.0) {
                        push(&stack, k + 2, i - 1, 0, _jmp_buf, o);
                        push(&stack, k, 0, 1, _jmp_buf, o);
                     } else {
                        push(&stack, k + 2, i - 1, 0, _jmp_buf, o);
                     }
                     break;
                  }
               }
            }
         }
      }
      if(top->mtrx==0) {
         bp[i - 1] = j;
         bp[j - 1] = i;
         SH[0] = -1.0;
         SH[1] = _INFINITY;
         if (equal(dpt[i][j].s, stackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq1[j]][numSeq1[j-1]] + dpt[i + 1][j - 1].s) &&
             equal(dpt[i][j].h, stackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq1[j]][numSeq1[j-1]] + dpt[i + 1][j - 1].h)) {
            push(&stack, i + 1, j - 1, 0, _jmp_buf, o);
         }
         else{
            //calc_hairpin(i, j, SH, 1, entropyDPT, enthalpyDPT, RC, numSeq1, oligo1_len); // 1 means that we use this method in traceback
            int d, done;
            for (done = 0, d = j - i - 3; d >= min_hrpn_loop + 1 && d >= j - i - 2 - maxLoop && !done; --d)
              for (ii = i + 1; ii < j - d; ++ii) {
                 jj = d + ii;
                 SH[0] = -1.0;
                 SH[1] = _INFINITY;
                 calc_bulge_internal_monomer(i, j, ii, jj, SH, 1, maxLoop, dpt, RC, numSeq1);
                 if (equal(dpt[i][j].s, SH[0] + dpt[ii][jj].s) &&
                     equal(dpt[i][j].h, SH[1] + dpt[ii][jj].h)) {
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

static int
thal_check_errors(const unsigned char *oligo_f, const unsigned char *oligo_r, int *len_f, int *len_r, const thal_args *a, thal_results *o){
   if (NULL == o) return 1; /* Leave it to the caller to crash */
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
   if((a->type != thal_any) && (a->type != thal_end1) && (a->type != thal_end2) && (a->type != thal_hairpin)){
      strcpy(o->msg, "Illegal type");
      return 1;
   }
   /* saltCorrectS() takes the logarithm of this value */
   if (a->mv + 120.0 * sqrt(fmax(0.0, a->dv - (a->dv <= 0 ? a->dv : a->dntp))) <= 0.0) {
      strcpy(o->msg, "Salt concentration too low for the salt correction");
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


static struct dpt_entry **
allocate_DPT(int oligo1_len, int oligo2_len, jmp_buf _jmp_buf, thal_results *o){
   //Add one to each dimension due to the way the DPT is indexed
   //Row i=0 and column j=0 are never used, but the wasted memory is negligible
   struct dpt_entry *dpt = (struct dpt_entry *) safe_malloc(sizeof(struct dpt_entry) * (oligo1_len + 1) * (oligo2_len +1), _jmp_buf, o);
   struct dpt_entry **rows = (struct dpt_entry **) safe_malloc(sizeof(struct dpt_entry *) * (oligo1_len + 1), _jmp_buf, o);
   for(int i = 0; i < oligo1_len+1; i++){
      rows[i] = &dpt[i * (oligo2_len +1)];
   }
   /* Pre-zero with a sentinel so the dump helpers can distinguish
      unvisited cells from visited ones. Production behavior unchanged
      because every cell that's read is also written first. */
   for (int i = 0; i <= oligo1_len; ++i) {
      for (int j = 0; j <= oligo2_len; ++j) {
         rows[i][j].h = NAN;
         rows[i][j].s = NAN;
         rows[i][j].tb_i = -2;  /* -2 = unvisited; -1 = visited, no predecessor */
         rows[i][j].tb_j = -2;
      }
   }
   return rows;
}

static void
free_DPT(struct dpt_entry **dpt){
   free(dpt[0]);
   free(dpt);
}

//=====================================================================================
//Functions for string manipulation
//=====================================================================================

unsigned char
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


/* get_thermodynamic_values, destroy_thal_structures, and all parameter
   file loaders live in thal_params.c. */



