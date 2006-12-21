/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006
Whitehead Institute for Biomedical Research, Steve Rozen
(http://jura.wi.mit.edu/rozen), and Helen Skaletsky
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

   * Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the
distribution.
   * Neither the names of the copyright holders nor contributors may
be used to endorse or promote products derived from this software
without specific prior written permission.

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

#include <limits.h>
#include <math.h>
#include <string.h>
#include "oligotm.h"
#include "primer3_release.h"

#define A_CHAR 'A'
#define G_CHAR 'G'
#define T_CHAR 'T'
#define C_CHAR 'C'
#define N_CHAR 'N'

#define CATID5(A,B,C,D,E) A##B##C##D##E
#define CATID2(A,B) A##B

#define DO_PAIR(LAST,THIS)          \
  if (CATID2(THIS,_CHAR) == c) {    \
     dh += CATID5(H,_,LAST,_,THIS); \
     ds += CATID5(S,_,LAST,_,THIS); \
     goto CATID2(THIS,_STATE);      \
  }

#define STATE(LAST)     \
   CATID2(LAST,_STATE): \
   c = *s; s++;         \
   DO_PAIR(LAST,A)      \
   else DO_PAIR(LAST,T) \
   else DO_PAIR(LAST,G) \
   else DO_PAIR(LAST,C) \
   else DO_PAIR(LAST,N) \
   else if ('\0' == c)  \
             goto DONE; \
   else goto ERROR \

#define DO_PAIR2(LAST,THIS)          \
     if (CATID2(THIS,_CHAR) == c) { \
      dh += CATID5(DH,_,LAST,_,THIS); \
      ds += CATID5(DS,_,LAST,_,THIS); \
      goto CATID2(THIS,_STATE2);      \
   }
#define STATE2(LAST)     \
   CATID2(LAST,_STATE2): \
   c = *s; s++;         \
   DO_PAIR2(LAST,A)      \
   else DO_PAIR2(LAST,T) \
   else DO_PAIR2(LAST,G) \
   else DO_PAIR2(LAST,C) \
   else DO_PAIR2(LAST,N) \
   else if ('\0' == c)  \
   goto DONE; \
   else goto ERROR \

double 
oligotm(s, DNA_nM, K_mM, tm_santalucia, salt_corrections)
     const  char *s;
     double DNA_nM;
     double K_mM;
     int tm_santalucia;
     int salt_corrections;
{
   register int dh = 0, ds = 0;
   register char c;
   double delta_H, delta_S;
   int len, sym;
   const char* d = s;
   len = (strlen(s)-1);
   sym=Symmetry(s); /*Add symmetry correction if seq is symmetrical*/
   if(tm_santalucia==0) {
      ds=108;
   }
   else {
      if(sym==1) {
	 ds+=14;
      }
	 
      /** Terminal AT penalty **/
      
      if(strncmp("A", s, 1)==0 || strncmp("T", s, 1)==0 || strncmp("a", s, 1)==0 || strncmp("t", s, 1)==0)  {
	   ds += -41;
	   dh += -23;
	} else if(strncmp("C", s, 1)==0 || strncmp("G", s, 1)==0 || strncmp("c", s, 1)==0 || strncmp("g", s, 1)==0) {
	   ds += 28;
	   dh += -1;
	}
      s+=len;
      if(strncmp("T", s, 1)==0 || strncmp("A", s, 1)==0 || strncmp("t", s, 1)==0 || strncmp("a", s, 1)==0) {
	 ds += -41;
	 dh += -23;
      } else if (strncmp("C", s, 1)==0 || strncmp("G", s, 1)==0 || strncmp("c", s, 1)==0 || strncmp("g", s, 1)==0) {
	 ds += 28;
	 dh += -1;
      }
      s-=len;
   }
   /* Use a finite-state machine (DFA) to calucluate dh and ds for s. */
   c = *s; s++;
   if(tm_santalucia==0) {
      if ((c == 'A')||(c=='a')) goto A_STATE;
      else if ((c == 'G')||(c=='g')) goto G_STATE;
      else if ((c == 'T')||(c=='t')) goto T_STATE;
      else if ((c == 'C') ||(c=='c')) goto C_STATE;
      else if ((c == 'N') ||(c=='n')) goto N_STATE;
      else goto ERROR;
      STATE(A);
      STATE(T);
      STATE(G);
      STATE(C);
      STATE(N);
   } else {
      if ((c == 'A')||(c=='a')) goto A_STATE2;
      else if ((c == 'G')||(c=='g')) goto G_STATE2;
      else if ((c == 'T')||(c=='t')) goto T_STATE2;
      else if ((c == 'C') ||(c=='c')) goto C_STATE2;
      else if ((c == 'N') ||(c=='n')) goto N_STATE2;
      else goto ERROR;
      STATE2(A);
      STATE2(T);
      STATE2(G);
      STATE2(C);
      STATE2(N);
   }
   
   
 DONE:  /* dh and ds are now computed for the given sequence. */
    delta_H = dh * -100.0;  /* 
			     * Nearest-neighbor thermodynamic values for dh
			     * are given in 100 cal/mol of interaction.
			     */
    delta_S = ds * -0.1;     /*
			      * Nearest-neighbor thermodynamic values for ds
			      * are in in .1 cal/K per mol of interaction.
			      */

    /* Tm calculation: tm_santalucia
     * 
     * 0 - See Rychlik, Spencer, Rhoads, "Optimization of the annealing temperature for
     *     DNA amplification in vitro."  Nucleic Acids Research, vol 18, no 21, page 6409 (1990).
     *     Article free at 
     *     http://www.pubmedcentral.nih.gov/articlerender.fcgi?tool=pubmed&pubmedid=2243783 See eqn (ii)
     *     Until Primer3 version 1.0.1 (including the version 1.0.1)
     * 1 - SantaLucia, Biochemistry vol 95, 1998 for Tm calculations
     * 
     * Salt correction formulas: salt_corrections
     * 
     * 0 - Schildkraut, C, and Lifson, S. (1965) Biopolymers, 3, 195-208. (until Primer3 version 1.0.1 (including the version 1.0.1))
     * 1 - SantaLucia98: See SantaLucia, Biochemistry vol 95, 1998
     * 2 - Owczarzy04: See Owczarzy, et al., Biochemistry vol 43, 2004
     */ 
   double Tm=0;
   len=len+1;
   if (salt_corrections==0) {
      double correction=- 273.15 + 16.6 * log10(K_mM/1000.0);
      Tm = delta_H / (delta_S + 1.987 * log(DNA_nM/4000000000.0)) + correction;
   } else if (salt_corrections==1) {
      delta_S = delta_S + 0.368 * (len - 1) * log(K_mM / 1000.0 );
      if(sym==1) { /* primer is symmetrical */
	 Tm = delta_H / (delta_S + 1.987 * log(DNA_nM/1000000000.0)) - 273.15;
      } else {
	 Tm = delta_H / (delta_S + 1.987 * log(DNA_nM/4000000000.0)) - 273.15;
      }      
   } else if (salt_corrections==2) {
      double gcPercent=0;
      int i;
      for(i=0; i<=len && d != NULL && d != '\0';) {
	 if(*d == 'C' || *d == 'G' || *d == 'c' || *d == 'g') {
	    gcPercent++;
	 }
	 *d++;
	 i++;
      }      
      gcPercent = (double)gcPercent/((double)len);
      double correction = (((4.29 * gcPercent) - 3.95) * pow(10,-5) * log(K_mM / 1000.0)) + (9.40 * pow(10,-6) * (pow(log(K_mM / 1000.0),2)));
      if(sym==1) { /* primer is symmetrical */
	 Tm = (1/((1/(delta_H / (delta_S + 1.9872 * log(DNA_nM/1000000000.0)))) + correction))-273.15;
      } else {
	 Tm = (1/((1/(delta_H / (delta_S + 1.9872 * log(DNA_nM/4000000000.0)))) + correction))-273.15;
      }
            
	}
   return Tm;
 ERROR:  /* 
	  * length of s was less than 2 or there was an illegal character in
	  * s.
	  */
    return OLIGOTM_ERROR;
}
#undef DO_PAIR
#undef DO_PAIR2

#define DO_PAIR(LAST,THIS)          \
  if (CATID2(THIS,_CHAR) == c) {    \
     dg += CATID5(G,_,LAST,_,THIS); \
     goto CATID2(THIS,_STATE);      \
  }

#define DO_PAIR2(LAST,THIS)          \
     if (CATID2(THIS,_CHAR) == c) { \
     dg += CATID5(DG,_,LAST,_,THIS); \
     goto CATID2(THIS,_STATE2);      \
}

double 
oligodg(s, tm_santalucia)
    const char *s;       /* The sequence. */
    int tm_santalucia; 
{
   register int dg = 0;
   register char c;
   /* Use a finite-state machine (DFA) to calucluate dg s. */
   c = *s; s++;
   if(tm_santalucia!=0) {      
      dg=-1960; /* Initial dG */
      if(c == 'A' || c == 'T' || c == 'a' || c == 't')  {
	 dg += -50; /* terminal AT penalty */
      }
      if ((c == 'A')||(c=='a')) goto A_STATE2;
      else if ((c == 'G')||(c=='g')) goto G_STATE2;
      else if ((c == 'T')||(c=='t')) goto T_STATE2;
      else if ((c == 'C') ||(c=='c')) goto C_STATE2;
      else if ((c == 'N') ||(c=='n')) goto N_STATE2;
      else goto ERROR;
      STATE2(A);
      STATE2(T);
      STATE2(G);
      STATE2(C);
      STATE2(N);
     } else {
    if ((c == 'A') || (c == 'a')) goto A_STATE;
    else if ((c == 'G') || (c == 'g')) goto G_STATE;
    else if ((c == 'T') || (c == 't')) goto T_STATE;
    else if ((c == 'C') || (c == 'c')) goto C_STATE;
    else if ((c == 'N') || (c == 'n')) goto N_STATE;
    else goto ERROR;
    STATE(A);
    STATE(T);
    STATE(G);
    STATE(C);
    STATE(N);

     }
DONE:  /* dg is now computed for the given sequence. */
   if(tm_santalucia!=0) {
      int sym;
      --s; --s; c = *s;
      if(c == 'A' || c == 'T' || c == 'a' || c == 't')  {
	 dg += -50; /* terminal AT penalty */
      }
      sym=Symmetry(s);
      if(sym==1)   {
	 dg +=-430; /* symmetry correction for dG */
      }
   }
   return dg / 1000.0;

 ERROR:  /* 
	  * length of s was less than 2 or there was an illegal character in
	  * s.
	  */
    return OLIGOTM_ERROR;
}

double end_oligodg(s, len, tm_santalucia)
  const char *s;  
  int len; /* The number of characters to return. */
  int tm_santalucia;
{
  int x = strlen(s);
  return x < len ? oligodg(s,tm_santalucia) : oligodg(s + (x - len),tm_santalucia);
}

double seqtm(seq, dna_conc, salt_conc, nn_max_len, tm_santalucia, salt_corrections)
  const  char *seq;
  double dna_conc;
  double salt_conc;
  int    nn_max_len;
  int tm_santalucia;
  int salt_corrections;
{
  int len = strlen(seq);
  return (len > nn_max_len)
    ? long_seq_tm(seq, 0, len, salt_conc) : oligotm(seq, dna_conc, salt_conc, tm_santalucia, salt_corrections);
}

/* See oligotm.h for documentation on this function and the formula it
   uses. */
double
long_seq_tm(s, start, len, salt_conc)
  const char *s;
  int start, len;
  double salt_conc;
{
  int GC_count = 0;
  const char *p, *end;

  if(start + len > strlen(s) || start < 0 || len <= 0) return OLIGOTM_ERROR;
  end = &s[start + len];
  /* Length <= 0 is nonsensical. */
  for (p = &s[start]; p < end; p++) {
    if ('G' == *p || 'g' == *p || 'C' == *p || 'c' == *p)
      GC_count++;
  }

  return
    81.5
    + (16.6 * log10(salt_conc / 1000.0))
    + (41.0 * (((double) GC_count) / len))
    - (600.0 / len);

}

int Symmetry(const char* seq) { /* for testing if string is symmetrical*/ 
   register char s;
   register char e;
   const char *seq_end=seq;
   int i = 0;
   int seq_len=strlen(seq);
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
      if((s=='A' && e!='T') || (s=='T' && e!='A') || (e=='A' && s!='T') || (e=='T' && s!='A') || (s=='a' && e!='t') || (s=='t' && e!='a') || (e=='a' && s!='t') || (e=='t' && s!='a')) {
	 return 0;
      }
      if((s=='C' && e!='G') || (s=='G' && e!='C') || (e=='C' && s!='G') || (e=='G' && s!='C') || (s=='c' && e!='g') || (s=='g' && e!='c') || (e=='c' && s!='g') || (e=='g' && s!='c')) {
	 return 0;
      }
      seq++;
      seq_end--;
   }
   return 1;
}
