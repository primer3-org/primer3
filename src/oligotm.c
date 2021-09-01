/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007
Whitehead Institute for Biomedical Research, Steve Rozen
(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky
All rights reserved.

    This file is part of the oligotm library.

    The oligotm library is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    The oligotm library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the oligtm library (file gpl-2.0.txt in the source
    distribution); if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

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

#define T_KELVIN 273.15;
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

/*
 * Two tables of nearest-neighbor parameters for di-nucleotide
 * base pairs.
 *
 * These are included in this file because they are not needed by
 * clients (callers) of oligtm().
 */

/* Table 1 (old parameters):
 * See table 2 in the paper [Breslauer KJ, Frank R, Bl�cker H and
 * Marky LA (1986) "Predicting DNA duplex stability from the base
 * sequence" Proc Natl Acad Sci 83:4746-50
 * http://dx.doi.org/10.1073/pnas.83.11.3746]
 */

#define S_A_A 240
#define S_A_C 173
#define S_A_G 208
#define S_A_T 239
#define S_A_N 215

#define S_C_A 129
#define S_C_C 266
#define S_C_G 278
#define S_C_T 208
#define S_C_N 220

#define S_G_A 135
#define S_G_C 267
#define S_G_G 266
#define S_G_T 173
#define S_G_N 210

#define S_T_A 169
#define S_T_C 135
#define S_T_G 129
#define S_T_T 240
#define S_T_N 168

#define S_N_A 168
#define S_N_C 210
#define S_N_G 220
#define S_N_T 215
#define S_N_N 203


#define H_A_A  91
#define H_A_C  65
#define H_A_G  78
#define H_A_T  86
#define H_A_N  80

#define H_C_A  58
#define H_C_C 110
#define H_C_G 119
#define H_C_T  78
#define H_C_N  91

#define H_G_A  56
#define H_G_C 111
#define H_G_G 110
#define H_G_T  65
#define H_G_N  85

#define H_T_A  60
#define H_T_C  56
#define H_T_G  58
#define H_T_T  91
#define H_T_N  66

#define H_N_A  66
#define H_N_C  85
#define H_N_G  91
#define H_N_T  80
#define H_N_N  80

/* Delta G's of disruption * 1000. */
#define G_A_A  1900
#define G_A_C  1300
#define G_A_G  1600
#define G_A_T  1500
#define G_A_N  1575

#define G_C_A  1900
#define G_C_C  3100
#define G_C_G  3600
#define G_C_T  1600
#define G_C_N  2550

#define G_G_A  1600
#define G_G_C  3100
#define G_G_G  3100
#define G_G_T  1300
#define G_G_N  2275

#define G_T_A   900
#define G_T_C  1600
#define G_T_G  1900
#define G_T_T  1900
#define G_T_N  1575

#define G_N_A  1575
#define G_N_C  2275
#define G_N_G  2550
#define G_N_T  1575
#define G_N_N  1994

/* Table 2, new parameters:
 * Tables of nearest-neighbor thermodynamics for DNA bases, from the
 * paper [SantaLucia JR (1998) "A unified view of polymer, dumbbell
 * and oligonucleotide DNA nearest-neighbor thermodynamics", Proc Natl
 * Acad Sci 95:1460-65 http://dx.doi.org/10.1073/pnas.95.4.1460]
 */

#define DS_A_A 222
#define DS_A_C 224
#define DS_A_G 210
#define DS_A_T 204
#define DS_A_N 224

#define DS_C_A 227
#define DS_C_C 199
#define DS_C_G 272
#define DS_C_T 210
#define DS_C_N 272

#define DS_G_A 222
#define DS_G_C 244
#define DS_G_G 199
#define DS_G_T 224
#define DS_G_N 244

#define DS_T_A 213
#define DS_T_C 222
#define DS_T_G 227
#define DS_T_T 222
#define DS_T_N 227

#define DS_N_A 168
#define DS_N_C 210
#define DS_N_G 220
#define DS_N_T 215
#define DS_N_N 220


#define DH_A_A  79
#define DH_A_C  84
#define DH_A_G  78
#define DH_A_T  72
#define DH_A_N  72

#define DH_C_A  85
#define DH_C_C  80
#define DH_C_G 106
#define DH_C_T  78
#define DH_C_N  78

#define DH_G_A  82
#define DH_G_C  98
#define DH_G_G  80
#define DH_G_T  84
#define DH_G_N  80

#define DH_T_A  72
#define DH_T_C  82
#define DH_T_G  85
#define DH_T_T  79
#define DH_T_N  72

#define DH_N_A  72
#define DH_N_C  80
#define DH_N_G  78
#define DH_N_T  72
#define DH_N_N  72

/* Delta G's of disruption * 1000. */
#define DG_A_A  1000
#define DG_A_C  1440
#define DG_A_G  1280
#define DG_A_T  880
#define DG_A_N  880

#define DG_C_A  1450
#define DG_C_C  1840
#define DG_C_G  2170
#define DG_C_T  1280
#define DG_C_N  1450

#define DG_G_A  1300
#define DG_G_C  2240
#define DG_G_G  1840
#define DG_G_T  1440
#define DG_G_N  1300

#define DG_T_A   580
#define DG_T_C  1300
#define DG_T_G  1450
#define DG_T_T  1000
#define DG_T_N   580

#define DG_N_A   580
#define DG_N_C  1300
#define DG_N_G  1280
#define DG_N_T   880
#define DG_N_N   580

/* End of tables nearest-neighbor parameter. */


/* Calculate the melting temperature of oligo s.  See
   oligotm.h for documentation of arguments.
*/
double 
oligotm(const  char *s,
     double DNA_nM,
     double K_mM,
     double divalent_conc,
     double dntp_conc,
     tm_method_type  tm_method,
     salt_correction_type salt_corrections)
{
   register int dh = 0, ds = 0;
  register char c;
  double delta_H, delta_S;
  double Tm; /* Melting temperature */
  double correction;
  int len, sym;
  const char* d = s;
   if(divalent_to_monovalent(divalent_conc, dntp_conc) == OLIGOTM_ERROR)
     return OLIGOTM_ERROR;
   
  /** K_mM = K_mM + divalent_to_monovalent(divalent_conc, dntp_conc); **/
  if (tm_method != breslauer_auto
      && tm_method != santalucia_auto)
    return OLIGOTM_ERROR;
  if (salt_corrections != schildkraut
      && salt_corrections != santalucia
      && salt_corrections != owczarzy)
    return OLIGOTM_ERROR;

  len = (strlen(s)-1);

  sym = symmetry(s); /*Add symmetry correction if seq is symmetrical*/
  if( tm_method == breslauer_auto ) {
    ds=108;
  }
  else {
    if(sym == 1) {
      ds+=14;
    }

    /** Terminal AT penalty **/
      
    if(strncmp("A", s, 1)==0
       || strncmp("T", s, 1)==0)  {
      ds += -41;
      dh += -23;
    } else if (strncmp("C", s, 1)==0 
               || strncmp("G", s, 1)==0) {
      ds += 28;
      dh += -1;
    }
    s+=len;
    if(strncmp("T", s, 1)==0 
       || strncmp("A", s, 1)==0) {
      ds += -41;
      dh += -23;
    } else if (strncmp("C", s, 1)==0 
               || strncmp("G", s, 1)==0) {
      ds += 28;
      dh += -1;
    }
    s-=len;
  }
  /* Use a finite-state machine (DFA) to calucluate dh and ds for s. */
  c = *s; s++;
  if (tm_method == breslauer_auto) {
    if (c == 'A') goto A_STATE;
    else if (c == 'G') goto G_STATE;
    else if (c == 'T') goto T_STATE;
    else if (c == 'C') goto C_STATE;
    else if (c == 'N') goto N_STATE;
    else goto ERROR;
    STATE(A);
    STATE(T);
    STATE(G);
    STATE(C);
    STATE(N);
  } else {
    if (c == 'A') goto A_STATE2;
    else if (c == 'G') goto G_STATE2;
    else if (c == 'T') goto T_STATE2;
    else if (c == 'C') goto C_STATE2;
    else if (c == 'N') goto N_STATE2;
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
  Tm=0;  /* Melting temperature */
  len=len+1;

   /**********************************************/
  if (salt_corrections == schildkraut) {
     K_mM = K_mM + divalent_to_monovalent(divalent_conc, dntp_conc);
     correction = 16.6 * log10(K_mM/1000.0) - T_KELVIN;
     Tm = delta_H / (delta_S + 1.987 * log(DNA_nM/4000000000.0)) + correction;
  } else if (salt_corrections== santalucia) {
    K_mM = K_mM + divalent_to_monovalent(divalent_conc, dntp_conc);
    delta_S = delta_S + 0.368 * (len - 1) * log(K_mM / 1000.0 );
    if(sym == 1) { /* primer is symmetrical */
      /* Equation A */
      Tm = delta_H / (delta_S + 1.987 * log(DNA_nM/1000000000.0)) - T_KELVIN;
    } else {
      /* Equation B */
      Tm = delta_H / (delta_S + 1.987 * log(DNA_nM/4000000000.0)) - T_KELVIN;
    }      
  } else if (salt_corrections == owczarzy) {
     double gcPercent=0;
     double free_divalent; /* conc of divalent cations minus dNTP conc */
     int i;
     for(i = 0; i <= len && d != NULL && *d != '\0';) {
        if(*d == 'C' || *d == 'G') {
           gcPercent++;
        }  
        d++;
        i++;
     }
     gcPercent = (double)gcPercent/((double)len);
     /**** BEGIN: UPDATED SALT BY OWCZARZY *****/
      /* different salt corrections for monovalent (Owczarzy et al.,2004) 
      and divalent cations (Owczarzy et al.,2008)
      */
     /* competition bw magnesium and monovalent cations, see Owczarzy et al., 2008 Figure 9 and Equation 16 */
     
     static const double crossover_point = 0.22; /* depending on the value of div_monov_ratio respect
                                             to value of crossover_point Eq 16 (divalent corr, Owczarzy et al., 2008)
                                             or Eq 22 (monovalent corr, Owczarzy et al., 2004) should be used */
     double div_monov_ratio;
     if(dntp_conc >= divalent_conc) {
        free_divalent = 0.00000000001; /* to not to get log(0) */
     } else {
        free_divalent = (divalent_conc - dntp_conc)/1000.0;
     }
     static double a = 0,b = 0,c = 0,d = 0,e = 0,f = 0,g = 0;
      if(K_mM==0) {
         div_monov_ratio = 6.0;
      } else {
         div_monov_ratio = (sqrt(free_divalent))/(K_mM/1000); /* if conc of monov cations is provided
                                                                    a ratio is calculated to further calculate
                                                                    the _correct_ correction */
      }
     if (div_monov_ratio < crossover_point) {
        /* use only monovalent salt correction, Eq 22 (Owczarzy et al., 2004) */
        correction
          = (((4.29 * gcPercent) - 3.95) * pow(10,-5) * log(K_mM / 1000.0))
            + (9.40 * pow(10,-6) * (pow(log(K_mM / 1000.0),2)));
     } else {
        /* magnesium effects are dominant, Eq 16 (Owczarzy et al., 2008) is used */
        b =- 9.11 * pow(10,-6);
        c = 6.26 * pow(10,-5);
        e =- 4.82 * pow(10,-4);
        f = 5.25 * pow(10,-4);
        a = 3.92 * pow(10,-5);
        d = 1.42 * pow(10,-5);
        g = 8.31 * pow(10,-5);
        if(div_monov_ratio < 6.0) {
           /* in particular ratio of conc of monov and div cations
            *             some parameters of Eq 16 must be corrected (a,d,g) */
           a = 3.92 * pow(10,-5) * (0.843 - (0.352 * sqrt(K_mM/1000.0) * log(K_mM/1000.0)));
           d = 1.42 * pow(10,-5) * (1.279 - 4.03 * pow(10,-3) * log(K_mM/1000.0) - 8.03 * pow(10,-3) * pow(log(K_mM/1000.0),2));
           g = 8.31 * pow(10,-5) * (0.486 - 0.258 * log(K_mM/1000.0) + 5.25 * pow(10,-3) * pow(log(K_mM/1000.0),3));
        }
        
        correction = a + (b * log(free_divalent))
          + gcPercent * (c + (d * log(free_divalent)))
            + (1/(2 * (len - 1))) * (e + (f * log(free_divalent))
                                     + g * (pow((log(free_divalent)),2)));
     }
     /**** END: UPDATED SALT BY OWCZARZY *****/
     if (sym == 1) {
           /* primer is symmetrical */
        /* Equation A */
        Tm = 1/((1/(delta_H 
                    /
                    (delta_S + 1.9872 * log(DNA_nM/1000000000.0)))) + correction) - T_KELVIN;
     } else {
        /* Equation B */
        Tm = 1/((1/(delta_H
                    /
                    (delta_S + 1.9872 * log(DNA_nM/4000000000.0)))) + correction) - T_KELVIN;
     } 
  } /* END else if (salt_corrections == owczarzy) { */
   
   
   /***************************************/
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
oligodg(const char *s,      /* The sequence. */
    int tm_method) 
{
   register int dg = 0;
   register char c;

  if (tm_method != breslauer_auto
      && tm_method != santalucia_auto)
    return OLIGOTM_ERROR;

   /* Use a finite-state machine (DFA) to calucluate dg s. */
   c = *s; s++;
   if(tm_method != breslauer_auto) {      
      dg=-1960; /* Initial dG */
      if(c == 'A' || c == 'T')  {
         dg += -50; /* terminal AT penalty */
      }
      if (c == 'A') goto A_STATE2;
      else if (c == 'G') goto G_STATE2;
      else if (c == 'T') goto T_STATE2;
      else if (c == 'C') goto C_STATE2;
      else if (c == 'N') goto N_STATE2;
      else goto ERROR;
      STATE2(A);
      STATE2(T);
      STATE2(G);
      STATE2(C);
      STATE2(N);
     } else {
    if (c == 'A') goto A_STATE;
    else if (c == 'G') goto G_STATE;
    else if (c == 'T') goto T_STATE;
    else if (c == 'C') goto C_STATE;
    else if (c == 'N') goto N_STATE;
    else goto ERROR;
    STATE(A);
    STATE(T);
    STATE(G);
    STATE(C);
    STATE(N);

     }
DONE:  /* dg is now computed for the given sequence. */
   if(tm_method != breslauer_auto) {
      int sym;
      --s; --s; c = *s;
      if(c == 'A' || c == 'T')  {
         dg += -50; /* terminal AT penalty */
      }
      sym = symmetry(s);
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

double end_oligodg(const char *s,  
  size_t len, /* The number of characters to return. */
  int tm_method)
{
  int x = strlen(s);

  if (tm_method != breslauer_auto
      && tm_method != santalucia_auto)
    return OLIGOTM_ERROR;

  return 
    x < len 
    ? oligodg(s,tm_method) :
    oligodg(s + (x - len),tm_method);
}

/* See oligotm.h for documentation of arguments. */
double seqtm(const  char *seq,
  double dna_conc,
  double salt_conc,
  double divalent_conc,
  double dntp_conc,
  int    nn_max_len,
  tm_method_type tm_method,
  salt_correction_type salt_corrections)
{
  int len = strlen(seq);
   if (tm_method != breslauer_auto
      && tm_method != santalucia_auto)
    return OLIGOTM_ERROR;
  if (salt_corrections != schildkraut
      && salt_corrections != santalucia
      && salt_corrections != owczarzy)
    return OLIGOTM_ERROR;

  if (len > nn_max_len) {
          return long_seq_tm(seq, 0, len, salt_conc, divalent_conc, dntp_conc);
  } else {
          return oligotm(seq, dna_conc, salt_conc, 
                      divalent_conc, dntp_conc, tm_method, salt_corrections);
  }
}

/* See oligotm.h for documentation on this function and the formula it
   uses. */
double
long_seq_tm(const char *s,
  int start,
  int len,
  double salt_conc,
  double divalent_conc,
  double dntp_conc)
{
  int GC_count = 0;
  const char *p, *end;

  if (divalent_to_monovalent(divalent_conc, dntp_conc) == OLIGOTM_ERROR)
    return OLIGOTM_ERROR;
  
  salt_conc = salt_conc + divalent_to_monovalent(divalent_conc, dntp_conc);
  
  if ((unsigned) (start + len) > strlen(s) || start < 0 || len <= 0)
    return OLIGOTM_ERROR;
  end = &s[start + len];
  /* Length <= 0 is nonsensical. */
  for (p = &s[start]; p < end; p++) {
    if ('G' == *p || 'C' == *p)
      GC_count++;
  }
  
  return
    81.5
    + (16.6 * log10(salt_conc / 1000.0))
    + (41.0 * (((double) GC_count) / len))
    - (600.0 / len);
}

 /* Return 1 if string is symmetrical, 0 otherwise. */ 
int symmetry(const char* seq) {
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

/* Convert divalent salt concentration to monovalent */
double divalent_to_monovalent(double divalent, 
                              double dntp)
{
   if(divalent==0) dntp=0;
   if(divalent<0 || dntp<0) return OLIGOTM_ERROR;
   if(divalent<dntp) 
     /* According to theory, melting temperature does not depend on
        divalent cations */
     divalent=dntp;  
   return 120*(sqrt(divalent-dntp));
}
