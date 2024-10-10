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
/* #include <stdio.h> */ /* Only for testing */

#ifndef MAX_PRIMER_LENGTH
#define MAX_PRIMER_LENGTH 36
#endif
static const double t_kelvin = 273.15;
/*
 * Two tables of nearest-neighbor parameters for di-nucleotide
 * base pairs.
 *
 * These are included in this file because they are not needed by
 * clients (callers) of oligtm().
 */

/* Table 1 (old parameters):
 * See table 2 in the paper [Breslauer KJ, Frank R, Bloecker H and
 * Marky LA (1986) "Predicting DNA duplex stability from the base
 * sequence" Proc Natl Acad Sci 83:4746-50
 * http://dx.doi.org/10.1073/pnas.83.11.3746]
 */

/* dH *-100 cal/mol */
static const int Breslauer_1986_dH[5][5] =  {
     {91,  65,  78, 86, 80},  /* AA, AC, AG, AT, AN; */
     {58, 110, 119, 78, 91},  /* CA, CC, CG, CT, CN; */
     {56, 111, 110, 65, 85},  /* GA, GC, GG, GT, GN; */
     {60,  56,  58, 91, 66},  /* TA, TC, TG, TT, TN; */
     {66,  85,  91, 80, 80}}; /* NA, NC, NG, NT, NN; */

/* dS *-0.1 cal/k*mol */
static const int Breslauer_1986_dS[5][5] =  {
     {240, 173, 208, 239, 215},  /* AA, AC, AG, AT, AN; */
     {129, 266, 278, 208, 220},  /* CA, CC, CG, CT, CN; */
     {135, 267, 266, 173, 210},  /* GA, GC, GG, GT, GN; */
     {169, 135, 129, 240, 168},  /* TA, TC, TG, TT, TN; */
     {168, 210, 220, 215, 203}}; /* NA, NC, NG, NT, NN; */

/* dG *-0.001 cal/mol */
static const int Breslauer_1986_dG[5][5] =  {
     {1900, 1300, 1600, 1500, 1575},  /* AA, AC, AG, AT, AN; */
     {1900, 3100, 3600, 1600, 2550},  /* CA, CC, CG, CT, CN; */
     {1600, 3100, 3100, 1300, 2275},  /* GA, GC, GG, GT, GN; */
     { 900, 1600, 1900, 1900, 1575},  /* TA, TC, TG, TT, TN; */
     {1575, 2275, 2550, 1575, 1994}}; /* NA, NC, NG, NT, NN; */


/* Table 2, new parameters:
 * Tables of nearest-neighbor thermodynamics for DNA bases, from the
 * paper [SantaLucia JR (1998) "A unified view of polymer, dumbbell
 * and oligonucleotide DNA nearest-neighbor thermodynamics", Proc Natl
 * Acad Sci 95:1460-65 http://dx.doi.org/10.1073/pnas.95.4.1460]
 */

/* dH *-100 cal/mol */
static const int SantaLucia_1998_dH[5][5] =  {
     {79, 84,  78, 72, 72},  /* AA, AC, AG, AT, AN; */
     {85, 80, 106, 78, 78},  /* CA, CC, CG, CT, CN; */
     {82, 98,  80, 84, 80},  /* GA, GC, GG, GT, GN; */
     {72, 82,  85, 79, 72},  /* TA, TC, TG, TT, TN; */
     {72, 80,  78, 72, 72}}; /* NA, NC, NG, NT, NN; */

/* dS *-0.1 cal/k*mol */
static const int SantaLucia_1998_dS[5][5] =  {
     {222, 224, 210, 204, 224},  /* AA, AC, AG, AT, AN; */
     {227, 199, 272, 210, 272},  /* CA, CC, CG, CT, CN; */
     {222, 244, 199, 224, 244},  /* GA, GC, GG, GT, GN; */
     {213, 222, 227, 222, 227},  /* TA, TC, TG, TT, TN; */
     {168, 210, 220, 215, 220}}; /* NA, NC, NG, NT, NN; */

/* dG *-0.001 cal/mol */
static const int SantaLucia_1998_dG[5][5] =  {
     {1000, 1440, 1280,  880,  880},  /* AA, AC, AG, AT, AN; */
     {1450, 1840, 2170, 1280, 1450},  /* CA, CC, CG, CT, CN; */
     {1300, 2240, 1840, 1440, 1300},  /* GA, GC, GG, GT, GN; */
     { 580, 1300, 1450, 1000,  580},  /* TA, TC, TG, TT, TN; */
     { 580, 1300, 1280,  880,  580}}; /* NA, NC, NG, NT, NN; */


/* Table 3, updated parameters:
 * Tables of nearest-neighbor thermodynamics for DNA bases, from the
 * paper [SantaLucia JR and Hicks (2004) "The thermodynamics of DNA
 * Structuralmotifs", PAnnu. Rev. Biophys. Biomol. Struct. 2004. 33:
 * 415–40 http://dx.doi.org/10.1146/annurev.biophys.32.110601.141800]
 */

/* dH *-100 cal/mol */
static const int SantaLucia_2004_dH[5][5] =  {
     {76, 84,  78, 72, 72},  /* AA, AC, AG, AT, AN; */
     {85, 80, 106, 78, 78},  /* CA, CC, CG, CT, CN; */
     {82, 98,  80, 84, 80},  /* GA, GC, GG, GT, GN; */
     {72, 82,  85, 76, 72},  /* TA, TC, TG, TT, TN; */
     {72, 80,  78, 72, 72}}; /* NA, NC, NG, NT, NN; */

/* dS *-0.1 cal/k*mol */
static const int SantaLucia_2004_dS[5][5] =  {
     {213, 224, 210, 204, 224},  /* AA, AC, AG, AT, AN; */
     {227, 199, 272, 210, 272},  /* CA, CC, CG, CT, CN; */
     {222, 244, 199, 224, 244},  /* GA, GC, GG, GT, GN; */
     {213, 222, 227, 213, 227},  /* TA, TC, TG, TT, TN; */
     {168, 210, 220, 215, 220}}; /* NA, NC, NG, NT, NN; */

/* dG *-0.001 cal/mol */
/* unmodified to SantaLucia_1998_dG */

/* Calculate the melting temperature of oligo s.  See
   oligotm.h for documentation of arguments.
*/
tm_ret
oligotm(const  char *s,
        double DNA_nM,
        double K_mM,
        double divalent_conc,
        double dntp_conc,
        double dmso_conc,
        double dmso_fact,
        double formamide_conc,
        tm_method_type  tm_method,
        salt_correction_type salt_corrections,
        double annealing_temp)
{
  int oligo_int[MAX_PRIMER_LENGTH+1];
  register int dh = 0, ds = 0;
  double delta_H = 0.0;
  double delta_S = 0.0;
  double ddG, ka;
  tm_ret ret;
  ret.bound = OLIGOTM_ERROR;
  double correction;
  int len, sym, i;
  int GC_count = 0;

  if(divalent_to_monovalent(divalent_conc, dntp_conc) == OLIGOTM_ERROR) {
    ret.Tm = OLIGOTM_ERROR;
    return ret;
  }
  /** K_mM = K_mM + divalent_to_monovalent(divalent_conc, dntp_conc); **/
  if (tm_method != breslauer_auto && tm_method != santalucia_auto
      && tm_method != santalucia_2004) {
    ret.Tm = OLIGOTM_ERROR;
    return ret;
  }
  if (salt_corrections != schildkraut
      && salt_corrections != santalucia
      && salt_corrections != owczarzy) {
    ret.Tm = OLIGOTM_ERROR;
    return ret;
  }

  /* Translate to int array */
  len = (strlen(s)-1);
  for (i = 0; i < len + 1; i++) {
    switch (s[i]) {
      case 'A':
        oligo_int[i] = 0;
        break;
      case 'C':
        oligo_int[i] = 1;
        GC_count++;
        break;
      case 'G':
        oligo_int[i] = 2;
        GC_count++;
        break;
      case 'T':
        oligo_int[i] = 3;
        break;
      case 'N':
        oligo_int[i] = 4;
        break;
      default:
        ret.Tm = OLIGOTM_ERROR;
        return ret;
    }
  }

  sym = symmetry(s); /*Add symmetry correction if seq is symmetrical*/
  if(tm_method == santalucia_auto) {
    if(sym == 1) {
      ds += 14;
    }

    /** Terminal penalty **/
    switch (s[0]) {
      case 'A':
      case 'T':
        ds += -41;
        dh += -23;
        break;
      case 'C':
      case 'G':
        ds += 28;
        dh += -1;
        break;
    }
    /* Sum the pairs up */
    for (i = 0; i < len; i++) {
      ds += SantaLucia_1998_dS[oligo_int[i]][oligo_int[i+1]];
      dh += SantaLucia_1998_dH[oligo_int[i]][oligo_int[i+1]];
    }
    /** Terminal penalty **/
    switch (s[len]) {
      case 'A':
      case 'T':
        ds += -41;
        dh += -23;
        break;
      case 'C':
      case 'G':
        ds += 28;
        dh += -1;
        break;
    }
  } else if(tm_method == santalucia_2004) {
    ds += 57;
    dh += -2;

    if(sym == 1) {
      ds += 14;
    }

    /** Terminal penalty **/
    switch (s[0]) {
      case 'A':
      case 'T':
        ds += -69;
        dh += -22;
        break;
    }
    /* Sum the pairs up */
    for (i = 0; i < len; i++) {
      ds += SantaLucia_2004_dS[oligo_int[i]][oligo_int[i+1]];
      dh += SantaLucia_2004_dH[oligo_int[i]][oligo_int[i+1]];
    }
    /** Terminal penalty **/
    switch (s[len]) {
      case 'A':
      case 'T':
        ds += -69;
        dh += -22;
        break;
    }
  } else {
    ds = 108;
    /* Sum the pairs up */
    for (i = 0; i < len; i++) {
      ds += Breslauer_1986_dS[oligo_int[i]][oligo_int[i+1]];
      dh += Breslauer_1986_dH[oligo_int[i]][oligo_int[i+1]];
    }
  }

  delta_H = dh * -100.0;  /*
                           * Nearest-neighbor thermodynamic values for dh
                           * are given in 100 cal/mol of interaction.
                           */
  delta_S = ds * -0.1;     /*
                            * Nearest-neighbor thermodynamic values for ds
                            * are in in .1 cal/K per mol of interaction.
                            */
  ret.Tm=0;  /* Melting temperature */
  len=len+1;

  /* For testing */
  /* printf("dH: %.2f kcal/mol\ndS: %.2f cal/k*mol\n", delta_H / 1000.0, delta_S); */

   /**********************************************/
  if (salt_corrections == schildkraut) {
    K_mM = K_mM + divalent_to_monovalent(divalent_conc, dntp_conc);
    correction = 16.6 * log10(K_mM/1000.0);
    ret.Tm = delta_H / (delta_S + 1.987 * log(DNA_nM/4000000000.0)) + correction - t_kelvin;
    ret.Tm -= dmso_conc * dmso_fact;
    ret.Tm += (0.453 * ((double) GC_count) / len - 2.88) * formamide_conc;
    if (annealing_temp > 0.0) {
      ddG = delta_H - (annealing_temp - correction + t_kelvin) * delta_S;
      ka = exp(-ddG / (1.987 * (annealing_temp - correction + t_kelvin)));
      ret.bound = (1 / (1 + sqrt(1/((DNA_nM/4000000000.0) * ka)))) * 100;
    }
  } else if (salt_corrections == santalucia) {
    K_mM = K_mM + divalent_to_monovalent(divalent_conc, dntp_conc);
    delta_S = delta_S + 0.368 * (len - 1) * log(K_mM / 1000.0 );
    if(sym == 1) { /* primer is symmetrical */
      /* Equation A */
      ret.Tm = delta_H / (delta_S + 1.987 * log(DNA_nM/1000000000.0)) - t_kelvin;
      ret.Tm -= dmso_conc * dmso_fact;
      ret.Tm += (0.453 * ((double) GC_count) / len - 2.88) * formamide_conc;
      if (annealing_temp > 0.0) {
        ddG = delta_H - (annealing_temp + t_kelvin) * delta_S;
        ka = exp(-ddG / (1.987 * (annealing_temp + t_kelvin)));
        ret.bound = (1 / (1 + sqrt(1/((DNA_nM/1000000000.0) * ka)))) * 100;
      }
    } else {
      /* Equation B */
      ret.Tm = delta_H / (delta_S + 1.987 * log(DNA_nM/4000000000.0)) - t_kelvin;
      ret.Tm -= dmso_conc * dmso_fact;
      ret.Tm += (0.453 * ((double) GC_count) / len - 2.88) * formamide_conc;
      if (annealing_temp > 0.0) {
        ddG = delta_H - (annealing_temp + t_kelvin) * delta_S;
        ka = exp(-ddG / (1.987 * (annealing_temp + t_kelvin)));
        ret.bound = (1 / (1 + sqrt(1/((DNA_nM/4000000000.0) * ka)))) * 100;
      }
    }
  } else if (salt_corrections == owczarzy) {
    double gcFract;
    double free_divalent; /* conc of divalent cations minus dNTP conc */
    gcFract = (double)GC_count/((double)len);
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
      correction = (((4.29 * gcFract) - 3.95) * pow(10,-5) * log(K_mM / 1000.0))
                   + (9.40 * pow(10,-6) * (pow(log(K_mM / 1000.0),2)));
    } else {
      /* magnesium effects are dominant, Eq 16 (Owczarzy et al., 2008) is used */
      a =  3.92 * pow(10,-5);
      b = -9.11 * pow(10,-6);
      c =  6.26 * pow(10,-5);
      d =  1.42 * pow(10,-5);
      e = -4.82 * pow(10,-4);
      f =  5.25 * pow(10,-4);
      g =  8.31 * pow(10,-5);
      if(div_monov_ratio < 6.0) {
        /* in particular ratio of conc of monov and div cations
         *             some parameters of Eq 16 must be corrected (a,d,g) */
        a = 3.92 * pow(10,-5) * (0.843 - (0.352 * sqrt(K_mM/1000.0) * log(K_mM/1000.0)));
        d = 1.42 * pow(10,-5) * (1.279 - 4.03 * pow(10,-3) * log(K_mM/1000.0) - 8.03 * pow(10,-3) * pow(log(K_mM/1000.0),2));
        g = 8.31 * pow(10,-5) * (0.486 - 0.258 * log(K_mM/1000.0) + 5.25 * pow(10,-3) * pow(log(K_mM/1000.0),3));
      }

      correction = a + (b * log(free_divalent))
                     + gcFract * (c + (d * log(free_divalent)))
                     + (1/(2 * (((double) len)- 1))) * (e + (f * log(free_divalent)) + g * (pow((log(free_divalent)),2)));
    }
    /**** END: UPDATED SALT BY OWCZARZY *****/
    if (sym == 1) {
      /* primer is symmetrical */
      /* Equation A */
      ret.Tm = 1/((1/(delta_H
                    /
                    (delta_S + 1.9872 * log(DNA_nM/1000000000.0)))) + correction) - t_kelvin;
    } else {
      /* Equation B */
      ret.Tm = 1/((1/(delta_H
                      /
                      (delta_S + 1.9872 * log(DNA_nM/4000000000.0)))) + correction) - t_kelvin;
    }
    ret.Tm -= dmso_conc * dmso_fact;
    ret.Tm += (0.453 * ((double) GC_count) / ((double) len) - 2.88) * formamide_conc;
  } /* END else if (salt_corrections == owczarzy) { */
  /***************************************/
  return ret;
}

double
oligodg(const char *s,      /* The sequence. */
        int tm_method)
{
  register int dg = 0;
  int oligo_int[MAX_PRIMER_LENGTH+1];
  int len, i, sym;

  if (tm_method != breslauer_auto
      && tm_method != santalucia_auto)
    return OLIGOTM_ERROR;

  /* Translate to int array */
  len = (strlen(s)-1);
  for (i = 0; i < len + 1; i++) {
    switch (s[i]) {
      case 'A':
        oligo_int[i] = 0;
        break;
      case 'C':
        oligo_int[i] = 1;
        break;
      case 'G':
        oligo_int[i] = 2;
        break;
      case 'T':
        oligo_int[i] = 3;
        break;
      case 'N':
        oligo_int[i] = 4;
        break;
      default:
        return OLIGOTM_ERROR;
    }
  }

  if ((tm_method == santalucia_auto) || (tm_method == santalucia_2004)) {
    dg  = -1960; /* Initial dG */

    sym = symmetry(s);
    if(sym == 1) {
      dg += -430; /* symmetry correction for dG */
    }

    /** Terminal penalty **/
    switch (s[0]) {
      case 'A':
      case 'T':
        dg += -50; /* terminal AT penalty */
        break;
    }
     /* Sum the pairs up */
    for (i = 0; i < len; i++) {
      dg += SantaLucia_1998_dG[oligo_int[i]][oligo_int[i+1]];
    }
    /** Terminal penalty **/
    switch (s[len]) {
      case 'A':
      case 'T':
        dg += -50; /* terminal AT penalty */
        break;
    }
  } else {
    /* Sum the pairs up */
    for (i = 0; i < len; i++) {
      dg += Breslauer_1986_dG[oligo_int[i]][oligo_int[i+1]];
    }
  }

  return dg / 1000.0;
}

double end_oligodg(const char *s,
  int len, /* The number of characters to return. */
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
tm_ret seqtm(const  char *seq,
             double dna_conc,
             double salt_conc,
             double divalent_conc,
             double dntp_conc,
             double dmso_conc,
             double dmso_fact,
             double formamide_conc,
             int    nn_max_len,
             tm_method_type tm_method,
             salt_correction_type salt_corrections,
             double annealing_temp)
{
  int len = strlen(seq);
  tm_ret ret;
  ret.bound = OLIGOTM_ERROR;
  ret.Tm = OLIGOTM_ERROR;
  if (tm_method != breslauer_auto
      && tm_method != santalucia_auto
      && tm_method != santalucia_2004)
    return ret;
  if (salt_corrections != schildkraut
      && salt_corrections != santalucia
      && salt_corrections != owczarzy)
    return ret;

  if (len > nn_max_len) {
    return long_seq_tm(seq, 0, len, salt_conc, divalent_conc, dntp_conc,
                       dmso_conc, dmso_fact, formamide_conc);
  } else {
    return oligotm(seq, dna_conc, salt_conc, divalent_conc, dntp_conc, dmso_conc,
                   dmso_fact, formamide_conc, tm_method, salt_corrections, annealing_temp);
  }
}

/* See oligotm.h for documentation on this function and the formula it
   uses. */
tm_ret
long_seq_tm(const char *s,
            int start,
            int len,
            double salt_conc,
            double divalent_conc,
            double dntp_conc,
            double dmso_conc,
            double dmso_fact,
            double formamide_conc)
{
  int GC_count = 0;
  const char *p, *end;
  tm_ret ret;
  ret.bound = OLIGOTM_ERROR;
  ret.Tm = OLIGOTM_ERROR;

  if (divalent_to_monovalent(divalent_conc, dntp_conc) == OLIGOTM_ERROR)
    return ret;

  salt_conc = salt_conc + divalent_to_monovalent(divalent_conc, dntp_conc);

  if ((unsigned) (start + len) > strlen(s) || start < 0 || len <= 0)
    return ret;
  end = &s[start + len];
  /* Length <= 0 is nonsensical. */
  for (p = &s[start]; p < end; p++) {
    if ('G' == *p || 'C' == *p)
      GC_count++;
  }

  ret.Tm = 81.5 - dmso_conc * dmso_fact
                + (0.453 * ((double) GC_count) / ((double) len) - 2.88) * formamide_conc
                + (16.6 * log10(salt_conc / 1000.0))
                + (41.0 * (((double) GC_count) / ((double) len)))
                - (600.0 / ((double) len));

  return ret;
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
