/*
 * Copyright (c) 1996, Whitehead Institute for Biomedical Research. All rights
 * reserved.  Please see full software use agreement in primer3_main.c or by
 * executing primer3 with -h.
 */

#include <limits.h>
#include <math.h>
#include <string.h>
#include "oligotm.h"
#include "primer3_release.h"

/* 
 * Tables of nearest-neighbor thermodynamics for DNA bases.  See Breslauer,
 * Frank, Bloecker, and Markey, Proc. Natl. Acad. Sci. USA, vol 83, page 3748,
 * table 2.
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

double 
oligotm(s, DNA_nM, K_mM)
     const  char *s;
     double DNA_nM;
     double K_mM;
{
    register int dh = 0, ds = 108;
    register char c;
    double delta_H, delta_S;

    /* Use a finite-state machine (DFA) to calucluate dh and ds for s. */
    c = *s; s++;
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

 DONE:  /* dh and ds are now computed for the given sequence. */
    delta_H = dh * -100.0;  /* 
			     * Nearest-neighbor thermodynamic values for dh
			     * are given in 100 cal/mol of interaction.
			     */
    delta_S = ds * -0.1;     /*
			      * Nearest-neighbor thermodynamic values for ds
			      * are in in .1 cal/K per mol of interaction.
			      */

    /* 
     * See Rychlik, Spencer, Roads, Nucleic Acids Research, vol 18, no 21,
     * page 6410, eqn (ii).
     */
    return delta_H / (delta_S + 1.987 * log(DNA_nM/4000000000.0))
	- 273.15 + 16.6 * log10(K_mM/1000.0);

 ERROR:  /* 
	  * length of s was less than 2 or there was an illegal character in
	  * s.
	  */
    return OLIGOTM_ERROR;
}
#undef DO_PAIR

#define DO_PAIR(LAST,THIS)          \
  if (CATID2(THIS,_CHAR) == c) {    \
     dg += CATID5(G,_,LAST,_,THIS); \
     goto CATID2(THIS,_STATE);      \
  }

double 
oligodg(s)
    const char *s;       /* The sequence. */
{
    register int dg = 0;
    register char c;

    /* Use a finite-state machine (DFA) to calucluate dg s. */
    c = *s; s++;
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

 DONE:  /* dg is now computed for the given sequence. */
    return dg / 1000.0;

 ERROR:  /* 
	  * length of s was less than 2 or there was an illegal character in
	  * s.
	  */
    return OLIGOTM_ERROR;
}

double end_oligodg(s, len)
  const char *s;
  int len; /* The number of characters to return. */
{
  int x = strlen(s);
  return x < len ? oligodg(s) : oligodg(s + (x - len));
}

double seqtm(seq, dna_conc, salt_conc, nn_max_len)
  const  char *seq;
  double dna_conc;
  double salt_conc;
  int    nn_max_len;
{
  int len = strlen(seq);
  return (len > nn_max_len)
    ? long_seq_tm(seq, 0, len, salt_conc) : oligotm(seq, dna_conc, salt_conc);
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
  const char *p, *end = &s[len];
  if (len <= 0) return OLIGOTM_ERROR;
  /* Length <= 0 is nonsensical. */
  for (p = &s[0]; p < end; p++) {
    if ('G' == *p || 'g' == *p || 'C' == *p || 'c' == *p)
      GC_count++;
  }

  return
    81.5 
    + (16.6 * log10(salt_conc / 1000.0))
    + (41.0 * (((double) GC_count) / len))
    - (600.0 / len);

}
