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
     
#ifndef _OLIGO_TM
#define _OLIGO_TM 1
#define OLIGOTM_ERROR -999999.9999

/* Return the delta G of the last len bases of oligo if oligo is at least len
   bases long; otherwise return the delta G of oligo. */
double end_oligodg(const char *oligo, int len, int tm_santalucia);

/* Calculate the melting temperature of substr(seq, start, length) using the
   formula from Bolton and McCarthy, PNAS 84:1390 (1962) as presented in
   Sambrook, Fritsch and Maniatis, Molecular Cloning, p 11.46 (1989, CSHL
   Press).

   Tm = 81.5 + 16.6(log10([Na+])) + .41*(%GC) - 600/length

   Where [Na+] is the molar sodium concentration, (%GC) is the percent of Gs
   and Cs in the sequence, and length is the length of the sequence.

   A similar formula is used by the prime primer selection program in GCG
   (http://www.gcg.com), which instead uses 675.0 / length in the last term
   (after F. Baldino, Jr, M.-F. Chesselet, and M.E.  Lewis, Methods in
   Enzymology 168:766 (1989) eqn (1) on page 766 without the mismatch and
   formamide terms).  The formulas here and in Baldino et al. assume Na+ rather
   than K+.  According to J.G. Wetmur, Critical Reviews in BioChem. and
   Mol. Bio. 26:227 (1991) 50 mM K+ should be equivalent in these formulae to .2
   M Na+.

   This function takes salt_conc to be the millimolar (mM) concentration,
   since mM is the usual units in PCR applications.

 */
double long_seq_tm(const char *seq, int start, int length, double salt_conc);

/* Return the melting temperature of the given oligo calculated as specified by user. 
 
 if tm_santalucia==1 then the table of nearest-neighbor thermodynamic parameters 
 and method for Tm calculation suggested by "SantaLucia J Jr (1998) A unified view 
 of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. 
 Proc Natl Acad Sci USA. 95(4):1460-5" is used. By default. Added by T. Koressaar
 
 if tm_santalucia==0 then method for Tm calculations suggested by "Rychlik, Spencer, 
 Roads, Nucleic Acids Research, vol 18, no 21, page 6410" and the table of 
 nearest-neighbor thermodynamic parameters suggested by "Breslauer, Frank, Bloecker, 
 and Markey,Proc. Natl. Acad. Sci. USA, vol 83, page 3748" is used. 
 This is the method and the table Primer3 used until version 1.0.1 (including the version 1.0.1)
 
 if salt_corrections==0 then formula for salt correction suggested by Schildkraut and Lifson is used.
 "Schildkraut, C, and Lifson, S. (1965) Dependence of the melting temperature of DNA on salt concentration. 
 Biopolymers, 3, 195-208." This is the formula Primer3 used until version 1.0.1 (including the version 1.0.1)

 if salt_corrections==1 then formula for salt correction suggested by SantaLucia is used.
 "SantaLucia JR. (1998). A unified view of polymer, dumbbell and oligonucleotide DNA 
 nearest-neighbor thermodynamics. Proc. Natl. Acad. Sci., 95, 1460-65." Used by default. Added by T.Koressaar
 
 if salt_corrections==2 then formula for salt correction suggested by Owczarzy et al. is used
 "Owczarzy R, You Y, Moreira BG, Manthey JA, Huang L, Behlke MA and Walder JA. (2004) 
 Effects of Sodium Ions on DNA Duplex Oligomers: Improved Predictions of Melting Temperatures. 
 Biochemistry, 43, 3537-54." Added by T.Koressaar
 
 */

double oligotm(const  char *seq, /* The sequence. */
               double dna_conc,  /* DNA concentration (nanomolar). */
               double salt_conc,  /* Salt concentration (millimolar). */
	       int tm_santalucia, /* the same as in function double seqtm(). See forward */
	       int salt_corrections /* the same as in function double seqtm(). See forward */
	       );

/* Return the delta G of disruption of oligo using the nearest neighbor model;
   seq should be relatively short, given the characteristics of the nearest
   neighbor model. */
double oligodg(const char *oligo, int tm_santalucia);

/* Return the melting temperature of a given sequence, 'seq'. */
double seqtm(const  char *seq,  /* The sequence. */
             double dna_conc,   /* DNA concentration (nanomolar). */
             double salt_conc,  /* Salt concentration (millimolar). */
             int    nn_max_len,  /* The maximum sequence length for using the
				   nearest neighbor model (as implemented
				   in oligotm. For sequences longer
				   than this use the "GC%" formula implemented
				   in long_seq_tm. */
	     int tm_santalucia, /* if =1 then the table of thermodynamic parameters of SantaLucia
				 1998 (also the method for calculating Tm) is used for oligo 
				 melting temperature calculation; else the table of thermodynamic 
				 parameters of Breslauer et al., 1986 is used for Tm calculation 
				 (and the method suggested by Rychlik et al., 1990) */
	     int salt_corrections /* The salt correction formula for correcting the melting temperature 
				    calculated with nearest-neighbor model
				    0 - salt correction formula suggested by Schildkraut and Lifson 
				    (Schildkraut and Lifson, 1965)
				    1 - salt correction formula suggested by SantaLucia (SantaLucia, 1998)
				    2 - salt correction formula suggested by Owczarzy et al., (Owczarzy et al., 2004)
				   */
	     );


/* Returns 1 if the sequence is self-complementary or symmetrical; 0 otherwise*/
int Symmetry(const char *seq);

/*
 * Tables of nearest-neighbor thermodynamics for DNA bases.
 * See Breslauer, Frank, Blocker, and Markey,
 * "Predicting DNA duplex stability from the base sequence."
 * Proc. Natl. Acad. Sci. USA, vol 83, page 3746 (1986).
 * Article free at
 * http://www.pubmedcentral.nih.gov/picrender.fcgi?artid=323600&blobtype=pdf
 * See table 2.
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

/*
 * Tables of nearest-neighbor thermodynamics for DNA bases.
 * See SantaLucia, Proc. Natl. Sci. USA, vol 95, pages 1460-1465
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

#endif

