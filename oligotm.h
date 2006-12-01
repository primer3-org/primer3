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
double end_oligodg(const char *oligo, int len);

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

/* Return the melting temperature of the given oligo, as calculated using eqn
    (ii) in Rychlik, Spencer, Roads, Nucleic Acids Research, vol 18, no 21, page
    6410, with tables of nearest-neighbor thermodynamics for DNA bases as
    provided in Breslauer, Frank, Bloecker, and Markey,
    Proc. Natl. Acad. Sci. USA, vol 83, page 3748. */
double oligotm(const  char *seq, /* The sequence. */
               double dna_conc,  /* DNA concentration (nanomolar). */
               double salt_conc  /* Salt concentration (millimolar). */
	       );

/* Return the delta G of disruption of oligo using the nearest neighbor model;
   seq should be relatively short, given the characteristics of the nearest
   neighbor model. */
double oligodg(const char *oligo);

/* Return the melting temperature of a given sequence, 'seq'. */
double seqtm(const  char *seq,  /* The sequence. */
             double dna_conc,   /* DNA concentration (nanomolar). */
             double salt_conc,  /* Salt concentration (millimolar). */
             int    nn_max_len  /* The maximum sequence length for using the
				   nearest neighbor model (as implemented
				   in oligotm. For sequences longer
				   than this use the "GC%" formula implemented
				   in long_seq_tm. */
	     );
#endif

