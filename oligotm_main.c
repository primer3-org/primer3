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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "oligotm.h"

/* Print the melting temperature of an oligo on stdout. */
int
main(argc, argv)
    int argc;
    const char **argv;
{
  double tm;

  char *msg = "USAGE: %s OPTIONS oligo\n"
    "OPTIONS:\n"
    "-k  salt_conc     - concentration of monovalent cations in mM, by default 50mM\n"
    "-d  dna_conc      - concentration of DNA strands in nM, by default 50nM\n"
    "-ts tm_santalucia - option for the table of thermodynamic parameters and\n"
    "                     for the method of melting temperature calculation:\n"
    "                     1 - Breslauer et al., 1986 and Rychlik et al., 1990\n"
    "                         (used by primer3 until version 1.1.0)\n"
    "                     2 - SantaLucia 1998, by default\n"
    "-sc salt_correct - option for salt correction formula for the melting temperature calculation\n"
    "                     1 - Schildkraut and Lifson 1965 (used by primer3 until version 1.1.0)\n"
    "                     2 - SantaLucia 1998, by default\n"
    "                     3 - Owczarzy et al., 2004\n\n"
    "-i               - prints references to publications which were used for thermodynamic calculations\n"
    "oligo            - oligonucleotide sequence of between 2 and 36 bases\n"
    "                   (bases must be uppercase)\n"
    "Prints oligo's melting temperature on stdout.\n";
   
  char *info = "1. Breslauer KJ, Frank R, Blöcker H and Marky LA. (1986) Predicting DNA duplex stability from the base sequence. Proc. Natl. Acad. Sci., 83, 4746-50.\n\n"
    "2. Rychlik W, Spencer WJ and Rhoads RE. (1990) Optimization of the annealing temperature for DNA amplification in vitro. Nucleic Acids Res., 118, 6409-12.\n\n"
    "3. SantaLucia JR. (1998). A unified view of polymer, dumbbell and oligonucleotide DNA nearest-neighbor thermodynamics. Proc. Natl. Acad. Sci., 95, 1460-65.\n\n"
    "4. Schildkraut, C, and Lifson, S. (1965) Dependence of the melting temperature of DNA on salt concentration. Biopolymers, 3, 195-208.\n\n"
    "5. Owczarzy R, You Y, Moreira BG, Manthey JA, Huang L, Behlke MA and Walder JA. (2004) Effects of Sodium Ions on DNA Duplex Oligomers: Improved Predictions of Melting Temperatures. Biochemistry, 43, 3537-54.\n";

  char *endptr;
  long k = 50, d = 50;
  int tm_santalucia=2, salt_corrections=2;
  int i;
   
  if (argc < 2 || argc > 10) {
    fprintf(stderr, msg, argv[0]);       
    return -1;
  }

  for (i=1; i < argc; ++i) {
    if (!strncmp("-k", argv[i], 2)) {
      k = strtol(argv[i+1], &endptr, 10);
      if ('\0' != *endptr) {
	fprintf(stderr, msg, argv[0]);
	exit(-1);
      }
      i++;
    } else if (!strncmp("-d", argv[i], 2)) {
      d = strtol(argv[i+1], &endptr, 10);
      if ('\0' != *endptr) {
	fprintf(stderr, msg, argv[0]);
	exit(-1);
      }
      i++;
    } else if (!strncmp("-ts", argv[i], 3)) { /* added by T.Koressaar */
      tm_santalucia = atoi(argv[i+1]);
      if (!tm_santalucia || '0' == tm_santalucia) {	      
	fprintf(stderr, msg, argv[0]);
	exit(-1);
      }
      i++;
    } else if (!strncmp("-sc", argv[i], 3)) { /* added by T.Koressaar */
      salt_corrections = atoi(argv[i+1]);
      if (!salt_corrections || '0' == salt_corrections) {
	fprintf(stderr, msg, argv[0]);
	exit(-1);
      }
      i++;
    } else if (!strncmp("-i", argv[i], 2)) {
      fprintf(stderr, info, argv[0]);
      exit(-1);
    } else if (!strncmp("-", argv[i], 1)) {
      /* Unknown option. */
      fprintf(stderr, msg, argv[0]);
      exit(-1);
    } else
      break;		/* all args processed. go on to sequences. */
  }
   
  /* added by T.Koressaar, start */
  if (tm_santalucia==1) {
    tm_santalucia=0;
  } else if (tm_santalucia==2) {
    tm_santalucia=1;
  } else {
    tm_santalucia=1;
  }
   
  if (salt_corrections==1) {
    salt_corrections=0;
  } else if (salt_corrections==2) {
    salt_corrections=1;
  } else if (salt_corrections==3) {
    salt_corrections=2;
  } else {
    salt_corrections=1;
  }
  /* added by T.Koressaar, end */
  if(!argv[i]) { /* if no oligonucleotide sequence is specified */
    fprintf(stderr, msg, argv[0]);
    exit(-1);
  }
   
  tm = oligotm(argv[i], d, k, tm_santalucia, salt_corrections);
  if (OLIGOTM_ERROR == tm) {
    fprintf(stderr,
	    "%s: length of %s is less than 2 or it contains an illegal character\n",
	    argv[0], argv[i]);
    return -1;
  }
  fprintf(stdout, "%f\n", tm);
  return 0;
}
