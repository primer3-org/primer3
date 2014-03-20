/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008
Whitehead Institute for Biomedical Research, Steve Rozen
(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky.
All rights reserved.

    This file is part of the primer3 suite and libraries.

    The primer3 suite and libraries are free software;
    you can redistribute them and/or modify them under the terms
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
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "oligotm.h"

/* Print the melting temperature of an oligo on stdout. */
/* This program provides a command line interface to
   the function oligotm() in oligtm.c
*/
int
main(int argc, char **argv)
{
  double tm;

  const char *msg = "USAGE: %s OPTIONS oligo\n"
    "\n"
    "where oligo is a DNA sequence of between 2 and 36 bases\n"
    "\n"
    "and\n"
     "\n"
     "OPTIONS can include any of the the following:\n"
     "\n"
     "-mv monovalent_conc - concentration of monovalent cations in mM, by default 50mM\n"
     "\n"
     "-dv divalent_conc   - concentration of divalent cations in mM, by default 0mM\n"
     "\n"
     "-n  dNTP_conc       - concentration of deoxynycleotide triphosphate in mM, by default 0mM\n"
     "\n"
     "-d  dna_conc        - concentration of DNA strands in nM, by default 50nM\n"
     "\n"

    "-tp [0|1]     - Specifies the table of thermodynamic parameters and\n"
    "                the method of melting temperature calculation:\n"
    "                 0  Breslauer et al., 1986 and Rychlik et al., 1990\n"
    "                    (used by primer3 up to and including release 1.1.0).\n"
    "                    This is the default, but _not_ the recommended value.\n"
    "                 1  Use nearest neighbor parameters from SantaLucia 1998\n"
    "                    *THIS IS THE RECOMMENDED VALUE*\n"
    "\n"
    "-sc [0..2]    - Specifies salt correction formula for the melting \n"
    "                 temperature calculation\n"
    "                  0  Schildkraut and Lifson 1965, used by primer3 up to \n"
    "                     and including release 1.1.0.\n"
    "                     This is the default but _not_ the recommended value.\n"
    "                  1  SantaLucia 1998\n"
    "                     *THIS IS THE RECOMMENDED VAULE*\n"
    "                  2  Owczarzy et al., 2004\n\n"
    "\n\n"
    "Prints oligo's melting temperature on stdout.\n";
   
   const char *copyright = 
"Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006\n"
"Whitehead Institute for Biomedical Research, Steve Rozen\n"
"(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky\n"
"All rights reserved.\n"
"\n"
"    This file is part of the oligotm library.\n"
"\n"
"    The oligotm library is free software; you can redistribute it and/or modify\n"
"    it under the terms of the GNU General Public License as published by\n"
"    the Free Software Foundation; either version 2 of the License, or\n"
"    (at your option) any later version.\n"
"\n"
"    The oligotm library is distributed in the hope that it will be useful,\n"
"    but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"    GNU General Public License for more details.\n"
"\n"
"    You should have received a copy of the GNU General Public License\n"
"    along with the oligtm library (file gpl-2.0.txt in the source\n"
"    distribution).  If not, see http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt;\n"
"    or write to the Free Software Foundation, Inc.,\n"
"    51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA\n";

   char *endptr, *seq;
   double mv = 50, d = 50;
   double dv = 0, n = 0;
   int tm_santalucia=0, salt_corrections=0;
   int i, j, len;
   if (argc < 2 || argc > 14) {
     fprintf(stderr, msg, argv[0]);       
     fprintf(stderr, "%s", copyright);
     return -1;
   }

   for (i=1; i < argc; ++i) {
     if (!strncmp("-mv", argv[i], 3)) { /* conc of monovalent cations */
       if (i+1 >= argc) {
	 /* Missing value */
	 fprintf(stderr, msg, argv[0]);
	 exit(-1);
       }
       mv = strtod(argv[i+1], &endptr);
       if ('\0' != *endptr) {
	 fprintf(stderr, msg, argv[0]);
	 exit(-1);
       }
       i++;
     } else if (!strncmp("-dv", argv[i], 3)) { /* conc of divalent cations; added by T.Koressaar */
       if (i+1 >= argc) {
	 /* Missing value */
	 fprintf(stderr, msg, argv[0]);
	 exit(-1);
       }
       dv = strtod(argv[i+1], &endptr);
       if('\0' != *endptr) {
	 fprintf(stderr, msg, argv[0]);
	 exit(-1);
       }
       i++;
     } else if (!strncmp("-n", argv[i], 2)) { /* conc of dNTP; added by T.Koressaar */
       if (i+1 >= argc) {
	 /* Missing value */
	 fprintf(stderr, msg, argv[0]);
	 exit(-1);
       }
       n = strtod(argv[i+1], &endptr);
       if('\0' != *endptr) {
	 fprintf(stderr, msg, argv[0]);
	 exit(-1);
       }
       i++;
     } else if (!strncmp("-d", argv[i], 2)) {
       if (i+1 >= argc) {
	 /* Missing value */
	 fprintf(stderr, msg, argv[0]);
	 exit(-1);
       }
       d = strtod(argv[i+1], &endptr);
       if ('\0' != *endptr) {
	 fprintf(stderr, msg, argv[0]);
	 exit(-1);
       }
       i++;
     } else if (!strncmp("-tp", argv[i], 3)) { /* added by T.Koressaar */
       if (i+1 >= argc) {
	 /* Missing value */
	 fprintf(stderr, msg, argv[0]);
	 exit(-1);
       }
       tm_santalucia = (int)strtol(argv[i+1], &endptr, 10);
       if ('\0' != *endptr || tm_santalucia<0 || tm_santalucia>1) {	  
	 fprintf(stderr, msg, argv[0]);
	 exit(-1);
       }
       i++;
     } else if (!strncmp("-sc", argv[i], 3)) { /* added by T.Koressaar */
       if (i+1 >= argc) {
	 /* Missing value */
	 fprintf(stderr, msg, argv[0]);
	 exit(-1);
       }
       salt_corrections = (int)strtol(argv[i+1], &endptr, 10);
       if ('\0' != *endptr || salt_corrections<0 || salt_corrections>2) {
	 fprintf(stderr, msg, argv[0]);
	 exit(-1);
       }
       i++;
     } else if (!strncmp("-", argv[i], 1)) {
       /* Unknown option. */
       fprintf(stderr, msg, argv[0]);
       exit(-1);
     } else
       break;		/* all args processed. go on to sequences. */
   }
   
  if(!argv[i]) { /* if no oligonucleotide sequence is specified */
    fprintf(stderr, msg, argv[0]);
    exit(-1);
  }
  /* input sequence to uppercase */
  seq = argv[i];
  len=strlen(seq);
  for(j=0;j<len;j++) seq[j]=toupper(seq[j]);
   
  tm = oligotm(seq, d, mv, dv, n, (tm_method_type) tm_santalucia, (salt_correction_type) salt_corrections);
  if (OLIGOTM_ERROR == tm) {
    fprintf(stderr,
	    "%s ERROR: length of sequence %s is less than 2 or\n"
	    "             the sequence contains an illegal character or\n" 
	    "             you have specified incorrect value for concentration of divalent cations or\n"
	    "             you have specified incorrect value for concentration of dNTPs\n",
	    argv[0], argv[i]);
    return -1;
  }
  fprintf(stdout, "%f\n", tm);
  return 0;
}
