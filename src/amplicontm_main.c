/*
    Amplicon3 calculates melting temperatures for amplicons.
    Copyright (c) 2021 Andreas Untergasser

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    ----------------------------------------------------------------

    Amplicon3 is based on MeltPolymer.c from the DECIPHER package
    created by Erik Wright. DECIPHER is GPL-3 licensed.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "amplicontm.h"

/* Reads the whole string s as a finite number into *out. Returns 1 on
   success, 0 for an empty string, trailing characters, nan or inf. */
static int
parse_finite(const char *s, double *out)
{
  char *end;
  *out = strtod(s, &end);
  if (end == s || *end != '\0' || !isfinite(*out)) {
    return 0;
  }
  return 1;
}

/* Print the melting temperature of an amplicon on stdout. */
/* This program provides a command line interface to
   the function amplicontm() in amplicontm.c
*/
int
main(int argc, char **argv)
{
  amplicon_result ret;
  int no_melt_point = 0;

  const char *msg = "USAGE: %s OPTIONS amplicon\n"
    "\n"
    "where amplicon is a DNA sequence of at least 20 bases.\n"
    "Characters other than A, C, G and T (also N and IUPAC codes)\n"
    "are removed from the sequence before the calculation.\n"
    "\n"
    "and\n"
    "\n"
    "OPTIONS can include any of the the following:\n"
    "\n"
    "-mv monovalent_conc     - concentration of monovalent cations in mM, by default 50mM\n"
    "\n"
    "-dv divalent_conc       - concentration of divalent cations in mM, by default 1.5mM\n"
    "\n"
    "-n  dNTP_conc           - concentration of deoxynycleotide triphosphate in mM, by default 0.6mM\n"
    "\n"
    "-dmso dmso_conc         - concentration of DMSO in percent, by default 0%%\n"
    "\n"
    "-dmso_fact dmso_fact    - DMSO correction factor, by default 0.6\n"
    "\n"
    "-formamid formamid_conc - concentration of formamide in mol/l, by default 0 mol/l\n"
    "\n"
    "-fs find_salt_conc  - provide a measured Tm (>= 0) and obtain a suitable salt concentration (experimental)\n"
    "\n"

    "-tp [0|1]     - Specifies the table of thermodynamic parameters and\n"
    "                the method of melting temperature calculation:\n"
    "                 0  Breslauer et al., 1986 and Rychlik et al., 1990\n"
    "                    (not supported by amplicon3)\n"
    "                 1  Use nearest neighbor parameters from SantaLucia 1998\n"
    "                    *This is the default and recommended value*\n"
    "\n"
    "-sc [0..2]    - Specifies salt correction formula for the melting \n"
    "                 temperature calculation\n"
    "                  0  Schildkraut and Lifson 1965, used by primer3 up to \n"
    "                     and including release 1.1.0.\n"
    "                  1  SantaLucia 1998\n"
    "                     *This is the default and recommended value*\n"
    "                  2  Owczarzy et al., 2004 (not supported by amplicon3)\n"
    "\n"
    "-mf [0..1]    - Specifies the formula for the melting temperature\n"
    "                 calculation\n"
    "                  0  Bolton and McCarthy 1962, used by primer3 up to \n"
    "                     and including release 2.4.2.\n"
    "                  1  Wright and Vetsigian 2016\n"
    "                     *This is the default and recommended value*\n"
    "\n"
    "-o [0..2]     - Specifies the output format on stdout\n"
    "                  0  A table of helicities.\n"
    "                  1  A short list of the calculated parameters \n"
    "                     *This is the default*\n"
    "                  2  A long list of the calculated parameters including \n"
    "                     the melting curve and the derivative curve \n"
    "\n";

   const char *copyright = 
"    Amplicon3 calculates melting temperatures for amplicons.\n"
"    Copyright (c) 2021 Andreas Untergasser\n"
"\n"
"    This program is free software: you can redistribute it and/or modify\n"
"    it under the terms of the GNU General Public License as published by\n"
"    the Free Software Foundation, either version 3 of the License, or\n"
"    (at your option) any later version.\n"
"\n"
"    This program is distributed in the hope that it will be useful,\n"
"    but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"    GNU General Public License for more details.\n"
"\n"
"    You should have received a copy of the GNU General Public License\n"
"    along with this program.  If not, see <https://www.gnu.org/licenses/>.\n"
"\n"
"    ----------------------------------------------------------------\n"
"\n"
"    Amplicon3 is based on MeltPolymer.c from the DECIPHER package\n"
"    created by Erik Wright. DECIPHER is GPL-3 licensed.\n"
"\n";

   char *endptr, *seq;
   double mv = 50.0, dv = 1.5, dntp = 0.6, fs_temp = -10.0;
   double dmso = 0.0, dmso_fact = 0.6, formamid = 0.0;

   int tm_parameters=1, salt_corrections=1, tm_formula=1, output=1;
   int i, j, k;
   if (argc < 2 || argc > 24) {
     fprintf(stderr, msg, argv[0]);       
     fprintf(stderr, "%s", copyright);
     return -1;
   }

   for (i=1; i < argc; ++i) {
     if (strcmp(argv[i], "-mv") == 0) { /* conc of monovalent cations */
       if (i+1 >= argc) {
         /* Missing value */
         fprintf(stderr, msg, argv[0]);
         exit(-1);
       }
       if (!parse_finite(argv[i+1], &mv) || mv < 0.0) {
         fprintf(stderr, msg, argv[0]);
         exit(-1);
       }
       i++;
     } else if (strcmp(argv[i], "-dv") == 0) { /* conc of divalent cations */
       if (i+1 >= argc) {
         /* Missing value */
         fprintf(stderr, msg, argv[0]);
         exit(-1);
       }
       if (!parse_finite(argv[i+1], &dv) || dv < 0.0) {
         fprintf(stderr, msg, argv[0]);
         exit(-1);
       }
       i++;
     } else if (strcmp(argv[i], "-n") == 0) { /* conc of dNTP */
       if (i+1 >= argc) {
         /* Missing value */
         fprintf(stderr, msg, argv[0]);
         exit(-1);
       }
       if (!parse_finite(argv[i+1], &dntp) || dntp < 0.0) {
         fprintf(stderr, msg, argv[0]);
         exit(-1);
       }
       i++;
     } else if (strcmp(argv[i], "-dmso_fact") == 0) { /* DMSO correction factor */
       if (i+1 >= argc) {
         /* Missing value */
         fprintf(stderr, msg, argv[0]);
         exit(-1);
       }
       if (!parse_finite(argv[i+1], &dmso_fact) || dmso_fact < 0.0) {
         fprintf(stderr, msg, argv[0]);
         exit(-1);
       }
       i++;
     } else if (strcmp(argv[i], "-dmso") == 0) { /* concentration of DMSO */
       if (i+1 >= argc) {
         /* Missing value */
         fprintf(stderr, msg, argv[0]);
         exit(-1);
       }
       if (!parse_finite(argv[i+1], &dmso) || dmso < 0.0 || dmso > 100.0) {
         fprintf(stderr, msg, argv[0]);
         exit(-1);
       }
       i++;
     } else if (strcmp(argv[i], "-formamid") == 0) { /* concentration of formamid */
       if (i+1 >= argc) {
         /* Missing value */
         fprintf(stderr, msg, argv[0]);
         exit(-1);
       }
       /* Pure formamide is about 25.1 mol/l */
       if (!parse_finite(argv[i+1], &formamid) || formamid < 0.0 || formamid > 25.1) {
         fprintf(stderr, msg, argv[0]);
         exit(-1);
       }
       i++;
     } else if (strcmp(argv[i], "-tp") == 0) { /* parameters for melting temperature calculation */
       if (i+1 >= argc) {
         /* Missing value */
         fprintf(stderr, msg, argv[0]);
         exit(-1);
       }
       tm_parameters = (int)strtol(argv[i+1], &endptr, 10);
       if ('\0' != *endptr || tm_parameters<0 || tm_parameters>1) {
         fprintf(stderr, msg, argv[0]);
         exit(-1);
       }
       i++;
     } else if (strcmp(argv[i], "-sc") == 0) { /* method of salt correction */
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
     } else if (strcmp(argv[i], "-mf") == 0) { /* formula for melting temperature calculation */
       if (i+1 >= argc) {
         /* Missing value */
         fprintf(stderr, msg, argv[0]);
         exit(-1);
       }
       tm_formula = (int)strtol(argv[i+1], &endptr, 10);
       if ('\0' != *endptr || tm_formula<0 || tm_formula>1) {
         fprintf(stderr, msg, argv[0]);
         exit(-1);
       }
       i++;
     } else if (strcmp(argv[i], "-fs") == 0) { /* find salt concentation */
       if (i+1 >= argc) {
         /* Missing value */
         fprintf(stderr, msg, argv[0]);
         exit(-1);
       }
       if (!parse_finite(argv[i+1], &fs_temp) || fs_temp < 0.0) {
         fprintf(stderr, msg, argv[0]);
         exit(-1);
       }
       i++;
     } else if (strcmp(argv[i], "-o") == 0) { /* output helicities table or list */
       if (i+1 >= argc) {
         /* Missing value */
         fprintf(stderr, msg, argv[0]);
         exit(-1);
       }
       output = (int)strtol(argv[i+1], &endptr, 10);
       if ('\0' != *endptr || output<0 || output>2) {
         fprintf(stderr, msg, argv[0]);
         exit(-1);
       }
       i++;
     } else if (!strncmp("-", argv[i], 1)) {
       /* Unknown option. */
       fprintf(stderr, msg, argv[0]);
       exit(-1);
     } else
       break;                /* all args processed. go on to sequences. */
   }
   
  if(!argv[i]) { /* if no oligonucleotide sequence is specified */
    fprintf(stderr, msg, argv[0]);
    exit(-1);
  }
  /* input sequence to uppercase */
  seq = argv[i];

  if (fs_temp < 0.0) {
    ret = amplicontm(seq,
                     mv,
                     dv,
                     dntp,
                     dmso,
                     dmso_fact,
                     formamid,
                     (amp_tm_parameters_type) tm_parameters,
                     (amp_salt_correction_type) salt_corrections,
                     (amp_tm_method_type) tm_formula,
                     output);
  } else {
    ret = ampliconfindsalt(seq,
                           fs_temp,
                           dmso,
                           dmso_fact,
                           formamid,
                           (amp_tm_parameters_type) tm_parameters,
                           (amp_salt_correction_type) salt_corrections,
                           (amp_tm_method_type) tm_formula,
                           output);
  }

  /* With -o 2 the curves are printed even without a melting point,
     they show why there is none; the program still fails at the end */
  if (ret.error == 3 && output == 2 &&
      ret.temp != NULL && ret.melt != NULL && ret.deriv != NULL) {
      no_melt_point = 1;
      ret.error = 0;
  }

  if (ret.error == 1) {
      fprintf(stderr, "%s ERROR: out of memory\n", argv[0]);
      return -1;
  } else if (ret.error == 2) {
      fprintf(stderr, "%s ERROR: input is corrupted. Sodium equivalent concentration must be at least 0.01M.\n", argv[0]);
      return -1;
  } else if (ret.error == 3) {
      fprintf(stderr, "%s ERROR: no melting temperature found between 60 and 95 C.\n", argv[0]);
      return -1;
  } else if (ret.error == 4) {
      fprintf(stderr, "%s ERROR: -sc 2 salt_corrections==owczarzy is not supported.\n", argv[0]);
      return -1;
  } else if (ret.error == 5) {
      fprintf(stderr, "%s ERROR: -tp 0 breslauer parameters for melting temperature calculation are not supported.\n", argv[0]);
      return -1;
  } else if (ret.error == 6) {
      fprintf(stderr, "%s ERROR: the sequence contains no A, C, G or T.\n", argv[0]);
      return -1;
  } else if (ret.error == 7) {
      fprintf(stderr, "%s ERROR: the sequence has fewer than %d A, C, G or T.\n", argv[0], AMP_MIN_LENGTH);
      return -1;
  } else if (ret.error == 8) {
      fprintf(stderr, "%s ERROR: unknown parameters for melting temperature calculation or salt correction.\n", argv[0]);
      return -1;
  } else if (ret.error == 9) {
      fprintf(stderr, "%s ERROR: concentrations of salt, divalent cations and dNTPs must not be negative.\n", argv[0]);
      return -1;
  } else if (ret.error == 11) {
      fprintf(stderr, "%s ERROR: the melting temperature is below %.1f C even at the highest sodium concentration of about 1 M.\n", argv[0], fs_temp);
      return -1;
  } else if (ret.error == 10) {
      fprintf(stderr, "%s ERROR: the melting temperature is above %.1f C even at the lowest sodium concentration of 10 mM.\n", argv[0], fs_temp);
      return -1;
  } else if ((output != 1) && (tm_formula == 0)) {
      fprintf(stderr, "%s ERROR: -mf 0 requires -o 1.\n", argv[0]);
      return -1;
  } else {
      if (output == 0) {
           printf("\t");
           for (j = 0; j < ret.seq_len; j++) {
               if (j == ret.seq_len - 1) {
                   printf("%d\n", j + 1);
               } else {
                   printf("%d\t", j + 1);
               }
           }
           printf("\t");
           for (j = 0; j < ret.seq_len; j++) {
               if (j == ret.seq_len - 1) {
                   printf("%c\n", ret.seq[j]);
               } else {
                   printf("%c\t", ret.seq[j]);
               }
           }
           for (j = 0; j < ret.temp_len; j++) {
               printf("%.1f\t", ret.temp[j]);
               for (k = 0; k < ret.seq_len; k++) {
                   if (k == ret.seq_len - 1) {
                       printf("%.6f\n", *(ret.pos_prob + j + k * ret.temp_len));
                   } else {
                       printf("%.6f\t", *(ret.pos_prob + j + k * ret.temp_len));
                   }
               }
           }
      } else {
           printf("AMPLICON_MONOVALENT=%.1f\n", ret.mv);
           printf("AMPLICON_DIVALENT=%.1f\n", ret.dv);
           printf("AMPLICON_DNTPS=%.1f\n", ret.dntp);
           printf("AMPLICON_DMSO=%.1f\n", ret.dmso);
           printf("AMPLICON_DMSO_CORRECTION=%.1f\n", ret.dmso_fact);
           printf("AMPLICON_FORMAMID=%.1f\n", ret.formamid);
           printf("AMPLICON_PRODUCT_SIZE=%d\n", ret.seq_len);
           if (output == 2) {
               printf("AMPLICON_SEQUENCE=%s\n", ret.seq);
           }
           printf("AMPLICON_GC_PERCENT=%.1f\n", ret.seq_gc);
           if (ret.melt_len > 0) {
               printf("AMPLICON_MELTPOINTS=");
               for (j = 0; j < ret.melt_len; j++) {
                   if (j == ret.melt_len - 1) {
                       printf("%.1f\n", ret.melt_points[j]);
                   } else {
                       printf("%.1f,", ret.melt_points[j]);
                   }
               }
           }
           if (output == 2) {
               printf("AMPLICON_TEMPERATURES=");
               for (j = 0; j < ret.temp_len; j++) {
                   if (j == ret.temp_len - 1) {
                       printf("%.1f\n", ret.temp[j]);
                   } else {
                       printf("%.1f,", ret.temp[j]);
                   }
               }
               printf("AMPLICON_MELT_CURVE=");
               for (j = 0; j < ret.temp_len; j++) {
                   if (j == ret.temp_len - 1) {
                       printf("%.6e\n", ret.melt[j]);
                   } else {
                       printf("%.6e,", ret.melt[j]);
                   }
               }
               printf("AMPLICON_DERIVATIVE_CURVE=");
               for (j = 0; j < ret.temp_len; j++) {
                   if (j == ret.temp_len - 1) {
                       printf("%.6e\n", ret.deriv[j]);
                   } else {
                       printf("%.6e,", ret.deriv[j]);
                   }
               }
           }
      }
  }

  free_amplicon_result(&ret);
  if (no_melt_point == 1) {
      fprintf(stderr, "%s ERROR: no melting temperature found between 60 and 95 C.\n", argv[0]);
      return -1;
  }
  return 0;
}
