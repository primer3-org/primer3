/*
 Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008,2009
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
#include <sys/types.h>
#include <sys/stat.h>

#include "thal.h"
#include "thal_parameters.h"

/* Check on which OS we compile */
#if defined(_WIN32) || defined(WIN32) || defined (__WIN32__) || defined(__CYGWIN__) || defined(__MINGW32__)
#define OS_WIN
#endif

#define DEBUG

char *endptr; /* reading input */
int i; /* index */
const unsigned char *oligo1, *oligo2; /* inserted oligo sequences */
char *path = NULL; /* path to the parameter files */

/* Beginning of main */
int main(int argc, char** argv) 
{   
   const char* usage;
   int tmp_ret = 0;
   thal_args a;
   thal_results o;
   set_thal_default_args(&a);
   thal_mode mode = THL_GENERAL; /* by default print only melting temperature, 
                                    do not draw structure or print any additional parameters */
   int thal_debug = 0;
   int thal_only = 0;
   int batch_mode = 0;
   if((mode == THL_DEBUG) || (mode == THL_DEBUG_F)) {
#undef DEBUG
   } else {
#define DEBUG
   }

   usage = "USAGE: %s OPTIONS oligo\n"
     "-mv monovalent_conc  - concentration of monovalent cations in mM, by default 50 mM\n"
     "\n"
     "-dv divalent_conc    - concentration of divalent cations in mM, by default 0 mM\n"
     "\n"
     "-n  dNTP_conc        - concentration of deoxynycleotide triphosphate in mM, by default 0 mM\n"
     "\n"
     "-d  dna_conc         - concentration of DNA strands in nM, by default 50 nM\n"
     "\n"
     "-a  mode             - alignment type, END1, END2, ANY and HAIRPIN, by default ANY (when duplex)\n"
     "\n"
     "-t  temp             - temperature at which duplex is calculated, by default 37C\n"
     "\n"
     "-r                   - causes the alignment NOT to be displayed on stderr, _only_ Tm is printed\n"
     "\n"
     "-maxloop size        - the maximum size of secondary structures loops.\n" 
     "                       Default is 30 (this is maximum allowed length, currently).\n"
     "\n"
     "-path <path>         - the path to the thermodynamic parameter files\n"
     "\n"
     "-b                   - batch mode, read primers from stdin\n"
     "                       Input format is: oligo1 [oligo2]\n"
     "\n"
     "-s1 DNA_oligomer\n"
     "\n"
     "-s2 DNA_oligomer\n"
     "\n";
   if(argc < 2) {
#ifdef DEBUG
      fprintf(stderr, usage, argv[0]);
#endif
      return -1;
   }
   /* BEGIN: READ the INPUT */
   for(i = 1; i < argc; ++i) {
      if (!strncmp("-mv", argv[i], 3)) { /* conc of monovalent cations */
         if(argv[i+1]==NULL) {
#ifdef DEBUG
            fprintf(stderr, usage, argv[0]);
#endif
            exit(-1);
         }
         a.mv = strtod(argv[i+1], &endptr);
         if ('\0' != *endptr || a.mv < 0.0) {
#ifdef DEBUG
            fprintf(stderr, usage, argv[0]);
#endif
            exit(-1);
         }
         i++;
      } else if (!strncmp("-dv", argv[i], 3)) { /* conc of divalent cations */
         if(argv[i+1]==NULL) {
#ifdef DEBUG
            fprintf(stderr, usage, argv[0]);
#endif
            exit(-1);
         }
         a.dv = strtod(argv[i+1], &endptr);
         if('\0' != *endptr || a.dv < 0.0) {
#ifdef DEBUG
            fprintf(stderr, usage, argv[0]);
#endif
            exit(-1);
         }
         i++;
      } else if (!strcmp("-path", argv[i])) {
        if(argv[i+1]==NULL) {
#ifdef DEBUG
            fprintf(stderr, usage, argv[0]);
#endif
            exit(-1);
         }
         path = (char*)argv[i+1];
         i++;
      } else if (!strncmp("-s1", argv[i], 3)) { /* first sequence in 5'->3' direction */
         if(argv[i+1]==NULL) {
#ifdef DEBUG
            fprintf(stderr, usage, argv[0]);
#endif
            exit(-1);
         }
         oligo1 = (const unsigned char*)argv[i+1];
         i++;         
      } else if (!strncmp("-s2", argv[i], 3)) { /* second sequence in 5'->3' direction */
         if(argv[i+1]==NULL) {
#ifdef DEBUG
            fprintf(stderr, usage, argv[0]);
#endif
            exit(-1);
         }
         oligo2 = (const unsigned char*)argv[i+1];
         i++;
      } else if (!strncmp("-a", argv[i], 2)) {          /* annealing type END1, END2, ANY, considered only when duplexis; 
                                                  by default ANY  */
         if(argv[i+1]==NULL) {
#ifdef DEBUG
            fprintf(stderr, usage, argv[0]);
#endif
            exit(-1);
         }
         if(strcmp(argv[i+1],"END1")==0) {
            a.type = thal_end1;
         } else if(strcmp(argv[i+1],"END2")==0) {
            a.type = thal_end2;
         } else if(strcmp(argv[i+1],"HAIRPIN")==0) {
            a.type = thal_hairpin;
            a.dimer = 0;
         } else if (strcmp(argv[i+1], "ANY")==0) {
               a.type = thal_any; /* ANY */  
         } else {
#ifdef DEBUG
            fprintf(stderr, usage, argv[0]);
#endif
            exit(-1);
         }
         i++;
      } else if (!strncmp("-d", argv[i], 2)) { /* dna conc */
         if(argv[i+1]==NULL) {
#ifdef DEBUG
            fprintf(stderr, usage, argv[0]);
#endif
            exit(-1);
         }
         a.dna_conc = strtod(argv[i+1], &endptr);
         if('\0' != *endptr || a.dna_conc < 0 || a.dna_conc == 0) {
#ifdef DEBUG
            fprintf(stderr, usage, argv[0]);
#endif
            exit(-1);
         }
         i++;
      } else if (!strncmp("-r", argv[i], 2)) { /* only temp is calculated */
         thal_only = 1;
      } else if (!strncmp("-b", argv[i], 2)) { /* batch mode */
         batch_mode = 1;
      } else if (!strncmp("-t", argv[i], 2)) { /* temperature at which sec str are calculated */
         if(argv[i+1]==NULL) {
#ifdef DEBUG
            fprintf(stderr, usage, argv[0]);
#endif
            exit(-1);
         }
         a.temp = strtod(argv[i+1], &endptr) + ABSOLUTE_ZERO;
         if('\0' != *endptr) {
#ifdef DEBUG
            fprintf(stderr, usage, argv[0]);
#endif
            exit(-1);
         }
         i++;
      } else if (!strncmp("-n", argv[i], 2)) { /* concentration of dNTPs */
         if(argv[i+1]==NULL) {
#ifdef DEBUG
            fprintf(stderr, usage, argv[0]);
#endif
            exit(-1);
         }
         a.dntp = strtod(argv[i+1], &endptr);
         if('\0' != *endptr || a.dntp < 0.0) {
#ifdef DEBUG
            fprintf(stderr, usage, argv[0]);
#endif
            exit(-1);
         }
         i++;
      } else if (!strncmp("-maxloop", argv[i], 8)) { /* maximum size of loop calculated; 
                                                      this value can not be larger than 30 */
         if(argv[i+1]==NULL) {
#ifdef DEBUG
            fprintf(stderr, usage, argv[0]);
#endif
            exit(-1);
         }
         a.maxLoop = (int) (strtod(argv[i+1], &endptr));
                 
         if(a.maxLoop > MAX_LOOP ) {
            a.maxLoop = MAX_LOOP;
#ifdef DEBUG
            fputs("Warning: the maximum size of secondary structures loop is set to default (30)\n", stderr);
#endif
         }  else if(a.maxLoop < MIN_LOOP) {         
            a.maxLoop = MIN_LOOP;
#ifdef DEBUG
            fputs("Warning: the maximum size of secondary structures loop was set to minimum size of allowed loop length (0)\n", stderr);
#endif
         } 
         if('\0' != *endptr || a.maxLoop < 0) {
#ifdef DEBUG
            fprintf(stderr, usage, argv[0]);
#endif
            exit(-1);
         }
         i++;
      } else if(!strncmp("-", argv[i], 1)) { /* Unknown option. */
#ifdef DEBUG
         fprintf(stderr, usage, argv[0]);
#endif
         exit(-1);
      } else {
         break;
      }
   }
   /* END reading INPUT */
   
   if (!batch_mode) {
      /* check the input correctness */
      if(a.dimer!=0 && (oligo2==NULL || oligo1==NULL)) { /* if user wants to calculate structure 
                                                          of dimer then two sequences must be defined*/
         fprintf(stderr, usage, argv[0]);
         exit(-1);
      }
      if(a.dimer==0 && (oligo2==NULL && oligo1==NULL)) { /* if user wants to calculate structure
                                                          of monomer then only one sequence must be defined */
         fprintf(stderr, usage, argv[0]);
         exit(-1);
      }
   }
   /* read default thermodynamic parameters */
   thal_parameters thermodynamic_parameters;
   thal_set_null_parameters(&thermodynamic_parameters);

   if (path != NULL) {
      thal_load_parameters(path, &thermodynamic_parameters, &o);
   } else {
      set_default_thal_parameters(&thermodynamic_parameters);
   }
   get_thermodynamic_values(&thermodynamic_parameters, &o);

   if (tmp_ret) {
     fprintf(stderr, "%s\n", o.msg);
     exit(-1);
   }

   /* Set the correct mode */
   if (thal_only) {
     if (thal_debug) {
       mode = THL_DEBUG_F;
     } else {
       mode = THL_FAST;
     }
   } else {
     if (thal_debug) {
       mode = THL_DEBUG;
     } else {
       mode = THL_GENERAL;
     }
   }

   do {
      static char line[200], primer1[100], primer2[100];
      if (batch_mode) {
         oligo1 = oligo2 = NULL;
         char *result = fgets(line, 200, stdin);
         if (result == NULL) break;
         line[strcspn(line, "\n")] = 0;
         int num = sscanf(line, "%s%s", primer1, primer2);
         if (num >= 1) oligo1 = (unsigned char *) primer1;
         if (num >= 2) oligo2 = (unsigned char *) primer2;
         if (num < 1) break;
         printf("BEGIN %s\n", line);
      }
      /* execute thermodynamical alignment */
      if(a.dimer==0 && oligo1!=NULL){
         thal(oligo1,oligo1,&a,mode,&o);
      } else if(a.dimer==0 && oligo1==NULL && oligo2!=NULL) {
         thal(oligo2,oligo2,&a,mode,&o);
      } else {
         thal(oligo1,oligo2,&a,mode,&o);
      }
      /* encountered error during thermodynamical calc */
      if (o.temp == THAL_ERROR_SCORE) {
         tmp_ret = fprintf(stderr, "Error: %s\n", o.msg);
         exit(-1);
      }
      if((mode == THL_FAST) || (mode == THL_DEBUG_F))
        printf("%f\n",o.temp);
      if (batch_mode) {
         printf("END %s\n", line);
         fflush(stdout);
      }
   } while (batch_mode);
   /* cleanup */
   destroy_thal_structures();
   thal_free_parameters(&thermodynamic_parameters);
   return EXIT_SUCCESS;
}
