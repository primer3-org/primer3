/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007
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

#include <errno.h>
#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <string.h>
#include <stdlib.h> /* free */
#include <unistd.h>
#include <float.h>
#include "primer3_release.h"
#include "format_output.h"
#include "dpal.h"
#include "oligotm.h"
#include "libprimer3.h"
#include "read_boulder.h"
#include "print_boulder.h"

#ifndef MAX_PRIMER_LENGTH
#error "Define MAX_PRIMER_LENGTH in Makefile..."
  /* to ensure that MAX_PRIMER_LENGTH <= DPAL_MAX_ALIGN. */
#endif
#if (MAX_PRIMER_LENGTH > DPAL_MAX_ALIGN) 
#error "MAX_PRIMER_LENGTH must be <= DPAL_MAX_ALIGN"
#endif

static void   print_usage();
static void   sig_handler(int);

/* Global static variables. */
static const char * copyright[] = {
"",
"Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007",
"Whitehead Institute for Biomedical Research, Steve Rozen",
"(http://jura.wi.mit.edu/rozen), and Helen Skaletsky",
"All rights reserved.",
"",
"Redistribution and use in source and binary forms, with or without",
"modification, are permitted provided that the following conditions are",
"met:",
"",
"   * Redistributions of source code must retain the above copyright",
"notice, this list of conditions and the following disclaimer.",
"   * Redistributions in binary form must reproduce the above",
"copyright notice, this list of conditions and the following disclaimer",
"in the documentation and/or other materials provided with the",
"distribution.",
"   * Neither the names of the copyright holders nor contributors may",
"be used to endorse or promote products derived from this software",
"without specific prior written permission.",
"",
"THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS",
"\"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT",
"LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR",
"A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT",
"OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,",
"SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT",
"LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,",
"DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY",
"THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT",
"(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE",
"OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.",
NULL
};

/* Other global variables. */
const char *pr_program_name;
int pr_program_name_len;

typedef struct oligo_array {
  int len;
  primer_rec *data;
} oligo_array;

int
main(argc,argv)
    int argc;
    char *argv[]; 
{ 
  program_args prog_args;
  primer_args *global_pa;
  seq_args *sa;
  primer3_state *p3state;

  int input_found=0;

  pr_program_name = argv[0];
  pr_program_name_len = strlen(argv[0]);
  /* 
   * We set up some signal handlers in case someone starts up the program
   * from the command line, wonders why nothing is happening, and then kills
   * the program.
   */
  signal(SIGINT, sig_handler);
  signal(SIGTERM, sig_handler);

  /* 
   * We allocate the following structures on the heap rather than on the
   * stack in order to take advantage of testcenter memory access checking
   * during testing.
   */
  global_pa = pr_safe_malloc(sizeof(*global_pa));

  pr_set_default_global_args(global_pa);
  memset(&prog_args, 0, sizeof(prog_args));

  while (--argc > 0) {
    argv++;
    if (!strcmp(*argv, "-format_output"))
      prog_args.format_output = 1;
    else if (!strcmp(*argv, "-2x_compat")) {
      /* prog_args.twox_compat = 1; */
      pr_append_new_chunk(&global_pa->glob_err, 
			  "flag -2x_compat is no longer supported");
    } else if (!strcmp(*argv, "-strict_tags"))
      prog_args.strict_tags = 1;
    else  {
      print_usage();
      exit(-1);
    }
  }

  /* 
   * Read the data from input stream record by record and process it if
   * there are no errors.
   */
  while (1) {
    sa = pr_safe_malloc(sizeof(*sa));

    if (read_record(&prog_args, global_pa, sa) <= 0) {
      free(sa);  /* Free ok, because there are no pointers
		    from sa at this point. */
      break;
    }

    input_found = 1;
    p3state = create_primer3_state();

    if (NULL == sa->error.data && NULL == global_pa->glob_err.data) {
      choose_primers(p3state, global_pa, sa);
    }

    if (NULL != global_pa->glob_err.data) 
      pr_append_new_chunk(&sa->error, global_pa->glob_err.data);

    if (pick_pcr_primers == global_pa->primer_task
	|| pick_pcr_primers_and_hyb_probe == global_pa->primer_task) {
      if (prog_args.format_output) {
	format_pairs(stdout, global_pa, sa, &p3state->best_pairs);
      }
      else {
	boulder_print_pairs(&prog_args, global_pa, sa, &p3state->best_pairs);
      }
    } else if(global_pa->primer_task == pick_left_only) {
      if (prog_args.format_output) 
	format_oligos(stdout, global_pa, sa, p3state->f, p3state->n_f, OT_LEFT);
      else boulder_print_oligos(global_pa, sa, p3state->n_f, OT_LEFT, p3state->f, p3state->r, p3state->mid);
    } else if(global_pa->primer_task == pick_right_only) {
      if (prog_args.format_output) 
	format_oligos(stdout, global_pa, sa, p3state->r, p3state->n_r, OT_RIGHT);
      else boulder_print_oligos(global_pa, sa, p3state->n_r, OT_RIGHT, p3state->f, p3state->r, p3state->mid);
    }
    else if(global_pa->primer_task == pick_hyb_probe_only) {
      if(prog_args.format_output) 
	format_oligos(stdout, global_pa, sa, p3state->mid, p3state->n_m, OT_INTL);
      else boulder_print_oligos(global_pa, sa, p3state->n_m, OT_INTL, p3state->f, p3state->r, p3state->mid);
    }

    if (NULL != global_pa->glob_err.data) {
      fprintf(stderr, "%s: %s\n", pr_program_name, global_pa->glob_err.data);
      free(global_pa->glob_err.data);
      exit(-4);
    }

    destroy_primer3_state(p3state);
    destroy_seq_args(sa);
  }

  /* To avoid being distracted when looking for leaks: */
  free(global_pa);
    
  if (0 == input_found) {
    print_usage();
    exit(-3);
  }
  return 0;
}

static void
print_usage()
{
    const char **p;
    p = &copyright[0];
    while (NULL != *p) fprintf(stderr, "%s\n", *p++);
    fprintf(stderr, 
	    "\n\nUSAGE: %s %s %s\n", pr_program_name,
	    "[-format_output]",
	    "[-strict_tags]");
    fprintf(stderr, "This is primer3 (%s)\n", pr_release);
    fprintf(stderr, "Input must be provided on standard input.\n");
    fprintf(stderr, "For example:\n");
    fprintf(stderr, "$ primer3_core < my_input_file\n");
}

static void
sig_handler(signal)
    int signal;
{
    print_usage();
    fprintf(stderr, "%s: received signal %d\n", pr_program_name, signal);
    exit(signal);
}
