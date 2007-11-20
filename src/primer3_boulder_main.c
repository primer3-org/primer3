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

#include <signal.h>
#include <string.h> /* strlen, memset, strcmp */
#include <stdlib.h> /* free */
#include "format_output.h"
#include "libprimer3.h"
#include "read_boulder.h"
#include "print_boulder.h"

static void   print_usage();
static void   sig_handler(int);

/* FIX ME -- deal with program name. */

/* Other global variables. */
static const char * pr_release = "primer3 release 1.1.2";
static const char *pr_program_name;
static int pr_program_name_len;

int
main(argc,argv)
    int argc;
    char *argv[]; 
{ 
  program_args prog_args;
  primer_args *global_pa;
  seq_args *sa;
  p3retval *retval = NULL;

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
  if (!(global_pa = malloc(sizeof(*global_pa)))) {
    exit(-2); /* Out of memory. */
  }

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

    /* Values for sa (seq_args *) are _not_ retained across different
       input records. */
    if (!(sa = malloc(sizeof(*sa)))) {
      exit(-2); /* Out of memory. */
    }

    if (read_record(&prog_args, global_pa, sa) <= 0) {
      free(sa);  /* free(s) is ok, because a return of 0 
		    indicates end-of-file, in which
		    case there are no pointers from sa .*/
      break;
    }

    /* We need to create the retval even if we are not going to call
       choose_primers because of user errors discovered in
       read_record().  This in turn is because (1) we count on
       boulder_print_pairs to print out the error tag and the final =,
       and (2) we count on format_output to print the error when
       prog_args.format_output is set.  We clean this up if we
       continue to provide boulder IO output and
       'format_{pairs,oligos}' for more than the next couple of
       releases. */
    if (!(retval = create_p3retval())) {
      exit(-2);
    }

    input_found = 1;

    /* Theoretically it would be possible to correct
       errorneous global input in a subsequent record, but
       way too complicated for the payoff. */
    if (NULL == sa->error.data && NULL == global_pa->glob_err.data) {
      retval = choose_primers(retval, global_pa, sa);
      if (NULL == retval) exit(-2);
    }

    if (NULL != global_pa->glob_err.data) 
      pr_append_new_chunk(&sa->error, global_pa->glob_err.data);

    if (pick_pcr_primers == global_pa->primer_task
	|| pick_pcr_primers_and_hyb_probe == global_pa->primer_task) {
      if (prog_args.format_output) {
	format_pairs(stdout, global_pa, sa, 
		     &retval->best_pairs, pr_release);
      }
      else {
	boulder_print_pairs(&prog_args, global_pa, sa, &retval->best_pairs);
      }
    } else if(global_pa->primer_task == pick_left_only) {
      if (prog_args.format_output) 
	format_oligos(stdout, global_pa, sa, retval->f,
		      retval->n_f, OT_LEFT, pr_release);
      else 
	boulder_print_oligos(global_pa, sa, retval->n_f, 
			     OT_LEFT, retval->f, retval->r, retval->mid);
    } else if(global_pa->primer_task == pick_right_only) {
      if (prog_args.format_output) 
	format_oligos(stdout, global_pa, sa, retval->r,
		      retval->n_r, OT_RIGHT,
		      pr_release);
      else 
	boulder_print_oligos(global_pa, sa, retval->n_r, OT_RIGHT,
				retval->f, retval->r, retval->mid);
    }
    else if(global_pa->primer_task == pick_hyb_probe_only) {
      if(prog_args.format_output) 
	format_oligos(stdout, global_pa, sa, retval->mid,
		      retval->n_m, OT_INTL, pr_release);
      else 
	boulder_print_oligos(global_pa, sa, retval->n_m, 
			     OT_INTL, retval->f, retval->r, retval->mid);
    }

    if (NULL != global_pa->glob_err.data) {
      fprintf(stderr, "%s: %s\n", pr_program_name, global_pa->glob_err.data);
      free(global_pa->glob_err.data);
      exit(-4);
    }

    destroy_p3retval(retval);
    destroy_seq_args(sa);
  }

  /* To avoid being distracted when looking for leaks: */
  destroy_seq_lib(global_pa->p_args.repeat_lib);
  destroy_seq_lib(global_pa->o_args.repeat_lib);
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
    const char **p = libprimer3_copyright();
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
