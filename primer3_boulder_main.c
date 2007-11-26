/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007
Whitehead Institute for Biomedical Research, Steve Rozen
(http://jura.wi.mit.edu/rozen), and Helen Skaletsky
All rights reserved.

    This file is part of primer3 and the libprimer3 library.

    Primer3 and the libprimer3 library are free software;
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

#include <signal.h>
#include <string.h> /* strlen(), memset(), strcmp() */
#include <stdlib.h> /* free() */
#include "format_output.h"
#include "libprimer3.h"
#include "read_boulder.h"
#include "print_boulder.h"

static void   print_usage();
static void   sig_handler(int);

/* Other global variables. */
static const char * pr_release = "primer3 release 1.1.2";
static const char *pr_program_name;

int
main(argc,argv)
    int argc;
    char *argv[]; 
{ 
  program_args prog_args;
  int         format_output = 0;
  primer_args *global_pa;
  seq_args *sa;
  pr_append_str *fatal_parse_err = NULL;
  pr_append_str *nonfatal_parse_err = NULL;
  pr_append_str *combined_retval_err = NULL;
  p3retval *retval = NULL;
  oligo_type oligot = OT_LEFT; /* Silence warning */
  int num_oligo = 0;
  primer_rec *oligo = NULL;
 
  int input_found=0;

  pr_program_name = argv[0];
  p3_set_program_name(pr_program_name);

  /* 
   * We set up some signal handlers in case someone starts up the program
   * from the command line, wonders why nothing is happening, and then kills
   * the program.
   */
  signal(SIGINT, sig_handler);
  signal(SIGTERM, sig_handler);

  /* 
   * We allocate the following structures on the heap rather than on the
   * stack in order to take advantage of memory access checking
   * during testing.
   */
  if (!(global_pa = malloc(sizeof(*global_pa)))) {
    exit(-2); /* Out of memory. */
  }

  pr_set_default_global_args(global_pa);
  memset(&prog_args, 0, sizeof(prog_args));

  while (--argc > 0) {
    argv++;
    if (!strcmp(*argv, "-format_output")) {
      prog_args.format_output = 1;
      format_output = 1;
    }  else if (!strcmp(*argv, "-2x_compat")) {
      printf( "PRIMER_ERROR=flag -2x_compat is no longer supported\n=\n");
      exit (-1);
    } else if (!strcmp(*argv, "-strict_tags"))
      prog_args.strict_tags = 1;
    else  {
      print_usage();
      exit(-1);
    }
  }

  fatal_parse_err    = create_pr_append_str();
  nonfatal_parse_err = create_pr_append_str();
  combined_retval_err = create_pr_append_str();
  if (NULL == fatal_parse_err 
      || NULL == nonfatal_parse_err
      || NULL == combined_retval_err) {
    exit(-2); /* Out of memory */
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

    pr_set_empty(fatal_parse_err);
    pr_set_empty(nonfatal_parse_err);
    retval = NULL;

    if (read_record(&prog_args, !format_output, global_pa, sa, 
		    fatal_parse_err, nonfatal_parse_err)
	<= 0) {
      break;
    }
    input_found = 1;


    if (fatal_parse_err->data != NULL) {
      if (format_output) {
	format_error(stdout, sa->sequence_name, fatal_parse_err->data);
      } else {
	boulder_print_error(fatal_parse_err->data);
      }
      fprintf(stderr, "%s: %s\n", 
	      pr_program_name, fatal_parse_err->data);
      exit(-4);
    }

    /* FIX ME -- read in mispriming libraries here */


    if (!pr_is_empty(nonfatal_parse_err)) {
      if (format_output) {
	format_error(stdout, sa->sequence_name, 
		   nonfatal_parse_err->data);
      } else {
	boulder_print_error(nonfatal_parse_err->data);
      }
      goto finish_loop;
    }

    /* FIX ME create retval inside choose_primers */
    if (!(retval = create_p3retval())) {
      exit(-2);
    }
    retval = choose_primers(retval, global_pa, sa);
    if (NULL == retval) exit(-2); /* Out of memory. */

    if (!pr_is_empty(&retval->glob_err)
	||
	!pr_is_empty(&retval->per_sequence_err)
	/* ||
	   !pr_is_empty(&sa->error) */ ) {
      pr_append_new_chunk(combined_retval_err, 
			  retval->glob_err.data);
      pr_append_new_chunk(combined_retval_err, 
			  retval->per_sequence_err.data);
      /* pr_append_new_chunk(combined_retval_err, 
	 sa->error.data); */
      if (format_output) {
	format_error(stdout, sa->sequence_name,
		     combined_retval_err->data);
      } else {
	boulder_print_error(combined_retval_err->data);
      }
      goto finish_loop;
    }

    /* PR_ASSERT(pr_is_empty(&sa->error)) */
    PR_ASSERT(pr_is_empty(&retval->glob_err))
    PR_ASSERT(pr_is_empty(&retval->per_sequence_err))

    if (global_pa->pick_left_primer && global_pa->pick_right_primer) {
      if (format_output) {
	format_pairs(stdout, global_pa, sa, 
		     &retval->best_pairs, 
		     pr_release);
      }
      else {
	boulder_print_pairs(&prog_args, global_pa, sa,
			    &retval->best_pairs);
      }
    } else {
      if (global_pa->pick_left_primer) {
	oligot = OT_LEFT;
	oligo = retval->f;
	num_oligo = retval->n_f;
      } else if (global_pa->pick_right_primer) {
	oligot = OT_RIGHT;
	oligo = retval->r;
	num_oligo = retval->n_r;
      } else if (global_pa->pick_internal_oligo) {
	oligot = OT_INTL;
	oligo = retval->mid;
	num_oligo = retval->n_m;
      } else {
	fprintf(stderr, "%s: fatal programming error\n", pr_program_name);
	abort();
      }

      if (format_output) {
	format_oligos(stdout, global_pa, sa, oligo,
		      num_oligo, oligot, pr_release);
      } else {
	boulder_print_oligos(global_pa, sa, num_oligo,
			     oligot, oligo);
      }
    }


  finish_loop:
    if (NULL != retval) {
      if (NULL != retval->glob_err.data) {
	fprintf(stderr, "%s: %s\n", pr_program_name, retval->glob_err.data);
	exit(-4);
      }
    }
    destroy_p3retval(retval); /* This works even if retval is NULL */
    destroy_seq_args(sa);
  } /* while (1) (processing boulder io records) */

  /* To avoid being distracted when looking for leaks: */
  destroy_seq_lib(global_pa->p_args.repeat_lib);
  destroy_seq_lib(global_pa->o_args.repeat_lib);
  free(global_pa);
  destroy_pr_append_str(fatal_parse_err);
  destroy_pr_append_str(nonfatal_parse_err);
  destroy_seq_args(sa);
    
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
