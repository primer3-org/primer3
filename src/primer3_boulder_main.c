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
#include <ctype.h>
#include <string.h> /* strcmp() */
#include <stdlib.h> /* free() */
#include "format_output.h"
#include "libprimer3.h"
#include "read_boulder.h"
#include "print_boulder.h"

/* Some function prototypes */
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
  /* Setup the input data structures handlers */
  int format_output = 0;
  int strict_tags = 0;
  int io_version = 0;

  p3_global_settings *global_pa;
  seq_args *sa;

  /* Setup the error structures handlers */
  pr_append_str *fatal_parse_err = NULL;
  pr_append_str *nonfatal_parse_err = NULL;
  pr_append_str *adjust_seq_args_warnings = NULL;
  
  /* Setup the output data structure handlers */
  p3retval *retval = NULL;
  int input_found=0;

  /* Get the program name for correct error messages */
  pr_program_name = argv[0];
  p3_set_program_name(pr_program_name);

  /* We set up some signal handlers in case someone starts up the program
   * from the command line, wonders why nothing is happening, and then kills
   * the program. */
  signal(SIGINT, sig_handler);
  signal(SIGTERM, sig_handler);

  /* Allocate the space for global settings and fill in default parameters */
  global_pa = p3_create_global_settings();
  if (!global_pa) {
    exit(-2); /* Out of memory. */
  }
  
  /* Read in the flags provided with the program call */
  while (--argc > 0) {
    argv++;
    if (!strcmp(*argv, "-format_output")) {
      format_output = 1;
    } else if (!strcmp(*argv, "-2x_compat")) {
      printf( "PRIMER_ERROR=flag -2x_compat is no longer supported\n=\n");
      exit (-1);
    } else if (!strncmp(*argv, "-io_version=", 10)) {
      /* This reads in the version number required for extended io functions */
      /* There may be a better way, but it works */
      char tag2int[20];
      strncpy (tag2int,*argv,19);
      int counter = 12;
      while (isdigit(tag2int[counter])) {
	if (isdigit(tag2int[counter])) {
	  io_version=10*io_version+(tag2int[counter] - '0');
	}
	if ( counter > 20 ) {
	  break; /* Just to be save */
	}
	counter++;
      }
    } else if (!strcmp(*argv, "-strict_tags")) {
      strict_tags = 1;
    } else  {
      print_usage();
      exit(-1);
    }
  }

  /* Allocate the space for empty error messages */
  fatal_parse_err     = create_pr_append_str();
  nonfatal_parse_err  = create_pr_append_str();
  adjust_seq_args_warnings       = create_pr_append_str();
  if (NULL == fatal_parse_err 
      || NULL == nonfatal_parse_err
      || NULL == adjust_seq_args_warnings) {
    exit(-2); /* Out of memory */
  }

  /* Read the data from input stream record by record and process it if
   * there are no errors. This is were the work is done! */
  while (1) {

    /* Create and initialize a seq_args data structure. sa (seq_args *) is 
     * initialized here because Values are _not_ retained across different
     * input records. */
    if (!(sa = create_seq_arg())) {
      exit(-2);
    }

    /* Reset all errors handlers and the return structure */
    pr_set_empty(fatal_parse_err);
    pr_set_empty(nonfatal_parse_err);
    pr_set_empty(adjust_seq_args_warnings);
    retval = NULL;

    /* Read data from stdin until a "=" line occurs.  Assign parameter
     * values for primer picking to pa and sa. Perform initial data
     * checking. */
    if (read_record(&strict_tags, &io_version, !format_output, global_pa, sa, 
		    fatal_parse_err, nonfatal_parse_err) <= 0) {
      break; /* leave the program loop and complain later */
    }
    
    input_found = 1;
    if (global_pa->primer_task == pick_pcr_primers_and_hyb_probe) {
      PR_ASSERT(global_pa->pick_internal_oligo);
    }

    /* If there are fatal errors, write the proper message and exit */
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

    /* POSSIBLE CHANGE -- read in mispriming libraries here? */

    /* Modify some of the arguments */
    p3_adjust_seq_args(global_pa, sa, nonfatal_parse_err, adjust_seq_args_warnings);

    /* If there are nonfatal errors, write the proper message
     * and skip to the end of the loop */
    if (!pr_is_empty(nonfatal_parse_err)) {
      if (format_output) {
	format_error(stdout, sa->sequence_name, 
		     nonfatal_parse_err->data);
      } else {
	boulder_print_error(nonfatal_parse_err->data);
      }
      goto loop_wrap_up;
    }
    
    /* Pick the primers - the central function */
    retval = choose_primers(global_pa, sa);
    if (NULL == retval) exit(-2); /* Out of memory. */

    if (pr_is_empty(&retval->glob_err)
	&& pr_is_empty(&retval->per_sequence_err)) {
      /* We need to test for error first, because
	 p3_print_oligo_lists does not handle
	 partial outputs in retval.  FIX ME,
	 move test inside  */
      /* Create files with left, right, and internal oligos. */
      if (global_pa->file_flag) {  /* FIX ME, get file_flag out of global settings */
	p3_print_oligo_lists(retval, sa, global_pa,
			     &retval->per_sequence_err); 
      }
    }

    if (format_output) {
      print_format_output(stdout, &io_version, global_pa, 
			  sa, retval, pr_release, adjust_seq_args_warnings);
    } else {
      /* Use boulder output */
      boulder_print(&io_version, global_pa, sa, retval, adjust_seq_args_warnings);
    }

  loop_wrap_up: /* Here the failed loops join in again */
    if (NULL != retval) {
      /* Check for errors and print them */
      if (NULL != retval->glob_err.data) {
	fprintf(stderr, "%s: %s\n", pr_program_name, retval->glob_err.data);
	exit(-4);
      }
    }
    /* Delete the data structures out of the memory */
    destroy_p3retval(retval); /* This works even if retval is NULL */
    destroy_seq_args(sa);
  } 
  /* while (1) (processing boulder io records) 
   * End of the primary working loop */

  /* To avoid being distracted when looking for leaks: */
  p3_destroy_global_settings(global_pa);
  global_pa = NULL;
  destroy_pr_append_str(fatal_parse_err);
  destroy_pr_append_str(nonfatal_parse_err);
  destroy_pr_append_str(adjust_seq_args_warnings);
  destroy_seq_args(sa);
  
  /* If it could not read input complain and die */
  if (0 == input_found) {
    print_usage();
    exit(-3);
  }
  return 0;
}

/* Print out copyright and a short usage message*/
static void
print_usage()
{
  fprintf(stderr, primer3_copyright());

  fprintf(stderr, "\n\nUSAGE: %s %s %s %s\n", pr_program_name,
	  "[-format_output]", "[-io_version=xxx]", "[-strict_tags]");
  fprintf(stderr, "This is primer3 (%s)\n", pr_release);
  fprintf(stderr, "Input must be provided on standard input.\n");
  fprintf(stderr, "For example:\n");
  fprintf(stderr, "$ primer3_core < my_input_file\n");
}

/* Print out copyright, a short usage message and the signal */
static void
sig_handler(int signal) {
    print_usage();
    fprintf(stderr, "%s: received signal %d\n", pr_program_name, signal);
    exit(signal);
}
