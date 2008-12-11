/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008
Whitehead Institute for Biomedical Research, Steve Rozen
(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky
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

#define FILE_NAME_SIZE 80

/* Some function prototypes */
static void   print_usage();
static void   sig_handler(int);

/* Other global variables. */
static const char * pr_release = "primer3 release 2.0.0";
static const char *pr_program_name;

int
main(int argc, char *argv[]) { 
  /* Setup the input data structures handlers */
  int format_output = 0;
  int strict_tags = 0;
  int dump_args = 0 ; /* set to 1 if dumping arguments to choose_primers */
  int io_version = 4;

  /* Some space for file names */
  char *tmp_file_name = NULL;
  char p3_all_file[FILE_NAME_SIZE];
  char p3_settings_file[FILE_NAME_SIZE];

  p3_global_settings *global_pa;
  seq_args *sarg;
  read_boulder_record_results read_boulder_record_res = {0,0};

  pr_append_str fatal_parse_err;
  pr_append_str nonfatal_parse_err;
  
  /* Retval will point to the return value from choose_primers(). */
  p3retval *retval = NULL;
  int input_found=0;

  p3_all_file[0] = '\0';
  p3_settings_file[0] = '\0';

  init_pr_append_str(&fatal_parse_err);
  init_pr_append_str(&nonfatal_parse_err);


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
  if (dump_args) global_pa->dump = 1 ;
  
  /* Read in the flags provided with the program call */
  while (--argc > 0) {
    argv++;
    if (!strcmp(*argv, "-format_output")) {
      format_output = 1;
    } else if (!strcmp(*argv, "-about")) {
          printf( "%s\n", pr_release);
          exit (0);
    } else if (!strcmp(*argv, "-2x_compat")) {
          printf( "PRIMER_ERROR=flag -2x_compat is no longer supported\n=\n");
          exit (-1);
    } else if (!strcmp(*argv, "-io_version=3")) {
          io_version = 3;
    } else if (!strcmp(*argv, "-io_version=4")) {
          io_version = 4;
    } else if (!strncmp(*argv, "-p3_settings_file=", 18)) {
      tmp_file_name = strchr(*argv,'=') + 1;
      strncpy (p3_settings_file,tmp_file_name,FILE_NAME_SIZE-1);
    } else if (!strcmp(*argv, "-strict_tags")) {
      strict_tags = 1;
    } else  {
      print_usage();
      exit(-1);
    }
  }

  /* Settings files have to be read in just below, and
     the functions need a temporary sarg */
  if (!(sarg = create_seq_arg())) {
    exit(-2);
  }

  /* Read data from the settings file until a "=" line occurs.  Assign parameter
   * values for primer picking to pa and sa. */
  if (p3_settings_file[0] != '\0') {
    read_p3_file(p3_settings_file, settings, global_pa, 
                 sarg, &fatal_parse_err, &nonfatal_parse_err,
                 &read_boulder_record_res);
  }

  /* We also need to print out errors here because the loop erases all
   *  errors at start */
  /* If there are fatal errors, write the proper message and exit */
  if (fatal_parse_err.data != NULL) {
    if (format_output) {
        format_error(stdout, sarg->sequence_name, fatal_parse_err.data);
    } else {
        print_boulder_error(fatal_parse_err.data);
    }
    fprintf(stderr, "%s: %s\n", 
              pr_program_name, fatal_parse_err.data);
    destroy_seq_args(sarg);
    exit(-4);
  }

  /* If there are nonfatal errors, write the proper message
   * and skip to the end of the loop */
  if (!pr_is_empty(&nonfatal_parse_err)) {
    if (format_output) {
        format_error(stdout, sarg->sequence_name, 
                     nonfatal_parse_err.data);
    } else {
        print_boulder_error(nonfatal_parse_err.data);
    }
  }

  /* The temporary sarg is not needed any more */
  destroy_seq_args(sarg);
  sarg = NULL;
  
  /* Read the data from input stream record by record and process it if
   * there are no errors. This is where the work is done. */
  while (1) {
    /* Create and initialize a seq_args data structure. sa (seq_args *) is 
     * initialized here because Values are _not_ retained across different
     * input records. */
    if (!(sarg = create_seq_arg())) {
      exit(-2);
    }

    /* Reset all errors handlers and the return structure */
    pr_set_empty(&fatal_parse_err);
    pr_set_empty(&nonfatal_parse_err);
    retval = NULL;

    /* See read_boulder.h for documentation on read_boulder_record().*/
    if (!read_boulder_record(stdin, 
                            &strict_tags, 
                            &io_version,
                            !format_output, 
                            all_parameters,
                            global_pa, 
                            sarg, 
                            &fatal_parse_err, 
                            &nonfatal_parse_err,
                            &read_boulder_record_res)) {
      break; /* There were no more boulder records */
    }
    
    input_found = 1;
    if ((global_pa->primer_task == pick_detection_primers) 
                && (global_pa->pick_internal_oligo == 1)){
      PR_ASSERT(global_pa->pick_internal_oligo);
    }

    /* If there are fatal errors, write the proper message and exit */
    if (fatal_parse_err.data != NULL) {
      if (format_output) {
        format_error(stdout, sarg->sequence_name, fatal_parse_err.data);
      } else {
        print_boulder_error(fatal_parse_err.data);
      }
      fprintf(stderr, "%s: %s\n", 
              pr_program_name, fatal_parse_err.data);
      destroy_p3retval(retval);
      destroy_seq_args(sarg);
      exit(-4);
    }

    /* If there are nonfatal errors, write the proper message
     * and skip to the end of the loop */
    if (!pr_is_empty(&nonfatal_parse_err)) {
      if (format_output) {
        format_error(stdout, sarg->sequence_name, 
                     nonfatal_parse_err.data);
      } else {
        print_boulder_error(nonfatal_parse_err.data);
      }
      goto loop_wrap_up;
    }
    
    if (read_boulder_record_res.file_flag && sarg->sequence_name == NULL) {
      /* We will not have a base name for the files */
      if (format_output) {
        format_error(stdout, NULL, 
                     "Need PRIMER_SEQUENCE_ID if PRIMER_FILE_FLAG is not 0");
      } else {
        print_boulder_error("Need PRIMER_SEQUENCE_ID if PRIMER_FILE_FLAG is not 0");
      }
      goto loop_wrap_up;
    }

    /* Pick the primers - the central function */
    p3_set_gs_primer_file_flag(global_pa, 
                               read_boulder_record_res.file_flag);
    retval = choose_primers(global_pa, sarg);  
    if (NULL == retval) exit(-2); /* Out of memory. */
    
    /* This is old code to make it compatible with version 3 input.
       In future versions it can be deleted!
       If it was necessary to use a left_input, right_input,
       or internal_oligo_input primer that was
       unacceptable, then add warnings. */

    if (global_pa->pick_anyway && (io_version == 3
                || format_output)) {
      if (sarg->left_input) {
        add_must_use_warnings(&retval->warnings,
                              "Left primer", &retval->fwd.expl);
      }
      if (sarg->right_input) {
        add_must_use_warnings(&retval->warnings,
                              "Right primer", &retval->rev.expl);
      }
      if (sarg->internal_input) {
        add_must_use_warnings(&retval->warnings,
                              "Hybridization probe", &retval->intl.expl);
      }
    }
    /* End of the old code for compartibility. */

    if (pr_is_empty(&retval->glob_err)
        && pr_is_empty(&retval->per_sequence_err)) {
      /* We need to test for errors before we call
         p3_print_oligo_lists. This function only works on retval as
         returned when there were no errors. */
      if (read_boulder_record_res.file_flag) {
        /* Create files with left, right, and internal oligos. */
        p3_print_oligo_lists(retval, sarg, global_pa,
                             &retval->per_sequence_err,
                             sarg->sequence_name);
      }
    }

    if (format_output) {
      print_format_output(stdout, &io_version, global_pa, 
                          sarg, retval, pr_release,
                          read_boulder_record_res.explain_flag);
    } else {
      /* Use boulder output */
      print_boulder(/* & */io_version, global_pa, sarg, retval, 
                    read_boulder_record_res.explain_flag);
    }

  loop_wrap_up: /* Here the failed loops join in again */
    if (NULL != retval) {
      /* Check for errors and print them */
      if (NULL != retval->glob_err.data) {
        fprintf(stderr, "%s: %s\n", pr_program_name, retval->glob_err.data);
        destroy_p3retval(retval);
        destroy_seq_args(sarg);
        exit(-4);
      }
    }

    destroy_p3retval(retval); /* This works even if retval is NULL */
    retval = NULL;
    destroy_seq_args(sarg);
    sarg = NULL;

  }   /* while (1) (processing boulder io records) ...
         End of the primary working loop */

  /* To avoid being distracted when looking for leaks: */
  p3_destroy_global_settings(global_pa);
  global_pa = NULL;
  destroy_seq_args(sarg);
  destroy_pr_append_str_data(&nonfatal_parse_err);
  destroy_pr_append_str_data(&fatal_parse_err);
  
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

  fprintf(stderr, "\n\nUSAGE: %s %s %s %s %s\n", pr_program_name,
          "[-format_output]", "[-io_version=3|-io_version=4]", "[-p3_settings_file=<file_path>]", "[-strict_tags]");
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
