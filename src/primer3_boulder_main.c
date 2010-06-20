/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008,2009,2010
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
#include <getopt.h> /* getopt() */
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "thal.h"
#include "format_output.h"
#include "libprimer3.h"
#include "read_boulder.h"
#include "print_boulder.h"

/* Check on which OS we compile */
#if defined(_WIN32) || defined(WIN32) || defined (__WIN32__) || defined(__CYGWIN__) || defined(__MINGW32__)
#define OS_WIN
#endif

#define FILE_NAME_SIZE 80

/* Some function prototypes */
static void   print_usage();
static void   sig_handler(int);
static void   read_thermodynamic_parameters();

/* Other global variables. */
static const char *pr_release;
static const char *pr_program_name;

int
main(int argc, char *argv[]) 
{ 
  /* Setup the input data structures handlers */
  int format_output = 0;
  int strict_tags = 0;
  int echo_settings = 0;
  int dump_args = 0 ; /* set to 1 if dumping arguments to choose_primers */
  int io_version = 4;

  /* Some space for file names */
  char p3_settings_file[FILE_NAME_SIZE];

  p3_global_settings *global_pa;
  seq_args *sarg;
  read_boulder_record_results read_boulder_record_res = {0,0};

  pr_append_str fatal_parse_err;
  pr_append_str nonfatal_parse_err;
  pr_append_str warnings;

  /* Some variables needed by getopt */
  int opt, option_index = 0;
  struct option long_options[] = {
    {"about", no_argument, 0, 'a'},
    {"format_output", no_argument, &format_output, 1},
    {"strict_tags", no_argument, &strict_tags, 1},
    {"p3_settings_file", required_argument, 0, 'p'},
    {"echo_settings_file", no_argument, &echo_settings, 1},
    {"io_version", required_argument, 0, 'i'},
    {"2x_compat", no_argument, 0, '2'},
    {"output", required_argument, 0, 'o'},
    {"error", required_argument, 0, 'e'},
    {0, 0, 0, 0}
  };
  int about = 0, output = 0, error = 0, compat = 0, invalid_flag = 0;
  char output_file[FILE_NAME_SIZE], error_file[FILE_NAME_SIZE]; 
 
  /* Retval will point to the return value from choose_primers(). */
  p3retval *retval = NULL;
  int input_found=0;

  p3_settings_file[0] = '\0';

  init_pr_append_str(&fatal_parse_err);
  init_pr_append_str(&nonfatal_parse_err);
  init_pr_append_str(&warnings);

  /* Get the program name for correct error messages */
  pr_release = libprimer3_release();
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
  opterr = 0;
  while ((opt = getopt_long_only(argc, argv, "", long_options, &option_index)) != -1) {
    switch (opt) {
    case 'a':
      about = 1;
      break;
    case 'p':
      strncpy(p3_settings_file, optarg, FILE_NAME_SIZE - 1);
      break;
    case 'i':
      if (!strcmp(optarg, "3"))
        io_version = 3;
      else if (!strcmp(optarg, "4"))
        io_version = 4;
      else
        io_version = -1;
      break;
    case '2':
      compat = 1;
      break;
    case 'o':
      output = 1;
      strncpy(output_file, optarg, FILE_NAME_SIZE - 1);
      break;
    case 'e':
      error = 1;
      strncpy(error_file, optarg, FILE_NAME_SIZE - 1);
      break;
    case '?':
      invalid_flag = 1;
      break;
    }
  }
  /* Open the output and error files specified */
  if (error == 1) {
    /* reassign stderr */
    if (freopen(error_file, "w", stderr) == NULL) {
      fprintf(stderr, "Error creating file %s\n", error_file);
      exit(-1);
    }
  }
  if (output == 1) {
    /* reassign stdout */
    if (freopen(output_file, "w", stdout) == NULL) {
      fprintf(stderr, "Error creating file %s\n", output_file);
      exit(-1);
    }
  }
  /* We do any printing after redirecting stdout and stderr */
  if (about == 1) {
    printf("%s\n", pr_release);
    exit(0);
  }
  if ((io_version == -1) || (invalid_flag == 1)) {
    print_usage();
    exit(-1);
  }
  if (compat == 1) {
    printf("PRIMER_ERROR=flag -2x_compat is no longer supported\n=\n");
    exit(-1);
  }
  /* Check if an input file has been specified */
  if (optind < argc) {
    if (optind + 1 != argc) {
      print_usage();
      exit(-1);
    }
    if (freopen(argv[optind], "r", stdin) == NULL) {
      fprintf(stderr, "Error opening file %s\n", argv[optind]);
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
    read_p3_file(p3_settings_file, settings, echo_settings && !format_output,
		 global_pa, sarg, &fatal_parse_err, 
		 &nonfatal_parse_err, &warnings, &read_boulder_record_res);
    /* Check if the thermodynamical alignment flag was given */ 
    if (global_pa->thermodynamic_alignment == 1)
      read_thermodynamic_parameters();
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
    pr_set_empty(&warnings);
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
			     &warnings,
			     &read_boulder_record_res)) {
      break; /* There were no more boulder records */
    }

    /* Check if the thermodynamical alignment flag was given and the
       path to the parameter files changed - we need to reread them */
    if ((global_pa->thermodynamic_alignment == 1)
	&& (thermodynamic_path_changed == 1))
      read_thermodynamic_parameters();
    
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

    /* Print any warnings and continue processing */
    if (!pr_is_empty(&warnings)) {
      if (format_output) {
        format_warning(stdout, sarg->sequence_name, 
		       warnings.data);
      } else {
        print_boulder_warning(warnings.data);
      }
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
    /* End of the old code for compatibility. */

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
  if (global_pa->thermodynamic_alignment == 1)
    destroy_thal_structures();
  p3_destroy_global_settings(global_pa);
  global_pa = NULL;
  destroy_seq_args(sarg);
  destroy_pr_append_str_data(&nonfatal_parse_err);
  destroy_pr_append_str_data(&fatal_parse_err);
  destroy_pr_append_str_data(&warnings);
  destroy_dpal_thal_arg_holder();
  free(thermodynamic_params_path);
  /* If it could not read input, then complain and die */
  if (0 == input_found) {
    print_usage();
    exit(-3);
  }
  return 0;
}

/* Reads the thermodynamic parameters if the thermodynamic alignment
   tag was set to 1 */
static void 
read_thermodynamic_parameters()
{
  thal_results o;
  /* if the path to the parameter files did not change, we do not want to read again */
  if (thermodynamic_path_changed == 0) return;
  /* check that the path to the parameters folder was given */
  if (thermodynamic_params_path == NULL) {
#ifdef OS_WIN
    /* in windows check for .\\primer3_config */
    struct stat st;
    if ((stat(".\\primer3_config", &st) == 0) && S_ISDIR(st.st_mode)) {
      thermodynamic_params_path = 
	(char*) malloc(strlen(".\\primer3_config\\") * sizeof(char) + 1);
      if (NULL == thermodynamic_params_path) exit (-2); /* Out of memory */
      strcpy(thermodynamic_params_path, ".\\primer3_config\\");
    } else {
      /* no default directory found, error */
      printf("PRIMER_ERROR=thermodynamic approach chosen, but path to thermodynamic parameters not specified\n=\n");
      exit(-1);
    }
#else
    /* in linux, check for ./primer3_config and /opt/primer3_config */
    struct stat st;
    if ((stat("./primer3_config", &st) == 0) && S_ISDIR(st.st_mode)) {
      thermodynamic_params_path = 
	(char*) malloc(strlen("./primer3_config/") * sizeof(char) + 1);
      if (NULL == thermodynamic_params_path) exit (-2); /* Out of memory */
      strcpy(thermodynamic_params_path, "./primer3_config/");
    } else if ((stat("/opt/primer3_config", &st) == 0)  && S_ISDIR(st.st_mode)) {
      thermodynamic_params_path = 
	(char*) malloc(strlen("/opt/primer3_config/") * sizeof(char) + 1);
      if (NULL == thermodynamic_params_path) exit (-2); /* Out of memory */
      strcpy(thermodynamic_params_path, "/opt/primer3_config/");
    } else {
      /* no default directory found, error */
      printf("PRIMER_ERROR=thermodynamic approach chosen, but path to thermodynamic parameters not specified\n=\n");	
      exit(-1);
    }
#endif
  }
  /* read in the thermodynamic parameters */
  if (get_thermodynamic_values(thermodynamic_params_path, &o)) {
    fprintf(stderr, "%s\n", o.msg);
    exit(-1);
  }
  /* mark that the last given path was used for reading the parameters */
  thermodynamic_path_changed = 0;
}

/* Print out copyright and a short usage message*/
static void
print_usage()
{
  fprintf(stderr, "%s", primer3_copyright());

  fprintf(stderr, "\n\nUSAGE: %s %s %s %s %s %s %s %s %s\n", pr_program_name,
          "[-format_output]", "[-io_version=3|-io_version=4]", "[-p3_settings_file=<file_path>]", "[-echo_settings_file]",
	  "[-strict_tags]", "[-output=<file_path>]", "[-error=<file_path>]", "[input_file]");
  fprintf(stderr, "This is primer3 (%s)\n", pr_release);
  fprintf(stderr, "Input can also be provided on standard input.\n");
  fprintf(stderr, "For example:\n");
  fprintf(stderr, "$ primer3_core < my_input_file\n");
}

/* Print out copyright, a short usage message and the signal */
static void
sig_handler(int signal) 
{
    print_usage();
    fprintf(stderr, "%s: received signal %d\n", pr_program_name, signal);
    exit(signal);
}
