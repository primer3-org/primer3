/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008,2009,
              2010,2011,2012
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
#include "libprimer3.h"
#include "format_output.h"
#include "read_boulder.h"
#include "print_boulder.h"

/* Check on which OS we compile */
#if defined(_WIN32) || defined(WIN32) || defined (__WIN32__) || defined(__CYGWIN__) || defined(__MINGW32__)
#define OS_WIN
#endif

/* Some function prototypes */
static void   print_usage();
static void   sig_handler(int);
#if !defined(OS_WIN)
static void   validate_kmer_lists_path();
#endif

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
  int io_version = 4;
  int default_version = 2;
  int dump_args = 0;

  p3_global_settings *global_pa;
  seq_args *sarg;
  read_boulder_record_results read_boulder_record_res = {0,0};

  pr_append_str p3_settings_path;
  pr_append_str output_path;
  pr_append_str error_path;
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
    {"default_version", required_argument, 0, 'd'},
    {"Dump_args", no_argument, 0, 'D'},
    {"2x_compat", no_argument, 0, '2'},
    {"output", required_argument, 0, 'o'},
    {"error", required_argument, 0, 'e'},
    {0, 0, 0, 0}
  };
  int about = 0, compat = 0, invalid_flag = 0;

  /* Retval will point to the return value from choose_primers(). */
  p3retval *retval = NULL;
  int input_found=0;

  init_pr_append_str(&fatal_parse_err);
  init_pr_append_str(&nonfatal_parse_err);
  init_pr_append_str(&warnings);
  init_pr_append_str(&p3_settings_path);
  init_pr_append_str(&output_path);
  init_pr_append_str(&error_path);

  /* Get the program name for correct error messages */
  pr_release = libprimer3_release();
  pr_program_name = argv[0];
  p3_set_program_name(pr_program_name);

  /* We set up some signal handlers in case someone starts up the program
   * from the command line, wonders why nothing is happening, and then kills
   * the program. */
  signal(SIGINT, sig_handler);
  signal(SIGTERM, sig_handler);

  /* Read in the flags provided with the program call */
  opterr = 0;
  while ((opt = getopt_long_only(argc, argv, "", long_options, &option_index)) != -1) {
    switch (opt) {
    case 'a':
      about = 1;
      break;
    case 'p':
      if (pr_append_external(&p3_settings_path, optarg)) {
        exit(-2); /* Out of memory. */
      }
      break;
    case 'i':
      if (!strcmp(optarg, "4"))
        io_version = 4;
      else
        io_version = -1;
      break;

    case 'd': /* default_version */
      if (!strcmp(optarg, "1"))
        default_version = 1;
      else if (!strcmp(optarg, "2"))
        default_version = 2;
      else
        default_version = -1;
      break;

    case 'D':  /* Undocumented flag for testing; causes 
                  values of arguments to be echoed to
                  stdout. */
      dump_args = 1;
      break;
    case '2':
      compat = 1;
      break;
    case 'o':
      if (pr_append_external(&output_path, optarg)) {
        exit(-2); /* Out of memory. */
      }

      break;
    case 'e':
      if (pr_append_external(&error_path, optarg)) {
        exit(-2); /* Out of memory. */
      }

      break;
    case '?':
      invalid_flag = 1;
      break;
    }
  }
  /* Open the output and error files specified */

  if (!pr_is_empty(&error_path)) {
    /* reassign stderr */
    if (freopen(pr_append_str_chars(&error_path), "w", stderr)
        == NULL) {
      fprintf(stderr, "Error creating file %s\n",
              pr_append_str_chars(&error_path));
      exit(-1);
    }
    destroy_pr_append_str_data(&error_path);
  }

  if (!pr_is_empty(&output_path)) {
    /* reassign stdout */
    if (freopen(pr_append_str_chars(&output_path), "w", stdout)
        == NULL) {
      fprintf(stderr, "Error creating file %s\n",
              pr_append_str_chars(&output_path));
      exit(-1);
    }
    destroy_pr_append_str_data(&output_path);
  }
  /* We do any printing after redirecting stdout and stderr */
  if (about == 1) {
    printf("%s\n", pr_release);
    exit(0);
  }
  if ((io_version == -1) || (invalid_flag == 1) || (default_version == -1)) {
    print_usage();
    exit(-1);
  }
  if (compat == 1) {
    printf("PRIMER_ERROR=flag -2x_compat is no longer supported\n=\n");
    exit(-1);
  }

  /* Check if an input file has been specified on the command line. */
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

  /* Allocate the space for global settings and fill in default parameters */
  if (default_version == 1)
    global_pa = p3_create_global_settings_default_version_1();
  else if (default_version == 2) 
    global_pa = p3_create_global_settings();
  else {
    print_usage();
    exit(-1);
  }
  /* Load default thal parameters */
  thal_results o;
  if (get_thermodynamic_values(&global_pa->thermodynamic_parameters, &o)) {
    fprintf(stderr, "%s\n", o.msg);
    exit(-1);
  }

  if (!global_pa) {
    exit(-2); /* Out of memory. */
  }
  global_pa->dump = dump_args ;

  /* Settings files have to be read in just below, and
     the functions need a temporary sarg */
  if (!(sarg = create_seq_arg())) {
    exit(-2);
  }

  /* Read data from the settings file until a "=" line occurs.  Assign
     parameter values for primer picking to pa and sa. */
  if (!pr_is_empty(&p3_settings_path)) {
    read_p3_file(pr_append_str_chars(&p3_settings_path),
                 settings,
                 echo_settings && !format_output,
                 strict_tags,
                 global_pa, sarg, &fatal_parse_err,
                 &nonfatal_parse_err, &warnings, &read_boulder_record_res);
    destroy_pr_append_str_data(&p3_settings_path);
    /* Check if masking template flag was given */
#if !defined(OS_WIN)
    if (global_pa->mask_template == 1)
       validate_kmer_lists_path();
#endif
   }     
  /* We also need to print out errors here because the loop erases all
     errors at start. If there are fatal errors, write the proper
     message and exit */
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
    
#if !defined(OS_WIN)
    if(global_pa->mask_template){
           global_pa->lowercase_masking=global_pa->mask_template;
     }        

     /* Check if template masking flag was given */
     if (global_pa->mask_template == 1)
        validate_kmer_lists_path();

    /* Check that we found the kmer lists in case masking flag was set to 1. */
    if (global_pa->mask_template == 1 && kmer_lists_path == NULL){
        printf("PRIMER_ERROR=masking template chosen, but path to kmer lists not specified\n=\n");
        exit(-1);
    }
    /* Set up some masking parameters */
    /* edited by M. Lepamets */
   if (global_pa->mask_template == 1) {
      global_pa->mp.window_size = DEFAULT_WORD_LEN_2;
                                          
      if (global_pa->pick_right_primer == 0) global_pa->mp.mdir = fwd;
      else if (global_pa->pick_left_primer == 0) global_pa->mp.mdir = rev;
      /* Check if masking parameters (k-mer list usage) have changed */
      if (global_pa->masking_parameters_changed == 1) {
        delete_formula_parameters (global_pa->mp.fp, global_pa->mp.nlists);
        global_pa->mp.fp = create_default_formula_parameters (global_pa->mp.list_prefix, kmer_lists_path, &fatal_parse_err);
        global_pa->masking_parameters_changed = 0;
      }
    }
#endif
    input_found = 1;
    if ((global_pa->primer_task == generic)
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

    /* If it was necessary to use a left_input, right_input,
       or internal_oligo_input primer that was
       unacceptable, then add warnings. */
    if (global_pa->pick_anyway && format_output) {
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
      print_boulder(io_version, global_pa, sarg, retval,
                    read_boulder_record_res.explain_flag);
    }

  loop_wrap_up: /* Here the failed loops join in again */
    if (NULL != retval) {
      /* Check for errors and print them */
      if (NULL != retval->glob_err.data) {
        fprintf(stderr, "%s: %s\n", pr_program_name, retval->glob_err.data);
        destroy_secundary_structures(global_pa, retval);
        destroy_p3retval(retval);
        destroy_seq_args(sarg);
        exit(-4);
      }
    }
    destroy_secundary_structures(global_pa, retval); /* This works even if retval is NULL */
    destroy_p3retval(retval); /* This works even if retval is NULL */
    retval = NULL;
    destroy_seq_args(sarg);
    sarg = NULL;

  }   /* while (1) (processing boulder io records) ...
         End of the primary working loop */

  /* To avoid being distracted when looking for leaks: */
  destroy_thal_structures();

#if !defined(OS_WIN)  
  if(global_pa->mask_template == 1){
    delete_formula_parameters (global_pa->mp.fp, global_pa->mp.nlists);
    /* free(global_pa->mp.list_prefix); */
  }
#endif
     
  p3_destroy_global_settings(global_pa);
  global_pa = NULL;
  destroy_seq_args(sarg);
  destroy_pr_append_str_data(&nonfatal_parse_err);
  destroy_pr_append_str_data(&fatal_parse_err);
  destroy_pr_append_str_data(&warnings);
  destroy_dpal_thal_arg_holder();
#if !defined(OS_WIN)  
  free(kmer_lists_path);
#endif
  /* If it could not read input, then complain and die */
  if (0 == input_found) {
    print_usage();
    exit(-3);
  }
  return 0;
}

#if 0
/* in windows check for ..\\kmer_lists */
static void validate_kmer_lists_path(){
   if (kmer_lists_path == NULL) {
      struct stat st;
      if ((stat("..\\kmer_lists", &st) == 0) && S_ISDIR(st.st_mode)) {
         kmer_lists_path =
           (char*) malloc(strlen("..\\kmer_lists\\") * sizeof(char) + 1);
         if (NULL == kmer_lists_path) exit (-2); /* Out of memory */
         strcpy(kmer_lists_path, "..\\kmer_lists\\");
      } else {
         /* no default directory found */
         return;
      }
   } 
}
#endif

#if !defined(OS_WIN)
/* in linux, check for ../kmer_lists and /opt/kmer_lists */
static void validate_kmer_lists_path(){
   if (kmer_lists_path == NULL) {
      struct stat st;
      if ((stat("../kmer_lists", &st) == 0) && S_ISDIR(st.st_mode)) {
         kmer_lists_path =
           (char*) malloc(strlen("../kmer_lists/") * sizeof(char) + 1);
         if (NULL == kmer_lists_path) exit (-2); /* Out of memory */
         strcpy(kmer_lists_path, "../kmer_lists/");
      } else if ((stat("/opt/kmer_lists", &st) == 0)  && S_ISDIR(st.st_mode)) {
         kmer_lists_path =
           (char*) malloc(strlen("/opt/kmer_lists/") * sizeof(char) + 1);
         if (NULL == kmer_lists_path) exit (-2); /* Out of memory */
         strcpy(kmer_lists_path, "/opt/kmer_lists/");
      } else {
         /* no default directory found */
         return;
      }
   }
}
#endif
      
   
/* Print out copyright and a short usage message*/
static void
print_usage()
{
  fprintf(stderr, "%s", primer3_copyright());

  fprintf(stderr, "\n\nUSAGE: %s %s %s %s %s %s %s %s %s %s\n", pr_program_name,
          "[--format_output]", 
          "[--default_version=1|--default_version=2]",
          "[--io_version=4]", 
          "[--p3_settings_file=<file_path>]",
          "[--echo_settings_file]",
          "[--strict_tags]", 
          "[--output=<file_path>]", 
          "[--error=<file_path>]",
          "[input_file]");
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
