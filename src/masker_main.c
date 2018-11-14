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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "libprimer3.h"

/* maximum number of variables allowed in the probability formula */
#define MAX_VARIABLES 100

/* flags for choosing between hard and soft masking */
#define HARD_MASKING 0
#define SOFT_MASKING 1


static void print_help (int exit_value);
void print_parameters (masker_parameters *mp);

const char *pr_programme_name = "primer3_masker";

int
main (int argc, const char *argv[])
{
  unsigned int idx;
  char *end;
  int debug = 0;

  /* sequence file name specified by the user, otherwise STDIN will be used */
  const char *sequence_file_name = NULL;
  const char *lists_file_name = NULL;
  const char *kmer_lists_path = NULL;

  /* data structure for all k-mer lists used for masking */
  unsigned int nlists = 0;
  unsigned int nlist_parameters = 0;
  unsigned int npos = 0;
  unsigned int list_pos[MAX_VARIABLES], list_components[MAX_VARIABLES];

  masker_parameters mp = {};
  parameters_builder pbuilder = {};
  input_sequence *input_seq = NULL;

  pr_append_str parse_err;
  pr_append_str warnings;

  init_pr_append_str (&parse_err);
  init_pr_append_str (&warnings);

  /* fill mp with default parameters */
  mp.mdir = DEFAULT_MASKING_DIRECTION;
  mp.failure_rate = DEFAULT_FAILURE_RATE;
  mp.abs_cutoff = DEFAULT_ABS_CUTOFF;
  mp.nucl_masked_in_5p_direction = DEFAULT_M5P;
  mp.nucl_masked_in_3p_direction = DEFAULT_M3P;
  mp.print_sequence = PRINT_SEQUENCE;
  mp.do_soft_masking = HARD_MASKING;
  mp.masking_char = DEFAULT_MASK_CHAR;
  mp.list_prefix = DEFAULT_LIST_FILE_PREFIX;
  kmer_lists_path = "../kmer_lists/";
  if(argc < 2){
   print_help (0);
  }

  /* parsing and checking the commandline arguments */
  for (idx = 1; (int)idx < argc; idx++) {

      if (!strcmp (argv[idx], "-h") || !strcmp (argv[idx], "--help") || !strcmp (argv[idx], "-?") || (argc < 2)) {
      print_help (0);

    } else if ((int)idx == argc - 1 && argv[idx][0] != '-') {
      sequence_file_name = argv[idx];
      FILE *file = fopen(sequence_file_name, "r");
      if (file){
        fclose(file);
      } else {
        print_help(1);
      }

    } else if (!strcmp (argv[idx], "-lf") || !strcmp (argv[idx], "--lists_file")) {
      /* lists specified from the file */
      if (!argv[idx + 1] || argv[idx + 1][0] == '-') {
        pr_append_new_chunk_external (&warnings, "No lists file specified.");
        idx += 1;
        print_help(1);
      }
      lists_file_name = argv[idx + 1];
      idx += 1;

    } else if (!strcmp (argv[idx], "-l") || !strcmp (argv[idx], "--list")) {
      /* lists specified from the commandline */
      if (nlist_parameters == MAX_VARIABLES) {
        pr_append_new_chunk_external (&parse_err, "Maximum number of list variables reached.");
        print_help(1);
      }

      if (!argv[idx + 1]) {
        pr_append_new_chunk_external (&warnings, "No list name specified with -l parameter!.");
        print_help(1);
      }

      /* get the positions of list files */
      list_pos[nlist_parameters] = idx;
      while (argv[idx + 1]) {
        if (argv[idx + 1][0] == '-') {
          if (!argv[idx + 1][1]) break;
            strtod (argv[idx + 1] + 1, &end);
          if (*end != 0) break;
        } else if (idx + 1 != list_pos[nlist_parameters] + 1 && idx + 1 != list_pos[nlist_parameters] + 4) {
          strtod (argv[idx + 1], &end);
          if (*end != 0) break;
        } else if (idx + 1 == list_pos[nlist_parameters] + 4 && strcmp (argv[idx + 1], "sq")) {
          break;
        }
        idx += 1;
      }
      list_components[nlist_parameters] = idx - list_pos[nlist_parameters];
      nlist_parameters += 1;
      npos += 1;

    } else if (!strcmp (argv[idx], "-lp") || !strcmp (argv[idx], "--list_prefix")) {
      /* lists specified by the prefix */
      if (!argv[idx + 1] || argv[idx + 1][0] == '-') {
        pr_append_new_chunk_external (&warnings, "No list prefix specified! Using the default value.");
        idx += 1;
        print_help(1);
      }
      mp.list_prefix = (char *)argv[idx + 1];
      idx += 1;
    } else if (!strcmp (argv[idx], "-lh") || !strcmp (argv[idx], "--kmer_lists_path")) {
      /* kmer lists path */
      if (!argv[idx + 1] || argv[idx + 1][0] == '-') {
         pr_append_new_chunk_external (&warnings, "No list prefix specified! Using the default value.");
         idx += 1;
         print_help(1);
      }
         kmer_lists_path = (char *)argv[idx + 1];
         idx += 1;
    } else if (!strcmp (argv[idx], "-p") || !strcmp (argv[idx], "--probability_cutoff")) {
      if (!argv[idx + 1]) {
        pr_append_new_chunk_external (&warnings, "No cutoff value specified! Using the default value.");
        idx += 1;
        continue;
      }
      mp.failure_rate = strtod (argv[idx + 1], &end);
      mp.abs_cutoff = 0;
      if (*end != 0 || mp.failure_rate < 0 || mp.failure_rate > 1) {
        pr_append_new_chunk_external (&parse_err, "Invalid cutoff value: ");
        pr_append_external (&parse_err, argv[idx + 1]);
        /*free(parse_err.data);*/
        break;
      }
      idx += 1;

    } else if (!strcmp (argv[idx], "-a") || !strcmp (argv[idx], "--absolute_value_cutoff")) {
      if (!argv[idx + 1]) {
        pr_append_new_chunk_external (&warnings, "No absolute cutoff value specified! Using the default value.");
        idx += 1;
        break;
      }
      mp.abs_cutoff = strtod (argv[idx + 1], &end);
      mp.failure_rate = 0.0;
      if (*end != 0 || mp.abs_cutoff < 0) {
        pr_append_new_chunk_external (&parse_err, "Invalid absolute cutoff value: ");
        pr_append_external (&parse_err, argv[idx + 1]);
        break;
      }
      idx += 1;

    } else if (!strcmp (argv[idx], "-m5") || !strcmp (argv[idx], "--mask_5p")) {
      if (!argv[idx + 1]) {
        pr_append_new_chunk_external (&warnings, "Number of nucleotides masked in 5' direction not specified! Using the default value.");
        idx += 1;
        continue;
      }
      mp.nucl_masked_in_5p_direction = strtod (argv[idx + 1], &end);
      if (*end != 0) {
        pr_append_new_chunk_external (&parse_err, "Invalid number of nucleotides masked in 5' direction: ");
        pr_append_external (&parse_err, argv[idx + 1]);
        print_help(1);
      }
      idx += 1;

    } else if (!strcmp (argv[idx], "-m3") || !strcmp (argv[idx], "--mask_3p")) {
      if (!argv[idx + 1]) {
        pr_append_new_chunk_external (&warnings, "Number of nucleotides masked in 3' direction not specified! Using the default value.");
        idx += 1;
        continue;
      }
      mp.nucl_masked_in_3p_direction = strtod (argv[idx + 1], &end);
      if (*end != 0) {
        pr_append_new_chunk_external (&parse_err, "Invalid number of nucleotides masked in 3' direction: ");
        pr_append_external (&parse_err, argv[idx + 1]);
        print_help(1);
      }
      idx += 1;

    } else if (!strcmp (argv[idx], "-c") || !strcmp (argv[idx], "--masking_char")) {
      if (!argv[idx + 1 || argv[idx + 1][0] == '-']) {
        pr_append_new_chunk_external (&warnings, "Character for masking not specified! Using the default value.");
        idx += 1;
        continue;
      }
      mp.masking_char = argv[idx + 1][0];
      if (strlen(argv[idx + 1]) > 1 || mp.masking_char < 33 || mp.masking_char > 126) {
        pr_append_new_chunk_external (&parse_err, "Invalid character for masking: ");
        pr_append_external (&parse_err, argv[idx + 1]);
        print_help(1);
      }
      idx += 1;

    } else if (!strcmp (argv[idx], "-d") || !strcmp (argv[idx], "--masking_direction")) {
      if (!argv[idx + 1] || argv[idx + 1][0] == '-') {
        pr_append_new_chunk_external (&warnings, "Masking direction not specified! Masking both strands by default.");
      } else if (!strcmp (argv[idx + 1], "both")) {
        mp.mdir = both_on_same;
      } else if (!strcmp (argv[idx + 1], "fwd")) {
        mp.mdir = fwd;
      } else if (!strcmp (argv[idx + 1], "rev")) {
        mp.mdir = rev;
      } else {
        pr_append_new_chunk_external (&warnings, "Unknown masking direction: ");
        pr_append_external (&warnings, argv[idx + 1]);
        pr_append_external (&warnings, ". Masking both strands by default.");
      }
      idx += 1;

    } else if (!strcmp (argv[idx], "-s") || !strcmp (argv[idx], "--soft_mask")) {
      mp.do_soft_masking = SOFT_MASKING;

    } else if (!strcmp (argv[idx], "-D")) {
      debug += 1;

    } else {
      pr_append_new_chunk_external (&parse_err, "Unknown parameter: ");
      pr_append_external (&parse_err, argv[idx]);
      print_help(1);
    }
  }
        /*FILE *file = fopen(sequence_file_name, "r");
        if (file){
           fclose(file);
        } else {
            print_help(1);
  }
         */
  input_seq = create_input_sequence_from_file_name (sequence_file_name, &parse_err);

  if (parse_err.data != NULL) {
    fprintf(stderr, "%s -> parsing commandline: ERROR: %s\n", pr_programme_name, parse_err.data);
    exit(-1);
  }
  if (warnings.data != NULL) {
    fprintf(stderr, "%s -> parsing commandline: WARNING: %s\n", pr_programme_name, warnings.data);
  }

  if (lists_file_name) {
    /* if lists are given in a text file */
    mp.fp = read_formula_parameters_from_file (lists_file_name, &nlist_parameters, &pbuilder, &mp.formula_intercept, &parse_err);
    nlists = pbuilder.nfp;
  }

  if (npos != 0) {
    /* if lists are given by commandline arguments (can be added to the ones given in a file) */
    pbuilder.fp_array = mp.fp;

    for (idx = 0; idx < npos; idx++) {
      unsigned int pos = list_pos[idx] + 1;
      char *values[4];
      unsigned int nvalues = list_components[idx];
      char *end;
      memcpy (&values, &argv[pos], nvalues * sizeof(*argv));

      if (nvalues == 1) {
        double ic;
        double neg = 1.0;
        if (values[0][0] == '-') {
          values[0] += 1;
          neg = -1.0;
        }
        ic = strtod (values[0], &end);
        if (*end == 0) {
          mp.formula_intercept = ic * neg;
          continue;
        }
      }
      add_variable_to_formula_parameters (values, nvalues, &pbuilder, &parse_err);
    }
    nlists = pbuilder.nfp;
    mp.fp = pbuilder.fp_array;

  } else if (nlists == 0 && !lists_file_name) {
    /* if there are no lists specified use the default formula */
    mp.fp = create_default_formula_parameters (mp.list_prefix, kmer_lists_path, &parse_err);
    mp.formula_intercept = DEFAULT_INTERCEPT;
    nlists = DEFAULT_NLISTS;
    nlist_parameters = DEFAULT_NLIST_PARAMETERS;
  }

  mp.nlists = nlists;

  if (mp.abs_cutoff > 0 && nlist_parameters != 1) {
    fprintf (stderr, "Error: Absolute value cutoff works with one list and one k-mer frequency parameter only. Currently you are using %u lists and %u parameters.\n", nlists, nlist_parameters);
    print_help(1);
  }

  if (parse_err.data != NULL) {
    fprintf(stderr, "%s -> building formula: ERROR: %s\n", pr_programme_name, parse_err.data);
    delete_input_sequence (input_seq);
    exit(-1);
  }

  if (debug > 0) print_parameters (&mp);

  read_and_mask_sequence (input_seq, NULL, &mp, &parse_err, debug);

  if (parse_err.data != NULL) {
    fprintf(stderr, "%s -> masking sequence: ERROR: %s\n", pr_programme_name, parse_err.data);
    delete_input_sequence (input_seq);
    delete_formula_parameters (mp.fp, nlists);
    exit(-1);
  }

  destroy_pr_append_str_data (&warnings);
  destroy_pr_append_str_data (&parse_err);
  delete_input_sequence (input_seq);
  delete_formula_parameters (mp.fp, nlists);
  free(pbuilder.used_lists);

  if (debug > 0) fprintf (stderr, "Done!\n");
  return 0;
}

void
print_parameters (masker_parameters *mp)
{
  unsigned int i;
  fprintf (stderr, "Current masker parameters:\n");
  fprintf (stderr, "    masking_direction = %s\n", mp->mdir == both_on_same ? "both" : (mp->mdir == fwd ? "fwd" : "rev"));
  fprintf (stderr, "    failure_rate = %f\n", mp->failure_rate);
  fprintf (stderr, "    abs_cutoff = %u\n", mp->abs_cutoff);
  fprintf (stderr, "    m5p = %u\n", mp->nucl_masked_in_5p_direction);
  fprintf (stderr, "    m3p = %u\n", mp->nucl_masked_in_3p_direction);
  fprintf (stderr, "    print_sequence = %d\n", mp->print_sequence);
  fprintf (stderr, "    do_soft_masking = %d\n", mp->do_soft_masking);
  fprintf (stderr, "    masking_char = %c\n", mp->masking_char);
  fprintf (stderr, "    nlists = %u\n", mp->nlists);
  fprintf (stderr, "    intercept = %f\n", mp->formula_intercept);
  for (i = 0; i < mp->nlists; i++) {
    fprintf (stderr, "    LIST: %u\n", i);
    formula_parameters *fp = mp->fp[i];
    if (!fp) continue;
    fprintf (stderr, "        name = %s\n", fp->list_file_name);
    fprintf (stderr, "        oligo_length = %u\n", fp->oligo_length);
    fprintf (stderr, "        words_in_list = %llu\n", fp->words_in_list);
    fprintf (stderr, "        mm0 %f mm1 %f mm2 %f mm0_2 %f mm1_2 %f mm2_2 %f\n", fp->mm0, fp->mm1, fp->mm2, fp->mm0_2, fp->mm1_2, fp->mm2_2);
  }
  return;

}

static void
print_help (int exit_value)
{
  fprintf (stdout, "Usage: ./primer3_masker [OPTIONS] <INPUTFILE>\n");
  fprintf (stdout, "Options:\n");
  fprintf (stdout, "    -h, --help                   - print this usage screen and exit\n\n");
  fprintf (stdout, "    -p, --probability_cutoff     - masking cutoff [0, 1] (default: >=0.1)\n");
  fprintf (stdout, "    -lh, --kmer_lists_path       - path to the kmer list files (default: ../kmer_lists/)\n");
  fprintf (stdout, "    -lp, --list_prefix           - prefix of the k-mer lists to use with default model (default: homo_sapiens)\n\n");
  fprintf (stdout, "    -a, --absolute_value_cutoff  - masking cutoff based on k-mer count; requires a single list name, defined with -l\n");
  fprintf (stdout, "    -l, --list                   - define a single k-mer list for masking with absolute cutoff option -a\n\n");
  fprintf (stdout, "    -m5, --mask_5p               - nucleotides to mask in 5' direction (default: 1)\n");
  fprintf (stdout, "    -m3, --mask_3p               - nucleotides to mask in 3' direction (default: 0)\n");
  fprintf (stdout, "    -c, --masking_char           - character used for masking (default: N)\n");
  fprintf (stdout, "    -s, --soft_mask              - use soft masking (default: false)\n");
  fprintf (stdout, "    -d, --masking_direction      - a strand to mask (fwd, rev, both) (default: both)\n\n");
  exit (exit_value);
}
