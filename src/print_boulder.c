/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008
Whitehead Institute for Biomedical Research, Steve Rozen
(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky
All rights reserved.

    This file is part of primer3 and the primer3 suite.

    Primer3 and the primer3 suite are free software;
    you can redistribute them and/or modify them under the terms
    of the GNU General Public License as published by the Free
    Software Foundation; either version 2 of the License, or (at
    your option) any later version.

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this file (file gpl-2.0.txt in the source
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

#include <string.h> /* strcpy */
#include "print_boulder.h"

static const char *pr_program_name = "Program name is probably primer3_core";

static void   print_all_explain(const p3_global_settings *,
                                const seq_args *, 
                                const p3retval *, 
                                int   io_version);

/* Print the data for chosen primer pairs to stdout in "boulderio" format. */
void
print_boulder(int io_version,
              const p3_global_settings *pa,
              const seq_args *sa,
              const p3retval *retval,
              int   explain_flag) {
  /* The pointers to warning tag */
  char *warning;

  /* A place to put a string containing all error messages */
  pr_append_str *combined_retval_err = NULL;

  /* A small spacer; WARNING this is a fixed size
     buffer, but plenty bigger than
     log(2^64, 10), the longest character
     string that is needed for a 64 bit integer. */
  char suffix [100];

  /* Pointers for the primer set just printing */
  primer_rec *fwd, *rev, *intl;
 
  /* Variables only used for Primer Lists */
  int num_fwd, num_rev, num_int, num_pair, num_print;
  int print_fwd = 0;
  int print_rev = 0;
  int print_int = 0;
    
  /* Switches for printing this primer */
  int go_fwd = 0;
  int go_rev = 0;
  int go_int = 0;
    
  /* The number of loop cycles */
  int loop_max;
    
  /* That links to the included region */
  int i, incl_s = sa->incl_s;
    
  /* This deals with the renaming of the internal oligo */
  const char *new_oligo_name = "INTERNAL";
  char *int_oligo = (char*) new_oligo_name;

  /* Check: are all pointers linked to something*/
  PR_ASSERT(NULL != pa);
  PR_ASSERT(NULL != sa);
        
  /* Check if there are warnings and print them */
  if ((warning = p3_get_rv_and_gs_warnings(retval, pa)) != NULL) { 
    printf("PRIMER_WARNING=%s\n", warning);
    free(warning);
  }

  /* Check if a settings file was read an print its id 
   * Not needed anymore, we do this when we read the setting file. */
  /*
  if (pa->settings_file_id != NULL) { 
    printf("P3_FILE_ID=%s\n", pa->settings_file_id);
    free(warning);
  }
  */

  combined_retval_err = create_pr_append_str();
  if (NULL == combined_retval_err) exit(-2); /* Out of memory */

  if (pr_append_new_chunk_external(combined_retval_err, 
                                   retval->glob_err.data))
    exit(-2);

  if (pr_append_new_chunk_external(combined_retval_err, 
                                   retval->per_sequence_err.data)) 
    exit(-2);

  /* Check if there are errors, print and return */
  if (!pr_is_empty(combined_retval_err)) {
    print_boulder_error(pr_append_str_chars(combined_retval_err));
    destroy_pr_append_str(combined_retval_err);
    return;
  }
  destroy_pr_append_str(combined_retval_err);

  /* Get how many primers are in the array */
  num_fwd = retval->fwd.num_elem;
  num_rev = retval->rev.num_elem;
  num_int = retval->intl.num_elem;
  num_pair = retval->best_pairs.num_pairs;

  /* Prints out statistics about the primers */
  if (explain_flag) print_all_explain(pa, sa, retval, io_version);
    
  /* Print out the stop codon if a reading frame was specified */
  if (!PR_START_CODON_POS_IS_NULL(sa))
    printf("PRIMER_STOP_CODON_POSITION=%d\n", retval->stop_codon_pos);
    
  /* How often has the loop to be done? */
  if (retval->output_type == primer_list) {
    /* For Primer Lists: Figure out how many primers are in
     * the array that can be printed. If more than needed,
     * set it to the number requested. */

    /* Get how may primers should be printed */
    num_print = pa->num_return;
    /* Set how many primers will be printed */
    print_fwd = (num_print < num_fwd) ? num_print : num_fwd;
    print_rev = (num_print < num_rev) ? num_print : num_rev;
    print_int = (num_print < num_int) ? num_print : num_int;
    /* Get which list has to print most primers */
    loop_max = 0;
    if (loop_max < print_fwd) {
      loop_max = print_fwd;
    }
    if (loop_max < print_rev) {
      loop_max = print_rev;
    }
    if (loop_max < print_int) {
      loop_max = print_int;
    }
    /* Now the vars are there how often we have to go
     * through the loop and how many of each primer can
     * be printed. */
    num_pair = 0;
  } else {
    loop_max = num_pair;
    /* Set how many primers will be printed */
    print_fwd = num_pair;
    print_rev = num_pair;
    if (num_int != 0) {
      print_int = num_pair;
    }
  }

  if (io_version == 4) {
    printf("PRIMER_LEFT_NUM_RETURNED=%d\n", print_fwd);
    printf("PRIMER_RIGHT_NUM_RETURNED=%d\n",  print_rev);
    printf("PRIMER_%s_NUM_RETURNED=%d\n", int_oligo, print_int);
    printf("PRIMER_PAIR_NUM_RETURNED=%d\n", num_pair);
  }
  
  /* --------------------------------------- */
  /* Start of the loop printing all pairs or primers or oligos */
  for(i=0; i<loop_max; i++) {
    /* What needs to be printed */
    /* The conditions for primer lists */

    if (retval->output_type == primer_list) {
      /* Attach the selected primers to the pointers */
      fwd = &retval->fwd.oligo[i];
      rev = &retval->rev.oligo[i];
      intl = &retval->intl.oligo[i];
      /* Do fwd oligos have to be printed? */
      if ((pa->pick_left_primer) && (i < print_fwd)) {
        go_fwd = 1;
      } else {
        go_fwd = 0;
      }
      /* Do rev oligos have to be printed? */
      if ((pa->pick_right_primer) && (i < print_rev)) {
        go_rev = 1;
      } else {
        go_rev = 0;
      }
      /* Do int oligos have to be printed? */
      if ((pa->pick_internal_oligo) && (i < print_int)) {
        go_int = 1;
      } else {
        go_int = 0;
      }
    }  else {
      /* We will print primer pairs or pairs plus internal oligos */
      /* Get pointers to the primer_rec's that we will print */
      fwd  = retval->best_pairs.pairs[i].left;
      rev  = retval->best_pairs.pairs[i].right;
      intl = retval->best_pairs.pairs[i].intl;
      /* Pairs must have fwd and rev primers */
      go_fwd = 1;
      go_rev = 1;
      /* Do hyb oligos have to be printed? */
      if (pa->pick_internal_oligo == 1) {
        go_int = 1;
      } else {
        go_int = 0;
      }
    }
      
    /* Get the number for pimer counting in suffix[0] */
    sprintf(suffix, "_%d", i);
    
    /* Print out the Pair Penalties */
    if (retval->output_type == primer_pairs) {
      printf("PRIMER_PAIR%s_PENALTY=%f\n", suffix,
	     retval->best_pairs.pairs[i].pair_quality);
    }

    /* Print single primer penalty */
    if (go_fwd == 1)
      printf("PRIMER_LEFT%s_PENALTY=%f\n", suffix, fwd->quality);
    if (go_rev == 1)
      printf("PRIMER_RIGHT%s_PENALTY=%f\n", suffix, rev->quality);
    if (go_int == 1)
      printf("PRIMER_%s%s_PENALTY=%f\n", int_oligo, suffix, intl->quality);

    /* Print the oligo_problems */
    if (io_version == 4) {
      if (go_fwd == 1 && p3_ol_has_any_problem(fwd))
        printf("PRIMER_LEFT%s_PROBLEMS=%s\n", suffix, p3_get_ol_problem_string(fwd));
      if (go_rev == 1 && p3_ol_has_any_problem(rev))
        printf("PRIMER_RIGHT%s_PROBLEMS=%s\n", suffix, p3_get_ol_problem_string(rev));
      if (go_int == 1 && p3_ol_has_any_problem(intl))
        printf("PRIMER_%s%s_PROBLEMS=%s\n", int_oligo, suffix, p3_get_ol_problem_string(intl));
    }

    /* Print primer sequences. */
    if (go_fwd == 1)
      printf("PRIMER_LEFT%s_SEQUENCE=%s\n", suffix,
             pr_oligo_sequence(sa, fwd));
    if (go_rev == 1)
      printf("PRIMER_RIGHT%s_SEQUENCE=%s\n", suffix,
             pr_oligo_rev_c_sequence(sa, rev));
    if(go_int == 1)
      printf("PRIMER_%s%s_SEQUENCE=%s\n", int_oligo, suffix,
             pr_oligo_sequence(sa,intl));
        
    /* Print primer start and length */
    if (go_fwd == 1)
      printf("PRIMER_LEFT%s=%d,%d\n", suffix,
             fwd->start + incl_s + pa->first_base_index,
             fwd->length);
    if (go_rev == 1)
      printf("PRIMER_RIGHT%s=%d,%d\n", suffix,
             rev->start + incl_s + pa->first_base_index,
             rev->length);
    if (go_int == 1)
      printf("PRIMER_%s%s=%d,%d\n", int_oligo, suffix,
             intl->start + incl_s + pa->first_base_index,
             intl->length);

    /* Print primer Tm */
    if (go_fwd == 1)
      printf("PRIMER_LEFT%s_TM=%.3f\n", suffix, fwd->temp);
    if (go_rev == 1)
      printf("PRIMER_RIGHT%s_TM=%.3f\n", suffix, rev->temp);
    if (go_int == 1)
      printf("PRIMER_%s%s_TM=%.3f\n", int_oligo, suffix, intl->temp);

    /* Print primer GC content */
    if (go_fwd == 1)
      printf("PRIMER_LEFT%s_GC_PERCENT=%.3f\n", suffix, fwd->gc_content);
    if (go_rev == 1)
      printf("PRIMER_RIGHT%s_GC_PERCENT=%.3f\n", suffix, rev->gc_content);
    if (go_int == 1)
      printf("PRIMER_%s%s_GC_PERCENT=%.3f\n", int_oligo, suffix,
             intl->gc_content);

    /* Print primer self_any */
    if (go_fwd == 1 && pa->thermodynamic_oligo_alignment==0)
      printf("PRIMER_LEFT%s_SELF_ANY=%.2f\n", suffix,
             fwd->self_any);
    if (go_rev == 1 && pa->thermodynamic_oligo_alignment==0)
      printf("PRIMER_RIGHT%s_SELF_ANY=%.2f\n", suffix,
             rev->self_any);
    if (go_int == 1 && pa->thermodynamic_oligo_alignment==0)
      printf("PRIMER_%s%s_SELF_ANY=%.2f\n", int_oligo, suffix,
             intl->self_any);
     if (go_int == 1 && pa->thermodynamic_oligo_alignment==1)
       printf("PRIMER_%s%s_SELF_ANY_TH=%.2f\n", int_oligo, suffix,
	      intl->self_any);
    /* Print primer self_any thermodynamical approach */
     if (go_fwd == 1 && pa->thermodynamic_oligo_alignment==1)
       printf("PRIMER_LEFT%s_SELF_ANY_TH=%.2f\n", suffix,
	      fwd->self_any);
     if (go_rev == 1 && pa->thermodynamic_oligo_alignment==1)
       printf("PRIMER_RIGHT%s_SELF_ANY_TH=%.2f\n", suffix,
	      rev->self_any);
    /* Print primer self_end*/
    if (go_fwd == 1 && pa->thermodynamic_oligo_alignment==0)
      printf("PRIMER_LEFT%s_SELF_END=%.2f\n", suffix,
             fwd->self_end);
    if (go_rev == 1 && pa->thermodynamic_oligo_alignment==0)
      printf("PRIMER_RIGHT%s_SELF_END=%.2f\n", suffix,
             rev->self_end);
    if (go_int == 1 && pa->thermodynamic_oligo_alignment==0)
      printf("PRIMER_%s%s_SELF_END=%.2f\n", int_oligo, suffix,
             intl->self_end);
     if (go_int == 1 && pa->thermodynamic_oligo_alignment==1)
       printf("PRIMER_%s%s_SELF_END_TH=%.2f\n", int_oligo, suffix,
	      intl->self_end);
     /* Print primer self_end thermodynamical approach */
     if (go_fwd == 1 && pa->thermodynamic_oligo_alignment==1)
       printf("PRIMER_LEFT%s_SELF_END_TH=%.2f\n", suffix,
	      fwd->self_end);
     if (go_rev == 1 && pa->thermodynamic_oligo_alignment==1)
       printf("PRIMER_RIGHT%s_SELF_END_TH=%.2f\n", suffix,
	      rev->self_end);
     /* Print primer hairpin */
     if (go_fwd == 1 && pa->thermodynamic_oligo_alignment==1)
       printf("PRIMER_LEFT%s_HAIRPIN_TH=%.2f\n", suffix,
	      fwd->hairpin_th);
     if (go_rev == 1 && pa->thermodynamic_oligo_alignment==1)
       printf("PRIMER_RIGHT%s_HAIRPIN_TH=%.2f\n", suffix,
	      rev->hairpin_th);
     if (go_int == 1 && pa->thermodynamic_oligo_alignment==1)
       printf("PRIMER_%s%s_HAIRPIN_TH=%.2f\n", int_oligo, suffix,
	      intl->hairpin_th);
     /*Print out primer mispriming scores */
     if (seq_lib_num_seq(pa->p_args.repeat_lib) > 0) {
       if (go_fwd == 1)
	 printf("PRIMER_LEFT%s_LIBRARY_MISPRIMING=%.2f, %s\n", suffix,
		fwd->repeat_sim.score[fwd->repeat_sim.max],
		fwd->repeat_sim.name);
       if (go_rev == 1)
	 printf("PRIMER_RIGHT%s_LIBRARY_MISPRIMING=%.2f, %s\n", suffix,
		rev->repeat_sim.score[rev->repeat_sim.max],
		rev->repeat_sim.name);
       if (retval->output_type == primer_pairs)
	 printf("PRIMER_PAIR%s_LIBRARY_MISPRIMING=%.2f, %s\n", suffix,
		retval->best_pairs.pairs[i].repeat_sim,
		retval->best_pairs.pairs[i].rep_name);
     }
		
     /* Print out internal oligo mispriming scores */
     if (go_int == 1 && seq_lib_num_seq(pa->o_args.repeat_lib) > 0)
       printf("PRIMER_%s%s_LIBRARY_MISHYB=%.2f, %s\n", int_oligo, suffix,
	      intl->repeat_sim.score[intl->repeat_sim.max],
	      intl->repeat_sim.name);

    /* If a sequence quality was provided, print it*/
    if (NULL != sa->quality){
      if (go_fwd == 1)
        printf("PRIMER_LEFT%s_MIN_SEQ_QUALITY=%d\n", suffix,
               fwd->seq_quality);
      if (go_rev == 1) 
        printf("PRIMER_RIGHT%s_MIN_SEQ_QUALITY=%d\n", suffix,
               rev->seq_quality);
      if (go_int == 1 && (retval->output_type == primer_list)) 
        printf("PRIMER_%s%s_MIN_SEQ_QUALITY=%d\n", int_oligo, suffix,
               intl->seq_quality);
      /* Has to be here and in primer pairs for backward compatibility */
    }
        
    /* Print position penalty, this is for backward compatibility */
    if (!_PR_DEFAULT_POSITION_PENALTIES(pa)
        || !PR_START_CODON_POS_IS_NULL(sa)){
      printf("PRIMER_LEFT%s_POSITION_PENALTY=%f\n", suffix,
             fwd->position_penalty);
      printf("PRIMER_RIGHT%s_POSITION_PENALTY=%f\n", suffix,
             rev->position_penalty);
    }
        
    /* Print primer end stability */
    if (go_fwd == 1)
      printf("PRIMER_LEFT%s_END_STABILITY=%.4f\n",
             suffix, fwd->end_stability);
    if (go_rev == 1)
      printf("PRIMER_RIGHT%s_END_STABILITY=%.4f\n",
             suffix, rev->end_stability);

    /* Print primer template mispriming */
    if ( (pa->thermodynamic_template_alignment == 0) && (go_fwd == 1) && 
         (oligo_max_template_mispriming(fwd) != ALIGN_SCORE_UNDEF))
      printf("PRIMER_LEFT%s_TEMPLATE_MISPRIMING=%.4f\n", suffix,
             oligo_max_template_mispriming(fwd));
    if ( (pa->thermodynamic_template_alignment == 0) && (go_rev == 1) && 
         (oligo_max_template_mispriming(rev) != ALIGN_SCORE_UNDEF))
      printf("PRIMER_RIGHT%s_TEMPLATE_MISPRIMING=%.4f\n", suffix,
             oligo_max_template_mispriming(rev));

     /* Print primer template mispriming, thermodynamical approach*/
     if ( (pa->thermodynamic_template_alignment == 1) && (go_fwd == 1) &&
	  (oligo_max_template_mispriming_thermod(fwd) != ALIGN_SCORE_UNDEF)) {
       printf("PRIMER_LEFT%s_TEMPLATE_MISPRIMING_TH=%.4f\n", suffix,
	      oligo_max_template_mispriming_thermod(fwd));
     }

     if ( (pa->thermodynamic_template_alignment == 1) && (go_rev == 1) &&
	  (oligo_max_template_mispriming_thermod(rev) != ALIGN_SCORE_UNDEF)) {
       printf("PRIMER_RIGHT%s_TEMPLATE_MISPRIMING_TH=%.4f\n", suffix,
	      oligo_max_template_mispriming_thermod(rev));
#if 0
       printf("DEBUG_PRIMER_RIGHT%s_TEMPLATE_MISPRIMING_TOP_TH=%.4f\n", suffix,
	      rev->template_mispriming);
       printf("DEBUG_PRIMER_RIGHT%s_TEMPLATE_MISPRIMING_R_TH=%.4f\n", suffix,
	      rev->template_mispriming_r);
#endif
     }
     /************************************************************************************/
     /* Print the pair parameters*/
    if (retval->output_type == primer_pairs) {
      if (go_int == 1 && NULL != sa->quality) /* FIX ME - Uptate the tests */
        printf("PRIMER_%s%s_MIN_SEQ_QUALITY=%d\n", int_oligo,
               suffix, intl->seq_quality);
      /* Print pair comp_any */
       if(pa->thermodynamic_oligo_alignment==0)
	 printf("PRIMER_PAIR%s_COMPL_ANY=%.2f\n", suffix,
		retval->best_pairs.pairs[i].compl_any);
       if(pa->thermodynamic_oligo_alignment==1)
	 printf("PRIMER_PAIR%s_COMPL_ANY_TH=%.2f\n", suffix,
		retval->best_pairs.pairs[i].compl_any);
       /* Print pair comp_end */
       if(pa->thermodynamic_oligo_alignment==0)
	 printf("PRIMER_PAIR%s_COMPL_END=%.2f\n", suffix,
		retval->best_pairs.pairs[i].compl_end);
       if(pa->thermodynamic_oligo_alignment==1)
	 printf("PRIMER_PAIR%s_COMPL_END_TH=%.2f\n", suffix,
		retval->best_pairs.pairs[i].compl_end);
       /* Print product size */
       printf("PRIMER_PAIR%s_PRODUCT_SIZE=%d\n", suffix,
	      retval->best_pairs.pairs[i].product_size);
       /* Print the product Tm if a Tm range is defined */
       if (pa->product_max_tm != PR_DEFAULT_PRODUCT_MAX_TM ||
	   pa->product_min_tm != PR_DEFAULT_PRODUCT_MIN_TM) {
	 printf("PRIMER_PAIR%s_PRODUCT_TM=%.4f\n", suffix,
		retval->best_pairs.pairs[i].product_tm);
	 
	 printf("PRIMER_PAIR%s_PRODUCT_TM_OLIGO_TM_DIFF=%.4f\n", suffix,
		retval->best_pairs.pairs[i].product_tm_oligo_tm_diff);
	 
	 printf("PRIMER_PAIR%s_T_OPT_A=%.4f\n", suffix,
		retval->best_pairs.pairs[i].t_opt_a);
       }
      
      /* Print the primer pair template mispriming */
      if ((pa->thermodynamic_template_alignment == 0) && (retval->best_pairs.pairs[i].template_mispriming != ALIGN_SCORE_UNDEF))
        printf("PRIMER_PAIR%s_TEMPLATE_MISPRIMING=%.2f\n", suffix,
               retval->best_pairs.pairs[i].template_mispriming);
       /* Print the primer pair template mispriming. Thermodynamic approach.  */
       if ((pa->thermodynamic_template_alignment == 1) && (retval->best_pairs.pairs[i].template_mispriming != ALIGN_SCORE_UNDEF))
	 printf("PRIMER_PAIR%s_TEMPLATE_MISPRIMING_TH=%.2f\n", suffix,
		retval->best_pairs.pairs[i].template_mispriming);

    } /* End of print parameters of primer pairs */
        
  } /* End of the big loop printing all data */
    
    /* End the print with newline and flush all buffers */
  printf("=\n");
  if (fflush(stdout) == EOF) {
    perror("fflush(stdout) failed");
    exit(-1);
  }
}

void
print_boulder_error(const char *err) 
{
  printf("PRIMER_ERROR=%s\n=\n", err);
  if (fflush(stdout) == EOF) {
    perror("fflush(stdout) failed");
    exit(-1);
  }
}

void
print_boulder_warning(const char *err) 
{
  printf("PRIMER_WARNING=%s\n", err);
  if (fflush(stdout) == EOF) {
    perror("fflush(stdout) failed");
    exit(-1);
  }
}

static void
print_all_explain(const p3_global_settings *pa,
                  const seq_args *sa, 
                  const p3retval *retval,
                  int io_version)
{
  if (pa->pick_left_primer == 1
      && !(pa->pick_anyway && sa->left_input))
    printf("PRIMER_LEFT_EXPLAIN=%s\n",
           p3_get_oligo_array_explain_string(p3_get_rv_fwd(retval)));

  if (pa->pick_right_primer == 1
      && !(pa->pick_anyway && sa->right_input))
    printf("PRIMER_RIGHT_EXPLAIN=%s\n",
           p3_get_oligo_array_explain_string(p3_get_rv_rev(retval)));

  if ( pa->pick_internal_oligo == 1
      && !(pa->pick_anyway && sa->internal_input)) 
    printf("PRIMER_INTERNAL_EXPLAIN=%s\n",
           p3_get_oligo_array_explain_string(p3_get_rv_intl(retval)));

  if (pa->pick_right_primer == 1 
      && pa->pick_left_primer == 1) {
    printf("PRIMER_PAIR_EXPLAIN=%s\n", 
           p3_get_pair_array_explain_string(p3_get_rv_best_pairs(retval)));
  }
}
