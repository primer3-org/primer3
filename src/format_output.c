/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008,2009,2010
Whitehead Institute for Biomedical Research, Steve Rozen
(http://purl.com/STEVEROZEN), Andreas Untergasser and Helen Skaletsky
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
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANYT
HEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>
#include <string.h>
#include "format_output.h"
#include "libprimer3.h"

#include "p3_exit_codes.h"
#define FORWARD 1
#define INTERNAL 0
#define REVERSE -1

static const char *pr_program_name = "Program name is probably primer3_core";

static void format_pairs(FILE *f, 
                         const p3_global_settings *pa,
                         const seq_args *sa, 
                         const p3retval *retval, 
                         const pair_array_t *best_pairs,
                         const char *pr_release,
                         const pr_append_str *combined_retval_err,
                         int explain_flag);

static void format_oligos(FILE *, 
                          const p3_global_settings *, 
                          const seq_args *, 
                          const p3retval *retval, 
                          const char*,
                          const pr_append_str *combined_retval_err,
                          int explain_flag);

static int lib_sim_specified(const p3_global_settings *);
static void print_explain(FILE *, const p3_global_settings *, const seq_args *,
                          const p3retval *retval, int, const char *);
static void print_pair_info(FILE *, const primer_pair *, const seq_args *,
                            const p3_global_settings *);
static void print_oligo(FILE *, const char *, const seq_args *,
                        const primer_rec *, int, const p3_global_settings *, 
                        const seq_lib*, int);
static void print_oligo_header(FILE *, const char *, const int, const int);
static void print_pair_array(FILE *, const char*, int,
                             const interval_array_t, 
                             const p3_global_settings*, const seq_args*);
static void print_2_pair_array(FILE *, const char*, int,
                               const interval_array_t, 
                               const interval_array_t, 
                               const p3_global_settings*, const seq_args*);
static void print_rest(FILE *, const p3_global_settings *, 
                       const seq_args *,  const pair_array_t *);
static int  print_seq(FILE *, const p3_global_settings *, const seq_args *, 
                      const p3retval *retval, const primer_rec *left,
                      const primer_rec *intl, const primer_rec *right,
                      const pair_array_t *, int);
static void print_seq_lines(FILE *, const char *s, const char *n, int, int,
                            int, const p3_global_settings *);
static void print_stat_line(FILE *, const char *, oligo_stats s, int, int, int);
static void print_summary(FILE *, const p3_global_settings *, 
                          const seq_args *, const pair_array_t *, int);

void
print_format_output(FILE *f,
                    const int *io_version,
                    const p3_global_settings *pa,
                    const seq_args *sa,
                    const p3retval *retval,
                    const char *pr_release,
                    int   explain_flag)
{  
  
  /* A place to put a string containing all error messages */
  pr_append_str *combined_retval_err = NULL;

  combined_retval_err = create_pr_append_str();
  if (NULL == combined_retval_err) {
    exit(P3_EXIT_OUT_OF_MEMORY); /* Out of memory */
  }

  if (pr_append_new_chunk_external(combined_retval_err, 
                                   retval->glob_err.data))
    exit(P3_EXIT_OUT_OF_MEMORY);

  if (pr_append_new_chunk_external(combined_retval_err, 
                                   retval->per_sequence_err.data)) 
    exit(P3_EXIT_OUT_OF_MEMORY);

  /* Print as primer pairs */
  if (retval->output_type == primer_pairs) {
    format_pairs(f, pa, sa, retval, &retval->best_pairs, 
                 pr_release, combined_retval_err, explain_flag);
    
    /* Print as primer list */
  } else {
    format_oligos(f, pa, sa, retval, pr_release,
                  combined_retval_err, explain_flag);
  }

  destroy_pr_append_str(combined_retval_err);
}

static void
format_pairs(FILE *f,
             const p3_global_settings *pa,
             const seq_args *sa,
             const p3retval *retval,
             const pair_array_t *best_pairs,
             const char *pr_release,
             const pr_append_str *combined_retval_err,
             int explain_flag)
{
  char *warning;
  int print_lib_sim = lib_sim_specified(pa);

  PR_ASSERT(NULL != f);
  PR_ASSERT(NULL != pa);
  PR_ASSERT(NULL != sa);
  
  /* If there are errors, print them and return */
  if (!pr_is_empty(combined_retval_err)) {
    format_error(f, sa->sequence_name, 
                 pr_append_str_chars(combined_retval_err));
    return;
    }
  
  /* Print the sequence name if it is provided */
  if (NULL != sa->sequence_name)
    fprintf(f, "PRIMER PICKING RESULTS FOR %s\n\n", sa->sequence_name);
  
  /* Print if a mispriming libraby was used and which one */
  if (pa->p_args.repeat_lib != NULL)
    fprintf(f, "Using mispriming library %s\n",
            pa->p_args.repeat_lib->repeat_file);
  else
    fprintf(f, "No mispriming library specified\n");

  /* Print if a mispriming libraby for the internal oligo 
   * was used and which one */
  if ( pa->pick_internal_oligo == 1 ) {
    if (pa->o_args.repeat_lib != NULL)
      fprintf(f, "Using internal oligo mishyb library %s\n",
              pa->o_args.repeat_lib->repeat_file);
    else
      fprintf(f, "No internal oligo mishyb library specified\n");
  }

  /* Does the sequence start at position 0 or 1 ? */
  fprintf(f, "Using %d-based sequence positions\n",
          pa->first_base_index);
  
  /* Complain if no primers are in the array */
  if (best_pairs->num_pairs == 0) fprintf(f, "NO PRIMERS FOUND\n\n");
  
  /* Print out the warings */
  if ((warning = p3_get_rv_and_gs_warnings(retval, pa)) != NULL) {
    fprintf(f, "WARNING: %s\n\n", warning);
    free(warning);
  }
  
  /* Print the results for the best pair */
  print_summary(f, pa, sa, best_pairs, 0);
  fprintf(f, "\n");

  /* Print nicely out the sequence with the best pair */
  if (print_seq(f, pa, sa, retval, NULL, NULL, NULL, best_pairs, 0)) {
    exit(P3_EXIT_OUT_OF_MEMORY); /* ENOMEM */
  }
  
  /* Print out the alternative pairs */
  if (best_pairs->num_pairs > 1 ) print_rest(f, pa, sa, best_pairs);
  
  /* Print the primer picking statistics */
  if (explain_flag)
    print_explain(f, pa, sa, retval, print_lib_sim, pr_release);
  
  /* Flush the buffers and return */
  fprintf(f, "\n\n");
  if (fflush(f) == EOF) {
    perror("fflush(f) failed");
    exit(P3_EXIT_OUTPUT_ERROR);
  }

}

/* Prints out the results of a primer pair */
static void
print_summary(FILE *f, const p3_global_settings *pa,
              const seq_args *sa, const pair_array_t *best_pairs,
              int num)
{
    int seq_len = strlen(sa->sequence);
    int print_lib_sim = lib_sim_specified(pa);
    primer_pair *p;
    p = best_pairs->pairs + num;
    if (best_pairs->num_pairs > 0) {
        /* 
         * If the following format changes, also change the format in
         * print_oligo.
         */
        print_oligo_header(f, "OLIGO", print_lib_sim, pa->thermodynamic_oligo_alignment);
        print_oligo(f, "LEFT PRIMER", sa, p->left, FORWARD, pa,
                    pa->p_args.repeat_lib,
                    print_lib_sim);
        print_oligo(f, "RIGHT PRIMER", sa, p->right, REVERSE, pa,
                    pa->p_args.repeat_lib,
                    print_lib_sim);
        if ( pa->pick_internal_oligo == 1 )
            print_oligo(f, "INTERNAL OLIGO", sa, p->intl, INTERNAL,
                        pa, pa->o_args.repeat_lib,
                        print_lib_sim);
    }
    fprintf(f, "SEQUENCE SIZE: %d\n", seq_len);
    fprintf(f, "INCLUDED REGION SIZE: %d\n\n", sa->incl_l);

    if (best_pairs->num_pairs > 0) print_pair_info(f, p, sa, pa);
    print_pair_array(f, "TARGETS", sa->tar2.count, sa->tar2.pairs, pa, sa);
    print_pair_array(f, "EXCLUDED REGIONS", sa->excl2.count, sa->excl2.pairs, pa, sa);
    print_pair_array(f, "INTERNAL OLIGO EXCLUDED REGIONS",
                     sa->excl_internal2.count, sa->excl_internal2.pairs, pa, sa);
    print_2_pair_array(f, "PAIR_OK_REGIONS", sa->ok_regions.count, sa->ok_regions.left_pairs, 
                     sa->ok_regions.right_pairs, pa, sa);
}

/* Print column headers for lines printed by print_oligo(). */
static void
print_oligo_header(FILE *f, const char *s, const int print_lib_sim,
                   const int thermodynamic_oligo_alignment)
{
   if(thermodynamic_oligo_alignment==0)
     fprintf(f,
             "%-16s start  len      tm     gc%%   any    3' %sseq\n",
             s, print_lib_sim ? "  rep " : "");
   else 
     fprintf(f,
             "%-16s start  len      tm     gc%%  any_th  3'_th hairpin %sseq\n",
             s, print_lib_sim ? "  rep " : "");
}

/* Print the line with the parameters */
static void
print_oligo(FILE *f,
            const char *title,
            const seq_args *sa,
            const primer_rec *o,
            int dir,
            const p3_global_settings *pa,
            const seq_lib *seqlib,
            int print_lib_sim)
{
    const char *format1;
    char *seq;
    /* The primers with their overhangs, as in the Boulder-IO output */
    if (dir == FORWARD) {
      seq = pr_oligo_overhang_sequence(sa, o);
    }
    if (dir == REVERSE) {
      seq = pr_oligo_rev_c_overhang_sequence(sa, o);
    }
    if (dir == INTERNAL) {
      if (o->rev_comp > 0) {
        seq = pr_internal_rev_c_sequence(sa, o);
      } else {
        seq = pr_oligo_sequence(sa, o);
      }
    }
    if(pa->thermodynamic_oligo_alignment==1) {
      format1 = "%-16s %5d %4d %7.2f %7.2f   %5.2f  %5.2f ";
      fprintf(f, format1,
              title, o->start + sa->incl_s + pa->first_base_index,
              o->length, o->temp, o->gc_content, o->self_any,
              o->self_end);
      fprintf(f, "  %5.2f ",  o->hairpin_th);
    } else {
      format1  = "%-16s %5d %4d %7.2f %7.2f %5.2f %5.2f ";
      fprintf(f, format1,
              title, o->start + sa->incl_s + pa->first_base_index,
              o->length, o->temp, o->gc_content, o->self_any,
              o->self_end);
    }
    if (print_lib_sim) {
        if (seqlib != NULL && o->repeat_sim.score != NULL)
            fprintf(f, "%5.2f ",  o->repeat_sim.score[o->repeat_sim.max]);
        else 
            fprintf(f, "%5s ", "");
    }
    fprintf(f, "%s\n", seq);
    if (PR_DEFAULT_INSIDE_PENALTY != pa->inside_penalty
        || PR_DEFAULT_OUTSIDE_PENALTY != pa->outside_penalty)
      fprintf(f, "POSITION PENALTY, QUALITY: %f, %f\n",
              o->position_penalty, o->quality);
}

static void
print_pair_array(FILE *f, const char* title, int num,
                 const interval_array_t array,
                 const p3_global_settings *pa,
                 const seq_args *sa)
{
    int j;
    if (num > 0) {
        fprintf(f, "%s (start, len)*:", title);
        for (j = 0; j < num; j++)
            fprintf(f, " %d,%d", 
                    array[j][0] + pa->first_base_index + sa->incl_s,
                    array[j][1]);
        fprintf(f, "\n");
    }
}

static void
print_2_pair_array(FILE *f, const char* title, int num,
                   const interval_array_t left_array,
                   const interval_array_t right_array,
                   const p3_global_settings *pa,
                   const seq_args *sa)
{
    int j;
    if (num > 0) {
        fprintf(f, "%s (left_start, left_len, right_start, right_len)*:", title);
        for (j = 0; j < num; j++) {
          if ((left_array[j][0] == -1) && (left_array[j][1] == -1))
            fprintf(f, " ,,");
          else fprintf(f, " %d,%d,", 
                       left_array[j][0] + pa->first_base_index + sa->incl_s,
                       left_array[j][1]);
          if ((right_array[j][0] == -1) && (right_array[j][1] == -1))
            fprintf(f, ",");
          else fprintf(f, "%d,%d", 
                       right_array[j][0] + pa->first_base_index + sa->incl_s,
                       right_array[j][1]);
        }
        fprintf(f, "\n");
    }
}

#define VECTOR           (1<<0)
#define LEFT_OLIGO       (1<<1)
#define RIGHT_OLIGO      (1<<2)
#define INTL_OLIGO       (1<<3)
#define TARGET           (1<<4)
#define EXCL_REGION      (1<<5)
#define INTL_EXCL_REGION (1<<6)

/* Prints out the asci picture of the sequence */
/* Return 1 on ENOMEM. Otherwise return 0. retval must not be NULL. */
/* left, intl and right are the oligos to draw when there are no primer
   pairs and the task is not pick_sequencing_primers; each may be NULL. */
static int
print_seq(FILE *f,
    const p3_global_settings *pa,
    const seq_args *sa,
    const p3retval *retval,
    const primer_rec *left,
    const primer_rec *intl,
    const primer_rec *right,
    const pair_array_t *best_pairs,
    int num)  /* The number of primer pair to print. */
{
    primer_rec *h2 = NULL;
    primer_rec *h3 = NULL;
    int len, i, j, start;
    int something_found = 0, vector_found = 0;
    int *notes;
    char *notestr;
    primer_pair *p;
    p = NULL;
    if(retval->output_type == primer_pairs) {
      p = best_pairs->pairs + num;
    }
    len = strlen(sa->sequence);
    if (!(notes = (int*) malloc(sizeof(*notes) * len))) return 1;
    memset(notes, 0, sizeof(*notes) * len);
    if (!(notestr = (char*) malloc(len + 1))) { free(notes); return 1; }
    memset(notestr, ' ', len);
    notestr[len] = '\0';

    for (i = 0; i < len; i++) {
        if (i < sa->incl_s || i >= sa->incl_s + sa->incl_l)
            notes[i] |= VECTOR;

        if (retval->output_type == primer_pairs && 
            best_pairs->num_pairs > 0) {
            if (i >= p->left->start + sa->incl_s
                && i < p->left->start + p->left->length + sa->incl_s)
                notes[i] |= LEFT_OLIGO;
            if (i >= p->right->start - p->right->length + 1 + sa->incl_s
                && i <= p->right->start + sa->incl_s)
                notes[i] |= RIGHT_OLIGO;
            if ( pa->pick_internal_oligo == 1
                && i >= p->intl->start + sa->incl_s 
                && i < p->intl->start + p->intl->length + sa->incl_s)
                notes[i] |= INTL_OLIGO;
        }
        else if (pa->primer_task != pick_sequencing_primers) {
            if (left != NULL
                && i >= left->start + sa->incl_s
                && i < left->start + left->length + sa->incl_s) {
                notes[i] |= LEFT_OLIGO;
            }
            if (right != NULL
                && i >= right->start - right->length + 1 + sa->incl_s
                && i <= right->start + sa->incl_s) {
                notes[i] |= RIGHT_OLIGO;
            }
            if (intl != NULL
                && i >= intl->start + sa->incl_s
                && i < intl->start + intl->length + sa->incl_s) {
                notes[i] |= INTL_OLIGO;
            }
        } else {
                if (pa->pick_right_primer && retval->rev.num_elem > 0){
                  h2 = retval->rev.oligo;
              for (j = 0; j < pa->num_return; j++) {
                    if(j > retval->rev.num_elem -1) break;
                    h3 = h2 + j;
                if(i >= h3->start - h3->length + 1 + sa->incl_s
                        && i <= h3->start + sa->incl_s)
                   notes[i] |= RIGHT_OLIGO;
                  }
            }           
                if (pa->pick_left_primer && retval->fwd.num_elem > 0){
                  h2 = retval->fwd.oligo;
              for (j = 0; j < pa->num_return; j++) {
                    if(j > retval->fwd.num_elem -1) break;
                    h3 = h2 + j;
                if(i < h3->start + h3->length + sa->incl_s &&
                   i >= h3->start + sa->incl_s) {
                   notes[i] |= LEFT_OLIGO;
                }
                  }
            }               
        }

        for (j = 0; j < sa->tar2.count; j++) {
            start = sa->tar2.pairs[j][0] + sa->incl_s;
            if (i >= start && i < start + sa->tar2.pairs[j][1])
                notes[i] |= TARGET;
        }
        for (j = 0; j < sa->excl2.count; j++) {
            start = sa->excl2.pairs[j][0] + sa->incl_s;
            if (i >= start && i < start + sa->excl2.pairs[j][1])
                notes[i] |= EXCL_REGION;
        }
        for (j = 0; j < sa->excl_internal2.count; j++) {
            start = sa->excl_internal2.pairs[j][0] + sa->incl_s;
            if (i >= start && i < start + sa->excl_internal2.pairs[j][1])
                notes[i] |= INTL_EXCL_REGION;
        }
    }
    for (i = 0; i < len; i++) {
        if (notes[i] & VECTOR) {
            vector_found = 1;
            notestr[i] = '.';
        }
        else if (notes[i] & EXCL_REGION)
            notestr[i] = 'X';
        else if (notes[i] & INTL_EXCL_REGION)
            notestr[i] = 'x';
        else if ((pa->primer_task == pick_sequencing_primers)
                && (notes[i] & LEFT_OLIGO)
                && (notes[i] & RIGHT_OLIGO))
            notestr[i] = '^';
        else if ((pa->primer_task == pick_sequencing_primers)
                && (notes[i] & LEFT_OLIGO))
            notestr[i] = '>';
        else if ((pa->primer_task == pick_sequencing_primers)
                && (notes[i] & RIGHT_OLIGO))
            notestr[i] = '<';
        else if ((notes[i] & TARGET) && (notes[i] & LEFT_OLIGO))
            notestr[i] = ')';
        else if ((notes[i] & TARGET) && (notes[i] & RIGHT_OLIGO))
            notestr[i] = '(';
        else if (notes[i] & TARGET)
            notestr[i] = '*';
        else if (notes[i] & LEFT_OLIGO)
            notestr[i] = '>';
        else if (notes[i] & RIGHT_OLIGO)
            notestr[i] = '<';
        else if (notes[i] & INTL_OLIGO)
            notestr[i] = '^';

        if (notes[i] != 0) something_found = 1;
    }

    print_seq_lines(f, sa->sequence, notestr, len, 60, something_found, pa);

    if (something_found)
        fprintf(f, "KEYS (in order of precedence):\n");

    if (vector_found)
        fprintf(f, "...... vector sequence\n");

    if (sa->excl2.count > 0)
        fprintf(f, "XXXXXX excluded region\n");

    if (pa->pick_internal_oligo ==1 && sa->excl_internal2.count > 0)
        fprintf(f, "xxxxxx excluded region for internal oligo\n");

    if (sa->tar2.count > 0)
        fprintf(f, "****** target\n");

    if (retval->output_type == primer_pairs &&
        best_pairs->num_pairs > 0) {
           fprintf(f, ">>>>>> left primer\n");
           fprintf(f, "<<<<<< right primer\n");
           if ( pa->pick_internal_oligo == 1 )
              fprintf(f, "^^^^^^ internal oligo\n");
    } else if (pa->primer_task == pick_sequencing_primers) {
        fprintf(f, ">>>>>> left primer\n");
        fprintf(f, "<<<<<< right primer\n");
        fprintf(f, "^^^^^^ left primer / right primer overlap\n");
            
    }
    else {
        if (left != NULL) {
            fprintf(f, ">>>>>> left primer\n");
        }
        if (right != NULL) {
            fprintf(f, "<<<<<< right primer\n");
        }
        if (intl != NULL) {
            fprintf(f, "^^^^^^ internal oligo\n");
        }
    }

    if (something_found) fputc('\n', f);
    free(notes);
    free(notestr);
    return 0;
}

static void
print_seq_lines(FILE *f, const char *s, const char *n, int seq_size, 
                int line_size, int something_found,
                const p3_global_settings *pa)
{
    int i = 0;
    while (seq_size > line_size) {
        fprintf(f, "%5d ", i + pa->first_base_index);
        fwrite(s, sizeof(*s), line_size, f);
        fputc('\n', f);
        if (something_found) {
            fprintf(f, "      ");
            fwrite(n, sizeof(*n), line_size, f);
            fprintf(f, "\n\n");
        }
        seq_size -= line_size;
        s += line_size;
        n += line_size;
        i += line_size;
    }
    if (something_found)
        fprintf(f, "%5d %s\n      %s\n\n", i + pa->first_base_index, s, n);
    else
        fprintf(f, "%5d %s\n\n", i + pa->first_base_index, s);
}

static void
print_pair_info(FILE *f, const primer_pair *p, const seq_args *sa,
                const p3_global_settings *pa)
{
  /* The product includes the overhangs, as in the Boulder-IO output */
  int product_size = p->product_size;
  if (NULL != sa->overhang_left) {
    product_size += strlen(sa->overhang_left);
  }
  if (NULL != sa->overhang_right) {
    product_size += strlen(sa->overhang_right);
  }
  fprintf(f, "PRODUCT SIZE: %d, ", product_size);
   if(pa->thermodynamic_oligo_alignment==0)
     fprintf(f, "PAIR ANY COMPL: %.2f, PAIR 3' COMPL: %.2f\n",
             p->compl_any, p->compl_end);
   else 
     fprintf(f, "PAIR ANY_TH COMPL: %.2f, PAIR 3'_TH COMPL: %.2f\n",
             p->compl_any, p->compl_end);
  if (pa->product_max_tm != PR_DEFAULT_PRODUCT_MAX_TM
      || pa->product_min_tm != PR_DEFAULT_PRODUCT_MIN_TM) {
    fprintf(f, "PRODUCT Tm: %.4f, ", p->product_tm);
    fprintf(f, "PRODUCT Tm - min(OLIGO Tm): %.4f\n",
            p->product_tm_oligo_tm_diff);
  }
}

static void
print_rest(FILE *f, const p3_global_settings *pa,
           const seq_args *sa,
           const pair_array_t *best_pairs)
{
    int i;
    int print_lib_sim = lib_sim_specified(pa);

    fprintf(f, "ADDITIONAL OLIGOS\n");
    fprintf(f, "   "); print_oligo_header(f, "", print_lib_sim, pa->thermodynamic_oligo_alignment);
    for (i = 1; i < best_pairs->num_pairs; i++) {
        fprintf(f, "\n%2d ", i);
        print_oligo(f, "LEFT PRIMER", sa, best_pairs->pairs[i].left, FORWARD,
                    pa, pa->p_args.repeat_lib, print_lib_sim);
        fprintf(f, "   ");
        print_oligo(f, "RIGHT PRIMER", sa, best_pairs->pairs[i].right, REVERSE,
                    pa, pa->p_args.repeat_lib, print_lib_sim);
        if ( pa->pick_internal_oligo == 1) {
            fprintf(f, "   ");
            print_oligo(f, "INTERNAL OLIGO", sa, best_pairs->pairs[i].intl,
                        INTERNAL, pa, pa->o_args.repeat_lib, print_lib_sim);
        }
        if (best_pairs->pairs[i].product_size > 0) {
            fprintf(f, "   ");
            print_pair_info(f, &best_pairs->pairs[i], sa, pa);
        }
    }
}

/* Print out the statistics of primer picking */
/* This function does _not_ print out the no_orf statistic. */
static void
print_explain(FILE *f,
              const p3_global_settings *pa,
              const seq_args *sa,
              const p3retval *retval,
              int print_lib_sim,
              const char *pr_release)
{
   
   
  const char *format;  /* Format string for the table headers */
  /* Template masking uses the lowercase masking check */
  int lowercase_masking = pa->lowercase_masking || pa->mask_template;
  if(pa->thermodynamic_oligo_alignment==0)
    format    = "%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s";
  else format = "%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%7s%6s";

  fprintf(f, "\nStatistics\n");
   
  if (!pa->pick_anyway || !(
                            (pa->pick_left_primer == 1 && pa->pick_internal_oligo == 0
                             && pa->pick_right_primer == 1
                             && sa->left_input && sa->right_input)
                            || (pa->pick_left_primer == 1 && pa->pick_internal_oligo == 1
                                && pa->pick_right_primer == 1
                                && sa->left_input && sa->right_input && sa->internal_input)
                            || (pa->pick_left_primer == 1 && pa->pick_internal_oligo == 0
                                && pa->pick_right_primer == 0
                                && sa->left_input)
                            || (pa->pick_left_primer == 0 && pa->pick_internal_oligo == 0
                                && pa->pick_right_primer == 1
                                && sa->right_input)
                            || (pa->pick_left_primer == 0 && pa->pick_internal_oligo == 1
                                && pa->pick_right_primer == 0
                                && sa->internal_input))) {
      /* print first row of header */
      if(pa->thermodynamic_oligo_alignment==0)
        fprintf(f, format, "", "con", "too",  "in",  "in",  "not", "",    "no",
                "tm",  "tm",  "  high", "high");
      else 
        fprintf(f, format, "", "con", "too",  "in",  "in",  "not", "",    "no",
                "tm",  "tm",  "   high", "high");
      if(pa->thermodynamic_oligo_alignment==1) 
        fprintf(f, "%6s", "high");
      if(print_lib_sim) 
        fprintf(f, "%6s", "high");
      fprintf(f, "%6s%6s", "", "high");
      if(lowercase_masking) 
        fprintf(f, "%6s", "lower");
      fprintf(f, "%6s\n", "");
      
      /* print second row of header */
      if(pa->thermodynamic_oligo_alignment==0) {
        fprintf(f, format, "", "sid", "many", "tar", "excl", "ok", "bad","GC",
                "too", "too", "any",  "3'");
      } else {
        fprintf(f, format, "", "sid", "many", "tar", "excl", "ok", "bad","GC",
                "too", "too", " any_th",  "3'_th");
        fprintf(f, "%6s", "hair-");
      }
           
      if(print_lib_sim) 
        fprintf(f, "%6s", "lib");
      fprintf(f, "%6s%6s", "poly", "end");
      if(lowercase_masking) 
        fprintf(f, "%6s", " case");
      fprintf(f, "%6s\n", "");
      
      /* print third row of header */
      if(pa->thermodynamic_oligo_alignment==0) 
        fprintf(f, format, "", "ered","Ns",   "get", "reg",  "reg", "GC%", "clamp",
                "low", "high"," compl", "compl");
      else
        fprintf(f, format, "", "ered","Ns",   "get", "reg",  "reg", "GC%", "clamp",
                "low", "high","  compl", "compl");
      if(pa->thermodynamic_oligo_alignment==1) 
        fprintf(f, "%6s", " pin");
      if(print_lib_sim) 
        fprintf(f, "%6s", "sim");
      fprintf(f, "%6s%6s", "X", "stab");
      if(lowercase_masking) 
        fprintf(f, "%6s", " end");
      if(lowercase_masking)
        fprintf(f, "%6s\n", "ok  ");
      else
        fprintf(f, "%6s\n", "ok");
  }
  if (pa->pick_left_primer == 1
      && !(pa->pick_anyway && sa->left_input))
    print_stat_line(f, "Left", retval->fwd.expl, 
                    print_lib_sim, lowercase_masking, pa->thermodynamic_oligo_alignment);
  
  if (pa->pick_right_primer == 1
      && !(pa->pick_anyway && sa->right_input))
    print_stat_line(f, "Right", retval->rev.expl,
                    print_lib_sim, lowercase_masking, pa->thermodynamic_oligo_alignment);
  
  if (pa->pick_internal_oligo == 1
      && !(pa->pick_anyway && sa->internal_input))
    print_stat_line(f, "Intl", retval->intl.expl, 
                    print_lib_sim, lowercase_masking, pa->thermodynamic_oligo_alignment);
  
  if (pa->pick_left_primer == 1 && pa->pick_right_primer == 1) {
    fprintf(f, "Pair Stats:\n%s\n",
            p3_get_pair_array_explain_string(p3_get_rv_best_pairs(retval)));
  }
  fprintf(f, "%s\n", pr_release);
}

static void
print_stat_line(FILE *f, const char *t, oligo_stats s,
                int print_lib_sim, int lowercase_masking,
                int thermodynamic_oligo_alignment)
{
   const char *format = "%-6s%6d%6d%6d%6d%6d%6d%6d%6d%6d";
   if (!strcmp(t, "Left"))
     fprintf(f,format, t, s.considered, s.ns, s.target, s.excluded,
             s.not_in_any_left_ok_region, s.gc, s.gc_clamp, s.temp_min, s.temp_max);
   else if (!strcmp(t, "Right"))
          fprintf(f,format, t, s.considered, s.ns, s.target, s.excluded,
                  s.not_in_any_right_ok_region, s.gc, s.gc_clamp, s.temp_min, s.temp_max);
   else fprintf(f,format, t, s.considered, s.ns, s.target, s.excluded,
                0, s.gc, s.gc_clamp, s.temp_min, s.temp_max);
   if(thermodynamic_oligo_alignment) {
      fprintf(f, " %6d%6d%6d", s.compl_any, s.compl_end, s.hairpin_th);
   } else {
      fprintf(f, "%6d%6d", s.compl_any, s.compl_end);
   }
   if (print_lib_sim) fprintf(f, "%6d",s.repeat_score);
   fprintf(f, "%6d%6d",s.poly_x, s.stability);
   if (lowercase_masking) fprintf(f, "%6d",s.gmasked);
   fprintf(f, "%6d\n", s.ok);
}


/* 
 * Return true iff a check for library similarity has been specified for
 * either the primer pair or the internal oligo.
 */
static int
lib_sim_specified(const p3_global_settings *pa) {
  return (pa->p_args.repeat_lib || pa->o_args.repeat_lib);
}

void
format_error(FILE *f, const char* seq_name, const char *err)
{
  if (NULL != seq_name)
    fprintf(f, "PRIMER PICKING RESULTS FOR %s\n\n", seq_name);
  if (err != NULL) 
    fprintf(f, "INPUT PROBLEM: %s\n\n", err);
}

void
format_warning(FILE *f, const char* seq_name, const char *err)
{
  if (NULL != seq_name)
    fprintf(f, "WARNINGS FOR %s\n\n", seq_name);
  if (err != NULL) 
    fprintf(f, "WARNING: %s\n\n", err);
}

/* Print the rows of oligos after the first rest_count ones of one
   oligo list, at most pa->num_return rows in all. Prints nothing if
   there is no such row. If rest_count is 1, the first table printed
   gets the heading "ADDITIONAL OLIGOS"; *additional_printed records
   that it was printed. If rest_count is 0 and newline_after is set,
   the table is followed by "\n ". */
static void
print_more_oligos(FILE *f,
                  const char *title,
                  const seq_args *sa,
                  const primer_rec *oligos,
                  int num_elem,
                  int type,
                  const p3_global_settings *pa,
                  const seq_lib *lib,
                  int print_lib_sim,
                  int rest_count,
                  int *additional_printed,
                  int newline_after)
{
  int i;
  int last = num_elem;

  if (pa->num_return < last) {
    last = pa->num_return;
  }
  if (last <= rest_count) {
    return;
  }
  if (rest_count == 1 && *additional_printed == 0) {
    fprintf(f, "ADDITIONAL OLIGOS\n");
    *additional_printed = 1;
  }
  fprintf(f, "   ");
  print_oligo_header(f, "", print_lib_sim, pa->thermodynamic_oligo_alignment);
  for (i = rest_count; i < last; i++) {
    fprintf(f, "%2d ", i + 1 - rest_count);
    print_oligo(f, title, sa, oligos + i, type, pa, lib, print_lib_sim);
  }
  if (rest_count == 0 && newline_after) {
    fprintf(f, "\n ");
  }
}

/* Format and print out one oligo */
static void 
format_oligos(FILE *f,
              const p3_global_settings *pa,
              const seq_args    *sa,
              const p3retval *retval,
              const char* pr_release,
              const pr_append_str *combined_retval_err,
              int explain_flag)
{
  char *warning;
  int print_lib_sim = lib_sim_specified(pa);
  int print_primers = 0;
  const primer_rec *best_left = NULL;
  const primer_rec *best_intl = NULL;
  const primer_rec *best_right = NULL;
  pair_array_t *best_pairs;
  int rest_count = 0;
  int additional_printed = 0;

  PR_ASSERT(NULL != f);
  PR_ASSERT(NULL != pa);
  PR_ASSERT(NULL != sa);

  best_pairs = NULL;

  if (!pr_is_empty(combined_retval_err)) {
    format_error(f, sa->sequence_name, pr_append_str_chars(combined_retval_err));
    return;
  }

  if (NULL != sa->sequence_name)
    fprintf(f, "PRIMER PICKING RESULTS FOR %s\n\n", sa->sequence_name);
  if (pa->pick_left_primer || pa->pick_right_primer) {
    if (pa->p_args.repeat_lib != NULL)
      fprintf(f, "Using mispriming library %s\n",
              pa->p_args.repeat_lib->repeat_file);
    else
      fprintf(f, "No mispriming library specified\n");
  } 
  if (pa->pick_internal_oligo) {
    if (pa->o_args.repeat_lib != NULL)
      fprintf(f, "Using internal oligo mishyb library %s\n",
              pa->o_args.repeat_lib->repeat_file);
    else
      fprintf(f, "No internal oligo mishyb library specified\n");
  }
  fprintf(f, "Using %d-based sequence positions\n",
          pa->first_base_index);

  if (pa->pick_left_primer) {
    if (retval->fwd.num_elem == 0){
        fprintf(f, "NO LEFT PRIMER FOUND\n\n");
    } else {
        print_primers = 1;
    }
  }
  if (pa->pick_internal_oligo) {
    if (retval->intl.num_elem == 0){
      fprintf(f, "NO INTERNAL OLIGO FOUND\n\n");
    } else {
      print_primers = 1;
    }
  }
  if (pa->pick_right_primer) {
    if (retval->rev.num_elem == 0){
      fprintf(f, "NO RIGHT PRIMER FOUND\n\n");
    } else {
      print_primers = 1;
    }
  }
  if ((warning = p3_get_rv_and_gs_warnings(retval, pa)) != NULL) {
    fprintf(f, "WARNING: %s\n\n", warning);
    free(warning);
  }
  if ((pa->primer_task != pick_primer_list )
         && (pa->primer_task != pick_sequencing_primers)) { 
    if (print_primers == 1) { 
      print_oligo_header(f, "OLIGO", print_lib_sim, pa->thermodynamic_oligo_alignment);
    }
    /* Print out the first line with the best primers */
    if ((pa->pick_left_primer) && (retval != NULL )
               && (retval->fwd.num_elem > 0)){
      print_oligo(f, "LEFT_PRIMER", sa, retval->fwd.oligo, FORWARD, 
                  pa, pa->p_args.repeat_lib, print_lib_sim);
      best_left = retval->fwd.oligo;
      rest_count = 1;
    }
    if ((pa->pick_internal_oligo) && (retval != NULL )
             && (retval->intl.num_elem > 0)){
      print_oligo(f, "INTERNAL_OLIGO", sa, retval->intl.oligo, INTERNAL, 
                  pa, pa->o_args.repeat_lib, print_lib_sim);
      best_intl = retval->intl.oligo;
      rest_count = 1;
    }
    if ((pa->pick_right_primer) && (retval != NULL )
             && (retval->rev.num_elem > 0)) {
      print_oligo(f, "RIGHT_PRIMER", sa, retval->rev.oligo, REVERSE, 
                  pa, pa->p_args.repeat_lib, print_lib_sim);
          best_right = retval->rev.oligo;
          rest_count = 1;
    }
  }
  if(print_primers == 1) { 
    fprintf(f, "SEQUENCE SIZE: %ld\n", (long int) strlen(sa->sequence));
    fprintf(f, "INCLUDED REGION SIZE: %d\n\n", sa->incl_l);
    print_pair_array(f, "TARGETS", sa->tar2.count, sa->tar2.pairs, pa, sa);
    print_pair_array(f, "EXCLUDED REGIONS", sa->excl2.count, sa->excl2.pairs, pa, sa);
    print_pair_array(f, "INTERNAL OLIGO EXCLUDED REGIONS",
            sa->excl_internal2.count, sa->excl_internal2.pairs, pa, sa);
    print_2_pair_array(f, "PAIR_OK_REGIONS", sa->ok_regions.count, sa->ok_regions.left_pairs, 
                       sa->ok_regions.right_pairs, pa, sa);
  }
  if (pa->primer_task != pick_primer_list ) {
    if (print_seq(f, pa, sa, retval, best_left, best_intl, best_right,
                  best_pairs, 0)) {
      exit(P3_EXIT_OUT_OF_MEMORY); /* ENOMEM */
    }
  }
  fprintf(f, "\n");
  /* Print out the other primers */
  if ((pa->pick_left_primer) && (retval != NULL)) {
    print_more_oligos(f, "LEFT_PRIMER", sa, retval->fwd.oligo,
                      retval->fwd.num_elem, FORWARD, pa,
                      pa->p_args.repeat_lib, print_lib_sim, rest_count,
                      &additional_printed, 1);
  }
  if ((pa->pick_internal_oligo) && (retval != NULL)) {
    print_more_oligos(f, "INTERNAL_OLIGO", sa, retval->intl.oligo,
                      retval->intl.num_elem, INTERNAL, pa,
                      pa->o_args.repeat_lib, print_lib_sim, rest_count,
                      &additional_printed, 1);
  }
  if ((pa->pick_right_primer) && (retval != NULL)) {
    print_more_oligos(f, "RIGHT_PRIMER", sa, retval->rev.oligo,
                      retval->rev.num_elem, REVERSE, pa,
                      pa->p_args.repeat_lib, print_lib_sim, rest_count,
                      &additional_printed, 0);
  }
  if (explain_flag) 
    print_explain(f, pa, sa, retval, print_lib_sim, pr_release);
  fprintf(f, "\n\n");
  if (fflush(f) == EOF) {
    perror("fflush(f) failed");
    exit(P3_EXIT_OUTPUT_ERROR);
  }
}

