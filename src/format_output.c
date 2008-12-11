/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008
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
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>
#include <string.h>
#include "format_output.h"
#include "libprimer3.h"

#define FORWARD 1
#define REVERSE -1

static char *pr_program_name = "Program name is probably primer3_core";

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
static void print_pair_info(FILE *, const primer_pair *,
                            const p3_global_settings *);
static void print_oligo(FILE *, const char *, const seq_args *,
                        const primer_rec *, int, const p3_global_settings *, 
                        const seq_lib*, int);
static void print_oligo_header(FILE *, const char *, const int);
static void print_pair_array(FILE *, const char*, int,
                             const interval_array_t, 
                             const p3_global_settings*, const seq_args*);
static void print_rest(FILE *, const p3_global_settings *, 
                       const seq_args *,  const pair_array_t *);
static int  print_seq(FILE *, const p3_global_settings *, const seq_args *, 
                      const p3retval *retval, primer_rec *h,
                      const pair_array_t *, int);
static void print_seq_lines(FILE *, const char *s, const char *n, int, int,
                            int, const p3_global_settings *);
static void print_stat_line(FILE *, const char *, oligo_stats s, int, int);
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
  if (NULL == combined_retval_err) exit(-2); /* Out of memory */

  if (pr_append_new_chunk_external(combined_retval_err, 
                                   retval->glob_err.data))
    exit(-2);

  if (pr_append_new_chunk_external(combined_retval_err, 
                                   retval->per_sequence_err.data)) 
    exit(-2);

  /* Print as primer pairs */
  if (retval->output_type == primer_pairs) {
    format_pairs(f, pa, sa, retval, &retval->best_pairs, 
                 pr_release, combined_retval_err, explain_flag);
    
    /* Print as primer list */
  } else {
    format_oligos(stdout, pa, sa, retval, pr_release,
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
  primer_rec *h = NULL;

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
  if (print_seq(f, pa, sa, retval, h, best_pairs, 0)) exit(-2); /* ENOMEM */
  
  /* Print out the alternative pairs */
  if (best_pairs->num_pairs > 1 ) print_rest(f, pa, sa, best_pairs);
  
  /* Print the primer picking statistics */
  if (explain_flag)
    print_explain(f, pa, sa, retval, print_lib_sim, pr_release);
  
  /* Flush the buffers and return */
  fprintf(f, "\n\n");
  if (fflush(f) == EOF) {
    perror("fflush(f) failed");
    exit(-1);
  }

}

/* Prints out the results of a primer pair */
static void
print_summary(f, pa, sa, best_pairs, num)
    FILE *f;
    const p3_global_settings *pa;
    const seq_args *sa;
    const pair_array_t *best_pairs;
    int num;
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
        print_oligo_header(f, "OLIGO", print_lib_sim);
        print_oligo(f, "LEFT PRIMER", sa, p->left, FORWARD, pa,
                    pa->p_args.repeat_lib,
                    print_lib_sim);
        print_oligo(f, "RIGHT PRIMER", sa, p->right, REVERSE, pa,
                    pa->p_args.repeat_lib,
                    print_lib_sim);
        if ( pa->pick_internal_oligo == 1 )
            print_oligo(f, "INTERNAL OLIGO", sa, p->intl, FORWARD,
                        pa, pa->o_args.repeat_lib,
                        print_lib_sim);
    }
    fprintf(f, "SEQUENCE SIZE: %d\n", seq_len);
    fprintf(f, "INCLUDED REGION SIZE: %d\n\n", sa->incl_l);

    if (best_pairs->num_pairs > 0) print_pair_info(f, p, pa);
    print_pair_array(f, "TARGETS", sa->tar2.count, sa->tar2.pairs, pa, sa);
    print_pair_array(f, "EXCLUDED REGIONS", sa->excl2.count, sa->excl2.pairs, pa, sa);
    print_pair_array(f, "INTERNAL OLIGO EXCLUDED REGIONS",
                     sa->excl_internal2.count, sa->excl_internal2.pairs, pa, sa);
}

/* Print column headers for lines printed by print_oligo(). */
static void
print_oligo_header(f, s, print_lib_sim)
    FILE *f;
    const char *s;
    const int print_lib_sim;
{
    fprintf(f,
            "%-16s start  len      tm     gc%%   any    3' %sseq\n",
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
    const char *format1 = "%-16s %5d %4d %7.2f %7.2f %5.2f %5.2f ";
    char *seq = (FORWARD == dir) 
        ? pr_oligo_sequence(sa, o) : pr_oligo_rev_c_sequence(sa, o);

    fprintf(f, format1,
            title, o->start + sa->incl_s + pa->first_base_index,
            o->length, o->temp, o->gc_content, 0.01 * o->self_any,
            0.01 * o->self_end);

    if (print_lib_sim) {
        if (seqlib != NULL) 
            fprintf(f, "%5.2f ",  0.01 * o->repeat_sim.score[o->repeat_sim.max]);
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
print_pair_array(f, title, num, array, pa, sa)
    FILE *f;
    const char* title;
    int num;
    const interval_array_t array;
    const p3_global_settings *pa;
    const seq_args *sa;
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

#define VECTOR           (1<<0)
#define LEFT_OLIGO       (1<<1)
#define RIGHT_OLIGO      (1<<2)
#define INTL_OLIGO       (1<<3)
#define TARGET           (1<<4)
#define EXCL_REGION      (1<<5)
#define INTL_EXCL_REGION (1<<6)

/* Prints out the asci picture of the sequence */
/* Return 1 on ENOMEM. Otherwise return 0. */
static int
print_seq(FILE *f,
    const p3_global_settings *pa,
    const seq_args *sa,
    const p3retval *retval,
    primer_rec *h,
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
    if (!(notes = malloc(sizeof(*notes) * len))) return 1;
    memset(notes, 0, sizeof(*notes) * len);
    if (!(notestr = malloc(len + 1))) return 1;
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
        else if (h != NULL) {
            if(pa->pick_left_primer == 1 &&
               i < h->start + h->length + sa->incl_s &&
               i >= h->start + sa->incl_s)
               notes[i] |= LEFT_OLIGO;
            else if(pa->pick_right_primer == 1 &&
               i >= h->start - h->length + 1 + sa->incl_s
               && i <= h->start + sa->incl_s)
               notes[i] |= RIGHT_OLIGO;
            else if(pa->pick_internal_oligo == 1 &&
                 i >= h->start + sa->incl_s                &&
                 i < h->start + h->length + sa->incl_s)
               notes[i] |= INTL_OLIGO;
        } else if (pa->primer_task == pick_sequencing_primers) {
    	    if (pa->pick_right_primer && &retval->rev != NULL 
    	             && retval->rev.num_elem > 0){
    	      h2 = retval->rev.oligo;
              for (j = 0; j < pa->num_return; j++) {
    	        if(j > retval->rev.num_elem -1) break;
    	        h3 = h2 + j;
                if(i >= h3->start - h3->length + 1 + sa->incl_s
                        && i <= h3->start + sa->incl_s)
                   notes[i] |= RIGHT_OLIGO;
    	      }
            }   	
    	    if (pa->pick_left_primer && &retval->fwd != NULL 
    	             && retval->fwd.num_elem > 0){
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
    else if (pa->pick_left_primer == 1 && h != NULL)
           fprintf(f, ">>>>>> left primer\n");
    else if (pa->pick_right_primer == 1 && h != NULL)
           fprintf(f, "<<<<<< right primer\n");
    else if (pa->pick_internal_oligo == 1 && h != NULL)
           fprintf(f, "^^^^^^ internal oligo\n");

    if (something_found) fputc('\n', f);
    free(notes);
    free(notestr);
    return 0;
}

static void
print_seq_lines(f, s, n, seq_size, line_size, something_found, pa)
    FILE *f;
    const char *s, *n;
    int seq_size, line_size, something_found;
    const p3_global_settings *pa;
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
print_pair_info(f, p, pa)
    FILE *f;
    const primer_pair *p;
    const p3_global_settings *pa;
{
  fprintf(f, "PRODUCT SIZE: %d, ", p->product_size);
  fprintf(f, "PAIR ANY COMPL: %.2f, PAIR 3' COMPL: %.2f\n",
          0.01 * p->compl_any, 0.01 * p->compl_end);

  if (pa->product_max_tm != PR_DEFAULT_PRODUCT_MAX_TM
      || pa->product_min_tm != PR_DEFAULT_PRODUCT_MIN_TM) {
    printf("PRODUCT Tm: %.4f, ", p->product_tm);
    printf("PRODUCT Tm - min(OLIGO Tm): %.4f\n",
           p->product_tm_oligo_tm_diff);
  }
}

static void
print_rest(f, pa, sa, best_pairs)
    FILE *f;
    const p3_global_settings *pa;
    const seq_args *sa;
    const pair_array_t *best_pairs;
{
    int i;
    int print_lib_sim = lib_sim_specified(pa);

    fprintf(f, "ADDITIONAL OLIGOS\n");
    fprintf(f, "   "); print_oligo_header(f, "", print_lib_sim);
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
                        FORWARD, pa, pa->o_args.repeat_lib, print_lib_sim);
        }
        if (best_pairs->pairs[i].product_size > 0) {
            fprintf(f, "   ");
            print_pair_info(f, &best_pairs->pairs[i], pa);
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


  char *format;

  if (print_lib_sim) {
    if (pa->lowercase_masking) {
      format = "%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s\n";
    } else {
      format = "%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s\n";
    }
  } else {
    if (pa->lowercase_masking) {
      format = "%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s\n";
    } else {
      format = "%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s\n";
    }
  }

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

    if (print_lib_sim) {
      if (pa->lowercase_masking) {
        fprintf(f, format,
                "", "con", "too",  "in",  "in",  "",    "no",
                "tm",  "tm",  "high", "high", "high",
                "", "high", " lower", "");
        fprintf(f, format,
                "", "sid", "many", "tar", "excl", "bad","GC",
                "too", "too", "any",  "3'", "lib",
                "poly", "end", " case", "");
        fprintf(f, format,
                "", "ered","Ns",   "get", "reg",  "GC%", "clamp",
                "low", "high","compl", "compl", "sim",
                "X",  "stab", " end", "ok  ");
      } else {
        fprintf(f, format,
                "", "con", "too",  "in",  "in",  "",    "no",
                "tm",  "tm",  "high", "high", "high",
                "", "high", "");
        fprintf(f, format,
                "", "sid", "many", "tar", "excl", "bad","GC",
                "too", "too", "any",  "3'", "lib",
                "poly", "end", "");
        fprintf(f, format,
                "", "ered","Ns",   "get", "reg",  "GC%", "clamp",
                "low", "high","compl", "compl", "sim",
                "X",  "stab", "ok");
      }
    } else {
      if (pa->lowercase_masking) {
        fprintf(f, format,
                "", "con", "too",  "in",  "in",  "",    "no",
                "tm",  "tm",  "high", "high",
                "", "high", " lower", "");
        fprintf(f, format,
                "", "sid", "many", "tar", "excl", "bad","GC",
                "too", "too", "any",  "3'",
                "poly", "end", " case", "");
        fprintf(f, format,
                "", "ered","Ns",   "get", "reg",  "GC%", "clamp",
                "low", "high","compl", "compl",
                "X", "stab", " end", "ok  ");
      } else {
        fprintf(f, format,
                "", "con", "too",  "in",  "in",  "",    "no",
                "tm",  "tm",  "high", "high",
                "", "high", "");
        fprintf(f, format,
                "", "sid", "many", "tar", "excl", "bad","GC",
                "too", "too", "any",  "3'",
                "poly", "end", "");
        fprintf(f, format,
                "", "ered","Ns",   "get", "reg",  "GC%", "clamp",
                "low", "high","compl", "compl",
                "X", "stab", "ok");
      }
    }

  }

  if (pa->pick_left_primer == 1
      && !(pa->pick_anyway && sa->left_input))
    print_stat_line(f, "Left", retval->fwd.expl, 
                    print_lib_sim, pa->lowercase_masking);

  if (pa->pick_right_primer == 1
      && !(pa->pick_anyway && sa->right_input))
    print_stat_line(f, "Right", retval->rev.expl,
                    print_lib_sim, pa->lowercase_masking);

  if (pa->pick_internal_oligo == 1
      && !(pa->pick_anyway && sa->internal_input))
    print_stat_line(f, "Intl", retval->intl.expl, 
                    print_lib_sim, pa->lowercase_masking);

  if (pa->pick_left_primer == 1 && pa->pick_right_primer == 1) {
    fprintf(f, "Pair Stats:\n%s\n",
            p3_get_pair_array_explain_string(p3_get_rv_best_pairs(retval)));
  }
  fprintf(f, "%s\n", pr_release);
}

static void
print_stat_line(f, t, s, print_lib_sim, lowercase_masking)
    FILE *f;
    const char *t;
    oligo_stats s;
    int print_lib_sim;
    int lowercase_masking;
{
  /* Modified by T. Koressaar to output statistics
     for lowercase masking. */
  if (print_lib_sim)
    if (lowercase_masking) {
      fprintf(f,
              "%-6s%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n",
              t, s.considered, s.ns, s.target, s.excluded,
              s.gc, s.gc_clamp, s.temp_min, s.temp_max,
              s.compl_any, s.compl_end, s.repeat_score,
              s.poly_x, s.stability, s.gmasked, s.ok);
    } else {
      fprintf(f,
              "%-6s%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n",
              t, s.considered, s.ns, s.target, s.excluded,
              s.gc, s.gc_clamp, s.temp_min, s.temp_max,
              s.compl_any, s.compl_end, s.repeat_score,
              s.poly_x, s.stability, s.ok);
    } else {
      if (lowercase_masking) {
        fprintf(f,
                "%-6s%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n",
                t, s.considered, s.ns, s.target, s.excluded,
                s.gc, s.gc_clamp, s.temp_min, s.temp_max,
                s.compl_any, s.compl_end, s.poly_x, s.stability, s.gmasked, s.ok);
      } else {
        fprintf(f,
                "%-6s%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n",
                t, s.considered, s.ns, s.target, s.excluded,
                s.gc, s.gc_clamp, s.temp_min, s.temp_max,
                s.compl_any, s.compl_end, s.poly_x, s.stability, s.ok);
      }
    }
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
  int i;
  int print_primers = 0;
  primer_rec  *h = NULL;
  pair_array_t *best_pairs;
  primer_rec *p;
  int rest_count = 0;

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
      print_oligo_header(f, "OLIGO", print_lib_sim);
    }
    /* Print out the first line with the best primers */
    if ((pa->pick_left_primer) && (&retval->fwd != NULL )
               && (retval->fwd.num_elem > 0)){
      print_oligo(f, "LEFT_PRIMER", sa, retval->fwd.oligo, FORWARD, 
                  pa, pa->p_args.repeat_lib, print_lib_sim);
      h = retval->fwd.oligo;
      rest_count = 1;
    }
    if ((pa->pick_internal_oligo) && (&retval->intl != NULL )
             && (retval->intl.num_elem > 0)){
      print_oligo(f, "INTERNAL_OLIGO", sa, retval->intl.oligo, FORWARD, 
                  pa, pa->p_args.repeat_lib, print_lib_sim);
      h = retval->intl.oligo;
      rest_count = 1;
    }
    if ((pa->pick_right_primer) && (&retval->rev != NULL )
             && (retval->rev.num_elem > 0)) {
      print_oligo(f, "RIGHT_PRIMER", sa, retval->rev.oligo, REVERSE, 
                  pa, pa->p_args.repeat_lib, print_lib_sim);
	  h = retval->rev.oligo;
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
  }
  if (pa->primer_task != pick_primer_list ) {
    if (print_seq(f, pa, sa, retval, h, best_pairs, 0)) exit(-2); /* ENOMEM */
  }
  fprintf(f, "\n");
  /* Print out the other primers */
  if ((pa->pick_left_primer) && (&retval->fwd != NULL )
             && (retval->fwd.num_elem > rest_count)){
    int n = retval->fwd.num_elem;
    h = retval->fwd.oligo;
    if (rest_count == 1) {  
      fprintf(f, "ADDITIONAL OLIGOS\n");
    }
    fprintf(f, "   "); print_oligo_header(f, "", print_lib_sim);
    for (i = rest_count; i < pa->num_return; i++) {
      if(i > n-1) break;
      p = h + i;
      fprintf(f, "%2d ", i + 1 - rest_count);
      print_oligo(f, "LEFT_PRIMER", sa, p, FORWARD, pa,
                  pa->p_args.repeat_lib, print_lib_sim);
    }
    if (rest_count == 0) {  
      fprintf(f, "\n ");
    }
  }
  if ((pa->pick_internal_oligo) && (&retval->intl != NULL )
           && (retval->intl.num_elem > rest_count)){
    int n = retval->intl.num_elem;
    h = retval->intl.oligo;  
    if (rest_count == 1) {  
      fprintf(f, "ADDITIONAL OLIGOS\n");
    }
    fprintf(f, "   "); print_oligo_header(f, "", print_lib_sim);
    for (i = rest_count; i < pa->num_return; i++) {
      if(i > n-1) break;
      p = h + i;
      fprintf(f, "%2d ", i + 1 - rest_count);
      print_oligo(f, "INTERNAL_OLIGO", sa, p, FORWARD, pa,
                  pa->p_args.repeat_lib, print_lib_sim);
      }
    if (rest_count == 0) {  
      fprintf(f, "\n ");
    }
  }
  if ((pa->pick_right_primer) && (&retval->rev != NULL )
           && (retval->rev.num_elem > rest_count)) {
    int n = retval->rev.num_elem;
    h = retval->rev.oligo; 
    if (rest_count == 1) {  
      fprintf(f, "ADDITIONAL OLIGOS\n");
    }
    fprintf(f, "   "); print_oligo_header(f, "", print_lib_sim);
    for (i = rest_count; i < pa->num_return; i++) {
      if(i > n-1) break;
      p = h + i;
      fprintf(f, "%2d ", i + 1 - rest_count);
      print_oligo(f, "RIGHT_PRIMER", sa, p, REVERSE, pa, 
                  pa->p_args.repeat_lib, print_lib_sim);
     }
  }
  if (explain_flag) 
    print_explain(f, pa, sa, retval, print_lib_sim, pr_release);
  fprintf(f, "\n\n");
  if (fflush(f) == EOF) {
    perror("fflush(f) failed");
    exit(-1);
  }
}

