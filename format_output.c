/*
 * Copyright (c) 1997 Whitehead Institute for Biomedical Research. All rights
 * reserved.  Please see full software use agreement in primer3_main.c or by
 * executing primer3 with -h.
 */

#include <stdio.h>
#include <string.h>
#include "format_output.h"
#include "primer3_release.h"

#define FORWARD 1
#define REVERSE -1

static int lib_sim_specified(const primer_args *);
static void print_explain(FILE *, const primer_args *,
			  const seq_args *, int);
static void print_pair_info(FILE *, const primer_pair *,
			    const primer_args *);
static void print_oligo(FILE *, const char *, const seq_args *,
			const primer_rec *, int, const primer_args *, 
			const seq_lib, int);
static void print_oligo_header(FILE *, const char *, const int);
static void print_pair_array(FILE *, const char*, int,
			     const interval_array_t, 
			     const primer_args*, const seq_args*);
static void print_rest(FILE *, const primer_args *, 
		       const seq_args *,  const pair_array_t *);
static int print_seq(FILE *, const primer_args *, const seq_args *, 
		     primer_rec *h, const pair_array_t *, int);
static void print_seq_lines(FILE *, const char *s, const char *n, int, int,
			    int, const primer_args *);
static void print_stat_line(FILE *, const char *, oligo_stats s, int);
static void print_summary(FILE *, const primer_args *, 
			  const seq_args *, const pair_array_t *, int);
static void print_oligo_summary(FILE *, const primer_args *, 
			  const seq_args *, primer_rec *, 
			  oligo_type, int);

/*
 * ==========================================================================
 * External APIs
 * ==========================================================================
 */

/*
 * Returns 0 for success
 *         1 for failure
 */
int
format_pairs(FILE *f,
	     const primer_args *pa,
	     const seq_args *sa,
	     const pair_array_t *best_pairs)
{
    char *warning;
    int print_lib_sim = lib_sim_specified(pa);
    primer_rec *h;

    PR_ASSERT(NULL != f);
    PR_ASSERT(NULL != pa);
    PR_ASSERT(NULL != sa);

    h = NULL;
    if (NULL != sa->sequence_name)
	fprintf(f, "PRIMER PICKING RESULTS FOR %s\n\n", sa->sequence_name);

    if (sa->error.data != NULL) 
	fprintf(f, "INPUT PROBLEM: %s\n\n", sa->error.data);
    else {
	if (pa->repeat_lib.repeat_file != NULL)
	    fprintf(f, "Using mispriming library %s\n",
		    pa->repeat_lib.repeat_file);
	else
	    fprintf(f, "No mispriming library specified\n");

	if ( pa->primer_task == 1) {
	  if (pa->io_mishyb_library.repeat_file != NULL)
	    fprintf(f, "Using internal oligo mishyb library %s\n",
		    pa->io_mishyb_library.repeat_file);
	  else
	    fprintf(f, "No internal oligo mishyb library specified\n");
	}

	fprintf(f, "Using %d-based sequence positions\n",
		pa->first_base_index);
	if (best_pairs->num_pairs == 0) fprintf(f, "NO PRIMERS FOUND\n\n");
	if ((warning = pr_gather_warnings(sa, pa)) != NULL) {
	  fprintf(f, "WARNING: %s\n\n", warning);
	  free(warning);
	}
	print_summary(f, pa, sa, best_pairs, 0);
	fprintf(f, "\n");

	if (print_seq(f, pa, sa, h, best_pairs, 0))
	    return 1;
	if (best_pairs->num_pairs > 1 ) print_rest(f, pa, sa, best_pairs);
	if (pa->explain_flag) print_explain(f, pa, sa, print_lib_sim);
	fprintf(f, "\n\n");
	if (fflush(f) == EOF) {
	  perror("fflush(f) failed");
	  return 1;
	}
    }
    return 0;
}

/*
 * Returns 0 for success
 *         1 for failure
 */
int
format_oligos(FILE *f,
	      const primer_args *pa,
	      const seq_args    *sa,
	      primer_rec  *h,
	      int n,
	      oligo_type l)
{
  char *warning;
  int print_lib_sim = lib_sim_specified(pa);
  int i;
  pair_array_t *best_pairs;
  primer_rec *p;
  char type[20];

  PR_ASSERT(NULL != f);
  PR_ASSERT(NULL != pa);
  PR_ASSERT(NULL != sa);

  best_pairs = NULL;
  if (NULL != sa->sequence_name)
    fprintf(f, "PRIMER PICKING RESULTS FOR %s\n\n", sa->sequence_name);

  if (sa->error.data != NULL) 
    fprintf(f, "INPUT PROBLEM: %s\n\n", sa->error.data);
  else {
    if (l != OT_INTL ) {
      if (pa->repeat_lib.repeat_file != NULL)
	fprintf(f, "Using mispriming library %s\n",
		pa->repeat_lib.repeat_file);
      else
	fprintf(f, "No mispriming library specified\n");
    } else {
      if ( pa->primer_task == 1) {
	if (pa->io_mishyb_library.repeat_file != NULL)
	  fprintf(f, "Using internal oligo mishyb library %s\n",
		  pa->io_mishyb_library.repeat_file);
	else
	  fprintf(f, "No internal oligo mishyb library specified\n");
      }
    }
  }

  if(l == OT_LEFT) strcpy(type, "LEFT_PRIMER");
  else if(l == OT_RIGHT) strcpy(type, "RIGHT_PRIMER");
  else strcpy(type, "INTERNAL_OLIGO");

  fprintf(f, "Using %d-based sequence positions\n",
	  pa->first_base_index);
  if (n == 0) fprintf(f, "NO OLIGOS FOUND\n\n");
  if ((warning = pr_gather_warnings(sa, pa)) != NULL) {
    fprintf(f, "WARNING: %s\n\n", warning);
    free(warning);
  }

  if(n > 0) print_oligo_summary(f, pa, sa, h, l, 0);
  else h = NULL;
  if (print_seq(f, pa, sa, h, best_pairs, 0))
      return 1;
  fprintf(f, "\n");
  if(n > 1) {
     fprintf(f, "ADDITIONAL OLIGOS\n");
     fprintf(f, "   "); print_oligo_header(f, "", print_lib_sim);
     for(i = 1; i < pa->num_return; i++) {
	if(i > n-1) break;
	p = h + i;
	fprintf(f, "%2d ", i);
	if(OT_LEFT == l || OT_INTL == l)
	       print_oligo(f, type, sa, p, FORWARD,pa, pa->repeat_lib, print_lib_sim);
        else   print_oligo(f, type, sa, p, REVERSE,pa, pa->repeat_lib, print_lib_sim);
     }
   }
   if(pa->explain_flag) print_explain(f, pa, sa, print_lib_sim);
   fprintf(f, "\n\n");
   if (fflush(f) == EOF) {
      perror("fflush(f) failed");
      return 1;
   }

   return 0;
}

/*
 * ==========================================================================
 * Internal functions - mostly children of the read_record function.
 * These should all be static.
 * ==========================================================================
 */

static void
print_summary(f, pa, sa, best_pairs, num)
    FILE *f;
    const primer_args *pa;
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
	print_oligo(f, "LEFT PRIMER", sa, p->left, FORWARD, pa, pa->repeat_lib,
		    print_lib_sim);
	print_oligo(f, "RIGHT PRIMER", sa, p->right, REVERSE, pa, pa->repeat_lib,
		    print_lib_sim);
	if ( pa->primer_task == 1)
	    print_oligo(f, "INTERNAL OLIGO", sa, p->intl, FORWARD, pa, pa->io_mishyb_library,
			print_lib_sim);
    }
    fprintf(f, "SEQUENCE SIZE: %d\n", seq_len);
    fprintf(f, "INCLUDED REGION SIZE: %d\n\n", sa->incl_l);

    if (best_pairs->num_pairs > 0) print_pair_info(f, p, pa);
    print_pair_array(f, "TARGETS", sa->num_targets, sa->tar, pa, sa);
    print_pair_array(f, "EXCLUDED REGIONS", sa->num_excl, sa->excl, pa, sa);
    print_pair_array(f, "INTERNAL OLIGO EXCLUDED REGIONS",
		     sa->num_internal_excl, sa->excl_internal, pa, sa);
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

static void
print_oligo(f, title, sa, o, dir, pa, seqlib, print_lib_sim)
    FILE *f;
    const char *title;
    const seq_args *sa;
    const primer_rec *o;
    int dir;
    const primer_args *pa;
    const seq_lib seqlib;
    int print_lib_sim;
{
    const char *format1 = "%-16s %5d %4d %7.2f %7.2f %5.2f %5.2f ";
    char *seq = (FORWARD == dir) 
	? pr_oligo_sequence(sa, o) : pr_oligo_rev_c_sequence(sa, o);

    fprintf(f, format1,
	    title, o->start + sa->incl_s + pa->first_base_index,
	    o->length, o->temp, o->gc_content, 0.01 * o->self_any,
	    0.01 * o->self_end);

    if (print_lib_sim) {
	if (seqlib.repeat_file) 
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
    const primer_args *pa;
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

static int
print_seq(f, pa, sa, h, best_pairs, num)
    FILE *f;
    const primer_args *pa;
    const seq_args *sa;
    primer_rec *h;
    const pair_array_t *best_pairs;
    int num;  /* The number of primer pair to print. */
{
    int len, i, j, start;
    int something_found = 0, vector_found = 0;
    int *notes;
    char *notestr;
    primer_pair *p;
    p = NULL;
    if(pa->primer_task == pick_pcr_primers ||
       pa->primer_task == pick_pcr_primers_and_hyb_probe)
				      p = best_pairs->pairs + num;
    len = strlen(sa->sequence);
    if (NULL == (notes = malloc(sizeof(*notes) * len)))
	return 1;
    memset(notes, 0, sizeof(*notes) * len);
    if (NULL == (notestr = malloc(len + 1)))
	return 1;
    memset(notestr, ' ', len);
    notestr[len] = '\0';

    for (i = 0; i < len; i++) {
	if (i < sa->incl_s || i >= sa->incl_s + sa->incl_l)
	    notes[i] |= VECTOR;

	if ((pa->primer_task == pick_pcr_primers ||
	    pa->primer_task == pick_pcr_primers_and_hyb_probe) &&
	    best_pairs->num_pairs > 0) {
	    if (i >= p->left->start + sa->incl_s
		&& i < p->left->start + p->left->length + sa->incl_s)
		notes[i] |= LEFT_OLIGO;
	    if (i >= p->right->start - p->right->length + 1 + sa->incl_s
		&& i <= p->right->start + sa->incl_s)
		notes[i] |= RIGHT_OLIGO;
	    if ( pa->primer_task == 1
		&& i >= p->intl->start + sa->incl_s 
		&& i < p->intl->start + p->intl->length + sa->incl_s)
		notes[i] |= INTL_OLIGO;
	}
	else if (h != NULL) {
	    if(pa->primer_task == pick_left_only &&
	       i < h->start + h->length + sa->incl_s &&
	       i >= h->start + sa->incl_s)
	       notes[i] |= LEFT_OLIGO;
            else if(pa->primer_task == pick_right_only &&
	       i >= h->start - h->length + 1 + sa->incl_s
	       && i <= h->start + sa->incl_s)
	       notes[i] |= RIGHT_OLIGO;
            else if(pa->primer_task == pick_hyb_probe_only &&
	         i >= h->start + sa->incl_s                &&
	         i < h->start + h->length + sa->incl_s)
	       notes[i] |= INTL_OLIGO;
        }

	for (j = 0; j < sa->num_targets; j++) {
	    start = sa->tar[j][0] + sa->incl_s;
	    if (i >= start && i < start + sa->tar[j][1])
		notes[i] |= TARGET;
	}
	for (j = 0; j < sa->num_excl; j++) {
	    start = sa->excl[j][0] + sa->incl_s;
	    if (i >= start && i < start + sa->excl[j][1])
		notes[i] |= EXCL_REGION;
	}
	for (j = 0; j < sa->num_internal_excl; j++) {
	    start = sa->excl_internal[j][0] + sa->incl_s;
	    if (i >= start && i < start + sa->excl_internal[j][1])
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

    if (sa->num_excl > 0)
	fprintf(f, "XXXXXX excluded region\n");

    if (pa->primer_task == 1
	&& sa->num_internal_excl > 0)
	fprintf(f, "xxxxxx excluded region for internal oligo\n");

    if (sa->num_targets > 0)
	fprintf(f, "****** target\n");

    if ((pa->primer_task == pick_pcr_primers ||
	pa->primer_task == pick_pcr_primers_and_hyb_probe) &&
	best_pairs->num_pairs > 0) {
	   fprintf(f, ">>>>>> left primer\n");
	   fprintf(f, "<<<<<< right primer\n");
	   if ( pa->primer_task == 1)
	      fprintf(f, "^^^^^^ internal oligo\n");
    }
    else if (pa->primer_task == pick_left_only && h != NULL)
	   fprintf(f, ">>>>>> left primer\n");
    else if (pa->primer_task == pick_right_only && h != NULL)
	   fprintf(f, "<<<<<< right primer\n");
    else if (pa->primer_task == pick_hyb_probe_only && h != NULL)
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
    const primer_args *pa;
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
    const primer_args *pa;
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
    const primer_args *pa;
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
		    pa, pa->repeat_lib, print_lib_sim);
        fprintf(f, "   ");
	print_oligo(f, "RIGHT PRIMER", sa, best_pairs->pairs[i].right, REVERSE,
		    pa, pa->repeat_lib, print_lib_sim);
	if ( pa->primer_task == 1) {
            fprintf(f, "   ");
	    print_oligo(f, "INTERNAL OLIGO", sa, best_pairs->pairs[i].intl,
			FORWARD, pa, pa->io_mishyb_library, print_lib_sim);
	}
        if (best_pairs->pairs[i].product_size > 0) {
	    fprintf(f, "   ");
	    print_pair_info(f, &best_pairs->pairs[i], pa);
	}
    }
}

/* This function does _not_ print out the no_orf statistic. */
static void
print_explain(f, pa, sa, print_lib_sim)
    FILE *f;
    const primer_args *pa;
    const seq_args *sa;
    int print_lib_sim;
{
  const pair_stats *x;
  const char *format = print_lib_sim
    ? "%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s\n"
    : "%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s\n";

  fprintf(f, "\nStatistics\n");

  if (!pa->pick_anyway
      || !((pick_pcr_primers == pa->primer_task 
	   && sa->left_input && sa->right_input)
	  || (pick_pcr_primers_and_hyb_probe == pa->primer_task
	      && sa->left_input && sa->right_input && sa->internal_input)
	  || (pick_left_only == pa->primer_task
	      && sa->left_input)
	  || (pick_right_only == pa->primer_task
	      && sa->right_input)
	  || (pick_hyb_probe_only == pa->primer_task
	      && sa->internal_input))) {

    if (print_lib_sim) {
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

  if ((pick_pcr_primers == pa->primer_task
       || pick_left_only == pa->primer_task
       || pick_pcr_primers_and_hyb_probe == pa->primer_task)
      && !(pa->pick_anyway && sa->left_input))
    print_stat_line(f, "Left", sa->left_expl, print_lib_sim);

  if ((pick_pcr_primers == pa->primer_task
       || pick_right_only  == pa->primer_task
       || pick_pcr_primers_and_hyb_probe == pa->primer_task)
      && !(pa->pick_anyway && sa->right_input))
    print_stat_line(f, "Right", sa->right_expl, print_lib_sim);

  if ((pick_pcr_primers_and_hyb_probe == pa->primer_task
       || pick_hyb_probe_only == pa->primer_task)
      && !(pa->pick_anyway && sa->internal_input))
    print_stat_line(f, "Intl", sa->intl_expl, print_lib_sim);

  if (pick_pcr_primers == pa->primer_task
      || pick_pcr_primers_and_hyb_probe == pa->primer_task) {
    fprintf(f, "Pair Stats:\n");
    x = &sa->pair_expl;
    pr_print_pair_explain(f, sa);
  }
  fprintf(f, "%s\n", pr_release());
}

static void
print_stat_line(f, t, s, print_lib_sim)
    FILE *f;
    const char *t;
    oligo_stats s;
    int print_lib_sim;
{
    const char *format = print_lib_sim
	? "%-6s%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n"
	: "%-6s%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n";

    if (print_lib_sim)
	fprintf(f, format,
		t, s.considered, s.ns, s.target, s.excluded,
		s.gc, s.gc_clamp, s.temp_min, s.temp_max,
		s.compl_any, s.compl_end, s.repeat, s.poly_x, s.stability, s.ok);
    else 
	fprintf(f, format,
		t, s.considered, s.ns, s.target, s.excluded,
		s.gc, s.gc_clamp, s.temp_min, s.temp_max,
		s.compl_any, s.compl_end, s.poly_x, s.stability, s.ok);
}

/* 
 * Return true iff a check for library similarity has been specified for
 * either the primer pair or the internal oligo.
 */

static int
lib_sim_specified(pa)
  const primer_args *pa;
{
  return (pa->repeat_lib.repeat_file || pa->io_mishyb_library.repeat_file);
}

static void
print_oligo_summary(f, pa, sa, h, l, num)
    FILE *f;
    const primer_args *pa;
    const seq_args *sa;
    primer_rec *h;
    oligo_type l;
    int num;
{
    int seq_len = strlen(sa->sequence);
    int print_lib_sim = lib_sim_specified(pa);
    primer_rec *p;
    char type[20];

    if(l == OT_LEFT) strcpy(type, "LEFT_PRIMER");
    else if(l == OT_RIGHT) strcpy(type, "RIGHT_PRIMER");
    else strcpy(type, "INTERNAL_OLIGO");

    p = h + num;
	/* 
	 * If the following format changes, also change the format in
	 * print_oligo.
	 */
    print_oligo_header(f, "OLIGO", print_lib_sim);
    if(OT_LEFT == l || OT_INTL == l)
	    print_oligo(f, type, sa, p, FORWARD, pa, pa->repeat_lib,
		    print_lib_sim);
    else print_oligo(f, type, sa, p, REVERSE, pa, pa->repeat_lib,
		    print_lib_sim);
    
    fprintf(f, "SEQUENCE SIZE: %d\n", seq_len);
    fprintf(f, "INCLUDED REGION SIZE: %d\n\n", sa->incl_l);

    print_pair_array(f, "TARGETS", sa->num_targets, sa->tar, pa, sa);
    print_pair_array(f, "EXCLUDED REGIONS", sa->num_excl, sa->excl, pa, sa);
    print_pair_array(f, "INTERNAL OLIGO EXCLUDED REGIONS",
		     sa->num_internal_excl, sa->excl_internal, pa, sa);
}

