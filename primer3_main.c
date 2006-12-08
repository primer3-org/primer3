/* 
 * Copyright (c) 1996, 1997, 1998 Whitehead Institute for Biomedical Research. All
 * rights reserved.  Please see full software use agreement in primer3_main.c or
 * by executing primer3 with -h.
 */

#include <stdio.h>
#include <signal.h>
#include <unistd.h>
#include "primer3.h"
#include "boulder_input.h"
#include "primer3_release.h"
#include "format_output.h"

static void   sig_handler(int);
static void   print_usage();
static void   boulder_print_pairs(const program_args *, const primer_args *, 
				  const seq_args *, const pair_array_t *);
static void   boulder_print_oligos(const primer_args *, 
				   const seq_args *, int, oligo_type,
				   primer_state *);
static void   print_all_explain(const primer_args *, const seq_args *);
static void   print_explain(const oligo_stats *, oligo_type);
static void *pr_safe_malloc(size_t); /* A fail-stop wrapper for malloc. */

/* Global variables - set in main */
static const char *pr_program_name;
static int pr_program_name_len;

typedef struct oligo_array {
  int len;
  primer_rec *data;
} oligo_array;

#define FREE_STUFF { \
    free(sa); 						      \
    free_seq_lib(&pa->repeat_lib);                            \
    free_seq_lib(&pa->io_mishyb_library); free(pa);           \
}

int
main(argc,argv)
    int argc;
    char *argv[]; 
{ 
    program_args prog_args;
    primer_args *pa;
    seq_args *sa;
    int input_found=0;
    primer_state *state;

    pr_program_name = argv[0];
    pr_program_name_len = strlen(argv[0]);
    /* 
     * We set up some signal handlers in case someone starts up the program
     * from the command line, wonders why nothing is happening, and then kills
     * the program.
     */
    signal(SIGINT, sig_handler);
    signal(SIGTERM, sig_handler);

    /* 
     * We allocate the following structures on the heap rather than on the
     * stack in order to take advantage of testcenter memory access checking
     * during testing.
     */
    pa = pr_safe_malloc(sizeof(*pa));
    sa = pr_safe_malloc(sizeof(*sa));

    set_default_global_primer_args(pa);

    memset(&prog_args, 0, sizeof(prog_args));
    while (--argc > 0) {
	argv++;
	if (!strcmp(*argv, "-format_output"))
	    prog_args.format_output = 1;
	else if (!strcmp(*argv, "-2x_compat"))
	    prog_args.twox_compat = 1;
	else if (!strcmp(*argv, "-strict_tags"))
	    prog_args.strict_tags = 1;
	else  {
	    print_usage();
	    FREE_STUFF;
	    exit(-1);
	}
    }

    state = primer3_create();
    if (!state)
	exit(-2);

    /* 
     * Read the data from input stream record by record and process it if
     * there are no errors.
     */
    while ((read_record(&prog_args, pa, sa)) > 0) {
	input_found = 1;
	state->n_f = state->n_m = state->n_r = 0;
	if (NULL == sa->error.data && NULL == pa->glob_err.data) {
	    primer3_choose(state, pa, sa);
        }
	if (NULL != pa->glob_err.data) 
	    pr_append_new_chunk(&sa->error, pa->glob_err.data);

	switch(pa->primer_task) {
	case pick_pcr_primers:
	case pick_pcr_primers_and_hyb_probe:
	    if (prog_args.format_output)
		format_pairs(stdout, pa, sa, &state->best_pairs);
	    else
		boulder_print_pairs(&prog_args, pa, sa, &state->best_pairs);
	    break;

	case pick_left_only:
	    if (prog_args.format_output) 
		format_oligos(stdout, pa, sa, state->f, state->n_f, OT_LEFT);
	    else
		boulder_print_oligos(pa, sa, state->n_f, OT_LEFT, state);
	    break;

	case pick_right_only:
	    if (prog_args.format_output) 
		format_oligos(stdout, pa, sa, state->r, state->n_r, OT_RIGHT);
	    else
		boulder_print_oligos(pa, sa, state->n_r, OT_RIGHT, state);
	    break;
	    
	case pick_hyb_probe_only:
	    if(prog_args.format_output) 
		format_oligos(stdout, pa, sa, state->mid, state->n_m, OT_INTL);
	    else
		boulder_print_oligos(pa, sa, state->n_m, OT_INTL, state);
	    break;

	default:
	    fprintf(stderr, "Unknown primer_task %d\n", pa->primer_task);
        }

	free_record(sa);

	if (NULL != pa->glob_err.data) {
	    fprintf(stderr, "%s: %s\n", pr_program_name, pa->glob_err.data);
	    free(pa->glob_err.data);
	    FREE_STUFF;
	    exit(-4);
	}
    }
    FREE_STUFF; 

    primer3_destroy(state);
    if (0 == input_found) {
	print_usage();
	exit(-3);
    }
    return 0;
}
#undef FREE_STUFF

static void
sig_handler(signal)
    int signal;
{
    print_usage();
    fprintf(stderr, "%s: received signal %d\n", pr_program_name, signal);
    exit(signal);
}

static void
print_usage()
{
    fprintf(stderr, "%s", pr_copyright());
    fprintf(stderr, 
	    "\n\nUSAGE: %s %s %s\n", pr_program_name,
	    "[-format_output] [-2x_compat]",
	    "[-strict_tags]");
    fprintf(stderr, "This is primer3 (%s)\n", pr_release());
    fprintf(stderr, "Input must be provided on standard input.\n");
}

/* Print the data for chosen primer pairs to stdout in "boulderio" format. */
static void
boulder_print_pairs(prog_args, pa, sa, best_pairs)
    const program_args *prog_args;
    const primer_args *pa;
    const seq_args *sa;
    const pair_array_t *best_pairs;
{
    const char *left_tag, *right_tag, *intl_tag, *prod_size_tag;
    char *warning;
    char suffix [3];
    primer_rec *fwd, *rev, *intl;
    int i, incl_s = sa->incl_s;

    PR_ASSERT(NULL != pa);
    PR_ASSERT(NULL != sa);
    PR_ASSERT(NULL != prog_args);

    if (prog_args->twox_compat) {
	left_tag = "FORWARD_PRIMER";
	right_tag = "REVERSE_PRIMER";
	intl_tag = "MIDDLE_OLIGO";
	prod_size_tag = "PRODUCT_SIZE";
    } else {
	left_tag = "PRIMER_LEFT";
	right_tag = "PRIMER_RIGHT";
	intl_tag = "PRIMER_INTERNAL_OLIGO";
	prod_size_tag = "PRIMER_PRODUCT_SIZE";
    }

    if ((warning = pr_gather_warnings(sa, pa)) != NULL) {
      
	printf("PRIMER_WARNING=%s\n", warning);
	free(warning);
    }
    
    if (sa->error.data != NULL) {
	printf("PRIMER_ERROR=%s\n=\n", sa->error.data);
	if (fflush(stdout) == EOF) {
	    perror("fflush(stdout) failed");
	    exit(-1);
	}
	return;
    }

    if (pa->explain_flag) print_all_explain(pa, sa);

    if (!PR_START_CODON_POS_IS_NULL(sa))
      printf("PRIMER_STOP_CODON_POSITION=%d\n", sa->stop_codon_pos);

    for(i=0; i<best_pairs->num_pairs; i++) {
	fwd = best_pairs->pairs[i].left;
	rev = best_pairs->pairs[i].right;
	intl = best_pairs->pairs[i].intl;

	if (i == 0) suffix[0] = '\0';
	else sprintf(suffix, "_%d", i);

	printf("PRIMER_PAIR_PENALTY%s=%.4f\n", suffix,
	       best_pairs->pairs[i].pair_quality);
	printf("PRIMER_LEFT%s_PENALTY=%f\n", suffix,
	       fwd->quality);
	printf("PRIMER_RIGHT%s_PENALTY=%f\n", suffix,
		   rev->quality);
	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe)
	  printf("PRIMER_INTERNAL_OLIGO%s_PENALTY=%f\n", suffix,
		 intl->quality);

        /* Print sequences. */
	printf("PRIMER_LEFT%s_SEQUENCE=%s\n", suffix,
	       pr_oligo_sequence(sa, fwd));
	printf("PRIMER_RIGHT%s_SEQUENCE=%s\n", suffix,
	       pr_oligo_rev_c_sequence(sa, rev));
	if( pa->primer_task == pick_pcr_primers_and_hyb_probe)
#if USE_OLD_FORMAT_MISTAKE
	    printf("PRIMER_INTERNAL%s_OLIGO_SEQUENCE=%s\n", suffix,
		   pr_oligo_sequence(sa,intl));
#else
	    printf("PRIMER_INTERNAL_OLIGO%s_SEQUENCE=%s\n", suffix,
		   pr_oligo_sequence(sa,intl));
#endif

	printf("%s%s=%d,%d\n", left_tag, suffix,
	       fwd->start + incl_s + pa->first_base_index,
	       fwd->length);
	printf("%s%s=%d,%d\n", right_tag, suffix,
	       rev->start + incl_s + pa->first_base_index,
	       rev->length);
	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe)
	    printf("%s%s=%d,%d\n", intl_tag, suffix,
		   intl->start + incl_s + pa->first_base_index,
		   intl->length);

	printf("PRIMER_LEFT%s_TM=%.3f\n", suffix, fwd->temp);
	printf("PRIMER_RIGHT%s_TM=%.3f\n", suffix, rev->temp);
	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe)
	    printf("PRIMER_INTERNAL_OLIGO%s_TM=%.3f\n",suffix, intl->temp);


	printf("PRIMER_LEFT%s_GC_PERCENT=%.3f\n", suffix, fwd->gc_content);
	printf("PRIMER_RIGHT%s_GC_PERCENT=%.3f\n", suffix, rev->gc_content);
	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe)
	    printf("PRIMER_INTERNAL_OLIGO%s_GC_PERCENT=%.3f\n",suffix,
		   intl->gc_content);

	printf("PRIMER_LEFT%s_SELF_ANY=%.2f\n", suffix,
	       fwd->self_any / PR_ALIGN_SCORE_PRECISION);
	printf("PRIMER_RIGHT%s_SELF_ANY=%.2f\n", suffix,
	       rev->self_any / PR_ALIGN_SCORE_PRECISION);
	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe)
	    printf("PRIMER_INTERNAL_OLIGO%s_SELF_ANY=%.2f\n", suffix,
		   intl->self_any / PR_ALIGN_SCORE_PRECISION);
	printf("PRIMER_LEFT%s_SELF_END=%.2f\n", suffix,
	       fwd->self_end / PR_ALIGN_SCORE_PRECISION);
	printf("PRIMER_RIGHT%s_SELF_END=%.2f\n",
	       suffix,rev->self_end / PR_ALIGN_SCORE_PRECISION);
	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe)
	    printf("PRIMER_INTERNAL_OLIGO%s_SELF_END=%.2f\n", suffix,
		   intl->self_end / PR_ALIGN_SCORE_PRECISION);
        if (pa->repeat_lib.seq_num > 0) {
	    printf("PRIMER_LEFT%s_MISPRIMING_SCORE=%.2f, %s\n", suffix,
		   fwd->repeat_sim.score[fwd->repeat_sim.max] / PR_ALIGN_SCORE_PRECISION,
		   fwd->repeat_sim.name);
            printf("PRIMER_RIGHT%s_MISPRIMING_SCORE=%.2f, %s\n", suffix,
		   rev->repeat_sim.score[rev->repeat_sim.max] / PR_ALIGN_SCORE_PRECISION,
                   rev->repeat_sim.name);
            printf("PRIMER_PAIR%s_MISPRIMING_SCORE=%.2f, %s\n", suffix,
		   best_pairs->pairs[i].repeat_sim / PR_ALIGN_SCORE_PRECISION,
		   best_pairs->pairs[i].rep_name);
        }
	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe
	    && pa->io_mishyb_library.seq_num > 0)
	    printf("PRIMER_INTERNAL_OLIGO%s_MISHYB_SCORE=%.2f, %s\n", suffix,
		   intl->repeat_sim.score[intl->repeat_sim.max]
		   / PR_ALIGN_SCORE_PRECISION,
		   intl->repeat_sim.name);
	if (NULL != sa->quality){
	   printf("PRIMER_LEFT%s_MIN_SEQ_QUALITY=%d\n", suffix,
		   fwd->seq_quality);
           printf("PRIMER_RIGHT%s_MIN_SEQ_QUALITY=%d\n", suffix,
		   rev->seq_quality);
        }

	if (!DEFAULT_POSITION_PENALTIES(pa)
            || !PR_START_CODON_POS_IS_NULL(sa)) {
	  printf("PRIMER_LEFT%s_POSITION_PENALTY=%f\n", suffix,
		 fwd->position_penalty);
	  printf("PRIMER_RIGHT%s_POSITION_PENALTY=%f\n", suffix,
		 rev->position_penalty);
	}

	printf("PRIMER_LEFT%s_END_STABILITY=%.4f\n",
	       suffix, fwd->end_stability);
	printf("PRIMER_RIGHT%s_END_STABILITY=%.4f\n",
	       suffix, rev->end_stability);

	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe
	    && NULL != sa->quality)
	   printf("PRIMER_INTERNAL_OLIGO%s_MIN_SEQ_QUALITY=%d\n",
		   suffix, intl->seq_quality);
	printf("PRIMER_PAIR%s_COMPL_ANY=%.2f\n", suffix,
	       best_pairs->pairs[i].compl_any / PR_ALIGN_SCORE_PRECISION);
	printf("PRIMER_PAIR%s_COMPL_END=%.2f\n", suffix,
	       best_pairs->pairs[i].compl_end  / PR_ALIGN_SCORE_PRECISION);

	printf("%s%s=%d\n", prod_size_tag, suffix,
				     best_pairs->pairs[i].product_size);

	if (pa->product_max_tm != PR_DEFAULT_PRODUCT_MAX_TM
	    || pa->product_min_tm != PR_DEFAULT_PRODUCT_MIN_TM) {
	  printf("PRIMER_PRODUCT_TM%s=%.4f\n", suffix,
		 best_pairs->pairs[i].product_tm);

	  printf("PRIMER_PRODUCT_TM_OLIGO_TM_DIFF%s=%.4f\n", suffix,
	      best_pairs->pairs[i].product_tm_oligo_tm_diff);

	   printf("PRIMER_PAIR_T_OPT_A%s=%.4f\n", suffix,
	      best_pairs->pairs[i].t_opt_a);
	}


    }
    printf("=\n");
    if (fflush(stdout) == EOF) {
	perror("fflush(stdout) failed");
	exit(-1);
    }
}

/*
 * Note this accesses the internals of the primer_state structure, which
 * implies this function should be in the library or we need to have a better
 * return mechanism.
 */
static void 
boulder_print_oligos(pa, sa, n, l, state)
    const primer_args *pa;
    const seq_args *sa;
    int n;
    oligo_type l;
    primer_state *state;
{
    char *warning;
    int i, j;
    char suffix [3], type[256];
    /* type must be larger than the length of "PRIMER_INTERNAL_OLIGO". */

    primer_rec *oligo;
    int incl_s = sa->incl_s;


    if ((warning = pr_gather_warnings(sa, pa)) != NULL) {
	printf("PRIMER_WARNING=%s\n", warning);
	free(warning);
    }
    if (sa->error.data != NULL) {
	printf("PRIMER_ERROR=%s\n=\n", sa->error.data);
	if (fflush(stdout) == EOF) {
	    perror("fflush(stdout) failed");
	    exit(-1);
        }
	return;
    }

    if(l == OT_LEFT) strcpy(type, "PRIMER_LEFT");
    else if(l == OT_RIGHT) strcpy(type, "PRIMER_RIGHT");
    else strcpy(type, "PRIMER_INTERNAL_OLIGO");

    if (pa->explain_flag) print_all_explain(pa, sa);

    i = 0;
    j = (pa->num_return < n) ? pa->num_return : n;
    if(l == OT_LEFT) oligo = state->f;
    else if (l == OT_RIGHT) oligo = state->r;
    else  oligo = state->mid;

    while(i < j) {
	if (i == 0) suffix[0] = '\0';
	else sprintf(suffix, "_%d", i);

	printf("%s%s_PENALTY=%.4f\n", type, suffix, oligo[i].quality);
	if(l == OT_RIGHT)
	   printf("%s%s_SEQUENCE=%s\n", type, suffix,
		  pr_oligo_rev_c_sequence(sa, &oligo[i]));
        else printf("%s%s_SEQUENCE=%s\n", type, suffix,
		    pr_oligo_sequence(sa, &oligo[i])); 
	printf("%s%s=%d,%d\n", type, suffix, 
		       oligo[i].start + incl_s + pa->first_base_index,
		       oligo[i].length);
	printf("%s%s_TM=%.3f\n", type, suffix, oligo[i].temp);
	printf("%s%s_GC_PERCENT=%.3f\n", type, suffix, oligo[i].gc_content);
        printf("%s%s_SELF_ANY=%.2f\n", type, suffix,
		       oligo[i].self_any / PR_ALIGN_SCORE_PRECISION);
        printf("%s%s_SELF_END=%.2f\n", type, suffix,
		       oligo[i].self_end / PR_ALIGN_SCORE_PRECISION);
        if ((l == OT_LEFT || l == OT_RIGHT) && pa->repeat_lib.seq_num > 0 ) 
	    printf("%s%s_MISPRIMING_SCORE=%.2f, %s\n", type, suffix,
		    oligo[i].repeat_sim.score[oligo[i].repeat_sim.max] /PR_ALIGN_SCORE_PRECISION,
		    oligo[i].repeat_sim.name);
        if (l == OT_INTL && pa->io_mishyb_library.seq_num > 0)
	    printf("%s%s_MISHYB_SCORE=%.2f,%s\n", type, suffix,
	           oligo[i].repeat_sim.score[oligo[i].repeat_sim.max]/ PR_ALIGN_SCORE_PRECISION,
                   oligo[i].repeat_sim.name);
	if (NULL != sa->quality)printf("%s%s_MIN_SEQ_QUALITY=%d\n",
		   type, suffix, oligo[i].seq_quality);
        if (PR_DEFAULT_INSIDE_PENALTY != pa->inside_penalty
	    || PR_DEFAULT_OUTSIDE_PENALTY != pa->outside_penalty != 0.0)
	    printf("%s%s_POSITION_PENALTY=%f\n", type, suffix,
		    oligo[i].position_penalty);
	if(l == OT_LEFT || l == OT_RIGHT)
	    printf("%s%s_END_STABILITY=%.4f\n", type, suffix,
		    oligo[i].end_stability);

        i++;
    }
    printf("=\n");
}

static void
print_all_explain(pa, sa)
    const primer_args *pa;
    const seq_args *sa;
{
  if (pa->explain_flag) {
    if (pa->primer_task != pick_right_only
	&& pa->primer_task != pick_hyb_probe_only
	&& !(pa->pick_anyway && sa->left_input))
      print_explain(&sa->left_expl,OT_LEFT);

    if (pa->primer_task != pick_left_only 
	&& pa->primer_task != pick_hyb_probe_only
	&& !(pa->pick_anyway && sa->right_input))
      print_explain(&sa->right_expl,OT_RIGHT);

    if ((pa->primer_task == pick_hyb_probe_only
	 || pa->primer_task == pick_pcr_primers_and_hyb_probe)
	&& !(pa->pick_anyway && sa->internal_input)) 
      print_explain(&sa->intl_expl, OT_INTL);

    if (pa->primer_task  == pick_pcr_primers
	|| pa->primer_task == pick_pcr_primers_and_hyb_probe) {
      printf("PRIMER_PAIR_EXPLAIN=");
      pr_print_pair_explain(stdout, sa);
    }
  }
}

static void
print_explain(stat, l)
    const oligo_stats *stat;
    oligo_type l;
{
    if(OT_LEFT == l)printf("PRIMER_LEFT_EXPLAIN=");
    else if(OT_RIGHT == l)printf("PRIMER_RIGHT_EXPLAIN=");
	 else printf("PRIMER_INTERNAL_OLIGO_EXPLAIN=");

    printf("considered %d", stat->considered);
    if (stat->no_orf) printf(", would not amplify any of the ORF %d", stat->no_orf);
    if(stat->ns)printf(", too many Ns %d", stat->ns);
    if(stat->target)printf(", overlap target %d", stat->target);
    if(stat->excluded)printf(", overlap excluded region %d", stat->excluded);
    if(stat->gc)printf(", GC content failed %d", stat->gc);
    if(stat->gc_clamp)printf(", GC clamp failed %d", stat->gc_clamp);
    if(stat->temp_min)printf(", low tm %d", stat->temp_min);
    if(stat->temp_max)printf(", high tm %d", stat->temp_max);
    if(stat->compl_any)printf(", high any compl %d", stat->compl_any);
    if(stat->compl_end)printf(", high end compl %d", stat->compl_end);
    if(stat->repeat)printf(", high repeat similarity %d", stat->repeat);
    if(stat->poly_x)printf(", long poly-x seq %d", stat->poly_x);
    if(stat->seq_quality)printf(",low sequence quality %d", stat->seq_quality);
    if (stat->stability) printf(",high 3' stability %d", stat->stability);

    printf(", ok %d\n", stat->ok);
}

/* 
 * Panic messages for when the program runs out of memory.  pr_program_name and
 * pr_program_name_len must be set at the beginning of main.
 */
#define OOM_MESSAGE      ": out of memory\n"
#define OOM_MESSAGE_LEN  16
#define OOM_ERROR \
    write(2, pr_program_name, pr_program_name_len), \
    write(2, OOM_MESSAGE, OOM_MESSAGE_LEN),   \
    exit(-2)

static void *
pr_safe_malloc(size_t x)
{
    void *r = malloc(x);
    if (NULL == r) OOM_ERROR;
    return r;
}

