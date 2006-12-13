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

