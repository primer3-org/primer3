/*
 * Copyright (c) 1996, Whitehead Institute for Biomedical Research. All rights
 * reserved.  Please see full software use agreement in primer3_main.c or by
 * executing primer3 with -h.
 */

#include "dpal.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

int
main(argc, argv)
    int argc;
    const char**argv;
{
    dpal_args a;
    dpal_results r;
    const char *s1, *s2;
    int i;
    int print_align_end = 0; /* 
                              * Print align_end_1 and align_end_2 from
                              * dpal_results.
			      */
    int use_ambiguity_codes = 0;
    int use_h_matrix = 0;
    char mode;
    char *endptr;
    const char *msg = 
      "Usage: %s [-g <gval>] [-l <lval>] [-m <mval>] [-h] [-p] [-s] [-e] <seq1> <seq2> <mode>\n"
      "where\n"
      "<gval> and <lval> are floats (.01 precision) specifying penalties for\n"
      "    creating or widening a gap respectively (<gval> and <lval> are\n"
      "    subtracted from the output score).\n"
      "-a causes the scoring matrix to be modified by dpal_set_ambiguity_codes.\n"
      "-e causes the end postion of the alignment in both sequences to\n"
      "   be printed.  Do not confuse with the 'e' <mode>\n"
      "-h use a different scoring matrix: G and C matches = 3, A and T = 2\n"
      "-p causes the alignment to be displayed.\n"
      "-s causes _only_ the score to printed.\n"
      "<mval> is the maximum allowed gap (default is 3).\n"
      "<seq1> and <seq2> are the sequences to be aligned.\n"
      "<mode> is one of g, G, l, or L specifying a global,\n"
      "       global end-anchored, local, or local end-achored\n"
      "       alignment respectively.  For backward compatibility\n"
      "       e is equivalent to G.\n";

    if (argc < 4) {
	fprintf(stderr, msg, argv[0]);
	exit(-1);
    }
    dpal_set_default_nt_args(&a);
    for (i=1; i < argc; ++i) {
	if (!strncmp("-p", argv[i], 2)) {
	    a.debug = 1;
	} else if (!strncmp("-l", argv[i], 2)) {
	    a.gapl = strtod(argv[i+1],(char **)NULL) * -100;
	    i++;
        } else if (!strncmp("-e", argv[i], 2)) {
	    print_align_end = 1;
	} else if (!strncmp("-a", argv[i], 2)) {
	    use_ambiguity_codes = 1;
	} else if (!strncmp("-h", argv[i], 2)) {
	    use_h_matrix = 1;
	} else if (!strncmp("-g", argv[i], 2)) {
	    a.gap = strtod(argv[i+1],(char **)NULL) * -100;
	    i++;
	} else if (!strncmp("-m", argv[i], 2)) {
	    a.max_gap = strtol(argv[i+1], &endptr, 10);
	    if ('\0' != *endptr) {
		fprintf(stderr, msg, argv[0]);
		exit(-1);
	    }
	    i++;
	} else if (!strncmp("-s", argv[i], 2)) {
	    a.score_only = 1;
	} else if (!strncmp("-e", argv[i], 2)) {
	    print_align_end = 1;
	} else if (!strncmp("-f1", argv[i], 3)) {
	    a.force_generic = 1;
	} else if (!strncmp("-f2", argv[i], 3)) {
	    a.force_long_generic = 1;
	} else if (!strncmp("-f3", argv[i], 3)) {
	    a.force_long_maxgap1 = 1;
	} else if (!strncmp("-", argv[i], 1)) {
	    /* Unknown option. */
	    fprintf(stderr, msg, argv[0]);
	    exit(-1);
	} else
	    break;		/* all args processed. go on to sequences. */
    }
    if (use_h_matrix) dpal_set_h_nt_matrix(&a);
    if (use_ambiguity_codes) dpal_set_ambiguity_code_matrix(&a);
    if (a.score_only && a.debug) {
	fprintf(stderr, msg, argv[0]);
	exit(-1);
    }
    a.flag = -1;
    
    if ((argc - i) != 3){
      /* Not enough remaining arguments. */
      fprintf(stderr, msg, argv[0]);
      exit(-1);
    }
    s1 = argv[i];
    s2 = argv[i+1];
    mode = *argv[i+2];
    if ('l' == mode)
	a.flag = DPAL_LOCAL;
    else if ('e' == mode || 'G' == mode)
	a.flag = DPAL_GLOBAL_END;
    else if ('g' == mode)
	a.flag = DPAL_GLOBAL;
    else if ('L' == mode)
	a.flag = DPAL_LOCAL_END;
    else {
	fprintf(stderr, msg, argv[0]);
	exit(-1);
    }
    if(print_align_end == 1) a.force_long_generic = 1;
    dpal((const unsigned char *)s1, (const unsigned char *)s2, &a, &r);
    printf("|%s|  |%s| %c ", s1, s2, mode); 
    if (a.score_only) {
	printf("%.2f\n", 0.01 * r.score);
	if (print_align_end) {
	  if(r.align_end_1 >= 0) printf("align_end_1=%d ",r.align_end_1);
	  if(r.align_end_2 >= 0) printf("align_end_2=%d\n ",r.align_end_2);
        }
    }
    else {
	printf("score=%.2f len=%d ", (0.01 * r.score), r.path_length);
	if (print_align_end) {
	  if(r.align_end_1 >= 0) printf("align_end_1=%d ",r.align_end_1);
	  if(r.align_end_2 >= 0) printf("align_end_2=%d ",r.align_end_2);
	}
	for (i=0; i<r.path_length; i++)
	    printf("|%d,%d", r.path[i][0],r.path[i][1]);
	printf("|\n");
    }
    return 0;
}
