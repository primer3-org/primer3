/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006
Whitehead Institute for Biomedical Research, Steve Rozen
(http://jura.wi.mit.edu/rozen), and Helen Skaletsky
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

   * Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the
distribution.
   * Neither the names of the copyright holders nor contributors may
be used to endorse or promote products derived from this software
without specific prior written permission.

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
#include <stdlib.h>
#include <string.h>
#include "oligotm.h"

/* Print the melting tm of an oligo on stdout. */
int
main(argc, argv)
    int argc;
    const char **argv;
{
    double tm;
    char *msg = "Usage: %s [-k salt-conc] [-d template-conc] oligo\n\n"
                "where oligo is an oligonucleotide sequence of between 2 and 36 bases.\n"
		"(Bases in oligo must be uppercase.)\n"
                "salt-conc is mM salt (usually K) concentration; defaults to 50mM\n"
                "template-conc is nM template concentration\n; defaults to 50nM\n"
		"Prints oligo's melting temperature on stdout.\n";
    char *endptr;
    long k = 50, d = 50;
    int i;

    if (argc < 2 || argc > 6) {
	fprintf(stderr, msg, argv[0]);
        return -1;
    }

    for (i=1; i < argc; ++i) {
	if (!strncmp("-k", argv[i], 2)) {
	    k = strtol(argv[i+1], &endptr, 10);
	    if ('\0' != *endptr) {
		fprintf(stderr, msg, argv[0]);
		exit(-1);
	    }
	    i++;
	} else if (!strncmp("-d", argv[i], 2)) {
	    d = strtol(argv[i+1], &endptr, 10);
	    if ('\0' != *endptr) {
		fprintf(stderr, msg, argv[0]);
		exit(-1);
	    }
	    i++;
	} else if (!strncmp("-", argv[i], 1)) {
	    /* Unknown option. */
	    fprintf(stderr, msg, argv[0]);
	    exit(-1);
	} else
	    break;		/* all args processed. go on to sequences. */
    }

    tm = oligotm(argv[i], d, k);
    if (OLIGOTM_ERROR == tm) {
	fprintf(stderr,
		"%s: length of %s is less than 2 or it contains an illegal character\n",
		argv[0], argv[i]);
        return -1;
    }
    fprintf(stdout, "%f\n", tm);
    return 0;
}




