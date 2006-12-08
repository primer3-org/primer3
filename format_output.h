#ifndef PR_FORMAT_OUTPUT_H
#define PR_FORMAT_OUTPUT_H 1

#include "primer3.h"

/* 
 * Format the pair in p (plus the middle oligo if appropriate)
 * on f.  If p->product_size = 0 then primer choice failed.
 *
 * Returns 0 for success
 *         1 for failure
 */
int format_pairs(FILE *f, const primer_args *pa,
		 const seq_args *sa, 
		 const pair_array_t *);

int format_oligos(FILE *, const primer_args *, const seq_args *, 
		  primer_rec *, int, oligo_type);
#endif
