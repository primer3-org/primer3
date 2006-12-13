#ifndef BOULDER_INPUT_H
#define BOULDER_INPUT_H 1
#include "primer3.h"

typedef struct program_args {
    char format_output;
    char twox_compat;
    char strict_tags;
} program_args;

int read_record(const program_args *, primer_args *, /*primer_args *,*/
		seq_args *);
void free_record(seq_args *sa);
void   free_seq_lib( seq_lib *);

/* Boulder format utility printing functions */
void   boulder_print_pairs(const program_args *, const primer_args *, const seq_args *, const pair_array_t *);
void   boulder_print_oligos(const primer_args *, 
				   const seq_args *, int, oligo_type,
				   primer_state *);
void   print_all_explain(const primer_args *, const seq_args *);
void   print_explain(const oligo_stats *, oligo_type);




#endif



