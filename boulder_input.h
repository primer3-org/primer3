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
#endif



