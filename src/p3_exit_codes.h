/* Exit status codes shared by primer3_core and the other programs
   (ntthal, ntdpal, oligotm, amplicon3_core, primer3_masker).
   Each code means the same in every program. They are documented in
   primer3_manual.htm, section "EXIT STATUS CODES"; changing a value
   changes the interface for callers. */
#ifndef P3_EXIT_CODES_H
#define P3_EXIT_CODES_H

/* Success */
#define P3_EXIT_OK            0
/* Invalid command line: unknown option, missing or invalid option
   value, wrong number of arguments, options that do not fit together */
#define P3_EXIT_USAGE         1
/* Out of memory */
#define P3_EXIT_OUT_OF_MEMORY 2
/* Empty input: no input record was found (primer3_core) */
#define P3_EXIT_EMPTY_INPUT   3
/* The input cannot be read or used: input file, settings, parameter
   or k-mer list files that cannot be read, an error in a global input
   tag, an invalid sequence, or a calculation that fails for the
   given input */
#define P3_EXIT_INPUT_ERROR   4
/* The output cannot be written: an output or error file cannot be
   created, or writing (flushing) the output fails */
#define P3_EXIT_OUTPUT_ERROR  5

#endif
