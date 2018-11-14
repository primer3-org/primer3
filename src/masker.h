#ifndef MASKER_H
#define MASKER_H


#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ALPHABET defines which characters are considered as nucleotides */
#define ALPHABET "ACGTUacgtu"

#define FWD 1
#define REV 2
#define BOTH 3

#define MAX_BUFFER_SIZE 5000
#define MAX_SPLITS 10

#define NUCLEOTIDE 0
#define WHITESPACE 1
#define MASKED_CHAR 2

#define FLUSH_ALL 1
#define KEEP_UNCERTAIN_REGION 0

#define DEFAULT_COEF 1.0

/* default formula parameter coefficients */
#define DEFAULT_NLISTS 2
#define DEFAULT_NLIST_PARAMETERS 2
#define DEFAULT_WORD_LEN_1 11
#define DEFAULT_WORD_LEN_2 16
#define DEFAULT_COEF_1 0.1772
#define DEFAULT_COEF_2 0.239
#define DEFAULT_INTERCEPT -4.336

/* default masking parameter values */
#define PRINT_SEQUENCE 1 /* used only for commandline tool */
#define DEFAULT_FAILURE_RATE 0.1
#define DEFAULT_ABS_CUTOFF 0 
#define DEFAULT_MASK_CHAR 'N'
#define DEFAULT_MASKING_DIRECTION both_on_same /* used only for commandline tool */
#define DEFAULT_M5P 1
#define DEFAULT_M3P 0
#define DEFAULT_LIST_FILE_PREFIX (char *)"homo_sapiens"

/* macros for moving forward and backward in masking buffer */
#define TAKE_STEP_BACK(i) \
  (i == 0) ? (MAX_BUFFER_SIZE - 1) : (i - 1)
  
#define TAKE_STEP_FORWARD(i) \
  (i == MAX_BUFFER_SIZE - 1) ? 0 : (i + 1)

struct pr_append_str;
  
/*
 * masking_direction: One can mask the forward, reverse or both strands.
 * both_on_same outputs one sequence that contains masked nucleotides according to both strands.
 * both_separately outputs two masked sequences separately.
 */
typedef enum masking_direction {
  both_on_same = 0,
  both_separately = 1,
  fwd = 2,
  rev = 3
} masking_direction;

/*
 * input_sequence: A wrapper for input sequence. Input can be either a stream or a 
 * string variable 
 */
typedef struct input_sequence {
  FILE *sequence_file;
  const char *sequence_string;
  size_t input_size;
  size_t current_pos;
} input_sequence;

/*
 * output_sequence: A wrapper for output sequence in case the output is given in a 
 * string variable
 */
typedef struct output_sequence {
  char *sequence;
  unsigned int pos;
  
  /* These will be used instead of char *sequence
   * in case of both_separately masking direction */
  char *sequence_fwd;
  char *sequence_rev;
} output_sequence;

/*
 * formula_parameters: All parameters concerning one input k-mer list 
 */
typedef struct formula_parameters {
  /* If the list is created with GenomeTester4, 
   * 210 char should be enough to contain the full list name */
  char list_file_name[210];
  unsigned int oligo_length;
  
  /* binary mask is used for cutting the k-mer into the size of oligo_length */
  unsigned long long binary_mask;
  
  /* number of unique k-mers in the given k-mer list */
  unsigned long long words_in_list;
  
  /* pointer to k-mer list */
  const char *word_list;
  const char *pointer;
  size_t size;
  
  /* coefficients for all possible masking formula variables (and their squares) 
   * concerning this k-mer list. 
   * If certain variables are not used, their coefficiest are equal to 0 */
  double mm0;
  double mm1;
  double mm2;
  double mm0_2;
  double mm1_2;
  double mm2_2;
} formula_parameters;

/*
 * parameters_builder: Helping structure for reading the formula parameters
 */
typedef struct parameters_builder {
  formula_parameters **fp_array;
  double intercept;
  char **used_lists;
  unsigned int nslots;
  unsigned int nfp;
} parameters_builder;

/*
 * masker_parameters: Masker parameters contain all the user specified 
 * (otherwise default) parameter values 
 */
typedef struct masker_parameters {
  /* strand to mask */
  masking_direction mdir;
  
  /* primer failure rate cutoff used in primer design, 
   * potential locations in a sequence for primers with PCR 
   * failure rate over the given cutoff are masked
   *(see function calculate_scores() from masker.c) */
  double failure_rate;
  
  /* absolute value cutoff, this can be used for masking all the k-mers in a sequence
   * that have the frequency over abs_cutoff in a k-mer list */
  unsigned int abs_cutoff;
  
  /* number of nucleotides masked in 5' and 3' direction with respect
   * to the 3' end of a primer */
  int nucl_masked_in_5p_direction;
  int nucl_masked_in_3p_direction;
  
  /* If masker is used as a separate application then always print_sequence=1, 
   * i.e the output is sent to stdout.
   * If print_sequence=0 the output is written in a string variable and can be forwarded
   * to the next function */
  int print_sequence;
  
  /* if do_soft_masking=1, masked nucleotides and converted to lower-case, else 
   * masked nucleotide are converted to masking_char ('N' by default) */
  int do_soft_masking;
  char masking_char;
  
  /* size of the masking window */
  int window_size;
  
  /* number of k-mer lists used in the masking formula */
  unsigned int nlists;
  /* k-mer lists and all their parameters which are used in the masking formula */
  char *list_prefix;
  formula_parameters **fp;
  double formula_intercept;
} masker_parameters;

/* 
 * masking_buffer: A round buffer for storing all the characters of the input sequence
 * as well as the masking and output information 
 */
typedef struct masking_buffer {
  /* the main array containing the characters that are read
   * from the sequence input */
  char buffer[MAX_BUFFER_SIZE];
  
  /* indicates the positions in the buffer that do not contain a
   * nucleotide character */
  int non_nucleotide_positions[MAX_BUFFER_SIZE];
  
  /* indicates the nucleotide positions in the buffer that
   * will be masked in he output */
  int mask_positions_fwd[MAX_BUFFER_SIZE];
  int mask_positions_rev[MAX_BUFFER_SIZE];
  
  /* indices for reading from the buffer and 
   writing to the buffer */
  unsigned int ri;
  unsigned int wi;
  
  /* index of the 5' end of the k-mer last entered plus the 
   * number of nucleotides masked in 3' direction. This index is used
   * to determine the first positions which can still be altered and therefore cannot be 
   * flushed out of the buffer yet */
  unsigned int ei;
  
  /* mi indicates the number of nucleotides that are
   * yet to be masked in 3' direction */
  unsigned int mi;
} masking_buffer;

/*
 * oligo_pair: contains k-mer pairs and their calculated scores
 */
typedef struct oligo_pair {
  /* fwd and rev k-mers */
  unsigned long long fwd;
  unsigned long long rev;
  
  /* a k-mer score is a failure rate of such a primer that 
   *contains a given k-mer in its 3' end */
  double score_fwd;
  double score_rev;
  
  /* abs_score is a number of given k-mers in a k-mer list */
  unsigned int abs_score;
} oligo_pair;

/*
 * oligo_counts: contains all k-mer frequencies (for both forward and reverse strand) 
 * retrieved from the k-mer lists 
 */
typedef struct oligo_counts {
  unsigned int oligo_length;
  unsigned int count_mm0_fwd;
  unsigned int count_mm1_fwd;
  unsigned int count_mm2_fwd;
  unsigned int count_mm0_rev;
  unsigned int count_mm1_rev;
  unsigned int count_mm2_rev;
} oligo_counts;

/*
 * 
 */
input_sequence *create_input_sequence_from_file_name (const char *input_file_name, pr_append_str *parse_err);
input_sequence *create_input_sequence_from_string (char *input_string, pr_append_str *parse_err);
void delete_input_sequence (input_sequence *input_seq);

int get_next_char_from_input (input_sequence *input_seq, unsigned long long *current_pos);
char *get_header_name_from_input (input_sequence *input_seq, unsigned long long header_pos, unsigned long long current_pos, pr_append_str *parse_err);

/*
 * 
 */
formula_parameters *create_formula_parameters_from_list_file_name (const char *list_file_name, pr_append_str *parse_err);
formula_parameters *create_formula_parameters_from_list_file_prefix (const char *list_name_prefix, const char *kmer_lists_path, unsigned int word_length, pr_append_str *parse_err);
formula_parameters **create_default_formula_parameters (const char *list_name_prefix, const char *kmer_lists_path, pr_append_str *parse_err);
formula_parameters **read_formula_parameters_from_file (const char *lists_file_name, unsigned int *nlist_parameters, parameters_builder *pbuilder, double *intercept, pr_append_str *parse_err);
int add_variable_to_formula_parameters (char **list_values, unsigned int nvalues, parameters_builder *pbuilder, pr_append_str *parse_err);
void delete_formula_parameters (formula_parameters **fp, unsigned int nlists);

/*
 * 
 */
output_sequence *create_output_sequence (unsigned long long seq_len, masking_direction mdir, pr_append_str *parse_err);
void delete_output_sequence (output_sequence *output_seq);

void write_header_to_output (output_sequence *output_seq, char *header_name, const masker_parameters *mp, pr_append_str *parse_err);
void write_char_to_output (output_sequence *output_seq, char c, char c_other, const masker_parameters *mp, pr_append_str *parse_err);

/*
 * 
 */
masking_buffer *create_masking_buffer (unsigned int word_length, pr_append_str *parse_err);
void initialize_masking_buffer (masking_buffer *mbuffer, unsigned int word_length);
void delete_masking_buffer (masking_buffer *mbuffer);
void add_char_to_buffer (char c, masking_buffer *mbuffer, int char_type);
void empty_buffer (output_sequence *output_seq, const masker_parameters *mp, masking_buffer *mbuffer, int flush_all, pr_append_str *parse_err);

/*
 * 
 */
void get_oligo_frequencies (oligo_counts *oc, formula_parameters *fp, unsigned long long word, unsigned int mm, int strand);

/*
 * 
 */
void calculate_scores (oligo_pair *h, const masker_parameters *mp, unsigned int word_length);
void mask_oligo_region (oligo_pair *h, const masker_parameters *mp, masking_buffer *mbuffer, unsigned int word_length, int debug);
void read_and_mask_sequence (input_sequence *input_seq, output_sequence *output_seq, const masker_parameters *mp, pr_append_str *parse_err, int debug);


const char *mmap_by_filename (const char *filename, size_t *size);
unsigned long long get_reverse_complement (unsigned long long word, unsigned int word_length);
/*
 * 
 */
unsigned long long create_binary_mask (unsigned int word_length);

/* 
 * converts nucleotide character into a 2-bit integer as follows:
 * 'A' -> 0; 'C' -> 1; 'G' -> 2; 'T'/'U' -> 3
 */
unsigned long long get_nucl_value (char nucl);

unsigned long long string_to_word (const char *s, unsigned int string_length, unsigned int word_length);
char *word_to_string (unsigned long long word, unsigned int wordlength);

void strip_string (char string[]);
char **split_string (char string[], char c, unsigned int *nchunks);


#endif









