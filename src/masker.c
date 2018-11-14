#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <math.h>
#include <unistd.h>
#include <errno.h>

#include "libprimer3.h"

unsigned int glistmaker_code_match = 'G' << 24 | 'T' << 16 | '4' << 8 | 'C';

/*==================================================================================
 * 
 * Input wrapper functions:
 * 
 *==================================================================================*/ 

input_sequence *
create_input_sequence_from_file_name (const char *input_file_name, pr_append_str *parse_err) 
{
  input_sequence *input_seq = (input_sequence *) malloc (sizeof(input_sequence));
  memset (input_seq, 0, sizeof(input_sequence));
  if (!input_file_name) input_seq->sequence_file = stdin;
  else input_seq->sequence_file = fopen(input_file_name, "r");
  if (!input_seq->sequence_file) {
    pr_append_new_chunk_external (parse_err, "Input file not found: ");
    pr_append_external (parse_err, input_file_name);
    return NULL;
  }
  return input_seq;
}

input_sequence *
create_input_sequence_from_string (char *input_string, pr_append_str *parse_err) 
{
  input_sequence *input_seq = (input_sequence *) malloc (sizeof(input_sequence));
  if (!input_seq) {
    pr_append_new_chunk_external (parse_err, "Memory allocation for input sequence failed!");
    return input_seq;
  }
  memset (input_seq, 0, sizeof(input_sequence));
  input_seq->sequence_string = input_string;
  input_seq->input_size = strlen(input_string);
  input_seq->current_pos = 0;
  return input_seq;
}

int
get_next_char_from_input (input_sequence *input_seq, unsigned long long *current_pos) 
{
  int c = 0;
  if (input_seq->sequence_file) {
    *current_pos = ftell(input_seq->sequence_file);
    c = fgetc (input_seq->sequence_file);
  } else if (input_seq->sequence_string && input_seq->input_size > 0) {
    if (input_seq->current_pos == input_seq->input_size) return -1;
    *current_pos = input_seq->current_pos;
    c = (int) input_seq->sequence_string[input_seq->current_pos];
    input_seq->current_pos += 1;
  }
  return c;
}

char *
get_header_name_from_input (input_sequence *input_seq, unsigned long long header_pos, unsigned long long current_pos, pr_append_str *parse_err)
{
  void *v = NULL;
  char *header_name = (char *) malloc (sizeof(char) * (current_pos - header_pos + 2));
  if (!header_name) {
    pr_append_new_chunk_external (parse_err, "Memory allocation for header name failed!");
    free(header_name);
    return NULL;
  }
  if (input_seq->sequence_file) {
    fseek (input_seq->sequence_file, header_pos, SEEK_SET);
    v = fgets (header_name, current_pos - header_pos + 2, input_seq->sequence_file);
  } else if (input_seq->sequence_string && input_seq->input_size > 0) {
    v = memcpy (header_name, input_seq->sequence_string + header_pos, current_pos - header_pos + 1);
  }
  if (!v) {
    pr_append_new_chunk_external (parse_err, "Reading header name failed!");
    free(header_name);
    return NULL;
  }
  return header_name;
}

void
delete_input_sequence (input_sequence *input_seq)
{
  if (!input_seq) return;
  if (input_seq->sequence_file && input_seq->sequence_file != stdin) {
    fclose (input_seq->sequence_file);
  }
  if (input_seq) free ((void *) input_seq);
  return;
}

/*==================================================================================
 * 
 * Formula parameters functions:
 * 
 *==================================================================================*/

formula_parameters *
create_formula_parameters_from_list_file_name (const char *list_file_name, pr_append_str *parse_err)
{
  const char *data;
  size_t size;
  unsigned long long header_size;
  unsigned int magic;
  formula_parameters *fp = (formula_parameters *) malloc (sizeof (formula_parameters));
  if (!fp) {
    pr_append_new_chunk_external (parse_err, "Memory allocation for formula parameters failed!");
    return fp;
  }
  memset (fp, 0, sizeof(formula_parameters));
  strcpy(fp->list_file_name, list_file_name);
  data = mmap_by_filename (fp->list_file_name, &size);
  if (!data) {
    pr_append_new_chunk_external (parse_err, "List file not found: ");
    pr_append_external (parse_err, fp->list_file_name);
    pr_append_external (parse_err, ". Lists can be specified by names or prefixes from the commandline or text file.");
    return NULL;
  }
  memcpy (&magic, data, sizeof(unsigned int));
  if (magic != glistmaker_code_match) {
    pr_append_new_chunk_external (parse_err, "Given file is not a list file: ");
    pr_append_external (parse_err, fp->list_file_name);
    return NULL;
  }
  
  memcpy (&fp->oligo_length, data + 12, sizeof(unsigned int));
  memcpy (&fp->words_in_list, data + 16, sizeof(unsigned int));
  memcpy (&header_size, data + 32, sizeof(unsigned long long));
  if (!fp->words_in_list){
     pr_append_new_chunk_external (parse_err, "List file contains no kmers: ");
           pr_append_external (parse_err, fp->list_file_name);
     return NULL;
  }
  fp->word_list = data + header_size;
        fp->pointer = data;
  fp->size = size;
  fp->binary_mask = create_binary_mask (fp->oligo_length);
  return fp;
}

formula_parameters *
create_formula_parameters_from_list_file_prefix (const char *list_name_prefix, const char *kmer_lists_path, unsigned int word_length, pr_append_str *parse_err)
{
  char list_file_name[300];
  formula_parameters *fp;
  sprintf(list_file_name, "%s%s_%u.list", kmer_lists_path, list_name_prefix, word_length);
  if(0 != access(list_file_name,0)){
    pr_append_new_chunk_external (parse_err, "Cannot find list file");
    return NULL;
  }
  fp = create_formula_parameters_from_list_file_name (list_file_name, parse_err);
  return fp;
}

int
add_variable_to_formula_parameters (char **list_values, unsigned int nvalues, parameters_builder *pbuilder, pr_append_str *parse_err)
{
  unsigned int i = 0;
  char *list_name = NULL;
  formula_parameters *fp = NULL;
  int add_parameters_to_existing_list = 0;
  unsigned int add_position;
  
  list_name = list_values[0];
  for (i = 0; i < pbuilder->nfp; i++) {
    if (!strcmp (list_name, pbuilder->used_lists[i])) {
      add_parameters_to_existing_list = 1;
      add_position = i;
      break;
    }
  }
  
  if (!add_parameters_to_existing_list) {
    fp = create_formula_parameters_from_list_file_name (list_name, parse_err);
    if (!fp) return 1;
    if (pbuilder->nfp >= pbuilder->nslots) {
      pbuilder->nslots = (pbuilder->nslots + 1) * 2;
      pbuilder->used_lists = (char **) realloc (pbuilder->used_lists, pbuilder->nslots * sizeof(char *));
      pbuilder->fp_array = (formula_parameters **) realloc (pbuilder->fp_array, pbuilder->nslots * sizeof(formula_parameters *));
      if (!pbuilder->used_lists || !pbuilder->fp_array) {
        pr_append_new_chunk_external (parse_err, "Memory allocation for parameters builder failed!");
        free(pbuilder->used_lists);
        free(pbuilder->fp_array);
        return 1;
      }
    }
    
    pbuilder->used_lists[pbuilder->nfp] = list_name;
    pbuilder->fp_array[pbuilder->nfp] = fp;
    add_position = pbuilder->nfp;
    pbuilder->nfp += 1;
  }
  
  if (fp || add_parameters_to_existing_list) {
    int squared = 0;
    unsigned int mm = 0;
    double coef = DEFAULT_COEF;
    char *end = NULL;
    if (nvalues > 1) {
      coef = (list_values[1][0] == '-') ? strtod (list_values[1] + 1, &end) * (-1.0) : strtod (list_values[1], &end);
      if (*end != 0) {
        pr_append_new_chunk_external (parse_err, "Invalid coefficient value: ");
        pr_append_external (parse_err, list_values[1]);
        return 2;
      }
    }
    if (nvalues > 2) {
      mm = strtol (list_values[2], &end, 10);
      if (*end != 0 || mm > 2) {
        pr_append_new_chunk_external (parse_err, "Invalid mismatches value specified: ");
        pr_append_external (parse_err, list_values[2]);
        pr_append_external (parse_err, ". Must be a positive integer less than 2.");
        return 3;
      }
    }
    if (nvalues > 3) {
      if (!strcmp(list_values[3], "sq")) squared = 1;
    }
    switch (mm) {
      case 0:
        if (squared) pbuilder->fp_array[add_position]->mm0_2 = coef;
        else pbuilder->fp_array[add_position]->mm0 = coef;
        break;
      case 1:
        if (squared) pbuilder->fp_array[add_position]->mm1_2 = coef;
        else pbuilder->fp_array[add_position]->mm1 = coef;
        break;
      case 2:
        if (squared) pbuilder->fp_array[add_position]->mm2_2 = coef;
        else pbuilder->fp_array[add_position]->mm2 = coef;
        break;
    }
    
  }  
  return 0;
}

formula_parameters **
create_default_formula_parameters (const char *list_name_prefix, const char *kmer_lists_path, pr_append_str *parse_err)
{
  formula_parameters *fp1 = create_formula_parameters_from_list_file_prefix (list_name_prefix, kmer_lists_path, DEFAULT_WORD_LEN_1, parse_err);
  formula_parameters *fp2 = create_formula_parameters_from_list_file_prefix (list_name_prefix, kmer_lists_path, DEFAULT_WORD_LEN_2, parse_err);
  if (!fp1 || !fp2) return NULL; 
  formula_parameters **fp = (formula_parameters **) malloc (2 * sizeof (formula_parameters *));
  if (!fp) {
    pr_append_new_chunk_external (parse_err, "Memory allocation for formula parameters failed!");
    return fp;
  }
  
  fp[0] = fp1;
  fp[1] = fp2;
    
  fp1->mm0 = DEFAULT_COEF_1;
  fp2->mm0 = DEFAULT_COEF_2;
  
  return fp;
}

formula_parameters ** 
read_formula_parameters_from_file (const char *lists_file_name, unsigned int *nlist_parameters, parameters_builder *pbuilder, double *intercept, pr_append_str *parse_err)
{
  FILE *lists = fopen (lists_file_name, "r");
  char *line = NULL;
  size_t line_length = 0;
  int read_chars;
  int v;
  
  if (!lists) {
    pr_append_new_chunk_external (parse_err, "File not found: ");
    pr_append_external (parse_err, lists_file_name);
    return NULL;
  }

  while ((read_chars = (int) getline(&line, &line_length, lists)) > 1) {
    char **values = NULL;
    unsigned int nvalues = 0;
    
    line[read_chars] = '\0';
    strip_string(line);
    values = split_string(line, ' ', &nvalues);
    if (nvalues == 1) {
      double ic;
      double neg = 1.0;
      char *end = NULL;
      if (values[0][0] == '-') {
        values[0] += 1;
        neg = -1.0;
      }
      ic = strtod (values[0], &end);
      if (*end == 0) {
        *intercept = ic * neg;
        continue;
      }
    }
    
    v = add_variable_to_formula_parameters (values, nvalues, pbuilder, parse_err);
    if (v) {
        free(pbuilder->used_lists);
        free(pbuilder->fp_array);
        return NULL;
    }
    *nlist_parameters += 1;
  }
  return pbuilder->fp_array;
}

void
delete_formula_parameters (formula_parameters **fp, unsigned int nlists)
{
  unsigned int i;
  if (!fp) return;
  for (i = 0; i < nlists; i++) {
    if (fp[i]->pointer) munmap ((void *) fp[i]->pointer, fp[i]->size);
    if (fp[i]) free ((void *) fp[i]);
  }
  if (fp) free ((void *) fp);
  return;
}

/*==================================================================================
 * 
 * Output function:
 * 
 *==================================================================================*/


output_sequence *
create_output_sequence (unsigned long long seq_len, masking_direction mdir, pr_append_str *parse_err) 
{
  output_sequence *output_seq = (output_sequence *) malloc (sizeof (output_sequence));
  memset (output_seq, 0, sizeof (output_sequence));
  if (!output_seq) {
    pr_append_new_chunk_external (parse_err, "Memory allocation for output sequence failed!");
    return output_seq;
  }
  if (mdir == both_separately) {
    output_seq->sequence_fwd = (char *) malloc (seq_len + 1);
    memset (output_seq->sequence_fwd, 0, seq_len + 1);
    output_seq->sequence_rev = (char *) malloc (seq_len + 1);
    memset (output_seq->sequence_rev, 0, seq_len + 1);
  } else {
    output_seq->sequence = (char *) malloc (seq_len + 1);
    memset (output_seq->sequence, 0, seq_len + 1);
  }
  if (!output_seq->sequence_fwd && !output_seq->sequence_rev && !output_seq->sequence) {
    pr_append_new_chunk_external (parse_err, "Memory allocation for output sequence failed!");
    return NULL;    
  }
  output_seq->pos = 0;
  return output_seq;
}

void 
delete_output_sequence (output_sequence *output_seq) 
{
  if (!output_seq) return;
  if (output_seq->sequence) free ((void *) output_seq->sequence);
  if (output_seq->sequence_fwd) free ((void *) output_seq->sequence_fwd);
  if (output_seq->sequence_rev) free ((void *) output_seq->sequence_rev);
  if (output_seq) free ((void* ) output_seq);
  return;
}

void 
write_header_to_output (output_sequence *output_seq, char *header_name, const masker_parameters *mp, pr_append_str *parse_err)
{
  void *v = NULL;
  if (mp->print_sequence) 
    fprintf (stdout, "%s", header_name);
  else if (output_seq) {
    if (mp->mdir == both_separately) {
      v = memcpy (output_seq->sequence_fwd + output_seq->pos, header_name, strlen(header_name));
      if (v) v = memcpy (output_seq->sequence_rev + output_seq->pos, header_name, strlen(header_name));
    } else {
      v = memcpy (output_seq->sequence + output_seq->pos, header_name, strlen(header_name));
    }  
    if (!v) {
      pr_append_new_chunk_external (parse_err, "Writing header to output failed!");
      return;
    }
    output_seq->pos += strlen(header_name);
  }
  return;
}

void 
write_char_to_output (output_sequence *output_seq, char c, char c_other, const masker_parameters *mp, pr_append_str *parse_err)
{
  if (mp->print_sequence) {
    fprintf (stdout, "%c", c);
  } else if (output_seq) {
    if (mp->mdir == both_separately) {
      output_seq->sequence_fwd[output_seq->pos] = c;
      output_seq->sequence_rev[output_seq->pos] = c_other;
    } else {
      output_seq->sequence[output_seq->pos] = c;
    }
    output_seq->pos += 1;
  }
  return;
}



/*==================================================================================
 * 
 * Masking buffer functions:
 * 
 *==================================================================================*/

masking_buffer *
create_masking_buffer (unsigned int word_length, pr_append_str *parse_err)
{
  masking_buffer *mbuffer = (masking_buffer *) malloc (sizeof(masking_buffer));
  if (!mbuffer) {
    pr_append_new_chunk_external (parse_err, "Memory allocation for masking buffer failed!");
    return mbuffer;
  }  
  initialize_masking_buffer (mbuffer, word_length);
  return mbuffer;
}

void 
initialize_masking_buffer (masking_buffer *mbuffer, unsigned int word_length)
{
  memset (mbuffer, 0, sizeof (masking_buffer));
  mbuffer->ei = MAX_BUFFER_SIZE - word_length + 1;
  return;
}

void 
delete_masking_buffer (masking_buffer *mbuffer)
{
  if (mbuffer) {
    free ((void *) mbuffer);
  }
  return;
}

void
add_char_to_buffer (char c, masking_buffer *mbuffer, int char_type)
{
  mbuffer->buffer[mbuffer->wi] = c;
  mbuffer->mask_positions_fwd[mbuffer->wi] = 0;
  mbuffer->mask_positions_rev[mbuffer->wi] = 0;
  mbuffer->non_nucleotide_positions[mbuffer->wi] = 0;
  
  if (char_type != WHITESPACE) {
    if (mbuffer->mi > 0) {
      mbuffer->mask_positions_fwd[mbuffer->wi] = 1;
      mbuffer->mi -= 1;
    } else if (char_type == MASKED_CHAR) {
      mbuffer->mask_positions_rev[mbuffer->wi] = 1;
      mbuffer->mask_positions_fwd[mbuffer->wi] = 1;
    }
    while (mbuffer->non_nucleotide_positions[mbuffer->ei] && !(mbuffer->mask_positions_fwd[mbuffer->ei])) {
      mbuffer->ei = TAKE_STEP_FORWARD(mbuffer->ei);
    }
    /* when it finds the first non-whitespace character it increases the ei one more time */
    mbuffer->ei = TAKE_STEP_FORWARD(mbuffer->ei);
  } 
  if (char_type == WHITESPACE || char_type == MASKED_CHAR) {
    mbuffer->non_nucleotide_positions[mbuffer->wi] = 1;
  }
  
  mbuffer->wi = TAKE_STEP_FORWARD(mbuffer->wi);
  return;
}

/* next character is not nucleotide, or if nucleotide then already masked */
#define COND0(i) (mbuffer->non_nucleotide_positions[i] || \
  (mp->do_soft_masking && mbuffer->buffer[i] >= 'a'))

/* output contains only one sequence and the next character is masked or
 * output contains two sequences and the next character is masked on the forward strand */
#define COND1(i) ((mp->mdir != both_separately && (mbuffer->mask_positions_fwd[i] || \
  mbuffer->mask_positions_rev[i])) || (mp->mdir == both_separately && mbuffer->mask_positions_fwd[i]))
  
/* output contains two sequences and the next character is masked on the reverse strand */
#define COND2(i) (mp->mdir == both_separately && mbuffer->mask_positions_rev[i])
  
void
empty_buffer (output_sequence *output_seq, const masker_parameters *mp, masking_buffer *mbuffer, int flush_all, pr_append_str *parse_err)
{
  unsigned int end = mbuffer->ei;
  if (flush_all) end = mbuffer->wi;
  while (mbuffer->ri != end) {
    if (COND0(mbuffer->ri)) {
      write_char_to_output (output_seq, mbuffer->buffer[mbuffer->ri], mbuffer->buffer[mbuffer->ri], mp, parse_err);
    } else {
      if (mp->do_soft_masking) {
        write_char_to_output (output_seq, COND1(mbuffer->ri) ? mbuffer->buffer[mbuffer->ri] + 32 : mbuffer->buffer[mbuffer->ri], 
                  COND2(mbuffer->ri) ? mbuffer->buffer[mbuffer->ri] + 32 : mbuffer->buffer[mbuffer->ri], mp, parse_err);
      } else {
        write_char_to_output (output_seq, COND1(mbuffer->ri) ? mp->masking_char : mbuffer->buffer[mbuffer->ri], 
                  COND2(mbuffer->ri) ? mp->masking_char : mbuffer->buffer[mbuffer->ri], mp, parse_err);        
      }
    }
    mbuffer->ri = TAKE_STEP_FORWARD(mbuffer->ri);
  }
  return;
}

#undef COND0
#undef COND1
#undef COND2

/*==================================================================================
 * 
 * K-mer searching functions:
 * 
 *==================================================================================*/


unsigned int 
binary_search (formula_parameters *fp, unsigned long long word)
{
  unsigned long long current_word, low, high, mid;
  unsigned int freq;
  low = 0;
  high = fp->words_in_list - 1;
  mid = (low + high) / 2;
  while (low <= high) {
    current_word = *((unsigned long long *) (fp->word_list + mid * (sizeof (unsigned long long) + sizeof (unsigned int))));
    if (current_word < word) {
      low = mid + 1;
    } else if (current_word > word) {
      if (mid == 0) break;
      high = mid - 1;
    } else {
      freq = *((unsigned int *) (fp->word_list + mid * (sizeof (unsigned long long) + sizeof (unsigned int)) + sizeof (unsigned long long)));
      return freq;
    }
    mid = (low + high) / 2;
  }
  return 0;
}


unsigned int 
get_frequency_of_canonical_oligo (formula_parameters *fp, unsigned long long word)
{
  unsigned int freq_fwd = 0, freq_rev = 0;
  freq_fwd = binary_search (fp, word);
  if (!freq_fwd) {
    freq_rev = binary_search (fp, get_reverse_complement (word, fp->oligo_length));
    
    //fprintf (stderr, "rev %u\n", freq_rev);
    if(freq_rev==0){ /*heuristics is used, by default, for speed and memory issues, 
          kmer lists are generated with kmers freq >1*/
       freq_rev=1;
    }
    return freq_rev;
  }
  
  //fprintf (stderr, "fwd %u\n", freq_fwd);
  if(freq_fwd==0){ /*heuristics is used, by default, for speed and memory issues,
               kmer lists are generated with kmers freq >1*/
         freq_fwd=1;
        }
  return freq_fwd;
}

void
get_oligo_frequencies (oligo_counts *oc, formula_parameters *fp, unsigned long long word, unsigned int mm, int strand)
{
  unsigned int i, j;
  unsigned int count_0mm = 0;
  unsigned int count_1mm = 0;
  unsigned int count_2mm = 0;
  unsigned int mismatch;
  
  word &= fp->binary_mask;
  count_0mm = get_frequency_of_canonical_oligo (fp, word);
  
  if (mm > 0) {
    for (i = 0; i < fp->oligo_length; i++) {
      for (mismatch = 1; mismatch < 4; mismatch++) {
        unsigned long long mask = mismatch << (2 * i);
        count_1mm += get_frequency_of_canonical_oligo (fp, word ^ mask);
        if (mm > 1) {
          for (j = i + 1; j < fp->oligo_length; j++) {
            unsigned long long mask2 = mismatch << (2 * j);
            count_2mm += get_frequency_of_canonical_oligo (fp, word ^ mask ^ mask2);
          }
        }
      }
    }
    
  }
  count_1mm += count_0mm;
  count_2mm += count_1mm;
  
  if (strand != REV) {
    oc->count_mm0_fwd = count_0mm;
    oc->count_mm1_fwd = count_1mm;
    oc->count_mm2_fwd = count_2mm;
  }
  if (strand != FWD) {
    oc->count_mm0_rev = count_0mm;
    oc->count_mm1_rev = count_1mm;
    oc->count_mm2_rev = count_2mm;    
  }
  return;
}




/*==================================================================================
 * 
 * Masking functions:
 * 
 *==================================================================================*/

void
calculate_scores (oligo_pair *h, const masker_parameters *mp, unsigned int word_length)
{
  formula_parameters **fp_array = mp->fp; // fp[0] on pointer!!
  unsigned int nlists = mp->nlists;
  unsigned int i;
  
  for (i = 0; i < nlists; i++) {
    oligo_counts oc = {0};
    formula_parameters *fp = fp_array[i];
    unsigned int mm = (fp->mm2 || fp->mm2_2) ? 2 : ((fp->mm1 || fp->mm1_2) ? 1 : 0);
    
    
    if ((mp->mdir == both_on_same || mp->mdir == both_separately || mp->abs_cutoff) && word_length == fp->oligo_length) {
      double score = 0.0;
      unsigned int abs_score = 0;
      
      if (mp->mdir == rev) get_oligo_frequencies (&oc, fp, h->rev, mm, BOTH);
      else get_oligo_frequencies (&oc, fp, h->fwd, mm, BOTH);
      
      if (oc.count_mm0_fwd || oc.count_mm0_rev) {
        unsigned int count = oc.count_mm0_fwd ? oc.count_mm0_fwd : oc.count_mm0_rev;
        score += fp->mm0 * log(oc.count_mm0_fwd) + fp->mm0_2 * log(oc.count_mm0_fwd) * log(oc.count_mm0_fwd);
        abs_score = count;
      }
      if (oc.count_mm1_fwd) {
        score += fp->mm1 * log(oc.count_mm1_fwd) + fp->mm1_2 * log(oc.count_mm1_fwd) * log(oc.count_mm1_fwd);
        abs_score = oc.count_mm1_fwd;
      }  
      if (oc.count_mm2_fwd) {
        score += fp->mm2 * log(oc.count_mm2_fwd) + fp->mm2_2 * log(oc.count_mm2_fwd) * log(oc.count_mm2_fwd);
        abs_score = oc.count_mm2_fwd;
      }
      if (mp->abs_cutoff) {
        h->abs_score = abs_score;
      } else {
        h->score_fwd += score;
        h->score_rev += score;
      }
    } else {
      if (mp->mdir != rev) {
        get_oligo_frequencies (&oc, fp, h->fwd, mm, FWD);
        
        //fprintf (stderr, "oc.count.fwd %u\n", oc.count_mm0_fwd);
        
        if (oc.count_mm0_fwd) h->score_fwd += fp->mm0 * log(oc.count_mm0_fwd) + fp->mm0_2 * log(oc.count_mm0_fwd) * log(oc.count_mm0_fwd);
        if (oc.count_mm1_fwd) h->score_fwd += fp->mm1 * log(oc.count_mm1_fwd) + fp->mm1_2 * log(oc.count_mm1_fwd) * log(oc.count_mm1_fwd);
        if (oc.count_mm2_fwd) h->score_fwd += fp->mm2 * log(oc.count_mm2_fwd) + fp->mm2_2 * log(oc.count_mm2_fwd) * log(oc.count_mm2_fwd);
        
        //fprintf (stderr, "1: sõnad: %s %s, skoorid %f %f, abs_skoor %u\n", word_to_string(h->fwd, word_length), word_to_string(h->rev, word_length), h->score_fwd, h->score_rev, h->abs_score);
      }
      if (mp->mdir != fwd) {
        get_oligo_frequencies (&oc, fp, h->rev, mm, REV);
        if (oc.count_mm0_rev) h->score_rev += fp->mm0 * log(oc.count_mm0_rev) + fp->mm0_2 * log(oc.count_mm0_rev) * log(oc.count_mm0_rev);
        if (oc.count_mm1_rev) h->score_rev += fp->mm1 * log(oc.count_mm1_rev) + fp->mm1_2 * log(oc.count_mm1_rev) * log(oc.count_mm1_rev);
        if (oc.count_mm2_rev) h->score_rev += fp->mm2 * log(oc.count_mm2_rev) + fp->mm2_2 * log(oc.count_mm2_rev) * log(oc.count_mm2_rev);
      }
    }
  }
  if (h->score_fwd) h->score_fwd = exp(h->score_fwd + mp->formula_intercept) / (1 + exp(h->score_fwd + mp->formula_intercept));
  if (h->score_rev) h->score_rev = exp(h->score_rev + mp->formula_intercept) / (1 + exp(h->score_rev + mp->formula_intercept));
  
  //fprintf (stderr, "2: sõnad: %s %s, skoorid %f %f, abs_skoor %u\n", word_to_string(h->fwd, word_length), word_to_string(h->rev, word_length), h->score_fwd, h->score_rev, h->abs_score);
  
  return;
}

void 
mask_oligo_region (oligo_pair *h, const masker_parameters *mp, masking_buffer *mbuffer, unsigned int word_length, int debug)
{
  calculate_scores (h, mp, word_length);
  
  if (debug > 1) {
    fprintf (stderr, "score-fwd: %f score-rev: %f\n", h->score_fwd, h->score_rev);
  }
  
  if (mp->mdir != rev && ((mp->failure_rate && h->score_fwd > mp->failure_rate) || 
    (mp->abs_cutoff && h->abs_score >= mp->abs_cutoff))) {
    int masked = 0, i = TAKE_STEP_BACK(mbuffer->wi);
    while (masked < mp->nucl_masked_in_5p_direction) {
      if (!mbuffer->non_nucleotide_positions[i] && !mbuffer->mask_positions_fwd[i]) {
        mbuffer->mask_positions_fwd[i] = 1;
        masked += 1;
      } else if (mbuffer->mask_positions_fwd[i]) {
        masked += 1;
      }
      i = TAKE_STEP_BACK(i);
    }
    mbuffer->mi = mp->nucl_masked_in_3p_direction;
  }
  
  if (mp->mdir != fwd && ((mp->failure_rate && h->score_rev > mp->failure_rate) || 
    (mp->abs_cutoff && h->abs_score >= mp->abs_cutoff))) {
    int masked = 0, i = TAKE_STEP_BACK(mbuffer->ei);    
    while (masked < mp->nucl_masked_in_5p_direction + mp->nucl_masked_in_3p_direction) {
      if (!mbuffer->non_nucleotide_positions[i] && !mbuffer->mask_positions_rev[i]) {
        mbuffer->mask_positions_rev[i] = 1;
        masked += 1;
      } else if (mbuffer->mask_positions_rev[i]) {
        masked += 1;
      }
      i = TAKE_STEP_FORWARD(i);
    }
    
  }
  return;
}

void
read_and_mask_sequence (input_sequence *input_seq, output_sequence *output_seq, const masker_parameters *mp, 
      pr_append_str *parse_err, int debug)
{
  int is_header = 0, init_round = 1;
  unsigned long long word_fwd = 0, word_rev = 0, nucl_value;
  unsigned int current_length = 0;
  unsigned int word_length = 0;
  unsigned long long binary_mask = 0;
  masking_buffer *mbuffer;
  unsigned long long header_pos = 0, current_pos = 0;
  unsigned int i;
  
  /* size of window is the length of the longest k-mers that we use */
  for (i = 0; i < mp->nlists; i++) {
    if (mp->fp[i]->oligo_length > word_length) {
      binary_mask = mp->fp[i]->binary_mask;
      word_length = mp->fp[i]->oligo_length;
    }
  }
  
  /* buffer for masking*/
  mbuffer = create_masking_buffer (word_length + mp->nucl_masked_in_3p_direction, parse_err);
  
  while (1) {
    oligo_pair h = {0};
    int c = get_next_char_from_input (input_seq, &current_pos);
    /* EOF or end of string */
    if (c < 0) break;
    
    if (debug > 1) {
      fprintf (stderr, "pos: %llu, input: %c\n", current_pos, c);
    }
    
    /* In case it is a FASTA file we need to consider 
     * the headers separately */
    if (c == '>') {
      header_pos = current_pos;
      word_fwd = word_rev = 0;
      current_length = 0;
      is_header = 1;
    }
    /* continue until the end of header */
    if (is_header) {
      if (c == 10 || c == 13) {
        char *header_name = get_header_name_from_input (input_seq, header_pos, current_pos, parse_err);
        empty_buffer (output_seq, mp, mbuffer, FLUSH_ALL, parse_err);
        write_header_to_output (output_seq, header_name, mp, parse_err);
        initialize_masking_buffer (mbuffer, word_length + mp->nucl_masked_in_3p_direction);
        init_round = 1;
        is_header = 0;
        free(header_name);
      }
    } else {
      if (!init_round && mbuffer->wi == mbuffer->ri) {
        empty_buffer (output_seq, mp, mbuffer, KEEP_UNCERTAIN_REGION, parse_err);
      }
      init_round = 0;
      
      if (!strchr(ALPHABET, c) && c > ' ') {
        add_char_to_buffer (c, mbuffer, MASKED_CHAR);
        word_fwd = word_rev = 0;
        current_length = 0;
        continue;
      } else if (c <= ' ') {
        add_char_to_buffer (c, mbuffer, WHITESPACE);
        continue;
      } else {
        add_char_to_buffer (c, mbuffer, NUCLEOTIDE);
      }
      
      /* add next nucleotide to the word */
      nucl_value = get_nucl_value (c);
      if (mp->mdir != rev) {
        word_fwd <<= 2;
        word_fwd |= nucl_value;
      }
      if (mp->mdir != fwd) {
        word_rev >>= 2;
        word_rev |= ((~nucl_value & 3) << ((word_length - 1) * 2));
      }
      current_length += 1;
      
      if (current_length > word_length) {
        word_fwd &= binary_mask;
        word_rev &= binary_mask;
        current_length = word_length;
      }
      if (current_length == word_length) {
        h.fwd = word_fwd;
        h.rev = word_rev;
        
        if (debug > 1) {
          fprintf (stderr, "%llu %llu\n", h.fwd, h.rev);
        }
        mask_oligo_region (&h, mp, mbuffer, word_length, debug);
      }
    }
  }
  
  empty_buffer (output_seq, mp, mbuffer, FLUSH_ALL, parse_err);
  delete_masking_buffer (mbuffer);
  return;
} 

/*==================================================================================
 * 
 * Helping functions:
 * 
 *==================================================================================*/

const char * 
mmap_by_filename (const char *filename, size_t *size)
{
  struct stat st;
  int status, handle;
  const char *data;

  status = stat (filename, &st);
  if (status < 0) {
    return NULL;
  }

  handle = open (filename, O_RDONLY);
  if (handle < 0) {
    return NULL;
  }

  data = (const char *) mmap (NULL, st.st_size, PROT_READ, MAP_PRIVATE, handle, 0);
  if (data == (const char *) -1) {
    return NULL;
  } else {
    *size = (size_t)st.st_size;
  }

  close (handle);
  return data;
}

unsigned long long 
create_binary_mask (unsigned int word_length)
{
  unsigned int i;
  unsigned long long mask = 0L;

  for (i = 0; i < 2 * word_length; i++) {
    mask = (mask << 1) | 1;
  }
  return mask;
}


unsigned long long 
get_nucl_value (char nucl)
{
  static unsigned long long bit1 = 1 << 2;
  static unsigned long long bit2 = 3 << 1;
  if (nucl & bit1) {
    return ((nucl >> 4) | 2) & 3;
  }
  return (nucl & bit2) >> 1;
}

unsigned long long 
get_reverse_complement (unsigned long long word, unsigned int word_length)
{
  unsigned int i;
  unsigned long long mask, v, revcompl = 0L;

  word = ~word;
  mask = 3;
  for (i = 0; i < word_length; i++) {
    v = word & mask;
    revcompl <<= 2;
    revcompl |= v;
    word >>= 2;
  }
  return revcompl;
}

unsigned long long 
string_to_word (const char *s, unsigned int string_length, unsigned int word_length)
{
  unsigned int i;
  unsigned long long word = 0L;

  for (i = string_length - word_length; i < string_length; i++) {
    word <<= 2;
    word |= get_nucl_value (s[i]);
  }
  return word;
}

char * 
word_to_string (unsigned long long word, unsigned int wordlength)
{
  char *s = (char *) malloc (wordlength + 1);
  unsigned int i, temp;

  for (i = 0; i < wordlength; i++) {
    temp = word & 3;
    s[wordlength - i - 1] = ALPHABET[temp];
    word >>= 2;
  }
  s[wordlength] = 0;
  return s;
}

char **
split_string (char string[], char c, unsigned int *nchunks)
{
  char **string_chunks  = (char **) malloc (MAX_SPLITS * sizeof(char *));
  char *p;
  char tmp[100];
  unsigned int i = 0, l;
  
  while ((p = strchr (string, c))) {
    /* length of substring to copy */
    l = (int) (p - string);
    
    /* add substring to the output */
    if (l > 0) {
      memcpy (tmp, string, l);
      tmp[l] = '\0';
      string_chunks[i] = (char *) malloc ((l + 1) * sizeof(char));
      strcpy (string_chunks[i], tmp);
      i += 1;
      *nchunks += 1;
    }
    string = p + 1;
  }
  /* length of the last substring */
  l = strlen(string);
  
  if (l > 0) {
  
    memcpy (tmp, string, l);
    tmp[l] = '\0';
    /* add last substring to output */
    string_chunks[i] = (char *) malloc ((l + 1) * sizeof(char));
    strcpy (string_chunks[i], tmp);
    *nchunks += 1;
  }
  
  return string_chunks;
}


void
strip_string (char string[])
{
  size_t l = strlen(string);
  if (string[l - 1] == '\n') string[l - 1] = '\0';
  return;
  
}


