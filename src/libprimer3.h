/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008
Whitehead Institute for Biomedical Research, Steve Rozen
(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky
All rights reserved.

    This file is part of primer3 and the libprimer3 library.

    Primer3 and the libprimer3 library are free software;
    you can redistribute them and/or modify them under the terms
    of the GNU General Public License as published by the Free
    Software Foundation; either version 2 of the License, or (at
    your option) any later version.

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this file (file gpl-2.0.txt in the source
    distribution); if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

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

#ifndef LIBPRIMER3_H
#define LIBPRIMER3_H

#include <setjmp.h>
#include <stdio.h> /* FILE */
#include <stdlib.h>
#include <limits.h> /* SHRT_MIN, ULONG_MAX */

#include "oligotm.h"

#ifdef __cplusplus
  extern "C" {
#endif

/* ALIGN_SCORE_UNDEF is used only libprimer3 and clients, not in dpal */
#define ALIGN_SCORE_UNDEF             SHRT_MIN

/* These next 5 are exposed for format_output.c -- probabaly should be reviewed. */
#define PR_INFINITE_POSITION_PENALTY -1.0
#define PR_DEFAULT_INSIDE_PENALTY     PR_INFINITE_POSITION_PENALTY
#define PR_DEFAULT_OUTSIDE_PENALTY    0.0
#define PR_DEFAULT_PRODUCT_MAX_TM     1000000.0
#define PR_DEFAULT_PRODUCT_MIN_TM     -1000000.0

/*  Exposed in the read_boulder input routine.... */
#define PR_NULL_START_CODON_POS       -1000000
#define PR_DEFAULT_START_CODON_POS    PR_NULL_START_CODON_POS
#define PR_START_CODON_POS_IS_NULL(SA) ((SA)->start_codon_pos <= PR_NULL_START_CODON_POS)

#define _PR_DEFAULT_POSITION_PENALTIES(PA) \
    (PR_DEFAULT_INSIDE_PENALTY == pa->inside_penalty \
     && PR_DEFAULT_OUTSIDE_PENALTY == pa->outside_penalty)

#define PR_ALIGN_SCORE_PRECISION 100.0

#define MACRO_STRING(X) #X
/* pr_progam_name must be set in main(). */
#define PR_ASSERT(COND)                                  \
if (!(COND)) {                                           \
    fprintf(stderr, "%s:%s:%d, assertion (%s) failed\n", \
           pr_program_name, __FILE__, __LINE__,          \
           MACRO_STRING(COND));                          \
    abort();                                             \
}

/* ANDREAS, we need documentation on these */
/* Enum to define tasks primer3 can do */
typedef enum task { 
  pick_pcr_primers               = 0,
  /*  For backward compatibility, equivalent to pick_detection_primers */
  pick_pcr_primers_and_hyb_probe = 1,
  pick_left_only                 = 2,
  pick_right_only                = 3,
  pick_hyb_probe_only            = 4,
  pick_detection_primers         = 5,
  pick_cloning_primers           = 6,
  pick_discriminative_primers    = 7,    
  pick_sequencing_primers        = 8,
  pick_primer_list               = 9,
  check_primers                  = 10,
} task;

/* Enum explaining if output are pairs */
typedef enum p3_output_type {
  primer_pairs    = 0,
  primer_list     = 1,
} p3_output_type;

/* pr_append_str is an append-only string ADT. */
typedef struct pr_append_str {
  int storage_size;
  char *data;
} pr_append_str;

/* 
 * Arguments to the primer program as a whole.  Values for these arguments are
 * retained _across_ different input records.  (These are the so-called
 * "Global" arguments in the documentation.)
 */
typedef struct oligo_weights {
  double temp_gt;
  double temp_lt;
  double gc_content_gt;
  double gc_content_lt;
  double compl_any;
  double compl_end;
  double repeat_sim;
  double length_lt;
  double length_gt;
  double seq_quality;
  double end_quality;
  double pos_penalty;
  double end_stability;
  double num_ns;
  double template_mispriming;
} oligo_weights;

typedef struct pair_weights {
  double primer_quality;
  double io_quality;
  double diff_tm;
  double compl_any;
  double compl_end;
  double product_tm_lt;
  double product_tm_gt;
  double product_size_lt;
  double product_size_gt;
  double repeat_sim;
  double template_mispriming;
} pair_weights;

typedef struct sequencing_parameters {
  int lead;
  int spacing;
  int interval;
  int accuracy;
} sequencing_parameters;

#include "p3_seq_lib.h"

typedef struct args_for_one_oligo_or_primer {
  seq_lib       *repeat_lib;
  oligo_weights weights;
  double opt_tm;
  double min_tm;
  double max_tm;
  double opt_gc_content;
  double max_gc;
  double min_gc;

  /* Warning: also used for product Tm (TO DO: factor this out) */
  double salt_conc;

  /*
    DIVALENT_CONC and DNTP_CONC are both needed for enabling use of
    divalent cations for calculation of melting temperature of short
    and long oligos.  The formula for converting the divalent cations
    to monovalent cations is in the paper [Ahsen von N, Wittwer CT,
    Schutz E (2001) "Oligonucleotide Melting Temperatures under PCR
    Conditions: Nearest-Neighbor Corrections for Mg^2+,
    Deoxynucleotide Triphosphate, and Dimethyl Sulfoxide
    Concentrations with Comparision to Alternative Empirical
    Formulas", Clinical Chemistry 47:1956-61
    http://www.clinchem.org/cgi/content/full/47/11/1956] The default
    is 0.0.  (New in v. 1.1.0, added by Maido Remm and Triinu
    Koressaar.)
  */
  double divalent_conc;
  double dntp_conc;

  double dna_conc;
  int    num_ns_accepted;
  int    opt_size;
  int    min_size;
  int    max_size;
  int    max_poly_x;      /* 
                           * Maximum length of mononucleotide sequence in an
                           * oligo.
                           */

  int    min_end_quality;
  int    min_quality;       /* Minimum quality permitted for oligo sequence.*/

  short  max_self_any;  
  short  max_self_end;
  short  max_repeat_compl;   /* 
                              * Acceptable complementarity with repeat
                              * sequences.
                              */

  short  max_template_mispriming;
} args_for_one_oligo_or_primer;

/* Maxima needed for interface data structures. */
#define PR_MAX_INTERVAL_ARRAY 200 
/* 
 * Maximum number of input intervals
 * supported; used for targets, excluded
 * regions, product-size intervals, etc.
 */

typedef struct p3_global_settings {
  /* ================================================== */
  /* Arguments that control behavior of choose_primers() */
  task   primer_task;          /* 2 if left primer only, 3 if right primer only,
                                * 4 if internal oligo only.    */

  int    pick_left_primer;
  int    pick_right_primer;
  int    pick_internal_oligo;

  int    file_flag;   /* TO DO, See if this can be factored out. */

  int    first_base_index;  /* 
                             * The index of the first base in the input
                             * sequence.  This parameter is ignored within
                             * pr_choice; pr_choice's caller must assure that
                             * all indexes are 0-based.  However, this
                             * parameter should used by output routines to
                             * adjust base indexes.
                             */

  int    liberal_base;   /* 
                          * If non-0 then turn characters other than
                          * [ATGCNatgcn] into N.
                          */

  int    num_return; /* The number of best primer pairs to return. */

  int    pick_anyway;    /* Pick even if input primer or oligos
                            violate constraints. */

  int    lib_ambiguity_codes_consensus;
  /* If non-0, treat ambiguity codes in a mispriming/mishyb
     library as representing a consensus.  So, for example,
     S would match C or G.  N would match any nucleotide.
     It turns out that this _not_ what one normally wants,
     since many libraries contain strings of N, which then
     match every oligo (very bad).
  */

  int    quality_range_min;
  int    quality_range_max;

  /* ================================================== */
  /* Arguments for individual oligos and/or primers */
  args_for_one_oligo_or_primer p_args;
  args_for_one_oligo_or_primer o_args;

  /* 
     Added by T.Koressaar for updated table thermodynamics.  Specifies
     details of melting temperature calculation.  (New in v. 1.1.0,
     added by Maido Remm and Triinu Koressaar.)
      
     A value of 1 (recommended) directs primer3 to use the table of
     thermodynamic values and the method for melting temperature
     calculation suggested in the paper [SantaLucia JR (1998) "A
     unified view of polymer, dumbbell and oligonucleotide DNA
     nearest-neighbor thermodynamics", Proc Natl Acad Sci 95:1460-65
     http://dx.doi.org/10.1073/pnas.95.4.1460].
      
     A value of 0 directs primer3 to a backward compatible calculation
     (in other words, the only calculation availble in previous
     version of primer3).
      
     This backward compatible calculation uses the table of
     thermodynamic parameters in the paper [Breslauer KJ, Frank R,
     Bloecker H and Marky LA (1986) "Predicting DNA duplex stability
     from the base sequence" Proc Natl Acad Sci 83:4746-50
     http://dx.doi.org/10.1073/pnas.83.11.3746],
     and the method in the paper [Rychlik W, Spencer WJ and Rhoads
     RE (1990) "Optimization of the annealing temperature for DNA
     amplification in vitro", Nucleic Acids Res 18:6409-12
     http://www.pubmedcentral.nih.gov/articlerender.fcgi?tool=pubmed&pubmedid=2243783].
       
     The default value is 0 only for backward compatibility.
  */
  tm_method_type tm_santalucia;  

  /* 
     Added by T.Koressaar for salt correction for Tm calculation.  A
     value of 1 (recommended) directs primer3 to use the salt correction
     formula in the paper [SantaLucia JR (1998) "A unified view of polymer,
     dumbbell and oligonucleotide DNA nearest-neighbor thermodynamics",
     Proc Natl Acad Sci 95:1460-65
     http://dx.doi.org/10.1073/pnas.95.4.1460]

     A value of 0 directs primer3 to use the the salt correction
     formula in the paper [Schildkraut, C, and Lifson, S (1965)
     "Dependence of the melting temperature of DNA on salt
     concentration", Biopolymers 3:195-208 (not available on-line)].
     This was the formula used in previous version of primer3.

     A value of 2 directs primer3 to use the salt correction formula
     in the paper [Owczarzy R, You Y, Moreira BG, Manthey JA, Huang L,
     Behlke MA and Walder JA (2004) "Effects of sodium ions on DNA
     duplex oligomers: Improved predictions of melting temperatures",
     Biochemistry 43:3537-54 http://dx.doi.org/10.1021/bi034621r].

     The default is 0 only for backward compatibility.
  */
  salt_correction_type salt_corrections; 

  /* Arguments applicable to primers but not oligos */

  double max_end_stability;
  /* The maximum value allowed for the delta
   * G of disruption for the 5 3' bases of
   * a primer.
   */

  int    max_end_gc;
  /* The maximum number of allowed G or C
   * for the 5 3' bases of a primer.
   */

  int    gc_clamp;              /* Required number of GCs at 3' end. */

  /* ================================================== */
  /* Arguments related to primer and/or oligo
     location in the template. */

  int lowercase_masking; 
  /* 
     Added by T.Koressaar. Enables design of primers from lowercase
     masked template.  A value of 1 directs primer3 to reject primers
     overlapping lowercase a base exactly at the 3' end.

     This property relies on the assumption that masked features
     (e.g. repeats) can partly overlap primer, but they cannot overlap the
     3'-end of the primer.  In other words, lowercase bases at other
     positions in the primer are accepted, assuming that the masked
     features do not influence the primer performance if they do not
     overlap the 3'-end of primer.
  */

  /* Parameters used to calculate the position of sequencing primers */
  sequencing_parameters sequencing;

  double outside_penalty; /* Multiply this value times the number of NTs
                           * from the 3' end to the the (unique) target to
                           * get the 'position penalty'.
                           * Meaningless if there are multiple targets
                           * or if the primer cannot be part of a pair
                           * that spans the target.
                           */

  double inside_penalty;  /* Multiply this value times the number of NT
                           * positions by which the primer overlaps
                           * the (unique) target to the 'position penalty'.
                           * Meaningless if there are multiple targets
                           * or if the primer cannot be part of a pair
                           * that spans the target.
                           */

  /* ================================================== */
  /* Arguments for primer pairs and products. */

  /* Warning,use p3_empty_gs_product_size_range and
     p3_add_to_gs_product_size_range to set this next
     three slot. */
  int    pr_min[PR_MAX_INTERVAL_ARRAY]; /* Minimum product sizes. */
  int    pr_max[PR_MAX_INTERVAL_ARRAY]; /* Maximum product sizes. */
  int    num_intervals;         /* 
                                 * Number of product size intervals
                                 * (i.e. number of elements in pr_min and
                                 * pr_max)
                                 */

  int    product_opt_size;
  double product_max_tm;
  double product_min_tm;
  double product_opt_tm;
  short  pair_max_template_mispriming;
  short  pair_repeat_compl;
  short  pair_compl_any;
  short  pair_compl_end;

  /* Max difference between temperature of primer and temperature of
     product.  Cannot be calculated until product is known. */
  double max_diff_tm; 

  pair_weights  pr_pair_weights;

  /* Minimum number of base pairs between the 3' ends of any two left
     or any two right primers when returning num_return primer pairs.
     The objective is get 'truly different' primer pairs.

     Primers that end at e.g. 30 and 31 have a three-prime distance
     of 1.

     0 indicates a primer pair is ok if it has not already appeared in
     the output list (default behavior and behavior in previous
     releases). This is the most liberal behavior.

     n > 0 indicates that a primer pair is ok if:

     NOT (3' end of left primer closer than n to the 3' end of a left
     primer in an existing pair)

     AND

     NOT (3' end of right primer closer than n
     to the 3' end of right primer in an existing pair)
  */
  int    min_three_prime_distance; 
  
  /* The number of basepairs the primer has to overlap an overlap position. */
  int    pos_overlap_primer_end; 
 
  char *settings_file_id;
  int dump;  /* dump fields for global settings and seq args if dump == 1 */
} p3_global_settings;

typedef enum oligo_type { OT_LEFT = 0, OT_RIGHT = 1, OT_INTL = 2 }
  oligo_type;

typedef struct rep_sim {
  char *name;      /* Name of the sequence from given file in fasta
                    * format with maximum similarity to the oligo.
                    */
  short min;       /* 
                    * The minimum score in slot 'score' (below).
                    * (Used when the objective function involves
                    * minimization of mispriming possibilities.)
                    */
  short max;       /* The maximum score in slot 'score' (below). */
  short *score;    /* 
                    * Array of similarity (i.e. false-priming) scores,
                    * one for each entry in the 'repeat_lib' slot
                    * of the primargs struct. 
                    */
} rep_sim;

#if (ULONG_MAX < 4294967295UL)
CANNOT COMPILE FOR THIS SYSTEM (< 32 bits in an unsigned long it)
#endif
typedef struct oligo_problems {
  unsigned long int prob;
} oligo_problems;

typedef struct primer_rec {

  rep_sim repeat_sim;
  /* Name of the sequence from given file in fasta
   * format with maximum similarity to the oligo
   * and corresponding alignment score.
   */

  double temp; /* The oligo melting temperature calculated for the
                * primer. */

  double gc_content;

  double position_penalty; 
  /*
   * Penalty for distance from "ideal" position as specified
   * by inside_penalty and outside_penalty.
   */

  double quality;  /* Part of objective function due to this primer. */

  double end_stability;
  /* Delta G of disription of 5 3' bases. */

  int    start;    /* Position of the 5'-most base within the primer
                      WITH RESPECT TO THE seq_args FIELD
                      trimmed_seq. */

  int    seq_quality; /* Minimum quality score of bases included. */   
  int    seq_end_quality;  /* Minimum quality core of the 5 3' bases. */


  short  self_any; /* Self complementarity as local alignment * 100. */

  short  self_end; /* Self complementarity at 3' end * 100. */

  short  template_mispriming;
  /* Max 3' complementarity to any ectopic site in template
     on the given template strand. */

  short  template_mispriming_r;
  /* Max 3' complementarity to any ectopic site in the
     template on the reverse complement of the given template
     strand. */

  char   length;   /* Length of the oligo. */
  char   num_ns;   /* Number of Ns in the oligo. */

  char   must_use; /* Non-0 if the oligo must be used even if it is illegal. */

  oligo_problems problems;

} primer_rec;

const char *
p3_primer_rec_problems_to_string(const primer_rec *);

int p3_ol_is_ok(const primer_rec *);

/* 
 * The structure for a pair of primers and for a pair of
 * primers plus an internal oligo.
 */
typedef struct primer_pair {
  double pair_quality;  /* Penalty value of the primer pair */

  double diff_tm;       /* Absolute value of the difference between melting
                         * temperatures for left and right primers. 
                         */
   
  double product_tm;    /* Estimated melting temperature of the product. */

  double product_tm_oligo_tm_diff;
                        /* Difference in Tm between the primer with lowest Tm
                           the product Tm. */

  double t_opt_a;

  int    compl_any;     /* 
                         * Local complementarity score between left and right
                         * primers (* 100).
                         */

  int    compl_end;     /* 
                         * 3'-anchored global complementatory score between *
                         * left and right primers (* 100).
                         */

  int    template_mispriming;
                        /* Maximum total mispriming score of both primers
                           to ectopic sites in the template, on "same"
                           strand (* 100). */

  short  repeat_sim;    /* Maximum total similarity of both primers to the
                         * sequence from given file in fasta format.
                         */
  primer_rec *left;     /* Left primer. */
  primer_rec *right;    /* Right primer. */
  primer_rec *intl;     /* Internal oligo. */

  int    product_size;  /* Product size. */
  int    target;        /* 
                         * 1 if there is a target between the right and left
                         * primers.
                         */
  char   *rep_name;
} primer_pair;

typedef int interval_array_t[PR_MAX_INTERVAL_ARRAY][2];

typedef struct interval_array_t2 {
  int pairs[PR_MAX_INTERVAL_ARRAY][2];
  int count;
} interval_array_t2;

int 
interval_array_t2_count(const interval_array_t2 *array);

const int *
interval_array_t2_get_pair(const interval_array_t2 *array, int i);

typedef struct oligo_stats {
  int considered;          /* Total number of tested oligos of given type   */
  int ns;                  /* Number of oligos rejected because of Ns       */
  int target;              /* Overlapping targets.                          */
  int excluded;            /* Overlapping excluded regions.                 */
  int gc;                  /* Unacceptable GC content.                      */
  int gc_clamp;            /* Don't have required number of GCs at 3' end.  */
  int gc_end_high;         /* Too many G+Cs at the 3' end.                  */
  int temp_min;            /* Melting temperature below t_min.              */
  int temp_max;            /* Melting temperature more than t_max.          */
  int size_min;            /* Primer shorter than minimal size.             */
  int size_max;            /* Primer longer than minimal size.              */
  int compl_any;           /* Self-complementarity too high.                */
  int compl_end;           /* Self-complementarity at 3' end too high.      */
  int repeat_score;        /* Complementarity with repeat sequence too high.*/
  int poly_x;              /* Long mononucleotide sequence inside.          */
  int seq_quality;         /* Low quality of bases included.                */
  int stability;           /* Stability of 5 3' bases too high.             */
  int no_orf;              /* Would not amplify any of the specified ORF
                             (valid for left primers only).                 */
  int template_mispriming; /* Template mispriming score too high.           */
  int ok;                  /* Number of acceptable oligos.                  */
   int gmasked;            /* edited by T. Koressaar, number of gmasked oligo*/
} oligo_stats;

typedef struct pair_stats {
  int considered;          /* Total number of pairs or triples tested.      */
  int product;             /* Pairs providing incorrect product size.       */
  int target;              /* Pairs without any target between primers.     */
  int temp_diff;           /* Melting temperature difference too high.      */
  int compl_any;           /* Pairwise complementarity larger than allowed. */
  int compl_end;           /* The same for 3' end complementarity.          */
  int internal;            /* Internal oligo was not found.                 */
  int repeat_sim;          /* Complementarity with repeat sequence too high.*/
  int high_tm;             /* Product Tm too high.                          */
  int low_tm;              /* Product Tm too low.                           */
  int template_mispriming; /* Sum of template mispriming scores too high.   */

  /* Neither oligo in the pairs overlaps one of the "required sites".       */
  int does_not_overlap_a_required_point;

  /* One of the oligos in the pair overlaps an oligo in a better_pair:       */
  int overlaps_oligo_in_better_pair;

  int ok;                  /* Number that were ok.                          */
} pair_stats;

typedef struct pair_array_t {
  int         storage_size;
  int         num_pairs;
  primer_pair *pairs;
  pair_stats  expl;
} pair_array_t;

/*
 * Arguments relating to a single particular source sequence (for which
 * we will pick primer(s), etc.
 */
typedef struct seq_args {

                          /* The net next 3 slots are presented as
                           * indexes within the sequence slot, but
                           * they are recalculated to be indexes
                           * within trimmed_seq (i.e. within the
			   * "included region").
                           */

  interval_array_t2 tar2; /* The targets.  tar2->pairs[i][0] s the start
                           * of the ith target, tar[i][1] its length.  */

  interval_array_t2 excl2;/* The number of excluded regions. */

  interval_array_t2 excl_internal2; 
                          /* Number of excluded regions for internal
                             oligo; similar to excl2.*/
  
  int primer_overlap_pos[PR_MAX_INTERVAL_ARRAY]; /* List of overlap positions. */
  int primer_overlap_pos_count; /* The number of intron sites. */

  int incl_s;             /* The 0-based start of included region. */
  int incl_l;             /* 
                           * The length of the included region, which is
                           * also the length of the trimmed_seq field.
                           */
  int  start_codon_pos;   /* Index of first base of the start codon. */

  int  *quality;             /* Vector of quality scores. */
  int  n_quality;            /* Number of valid elements in 'quality' */
  int  quality_storage_size; /* Amount of storage quality points to. */

  char *sequence;         /* The template sequence itself as input, 
                             not trimmed, not up-cased. */
  char *sequence_name;    /* An identifier for the sequence. */
  char *sequence_file;    /* Another identifer for the sequence. */
  char *trimmed_seq;      /* The included region only, _UPCASED_. */

  /* Element add by T. Koressaar support lowercase masking: */
  char *trimmed_orig_seq; /* Trimmed version of the original,
                             mixed-case sequence. */

  char *upcased_seq;      /* Upper case version of sequence
                             (_not_ trimmed). */

  char *upcased_seq_r;    /* Upper case version of sequence, 
                             other strand (_not_ trimmed). */

  char *left_input;       /* A left primer to check or design around. */

  char *right_input;      /* A right primer to check or design around. */

  char *internal_input;   /* An internal oligo to check or design around. */
  
  int force_left_start;   /* The 0-based forced 5' start left primer. */
  int force_left_end;     /* The 0-based forced 3' end left primer. */
  int force_right_start;  /* The 0-based forced 5' start right primer. */
  int force_right_end;    /* The 0-based forced 3' end right primer. */

} seq_args;

/* oligo_array is used to store a list of oligos or primers */
typedef struct oligo_array {

 /* Array of oligo (primer) records. */
 primer_rec *oligo;

 /* Number of initialized elements */
 int num_elem;

 /* Storage lengths of oligo */
 int storage_size;

 /* Type of oligos in the array */
 oligo_type type;

 /* Primers statistics. */
 oligo_stats expl;

} oligo_array;

/*
 * The return value for for primer3. 
 * After use, free memory with destroy_p3retval().
 */
typedef struct p3retval {
        
  /* Arrays of oligo (primer) records. */
  oligo_array fwd, intl, rev;
        
  /* Array of best primer pairs */
  pair_array_t best_pairs;

  /* Enum to store type of output */
  p3_output_type output_type;

  /* Place for error messages */
  pr_append_str glob_err;
  pr_append_str per_sequence_err;
  pr_append_str warnings;

  /* 
   * An optional _output_, meaninful if a
   * start_codon_pos is "not null".  The position of
   * the intial base of the leftmost stop codon that
   * is to the right of sa->start_codon_pos.
   */
  int stop_codon_pos;

  int upstream_stop_codon;  /* TO DO needs docs */

} p3retval;

/* Deallocate a primer3 state */
void destroy_p3retval(p3retval *);


/* get elements of p3retval */
const pair_array_t *p3_get_rv_best_pairs(const p3retval *r);
const oligo_array *p3_get_rv_fwd(const p3retval *r);
const oligo_array *p3_get_rv_intl(const p3retval *r);
const oligo_array *p3_get_rv_rev(const p3retval *r);

/* TO DO -- needs more documentation
   It is the responsibility of caller to free the return value. */
char *p3_get_rv_and_gs_warnings(const p3retval *retval, 
                                      const p3_global_settings *pa);

/* Return a char * describing global errors (usually/always?)  errors
   caused by problems detected in the p3_global_settings * argument to
   choose primers. Returned storage is free'ed on calling
   destroy_p3retval(r).  Returns NULL if no error. */
const char *p3_get_rv_global_errors(const p3retval *r);

/* Return a char * describing per-sequence errors (usually/always?)
   errors caused by problems detected in the seq_args * argument to
   choose primers. Returned storage is free'ed on calling
   destroy_p3retval(r). Returns NULL if no error. */
const char *p3_get_rv_per_sequence_errors(const p3retval *r);

/* Return a char * describing warnings generated in
   choose_primers. Returned storage is free'ed on calling
   destroy_p3retval(r). Returns NULL if no warnings. */
const char *p3_get_rv_warnings(const p3retval *r);

p3_output_type p3_get_rv_output_type(const p3retval *r);
int            p3_get_rv_stop_codon_pos(p3retval *r);

/* Get elements of an oligo_array */
int p3_get_oa_n(const oligo_array *x);
const  primer_rec *p3_get_oa_i(const oligo_array *x, int i);

/* We still need accessors (getters) for the elements of
   primer_rec */

/* Functions for seq_args -- create, destroy, set, set slots */
seq_args *create_seq_arg();
void destroy_seq_args(seq_args *);
int p3_set_sa_sequence(seq_args *sargs, const char *sequence);
void p3_set_sa_primer_sequence_quality(seq_args *sargs, int quality);
int p3_set_sa_sequence_name(seq_args *sargs, const char* sequence_name);
int p3_set_sa_left_input(seq_args *sargs, const char *left_input);
int p3_set_sa_right_input(seq_args *sargs, const char *right_input);
int p3_set_sa_internal_input(seq_args *sargs, const char *internal_input);
void p3_set_sa_empty_quality(seq_args *sargs);
void p3_sa_add_to_quality_array(seq_args *sargs, int quality);



/* The following three functions return 0 on success,
   1 on error (no space for additional intervals). */
int p3_add_to_sa_tar2(seq_args *, int, int);
int p3_add_to_sa_excl2(seq_args *, int, int);
int p3_add_to_sa_excl_internal2(seq_args *, int, int);

const interval_array_t2 *p3_get_sa_tar2(const seq_args *sargs);
const interval_array_t2 *p3_get_sa_excl2(const seq_args *sargs);
const interval_array_t2 *p3_get_sa_excl_internal2(const seq_args *sargs);

/*
  use p3_add_to_interval_array(interval_array_t2 *interval_arr, int i1, int i2);
  to do the sets for tar2, excl2, nd excl_internal2
*/
int p3_add_to_interval_array(interval_array_t2 *interval_arr, int i1, int i2);

/*
  included region
*/
void p3_set_sa_incl_s(seq_args *sargs, int incl_s);
void p3_set_sa_incl_l(seq_args *sargs, int incl_l);

void p3_set_sa_n_quality(seq_args *sargs, int n_quality) ;
void p3_set_sa_start_codon_pos(seq_args *sargs, int start_codon_pos);
int p3_set_sa_sequence_file(seq_args *sargs, const char *sequence_file);
int p3_set_sa_trimmed_sequence(seq_args *sargs, const char *trimmed_sequence);
int p3_set_sa_trimmed_original_sequence(seq_args *sargs, const char *trimmed_original_sequence);
int p3_set_sa_upcased_sequence(seq_args *sargs, const char *upcased_sequencd);


/* ============================================================ */
/* Functions for p3_global_settings -- create, destroy,
   set slots, get slots */
/* ============================================================ */

p3_global_settings *p3_create_global_settings();
void p3_destroy_global_settings(p3_global_settings *);

void p3_set_gs_prmin (p3_global_settings * p , int val, int i);
void p3_set_gs_prmax (p3_global_settings * p , int val, int i);
void p3_set_gs_primer_opt_size(p3_global_settings * p , int val);
void p3_set_gs_primer_min_size(p3_global_settings * p , int val);
void p3_set_gs_primer_max_size(p3_global_settings * p , int val);
void p3_set_gs_primer_max_polyx(p3_global_settings * p , int val);
void p3_set_gs_primer_opt_tm(p3_global_settings * p , double product_opt_tm);
void p3_set_gs_primer_opt_gc_percent(p3_global_settings * p , double val);
void p3_set_gs_primer_min_tm(p3_global_settings * p , double product_min_tm);
void p3_set_gs_primer_max_tm(p3_global_settings * p , double product_max_tm);
void p3_set_gs_primer_max_diff_tm(p3_global_settings * p , double val);
void p3_set_gs_primer_tm_santalucia(p3_global_settings * p,
                                    tm_method_type val);
void p3_set_gs_primer_salt_corrections(p3_global_settings * p,
                                       salt_correction_type salt_corrections);
void p3_set_gs_primer_min_gc(p3_global_settings * p , double val);
void p3_set_gs_primer_max_gc(p3_global_settings * p , double val);
void p3_set_gs_primer_salt_conc(p3_global_settings * p , double val);
void p3_set_gs_primer_divalent_conc(p3_global_settings * p , double val);
void p3_set_gs_primer_dntp_conc(p3_global_settings * p , double val);
void p3_set_gs_primer_dna_conc(p3_global_settings * p , double val);
void p3_set_gs_primer_num_ns_accepted(p3_global_settings * p , int val);
void p3_set_gs_primer_product_opt_size(p3_global_settings * p , int val);
void p3_set_gs_primer_self_any(p3_global_settings * p , double val);
void p3_set_gs_primer_self_end(p3_global_settings * p , double val);
void p3_set_gs_primer_file_flag(p3_global_settings * p , int val);
void p3_set_gs_primer_pick_anyway(p3_global_settings * p , int val);
void p3_set_gs_primer_gc_clamp(p3_global_settings * p , int val);
void p3_set_gs_primer_explain_flag(p3_global_settings * p , int val);
void p3_set_gs_primer_liberal_base(p3_global_settings * p , int val);
void p3_set_gs_primer_first_base_index(p3_global_settings * p , int val);
void p3_set_gs_primer_num_return(p3_global_settings * p , int val);
void p3_set_gs_primer_min_quality(p3_global_settings * p , int val);
void p3_set_gs_primer_min_end_quality(p3_global_settings * p , int val);
void p3_set_gs_primer_quality_range_min(p3_global_settings * p , int val);
void p3_set_gs_primer_quality_range_max(p3_global_settings * p , int val);
void p3_set_gs_primer_product_max_tm(p3_global_settings * p , double val);
void p3_set_gs_primer_product_min_tm(p3_global_settings * p , double val);
void p3_set_gs_primer_product_opt_tm(p3_global_settings * p , double val);
void p3_set_gs_primer_task(p3_global_settings * p , char * primer_task);
void p3_set_gs_primer_pick_left_primer(p3_global_settings * p , int pick_left_primer);
void p3_set_gs_primer_pick_right_primer(p3_global_settings * p , int pick_right_primer);
void p3_set_gs_primer_pick_internal_oligo(p3_global_settings * p , int pick_internal_oligo);
void p3_set_gs_primer_internal_oligo_opt_size(p3_global_settings * p , int val);
void p3_set_gs_primer_internal_oligo_max_size(p3_global_settings * p , int val);
void p3_set_gs_primer_internal_oligo_min_size(p3_global_settings * p , int val);
void p3_set_gs_primer_internal_oligo_max_poly_x(p3_global_settings * p , int val);
void p3_set_gs_primer_internal_oligo_opt_tm(p3_global_settings * p , double val);
void p3_set_gs_primer_internal_oligo_max_tm(p3_global_settings * p , double val);
void p3_set_gs_primer_internal_oligo_min_tm(p3_global_settings * p , double val);
void p3_set_gs_primer_internal_oligo_min_gc(p3_global_settings * p , double val);
void p3_set_gs_primer_internal_oligo_max_gc(p3_global_settings * p , double val);
void p3_set_gs_primer_internal_oligo_salt_conc(p3_global_settings * p , double val);
void p3_set_gs_primer_internal_oligo_divalent_conc(p3_global_settings * p , double val);
void p3_set_gs_primer_internal_oligo_dntp_conc(p3_global_settings * p , double val);
void p3_set_gs_primer_internal_oligo_dna_conc(p3_global_settings * p , double val);
void p3_set_gs_primer_internal_oligo_num_ns(p3_global_settings * p , int val);
void p3_set_gs_primer_internal_oligo_min_quality(p3_global_settings * p , int val);
void p3_set_gs_primer_internal_oligo_self_any(p3_global_settings * p , double val);
void p3_set_gs_primer_internal_oligo_self_end(p3_global_settings * p , double val);
void p3_set_gs_primer_max_mispriming(p3_global_settings * p , double val);
void p3_set_gs_primer_internal_oligo_max_mishyb(p3_global_settings * p , double val);
void p3_set_gs_primer_pair_max_mispriming(p3_global_settings * p , double val);
void p3_set_gs_primer_max_template_mispriming(p3_global_settings * p , double val);
void p3_set_gs_primer_internal_oligo_max_template_mishyb(p3_global_settings * p , double val);
void p3_set_gs_primer_lib_ambiguity_codes_consensus(p3_global_settings * p , int val);
void p3_set_gs_primer_inside_penalty(p3_global_settings * p , double val);
void p3_set_gs_primer_outside_penalty(p3_global_settings * p , double val);
void p3_set_gs_primer_mispriming_library(p3_global_settings * p , char * val);
void p3_set_gs_primer_internal_oligo_mishyb_library(p3_global_settings * p , char * val);
void p3_set_gs_primer_max_end_stability(p3_global_settings * p , double val);
void p3_set_gs_primer_lowercase_masking(p3_global_settings * p , int val);
void p3_set_gs_primer_wt_tm_gt(p3_global_settings * p , double val);
void p3_set_gs_primer_wt_tm_lt(p3_global_settings * p , double val);
void p3_set_gs_primer_wt_gc_percent_gt(p3_global_settings * p , double val);
void p3_set_gs_primer_wt_gc_percent_lt(p3_global_settings * p , double val);
void p3_set_gs_primer_wt_size_lt(p3_global_settings * p , double val);
void p3_set_gs_primer_wt_size_gt(p3_global_settings * p , double val);
void p3_set_gs_primer_wt_compl_any(p3_global_settings * p , double val);
void p3_set_gs_primer_wt_compl_end(p3_global_settings * p , double val);
void p3_set_gs_primer_wt_num_ns(p3_global_settings * p , double val);
void p3_set_gs_primer_wt_rep_sim(p3_global_settings * p , double val);
void p3_set_gs_primer_wt_seq_qual(p3_global_settings * p , double val);
void p3_set_gs_primer_wt_end_qual(p3_global_settings * p , double val);
void p3_set_gs_primer_wt_pos_penalty(p3_global_settings * p , double val);
void p3_set_gs_primer_wt_end_stability(p3_global_settings * p , double val);
void p3_set_gs_primer_wt_template_mispriming(p3_global_settings * p , double val);
void p3_set_gs_primer_io_wt_tm_gt(p3_global_settings * p , double val);
void p3_set_gs_primer_io_wt_tm_lt(p3_global_settings * p , double val);
void p3_set_gs_primer_io_wt_gc_percent_gt(p3_global_settings * p , double val);
void p3_set_gs_primer_io_wt_gc_percent_lt(p3_global_settings * p , double val);
void p3_set_gs_primer_io_wt_size_lt(p3_global_settings * p , double val);
void p3_set_gs_primer_io_wt_size_gt(p3_global_settings * p , double val);
void p3_set_gs_primer_io_wt_wt_coml_any(p3_global_settings * p , double val);
void p3_set_gs_primer_io_wt_compl_end(p3_global_settings * p , double val);
void p3_set_gs_primer_io_wt_num_ns(p3_global_settings * p , double val);
void p3_set_gs_primer_io_wt_rep_sim(p3_global_settings * p , double val);
void p3_set_gs_primer_io_wt_seq_qual(p3_global_settings * p , double val);
void p3_set_gs_primer_io_wt_end_qual(p3_global_settings * p , double val);
void p3_set_gs_primer_io_wt_template_mishyb(p3_global_settings * p , double val);
void p3_set_gs_primer_pair_wt_pr_penalty(p3_global_settings * p , double val);
void p3_set_gs_primer_pair_wt_io_penalty(p3_global_settings * p , double val);
void p3_set_gs_primer_pair_wt_diff_tm(p3_global_settings * p , double val);
void p3_set_gs_primer_pair_wt_compl_any(p3_global_settings * p , double val);
void p3_set_gs_primer_pair_wt_compl_end(p3_global_settings * p , double val);
void p3_set_gs_primer_pair_wt_product_tm_lt(p3_global_settings * p , double val);
void p3_set_gs_primer_pair_wt_product_tm_gt(p3_global_settings * p , double val);
void p3_set_gs_primer_pair_wt_product_size_gt(p3_global_settings * p , double val);
void p3_set_gs_primer_pair_wt_product_size_lt(p3_global_settings * p , double val);
void p3_set_gs_primer_pair_wt_rep_sim(p3_global_settings * p , double val);
void p3_set_gs_primer_pair_wt_template_mispriming(p3_global_settings * p , double val);
void p3_set_gs_first_base_index(p3_global_settings * p , int first_base_index);
void p3_set_gs_liberal_base(p3_global_settings * p , int liberal_base);
void p3_set_gs_num_return(p3_global_settings * p , int num_return);
void p3_set_gs_pick_anyway(p3_global_settings * p , int pick_anyway);
void p3_set_gs_lib_ambiguity_codes_consensus(p3_global_settings * p , int lib_ambiguity_codes_consensus);
void p3_set_gs_quality_range_min(p3_global_settings * p , int quality_range_min);
void p3_set_gs_quality_range_max(p3_global_settings * p , int quality_range_max);

void p3_empty_gs_product_size_range(p3_global_settings *pgs);
/* Return 1 on error (product size range is full);
   otherwise return 0. */

int  p3_add_to_gs_product_size_range(p3_global_settings *, int, int);
args_for_one_oligo_or_primer *p3_get_global_settings_p_args(p3_global_settings * p);
args_for_one_oligo_or_primer *p3_get_global_settings_o_args(p3_global_settings * p);
int p3_set_afogop_seq_lib(args_for_one_oligo_or_primer *, seq_lib *);
int p3_set_afogop_opt_tm(args_for_one_oligo_or_primer *, double);

void p3_set_gs_max_end_stability(p3_global_settings * p , int max_end_stability);
void p3_set_gs_gc_clamp(p3_global_settings * p , int gc_clamp);
void p3_set_gs_lowercase_masking(p3_global_settings * p , int lowercase_masking);
void p3_set_gs_outside_penalty(p3_global_settings * p , double outside_penalty);
void p3_set_gs_inside_penalty(p3_global_settings * p , double inside_penalty);

/* DO NOT USE */
/* void p3_set_gs_num_intervals(p3_global_settings * p , int num_intervals); */

void p3_set_gs_pair_max_template_mispriming(p3_global_settings * p ,
                                            double pair_max_template_mispriming);

void p3_set_gs_pair_repeat_compl(p3_global_settings * p, double pair_repeat_compl); 

void p3_set_gs_pair_compl_any(p3_global_settings * p , double  pair_compl_any);
void p3_set_gs_pair_compl_end(p3_global_settings * p , double  pair_compl_end);


void p3_set_gs_min_three_prime_distance(p3_global_settings *p, int min_distance);

/* 
 * Choose individual primers or oligos, or primer pairs, or primer
 * pairs with internal oligos. On ENOMEM return NULL and set errno. 
 * Otherwise return retval (updated).  Errors are returned in 
 * in retval.
 */
p3retval *choose_primers(const p3_global_settings *pa, 
                         seq_args *sa);

/* For testing/debugging: print the values in pa and sa to stdout. */
void p3_print_args(const p3_global_settings *pa, seq_args *sa) ;

/* Print out the content of one primer array */
/* Return 1 on error, otherwise 0. */
int    p3_print_one_oligo_list(const seq_args *, 
                               int, const primer_rec[],
                               const oligo_type, const int, 
                               const int, FILE *);

char  *pr_oligo_sequence(const seq_args *, const primer_rec *);

char  *pr_oligo_rev_c_sequence(const seq_args *, const primer_rec *);

void  pr_set_default_global_args(p3_global_settings *);

/* Return NULL on ENOMEM */
pr_append_str *create_pr_append_str();
void          init_pr_append_str(pr_append_str *s);

void          pr_set_empty(pr_append_str *);
int           pr_is_empty(const pr_append_str *);
void          destroy_pr_append_str(pr_append_str *);
void          destroy_pr_append_str_data(pr_append_str *str);

/* Return 1 on ENOMEM, otherwise 0 */
int           pr_append_external(pr_append_str *, const char *);

int           pr_append_w_sep_external(pr_append_str *x, 
                                       const char *sep,
                                       const char *s);

int           pr_append_new_chunk_external(pr_append_str *, const char *);

const char *  pr_append_str_chars(const pr_append_str *x);

/* Warning, return pointers to static storage;  overwritten on each call. */
const char *p3_get_pair_array_explain_string(const pair_array_t *);
const char *p3_get_oligo_array_explain_string(const oligo_array *);

const char  *libprimer3_release(void);
const char *primer3_copyright(void);

/* An accessor function for a primer_rec *. */
short oligo_max_template_mispriming(const primer_rec *);

int   strcmp_nocase(const char *, const char *);

void
p3_set_program_name(const char *name);

/* 
 * Utility function for C clients -- will not overflow
 * buffer.  Warning: points to static storage that is over-written
 * on the next call.
 */
char* p3_read_line(FILE *file);

int        p3_ol_has_any_problem(const primer_rec *oligo);
const char *p3_get_ol_problem_string(const primer_rec *oligo);

int    p3_print_oligo_lists(const p3retval*, 
                            const seq_args *, 
                            const p3_global_settings *, 
                            pr_append_str *err,
                            const char* file_name);

/* Must be called from main to allow dependency on 
 * input version=3 in future version can it be 
 * deleted.                                       */
void add_must_use_warnings(pr_append_str *warning,
                      const char* text,
                      const oligo_stats *stats);



/* Reverse and complement the sequence seq and put the result in s.
   WARNING: It is up the caller to ensure that s points to enough
   space. */ 
void   p3_reverse_complement(const char *, char *);

#ifdef __cplusplus
    }
#endif

#endif
