/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007
Whitehead Institute for Biomedical Research, Steve Rozen
(http://jura.wi.mit.edu/rozen), and Helen Skaletsky
All rights reserved.

    This file is part of primer3.

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

#include <string.h> /* strcpy */
#include "print_boulder.h"

static char *pr_program_name = "Program name is probably primer3_core";

static void   print_all_explain(const primer_args *, const seq_args *);
static void   print_explain(const oligo_stats *, oligo_type);

/* Print the data for chosen primer pairs to stdout in "boulderio" format. */
void
boulder_print_pairs(prog_args, pa, sa, best_pairs)
    const program_args *prog_args;
    const primer_args *pa;
    const seq_args *sa;
    const pair_array_t *best_pairs;
{
    const char *left_tag, *right_tag, *intl_tag, *prod_size_tag;
    char *warning;
    char suffix [3];
    primer_rec *fwd, *rev, *intl;
    int i, incl_s = sa->incl_s;

    PR_ASSERT(NULL != pa);
    PR_ASSERT(NULL != sa);
    PR_ASSERT(NULL != prog_args);

    if (0 /* prog_args->twox_compat */) {
      /* twox_compat no longer supported. */
	left_tag = "FORWARD_PRIMER";
	right_tag = "REVERSE_PRIMER";
	intl_tag = "MIDDLE_OLIGO";
	prod_size_tag = "PRODUCT_SIZE";
    } else {
	left_tag = "PRIMER_LEFT";
	right_tag = "PRIMER_RIGHT";
	intl_tag = "PRIMER_INTERNAL_OLIGO";
	prod_size_tag = "PRIMER_PRODUCT_SIZE";
    }

    if ((warning = pr_gather_warnings(sa, pa)) != NULL) {
      
	printf("PRIMER_WARNING=%s\n", warning);
	free(warning);
    }
    
    if (sa->error.data != NULL) {
	printf("PRIMER_ERROR=%s\n=\n", sa->error.data);
	if (fflush(stdout) == EOF) {
	    perror("fflush(stdout) failed");
	    exit(-1);
	}
	return;
    }

    if (pa->explain_flag) print_all_explain(pa, sa);

    if (!PR_START_CODON_POS_IS_NULL(sa))
      printf("PRIMER_STOP_CODON_POSITION=%d\n", sa->stop_codon_pos);

    for(i=0; i<best_pairs->num_pairs; i++) {
	fwd = best_pairs->pairs[i].left;
	rev = best_pairs->pairs[i].right;
	intl = best_pairs->pairs[i].intl;

	if (i == 0) suffix[0] = '\0';
	else sprintf(suffix, "_%d", i);

	printf("PRIMER_PAIR_PENALTY%s=%.4f\n", suffix,
	       best_pairs->pairs[i].pair_quality);
	printf("PRIMER_LEFT%s_PENALTY=%f\n", suffix,
	       fwd->quality);
	printf("PRIMER_RIGHT%s_PENALTY=%f\n", suffix,
		   rev->quality);
	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe)
	  printf("PRIMER_INTERNAL_OLIGO%s_PENALTY=%f\n", suffix,
		 intl->quality);

        /* Print sequences. */
	printf("PRIMER_LEFT%s_SEQUENCE=%s\n", suffix,
	       pr_oligo_sequence(sa, fwd));
	printf("PRIMER_RIGHT%s_SEQUENCE=%s\n", suffix,
	       pr_oligo_rev_c_sequence(sa, rev));
	if( pa->primer_task == pick_pcr_primers_and_hyb_probe)
#if USE_OLD_FORMAT_MISTAKE
	    printf("PRIMER_INTERNAL%s_OLIGO_SEQUENCE=%s\n", suffix,
		   pr_oligo_sequence(sa,intl));
#else
	    printf("PRIMER_INTERNAL_OLIGO%s_SEQUENCE=%s\n", suffix,
		   pr_oligo_sequence(sa,intl));
#endif

	printf("%s%s=%d,%d\n", left_tag, suffix,
	       fwd->start + incl_s + pa->first_base_index,
	       fwd->length);
	printf("%s%s=%d,%d\n", right_tag, suffix,
	       rev->start + incl_s + pa->first_base_index,
	       rev->length);
	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe)
	    printf("%s%s=%d,%d\n", intl_tag, suffix,
		   intl->start + incl_s + pa->first_base_index,
		   intl->length);

	printf("PRIMER_LEFT%s_TM=%.3f\n", suffix, fwd->temp);
	printf("PRIMER_RIGHT%s_TM=%.3f\n", suffix, rev->temp);
	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe)
	    printf("PRIMER_INTERNAL_OLIGO%s_TM=%.3f\n",suffix, intl->temp);


	printf("PRIMER_LEFT%s_GC_PERCENT=%.3f\n", suffix, fwd->gc_content);
	printf("PRIMER_RIGHT%s_GC_PERCENT=%.3f\n", suffix, rev->gc_content);
	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe)
	    printf("PRIMER_INTERNAL_OLIGO%s_GC_PERCENT=%.3f\n",suffix,
		   intl->gc_content);

	printf("PRIMER_LEFT%s_SELF_ANY=%.2f\n", suffix,
	       fwd->self_any / PR_ALIGN_SCORE_PRECISION);
	printf("PRIMER_RIGHT%s_SELF_ANY=%.2f\n", suffix,
	       rev->self_any / PR_ALIGN_SCORE_PRECISION);
	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe)
	    printf("PRIMER_INTERNAL_OLIGO%s_SELF_ANY=%.2f\n", suffix,
		   intl->self_any / PR_ALIGN_SCORE_PRECISION);
	printf("PRIMER_LEFT%s_SELF_END=%.2f\n", suffix,
	       fwd->self_end / PR_ALIGN_SCORE_PRECISION);
	printf("PRIMER_RIGHT%s_SELF_END=%.2f\n",
	       suffix,rev->self_end / PR_ALIGN_SCORE_PRECISION);
	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe)
	    printf("PRIMER_INTERNAL_OLIGO%s_SELF_END=%.2f\n", suffix,
		   intl->self_end / PR_ALIGN_SCORE_PRECISION);
        if (seq_lib_num_seq(pa->p_args.repeat_lib) > 0) {
	    printf("PRIMER_LEFT%s_MISPRIMING_SCORE=%.2f, %s\n", suffix,
		   fwd->repeat_sim.score[fwd->repeat_sim.max] / PR_ALIGN_SCORE_PRECISION,
		   fwd->repeat_sim.name);
            printf("PRIMER_RIGHT%s_MISPRIMING_SCORE=%.2f, %s\n", suffix,
		   rev->repeat_sim.score[rev->repeat_sim.max] / PR_ALIGN_SCORE_PRECISION,
                   rev->repeat_sim.name);
            printf("PRIMER_PAIR%s_MISPRIMING_SCORE=%.2f, %s\n", suffix,
		   best_pairs->pairs[i].repeat_sim / PR_ALIGN_SCORE_PRECISION,
		   best_pairs->pairs[i].rep_name);
        }
	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe
	    && seq_lib_num_seq(pa->o_args.repeat_lib) > 0)
	    printf("PRIMER_INTERNAL_OLIGO%s_MISHYB_SCORE=%.2f, %s\n", suffix,
		   intl->repeat_sim.score[intl->repeat_sim.max]
		   / PR_ALIGN_SCORE_PRECISION,
		   intl->repeat_sim.name);
	if (NULL != sa->quality){
	   printf("PRIMER_LEFT%s_MIN_SEQ_QUALITY=%d\n", suffix,
		   fwd->seq_quality);
           printf("PRIMER_RIGHT%s_MIN_SEQ_QUALITY=%d\n", suffix,
		   rev->seq_quality);
        }

	if (!_PR_DEFAULT_POSITION_PENALTIES(pa)
            || !PR_START_CODON_POS_IS_NULL(sa)) {
	  printf("PRIMER_LEFT%s_POSITION_PENALTY=%f\n", suffix,
		 fwd->position_penalty);
	  printf("PRIMER_RIGHT%s_POSITION_PENALTY=%f\n", suffix,
		 rev->position_penalty);
	}

	printf("PRIMER_LEFT%s_END_STABILITY=%.4f\n",
	       suffix, fwd->end_stability);
	printf("PRIMER_RIGHT%s_END_STABILITY=%.4f\n",
	       suffix, rev->end_stability);

	if (oligo_max_template_mispriming(fwd) != ALIGN_SCORE_UNDEF)
	  printf("PRIMER_LEFT%s_TEMPLATE_MISPRIMING=%.4f\n", suffix,
		 oligo_max_template_mispriming(fwd)
		 / PR_ALIGN_SCORE_PRECISION);

	if (oligo_max_template_mispriming(rev) != ALIGN_SCORE_UNDEF)
	  printf("PRIMER_RIGHT%s_TEMPLATE_MISPRIMING=%.4f\n", suffix,
		 oligo_max_template_mispriming(rev)
		 / PR_ALIGN_SCORE_PRECISION);


	if ( pa->primer_task == pick_pcr_primers_and_hyb_probe
	    && NULL != sa->quality)
	   printf("PRIMER_INTERNAL_OLIGO%s_MIN_SEQ_QUALITY=%d\n",
		   suffix, intl->seq_quality);
	printf("PRIMER_PAIR%s_COMPL_ANY=%.2f\n", suffix,
	       best_pairs->pairs[i].compl_any / PR_ALIGN_SCORE_PRECISION);
	printf("PRIMER_PAIR%s_COMPL_END=%.2f\n", suffix,
	       best_pairs->pairs[i].compl_end  / PR_ALIGN_SCORE_PRECISION);

	printf("%s%s=%d\n", prod_size_tag, suffix,
				     best_pairs->pairs[i].product_size);

	if (pa->product_max_tm != PR_DEFAULT_PRODUCT_MAX_TM
	    || pa->product_min_tm != PR_DEFAULT_PRODUCT_MIN_TM) {
	  printf("PRIMER_PRODUCT_TM%s=%.4f\n", suffix,
		 best_pairs->pairs[i].product_tm);

	  printf("PRIMER_PRODUCT_TM_OLIGO_TM_DIFF%s=%.4f\n", suffix,
	      best_pairs->pairs[i].product_tm_oligo_tm_diff);

	   printf("PRIMER_PAIR%s_T_OPT_A=%.4f\n", suffix,
	      best_pairs->pairs[i].t_opt_a);
	}

	if (best_pairs->pairs[i].template_mispriming != ALIGN_SCORE_UNDEF)
	  printf("PRIMER_PAIR%s_TEMPLATE_MISPRIMING=%.2f\n", suffix,
		 best_pairs->pairs[i].template_mispriming 
		 / PR_ALIGN_SCORE_PRECISION);

    }
    printf("=\n");
    if (fflush(stdout) == EOF) {
	perror("fflush(stdout) failed");
	exit(-1);
    }
}

void 
boulder_print_oligos(pa, sa, n, l, f, r, mid)
    const primer_args *pa;
    const seq_args *sa;
    int n;
    oligo_type l;
    primer_rec *f;
    primer_rec *r;
    primer_rec *mid;
{
    char *warning;
    int i, j;
    char suffix [3], type[256];
    /* type must be larger than the length of "PRIMER_INTERNAL_OLIGO". */

    primer_rec *oligo;
    int incl_s = sa->incl_s;


    if ((warning = pr_gather_warnings(sa, pa)) != NULL) {
	printf("PRIMER_WARNING=%s\n", warning);
	free(warning);
    }
    if (sa->error.data != NULL) {
	printf("PRIMER_ERROR=%s\n=\n", sa->error.data);
	if (fflush(stdout) == EOF) {
	    perror("fflush(stdout) failed");
	    exit(-1);
        }
	return;
    }

    if(l == OT_LEFT) strcpy(type, "PRIMER_LEFT");
    else if(l == OT_RIGHT) strcpy(type, "PRIMER_RIGHT");
    else strcpy(type, "PRIMER_INTERNAL_OLIGO");

    if (pa->explain_flag) print_all_explain(pa, sa);

    i = 0;
    j = (pa->num_return < n) ? pa->num_return : n;
    if(l == OT_LEFT) oligo = f;
    else if (l == OT_RIGHT) oligo = r;
    else  oligo = mid;

    while(i < j) {
	if (i == 0) suffix[0] = '\0';
	else sprintf(suffix, "_%d", i);

	printf("%s%s_PENALTY=%.4f\n", type, suffix, oligo[i].quality);
	if(l == OT_RIGHT)
	   printf("%s%s_SEQUENCE=%s\n", type, suffix,
		  pr_oligo_rev_c_sequence(sa, &oligo[i]));
        else printf("%s%s_SEQUENCE=%s\n", type, suffix,
		    pr_oligo_sequence(sa, &oligo[i])); 
	printf("%s%s=%d,%d\n", type, suffix, 
		       oligo[i].start + incl_s + pa->first_base_index,
		       oligo[i].length);
	printf("%s%s_TM=%.3f\n", type, suffix, oligo[i].temp);
	printf("%s%s_GC_PERCENT=%.3f\n", type, suffix, oligo[i].gc_content);
        printf("%s%s_SELF_ANY=%.2f\n", type, suffix,
		       oligo[i].self_any / PR_ALIGN_SCORE_PRECISION);
        printf("%s%s_SELF_END=%.2f\n", type, suffix,
		       oligo[i].self_end / PR_ALIGN_SCORE_PRECISION);
        if ((l == OT_LEFT || l == OT_RIGHT) && seq_lib_num_seq(pa->p_args.repeat_lib) > 0 )
	    printf("%s%s_MISPRIMING_SCORE=%.2f, %s\n", type, suffix,
		    oligo[i].repeat_sim.score[oligo[i].repeat_sim.max] /PR_ALIGN_SCORE_PRECISION,
		    oligo[i].repeat_sim.name);
        if (l == OT_INTL && seq_lib_num_seq(pa->o_args.repeat_lib) > 0)
	    printf("%s%s_MISHYB_SCORE=%.2f,%s\n", type, suffix,
	           oligo[i].repeat_sim.score[oligo[i].repeat_sim.max]/ PR_ALIGN_SCORE_PRECISION,
                   oligo[i].repeat_sim.name);
	if (NULL != sa->quality)printf("%s%s_MIN_SEQ_QUALITY=%d\n",
		   type, suffix, oligo[i].seq_quality);
        if (PR_DEFAULT_INSIDE_PENALTY != pa->inside_penalty
	    || PR_DEFAULT_OUTSIDE_PENALTY != pa->outside_penalty != 0.0)
	    printf("%s%s_POSITION_PENALTY=%f\n", type, suffix,
		    oligo[i].position_penalty);
	if(l == OT_LEFT || l == OT_RIGHT)
	    printf("%s%s_END_STABILITY=%.4f\n", type, suffix,
		    oligo[i].end_stability);

	if (oligo_max_template_mispriming(&oligo[i]) != ALIGN_SCORE_UNDEF)
	  printf("%s%s_TEMPLATE_MISPRIMING=%.2f\n", type, suffix,
		 oligo_max_template_mispriming(&oligo[i])
		 / PR_ALIGN_SCORE_PRECISION);

        i++;
    }
    printf("=\n");
}

static void
print_all_explain(pa, sa)
    const primer_args *pa;
    const seq_args *sa;
{
  if (pa->explain_flag) {
    if (pa->primer_task != pick_right_only
	&& pa->primer_task != pick_hyb_probe_only
	&& !(pa->pick_anyway && sa->left_input))
      print_explain(&sa->left_expl,OT_LEFT);

    if (pa->primer_task != pick_left_only 
	&& pa->primer_task != pick_hyb_probe_only
	&& !(pa->pick_anyway && sa->right_input))
      print_explain(&sa->right_expl,OT_RIGHT);

    if ((pa->primer_task == pick_hyb_probe_only
	 || pa->primer_task == pick_pcr_primers_and_hyb_probe)
	&& !(pa->pick_anyway && sa->internal_input)) 
      print_explain(&sa->intl_expl, OT_INTL);

    if (pa->primer_task  == pick_pcr_primers
	|| pa->primer_task == pick_pcr_primers_and_hyb_probe) {
      printf("PRIMER_PAIR_EXPLAIN=");
      pr_print_pair_explain(stdout, sa);
    }
  }
}

static void
print_explain(stat, l)
    const oligo_stats *stat;
    oligo_type l;
{
    if(OT_LEFT == l)printf("PRIMER_LEFT_EXPLAIN=");
    else if(OT_RIGHT == l)printf("PRIMER_RIGHT_EXPLAIN=");
	 else printf("PRIMER_INTERNAL_OLIGO_EXPLAIN=");

    printf("considered %d", stat->considered);
    if (stat->no_orf) printf(", would not amplify any of the ORF %d", stat->no_orf);
    if(stat->ns)printf(", too many Ns %d", stat->ns);
    if(stat->target)printf(", overlap target %d", stat->target);
    if(stat->excluded)printf(", overlap excluded region %d", stat->excluded);
    if(stat->gc)printf(", GC content failed %d", stat->gc);
    if(stat->gc_clamp)printf(", GC clamp failed %d", stat->gc_clamp);
    if(stat->temp_min)printf(", low tm %d", stat->temp_min);
    if(stat->temp_max)printf(", high tm %d", stat->temp_max);
    if(stat->compl_any)printf(", high any compl %d", stat->compl_any);
    if(stat->compl_end)printf(", high end compl %d", stat->compl_end);
    if(stat->repeat_score) printf(", high repeat similarity %d", stat->repeat_score);
    if(stat->poly_x)printf(", long poly-x seq %d", stat->poly_x);
    if(stat->seq_quality)printf(",low sequence quality %d", stat->seq_quality);
    if (stat->stability) printf(",high 3' stability %d", stat->stability);
    if (stat->template_mispriming) printf(",high template mispriming score %d",
					  stat->template_mispriming);
    /* edited by T. Koressaar for lowercase masking */
    if(stat->gmasked) printf(",lowercase masking of 3' end %d",stat->gmasked);
   
    printf(", ok %d\n", stat->ok);
}

