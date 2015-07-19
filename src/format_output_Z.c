/*
 Copyright Notice and Disclaimer for Edesign

 Copyright (c) 2013,2014,2015 RIKEN and K.K.DNAFORM. All Rights Reserved
 The Edesign is based on the Primer3 program (version 2.3.4) of the Whitehead Institute (http://primer3.ut.ee/).
 
       This file is part of Edesign software.

       This software is free software;
       you can redistribute it and/or modify it under the terms
       of the GNU General Public License as published by the Free
       Software Foundation; either version 2 of the License, or (at
       your option) any later version.

       This software is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with this software (file gpl-2.0.txt in the source
       distribution); if not, write to the Free Software
       Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON A THEORY
 OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 */

#include <stdio.h>
#include <string.h>
#include "format_output_Z.h"
#include "libprimer3_Z.h"

#define FORWARD 1
#define REVERSE -1

/* Oligo titiles used in the output page - Y.Kimura */
/* Normal oligo titiles */
#define TITLE_LEFT        "LEFT Primer"
#define TITLE_RIGHT       "RIGHT Primer"
#define TITLE_INTL_FW     "INTERNAL Probe (Forward)"
#define TITLE_INTL_RV     "INTERNAL Probe (Reverse)"
/* Modified oligo titiles */
#define TITLE_MOD_LEFT    "LEFT Eprimer"
#define TITLE_MOD_RIGHT   "RIGHT Eprimer"
#define TITLE_MOD_INTL_FW "INTERNAL Eprobe (Forward)"
#define TITLE_MOD_INTL_RV "INTERNAL Eprobe (Reverse)"

static void format_pairs(FILE *f, 
                         const p3_global_settings *pa,
                         const seq_args *sa, 
                         const p3retval *retval, 
                         const pair_array_t *best_pairs,
                         const char *pr_release,
                         const pr_append_str *combined_retval_err,
                         int explain_flag);

static void format_oligos(FILE *, 
                          const p3_global_settings *, 
                          const seq_args *, 
                          const p3retval *retval, 
                          const char*,
                          const pr_append_str *combined_retval_err,
                          int explain_flag);

static int lib_sim_specified(const p3_global_settings *);
static void print_explain(FILE *, const p3_global_settings *, const seq_args *,
                          const p3retval *retval, int, const char *);
static void print_pair_info(FILE *, const primer_pair *,
                            const p3_global_settings *);
/*
static void print_oligo(FILE *, const char *, const seq_args *,
                        const primer_rec *, int, const p3_global_settings *, 
                        const seq_lib*, int);
*/

/* Function of printing result including Z code
 * Created - 20120911 N.Kasahara
 * Modified -  20140109 Y.Kimura 
 */
static void print_oligo_Z(FILE *, const char *, const seq_args *,
                        const primer_rec *, oligo_type, const p3_global_settings *, 
                        const seq_lib*, int);

static void print_oligo_header(FILE *, const char *, const int, const int, const int);
static void print_pair_array(FILE *, const char*, int,
                             const interval_array_t, 
                             const p3_global_settings*, const seq_args*);


static void print_2_pair_array(FILE *, const char*, int,
			       const interval_array_t, 
			       const interval_array_t, 
			       const p3_global_settings*, const seq_args*);
static void print_rest(FILE *, const p3_global_settings *, 
                       const seq_args *,  const pair_array_t *);

static int  print_seq(FILE *, const p3_global_settings *, const seq_args *, 
                      const p3retval *retval, primer_rec *h,
                      const pair_array_t *, int);

static void print_seq_lines(FILE *, char *s, char *n[2], int, int,
                            int something_found[2], const p3_global_settings *);

static void print_stat_line(FILE *, const char *, oligo_stats s, int, int, int);
static void print_summary(FILE *, const p3_global_settings *, 
                          const seq_args *, const pair_array_t *, int);

void
print_format_output(FILE *f,
                    const int *io_version,
                    const p3_global_settings *pa,
                    const seq_args *sa,
                    const p3retval *retval,
                    const char *pr_release,
                    int   explain_flag)
{  
  
  /* A place to put a string containing all error messages */
  pr_append_str *combined_retval_err = NULL;

  combined_retval_err = create_pr_append_str();
  if (NULL == combined_retval_err) exit(-2); /* Out of memory */

  if (pr_append_new_chunk_external(combined_retval_err, 
                                   retval->glob_err.data))
    exit(-2);

  if (pr_append_new_chunk_external(combined_retval_err, 
                                   retval->per_sequence_err.data)) 
    exit(-2);

  /* Print as primer pairs */
  if (retval->output_type == primer_pairs) {
    format_pairs(f, pa, sa, retval, &retval->best_pairs, 
                 pr_release, combined_retval_err, explain_flag);
    
    /* Print as primer list */
  } else {
    format_oligos(stdout, pa, sa, retval, pr_release,
                  combined_retval_err, explain_flag);
  }

  destroy_pr_append_str(combined_retval_err);
}

static void
format_pairs(FILE *f,
             const p3_global_settings *pa,
             const seq_args *sa,
             const p3retval *retval,
             const pair_array_t *best_pairs,
             const char *pr_release,
             const pr_append_str *combined_retval_err,
             int explain_flag)
{
  char *warning;
  int print_lib_sim = lib_sim_specified(pa);
  primer_rec *h = NULL;

  PR_ASSERT(NULL != f);
  PR_ASSERT(NULL != pa);
  PR_ASSERT(NULL != sa);
  
  /* If there are errors, print them and return */
  if (!pr_is_empty(combined_retval_err)) {
    format_error(f, sa->sequence_name, 
                 pr_append_str_chars(combined_retval_err));
    return;
    }
  
  /* Print the sequence name if it is provided */
  if (NULL != sa->sequence_name)
    fprintf(f, "PRIMER PICKING RESULTS FOR %s\n\n", sa->sequence_name);
  
  /* Print if a mispriming libraby was used and which one */
  if (pa->p_args.repeat_lib != NULL)
    fprintf(f, "Using mispriming library %s\n",
            pa->p_args.repeat_lib->repeat_file);
  else
    fprintf(f, "No mispriming library specified\n");

  /* Print if a mispriming libraby for the internal probe 
   * was used and which one */
  if ( pa->pick_internal_oligo == 1 ) {
    if (pa->o_args.repeat_lib != NULL)
      fprintf(f, "Using internal probe mishyb library %s\n",
              pa->o_args.repeat_lib->repeat_file);
    else
      fprintf(f, "No internal probe mishyb library specified\n");
  }

  /* Does the sequence start at position 0 or 1 ? */
  fprintf(f, "Using %d-based sequence positions\n",
          pa->first_base_index);
  
  /* Complain if no primers are in the array */
  if (best_pairs->num_pairs == 0) fprintf(f, "NO PRIMERS FOUND\n\n");
  
  /* Print out the warings */
  if ((warning = p3_get_rv_and_gs_warnings(retval, pa)) != NULL) {
    fprintf(f, "\nWARNING: %s\n", warning);
    free(warning);
  }
  
  /* Print the results for the best pair */
  print_summary(f, pa, sa, best_pairs, 0);
  fprintf(f, "\n");

  /* Print nicely out the sequence with the best pair */
  if (print_seq(f, pa, sa, retval, h, best_pairs, 0)) exit(-2); /* ENOMEM */
  
  /* Print out the alternative pairs */
  if (best_pairs->num_pairs > 1 ) print_rest(f, pa, sa, best_pairs);
  
  /* Print the primer picking statistics */
  if (explain_flag)
    print_explain(f, pa, sa, retval, print_lib_sim, pr_release);
  
  /* Flush the buffers and return */
  fprintf(f, "\n\n");
  if (fflush(f) == EOF) {
    perror("fflush(f) failed");
    exit(-1);
  }

}

/* Prints out the results of a primer pair */
/* add branching condition for printing internal_oligo_direction and modification 
 * Created - 20121204 N.Kasahara
 * Modified - 20140109 Y.Kimura
 */
static void
print_summary(FILE *f, const p3_global_settings *pa,
              const seq_args *sa, const pair_array_t *best_pairs,
              int num)
{
	int seq_len = strlen(sa->sequence);
	int print_lib_sim = lib_sim_specified(pa);
	primer_pair *p;
	p = best_pairs->pairs + num;


	if (best_pairs->num_pairs > 0) {
		/* 
		 * If the following format changes, also change the format in print_oligo.
		 */
		print_oligo_header(f, "\nOLIGO ", print_lib_sim, pa->thermodynamic_alignment, sa->tar2.genotyping);
		if(pa->modify_left_primer==1){
			print_oligo_Z(f, TITLE_MOD_LEFT, sa, p->left, OT_LEFT, pa, pa->p_args.repeat_lib, print_lib_sim);
		} else {
			print_oligo_Z(f, TITLE_LEFT, sa, p->left, OT_LEFT, pa, pa->p_args.repeat_lib, print_lib_sim);
		}

		if(pa->modify_right_primer==1){
			print_oligo_Z(f, TITLE_MOD_RIGHT, sa, p->right, OT_RIGHT, pa, pa->p_args.repeat_lib, print_lib_sim);
		} else {
			print_oligo_Z(f, TITLE_RIGHT, sa, p->right, OT_RIGHT, pa, pa->p_args.repeat_lib, print_lib_sim);
		}

		if ( pa->pick_internal_oligo == 1 ){
			if(pa->modify_internal_oligo==1){
			/* Modified internal probe */
				if((p->intl)->oligo_dir==0){
					print_oligo_Z(f, TITLE_MOD_INTL_FW, sa, p->intl, OT_INTL, pa, pa->o_args.repeat_lib, print_lib_sim);
				} else if((p->intl)->oligo_dir==1){
					print_oligo_Z(f, TITLE_MOD_INTL_RV, sa, p->intl, OT_INTL, pa, pa->o_args.repeat_lib, print_lib_sim);
				}
			} else {
			/* Internal probe */
				if((p->intl)->oligo_dir==0){
					print_oligo_Z(f, TITLE_INTL_FW, sa, p->intl, OT_INTL, pa, pa->o_args.repeat_lib, print_lib_sim);
				} else if((p->intl)->oligo_dir==1){
					print_oligo_Z(f, TITLE_INTL_RV, sa, p->intl, OT_INTL, pa, pa->o_args.repeat_lib, print_lib_sim);
				}
			}
		}
	}
	if (best_pairs->num_pairs > 0) { 
		print_pair_info(f, p, pa);
		/* Print penalty score - Y.Kimura */
		fprintf(f, "penalty score: %.4f\n", p->pair_quality);
	}
	fprintf(f, "\n");
	fprintf(f, "SEQUENCE SIZE: %d\n", seq_len);
	fprintf(f, "INCLUDED REGION SIZE: %d\n", sa->incl_l);
	if(sa->tar2.genotyping==0){
		print_pair_array(f, "TARGETS", sa->tar2.count, sa->tar2.pairs, pa, sa);
	}else if(sa->tar2.genotyping==1){
		fprintf(f,"TARGET VARIANT * (position,size):");
		if(sa->tar2.char_count==sa->tar2.VAR_sequence_number){
			fprintf(f," %d,%d",sa->tar2.VAR_start[0],sa->tar2.char_num[sa->tar2.VAR_sequence_number]-1);
		} else {
			fprintf(f," %d,%d",sa->tar2.VAR_start[0],sa->tar2.char_num[sa->tar2.VAR_sequence_number]);
		}
		fprintf(f,"\n\n");
	}
	print_pair_array(f, "EXCLUDED REGIONS", sa->excl2.count, sa->excl2.pairs, pa, sa);
	print_pair_array(f, "INTERNAL PROBE EXCLUDED REGIONS", sa->excl_internal2.count, sa->excl_internal2.pairs, pa, sa);
	print_2_pair_array(f, "PAIR_OK_REGIONS", sa->ok_regions.count, sa->ok_regions.left_pairs, sa->ok_regions.right_pairs, pa, sa);
}

/* Print column headers for lines printed by print_oligo(). */
static void
print_oligo_header(FILE *f, const char *s, const int print_lib_sim, const int thermodynamic_alignment, const int genotyping)
{
	if(genotyping==1){
		/* add tm against Variation : tm(VAR) */
		if(thermodynamic_alignment==0)
			fprintf(f, "%-25s  start  len      tm  tm(VAR)     gc%%   any    3' hairpin %sseq\n", s, print_lib_sim ? "  rep " : "");
		else 
			fprintf(f, "%-25s  start  len      tm  tm(VAR)     gc%%  any_th  3'_th hairpin %sseq\n", s, print_lib_sim ? "  rep " : "");
	}else{
		if(thermodynamic_alignment==0)
			fprintf(f, "%-25s  start  len      tm     gc%%   any    3' hairpin %sseq\n", s, print_lib_sim ? "  rep " : "");
		else 
			fprintf(f, "%-25s  start  len      tm     gc%%  any_th  3'_th hairpin %sseq\n", s, print_lib_sim ? "  rep " : "");
	}
}

/* Print the line with the parameters */
/* Function of printing result including Z
 * Created - 20120911 N.Kasahara
 * Modified -  20140109 Y.Kimura
 */
static void
print_oligo_Z(FILE *f,
            const char *title,
            const seq_args *sa,
            const primer_rec *o,
            oligo_type otype,
            const p3_global_settings *pa,
            const seq_lib *seqlib,
            int print_lib_sim)
{
    const char *format1;
    const char* seq = pr_oligo_sequence_Z(pa,otype,sa,o);
//    char *seq = (OT_LEFT == otype || (OT_INTL == otype && o->oligo_dir == 0) )
//        ? pr_oligo_sequence(sa, o) : pr_oligo_rev_c_sequence(sa, o);
//
//	int len=0;
//	len=strlen(seq);
//
//	if(otype==OT_LEFT){
//		/* LEFT Primer */
//		if(pa->znum_left_primer==1 && pa->modify_left_primer==1){
//			seq[pa->zpos_left_primer]='Z';
//		} else if(pa->znum_left_primer!=1 && pa->modify_left_primer==1){
//			seq[o->seq_z]='Z';
//		}
//	}else if(otype==OT_RIGHT){
//		/* RIGHT Primer */
//		if(pa->znum_right_primer==1 && pa->modify_right_primer==1){
//			seq[pa->zpos_right_primer]='Z';
//		}else if(pa->znum_right_primer!=1 && pa->modify_right_primer==1){
//			seq[len-o->seq_z-1]='Z';
//		}
//	}else if(otype==OT_INTL){
//		if(o->oligo_dir==0){
//		/* Forward Internal probe */
//			if(pa->znum_internal_oligo==1 && pa->modify_internal_oligo==1){
//				seq[pa->zpos_internal_oligo]='Z';
//			} else if(pa->znum_internal_oligo!=1 && pa->modify_internal_oligo==1){
//				seq[o->seq_z]='Z';
//			}
//		} else if(o->oligo_dir==1){
//		/* Reverse Internal probe */
//			if(pa->znum_internal_oligo_REV==1 && pa->modify_internal_oligo==1){
//				seq[pa->zpos_internal_oligo_REV]='Z';
//			} else if(pa->znum_internal_oligo_REV!=1 && pa->modify_internal_oligo==1){
//				seq[len-o->seq_z-1]='Z';
//			}
//		}
//	}

	if(sa->tar2.genotyping==1){
		if(pa->thermodynamic_alignment==1) {
			format1 = "%-25s %5d %4d %7.2f  %7.2f %7.2f   %5.2f  %5.2f ";
		} else {
			format1 = "%-25s %5d %4d %7.2f  %7.2f %7.2f %5.2f %5.2f ";
		}
		fprintf(f, format1, title, o->start + sa->incl_s + pa->first_base_index, o->length, o->temp, o->temp_var, o->gc_content, o->self_any, o->self_end);
		fprintf(f, "  %5.2f ",  o->hairpin_th);
	}else{
		if(pa->thermodynamic_alignment==1) {
			format1 = "%-25s %5d %4d %7.2f %7.2f   %5.2f  %5.2f ";
		} else {
			format1 = "%-25s %5d %4d %7.2f %7.2f %5.2f %5.2f ";
		}
		fprintf(f, format1, title, o->start + sa->incl_s + pa->first_base_index, o->length, o->temp, o->gc_content, o->self_any, o->self_end);
		fprintf(f, "  %5.2f ",  o->hairpin_th);
	}
	if (print_lib_sim) {
		if (seqlib != NULL) fprintf(f, "%5.2f ",  o->repeat_sim.score[o->repeat_sim.max]);
		  else fprintf(f, "%5s ", "");
	}
	fprintf(f, "%s\n", seq);
	if (PR_DEFAULT_INSIDE_PENALTY != pa->inside_penalty || PR_DEFAULT_OUTSIDE_PENALTY != pa->outside_penalty)
	  fprintf(f, "POSITION PENALTY, QUALITY: %f, %f\n", o->position_penalty, o->quality);
}

static void
print_pair_array(FILE *f, const char* title, int num,
                 const interval_array_t array,
                 const p3_global_settings *pa,
                 const seq_args *sa)
{
    int j;
    if (num > 0) {
        fprintf(f, "%s (position,size):", title);
        for (j = 0; j < num; j++)
            fprintf(f, " %d,%d", 
                    array[j][0] + pa->first_base_index + sa->incl_s,
                    array[j][1]);
        fprintf(f, "\n");
    }
}


static void
print_2_pair_array(FILE *f, const char* title, int num,
		   const interval_array_t left_array,
		   const interval_array_t right_array,
		   const p3_global_settings *pa,
		   const seq_args *sa)
{
	int j;
	if (num > 0) {
		fprintf(f, "%s (left_pos, left_len, right_pos, right_len):", title);
		for (j = 0; j < num; j++) {
			if ((left_array[j][0] == -1) && (left_array[j][1] == -1))
				fprintf(f, " ,,");
			else 
				fprintf(f, " %d,%d,", left_array[j][0] + pa->first_base_index + sa->incl_s, left_array[j][1]);
			if ((right_array[j][0] == -1) && (right_array[j][1] == -1))
				fprintf(f, ",");
			 else 
			 	fprintf(f, "%d,%d", right_array[j][0] + pa->first_base_index + sa->incl_s, right_array[j][1]);
		}
		fprintf(f, "\n");
	}
}

#define VECTOR           (1<<0)
#define LEFT_OLIGO       (1<<1)
#define RIGHT_OLIGO      (1<<2)
#define INTL_OLIGO       (1<<3)
#define TARGET           (1<<4)
#define EXCL_REGION      (1<<5)
#define INTL_EXCL_REGION (1<<6)

/* Prints out the asci picture of the sequence */
/* Return 1 on ENOMEM. Otherwise return 0. */
static int
print_seq(FILE *f,
    const p3_global_settings *pa,
    const seq_args *sa,
    const p3retval *retval,
    primer_rec *h,
    const pair_array_t *best_pairs,
    int num)  /* The number of primer pair to print. */
{
	primer_rec *h2 = NULL;
	primer_rec *h3 = NULL;
	int len, i, j, start;
	int something_found[2] = {0,0};
	int vector_found = 0;
	int *notes;
	char *notestr[2];
	primer_pair *p;
	p = NULL;

	/* the number of variation sequences - 20121207 N.Kasahara */
	int var_num;

	if(retval->output_type == primer_pairs) {
		p = best_pairs->pairs + num;
	}
	len = strlen(sa->sequence);
	if (!(notes = (int*) malloc(sizeof(*notes) * len))) return 1;
	memset(notes, 0, sizeof(*notes) * len);
	if (!(notestr[0] = (char*) malloc(len + 1))) { free(notes); return 1; }
	if (!(notestr[1] = (char*) malloc(len + 1))) { free(notes); free(notestr[0]); return 1; }

	for (i = 0; i < len; i++) {
		if (i < sa->incl_s || i >= sa->incl_s + sa->incl_l)
			notes[i] |= VECTOR;

		if (retval->output_type == primer_pairs && 
			best_pairs->num_pairs > 0) {
			if (i >= p->left->start + sa->incl_s
			    && i < p->left->start + p->left->length + sa->incl_s)
				notes[i] |= LEFT_OLIGO;
			if (i >= p->right->start - p->right->length + 1 + sa->incl_s
			    && i <= p->right->start + sa->incl_s){
				notes[i] |= RIGHT_OLIGO;
			}

			if(pa->pick_internal_oligo==1){
				/* forward internal oligo */
				if (p->intl->oligo_dir==0 &&  i >= p->intl->start + sa->incl_s 
				    && i < p->intl->start + p->intl->length + sa->incl_s){
					notes[i] |= INTL_OLIGO;
				}
				/* reverse internal oligo */
				if(p->intl->oligo_dir==1 && i>=p->intl->start-p->intl->length+1+sa->incl_s
				    && i<=p->intl->start+sa->incl_s){
					notes[i]|=INTL_OLIGO;
				}
			}
		}
		else if (h != NULL) {

			if(pa->pick_left_primer == 1 &&
			    i < retval->fwd.oligo->start + retval->fwd.oligo->length + sa->incl_s &&
			    i >= retval->fwd.oligo->start + sa->incl_s)
				notes[i] |= LEFT_OLIGO;
			if(pa->pick_right_primer == 1 &&
			    i>=retval->rev.oligo->start - retval->rev.oligo->length + 1 + sa->incl_s &&
			    i<=retval->rev.oligo->start + sa->incl_s){
				notes[i]|=RIGHT_OLIGO;
			}

			/* forward internal oligo */
			if(pa->pick_internal_oligo == 1 && retval->intl.oligo->oligo_dir == 0 &&
			    i >= retval->intl.oligo->start + sa->incl_s                &&
			    i < retval->intl.oligo->start + retval->intl.oligo->length + sa->incl_s){
				notes[i] |= INTL_OLIGO;
			}
			/* reverse internal oligo */
			else if(pa->pick_internal_oligo == 1 && retval->intl.oligo->oligo_dir == 1 &&
			    i >= retval->intl.oligo->start-retval->intl.oligo->length + 1 + sa->incl_s &&
			    i <= retval->intl.oligo->start+sa->incl_s){
				notes[i]|=INTL_OLIGO;
			}
		} else if (pa->primer_task == pick_sequencing_primers) {
			if (pa->pick_right_primer && &retval->rev != NULL && retval->rev.num_elem > 0){
				h2 = retval->rev.oligo;
				for (j = 0; j < pa->num_return; j++) {
					if(j > retval->rev.num_elem -1) break;
					h3 = h2 + j;
					if(i >= h3->start - h3->length + 1 + sa->incl_s
					    && i <= h3->start + sa->incl_s) notes[i] |= RIGHT_OLIGO;
				}
			}
			if (pa->pick_left_primer && &retval->fwd != NULL && retval->fwd.num_elem > 0){
				h2 = retval->fwd.oligo;
				for (j = 0; j < pa->num_return; j++) {
					if(j > retval->fwd.num_elem -1) break;
					h3 = h2 + j;
					if(i < h3->start + h3->length + sa->incl_s && i >= h3->start + sa->incl_s) {
						notes[i] |= LEFT_OLIGO;
					}
				}
			}
		}

		/* Adjusting printing difference TARGET and Variation - 20121206 N.Kasahara */
		if(sa->tar2.genotyping==0){
			for (j = 0; j < sa->tar2.count; j++) {
				start = sa->tar2.pairs[j][0] + sa->incl_s;
				if (i >= start && i < start + sa->tar2.pairs[j][1])
					notes[i] |= TARGET;
			}
		}
		if(sa->tar2.genotyping==1){
			var_num=sa->tar2.VAR_sequence_number;
			if(retval->intl.oligo->oligo_dir==0){
				j=sa->tar2.VAR_sequence_number;
				start=sa->tar2.VAR_start[0]-1;
				if(j==0){
					if(i>=start && i<start+sa->tar2.char_num[var_num]){
						notes[i]|=TARGET;
					}
				}else if (j==sa->tar2.char_count){
					if(i>=start && i<start+sa->tar2.char_num[var_num]-1){
						notes[i]|=TARGET;
					}
				}else if(j!=0 && j!=sa->tar2.char_count){
					if(i>=start && i<start+sa->tar2.char_num[var_num]){
						notes[i]|=TARGET;
					}
				}
			} else if(retval->intl.oligo->oligo_dir==1){
				j=sa->tar2.VAR_sequence_number;
				start=sa->tar2.VAR_start[0]-1;
				if(j==0){
					if(i>=start && i<start+sa->tar2.char_num[var_num]){
						notes[i]|=TARGET;
					}
				}else if (j==sa->tar2.char_count){
					if(i>=start && i<start+sa->tar2.char_num[var_num]-1){
						notes[i]|=TARGET;
					}
				}else if(j!=0 && j!=sa->tar2.char_count){
					if(i>=start && i<start+sa->tar2.char_num[var_num]){
						notes[i]|=TARGET;
					}
				}
			}
		}
		for (j = 0; j < sa->excl2.count; j++) {
			start = sa->excl2.pairs[j][0] + sa->incl_s;
			if (i >= start && i < start + sa->excl2.pairs[j][1])
				notes[i] |= EXCL_REGION;
		}
		for (j = 0; j < sa->excl_internal2.count; j++) {
			start = sa->excl_internal2.pairs[j][0] + sa->incl_s;
			if (i >= start && i < start + sa->excl_internal2.pairs[j][1])
				notes[i] |= INTL_EXCL_REGION;
		}
	}
	for (i = 0; i < len; i++) {
		// Initialization
		notestr[0][i] = notestr[1][i] = ' ';

		// Upper Notes
		if (notes[i] & VECTOR) {
			vector_found = 1;
			notestr[0][i] = '.';
		} else if (notes[i] & EXCL_REGION) {
			notestr[0][i] = 'X';
		} else if (pa->pick_left_primer==1 && (notes[i] & LEFT_OLIGO)){
			notestr[0][i] = '>';
		}else if (pa->pick_right_primer==1 && (notes[i] & RIGHT_OLIGO)){
			notestr[0][i] = '<';
		} else if ((pa->primer_task == pick_sequencing_primers) && (notes[i] & LEFT_OLIGO)) {
			notestr[0][i] = '>';
		} else if ((pa->primer_task == pick_sequencing_primers) && (notes[i] & RIGHT_OLIGO)) {
			notestr[0][i] = '<';
		}

		// Lower Notes
		if (notes[i] & INTL_EXCL_REGION){
			notestr[1][i] = 'x';
		} else if ((pa->primer_task == pick_sequencing_primers) && (notes[i] & LEFT_OLIGO) && (notes[i] & RIGHT_OLIGO)){
			notestr[1][i] = '^';
		} else if (pa->pick_internal_oligo==1 && (notes[i] & INTL_OLIGO)){
			notestr[1][i] = '^';
		}

		// Target
		if (notes[i] & TARGET){
			notestr[1][i] = '*';
		}

	}
	for(j=0; j<2; ++j){
		for (i = 0; i < len; i++) {
			if (notestr[j][i] != ' '){
				something_found[j] = 1;
				break;
			}
		}
	}

	print_seq_lines(f, sa->sequence, notestr, len, 60, something_found, pa);

//    if (something_found)
	if (something_found[0] || something_found[1])
		fprintf(f, "KEYS (in order of precedence):\n");

	if (vector_found)
		fprintf(f, "...... vector sequence\n");

	if (sa->excl2.count > 0)
		fprintf(f, "XXXXXX excluded region\n");

	if (pa->pick_internal_oligo ==1 && sa->excl_internal2.count > 0)
		fprintf(f, "xxxxxx excluded region for internal probe\n");

	if (sa->tar2.count > 0 && sa->tar2.genotyping==0) 
		fprintf(f, "****** target\n");
	else if(sa->tar2.count>0 && sa->tar2.genotyping==1) 
		fprintf(f, "****** TARGET VARIANT\n");

	if (retval->output_type == primer_pairs && best_pairs->num_pairs > 0) {
		fprintf(f, ">>>>>> left primer\n");
		fprintf(f, "<<<<<< right primer\n");
		if ( pa->pick_internal_oligo == 1 ) 
			fprintf(f, "^^^^^^ internal probe\n");
	} else if (pa->primer_task == pick_sequencing_primers) {
		fprintf(f, ">>>>>> left primer\n");
		fprintf(f, "<<<<<< right primer\n");
		fprintf(f, "^^^^^^ left primer / right primer overlap\n");
	}

	if (pa->pick_left_primer == 1 && h != NULL) 
		fprintf(f, ">>>>>> left primer\n");
	if (pa->pick_right_primer == 1 && h != NULL) 
		fprintf(f, "<<<<<< right primer\n");
	if (pa->pick_internal_oligo == 1 && h != NULL) 
		fprintf(f, "^^^^^^ internal probe\n");

	if (something_found[0] | something_found[1]) fputc('\n', f);
	free(notes);
	free(notestr[0]);
	free(notestr[1]);
	return 0;
}


static void
print_seq_lines(FILE *f, char *s, char *n[2], int seq_size,
                int line_size, int something_found[2],
                const p3_global_settings *pa)
{
    char* n0 = n[0];
    char* n1 = n[1];
    int i = 0;

    while (seq_size > line_size) {
        /* The first line */
        if (something_found[0]) {
            fprintf(f, "      ");
            fwrite(n0, sizeof(*n0), line_size, f);
            fprintf(f, "\n");
        }

        /* The second line */
        fprintf(f, "%5d ", i + pa->first_base_index);
        fwrite(s, sizeof(*s), line_size, f);
        fputc('\n', f);

        /* The third line */
        if (something_found[1]) {
            fprintf(f, "      ");
            fwrite(n1, sizeof(*n1), line_size, f);
            fprintf(f, "\n");
        }

        /* spacing */
        fputc('\n', f);

        seq_size -= line_size;
        s += line_size;
        n0 += line_size;
        n1 += line_size;
        i += line_size;
    }

    /* Tails */
    /* The first line */
    if (something_found[0]){
        *(n0+seq_size) = NULL; /* add NULL to the end of the sequence - 20130117 N.Kasahara */
        fprintf(f, "      %s\n", n0);
    }
    /* The second line */
    fprintf(f, "%5d %s\n", i + pa->first_base_index, s);

    /* The third line */
    if (something_found[1]){
        *(n1+seq_size) = NULL; /* add NULL to the end of the sequence - 20130117 N.Kasahara */
        fprintf(f, "      %s\n\n\n", n1);
    }
}

static void
print_pair_info(FILE *f, const primer_pair *p, const p3_global_settings *pa)
{
	double pi_pair_compl_any;
	double pi_pair_compl_end;

	fprintf(f, "product size: %d,  ", p->product_size);
	if(pa->thermodynamic_alignment==0)
		fprintf(f, "Primer pair compl any: %.2f, 3': %.2f", p->compl_any, p->compl_end);
	else 
		fprintf(f, "Primer pair compl any_th: %.2f, 3'_th: %.2f", p->compl_any, p->compl_end);

	if(pa->pick_internal_oligo==1){
		/* output Primer-Probe pair complementarity is larger one between LEFT to INTERNAL and RIGHT to INTERNAL - 20150707 Y.Kimura */
		if(p->li_pair_compl_any > p->ri_pair_compl_any){
			pi_pair_compl_any = p->li_pair_compl_any;
		}else{
			pi_pair_compl_any = p->ri_pair_compl_any;
		}
		if(p->li_pair_compl_end > p->ri_pair_compl_end){
			pi_pair_compl_end = p->li_pair_compl_end;
		}else{
			pi_pair_compl_end = p->ri_pair_compl_end;
		}
		fprintf(f, ",  LEFT Primer-Probe compl ");
		if(pa->thermodynamic_alignment==0)
			fprintf(f, "any: %.2f, 3': %.2f", p->li_pair_compl_any, p->li_pair_compl_end);
			//fprintf(f, "any: %.2f, 3': %.2f\n", pi_pair_compl_any, pi_pair_compl_end);
		else 
			fprintf(f, "any_th: %.2f, 3'_th: %.2f", p->li_pair_compl_any, p->li_pair_compl_end);
			//fprintf(f, "any_th: %.2f, 3'_th: %.2f\n", pi_pair_compl_any, pi_pair_compl_end);
		fprintf(f, ",  RIGHT Primer-Probe compl ");
		if(pa->thermodynamic_alignment==0)
			fprintf(f, "any: %.2f, 3': %.2f\n", p->ri_pair_compl_any, p->ri_pair_compl_end);
		else 
			fprintf(f, "any_th: %.2f, 3'_th: %.2f\n", p->ri_pair_compl_any, p->ri_pair_compl_end);
	}else{
		fprintf(f, "\n");
	}

	if (pa->product_max_tm != PR_DEFAULT_PRODUCT_MAX_TM || pa->product_min_tm != PR_DEFAULT_PRODUCT_MIN_TM) {
		printf("PRODUCT Tm: %.4f, ", p->product_tm);
		printf("PRODUCT Tm - min(OLIGO Tm): %.4f\n", p->product_tm_oligo_tm_diff);
	}
}


static void
print_rest(FILE *f, const p3_global_settings *pa,
           const seq_args *sa,
           const pair_array_t *best_pairs)
{
    int i;
    int print_lib_sim = lib_sim_specified(pa);


    fprintf(f, "ADDITIONAL OLIGOS\n");
    fprintf(f, "  "); print_oligo_header(f, "", print_lib_sim, pa->thermodynamic_alignment, sa->tar2.genotyping);
    for (i = 1; i < best_pairs->num_pairs; i++) {
        fprintf(f, "\n%2d ", i);

		if(pa->modify_left_primer==1){
			print_oligo_Z(f, TITLE_MOD_LEFT, sa, best_pairs->pairs[i].left, OT_LEFT, pa, pa->p_args.repeat_lib, print_lib_sim);
		} else {
			print_oligo_Z(f, TITLE_LEFT, sa, best_pairs->pairs[i].left, OT_LEFT, pa, pa->p_args.repeat_lib, print_lib_sim);
		}

		fprintf(f, "   ");
		if(pa->modify_right_primer==1){
			print_oligo_Z(f, TITLE_MOD_RIGHT, sa, best_pairs->pairs[i].right, OT_RIGHT, pa, pa->p_args.repeat_lib, print_lib_sim);
		} else {
			print_oligo_Z(f, TITLE_RIGHT, sa, best_pairs->pairs[i].right, OT_RIGHT, pa, pa->p_args.repeat_lib, print_lib_sim);
		}

		if(pa->pick_internal_oligo==1){
			fprintf(f, "   ");
			if(pa->modify_internal_oligo==1){
			/* Modified internal oligo */
				/* printing internal_oligo_direction - 20121204 N.Kasahara */
				if((best_pairs->pairs[i].intl)->oligo_dir==0){
					print_oligo_Z(f, TITLE_MOD_INTL_FW,sa, best_pairs->pairs[i].intl, OT_INTL, pa, pa->o_args.repeat_lib, print_lib_sim);
				} else if((best_pairs->pairs[i].intl)->oligo_dir==1){
					print_oligo_Z(f, TITLE_MOD_INTL_RV,sa, best_pairs->pairs[i].intl, OT_INTL, pa, pa->o_args.repeat_lib, print_lib_sim);
				}
			} else {
			/* Normal internal oligo */
				/* printing internal_oligo_direction - 20121204 N.Kasahara */
				if((best_pairs->pairs[i].intl)->oligo_dir==0){
					print_oligo_Z(f, TITLE_INTL_FW, sa, best_pairs->pairs[i].intl, OT_INTL, pa, pa->o_args.repeat_lib, print_lib_sim);
				} else if((best_pairs->pairs[i].intl)->oligo_dir==1){
					print_oligo_Z(f, TITLE_INTL_RV, sa, best_pairs->pairs[i].intl, OT_INTL, pa, pa->o_args.repeat_lib, print_lib_sim);
				}
			}
		}
		if (best_pairs->pairs[i].product_size > 0) {
			fprintf(f, "   ");
			print_pair_info(f, &best_pairs->pairs[i], pa);
			/* Print penalty score - Y.Kimura */
			fprintf(f, "   ");
			fprintf(f, "penalty score: %.4f\n", best_pairs->pairs[i].pair_quality);
		}
	}
}

/* Print out the statistics of primer picking */
/* This function does _not_ print out the no_orf statistic. */
static void
print_explain(FILE *f,
	      const p3_global_settings *pa,
	      const seq_args *sa,
	      const p3retval *retval,
	      int print_lib_sim,
	      const char *pr_release)
{
   
   
  const char *format;  /* Format string for the table headers */
  if(pa->thermodynamic_alignment==0)
    format    = "%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s";
  else format = "%6s%6s%6s%6s%6s%6s%6s%6s%6s%6s%7s%6s";
  
  fprintf(f, "\nStatistics\n");
   
  if (!pa->pick_anyway || !(
			    (pa->pick_left_primer == 1 && pa->pick_internal_oligo == 0
			     && pa->pick_right_primer == 1
			     && sa->left_input && sa->right_input)
			    || (pa->pick_left_primer == 1 && pa->pick_internal_oligo == 1
				&& pa->pick_right_primer == 1
				&& sa->left_input && sa->right_input && sa->internal_input)
			    || (pa->pick_left_primer == 1 && pa->pick_internal_oligo == 0
				&& pa->pick_right_primer == 0
				&& sa->left_input)
			    || (pa->pick_left_primer == 0 && pa->pick_internal_oligo == 0
				&& pa->pick_right_primer == 1
				&& sa->right_input)
			    || (pa->pick_left_primer == 0 && pa->pick_internal_oligo == 1
				&& pa->pick_right_primer == 0
				&& sa->internal_input))) {
      /* print first row of header */
      if(pa->thermodynamic_alignment==0)
		fprintf(f, format, "", "con", "too",  "in",  "in",  "not", "", "no", "tm",  "tm",  "   high", "high");
      else 
		fprintf(f, format, "", "con", "too",  "in",  "in",  "not", "", "no", "tm",  "tm",  "   high", "high");
//      if(pa->thermodynamic_alignment==1) 
		fprintf(f, "%6s", "high");
      if(print_lib_sim) 
		fprintf(f, "%6s", "high");
      fprintf(f, "%6s%6s", "", "high");
      if(pa->lowercase_masking) 
		fprintf(f, "%6s", "lower");
      fprintf(f, "%6s\n", "");
      
      /* print second row of header */
      if(pa->thermodynamic_alignment==0) {
		fprintf(f, format, "", "sid", "many", "tar", "excl", "ok", "bad","GC", "too", "too", "    any",  "3'");
		fprintf(f, "%6s", "hair-");
      } else {
		fprintf(f, format, "", "sid", "many", "tar", "excl", "ok", "bad","GC", "too", "too", " any_th",  "3'_th");
		fprintf(f, "%6s", "hair-");
      }
	   
      if(print_lib_sim) 
	fprintf(f, "%6s", "lib");
      fprintf(f, "%6s%6s", "poly", "end");
      if(pa->lowercase_masking) 
	fprintf(f, "%6s", " case");
      fprintf(f, "%6s\n", "");
      
      /* print third row of header */
      if(pa->thermodynamic_alignment==0) 
	fprintf(f, format, "", "ered","Ns",   "get", "reg",  "reg", "GC%", "clamp", "low", "high","  compl", "compl");
      else
	fprintf(f, format, "", "ered","Ns",   "get", "reg",  "reg", "GC%", "clamp", "low", "high","  compl", "compl");
//      if(pa->thermodynamic_alignment==1) 
	fprintf(f, "%6s", " pin");
      if(print_lib_sim) 
	fprintf(f, "%6s", "sim");
      fprintf(f, "%6s%6s", "X", "stab");
      if(pa->lowercase_masking) 
	fprintf(f, "%6s", " end");
      if(pa->lowercase_masking)
	fprintf(f, "%6s\n", "ok  ");
      else
	fprintf(f, "%6s\n", "ok");
  }
  if (pa->pick_left_primer == 1
      && !(pa->pick_anyway && sa->left_input))
    print_stat_line(f, "Left", retval->fwd.expl, 
		    print_lib_sim, pa->lowercase_masking, pa->thermodynamic_alignment);
  
  if (pa->pick_right_primer == 1
      && !(pa->pick_anyway && sa->right_input))
    print_stat_line(f, "Right", retval->rev.expl,
		    print_lib_sim, pa->lowercase_masking, pa->thermodynamic_alignment);
  
  if (pa->pick_internal_oligo == 1
      && !(pa->pick_anyway && sa->internal_input))
    print_stat_line(f, "Intl", retval->intl.expl, 
		    print_lib_sim, pa->lowercase_masking, pa->thermodynamic_alignment);
  
  if (pa->pick_left_primer == 1 && pa->pick_right_primer == 1) {
    fprintf(f, "Pair Stats:\n%s\n",
	    p3_get_pair_array_explain_string(p3_get_rv_best_pairs(retval)));
  }
}

static void
print_stat_line(FILE *f, const char *t, oligo_stats s,
                int print_lib_sim, int lowercase_masking,
                int thermodynamic_alignment)
{
   const char *format = "%-6s%6d%6d%6d%6d%6d%6d%6d%6d%6d";
   if (!strcmp(t, "Left"))
     fprintf(f,format, t, s.considered, s.ns, s.target, s.excluded,
	     s.not_in_any_left_ok_region, s.gc, s.gc_clamp, s.temp_min, s.temp_max);
   else if (!strcmp(t, "Right"))
	  fprintf(f,format, t, s.considered, s.ns, s.target, s.excluded,
		  s.not_in_any_right_ok_region, s.gc, s.gc_clamp, s.temp_min, s.temp_max);
   else fprintf(f,format, t, s.considered, s.ns, s.target, s.excluded,
		0, s.gc, s.gc_clamp, s.temp_min, s.temp_max);
   if(thermodynamic_alignment) {
      fprintf(f, " %6d%6d%6d", s.compl_any, s.compl_end, s.hairpin_th);
   } else {
      fprintf(f, " %6d%6d%6d", s.compl_any, s.compl_end, s.hairpin_th);
      //fprintf(f, "%6d%6d", s.compl_any, s.compl_end);
   }
   if (print_lib_sim) fprintf(f, "%6d",s.repeat_score);
   fprintf(f, "%6d%6d",s.poly_x, s.stability);
   if (lowercase_masking) fprintf(f, "%6d",s.gmasked);
   fprintf(f, "%6d\n", s.ok);
}


/* 
 * Return true iff a check for library similarity has been specified for
 * either the primer pair or the internal probe.
 */
static int
lib_sim_specified(const p3_global_settings *pa) {
  return (pa->p_args.repeat_lib || pa->o_args.repeat_lib);
}

void
format_error(FILE *f, const char* seq_name, const char *err)
{
  if (NULL != seq_name)
    fprintf(f, "PRIMER PICKING RESULTS FOR %s\n\n", seq_name);
  if (err != NULL) 
    fprintf(f, "INPUT PROBLEM: %s\n\n", err);
}

void
format_warning(FILE *f, const char* seq_name, const char *err)
{
  if (NULL != seq_name)
    fprintf(f, "WARNINGS FOR %s\n\n", seq_name);
  if (err != NULL) 
    fprintf(f, "INPUT PROBLEM: %s\n\n", err);
}

/* Format and print out one oligo */
/* add function for printing internal_oligo_direction and modification
   20121204 - N.Kasahara
   20140109 - Y.Kimura   */

static void 
format_oligos(FILE *f,
              const p3_global_settings *pa,
              const seq_args    *sa,
              const p3retval *retval,
              const char* pr_release,
              const pr_append_str *combined_retval_err,
              int explain_flag)
{
	char *warning;
	int print_lib_sim = lib_sim_specified(pa);
	int i;
	int print_primers = 0;
	primer_rec  *h = NULL;
	pair_array_t *best_pairs;
	primer_rec *p;
	int rest_count = 0;

	PR_ASSERT(NULL != f);
	PR_ASSERT(NULL != pa);
	PR_ASSERT(NULL != sa);

	best_pairs = NULL;

	if (!pr_is_empty(combined_retval_err)) {
		format_error(f, sa->sequence_name, pr_append_str_chars(combined_retval_err));
		return;
	}

	if (NULL != sa->sequence_name)
		fprintf(f, "PRIMER PICKING RESULTS FOR %s\n\n", sa->sequence_name);
	if (pa->pick_left_primer || pa->pick_right_primer) {
		if (pa->p_args.repeat_lib != NULL)
			fprintf(f, "Using mispriming library %s\n",
			        pa->p_args.repeat_lib->repeat_file);
		else
			fprintf(f, "No mispriming library specified\n");
	} 
	if (pa->pick_internal_oligo) {
		if (pa->o_args.repeat_lib != NULL)
			fprintf(f, "Using internal probe mishyb library %s\n", pa->o_args.repeat_lib->repeat_file);
		else
			fprintf(f, "No internal probe mishyb library specified\n");
	}
	fprintf(f, "Using %d-based sequence positions\n", pa->first_base_index);

	if (pa->pick_left_primer) {
		if (retval->fwd.num_elem == 0){
			fprintf(f, "NO LEFT PRIMER FOUND\n\n");
		} else {
			print_primers = 1;
		}
	}
	if (pa->pick_internal_oligo) {
		if (retval->intl.num_elem == 0){
			fprintf(f, "NO INTERNAL PROBE FOUND\n\n");
		} else {
			print_primers = 1;
		}
	}
	if (pa->pick_right_primer) {
		if (retval->rev.num_elem == 0){
			fprintf(f, "NO RIGHT PRIMER FOUND\n\n");
		} else {
			print_primers = 1;
	}
	}
	if ((warning = p3_get_rv_and_gs_warnings(retval, pa)) != NULL) {
		fprintf(f, "\nWARNING: %s\n", warning);
		free(warning);
	}
	if ((pa->primer_task != pick_primer_list ) && (pa->primer_task != pick_sequencing_primers)) { 
		if (print_primers == 1) { 
			print_oligo_header(f, "OLIGO", print_lib_sim, pa->thermodynamic_alignment, sa->tar2.genotyping);
		}

		/* Print out the first line with the best primers */
		if ((pa->pick_left_primer) && (&retval->fwd != NULL ) && (retval->fwd.num_elem > 0)){
			if(pa->modify_left_primer==1){
				print_oligo_Z(f, TITLE_MOD_LEFT, sa, retval->fwd.oligo, OT_LEFT, pa, pa->p_args.repeat_lib, print_lib_sim);
			} else {
				print_oligo_Z(f, TITLE_LEFT, sa, retval->fwd.oligo, OT_LEFT, pa, pa->p_args.repeat_lib, print_lib_sim);
			}
			h = retval->fwd.oligo;
			rest_count = 1;
		}

		if ((pa->pick_internal_oligo) && (&retval->intl != NULL ) && (retval->intl.num_elem > 0)){
			if(pa->modify_internal_oligo==1){
				if((retval->intl.oligo)->oligo_dir==0){
					print_oligo_Z(f, TITLE_MOD_INTL_FW, sa, retval->intl.oligo, OT_INTL, pa, pa->p_args.repeat_lib, print_lib_sim);
				}else if((retval->intl.oligo)->oligo_dir==1){
					print_oligo_Z(f, TITLE_MOD_INTL_RV, sa, retval->intl.oligo, OT_INTL, pa, pa->p_args.repeat_lib, print_lib_sim);
				}
			} else {
				if((retval->intl.oligo)->oligo_dir==0){
					print_oligo_Z(f, TITLE_INTL_FW, sa, retval->intl.oligo, OT_INTL, pa, pa->p_args.repeat_lib, print_lib_sim);
				}else if((retval->intl.oligo)->oligo_dir==1){
					print_oligo_Z(f, TITLE_INTL_RV, sa, retval->intl.oligo, OT_INTL, pa, pa->p_args.repeat_lib, print_lib_sim);
				}
			}
			h = retval->intl.oligo;
			rest_count = 1;
		}

		if ((pa->pick_right_primer) && (&retval->rev != NULL ) && (retval->rev.num_elem > 0)){
			if(pa->modify_right_primer==1){
				print_oligo_Z(f, TITLE_MOD_RIGHT, sa, retval->rev.oligo, OT_RIGHT, pa, pa->p_args.repeat_lib, print_lib_sim);
			} else {
				print_oligo_Z(f, TITLE_RIGHT, sa, retval->rev.oligo, OT_RIGHT, pa, pa->p_args.repeat_lib, print_lib_sim);
			}
			h = retval->rev.oligo;
			rest_count = 1;
		}
	}

	if(print_primers == 1) { 
		fprintf(f, "SEQUENCE SIZE: %ld\n", (long int) strlen(sa->sequence));
		fprintf(f, "INCLUDED REGION SIZE: %d\n\n", sa->incl_l);

		if(sa->tar2.genotyping==0){
			print_pair_array(f, "TARGETS", sa->tar2.count, sa->tar2.pairs, pa, sa);
		} else if(sa->tar2.genotyping==1){
			fprintf(f,"TARGET VARIANT * (position,size):");
			if(sa->tar2.char_count==sa->tar2.VAR_sequence_number){
				fprintf(f," %d,%d",sa->tar2.VAR_start[0],sa->tar2.char_num[sa->tar2.VAR_sequence_number]-1);
			} else {
				fprintf(f," %d,%d",sa->tar2.VAR_start[0],sa->tar2.char_num[sa->tar2.VAR_sequence_number]);
			}
			fprintf(f,"\n\n");
		}

		print_pair_array(f, "EXCLUDED REGIONS", sa->excl2.count, sa->excl2.pairs, pa, sa);
		print_pair_array(f, "INTERNAL PROBE EXCLUDED REGIONS", sa->excl_internal2.count, sa->excl_internal2.pairs, pa, sa);
		print_2_pair_array(f, "PAIR_OK_REGIONS", sa->ok_regions.count, sa->ok_regions.left_pairs, sa->ok_regions.right_pairs, pa, sa);
	}

	if (pa->primer_task != pick_primer_list ) {
		if(pa->pick_internal_oligo==1){
			if(print_seq(f,pa,sa,retval,retval->intl.oligo,best_pairs,0)) exit(-2);
		}
		else if(pa->pick_left_primer==1){
			if(print_seq(f,pa,sa,retval,retval->fwd.oligo,best_pairs,0)) exit(-2);
		}
		else if(pa->pick_right_primer==1){
			if(print_seq(f,pa,sa,retval,retval->rev.oligo,best_pairs,0)) exit(-2);
		} 
		else if (print_seq(f, pa, sa, retval, h, best_pairs, 0)) exit(-2); /* ENOMEM */
	}
	fprintf(f, "\n");

	/* Print out the other primers */
	if ((pa->pick_left_primer) && (&retval->fwd != NULL ) && (retval->fwd.num_elem > rest_count)){
		int n = retval->fwd.num_elem;
		h = retval->fwd.oligo;
		if (rest_count == 1) {  
			fprintf(f, "ADDITIONAL OLIGOS\n");
		}
		fprintf(f, "   "); print_oligo_header(f, "", print_lib_sim, pa->thermodynamic_alignment, sa->tar2.genotyping);
		for (i = rest_count; i < pa->num_return; i++) {
			if(i > n-1) break;
			p = h + i;
			fprintf(f, "%2d ", i + 1 - rest_count);
			if(pa->modify_left_primer==1){
				print_oligo_Z(f, TITLE_MOD_LEFT, sa, p, OT_LEFT, pa, pa->p_args.repeat_lib, print_lib_sim);
			} else {
				print_oligo_Z(f, TITLE_LEFT, sa, p, OT_LEFT, pa, pa->p_args.repeat_lib, print_lib_sim);
			}
		}
		if (rest_count == 0) {  
			fprintf(f, "\n ");
		}
	}

	if ((pa->pick_internal_oligo) && (&retval->intl != NULL ) && (retval->intl.num_elem > rest_count)){
		int n = retval->intl.num_elem;
		h = retval->intl.oligo;  

		if (rest_count == 1) {  
			fprintf(f, "ADDITIONAL OLIGOS\n");
		}
		fprintf(f, "   "); print_oligo_header(f, "", print_lib_sim, pa->thermodynamic_alignment, sa->tar2.genotyping);
		for (i = rest_count; i < pa->num_return; i++) {
			if(i > n-1) break;
			p = h + i;
			fprintf(f, "%2d ", i + 1 - rest_count);

			if(pa->modify_internal_oligo==1){
				if(p->oligo_dir==0){
					print_oligo_Z(f, TITLE_MOD_INTL_FW, sa, p, OT_INTL, pa, pa->p_args.repeat_lib, print_lib_sim);
				}else if(p->oligo_dir==1){
					print_oligo_Z(f, TITLE_MOD_INTL_RV, sa, p, OT_INTL, pa, pa->p_args.repeat_lib, print_lib_sim);
				}
			} else {
				if(p->oligo_dir==0){
					print_oligo_Z(f, TITLE_INTL_FW, sa, p, OT_INTL, pa, pa->p_args.repeat_lib, print_lib_sim);
				}else if(p->oligo_dir==1){
					print_oligo_Z(f, TITLE_INTL_RV, sa, p, OT_INTL, pa, pa->p_args.repeat_lib, print_lib_sim);
				}
			}
		}
		if (rest_count == 0) {  
			fprintf(f, "\n ");
		}
	}

	if ((pa->pick_right_primer) && (&retval->rev != NULL ) && (retval->rev.num_elem > rest_count)) {
		int n = retval->rev.num_elem;
		h = retval->rev.oligo; 
		if (rest_count == 1) {  
			fprintf(f, "ADDITIONAL OLIGOS\n");
		}
		fprintf(f, "   "); print_oligo_header(f, "", print_lib_sim, pa->thermodynamic_alignment, sa->tar2.genotyping);
		for (i = rest_count; i < pa->num_return; i++) {
			if(i > n-1) break;
			p = h + i;
			fprintf(f, "%2d ", i + 1 - rest_count);
			if(pa->modify_right_primer==1){
				print_oligo_Z(f, TITLE_MOD_RIGHT, sa, p, OT_RIGHT, pa, pa->p_args.repeat_lib, print_lib_sim);
			} else {
				print_oligo_Z(f, TITLE_RIGHT, sa, p, OT_RIGHT, pa, pa->p_args.repeat_lib, print_lib_sim);
			}
		}
	}
	if (explain_flag) 
		print_explain(f, pa, sa, retval, print_lib_sim, pr_release);
		fprintf(f, "\n\n");
	if (fflush(f) == EOF) {
		perror("fflush(f) failed");
		exit(-1);
	}
}

