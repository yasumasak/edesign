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

#include <limits.h>
#include <stdlib.h>  /* strtod, strtol,... */
#include <ctype.h> /* toupper, isspace */
#include <string.h> /* memset, strlen,  strcmp, ... */
#include <float.h>
#include "read_boulder_Z.h"
#include "modification_Z.h"

#define INIT_BUF_SIZE 1024
#define INIT_LIB_SIZE  500
#define PR_MAX_LIBRARY_WT 100.0

char *thermodynamic_params_path = NULL;
int   thermodynamic_path_changed = 1;

/* Static functions. */
static void *_rb_safe_malloc(size_t x);
static void  out_of_memory_error();

static void   parse_double(const char *, const char *, double *,
                           pr_append_str *);

static void   parse_int(const char *, const char *, int *, pr_append_str *);

static const char *parse_int_pair(const char *, const char *,
                                  char, int *, int *,
                                  pr_append_str *);

static char  *parse_2_int_pair(const char*, char*, char,
                               char, int*, int*, int*, int*,
                               pr_append_str*); 

static void   parse_interval_list(const char *tag_name,
                                  const char *datum,
                                  interval_array_t2 *interval_arr,
                                  pr_append_str *err);

/* Get variation position to new elements of data structure
   - 20121217 N.Kasahara */
static void   parse_interval_list_VAR(const char *tag_name,
                                    const char  *datum,
                                    interval_array_t2 *interval_arr,
                                    pr_append_str *err);


static void   parse_2_interval_list(const char *tag_name,
                                    char *datum,
                                    interval_array_t4 *interval_arr,
                                    pr_append_str *err);   

static int    parse_intron_list(char *, int *, int *);

static void   parse_product_size(const char *, char *, p3_global_settings *,
                                 pr_append_str *);

static void   pr_append(pr_append_str *, const char *);
static void   pr_append_new_chunk(pr_append_str *x, const char *s);

static void   tag_syntax_error(const char *, const char *,  pr_append_str *);
static int    parse_seq_quality(char *, seq_args *);

/* 
 * Hack to support old SunOS headers.  (We do not try to declare _all_
 * undeclared functions; only those with non-int return types.)
 */
#ifndef __cplusplus
extern double strtod();
#endif

/* 
 * See read_boulder.h for description.
 */
#define COMPARE(TAG) (!strncmp(s, TAG, tag_len) \
                      && ('=' == s[tag_len] || ' ' == s[tag_len]) \
                      && '\0' == TAG[tag_len])

#define COMPARE_AND_MALLOC(TAG,T)                  \
   if (COMPARE(TAG)) {                             \
       if (T) {                                    \
           pr_append_new_chunk(parse_err,          \
                               "Duplicate tag: "); \
           pr_append(parse_err, TAG);              \
       } else {                                    \
           T = (char*) _rb_safe_malloc(datum_len + 1); \
           strcpy(T, datum);                       \
       }                                           \
       continue;                                   \
   }

#define COMPARE_FLOAT(TAG,T)                      \
   if (COMPARE(TAG)) {                            \
       parse_double(TAG, datum, &(T), parse_err); \
       continue;                                  \
   }

/* This macro added 2011 08 28 to allow use of 
   of the 'setting' functions, such as 
   p3_set_gs_primer_self_end, in setting
   arguments. At this point, this macro is
   only used twice, and there are no analogous
   macros for integers,etc.
*/
#define COMPARE_FLOAT_USE_FN(TAG, set_function)        \
   if (COMPARE(TAG)) {                                 \
     parse_double(TAG, datum, &tmp_double, parse_err); \
     set_function(pa, tmp_double);                     \
     continue;                                         \
   }

#define COMPARE_INT(TAG,T)                     \
   if (COMPARE(TAG)) {                         \
       parse_int(TAG, datum, &(T), parse_err); \
       continue;                               \
   }

#define COMPARE_INTERVAL_LIST(TAG, PLACE)                  \
   if (COMPARE(TAG)) {                                     \
       parse_interval_list(TAG, datum, PLACE, parse_err);  \
       continue;                                           \
   }

/* COMPARE_INTERVAL_LIST_VAR function
   20140830 -  Y.Kimura */
#define COMPARE_INTERVAL_LIST_VAR(TAG, PLACE)                  \
   if (COMPARE(TAG)) {                                         \
       parse_interval_list_VAR(TAG, datum, PLACE, parse_err);  \
   }

#define COMPARE_2_INTERVAL_LIST(TAG, PLACE)                \
   if (COMPARE(TAG)) {                                     \
       parse_2_interval_list(TAG, datum, PLACE, parse_err);\
       continue;                                           \
   }

/* parameter size definition - N.Kasahara*/
#define INIT_BUF_SIZE 1024

/* c-basic-offset in emacs is set to 2 */
/* 
 * See read_boulder.h for description.
 */
int
read_boulder_record(FILE *file_input,
                    const int *strict_tags,
                    const int *io_version,
                    int   echo_output, /* should be echo_input */
                    const p3_file_type file_type,
                    p3_global_settings *pa, 
                    seq_args *sa, 
                    pr_append_str *glob_err,  /* Really should be called fatal_parse_err */
                    pr_append_str *nonfatal_parse_err,
                    pr_append_str *warnings,
                    read_boulder_record_results *res) 
{
	int line_len;
	int tag_len, datum_len;
	int data_found = 0;
	int pick_internal_oligo = 2;
	char *s, *n, *datum, *task_tmp = NULL;
	const char *p;
	pr_append_str *parse_err;
	pr_append_str *non_fatal_err;
	char *repeat_file_path = NULL, *int_repeat_file_path = NULL;
	int tmp_int;
	/* int min_3_prime = 0, min_5_prime = 0; Removed 10/20/2010 */
	int min_3_prime_distance_global = 0;    /* needed to check if both global and specific*/
	int min_3_prime_distance_specific = 0;  /* are given in same boulder record */
	int min_three_prime_distance;           /* holder for the value of this tag */
	double tmp_double;

	/* variables for checking Z - 20120910 N.Kasahara */
	int i, j, k, m, h = 0;
	int zcount = 0; /* numbers of appearing Z in sa->***input sequence */
	int length = 0; /* length of sa->***input sequence */
	int seq_len = 0; /* length of sa->sequence */

	/* variables for variation type - 20121205 N.Kasahara */
	int start=0;
	non_fatal_err = nonfatal_parse_err;

	sa->left_input = NULL;
	sa->right_input = NULL;
	sa->internal_input = NULL;

	while (((s = p3_read_line(file_input)) != NULL) && (strcmp(s,"="))) {
		/* If we are reading from a settings file, then skip every 
		   line except those begining "PRIMER_" or "P3_FILE_ID".
		   Hint: strncomp returns 0 if both strings are equal */
		if (file_type == settings 
		    && strncmp(s, "PRIMER_", 7) /* "s does not begin with 'PRIMER_'" */
		    && strncmp(s, "P3_FILE_ID", 10) /* "s does not begin with 'P3_FILE_ID'" */
		 ) {
			continue;
		}
		/* Silently ignore all primer3plus tags */
		/* Removed 10/20/10 because P3P_ tags
		   are already ignored in settings files. */
		/* if (!(strncmp(s, "P3P_", 4))) {
		  continue;
		  } */

		data_found = 1;
		/* Print out the input */
		if (echo_output) printf("%s\n", s);
		line_len = strlen(s);

		/* If the line has an "=" read the tag in the right place */
		if ((n=strchr(s,'=')) == NULL) {
		  /* The input line is illegal because it has no
		   * "=" in it, but we still will read to the end
		   * of the record. */
			pr_append_new_chunk(glob_err, "Input line with no '=': ");
			pr_append(glob_err, s);
			continue;
		} 
		/* Read in the new tags used from primer3 version 2.0 on */    
		else {
			/* Get the tag and the value pointers */
			tag_len = n - s;
			datum = n + 1;
			datum_len = line_len - tag_len - 1;

			/* Process "Sequence" (i.e. Per-Record) Arguments". */
			parse_err = non_fatal_err;

			/* COMPARE_AND_MALLOC("SEQUENCE", sa->sequence); */
			if (COMPARE("SEQUENCE_TEMPLATE")) {   /* NEW WAY */
				if (/* p3_get_seq_arg_sequence(sa) */ sa->sequence) {
					pr_append_new_chunk(parse_err, "Duplicate tag: ");
					pr_append(parse_err, "SEQUENCE_TEMPLATE"); 
				} else {
					if (p3_set_sa_sequence(sa, datum)) exit(-2);
				}
				continue;
			}
			if (COMPARE("SEQUENCE_QUALITY")) {
				if ((sa->n_quality = parse_seq_quality(datum, sa)) == 0) {
					pr_append_new_chunk(parse_err, "Error in sequence quality data");
				}
				continue;
			}
			COMPARE_AND_MALLOC("SEQUENCE_ID", sa->sequence_name);
			COMPARE_AND_MALLOC("SEQUENCE_PRIMER", sa->left_input);
			COMPARE_AND_MALLOC("SEQUENCE_PRIMER_REVCOMP", sa->right_input);
			COMPARE_AND_MALLOC("SEQUENCE_INTERNAL_OLIGO", sa->internal_input);
			COMPARE_2_INTERVAL_LIST("SEQUENCE_PRIMER_PAIR_OK_REGION_LIST", &sa->ok_regions);

			/* change for reading Variation position information - 20121205 N.Kasahara */
			/* add function of getting Variation position information - 20121217 - N.Kasahara */
			/* rewirte SEQUENCE_TARGET (add COMPARE_INTERVAL_LIST_VAR function) - 20140830 Y.Kimura*/
			//      COMPARE_INTERVAL_LIST("SEQUENCE_TARGET", &sa->tar2);
			if(strstr(s,"SEQUENCE_TARGET")!=NULL){
				if(sa->tar2.genotyping==1){
					COMPARE_INTERVAL_LIST_VAR("SEQUENCE_TARGET", &sa->tar2);
					parse_interval_list_VAR("SEQUENCE_TARGET",datum,&sa->tar2,parse_err);
					sa->tar2.pairs[0][0]=sa->tar2.VAR_pairs[0][0];
					sa->tar2.pairs[0][1]=sa->tar2.VAR_pairs[0][1];
					sa->tar2.count=sa->tar2.char_count+1;
					continue;
				}else {
					COMPARE_INTERVAL_LIST("SEQUENCE_TARGET", &sa->tar2);
				}
			}

			COMPARE_INTERVAL_LIST("SEQUENCE_EXCLUDED_REGION", &sa->excl2);
			COMPARE_INTERVAL_LIST("SEQUENCE_INTERNAL_EXCLUDED_REGION", &sa->excl_internal2);

			if (COMPARE("SEQUENCE_OVERLAP_JUNCTION_LIST")) {
				if (parse_intron_list(datum, sa->primer_overlap_junctions, 
				                      &sa->primer_overlap_junctions_count) == 0) {
					pr_append_new_chunk(parse_err, "Error in SEQUENCE_PRIMER_OVERLAP_JUNCTION_LIST");
				}
				continue;
			}

			if (COMPARE("SEQUENCE_INCLUDED_REGION")) {
				p = parse_int_pair("SEQUENCE_INCLUDED_REGION", datum, ',',
				                   &sa->incl_s, &sa->incl_l, parse_err);
				if (NULL == p) /* An error; the message is already * in parse_err.  */
					continue;
				while (' ' == *p || '\t' == *p) p++;
				if (*p != '\n' && *p != '\0')
					tag_syntax_error("SEQUENCE_INCLUDED_REGION", datum, parse_err);
				continue;
			}
			COMPARE_INT("SEQUENCE_START_CODON_POSITION", sa->start_codon_pos);
			COMPARE_INT("SEQUENCE_FORCE_LEFT_START", sa->force_left_start);
			COMPARE_INT("SEQUENCE_FORCE_LEFT_END", sa->force_left_end);
			COMPARE_INT("SEQUENCE_FORCE_RIGHT_START", sa->force_right_start);
			COMPARE_INT("SEQUENCE_FORCE_RIGHT_END", sa->force_right_end);
			/* Process "Global" Arguments (those that persist between boulder
			 * records). */
			parse_err = glob_err;  /* These errors are considered fatal. */
			if (COMPARE("PRIMER_PRODUCT_SIZE_RANGE")) {
				parse_product_size("PRIMER_PRODUCT_SIZE_RANGE", datum, pa, parse_err);
				continue;
			}
			COMPARE_INT("PRIMER_OPT_SIZE", pa->p_args.opt_size);
			COMPARE_INT("PRIMER_MIN_SIZE", pa->p_args.min_size);
			COMPARE_INT("PRIMER_MAX_SIZE", pa->p_args.max_size);
			COMPARE_INT("PRIMER_MAX_POLY_X", pa->p_args.max_poly_x);
			COMPARE_FLOAT("PRIMER_OPT_TM", pa->p_args.opt_tm);
			COMPARE_FLOAT("PRIMER_OPT_GC_PERCENT", pa->p_args.opt_gc_content);
			COMPARE_FLOAT("PRIMER_MIN_TM", pa->p_args.min_tm);
			COMPARE_FLOAT("PRIMER_MAX_TM", pa->p_args.max_tm);
			COMPARE_FLOAT("PRIMER_PAIR_MAX_DIFF_TM", pa->max_diff_tm);
			if (COMPARE("PRIMER_TM_FORMULA")) {
				parse_int("PRIMER_TM_FORMULA", datum, &tmp_int, parse_err);
				pa->tm_santalucia = (tm_method_type) tmp_int;    /* added by T.Koressaar */
				continue;
			}
			if (COMPARE("PRIMER_SALT_CORRECTIONS")) {
				parse_int("PRIMER_SALT_CORRECTIONS", datum, &tmp_int, parse_err);
				pa->salt_corrections = (salt_correction_type) tmp_int; /* added by T.Koressaar */
				continue;
			}
			COMPARE_FLOAT("PRIMER_MIN_GC", pa->p_args.min_gc);
			COMPARE_FLOAT("PRIMER_MAX_GC", pa->p_args.max_gc);
			COMPARE_FLOAT("PRIMER_SALT_MONOVALENT", pa->p_args.salt_conc);
			COMPARE_FLOAT("PRIMER_SALT_DIVALENT", pa->p_args.divalent_conc);
			COMPARE_FLOAT("PRIMER_DNTP_CONC", pa->p_args.dntp_conc);
			COMPARE_FLOAT("PRIMER_DNA_CONC", pa->p_args.dna_conc);
			COMPARE_INT("PRIMER_MAX_NS_ACCEPTED", pa->p_args.num_ns_accepted);
			COMPARE_INT("PRIMER_PRODUCT_OPT_SIZE", pa->product_opt_size);
			COMPARE_FLOAT("PRIMER_MAX_SELF_ANY", pa->p_args.max_self_any);

			/* COMPARE_FLOAT("PRIMER_MAX_SELF_END", pa->p_args.max_self_end); */
			/* NEW */ COMPARE_FLOAT_USE_FN("PRIMER_MAX_SELF_END", p3_set_gs_primer_self_end)
			/* if (COMPARE("PRIMER_MAX_SELF_END")) {
			  parse_double("PRIMER_MAX_SELF_END", datum, &tmp_double, parse_err);
			  p3_set_gs_primer_self_end(pa, tmp_double);
			  continue;
			  } */

			COMPARE_FLOAT("PRIMER_MAX_SELF_ANY_TH", pa->p_args.max_self_any_th);
			COMPARE_FLOAT("PRIMER_MAX_SELF_END_TH", pa->p_args.max_self_end_th);
			COMPARE_FLOAT("PRIMER_MAX_HAIRPIN_TH", pa->p_args.max_hairpin_th);   
			COMPARE_FLOAT("PRIMER_PAIR_MAX_COMPL_ANY", pa->pair_compl_any);
			COMPARE_FLOAT("PRIMER_PAIR_MAX_COMPL_END", pa->pair_compl_end);
			COMPARE_FLOAT("PRIMER_PAIR_MAX_COMPL_ANY_TH", pa->pair_compl_any_th);
			COMPARE_FLOAT("PRIMER_PAIR_MAX_COMPL_END_TH", pa->pair_compl_end_th);
			COMPARE_INT("P3_FILE_FLAG", res->file_flag);
			COMPARE_INT("PRIMER_PICK_ANYWAY", pa->pick_anyway);
			COMPARE_INT("PRIMER_GC_CLAMP", pa->gc_clamp);
			COMPARE_INT("PRIMER_MAX_END_GC", pa->max_end_gc);
			COMPARE_INT("PRIMER_EXPLAIN_FLAG", res->explain_flag);
			COMPARE_INT("PRIMER_LIBERAL_BASE", pa->liberal_base);
			COMPARE_INT("PRIMER_FIRST_BASE_INDEX", pa->first_base_index);
			COMPARE_INT("PRIMER_NUM_RETURN", pa->num_return);
			COMPARE_INT("PRIMER_MIN_QUALITY", pa->p_args.min_quality);
			COMPARE_INT("PRIMER_MIN_END_QUALITY", pa->p_args.min_end_quality);
			if (COMPARE("PRIMER_MIN_THREE_PRIME_DISTANCE")) {
				parse_int("PRIMER_MIN_THREE_PRIME_DISTANCE", datum, &(min_three_prime_distance), parse_err);
				/* check if specific tag also specified - error in this case */
				if (min_3_prime_distance_specific == 1) {
					pr_append_new_chunk(glob_err,
					                    "Both PRIMER_MIN_THREE_PRIME_DISTANCE and PRIMER_{LEFT/RIGHT}_MIN_THREE_PRIME_DISTANCE specified");
				} else {
					min_3_prime_distance_global = 1;
					/* set up individual flags */
					pa->min_left_three_prime_distance = min_three_prime_distance;
					pa->min_right_three_prime_distance = min_three_prime_distance;
				}
				continue;
			}
			if (COMPARE("PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE")) {
				parse_int("PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE", datum, &(pa->min_left_three_prime_distance), parse_err);
				/* check if global tag also specified - error in this case */
				if (min_3_prime_distance_global == 1) {
					pr_append_new_chunk(glob_err,
					                    "Both PRIMER_MIN_THREE_PRIME_DISTANCE and PRIMER_{LEFT/RIGHT}_MIN_THREE_PRIME_DISTANCE specified");
				} else {
					min_3_prime_distance_specific = 1;
				}
				continue;
			}
			if (COMPARE("PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE")) {
				parse_int("PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE", datum, &(pa->min_right_three_prime_distance), parse_err);
				/* check if global tag also specified - error in this case */
				if (min_3_prime_distance_global == 1) {
				pr_append_new_chunk(glob_err,
				                    "Both PRIMER_MIN_THREE_PRIME_DISTANCE and PRIMER_{LEFT/RIGHT}_MIN_THREE_PRIME_DISTANCE specified");
				} else {
					min_3_prime_distance_specific = 1;
				}
				continue;
			}
			if (file_type == settings) {
			  if (COMPARE("P3_FILE_ID")) continue;
			}
			COMPARE_INT("PRIMER_QUALITY_RANGE_MIN", pa->quality_range_min);
			COMPARE_INT("PRIMER_QUALITY_RANGE_MAX", pa->quality_range_max);
			COMPARE_FLOAT("PRIMER_PRODUCT_MAX_TM", pa->product_max_tm);
			COMPARE_FLOAT("PRIMER_PRODUCT_MIN_TM", pa->product_min_tm);
			COMPARE_FLOAT("PRIMER_PRODUCT_OPT_TM", pa->product_opt_tm);
			COMPARE_INT("PRIMER_SEQUENCING_LEAD", pa->sequencing.lead);
			COMPARE_INT("PRIMER_SEQUENCING_SPACING", pa->sequencing.spacing);
			COMPARE_INT("PRIMER_SEQUENCING_INTERVAL", pa->sequencing.interval);
			COMPARE_INT("PRIMER_SEQUENCING_ACCURACY", pa->sequencing.accuracy);
			if (COMPARE("PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION")) {
				parse_int("PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION", datum, &pa->min_5_prime_overlap_of_junction, parse_err);
				/* min_5_prime = 1; Removed 10/20/2010 */
				continue;
			}
			if (COMPARE("PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION")) {
				parse_int("PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION", datum, &pa->min_3_prime_overlap_of_junction, parse_err);
				/* min_3_prime = 1; Removed 10/20/2010 */
				continue;
			}
			COMPARE_AND_MALLOC("PRIMER_TASK", task_tmp);
			COMPARE_INT("PRIMER_PICK_RIGHT_PRIMER", pa->pick_right_primer);
			COMPARE_INT("PRIMER_PICK_INTERNAL_OLIGO", pa->pick_internal_oligo);
			COMPARE_INT("PRIMER_PICK_LEFT_PRIMER", pa->pick_left_primer);

			/* parameters for using Variation information - 20121205 N.Kasahara */
			/*  modified Y.Kimura */
			COMPARE_INT("PRIMER_INTERNAL_GENOTYPING",sa->tar2.genotyping);
			COMPARE_INT("PRIMER_INTERNAL_MOD_NEAR_VAR",pa->mod_near_var);
			/*  Threshold for distance between modified nucleotide and variant - Y.Kimura */
			COMPARE_INT("PRIMER_INTERNAL_MOD_VAR_DISTANCE",pa->thre_distance_mod_var);
			/*  Terminal region of internal probe excluded for target variant  - Y.Kimura */
			COMPARE_INT("PRIMER_INTERNAL_TERMINAL_EXCLUDED_VAR",pa->excl_terminal_intl_var);
			/*  Trminal region of internal probe excluded for modification - Y.Kimura */
			COMPARE_INT("PRIMER_INTERNAL_5_PRIME_TERMINAL_EXCLUDED_MOD",pa->excl_5_prime_intl_mod);
			COMPARE_INT("PRIMER_INTERNAL_3_PRIME_TERMINAL_EXCLUDED_MOD",pa->excl_3_prime_intl_mod);
			/*  Terminal region of primer excluded for modification - Y.Kimura */
			COMPARE_INT("PRIMER_5_PRIME_TERMINAL_EXCLUDED_MOD",pa->excl_5_prime_primer_mod);
			COMPARE_INT("PRIMER_3_PRIME_TERMINAL_EXCLUDED_MOD",pa->excl_3_prime_primer_mod);

			/* parameters for treating Z - 20120903 N.Kasahara */
			/*  modified Y.Kimura */
			COMPARE_INT("PRIMER_MODIFY_RIGHT_PRIMER",pa->modify_right_primer);
			COMPARE_INT("PRIMER_MODIFY_LEFT_PRIMER",pa->modify_left_primer);
			COMPARE_INT("PRIMER_INTERNAL_MODIFY_OLIGO",pa->modify_internal_oligo);

			COMPARE_INT("PRIMER_INTERNAL_OPT_SIZE", pa->o_args.opt_size);
			COMPARE_INT("PRIMER_INTERNAL_MAX_SIZE", pa->o_args.max_size);
			COMPARE_INT("PRIMER_INTERNAL_MIN_SIZE", pa->o_args.min_size);
			COMPARE_INT("PRIMER_INTERNAL_MAX_POLY_X", pa->o_args.max_poly_x);
			COMPARE_FLOAT("PRIMER_INTERNAL_OPT_TM", pa->o_args.opt_tm);
			COMPARE_FLOAT("PRIMER_INTERNAL_OPT_GC_PERCENT", pa->o_args.opt_gc_content);
			COMPARE_FLOAT("PRIMER_INTERNAL_MAX_TM", pa->o_args.max_tm);
			COMPARE_FLOAT("PRIMER_INTERNAL_MIN_TM", pa->o_args.min_tm);
			/* add PRIMER_INTERNAL_MIN_TM_VAR     20141103 - Y.Kimura*/
			COMPARE_FLOAT("PRIMER_INTERNAL_MIN_TM_VAR", pa->o_args.min_tm_var);
			COMPARE_FLOAT("PRIMER_INTERNAL_MIN_GC", pa->o_args.min_gc);
			COMPARE_FLOAT("PRIMER_INTERNAL_MAX_GC", pa->o_args.max_gc);
			COMPARE_FLOAT("PRIMER_INTERNAL_SALT_MONOVALENT", pa->o_args.salt_conc);
			COMPARE_FLOAT("PRIMER_INTERNAL_SALT_DIVALENT", pa->o_args.divalent_conc);
			COMPARE_FLOAT("PRIMER_INTERNAL_DNTP_CONC", pa->o_args.dntp_conc);
			COMPARE_FLOAT("PRIMER_INTERNAL_DNA_CONC", pa->o_args.dna_conc);
			COMPARE_INT("PRIMER_INTERNAL_MAX_NS_ACCEPTED", pa->o_args.num_ns_accepted);
			COMPARE_INT("PRIMER_INTERNAL_MIN_QUALITY", pa->o_args.min_quality);
			COMPARE_FLOAT("PRIMER_INTERNAL_MAX_SELF_ANY", pa->o_args.max_self_any);

			/* COMPARE_FLOAT("PRIMER_INTERNAL_MAX_SELF_END", pa->o_args.max_self_end); */
			/* NEW  */ COMPARE_FLOAT_USE_FN("PRIMER_INTERNAL_MAX_SELF_END", p3_set_gs_primer_internal_oligo_self_end);
			/* if (COMPARE("PRIMER_INTERNAL_MAX_SELF_END")) { 
			 parse_double("PRIMER_INTERNAL_MAX_SELF_END", datum, &tmp_double, parse_err);
			 p3_set_gs_primer_internal_oligo_self_end(pa, tmp_double);
			 continue;
			 } */

			COMPARE_FLOAT("PRIMER_INTERNAL_MAX_SELF_ANY_TH", pa->o_args.max_self_any_th);
			COMPARE_FLOAT("PRIMER_INTERNAL_MAX_SELF_END_TH", pa->o_args.max_self_end_th);
			COMPARE_FLOAT("PRIMER_INTERNAL_MAX_HAIRPIN_TH", pa->o_args.max_hairpin_th);
			COMPARE_FLOAT("PRIMER_MAX_LIBRARY_MISPRIMING", pa->p_args.max_repeat_compl);
			COMPARE_FLOAT("PRIMER_INTERNAL_MAX_LIBRARY_MISHYB", pa->o_args.max_repeat_compl);

			/* parameters for internal_oligo template_mishybridization - N.Kasahara, Y.Kimura */
			COMPARE_FLOAT("PRIMER_INTERNAL_MAX_TEMPLATE_MISHYB",pa->o_args.max_template_mishyb);
			COMPARE_FLOAT("PRIMER_INTERNAL_MAX_TEMPLATE_MISHYB_TH",pa->o_args.max_template_mishyb_th);

			COMPARE_FLOAT("PRIMER_PAIR_MAX_LIBRARY_MISPRIMING", pa->pair_repeat_compl);
			/* Mispriming / mishybing in the template. */
			COMPARE_FLOAT("PRIMER_MAX_TEMPLATE_MISPRIMING", pa->p_args.max_template_mispriming);
			COMPARE_FLOAT("PRIMER_MAX_TEMPLATE_MISPRIMING_TH", pa->p_args.max_template_mispriming_th);
			COMPARE_FLOAT("PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING", pa->pair_max_template_mispriming);
			COMPARE_FLOAT("PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH", pa->pair_max_template_mispriming_th);

			 /* Control interpretation of ambiguity codes in mispriming
			    and mishyb libraries. */
			COMPARE_INT("PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS", pa->lib_ambiguity_codes_consensus);
			COMPARE_FLOAT("PRIMER_INSIDE_PENALTY", pa->inside_penalty);
			COMPARE_FLOAT("PRIMER_OUTSIDE_PENALTY", pa->outside_penalty);
			if (COMPARE("PRIMER_MISPRIMING_LIBRARY")) {
				if (repeat_file_path != NULL) {
					pr_append_new_chunk(glob_err,
					                    "Duplicate PRIMER_MISPRIMING_LIBRARY tag");
					free(repeat_file_path);
					repeat_file_path = NULL;
				} else {
					repeat_file_path = (char*) _rb_safe_malloc(strlen(datum) + 1);
					strcpy(repeat_file_path, datum);
				}
				continue;
			}
			if (COMPARE("PRIMER_INTERNAL_MISHYB_LIBRARY")) {
				if (int_repeat_file_path != NULL) {
					pr_append_new_chunk(glob_err,
					                    "Duplicate PRIMER_INTERNAL_MISHYB_LIBRARY tag");
					free(int_repeat_file_path);
					int_repeat_file_path = NULL;
				} else {
					int_repeat_file_path = (char*) _rb_safe_malloc(strlen(datum) + 1);
					strcpy(int_repeat_file_path, datum);
				}
				continue;
			}
			if (COMPARE("P3_COMMENT")) continue;
			COMPARE_FLOAT("PRIMER_MAX_END_STABILITY", pa->max_end_stability);

			COMPARE_INT("PRIMER_LOWERCASE_MASKING",
			            pa->lowercase_masking); 
			/* added by T. Koressaar */
			COMPARE_INT("PRIMER_THERMODYNAMIC_ALIGNMENT", pa->thermodynamic_alignment);
			if (COMPARE("PRIMER_THERMODYNAMIC_PARAMETERS_PATH")) {
				if (thermodynamic_params_path == NULL) {
					thermodynamic_params_path = (char*) _rb_safe_malloc(datum_len + 1);
					strcpy(thermodynamic_params_path, datum);
					thermodynamic_path_changed = 1;
				}
				/* check if path changes */
				else if (strcmp(thermodynamic_params_path, datum)) {
					free(thermodynamic_params_path);
					thermodynamic_params_path = (char*) _rb_safe_malloc(datum_len + 1); 
					strcpy(thermodynamic_params_path, datum);
					thermodynamic_path_changed = 1;
				}
				continue;
			}
			/* weights for objective functions  */
			/* CHANGE TEMP/temp -> TM/tm */
			COMPARE_FLOAT("PRIMER_WT_TM_GT", pa->p_args.weights.temp_gt);
			COMPARE_FLOAT("PRIMER_WT_TM_LT", pa->p_args.weights.temp_lt);
			COMPARE_FLOAT("PRIMER_WT_GC_PERCENT_GT", pa->p_args.weights.gc_content_gt);
			COMPARE_FLOAT("PRIMER_WT_GC_PERCENT_LT", pa->p_args.weights.gc_content_lt);
			COMPARE_FLOAT("PRIMER_WT_SIZE_LT", pa->p_args.weights.length_lt);
			COMPARE_FLOAT("PRIMER_WT_SIZE_GT", pa->p_args.weights.length_gt);
			COMPARE_FLOAT("PRIMER_WT_SELF_ANY", pa->p_args.weights.compl_any);
			COMPARE_FLOAT("PRIMER_WT_SELF_END", pa->p_args.weights.compl_end);
			COMPARE_FLOAT("PRIMER_WT_SELF_ANY_TH", pa->p_args.weights.compl_any_th);
			COMPARE_FLOAT("PRIMER_WT_SELF_END_TH", pa->p_args.weights.compl_end_th);
			COMPARE_FLOAT("PRIMER_WT_HAIRPIN_TH", pa->p_args.weights.hairpin_th);
			COMPARE_FLOAT("PRIMER_WT_NUM_NS", pa->p_args.weights.num_ns);
			COMPARE_FLOAT("PRIMER_WT_LIBRARY_MISPRIMING", pa->p_args.weights.repeat_sim);
			COMPARE_FLOAT("PRIMER_WT_SEQ_QUAL", pa->p_args.weights.seq_quality);
			COMPARE_FLOAT("PRIMER_WT_END_QUAL", pa->p_args.weights.end_quality);
			COMPARE_FLOAT("PRIMER_WT_POS_PENALTY", pa->p_args.weights.pos_penalty);
			COMPARE_FLOAT("PRIMER_WT_END_STABILITY", pa->p_args.weights.end_stability);
			COMPARE_FLOAT("PRIMER_WT_TEMPLATE_MISPRIMING", pa->p_args.weights.template_mispriming);
			COMPARE_FLOAT("PRIMER_WT_TEMPLATE_MISPRIMING_TH", pa->p_args.weights.template_mispriming_th);
			COMPARE_FLOAT("PRIMER_INTERNAL_WT_TM_GT", pa->o_args.weights.temp_gt);
			COMPARE_FLOAT("PRIMER_INTERNAL_WT_TM_LT", pa->o_args.weights.temp_lt);
			COMPARE_FLOAT("PRIMER_INTERNAL_WT_GC_PERCENT_GT", pa->o_args.weights.gc_content_gt);
			COMPARE_FLOAT("PRIMER_INTERNAL_WT_GC_PERCENT_LT", pa->o_args.weights.gc_content_lt);
			COMPARE_FLOAT("PRIMER_INTERNAL_WT_SIZE_LT", pa->o_args.weights.length_lt);
			COMPARE_FLOAT("PRIMER_INTERNAL_WT_SIZE_GT", pa->o_args.weights.length_gt);
			COMPARE_FLOAT("PRIMER_INTERNAL_WT_SELF_ANY", pa->o_args.weights.compl_any);
			COMPARE_FLOAT("PRIMER_INTERNAL_WT_SELF_END", pa->o_args.weights.compl_end);
			COMPARE_FLOAT("PRIMER_INTERNAL_WT_SELF_ANY_TH", pa->o_args.weights.compl_any_th);
			COMPARE_FLOAT("PRIMER_INTERNAL_WT_SELF_END_TH", pa->o_args.weights.compl_end_th);
			COMPARE_FLOAT("PRIMER_INTERNAL_WT_HAIRPIN_TH", pa->o_args.weights.hairpin_th);
			COMPARE_FLOAT("PRIMER_INTERNAL_WT_NUM_NS", pa->o_args.weights.num_ns);
			COMPARE_FLOAT("PRIMER_INTERNAL_WT_LIBRARY_MISHYB", pa->o_args.weights.repeat_sim);
			COMPARE_FLOAT("PRIMER_INTERNAL_WT_SEQ_QUAL", pa->o_args.weights.seq_quality);
			COMPARE_FLOAT("PRIMER_INTERNAL_WT_END_QUAL", pa->o_args.weights.end_quality);

			/* PRIMER_INTERNAL_WT_DIFF_TM_VAR     2014/11/03 - Y.Kimura */
			COMPARE_FLOAT("PRIMER_INTERNAL_WT_DIFF_TM_VAR", pa->o_args.weights.diff_temp_var);
			/* PRIMER_INTERNAL_WT_MOD_POS     2014/11/03 - Y.Kimura */
			COMPARE_FLOAT("PRIMER_INTERNAL_WT_MOD_POS", pa->o_args.weights.mod_pos);

			/* parameters for internal probe template_mishybridazation  - Y.Kimura */
			COMPARE_FLOAT("PRIMER_INTERNAL_WT_TEMPLATE_MISHYB", pa->o_args.weights.template_mishyb);
			COMPARE_FLOAT("PRIMER_INTERNAL_WT_TEMPLATE_MISHYB_TH", pa->o_args.weights.template_mishyb_th);

			COMPARE_FLOAT("PRIMER_PAIR_WT_PR_PENALTY", pa->pr_pair_weights.primer_quality);
			COMPARE_FLOAT("PRIMER_PAIR_WT_IO_PENALTY", pa->pr_pair_weights.io_quality);
			COMPARE_FLOAT("PRIMER_PAIR_WT_DIFF_TM", pa->pr_pair_weights.diff_tm);
			COMPARE_FLOAT("PRIMER_PAIR_WT_COMPL_ANY", pa->pr_pair_weights.compl_any);
			COMPARE_FLOAT("PRIMER_PAIR_WT_COMPL_END", pa->pr_pair_weights.compl_end);
			COMPARE_FLOAT("PRIMER_PAIR_WT_COMPL_ANY_TH", pa->pr_pair_weights.compl_any_th);
			COMPARE_FLOAT("PRIMER_PAIR_WT_COMPL_END_TH", pa->pr_pair_weights.compl_end_th);
			COMPARE_FLOAT("PRIMER_PAIR_WT_PRODUCT_TM_LT", pa->pr_pair_weights.product_tm_lt);
			COMPARE_FLOAT("PRIMER_PAIR_WT_PRODUCT_TM_GT", pa->pr_pair_weights.product_tm_gt);
			COMPARE_FLOAT("PRIMER_PAIR_WT_PRODUCT_SIZE_GT", pa->pr_pair_weights.product_size_gt);
			COMPARE_FLOAT("PRIMER_PAIR_WT_PRODUCT_SIZE_LT", pa->pr_pair_weights.product_size_lt);
			COMPARE_FLOAT("PRIMER_PAIR_WT_LIBRARY_MISPRIMING", pa->pr_pair_weights.repeat_sim);
			COMPARE_FLOAT("PRIMER_PAIR_WT_TEMPLATE_MISPRIMING", pa->pr_pair_weights.template_mispriming);
			COMPARE_FLOAT("PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH", pa->pr_pair_weights.template_mispriming_th);

			/* flag of internal_oligo_direction - 20121009 N.Kasahara */
			if(COMPARE("PRIMER_INTERNAL_OLIGO_DIRECTION")){
				pa->internal_oligo_direction=atoi(datum);
				if(pa->internal_oligo_direction==2 && sa->internal_input!=NULL){
					pr_append_new_chunk(glob_err,
					"in case of inputting internal probe sequence, can't assign internal probe direction to be ANY");
					return(1);
				}
				continue;
			}
		}
		/* End of reading the tags in the right place */

		/*  Complain about unrecognized tags */
		if (*strict_tags == 1) {
			pr_append_new_chunk(glob_err, "Unrecognized tag: ");
			pr_append(glob_err, s);
			fprintf(stderr, "Unrecognized tag: %s\n", s);
		}
	}  /* while ((s = p3_read_line(stdin)) != NULL && strcmp(s,"=")) { */

	/* check input left primer sequence - 20120911 N.Kasahara */
	if(sa->left_input !=NULL){
		length=strlen(sa->left_input);
	
		/* memory allocate for keeping input sequence - 20130305 N.Kasahara */
		if((sa->left_input_original =(char *)malloc(sizeof(char)*(length+1)))==NULL){
			(void)printf("memory allocation error of sa->left_input_original\n");
			return(0);
		}
		strcpy(sa->left_input_original,sa->left_input);

		for(i=0,zcount=0;i<length;i++){
			if(*(sa->left_input+i)=='Z' ||*(sa->left_input+i)=='z'){
				pa->zpos_left_primer=i;
				zcount++;
				pa->znum_left_primer=zcount;
			}
			if(isalpha(*(sa->left_input+i))==0) {
				*(sa->left_input+i)=0x00;
				break;
			}
		}
		for(i=0,zcount=0;i<length;i++){
			//if((*(sa->left_input+i)==MODIFIABLE_BASE_UP || *(sa->left_input+i)==MODIFIABLE_BASE_LW) && pa->znum_left_primer==0){
			if(IS_MODIFIABLE_BASE(*(sa->left_input+i)) && pa->znum_left_primer==0){
				zcount++;
			}
		}
		if(pa->znum_left_primer==1){
			if(*(sa->left_input+pa->zpos_left_primer)=='Z'){
				*(sa->left_input+pa->zpos_left_primer)=MODIFIABLE_BASE_UP;
			}else if(*(sa->left_input+pa->zpos_left_primer)=='z'){
				*(sa->left_input+pa->zpos_left_primer)=MODIFIABLE_BASE_LW;
			}
		}
		if(pa->znum_left_primer>=2){
			pr_append_new_chunk(glob_err,"Too Many Z in LEFT PRIMER sequence\n");
		}
		if(pa->znum_left_primer!=0 && pa->modify_left_primer==0 ) {
			pr_append_new_chunk(glob_err,"Z in LEFT PRIMER sequence while modification option is not set\n");
		}
		/* WARNING message in case of no available base for modification - 20130222 N.Kasahara */
		if(zcount==0 && pa->modify_left_primer==1 && pa->znum_left_primer==0){
			pr_append_new_chunk(glob_err,"No available base for modification in LEFT PRIMER sequence\n");
		}	

		/* position check of Z - 20130301 N.Kasahara */
		/* 20141102 Y.Kimura modified */
		if(pa->znum_left_primer==1){
			if(pa->zpos_left_primer < pa->excl_5_prime_primer_mod){
				pr_append_new_chunk(glob_err,"Modified nucleotide (Z) is set in the excluded region of LEFT PRIMER 5'end\n");
			}
			if(pa->zpos_left_primer >= (length - pa->excl_3_prime_primer_mod)){
				pr_append_new_chunk(glob_err,"Modified nucleotide (Z) is set in the excluded region of LEFT PRIMER 3'end\n");
			}
		}
	}

	/* check input right primer sequence - 20120911 N.Kasahara */
	if(sa->right_input!=NULL){
		length=strlen(sa->right_input);
		
		/*	memory allocate for keeping input sequence	*/
		/*		20130305 - 	N.Kasahara		*/
		if((sa->right_input_original =(char *)malloc(sizeof(char)*(length+1)))==NULL){
			(void)printf("memory allocation error of sa->right_input_original\n");
			return(0);
		}
		strcpy(sa->right_input_original,sa->right_input);

		for(i=0,zcount=0;i<length;i++){
			if(*(sa->right_input+i)=='Z'||*(sa->right_input+i)=='z'){
				pa->zpos_right_primer=i;
				zcount++;
				pa->znum_right_primer=zcount;
			}
			if(isalpha(*(sa->right_input+i))==0) {
				*(sa->right_input+i)=0x00;
				break;
			}
		}
		for(i=0,zcount=0;i<length;i++){
			//if((*(sa->right_input+i)==MODIFIABLE_BASE_UP || *(sa->right_input+i)==MODIFIABLE_BASE_LW) && pa->znum_right_primer==0){
			if(IS_MODIFIABLE_BASE(*(sa->right_input+i)) && pa->znum_right_primer==0){
				zcount++;
			}
		}
		if(pa->znum_right_primer==1){
			if(*(sa->right_input+pa->zpos_right_primer)=='Z'){
				*(sa->right_input+pa->zpos_right_primer)=MODIFIABLE_BASE_UP;
			}else if(*(sa->right_input+pa->zpos_right_primer)=='z'){
				*(sa->right_input+pa->zpos_right_primer)=MODIFIABLE_BASE_LW;
			}
		}
		if(pa->znum_right_primer>=2){
			pr_append_new_chunk(glob_err,"Too Many Z in RIGHT PRIMER sequence\n");
		}
		if(pa->znum_right_primer!=0 && pa->modify_right_primer==0){
			pr_append_new_chunk(glob_err,"Z in RIGHT PRIMER sequence while modification option it is not set\n");
		}
		/* WARNING message in case of no available base for modification - 20130222 N.Kasahara */
		if(zcount==0 && pa->modify_right_primer==1 && pa->znum_right_primer==0){
			pr_append_new_chunk(glob_err,"No available base for modification in RIGHT PRIMER sequence\n");
		}	

		/* position check of Z - 20130301 N.Kasahara */
		/* 20141102 Y.Kimura modified */
		if(pa->znum_right_primer==1){
			if(pa->zpos_right_primer < pa->excl_5_prime_primer_mod){
				pr_append_new_chunk(glob_err,"Modified nucleotide (Z) is set in the excluded region of RIGHT PRIMER 5'end\n");
			}
			if(pa->zpos_right_primer >= (length - pa->excl_3_prime_primer_mod)){
				pr_append_new_chunk(glob_err,"Modified nucleotide (Z) is set in the excluded region of RIGHT PRIMER 3'end\n");
			}
		}
	}

	/* check input internal oligo sequence - 20120911 N.Kasahara */
	if(sa->internal_input!=NULL && pa->internal_oligo_direction==0){
		length=strlen(sa->internal_input);

		/*	memory allocate for keeping input sequence	*/
		/*		20130305 - 	N.Kasahara		*/
		if((sa->internal_input_original =(char *)malloc(sizeof(char)*(length+1)))==NULL){
			(void)printf("memory allocation error of sa->internal_input_original\n");
			return(0);
		}
		strcpy(sa->internal_input_original,sa->internal_input);

		for(i=0,zcount=0;i<length;i++){
			if(*(sa->internal_input+i)=='Z' || *(sa->internal_input+i)=='z'){
				pa->zpos_internal_oligo=i;
				zcount++;
				pa->znum_internal_oligo=zcount;
			}
			if(isalpha(*(sa->internal_input+i))==0) {
				*(sa->internal_input+i)=0x00;
				break;
			}
		}
		for(i=0,zcount=0;i<length;i++){
			//if((*(sa->internal_input+i)==MODIFIABLE_BASE_UP||*(sa->internal_input+i)==MODIFIABLE_BASE_LW) && pa->znum_internal_oligo==0){
			if(IS_MODIFIABLE_BASE(*(sa->internal_input+i)) && pa->znum_internal_oligo==0){
				zcount++;
			}
		}
		if(pa->znum_internal_oligo == 1){
			if(*(sa->internal_input + pa->zpos_internal_oligo) == 'Z'){
				*(sa->internal_input + pa->zpos_internal_oligo) = MODIFIABLE_BASE_UP;
			}else if(*(sa->internal_input + pa->zpos_internal_oligo) == 'z'){
				*(sa->internal_input + pa->zpos_internal_oligo) = MODIFIABLE_BASE_LW;
			}
		}
		if(pa->znum_internal_oligo>=2){
			pr_append_new_chunk(glob_err,"Too Many Z in INTERNAL PROBE sequence\n");
		}
		if(pa->znum_internal_oligo!=0 && pa->modify_internal_oligo==0){
			pr_append_new_chunk(glob_err,"Z in INTERNAL PROBE sequence while modification option is not set\n");
		}
		/* WARNING message in case of no available base for modification - 20130222 N.Kasahara */
		if(zcount==0 && pa->modify_internal_oligo==1 && pa->znum_internal_oligo==0){
			pr_append_new_chunk(glob_err,"No available base for modification in INTERNAL PROBE sequence\n");
		}

		/* position check of Z - 20130301 N.Kasahara */
		/* 20141102 Y.Kimura modified */
		if(pa->znum_internal_oligo == 1){
			if(pa->zpos_internal_oligo < pa->excl_5_prime_intl_mod){
				pr_append_new_chunk(glob_err,"Modified nucleotide (Z) is set in the excluded region of INTERNAL PROBE 5'end\n");
			}
			if(pa->zpos_internal_oligo >= (length - pa->excl_3_prime_intl_mod)){
				pr_append_new_chunk(glob_err,"Modified nucleotide (Z) is set in the excluded region of INTERNAL PROBE 3'end\n");
			}
		}

	}else if(sa->internal_input!=NULL && pa->internal_oligo_direction==1){
		length=strlen(sa->internal_input);

		/*	memory allocate for keeping input sequence	*/
		/*		20130305 - 	N.Kasahara		*/
		if((sa->internal_input_original =(char *)malloc(sizeof(char)*(length+1)))==NULL){
			(void)printf("memory allocation error of sa->internal_input_original\n");
			return(0);
		}
		strcpy(sa->internal_input_original,sa->internal_input);

		for(i=0,zcount=0;i<length;i++){
			if(*(sa->internal_input+i)=='Z' || *(sa->internal_input+i)=='z'){
				pa->zpos_internal_oligo_REV=i;
				zcount++;
				pa->znum_internal_oligo_REV=zcount;
			}
			if(isalpha(*(sa->internal_input+i))==0){
				*(sa->internal_input+i)=0x00;
				break;
			}
		}
		for(i=0,zcount=0;i<length;i++){
			//if((*(sa->internal_input+i)==MODIFIABLE_BASE_UP||*(sa->internal_input+i)==MODIFIABLE_BASE_LW) && pa->znum_internal_oligo_REV==0){
			if(IS_MODIFIABLE_BASE(*(sa->internal_input+i)) && pa->znum_internal_oligo_REV==0){
				zcount++;
			}
		}
		if(pa->znum_internal_oligo_REV==1){
			if(*(sa->internal_input + pa->zpos_internal_oligo_REV) == 'Z'){
				*(sa->internal_input + pa->zpos_internal_oligo_REV) = MODIFIABLE_BASE_UP;
			}else if(*(sa->internal_input + pa->zpos_internal_oligo_REV) == 'z'){
				*(sa->internal_input + pa->zpos_internal_oligo_REV) = MODIFIABLE_BASE_LW;
			}
		}
		if(pa->znum_internal_oligo_REV>=2){
			pr_append_new_chunk(glob_err,"Too Many Z in INTERNAL PROBE sequence\n");
		}
		if(pa->znum_internal_oligo_REV!=0 && pa->modify_internal_oligo==0){
			pr_append_new_chunk(glob_err,"Z in INTERNAL PROBE sequence while modification option is not set\n");
		}
		/* WARNING message in case of no available base for modification - 20130222 N.Kasahara */
		if(zcount==0 && pa->modify_internal_oligo==1 && pa->znum_internal_oligo_REV==0){
			pr_append_new_chunk(glob_err,"No available base for modification in INTERNAL PROBE sequence\n");
		}

		/* position check of Z - 20130301 N.Kasahara */
		/* 20141102 Y.Kimura modified */
		if(pa->znum_internal_oligo_REV==1){
			if(pa->zpos_internal_oligo_REV < pa->excl_5_prime_intl_mod){
				pr_append_new_chunk(glob_err,"Modified nucleotide (Z) is set in the excluded region of INTERNAL PROBE 5'end\n");
			}
			if(pa->zpos_internal_oligo_REV >= (length - pa->excl_3_prime_intl_mod)){
				pr_append_new_chunk(glob_err,"Modified nucleotide (Z) is set in the excluded region of INTERNAL PROBE 3'end\n");
			}
		}

	}else if(sa->internal_input!=NULL && pa->internal_oligo_direction==2){
		pr_append_new_chunk(glob_err,"To use own internal probe sequence, Internal Probe Direction cannot be set as ANY\n");
	}

	/* check variation position and GENOTYPING flag - 20121206 N.Kasahara */
	if(sa->tar2.genotyping==1 && sa->tar2.count==0){
		pr_append_new_chunk(glob_err,"Target Variation position is NOT set.\n");
	}

	/* change for making sequences for each variation pattern - 20121205 - N.Kasahara */
	/* 20140926 modified by A.Sekiguchi to count '/' even if genotyping==0 */
	if(sa->sequence!=NULL){
		seq_len=strlen(sa->sequence);
		for(i=0,zcount=0,sa->tar2.char_count=0;i<seq_len;i++){
			if(*(sa->sequence+i)=='/'){
				sa->tar2.char_count++;
			}
		}	/*	count the number of '/' in template sequence (sa->sequence)	*/
	}

	if(sa->sequence!=NULL && sa->tar2.genotyping!=0 && sa->tar2.count!=0){
		/*	memory allocation for sequences of each Variation pattern	*/
		sa->sequence_all_var=(char *)malloc(sizeof(char)*((sa->tar2.char_count+1)*seq_len+1));
		sa->tar2.char_position=(int *)malloc(sizeof(int)*((sa->tar2.char_count+1)));
		sa->tar2.char_num=(int *)malloc(sizeof(int)*((sa->tar2.char_count+1)));
		sa->tar2.var_len=(int *)malloc(sizeof(int)*((sa->tar2.char_count+1)));

		*(sa->sequence_all_var)=NULL;
		
		for(i=0,j=0,k=0;j<=sa->tar2.char_count && i<seq_len;i++){
			if(*(sa->sequence+i)=='/'){
				sa->tar2.char_position[j]=i;
				j++;
			}else if(i==(sa->tar2.VAR_pairs[k][0]+sa->tar2.VAR_pairs[k][1])-1){
				sa->tar2.char_position[j]=i;
				j++;k++;
			}
		}
		for(i=0;i<=sa->tar2.count;i++){
			if(i==0) sa->tar2.VAR_start[i]=sa->tar2.VAR_pairs[i][0];
			else     sa->tar2.VAR_start[i]=(sa->tar2.char_position[i-1])+2;
		}

		/* loop for each [ ] */
		for(k=0;k<sa->tar2.count;k++){
			start=sa->tar2.VAR_pairs[k][0];
			length=sa->tar2.VAR_pairs[k][1];

			for(i=0;i<=sa->tar2.char_count;i++){
				if(i==0){
					sa->tar2.char_num[i]=sa->tar2.char_position[i]-start+1;
				} else if(i!=0 && i<sa->tar2.char_count){
					sa->tar2.char_num[i]=sa->tar2.char_position[i]-sa->tar2.char_position[i-1]-1;
				} else if(i!=0 && i==sa->tar2.char_count){
					sa->tar2.char_num[i]=start+length-sa->tar2.char_position[i-1]-1;
				}
			}
			*(sa->sequence_all_var)=NULL;
			for(m=0,j=0,h=0;m<=(sa->tar2.char_count);m++){
				if(m>0) h=sa->tar2.var_len[m-1];
				for(i=0;i<start-1;i++,j++){
					*(sa->sequence_all_var+j)=*(sa->sequence+i);
				}
				if(m==0){
					for(;i<sa->tar2.char_position[m];i++,j++){
						*(sa->sequence_all_var+j)=*(sa->sequence+i);
					}
				} else if(m!=0){
					for(i=sa->tar2.char_position[m-1]+1;i<sa->tar2.char_position[m];i++,j++){
						*(sa->sequence_all_var+j)=*(sa->sequence+i);
					}
				}

				for(i=start+length-1;i<seq_len;i++,j++){
					*(sa->sequence_all_var+j)=*(sa->sequence+i);
				}
				
				sa->tar2.var_len[m]=j;
			}
		}
	}

	/* 20121206 - N.Kasahara	*/
	for(k=1;k<=sa->tar2.count;k++){
		sa->tar2.VAR_pairs[k][0]=sa->tar2.VAR_pairs[0][0];
		sa->tar2.VAR_pairs[k][1]=sa->tar2.VAR_pairs[0][1];
	}

	/* set flag of pick_primers for special case (temporary processing) - 20121210 N.Kasahara */
	if((sa->tar2.genotyping==1 && pa->pick_left_primer==1 && pa->pick_right_primer==0) || 
		(sa->tar2.genotyping==1 && pa->pick_left_primer==0 && pa->pick_right_primer==1)){
		/* WARNING message - 20121211 N.Kasahara */
		(void)printf("WARNING : if pick only one primer (left or right), change picking both primers\n");
		pa->pick_left_primer=1;
		pa->pick_right_primer=1;
	}

	/* Check if the record was terminated by "=" */
	if (NULL == s) { /* End of file. */
		if (data_found) {
			 pr_append_new_chunk(glob_err, 
			                     "Final record not terminated by '='");
			return 1;
		} else  return 0;
	}

	/* Figure out the right settings for the tasks*/
	if (task_tmp != NULL) {

		if (!strcmp_nocase(task_tmp, "pick_pcr_primers")) {
			pa->primer_task = generic;
			pa->pick_left_primer = 1;
			pa->pick_right_primer = 1;
			pa->pick_internal_oligo = 0;
		} else if (!strcmp_nocase(task_tmp, "pick_pcr_primers_and_hyb_probe")) {
			pa->primer_task = generic; 
			pa->pick_left_primer = 1;
			pa->pick_right_primer = 1;
			pa->pick_internal_oligo = 1;
		} else if (!strcmp_nocase(task_tmp, "pick_left_only")) {
			pa->primer_task = generic;
			pa->pick_left_primer = 1;
			pa->pick_right_primer = 0;
			pa->pick_internal_oligo = 0;
		} else if (!strcmp_nocase(task_tmp, "pick_right_only")) {
			pa->primer_task = generic;
			pa->pick_left_primer = 0;
			pa->pick_right_primer = 1;
			pa->pick_internal_oligo = 0;
		} else if (!strcmp_nocase(task_tmp, "pick_hyb_probe_only")) {
			pa->primer_task = generic;
			pa->pick_left_primer = 0;
			pa->pick_right_primer = 0;
			pa->pick_internal_oligo = 1;
		} else if (*io_version == 3) {
			pr_append_new_chunk(glob_err, "Unrecognized PRIMER_TASK");
		} else if (!strcmp_nocase(task_tmp, "generic")) {
			pa->primer_task = generic;
		} else if (!strcmp_nocase(task_tmp, "pick_detection_primers")) {
			pa->primer_task = generic; /* Deliberate duplication for backward compatibility. */
		} else if (!strcmp_nocase(task_tmp, "pick_cloning_primers")) {
			pa->primer_task = pick_cloning_primers;
		} else if (!strcmp_nocase(task_tmp, "pick_discriminative_primers")) {
			pa->primer_task = pick_discriminative_primers;
		} else if (!strcmp_nocase(task_tmp, "pick_sequencing_primers")) {
			pa->primer_task = pick_sequencing_primers;
		} else if (!strcmp_nocase(task_tmp, "pick_primer_list")) {
			pa->primer_task = pick_primer_list;
		} else if (!strcmp_nocase(task_tmp, "check_primers")) {
			pa->primer_task = check_primers;
			/* check_primers sets the picking flags itself */
			pa->pick_left_primer = 0;
			pa->pick_right_primer = 0;
			pa->pick_internal_oligo = 0;
			if (sa->left_input){
				pa->pick_left_primer = 1;
			}
			if (sa->right_input){
				pa->pick_right_primer = 1;
			}
			if (sa->internal_input){
				pa->pick_internal_oligo = 1;
			}

			/* checking Z in sa->****input - 201200910 N.Kasahara */
			if(sa->left_input){
				length=strlen(sa->left_input);
				for(i=0;i<length;i++){
					if(*(sa->left_input+i)=='Z' || *(sa->left_input+i)=='z'){
						zcount++;
					}
				}
				if(zcount==1)	pa->modify_left_primer=1;
				else if (zcount>1) pr_append_new_chunk(glob_err,"Too Many Z in sa->left_input\n");
			}
			if(sa->right_input){
				length=strlen(sa->right_input);
				for(i=0;i<length;i++){
					if(*(sa->right_input+i)=='Z' || *(sa->right_input+i)=='z'){
						zcount++;
					}
				}
				if(zcount==1) pa->modify_left_primer=1;
				else if(zcount>1) pr_append_new_chunk(glob_err,"Too Many Z in sa->right_input\n");
			}
			if(sa->internal_input){
				length=strlen(sa->internal_input);
				for(i=0;i<length;i++){
					if(*(sa->internal_input+i)=='Z' || *(sa->internal_input+i)=='z'){
						zcount++;
					}
				}
				if(zcount==1) pa->modify_internal_oligo=1;
				else if(zcount>1) pr_append_new_chunk(glob_err,"Too Many Z in sa->internal_input\n");
			}


			/* add check function of consistency between pa->pick_internal_oligo and sa->tar2.genotyping - 20121210 N.Kasahara */
			if(pa->pick_internal_oligo==0 && sa->tar2.genotyping==1){
				(void)printf("WARNING : check pick_internal_oligo, if check Genotyping\n\n");
				pa->pick_internal_oligo=1;	/*	change flag 	*/
			}
		} else pr_append_new_chunk(glob_err, "Unrecognized PRIMER_TASK");
		free(task_tmp);
	}

	/* 
	 * WARNING: read_and_create_seq_lib uses p3_read_line, so repeat
	 * library files cannot be read inside the 
	 * while ((s = p3_read_line(stdin))...) loop above.
	 *
	 * FIX ME, in fact the reading
	 * of the library contents probably belongs inside
	 * primer3_boulder_main.c or libprimer3.c.
	 */

	/* Read in the repeat libraries */
	if (NULL != repeat_file_path) {
		destroy_seq_lib(pa->p_args.repeat_lib);
		if ('\0' == *repeat_file_path) {
			/* Input now specifies no repeat library. */
			pa->p_args.repeat_lib = NULL;
		}
		else {
			pa->p_args.repeat_lib = read_and_create_seq_lib(repeat_file_path, "mispriming library");
			if(pa->p_args.repeat_lib->error.data != NULL) {
				pr_append_new_chunk(glob_err, pa->p_args.repeat_lib->error.data);
			}
		}
		free(repeat_file_path);
		repeat_file_path = NULL;
	}

	/* Read in the repeat libraries for internal oligo */
	if (NULL != int_repeat_file_path) {
		destroy_seq_lib(pa->o_args.repeat_lib);
		if ('\0' == *int_repeat_file_path) {
			/* Input now specifies no mishybridization library. */
			pa->o_args.repeat_lib = NULL;
		}
		else {
			pa->o_args.repeat_lib = 
			  read_and_create_seq_lib(int_repeat_file_path,
			                          "internal oligo mishyb library");
			if(pa->o_args.repeat_lib->error.data != NULL) {
				pr_append_new_chunk(glob_err, pa->o_args.repeat_lib->error.data);
			}
		}
		free(int_repeat_file_path);
		int_repeat_file_path = NULL;
	}

	/* Fix very old tags for backward compatibility */
	if (*io_version == 3) {

		/* This next belongs here rather than libprimer3, because it deals
		   with potential incompatibility between new tags and old tags
		   (which we have kept for backward compatibility).  */
		if (pa->primer_task == pick_pcr_primers_and_hyb_probe) {
			PR_ASSERT(pa->pick_internal_oligo);
		}
		/* Give a error if the tasks don't match */
		if((pick_internal_oligo == 1 || pick_internal_oligo == 0) &&
		   (pa->primer_task == pick_left_only || 
		    pa->primer_task == pick_right_only ||
		    pa->primer_task == pick_hyb_probe_only)) {
			pr_append_new_chunk(glob_err, "Contradiction in primer_task definition");
			if (pa->primer_task == pick_pcr_primers_and_hyb_probe) {
				PR_ASSERT(pa->pick_internal_oligo);
			}
		} else if (pick_internal_oligo == 1) {
			pa->pick_left_primer = 1;
			pa->pick_right_primer = 1;
			pa->pick_internal_oligo = 1;
			if (pa->primer_task == pick_pcr_primers_and_hyb_probe) {
				PR_ASSERT(pa->pick_internal_oligo);
			}
		} else if (pick_internal_oligo == 0) {
			pa->pick_left_primer = 1;
			pa->pick_right_primer = 1;
			pa->pick_internal_oligo = 0;
			if (pa->primer_task == pick_pcr_primers_and_hyb_probe) {
				PR_ASSERT(pa->pick_internal_oligo);
			}
		}
		if (pa->primer_task == pick_pcr_primers_and_hyb_probe) {
			PR_ASSERT(pa->pick_internal_oligo);
		}
	}

	if (pa->primer_task == pick_pcr_primers_and_hyb_probe) {
		PR_ASSERT(pa->pick_internal_oligo);
	}

	/* Removed 10/20/2010 as excessively compulsive, especially in the context
	   of taking input from web pages.
	   if ((min_3_prime || min_5_prime) && (sa->primer_overlap_junctions_count == 0)) {
	   pr_append_new_chunk(warnings,
                         "SEQUENCE_OVERLAP_JUNCTION_LIST not given, but "
	                     "PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION or "
	                     "PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION specified");
	   } */

	return 1;
}
#undef COMPARE
#undef COMPARE_AND_MALLOC
#undef COMPARE_INT
#undef COMPARE_FLOAT
#undef COMPARE_INTERVAL_LIST

static void
pr_append2(pr_append_str *err,
           const char* s1, 
           const char* s2) {
	pr_append_new_chunk(err, s1);
	pr_append(err, s2);
}

int 
read_p3_file(const char *file_name,
             const p3_file_type expected_file_type,
             int echo_output,
             int strict_tags,
             p3_global_settings *pa, 
             seq_args *sa,
             pr_append_str *fatal_err,
             pr_append_str *nonfatal_err,
             pr_append_str *warnings,
             read_boulder_record_results *read_boulder_record_res) 
{
	FILE *file;
	int ret_par = 0;
	int io_version = 4;
	char *line1;
	char *line2;
	char *line3;
	p3_file_type file_type = all_parameters;

	/* Check if a file name was provided */
	PR_ASSERT(NULL != file_name);

	/* Open the file */
	if (!(file = fopen(file_name,"r"))) {
		pr_append2(fatal_err, "Cannot open ", file_name);
		return ret_par;
	}

	/* Parse and interpret the first 3 lines */

	/* Line 1 */
	line1 = p3_read_line(file);
	if (!line1) {
		pr_append2(fatal_err, "Settings file is empty: ", file_name);
		return ret_par;
	}

	/* Line 2 */
	line2 = p3_read_line(file);
	if (!line2) {
		pr_append2(fatal_err, 
		     "Incorrect file format (too few lines) in ", 
		     file_name);
		return ret_par;
	}
	if (!strcmp(line2,"P3_FILE_TYPE=all_parameters")) {
		file_type = all_parameters;
	} else if (!strcmp(line2,"P3_FILE_TYPE=sequence")) {
		file_type = sequence;
	} else if (!strcmp(line2,"P3_FILE_TYPE=settings")) {
		file_type = settings;
	} else {
		pr_append2(fatal_err, "Unknown file type in at line 2 (line2='", line2);
		pr_append(fatal_err, "') in ");
		pr_append(fatal_err, file_name);
		return ret_par;
	}
	if (echo_output) {
		printf("P3_SETTINGS_FILE_USED=%s\n", file_name);
		printf("%s\n", line2);
	}

	/* Line 3 */
	line3 = p3_read_line(file);
	if (!line3) {
		pr_append2(fatal_err, "Incorrect file format (too few lines) in ", file_name);
		return ret_par;
	}
	if (strcmp(line3, "")) {
		pr_append2(fatal_err, "Line 3 must be empty in ", file_name);
		return ret_par;
	}

	/* Check if the file type matches the expected type */
	if (file_type != expected_file_type){
		pr_append_new_chunk(nonfatal_err, "Unexpected P3 file type parsed");
	}

	/* read the file */
	ret_par = read_boulder_record(file, &strict_tags, &io_version, 
                                  echo_output, expected_file_type,
                                  pa, sa, fatal_err, 
                                  nonfatal_err, warnings, 
                                  read_boulder_record_res);

	if (echo_output) printf("P3_SETTINGS_FILE_END=\n");
	if (file) fclose(file);
	return ret_par;
}

static void
tag_syntax_error(const char *tag_name, const char *datum,
                 pr_append_str *err)
{
	pr_append_new_chunk(err, "Illegal ");
	pr_append(err, tag_name);
	pr_append(err, " value: ");
	pr_append(err, datum);
}

static void
parse_double(const char *tag_name, const char *datum,
             double *out, pr_append_str *err)
{
	char *nptr;
	*out = strtod(datum, &nptr);
	if (nptr == datum) {
		/* Empty string or complete junk. */
		tag_syntax_error(tag_name, datum, err);
		*out = 0.0;
		return;
	}
	/* Look for trailing junk. */
	while (*nptr != '\n' && *nptr != '\0') {
		if (*nptr != ' ' && *nptr != '\t') {
			tag_syntax_error(tag_name, datum, err);
			break;
		}
		nptr++;
	}
}

static void
parse_int(const char *tag_name, const char *datum,
          int *out, pr_append_str *err)
{
	char *nptr;
	long tlong;
	tlong = strtol(datum, &nptr, 10);
	if (tlong > INT_MAX || tlong < INT_MIN) {
		tag_syntax_error(tag_name, datum, err);
		pr_append(err, " (value too large or too small)");
		return;
	}
	*out = tlong;
	if (nptr == datum) {
		/* Empty string or complete junk. */
		tag_syntax_error(tag_name, datum, err);
		return;
	}
	/* Look for trailing junk. */
	while (*nptr != '\n' && *nptr != '\0') {
		if (*nptr != ' ' && *nptr != '\t') {
			tag_syntax_error(tag_name, datum, err);
			break;
		}
		nptr++;
	}
}

/* 
 * For correct input, return a pointer to the first non-tab, non-space
 * character after the second integer, and place the integers in out1 and
 * out2.  On incorrect input, return NULL;
 */
static const char *
parse_int_pair(const char    *tag_name, const char *datum,
               char          sep,              /* The separator, e.g. ',' or '-'. */
               int           *out1, int *out2, /* The 2 integers. */
               pr_append_str *err)             /* Error messages. */
{
	char *nptr, *tmp;
	long tlong;

	if(datum[0]==0) return NULL;	/* in case of *datum==NULL - 20130215 N.Kasahara */

	tlong = strtol(datum, &nptr, 10);
	if (tlong > INT_MAX || tlong < INT_MIN) {
		tag_syntax_error(tag_name, datum, err);
		pr_append(err, " (value too large or too small)");
		return NULL;
	}
	*out1 = tlong;
	if (nptr == datum) {
		tag_syntax_error(tag_name, datum, err);
		return NULL;
	}
	while (' ' == *nptr || '\t' == *nptr) nptr++;
	if (sep != *nptr) {
		tag_syntax_error(tag_name, datum, err);
		return NULL;
	}
	nptr++; /* Advance past separator. */
	while (' ' == *nptr || '\t' == *nptr) nptr++;
	tmp = nptr;
	tlong = strtol(tmp, &nptr, 10);
	if (tlong > INT_MAX || tlong < INT_MIN) {
		tag_syntax_error(tag_name, datum, err);
		pr_append(err, " (value too large or too small)");
		return NULL;
	}
	*out2 = tlong;
	if (nptr == tmp) {
		tag_syntax_error(tag_name, datum, err);
		return NULL;
	}
	while (' ' == *nptr || '\t' == *nptr) nptr++;

	/* A hack to live with the old TARGET syntax. */
	if (',' == *nptr && !strcmp(tag_name, "TARGET")) {
		/* Skip the old-fashioned "description". */
		while (' ' != *nptr && '\t' != *nptr && '\0' != *nptr && '\n' != *nptr) nptr++;
		/* Advance to non-space, non-tab. */
		while (' ' == *nptr || '\t' == *nptr) nptr++;
	}
	return nptr;
}


static void
parse_interval_list(const char *tag_name,
                    const char *datum,
                    interval_array_t2 *interval_arr,
                    pr_append_str *err)
{
	const char *p = datum;
	int i1, i2;
	int ret = 0;
	
	if (p == NULL) return; /* in case of *datum==NULL  - 20130215 N.Kasahara  */

	while (' ' == *p || '\t' == *p) p++;
	while (*p != '\0' && *p != '\n') {
		p = parse_int_pair(tag_name, p, ',', &i1, &i2, err);
		if (NULL == p) return;
		ret = p3_add_to_interval_array(interval_arr, i1, i2);
		if (ret) {
			pr_append_new_chunk(err, "Too many elements for tag ");
			pr_append(err, tag_name); return;
		}
	}
}

/* Get variation position to new elements of data structure
   - 20121217 N.Kasahara */
static void
parse_interval_list_VAR(const char *tag_name,
                        const char *datum,
                        interval_array_t2 *interval_arr,
                        pr_append_str *err)
{
	const char *p = datum;
	int i1, i2;
	int ret = 0;
	
	if (p == NULL) return; /* in case of *datum==NULL  - 20130215 N.Kasahara  */

	while (' ' == *p || '\t' == *p) p++;
	while (*p != '\0' && *p != '\n') {
		p = parse_int_pair(tag_name, p, ',', &i1, &i2, err);
		if (p == NULL) return;
		ret = p3_add_to_interval_array_VAR(interval_arr, i1, i2);
		if (ret) {
			pr_append_new_chunk(err, "Too many elements for tag (VAR)");
			pr_append(err,tag_name);
			return;
		}
	}
}



/*
 * For correct input, return a pointer to the first non-tab, non-space
 * character after the forth integer and after the separator sep2, and
 * place the integers in out1, out2, out3 and out4. On incorrect input,
 * return NULL; If any of the 4 integers is not specified, the
 * corresponding output value will be -1.
 */
static char *
parse_2_int_pair(const char    *tag_name, char *datum,
                 char          sep,              /* The separator between 2 numbers, e.g. ',' or '-'. */
                 char          sep2,             /* Separator between 2 intervals */
                 int           *out1, int *out2, 
                 int           *out3, int *out4, /* The 4 integers. */
                 pr_append_str *err)             /* Error messages. */
{
	char *nptr, *tmp;
	long tlong;

	nptr = datum;

	if (*nptr == sep) {
		*out1 = -1;
	} else {
		tlong = strtol(datum, &nptr, 10);
		if (tlong > INT_MAX || tlong < INT_MIN) {
			tag_syntax_error(tag_name, datum, err);
			pr_append(err, " (value too large or too small)");
			return NULL;
		}
		*out1 = tlong;
		if (nptr == datum) {
			tag_syntax_error(tag_name, datum, err);
			return NULL;
		}
	}
	while (' ' == *nptr || '\t' == *nptr) nptr++;
	if (sep != *nptr) {
			tag_syntax_error(tag_name, datum, err);
			return NULL;
	}
	nptr++; /* Advance past separator. */
	while (' ' == *nptr || '\t' == *nptr) nptr++;
	tmp = nptr;
	if (*nptr == sep) {
		*out2 = -1;
	} else {
		tlong = strtol(tmp, &nptr, 10);
		if (tlong > INT_MAX || tlong < INT_MIN) {
			tag_syntax_error(tag_name, datum, err);
			pr_append(err, " (value too large or too small)");
			return NULL;
		}
		*out2 = tlong;
		if (nptr == tmp) {
	tag_syntax_error(tag_name, datum, err);
			return NULL;
		}
	}
	while (' ' == *nptr || '\t' == *nptr) nptr++;
	if (sep != *nptr) {
			tag_syntax_error(tag_name, datum, err);
			return NULL;
	}
	nptr++; /* Advance past separator. */
	while (' ' == *nptr || '\t' == *nptr) nptr++;
	tmp = nptr;
	if (*nptr == sep) {
		*out3 = -1;
	} else {
		tlong = strtol(tmp, &nptr, 10);
		if (tlong > INT_MAX || tlong < INT_MIN) {
			tag_syntax_error(tag_name, datum, err);
			pr_append(err, " (value too large or too small)");
			return NULL;
		}
		*out3 = tlong;
		if (nptr == tmp) {
	tag_syntax_error(tag_name, datum, err);
			return NULL;
		}
	}
	while (' ' == *nptr || '\t' == *nptr) nptr++;
	if (sep != *nptr) {
			tag_syntax_error(tag_name, datum, err);
			return NULL;
	}
	nptr++; /* Advance past separator. */
	while (' ' == *nptr || '\t' == *nptr) nptr++;
	tmp = nptr;
	if ((*nptr == '\n') || (*nptr == '\0') || (*nptr == sep2)) {
		*out4 = -1;
	} else {
		tlong = strtol(tmp, &nptr, 10);
		if (tlong > INT_MAX || tlong < INT_MIN) {
			tag_syntax_error(tag_name, datum, err);
			pr_append(err, " (value too large or too small)");
			return NULL;
		}
		*out4 = tlong;
		if (nptr == tmp) {
	tag_syntax_error(tag_name, datum, err);
			return NULL;
		}
	}
	while (' ' == *nptr || '\t' == *nptr) nptr++;
	/* Should have reached an interval separator or end of line/string */
	if ((*nptr != sep2) && (*nptr != '\0') && (*nptr != '\n')) {
		tag_syntax_error(tag_name, datum, err);
		return NULL;
	}
	if (*nptr == sep2) {
		/* Advance past the interval separator */
		nptr++;
	}
	while (' ' == *nptr || '\t' == *nptr) nptr++;

	return nptr;
}

static void
parse_2_interval_list(const char *tag_name,
                      char *datum,
                      interval_array_t4 *interval_arr,
                      pr_append_str *err)
{
	char *p = datum;
	int i1, i2, i3, i4;
	int ret = 0;

	interval_arr->count = 0;
	interval_arr->any_pair = 0;
	interval_arr->any_left = interval_arr->any_right = 0;

	while (' ' == *p || '\t' == *p) p++;
	while (*p != '\0' && *p != '\n') {
		p = parse_2_int_pair(tag_name, p, ',', ';', &i1, &i2, &i3, &i4, err);
		if (NULL == p) return;
		ret = p3_add_to_2_interval_array(interval_arr, i1, i2, i3, i4);
		if (ret == 1) {
			pr_append_new_chunk(err, "Too many elements for tag ");
			pr_append(err, tag_name); return;
		} else if (ret == 2) {
			pr_append_new_chunk(err, "Invalid range at tag  ");
			pr_append(err, tag_name); return;
		}
	}
}


static int
parse_intron_list(char *s,
                  int *list,
                  int *count) 
{
	long t;
	char *p, *q;

	*count = 0;

	p = q = s;
	
	if(s[0]=='\0') return(INIT_BUF_SIZE);	/* in case of string *datum==NULL - 20130215 N.Kasahara */

	while (*q != '\0' && *q != '\n') {
		t = strtol(p, &q, 10);
		if (q == p) {
			while (*q != '\0') {
				if (!isspace(*q)) {
					*count = 0;
					return 0; 
				}
				q++;
			}
			return *count;
		}
		if (t > INT_MAX || t < INT_MIN) {
			return 0;
		}
		list[*count] = t;
		(*count)++;

		p = q;
	}
	return *count;
}

static void
parse_product_size(const char *tag_name, char *in,
                   p3_global_settings *pa,
                   pr_append_str *err)
{
	char *q, *s = in;
	const char *p;
	int i;
	/* 
	 * Handle possible double quotes around the value.
	 * (This handling is needed for backward compatibility with v2.)
	 */
	if ('"' == *s)  {
		s++;
		in++;
		q = strchr(s, '"');
		if (NULL == q) {
			pr_append_new_chunk(err, tag_name);
			pr_append(err, " begins but does not end with a quote");
			return;
		}
		/* Ignore the " and everything after it. */
		*q = '\0';
	}
	p = in;
	while (' ' == *p || '\t' == *p) p++;
	i = 0;
	while (*p != '\0' && *p != '\n') {
		if (i >= PR_MAX_INTERVAL_ARRAY) {
			pr_append_new_chunk(err, "Too many values for ");
			pr_append(err, tag_name);
			return;
		}
		p = parse_int_pair(tag_name, p, '-',
		                   &pa->pr_min[i], &pa->pr_max[i], err);
		if (NULL == p) return;
		i++;
	}
	pa->num_intervals = i;
}

/* Would be worthwhile to test this a bit more off-by-one errors on
   length of the quality score */

/* This function returns the number of elements in
   the sequence quality vector, and updates
   sargs->quality,  sargs->n_quality,
   and sargs->quality_storage_size */
static int
parse_seq_quality(char *s,
                  seq_args *sargs) 
{
	long t;
	char *p, *q;

	p3_set_sa_empty_quality(sargs);

	p = q = s;

	while (*q != '\0' && *q != '\n') {
		t = strtol(p, &q, 10);
		if (q == p) {
			while (*q != '\0') {
				if (!isspace(*q)) {
					p3_set_sa_empty_quality(sargs);
					return 0; 
				}
				q++;
			}
			return sargs->n_quality;
		}
		p3_sa_add_to_quality_array(sargs, t);

		p = q;
	}
	return sargs->n_quality;
}

/* =========================================================== */
/* Fail-stop wrapper for memory allocation.                    */
/* =========================================================== */
/* 
 * Panic messages for when the program runs out of memory.
 */

static void *
_rb_safe_malloc(size_t x)
{
	void *r = malloc(x);
	if (NULL == r)
			out_of_memory_error();
	return r;
}

static void
pr_append(pr_append_str *x,
          const char *s)
{
	if (pr_append_external(x, s))
			out_of_memory_error();
}

static void
pr_append_new_chunk(pr_append_str *x,
                    const char *s)
{
	if (pr_append_new_chunk_external(x, s))
		out_of_memory_error();
}

static void
out_of_memory_error() 
{
	fprintf(stderr, "out of memory in read_boulder\n");
	exit(-2);
}

/* End of fail-stop wrappers. */
/* =========================================================== */
