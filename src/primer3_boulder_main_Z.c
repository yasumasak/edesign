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

#include <signal.h>
#include <ctype.h>
#include <string.h> /* strcmp() */
#include <stdlib.h> /* free() */
#include <getopt.h> /* getopt() */
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "thal_Z.h"
#include "format_output_Z.h"
#include "libprimer3_Z.h"
#include "read_boulder_Z.h"
#include "print_boulder_Z.h"

/* Check on which OS we compile */
#if defined(_WIN32) || defined(WIN32) || defined (__WIN32__) || defined(__CYGWIN__) || defined(__MINGW32__)
#define OS_WIN
#endif

/* Some function prototypes */
static void   print_usage();
static void   sig_handler(int);
static void   read_thermodynamic_parameters();

/* Other global variables. */
static const char *pr_release;
static const char *pr_program_name;

int
main(int argc, char *argv[])
{
  /* Setup the input data structures handlers */
  int format_output = 0;
  int strict_tags = 0;
  int echo_settings = 0;
  int io_version = 4;
  int default_version = 2;
  int dump_args = 0;
  char *arg_config_path = NULL;

  p3_global_settings *global_pa;
  seq_args *sarg;
  read_boulder_record_results read_boulder_record_res = {0,0};

  pr_append_str p3_settings_path;
  pr_append_str output_path;
  pr_append_str error_path;
  pr_append_str fatal_parse_err;
  pr_append_str nonfatal_parse_err;
  pr_append_str warnings;

  /* Some variables needed by getopt */
  int opt, option_index = 0;
  struct option long_options[] = {
    {"about", no_argument, 0, 'a'},
    {"format_output", no_argument, &format_output, 1},
//    {"strict_tags", no_argument, &strict_tags, 1},
    {"p3_settings_file", required_argument, 0, 'p'},
    {"echo_settings_file", no_argument, &echo_settings, 1},
//    {"io_version", required_argument, 0, 'i'},
//    {"default_version", required_argument, 0, 'd'},
    {"config_path", required_argument, 0, 'c'},
    {"Dump_args", no_argument, 0, 'D'},
    {"2x_compat", no_argument, 0, '2'},
    {"output", required_argument, 0, 'o'},
    {"error", required_argument, 0, 'e'},
    {0, 0, 0, 0}
  };
  int about = 0, compat = 0, invalid_flag = 0;

  /* Retval will point to the return value from choose_primers(). */
  p3retval *retval = NULL;
  int input_found=0;

	/*  variables for variation sequence		*/
	int	length;
	int	m,i,j;

  init_pr_append_str(&fatal_parse_err);
  init_pr_append_str(&nonfatal_parse_err);
  init_pr_append_str(&warnings);
  init_pr_append_str(&p3_settings_path);
  init_pr_append_str(&output_path);
  init_pr_append_str(&error_path);

  /* Get the program name for correct error messages */
  pr_release = libprimer3_release();
  pr_program_name = argv[0];
  p3_set_program_name(pr_program_name);

  /* We set up some signal handlers in case someone starts up the program
   * from the command line, wonders why nothing is happening, and then kills
   * the program. */
  signal(SIGINT, sig_handler);
  signal(SIGTERM, sig_handler);

  /* Read in the flags provided with the program call */
  opterr = 0;
  while ((opt = getopt_long_only(argc, argv, "", long_options, &option_index)) != -1) {
    switch (opt) {
    case 'a':
      about = 1;
      break;
    case 'p':
      if (pr_append_external(&p3_settings_path, optarg)) {
        exit(-2); /* Out of memory. */
      }
      break;
    case 'c':
      arg_config_path = (char*) malloc(strlen(optarg) * sizeof(char) + 1);
      if (NULL == arg_config_path) exit (-2); /* Out of memory */
      strcpy(arg_config_path, optarg);
      break;

/* No io_version selection, io_version = 4
    case 'i':
      if (!strcmp(optarg, "3"))
        io_version = 3;
      else if (!strcmp(optarg, "4"))
        io_version = 4;
      else
        io_version = -1;
      break;
 */

/* No default_version selection, io_version = 2
    case 'd':
      if (!strcmp(optarg, "1"))
        default_version = 1;
      else if (!strcmp(optarg, "2"))
        default_version = 2;
      else
        default_version = -1;
      break;
 */

    case 'D':  /* Undocumented flag for testing; causes 
		  values of arguments to be echoed to
		  stdout. */
      dump_args = 1;
      break;
    case '2':
      compat = 1;
      break;
    case 'o':
      if (pr_append_external(&output_path, optarg)) {
	exit(-2); /* Out of memory. */
      }

      break;
    case 'e':
      if (pr_append_external(&error_path, optarg)) {
	exit(-2); /* Out of memory. */
      }

      break;
    case '?':
      invalid_flag = 1;
      break;
    }
  }
  /* Open the output and error files specified */

  if (!pr_is_empty(&error_path)) {
    /* reassign stderr */
    if (freopen(pr_append_str_chars(&error_path), "w", stderr)
	== NULL) {
      fprintf(stderr, "Error creating file %s\n",
	      pr_append_str_chars(&error_path));
      exit(-1);
    }
    destroy_pr_append_str_data(&error_path);
  }

  if (!pr_is_empty(&output_path)) {
    /* reassign stdout */
    if (freopen(pr_append_str_chars(&output_path), "w", stdout)
	== NULL) {
      fprintf(stderr, "Error creating file %s\n",
	      pr_append_str_chars(&output_path));
      exit(-1);
    }
    destroy_pr_append_str_data(&output_path);
  }
  /* We do any printing after redirecting stdout and stderr */
  if (about == 1) {
    printf("%s\n", pr_release);
    exit(0);
  }
  //if ((io_version == -1) || (invalid_flag == 1) || (default_version == -1)) {
  if ((invalid_flag == 1) || (default_version == -1)) {
    print_usage();
    exit(-1);
  }
  if (compat == 1) {
    printf("PRIMER_ERROR=flag -2x_compat is no longer supported\n=\n");
    exit(-1);
  }

  /* Check if an input file has been specified on the command line. */
  if (optind < argc) {
    if (optind + 1 != argc) {
      print_usage();
      exit(-1);
    }
    if (freopen(argv[optind], "r", stdin) == NULL) {
      fprintf(stderr, "Error opening file %s\n", argv[optind]);
      exit(-1);
    }
  }

  /* Allocate the space for global settings and fill in default parameters */
  if (default_version == 1)
    global_pa = p3_create_global_settings_default_version_1();
  else if (default_version == 2) 
    global_pa = p3_create_global_settings();
  else {
    print_usage();
    exit(-1);
  }
    
  if (!global_pa) {
    exit(-2); /* Out of memory. */
  }
  global_pa->dump = dump_args ;

  /* Settings files have to be read in just below, and
     the functions need a temporary sarg */
  if (!(sarg = create_seq_arg())) {
    exit(-2);
  }

  /* Read data from the settings file until a "=" line occurs.  Assign
     parameter values for primer picking to pa and sa. */
  if (!pr_is_empty(&p3_settings_path)) {
    read_p3_file(pr_append_str_chars(&p3_settings_path),
		 settings,
		 echo_settings && !format_output,
		 strict_tags,
		 global_pa, sarg, &fatal_parse_err,
		 &nonfatal_parse_err, &warnings, &read_boulder_record_res);
    destroy_pr_append_str_data(&p3_settings_path);
    /* Check if the thermodynamical alignment flag was given */
    //if (global_pa->thermodynamic_alignment == 1)
      read_thermodynamic_parameters();
  }

  /* We also need to print out errors here because the loop erases all
     errors at start. If there are fatal errors, write the proper
     message and exit */
  if (fatal_parse_err.data != NULL) {
    if (format_output) {
        format_error(stdout, sarg->sequence_name, fatal_parse_err.data);
    } else {
        print_boulder_error(fatal_parse_err.data);
    }
    fprintf(stderr, "%s: %s\n",
              pr_program_name, fatal_parse_err.data);
    destroy_seq_args(sarg);
    exit(-4);
  }

  /* If there are nonfatal errors, write the proper message
   * and skip to the end of the loop */
  if (!pr_is_empty(&nonfatal_parse_err)) {
    if (format_output) {
        format_error(stdout, sarg->sequence_name,
                     nonfatal_parse_err.data);
    } else {
        print_boulder_error(nonfatal_parse_err.data);
    }
  }

  /* The temporary sarg is not needed any more */
  destroy_seq_args(sarg);
  sarg = NULL;

  /* Read the data from input stream record by record and process it if
   * there are no errors. This is where the work is done. */
  while (1) {
    /* Create and initialize a seq_args data structure. sa (seq_args *) is
     * initialized here because Values are _not_ retained across different
     * input records. */
    if (!(sarg = create_seq_arg())) {
      exit(-2);
    }

    /* Reset all errors handlers and the return structure */
    pr_set_empty(&fatal_parse_err);
    pr_set_empty(&nonfatal_parse_err);
    pr_set_empty(&warnings);
    retval = NULL;



    /* See read_boulder.h for documentation on read_boulder_record().*/
    if (!read_boulder_record(stdin,
			     &strict_tags,
			     &io_version,
			     !format_output,
			     all_parameters,
			     global_pa,
			     sarg,
			     &fatal_parse_err,
			     &nonfatal_parse_err,
			     &warnings,
			     &read_boulder_record_res)) {
      break; /* There were no more boulder records */
    }
		

	/* Error action when genotyping mode is set but target variant is not properly set - 20121206 N.Kasahara */
	if(sarg->tar2.genotyping==1 && sarg->tar2.char_count==0){
		printf("No target variant is correctly set, while \"Genotyping by Internal Probe\" is checked.\n\n");
		printf("Please specify your target variant by using \"[\", \"]\" and \"/\" in the target sequence.\n");
		printf("Different variants should be separated by \"/\".\n");
		exit(-1);
	}

	/* Error action when target variant set but genotyping mode is not set - 20140926 A.Sekiguchi */
	if(sarg->tar2.genotyping==0 && sarg->tar2.char_count>0){
		printf("The target variant is specified, while the \"Genotyping by Internal Probe\" option is unchecked.\n\n");
		printf("Please check the genotyping option when you put a sequence with variant symobls, such as \"[\", \"]\" and \"/\".\n");
		exit(-1);
	}

	/* calculating each Variation sequence
	 * 20121205 - N.Kasahara
	 * 20140108, 20140827 - Y.Kimura modified */
	if(sarg->tar2.genotyping == 1 && sarg->tar2.char_count != 0){
		length = strlen(sarg->sequence);
		sarg->sequence_original = (char *)malloc(length + 1);
		sarg->sequence_var = (char *)malloc(length + 1);
		strcpy(sarg->sequence_original, sarg->sequence);
		strcpy(sarg->sequence_var, sarg->sequence);

		sarg->incl_s = 0;
		sarg->incl_l = -1;

		/* target sequence */
		m = 0;
		(void)memset(sarg->sequence, 0, sizeof(*(sarg->sequence)));
		j = 0;
		for(i = 0; j < sarg->tar2.var_len[m]; i++, j++){
			*(sarg->sequence + i) = *(sarg->sequence_all_var + j);
		}
		*(sarg->sequence + i) = '\0';
		sarg->tar2.VAR_sequence_number = m;

		/* variant sequence */
		m = 1;
		(void)memset(sarg->sequence_var, 0, sizeof(*(sarg->sequence_var)));
		j = sarg->tar2.var_len[m - 1];
		for(i = 0; j < sarg->tar2.var_len[m]; i++, j++){
			*(sarg->sequence_var + i) = *(sarg->sequence_all_var + j);
		}
		*(sarg->sequence_var + i) = '\0';

		/* show target sequence */
		//printf("TARGET SEQUENCE: \n%s\n\n",sarg->sequence);
		//printf("VARIANT SEQUENCE:\n%s\n\n",sarg->sequence_var);
	}else{
		/* show target sequence */
		//printf("TARGET SEQUENCE:\n%s\n\n",sarg->sequence);
	}

	/* Set thermodynamic parameters path by the argument option */
	if (arg_config_path != NULL) {
		thermodynamic_params_path = (char*) malloc(strlen(arg_config_path) * sizeof(char) + 1);
		if (NULL == thermodynamic_params_path) exit (-2); /* Out of memory */
		strcpy(thermodynamic_params_path, arg_config_path);
	}

    /* Check if the thermodynamical alignment flag was given and the
       path to the parameter files changed - we need to reread them */
//    if ((global_pa->thermodynamic_alignment == 1) && (thermodynamic_path_changed == 1))
    if (thermodynamic_path_changed == 1)
      read_thermodynamic_parameters();

    /* Check that we found the thermodynamic parameters in case the flag was set to 1. */
//    if ((global_pa->thermodynamic_alignment == 1) && (thermodynamic_params_path == NULL)) {
    if (thermodynamic_params_path == NULL) {
      /* no parameter directory found, error */
      printf("PRIMER_ERROR=thermodynamic approach chosen, but path to thermodynamic parameters not specified\n=\n");
      exit(-1);
    }

    input_found = 1;
    if ((global_pa->primer_task == generic)
                && (global_pa->pick_internal_oligo == 1)){
      PR_ASSERT(global_pa->pick_internal_oligo);
    }


    /* If there are fatal errors, write the proper message and exit */
    if (fatal_parse_err.data != NULL) {
      if (format_output) {
        format_error(stdout, sarg->sequence_name, fatal_parse_err.data);
      } else {
        print_boulder_error(fatal_parse_err.data);
      }
      fprintf(stderr, "%s: %s\n",
              pr_program_name, fatal_parse_err.data);
      destroy_p3retval(retval);
      destroy_seq_args(sarg);
      exit(-4);
    }


    /* If there are nonfatal errors, write the proper message
     * and skip to the end of the loop */
    if (!pr_is_empty(&nonfatal_parse_err)) {
      if (format_output) {
        format_error(stdout, sarg->sequence_name,
                     nonfatal_parse_err.data);
      } else {
        print_boulder_error(nonfatal_parse_err.data);
      }
      goto loop_wrap_up;
    }

    /* Print any warnings and continue processing */
    if (!pr_is_empty(&warnings)) {
      if (format_output) {
        format_warning(stdout, sarg->sequence_name,
		       warnings.data);
      } else {
        print_boulder_warning(warnings.data);
      }
    }

    if (read_boulder_record_res.file_flag && sarg->sequence_name == NULL) {
      /* We will not have a base name for the files */
      if (format_output) {
        format_error(stdout, NULL,
                     "Need PRIMER_SEQUENCE_ID if PRIMER_FILE_FLAG is not 0");
      } else {
        print_boulder_error("Need PRIMER_SEQUENCE_ID if PRIMER_FILE_FLAG is not 0");
      }
      goto loop_wrap_up;	
    }

    /* Pick the primers - the central function */
    p3_set_gs_primer_file_flag(global_pa,
                               read_boulder_record_res.file_flag);
    retval = choose_primers(global_pa, sarg);


    if (NULL == retval) exit(-2); /* Out of memory. */

    /* This is old code to make it compatible with version 3 input.
       In future versions it can be deleted!
       If it was necessary to use a left_input, right_input,
       or internal_oligo_input primer that was
       unacceptable, then add warnings. */

    if (global_pa->pick_anyway && (io_version == 3
                || format_output)) {
      if (sarg->left_input) {
        add_must_use_warnings(&retval->warnings,
                              "Left primer", &retval->fwd.expl);
      }
      if (sarg->right_input) {
        add_must_use_warnings(&retval->warnings,
                              "Right primer", &retval->rev.expl);
      }
      if (sarg->internal_input) {
        add_must_use_warnings(&retval->warnings,
                              "Internal probe", &retval->intl.expl);
      }
    }
    /* End of the old code for compatibility. */

    if (pr_is_empty(&retval->glob_err)
        && pr_is_empty(&retval->per_sequence_err)) {
      /* We need to test for errors before we call
         p3_print_oligo_lists. This function only works on retval as
         returned when there were no errors. */
      if (read_boulder_record_res.file_flag) {
        /* Create files with left, right, and internal oligos. */
        p3_print_oligo_lists(retval, sarg, global_pa,
                             &retval->per_sequence_err,
                             sarg->sequence_name);
      }
    }

    if (format_output) {
      print_format_output(stdout, &io_version, global_pa,
                          sarg, retval, pr_release,
                          read_boulder_record_res.explain_flag);
    } else {
      /* Use boulder output */
      print_boulder(io_version, global_pa, sarg, retval,
                    read_boulder_record_res.explain_flag);
    }

  loop_wrap_up: /* Here the failed loops join in again */
    if (NULL != retval) {
      /* Check for errors and print them */
      if (NULL != retval->glob_err.data) {
        fprintf(stderr, "%s: %s\n", pr_program_name, retval->glob_err.data);
        destroy_p3retval(retval);
        destroy_seq_args(sarg);
        exit(-4);
      }
    }
//    destroy_p3retval(retval); /* This works even if retval is NULL */
//    retval = NULL;
//    destroy_seq_args(sarg);
//    sarg = NULL;

  }   /* while (1) (processing boulder io records) ...
         End of the primary working loop */

  /* To avoid being distracted when looking for leaks: */
  //if (global_pa->thermodynamic_alignment == 1)
    destroy_thal_structures();
  p3_destroy_global_settings(global_pa);
  global_pa = NULL;
  destroy_seq_args(sarg);
  destroy_pr_append_str_data(&nonfatal_parse_err);
  destroy_pr_append_str_data(&fatal_parse_err);
  destroy_pr_append_str_data(&warnings);
  destroy_dpal_thal_arg_holder();
  free(thermodynamic_params_path);
  /* If it could not read input, then complain and die */
  if (0 == input_found) {
    print_usage();
    exit(-3);
  }
  return 0;
}

/* Reads the thermodynamic parameters if the thermodynamic alignment
   tag was set to 1 */
static void
read_thermodynamic_parameters()
{
  thal_results o;
  /* if the path to the parameter files did not change,
     we do not want to read again */
  if (thermodynamic_path_changed == 0) return;
  /* check that the path to the parameters folder was given */
  if (thermodynamic_params_path == NULL) {

#ifdef OS_WIN
    /* in windows check for .\\primer3_config_Z */
    struct stat st;
    if ((stat(".\\primer3_config_Z", &st) == 0) && S_ISDIR(st.st_mode)) {
      thermodynamic_params_path =
      (char*) malloc(strlen(".\\primer3_config_Z\\") * sizeof(char) + 1);
      if (NULL == thermodynamic_params_path) exit (-2); /* Out of memory */
      strcpy(thermodynamic_params_path, ".\\primer3_config_Z\\");
    } else {
      /* no default directory found */
      return;
    }
#else
    /* in linux, check for ./primer3_config_Z and /opt/primer3_config_Z */
    struct stat st;
    if ((stat("./primer3_config_Z", &st) == 0) && S_ISDIR(st.st_mode)) {
      thermodynamic_params_path =
      (char*) malloc(strlen("./primer3_config_Z/") * sizeof(char) + 1);
      if (NULL == thermodynamic_params_path) exit (-2); /* Out of memory */
      strcpy(thermodynamic_params_path, "./primer3_config_Z/");
    } else if ((stat("/opt/primer3_config_Z", &st) == 0)  && S_ISDIR(st.st_mode)) {
      thermodynamic_params_path =
      (char*) malloc(strlen("/opt/primer3_config_Z/") * sizeof(char) + 1);
      if (NULL == thermodynamic_params_path) exit (-2); /* Out of memory */
      strcpy(thermodynamic_params_path, "/opt/primer3_config_Z/");
    } else {
      /* no default directory found */
      return;
    }
#endif

  }
  /* read in the thermodynamic parameters */
  if (get_thermodynamic_values(thermodynamic_params_path, &o)) {
    fprintf(stderr, "%s\n", o.msg);
    exit(-1);
  }
  /* mark that the last given path was used for reading the parameters */
  thermodynamic_path_changed = 0;
}

/* Print out copyright and a short usage message*/
static void
print_usage()
{
  fprintf(stderr, "%s", edesign_copyright());

  fprintf(stderr, "\n\nUSAGE: %s %s %s %s %s %s %s %s\n", pr_program_name,
		"[-format_output]", 
		"[-p3_settings_file=<file_path>]",
		"[-echo_settings_file]",
		"[-config_path=<directory_path>]",
		"[-output=<file_path>]", 
		"[-error=<file_path>]",
		"[input_file]");
  fprintf(stderr, "\n"
"  -about          \n"
"       show the release number.\n"
"  -format_output  \n"
"       generate user-oriented (rather than program-oriented) output.\n"
"  -p3_settings_file=<file_path>\n"
"       specifies a settings file. Edesign reads the global parameters from this\n"
"       file first. Tags appearing in the settings file override default Edesign\n"
"       settings. The values set in the settings files can be also overridden by\n"
"       tags in the input file.\n"
"  -config_path=<directory_path>\n"
"       specifies the config directory for the thermodynamic parameters.\n"
"       This path overrides the path in the input file. The path should end with \"/\".\n"
"  -echo_settings_file\n"
"       print the input tags specified in the given settings file. \n"
"       If no settings file is provided or if the -format_output option is given, \n"
"       this argument will have no effect.\n"
"  -output=<file_path>\n"
"       the file where the output should be written. If omitted, stdout is used.\n"
"  -error=<file_path>\n"
"       the file where the error messages should be written. If omitted, stderr is used.\n\n");
  fprintf(stderr, "This is Edesign (%s)\n", pr_release);
  fprintf(stderr, "Input can also be provided on standard input.\n\n");
}

/* Print out copyright, a short usage message and the signal */
static void
sig_handler(int signal)
{
    print_usage();
    fprintf(stderr, "%s: received signal %d\n", pr_program_name, signal);
    exit(signal);
}

/* Set the program name */
void
p3_set_program_name(const char *pname)
{
   pr_program_name = pname;
}

