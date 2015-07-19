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

#ifndef PR_FORMAT_OUTPUT_H
#define PR_FORMAT_OUTPUT_H 1

#include "libprimer3_Z.h"
/* 
 * Format the pair in p (plus the middle oligo if appropriate)
 * on f.  If p->product_size = 0 then primer choice failed.
 * Both functions may exit on error, using the conventions
 * in primer3_boulder_main.c
 */

void print_format_output(FILE *f, 
                         const int *io_version, 
                         const p3_global_settings *pa,
                         const seq_args *sa, 
                         const p3retval *, 
                         const char *,
                         int   explain_flag);


void format_error(FILE *f, const char* seq_name, const char *err);
void format_warning(FILE *f, const char* seq_name, const char *err);

#endif
