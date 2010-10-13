/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008,2009,2010
Whitehead Institute for Biomedical Research, Steve Rozen
(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky
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

#ifndef BOULDER_INPUT_H
#define BOULDER_INPUT_H 1
#include "libprimer3.h"

typedef enum p3_file_type {
        all_parameters    = 0,
        sequence          = 1,
        settings          = 2,
} p3_file_type;

typedef struct read_boulder_record_results {
  int explain_flag;
  int file_flag;
} read_boulder_record_results;

extern char *thermodynamic_params_path; /* path to thermodynamic parameter files */
extern int   thermodynamic_path_changed;/* if this is set to 1, we need to re-read the thermodynamic parameters from new path */

/* 
 * Read data from file_input until a "=" line occurs.  Assign
 * parameter values for primer picking to pa and sarg. Perform initial
 * data checking. Return 0 if no records or no _more_ records were
 * found, and 1 otherwise.  If nonfatal_err->data is not NULL or
 * fatal_err->data is not NULL then the data is erroneous and should
 * not be processed. Echo the input lines to stdout.
 */
int read_boulder_record(FILE *file_input,
                        const int *strict_tags,
                        const int * io_version,
                        int echo_output,
                        const p3_file_type read_file_type,
                        p3_global_settings *pa, 
                        seq_args *sarg,
                        pr_append_str *fatal_err,
                        pr_append_str *nonfatal_err,
			pr_append_str *warnings,
                        read_boulder_record_results *);

/* Return null on error. */
/* pr_append_str is an append-only string ADT. */
int read_p3_file(const char *file_name,
                 const p3_file_type file_type,
		 int echo_output,
		 int strict_tags,
                 p3_global_settings *pa, 
                 seq_args *sarg,
                 pr_append_str *fatal_err,
                 pr_append_str *nonfatal_err,
		 pr_append_str *warnings,
                 read_boulder_record_results *read_boulder_record_res);

#endif
