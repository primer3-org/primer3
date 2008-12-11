/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008
Whitehead Institute for Biomedical Research, Steve Rozen
(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky
All rights reserved.

    This file is part of primer3 and primer3 software suite.

    Primer3 and the primer3 software suite are free software;
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

#ifndef P3_SEQ_LIB_H
#define P3_SEQ_LIB_H

/* WARNING -- #include libprimer3.h, not this file */

/* The seq_lib struct represents a library of sequences. */
/* Clients do not need to know the details of this structure. */
typedef struct seq_lib {
  char **names;         /* An array of sequence names. */
  char **seqs;          /* An array of sequences. */
  char **rev_compl_seqs;/* An array of reversed-complemented sequences.
			   x->rev_compl_seqs[i] is the reverse complement
			   of x->seqs[i], which lets us keep track of pairwise
			   mispriming.  See reverse_complement_seq_lib(). */
  double *weight;       /* An array of weights. */
  char   *repeat_file;  /* The path of the file containing the library. */
  pr_append_str error;  /* Global error message if any.  */
  pr_append_str warning;/* Warning message. */

  /* The number of names, sequences, and weights. */
  int seq_num;

  /* The size of storage allocated for names, sequences, and weights */
  int storage_size;  
} seq_lib;

/* ======================================================= */
/* Functions for creating and destroying a seq_lib object. */
/* ======================================================= */

/*  
 * Reads any file in fasta format and returns a newly allocated
 * seq_lib, lib.  Sets lib.error to a non-empty string on any error
 * other than ENOMEM.  Returns NULL on ENOMEM.
 */
seq_lib *
read_and_create_seq_lib(const char *filename, const char* errfrag);

void
destroy_seq_lib(seq_lib *lib);

/* number of sequences in a seq_lib* */
int
seq_lib_num_seq(const seq_lib* lib);

char *
seq_lib_warning_data(const seq_lib *lib);

int add_seq_and_rev_comp_to_seq_lib(seq_lib *sl,
				    char *seq, 
				    char *seq_id_plus, 
				    const char *errfrag);


#endif

