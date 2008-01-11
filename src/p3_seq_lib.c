/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008
Whitehead Institute for Biomedical Research, Steve Rozen
(http://jura.wi.mit.edu/rozen), and Helen Skaletsky
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

#include <ctype.h>       /* toupper() ... */
#include <string.h>      /* memset(), ... */
#include <setjmp.h>      /* longjmp(), ... */
#include "libprimer3.h"  /* Must include libprimer3.h, which
			    includes p3_seq_lib.h */

static void *p3sl_safe_malloc(size_t x);
static void *p3sl_safe_realloc(void *p, size_t x);
static void  p3sl_append_new_chunk(pr_append_str *x, const char *s);
static void  p3sl_append(pr_append_str *x, const char *s);


static jmp_buf _jmp_buf;

#define P3SL_INIT_BUF_SIZE 1024

#define INIT_LIB_SIZE  500
#define PR_MAX_LIBRARY_WT 100.0

static double parse_seq_name(char *s);
static char   upcase_and_check_char(char *s);
static void   reverse_complement_seq_lib(seq_lib  *lib);

/* See comments in p3_seq_lib.h */
seq_lib *
old_read_and_create_seq_lib(const char * filename, const char *errfrag)
{
    char  *p;
    FILE *file;
    int i, m, k;
    size_t j, n;
    char buf[2];
    char offender = '\0', tmp;
    seq_lib *lib; 

    if (setjmp(_jmp_buf) != 0)
      return NULL; /* If we get here, there was an error in
		      p3sl_safe_malloc or p3sl_safe_realloc. */

    lib =  p3sl_safe_malloc(sizeof(* lib));

    memset(lib, 0, sizeof(*lib));

    lib->repeat_file = p3sl_safe_malloc(strlen(filename) + 1);
    strcpy(lib->repeat_file, filename);

    if((file = fopen(lib->repeat_file,"r")) == NULL) {
	p3sl_append_new_chunk(&lib->error,
			    "Cannot open ");
	goto ERROR;
    }

    j = P3SL_INIT_BUF_SIZE;
    n = INIT_LIB_SIZE;
    lib->names = p3sl_safe_malloc(INIT_LIB_SIZE*sizeof(*lib->names));
    lib->seqs  = p3sl_safe_malloc(INIT_LIB_SIZE*sizeof(*lib->seqs));
    lib->weight= p3sl_safe_malloc(INIT_LIB_SIZE*sizeof(*lib->weight));
    lib->seq_num = 0;
    lib->storage_size = INIT_LIB_SIZE;

    /* Read in the file */
    i = -1;  m = 0; k = 0;
    while((p = p3_read_line(file))) {
	if(*p == '>'){
	    i++;
	    if(i >= n) {
		n += INIT_LIB_SIZE;
		lib->names = p3sl_safe_realloc(lib->names,n*sizeof(*lib->names));
		lib->seqs  = p3sl_safe_realloc(lib->seqs ,n*sizeof(*lib->seqs));
		lib->weight= p3sl_safe_realloc(lib->weight,
					     n*sizeof(*lib->weight));
	    }
	    p++;
	    lib->names[i] = p3sl_safe_malloc(strlen(p) + 1);
	    strcpy(lib->names[i],p);
	    lib->weight[i] = parse_seq_name(lib->names[i]);
	    lib->seqs[i] = p3sl_safe_malloc(P3SL_INIT_BUF_SIZE);
	    lib->seqs[i][0] = '\0';
	    lib->seq_num = i+1;
	    if(lib->weight[i] < 0) {
		p3sl_append_new_chunk(&lib->error, "Illegal weight in ");
		goto ERROR;
	    }
	    j = P3SL_INIT_BUF_SIZE;
	    k = 0;
	    if(i > 0) {
		/* i has already been incremented, so need to use i-1 */
		if(strlen(lib->seqs[i-1]) == 0) {
		    p3sl_append_new_chunk(&lib->error, "Empty sequence in ");
		    goto ERROR;
		}
		tmp = upcase_and_check_char(lib->seqs[i-1]);
		m += tmp;
		if (tmp && '\0' == offender) offender = tmp;
	    }
	    p--;
	}
	else {
	    if(i < 0){ 
		p3sl_append_new_chunk(&lib->error,
				    "Missing id line (expected '>') in ");
		goto ERROR;
	    } else {
		if(k+strlen(p) > j-2){
		    while(j-2 < k+ strlen(p))j += P3SL_INIT_BUF_SIZE;
		    lib->seqs[i] = p3sl_safe_realloc(lib->seqs[i], j);

		}
		strcat(lib->seqs[i], p);
		k += strlen(p);
	    }
	}
    }
    if(i < 0) {
	p3sl_append_new_chunk(&lib->error, "Empty ");
	goto ERROR;
    }
    else if(strlen(lib->seqs[i]) < 3) {
	p3sl_append_new_chunk(&lib->error, "Sequence length < 3 in ");
	goto ERROR;
    }
    tmp = upcase_and_check_char(lib->seqs[i]);
    m += tmp;
    if (tmp && '\0' == offender) offender = tmp;
    if (offender) {
	p3sl_append_new_chunk(&lib->warning,
			    "Unrecognized character (");
	buf[0] = offender;
	buf[1] = '\0';
	p3sl_append(&lib->warning, buf);
	p3sl_append(&lib->warning, ") in ");
	p3sl_append(&lib->warning, errfrag);
	p3sl_append(&lib->warning, " ");
	p3sl_append(&lib->warning, lib->repeat_file);
    }
    if (file) fclose(file);
    reverse_complement_seq_lib(lib);
    return lib;

 ERROR:
    p3sl_append(&lib->error, errfrag);
    p3sl_append(&lib->error, " ");
    p3sl_append(&lib->error, lib->repeat_file);
    if (file) fclose(file);
    return lib;
}

/* FIX ME, add library with weird characters to test the warning */

int
add_seq_to_seq_lib(seq_lib *sl,
		   char *seq, 
		   char *seq_id_plus, 
		   const char *errfrag) {

  int  i = sl->seq_num;
  int  ss = sl->storage_size;
  char offender;
  char buf[2];
  
  /* We need to allocate more storage */
  if (i >= ss) {
    ss += INIT_LIB_SIZE;
    sl->storage_size = ss;
    sl->names = p3sl_safe_realloc(sl->names, ss*sizeof(*sl->names));
    sl->seqs  = p3sl_safe_realloc(sl->seqs , ss*sizeof(*sl->seqs));
    sl->weight= p3sl_safe_realloc(sl->weight,
				   ss*sizeof(*sl->weight));
  }
  sl->seq_num = i + 1;

  sl->names[i] = p3sl_safe_malloc(strlen(seq_id_plus) + 1);
  strcpy(sl->names[i], seq_id_plus);
  sl->weight[i] = parse_seq_name(sl->names[i]);
  if(sl->weight[i] < 0) {
    p3sl_append_new_chunk(&sl->error, "Illegal weight");
    return 1;
  }

  sl->seqs[i] = p3sl_safe_malloc(strlen(seq) + 1);
  strcpy(sl->seqs[i], seq);
  if(strlen(sl->seqs[i]) == 0) {
    p3sl_append_new_chunk(&sl->error, "Empty sequence in ");
    return 1;
  }

  offender = upcase_and_check_char(sl->seqs[i]);
  if ('\0' != offender) {
    buf[0] = offender;
    buf[1] = '\0';
    p3sl_append(&sl->warning, buf);
    p3sl_append(&sl->warning, ") in ");
    p3sl_append(&sl->warning, errfrag);
    p3sl_append(&sl->warning, " ");
  }
    
  return 0;
}

seq_lib *
create_empty_seq_lib() {
  seq_lib *lib;

  if (setjmp(_jmp_buf) != 0)
    return NULL; /* If we get here, there was an error in
		    p3sl_safe_malloc or p3sl_safe_realloc. */

  lib =  p3sl_safe_malloc(sizeof(* lib));
  
  memset(lib, 0, sizeof(*lib));
  lib->repeat_file = NULL;
  lib->names = p3sl_safe_malloc(INIT_LIB_SIZE*sizeof(*lib->names));
  lib->seqs  = p3sl_safe_malloc(INIT_LIB_SIZE*sizeof(*lib->seqs));
  lib->weight= p3sl_safe_malloc(INIT_LIB_SIZE*sizeof(*lib->weight));
  lib->seq_num = 0;
  lib->storage_size = INIT_LIB_SIZE;
  return lib;
}

seq_lib *
read_and_create_seq_lib(const char * filename, const char *errfrag) {
    char  *p;
    FILE *file;
    char *seq_id_plus = NULL;
    char *seq = NULL;;
    size_t seq_storage_size;
    size_t seq_len;

    seq_lib *lib = create_empty_seq_lib();
    if (NULL == lib) return NULL; /* ENOMEM */

    if (setjmp(_jmp_buf) != 0)
      return NULL; /* If we get here, there was an error in
		      p3sl_safe_malloc or p3sl_safe_realloc. */


    lib->repeat_file = p3sl_safe_malloc(strlen(filename) + 1);
    strcpy(lib->repeat_file, filename);

    if((file = fopen(lib->repeat_file,"r")) == NULL) {
	p3sl_append_new_chunk(&lib->error,
			    "Cannot open ");
	goto ERROR;
    }

    seq = p3sl_safe_malloc(P3SL_INIT_BUF_SIZE);
    seq_storage_size = P3SL_INIT_BUF_SIZE;
    seq_len = 0;
    *seq = '\0';

    /* Read in the lines from the file */
    while(1) {
      p = p3_read_line(file);

      if (NULL == p) {
	/* End of file */
	if (seq_id_plus != NULL) {
	    if (seq_len == 0) {
	      p3sl_append_new_chunk(&lib->error,
				    "Empty sequence in ");
	      goto ERROR;
	    } else {
	      if (add_seq_to_seq_lib(lib, seq, seq_id_plus, errfrag)) {
		p3sl_append(&lib->error, " in ");
		goto ERROR;
	      }
	    }
	    free(seq_id_plus);
	    seq_id_plus = NULL;
	}
	break;
      }

      if ('>' == *p) {
	/* We found an id line */
	/* There are two possibilities */
	if (seq_id_plus == NULL) {
	  /* 1. This is the first id line in the file,
	     in which case seq_id_plus == NULL 
	     (and seq_len == 0). */
	  seq_id_plus = p3sl_safe_malloc(strlen(p) + 1);
	  p++;  /* skip past the '>' */
	  strcpy(seq_id_plus, p); 
	} else {
	  /* 2. This is NOT the first id line in the
	     file, in which case seq_id_plus != NULL */
	  if (seq_id_plus != NULL) {
	    if (seq_len == 0) {
	      p3sl_append_new_chunk(&lib->error,
				    "Empty sequence in ");
	      goto ERROR;
	    } else {
	      if (add_seq_to_seq_lib(lib, seq, seq_id_plus, errfrag)) {
		p3sl_append(&lib->error, " in ");
		goto ERROR;
	      }
	      /*"emtpy" the buffer */
	      seq_len = 0;
	      *seq = '\0';
	    }
	    free(seq_id_plus);
	    seq_id_plus = p3sl_safe_malloc(strlen(p));
	    p++;  /* skip past the '>' */
	    strcpy(seq_id_plus, p); 
	  }
	}
      } else {
	/* A sequence line */
	if (seq_id_plus == NULL) {
	  p3sl_append_new_chunk(&lib->error,
				    "Missing id line (expected '>') in ");
	  goto ERROR;
	}
	while ( strlen(p)+ seq_len + 1 > seq_storage_size ) {
	  seq_storage_size *= 2;
	  seq = p3sl_safe_realloc(seq, seq_storage_size);
	}
	strcat(seq, p);
	seq_len += strlen(p);
      }
    }

    if (lib->seq_num == 0) {
      	p3sl_append_new_chunk(&lib->error, "Empty ");
	goto ERROR;
    }
      
    reverse_complement_seq_lib(lib);

    if (file) fclose(file);
    if (seq != NULL) free(seq);
    if (seq_id_plus != NULL) free(seq_id_plus);
    return lib;

 ERROR:
    if (seq != NULL) free(seq);
    if (seq_id_plus != NULL) free(seq_id_plus);
    p3sl_append(&lib->error, errfrag);
    p3sl_append(&lib->error, " ");
    p3sl_append(&lib->error, lib->repeat_file);
    if (file) fclose(file);
    return lib;
}


/* 
 * Free exogenous storage associated with a seq_lib (but not the seq_lib
 * itself).  Silently ignore NULL p.  Set *p to 0 bytes.
 */
void
destroy_seq_lib(seq_lib *p)
{
    int i;
    if (NULL == p) return;

    if ( NULL != p->repeat_file) free(p->repeat_file);
    if (NULL != p->seqs) { 
	for(i = 0; i < p->seq_num; i++)
	    if (NULL != p->seqs[i]) free(p->seqs[i]);
	free(p->seqs);
    }
    if (NULL != p->names) {
	for(i = 0; i < p->seq_num; i++)
	    if (NULL != p->names[i]) free(p->names[i]);
	free(p->names);
    }
    if (NULL != p->weight) free(p->weight);
    if (NULL != p->error.data) free(p->error.data);
    if (NULL != p->warning.data) free(p->warning.data);
    if (NULL != p->rev_compl_seqs) free(p->rev_compl_seqs);
    free(p);
}

static void
reverse_complement_seq_lib(seq_lib  *lib)
{
    int i, n, k;
    if((n = lib->seq_num) == 0) return;
    else {
	lib->names = p3sl_safe_realloc(lib->names, 2*n*sizeof(*lib->names));
	lib->seqs = p3sl_safe_realloc(lib->seqs, 2*n*sizeof(*lib->seqs));
	lib->weight = p3sl_safe_realloc(lib->weight, 2*n*sizeof(*lib->weight));
	lib->rev_compl_seqs = p3sl_safe_malloc(2*n*sizeof(*lib->seqs));

	lib->seq_num *= 2;
	for(i=n; i<lib->seq_num; i++){
	    k = strlen(lib->names[i-n]);
	    lib->names[i] = p3sl_safe_malloc(k + 9);
	    strcpy(lib->names[i], "reverse ");
	    strcat(lib->names[i], lib->names[i-n]);
	    lib->seqs[i] = p3sl_safe_malloc(strlen(lib->seqs[i-n]) + 1);
	    p3_reverse_complement(lib->seqs[i-n], lib->seqs[i]);
	    lib->weight[i] = lib->weight[i-n];
	    lib->rev_compl_seqs[i-n] = lib->seqs[i];
	    lib->rev_compl_seqs[i] = lib->seqs[i-n];
       }
    }
    return;
}

int
seq_lib_num_seq(const seq_lib* lib) {
  if (NULL == lib) return 0;
  else return lib->seq_num;
}

char *
seq_lib_warning_data(const seq_lib *lib) {
  if (NULL == lib) return NULL;
  else return lib->warning.data;
}

static double
parse_seq_name(char *s)
{
    char *p, *q;
    double n;

    p = s;
    while( *p != '*' && *p != '\0' ) p++;
    if (*p == '\0' ) return 1;
    else {
	 p++;
	 n = strtod( p, &q );
	 if( q == p ) return -1;
    }
    if(n > PR_MAX_LIBRARY_WT) return -1;

    return n;
}

/* 
 * Removes spaces and "end-of-line" characters
 * from the sequence, replaces all other
 * characters except A, T, G, C and IUB/IUPAC
 * codes with N.  Returns 0 if there were no such
 * replacements and the first non-ACGT, non-IUB
 * character otherwise. 
 */
static char
upcase_and_check_char(char *s)
{
    int i, j, n, m;

    j = 0; m = 0;
    n = strlen(s);
    for(i=0; i<n; i++){
      
	switch(s[i])
	{
	case 'a' : s[i-j] = 'A'; break;
	case 'g' : s[i-j] = 'G'; break;
	case 'c' : s[i-j] = 'C'; break;
	case 't' : s[i-j] = 'T'; break;
	case 'n' : s[i-j] = 'N'; break;
	case 'A' : s[i-j] = 'A'; break;
	case 'G' : s[i-j] = 'G'; break;
	case 'C' : s[i-j] = 'C'; break;
	case 'T' : s[i-j] = 'T'; break;
	case 'N' : s[i-j] = 'N'; break;

        case 'b' : case 'B': 
        case 'd' : case 'D':
        case 'h' : case 'H':
        case 'v' : case 'V':
        case 'r' : case 'R':
        case 'y' : case 'Y':
        case 'k' : case 'K':
        case 'm' : case 'M':
	case 's' : case 'S':
	case 'w' : case 'W':
	  s[i-j] = toupper(s[i]); break;

	case '\n': j++;          break;
	case ' ' : j++;          break;
	case '\t': j++;          break;
	case '\r': j++;          break;
	default  : if (!m) m = s[i]; s[i-j] = 'N'; 
	}
    }
    s[n-j] = '\0';
    return m;
}

#undef INIT_LIB_SIZE
#undef PR_MAX_LIBRARY_WT

/* =========================================================== */
/* Malloc and realloc wrappers that longjmp() on failure       */
/* =========================================================== */
static void *
p3sl_safe_malloc(size_t x)
{
    void *r = malloc(x);
    if (NULL == r) longjmp(_jmp_buf, 1);
    return r;
}

static void *
p3sl_safe_realloc(void *p, size_t x)
{
    void *r = realloc(p, x);
    if (NULL == r) longjmp(_jmp_buf, 1);
    return r;
}

static void
p3sl_append_new_chunk(pr_append_str *x, const char *s)
{
  int r = pr_append_new_chunk_external(x, s);
  if (r) longjmp(_jmp_buf, 1);
}

static void
p3sl_append(pr_append_str *x, const char *s)
{
  int r = pr_append_external(x, s);
  if (r) longjmp(_jmp_buf, 1);
}
