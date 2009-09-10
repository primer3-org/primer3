/*
Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008
Whitehead Institute for Biomedical Research, Steve Rozen
(http://purl.com/STEVEROZEN/), Andreas Untergasser and Helen Skaletsky.
All rights reserved.

    This file is part of the primer3 suite and libraries.

    The primer3 suite and libraries are free software;
    you can redistribute them and/or modify them under the terms
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
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/* 

Test main for function long_seq_tm() in oligotm.c

Usage is e.g 

long_seq_tm_test AAAAGGGCCCCCCCCTTTTTTTTTTT 3 20

In this example 3 is the 0-based start of the
substring to use and 30 is the length of the
substring to use, so Tm is calculated on
GGGCCCCCCCCTTTTTTTTT (= 52.452902).

For testing compare to independent implementation
in ../test/long_seq_tm_test.pl.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "oligotm.h"
  /* double long_seq_tm(const char* s, int start, int len, double salt_conc, double divalent_conc, double dntp_conc); */

int
main(int argc, const char**argv)
{
  const char *s;
  double salt_conc = 50;
  double divalent_conc = 0;
  double dntp_conc = 0;
  double tm;
  int start, len;
  char *endptr;

  s = argv[1]; 
  if (0 == s) {
    fprintf(stderr, "\n%s: incorrect arguments.\n", argv[0]);
    fprintf(stderr, "See file long_seq_tm_test_main.c for usage.\n\n");
    exit(-1);
  }
  start = strtol(argv[2], &endptr, 10); 
  len = strtol(argv[3], &endptr, 10);
  printf("s=%s, start=%d, length=%d\n", s, start, len);
  
  tm = long_seq_tm(s, start, len, salt_conc, divalent_conc, dntp_conc);
  printf("tm = %f\n", tm);
  return 0;
}
