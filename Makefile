#  Copyright (c) 1996
#         Whitehead Institute for Biomedical Research. All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1.      Redistributions must reproduce the above copyright notice, this
# list of conditions and the following disclaimer in the  documentation
# and/or other materials provided with the distribution.  Redistributions of
# source code must also reproduce this information in the source code itself.
# 
# 2.      If the program is modified, redistributions must include a notice
# (in the same places as above) indicating that the redistributed program is
# not identical to the version distributed by Whitehead Institute.
# 
# 3.      All advertising materials mentioning features or use of this
# software  must display the following acknowledgment:
#         This product includes software developed by the
#         Whitehead Institute for Biomedical Research.
# 
# 4.      The name of the Whitehead Institute may not be used to endorse or
# promote products derived from this software  without specific prior written
# permission.
# 
# We also request that use of this software be cited in publications as 
# 
# Steve Rozen, Helen J. Skaletsky (1996)
#    Primer3. Code available at
#    http://www-genome.wi.mit.edu/genome_software/other/primer3.html
# 
# THIS SOFTWARE IS PROVIDED BY THE WHITEHEAD INSTITUTE ``AS IS'' AND  ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE  IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  ARE
# DISCLAIMED. IN NO EVENT SHALL THE WHITEHEAD INSTITUTE BE LIABLE  FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL  DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS  OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)  HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.

# NOTE: The distribution does not contain ntdpal_main.c or oligotm_main.c

MAX_PRIMER_LENGTH = 36

LDLIBS = -lm
CC      = gcc
O_OPTS  = 
CC_OPTS = -g -fPIC
P_DEFINES = -DDPAL_MAX_ALIGN=$(MAX_PRIMER_LENGTH) -DMAX_PRIMER_LENGTH=$(MAX_PRIMER_LENGTH)

CFLAGS  = $(CC_OPTS) $(O_OPTS)
LDFLAGS = -g
#LIBOPTS ='-static'
PRIMER_EXE = primer3_core
PRIMER_LIB = libprimer3.a
RANLIB = ranlib
#RANLIB = @echo

PRIMER_LIBOBJECTS=primer3_lib.o\
	       primer3_release.o\
               oligotm.o\
               dpal_primer.o\
               format_output.o\
               boulder_input.o

PRIMER_OBJECTS=primer3_main.o\
	       $(PRIMER_LIB)

EXES=$(PRIMER_EXE) ntdpal oligotm

all: $(PRIMER_EXE)

clean:
	-rm *.o $(EXES) *~ *.a

$(PRIMER_LIB): $(PRIMER_LIBOBJECTS)
	ar rv $@ $(PRIMER_LIBOBJECTS)
	$(RANLIB) $@
	gcc -shared -Wl,-soname,libprimer.so.1 -o libprimer.so.1.0.1 $(PRIMER_LIBOBJECTS)

$(PRIMER_EXE): $(PRIMER_OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $(PRIMER_OBJECTS) $(LIBOPTS) $(LDLIBS)

# For use with the "testcenter" testing program (CenterLine Software
# Inc, http://www.centerline.com)
$(PRIMER_EXE).tc: $(PRIMER_OBJECTS) /usr/lib/debug/malloc.o
	proof $(CC) $(CFLAGS) -o $@ $(PRIMER_OBJECTS) $(LIBOPTS) $(LDLIBS)

ntdpal: ntdpal_main.o $(PRIMER_LIB)
	$(CC) $(LDFLAGS) -o $@ ntdpal_main.o $(PRIMER_LIB)

oligotm: oligotm_main.o $(PRIMER_LIB)
	$(CC) $(CFLAGS) -o $@ oligotm_main.o $(PRIMER_LIB) $(LIBOPTS) $(LDLIBS)

boulder_input.o: boulder_input.c boulder_input.h primer3.h primer3_release.h dpal.h
	$(CC) -c $(CFLAGS) $(P_DEFINES) -o $@ boulder_input.c

dpal.o: dpal.c dpal.h primer3_release.h
	$(CC) -c $(CFLAGS) -o $@ dpal.c

dpal_primer.o: dpal.c dpal.h primer3_release.h
	$(CC) -c $(CFLAGS) $(P_DEFINES) -o $@ dpal.c

format_output.o: format_output.c primer3_release.h format_output.h primer3.h dpal.h
	$(CC) -c $(CFLAGS) $(P_DEFINES) -o $@ format_output.c

ntdpal_main.o: ntdpal_main.c dpal.h
	$(CC) -c $(CC_OPTS) -o $@ ntdpal_main.c
# We use CC_OPTS above rather than CFLAGS because
# gcc 2.7.2 crashes while compiling ntdpal_main.c with -O2


oligotm.o: oligotm.c oligotm.h primer3_release.h

oligotm_main.o: oligotm_main.c oligotm.h

primer3_main.o: $(PRIMER_LIB)

primer3_lib.o: primer3_lib.c primer3.h primer3_release.h dpal.h oligotm.h format_output.h
	$(CC) -c $(CFLAGS) $(P_DEFINES) primer3_lib.c

primer_test: $(PRIMER_EXE)
	cd ../test; primer_test.pl

backup:
	tar cvf backup.tar Makefile *.[ch] FUNCTIONS
	gzip backup.tar
