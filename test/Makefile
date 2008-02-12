# ======================================================================
# (c) Copyright 1996,1997,1998,1999,2000,2001,2004,2006,2007 Whitehead
# Institute for Biomedical Research, Steve Rozen, and Helen Skaletsky
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
#    * Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above
# copyright notice, this list of conditions and the following disclaimer
# in the documentation and/or other materials provided with the
# distribution.
#    * Neither the names of the copyright holders nor contributors may
# be used to endorse or promote products derived from this software
# without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# ======================================================================

VALGRIND_FLAG = 

test: makeexes testtm testdpal testcore
	-echo 'TESTING ALL:'

testcore:
	-echo 'Testing primer3_core'; perl p3test.pl $(VALGRIND_FLAG)

testdpal:
	-echo 'Testing alignment calcs'; perl dpal_test.pl $(VALGRIND_FLAG); 

testtm:
	-echo 'Testing Tm calcs'; perl oligotm_test.pl $(VALGRIND_FLAG);

makeexes:
	cd ../src; make

clean:
	-rm -rf primer_list_tmp/*
	-rm -rf primer1_list_tmp/*
	-rm -f primer_global_err/*.tmp primer_global_err/*.tmp2 primer_global_err/*.bak
	-rm -f *.tmp *.bak *~
	-rm -f hyb_probe_only.for hyb_probe_only.int 
	-rm -f right_only.for right_only.rev left_only.for
	-rm -f pcr_primers_and_io.int pcr_primers_and_io.rev pcr_primers_and_io.for
	-rm -f pcr_primers.rev pcr_primers.for
	-rm -f oligotm.*.valg  # valgrind output from oligotm_test.pl --valgrind
	-rm -f ntdpal.*.valg   # valgrind output from dpal_test.pl --valgrind
	-rm -f *.vg */*.vg     # valgrind output from p3test{,z}.pl --valgrind
