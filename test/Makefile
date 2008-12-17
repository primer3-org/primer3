# ======================================================================
# (c) Copyright 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008
# Whitehead Institute for Biomedical Research, Steve Rozen, 
# Andreas Untergasser and Helen Skaletsky
# All rights reserved.
# 
#   This file is part of primer3, the libprimer3 library, the oligotm 
#   library and the dpal library.
#
#   Primer3 and the libraries above are free software; you can
#   redistribute them and/or modify them under the terms of the GNU
#   General Public License as published by the Free Software Foundation;
#   either version 2 of the License, or (at your option) any later
#   version.
#
#   This software is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this file (file gpl-2.0.txt in the source distribution); if
#   not, write to the Free Software Foundation, Inc., 51 Franklin St,
#   Fifth Floor, Boston, MA 02110-1301 USA
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

# See ../src/Makefile for detailed instructions on how to use VAL
# The only valid non-empty value for VAL is --valgrind
VAL =


test: makeexes testtm testdpal testcore test_io_v3

test_io_v3:
	cd io_version_3; make test

testcore:
	echo; echo; echo 'TESTING primer3_core'; echo; perl p3test.pl $(VAL)

testdpal:
	echo; echo; echo 'TESTING ALIGNMENT CALCS'; echo; perl dpal_test.pl $(VAL); 

testtm:
	echo; echo; echo 'TESTING Tm CALCS'; echo; perl oligotm_test.pl $(VAL);

makeexes:
	cd ../src; make

clean:
	echo; echo; echo 'CLEAN UP TEST DIRECTORIES'; echo; echo;
	-rm -rf primer_list_tmp/*
	-rm -rf primer1_list_tmp/*
	-rm -f primer_global_err/*.tmp primer_global_err/*.tmp2 primer_global_err/*.bak
	-rm -f *.tmp *.bak *~
	-rm -f hyb_probe_only.for hyb_probe_only.int 
	-rm -f right_only.for right_only.rev left_only.for
	-rm -f pcr_primers_and_io.int pcr_primers_and_io.rev pcr_primers_and_io.for
	-rm -f pcr_primers.rev pcr_primers.for
	-rm -f syntax[0-1][0-9].for syntax[0-1][0-9].int syntax[0-1][0-9].rev
	-rm -f oligotm.*.valg  # valgrind output from oligotm_test.pl --valgrind
	-rm -f ntdpal.*.valg   # valgrind output from dpal_test.pl --valgrind
	-rm -f *.vg */*.vg     # valgrind output from p3test{,z}.pl --valgrind
	-rm -f *.vg.core.*     # core files from crashed valgrind test
	-rm -f ./core # in case we crashed and left a core file
	cd io_version_3; make clean
