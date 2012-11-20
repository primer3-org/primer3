# ======================================================================
# (c) Copyright 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008,2009,
#       2010,2011
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

WINMAKE=mingw32-make

# See ../src/Makefile for detailed instructions on how to use TESTOPTS
# The only valid non-empty values for TESTOPTS are --valgrind and --windows
TESTOPTS =


test: makeexes testcmdline testtm testdpal testthal testcore
	echo; echo; echo "Testing completed. You must check the test output for the string FAILED"; echo; echo;

testcmdline:
ifeq ($(TESTOPTS),--windows)
	perl cmdline_test.pl $(TESTOPTS)
else 
	perl cmdline_test.pl $(TESTOPTS)
endif

testcore:
ifeq ($(TESTOPTS),--windows)
	echo 'TESTING primer3_core' & perl p3test.pl $(TESTOPTS)
else
	echo; echo; echo 'TESTING primer3_core'; echo; perl p3test.pl $(TESTOPTS)
endif

testdpal:
ifeq ($(TESTOPTS),--windows)
	echo 'TESTING ALIGNMENT CALCS' & perl dpal_test.pl $(TESTOPTS)
else
	echo; echo; echo 'TESTING ALIGNMENT CALCS'; echo; perl dpal_test.pl $(TESTOPTS); 
endif

testthal:
ifeq ($(TESTOPTS),--windows)
	echo 'TESTING THERMODYNAMIC ALIGNMENT CALCS' & perl thal_test.pl $(TESTOPTS)
else
	echo; echo; echo 'TESTING THERMODYNAMIC ALIGNMENT CALCS'; echo; perl thal_test.pl $(TESTOPTS);
endif

testtm:
ifeq ($(TESTOPTS),--windows)
	perl oligotm_test.pl $(TESTOPTS)
else
	perl oligotm_test.pl $(TESTOPTS)
endif

makeexes:
ifeq ($(TESTOPTS),--windows)
	cd ..\src & $(WINMAKE)
else
	cd ../src; make
endif

clean:
ifeq ($(TESTOPTS),--windows)
	echo 'CLEAN UP TEST DIRECTORIES'
	del /Q /F primer1_list_tmp\*
	del /Q /F primer1_th_list_tmp\*
	del /Q /F primer_list_tmp\*
	del /Q /F th-w-other-tasks_list_tmp\*
	del /Q /F primer_global_err\*.tmp primer_global_err\*.tmp2 primer_global_err\*.bak
	del /Q /F *.tmp *.bak *.tmp.err *~
	del /Q /F hyb_probe_only.for hyb_probe_only.int 
	del /Q /F right_only.for right_only.rev left_only.for
	del /Q /F pcr_primers_and_io.int pcr_primers_and_io.rev pcr_primers_and_io.for
	del /Q /F pcr_primers.rev pcr_primers.for
	del /Q /F pick_primer_list_and_file_flag_th_*.for pick_primer_list_and_file_flag_th_*.rev
	del /Q /F pick_primer_list_and_file_flag_th_*.int
	del /Q /F syntax*.for syntax*.int syntax*.rev
else
	echo; echo; echo 'CLEAN UP TEST DIRECTORIES'; echo; echo;
	-rm -rf primer_list_tmp/*
	-rm -rf primer1_list_tmp/*
	-rm -rf primer1_th_list_tmp/*
	-rm -rf th-w-other-tasks_list_tmp/*
	-rm -f primer_global_err/*.tmp primer_global_err/*.tmp2 primer_global_err/*.bak
	-rm -f *.tmp *.bak *.tmp.err *~
	-rm -f hyb_probe_only.for hyb_probe_only.int 
	-rm -f right_only.for right_only.rev left_only.for
	-rm -f pcr_primers_and_io.int pcr_primers_and_io.rev pcr_primers_and_io.for
	-rm -f pcr_primers.rev pcr_primers.for
	-rm -f pick_primer_list_and_file_flag_th_[1-4].for pick_primer_list_and_file_flag_th_[1-4].rev
	-rm -f pick_primer_list_and_file_flag_th_[1-4].int
	-rm -f syntax[0-1][0-9].for syntax[0-1][0-9].int syntax[0-1][0-9].rev
	-rm -f oligotm.*.valg    # valgrind output from oligotm_test.pl --valgrind
	-rm -f ntdpal.*.valg     # valgrind output from dpal_test.pl --valgrind
	-rm -f ntthal.*.valg     # valgrind output from thal_test.pl --valgrind
	-rm -f oligotm.*.valg.*  # valgrind output from oligotm_test.pl --valgrind
	-rm -f ntdpal.*.valg.*   # valgrind output from dpal_test.pl --valgrind
	-rm -f ntthal.*.valg.*   # valgrind output from thal_test.pl --valgrind
	-rm -f *.vg */*.vg       # valgrind output from p3test{,z}.pl --valgrind
	-rm -f *.vg.core.*       # core files from crashed valgrind test
	-rm -f ./core # in case we crashed and left a core file
endif
