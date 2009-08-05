#!/bin/sh
# ======================================================================
# (c) Copyright 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008,2009
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

# ======================================================================
# CITING PRIMER3
#
# Steve Rozen and Helen J. Skaletsky (2000) Primer3 on the WWW for
# general users and for biologist programmers. In: Krawetz S, Misener S
# (eds) Bioinformatics Methods and Protocols: Methods in Molecular
# Biology. Humana Press, Totowa, NJ, pp 365-386.  Source code available
# from https://sourceforge.net/projects/primer3/
# ======================================================================
#
# This script has to be executed _before_ Makefile is executed. This script 
# determines and thereafter specifies to compiler the location of 
# files of thermodyamical parameters which are needed for calculating
# primers/oligos secondary structures by thermodynamical approach. This script
# can be used for Primer3 version 2.1.0 and newer. If problems occur while
# executing this script then you have probably wrong version of Makefile


# Alternative way for executing this script is to find a row in file 
# Makefile where is written "CC_TH_OPTS = -DPARAMFILES=" and to write
# the full path of location of themrodynamic tables. The format must 
# be like this CC_TH_OPTS = -DPARAMFILES=\"/your/specific/path/primer3/src/tables/\"
# Do not forget the backslashes and quote marks.

srcdir=`pwd`
paramfiles=${srcdir}/tables

makefile="./Makefile"
makefiletmp="./Makefiletmp"

if [ -e $makefiletmp ];
then
rm -rf "${makefiletmp}" 2> /dev/null 
fi

difftail='CC_TH_OPTS = -DPARAMFILES';
diffcount=4;
matchlen=25;

linenr_head=`grep "CC_TH_OPTS = -DPARAMFILES" -n $makefile | cut -f 1 -d \:`
linenr_head2=`expr $linenr_head - 1`
`head \-n $linenr_head2 $makefile > $makefiletmp`

echo "CC_TH_OPTS = -DPARAMFILES=\\\"$paramfiles/\\\"" >> $makefiletmp

linenr_total=`wc \-l $makefile | cut \-f 1 \-d \  `;
linenr_tail=`expr $linenr_total - $linenr_head2 - 1`;
`tail -n $linenr_tail $makefile >> $makefiletmp`

difftail2=`diff \-n $makefiletmp $makefile | tail \-n 1`
matchlen2=`expr match "$difftail2" "$difftail"`
if [ "$matchlen" -ge "$matchlen2" ]; then
  diffcount2=`diff $makefiletmp $makefile | wc \-l`
   if [ "$diffcount" -ge $diffcount2 ]; then
     `mv \-f $makefiletmp $makefile`
     echo "Succeeded.."
   else
    echo "Error: Number of lines different. Read Primer3 manual.\n"
  fi
else 
  echo "Error: Difference in content of wrong line. Read Primer3 manual.\n"
fi
