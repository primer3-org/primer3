#!/usr/bin/perl

# ======================================================================
# (c) Copyright 1996,1997,1998,1999,2000,2001,2004,2006,2007,2008,2010,
#  2011,2012,2016
# Whitehead Institute for Biomedical Research, Steve Rozen, 
# Andreas Untergasser and Helen Skaletsky
# All rights reserved.
# 
#   This file is part of the primer3 suite.
#
#   The primer3 suite is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License as
#   published by the Free Software Foundation; either version 2 of the
#   License, or (at your option) any later version.
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

use warnings 'all';
use strict;

if (-d "test") {
    if (-d "test/primer_list_tmp") {
        print "Folder excists: test/primer_list_tmp\n";
    } else {
        mkdir "test/primer_list_tmp";
        print "Created folder: test/primer_list_tmp\n";
    }
    if (-d "test/primer1_list_tmp") {
        print "Folder excists: test/primer1_list_tmp\n";
    } else {
        mkdir "test/primer1_list_tmp";
        print "Created folder: test/primer1_list_tmp\n";
    }
    if (-d "test/primer1_th_list_tmp") {
        print "Folder excists: test/primer1_th_list_tmp\n";
    } else {
        mkdir "test/primer1_th_list_tmp";
        print "Created folder: test/primer1_th_list_tmp\n";
    }
    if (-d "test/th-w-other-tasks_list_tmp") {
        print "Folder excists: test/th-w-other-tasks_list_tmp\n";
    } else {
        mkdir "test/th-w-other-tasks_list_tmp";
        print "Created folder: test/th-w-other-tasks_list_tmp\n";
    }

}
else {
    print "Run this script in the Primer3 root folder\n";
}


