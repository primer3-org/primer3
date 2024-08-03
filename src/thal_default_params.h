/*
This file is created by thal_default_params_create.py. Modify that script, not this file.
requires python3, python3-numpy
run "python3 thal_default_params_create.py" to regenerate this file.

Globals initialize in this file:
const double _INFINITY;
static double atpS[5][5]; AT penalty 
static double atpH[5][5];  AT penalty 
static int numTriloops;  hairpin triloop penalties 
static int numTetraloops;  hairpin tetraloop penalties 
static double dangleEntropies3[5][5][5]; thermodynamic paramteres for 3' dangling ends 
static double dangleEnthalpies3[5][5][5];  ther params for 3' dangling ends 
static double dangleEntropies5[5][5][5];   ther params for 5' dangling ends 
static double dangleEnthalpies5[5][5][5];  ther params for 5' dangling ends 
static double stackEntropies[5][5][5][5];  ther params for perfect match pairs 
static double stackEnthalpies[5][5][5][5];  ther params for perfect match pairs 
static double stackint2Entropies[5][5][5][5]; ther params for perfect match and internal mm 
static double stackint2Enthalpies[5][5][5][5];  ther params for perfect match and internal mm
static double interiorLoopEntropies[30];  interior loop params according to length of the loop 
static double bulgeLoopEntropies[30];  bulge loop params according to length of the loop 
static double hairpinLoopEntropies[30];  hairpin loop params accordint to length of the loop 
static double interiorLoopEnthalpies[30];  same as interiorLoopEntropies but values of entropy 
static double bulgeLoopEnthalpies[30];  same as bulgeLoopEntropies but values of entropy 
static double hairpinLoopEnthalpies[30];  same as hairpinLoopEntropies but values of entropy 
static double tstackEntropies[5][5][5][5];  ther params for terminal mismatches 
static double tstackEnthalpies[5][5][5][5];  ther params for terminal mismatches 
static double tstack2Entropies[5][5][5][5];  ther params for internal terminal mismatches 
static double tstack2Enthalpies[5][5][5][5];  ther params for internal terminal mismatches 
static struct triloop* triloopEntropies;  ther penalties for given triloop seq-s 
static struct triloop* triloopEnthalpies;  ther penalties for given triloop seq-s 
static struct tetraloop* tetraloopEntropies;  ther penalties for given tetraloop seq-s 
static struct tetraloop* tetraloopEnthalpies;  ther penalties for given tetraloop seq-s 
*/

#include <math.h>
#include <stdio.h>
#include "thal.h"


# ifdef INTEGER
const double _INFINITY = 999999.0;
# else
# ifdef INFINITY
const double _INFINITY = INFINITY;
# else
const double _INFINITY = 1.0 / 0.0;
# endif
# endif

static double atpS[5][5] = {
	{0.00000000001, 0.00000000001, 0.00000000001, 6.9, 0.00000000001},
	{0.00000000001, 0.00000000001, 0.00000000001, 0.00000000001, 0.00000000001},
	{0.00000000001, 0.00000000001, 0.00000000001, 0.00000000001, 0.00000000001},
	{6.9, 0.00000000001, 0.00000000001, 0.00000000001, 0.00000000001},
	{0.00000000001, 0.00000000001, 0.00000000001, 0.00000000001, 0.00000000001}};

static double atpH[5][5] = {
	{0.0, 0.0, 0.0, 2200.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{2200.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0}};

static double stackEntropies[5][5][5][5] = {
	{{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -22.2, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -22.4, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -21.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-20.4, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}},

	{{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -22.7, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -19.9, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -27.2, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-21.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}},

	{{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -22.2, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -24.4, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -19.9, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-22.4, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}},

	{{{-1.0, -1.0, -1.0, -21.3, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -22.2, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -22.7, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-22.2, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}},

	{{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}}};

static double stackEnthalpies[5][5][5][5] = {
	{{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, -7900.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, -8400.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, -7800.0, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{-7200.0, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}},

	{{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, -8500.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, -8000.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, -10600.0, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{-7800.0, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}},

	{{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, -8200.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, -9800.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, -8000.0, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{-8400.0, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}},

	{{{_INFINITY, _INFINITY, _INFINITY, -7200.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, -8200.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, -8500.0, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{-7900.0, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}},

	{{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}}};

static double stackint2Entropies[5][5][5][5] = {
	{{{-1.0, -1.0, -1.0, 12.9, -1.0},
	{-1.0, -1.0, -1.0, 20.2, -1.0},
	{-1.0, -1.0, -1.0, 7.4, -1.0},
	{1.7, 4.6, -2.3, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -9.8, -1.0, -1.0},
	{-1.0, -1.0, -3.8, -1.0, -1.0},
	{-1.0, -1.0, 3.2, -1.0, -1.0},
	{14.6, -4.4, -1.0, 0.2, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -4.2, -1.0, -1.0, -1.0},
	{-1.0, -0.6, -1.0, -1.0, -1.0},
	{-1.0, -13.2, -1.0, -1.0, -1.0},
	{-2.3, -1.0, -9.5, 0.9, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{1.7, -1.0, -1.0, -1.0, -1.0},
	{14.6, -1.0, -1.0, -1.0, -1.0},
	{-2.3, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -6.2, -8.3, -10.8, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}},

	{{{-1.0, -1.0, -1.0, 8.0, -1.0},
	{-1.0, -1.0, -1.0, 16.4, -1.0},
	{-4.2, 3.7, -2.3, -1.0, -1.0},
	{-1.0, -1.0, -1.0, 0.7, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, 14.2, -1.0, -1.0},
	{-1.0, -1.0, 8.9, -1.0, -1.0},
	{-0.6, -7.2, -1.0, -4.5, -1.0},
	{-1.0, -1.0, 13.5, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, 3.7, -1.0, -1.0, -1.0},
	{-1.0, -7.2, -1.0, -1.0, -1.0},
	{-13.2, -1.0, -15.3, -11.7, -1.0},
	{-1.0, -6.1, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{4.6, -1.0, -1.0, -1.0, -1.0},
	{-4.4, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -6.1, -8.0, -15.8, -1.0},
	{-6.2, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}},

	{{{-1.0, -1.0, -1.0, 0.7, -1.0},
	{-9.8, 14.2, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, 3.6, -1.0},
	{-1.0, -1.0, -1.0, -5.3, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-3.8, 8.9, -1.0, 5.4, -1.0},
	{-1.0, -1.0, -15.8, -1.0, -1.0},
	{-1.0, -1.0, -12.3, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -2.3, -1.0, -1.0, -1.0},
	{3.2, -1.0, -15.8, 10.4, -1.0},
	{-1.0, -15.3, -1.0, -1.0, -1.0},
	{-1.0, -8.0, -1.0, 16.3, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-2.3, -1.0, -1.0, -1.0, -1.0},
	{-1.0, 13.5, -12.3, -8.4, -1.0},
	{-9.5, -1.0, -1.0, -1.0, -1.0},
	{-8.3, -1.0, 9.5, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}},

	{{{-1.0, 8.0, 0.7, -1.0, -1.0},
	{-1.0, -1.0, -1.0, 0.7, -1.0},
	{-1.0, -1.0, -1.0, -1.7, -1.0},
	{-1.0, -1.0, -1.0, -1.5, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{20.2, 16.4, -1.0, 0.7, -1.0},
	{-1.0, -1.0, 5.4, -1.0, -1.0},
	{-1.0, -1.0, 10.4, -1.0, -1.0},
	{-1.0, -1.0, -8.4, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{7.4, -1.0, 3.6, -1.7, -1.0},
	{-1.0, -4.5, -1.0, -1.0, -1.0},
	{-1.0, -11.7, -1.0, -6.2, -1.0},
	{-1.0, -15.8, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, 0.7, -5.3, -1.5, -1.0},
	{0.2, -1.0, -1.0, -1.0, -1.0},
	{0.9, -1.0, 16.3, -1.0, -1.0},
	{-10.8, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}},

	{{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}}};

static double stackint2Enthalpies[5][5][5][5] = {
	{{{_INFINITY, _INFINITY, _INFINITY, 4700.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, 7600.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, 3000.0, _INFINITY},
	{1200.0, 2300.0, -600.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, -2900.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, -700.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, 500.0, _INFINITY, _INFINITY},
	{5300.0, -10.0, _INFINITY, 700.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, -900.0, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, 600.0, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, -4000.0, _INFINITY, _INFINITY, _INFINITY},
	{-700.0, _INFINITY, -3100.0, 1000.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{1200.0, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{5300.0, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{-700.0, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, -1200.0, -2500.0, -2700.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}},

	{{{_INFINITY, _INFINITY, _INFINITY, 3400.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, 6100.0, _INFINITY},
	{-900.0, 1900.0, -700.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, 1000.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, 5200.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, 3600.0, _INFINITY, _INFINITY},
	{600.0, -1500.0, _INFINITY, -800.0, _INFINITY},
	{_INFINITY, _INFINITY, 5200.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, 1900.0, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, -1500.0, _INFINITY, _INFINITY, _INFINITY},
	{-4000.0, _INFINITY, -4900.0, -4100.0, _INFINITY},
	{_INFINITY, -1500.0, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{2300.0, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{-10.0, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, -1500.0, -2800.0, -5000.0, _INFINITY},
	{-1200.0, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}},

	{{{_INFINITY, _INFINITY, _INFINITY, 700.0, _INFINITY},
	{-2900.0, 5200.0, -600.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, 1600.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, -1300.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, -600.0, _INFINITY, _INFINITY},
	{-700.0, 3600.0, _INFINITY, 2300.0, _INFINITY},
	{_INFINITY, _INFINITY, -6000.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, -4400.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, -700.0, _INFINITY, _INFINITY, _INFINITY},
	{500.0, _INFINITY, -6000.0, 3300.0, _INFINITY},
	{_INFINITY, -4900.0, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, -2800.0, _INFINITY, 5800.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{-600.0, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, 5200.0, -4400.0, -2200.0, _INFINITY},
	{-3100.0, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{-2500.0, _INFINITY, 4100.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}},

	{{{_INFINITY, 3400.0, 700.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, 1200.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, -100.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, 200.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{7600.0, 6100.0, _INFINITY, 1200.0, _INFINITY},
	{_INFINITY, _INFINITY, 2300.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, 3300.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, -2200.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{3000.0, _INFINITY, 1600.0, -100.0, _INFINITY},
	{_INFINITY, -800.0, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, -4100.0, _INFINITY, -1400.0, _INFINITY},
	{_INFINITY, -5000.0, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, 1000.0, -1300.0, 200.0, _INFINITY},
	{700.0, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{1000.0, _INFINITY, 5800.0, _INFINITY, _INFINITY},
	{-2700.0, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}},

	{{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}}};

static double tstackEntropies[5][5][5][5] = {
	{{{-1.0, -1.0, -1.0, -6.3, 1e-11},
	{-1.0, -1.0, -1.0, -7.0, 1e-11},
	{-1.0, -1.0, -1.0, -5.8, 1e-11},
	{-7.8, -4.0, -4.4, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -22.5, -1.0, 1e-11},
	{-1.0, -1.0, -7.1, -1.0, 1e-11},
	{-1.0, -1.0, -11.4, -1.0, 1e-11},
	{-3.8, -0.5, -1.0, -1.7, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -10.7, -1.0, -1.0, 1e-11},
	{-1.0, -6.0, -1.0, -1.0, 1e-11},
	{-1.0, -15.5, -1.0, -1.0, 1e-11},
	{-5.9, -1.0, -2.1, -8.7, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-7.8, -1.0, -1.0, -1.0, 1e-11},
	{-3.8, -1.0, -1.0, -1.0, 1e-11},
	{-5.9, -1.0, -1.0, -1.0, 1e-11},
	{-1.0, -6.3, -9.4, -6.5, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}},

	{{{-1.0, -1.0, -1.0, -5.9, 1e-11},
	{-1.0, -1.0, -1.0, -1.3, 1e-11},
	{-10.7, -5.9, -9.6, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -1.2, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -13.8, -1.0, 1e-11},
	{-1.0, -1.0, -10.6, -1.0, 1e-11},
	{-6.0, -5.1, -1.0, -8.0, 1e-11},
	{-1.0, -1.0, -7.8, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -5.9, -1.0, -1.0, 1e-11},
	{-1.0, -5.1, -1.0, -1.0, 1e-11},
	{-15.5, -1.0, -9.5, -9.0, 1e-11},
	{-1.0, -10.6, -1.0, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-4.0, -1.0, -1.0, -1.0, 1e-11},
	{-0.5, -1.0, -1.0, -1.0, 1e-11},
	{-1.0, -10.6, -18.7, -16.9, 1e-11},
	{-6.3, -1.0, -1.0, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}},

	{{{-1.0, -1.0, -1.0, -4.7, 1e-11},
	{-22.5, -13.8, -11.1, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -2.7, 1e-11},
	{-1.0, -1.0, -1.0, -9.8, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -11.1, -1.0, 1e-11},
	{-7.1, -10.6, -1.0, -13.5, 1e-11},
	{-1.0, -1.0, -19.2, -1.0, 1e-11},
	{-1.0, -1.0, -16.1, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -9.6, -1.0, -1.0, 1e-11},
	{-11.4, -1.0, -19.2, -15.9, 1e-11},
	{-1.0, -9.5, -1.0, -1.0, 1e-11},
	{-1.0, -18.7, -1.0, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-4.4, -1.0, -1.0, -1.0, 1e-11},
	{-1.0, -7.8, -16.1, -21.2, 1e-11},
	{-2.1, -1.0, -1.0, -1.0, 1e-11},
	{-9.4, -1.0, -1.0, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}},

	{{{-6.3, -5.9, -4.7, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -6.3, 1e-11},
	{-1.0, -1.0, -1.0, -10.5, 1e-11},
	{-1.0, -1.0, -1.0, -8.9, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-7.0, -1.3, -1.0, -6.3, 1e-11},
	{-1.0, -1.0, -13.5, -1.0, 1e-11},
	{-1.0, -1.0, -15.9, -1.0, 1e-11},
	{-1.0, -1.0, -21.2, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-5.8, -1.0, -2.7, -10.5, 1e-11},
	{-1.0, -8.0, -1.0, -1.0, 1e-11},
	{-1.0, -9.0, -1.0, -1.0, 1e-11},
	{-1.0, -16.9, -1.0, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.2, -9.8, -8.9, 1e-11},
	{-1.7, -1.0, -1.0, -1.0, 1e-11},
	{-8.7, -1.0, -1.0, -1.0, 1e-11},
	{-6.5, -1.0, -1.0, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}},

	{{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}}};

static double tstackEnthalpies[5][5][5][5] = {
	{{{_INFINITY, _INFINITY, _INFINITY, -2500.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, -2700.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, -2400.0, 0.0},
	{-3100.0, -1600.0, -1900.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, -8000.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, -3200.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, -4600.0, _INFINITY, 0.0},
	{-1800.0, -100.0, _INFINITY, -900.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, -4300.0, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, -2700.0, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, -6000.0, _INFINITY, _INFINITY, 0.0},
	{-2500.0, _INFINITY, -1100.0, -3200.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{-3100.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{-1800.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{-2500.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, -2300.0, -3500.0, -2400.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}},

	{{{_INFINITY, _INFINITY, _INFINITY, -2300.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, -700.0, 0.0},
	{-4300.0, -2600.0, -3900.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, -700.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, -5000.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, -3900.0, _INFINITY, 0.0},
	{-2700.0, -2100.0, _INFINITY, -3200.0, 0.0},
	{_INFINITY, _INFINITY, -3000.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, -2600.0, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, -2100.0, _INFINITY, _INFINITY, 0.0},
	{-6000.0, _INFINITY, -3800.0, -3800.0, 0.0},
	{_INFINITY, -3900.0, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{-1600.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{-100.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, -3900.0, -6600.0, -6100.0, 0.0},
	{-2300.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}},

	{{{_INFINITY, _INFINITY, _INFINITY, -2000.0, 0.0},
	{-8000.0, -5000.0, -4300.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, -1100.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, -3600.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, -4300.0, _INFINITY, 0.0},
	{-3200.0, -3900.0, _INFINITY, -4900.0, 0.0},
	{_INFINITY, _INFINITY, -700.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, -5900.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, -3900.0, _INFINITY, _INFINITY, 0.0},
	{-4600.0, _INFINITY, -700.0, -5700.0, 0.0},
	{_INFINITY, -3800.0, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, -6600.0, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{-1900.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, -3000.0, -5900.0, -7400.0, 0.0},
	{-1100.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{-3500.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}},

	{{{-2500.0, -2300.0, -2000.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, -2500.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, -3900.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, -3200.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{-2700.0, -700.0, _INFINITY, -2500.0, 0.0},
	{_INFINITY, _INFINITY, -4900.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, -5700.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, -7400.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{-2400.0, _INFINITY, -1100.0, -3900.0, 0.0},
	{_INFINITY, -3200.0, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, -3800.0, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, -6100.0, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, -700.0, -3600.0, -3200.0, 0.0},
	{-900.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{-3200.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{-2400.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}},

	{{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}}};

static double tstack2Entropies[5][5][5][5] = {
	{{{-1.0, -1.0, -1.0, -6.3, 1e-11},
	{-1.0, -1.0, -1.0, -7.0, 1e-11},
	{-1.0, -1.0, -1.0, -5.8, 1e-11},
	{-7.8, -4.0, -4.4, -13.5, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -22.5, -1.0, 1e-11},
	{-1.0, -1.0, -7.1, -1.0, 1e-11},
	{-1.0, -1.0, -11.4, -1.0, 1e-11},
	{-3.8, -0.5, -16.1, -1.7, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -10.7, -1.0, -1.0, 1e-11},
	{-1.0, -6.0, -1.0, -1.0, 1e-11},
	{-1.0, -15.5, -1.0, -1.0, 1e-11},
	{-5.9, -16.1, -2.1, -8.7, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-7.8, -1.0, -1.0, -1.0, 1e-11},
	{-3.8, -1.0, -1.0, -1.0, 1e-11},
	{-5.9, -1.0, -1.0, -1.0, 1e-11},
	{-13.6, -6.3, -9.4, -6.5, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}},

	{{{-1.0, -1.0, -1.0, -5.9, 1e-11},
	{-1.0, -1.0, -1.0, -1.3, 1e-11},
	{-10.7, -5.9, -9.6, -16.1, 1e-11},
	{-1.0, -1.0, -1.0, -1.2, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -13.8, -1.0, 1e-11},
	{-1.0, -1.0, -10.6, -1.0, 1e-11},
	{-6.0, -5.1, -19.3, -8.0, 1e-11},
	{-1.0, -1.0, -7.8, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -5.9, -1.0, -1.0, 1e-11},
	{-1.0, -5.1, -1.0, -1.0, 1e-11},
	{-15.5, -19.3, -9.5, -9.0, 1e-11},
	{-1.0, -10.6, -1.0, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-4.0, -1.0, -1.0, -1.0, 1e-11},
	{-0.5, -1.0, -1.0, -1.0, 1e-11},
	{-16.1, -10.6, -18.7, -16.9, 1e-11},
	{-6.3, -1.0, -1.0, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}},

	{{{-1.0, -1.0, -1.0, -4.7, 1e-11},
	{-22.5, -13.8, -11.1, -16.1, 1e-11},
	{-1.0, -1.0, -1.0, -2.7, 1e-11},
	{-1.0, -1.0, -1.0, -9.8, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -11.1, -1.0, 1e-11},
	{-7.1, -10.6, -19.3, -13.5, 1e-11},
	{-1.0, -1.0, -19.2, -1.0, 1e-11},
	{-1.0, -1.0, -16.1, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -9.6, -1.0, -1.0, 1e-11},
	{-11.4, -19.3, -19.2, -15.9, 1e-11},
	{-1.0, -9.5, -1.0, -1.0, 1e-11},
	{-1.0, -18.7, -1.0, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-4.4, -1.0, -1.0, -1.0, 1e-11},
	{-16.1, -7.8, -16.1, -21.2, 1e-11},
	{-2.1, -1.0, -1.0, -1.0, 1e-11},
	{-9.4, -1.0, -1.0, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}},

	{{{-6.3, -5.9, -4.7, -14.2, 1e-11},
	{-1.0, -1.0, -1.0, -6.3, 1e-11},
	{-1.0, -1.0, -1.0, -10.5, 1e-11},
	{-1.0, -1.0, -1.0, -8.9, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-7.0, -1.3, -16.1, -6.3, 1e-11},
	{-1.0, -1.0, -13.5, -1.0, 1e-11},
	{-1.0, -1.0, -15.9, -1.0, 1e-11},
	{-1.0, -1.0, -21.2, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-5.8, -16.1, -2.7, -10.5, 1e-11},
	{-1.0, -8.0, -1.0, -1.0, 1e-11},
	{-1.0, -9.0, -1.0, -1.0, 1e-11},
	{-1.0, -16.9, -1.0, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-13.5, -1.2, -9.8, -8.9, 1e-11},
	{-1.7, -1.0, -1.0, -1.0, 1e-11},
	{-8.7, -1.0, -1.0, -1.0, 1e-11},
	{-6.5, -1.0, -1.0, -1.0, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{1e-11, 1e-11, 1e-11, 1e-11, 1e-11},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}},

	{{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}}};

static double tstack2Enthalpies[5][5][5][5] = {
	{{{_INFINITY, _INFINITY, _INFINITY, -2500.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, -2700.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, -2400.0, 0.0},
	{-3100.0, -1600.0, -1900.0, -5000.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, -8000.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, -3200.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, -4600.0, _INFINITY, 0.0},
	{-1800.0, -100.0, -6000.0, -900.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, -4300.0, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, -2700.0, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, -6000.0, _INFINITY, _INFINITY, 0.0},
	{-2500.0, -6000.0, -1100.0, -3200.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{-3100.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{-1800.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{-2500.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{-5000.0, -2300.0, -3500.0, -2400.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}},

	{{{_INFINITY, _INFINITY, _INFINITY, -2300.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, -700.0, 0.0},
	{-4300.0, -2600.0, -3900.0, -6000.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, -700.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, -5000.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, -3900.0, _INFINITY, 0.0},
	{-2700.0, -2100.0, -7000.0, -3200.0, 0.0},
	{_INFINITY, _INFINITY, -3000.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, -2600.0, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, -2100.0, _INFINITY, _INFINITY, 0.0},
	{-6000.0, -7000.0, -3800.0, -3800.0, 0.0},
	{_INFINITY, -3900.0, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{-1600.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{-100.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{-6000.0, -3900.0, -6600.0, -6100.0, 0.0},
	{-2300.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}},

	{{{_INFINITY, _INFINITY, _INFINITY, -2000.0, 0.0},
	{-8000.0, -5000.0, -4300.0, -6000.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, -1100.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, -3600.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, -4300.0, _INFINITY, 0.0},
	{-3200.0, -3900.0, -7000.0, -4900.0, 0.0},
	{_INFINITY, _INFINITY, -700.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, -5900.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, -3900.0, _INFINITY, _INFINITY, 0.0},
	{-4600.0, -7000.0, -700.0, -5700.0, 0.0},
	{_INFINITY, -3800.0, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, -6600.0, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{-1900.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{-6000.0, -3000.0, -5900.0, -7400.0, 0.0},
	{-1100.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{-3500.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}},

	{{{-2500.0, -2300.0, -2000.0, -5000.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, -2500.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, -3900.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, -3200.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{-2700.0, -700.0, -6000.0, -2500.0, 0.0},
	{_INFINITY, _INFINITY, -4900.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, -5700.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, -7400.0, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{-2400.0, -6000.0, -1100.0, -3900.0, 0.0},
	{_INFINITY, -3200.0, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, -3800.0, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, -6100.0, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{-5000.0, -700.0, -3600.0, -3200.0, 0.0},
	{-900.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{-3200.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{-2400.0, _INFINITY, _INFINITY, _INFINITY, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}},

	{{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}}};

static double dangleEntropies3[5][5][5] = {
	{{-1.0, -1.0, -1.0, -1.1, -1.0},
	{-1.0, -1.0, -1.0, 14.2, -1.0},
	{-1.0, -1.0, -1.0, -13.1, -1.0},
	{-1.0, -1.0, -1.0, -12.6, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -16.5, -1.0, -1.0},
	{-1.0, -1.0, -7.4, -1.0, -1.0},
	{-1.0, -1.0, -10.4, -1.0, -1.0},
	{-1.0, -1.0, -15.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -3.9, -1.0, -1.0, -1.0},
	{-1.0, -0.1, -1.0, -1.0, -1.0},
	{-1.0, -11.2, -1.0, -1.0, -1.0},
	{-1.0, -13.1, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-0.8, -1.0, -1.0, -1.0, -1.0},
	{14.9, -1.0, -1.0, -1.0, -1.0},
	{-3.6, -1.0, -1.0, -1.0, -1.0},
	{10.4, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}};

static double dangleEnthalpies3[5][5][5] = {
	{{_INFINITY, _INFINITY, _INFINITY, -500.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, 4700.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, -4100.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, -3800.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, -5900.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, -2600.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, -3200.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, -5200.0, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, -2100.0, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, -200.0, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, -3900.0, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, -4400.0, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{-700.0, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{4400.0, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{-1600.0, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{2900.0, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}};

static double dangleEntropies5[5][5][5] = {
	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-7.6, -13.0, -15.0, -0.5, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-10.0, -11.9, -10.9, -13.8, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-17.1, -12.6, -14.0, -10.9, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{2.3, 3.3, -1.6, -20.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}},

	{{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0},
	{-1.0, -1.0, -1.0, -1.0, -1.0}}};

static double dangleEnthalpies5[5][5][5] = {
	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{-2900.0, -4100.0, -4200.0, -200.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{-3700.0, -4000.0, -3900.0, -4900.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{-6300.0, -4400.0, -5100.0, -4000.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{200.0, 600.0, -1100.0, -6900.0, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}},

	{{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY},
	{_INFINITY, _INFINITY, _INFINITY, _INFINITY, _INFINITY}}};

static double interiorLoopEntropies[30] = {
	-1.0, -1.0, -10.31, -11.6, -12.89, 
	-14.18, -14.83, -15.47, -15.79, -15.79, 
	-16.26, -16.76, -17.15, -17.41, -17.74, 
	-18.05, -18.34, -18.7, -18.96, -19.02, 
	-19.25, -19.48, -19.7, -19.9, -20.31, 
	-20.5, -20.68, -20.86, -21.03, -21.28};

static double interiorLoopEnthalpies[30] = {
	_INFINITY, _INFINITY, 0.0, 0.0, 0.0, 
	0.0, 0.0, 0.0, 0.0, 0.0, 
	0.0, 0.0, 0.0, 0.0, 0.0, 
	0.0, 0.0, 0.0, 0.0, 0.0, 
	0.0, 0.0, 0.0, 0.0, 0.0, 
	0.0, 0.0, 0.0, 0.0, 0.0};

static double bulgeLoopEntropies[30] = {
	-12.89, -9.35, -9.99, -10.31, -10.64, 
	-11.28, -11.92, -12.57, -13.21, -13.86, 
	-14.32, -14.5, -14.89, -15.47, -15.81, 
	-16.12, -16.41, -16.76, -17.02, -17.08, 
	-17.32, -17.55, -17.76, -17.97, -18.05, 
	-18.24, -18.42, -18.6, -18.77, -19.02};

static double bulgeLoopEnthalpies[30] = {
	0.0, 0.0, 0.0, 0.0, 0.0, 
	0.0, 0.0, 0.0, 0.0, 0.0, 
	0.0, 0.0, 0.0, 0.0, 0.0, 
	0.0, 0.0, 0.0, 0.0, 0.0, 
	0.0, 0.0, 0.0, 0.0, 0.0, 
	0.0, 0.0, 0.0, 0.0, 0.0};

static double hairpinLoopEntropies[30] = {
	-1.0, -1.0, -11.28, -11.28, -10.64, 
	-12.89, -13.54, -13.86, -14.5, -14.83, 
	-15.29, -16.12, -16.5, -16.44, -16.77, 
	-17.08, -17.38, -17.73, -17.99, -18.37, 
	-18.61, -18.84, -19.05, -19.26, -19.66, 
	-19.85, -20.04, -20.21, -20.38, -20.31};

static double hairpinLoopEnthalpies[30] = {
	_INFINITY, _INFINITY, 0.0, 0.0, 0.0, 
	0.0, 0.0, 0.0, 0.0, 0.0, 
	0.0, 0.0, 0.0, 0.0, 0.0, 
	0.0, 0.0, 0.0, 0.0, 0.0, 
	0.0, 0.0, 0.0, 0.0, 0.0, 
	0.0, 0.0, 0.0, 0.0, 0.0};

static int numTriloops = 16;
static int numTetraloops = 77;
static struct triloop defaultTriloopEntropies[] = {
	{{0,2,0,0,3}, 0},
	{{0,2,1,0,3}, 0},
	{{0,2,2,0,3}, 0},
	{{0,2,3,0,3}, 0},
	{{1,2,0,0,2}, 0},
	{{1,2,1,0,2}, 0},
	{{1,2,2,0,2}, 0},
	{{1,2,3,0,2}, 0},
	{{2,2,0,0,1}, 0},
	{{2,2,1,0,1}, 0},
	{{2,2,2,0,1}, 0},
	{{2,2,3,0,1}, 0},
	{{3,2,0,0,0}, 0},
	{{3,2,1,0,0}, 0},
	{{3,2,2,0,0}, 0},
	{{3,2,3,0,0}, 0}};

static struct triloop defaultTriloopEnthalpies[] = {
	{{0,2,0,0,3}, -1500},
	{{0,2,1,0,3}, -1500},
	{{0,2,2,0,3}, -1500},
	{{0,2,3,0,3}, -1500},
	{{1,2,0,0,2}, -2000},
	{{1,2,1,0,2}, -2000},
	{{1,2,2,0,2}, -2000},
	{{1,2,3,0,2}, -2000},
	{{2,2,0,0,1}, -2000},
	{{2,2,1,0,1}, -2000},
	{{2,2,2,0,1}, -2000},
	{{2,2,3,0,1}, -2000},
	{{3,2,0,0,0}, -1500},
	{{3,2,1,0,0}, -1500},
	{{3,2,2,0,0}, -1500},
	{{3,2,3,0,0}, -1500}};

static struct tetraloop defaultTetraloopEntropies[] = {
	{{0,0,0,0,0,3}, -650},
	{{0,0,0,0,1,3}, 1610},
	{{0,0,0,1,0,3}, 1610},
	{{0,1,3,3,2,3}, 4190},
	{{0,2,0,0,0,3}, 1610},
	{{0,2,0,2,0,3}, 1610},
	{{0,2,0,3,0,3}, 1610},
	{{0,2,1,0,0,3}, 1610},
	{{0,2,1,2,0,3}, 1610},
	{{0,2,1,3,3,3}, 1610},
	{{0,2,2,0,0,3}, 1610},
	{{0,2,2,2,0,3}, 1610},
	{{0,2,2,2,2,3}, 640},
	{{0,2,3,0,0,3}, 1610},
	{{0,2,3,2,0,3}, 1610},
	{{0,2,3,3,1,3}, 1610},
	{{0,3,3,1,2,3}, 1610},
	{{0,3,3,3,2,3}, 1610},
	{{0,3,3,3,3,3}, 1610},
	{{1,0,0,0,0,2}, -1290},
	{{1,0,0,0,1,2}, 0},
	{{1,0,0,1,0,2}, 0},
	{{1,0,0,1,1,2}, 0},
	{{1,1,3,3,2,2}, 2570},
	{{1,2,0,0,0,2}, 0},
	{{1,2,0,2,0,2}, 0},
	{{1,2,0,3,0,2}, 0},
	{{1,2,1,0,0,2}, 0},
	{{1,2,1,2,0,2}, 0},
	{{1,2,1,3,3,2}, 0},
	{{1,2,2,0,0,2}, 0},
	{{1,2,2,2,0,2}, 0},
	{{1,2,2,2,2,2}, -970},
	{{1,2,3,0,0,2}, 0},
	{{1,2,3,2,0,2}, 0},
	{{1,2,3,3,1,2}, 0},
	{{1,3,3,1,2,2}, 0},
	{{1,3,3,3,2,2}, 0},
	{{1,3,3,3,3,2}, 0},
	{{2,0,0,0,0,1}, -3230},
	{{2,0,0,0,1,1}, 0},
	{{2,0,0,1,0,1}, 0},
	{{2,1,3,3,2,1}, 2570},
	{{2,2,0,0,0,1}, 0},
	{{2,2,0,2,0,1}, 0},
	{{2,2,0,3,0,1}, 0},
	{{2,2,1,0,0,1}, 0},
	{{2,2,1,2,0,1}, 0},
	{{2,2,1,3,3,1}, 0},
	{{2,2,2,0,0,1}, 0},
	{{2,2,2,2,0,1}, 0},
	{{2,2,2,2,2,1}, -970},
	{{2,2,3,0,0,1}, 0},
	{{2,2,3,2,0,1}, 0},
	{{2,2,3,3,1,1}, 0},
	{{2,3,3,1,2,1}, 0},
	{{2,3,3,3,2,1}, 0},
	{{2,3,3,3,3,1}, 0},
	{{3,0,0,0,0,0}, 320},
	{{3,0,0,0,1,0}, 1610},
	{{3,0,0,1,0,0}, 1610},
	{{3,1,3,3,2,0}, 4190},
	{{3,2,0,0,0,0}, 1610},
	{{3,2,0,2,0,0}, 1610},
	{{3,2,0,3,0,0}, 1610},
	{{3,2,1,0,0,0}, 1610},
	{{3,2,1,2,0,0}, 1610},
	{{3,2,1,3,3,0}, 1610},
	{{3,2,2,0,0,0}, 1610},
	{{3,2,2,2,0,0}, 1610},
	{{3,2,2,2,2,0}, 640},
	{{3,2,3,0,0,0}, 1610},
	{{3,2,3,2,0,0}, 1610},
	{{3,2,3,3,1,0}, 1610},
	{{3,3,3,1,2,0}, 1610},
	{{3,3,3,3,2,0}, 1610},
	{{3,3,3,3,3,0}, 1610}};

static struct tetraloop defaultTetraloopEnthalpies[] = {
	{{0,0,0,0,0,3}, 500},
	{{0,0,0,0,1,3}, 700},
	{{0,0,0,1,0,3}, 1000},
	{{0,1,3,3,2,3}, 0},
	{{0,2,0,0,0,3}, -1100},
	{{0,2,0,2,0,3}, -1100},
	{{0,2,0,3,0,3}, -1500},
	{{0,2,1,0,0,3}, -1600},
	{{0,2,1,2,0,3}, -1100},
	{{0,2,1,3,3,3}, 200},
	{{0,2,2,0,0,3}, -1100},
	{{0,2,2,2,0,3}, -1100},
	{{0,2,2,2,2,3}, 500},
	{{0,2,3,0,0,3}, -1600},
	{{0,2,3,2,0,3}, -1100},
	{{0,2,3,3,1,3}, 800},
	{{0,3,3,1,2,3}, -200},
	{{0,3,3,3,2,3}, 0},
	{{0,3,3,3,3,3}, -500},
	{{1,0,0,0,0,2}, 500},
	{{1,0,0,0,1,2}, 700},
	{{1,0,0,1,0,2}, 1000},
	{{1,0,0,1,1,2}, 0},
	{{1,1,3,3,2,2}, 0},
	{{1,2,0,0,0,2}, -1100},
	{{1,2,0,2,0,2}, -1100},
	{{1,2,0,3,0,2}, -1500},
	{{1,2,1,0,0,2}, -1600},
	{{1,2,1,2,0,2}, -1100},
	{{1,2,1,3,3,2}, 200},
	{{1,2,2,0,0,2}, -1100},
	{{1,2,2,2,0,2}, -1000},
	{{1,2,2,2,2,2}, 500},
	{{1,2,3,0,0,2}, -1600},
	{{1,2,3,2,0,2}, -1100},
	{{1,2,3,3,1,2}, 800},
	{{1,3,3,1,2,2}, -200},
	{{1,3,3,3,2,2}, 0},
	{{1,3,3,3,3,2}, -500},
	{{2,0,0,0,0,1}, 500},
	{{2,0,0,0,1,1}, 700},
	{{2,0,0,1,0,1}, 1000},
	{{2,1,3,3,2,1}, 0},
	{{2,2,0,0,0,1}, -1100},
	{{2,2,0,2,0,1}, -1100},
	{{2,2,0,3,0,1}, -1600},
	{{2,2,1,0,0,1}, -1600},
	{{2,2,1,2,0,1}, -1100},
	{{2,2,1,3,3,1}, 200},
	{{2,2,2,0,0,1}, -1100},
	{{2,2,2,2,0,1}, -1100},
	{{2,2,2,2,2,1}, 500},
	{{2,2,3,0,0,1}, -1600},
	{{2,2,3,2,0,1}, -1100},
	{{2,2,3,3,1,1}, 800},
	{{2,3,3,1,2,1}, -200},
	{{2,3,3,3,2,1}, 0},
	{{2,3,3,3,3,1}, -500},
	{{3,0,0,0,0,0}, 500},
	{{3,0,0,0,1,0}, 700},
	{{3,0,0,1,0,0}, 1000},
	{{3,1,3,3,2,0}, 0},
	{{3,2,0,0,0,0}, -1100},
	{{3,2,0,2,0,0}, -1100},
	{{3,2,0,3,0,0}, -1600},
	{{3,2,1,0,0,0}, -1600},
	{{3,2,1,2,0,0}, -1100},
	{{3,2,1,3,3,0}, 200},
	{{3,2,2,0,0,0}, -1100},
	{{3,2,2,2,0,0}, -1100},
	{{3,2,2,2,2,0}, 500},
	{{3,2,3,0,0,0}, -1600},
	{{3,2,3,2,0,0}, -1100},
	{{3,2,3,3,1,0}, 800},
	{{3,3,3,1,2,0}, -200},
	{{3,3,3,3,2,0}, 0},
	{{3,3,3,3,3,0}, -500}};

static struct triloop *triloopEntropies = defaultTriloopEntropies;
static struct triloop *triloopEnthalpies = defaultTriloopEnthalpies;
static struct tetraloop *tetraloopEntropies = defaultTetraloopEntropies;
static struct tetraloop *tetraloopEnthalpies = defaultTetraloopEnthalpies;
