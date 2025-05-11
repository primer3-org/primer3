import numpy as np
import re

def cleanString(data):
    # Remove comments and empty lines
    data = re.sub("#[^\n]*\n", "", data)
    data = re.sub("\n\t+", "\n", data)
    data = re.sub("\n +", "\n", data)
    data = re.sub("\n\n+", "\n", data)
    data = re.sub("^\n", "", data)
    data = re.sub("\n$", "", data)
    return data.split("\n")

def removeFristColumn(data):
    for i in range(0, len(data)):
        data[i] = re.sub("^[ACTGN_ \t]*", "", data[i])
    return data

outstr = "#include <math.h>\n"
outstr += "#include <stdio.h>\n"
outstr += "#include \"thal.h\"\n\n"

outstr += """
# ifdef INTEGER
const double _INFINITY = 999999.0;
# else
# ifdef INFINITY
const double _INFINITY = INFINITY;
# else
const double _INFINITY = 1.0 / 0.0;
# endif
# endif

"""

infile = open("./primer3_config/stack.ds", "r")
stack_ds_values = cleanString(infile.read())
stack_ds_values = removeFristColumn(stack_ds_values)
infile.close()
infile = open("./primer3_config/stack.dh", "r")
stack_dh_values = cleanString(infile.read())
stack_dh_values = removeFristColumn(stack_dh_values)
infile.close()
infile = open("./primer3_config/tstack_tm_inf.ds", "r")
tstack_ds_values = cleanString(infile.read())
tstack_ds_values = removeFristColumn(tstack_ds_values)
infile.close()
infile = open("./primer3_config/tstack.dh", "r")
tstack_dh_values = cleanString(infile.read())
tstack_dh_values = removeFristColumn(tstack_dh_values)
infile.close()
infile = open("./primer3_config/tstack2.ds", "r")
tstack2_ds_values = cleanString(infile.read())
tstack2_ds_values = removeFristColumn(tstack2_ds_values)
infile.close()
infile = open("./primer3_config/tstack2.dh", "r")
tstack2_dh_values = cleanString(infile.read())
tstack2_dh_values = removeFristColumn(tstack2_dh_values)
infile.close()
infile = open("./primer3_config/stackmm.ds", "r")
stackmm_ds_values = cleanString(infile.read())
stackmm_ds_values = removeFristColumn(stackmm_ds_values)
infile.close()
infile = open("./primer3_config/stackmm.dh", "r")
stackmm_dh_values = cleanString(infile.read())
stackmm_dh_values = removeFristColumn(stackmm_dh_values)
infile.close()
infile = open("./primer3_config/loops.dh", "r")
loops_dh_values = cleanString(infile.read())
infile.close()
infile = open("./primer3_config/loops.ds", "r")
loops_ds_values = cleanString(infile.read())
infile.close()
infile = open("./primer3_config/dangle.dh", "r")
dangle_dh_values = cleanString(infile.read())
dangle_dh_values = removeFristColumn(dangle_dh_values)
infile.close()
infile = open("./primer3_config/dangle.ds", "r")
dangle_ds_values = cleanString(infile.read())
dangle_ds_values = removeFristColumn(dangle_ds_values)
infile.close()
infile = open("./primer3_config/triloop.ds", "r")
triloop_ds_values = cleanString(infile.read())
infile.close()
infile = open("./primer3_config/triloop.dh", "r")
triloop_dh_values = cleanString(infile.read())
infile.close()
infile = open("./primer3_config/tetraloop.ds", "r")
tetraloop_ds_values = cleanString(infile.read())
infile.close()
infile = open("./primer3_config/tetraloop.dh", "r")
tetraloop_dh_values = cleanString(infile.read())
infile.close()


stack_ds = np.zeros((5,5,5,5), dtype=float)
stack_dh = np.zeros((5,5,5,5), dtype=float)
stackmm_ds = np.zeros((5,5,5,5), dtype=float)
stackmm_dh = np.zeros((5,5,5,5), dtype=float)
tstack_ds = np.zeros((5,5,5,5), dtype=float)
tstack_dh = np.zeros((5,5,5,5), dtype=float)
tstack2_ds = np.zeros((5,5,5,5), dtype=float)
tstack2_dh = np.zeros((5,5,5,5), dtype=float)
dangle3_ds = np.zeros((5,5,5), dtype=float)
dangle3_dh = np.zeros((5,5,5), dtype=float)
dangle5_ds = np.zeros((5,5,5), dtype=float)
dangle5_dh = np.zeros((5,5,5), dtype=float)


idx = 0
for i in range(5):
    for j in range(5):
        for k in range(5):
            for l in range(5):
                if 4 in [i, k]:
                    tstack_ds[i][j][k][l] = -1.0
                    tstack_dh[i][j][k][l] = float('inf')
                    tstack2_ds[i][j][k][l] = -1.0
                    tstack2_dh[i][j][k][l] = float('inf')
                elif 4 in [j, l]:
                    tstack_ds[i][j][k][l] = 0.00000000001
                    tstack_dh[i][j][k][l] = 0.0
                    tstack2_ds[i][j][k][l] = 0.00000000001
                    tstack2_dh[i][j][k][l] = 0.0
                else:
                    tstack_ds[i][j][k][l] = float(tstack_ds_values[idx])
                    tstack_dh[i][j][k][l] = float(tstack_dh_values[idx])
                    tstack2_ds[i][j][k][l] = float(tstack2_ds_values[idx])
                    tstack2_dh[i][j][k][l] = float(tstack2_dh_values[idx])
                    if float('inf') in [tstack_ds[i][j][k][l], tstack_dh[i][j][k][l]]:
                        tstack_ds[i][j][k][l] = -1.0
                        tstack_dh[i][j][k][l] = float('inf')
                    if float('inf') in [tstack2_ds[i][j][k][l], tstack2_dh[i][j][k][l]]:
                        tstack2_ds[i][j][k][l] = -1.0
                        tstack2_dh[i][j][k][l] = float('inf')

                if 4 in [i, j, k, l]:
                    stack_ds[i][j][k][l] = -1.0
                    stack_dh[i][j][k][l] = float('inf')
                    stackmm_ds[i][j][k][l] = -1.0
                    stackmm_dh[i][j][k][l] = float('inf')
                else:
                    stack_ds[i][j][k][l] = float(stack_ds_values[idx])
                    stack_dh[i][j][k][l] = float(stack_dh_values[idx])
                    stackmm_ds[i][j][k][l] = float(stackmm_ds_values[idx])
                    stackmm_dh[i][j][k][l] = float(stackmm_dh_values[idx])
                    idx += 1
                    if float('inf') in [stack_ds[i][j][k][l], stack_dh[i][j][k][l]]:
                        stack_ds[i][j][k][l] = -1.0
                        stack_dh[i][j][k][l] = float('inf')
                    if float('inf') in [stackmm_ds[i][j][k][l], stackmm_dh[i][j][k][l]]:
                        stackmm_ds[i][j][k][l] = -1.0
                        stackmm_dh[i][j][k][l] = float('inf')

idx = 0
for i in range(5):
    for j in range(5):
        for k in range(5):
            if 4 in [i, j, k]:
                dangle3_ds[i][k][j] = -1.0
                dangle3_dh[i][k][j] = float('inf')
                dangle5_ds[i][j][k] = -1.0
                dangle5_dh[i][j][k] = float('inf')
            else:
                dangle3_ds[i][k][j] = float(dangle_ds_values[idx])
                dangle3_dh[i][k][j] = float(dangle_dh_values[idx])
                dangle5_ds[i][j][k] = float(dangle_ds_values[idx+64])
                dangle5_dh[i][j][k] = float(dangle_dh_values[idx+64])
                idx += 1
                if float('inf') in [dangle3_ds[i][k][j], dangle3_dh[i][k][j]]:
                    dangle3_ds[i][k][j] = -1.0
                    dangle3_dh[i][k][j] = float('inf')
                if float('inf') in [dangle5_ds[i][j][k], dangle5_dh[i][j][k]]:
                    dangle5_ds[i][j][k] = -1.0
                    dangle5_dh[i][j][k] = float('inf')




stack_ds = stack_ds.flatten()
stack_dh = stack_dh.flatten()
tstack_ds = tstack_ds.flatten()
tstack_dh = tstack_dh.flatten()
tstack2_ds = tstack2_ds.flatten()
tstack2_dh = tstack2_dh.flatten()
stackmm_ds = stackmm_ds.flatten()
stackmm_dh = stackmm_dh.flatten()
dangle3_ds = dangle3_ds.flatten()
dangle3_dh = dangle3_dh.flatten()
dangle5_ds = dangle5_ds.flatten()
dangle5_dh = dangle5_dh.flatten()



matricies_5_5_5_5 = []
matricies_5_5_5_5.append(["static double stackEntropies[5][5][5][5]", stack_ds])
matricies_5_5_5_5.append(["static double stackEnthalpies[5][5][5][5]", stack_dh])
matricies_5_5_5_5.append(["static double stackint2Entropies[5][5][5][5]", stackmm_ds])
matricies_5_5_5_5.append(["static double stackint2Enthalpies[5][5][5][5]", stackmm_dh])
matricies_5_5_5_5.append(["static double tstackEntropies[5][5][5][5]", tstack_ds])
matricies_5_5_5_5.append(["static double tstackEnthalpies[5][5][5][5]", tstack_dh])
matricies_5_5_5_5.append(["static double tstack2Entropies[5][5][5][5]", tstack2_ds])
matricies_5_5_5_5.append(["static double tstack2Enthalpies[5][5][5][5]", tstack2_dh])

matricies_5_5_5 = []
matricies_5_5_5.append(["static double dangleEntropies3[5][5][5]", dangle3_ds])
matricies_5_5_5.append(["static double dangleEnthalpies3[5][5][5]", dangle3_dh])
matricies_5_5_5.append(["static double dangleEntropies5[5][5][5]", dangle5_ds])
matricies_5_5_5.append(["static double dangleEnthalpies5[5][5][5]", dangle5_dh])

outstr += "static double atpS[5][5] = {"
for i in range(5):
    outstr += "\n\t{"
    for j in range(5):
        if ((i == 0) and (j == 3)) or ((i == 3) and (j == 0)):
            outstr += "6.9"
        else:
            outstr += "0.00000000001"
        if j != 4:
            outstr += ", "
    outstr += "}"
    if i != 4:
        outstr += ","
outstr += "};\n\n"

outstr += "static double atpH[5][5] = {"
for i in range(5):
    outstr += "\n\t{"
    for j in range(5):
        if ((i == 0) and (j == 3)) or ((i == 3) and (j == 0)):
            outstr += "2200.0"
        else:
            outstr += "0.0"
        if j != 4:
            outstr += ", "
    outstr += "}"
    if i != 4:
        outstr += ","
outstr += "};\n\n"

for mat in matricies_5_5_5_5:
    outstr += mat[0] + " = {\n\t"
    idx = 0
    for i in range(5):
        outstr += "{"
        for j in range(5):
            outstr += "{"
            for k in range(5):
                outstr += "{"
                for l in range(5):
                    if mat[1][idx] == float('inf'):
                        outstr += "_INFINITY"
                    else:
                        outstr += str(mat[1][idx])
                    idx += 1
                    if l != 4:
                        outstr += ", "
                outstr += "}"
                if k != 4:
                    outstr += ",\n\t"
            outstr += "}"
            if j != 4:
                outstr += ",\n\n\t"
        outstr += "}"
        if i != 4:
            outstr += ",\n\n\t"
    outstr += "};\n\n"

for mat in matricies_5_5_5:
    outstr += mat[0] + " = {\n\t"
    idx = 0
    for i in range(5):
        outstr += "{"
        for j in range(5):
            outstr += "{"
            for k in range(5):
                if mat[1][idx] == float('inf'):
                    outstr += "_INFINITY"
                else:
                    outstr += str(mat[1][idx])
                idx += 1
                if k != 4:
                    outstr += ", "
            outstr += "}"
            if j != 4:
                outstr += ",\n\t"
        outstr += "}"
        if i != 4:
            outstr += ",\n\n\t"
    outstr += "};\n\n"

bulge_loop_ds = np.zeros(30, dtype=float)
interior_loop_ds = np.zeros(30, dtype=float)
hairpin_loop_ds = np.zeros(30, dtype=float)
bulge_loop_dh = np.zeros(30, dtype=float)
interior_loop_dh = np.zeros(30, dtype=float)
hairpin_loop_dh = np.zeros(30, dtype=float)

for line in loops_dh_values:
    values = line.split("\t")
    idx = int(values[0]) -1
    interior_loop_dh[idx] = float(values[1])
    bulge_loop_dh[idx] = float(values[2])
    hairpin_loop_dh[idx] = float(values[3])

for line in loops_ds_values:
    values = line.split("\t")
    idx = int(values[0]) -1
    interior_loop_ds[idx] = float(values[1])
    bulge_loop_ds[idx] = float(values[2])
    hairpin_loop_ds[idx] = float(values[3])

loops_30 = []
loops_30.append(["static double interiorLoopEntropies[30]", interior_loop_ds])
loops_30.append(["static double interiorLoopEnthalpies[30]", interior_loop_dh])
loops_30.append(["static double bulgeLoopEntropies[30]", bulge_loop_ds])
loops_30.append(["static double bulgeLoopEnthalpies[30]", bulge_loop_dh])
loops_30.append(["static double hairpinLoopEntropies[30]", hairpin_loop_ds])
loops_30.append(["static double hairpinLoopEnthalpies[30]", hairpin_loop_dh])

for loop in loops_30:
    outstr += loop[0] + " = {"
    for i in range(len(loop[1])):
        if i % 5 == 0:
            outstr += "\n\t"
        if loop[1][i] == float('inf'):
            outstr += "_INFINITY"
        else:
            outstr += str(loop[1][i])
        if i != len(loop[1]) - 1:
            outstr += ", "
    outstr += "};\n\n"

outstr += f"static int numTriloops = {len(triloop_dh_values)};\n"
outstr += f"static int numTetraloops = {len(tetraloop_dh_values)};\n"

bases = ['A', 'C', 'G', 'T']

outstr += "static struct triloop defaultTriloopEntropies[] = {\n"
for i in range(len(triloop_ds_values)):
    line = triloop_ds_values[i].split("\t")
    outstr += "\t{{"
    for j in range(len(line[0])):
        outstr += f"{bases.index(line[0][j])}"
        if j != len(line[0]) - 1:
            outstr += ","
    outstr += "}, " + line[1] + "}"
    if i != len(triloop_ds_values) - 1:
        outstr += ","
    else:
        outstr += "};\n"
    outstr += "\n"

outstr += "static struct triloop defaultTriloopEnthalpies[] = {\n"
for i in range(len(triloop_dh_values)):
    line = triloop_dh_values[i].split("\t")
    outstr += "\t{{"
    for j in range(len(line[0])):
        outstr += f"{bases.index(line[0][j])}"
        if j != len(line[0]) - 1:
            outstr += ","
    outstr += "}, " + line[1] + "}"
    if i != len(triloop_dh_values) - 1:
        outstr += ","
    else:
        outstr += "};\n"
    outstr += "\n"

outstr += "static struct tetraloop defaultTetraloopEntropies[] = {\n"
for i in range(len(tetraloop_ds_values)):
    line = tetraloop_ds_values[i].split("\t")
    outstr += "\t{{"
    for j in range(len(line[0])):
        outstr += f"{bases.index(line[0][j])}"
        if j != len(line[0]) - 1:
            outstr += ","
    outstr += "}, " + line[1] + "}"
    if i != len(tetraloop_ds_values) - 1:
        outstr += ","
    else:
        outstr += "};\n"
    outstr += "\n"

outstr += "static struct tetraloop defaultTetraloopEnthalpies[] = {\n"
for i in range(len(tetraloop_dh_values)):
    line = tetraloop_dh_values[i].split("\t")
    outstr += "\t{{"
    for j in range(len(line[0])):
        outstr += f"{bases.index(line[0][j])}"
        if j != len(line[0]) - 1:
            outstr += ","
    outstr += "}, " + line[1] + "}"
    if i != len(tetraloop_dh_values) - 1:
        outstr += ","
    else:
        outstr += "};\n"
    outstr += "\n"

outstr += "static struct triloop *triloopEntropies = defaultTriloopEntropies;\n"
outstr += "static struct triloop *triloopEnthalpies = defaultTriloopEnthalpies;\n"
outstr += "static struct tetraloop *tetraloopEntropies = defaultTetraloopEntropies;\n"
outstr += "static struct tetraloop *tetraloopEnthalpies = defaultTetraloopEnthalpies;\n"

comment_str = """/*
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
*/\n\n"""

outfile = open("thal_default_params.h", "w+")
outfile.write(comment_str)
outfile.write(outstr)
outfile.close()