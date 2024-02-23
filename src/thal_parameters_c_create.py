import numpy as np


"""TODO:
static struct triloop* triloopEntropies = NULL; /* ther penalties for given triloop seq-s */
static struct triloop* triloopEnthalpies = NULL; /* ther penalties for given triloop seq-s */
static struct tetraloop* tetraloopEntropies = NULL; /* ther penalties for given tetraloop seq-s */
static struct tetraloop* tetraloopEnthalpies = NULL; /* ther penalties for given tetraloop seq-s */
"""

outstr = "#include <math.h>\n"
outstr += "#include <stdio.h>\n\n"

infile = open("./primer3_config/stack.ds", "r")
stack_ds_values = infile.read().splitlines()
infile.close()
infile = open("./primer3_config/stack.dh", "r")
stack_dh_values = infile.read().splitlines()
infile.close()
infile = open("./primer3_config/tstack_tm_inf.ds", "r")
tstack_ds_values = infile.read().splitlines()
infile.close()
infile = open("./primer3_config/tstack.dh", "r")
tstack_dh_values = infile.read().splitlines()
infile.close()
infile = open("./primer3_config/tstack2.ds", "r")
tstack2_ds_values = infile.read().splitlines()
infile.close()
infile = open("./primer3_config/tstack2.dh", "r")
tstack2_dh_values = infile.read().splitlines()
infile.close()
infile = open("./primer3_config/stackmm.ds", "r")
stackmm_ds_values = infile.read().splitlines()
infile.close()
infile = open("./primer3_config/stackmm.dh", "r")
stackmm_dh_values = infile.read().splitlines()
infile.close()
infile = open("./primer3_config/loops.dh", "r")
loops_dh_values = infile.read().splitlines()
infile.close()
infile = open("./primer3_config/loops.ds", "r")
loops_ds_values = infile.read().splitlines()
infile.close()
infile = open("./primer3_config/dangle.dh", "r")
dangle_dh_values = infile.read().splitlines()
infile.close()
infile = open("./primer3_config/dangle.ds", "r")
dangle_ds_values = infile.read().splitlines()
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
                if 4 in [i, j, k, l]:
                    stack_ds[i][j][k][l] = -1.0
                    stack_dh[i][j][k][l] = float('inf')
                    tstack_ds[i][j][k][l] = -1.0
                    tstack_dh[i][j][k][l] = float('inf')
                    tstack2_ds[i][j][k][l] = -1.0
                    tstack2_dh[i][j][k][l] = float('inf')
                    stackmm_ds[i][j][k][l] = -1.0
                    stackmm_dh[i][j][k][l] = float('inf')
                else:
                    stack_ds[i][j][k][l] = float(stack_ds_values[idx])
                    stack_dh[i][j][k][l] = float(stack_dh_values[idx])
                    tstack_ds[i][j][k][l] = float(tstack_ds_values[idx])
                    tstack_dh[i][j][k][l] = float(tstack_dh_values[idx])
                    tstack2_ds[i][j][k][l] = float(tstack2_ds_values[idx])
                    tstack2_dh[i][j][k][l] = float(tstack2_dh_values[idx])
                    stackmm_ds[i][j][k][l] = float(stackmm_ds_values[idx])
                    stackmm_dh[i][j][k][l] = float(stackmm_dh_values[idx])
                    idx += 1
                    if float('inf') in [stack_ds[i][j][k][l], stack_dh[i][j][k][l]]:
                        stack_ds[i][j][k][l] = -1.0
                        stack_dh[i][j][k][l] = float('inf')
                    if float('inf') in [stackmm_ds[i][j][k][l], stackmm_dh[i][j][k][l]]:
                        stackmm_ds[i][j][k][l] = -1.0
                        stackmm_dh[i][j][k][l] = float('inf')
                    if float('inf') in [tstack_ds[i][j][k][l], tstack_dh[i][j][k][l]]:
                        tstack_ds[i][j][k][l] = -1.0
                        tstack_dh[i][j][k][l] = float('inf')
                    if float('inf') in [tstack2_ds[i][j][k][l], tstack2_dh[i][j][k][l]]:
                        tstack2_ds[i][j][k][l] = -1.0
                        tstack2_dh[i][j][k][l] = float('inf')

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



outstr += "double stackEntropies[5][5][5][5] = {"
for i in range(len(stack_ds)):
    if i % 5 == 0:
        outstr += "\n\t"
    if i % 25 == 0:
        outstr += "\n\t"
    if stack_ds[i] == float('inf'):
        outstr += "INFINITY"
    else:
        outstr += str(stack_ds[i])
    if i != len(stack_ds) - 1:
        outstr += ", "
outstr += "};\n\n"

outstr += "double stackEnthalpies[5][5][5][5] = {"
for i in range(len(stack_dh)):
    if i % 5 == 0:
        outstr += "\n\t"
    if i % 25 == 0:
        outstr += "\n\t"
    if stack_dh[i] == float('inf'):
        outstr += "INFINITY"
    else:
        outstr += str(stack_dh[i])
    if i != len(stack_dh) - 1:
        outstr += ", "
outstr += "};\n\n"

outstr += "double tstackEntropies[5][5][5][5] = {"
for i in range(len(tstack_ds)):
    if i % 5 == 0:
        outstr += "\n\t"
    if i % 25 == 0:
        outstr += "\n\t"
    if tstack_ds[i] == float('inf'):
        outstr += "INFINITY"
    else:
        outstr += str(tstack_ds[i])
    if i != len(tstack_ds) - 1:
        outstr += ", "
outstr += "};\n\n"

outstr += "double tstackEnthalpies[5][5][5][5] = {"
for i in range(len(tstack_dh)):
    if i % 5 == 0:
        outstr += "\n\t"
    if i % 25 == 0:
        outstr += "\n\t"
    if tstack_dh[i] == float('inf'):
        outstr += "INFINITY"
    else:
        outstr += str(tstack_dh[i])
    if i != len(tstack_dh) - 1:
        outstr += ", "
outstr += "};\n\n"

outstr += "double tstack2Entropies[5][5][5][5] = {"
for i in range(len(tstack2_ds)):
    if i % 5 == 0:
        outstr += "\n\t"
    if i % 25 == 0:
        outstr += "\n\t"
    if tstack2_ds[i] == float('inf'):
        outstr += "INFINITY"
    else:
        outstr += str(tstack2_ds[i])
    if i != len(tstack2_ds) - 1:
        outstr += ", "
outstr += "};\n\n"

outstr += "double tstack2Enthalpies[5][5][5][5] = {"
for i in range(len(tstack2_dh)):
    if i % 5 == 0:
        outstr += "\n\t"
    if i % 25 == 0:
        outstr += "\n\t"
    if tstack2_dh[i] == float('inf'):
        outstr += "INFINITY"
    else:
        outstr += str(tstack2_dh[i])
    if i != len(tstack2_dh) - 1:
        outstr += ", "
outstr += "};\n\n"

outstr += "double stackint2Entropies[5][5][5][5] = {"
for i in range(len(stackmm_ds)):
    if i % 5 == 0:
        outstr += "\n\t"
    if i % 25 == 0:
        outstr += "\n\t"
    if stackmm_ds[i] == float('inf'):
        outstr += "INFINITY"
    else:
        outstr += str(stackmm_ds[i])
    if i != len(stackmm_ds) - 1:
        outstr += ", "
outstr += "};\n\n"

outstr += "double stackint2Enthalpies[5][5][5][5] = {"
for i in range(len(stackmm_dh)):
    if i % 5 == 0:
        outstr += "\n\t"
    if i % 25 == 0:
        outstr += "\n\t"
    if stackmm_dh[i] == float('inf'):
        outstr += "INFINITY"
    else:
        outstr += str(stackmm_dh[i])
    if i != len(stackmm_dh) - 1:
        outstr += ", "
outstr += "};\n\n"

outstr += "double dangle3Entropies[5][5][5] = {"
for i in range(len(dangle3_ds)):
    if i % 5 == 0:
        outstr += "\n\t"
    if i % 25 == 0:
        outstr += "\n\t"
    if dangle3_ds[i] == float('inf'):
        outstr += "INFINITY"
    else:
        outstr += str(dangle3_ds[i])
    if i != len(dangle3_ds) - 1:
        outstr += ", "
outstr += "};\n\n"

outstr += "double dangle3Enthalpies[5][5][5] = {"
for i in range(len(dangle3_dh)):
    if i % 5 == 0:
        outstr += "\n\t"
    if i % 25 == 0:
        outstr += "\n\t"
    if dangle3_dh[i] == float('inf'):
        outstr += "INFINITY"
    else:
        outstr += str(dangle3_dh[i])
    if i != len(dangle3_dh) - 1:
        outstr += ", "
outstr += "};\n\n"

outstr += "double dangle5Entropies[5][5][5] = {"
for i in range(len(dangle5_ds)):
    if i % 5 == 0:
        outstr += "\n\t"
    if i % 25 == 0:
        outstr += "\n\t"
    if dangle5_ds[i] == float('inf'):
        outstr += "INFINITY"
    else:
        outstr += str(dangle5_ds[i])
    if i != len(dangle5_ds) - 1:
        outstr += ", "
outstr += "};\n\n"

outstr += "double dangle5Enthalpies[5][5][5] = {"
for i in range(len(dangle5_dh)):
    if i % 5 == 0:
        outstr += "\n\t"
    if i % 25 == 0:
        outstr += "\n\t"
    if dangle5_dh[i] == float('inf'):
        outstr += "INFINITY"
    else:
        outstr += str(dangle5_dh[i])
    if i != len(dangle5_dh) - 1:
        outstr += ", "
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

outstr += "double interiorLoopEntropies[30] = {"
for i in range(len(interior_loop_ds)):
    if i % 5 == 0:
        outstr += "\n\t"
    if interior_loop_ds[i] == float('inf'):
        outstr += "INFINITY"
    else:
        outstr += str(interior_loop_ds[i])
    if i != len(interior_loop_ds) - 1:
        outstr += ", "
outstr += "};\n\n"

outstr += "double bulgeLoopEntropies[30] = {"
for i in range(len(bulge_loop_ds)):
    if i % 5 == 0:
        outstr += "\n\t"
    if bulge_loop_ds[i] == float('inf'):
        outstr += "INFINITY"
    else:
        outstr += str(bulge_loop_ds[i])
    if i != len(bulge_loop_ds) - 1:
        outstr += ", "
outstr += "};\n\n"

outstr += "double hairpinLoopEntropies[30] = {"
for i in range(len(hairpin_loop_ds)):
    if i % 5 == 0:
        outstr += "\n\t"
    if hairpin_loop_ds[i] == float('inf'):
        outstr += "INFINITY"
    else:
        outstr += str(hairpin_loop_ds[i])
    if i != len(hairpin_loop_ds) - 1:
        outstr += ", "
outstr += "};\n\n"

outstr += "double interiorLoopEnthalpies[30] = {"
for i in range(len(interior_loop_dh)):
    if i % 5 == 0:
        outstr += "\n\t"
    if interior_loop_dh[i] == float('inf'):
        outstr += "INFINITY"
    else:
        outstr += str(interior_loop_dh[i])
    if i != len(interior_loop_dh) - 1:
        outstr += ", "
outstr += "};\n\n"

outstr += "double bulgeLoopEnthalpies[30] = {"
for i in range(len(bulge_loop_dh)):
    if i % 5 == 0:
        outstr += "\n\t"
    if bulge_loop_dh[i] == float('inf'):
        outstr += "INFINITY"
    else:
        outstr += str(bulge_loop_dh[i])
    if i != len(bulge_loop_dh) - 1:
        outstr += ", "
outstr += "};\n\n"

outstr += "double hairpinLoopEnthalpies[30] = {"
for i in range(len(hairpin_loop_dh)):
    if i % 5 == 0:
        outstr += "\n\t"
    if hairpin_loop_dh[i] == float('inf'):
        outstr += "INFINITY"
    else:
        outstr += str(hairpin_loop_dh[i])
    if i != len(hairpin_loop_dh) - 1:
        outstr += ", "
outstr += "};\n\n"


outfile = open("thal_parameters_new.c", "w+")
outfile.write(outstr)
outfile.close()