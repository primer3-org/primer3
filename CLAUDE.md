# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Primer3 picks PCR primers, hybridization probes ("internal oligos") and sequencing primers from a DNA template. The main program is `primer3_core`. The code is C plus one C++11 file (`src/libprimer3.cc`), built with plain Makefiles. Tests are Perl scripts that diff program output against checked-in golden files.

- `src/primer3_manual.htm` is the reference for every input and output tag. Output tags are parsed by other tools, so renaming one is a breaking change.
- The manual is maintained in a separate repository, and the copy here is only updated from there. Edits to `src/primer3_manual.htm` in this repository are proposals that the maintainer carries over to that repository. Still make them, so the proposed manual text sits next to the code change.
- User-visible changes get an entry at the top of `src/release_notes.txt`, under the unreleased version heading.

## Build

Everything builds from `src/`:

```sh
cd src
make                  # all executables and static libs (a few seconds)
make O_OPTS=-O0       # unoptimized, for gdb/valgrind
make clean_src        # cleans src only; `make clean` also deletes test temp files
```

- The build has no warnings under `-Wall`.
- `make` doesn't notice changed flags. Run `make clean_src` before rebuilding with different `O_OPTS` or `CC_OPTS`.
- The compiler writes the header dependencies of each object file into `src/*.d` (`DEPFLAGS = -MMD -MP` in `src/Makefile`), so `make` rebuilds every object that includes a changed header. Object rules list only their source file; don't add headers by hand. A program built from a `*_main.c` needs its own object rule, because compiling and linking in one step writes no usable `.d` file.
- Windows builds use `make TESTOPTS=--windows` (MinGW). That build leaves out the masker.

## Code style

In new or changed C/C++ code:
- Write every `if`, `else if` and `else` with braces, never as a one-liner such as `if (x) y = 1;`.
- Don't use the ternary operator `?:` to choose a value; write `if (...) { ... } else { ... }`.

Don't reformat untouched code. But when new code belongs under an existing braceless `if`, put it in that `if`: add the braces there rather than repeating the condition in a second block.

## Tests

Run the test scripts from `test/`. They find the binaries through `../src/`.

```sh
cd src && make test                                # build, then run every suite (~1.75 min)
cd test && perl run_tests.pl --only=thal,p3 --fast # some suites; --fast skips the slowest p3 tests
perl p3test.pl --onetest primer_task               # run one primer3_core test
perl thal_test.pl     # also: dpal_test.pl oligotm_test.pl masker_test.pl cmdline_test.pl a3test.pl
```

`make test` runs `test/run_tests.pl`:
- It runs every suite, also after one fails, and passes `--action` to `p3test.pl` and `a3test.pl`.
- At the end it prints the run time of each suite, every failed test and a `RESULT:` line, and exits 1 if anything failed. So read the summary instead of searching the output.
- Options: `--windows` (also skips the masker), `--valgrind`, `--valgrinda`, `--valgrindb`, `--fast` (only for `p3test.pl`) and `--only=SUITE,...` with the suites `cmdline oligotm dpal thal masker a3 p3`.
- `--valgrinda` and `--valgrindb` split the valgrind run as CI does: `--valgrinda` runs the small suites with `--valgrind` and `a3test.pl` and `p3test.pl` with `--valgrinda`; `--valgrindb` runs only `p3test.pl --valgrindb`.
- With `GITHUB_ACTIONS=true` it also writes `::error` annotations and a job summary. The workflows call each script in its own step; only the valgrind jobs use it, through `make valgrnda` and `make valgrndb`.
- `make test_separately` is the old chain of targets, which stops at the first suite that exits non-zero.

When you run a single script, search its output for `[FAILED]`:
- `p3test.pl` and `a3test.pl` exit 0 even when tests fail, unless you pass `--action`.
- The other scripts exit non-zero on failure.

The scripts run the binaries in `../src/`, so don't rebuild while a script is running. For experiments with other build flags, use a copy outside the repository, for example `git archive HEAD | tar -x -C <dir>`.

To run under valgrind, first build with `make O_OPTS=-O0`. Then, in `test/`, run `make valgrnda` and `make valgrndb` (how CI splits it; they call `run_tests.pl --valgrinda` and `--valgrindb`), or `make TESTOPTS=--valgrind` for everything. Each script checks only the valgrind logs it wrote itself: `*.vg` from `cmdline_test.pl`, `a3test.pl` and `p3test.pl`, and `*.valg` from the other scripts. `src/paranoid_tests.sh` runs the full -O2, -O3 and valgrind sequence.

For AddressSanitizer and UBSan, rebuild and then run the test scripts as usual:

```sh
cd src && make clean_src
make O_OPTS="-O1 -fno-omit-frame-pointer" \
     CC_OPTS="-g -fsanitize=address,undefined -fno-sanitize-recover=undefined" \
     LDFLAGS="-g -fsanitize=address,undefined" LDLIBS="-lm -fsanitize=address,undefined"
```

A sanitizer finding makes the program exit non-zero, which the scripts report as a failure. `p3test.pl` takes about 6 minutes in this build. Run `make clean_src && make` afterwards to get a normal build back.

### How p3test.pl works

For each test name `X` in `@TESTS`, the driver runs `primer3_core --strict_tags` on `X_input` (Boulder-IO), writes `X.tmp`, and compares it with `X_output`. The comparison is exact except for three things: the executable name, `primer3 release N.N.N` strings and path separators. A plain `diff` will therefore also show release-number lines that the driver ignores.

The test name selects a variant:
- `X_formatted` reuses `X_input` and adds `--format_output`.
- `*_settings` loads `../settings_files/X.txt` through `--p3_settings_file`. `*_load_set` loads `X.set`.
- `primer`, `primer1`, `primer1_th` and `th-w-other-tasks` also write primer list files to `X_list_tmp/`, which are compared with `X_list_last/`.

Flags differ from test to test:
- Most tests run with `--default_version=1` (the old defaults), but the CLI default is 2. The `%default_version2` hash lists the tests that use 2.
- To reproduce a test by hand, add `--verbose`, which prints the exact command.
- Without `--onetest`, the driver first runs the fatal-error cases. Each `primer_global_err/*.in` must exit non-zero, and its stdout and stderr are compared with `.out` and `.out2`.

To change expected output on purpose:
1. Inspect `diff X_output X.tmp`.
2. Copy the `.tmp` file over the golden file.
3. Commit the golden file together with the code change.

To add a test, add the `_input` and `_output` files and append the name to `@TESTS`. CI runs on Linux, macOS and Windows, so printed numbers must match on all three. That is why `thal.c` is compiled with `-ffloat-store`.

### How thal_test.pl, dpal_test.pl and masker_test.pl work

These three drivers share one design:
- Each line of the input file (`thal_input`, `thal_long_input`, `dpal_input`, `masker_input`) is the argument list for one call of the tool.
- The stdout of all calls is concatenated into one `.tmp` file and compared with one golden file. stderr is never compared; all three drivers write it to `.tmp.err`.
- A call that exits non-zero adds a line `EXIT=<code>` at its place in the `.tmp` file. Lines that print nothing and produce an `EXIT` line are deliberate error cases, such as missing or invalid option values. On Windows the code is shifted (`>> 8`) so that it matches Linux.
- A new input line therefore changes the golden file at that position. Regenerate it by copying the `.tmp` file, and check that every new `EXIT` line is intended.

In `thal_test.pl`, the interactive runs (`-i`) use the same input and golden files as the normal runs. When comparing, it rounds the numbers on `Calculated thermodynamical parameters` lines to one decimal.

## Architecture

### primer3_core pipeline

1. `primer3_boulder_main.c` parses the CLI and creates `p3_global_settings`, with defaults from `pr_set_default_global_args_1` or `_2`. It applies an optional settings file. Then it loops over Boulder-IO records on stdin: `TAG=VALUE` lines, with a line holding only `=` ending each record.
2. `read_boulder.c` parses each record.
   - Global settings (`p3_global_settings`) carry over between records. Per-sequence data (`seq_args`) is recreated for each record.
   - Tags are matched by a chain of `COMPARE_INT`, `COMPARE_FLOAT`, `COMPARE_AND_MALLOC` and similar macros. `--strict_tags` makes unknown tags an error.
3. `choose_primers(pa, sa)` in `libprimer3.cc` is the library entry point. It returns a `p3retval` and runs these steps:
   1. `_adjust_seq_args`.
   2. `_pr_data_control`, which validates the input.
   3. Build the candidate lists, either with `make_detection_primer_lists` and `make_internal_oligo_list`, or through the task-specific paths for `pick_primer_list`, `pick_sequencing_primers` and `check_primers`.
   4. `choose_pair_or_triple`.

   Out-of-memory errors and thal sequence-length errors exit through `setjmp`/`longjmp`. See the warning comment in `choose_primers`.
4. Output:
   - `print_boulder.c` writes Boulder-IO output.
   - `format_output.c` writes the human-readable `--format_output` report.
   - A new output field usually needs changes in both.

Primer settings are in `pa->p_args` and internal-oligo (probe) settings are in `pa->o_args`. Both are `args_for_one_oligo_or_primer`. Default version 2 is version 1 with a few overrides: SantaLucia Tm and salt correction, thermodynamic oligo alignment, and the divalent and dNTP concentrations.

### State, coordinates and errors

- **Per-record and global state.** `seq_args` is created fresh for every record. `p3_global_settings` persists into all later records, and so do the file-scope globals in `read_boulder.c` such as `kmer_lists_path`. Don't store values derived from one record's input in the global settings; derive them where they are used. When you change code that reads or writes global settings, test it with a multi-record input.
- **Validation order.** `choose_primers()` calls `_adjust_seq_args()` before `_pr_data_control()`, so code in `_adjust_seq_args()` sees input that hasn't been validated yet. Settings that must be valid before masking, such as `PRIMER_MASK_TEMPLATE`, are therefore also checked in `read_boulder.c`.
- **Unvalidated input aborts.** Values that pass parsing but are out of range often end in a `PR_ASSERT` abort deep in the code. For new numeric input, check the range, and 0 where the value is a divisor or a step. `parse_double()` rejects NaN and infinity (it reports a syntax error and stores 0); values read with plain `strtod()`, as in the standalone CLIs, still need that check.
- **Coordinates.** Input positions count from `PRIMER_FIRST_BASE_INDEX` (default 0). `_adjust_seq_args()` makes them 0-based and relative to the included region. After that, intervals, junction lists, force positions, `start_codon_pos` and oligo positions index `trimmed_seq`, and `incl_s` is the 0-based start of the region in the full `sequence`. Output adds `incl_s + first_base_index` back.
- **Sequence copies in `seq_args`:**
  - `sequence`: the template as given.
  - `trimmed_seq`: the included region; `_pr_data_control()` upper-cases it and turns unknown bases into N.
  - `trimmed_orig_seq`: the included region in its original case, used for lowercase masking.
  - `upcased_seq` and `upcased_seq_r`: the whole template and its reverse complement, used for template mispriming.
  - `trimmed_masked_seq` and `trimmed_masked_seq_r`: the k-mer-masked region for each strand.
- **Oligo records (`primer_rec`).** `start` is the 5′ base: the leftmost base of a left primer or internal oligo, and the rightmost base of a right primer. `rev_comp` marks the printed orientation. For an internal oligo it is set when `PRIMER_INTERNAL_REV_COMP_G_RICH` prints the reverse complement, so checks that depend on orientation must use the printed sequence.
- **Errors.**
  - In `read_boulder.c`, a bad `SEQUENCE_*` tag is a per-record error: that record prints `PRIMER_ERROR` and the next record runs. A bad global tag is fatal: primer3_core prints the error and exits with 4.
  - In `p3retval`, `glob_err` is fatal, `per_sequence_err` ends only the current record, and `warnings` are printed as `PRIMER_WARNING`.
- **Exit codes.** All programs end with the codes in `src/p3_exit_codes.h`, and each code means the same in every program. Chapter 28 "EXIT STATUS CODES" of the manual documents them for users.

  | Code | Constant | Meaning |
  |---|---|---|
  | 0 | `P3_EXIT_OK` | success |
  | 1 | `P3_EXIT_USAGE` | invalid command line: unknown option, missing or invalid option value, wrong number of arguments, options that do not fit together |
  | 2 | `P3_EXIT_OUT_OF_MEMORY` | out of memory |
  | 3 | `P3_EXIT_EMPTY_INPUT` | no input record (primer3_core) |
  | 4 | `P3_EXIT_INPUT_ERROR` | input that cannot be read or used: input file, settings, library, parameter or k-mer files, a bad global tag, an invalid sequence, a calculation that fails for the input |
  | 5 | `P3_EXIT_OUTPUT_ERROR` | output that cannot be written: `--output`/`--error` file, failed `fflush` |

  - Call `exit()` only with these constants, never with a literal or a negative value. The golden files and `cmdline_test.pl` check the values, and callers depend on them.
  - Library code (`dpal.c`, `thal*.c`, `oligotm.c`, `masker.c`, `libprimer3.cc`) never calls `exit()`. It returns an error, and the `*_main.c` program picks the code.
  - When a program checks the same condition as another one, use the same code. For example, a salt concentration too low for the salt correction is 1 in `oligotm`, `ntthal` and `amplicon3_core`.
- **Failure values.** Calculations signal failure with sentinel values instead of status codes: `OLIGOTM_ERROR` (−999999.9999) from the Tm and ΔG functions, `THAL_ERROR_SCORE` from thal, and `DPAL_ERROR_SCORE` from dpal. `ALIGN_SCORE_UNDEF` marks a score that hasn't been computed yet. Check for them before using a value.
- **Fixed limits.** `PR_MAX_INTERVAL_ARRAY` (200 entries per interval or junction list), `MAX_PRIMER_LENGTH` (36), `THAL_MAX_ALIGN` (60) and `DPAL_MAX_ALIGN` (1600) size fixed arrays. When filling them from input, use the bounds-checked helpers such as `p3_add_to_interval_array()` and `p3_sa_add_to_overlap_junctions_array()`.
- **Shared buffers and process-wide state.**
  - `p3_read_line()` returns a static buffer that the next call overwrites and that is freed at end of file, so copy what you need. Settings files and mispriming libraries are read with the same function, which is why libraries are loaded only after the record loop in `read_boulder_record()`.
  - The alignment argument holders and `_jmp_buf` in `libprimer3.cc`, and the thal parameter tables, are file- or process-wide, so the library is not reentrant.

### Adding an input tag

Follow the pattern of commit 010c4e8:
1. Add a field in `libprimer3.h`: to `seq_args` if it applies to one record, to `p3_global_settings` if it should persist.
2. Set its default in `pr_set_default_global_args_1`.
3. Parse it in `read_boulder.c`.
4. Validate it in `_pr_data_control()`, including the allowed range of enum-like integers.
5. Implement the logic in `libprimer3.cc`.
6. If the tag affects output, update `print_boulder.c` and `format_output.c`. In `print_boulder.c`, guard every per-oligo line with `go_fwd`, `go_rev` or `go_int`, because list output can leave the oligo arrays empty.
7. Add test input and output, an entry in the release notes, and a proposed manual entry in `src/primer3_manual.htm` (see Overview).

When you add a value to an existing enum, update every place that checks it. For Tm methods (`tm_method_type`) that means `oligotm()`, `oligodg()`, `end_oligodg()` and `seqtm()`, the input validation, and the `-tp` check in `oligotm_main.c`. Salt corrections (`salt_correction_type`) work the same way with `-sc`.

### Component libraries

Each library also has a standalone CLI in a `*_main.c` file.
- `oligotm.c` (`oligotm`): melting temperature and salt corrections.
  - `seqtm()` is the entry point for any length. It uses `oligotm()` up to `MAX_NN_TM_LENGTH` (36) bases and `long_seq_tm()` above that.
  - `oligotm()` and `oligodg()` handle at most `MAX_PRIMER_LENGTH` bases.
- `dpal.c` (`ntdpal`): classic dynamic-programming alignment.
  - Always used for repeat-library mispriming.
  - Used for self-complementarity and template mispriming only when thermodynamic alignment is off (`thermodynamic_oligo_alignment` / `thermodynamic_template_alignment`).
- `thal*.c` (`ntthal`): thermodynamic alignment. `thal_internal.h` is private to these files.
  - `thal.c`: the DP core.
  - `thal_params.c`: parameter tables and loading.
  - `thal_trace.c`: traceback and debug dumps.
  - `thal_draw.c`: ASCII structure drawing.
  - `thal_util.c`: allocators.
  - `thal()` sets `thal_results.sec_struct`: NULL unless it drew a structure (`THL_STRUCT`), also after an error. The caller frees it.
  - For two sequences, `align_end_2` counts from the 3′ end of the second sequence, which thal reads in reverse: the terminal base pair is at index `strlen(seq2) - align_end_2`.
  - `align_thermod()` in libprimer3.cc reports a negative Tm as 0. Structures are drawn only for values above 0.
- Template mispriming (thermodynamic): `thermod_template_targets()` gives the oligo, `target` (the strand the primer binds at its own site; its own site is left out, so the 5′ and 3′ parts are aligned separately) and `target_r` (the other strand, aligned whole). The score and the `TEMPLATE_MISPRIMING_STUCT` structures both use it, so keep them together.
- `p3_seq_lib.c`: mispriming and repeat libraries (FASTA).
- `masker.c` (`primer3_masker`): k-mer based template masking. It is compiled out on Windows (`OS_WIN`).
  - The k-mer lists are binary GenomeTester4 files named `<prefix>_<k>.list`, with k = 11 and 16 and the default prefix `homo_sapiens`. They are memory-mapped: a header, then 12-byte records (an 8-byte k-mer and a 4-byte count). `kmer_lists/readme.txt` explains that the checked-in lists are test files and where to get full ones.
  - primer3_core looks for the lists in `PRIMER_MASK_KMERLIST_PATH`, then `../kmer_lists/`, then `/opt/kmer_lists/`. That is why masking tests run from `test/`.
- `amplicontm.c` (`amplicon3_core`): amplicon melting curves (Amplicon3). It is **GPL-3** because it derives from DECIPHER. The rest of the code is GPL-2+.

Makefile quirk: `dpal.c` and each `thal*.c` are compiled twice.
- `X.o` is used by the standalone executables.
- `X_primer.o`, built with `$(P_DEFINES)`, goes into the libraries linked into `primer3_core`.

A new thal source file therefore needs both rules, plus entries in `LIBTHAL` and the `NTTHAL_EXE` link line. There are no shared-library targets.

### Deliberate behavior

Don't "fix" these:
- `SEQUENCE_ID` is used unchanged as the base name of the `P3_FILE_FLAG` list files (`<SEQUENCE_ID>.for`, `.rev`, `.int`). Directory parts, including `..` and absolute paths, are allowed on purpose; don't reject `/`, `\` or a leading `.`.
- Internal oligos are not checked for template mispriming (`PRIMER_INTERNAL_MAX_TEMPLATE_MISHYB` is rejected as unsupported).
- The pair template-mispriming value is the sum of two separate bindings, one of each primer on opposite strands, in the orientation that could amplify a product: the larger of `left->template_mispriming + right->template_mispriming` (both on the strand they bind at their own site) and `left->template_mispriming_r + right->template_mispriming_r` (both on the other strand). The positions of the two bindings are not compared. It is used for the penalty, but there is no pair structure; only `PRIMER_LEFT/RIGHT_n_TEMPLATE_MISPRIMING_STUCT` exist.
- dpal's global modes (`g`, `G` in `ntdpal`) are semi-global, not Needleman-Wunsch: unaligned ends cost nothing.
- In `test/masker_lists.txt` the first path starts with `../test/` twelve times on purpose: it tests a bug with long paths in lists files.

### Settings files

`settings_files/*.txt` begin with three fixed lines: `Primer3 File - http://primer3.org`, `P3_FILE_TYPE=settings` and an empty line. In a settings file only `PRIMER_*` and `P3_FILE_ID` lines are read; everything else is skipped. They are loaded with `--p3_settings_file`, which the `*_settings` tests in `p3test.pl` exercise.

### Thermodynamic parameters

`src/primer3_config/*.ds` and `*.dh` hold the nearest-neighbor tables: entropy in `.ds`, enthalpy in `.dh`.
- Values are read by position. The leading base-pair labels and `#` comments are only there to orient the reader, so reordering lines changes meaning.
- The tables are compiled in through `src/thal_default_params.h`, which is **generated**. Don't edit it by hand. Edit the config files instead, then run `python3 thal_default_params_create.py` from `src/` (needs numpy).
- Changing values alters thal and primer3 output, so update the affected golden files in the same commit (see d2844b5).
- `PRIMER_THERMODYNAMIC_PARAMETERS_PATH` (and `ntthal -path`) load a directory of these files at run time through `thal_load_parameters()` and `get_thermodynamic_values()`. The loaded values replace the process-wide tables for the rest of the run.
- `test/primer3_config/` is an older copy in which some values differ (for example `stackmm`). Tests that set `PRIMER_THERMODYNAMIC_PARAMETERS_PATH` load it. It is not the source of truth.

To debug thal:
- `THAL_DPT_DUMP_MD=-` prints Markdown snapshots of the DP matrix to stderr. `THAL_DPT_DUMP_MD=<path>` writes them to a file.
- `THAL_DPT_DUMP_JSON=<path|->` writes a JSON array instead; the closing bracket is written at process exit.
- The variables are read once per process, and the file collects all `thal()` calls of the run (a primer3_core run or several ntthal lines). Each call starts with `## thal() call N` (Markdown) or carries `"call": N` (JSON).

## CI

The workflows in `.github/workflows/`:
- Ubuntu, macOS and Windows run on pushes to `main`, on pull requests and by hand (`workflow_dispatch`). Each test script runs in its own step, with `--action` for `p3test.pl` and `a3test.pl`. Windows leaves out the masker.
- The "paranoid" workflow runs on pushes to `main` and by hand. Its first job tests at -O2, with AddressSanitizer and UBSan, and at -O3; two more jobs run `make valgrnda` and `make valgrndb` on an -O0 build.
- After a failed test step the later test steps still run, as long as their build succeeded, and the `src` and `test` folders are always uploaded as artifacts. The job still fails.
