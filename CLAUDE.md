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
- Windows builds use `make TESTOPTS=--windows` (MinGW). That build leaves out the masker.

## Tests

Run the test scripts from `test/`. They find the binaries through `../src/`.

```sh
cd src && make test                               # build, then run every suite (~1.5 min)
cd test && perl p3test.pl --onetest primer_task   # run one primer3_core test
perl p3test.pl --fast --action                    # skip the slowest tests
perl thal_test.pl     # also: dpal_test.pl oligotm_test.pl masker_test.pl cmdline_test.pl a3test.pl
```

Check results by searching the output for `[FAILED]`:
- `p3test.pl` and `a3test.pl` exit 0 even when tests fail, unless you pass `--action`.
- The other scripts exit non-zero on failure, which stops `make test` before the later suites run.
- A plain `FAILED` also matches the reminder line that `make test` prints at the end.

To run under valgrind, first build with `make O_OPTS=-O0`. Then, in `test/`, run `make valgrnda` and `make valgrndb` (how CI splits it), or `make TESTOPTS=--valgrind` for everything. `src/paranoid_tests.sh` runs the full -O2, -O3 and valgrind sequence.

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
- The stdout of all calls is concatenated into one `.tmp` file and compared with one golden file. stderr is never compared: thal and masker write it to `.tmp.err`, while dpal leaves it on the terminal.
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
- **Validation order.** `choose_primers()` calls `_adjust_seq_args()` before `_pr_data_control()`, so code in `_adjust_seq_args()` sees input that hasn't been validated yet.
- **Coordinates.** Input positions count from `PRIMER_FIRST_BASE_INDEX` (default 0). `_adjust_seq_args()` makes them 0-based and relative to the included region. After that, intervals, junction lists, force positions, `start_codon_pos` and oligo positions index `trimmed_seq`, and `incl_s` is the 0-based start of the region in the full `sequence`. Output adds `incl_s + first_base_index` back.
- **Sequence copies in `seq_args`:**
  - `sequence`: the template as given.
  - `trimmed_seq`: the included region; `_pr_data_control()` upper-cases it and turns unknown bases into N.
  - `trimmed_orig_seq`: the included region in its original case, used for lowercase masking.
  - `upcased_seq` and `upcased_seq_r`: the whole template and its reverse complement, used for template mispriming.
  - `trimmed_masked_seq` and `trimmed_masked_seq_r`: the k-mer-masked region for each strand.
- **Oligo records (`primer_rec`).** `start` is the 5′ base: the leftmost base of a left primer or internal oligo, and the rightmost base of a right primer. `rev_comp` marks the printed orientation. For an internal oligo it is set when `PRIMER_INTERNAL_REV_COMP_G_RICH` prints the reverse complement, so checks that depend on orientation must use the printed sequence.
- **Errors.**
  - In `read_boulder.c`, a bad `SEQUENCE_*` tag is a per-record error: that record prints `PRIMER_ERROR` and the next record runs. A bad global tag is fatal: primer3_core prints the error and exits with -4.
  - In `p3retval`, `glob_err` is fatal, `per_sequence_err` ends only the current record, and `warnings` are printed as `PRIMER_WARNING`.
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
- `p3_seq_lib.c`: mispriming and repeat libraries (FASTA).
- `masker.c` (`primer3_masker`): k-mer based template masking. It is compiled out on Windows (`OS_WIN`).
  - The k-mer lists are binary GenomeTester4 files named `<prefix>_<k>.list`, with k = 11 and 16 and the default prefix `homo_sapiens`. They are memory-mapped: a header, then 12-byte records (an 8-byte k-mer and a 4-byte count). `kmer_lists/readme.txt` explains that the checked-in lists are test files and where to get full ones.
  - primer3_core looks for the lists in `PRIMER_MASK_KMERLIST_PATH`, then `../kmer_lists/`, then `/opt/kmer_lists/`. That is why masking tests run from `test/`.
- `amplicontm.c` (`amplicon3_core`): amplicon melting curves (Amplicon3). It is **GPL-3** because it derives from DECIPHER. The rest of the code is GPL-2+.

Makefile quirk: `dpal.c` and each `thal*.c` are compiled twice.
- `X.o` is used by the standalone executables.
- `X_primer.o`, built with `$(P_DEFINES)`, goes into the libraries linked into `primer3_core`.

A new thal source file therefore needs both rules, plus entries in `LIBTHAL`, `LIBTHAL_DYN` and the `NTTHAL_EXE` link line.

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
- `THAL_DPT_DUMP_JSON=<path>` writes a JSON stream instead.

## CI

The workflows in `.github/workflows/` run only on pushes to `main`.
- On Ubuntu, macOS and Windows, each test script runs separately, with `--action` for `p3test.pl` and `a3test.pl`.
- A "paranoid" job repeats the tests at -O2 and -O3 and under valgrind.
