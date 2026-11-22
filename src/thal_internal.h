/* Private header shared among thal.c, thal_trace.c, and thal_draw.c.
   Not part of the public API — do not include from outside this directory. */
#ifndef THAL_INTERNAL_H
#define THAL_INTERNAL_H

#include <math.h>
#include <setjmp.h>
#include <stddef.h>

#include "thal.h"  /* thal_mode, thal_results, etc. */

/* DPT cell.  The DP matrix in fillMatrix_{dimer,monomer} is
   `struct dpt_entry **`. */
struct dpt_entry {
   double h;
   double s;
   int tb_i;  /* traceback predecessor i, -2 = unvisited, -1 = no predecessor */
   int tb_j;  /* traceback predecessor j */
};

/* Defined in thal.c / thal_default_params.h. */
extern const double TEMP_KELVIN;
extern const double ABSOLUTE_ZERO;
extern const double MinEntropy;
extern const double _INFINITY;

#ifdef INTEGER
#  define isFinite(x) (x < _INFINITY / 2)
#else
#  define isFinite(x) isfinite(x)
#endif

/* Allocators that longjmp through `_jmp_buf` on OOM after writing an
   error message to `o`. Defined in thal_util.c. */
void *safe_calloc(size_t m, size_t n, jmp_buf _jmp_buf, thal_results *o);
void *safe_malloc(size_t n, jmp_buf _jmp_buf, thal_results *o);
void *safe_realloc(void *ptr, size_t n, jmp_buf _jmp_buf, thal_results *o);

/* Dump API.  No-ops when neither THAL_DPT_DUMP_MD nor THAL_DPT_DUMP_JSON
   is set in the environment.

     THAL_DPT_DUMP_MD=<path|->    markdown snapshots (- = stderr)
     THAL_DPT_DUMP_JSON=<path|->  JSON array of snapshots (- = stderr)

   The files stay open for all thal() calls of the process; see
   thal_trace.c. dpt_dump_close() ends one thal() call. */
void dpt_dump_close(void);
void dpt_dump(struct dpt_entry **dpt, int len1, int len2,
              int is_hairpin, const char *label);
void send5_dump(double *send5, double *hend5, int len, const char *label);

/* Secondary-structure text-art rendering. Defined in thal_draw.c. */
char *drawDimer(int *ps1, int *ps2, const thal_mode mode, double t37,
                const unsigned char *oligo1, const unsigned char *oligo2,
                int oligo1_len, int oligo2_len,
                jmp_buf _jmp_buf, thal_results *o);
char *drawHairpin(int *bp, double mh, double ms, const thal_mode mode,
                  double temp, const unsigned char *oligo1,
                  const unsigned char *oligo2, double saltCorrection,
                  int oligo1_len, int oligo2_len,
                  jmp_buf _jmp_buf, thal_results *o);

/* DNA-base → small-integer encoding (A=0, C=1, G=2, T=3, anything else=4).
   Defined in thal.c; needed by thal_params.c's triloop/tetraloop loaders. */
unsigned char str2int(char c);

/* Parameter file → table loaders. Each fills the destination array(s)
   from the corresponding `tp->*_dh` / `tp->*_ds` strings.  Defined and
   dispatched in thal_params.c, which also owns the (default-initialised)
   storage for the destination arrays. */
void getStack(double stackEntropies[5][5][5][5],
              double stackEnthalpies[5][5][5][5],
              const thal_parameters *tp, jmp_buf _jmp_buf, thal_results *o);
void getStackint2(double stackint2Entropies[5][5][5][5],
                  double stackint2Enthalpies[5][5][5][5],
                  const thal_parameters *tp, jmp_buf _jmp_buf, thal_results *o);
void getDangle(double dangleEntropies3[5][5][5],
               double dangleEnthalpies3[5][5][5],
               double dangleEntropies5[5][5][5],
               double dangleEnthalpies5[5][5][5],
               const thal_parameters *tp, jmp_buf _jmp_buf, thal_results *o);
void getLoop(double hairpinLoopEntropies[30],
             double interiorLoopEntropies[30],
             double bulgeLoopEntropies[30],
             double hairpinLoopEnthalpies[30],
             double interiorLoopEnthalpies[30],
             double bulgeLoopEnthalpies[30],
             const thal_parameters *tp, jmp_buf _jmp_buf, thal_results *o);
void getTstack(double tstackEntropies[5][5][5][5],
               double tstackEnthalpies[5][5][5][5],
               const thal_parameters *tp, jmp_buf _jmp_buf, thal_results *o);
void getTstack2(double tstack2Entropies[5][5][5][5],
                double tstack2Enthalpies[5][5][5][5],
                const thal_parameters *tp, jmp_buf _jmp_buf, thal_results *o);
void getTriloop(struct triloop **triloopEntropies,
                struct triloop **triloopEnthalpies,
                int *num, const thal_parameters *tp,
                jmp_buf _jmp_buf, thal_results *o);
void getTetraloop(struct tetraloop **tetraloopEntropies,
                  struct tetraloop **tetraloopEnthalpies,
                  int *num, const thal_parameters *tp,
                  jmp_buf _jmp_buf, thal_results *o);
void tableStartATS(double atp_value, double atpS[5][5]);
void tableStartATH(double atp_value, double atpH[5][5]);

/* Thermodynamic parameter tables.  Storage is provided by thal_params.c
   (the sole TU that includes thal_default_params.h, which carries the
   default initial values).  Everyone else reads them through these
   extern declarations. */
extern double atpS[5][5];
extern double atpH[5][5];
extern double dangleEntropies3[5][5][5];
extern double dangleEnthalpies3[5][5][5];
extern double dangleEntropies5[5][5][5];
extern double dangleEnthalpies5[5][5][5];
extern double stackEntropies[5][5][5][5];
extern double stackEnthalpies[5][5][5][5];
extern double stackint2Entropies[5][5][5][5];
extern double stackint2Enthalpies[5][5][5][5];
extern double interiorLoopEntropies[30];
extern double bulgeLoopEntropies[30];
extern double hairpinLoopEntropies[30];
extern double interiorLoopEnthalpies[30];
extern double bulgeLoopEnthalpies[30];
extern double hairpinLoopEnthalpies[30];
extern double tstackEntropies[5][5][5][5];
extern double tstackEnthalpies[5][5][5][5];
extern double tstack2Entropies[5][5][5][5];
extern double tstack2Enthalpies[5][5][5][5];
extern int numTriloops;
extern int numTetraloops;
extern struct triloop *triloopEntropies;
extern struct triloop *triloopEnthalpies;
extern struct tetraloop *tetraloopEntropies;
extern struct tetraloop *tetraloopEnthalpies;
extern struct triloop defaultTriloopEntropies[];
extern struct triloop defaultTriloopEnthalpies[];
extern struct tetraloop defaultTetraloopEntropies[];
extern struct tetraloop defaultTetraloopEnthalpies[];

#endif /* THAL_INTERNAL_H */
