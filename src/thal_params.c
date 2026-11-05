/* thal_params.c — load thermodynamic parameter files into thal's table arrays.

   Split out of thal.c.  The default parameters live in thal_default_params.h
   (statically initialised) and are owned by thal.c; this file only provides
   the *loaders* that read user-supplied parameter files and overwrite those
   arrays.  The arrays themselves are passed in by pointer from
   get_thermodynamic_values() in thal.c. */

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "thal_internal.h"
/* Parameter tables (default-initialized) live in this TU only. */
#include "thal_default_params.h"

#define INIT_BUF_SIZE 1024

/* AT-pair penalty used by tableStartATS/ATH via get_thermodynamic_values. */
static const double AT_H = 2200.0;
static const double AT_S = 6.9;

/* ----- low-level line-reading helpers ----- */

/* Pop one '\n'-delimited line from `*str` (advancing *str past the
   newline) and return it as a heap-allocated null-terminated string.
   Returns NULL if no line remains. longjmps through `_jmp_buf` on
   allocation failure. */
static char *
th_read_str_line(char **str, jmp_buf _jmp_buf, thal_results *o)
{
   if (*str == NULL) {
      return NULL;
   }
   char *ptr = *str;
   char *ini = *str;
   while (1) {
      if ((*ptr == '\n') || (*ptr == '\0')) {
         char *ret = NULL;
         if (!(ret = (char *) malloc(sizeof(char) * (ptr - ini + 1)))) {
#ifdef DEBUG
            fputs("Error in malloc()\n", stderr);
#endif
            strcpy(o->msg, "Out of memory");
            errno = ENOMEM;
            longjmp(_jmp_buf, 1);
         }
         strncpy(ret, ini, (ptr - ini + 1));
         ret[ptr - ini] = '\0';

         if (*ptr == '\0') {
            *str = NULL;
         } else {
            ptr++;
            if (*ptr == '\0') {
               *str = NULL;
            } else {
               *str = ptr;
            }
         }
         if (ptr == ini) {
            if (ret != NULL) {
               free(ret);
            }
            return NULL;
         } else {
            return ret;
         }
      }
      ptr++;
   }
}

/* Read one number from *str. The number may be "inf" (mapped to _INFINITY).
   "inf" is handled here rather than via sscanf to keep Windows happy. */
static double
readDouble(char **str, jmp_buf _jmp_buf, thal_results *o)
{
   double result;
   char *line = th_read_str_line(str, _jmp_buf, o);
   char *p;
   /* NULL: the file has fewer values than the table needs, or it was
      not loaded at all. */
   if (line == NULL) {
      strcpy(o->msg, "Thermodynamic parameter file has too few values");
      longjmp(_jmp_buf, 1);
   }
   /* Skip leading spaces with p, so that line can still be freed. */
   p = line;
   while (isspace(*p)) p++;
   if (!strncmp(p, "inf", 3)) {
      free(line);
      return _INFINITY;
   }
   if (sscanf(p, "%lf", &result) != 1) {
      snprintf(o->msg, 255, "Illegal value in thermodynamic parameter file: %.200s", p);
      free(line);
      longjmp(_jmp_buf, 1);
   }
   free(line);
   return result;
}

/* Read a line with three doubles (interior / bulge / hairpin loop entries). */
static void
readLoop(char **str, double *v1, double *v2, double *v3, jmp_buf _jmp_buf, thal_results *o)
{
   char *line = th_read_str_line(str, _jmp_buf, o);
   char *p, *q;
   if (line == NULL) {
      strcpy(o->msg, "Thermodynamic loop parameter file has too few lines");
      longjmp(_jmp_buf, 1);
   }
   p = line;
   while (isspace(*p)) p++;
   while (isdigit(*p)) p++;
   while (isspace(*p)) p++;
   q = p;
   /* Stop at the end of the line; a line with fewer than three values
      would otherwise be read past its end. */
   while (*q != '\0' && !isspace(*q)) q++;
   if (*q == '\0') {
      strcpy(o->msg, "Thermodynamic loop parameter line has too few values");
      free(line);
      longjmp(_jmp_buf, 1);
   }
   *q = '\0'; q++;
   if (!strcmp(p, "inf")) *v1 = _INFINITY;
   else sscanf(p, "%lf", v1);
   while (isspace(*q)) q++;
   p = q;
   while (*p != '\0' && !isspace(*p)) p++;
   if (*p == '\0') {
      strcpy(o->msg, "Thermodynamic loop parameter line has too few values");
      free(line);
      longjmp(_jmp_buf, 1);
   }
   *p = '\0'; p++;
   if (!strcmp(q, "inf")) *v2 = _INFINITY;
   else sscanf(q, "%lf", v2);
   while (isspace(*p)) p++;
   q = p;
   while (!isspace(*q) && (*q != '\0')) q++;
   *q = '\0';
   if (!strcmp(p, "inf")) *v3 = _INFINITY;
   else sscanf(p, "%lf", v3);
   if (line != NULL) {
      free(line);
   }
}

/* Read a triloop (5-char) or tetraloop (6-char) entry: short string + double. */
static int
readTLoop(char **str, char *s, double *v, int triloop, jmp_buf _jmp_buf, thal_results *o)
{
   char *line = th_read_str_line(str, _jmp_buf, o);
   if (!line) return -1;
   char *p = line, *q;
   while (isspace(*p)) p++;
   q = p;
   while (isalpha(*q)) q++;
   *q = '\0'; q++;
   if (triloop) {
      strncpy(s, p, 5);
   } else {
      strncpy(s, p, 6);
   }
   while (isspace(*q)) q++;
   p = q;
   while (!isspace(*p) && (*p != '\0')) p++;
   *p = '\0';
   if (!strcmp(q, "inf")) *v = _INFINITY;
   else sscanf(q, "%lg", v);
   if (line != NULL) {
      free(line);
   }
   return 0;
}

/* Strip `#` line comments and collapse runs of blank / whitespace-only
   leading lines from an in-memory parameter file. Operates in place. */
static void
removeComments(char *data)
{
   /* Remove # comments */
   int i = 0;
   int w = 0;
   int comm = 0;
   while (1) {
      if (data[i] == '\0') {
         data[w] = '\0';
         break;
      }
      if (data[i] == '#') {
         comm = 1;
      }
      if (comm == 0) {
         data[w] = data[i];
         w++;
      }
      if (data[i] == '\n') {
         comm = 0;
      }
      i++;
   }
   /* Remove empty lines and spaces & tabs at the start of lines */
   i = 0;
   w = 0;
   int last_newline = 1;
   while (1) {
      if (data[i] == '\0') {
         data[w] = '\0';
         return;
      }
      if (last_newline == 0) {
         if (data[i] == '\n') {
            last_newline = 1;
         }
         data[w] = data[i];
         w++;
      } else {
         switch (data[i]) {
         case ' ':
         case '\t':
         case '\n':
            break;
         default:
            last_newline = 0;
            data[w] = data[i];
            w++;
         }
      }
      i++;
   }
}

/* Strip a leading identifier column (e.g. "AA_AT" prefixes on the
   stack/dangle/tstack tables) from an in-memory parameter file by
   dropping every leading character in the [ACTGN_ \t] set. Operates
   in place. */
static void
removeFristColumn(char *data)
{
   int i = 0;
   int w = 0;
   while (1) {
      if (data[i] == '\0') {
         data[w] = '\0';
         return;
      }
      switch (data[i]) {
      case 'A':
      case 'C':
      case 'T':
      case 'G':
      case 'N':
      case '_':
      case ' ':
      case '\t':
         break;
      default:
         data[w] = data[i];
         w++;
      }
      i++;
   }
}

/* Slurp <dirname>/<fname> into a heap-allocated null-terminated string. */
static char *
readParamFile(const char *dirname, const char *fname, jmp_buf _jmp_buf, thal_results *o)
{
   FILE *file;
   char *ret = NULL;
   char *paramdir = NULL;
   paramdir = (char *) safe_malloc(strlen(dirname) + strlen(fname) + 2, _jmp_buf, o);
   strcpy(paramdir, dirname);
#ifdef OS_WIN
   if (paramdir[strlen(paramdir) - 1] != '\\') {
      strcat(paramdir, "\\\0");
   }
#else
   if (paramdir[strlen(paramdir) - 1] != '/') {
      strcat(paramdir, "/\0");
   }
#endif
   strcat(paramdir, fname);
   if (!(file = fopen(paramdir, "r"))) {
      snprintf(o->msg, 255, "Unable to open file %s", paramdir);
      if (paramdir != NULL) {
         free(paramdir);
         paramdir = NULL;
      }
      longjmp(_jmp_buf, 1);
      return NULL;
   }
   if (paramdir != NULL) {
      free(paramdir);
      paramdir = NULL;
   }
   /* fgetc() returns an int: storing it in a char would make the EOF
      marker look like byte 0xFF and append it to the text. */
   int c;
   size_t i = 0;
   size_t ssz = INIT_BUF_SIZE;
   ret = (char *) safe_malloc(ssz, _jmp_buf, o);
   while ((c = fgetc(file)) != EOF) {
      /* Keep room for the terminating '\0'. */
      if (i + 1 >= ssz) {
         if (ssz >= INT_MAX / 2) {
            strcpy(o->msg, "Out of memory");
            free(ret);
            fclose(file);
            longjmp(_jmp_buf, 1);
         }
         ssz += INIT_BUF_SIZE;
         ret = (char *) safe_realloc(ret, ssz, _jmp_buf, o);
      }
      ret[i++] = (char) c;
   }
   ret[i] = '\0';
   fclose(file);
   removeComments(ret);
   return ret;
}

/* ----- top-level thal_parameters helpers (declared in thal.h) ----- */

int
thal_set_null_parameters(thal_parameters *a)
{
   a->dangle_dh = NULL;
   a->dangle_ds = NULL;
   a->loops_dh = NULL;
   a->loops_ds = NULL;
   a->stack_dh = NULL;
   a->stack_ds = NULL;
   a->stackmm_dh = NULL;
   a->stackmm_ds = NULL;
   a->tetraloop_dh = NULL;
   a->tetraloop_ds = NULL;
   a->triloop_dh = NULL;
   a->triloop_ds = NULL;
   a->tstack_tm_inf_ds = NULL;
   a->tstack_dh = NULL;
   a->tstack2_dh = NULL;
   a->tstack2_ds = NULL;
   return 0;
}

int
thal_free_parameters(thal_parameters *a)
{
   if (NULL != a->dangle_dh) { free(a->dangle_dh); a->dangle_dh = NULL; }
   if (NULL != a->dangle_ds) { free(a->dangle_ds); a->dangle_ds = NULL; }
   if (NULL != a->loops_dh) { free(a->loops_dh); a->loops_dh = NULL; }
   if (NULL != a->loops_ds) { free(a->loops_ds); a->loops_ds = NULL; }
   if (NULL != a->stack_dh) { free(a->stack_dh); a->stack_dh = NULL; }
   if (NULL != a->stack_ds) { free(a->stack_ds); a->stack_ds = NULL; }
   if (NULL != a->stackmm_dh) { free(a->stackmm_dh); a->stackmm_dh = NULL; }
   if (NULL != a->stackmm_ds) { free(a->stackmm_ds); a->stackmm_ds = NULL; }
   if (NULL != a->tetraloop_dh) { free(a->tetraloop_dh); a->tetraloop_dh = NULL; }
   if (NULL != a->tetraloop_ds) { free(a->tetraloop_ds); a->tetraloop_ds = NULL; }
   if (NULL != a->triloop_dh) { free(a->triloop_dh); a->triloop_dh = NULL; }
   if (NULL != a->triloop_ds) { free(a->triloop_ds); a->triloop_ds = NULL; }
   if (NULL != a->tstack_tm_inf_ds) { free(a->tstack_tm_inf_ds); a->tstack_tm_inf_ds = NULL; }
   if (NULL != a->tstack_dh) { free(a->tstack_dh); a->tstack_dh = NULL; }
   if (NULL != a->tstack2_dh) { free(a->tstack2_dh); a->tstack2_dh = NULL; }
   if (NULL != a->tstack2_ds) { free(a->tstack2_ds); a->tstack2_ds = NULL; }
   return 0;
}

int
thal_load_parameters(const char *path, thal_parameters *a, thal_results *o)
{
   jmp_buf _jmp_buf;
   thal_free_parameters(a);
   if (setjmp(_jmp_buf) != 0) {
      /* o->msg holds the reason. Don't leave half-loaded tables behind. */
      thal_free_parameters(a);
      return -1;
   }
   a->dangle_dh = readParamFile(path, "dangle.dh", _jmp_buf, o);
   removeFristColumn(a->dangle_dh);
   a->dangle_ds = readParamFile(path, "dangle.ds", _jmp_buf, o);
   removeFristColumn(a->dangle_ds);
   a->loops_dh = readParamFile(path, "loops.dh", _jmp_buf, o);
   a->loops_ds = readParamFile(path, "loops.ds", _jmp_buf, o);
   a->stack_dh = readParamFile(path, "stack.dh", _jmp_buf, o);
   removeFristColumn(a->stack_dh);
   a->stack_ds = readParamFile(path, "stack.ds", _jmp_buf, o);
   removeFristColumn(a->stack_ds);
   a->stackmm_dh = readParamFile(path, "stackmm.dh", _jmp_buf, o);
   removeFristColumn(a->stackmm_dh);
   a->stackmm_ds = readParamFile(path, "stackmm.ds", _jmp_buf, o);
   removeFristColumn(a->stackmm_ds);
   a->tetraloop_dh = readParamFile(path, "tetraloop.dh", _jmp_buf, o);
   a->tetraloop_ds = readParamFile(path, "tetraloop.ds", _jmp_buf, o);
   a->triloop_dh = readParamFile(path, "triloop.dh", _jmp_buf, o);
   a->triloop_ds = readParamFile(path, "triloop.ds", _jmp_buf, o);
   a->tstack_tm_inf_ds = readParamFile(path, "tstack_tm_inf.ds", _jmp_buf, o);
   removeFristColumn(a->tstack_tm_inf_ds);
   a->tstack_dh = readParamFile(path, "tstack.dh", _jmp_buf, o);
   removeFristColumn(a->tstack_dh);
   a->tstack2_dh = readParamFile(path, "tstack2.dh", _jmp_buf, o);
   removeFristColumn(a->tstack2_dh);
   a->tstack2_ds = readParamFile(path, "tstack2.ds", _jmp_buf, o);
   removeFristColumn(a->tstack2_ds);
   return 0;
}

/* ----- per-table loaders ----- */

void
getStack(double stackEntropies[5][5][5][5], double stackEnthalpies[5][5][5][5],
         const thal_parameters *tp, jmp_buf _jmp_buf, thal_results *o)
{
   int i, j, ii, jj;
   char *pt_ds = tp->stack_ds;
   char *pt_dh = tp->stack_dh;
   for (i = 0; i < 5; ++i) {
      for (ii = 0; ii < 5; ++ii) {
         for (j = 0; j < 5; ++j) {
            for (jj = 0; jj < 5; ++jj) {
               if (i == 4 || j == 4 || ii == 4 || jj == 4) {
                  stackEntropies[i][ii][j][jj] = -1.0;
                  stackEnthalpies[i][ii][j][jj] = _INFINITY;
               } else {
                  stackEntropies[i][ii][j][jj] = readDouble(&pt_ds, _jmp_buf, o);
                  stackEnthalpies[i][ii][j][jj] = readDouble(&pt_dh, _jmp_buf, o);
                  if (!isFinite(stackEntropies[i][ii][j][jj]) || !isFinite(stackEnthalpies[i][ii][j][jj])) {
                     stackEntropies[i][ii][j][jj] = -1.0;
                     stackEnthalpies[i][ii][j][jj] = _INFINITY;
                  }
               }
            }
         }
      }
   }
}

void
getStackint2(double stackint2Entropies[5][5][5][5], double stackint2Enthalpies[5][5][5][5],
             const thal_parameters *tp, jmp_buf _jmp_buf, thal_results *o)
{
   int i, j, ii, jj;
   char *pt_ds = tp->stackmm_ds;
   char *pt_dh = tp->stackmm_dh;
   for (i = 0; i < 5; ++i) {
      for (ii = 0; ii < 5; ++ii) {
         for (j = 0; j < 5; ++j) {
            for (jj = 0; jj < 5; ++jj) {
               if (i == 4 || j == 4 || ii == 4 || jj == 4) {
                  stackint2Entropies[i][ii][j][jj] = -1.0;
                  stackint2Enthalpies[i][ii][j][jj] = _INFINITY;
               } else {
                  stackint2Entropies[i][ii][j][jj] = readDouble(&pt_ds, _jmp_buf, o);
                  stackint2Enthalpies[i][ii][j][jj] = readDouble(&pt_dh, _jmp_buf, o);
                  if (!isFinite(stackint2Entropies[i][ii][j][jj]) || !isFinite(stackint2Enthalpies[i][ii][j][jj])) {
                     stackint2Entropies[i][ii][j][jj] = -1.0;
                     stackint2Enthalpies[i][ii][j][jj] = _INFINITY;
                  }
               }
            }
         }
      }
   }
}

void
getDangle(double dangleEntropies3[5][5][5], double dangleEnthalpies3[5][5][5],
          double dangleEntropies5[5][5][5], double dangleEnthalpies5[5][5][5],
          const thal_parameters *tp, jmp_buf _jmp_buf, thal_results *o)
{
   int i, j, k;
   char *pt_ds = tp->dangle_ds;
   char *pt_dh = tp->dangle_dh;
   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       for (k = 0; k < 5; ++k) {
          if (i == 4 || j == 4) {
             dangleEntropies3[i][k][j] = -1.0;
             dangleEnthalpies3[i][k][j] = _INFINITY;
          } else if (k == 4) {
             dangleEntropies3[i][k][j] = -1.0;
             dangleEnthalpies3[i][k][j] = _INFINITY;
          } else {
             dangleEntropies3[i][k][j] = readDouble(&pt_ds, _jmp_buf, o);
             dangleEnthalpies3[i][k][j] = readDouble(&pt_dh, _jmp_buf, o);
             if (!isFinite(dangleEntropies3[i][k][j]) || !isFinite(dangleEnthalpies3[i][k][j])) {
                dangleEntropies3[i][k][j] = -1.0;
                dangleEnthalpies3[i][k][j] = _INFINITY;
             }
          }
       }

   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       for (k = 0; k < 5; ++k) {
          if (i == 4 || j == 4) {
             dangleEntropies5[i][j][k] = -1.0;
             dangleEnthalpies5[i][j][k] = _INFINITY;
          } else if (k == 4) {
             dangleEntropies5[i][j][k] = -1.0;
             dangleEnthalpies5[i][j][k] = _INFINITY;
          } else {
             dangleEntropies5[i][j][k] = readDouble(&pt_ds, _jmp_buf, o);
             dangleEnthalpies5[i][j][k] = readDouble(&pt_dh, _jmp_buf, o);
             if (!isFinite(dangleEntropies5[i][j][k]) || !isFinite(dangleEnthalpies5[i][j][k])) {
                dangleEntropies5[i][j][k] = -1.0;
                dangleEnthalpies5[i][j][k] = _INFINITY;
             }
          }
       }
}

void
getLoop(double hairpinLoopEntropies[30], double interiorLoopEntropies[30], double bulgeLoopEntropies[30],
        double hairpinLoopEnthalpies[30], double interiorLoopEnthalpies[30], double bulgeLoopEnthalpies[30],
        const thal_parameters *tp, jmp_buf _jmp_buf, thal_results *o)
{
   int k;
   char *pt_ds = tp->loops_ds;
   char *pt_dh = tp->loops_dh;
   for (k = 0; k < 30; ++k) {
      readLoop(&pt_ds, &interiorLoopEntropies[k], &bulgeLoopEntropies[k], &hairpinLoopEntropies[k], _jmp_buf, o);
      readLoop(&pt_dh, &interiorLoopEnthalpies[k], &bulgeLoopEnthalpies[k], &hairpinLoopEnthalpies[k], _jmp_buf, o);
   }
}

void
getTstack(double tstackEntropies[5][5][5][5], double tstackEnthalpies[5][5][5][5],
          const thal_parameters *tp, jmp_buf _jmp_buf, thal_results *o)
{
   int i1, j1, i2, j2;
   char *pt_ds = tp->tstack_tm_inf_ds;
   char *pt_dh = tp->tstack_dh;
   for (i1 = 0; i1 < 5; ++i1)
     for (i2 = 0; i2 < 5; ++i2)
       for (j1 = 0; j1 < 5; ++j1)
         for (j2 = 0; j2 < 5; ++j2)
           if (i1 == 4 || j1 == 4) {
              tstackEnthalpies[i1][i2][j1][j2] = _INFINITY;
              tstackEntropies[i1][i2][j1][j2] = -1.0;
           } else if (i2 == 4 || j2 == 4) {
              tstackEntropies[i1][i2][j1][j2] = 0.00000000001;
              tstackEnthalpies[i1][i2][j1][j2] = 0.0;
           } else {
              tstackEntropies[i1][i2][j1][j2] = readDouble(&pt_ds, _jmp_buf, o);
              tstackEnthalpies[i1][i2][j1][j2] = readDouble(&pt_dh, _jmp_buf, o);
              if (!isFinite(tstackEntropies[i1][i2][j1][j2]) || !isFinite(tstackEnthalpies[i1][i2][j1][j2])) {
                 tstackEntropies[i1][i2][j1][j2] = -1.0;
                 tstackEnthalpies[i1][i2][j1][j2] = _INFINITY;
              }
           }
}

void
getTstack2(double tstack2Entropies[5][5][5][5], double tstack2Enthalpies[5][5][5][5],
           const thal_parameters *tp, jmp_buf _jmp_buf, thal_results *o)
{
   int i1, j1, i2, j2;
   char *pt_ds = tp->tstack2_ds;
   char *pt_dh = tp->tstack2_dh;
   for (i1 = 0; i1 < 5; ++i1)
     for (i2 = 0; i2 < 5; ++i2)
       for (j1 = 0; j1 < 5; ++j1)
         for (j2 = 0; j2 < 5; ++j2)
           if (i1 == 4 || j1 == 4) {
              tstack2Enthalpies[i1][i2][j1][j2] = _INFINITY;
              tstack2Entropies[i1][i2][j1][j2] = -1.0;
           } else if (i2 == 4 || j2 == 4) {
              tstack2Entropies[i1][i2][j1][j2] = 0.00000000001;
              tstack2Enthalpies[i1][i2][j1][j2] = 0.0;
           } else {
              tstack2Entropies[i1][i2][j1][j2] = readDouble(&pt_ds, _jmp_buf, o);
              tstack2Enthalpies[i1][i2][j1][j2] = readDouble(&pt_dh, _jmp_buf, o);
              if (!isFinite(tstack2Entropies[i1][i2][j1][j2]) || !isFinite(tstack2Enthalpies[i1][i2][j1][j2])) {
                 tstack2Entropies[i1][i2][j1][j2] = -1.0;
                 tstack2Enthalpies[i1][i2][j1][j2] = _INFINITY;
              }
           }
}

void
getTriloop(struct triloop **triloopEntropies, struct triloop **triloopEnthalpies,
           int *num, const thal_parameters *tp, jmp_buf _jmp_buf, thal_results *o)
{
   int i, size;
   double value;
   char *pt_ds = tp->triloop_ds;
   *num = 0;
   size = 16;
   if ((*triloopEntropies != NULL) && (*triloopEntropies != defaultTriloopEntropies)) {
      free(*triloopEntropies);
      *triloopEntropies = NULL;
   }
   *triloopEntropies = (struct triloop *) safe_calloc(16, sizeof(struct triloop), _jmp_buf, o);
   while (readTLoop(&pt_ds, (*triloopEntropies)[*num].loop, &value, 1, _jmp_buf, o) != -1) {
      for (i = 0; i < 5; ++i)
        (*triloopEntropies)[*num].loop[i] = str2int((*triloopEntropies)[*num].loop[i]);
      (*triloopEntropies)[*num].value = value;
      ++*num;
      if (*num == size) {
         size *= 2;
         *triloopEntropies = (struct triloop *) safe_realloc(*triloopEntropies, size * sizeof(struct triloop), _jmp_buf, o);
      }
   }
   *triloopEntropies = (struct triloop *) safe_realloc(*triloopEntropies, *num * sizeof(struct triloop), _jmp_buf, o);

   char *pt_dh = tp->triloop_dh;
   *num = 0;
   size = 16;
   if ((*triloopEnthalpies != NULL) && (*triloopEnthalpies != defaultTriloopEnthalpies)) {
      free(*triloopEnthalpies);
      *triloopEnthalpies = NULL;
   }
   *triloopEnthalpies = (struct triloop *) safe_calloc(16, sizeof(struct triloop), _jmp_buf, o);
   while (readTLoop(&pt_dh, (*triloopEnthalpies)[*num].loop, &value, 1, _jmp_buf, o) != -1) {
      for (i = 0; i < 5; ++i)
        (*triloopEnthalpies)[*num].loop[i] = str2int((*triloopEnthalpies)[*num].loop[i]);
      (*triloopEnthalpies)[*num].value = value;
      ++*num;
      if (*num == size) {
         size *= 2;
         *triloopEnthalpies = (struct triloop *) safe_realloc(*triloopEnthalpies, size * sizeof(struct triloop), _jmp_buf, o);
      }
   }
   *triloopEnthalpies = (struct triloop *) safe_realloc(*triloopEnthalpies, *num * sizeof(struct triloop), _jmp_buf, o);
}

void
getTetraloop(struct tetraloop **tetraloopEntropies, struct tetraloop **tetraloopEnthalpies,
             int *num, const thal_parameters *tp, jmp_buf _jmp_buf, thal_results *o)
{
   int i, size;
   double value;
   char *pt_ds = tp->tetraloop_ds;
   *num = 0;
   size = 16;
   if ((*tetraloopEntropies != NULL) && (*tetraloopEntropies != defaultTetraloopEntropies)) {
      free(*tetraloopEntropies);
      *tetraloopEntropies = NULL;
   }
   *tetraloopEntropies = (struct tetraloop *) safe_calloc(16, sizeof(struct tetraloop), _jmp_buf, o);
   while (readTLoop(&pt_ds, (*tetraloopEntropies)[*num].loop, &value, 0, _jmp_buf, o) != -1) {
      for (i = 0; i < 6; ++i)
        (*tetraloopEntropies)[*num].loop[i] = str2int((*tetraloopEntropies)[*num].loop[i]);
      (*tetraloopEntropies)[*num].value = value;
      ++*num;
      if (*num == size) {
         size *= 2;
         *tetraloopEntropies = (struct tetraloop *) safe_realloc(*tetraloopEntropies, size * sizeof(struct tetraloop), _jmp_buf, o);
      }
   }
   *tetraloopEntropies = (struct tetraloop *) safe_realloc(*tetraloopEntropies, *num * sizeof(struct tetraloop), _jmp_buf, o);

   char *pt_dh = tp->tetraloop_dh;
   *num = 0;
   size = 16;
   if ((*tetraloopEnthalpies != NULL) && (*tetraloopEnthalpies != defaultTetraloopEnthalpies)) {
      free(*tetraloopEnthalpies);
      *tetraloopEnthalpies = NULL;
   }
   *tetraloopEnthalpies = (struct tetraloop *) safe_calloc(16, sizeof(struct tetraloop), _jmp_buf, o);
   while (readTLoop(&pt_dh, (*tetraloopEnthalpies)[*num].loop, &value, 0, _jmp_buf, o) != -1) {
      for (i = 0; i < 6; ++i)
        (*tetraloopEnthalpies)[*num].loop[i] = str2int((*tetraloopEnthalpies)[*num].loop[i]);
      (*tetraloopEnthalpies)[*num].value = value;
      ++*num;
      if (*num == size) {
         size *= 2;
         *tetraloopEnthalpies = (struct tetraloop *) safe_realloc(*tetraloopEnthalpies, size * sizeof(struct tetraloop), _jmp_buf, o);
      }
   }
   *tetraloopEnthalpies = (struct tetraloop *) safe_realloc(*tetraloopEnthalpies, *num * sizeof(struct tetraloop), _jmp_buf, o);
}

void
tableStartATS(double atp_value, double atpS[5][5])
{
   int i, j;
   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       atpS[i][j] = 0.00000000001;
   atpS[0][3] = atpS[3][0] = atp_value;
}

void
tableStartATH(double atp_value, double atpH[5][5])
{
   int i, j;
   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       atpH[i][j] = 0.0;
   atpH[0][3] = atpH[3][0] = atp_value;
}

/* ----- top-level dispatch ----- */

/* Read the thermodynamic values (parameters) from the parameter strings
   in `tp` into the static parameter tables (defined in this TU via
   thal_default_params.h).  Return 0 on success and -1 on error. */
int
get_thermodynamic_values(const thal_parameters *tp, thal_results *o)
{
   jmp_buf _jmp_buf;
   if (setjmp(_jmp_buf) != 0) {
      return -1;
   }
   getStack(stackEntropies, stackEnthalpies, tp, _jmp_buf, o);
   getStackint2(stackint2Entropies, stackint2Enthalpies, tp, _jmp_buf, o);
   getDangle(dangleEntropies3, dangleEnthalpies3, dangleEntropies5, dangleEnthalpies5, tp, _jmp_buf, o);
   getLoop(hairpinLoopEntropies, interiorLoopEntropies, bulgeLoopEntropies, hairpinLoopEnthalpies,
           interiorLoopEnthalpies, bulgeLoopEnthalpies, tp, _jmp_buf, o);
   getTstack(tstackEntropies, tstackEnthalpies, tp, _jmp_buf, o);
   getTstack2(tstack2Entropies, tstack2Enthalpies, tp, _jmp_buf, o);
   getTriloop(&triloopEntropies, &triloopEnthalpies, &numTriloops, tp, _jmp_buf, o);
   getTetraloop(&tetraloopEntropies, &tetraloopEnthalpies, &numTetraloops, tp, _jmp_buf, o);
   tableStartATS(AT_S, atpS);
   tableStartATH(AT_H, atpH);
   return 0;
}

void
destroy_thal_structures(void)
{
   if ((triloopEntropies != NULL) && (triloopEntropies != defaultTriloopEntropies)) {
      free(triloopEntropies);
      triloopEntropies = NULL;
   }
   if ((triloopEnthalpies != NULL) && (triloopEnthalpies != defaultTriloopEnthalpies)) {
      free(triloopEnthalpies);
      triloopEnthalpies = NULL;
   }
   if ((tetraloopEntropies != NULL) && (tetraloopEntropies != defaultTetraloopEntropies)) {
      free(tetraloopEntropies);
      tetraloopEntropies = NULL;
   }
   if ((tetraloopEnthalpies != NULL) && (tetraloopEnthalpies != defaultTetraloopEnthalpies)) {
      free(tetraloopEnthalpies);
      tetraloopEnthalpies = NULL;
   }
}
