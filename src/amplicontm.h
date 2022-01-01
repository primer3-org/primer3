/*
    Amplicon3 calculates melting temperatures for amplicons.
    Copyright (c) 2021 Andreas Untergasser

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    ----------------------------------------------------------------

    Amplicon3 is based on MeltPolymer.c from the DECIPHER package
    created by Erik Wright. DECIPHER is GPL-3 licensed.
*/

#ifndef _AMPLICON_TM
#define _AMPLICON_TM 1

#ifdef __cplusplus
  extern "C" {
#endif

/* error codes: 0 - no errors
                1 - out of memory

*/

typedef struct amplicon_result {
    int error;     /* See description above. */
    double mv;     /* Salt concentration (millimolar). */
    double dv;     /* Concentration of divalent cations (millimolar) */
    double dntp;   /* Concentration of dNTPs (millimolar) */
    double dmso;   /* Concentration of DMSO (%) */
    double dmso_fact;  /* DMSO correction factor (default 0.6) */
    double formamid;   /* Concentration of formamid (millimolar) */
    int seq_len;   /* Length of the sequence. */
    char *seq;     /* Sequence with ATCG only and in upper case. */
    double seq_gc;    /* GC content of the sequence in % . */
    int temp_len;  /* Number of temperatures. */
    double *temp;  /* The list of temperatures. */
    double *pos_prob;  /* The table of positional probabilities [temp_len][seq_len]. */
    double *melt;  /* The melt curve as list [temp_len]. */
    double *deriv;  /* The derivative of the melt curve as list [temp_len].  */
    int melt_len;  /* Number of melting points found. */
    double melt_points[10];  /* The list of temperatures. */
} amplicon_result;

void free_amplicon_result(amplicon_result *res);

void amp_init_all(void **alloc_box, int max_count);

void *amp_malloc(void **alloc_box, int *alloc_count, size_t size);

void amp_free_all(void **alloc_box, int alloc_count);

void amp_zero_int(int *arr, int element_count);

void amp_zero_double(double *arr, int element_count);

typedef enum amp_tm_parameters_type {
        breslauer_amp      = 0,
        santalucia_amp     = 1,
} amp_tm_parameters_type;

typedef enum amp_salt_correction_type {
        schildkraut_amps    = 0,
        santalucia_amps     = 1,
        owczarzy_amps       = 2,
} amp_salt_correction_type;

/*
   If tm_parameters==santalucia_amp, then the table of
   nearest-neighbor thermodynamic parameters and method for Tm
   calculation in the paper [SantaLucia JR (1998) "A unified view of
   polymer, dumbbell and oligonucleotide DNA nearest-neighbor
   thermodynamics", Proc Natl Acad Sci 95:1460-65
   http://dx.doi.org/10.1073/pnas.95.4.1460] is used.
   *THIS IS THE RECOMMENDED VALUE*.

   If tm_parameters==breslauer_amp, then method for Tm
   calculations in the paper [Rychlik W, Spencer WJ and Rhoads RE
   (1990) "Optimization of the annealing temperature for DNA
   amplification in vitro", Nucleic Acids Res 18:6409-12
   http://www.pubmedcentral.nih.gov/articlerender.fcgi?tool=pubmed&pubmedid=2243783].
   and the thermodynamic parameters in the paper [Breslauer KJ, Frank
   R, Bl?cker H and Marky LA (1986) "Predicting DNA duplex stability
   from the base sequence" Proc Natl Acad Sci 83:4746-50
   http://dx.doi.org/10.1073/pnas.83.11.3746], are is used.  This is
   the method and the table that primer3 used up to and including
   version 1.0.1

   If salt_corrections==schildkraut_amp, then formula for
   salt correction in the paper [Schildkraut, C, and Lifson, S (1965)
   "Dependence of the melting temperature of DNA on salt
   concentration", Biopolymers 3:195-208 (not available on-line)] is
   used.  This is the formula that primer3 used up to and including
   version 1.0.1.

   If salt_corrections==santalucia_amp, then formula for
   salt correction suggested by the paper [SantaLucia JR (1998) "A
   unified view of polymer, dumbbell and oligonucleotide DNA
   nearest-neighbor thermodynamics", Proc Natl Acad Sci 95:1460-65
   http://dx.doi.org/10.1073/pnas.95.4.1460] is used.

   *THIS IS THE RECOMMENDED VALUE*.

   The salt_corrections==owczarzy_amp, according to formula for
   salt correction in the paper [Owczarzy, R., Moreira, B.G., You, Y.,
   Behlke, M.A., and Walder, J.A. (2008) "Predicting stability of DNA
   duplexes in solutions containing magnesium and monovalent cations",
   Biochemistry 47:5336-53 http://dx.doi.org/10.1021/bi702363u] is not
   supported.

 */

typedef enum amp_tm_method_type {
        bolton_amp      = 0,
        wright_amp      = 1,
} amp_tm_method_type;

/*
   If tm_parameters==bolton_amp, then the the formula from Bolton and McCarthy,
   PNAS 84:1390 (1962) as presented in Sambrook, Fritsch and Maniatis,
   Molecular Cloning, p 11.46 (1989, CSHL Press).

   Tm = 81.5 + 16.6*(log10([Na+] + sqrt([Mg++])*3.795)) + 0.41*(%GC) - 600/length

   Where [Na+] is the molar sodium concentration, [Mg++] is the molar
   magnesium concentration, (%GC) is the percent of Gs and Cs in the
   sequence, and length is the length of the sequence. Be aware that
   the tm_parameters and salt_corrections are ignored using this
   formula.

   If tm_parameters==wright_amp, then the the formula from Wright and
   Vetsigian (2016) "DesignSignatures: a tool for designing primers
   that yields amplicons with distinct signatures.", Bioinformatics,
   32:1565-7 http://doi.org/10.1093/bioinformatics/btw047 is used.
   *THIS IS THE RECOMMENDED VALUE*.

 */




amplicon_result amplicontm(const  char *inseq,     /* The sequence. */
                           double monovalent_conc, /* Salt concentration (millimolar). */
                           double divalent_conc,   /* Concentration of divalent cations (millimolar) */
                           double dntp_conc,       /* Concentration of dNTPs (millimolar) */
                           double dmso,            /* Concentration of DMSO (%) */
                           double dmso_fact,       /* DMSO correction factor (default 0.6) */
                           double formamid,        /* Concentration of formamid (millimolar) */
                           amp_tm_parameters_type tm_parameters,       /* See description above. */
                           amp_salt_correction_type salt_corrections,  /* See description above. */
                           amp_tm_method_type tm_formula,              /* See description above. */
                           int output    /* 0 calc positional probabilities, 1 calc values, 2 calc 1 + curves */
                           );

amplicon_result ampliconfindsalt(const  char *inseq,  /* The sequence. */
                                 double fs_temp,      /* Desired melting temperature. */
                                 double dmso,         /* Concentration of DMSO (%) */
                                 double dmso_fact,    /* DMSO correction factor (default 0.6) */
                                 double formamid,     /* Concentration of formamid (millimolar) */
                                 amp_tm_parameters_type tm_parameters,       /* See description above. */
                                 amp_salt_correction_type salt_corrections,  /* See description above. */
                                 amp_tm_method_type tm_formula,              /* See description above. */
                                 int output    /* See description above. */
                                 );

#ifdef __cplusplus
  }
#endif

#endif
