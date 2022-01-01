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

#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "amplicontm.h"

#include <stdio.h>


void free_amplicon_result(amplicon_result *res) {
    if (res->seq != NULL) {
        free(res->seq);
    }
    if (res->temp != NULL) {
        free(res->temp);
    }
    if (res->pos_prob != NULL) {
        free(res->pos_prob);
    }
    if (res->melt != NULL) {
        free(res->melt);
    }
    if (res->deriv != NULL) {
        free(res->deriv);
    }
}

void amp_init_all(void **alloc_box, int max_count) {
    for (int i = 0 ; i < max_count ; i++) {
        alloc_box[i] = NULL;
    }
}

void *amp_malloc(void **alloc_box, int *alloc_count, size_t size) {
    *alloc_count = (*alloc_count) + 1;
    alloc_box[*alloc_count] = malloc(size);
    return alloc_box[*alloc_count];
}

void amp_free_all(void **alloc_box, int alloc_count) {
    for (int i = 0 ; i < alloc_count ; i++) {
        free(alloc_box[i + 1]);
    }
}

void amp_zero_int(int *arr, int element_count) {
    if (arr == NULL) {
        return;
    }
    for(int i = 0 ; i < element_count ; i++) {
        arr[i] = 0;
    }
}

void amp_zero_double(double *arr, int element_count) {
    if (arr == NULL) {
        return;
    }
    for(int i = 0 ; i < element_count ; i++) {
        arr[i] = 0.0;
    }
}


/* Calculate the melting temperature of amplicon s.  See
   amplicontm.h and oligotm.h for documentation of arguments.
*/
amplicon_result amplicontm(const  char *inseq,
                           double mv,
                           double dv,
                           double dntp,
                           double dmso,
                           double dmso_fact,
                           double formamid,
                           amp_tm_parameters_type  tm_parameters,
                           amp_salt_correction_type salt_corrections,
                           amp_tm_method_type tm_formula,
                           int output) {
    int i, j;
    int max_count = 15;
    int alloc_count = 0;
    double schildkraut_corr = 0.0;
    double last_min = 0.0;
    void *alloc_box[max_count];
    amp_init_all(alloc_box, max_count);

    // Initialize return structure
    amplicon_result ret;
    ret.error = 0;
    ret.mv = mv;
    ret.dv = dv;
    ret.dntp = dntp;
    ret.dmso = dmso;
    ret.dmso_fact = dmso_fact;
    ret.formamid = formamid;
    ret.seq_len = 0;
    ret.seq = NULL;
    ret.seq_gc = 0.0;
    ret.temp_len = 0;
    ret.temp = NULL;
    ret.pos_prob = NULL;
    ret.melt = NULL;
    ret.deriv = NULL;
    ret.melt_len = 0;
    for (i = 0; i < 10; i++) {
      ret.melt_points[i] = -10.0;
    }

    mv /= 1000.0;
    dv /= 1000.0;
    dntp /= 1000.0;

    // Check Input
    if (tm_parameters != breslauer_amp && tm_parameters != santalucia_amp) {
        ret.error = 2;
        return ret;
    }
    if (tm_parameters == breslauer_amp) {
        ret.error = 5;
        return ret;
    }
    if (salt_corrections != schildkraut_amps && salt_corrections != santalucia_amps && salt_corrections != owczarzy_amps) {
        ret.error = 2;
        return ret;
    }
    if (salt_corrections == owczarzy_amps) {
        ret.error = 4;
        return ret;
    }
    if (mv < 0.01 || dv < 0.0) {
        ret.error = 2;
        return ret;
    }

    int orilen = strlen(inseq);
    int gc_count = 0;
    ret.seq = (char *) malloc(sizeof(char) * (orilen + 1));
    for (i = 0 ; i < orilen ; i++) {
        ret.seq[ret.seq_len]=toupper(inseq[i]);
        switch(ret.seq[ret.seq_len]) {
            case 'C':
            case 'G':
                gc_count++;
            case 'A':
            case 'T':
                ret.seq_len++;
                break;
        }
    }
    ret.seq[ret.seq_len] = '\0';
    ret.seq_gc = 100.0 * gc_count / ret.seq_len;

    if (tm_formula == bolton_amp) {
        ret.melt_points[0] = 81.5 - ret.dmso * ret.dmso_fact
                                  + (0.453 * ret.seq_gc / 100.0 - 2.88) * ret.formamid
                                  + 16.6 * (log10(mv + sqrt(dv) * 3.795))
                                  + 0.41 * ret.seq_gc
                                  - 600 / ret.seq_len;
        ret.melt_len = 1;
        return ret;
    }

    // Translate sequence to array positions
    int *seq = (int *) amp_malloc(alloc_box, &alloc_count, sizeof(int) * ret.seq_len);
    if (seq == NULL) {
        ret.error = 1;
        amp_free_all(alloc_box, alloc_count);
        return ret;
    }
    for (i = 0; i < ret.seq_len; i++) {
        switch(ret.seq[i]) {
            case 'A':
                seq[i] = 0;
                break;
            case 'C':
                seq[i] = 1;
                break;
            case 'G':
                seq[i] = 2;
                break;
            case 'T':
                seq[i] = 3;
                break;
        }
    }

    ret.temp_len = 950 - 600 + 1;  // Max Temp, Min Temp, incl. last
    ret.temp = (double *) malloc(sizeof(double) * ret.temp_len);
    if ((ret.seq == NULL) || (ret.temp == NULL)) {
        ret.error = 1;
        amp_free_all(alloc_box, alloc_count);
        return ret;
    }
    for (i = 0; i < ret.temp_len; i++) {
        ret.temp[i] = (600.0 + i) / 10;  // Min Temp
    }
    double *t = ret.temp;

    double *rans = NULL;
    if (output == 0) {
        // list element for each temperature * sequence position
        ret.pos_prob = (double *) malloc(sizeof(double) * ret.temp_len * ret.seq_len);
        if (ret.pos_prob == NULL) {
            ret.error = 1;
            amp_free_all(alloc_box, alloc_count);
            return ret;
        }
        amp_zero_double(ret.pos_prob, ret.temp_len * ret.seq_len);
        rans = ret.pos_prob;
    } else {
        // list element for each temperature
        ret.melt = (double *) malloc(sizeof(double) * ret.temp_len);
        ret.deriv = (double *) malloc(sizeof(double) * ret.temp_len);
        if ((ret.melt == NULL) || (ret.deriv == NULL)) {
            ret.error = 1;
            amp_free_all(alloc_box, alloc_count);
            return ret;
        }
        amp_zero_double(ret.melt, ret.temp_len);
        amp_zero_double(ret.deriv, ret.temp_len);
        rans = ret.melt;
    }

    // constants
    double alpha = 2.15; // Blossey & Carlon, 2003
    double sigma = 1.26e-4; // Blossey & Carlon, 2003
    double Beta = 1e-7; // between 10^-6 and 10^-8
    double rescale_G = 1e1; // rescaling threshold
    double rescale_F = 1e-1; // rescaling factor
    double rF1 = 1e1;

    // [A C G T]
    // [0 1 2 3]
    double dH[4][4];
    double dHini[4];
    double dH_010[4];
    double dS[4][4];
    double dSini[4];
    double dS_010[4];

    if (tm_parameters == breslauer_amp) {
        // NN parameters (Breslauer, 1986)
        dH[0][0] = -9.1;   // AA
        dH[0][1] = -6.5;   // AC
        dH[0][2] = -7.8;   // AG
        dH[0][3] = -8.6;   // AT
        dH[1][0] = -5.8;   // CA
        dH[1][1] = -11.0;  // CC
        dH[1][2] = -11.9;  // CG
        dH[1][3] = -7.8;   // CT
        dH[2][0] = -5.6;   // GA
        dH[2][1] = -11.1;  // GC
        dH[2][2] = -11.0;  // GG
        dH[2][3] = -6.5;   // GT
        dH[3][0] = -6.0;   // TA
        dH[3][1] = -5.6;   // TC
        dH[3][2] = -5.8;   // TG
        dH[3][3] = -9.1;   // TT

        dHini[0] = 2.3;    // A
        dHini[1] = 0.1;    // C
        dHini[2] = 0.1;    // G
        dHini[3] = 2.3;    // T

        dH_010[0] = -4.54; // A
        dH_010[1] = -4.54; // C
        dH_010[2] = -4.54; // G
        dH_010[3] = -4.54; // T

        dS[0][0] = -24.0;  // AA
        dS[0][1] = -17.3;  // AC
        dS[0][2] = -20.8;  // AG
        dS[0][3] = -23.9;  // AT
        dS[1][0] = -12.9;  // CA
        dS[1][1] = -26.6;  // CC
        dS[1][2] = -27.8;  // CG
        dS[1][3] = -20.8;  // CT
        dS[2][0] = -13.5;  // GA
        dS[2][1] = -26.7;  // GC
        dS[2][2] = -26.6;  // GG
        dS[2][3] = -17.3;  // GT
        dS[3][0] = -16.9;  // TA
        dS[3][1] = -13.5;  // TC
        dS[3][2] = -12.9;  // TG
        dS[3][3] = -24.0;  // TT

        dSini[0] =  4.1;   // A
        dSini[1] = -2.8;   // C
        dSini[2] = -2.8;   // G
        dSini[3] =  4.1;   // T

        dS_010[0] = -20.2; // A
        dS_010[1] = -20.2; // C
        dS_010[2] = -20.2; // G
        dS_010[3] = -20.2; // T
    } else {  // default santalucia_auto */
        // Unified NN parameters (SantaLucia, 1998)
        dH[0][0] = -7.9;   // AA
        dH[0][1] = -8.4;   // AC
        dH[0][2] = -7.8;   // AG
        dH[0][3] = -7.2;   // AT
        dH[1][0] = -8.5;   // CA
        dH[1][1] = -8.0;   // CC
        dH[1][2] = -10.6;  // CG
        dH[1][3] = -7.8;   // CT
        dH[2][0] = -8.2;   // GA
        dH[2][1] = -9.8;   // GC
        dH[2][2] = -8.0;   // GG
        dH[2][3] = -8.4;   // GT
        dH[3][0] = -7.2;   // TA
        dH[3][1] = -8.2;   // TC
        dH[3][2] = -8.5;   // TG
        dH[3][3] = -7.9;   // TT

        dHini[0] = 2.3;    // A
        dHini[1] = 0.1;    // C
        dHini[2] = 0.1;    // G
        dHini[3] = 2.3;    // T

        dH_010[0] = -4.54; // A
        dH_010[1] = -4.54; // C
        dH_010[2] = -4.54; // G
        dH_010[3] = -4.54; // T

        dS[0][0] = -22.2;  // AA
        dS[0][1] = -22.4;  // AC
        dS[0][2] = -21.0;  // AG
        dS[0][3] = -20.4;  // AT
        dS[1][0] = -22.7;  // CA
        dS[1][1] = -19.9;  // CC
        dS[1][2] = -27.2;  // CG
        dS[1][3] = -21.0;  // CT
        dS[2][0] = -22.2;  // GA
        dS[2][1] = -24.4;  // GC
        dS[2][2] = -19.9;  // GG
        dS[2][3] = -22.4;  // GT
        dS[3][0] = -21.3;  // TA
        dS[3][1] = -22.2;  // TC
        dS[3][2] = -22.7;  // TG
        dS[3][3] = -22.2;  // TT

        dSini[0] =  4.1;   // A
        dSini[1] = -2.8;   // C
        dSini[2] = -2.8;   // G
        dSini[3] =  4.1;   // T

        dS_010[0] = -20.2; // A
        dS_010[1] = -20.2; // C
        dS_010[2] = -20.2; // G
        dS_010[3] = -20.2; // T
    }

    if (salt_corrections != owczarzy_amps) {
        if ((dv > 0.0) && (dntp >= 0.0) && (dv > dntp)) {
            mv += sqrt(dv - dntp) * 3.795;
        }
    }

    double RT, s_11[4][4], s_end[4], s_010[4], Q_tot;
    double avg;
    double LEP = pow(5, -alpha);
    double ratio = pow(7, -alpha)/LEP;
    int rescale_i, seq_length, exp1;
    seq_length = ret.seq_len;

    // variables for interpolating plateaus when output!=0
    double slope;

    int *stack = (int *) amp_malloc(alloc_box, &alloc_count, sizeof(int) * ret.temp_len);
    amp_zero_int(stack, ret.temp_len);

    double *V_10_LR = (double *) amp_malloc(alloc_box, &alloc_count, sizeof(double) * (seq_length + 1));
    double *U_01_LR = (double *) amp_malloc(alloc_box, &alloc_count, sizeof(double) * seq_length);
    double *U_11_LR = (double *) amp_malloc(alloc_box, &alloc_count, sizeof(double) * seq_length);
    int *rescale = (int *) amp_malloc(alloc_box, &alloc_count, sizeof(int) * seq_length);
    double *V_10_RL = (double *) amp_malloc(alloc_box, &alloc_count, sizeof(double) * (seq_length + 1));
    double *U_01_RL = (double *) amp_malloc(alloc_box, &alloc_count, sizeof(double) * seq_length);
    double *U_11_RL = (double *) amp_malloc(alloc_box, &alloc_count, sizeof(double) * seq_length);

    if (stack == NULL || V_10_LR == NULL || U_01_LR == NULL || U_11_LR == NULL ||
        rescale == NULL || V_10_RL == NULL || U_01_RL == NULL || U_11_RL == NULL) {
        ret.error = 1;
        amp_free_all(alloc_box, alloc_count);
        return ret;
    }

    int pos = 2, k = 0, it = 0;
    while (pos > 1) { // each temperature
        if (output==0) {
            if (k == ret.temp_len)
                break;
        } else { // stack-based iteration
            if (it < 2) { // initialize values
                if (it!=0) { // second temp
                    if (ret.temp_len ==1 ) // only one temp
                        break;
                    k = ret.temp_len - 1; // last temp
                }
                it++;
            } else {
                if (it==2) {
                    if (ret.temp_len == 2) // only two temps
                        break;
                    // initialize values
                    stack[0] = 0;
                    stack[1] = ret.temp_len - 1;
                    it++;
                }
                stack[pos] = (stack[pos - 1] - stack[0])/2 + stack[0];
                if (stack[pos]==stack[0]) {
                    stack[0] = stack[pos - 1];
                    pos = pos - 1;
                    continue;
                }
                k = stack[pos];
            }
        }

        if (salt_corrections == schildkraut_amps) {
            schildkraut_corr = 16.6 * log10(mv);
            mv = 1.0;
        }

        RT = 1.9871*(273.15 + ret.dmso * ret.dmso_fact - schildkraut_corr
                            - (0.453 * ret.seq_gc / 100.0 - 2.88) * ret.formamid
                            + t[k]);
        // statistical weights
        for (i = 0; i < 4; i++) {
            for (j = 0; j < 4; j++) {
                // statistical weight = exp(-dG/RT)
                s_11[i][j] = exp((-1000*dH[i][j] +
                                  (273.15 + ret.dmso * ret.dmso_fact - schildkraut_corr
                                   - (0.453 * ret.seq_gc / 100.0 - 2.88) * ret.formamid + t[k])
                                  * (dS[i][j] + 0.368*log(mv))) / RT);
            }
        }
        for (i = 0; i < 4; i++) {
            s_end[i] = exp((-1000*dHini[i] +
                            (273.15 + ret.dmso * ret.dmso_fact - schildkraut_corr
                             - (0.453 * ret.seq_gc / 100.0 - 2.88) * ret.formamid + t[k])
                            * (dSini[i] + 0.368*log(mv))) / RT);
            s_010[i] = exp((-1000*dH_010[i] +
                             (273.15 + ret.dmso * ret.dmso_fact - schildkraut_corr
                              - (0.453 * ret.seq_gc / 100.0 - 2.88) * ret.formamid + t[k])
                             * (dS_010[i] + 0.368*log(mv))) / RT);
        }

        // initialize to zero
        amp_zero_double(V_10_LR, seq_length + 1);
        amp_zero_double(U_01_LR, seq_length);
        amp_zero_double(U_11_LR, seq_length);
        amp_zero_int(rescale, seq_length);

        // initialize according to Tostesen (2003)
        rescale_i = 0; // current rescaling exponent
        V_10_LR[0] = 1;
        V_10_LR[1] = Beta*s_010[seq[0]];
        U_01_LR[1] = Beta;
        U_11_LR[1] = Beta*s_end[seq[0]]*s_11[seq[0]][seq[1]];
        V_10_LR[2] = s_010[seq[1]]*U_01_LR[1] + s_end[seq[1]]*U_11_LR[1];
        Q_tot = V_10_LR[0] + V_10_LR[1] + V_10_LR[2];

        // recursive loop according to Tostesen (2003)
        for (i = 2; i < seq_length; i++) {
            U_01_LR[i] = Beta*V_10_LR[0];

            // loop entropy = SUM(sigma*(N + stiffness)^-alpha)
            if (i > 2)
                U_01_LR[i] += ratio*(U_01_LR[i - 1] - U_01_LR[i]);
            exp1 = rescale_i - rescale[i - 2];
            switch (exp1) {
                case 0:
                    U_01_LR[i] += V_10_LR[i - 1]*sigma*LEP;
                    break;
                case -1:
                    U_01_LR[i] += V_10_LR[i - 1]*rF1*sigma*LEP;
                    break;
                default:
                    U_01_LR[i] += V_10_LR[i - 1]*pow(rescale_F, exp1)*sigma*LEP;
                    break;
            }

            U_11_LR[i] = s_11[seq[i - 1]][seq[i]]*(U_01_LR[i - 1]*s_end[seq[i - 1]] + U_11_LR[i - 1]);
            V_10_LR[i + 1] = s_010[seq[i]]*U_01_LR[i] + s_end[seq[i]]*U_11_LR[i];
            Q_tot += V_10_LR[i + 1];

            // if the numbers are getting too large
            // then rescale current elements
            if (Q_tot > rescale_G && i < seq_length - 3) {
                rescale_i++;
                rescale[i] = rescale_i;

                Q_tot *= rescale_F;
                V_10_LR[i + 1] *= rescale_F;
                U_01_LR[i] *= rescale_F;
                U_11_LR[i] *= rescale_F;
                V_10_LR[0] *= rescale_F;
            } else {
                rescale[i] = rescale_i;
            }
        }

        // repeat the recursion from right-to-left
        // NOTE:  sequence complementation is unnecessary

        // initialize to zero
        amp_zero_double(V_10_RL, seq_length + 1);
        amp_zero_double(U_01_RL, seq_length);
        amp_zero_double(U_11_RL, seq_length);

        // initialize according to Tostesen (2003)
        V_10_RL[0] = 1;
        V_10_RL[1] = Beta*s_010[seq[seq_length - 1]];
        U_01_RL[1] = Beta;
        U_11_RL[1] = Beta*s_end[seq[seq_length - 1]]*s_11[seq[seq_length - 2]][seq[seq_length - 1]];
        V_10_RL[2] = s_010[seq[seq_length - 2]]*U_01_RL[1] + s_end[seq[seq_length - 2]]*U_11_RL[1];
        // NOTE:  Q_tot is only calculated one direction

        // recursive loop according to Tostesen (2003)
        for (i = 2; i < seq_length; i++) {
            U_01_RL[i] = Beta*V_10_RL[0];

            // loop entropy = SUM(sigma*(N + stiffness)^-alpha)
            if (i > 2)
                U_01_RL[i] += ratio*(U_01_RL[i - 1] - U_01_RL[i]);
            exp1 = rescale[seq_length - i - 1] - rescale[seq_length - i];
            switch (exp1) {
                case 0:
                    U_01_RL[i] += V_10_RL[i - 1]*sigma*LEP;
                    break;
                case -1:
                    U_01_RL[i] += V_10_RL[i - 1]*rF1*sigma*LEP;
                    break;
                default:
                    U_01_RL[i] += V_10_RL[i - 1]*pow(rescale_F, exp1)*sigma*LEP;
                    break;
            }

            U_11_RL[i] = s_11[seq[seq_length - i - 1]][seq[seq_length - i]]*(U_01_RL[i - 1]*s_end[seq[seq_length - i]] + U_11_RL[i - 1]);
            V_10_RL[i + 1] = s_010[seq[seq_length - i - 1]]*U_01_RL[i] + s_end[seq[seq_length - i - 1]]*U_11_RL[i];

            // rescale at the opposite indicies as LR
            if (rescale[seq_length - i - 1]!=rescale[seq_length - i]) {
                V_10_RL[i + 1] *= rescale_F;
                U_01_RL[i] *= rescale_F;
                U_11_RL[i] *= rescale_F;
                V_10_RL[0] *= rescale_F;
            }
        }
        //  return ret;

        // calculate p(i) = closed
        // (Tostesen, 2003) Eq. 11
        if (output==0) {
            for (i = 1; i < seq_length - 1; i++) {
                *(rans + k + i*ret.temp_len) = (U_01_LR[i]*s_010[seq[i]]*U_01_RL[seq_length - i - 1] + U_01_LR[i]*s_end[seq[i]]*U_11_RL[seq_length - i - 1] + U_11_LR[i]*s_end[seq[i]]*U_01_RL[seq_length - i - 1] + U_11_LR[i]*U_11_RL[seq_length - i - 1])/(Beta*Q_tot);

                // there is a minor approximation error LR vs. RL
                // which sometimes results in p(i) slightly > 1
                if (*(rans + k + i*ret.temp_len) > 1)
                    *(rans + k + i*ret.temp_len) = 1;
            }

            *(rans + k) = V_10_RL[seq_length]/Q_tot;
            if (*(rans + k) > 1)
                *(rans + k) = 1;
            *(rans + k + (seq_length - 1)*ret.temp_len) = V_10_LR[seq_length]/Q_tot;
            if (*(rans + k + (seq_length - 1)*ret.temp_len) > 1)
                *(rans + k + (seq_length - 1)*ret.temp_len) = 1;
            k++;
        } else {
            for (i = 1; i < seq_length - 1; i++) {
                avg = (U_01_LR[i]*s_010[seq[i]]*U_01_RL[seq_length - i - 1] + U_01_LR[i]*s_end[seq[i]]*U_11_RL[seq_length - i - 1] + U_11_LR[i]*s_end[seq[i]]*U_01_RL[seq_length - i - 1] + U_11_LR[i]*U_11_RL[seq_length - i - 1])/(Beta*Q_tot);

                // there is a minor approximation error LR vs. RL
                // which sometimes results in p(i) slightly > 1
                if (avg > 1) {
                    *(rans + k) += 1;
                } else {
                    *(rans + k) += avg;
                }
            }

            avg = V_10_RL[seq_length]/Q_tot;
            if (avg > 1) {
                *(rans + k) += 1;
            } else {
                *(rans + k) += avg;
            }
            avg = V_10_LR[seq_length]/Q_tot;
            if (avg > 1) {
                *(rans + k) += 1;
            } else {
                *(rans + k) += avg;
            }

            *(rans + k) /= seq_length; // average helicity

            if (it > 2) {
                // calculate slope
                slope = (*(rans + stack[0]) - *(rans + stack[pos]))/(t[stack[0]] - t[stack[pos]]);

                if (-1*slope < 1e-3) { // linearly interpolate points
                    for (i = stack[0] + 1; i < stack[pos]; i++)
                        *(rans + i) = slope*(t[i] - t[stack[0]]) + *(rans + stack[0]);
                    stack[0] = stack[pos];

                    stack[0] = stack[pos];
                    while ((stack[pos - 1] + 1L)==stack[pos]) {
                        stack[0] = stack[pos - 1];
                        pos--;
                    }
                } else {
                    pos++; // iterate
                }
            }
        }

    }

    if (output!=0) {
        // negative derivative
        // cubic spline interpolation
        double a_11, a_12, a_22, a_23, a_33;
        double b_1, b_2, b_3;
        double k1, temp;
        for (k = 1; k < ret.temp_len - 1; k++) {
            a_12 = 1/(t[k] - t[k - 1]);
            a_11 = 2*a_12;
            a_23 = 1/(t[k + 1] - t[k]);
            a_22 = 2*(a_12 + a_23);
            a_33 = 2*a_23;
            b_1 = 3*(*(rans + k) - *(rans + k - 1))/((t[k] - t[k - 1])*(t[k] - t[k - 1]));
            b_3 = 3*(*(rans + k + 1) - *(rans + k))/((t[k + 1] - t[k])*(t[k + 1] - t[k]));
            b_2 = b_1 + b_3;

            k1 = (a_33*(a_11*b_2 - a_12*b_1) - b_3*a_23*a_11)/(a_33*(a_11*a_22 - a_12*a_12) - a_23*a_11*a_23);
            if (k==1) {
                *(ret.deriv + k - 1) = -1*(b_1 - a_12*k1)/a_11; // -k0
            } else {
                *(ret.deriv + k - 1) = -1*temp; // previous -k1
            }
            if (k==(ret.temp_len - 2)) {
                *(ret.deriv + k) = -1*k1; // current -k1
                *(ret.deriv + k + 1) = -1*(b_3 - a_23*k1)/a_33; // -k2
                if (*(ret.deriv + k + 1) < 0) // sometimes k2 is positive
                    *(ret.deriv + k + 1) = 0;
            }
            temp = k1;
        }

        // Find the peaks
        last_min = *ret.deriv;
        for (k = ret.temp_len - 2; k > 0; k--) {
            if ((*(ret.deriv + k) - *(ret.deriv + k - 1) < 0.0000000001) &&
                (*(ret.deriv + k) - *(ret.deriv + k + 1) < 0.0000000001)) {
                last_min = *(ret.deriv + k);
            }
            if ((*(ret.deriv + k) - *(ret.deriv + k - 1) > 0.0000000001) &&
                (*(ret.deriv + k) - *(ret.deriv + k + 1) > 0.0000000001) &&
                (*(ret.deriv + k) - last_min > 0.001)) {
                if (ret.melt_len < 10) {
                    ret.melt_points[ret.melt_len] = *(t + k);
                    ret.melt_len++;
                }
            }
        }
    }

    amp_free_all(alloc_box, alloc_count);

    return ret;
}

amplicon_result ampliconfindsalt(const  char *inseq,
                                 double fs_temp,
                                 double dmso,
                                 double dmso_fact,
                                 double formamid,
                                 amp_tm_parameters_type tm_parameters,
                                 amp_salt_correction_type salt_corrections,
                                 amp_tm_method_type tm_formula,
                                 int output) {
    int i = 0;
    int cont = 1;
    double mv = 10.0;
    double mv_min = 10.0;
    double mv_max = 10.0;
    double step = 10.0;
    amplicon_result ret;

    fs_temp = ((int) (fs_temp * 10.0 + 0.5)) / 10.0;

    while ((cont == 1) && (i < 100)) {
        i++;
        ret = amplicontm(inseq, mv, 0.0, 0.0, dmso, dmso_fact, formamid,
                         tm_parameters, salt_corrections, tm_formula, 1);
        if (ret.error > 0) {
            free_amplicon_result(&ret);
            return ret;
        }
       if (ret.melt_len == 0) {
            ret.error = 3;
            free_amplicon_result(&ret);
            return ret;
        }

        if (fabs(ret.melt_points[0] - fs_temp) < 0.01) {
            break;
        }

        if (ret.melt_points[0] < fs_temp) {
            mv += step;
        } else {
            mv -= step;
            step /= 10.0;
            if (step < 0.5) {
                break;
            }
        }
        free_amplicon_result(&ret);
    }
    mv_min = mv;
    while (mv_min > 10) {
        mv_min -= 1.0;
        ret = amplicontm(inseq, mv_min, 0.0, 0.0, dmso, dmso_fact, formamid,
                         tm_parameters, salt_corrections, tm_formula, 1);
        if (ret.error > 0) {
            free_amplicon_result(&ret);
            return ret;
        }
        if (ret.melt_len == 0) {
            ret.error = 3;
            free_amplicon_result(&ret);
            return ret;
        }

        if (fabs(ret.melt_points[0] - fs_temp) > 0.01) {
            mv_min += 1.0;
            break;
        }
        free_amplicon_result(&ret);
    }
    mv_max = mv;
    while (mv_max < 1000) {
        mv_max += 1.0;
        ret = amplicontm(inseq, mv_max, 0.0, 0.0, dmso, dmso_fact, formamid,
                         tm_parameters, salt_corrections, tm_formula, 1);
        if (ret.error > 0) {
            free_amplicon_result(&ret);
            return ret;
        }
        if (ret.melt_len == 0) {
            ret.error = 3;
            free_amplicon_result(&ret);
            return ret;
        }

        if (fabs(ret.melt_points[0] - fs_temp) > 0.01) {
            mv_max -= 1.0;
            break;
        }
        free_amplicon_result(&ret);
    }
    mv = (mv_max + mv_min) / 2.0;
    return amplicontm(inseq, mv, 0.0, 0.0, dmso, dmso_fact, formamid,
                      tm_parameters, salt_corrections, tm_formula, output);
}
