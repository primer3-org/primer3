/****************************************************************************
 *                   Simulates dsDNA Melting for Sequences                  *
 *                           Author: Erik Wright                            *
 ****************************************************************************/

/*
 * Rdefines.h is needed for the SEXP typedef, for the error(), INTEGER(),
 * GET_DIM(), LOGICAL(), NEW_INTEGER(), PROTECT() and UNPROTECT() macros,
 * and for the NA_INTEGER constant symbol.
 */
#include <Rdefines.h>

/*
 * R_ext/Rdynload.h is needed for the R_CallMethodDef typedef and the
 * R_registerRoutines() prototype.
 */
#include <R_ext/Rdynload.h>

/* for R_CheckUserInterrupt */
#include <R_ext/Utils.h>

// for math functions
#include <math.h>

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

//ans_start <- .Call("meltPolymer", myDNAStringSet, temps, ions, 1L, PACKAGE="DECIPHER")
SEXP meltPolymer(SEXP x, SEXP temps, SEXP ions, SEXP output)
{
	double Na = asReal(ions); // salt (Molar)
	int l = length(temps); // number of temperatures
	double *t = REAL(temps);
	// type of output
	// 1 = list of positional helicities
	// 2 = melt curve
	// 3 = derivative of curve
	int o = asInteger(output);
	
	XStringSet_holder x_set;
	x_set = hold_XStringSet(x);
	int x_length = get_length_from_XStringSet_holder(&x_set);
	Chars_holder x_s;
	int s, i, j;
	
	SEXP ret, ans;
	double *rans;
	if (o==1) {
		// list element for each sequence
		// row for each temperature
		// column for each sequence position
		PROTECT(ret = allocVector(VECSXP, x_length));
	} else {
		// row for each temperature
		// column for each sequence
		PROTECT(ret = allocMatrix(REALSXP, l, x_length));
		rans = REAL(ret);
		for (i = 0; i < l*x_length; i++)
			*(rans + i) = 0;
	}
	
	// constants
	double alpha = 2.15; // Blossey & Carlon, 2003
	double sigma = 1.26e-4; // Blossey & Carlon, 2003
	double Beta = 1e-7; // between 10^-6 and 10^-8
	double rescale_G = 1e1; // rescaling threshold
	double rescale_F = 1e-1; // rescaling factor
	double rF1 = 1e1;
	
	// Unified NN parameters (SantaLucia, 1998)
	// [A C G T]
	double dH[4][4] = {
		-7.9,-8.4,-7.8,-7.2
		,-8.5,-8.0,-10.6,-7.8
		,-8.2,-9.8,-8.0,-8.4
		,-7.2,-8.2,-8.5,-7.9
	};
	double dHini[4] = {2.3,0.1,0.1,2.3}; // ending [A C G T]
	double dH_010[4] = {-4.54,-4.54,-4.54,-4.54}; // isolated [A C G T]
	
	double dS[4][4] = {
		-22.2,-22.4,-21.0,-20.4
		,-22.7,-19.9,-27.2,-21.0
		,-22.2,-24.4,-19.9,-22.4
		,-21.3,-22.2,-22.7,-22.2
	};
	double dSini[4] = {4.1,-2.8,-2.8,4.1};
	double dS_010[4] = {-20.2,-20.2,-20.2,-20.2};
	
	double RT, s_11[4][4], s_end[4], s_010[4], Q_tot;
	double avg;
	double LEP = pow(5, -alpha);
	double ratio = pow(7, -alpha)/LEP;
	int rescale_i, seq_length, exp1;
	
	// variables for interpolating plateaus when o!=1
	double slope;
	
	for (s = 0; s < x_length; s++) {
		x_s = get_elt_from_XStringSet_holder(&x_set, s);
		int *seq = Calloc(x_s.length, int); // initialized to zero
		seq_length = 0; // x_s.length without non-DNA characters
		for (i = 0; i < x_s.length; i++) {
			switch (x_s.ptr[i]) {
				case 1: // A
					seq[seq_length] = 0;
					break;
				case 2: // C
					seq[seq_length] = 1;
					break;
				case 3: // AC
					seq[seq_length] = 0;
					break;
				case 4: // G
					seq[seq_length] = 2;
					break;
				case 5: // AG
					seq[seq_length] = 0;
					break;
				case 6: // CG
					seq[seq_length] = 1;
					break;
				case 7: // ACG
					seq[seq_length] = 0;
					break;
				case 8: // T
					seq[seq_length] = 3;
					break;
				case 9: // AT
					seq[seq_length] = 0;
					break;
				case 10: // CT
					seq[seq_length] = 3;
					break;
				case 11: // ACT
					seq[seq_length] = 0;
					break;
				case 12: // GT
					seq[seq_length] = 3;
					break;
				case 13: // AGT
					seq[seq_length] = 0;
					break;
				case 14: // CGT
					seq[seq_length] = 3;
					break;
				case 15: // ACGT
					seq[seq_length] = 0;
					break;
				default: // other (+, -, .)
					seq_length--;
					break;
			}
			seq_length++;
		}
		
		int *stack = Calloc(l, int); // initialized to zero
		
		if (o==1) {
			PROTECT(ans = allocMatrix(REALSXP, l, seq_length)); // [temp][pos]
			rans = REAL(ans);
		}
		
		int pos = 2, k = 0, it = 0;
		while (pos > 1) { // each temperature
			if (o==1) {
				if (k==l)
					break;
			} else { // stack-based iteration
				if (it < 2) { // initialize values
					if (it!=0) { // second temp
						if (l==1) // only one temp
							break;
						k = l - 1; // last temp
					}
					it++;
				} else {
					if (it==2) {
						if (l==2) // only two temps
							break;
						
						// initialize values
						stack[0] = 0;
						stack[1] = l - 1;
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
			
			RT = 0.0019871*(273.15 + t[k]);
			
			// statistical weights
			for (i = 0; i < 4; i++) {
				for (j = 0; j < 4; j++) {
					// statistical weight = exp(-dG/RT)
					s_11[i][j] = exp(-1*(dH[i][j] - (273.15 + t[k])*(dS[i][j] + 0.368*log(Na))/1000)/RT);
				}
			}
			for (i = 0; i < 4; i++) {
				s_end[i] = exp(-1*(dHini[i] - (273.15 + t[k])*(dSini[i] + 0.368*log(Na))/1000)/RT);
				s_010[i] = exp(-1*(dH_010[i] - (273.15 + t[k])*(dS_010[i] + 0.368*log(Na))/1000)/RT);
			}
			
			// allocate memory
			double *V_10_LR = Calloc(seq_length + 1, double); // initialized to zero
			double *U_01_LR = Calloc(seq_length, double); // initialized to zero
			double *U_11_LR = Calloc(seq_length, double); // initialized to zero
			int *rescale = Calloc(seq_length, int); // initialized to zero
			
			// initialize according to Tøstesen (2003)
			rescale_i = 0; // current rescaling exponent
			V_10_LR[0] = 1;
			V_10_LR[1] = Beta*s_010[seq[0]];
			U_01_LR[1] = Beta;
			U_11_LR[1] = Beta*s_end[seq[0]]*s_11[seq[0]][seq[1]];
			V_10_LR[2] = s_010[seq[1]]*U_01_LR[1] + s_end[seq[1]]*U_11_LR[1];
			Q_tot = V_10_LR[0] + V_10_LR[1] + V_10_LR[2];
			
			// recursive loop according to Tøstesen (2003)
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
				//Rprintf("\ni=%d rescale=%d V_10=%f U_01=%f U_11=%f V_10=%f", i, rescale_i, V_10_LR[i + 1], U_01_LR[i], U_11_LR[i], V_10_LR[0]);
			}
			
			// repeat the recursion from right-to-left
			// NOTE:  sequence complementation is unnecessary
			
			// allocate memory
			double *V_10_RL = Calloc(seq_length + 1, double); // initialized to zero
			double *U_01_RL = Calloc(seq_length, double); // initialized to zero
			double *U_11_RL = Calloc(seq_length, double); // initialized to zero
			
			// initialize according to Tøstesen (2003)
			V_10_RL[0] = 1;
			V_10_RL[1] = Beta*s_010[seq[seq_length - 1]];
			U_01_RL[1] = Beta;
			U_11_RL[1] = Beta*s_end[seq[seq_length - 1]]*s_11[seq[seq_length - 2]][seq[seq_length - 1]];
			V_10_RL[2] = s_010[seq[seq_length - 2]]*U_01_RL[1] + s_end[seq[seq_length - 2]]*U_11_RL[1];
			// NOTE:  Q_tot is only calculated one direction
			
			// recursive loop according to Tøstesen (2003)
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
                //Rprintf("\ni=%d rescale_F=%d V_10=%f U_01=%f U_11=%f V_10=%f", i, rescale_F, V_10_RL[i + 1], U_01_RL[i], U_11_RL[i], V_10_RL[0]);
			}
			
			// calculate p(i) = closed
			// (Tøstesen, 2003) Eq. 11
			if (o==1) {
				for (i = 1; i < seq_length - 1; i++) {
					*(rans + k + i*l) = (U_01_LR[i]*s_010[seq[i]]*U_01_RL[seq_length - i - 1] + U_01_LR[i]*s_end[seq[i]]*U_11_RL[seq_length - i - 1] + U_11_LR[i]*s_end[seq[i]]*U_01_RL[seq_length - i - 1] + U_11_LR[i]*U_11_RL[seq_length - i - 1])/(Beta*Q_tot);
					
					// there is a minor approximation error LR vs. RL
					// which sometimes results in p(i) slightly > 1
					if (*(rans + k + i*l) > 1)
						*(rans + k + i*l) = 1;
				}
				
				*(rans + k) = V_10_RL[seq_length]/Q_tot;
				if (*(rans + k) > 1)
					*(rans + k) = 1;
				*(rans + k + (seq_length - 1)*l) = V_10_LR[seq_length]/Q_tot;
				if (*(rans + k + (seq_length - 1)*l) > 1)
					*(rans + k + (seq_length - 1)*l) = 1;
				k++;
			} else {
				for (i = 1; i < seq_length - 1; i++) {
					avg = (U_01_LR[i]*s_010[seq[i]]*U_01_RL[seq_length - i - 1] + U_01_LR[i]*s_end[seq[i]]*U_11_RL[seq_length - i - 1] + U_11_LR[i]*s_end[seq[i]]*U_01_RL[seq_length - i - 1] + U_11_LR[i]*U_11_RL[seq_length - i - 1])/(Beta*Q_tot);
					
					// there is a minor approximation error LR vs. RL
					// which sometimes results in p(i) slightly > 1
					if (avg > 1) {
						*(rans + k + l*s) += 1;
					} else {
						*(rans + k + l*s) += avg;
					}
				}
				
				avg = V_10_RL[seq_length]/Q_tot;
				if (avg > 1) {
					*(rans + k + l*s) += 1;
				} else {
					*(rans + k + l*s) += avg;
				}
				avg = V_10_LR[seq_length]/Q_tot;
				if (avg > 1) {
					*(rans + k + l*s) += 1;
				} else {
					*(rans + k + l*s) += avg;
				}
				
				*(rans + k + l*s) /= seq_length; // average helicity
				
				if (it > 2) {
					// calculate slope
					slope = (*(rans + stack[0] + l*s) - *(rans + stack[pos] + l*s))/(t[stack[0]] - t[stack[pos]]);
					
					if (-1*slope < 1e-3) { // linearly interpolate points
						for (i = stack[0] + 1; i < stack[pos]; i++)
							*(rans + i + l*s) = slope*(t[i] - t[stack[0]]) + *(rans + stack[0] + l*s);
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
			
			Free(V_10_LR);
			Free(U_01_LR);
			Free(U_11_LR);
			Free(rescale);
			Free(V_10_RL);
			Free(U_01_RL);
			Free(U_11_RL);
			
			R_CheckUserInterrupt();
		}
		
		Free(seq);
		Free(stack);
		
		if (o==1) {
			SET_VECTOR_ELT(ret, s, ans);
			UNPROTECT(1); // ans
		} else if (o==3) { // negative derivative
			// cubic spline interpolation
			double a_11, a_12, a_22, a_23, a_33;
			double b_1, b_2, b_3;
			double k1, temp;
			for (k = 1; k < l - 1; k++) {
				a_12 = 1/(t[k] - t[k - 1]);
				a_11 = 2*a_12;
				a_23 = 1/(t[k + 1] - t[k]);
				a_22 = 2*(a_12 + a_23);
				a_33 = 2*a_23;
				b_1 = 3*(*(rans + k + l*s) - *(rans + k - 1 + l*s))/((t[k] - t[k - 1])*(t[k] - t[k - 1]));
				b_3 = 3*(*(rans + k + 1 + l*s) - *(rans + k + l*s))/((t[k + 1] - t[k])*(t[k + 1] - t[k]));
				b_2 = b_1 + b_3;
				
				k1 = (a_33*(a_11*b_2 - a_12*b_1) - b_3*a_23*a_11)/(a_33*(a_11*a_22 - a_12*a_12) - a_23*a_11*a_23);
				if (k==1) {
					*(rans + k - 1 + l*s) = -1*(b_1 - a_12*k1)/a_11; // -k0
				} else {
					*(rans + k - 1 + l*s) = -1*temp; // previous -k1
				}
				if (k==(l - 2)) {
					*(rans + k + l*s) = -1*k1; // current -k1
					*(rans + k + 1 + l*s) = -1*(b_3 - a_23*k1)/a_33; // -k2
					if (*(rans + k + 1 + l*s) < 0) // sometimes k2 is positive
						*(rans + k + 1 + l*s) = 0;
				}
				temp = k1;
			}
		}
	}
	
	UNPROTECT(1); // ret
	
	return ret;
}
