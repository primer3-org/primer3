/* thal_draw.c — secondary-structure text-art renderers for ntthal.
   Split out of thal.c (this is only invoked when mode != THL_FAST,
   so it's orthogonal to the DP algorithm itself). */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "thal_internal.h"

/* String-builder helpers — only used by the drawDimer/drawHairpin output
   formatting below. */
static void
strcatc(char *str, char c)
{
   str[strlen(str) + 1] = 0;
   str[strlen(str)] = c;
}

static void
save_append_string(char **ret, int *space, thal_results *o,
                   const char *str, jmp_buf _jmp_buf)
{
   int xlen, slen;
   if (str == NULL) {
      return;
   }
   if (*ret == NULL) {
      *ret = (char *) safe_malloc(sizeof(char) * 500, _jmp_buf, o);
      *ret[0] = '\0';
      *space = 500;
   }
   xlen = strlen(*ret);
   slen = strlen(str);
   if (xlen + slen + 1 > *space) {
      *space += 4 * (slen + 1);
      *ret = (char *) safe_realloc(*ret, *space, _jmp_buf, o);
   }
   strcpy(*ret + xlen, str);
}

static void
save_append_char(char **ret, int *space, thal_results *o,
                 const char str, jmp_buf _jmp_buf)
{
   char fix[3];
   fix[0] = str;
   fix[1] = '\0';
   save_append_string(ret, space, o, fix, _jmp_buf);
}

char *
drawDimer(int *ps1, int *ps2, const thal_mode mode, double t37,
          const unsigned char *oligo1, const unsigned char *oligo2,
          int oligo1_len, int oligo2_len, jmp_buf _jmp_buf, thal_results *o)
{
   (void) t37; /* parameter is unused; kept for ABI compat with thal_internal.h */
   int  ret_space = 0;
   char *ret_ptr = NULL;
   int ret_nr, ret_pr_once;
   char ret_para[400];
   char* ret_str[4];
   int i, j, k, numSS1, numSS2;
   char* duplex[4];

   if (mode != THL_STRUCT) {
      printf("Calculated thermodynamical parameters for dimer:\tdS = %g\tdH = %g\tdG = %g\tt = %g\n",
            (double) o->ds, (double) o->dh, (double) o->dg, (double) o->temp);
   } else {
      snprintf(ret_para, 400, "Tm: %.1f&deg;C  dG: %.0f cal/mol  dH: %.0f cal/mol  dS: %.0f cal/mol*K\\n",
               (double) o->temp, (double) o->dg, (double) o->dh, (double) o->ds);
   }

   duplex[0] = (char*) safe_malloc(oligo1_len + oligo2_len + 1, _jmp_buf, o);
   duplex[1] = (char*) safe_malloc(oligo1_len + oligo2_len + 1, _jmp_buf, o);
   duplex[2] = (char*) safe_malloc(oligo1_len + oligo2_len + 1, _jmp_buf, o);
   duplex[3] = (char*) safe_malloc(oligo1_len + oligo2_len + 1, _jmp_buf, o);
   duplex[0][0] = duplex[1][0] = duplex[2][0] = duplex[3][0] = 0;

   i = 0;
   numSS1 = 0;
   while (ps1[i++] == 0) ++numSS1;
   j = 0;
   numSS2 = 0;
   while (ps2[j++] == 0) ++numSS2;

   if (numSS1 >= numSS2){
      for (i = 0; i < numSS1; ++i) {
         strcatc(duplex[0], oligo1[i]);
         strcatc(duplex[1], ' ');
         strcatc(duplex[2], ' ');
      }
      for (j = 0; j < numSS1 - numSS2; ++j) strcatc(duplex[3], ' ');
      for (j = 0; j < numSS2; ++j) strcatc(duplex[3], oligo2[j]);
   } else {
      for (j = 0; j < numSS2; ++j) {
         strcatc(duplex[3], oligo2[j]);
         strcatc(duplex[1], ' ');
         strcatc(duplex[2], ' ');
      }
      for (i = 0; i < numSS2 - numSS1; ++i)
        strcatc(duplex[0], ' ');
      for (i = 0; i < numSS1; ++i)
        strcatc(duplex[0], oligo1[i]);
   }
   i = numSS1 + 1;
   j = numSS2 + 1;

   while (i <= oligo1_len) {
      while (i <= oligo1_len && ps1[i - 1] != 0 && j <= oligo2_len && ps2[j - 1] != 0) {
         strcatc(duplex[0], ' ');
         strcatc(duplex[1], oligo1[i - 1]);
         strcatc(duplex[2], oligo2[j - 1]);
         strcatc(duplex[3], ' ');
         ++i;
         ++j;
      }
      numSS1 = 0;
      while (i <= oligo1_len && ps1[i - 1] == 0) {
         strcatc(duplex[0], oligo1[i - 1]);
         strcatc(duplex[1], ' ');
         ++numSS1;
         ++i;
      }
      numSS2 = 0;
      while (j <= oligo2_len && ps2[j - 1] == 0) {
         strcatc(duplex[2], ' ');
         strcatc(duplex[3], oligo2[j - 1]);
         ++numSS2;
         ++j;
      }
      if (numSS1 < numSS2)
        for (k = 0; k < numSS2 - numSS1; ++k) {
           strcatc(duplex[0], '-');
           strcatc(duplex[1], ' ');
        }
      else if (numSS1 > numSS2)
        for (k = 0; k < numSS1 - numSS2; ++k) {
           strcatc(duplex[2], ' ');
           strcatc(duplex[3], '-');
        }
   }
   if (mode == THL_GENERAL) {
     printf("SEQ\t");
     printf("%s\n", duplex[0]);
     printf("SEQ\t");
     printf("%s\n", duplex[1]);
     printf("STR\t");
     printf("%s\n", duplex[2]);
     printf("STR\t");
     printf("%s\n", duplex[3]);
   }
   if (mode == THL_STRUCT) {
     ret_str[3] = NULL;
     ret_str[0] = (char*) safe_malloc(oligo1_len + oligo2_len + 10, _jmp_buf, o);
     ret_str[1] = (char*) safe_malloc(oligo1_len + oligo2_len + 10, _jmp_buf, o);
     ret_str[2] = (char*) safe_malloc(oligo1_len + oligo2_len + 10, _jmp_buf, o);
     ret_str[0][0] = ret_str[1][0] = ret_str[2][0] = '\0';

     /* Join top primer */
     strcpy(ret_str[0], "   ");
     strcat(ret_str[0], duplex[0]);
     ret_nr = 0;
     while (duplex[1][ret_nr] != '\0') {
       if (duplex[1][ret_nr] == 'A' || duplex[1][ret_nr] == 'T' ||
           duplex[1][ret_nr] == 'C' || duplex[1][ret_nr] == 'G' ||
           duplex[1][ret_nr] == '-') {
         ret_str[0][ret_nr + 3] = duplex[1][ret_nr];
       }
       ret_nr++;
     }
     if (strlen(duplex[1]) > strlen(duplex[0])) {
       ret_str[0][strlen(duplex[1]) + 3] = '\0';
     }
     /* Clean Ends */
     ret_nr = strlen(ret_str[0]) - 1;
     while (ret_nr > 0 && (ret_str[0][ret_nr] == ' ' || ret_str[0][ret_nr] == '-')) {
       ret_str[0][ret_nr--] = '\0';
     }
     /* Write the 5' */
     ret_nr = 3;
     ret_pr_once = 1;
     while (ret_str[0][ret_nr] != '\0' && ret_pr_once == 1) {
       if (ret_str[0][ret_nr] == 'A' || ret_str[0][ret_nr] == 'T' ||
           ret_str[0][ret_nr] == 'C' || ret_str[0][ret_nr] == 'G' ||
           ret_str[0][ret_nr] == '-') {
         ret_str[0][ret_nr - 3] = '5';
         ret_str[0][ret_nr - 2] = '\'';
         ret_pr_once = 0;
       }
       ret_nr++;
     }

     /* Create the align tics */
     strcpy(ret_str[1], "     ");
     for (i = 0 ; i < strlen(duplex[1]) ; i++) {
       if (duplex[1][i] == 'A' || duplex[1][i] == 'T' ||
           duplex[1][i] == 'C' || duplex[1][i] == 'G' ) {
         ret_str[1][i + 3] = '|';
       } else {
         ret_str[1][i + 3] = ' ';
       }
       ret_str[1][i + 4] = '\0';
     }
     /* Clean Ends */
     ret_nr = strlen(ret_str[1]) - 1;
     while (ret_nr > 0 && ret_str[1][ret_nr] == ' ') {
       ret_str[1][ret_nr--] = '\0';
     }
     /* Join bottom primer */
     strcpy(ret_str[2], "   ");
     strcat(ret_str[2], duplex[2]);
     ret_nr = 0;
     while (duplex[3][ret_nr] != '\0') {
       if (duplex[3][ret_nr] == 'A' || duplex[3][ret_nr] == 'T' ||
           duplex[3][ret_nr] == 'C' || duplex[3][ret_nr] == 'G' ||
           duplex[3][ret_nr] == '-') {
         ret_str[2][ret_nr + 3] = duplex[3][ret_nr];
       }
       ret_nr++;
     }
     if (strlen(duplex[3]) > strlen(duplex[2])) {
       ret_str[2][strlen(duplex[3]) + 3] = '\0';
     }
     /* Clean Ends */
     ret_nr = strlen(ret_str[2]) - 1;
     while (ret_nr > 0 && (ret_str[2][ret_nr] == ' ' || ret_str[2][ret_nr] == '-')) {
       ret_str[2][ret_nr--] = '\0';
     }
     /* Write the 5' */
     ret_nr = 3;
     ret_pr_once = 1;
     while (ret_str[2][ret_nr] != '\0' && ret_pr_once == 1) {
       if (ret_str[2][ret_nr] == 'A' || ret_str[2][ret_nr] == 'T' ||
           ret_str[2][ret_nr] == 'C' || ret_str[2][ret_nr] == 'G' ||
           ret_str[2][ret_nr] == '-') {
         ret_str[2][ret_nr - 3] = '3';
         ret_str[2][ret_nr - 2] = '\'';
         ret_pr_once = 0;
       }
       ret_nr++;
     }

     save_append_string(&ret_str[3], &ret_space, o, ret_para, _jmp_buf);
     save_append_string(&ret_str[3], &ret_space, o, ret_str[0], _jmp_buf);
     save_append_string(&ret_str[3], &ret_space, o, " 3\'\\n", _jmp_buf);
     save_append_string(&ret_str[3], &ret_space, o, ret_str[1], _jmp_buf);
     save_append_string(&ret_str[3], &ret_space, o, "\\n", _jmp_buf);
     save_append_string(&ret_str[3], &ret_space, o, ret_str[2], _jmp_buf);
     save_append_string(&ret_str[3], &ret_space, o, " 5\'\\n", _jmp_buf);

     ret_ptr = (char *) safe_malloc(strlen(ret_str[3]) + 1, _jmp_buf, o);
     strcpy(ret_ptr, ret_str[3]);
     if (ret_str[3]) {
       free(ret_str[3]);
     }
     free(ret_str[0]);
     free(ret_str[1]);
     free(ret_str[2]);
   }
   free(duplex[0]);
   free(duplex[1]);
   free(duplex[2]);
   free(duplex[3]);

   return ret_ptr;
}

char *
drawHairpin(int *bp, double mh, double ms, const thal_mode mode, double temp,
            const unsigned char *oligo1, const unsigned char *oligo2,
            double saltCorrection, int oligo1_len, int oligo2_len,
            jmp_buf _jmp_buf, thal_results *o)
{
   (void) oligo2;
   (void) oligo2_len;
   int  ret_space = 0;
   char *ret_ptr;
   int ret_last_l, ret_first_r, ret_center, ret_left_end, ret_right_start, ret_left_len, ret_right_len;
   int ret_add_sp_l, ret_add_sp_r;
   char ret_center_char;
   char ret_para[400];
   char* ret_str;
   /* Plain text */
   int i, N;
   ret_ptr = NULL;
   N = 0;
   double mg, t;
   if (!isFinite(ms) || !isFinite(mh)) {
      if((mode != THL_FAST) && (mode != THL_DEBUG_F)) {
        if (mode != THL_STRUCT) {
          printf("0\tdS = %g\tdH = %g\tinf\tinf\n", (double) ms,(double) mh);
#ifdef DEBUG
          fputs("No temperature could be calculated\n",stderr);
#endif
        }
      } else {
         o->temp = 0.0; /* lets use generalization here */
         strcpy(o->msg, "No predicted sec struc for given seq\n");
      }
   } else {
      if((mode != THL_FAST) && (mode != THL_DEBUG_F)) {
         for (i = 1; i < oligo1_len; ++i) {
            if(bp[i-1] > 0) N++;
         }
      } else {
         for (i = 1; i < oligo1_len; ++i) {
            if(bp[i-1] > 0) N++;
         }
      }
      t = (mh / (ms + (((N/2)-1) * saltCorrection))) - ABSOLUTE_ZERO;
      mg = mh - (temp * (ms + (((N/2)-1) * saltCorrection)));
      ms = ms + (((N/2)-1) * saltCorrection);
      o->dg = mg;
      o->ds = ms;
      o->dh = mh;
      o->temp = (double) t;
      if((mode != THL_FAST) && (mode != THL_DEBUG_F)) {
         if (mode != THL_STRUCT) {
           printf("Calculated thermodynamical parameters for hairpin:\tdS = %g\tdH = %g\tdG = %g\tt = %g\n",
                  (double) ms, (double) mh, (double) mg, (double) t);
         } else {
           snprintf(ret_para, 400, "Tm: %.1f&deg;C  dG: %.0f cal/mol  dH: %.0f cal/mol  dS: %.0f cal/mol*K\\n",
                   (double) t, (double) mg, (double) mh, (double) ms);
         }
      } else {
         return NULL;
      }
   }
   /* plain-text output */
   char* asciiRow;
   asciiRow = (char*) safe_malloc(oligo1_len, _jmp_buf, o);
   for(i = 0; i < oligo1_len; ++i) asciiRow[i] = '0';
   for(i = 1; i < oligo1_len+1; ++i) {
      if(bp[i-1] == 0) {
         asciiRow[(i-1)] = '-';
      } else {
         if(bp[i-1] > (i-1)) {
            asciiRow[(bp[i-1]-1)]='\\';
         } else  {
            asciiRow[(bp[i-1]-1)]='/';
         }
      }
   }
   if ((mode == THL_GENERAL) || (mode == THL_DEBUG)) {
     printf("SEQ\t");
     for(i = 0; i < oligo1_len; ++i) printf("%c",asciiRow[i]);
     printf("\nSTR\t%s\n", oligo1);
   }
   if (mode == THL_STRUCT) {
     ret_str = NULL;

     save_append_string(&ret_str, &ret_space, o, ret_para, _jmp_buf);

     ret_last_l = -1;
     ret_first_r = -1;
     ret_center_char = '|';
     for(i = 0; i < oligo1_len; ++i) {
       if (asciiRow[i] == '/') {
         ret_last_l = i;
       }
       if ((ret_first_r == -1) && (asciiRow[i] == '\\')) {
         ret_first_r = i;
       }
     }
     ret_center = ret_first_r - ret_last_l;
     if (ret_center % 2 == 0) {
       /* ret_center is odd */
       ret_left_end = ret_last_l + (ret_first_r - ret_last_l) / 2 - 1;
       ret_center_char = (char) oligo1[ret_left_end + 1];
       ret_right_start = ret_left_end + 2;
     } else {
       /* ret_center is even */
       ret_left_end = ret_last_l + (ret_first_r - ret_last_l - 1) / 2;
       ret_right_start = ret_left_end + 1;
     }
     ret_left_len = ret_left_end + 1;
     ret_right_len = oligo1_len - ret_right_start;
     ret_add_sp_l = 0;
     ret_add_sp_r = 0;
     if (ret_left_len > ret_right_len) {
       ret_add_sp_r = ret_left_len - ret_right_len + 1;
     }
     if (ret_right_len > ret_left_len) {
       ret_add_sp_l = ret_right_len - ret_left_len;
     }
     for (i = 0 ; i < ret_add_sp_l ; i++) {
       save_append_char(&ret_str, &ret_space, o, ' ', _jmp_buf);
     }
     save_append_string(&ret_str, &ret_space, o, "5' ", _jmp_buf);
     for (i = 0 ; i < ret_left_len ; i++) {
       save_append_char(&ret_str, &ret_space, o, (char) oligo1[i], _jmp_buf);
     }
     save_append_string(&ret_str, &ret_space, o, "U+2510\\n   ", _jmp_buf);
     for (i = 0 ; i < ret_add_sp_l ; i++) {
       save_append_char(&ret_str, &ret_space, o, ' ', _jmp_buf);
     }
     for (i = 0 ; i < ret_left_len ; i++) {
       if (asciiRow[i] == '/') {
         save_append_char(&ret_str, &ret_space, o, '|', _jmp_buf);
       } else {
         save_append_char(&ret_str, &ret_space, o, ' ', _jmp_buf);
       }
     }
     if (ret_center_char == '|' ) {
       save_append_string(&ret_str, &ret_space, o, "U+2502", _jmp_buf);
     } else {
       save_append_char(&ret_str, &ret_space, o, ret_center_char, _jmp_buf);
     }
     save_append_string(&ret_str, &ret_space, o, "\\n", _jmp_buf);
     for (i = 0 ; i < ret_add_sp_r - 1 ; i++) {
       save_append_char(&ret_str, &ret_space, o, ' ', _jmp_buf);
     }
     save_append_string(&ret_str, &ret_space, o, "3' ", _jmp_buf);
     for (i = oligo1_len ; i > ret_right_start - 1; i--) {
       save_append_char(&ret_str, &ret_space, o, (char) oligo1[i], _jmp_buf);
     }
     save_append_string(&ret_str, &ret_space, o, "U+2518\\n", _jmp_buf);

     ret_ptr = (char *) safe_malloc(strlen(ret_str) + 1, _jmp_buf, o);
     strcpy(ret_ptr, ret_str);
     if (ret_str != NULL) {
       free(ret_str);
     }
   }
   free(asciiRow);
   return ret_ptr;
}
