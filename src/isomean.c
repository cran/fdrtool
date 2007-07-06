/* isomean.c  (2007-07-06)  
 *
 * Copyright 2007 Korbinian Strimmer
 *
 * ported from R code originally by Kaspar Rufibach / June 2004
 *
 * This file is part of the `fdrtool' library for R and related languages.
 * It is made available under the terms of the GNU General Public
 * License, version 2, or at your option, any later version,
 * incorporated herein by reference.
 *
 * This program is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more
 * details. 
 *
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
 * MA 02111-1307, USA 
 */


#include <R.h>


/* 
 * input:  y       measured values in a regression setting
 *         w       weights
 *         n       length of y vector (> 1)
 * output: ghat    vector containing estimated (isotonic) values
 */
void C_isomean(double* y, double* w, int* n, double* ghat)
{
  int c, j, nn; 
  double neu;
  double* gew;
  int* k;
  
  nn = *n; /* nn > 1 */
 
  /* allocate vector - error handling is done by R */
  k = (int *) Calloc((size_t) nn, int);
  gew = (double *) Calloc((size_t) nn, double);

  c = 0;
  k[c] = 0;
  gew[c] = w[0];
  ghat[c] = y[0];

  for (j=1; j < nn; j++)
  {		
    c = c+1;
    k[c] = j;
    gew[c] = w[j];
    ghat[c] = y[j];

    /* c is at least 1 as nn is > 1 */
    while (ghat[c-1] >= ghat[c])
    {
      neu = gew[c]+gew[c-1];
      ghat[c-1] = ghat[c-1]+(gew[c]/neu)*(ghat[c]-ghat[c-1]);
      gew[c-1] = neu;
      c = c-1;

      if (c==0) break;
    }
  }

  while (nn >= 1)
  {
    for (j=k[c]; j < nn; j++)
    {
      ghat[j] = ghat[c];
    }
    nn = k[c];
    c = c-1;
  }

  /* free vector */
  Free(k); 
  Free(gew); 
}

