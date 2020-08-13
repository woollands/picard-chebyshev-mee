/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     May 2017
*  LAST MODIFIED:    May 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Generates Chebyshev polynomials of the first kind
*
* INPUT:
*    s   -- sign on tau (-1 or 1)
*    N   -- Chebyshev polynomial order
*    M   -- Number of sample points
*    arg -- Recursive OR Trigonometric formulation
*
* OUTPUTS:
*    T   -- Chebyshev polynomials
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "lsq_chebyshev_fit.h"
#include "c_functions.h"
#include "const.h"

void chebyshev(double s, int N, int M, int arg, double* T){

  // Cosine Sample Points
  double tau[(M+1)];
  memset( tau, 0.0, ((M+1)*sizeof(double)));
  for (int i=0; i<=M; i++){
    tau[i] = s*cos(i*C_PI/M);
  }

  if (arg == 1){
    // Chebyshev Polynomials (Recursive Formulation)
    for (int j=1; j<=N+1; j++){
      for (int i=1; i<=M+1; i++){
        if (j == 1){
          T[ID2(i,j,M+1)] = 1.0;
        }
        if (j == 2){
          T[ID2(i,j,M+1)] = tau[i-1];
        }
        if (j > 2){
          T[ID2(i,j,M+1)] = 2*tau[i-1]*T[ID2(i,j-1,M+1)] - T[ID2(i,j-2,M+1)];
        }
      }
    }
  }

  if (arg == 2){
    // Chebyshev Polyniomials (Trigonometric Formulation)
    for (int i=0; i<=M; i++){
      for (int j=0; j<=N; j++){
        T[ID2(i+1,j+1,M+1)] = cos(j*acos(tau[i]));
      }
    }
  }

}
