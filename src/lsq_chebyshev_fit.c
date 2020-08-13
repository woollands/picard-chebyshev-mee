/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     May 2017
*  LAST MODIFIED:    May 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Builds Least Sqaures Operator and Chebyshev Matrix
*
* INPUT:
*    s  -- sign on tau (-1 or 1)
*    N  -- Chebyshev polynomial order
*    M  -- Number of sample points
*
* OUTPUTS:
*    A  -- Least Squares Operator
*    T  -- Chebyshev matrix
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "lsq_chebyshev_fit.h"
#include "chebyshev.h"
#include "c_functions.h"

void lsq_chebyshev_fit(double s, int N, int M, double* T, double* A){

  // Generate Chebyshev Polyniomials
  chebyshev(s,N,M,2,T);

  // Weight Matrix
  double W[(M+1)*(M+1)];
  memset( W, 0.0, ((M+1)*(M+1)*sizeof(double)));
  for (int i=1; i<=M+1; i++){
    for (int j=1; j<=M+1; j++){
      if (i == j){
        W[ID2(i,j,M+1)] = 1.0;
      }
    }
  }
  W[ID2(1,1,M+1)] = 0.5;
  W[ID2(M+1,M+1,M+1)] = 0.5;

  // V Matrix
  double V[(N+1)*(N+1)];
  memset( V, 0.0, ((N+1)*(N+1)*sizeof(double)));
  for (int i=1; i<=N+1; i++){
    for (int j=1; j<=N+1; j++){
      if (i == j){
        V[ID2(i,j,N+1)] = 2.0/M;
      }
    }
  }
  if (M == N){
    V[ID2(1,1,N+1)] = 1.0/M;
    V[ID2(N+1,N+1,N+1)] = 1.0/M;
  }
  if (M > N){
    V[ID2(1,1,N+1)] = 1.0/M;
  }

  // T Transpose
  double TT[(N+1)*(M+1)];
  memset( TT, 0.0, ((N+1)*(M+1)*sizeof(double)));
  for (int i=1; i<=N+1; i++){
    for (int j=1; j<=M+1; j++){
      TT[ID2(i,j,N+1)] = T[ID2(j,i,M+1)];
    }
  }

  // Least Squares Operator
  double TTW[(M+1)*(M+1)];
  memset( TTW, 0.0, ((M+1)*(M+1)*sizeof(double)));
  matmul(TT,W,TTW,N+1,M+1,M+1,N+1,M+1,N+1);
  matmul(V,TTW,A,N+1,N+1,M+1,N+1,N+1,N+1);

}
