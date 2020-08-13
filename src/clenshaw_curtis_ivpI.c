/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Oct 2019
*  LAST MODIFIED:    Oct 2019
*  AFFILIATION:      Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA
*  DESCRIPTION:      Generates constant matrices for first order Clensahe & Curtis Quadrature
*
* INPUT:
*    N   -- Chebyshev polynomial order
*    M   -- Number of sample points
*
* OUTPUTS:
*    T1_1  -- Chebyshev Matrix [(M+1)x(N+1)]
*    P1_1  -- Picard Iteration Operator [(N+1)xN]
*    Ta_1  -- Chebyshev Matrix [(M+1)xN]
*    A_1   -- Least Squares Operator [Nx(M+1)]
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "lsq_chebyshev_fit.h"
#include "c_functions.h"
#include "chebyshev.h"
#include "const.h"

void clenshaw_curtis_ivpI( int N, int M, double* T1, double* P1, double* Ta, double* A ){

  // Least Squares Operator (A) [(N-1)x(M+1)]
  lsq_chebyshev_fit(-1.0,N-1,M,Ta,A);

  // Compute Constants of Integration (i.e. evaluated T at tau = -1).
  double Lconst[(N+1)*(N+1)];
  memset( Lconst, 0.0, ((N+1)*(N+1)*sizeof(double)));
  for (int k=0; k<=N; k++){
    Lconst[ID2(1,k+1,N+1)] = cos(k*acos(-1)); // Const of Integration (CoI)
  }

  // for (int i=0; i<=N; i++){
  //   for (int j=0; j<=N; j++){
  //     printf("Lcont %f\t",Lconst[ID2(i+1,j+1,N+1)]);
  //   }
  //   printf("\n");
  // }

  // S Matrix [(N+1)x(N+1)]
  double temp1[(N+1)*(N+1)];
  memset( temp1, 0.0, ((N+1)*(N+1)*sizeof(double)));
  temp1[ID2(1,1,N+1)] = 1.0;
  for (int i=1; i<N+1; i++){
    for (int j=1; j<N+1; j++){
      if (i == j){
        temp1[ID2(i+1,j+1,N+1)] = 1.0/(2.0*i);
      }
    }
  }

  // for (int i=0; i<=N; i++){
  //   for (int j=0; j<=N; j++){
  //     printf("temp1 %f\t",temp1[ID2(i+1,j+1,N+1)]);
  //   }
  //   printf("\n");
  // }

  double temp2[N*N];
  memset( temp2, 0.0, N*N*sizeof(double));
  for (int i=1; i<=N; i++){
    for (int j=1; j<=N; j++){
      if (i == j){
        temp2[ID2(i,j,N)] = 1.0;
      }
    }
  }

  // for (int i=1; i<=N; i++){
  //   for (int j=1; j<=N; j++){
  //     printf("temp2 %f\t",temp2[ID2(i,j,N)]);
  //   }
  //   printf("\n");
  // }

  double temp3[N*N];
  memset( temp3, 0.0, (N*N*sizeof(double)));
  for (int i=1; i<=N; i++){
    for (int j=3; j<=N; j++){
      if (i == j-2){
        temp3[ID2(i,j,N)] = -1.0;
      }
    }
  }

  // for (int i=1; i<=N; i++){
  //   for (int j=1; j<=N; j++){
  //     printf("temp3 %f\t",temp3[ID2(i,j,N)]);
  //   }
  //   printf("\n");
  // }

  double temp4[(N+1)*N];
  memset( temp4, 0.0, ((N+1)*N*sizeof(double)));
  for (int i=2; i<=N+1; i++){
    for (int j=1; j<=N; j++){
      temp4[ID2(i,j,N+1)] = temp2[ID2(i-1,j,N)] + temp3[ID2(i-1,j,N)];
    }
  }
  //
  // for (int i=1; i<=N+1; i++){
  //   for (int j=1; j<=N; j++){
  //     printf("temp4 %f\t",temp4[ID2(i,j,N+1)]);
  //   }
  //   printf("\n");
  // }

  double S[(N+1)*N];
  memset( S, 0.0, (N*N*sizeof(double)));
  matmul(temp1,temp4,S,N+1,N+1,N,N+1,N+1,N+1);
  S[ID2(1,1,N)] = 0.25;
  S[ID2(2,1,N)] = 1.0;


  // for (int i=1; i<=N+1; i++){
  //   for (int j=1; j<=N; j++){
  //     printf("S %f\t",S[ID2(i,j,N+1)]);
  //   }
  //   printf("\n");
  // }

  // Picard Integration Operator (first order system)
  double temp5[(N+1)*(N+1)];
  memset( temp5, 0.0, ((N+1)*(N+1)*sizeof(double)));
  for (int i=1; i<=N+1; i++){
    for (int j=1; j<=N+1; j++){
      temp5[ID2(i,j,N+1)] = -Lconst[ID2(i,j,N+1)];
      if (i == j){
        temp5[ID2(i,j,N+1)] = temp5[ID2(i,j,N+1)] + 1.0;
      }
    }
  }

  // for (int i=1; i<=N+1; i++){
  //   for (int j=1; j<=N+1; j++){
  //     printf("temp5 %f\t",temp5[ID2(i,j,N+1)]);
  //   }
  //   printf("\n");
  // }

  matmul(temp5,S,P1,N+1,N+1,N,N+1,N+1,N+1);  // [(N+1)xN]

  // for (int i=1; i<=N+1; i++){
  //   for (int j=1; j<=N; j++){
  //     printf("P1 %f\t",P1[ID2(i,j,N+1)]);
  //   }
  //   printf("\n");
  // }

  // // Chebyshev Matrix (interpolate "position" coefficients)
  chebyshev(-1.0,N,M,2,T1);    // arg4 = 2 -> Trig Cheby Poly

  // for (int i=1; i<=M+1; i++){
  //   for (int j=1; j<=N+1; j++){
  //     printf("T1 %f\t",T1[ID2(i,j,M+1)]);
  //   }
  //   printf("\n");
  // }

}
