/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com) (based on a similar code by Brent Macomber)
*  DATE WRITTEN:     Oct 2019
*  LAST MODIFIED:    Oct 2019
*  AFFILIATION:      Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA
*  DESCRIPTION:      One time build & store constant matrices required for the Picard-Chebyshev first order integration (MEEs)
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "clenshaw_curtis_ivpI.h"
#include "c_functions.h"
#include "matrix_loader.h"

int N, M;

// Initialize Storage Arrays
double MAT_T1[(Nmax-Nmin+1)][(Nmax+1)*(Nmax+1)];
double MAT_P1[(Nmax-Nmin+1)][(Nmax+1)*(Nmax+1)];
double MAT_Ta[(Nmax-Nmin+1)][(Nmax+1)*(Nmax+1)];
double MAT_A[(Nmax-Nmin+1)][(Nmax+1)*(Nmax+1)];

int main(){

  printf("Building constant matrices...\n");

  for (int i=Nmin; i<=Nmax; i++){
    N = i;
    M = N;
    printf("N = %i\n",N);

    // Compute Clenshaw-Curtis Quadrature Constant Matrices
    double T1[(M+1)*(N+1)];
    memset( T1, 0.0, ((M+1)*(N+1)*sizeof(double)));
    double P1[(N+1)*N];
    memset( P1, 0.0, ((N+1)*N*sizeof(double)));
    double Ta[(M+1)*N];
    memset( Ta, 0.0, ((M+1)*N*sizeof(double)));
    double A[N*(M+1)];
    memset( A, 0.0, (N*(M+1)*sizeof(double)));
    clenshaw_curtis_ivpI(N,M,T1,P1,Ta,A);

    // Build & Store Arrays
    for (int j=1; j<=M+1; j++){
      for (int k=1; k<=N+1; k++){
        MAT_T1[i-Nmin][ID2(j,k,Nmax+1)] = T1[ID2(j,k,M+1)];
      }
    }
    for (int j=1; j<=N+1; j++){
      for (int k=1; k<=N; k++){
        MAT_P1[i-Nmin][ID2(j,k,Nmax+1)] = P1[ID2(j,k,N+1)];
      }
    }
    for (int j=1; j<=M+1; j++){
      for (int k=1; k<=N; k++){
        MAT_Ta[i-Nmin][ID2(j,k,Nmax+1)] = Ta[ID2(j,k,M+1)];
      }
    }
    for (int j=1; j<=N; j++){
      for (int k=1; k<=M+1; k++){
        MAT_A[i-Nmin][ID2(j,k,Nmax+1)] = A[ID2(j,k,N)];
      }
    }

  }

  // Open Files
  FILE* fT1 = fopen("../matrices/T1_matrices.bin","wb");
  FILE* fP1 = fopen("../matrices/P1_matrices.bin","wb");
  FILE* fTa = fopen("../matrices/Ta_matrices.bin","wb");
  FILE* fA  = fopen("../matrices/A_matrices.bin","wb");

  // Confirm Opening
  if ( !fT1 ){
    printf("Failure to open fT1 for binary write: CHECK PATH\n");
  }
  if ( !fT1 ){
    printf("Failure to open fT1 for binary write: CHECK PATH\n");
  }
  if ( !fTa ){
    printf("Failure to open fTa for binary write: CHECK PATH\n");
  }
  if ( !fA ){
    printf("Failure to open fA for binary write: CHECK PATH\n");
  }

  printf("Saving constant matrices... ");

  // Write Binary Data
  fwrite( MAT_T1, sizeof(double), (Nmax-Nmin+1)*(Nmax+1)*(Nmax+1), fT1 );
  fwrite( MAT_P1, sizeof(double), (Nmax-Nmin+1)*(Nmax+1)*(Nmax+1), fP1 );
  fwrite( MAT_Ta, sizeof(double), (Nmax-Nmin+1)*(Nmax+1)*(Nmax+1), fTa );
  fwrite( MAT_A, sizeof(double), (Nmax-Nmin+1)*(Nmax+1)*(Nmax+1), fA );

  // Close Files
  fclose( fT1 );
  fclose( fP1 );
  fclose( fTa );
  fclose( fA );

  printf("Complete!\n");

}
