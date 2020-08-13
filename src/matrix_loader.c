/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com) (based on a similar code by Brent Macomber)
*  DATE WRITTEN:     May 2017
*  LAST MODIFIED:    May 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Loads constant matrices required for the Picard Chebyshev first order integration (MEEs)
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "matrix_loader.h"

// Initialize Arrays
double arr_T1[(Nmax-Nmin+1)][(Nmax+1)*(Nmax+1)];
double arr_P1[(Nmax-Nmin+1)][(Nmax+1)*(Nmax+1)];
double arr_Ta[(Nmax-Nmin+1)][(Nmax+1)*(Nmax+1)];
double arr_A[(Nmax-Nmin+1)][(Nmax+1)*(Nmax+1)];

void matrix_loader(){

  // Open Files
  FILE* fT1 = fopen("../matrices/T1_matrices.bin","rb");
  FILE* fP1 = fopen("../matrices/P1_matrices.bin","rb");
  FILE* fTa = fopen("../matrices/Ta_matrices.bin","rb");
  FILE* fA  = fopen("../matrices/A_matrices.bin","rb");

  // Confirm Opening
  if ( !fT1 ){
    printf("Failure to open fT1 for binary write: CHECK PATH\n");
  }
  if ( !fP1 ){
    printf("Failure to open fP1 for binary write: CHECK PATH\n");
  }
  if ( !fTa ){
    printf("Failure to open fTa for binary write: CHECK PATH\n");
  }
  if ( !fA ){
    printf("Failure to open fA for binary write: CHECK PATH\n");
  }

  // Read Binary Data
  int c1, c2, c3, c4;
  c1 = fread( arr_T1, sizeof(double), (Nmax-Nmin+1)*(Nmax+1)*(Nmax+1), fT1 );
  c2 = fread( arr_P1, sizeof(double), (Nmax-Nmin+1)*(Nmax+1)*(Nmax+1), fP1 );
  c3 = fread( arr_Ta, sizeof(double), (Nmax-Nmin+1)*(Nmax+1)*(Nmax+1), fTa );
  c4 = fread( arr_A, sizeof(double), (Nmax-Nmin+1)*(Nmax+1)*(Nmax+1), fA );

  // Check that correct amount of data was read
  if ( c1 != (Nmax-Nmin+1)*(Nmax+1)*(Nmax+1) || c2 != (Nmax-Nmin+1)*(Nmax+1)*(Nmax+1) \
  || c3 != (Nmax-Nmin+1)*(Nmax+1)*(Nmax+1) || c4 != (Nmax-Nmin+1)*(Nmax+1)*(Nmax+1)){

    printf("Files contained incorrect data. :-(\n");
  }

  // Close Files
  fclose( fT1 );
  fclose( fP1 );
  fclose( fTa );
  fclose( fA );

}
