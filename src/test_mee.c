/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Oct 2019
*  LAST MODIFIED:    Oct 2019
*  AFFILIATION:      Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA
*  DESCRIPTION:      Propagate modified Equinoctial elements
*  REFERENCE:
*/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "const.h"
#include "classical2mee.h"
#include "matrix_loader.h"
#include "c_functions.h"
#include "mee2rv.h"
#include "ecef2eci.h"
#include "eci2ecef.h"
#include "EGM2008.h"
#include "perturbed_gravity.h"
#include "inertial2radial.h"

int main(){

  struct timeval start, end;
  gettimeofday(&start,NULL);

  int N       = 200;    // Chebyshev Order
  int M       = N;      // Number of Nodes
  double Deg  = 100.0;   // Gravity Degree

  double T1[(M+1)*(N+1)];     // [(M+1)x(N+1)]
  memset( T1, 0.0, ((M+1)*(N+1)*sizeof(double)));
  double P1[(N+1)*N];         // [(N+1)xN]
  memset( P1, 0.0, ((N+1)*N*sizeof(double)));
  double Ta[(M+1)*N];         // [(M+1)xN]
  memset( Ta, 0.0, ((M+1)*N*sizeof(double)));
  double A[N*(M+1)];          // [Nx(M+1)]
  memset( A, 0.0, (N*(M+1)*sizeof(double)));

  // Initial Conditions (elements)
  double a   = 8000.0;
  double e   = 0.1;
  double inc = 10.0*C_PI/180.0;
  double Om  = 0.0*C_PI/180.0;
  double w   = 0.0*C_PI/180.0;
  double nu  = 0.0*C_PI/180.0;
  double P   = 2.0*C_PI*sqrt(pow(a,3)/C_MU);

  // Classical to MEE
  double mee0[6] = {0.0};
  classical2mee(a,e,inc,Om,w,nu,mee0);

  // MCPI Initialization
  double orb  = 1.0;
  double t0   = 0.0;
  double tf   = orb*P;
  double W1   = (tf + t0)/2.0;
  double W2   = (tf - t0)/2.0;

  double tau[N+1];
  memset( tau, 0.0, ((N+1)*sizeof(double)));
  double times[N+1];
  memset( times, 0.0, ((N+1)*sizeof(double)));
  for (int i=0; i<=M; i++){
    tau[i] = -cos(i*C_PI/M);
    times[i]   = (W2*tau[i] + W1);
  }

  // LOAD PRECOMPUTED MATRICES
  matrix_loader();
  int idN;
  idN = (N-10);

  // Retrive Data from Storage Arrays
  double* temp1;
  double* temp2;
  double* temp3;
  double* temp4;
  temp1 = &arr_T1[idN][0];
  temp2 = &arr_P1[idN][0];
  temp3 = &arr_Ta[idN][0];
  temp4 = &arr_A[idN][0];

  // BUILD MATRICES
  for (int j=1; j<=M+1; j++){
    for (int k=1; k<=N+1; k++){
      T1[ID2(j,k,M+1)] = temp1[ID2(j,k,Nmax+1)];  // Chebyshev Matrix
    }
  }
  for (int j=1; j<=N+1; j++){
    for (int k=1; k<=N; k++){
      P1[ID2(j,k,N+1)] = temp2[ID2(j,k,Nmax+1)];  // Integration Operator
    }
  }
  for (int j=1; j<=M+1; j++){
    for (int k=1; k<=N; k++){
      Ta[ID2(j,k,M+1)] = temp3[ID2(j,k,Nmax+1)];  // Chebyshev Matrix
    }
  }
  for (int j=1; j<=N; j++){
    for (int k=1; k<=M+1; k++){
      A[ID2(j,k,N)] = temp4[ID2(j,k,Nmax+1)];  // Least Squares Operator
    }
  }

  double mee[(M+1)*6];
  memset( mee, 0.0, ((M+1)*6*sizeof(double)));
  // Warm Start
  for (int j=1; j<=6; j++){
    for (int i=1; i<=M+1; i++){
      if (j<=5){
        mee[ID2(i,j,M+1)] = mee0[j-1];
      }
      if (j==6){
        mee[ID2(i,6,M+1)] = orb*C_PI+orb*C_PI*tau[i-1];
      }
    }
  }

  // MCPI IVP
  double mcpi_tol = 1.0e-13;
  double Feval[2] = {0.0};
  int itr         = 0;
  double errX     = 1.0;

  double mee_new[(M+1)*6];
  memset( mee_new, 0.0, ((M+1)*6*sizeof(double)));
  double mee_dot[(M+1)*6];
  memset( mee_dot, 0.0, ((M+1)*6*sizeof(double)));

  // clock_t startTime = clock();

  while (errX > mcpi_tol){

    // clock_t startTime = clock();

    int i;
    // #pragma omp parallel for default(none) shared(W2,M,mee,mee_dot,Deg,times,errX,mcpi_tol,itr,Feval)
    //private(r_eci,v_eci,xECEF,vECEF,aECEF,aECI,TB,ad_eci,xrdl,yrdl,zrdl,a_r,a_t,a_h,cosL,sinL,sqpm,q,s_sq)
    for (i=1; i<=M+1; i++){

      double r_eci[3]  = {0.0};
      double v_eci[3]  = {0.0};
      double xECEF[3]  = {0.0};
      double vECEF[3]  = {0.0};
      double aECEF[3]  = {0.0};
      double aECI[3]   = {0.0};
      double TB[3]     = {0.0};
      double ad_eci[3] = {0.0};
      double xrdl[3]   = {0.0};
      double yrdl[3]   = {0.0};
      double zrdl[3]   = {0.0};
      double a_r  = 0.0;
      double a_t  = 0.0;
      double a_h  = 0.0;
      double cosL = 0.0;
      double sinL = 0.0;
      double sqpm = 0.0;
      double q    = 0.0;
      double s_sq = 0.0;

      // Convert MEE to Cartesian ECI
      mee2rv(mee[ID2(i,1,M+1)],mee[ID2(i,2,M+1)],mee[ID2(i,3,M+1)],mee[ID2(i,4,M+1)],mee[ID2(i,5,M+1)],mee[ID2(i,6,M+1)],r_eci,v_eci);

      // Convert from ECI to ECEF
      eci2ecef(times[i-1],r_eci,v_eci,xECEF,vECEF);

      // Compute SH gravity
      // EGM2008(xECEF, aECEF, Deg);
      perturbed_gravity(times[i-1],xECEF,errX,i,M,Deg,aECEF,mcpi_tol,&itr,Feval);

      // Convert from ECEF to ECI
      ecef2eci(times[i-1],aECEF,aECI);

      // Compute Two-Body Acceleration
      for (int j=0; j<=2; j++){
        TB[j] = -C_MU*r_eci[j]/pow(sqrt(pow(r_eci[0],2)+pow(r_eci[1],2)+pow(r_eci[2],2)),3);
      }

      // Compute Perturbed Acceleration
      for (int j=0; j<=2; j++){
        ad_eci[j] = aECI[j] - TB[j];
      }

      // Acceleration Components
      inertial2radial( r_eci, v_eci, xrdl, yrdl, zrdl );
      a_r = ad_eci[0]*xrdl[0]+ad_eci[1]*xrdl[1]+ad_eci[2]*xrdl[2];
      a_t = ad_eci[0]*yrdl[0]+ad_eci[1]*yrdl[1]+ad_eci[2]*yrdl[2];
      a_h = ad_eci[0]*zrdl[0]+ad_eci[1]*zrdl[1]+ad_eci[2]*zrdl[2];

      // Gauss' Variational Eqns (MEEs)
      cosL = cos(mee[ID2(i,6,M+1)]);
      sinL = sin(mee[ID2(i,6,M+1)]);
      sqpm = sqrt(mee[ID2(i,1,M+1)]/C_MU);
      q    = 1.0 + mee[ID2(i,2,M+1)]*cosL + mee[ID2(i,3,M+1)]*sinL;
      s_sq = 1.0 + pow(mee[ID2(i,4,M+1)],2) + pow(mee[ID2(i,5,M+1)],2);
      mee_dot[ID2(i,1,M+1)] = W2*sqpm*2.0*mee[ID2(i,1,M+1)]*a_t/q;
      mee_dot[ID2(i,2,M+1)] = W2*sqpm*(a_r*sinL + (((q+1.0)*cosL + mee[ID2(i,2,M+1)])*a_t)/q - (mee[ID2(i,3,M+1)]*(mee[ID2(i,4,M+1)]*sinL - mee[ID2(i,5,M+1)]*cosL)*a_h)/q);
      mee_dot[ID2(i,3,M+1)] = W2*sqpm*(-a_r*cosL + (((q+1.0)*sinL + mee[ID2(i,3,M+1)])*a_t)/q + (mee[ID2(i,2,M+1)]*(mee[ID2(i,4,M+1)]*sinL - mee[ID2(i,5,M+1)]*cosL)*a_h)/q);
      mee_dot[ID2(i,4,M+1)] = W2*sqpm*s_sq/2.0/q*cosL*a_h;
      mee_dot[ID2(i,5,M+1)] = W2*sqpm*s_sq/2.0/q*sinL*a_h;
      mee_dot[ID2(i,6,M+1)] = W2*(sqrt(C_MU*mee[ID2(i,1,M+1)])*pow(q/mee[ID2(i,1,M+1)],2) + sqpm*(mee[ID2(i,4,M+1)]*sinL - mee[ID2(i,5,M+1)]*cosL)/q*a_h);

    } // M+1 loop

    // clock_t endTime = clock();
    // float elapsedTime = ((float) (endTime - startTime))/CLOCKS_PER_SEC/(1.0);
    // printf("Elapsed time: %f s\n",elapsedTime);

    // Integrate
    double gamma[N*(M+1)];
    memset( gamma, 0.0, (N*(M+1)*sizeof(double)));
    double beta[N*(M+1)];
    memset( beta, 0.0, (N*(M+1)*sizeof(double)));
    matmul(A,mee_dot,gamma,N,M+1,6,N,M+1,N);    // LSQ Coefficients
    matmul(P1,gamma,beta,N+1,N,6,N+1,N,N+1);    // Integration Operator
    for (int j=1; j<=6; j++){
      beta[ID2(1,j,M+1)] = beta[ID2(1,j,M+1)] + mee0[j-1];
    }
    matmul(T1,beta,mee_new,M+1,N+1,6,M+1,N+1,M+1);  // Evaluate dynamics

    // Check Convergence (Non-dimensionalize semimajor axis)
    double err = 0.0;
    double max_err = 0.0;
    for (int pp=1; pp<=M+1; pp++){
      for (int j=1; j<=6; j++){
        if (j == 1){
            err = fabs(mee_new[ID2(pp,j,M+1)]/DU-mee[ID2(pp,j,M+1)]/DU); // Nondimensionalize p
        }
        if (j > 1){
            err = fabs(mee_new[ID2(pp,j,M+1)]-mee[ID2(pp,j,M+1)]);
        }
        if (err > max_err){
          max_err = err;
          // printf("max_err %f\n",err);
        }
      }
    }
    errX = max_err;
    printf("err %1.16E\n",errX);

    // Update
    memcpy(mee,mee_new,(M+1)*6*sizeof(double));

    // Iteration Counter
    itr = itr + 1;
    // printf("iter %i\n",itr);
    if (itr > 20){
      break;
    }

  } // Picard-Chebyshev while-loop

  // clock_t endTime = clock();
  // float elapsedTime = ((float) (endTime - startTime))/CLOCKS_PER_SEC/(1.0);
  // printf("Elapsed time: %f s\n",elapsedTime);

  gettimeofday(&end,NULL);
  double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
  printf("time: %f\n",delta);

  for (int i=0; i<=5; i++){
      printf("%1.15f\t",mee[ID2(M+1,i,M+1)]);
  }
  printf("\n");

} // Main
