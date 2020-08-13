/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     March 2020
*  LAST MODIFIED:    March 2020
*  AFFILIATION:      Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA
*  DESCRIPTION:
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
#include "const.h"
#include "c_functions.h"

void mee2rv( double* tspan, double* states0, double* statesT, double* costates0,
  double Thr, double eta, double c, double rho, double P, double tol, double mps_tol){

    # Persistent TVEC TVECprev

    // Initialize
    double DEL       = 1.0e-5;
    int cont         = 0;
    double err       = 0.0;
    double tcntfinal = 8;

    // Set Initial Segment Scheme
    double TVEC
    double TVECprev;
    int seg_guess = 5;
    if (rho == 1){
      for (int i=1; i<=seg_guess; i++){
        TVEC[i] = i*tspan[1]/seg_guess;
        TVECprev[i] = TVEC[i];
      }
    }

    int N;
    int flag;
    double W1, W2;
    // MPS loop
    while(err > mps_tol && count < 15){
      for (int tcnt=1; tcnt <=tcntfinal; tcnt++){

        // Initial Condition Vector
        for (int i=0; i<=13; i++){
          if (i<7){
            IC[i] = states0[i];
          }
          if (i>6){
            IC[i] = costates[i-7];
          }
        }
        if (tcnt > 1){
          IC[tcnt+6] = IC[tcnt+6]+DEL;
        }
      }

      N = 40; // Nodes per segment (need to adaptively tune this in future)

      flag = 0;
      if (count > 0 && tcnt == 1){
        TVECprev = TVEC;
        TVEC = TVECsave;

        if length(TVEC) == length(TVECprev){
          if (max(abs(TVEC-TVEC)) < 1.0e-10){
            flag = 1;
          }
        }

      }

      // Integrate Dynamics
      [time,traj,GAMMA] = mcpi_mee(flag,tcnt,count,TVEC,TVECprev,IC,N,Thr,eta,c,rho,tol);

      if (tcnt == 1){}
      // Compute Switch Function
      S   = switch_function(traj,eta,c);

      TVECsave = [];
      for (i = 2:length(S)){
        if (S(i-1) < 0 && S(i) > 0 || S(i-1) > 0 && S(i) < 0){}
        TVECsave = [TVECsave; time(i-1)+abs(time(i)-time(i-1))/2];
      }
    }
  }

  if (tcnt == 1){
    REF = traj;   % Reference trajectory
    REFtime = time;
  }
  if (tcnt > 1){
    Del(:,tcnt-1) = [traj(end,1:6) traj(end,14)]' - [REF(end,1:6) REF(end,14)]';
  }


}


}
