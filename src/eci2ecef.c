/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Jun 2016
*  LAST MODIFIED:    Jun 2016
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Convert states from inertial frame to body frame
*
* INPUT:
*    t   -- time (s)
*    X   -- ECI Position (km)
*    V   -- ECI Velocity (km/s)
*
* OUTPUTS:
*    xB  -- ECEF Position (km)
*    vB  -- ECEF Velocity (km/s)
*/

#include "eci2ecef.h"
#include "const.h"

void eci2ecef(double t, double* X, double* V, double* xB, double* vB){

  double th       = t*C_omega;
  double cos_th   = cos(th);
  double sin_th   = sin(th);

  // Convert to Position to Rotating Frame
  xB[0]   =  cos_th*X[0] + sin_th*X[1];
  xB[1]   = -sin_th*X[0] + cos_th*X[1];
  xB[2]   = X[2];

  // Convert Velocity to the Rotating Frame
  if (V[0]+V[1]+V[2] != 0.0){
    vB[0]   = V[0] + C_omega*X[1];
    vB[1]   = V[1] - C_omega*X[0];
    vB[2]   = V[2];
  }
  if (V[0]+V[1]+V[2] == 0.0){
    vB[0] = 0.0;
    vB[1] = 0.0;
    vB[2] = 0.0;
  }

  return;
}
