/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Jun 2016
*  LAST MODIFIED:    Jun 2016
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Convert states from body frame to inertial frame
*
* INPUT:
*    t   -- time (s)
*    aB  -- ECEF Position (km)
*
* OUTPUTS:
*    acc -- ECEF Acceleration (km/s^2)
*/

#include "ecef2eci.h"

void ecef2eci (double t,double* aB, double* acc){

  double th       = t*C_omega;
  double cos_th   = cos(th);
  double sin_th   = sin(th);

  // Convert to Acceleration to Inertial Frame
  acc[0]     = cos_th*aB[0] - sin_th*aB[1];
  acc[1]     = sin_th*aB[0] + cos_th*aB[1];
  acc[2]     = aB[2];

  return;
}
