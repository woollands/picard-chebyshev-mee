/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Oct 2019
*  LAST MODIFIED:    Oct 2019
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

void inertial2radial( double* r_eci, double* v_eci, double* xrdl, double* yrdl, double* zrdl ){

  // Radial in x- and z- direction
  double hvec[3] = {0.0};
  cross_product_3D( r_eci, v_eci, hvec );

  double h = 0.0;
  h = sqrt(pow(hvec[0],2)+pow(hvec[1],2)+pow(hvec[2],2));

  for (int i=0; i<=2; i++){
    xrdl[i] = r_eci[i]/sqrt(pow(r_eci[0],2)+pow(r_eci[1],2)+pow(r_eci[2],2));
    zrdl[i] = hvec[i]/h;
  }

  // Radial in y-direction
  cross_product_3D(zrdl,xrdl,yrdl);

}
