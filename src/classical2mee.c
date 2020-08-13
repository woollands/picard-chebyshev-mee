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

void classical2mee( double a, double e, double inc, double Om, double w, double nu, double* mee ){

  mee[0] = a*(1.0 - pow(e,2));
  mee[1] = e*cos(w + Om);
  mee[2] = e*sin(w + Om);
  mee[3] = tan(inc/2)*cos(Om);
  mee[4] = tan(inc/2)*sin(Om);
  mee[5] = Om + w + nu;

}
