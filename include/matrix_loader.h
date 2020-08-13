/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Oct 2019
*  LAST MODIFIED:    Oct 2019
*  AFFILIATION:      Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA
*  DESCRIPTION:      Header file
*/

#ifndef __ML__
#define __ML__

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include "const.h"

# define Nmin 10
# define Nmax 200

extern double arr_T1[(Nmax-Nmin+1)][(Nmax+1)*(Nmax+1)];
extern double arr_P1[(Nmax-Nmin+1)][(Nmax+1)*(Nmax+1)];
extern double arr_Ta[(Nmax-Nmin+1)][(Nmax+1)*(Nmax+1)];
extern double arr_A[(Nmax-Nmin+1)][(Nmax+1)*(Nmax+1)];

void matrix_loader();

#endif
