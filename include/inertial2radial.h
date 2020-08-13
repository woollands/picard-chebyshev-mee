/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Oct 2019
*  LAST MODIFIED:    Oct 2019
*  AFFILIATION:      Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA
*  DESCRIPTION:      Header file
*/

#ifndef __IR__
#define __IR__

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include "const.h"

void inertial2radial( double* r_eci, double* v_eci, double* xrdl, double* yrdl, double* zrdl );

#endif
