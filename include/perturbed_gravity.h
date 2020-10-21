/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Feb 2017
*  LAST MODIFIED:    Feb 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __PERT__
#define __PERT__

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include "const.h"

void perturbed_gravity(double t, double* Xo, double err, int ii, int N, double deg, double* G, double tol, int* itr, double* Feval);

void Grav_Approx(double t, double* X, double* dX, double* Feval);

void Grav_Full(double t, double* Xo, double* acc, double tol, double deg, double* Feval);

#endif
