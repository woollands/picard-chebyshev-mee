#ifndef __EGM__
#define __EGM__

/*! \mainpage EGM2008.h Calculates gravity according to the EGM 2008 Spherical Harmonic Model
 *  AUTHORS:          Austin Probe (abprobe88@gmail.com) and Brent Macomber (brentmacomber@gmail.com)
 *  DATE WRITTEN:     October 2014
 *  LAST MODIFIED:    May 2016
 *  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
 *  DESCRIPTION:      Computes Spherical Harmonic gravity for Orbit propagation using MCPI
 */




#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#include "const.h"

/* Macro for looking up indices using "1" based addressing ( a la fortran and matlab )
 and column-major array format [Note ld = Num of rows]
 */
#define IDX2F(i,j,ld) ((((j)-1)*(ld))+((i)-1))
#define IDX3F(i,j,k,ld1,ld2) ((((k)-1)*(ld1*ld2))+(((j)-1)*(ld1))+((i)-1))

#define Max_Degree 210

/**
  Round a / b to nearest higher integer value
 */
inline int iDivUp(int a, int b)
{
    return (a % b != 0) ? (a / b + 1) : (a / b);
}


/*!
 * \brief Matrix Multiplication
 * This is a simple matrix multiplication function
 *
 * \param[in] A Vector representation of matrix A
 * \param[in] B Vector representation of matrix B
 * \param[in] m Column dimension of A
 * \param[in] n Shared dimension of A and B
 * \param[in] q Row dimension of B
 * \param[out] OUT Matrix Output
 */
void matmulEGM(double* A, double* B, double* OUT, int m, int n, int q);


/*!
 * \brief Gravity Evaluation
 * This is the function that evaluates the spherical harmonic
 * series to provide acceloration
 *
 * \param[in] p 3-element position in ECEF
 * \param[in] Deg Degree and order of the series to be used
 * \param[in] delV_del_r_phi_lambda 9-element vector of gravity partials used to compute Jacobian fo$
 * \param[out] Gxyz Gravitational Acceloration Output
 */
void EGM2008( double* xECEF, double* aECEF, int DEG);


/*!
 * \brief Gravity Potential Evaluation
 * This is the function that evaluates the spherical harmonic
 * series to provide gravitational potential
 *
 * \param[in] p 3 element position in ECEF
 * \param[in] Deg Degree and order of the series to be used
 * \param[out] Pot Gravitational Potential Output
 */
void EGM2008Pot( double* xECEF, double* Pot, int DEG);


/*!
 * \brief Legendre Polyniomial Evaluation
 * This is the function that computes the normalized associated
 * legendre polynomials based on geocentric latitude
 *
 * \param[in] phi Geocentric latitude
 * \param[in] Deg Degree and order of the series to be used
 * \param[out] P associated Legendre polynomial matrix
 * \param[out] scaleFactor Legendre scale factor
 */
void loc_gravLegendre( double phi, double* scaleFactor, double* P, int DEG );


/*!
 * \brief Internal Gravitational Acceloration Evaluation
 * This is the function that computes the gravitational acceloration based on
 * the associated Legendre polynomials and the state
 *
 * \param[in] p Position vector in ECEF
 * \param[in] P associated Legendre polynomial matrix
 * \param[in] scaleFactor Legendre scale factor
 * \param[in] Deg Degree and order of the series to be used
 * \param[in] r Position vector in ECEF
 * \param[in] smlambda Trigonometric function of longitude
 * \param[in] smlambda Trigonometric function of longitude
 * \param[out] Gxyz Gravitational Acceloration Output
 */
void loc_gravityPCPF( double* xECEF, double* P, int DEG , double* smlambda, double* cmlambda, double r, double *scaleFactor, double* aECEF );


/*!
 * \brief Internal Gravity Potential Evaluation
 * This is the function that computes the gravitational acceloration based on
 * the associated Legendre polynomials and the state
 *
 * \param[in] p Position vector in ECEF
 * \param[in] P associated Legendre polynomial matrix
 * \param[in] scaleFactor Legendre scale factor
 * \param[in] Deg Degree and order of the series to be used
 * \param[in] r Position vector in ECEF
 * \param[in] smlambda Trigonometric function of longitude
 * \param[in] smlambda Trigonometric function of longitude
 * \param[out] Pot Gravitational Potential Output
 */
void loc_gravityPot( double* xECEF, double* P, int DEG , double* smlambda, double* cmlambda, double r, double *scaleFactor, double* Pot );


/*!
 * \brief Jacobi Integral
 * This is the function computes the Jacobi Integral based on position and
 * state vector
 *
 * This is a useful tool for checking the accuracy of conservitive
 * orbit propagation
 *
 * \param[in] solN State (position and velocity) vector in ECEF
 * \param[in] Deg Degree and order of the series to be used
 * \param[out] H Jacobi Integral Output
 */
void jacobiIntegral(double time, double* solN, double* H, int Deg);


#endif
