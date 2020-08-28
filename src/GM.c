/*! \mainpage EGM2008.c Calculates gravity according to the EGM 2008 Spherical Harmonic Model
*  AUTHORS:          Austin Probe (abprobe88@gmail.com) and Brent Macomber (brentmacomber@gmail.com)
*  DATE WRITTEN:     October 2014
*  LAST MODIFIED:    May 2016
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Computes Spherical Harmonic gravity for Orbit propagation using MCPI
*
*
*/


/*!
@file EGM2008.c
*/

/*! \file EGM2008.c EGM2008 Function File
* This file contains the fuctions for the EGM2008 Gravity Spherical
* harmonic libriary
*/

#include "EGM2008.h"
//#include "EGM2008.cc"
#include "MoonGrav.cc"
//#include "eci2ecef.h"
//#include "const.hpp"

// Declare Needed Variables
/*!
* This is the allocation for the maximum size associated
* Legendre polynomical matrix
*/
double P [(Max_Degree+3)*(Max_Degree+3)] = {0.0};
/*!
* This is the allocation for the maximum size associated
* Legendre polynomical matrix scale factor
*/
double scaleFactor [(Max_Degree+3)*(Max_Degree+3)] = {0.0};

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
void matmulEGM(double* A, double* B, double* OUT, int m, int n, int q)
{
	for(int i=0;i<m; i++){
		for(int j=0;j<q;j++){
			double sum = {0.0};
			for(int j1=0;j1<n;j1++)
			sum += A[IDX2F(i+1,j1+1,m)]*B[IDX2F(j1+1,j+1,n)];
			OUT[IDX2F(i+1,j+1,m)] = sum;
		}
	}
}

/*!
* \brief Gravity Evaluation
* This is the function that evaluates the spherical harmonic
* series to provide acceleration
*
* \param[in] p 3 element position in ECEF
* \param[in] Deg Degree and order of the serries to be used
* \param[out] Gxyz Gravitational Acceloration Output
*/
void GM_B( double* p, double* Gxyz, int DEG)
{

	double r             = {0.0};
	double phic          = {0.0};
	double lambda        = {0.0};
	double slambda       = {0.0};
	double clambda       = {0.0};
	double x = {0.0};
	double y = {0.0};
	double z = {0.0};
	double smlambda[Max_Degree+1] = {0.0};
	double cmlambda[Max_Degree+1] = {0.0};

	// determine radius of this thread's node
	x = p[0];
	y = p[1];
	z = p[2];

	// Compute geocentric radius
	r = pow( x*x + y*y + z*z , 0.5 );
	// Compute geocentric latitude
	phic  = asin( z / r );
	// Compute lambda
	lambda  = atan2( y, x );
	while (lambda<0)
	lambda = lambda+2*C_PI;
	while (lambda>=2*C_PI)
	lambda = lambda-2*C_PI;


	slambda = sin(lambda);
	clambda = cos(lambda);
	smlambda[0] = 0.0;
	cmlambda[0] = 1.0;
	smlambda[1] = slambda;
	cmlambda[1] = clambda;

	for(int m=2;m<DEG+1;m++){
		smlambda[m] = 2.0*clambda*smlambda[m-1] - smlambda[m-2];
		cmlambda[m] = 2.0*clambda*cmlambda[m-1] - cmlambda[m-2];
	}


	loc_gravLegendre( phic, scaleFactor, P , DEG);

	loc_gravityPCPF( p, P, DEG, smlambda, cmlambda, r, scaleFactor, Gxyz );


}

/*!
* \brief Legendre Polyniomial Evaluation
* This is the function that computes the normalized associated
* legendre polynomials based on geocentric latitude
*
* \param[in] phi Geocentric latitude
* \param[in] Deg Degree and order of the serries to be used
* \param[out] P associated Legendre polynomial matrix
* \param[out] scaleFactor Legendre scale factor
*/
void loc_gravLegendre( double phi, double* scaleFactor, double* P, int DEG )
{

	int k, p;
	double cphi = {0.0};
	double sphi = {0.0};


	cphi = cos(0.5*C_PI - phi);
	sphi = sin(0.5*C_PI - phi);
	// Seeds for recursion formula
	P[IDX2F(1,1,Max_Degree+3)] = 1.0;            // n = 0, m = 0;
	scaleFactor[IDX2F(1,1, Max_Degree+3)] = 0.0;
	P[IDX2F(2,1, Max_Degree+3)] = sqrt(3.0)*cphi ; // n = 1, m = 0;
	scaleFactor[IDX2F(2,1,Max_Degree+3)]  = 1.0;
	P[IDX2F(2,2,Max_Degree+3)] = sqrt(3.0)*sphi; // n = 1, m = 1;
	scaleFactor[IDX2F(2,2,Max_Degree+3)] = 0.0;

	// // New Method
	//     int nn = 2;
	//     int mm = 0;
	//     double m, n;
	//     int limit = (DEG+3)*(DEG+4)/2;
	//     for (int counter = 3; counter <= limit;counter++){
	//         k = nn + 1;
	//         p = mm + 1;
	//         n = coefMatrix[counter][0];
	//         m = coefMatrix[counter][1];
	//
	//         switch ( (int)(coefMatrix[counter][4]) )
	//        {
	//           case 1:
	//             P[IDX2F(k,k,Max_Degree+3)] = sqrt(2*n+1.0)/sqrt(2.0*n)*sphi*P[IDX2F(k-1,k-1, Max_Degree+3)];
	//             scaleFactor[IDX2F(k,k,Max_Degree+3)] = 0.0;
	//              nn++;
	//              mm = 0;
	//              break;
	//           case 2:
	//             P[IDX2F(k,p,Max_Degree+3)] = (sqrt(2*n+1)/n)*(sqrt(2*n-1.0)*cphi*P[IDX2F(k-1,p,Max_Degree+3)] - (n-1)/sqrt(2*n-3)* P[IDX2F(k-2,p, Max_Degree+3)] );
	//             scaleFactor[IDX2F(k,p,Max_Degree+3)] = sqrt( (n+1)*(n)/2);
	//             mm++;
	//             break;
	//           case 3:
	//             P[IDX2F(k,p,Max_Degree+3)] = sqrt(2*n+1)/(sqrt(n+m)*sqrt(n-m))*(sqrt(2*n-1.0)*cphi*P[IDX2F(k-1,p,Max_Degree+3)] - sqrt(n+m-1.0)*sqrt(n-m-1)/sqrt(2*n-3)*P[IDX2F(k-2,p,Max_Degree+3)] );
	//             scaleFactor[IDX2F(k,p,Max_Degree+3)] = sqrt( (n+m+1)*(n-m));
	//             mm++;
	//             break;
	//
	//        }
	//
	//     }

	// Old Method
	for (int nn = 2; nn <= DEG+2;nn++){
		double n = (double)nn;
		k = nn + 1;
		for(int mm=0; mm<=n;mm++) {
			double m = (double)mm;
			p = mm + 1;
			// Compute normalized associated legendre polynomials, P, via recursion relations
			// Scale Factor needed for normalization of dUdphi partial derivative
			if (n == m){
				P[IDX2F(k,k,Max_Degree+3)] = sqrt(2*n+1.0)/sqrt(2.0*n)*sphi*P[IDX2F(k-1,k-1, Max_Degree+3)];
				scaleFactor[IDX2F(k,k,Max_Degree+3)] = 0.0;
			}
			else if (m == 0){
				P[IDX2F(k,p,Max_Degree+3)] = (sqrt(2*n+1)/n)*(sqrt(2*n-1.0)*cphi*P[IDX2F(k-1,p,Max_Degree+3)] - (n-1)/sqrt(2*n-3)* P[IDX2F(k-2,p, Max_Degree+3)] );
				scaleFactor[IDX2F(k,p,Max_Degree+3)] = sqrt( (n+1)*(n)/2);
			}
			else {
				P[IDX2F(k,p,Max_Degree+3)] = sqrt(2*n+1)/(sqrt(n+m)*sqrt(n-m))*(sqrt(2*n-1.0)*cphi*P[IDX2F(k-1,p,Max_Degree+3)] - sqrt(n+m-1.0)*sqrt(n-m-1)/sqrt(2*n-3)*P[IDX2F(k-2,p,Max_Degree+3)] );
				scaleFactor[IDX2F(k,p,Max_Degree+3)] = sqrt( (n+m+1)*(n-m));
			}
		}
	}
}

/*!
* \brief Internal Gravitational Acceloration Evaluation
* This is the function that computes the gravitational acceloration based on
* the associated Legendre polynomials and the state
*
* \param[in] p Position vector in ECEF
* \param[in] P associated Legendre polynomial matrix
* \param[in] scaleFactor Legendre scale factor
* \param[in] Deg Degree and order of the serries to be used
* \param[in] r Position vector in ECEF
* \param[in] smlambda Trigonometric function of longitude
* \param[in] smlambda Trigonometric function of longitude
* \param[out] Gxyz Gravitational Acceloration Output
*/
void loc_gravityPCPF( double* p, double* P, int DEG, double* smlambda, double* cmlambda, double r, double* scaleFactor, double* Gxyz )
{

	int k, j;
	double radu;
	double rRatio  = {0.0};
	double rRatio_n  = {0.0};

	// initialize summation of gravity in radial coordinates
	double dUdrSumN   = {1.0};
	double dUdphiSumN    = {0.0};
	double dUdlambdaSumN = {0.0};

	double dUdrSumM   = {0.0};
	double dUdphiSumM = {0.0};
	double dUdlambdaSumM = {0.0};

	double dUdr    = {0.0};
	double dUdphi  = {0.0};
	double dUdlambda = {0.0};



	double x = p[0];
	double y = p[1];
	double z = p[2];
	radu = r;
	rRatio = C_Req/radu;
	rRatio_n = rRatio;
	// summation of gravity in radial coordinates




	//    // New Method
	//        int nn = 2;
	//        int mm = 0;
	//        double m, n;
	//        dUdrSumM      = 0.0;
	//        dUdphiSumM    = 0.0;
	//        dUdlambdaSumM = 0.0;
	//        int limit = (DEG+1)*(DEG+2)/2;
	//        for (int counter = 3; counter <= limit;counter++){
	//            k = nn + 1;
	//            j = mm + 1;
	//            n = coefMatrix[counter][0];
	//            m = coefMatrix[counter][1];

	//            switch ( (int)(coefMatrix[counter][4]) )
	//           {
	//              case 1:
	//                dUdrSumM      = dUdrSumM + P[IDX2F(k,j,Max_Degree+3)] *(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
	//                dUdphiSumM    = dUdphiSumM + ((P[IDX2F(k,j+1,Max_Degree+3)]*scaleFactor[IDX2F(k,j,Max_Degree+3)]) - z/(sqrt(x*x + y*y))*m*P[IDX2F(k,j,Max_Degree+3)])*(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
	//                dUdlambdaSumM = dUdlambdaSumM + m*P[IDX2F(k,j,Max_Degree+3)]*(S[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] - C[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
	//                rRatio_n = rRatio_n*rRatio;
	//                dUdrSumN      = dUdrSumN      + dUdrSumM*rRatio_n*k;
	//                dUdphiSumN    = dUdphiSumN    + dUdphiSumM*rRatio_n;
	//                dUdlambdaSumN = dUdlambdaSumN + dUdlambdaSumM*rRatio_n;
	//                 nn++;
	//                 mm = 0;
	//                 dUdrSumM      = 0.0;
	//                 dUdphiSumM    = 0.0;
	//                 dUdlambdaSumM = 0.0;
	//                 break;
	//              default:
	//                dUdrSumM      = dUdrSumM + P[IDX2F(k,j,Max_Degree+3)] *(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
	//                dUdphiSumM    = dUdphiSumM + ((P[IDX2F(k,j+1,Max_Degree+3)]*scaleFactor[IDX2F(k,j,Max_Degree+3)]) - z/(sqrt(x*x + y*y))*m*P[IDX2F(k,j,Max_Degree+3)])*(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
	//                dUdlambdaSumM = dUdlambdaSumM + m*P[IDX2F(k,j,Max_Degree+3)]*(S[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] - C[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
	//                mm++;

	//           }

	//        }

	// // Old Method
	// for (int n = 2; n <= DEG; n++) {
	// 	k = n+1;
	// 	rRatio_n = rRatio_n*rRatio;
	// 	dUdrSumM      = 0.0;
	// 	dUdphiSumM    = 0.0;
	// 	dUdlambdaSumM = 0.0;
	// 	for (int m = 0; m <= n; m++){
	// 		j = m+1;
	// 		dUdrSumM      = dUdrSumM + P[IDX2F(k,j,Max_Degree+3)] *(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
	// 		dUdphiSumM    = dUdphiSumM + ((P[IDX2F(k,j+1,Max_Degree+3)]*scaleFactor[IDX2F(k,j,Max_Degree+3)]) - z/(sqrt(x*x + y*y))*m*P[IDX2F(k,j,Max_Degree+3)])*(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
	// 		dUdlambdaSumM = dUdlambdaSumM + m*P[IDX2F(k,j,Max_Degree+3)]*(S[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] - C[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
	// 	}
	// 	dUdrSumN      = dUdrSumN      + dUdrSumM*rRatio_n*k;
	// 	dUdphiSumN    = dUdphiSumN    + dUdphiSumM*rRatio_n;
	// 	dUdlambdaSumN = dUdlambdaSumN + dUdlambdaSumM*rRatio_n;
	// }

	// +1 to each index

	int mindx[] = {2,3,4,5,6};
	int nindx[] = {0,0,0,0,0};
	int nterms = sizeof(mindx)/sizeof(mindx[0]);
	int m;



	for (int z = 0; z < nterms; z++){
		k = mindx[z]+1;
		j = nindx[z]+1;
		m = j-1;



		rRatio_n 			= pow(rRatio,k-1);
		dUdrSumN 		 	= dUdrSumN + rRatio_n*(P[IDX2F(k,j,Max_Degree+3)] *(coeffMoon[IDXCOEF(k-1,j-1)][2]*cmlambda[j-1] + coeffMoon[IDXCOEF(k-1,j-1)][3]*smlambda[j-1]));
		dUdphiSumN    = dUdphiSumN + rRatio_n*((P[IDX2F(k,j+1,Max_Degree+3)]*scaleFactor[IDX2F(k,j,Max_Degree+3)]) - z/(sqrt(x*x + y*y))*m*P[IDX2F(k,j,Max_Degree+3)])*(coeffMoon[IDXCOEF(k-1,j-1)][2]*cmlambda[j-1] + coeffMoon[IDXCOEF(k-1,j-1)][3]*smlambda[j-1]);
		dUdlambdaSumN = dUdlambdaSumN + rRatio_n*(m*P[IDX2F(k,j,Max_Degree+3)]*(coeffMoon[IDXCOEF(k-1,j-1)][3]*cmlambda[j-1] - coeffMoon[IDXCOEF(k-1,j-1)][2]*smlambda[j-1]));
	}


	// gravity in spherical coordinates
	dUdr      = -C_MU/(radu*radu)*dUdrSumN ;
	dUdphi    =  C_MU/radu*dUdphiSumN ;
	dUdlambda =  C_MU/radu*dUdlambdaSumN ;

	//gravity in ECEF coordinates
	Gxyz[0] = ((1.0/radu)*dUdr - (z/(radu*radu*sqrt(x*x + y*y)))*dUdphi)*x - (dUdlambda/(x*x + y*y))*y;
	Gxyz[1] = ((1.0/radu)*dUdr - (z/(radu*radu*sqrt(x*x + y*y)))*dUdphi)*y + (dUdlambda/(x*x + y*y))*x;
	Gxyz[2] = (1.0/radu)*dUdr*z + ((sqrt(x*x + y*y))/(radu*radu))*dUdphi;

	// special case for poles
	/*
	atPole = abs(atan2(p[IDX2F(i+1,3,M)],sqrt(p[IDX2F(i+1,1,M)]^2 + p[IDX2F(i+1,2,M)]^2)))==C_PI/2;
	if any(atPole){
	gx(atPole) = 0;
	gy(atPole) = 0;
	gz(atPole) = (1./r(atPole)).*dUdr(atPole).*p((atPole),3);
}
*/

}

/*!
* \brief Internal Gravity Potential Evaluation
* This is the function that computes the gravitational acceloration based on
* the associated Legendre polynomials and the state
*
* \param[in] p Position vector in ECEF
* \param[in] P associated Legendre polynomial matrix
* \param[in] scaleFactor Legendre scale factor
* \param[in] Deg Degree and order of the serries to be used
* \param[in] r Position vector in ECEF
* \param[in] smlambda Trigonometric function of longitude
* \param[in] smlambda Trigonometric function of longitude
* \param[out] Pot Gravitational Potential Output
*/
// void loc_gravityPot( double* p, double* P, int DEG, double* smlambda, double* cmlambda, double r, double* scaleFactor, double* Pot )
// {
//
// 	int k, j;
// 	double radu;
// 	double rRatio  = {0.0};
// 	double rRatio_n  = {0.0};
//
// 	// initialize summation of gravity in radial coordinates
// 	double dUdrSumN   = {1.0};
// 	double dUdphiSumN    = {0.0};
// 	double dUdlambdaSumN = {0.0};
//
// 	double dUdrSumM   = {0.0};
// 	double dUdphiSumM = {0.0};
// 	double dUdlambdaSumM = {0.0};
//
// 	double dUdr  = {0.0};
// 	double dUdphi  = {0.0};
// 	double dUdlambda = {0.0};
// 	double USum = {0.0};
//
//
// 	double x = p[0];
// 	double y = p[1];
// 	double z = p[2];
// 	int n;
// 	int m;
//
// 	radu = r;
// 	rRatio = C_Req/radu;
// 	rRatio_n = rRatio;
// 	// summation of gravity in radial coordinates
//
// 	for (n = 2; n <= DEG; n++) {
// 		k = n+1;
// 		rRatio_n = rRatio_n*rRatio;
// 		dUdrSumM      = 0.0;
// 		dUdphiSumM    = 0.0;
// 		dUdlambdaSumM = 0.0;
// 		for (m = 0; m <= n; m++){
// 			j = m+1;
// 			dUdrSumM      = dUdrSumM + P[IDX2F(k,j,Max_Degree+3)] *(coeffMoon[k-2][j-2]*cmlambda[j-1] + coeffMoon[k-2][j-2]*smlambda[j-1]);
// 			dUdphiSumM    = dUdphiSumM + ((P[IDX2F(k,j+1,Max_Degree+3)]*scaleFactor[IDX2F(k,j,Max_Degree+3)]) - z/(sqrt(x*x + y*y))*m*P[IDX2F(k,j,Max_Degree+3)])*(coeffMoon[k-2][j-2]*cmlambda[j-1] + coeffMoon[k-2][j-2]*smlambda[j-1]);
// 			dUdlambdaSumM = dUdlambdaSumM + m*P[IDX2F(k,j,Max_Degree+3)]*(coeffMoon[k-2][j-2]*cmlambda[j-1] - coeffMoon[k-2][j-2]*smlambda[j-1]);
// 			// printf("UsumM: %e\n ",dUdrSumM);
// 		}
// 		dUdrSumN      = dUdrSumN      + dUdrSumM*rRatio_n;
// 		dUdphiSumN    = dUdphiSumN    + dUdphiSumM*rRatio_n;
// 		dUdlambdaSumN = dUdlambdaSumN + dUdlambdaSumM*rRatio_n;
// 	}
// 	// gravity in spherical coordinates
// 	dUdr      =  C_MU/(radu)*dUdrSumN ;
// 	dUdphi    =  C_MU/radu*dUdphiSumN ;
// 	dUdlambda =  C_MU/radu*dUdlambdaSumN ;
//
// 	// gravity in potential
// 	// *Pot = MU * sqrt(dUdr*dUdr + dUdphi*dUdphi + dUdlambda*dUdlambda)/radu;
// 	*Pot =  sqrt(dUdr*dUdr);
//
//
// }
