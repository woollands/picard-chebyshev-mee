/* const.h
AUTHOR:           Robyn Woollands (robyn.woollands@gmail.com)
DATE WRITTEN:     May 2016
AFFILIATION:      Department of Aerospace Engineering, Texas A & M University, College Station, TX
DESCRIPTION:      Constants for Moon
*/
#ifndef _CONSTANTS_
#define _CONSTANTS_
#define C_PI 3.1415926535897932      // Pi
#define C_MU 0.4902800456866000E+04  // Gravitational Constant [km^3/s^2]
#define C_MUCan 1                    // Gravitational Constant Canonical Units
#define C_omega 7292115.0e-011       // Angular Speed of Earth [rad/s]
#define C_Req 0.1738000000000000E+04 // Equatorial Radius of Moon [km]
#define DU C_Req
#define TU sqrt(pow(DU,3)/C_MU)
#endif
