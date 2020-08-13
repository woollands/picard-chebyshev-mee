/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Feb 2016
*  LAST MODIFIED:    Feb 2016
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Convert r & v to orbital elements
*
* INPUT:
*    r     -- Initial position vector (km)
*    v     -- Initial velocity vector (km/s)
*    tol   -- Integration Time (sec)
*
* OUTPUTS:
*    elm   -- Keplerian orbit elements (p,a,e,inc,Om,w,f,E,M,s)
*
*    WHERE: p  -- Semilatus Rectum (km)
*           a  -- Semimajor Axis (km)
*           e  -- Eccentricity
*           i  -- Inclination (rad)
*           Om -- Right Ascension of Ascending Node (rad)
*           w  -- Argument of Perigee (rad)
*           f  -- True Anomaly (rad)
*           E  -- Eccentric Anomaly (rad)
*           M  -- Mean Anomaly (rad)
*           s  -- Special case location of perigee (rad)
*                    -- Longitude of Perigee
*                    -- Argument of Latitude
*                    -- True Longitude
*
* REFERENCES:
* Vallado (p. 125 , Algorithm 9)
*
* COMMENTS:
*
*/

#include "rv2elm.h"
#include "adaptive_picard_chebyshev.h"
#include "c_functions.h"

void rv2elm(double* r, double* v, double tol, double* elm){

  // Position & Velocity Magnitudes
  double R, V;
  R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
  V = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

  // Angular Momentum Vector
  double H;
  double h[3] = {0.0};
  cross_product_3D(r,v,h);
  H = sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2]);

  // Line of Nodes Vector
  double nvec[3] = {0.0};
  double vert[3] = {0.0};
  double n;
  vert[0] = 0.0;
  vert[1] = 0.0;
  vert[2] = 1.0;
  cross_product_3D(vert,h,nvec);
  n = sqrt(nvec[0]*nvec[0] + nvec[1]*nvec[1] + nvec[2]*nvec[2]);

  // Eccentricity Vector
  double evec[3] = {0.0};
  double rv, e;
  rv = r[0]*v[0] + r[1]*v[1] + r[2]*v[2];
  for (int i=0; i<=2; i++){
    evec[i] = 1/C_MU*((pow(V,2) - C_MU/R)*r[i] - rv*v[i]);
  }
  e = sqrt(evec[0]*evec[0] + evec[1]*evec[1] + evec[2]*evec[2]);

  // Energy
  double xi;
  xi      = (pow(V,2))/2 - C_MU/R;

  // Semimajor Axis (a) & Semillatus Rectum (p)
  double a, p;
  if (fabs(1-e) <= tol){
    a   = INFINITY;
    p   = (pow(H,2))/C_MU;
  }
  else if (fabs(1-e) > tol){
    a   = -C_MU/2/xi;
    p   = a * (1 - pow(e,2));
  }

  // Inclination
  double inc;
  inc     = acos(h[2]/H);

  // Right Ascension of Ascending Node
  double Om;
  if (fabs(inc) >= tol){
    Om      = acos(nvec[0]/n);
    if (nvec[2] < 0){
      Om = 2*C_PI - Om;
    }
  } else if (fabs(inc) < tol){
    Om = 0.0;
  }

  // Argument of Perigee
  double w;
  if (fabs(inc) >= tol){
    w = acos((nvec[0]*evec[0] + nvec[1]*evec[1] + nvec[2]*evec[2]) / n / e);
    if (evec[2] < 0){
      w = 2*C_PI - w;
    }
  } else if (fabs(inc) < tol){
    w = 0.0;
  }

  // True Anomaly
  double f;
  double temp;
  temp = (evec[0]*r[0] + evec[1]*r[1] + evec[2]*r[2]) / R / e;
  if (fabs(temp-1.0) <= tol){
    f = 0.0;
  }
  else if (fabs(temp-1.0) > tol){
    f = acos(temp);
  }
  // f = acos((evec[0]*r[0] + evec[1]*r[1] + evec[2]*r[2]) / R / e);
  if (rv < 0){
    f = 2*C_PI - f;
  }

  // Mean Anomaly & Eccentric Anomaly
  double E, M;
  E = 2*atan2(sqrt(1-e)*tan(f/2),sqrt(1+e));
  if (E < 0){
    E = 2*C_PI + E;
  }
  M = E - e*sin(E);
  if (M < 0){
    M = 2*C_PI + M;
  }

  // Special Cases
  // Elliptical Equatorial (ascending node undefined)
  double s;
  if (inc < tol && e >= tol){
    s = acos(evec[0] / e);
    if (evec[1] < 0)
    s = 2*C_PI - s;   // Longitude of Perigee
  }
  // Circular Inclined (perigee undefined)
  else if (inc >= tol && e < tol){
    s = acos((nvec[0]*r[0] + nvec[1]*r[1] + nvec[2]*r[2])/R/n);    // Argument of Latitude
    if (r[2] < 0){
      s = 2*C_PI - s;
    }
    // Circular Equatorial (perigee & ascending node undefined)
    else if (inc < tol && e < tol){
      s = acos(r[0]/R);
      if (r[1] < 0){
        s = 2*C_PI - s;    // True Longitude
      }
    }
  }

// OUTPUT
elm[0] = p;
elm[1] = a;
elm[2] = e;
elm[3] = inc;
elm[4] = Om;
elm[5] = w;
elm[6] = f;
elm[7] = E;
elm[8] = M;
elm[9] = s;

}
