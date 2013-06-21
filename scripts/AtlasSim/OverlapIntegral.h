//date: 23.11.12

/** @file
 * Implement various integration of the overlap of the two beams.
 */

#include "BeamProfile.h"

#include <cmath>
#include "Math/Functor.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/IntegratorOptions.h"


// Store beam properties
//----------------
BeamProfile *B1; ///< Define beam 1 properties
BeamProfile *B2; ///< Define beam 2 properties

// Store internal integrator and range of integration
//----------------
/** Define integrator */
ROOT::Math::IntegratorMultiDim *ig;

/** Define integrator options **/
ROOT::Math::IntegratorMultiDimOptions * ig_opt;

/** Define functor integrand */
ROOT::Math::Functor *integrandFunctor; ///< cached integrand
/** Define integrator 1D */
ROOT::Math::IntegratorOneDim *ig1d;
/** Define functor 1D integrand */
ROOT::Math::Functor1D *integrandFunctor1D; ///< cached integrand

/** Default lower limit for integral ranges for all coordinates */
double cfgIntegralRangeLow[BeamProfile::cNCoord];
/** Default upper limit for integral ranges for all coordinates */
double cfgIntegralRangeUp[BeamProfile::cNCoord];

/** Define lower limit for the range of integration */
double integralRangeLow[BeamProfile::cNCoord];
/** Define upper limit for the range of integration */
double integralRangeUp[BeamProfile::cNCoord];


// Integration over all spatial coordinates
//----------------
/** Define generic overlap integrand as function of all 4 coordinates. */
double OverlapIntegrand(const double x[BeamProfile::cNCoord]) {
  return B1->rho(x[BeamProfile::cX],x[BeamProfile::cY],x[BeamProfile::cZ],x[BeamProfile::cS])*
    B2->rho(x[BeamProfile::cX],x[BeamProfile::cY],x[BeamProfile::cZ],x[BeamProfile::cS]);
}

double OverlapIntegrand_meanx(const double x[BeamProfile::cNCoord]) {
  return B1->rho(x[BeamProfile::cX],x[BeamProfile::cY],x[BeamProfile::cZ],x[BeamProfile::cS])*
    B2->rho(x[BeamProfile::cX],x[BeamProfile::cY],x[BeamProfile::cZ],x[BeamProfile::cS]) * x[BeamProfile::cX];
}

double OverlapIntegrand_meanx2(const double x[BeamProfile::cNCoord]) {
  return B1->rho(x[BeamProfile::cX],x[BeamProfile::cY],x[BeamProfile::cZ],x[BeamProfile::cS])*
    B2->rho(x[BeamProfile::cX],x[BeamProfile::cY],x[BeamProfile::cZ],x[BeamProfile::cS]) * x[BeamProfile::cX] * x[BeamProfile::cX];
}

double OverlapIntegrand_meany(const double x[BeamProfile::cNCoord]) {
  return B1->rho(x[BeamProfile::cX],x[BeamProfile::cY],x[BeamProfile::cZ],x[BeamProfile::cS])*
    B2->rho(x[BeamProfile::cX],x[BeamProfile::cY],x[BeamProfile::cZ],x[BeamProfile::cS]) * x[BeamProfile::cY];
}

double OverlapIntegrand_meany2(const double x[BeamProfile::cNCoord]) {
  return B1->rho(x[BeamProfile::cX],x[BeamProfile::cY],x[BeamProfile::cZ],x[BeamProfile::cS])*
    B2->rho(x[BeamProfile::cX],x[BeamProfile::cY],x[BeamProfile::cZ],x[BeamProfile::cS]) * x[BeamProfile::cY] * x[BeamProfile::cY];
}

double OverlapIntegrand_meanz(const double x[BeamProfile::cNCoord]) {
  return B1->rho(x[BeamProfile::cX],x[BeamProfile::cY],x[BeamProfile::cZ],x[BeamProfile::cS])*
    B2->rho(x[BeamProfile::cX],x[BeamProfile::cY],x[BeamProfile::cZ],x[BeamProfile::cS]) * x[BeamProfile::cZ];
}

double OverlapIntegrand_meanz2(const double x[BeamProfile::cNCoord]) {
  return B1->rho(x[BeamProfile::cX],x[BeamProfile::cY],x[BeamProfile::cZ],x[BeamProfile::cS])*
    B2->rho(x[BeamProfile::cX],x[BeamProfile::cY],x[BeamProfile::cZ],x[BeamProfile::cS]) * x[BeamProfile::cZ] * x[BeamProfile::cZ];
}


// Integration over all but one spatial coordinates
//----------------
BeamProfile::Coord nm1_coord; ///< N-1 coordinate to do not integrate
double nm1_value; ///< Value of coordinate to not be integrated

/** Get Overlap Integrand value for N-1 coordinates 
 * Set nm1_coord and nm1_value accordingly before calling this function.
 * @param x is the value of the other coordinates (exlcuding nm1_coord)
 */
double OverlapIntegrandNm1(const double *x) {
  double newX[BeamProfile::cNCoord];
  unsigned int idx=0;
  for (int nc=0; nc<BeamProfile::cNCoord; nc++) 
    if (nc == nm1_coord) {
      newX[nc] = nm1_value;
    } else {
      newX[nc] = x[idx];
      idx++;
    }
  return OverlapIntegrand(newX);
}

// Integration over time coordinate
//----------------
double oi_time_x[3]; ///< Values of (x,y,z)
/** Implement integrand as function of BeamProfile::cS. */
double OverlapIntegrandTime(const double x)
{
  return B1->rho(oi_time_x[BeamProfile::cX],oi_time_x[BeamProfile::cY],oi_time_x[BeamProfile::cZ],x)*
    B2->rho(oi_time_x[BeamProfile::cX],oi_time_x[BeamProfile::cY],oi_time_x[BeamProfile::cZ],x);

  //asymmetry beam profiles - if (oi_time_x[0]-GetBeamOffsetX > 0) return one beam profile / otherwise the other
  // perhaps to begin with you would only have asymmetry in one dimension, there would then be only
  // four combinations possible

  //or you could have a function which resets all of the beam properties based on that if statement, might be slightly easier to implement


}

// Public constructors and accessors
//----------------

void InitOverlapBeams() {
  //create new beams
  B1 = new BeamProfile(BeamProfile::BT_1);
  B2 = new BeamProfile(BeamProfile::BT_2);
  ig = new ROOT::Math::IntegratorMultiDim(); //by default: ROOT::Math::IntegrationMultiDim::kADAPTIVE 
  ig1d=0; //initialize when necessary
  
  integrandFunctor = 0;
  integrandFunctor1D = 0;

  ig_opt = 0;

  for (int idx=0; idx<BeamProfile::cNCoord; ++idx) {
    integralRangeUp[idx] = 0.0;
    integralRangeLow[idx] = 0.0;
  }
  

  //qwe

/*   cfgIntegralRangeLow[BeamProfile::cX] = -5.0; */
/*   cfgIntegralRangeUp[BeamProfile::cX] = 5.0; */
/*   cfgIntegralRangeLow[BeamProfile::cY] = -5.0; */
/*   cfgIntegralRangeUp[BeamProfile::cY] = 5.0; */

  cfgIntegralRangeLow[BeamProfile::cX] = -1.5;
  cfgIntegralRangeUp[BeamProfile::cX] =   1.5;
  cfgIntegralRangeLow[BeamProfile::cY] = -1.5;
  cfgIntegralRangeUp[BeamProfile::cY] =   1.5;


/*   cfgIntegralRangeLow[BeamProfile::cX] = -1.0; */
/*   cfgIntegralRangeUp[BeamProfile::cX] =   1.0; */
/*   cfgIntegralRangeLow[BeamProfile::cY] = -1.0; */
/*   cfgIntegralRangeUp[BeamProfile::cY] =   1.0; */


/*   Cfgintegralrangelow[Beamprofile::cX] = -20; */
/*   cfgIntegralRangeUp[BeamProfile::cX] = 20; */
/*   cfgIntegralRangeLow[BeamProfile::cY] = -20; */
/*   cfgIntegralRangeUp[BeamProfile::cY] = 20; */

/*   cfgIntegralRangeLow[BeamProfile::cX] = -30; */
/*   cfgIntegralRangeUp[BeamProfile::cX] = 30; */
/*   cfgIntegralRangeLow[BeamProfile::cY] = -30; */
/*   cfgIntegralRangeUp[BeamProfile::cY] = 30; */

  cfgIntegralRangeLow[BeamProfile::cZ] = -300.0;
  cfgIntegralRangeUp[BeamProfile::cZ] = 300.0;
  cfgIntegralRangeLow[BeamProfile::cS] = -1000.0;
  cfgIntegralRangeUp[BeamProfile::cS] = 1000.0;




}

void DestroyOverlapBeams() {
  delete B1;
  delete B2;
  delete ig;
  delete ig1d;
}


// Main integration routines
//----------------

/** Define overlap integral as function of a given coordinate and beam separation (x_i, h) 
 * Integrate B1*B2 in all coordinates but c.
 * @param c coordinate to do *not* integrate out
 * @param x coordinate 'c' value
 * @param ch coordinate (X or Y) to displace beams
 * @param h beam separation
 * @param fastProc do not load functor, assume the previous one (if one exists)
 * @return overlap integral for a given (x, h) point
 */
double OverlapIntegral_xh(BeamProfile::Coord c, double x, BeamProfile::Coord ch, double h, bool fastProc=false)
{
  if ( (not fastProc) or (integrandFunctor == 0) or (ig1d == 0)) {
    //init integrand, set N-1 integrand parameters
    if (integrandFunctor)
      delete integrandFunctor;
    nm1_coord = c;
    integrandFunctor = new ROOT::Math::Functor(&OverlapIntegrandNm1, (unsigned int)BeamProfile::cNCoord-1);
    ig->SetFunction(*integrandFunctor);
    //define integration ranges
    unsigned int idx=0;
    for (int nc=0; nc<BeamProfile::cNCoord; nc++) {
      if (nc == c)
	continue;
      integralRangeUp[idx] = cfgIntegralRangeUp[nc];
      integralRangeLow[idx] = cfgIntegralRangeLow[nc];
      idx++;
    }
    for (int idx2=idx; idx2 < BeamProfile::cNCoord; idx2++)
      integralRangeUp[idx2] = integralRangeLow[idx2] = 0.0;
  } //end init functor and integral range

  nm1_value = x;
  if (ch == BeamProfile::cX) {
    B1->SetBeamOffsetX(h/2);
    B2->SetBeamOffsetX(-h/2);
  } else if (ch == BeamProfile::cY) {
    B1->SetBeamOffsetY(h/2);
    B2->SetBeamOffsetY(-h/2);
  } else {
    cerr << "Not implemented to shift beam in the coordinate: " << ch << endl;
    throw;
  }

  return ig->Integral(integralRangeLow, integralRangeUp);

}

/** Define overlap integral as function of (x,y,z) and beam separation (h).
 * Integrate only on time BeamProfile::cS
 * @param ch coordinate (X or Y) to displace beams
 * @param h beam separation
 * @param fastProc do not load functor, assume the previous one (if one exists)
 * @param x contains (x, y, z) position
 * @param offset offset on the other transverse coordinate 
 * @return overlap integral for a given (x, y, z, h) point
 */
double OverlapIntegral_xyz_h(BeamProfile::Coord ch, double h, const double x[3], double offset=0, double driftx1=0, double driftx2=0, double drifty1=0, double drifty2=0, bool fastProc=false)
{
  if ( (not fastProc) or (integrandFunctor1D == 0)) {
    //init integrand
    if (integrandFunctor1D)
      delete integrandFunctor1D;
    if (ig1d)
      delete ig1d;
    integrandFunctor1D = new ROOT::Math::Functor1D(&OverlapIntegrandTime);
    ig1d = new ROOT::Math::Integrator(*integrandFunctor1D, ROOT::Math::IntegrationOneDim::kADAPTIVE,1.E-12,1.E-12);
    //define integration ranges
    integralRangeUp[0] = cfgIntegralRangeUp[BeamProfile::cS];
    integralRangeLow[0] = cfgIntegralRangeLow[BeamProfile::cS];
  }

  for (int i=0; i<3; i++)
    oi_time_x[i] = x[i];

  //beam positions are given in the ATLAS coordinate system
  //transverse offset in LHC system

  if (ch == BeamProfile::cX) {
    //left - verified
/*     B1->SetBeamOffsetX(h/2); */
/*     B2->SetBeamOffsetX(-h/2); */
 
//    cout << "Actual Separation = " << h + driftx2 + driftx1 << endl;

//Adding possibility for the beam to drift more than the quoted separation

    B1->SetBeamOffsetX(-h/2 -driftx1);
    //    cout << "beam1offset" << -h/2 -driftx1 << endl;  
    B2->SetBeamOffsetX(+h/2 +driftx2);
    //    cout << "beam2offset" << h/2 + driftx2 << endl;   

    //changed 
/*     B1->SetBeamOffsetY(-offset/2); */
/*     B2->SetBeamOffsetY(+offset/2); */

    B1->SetBeamOffsetY(-offset/2);
    B2->SetBeamOffsetY(+offset/2);


  } else if (ch == BeamProfile::cY) {
    //changed 
/*      B1->SetBeamOffsetY(-h/2);  */
/*      B2->SetBeamOffsetY(+h/2); */

    cout << "Actual Separation = " << h + drifty2 + drifty1 << endl;

     B1->SetBeamOffsetY(-h/2 -drifty1); 
     cout << "beam1offset" << -h/2 -drifty1 << endl;  
     B2->SetBeamOffsetY(+h/2 +drifty2);
     cout << "beam2offset" << h/2 + drifty2 << endl;   

    //left
/*     B1->SetBeamOffsetX(offset/2); */
/*     B2->SetBeamOffsetX(-offset/2); */

    B1->SetBeamOffsetX(+offset/2);
    B2->SetBeamOffsetX(-offset/2);


  } else {
    cerr << "Not implemented to shift beam in the coordinate: " << ch << endl;
    throw;
  }

  return ig1d->Integral(integralRangeLow[0], integralRangeUp[0]);
}

/** Define overlap integral as function of a given coordinate and beam separation (x_i, h) 
 * Integrate B1*B2 in all coordinates
 * @param fastProc do not load functor, assume the previous one (if one exists)
 * @return overlap integral for current B1,B2 beams
 */
double OverlapIntegral(bool fastProc=false) {
    if ( (not fastProc) or (integrandFunctor == 0) or (ig1d == 0)) {
      //init integrand, set N-1 integrand parameters
      if (integrandFunctor)
	delete integrandFunctor;
      integrandFunctor = new ROOT::Math::Functor(&OverlapIntegrand, (unsigned int)BeamProfile::cNCoord);
      //      ig_opt = new ROOT::Math::IntegratorMultiDimOptions(); 

/*       double test = pow(10,-300); */

/*           ig_opt->SetWKSize(100000000);  */
/*           ig_opt->SetDefaultAbsTolerance(test);  */
     
      ig->SetFunction(*integrandFunctor);

      //define integration ranges
      unsigned int idx=0;
      for (int nc=0; nc<BeamProfile::cNCoord; nc++) {
	integralRangeUp[idx] = cfgIntegralRangeUp[nc];
	integralRangeLow[idx] = cfgIntegralRangeLow[nc];
	idx++;
      }
    }
    
    //does the integration  

    return ig->Integral(integralRangeLow, integralRangeUp);

}

double OverlapIntegral_mean(BeamProfile::Coord ch, bool fastProc=false) {
    if ( (not fastProc) or (integrandFunctor == 0) or (ig1d == 0)) {
      //init integrand, set N-1 integrand parameters
      if (integrandFunctor)
	delete integrandFunctor;

      if (ch == BeamProfile::cX) integrandFunctor = new ROOT::Math::Functor(&OverlapIntegrand_meanx, (unsigned int)BeamProfile::cNCoord);     
      if (ch == BeamProfile::cY) integrandFunctor = new ROOT::Math::Functor(&OverlapIntegrand_meany, (unsigned int)BeamProfile::cNCoord);     
      if (ch == BeamProfile::cZ) integrandFunctor = new ROOT::Math::Functor(&OverlapIntegrand_meanz, (unsigned int)BeamProfile::cNCoord);     
      ig->SetFunction(*integrandFunctor);

      //define integration ranges
      unsigned int idx=0;
      for (int nc=0; nc<BeamProfile::cNCoord; nc++) {
	integralRangeUp[idx] = cfgIntegralRangeUp[nc];
	integralRangeLow[idx] = cfgIntegralRangeLow[nc];
	idx++;
      }
    }
    
    return ig->Integral(integralRangeLow, integralRangeUp);
}

double OverlapIntegral_x2(BeamProfile::Coord ch, bool fastProc=false) {
    if ( (not fastProc) or (integrandFunctor == 0) or (ig1d == 0)) {
      //init integrand, set N-1 integrand parameters
      if (integrandFunctor)
	delete integrandFunctor;

      if (ch == BeamProfile::cX) integrandFunctor = new ROOT::Math::Functor(&OverlapIntegrand_meanx2, (unsigned int)BeamProfile::cNCoord);     
      if (ch == BeamProfile::cY) integrandFunctor = new ROOT::Math::Functor(&OverlapIntegrand_meany2, (unsigned int)BeamProfile::cNCoord);     
      if (ch == BeamProfile::cZ) integrandFunctor = new ROOT::Math::Functor(&OverlapIntegrand_meanz2, (unsigned int)BeamProfile::cNCoord);     
      ig->SetFunction(*integrandFunctor);

      //define integration ranges
      unsigned int idx=0;
      for (int nc=0; nc<BeamProfile::cNCoord; nc++) {
	integralRangeUp[idx] = cfgIntegralRangeUp[nc];
	integralRangeLow[idx] = cfgIntegralRangeLow[nc];
	idx++;
      }
    }
    
    return ig->Integral(integralRangeLow, integralRangeUp);
}

