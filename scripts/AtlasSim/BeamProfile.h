//Date: 04.10.12

/** @file
 * Implements BeamProfile class to store beam properties and get beam density
 */

#include <iostream>

#include "TMath.h"
#include "TVectorD.h"
#include "TMatrixD.h"

#include "TF2.h"
#include "TPad.h"

/** Implements beam density function.
 * Coordinates: x,y,z,s where s=ct parametrizes the time.
 * Effects considered:
 *  - full covariance matrix for transverse beam size (x,y)
 *  - crossing angles (x-z, y-z)
 *  - beam displacement in the transverse plane (x,y)
 *  - hourglass effect on transverse beam size
 */
class BeamProfile {
 public:
  /** Type of beam: 1 or 2.
   * Determines the direction of the evolution.
   */
  enum BeamType {
    BT_1,
    BT_2
  };
  /** Alias for coordinate indexes. */
  enum Coord {
    cX=0, ///< space: x
    cY=1, ///< space: y
    cZ=2, ///< space: z
    cS=3, ///< time: s=ct
    cNCoord  ///< number of coordinates
  };
 protected:
  /// Beam evolution
  BeamType m_beamType;
  /// Full Covariance matrix 'a' in (x,y,z)
  TMatrixD m_cov;
  /// Full Covariance matrix 'b 'in (x,y,z)
  TMatrixD m_cov_b;
  /// Full Covariance matrix 'c 'in (x,y,z)
  TMatrixD m_cov_c;
  /// Weighting given to Gaussian 'a'
  double m_weighta;
  /// Weighting given to Gaussian 'b'
  double m_weightb;
  /// Inidividual beam sizes: x
  double m_sigma_x;
  double m_sigma_x_b;
  double m_sigma_x_c;
  /// Inidividual beam sizes: y
  double m_sigma_y;
  double m_sigma_y_b;
  double m_sigma_y_c;
  /// Inidividual beam sizes: z
  double m_sigma_z;
  double m_sigma_z_b;
  double m_sigma_z_c;
  /// Correlation coefficient
  double m_rhoxy;
  double m_rhoxy_b;
  double m_rhoxy_c;
  /// Half-crossing angle between beams: x-z
  double m_xingA_xz;
  /// Half-crossing angle between beams: y-z
  double m_xingA_yz;
  /// Beam transverse separation in horizontal plane
  double m_bsep_x;
  /// Beam transverse separation in vertical plane
  double m_bsep_y;
  /// Beta star for horizontal (x) direction
  double m_betas_x;
  /// Beta star for vertical (y) direction
  double m_betas_y;

  double overlapval;

 public:
  BeamProfile(BeamType bt);
  ~BeamProfile();

  //==========
  // Settings
  //==========

  /// Print summary of settings.
  void Print();

  /// Expert: Set full covariance matrix 
  void SetCovMatrix(TMatrixD cov);
  void SetCovMatrix_b(TMatrixD cov);

  ///Set weight of Gaussian 'a'
  void SetWeightGausA(double weight);
  ///Set weight of Gaussian 'b'
  void SetWeightGausB(double weight);

  /** Set transverse part of the covariance matrix.
   * transverse beam sizes and correlation only (assumes uncorrelated longitudinal).
   */
  void SetBeamXY(double sigma_x, double sigma_y, double rho=0);
  void SetBeamXY_b(double sigma_x, double sigma_y, double rho=0);
  void SetBeamXY_c(double sigma_x, double sigma_y, double rho=0);

  /** Set longitudinal part of the covariance matrix.
   * Longitudinal beam size only (assumes uncorrelated transverse).
   */
  void SetBeamZ(double sigma_z);
  void SetBeamZ_b(double sigma_z);
  void SetBeamZ_c(double sigma_z);

  /** Set the covariance matrix assuming uncorrelated transverse and longitudinal components. */
  //Gaussian 'a'
  void SetBeamXYZ(double sigma_x, double sigma_y, double sigma_z,  double rhoxy=0);
  //Gaussian 'b'
  void SetBeamXYZ_b(double sigma_x, double sigma_y, double sigma_z, double rhoxy=0);

  /** Set crossing angle (half of the beam-beam angle) in x-z plane. */
  void SetXingAngleXZ(double phi);

  /** Set crossing angle (half of the beam-beam angle) in y-z plane. */
  void SetXingAngleYZ(double phi);

  /** Set beam offset (separation) along horizontal plane. */
  void SetBeamOffsetX(double bsep);

  /** Set beam offset (separation) along vertical plane. */
  void SetBeamOffsetY(double bsep);

  /** Set beta* for horizontal direction (x) */
  void SetBetaStarX(double bs);

  /** Set beta* for vertical direction (y) */
  void SetBetaStarY(double bs);

  /** Set the same beta* for horizontal and vertical directions */
  void SetBetaStar(double bs);


  void SetOverlapVal(double val);

  //==========
  // Accessors for beam density
  //==========

  /** Retrieve density function at a given position,time. */
  //total
  double rho(double x, double y, double z, double s);
  //Gaussian 'a' only
  double rhoa(double x, double y, double z, double s);
  //Gaussian 'b' only
  double rhob(double x, double y, double z, double s);
  //Gaussian 'c' only
  double rhoc(double x, double y, double z, double s);



  /** Retrieve transverse density function.
   * Assumes (x,y) and (z,s) to be uncorrelated.
   */
  //double rho_T(double x, double y);

  /** Retrieve longitudinal density function.
   * Assumes (x,y) and (z,s) to be uncorrelated.
   */
  //double rho_L(double z, double s);

  //==========
  // Accessors for beam properties
  //==========
  double GetRhoXY() {return m_rhoxy;};
  double GetRhoXYb() {return m_rhoxy_b;};
  double GetSigmaX() {return m_sigma_x;}
  double GetSigmaXb() {return m_sigma_x_b;}
  double GetSigmaY() {return m_sigma_y;}
  double GetSigmaYb() {return m_sigma_y_b;}
  double GetSigmaZ() {return m_sigma_z;}
  double GetSigmaZb() {return m_sigma_z_b;}

  double GetWeightN() {return m_weighta;}
  double GetWeightNb() {return m_weightb;}

  TMatrixD GetCovMatrix() {return m_cov;}
  TMatrixD GetCovMatrix_b() {return m_cov_b;}

  //==========
  // Display methods
  //==========
  /// Draw XY profile
  //  void DrawXY(double min, double max, double z0=0, double s0=0, TPad *p=0);

  // Provide TF2 for drawing beam profile in cX,cY at given par[0]=cZ, par[1]=cS
  double DrawXY_TF2(double *x, double *par);


};


using namespace std;
using namespace TMath;

BeamProfile::BeamProfile(BeamType bt)
{
  m_beamType = bt;
  m_cov.ResizeTo(3,3);
  m_cov_b.ResizeTo(3,3);
  m_cov_c.ResizeTo(3,3);
  for (unsigned int cr=0; cr < cS; cr++){
    for (unsigned int cc=0; cc < cS; cc++){
      m_cov[cr][cc] = 0.0;
      m_cov_b[cr][cc] = 0.0;
      m_cov_c[cr][cc] = 0.0;
    }
  }
  m_xingA_xz = m_xingA_yz = 0.0;
  m_bsep_x = m_bsep_y = 0.0;
  m_betas_x = m_betas_y = 0.0;  // beta*=0 has the special meaning of disabling the effect
  m_sigma_x = m_sigma_y = m_sigma_z = m_sigma_x_b = m_sigma_y_b = m_sigma_z_b = 0;
  m_rhoxy = m_rhoxy_b = m_rhoxy_c = 0;
}

BeamProfile::~BeamProfile()
{

}

void BeamProfile::Print()
{
  cout << "===== BEAM SETTINGS SUMMARY =====" << endl;
  cout << "Beam Type: ";
  switch (m_beamType) {
  case BT_1:
    cout << "Beam 1";
    break;
  case BT_2:
    cout << "Beam 2";
    break;
  default:
    cout << m_beamType;
  }
  cout << endl;

  cout << "Covariance matrix 'a':" << endl;
  m_cov.Print();
  cout << endl;

  cout << "Covariance matrix 'b':" << endl;
  m_cov_b.Print();
  cout << endl;

  cout << "Half-Crossing angle XZ, YZ: " << m_xingA_xz << ", " << m_xingA_yz << endl;
  cout << "Beam offset (separation) in X, Y: " << m_bsep_x << ", " << m_bsep_y << endl;
  cout << "Beta* in X, Y (0=Infinity): " << m_betas_x << ", " << m_betas_y << endl;

  cout << "=================================" << endl;
  cout << endl;
}


void BeamProfile::SetCovMatrix(TMatrixD cov)
{  
  m_cov = cov;
  //set individual sigmas accordingly
  m_sigma_x = TMath::Sqrt(m_cov[0][0]);
  m_sigma_y = TMath::Sqrt(m_cov[1][1]);
  m_sigma_z = TMath::Sqrt(m_cov[2][2]);
  m_rhoxy = m_cov[0][1] / m_sigma_x / m_sigma_y;
}

void BeamProfile::SetCovMatrix_b(TMatrixD cov)
{
  m_cov_b = cov;
  //set individual sigmas accordingly
  m_sigma_x_b = TMath::Sqrt(m_cov_b[0][0]);
  m_sigma_y_b = TMath::Sqrt(m_cov_b[1][1]);
  m_sigma_z_b = TMath::Sqrt(m_cov_b[2][2]);
  m_rhoxy_b = m_cov_b[0][1] / m_sigma_x_b / m_sigma_y_b;
}

void BeamProfile::SetWeightGausA(double weight)
{
  m_weighta = weight;
}

void BeamProfile::SetWeightGausB(double weight)
{
  m_weightb = weight;
}

void BeamProfile::SetBeamXY(double sigma_x, double sigma_y, double rho)
{
  m_cov[cX][cX] = sigma_x*sigma_x;
  m_cov[cX][cY] = rho*sigma_x*sigma_y;
  m_cov[cY][cX] = rho*sigma_x*sigma_y;
  m_cov[cY][cY] = sigma_y*sigma_y;
  m_sigma_x = sigma_x;
  m_sigma_y = sigma_y;
  m_rhoxy = rho;
}

void BeamProfile::SetBeamXY_b(double sigma_x, double sigma_y, double rho)
{
  m_cov_b[cX][cX] = sigma_x*sigma_x;
  m_cov_b[cX][cY] = rho*sigma_x*sigma_y;
  m_cov_b[cY][cX] = rho*sigma_x*sigma_y;
  m_cov_b[cY][cY] = sigma_y*sigma_y;

  m_sigma_x_b = sigma_x;
  m_sigma_y_b = sigma_y;

  m_rhoxy_b = rho;
}

void BeamProfile::SetBeamXY_c(double sigma_x, double sigma_y, double rho)
{
  m_cov_c[cX][cX] = sigma_x*sigma_x;
  m_cov_c[cX][cY] = rho*sigma_x*sigma_y;
  m_cov_c[cY][cX] = rho*sigma_x*sigma_y;
  m_cov_c[cY][cY] = sigma_y*sigma_y;
  m_sigma_x_c = sigma_x;
  m_sigma_y_c = sigma_y;
  m_rhoxy_c = rho;
}

void BeamProfile::SetBeamZ(double sigma_z)
{
  m_cov[cZ][cZ] = sigma_z*sigma_z;
  m_sigma_z = sigma_z;
}

void BeamProfile::SetBeamZ_b(double sigma_z)
{
  m_cov_b[cZ][cZ] = sigma_z*sigma_z;
  m_sigma_z_b = sigma_z;
}

void BeamProfile::SetBeamZ_c(double sigma_z)
{
  m_cov_c[cZ][cZ] = sigma_z*sigma_z;
  m_sigma_z_c = sigma_z;
}

void BeamProfile::SetXingAngleXZ(double phi)
{
  m_xingA_xz = phi;
}

void BeamProfile::SetXingAngleYZ(double phi)
{
  m_xingA_yz = phi;
}

void BeamProfile::SetBeamOffsetX(double bsep)
{
  m_bsep_x = bsep;
}

void BeamProfile::SetBeamOffsetY(double bsep)
{
  m_bsep_y = bsep;
}


void BeamProfile::SetBetaStarX(double bs)
{
  m_betas_x = bs;
}

void BeamProfile::SetBetaStarY(double bs)
{
  m_betas_y = bs;
}

void BeamProfile::SetBetaStar(double bs)
{
  m_betas_x = m_betas_y = bs;
}

void BeamProfile::SetBeamXYZ(double sigma_x, double sigma_y, double sigma_z, double rhoxy)
{
  SetBeamXY(sigma_x, sigma_y, rhoxy);
  SetBeamZ(sigma_z);
}
void BeamProfile::SetBeamXYZ_b(double sigma_x, double sigma_y, double sigma_z, double rhoxy)
{
  SetBeamXY_b(sigma_x, sigma_y, rhoxy);
  SetBeamZ_b(sigma_z);
}

void BeamProfile::SetOverlapVal(double val)
{
  overlapval = val;
}

double BeamProfile::rhoa(double x, double y, double z, double s)
{
  double sign=1;
  if (m_beamType == BT_1)
    sign=-1; //opposite to the z axis
  TVectorD p(3);
  p[cX] = m_bsep_x + x*Cos(m_xingA_xz) - sign*z*Sin(m_xingA_xz);
  p[cY] = m_bsep_y + y*Cos(m_xingA_yz) - sign*z*Sin(m_xingA_yz);
  p[cZ] = z*Cos(m_xingA_xz)*Cos(m_xingA_yz) + sign*x*Sin(m_xingA_xz) + sign*y*Sin(m_xingA_yz) - sign*s;

  TMatrixD Cov(m_cov);

  if (m_betas_x != 0 || m_betas_y != 0) {
    double betas_x = m_betas_x != 0 ? Sqrt(1 + Power(z/betas_x, 2)) : 1.0;
    double betas_y = m_betas_y != 0 ? Sqrt(1 + Power(z/betas_y, 2)) : 1.0;
    Cov[cX][cX] *= betas_x*betas_x;
    Cov[cX][cY] *= betas_x*betas_y;
    Cov[cY][cX] *= betas_y*betas_x;
    Cov[cY][cY] *= betas_y*betas_y;
  }
  double rho;
  double detCov = Cov.Determinant();
  if (detCov == 0) {
    cerr << "Covariance matrix is singular! Returning rho=0" << endl;
    return 0.0;
  }
  rho = 1./( Power(2*Pi(),3./2) * Sqrt(Abs(detCov)) ); //normalization
  rho = rho*Exp(-1./2 * p*(Cov.Invert()*p)); //Gaussian 

  return rho;
}

double BeamProfile::rhob(double x, double y, double z, double s)
{
  double sign=1;
  if (m_beamType == BT_1)
    sign=-1; //opposite to the z axis
  TVectorD p(3);
  p[cX] = m_bsep_x + x*Cos(m_xingA_xz) - sign*z*Sin(m_xingA_xz);
  p[cY] = m_bsep_y + y*Cos(m_xingA_yz) - sign*z*Sin(m_xingA_yz);
  p[cZ] = z*Cos(m_xingA_xz)*Cos(m_xingA_yz) + sign*x*Sin(m_xingA_xz) + sign*y*Sin(m_xingA_yz) - sign*s;

  TMatrixD Cov(m_cov_b);

  if (m_betas_x != 0 || m_betas_y != 0) {
    double betas_x = m_betas_x != 0 ? Sqrt(1 + Power(z/betas_x, 2)) : 1.0;
    double betas_y = m_betas_y != 0 ? Sqrt(1 + Power(z/betas_y, 2)) : 1.0;
    Cov[cX][cX] *= betas_x*betas_x;
    Cov[cX][cY] *= betas_x*betas_y;
    Cov[cY][cX] *= betas_y*betas_x;
    Cov[cY][cY] *= betas_y*betas_y;
  }
  double rho;
  double detCov = Cov.Determinant();
  if (detCov == 0) {
    cerr << "Covariance matrix is singular! Returning rho=0 (b)" << endl;
    return 0.0;
  }

  rho = 1./( Power(2*Pi(),3./2) * Sqrt(Abs(detCov)) ); //normalization
  rho = rho*Exp(-1./2 * p*(Cov.Invert()*p)); //Gaussian 

  return rho;
}

double BeamProfile::rhoc(double x, double y, double z, double s)
{
  double sign=1;
  if (m_beamType == BT_1)
    sign=-1; //opposite to the z axis
  TVectorD p(3);
  p[cX] = m_bsep_x + x*Cos(m_xingA_xz) - sign*z*Sin(m_xingA_xz);
  p[cY] = m_bsep_y + y*Cos(m_xingA_yz) - sign*z*Sin(m_xingA_yz);
  p[cZ] = z*Cos(m_xingA_xz)*Cos(m_xingA_yz) + sign*x*Sin(m_xingA_xz) + sign*y*Sin(m_xingA_yz) - sign*s;

  TMatrixD Cov(m_cov_c);

  if (m_betas_x != 0 || m_betas_y != 0) {
    double betas_x = m_betas_x != 0 ? Sqrt(1 + Power(z/betas_x, 2)) : 1.0;
    double betas_y = m_betas_y != 0 ? Sqrt(1 + Power(z/betas_y, 2)) : 1.0;
    Cov[cX][cX] *= betas_x*betas_x;
    Cov[cX][cY] *= betas_x*betas_y;
    Cov[cY][cX] *= betas_y*betas_x;
    Cov[cY][cY] *= betas_y*betas_y;
  }
  double rho;
  double detCov = Cov.Determinant();
  if (detCov == 0) {
    cerr << "Covariance matrix is singular! Returning rho=0" << endl;
    return 0.0;
  }

  rho = 1./( Power(2*Pi(),3./2) * Sqrt(Abs(detCov)) ); //normalization
  rho = rho*Exp(-1./2 * p*(Cov.Invert()*p)); //Gaussian 

  return rho;
}


//For triple Gaussian
double BeamProfile::rho(double x, double y, double z, double s)
{
  double rho;
  rho = (m_weighta)*rhoa(x,y,z,s) + (m_weightb)*rhob(x,y,z,s) + (1-m_weighta-m_weightb)*rhoc(x,y,z,s);
  return rho;
}

//For double Gaussian
/* double BeamProfile::rho(double x, double y, double z, double s) */
/* { */
/*   double rho; */
/*   rho = (m_weighta)*rhoa(x,y,z,s) + (1-m_weighta)*rhob(x,y,z,s); */
/*   return rho; */
/* } */


/*
  void BeamProfile::DrawXY(double min, double max, double z0, double s0, TPad *p)
  {
  TF2 * f = new TF2("f", (Double_t (*)(Double_t*, Double_t*))(this->DrawXY_TF2), min, max, min, max, 2);
  f->SetParameter(0, z0);
  f->SetParameter(1, s0);
  if (p) p->cd();
  p->SetGridx();
  p->SetGridy();
  f->Draw();
  }
*/

double BeamProfile::DrawXY_TF2(double *x, double *par)
{
  return rho(x[0], x[1], par[0], par[1]);
}
