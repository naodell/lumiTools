//date 04.10.12

/** @file
 * Implement studies on beam properties
 */

#include "OverlapIntegral.h"
#include "Math/Derivator.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

#include "TApplication.h"
#include "TF1.h"
#include "TF12.h"
#include "TF2.h"
#include "TF3.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TStyle.h"
#include "TVirtualFFT.h"
#include "TTree.h"
#include "TNtupleD.h"
#include "TBranch.h"

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sstream>
#include <fstream>

#include "AtlasStyle.C"
#include "AtlasLabels.C"

// Global Settings
//----------------
//const int debug=100;
const int debug=5;


// Study z movement with displacement h
//----------------

// Beam overlap density as function of z and h
BeamProfile::Coord bs_xyz_dh_bsep_coord = BeamProfile::cY;
double bs_L_xyz_h_par_h = 0.0; ///< Set parameter h
double bs_L_xyz_h_par_driftx1_gradient = 0.0; ///< Set x-drift of beam 1
double bs_L_xyz_h_par_driftx1 = 0.0; ///< Set x-drift of beam 1
double bs_L_xyz_h_par_drifty1_gradient = 0.0; ///< Set y-drift of beam 1
double bs_L_xyz_h_par_drifty1 = 0.0; ///< Set y-drift of beam 1

double bs_L_xyz_h_par_driftx2_gradient = 0.0; ///< Set x-drift of beam 2
double bs_L_xyz_h_par_driftx2 = 0.0; ///< Set x-drift of beam 2
double bs_L_xyz_h_par_drifty2_gradient = 0.0; ///< Set y-drift of beam 2
double bs_L_xyz_h_par_drifty2 = 0.0; ///< Set y-drift of beam 2

double bs_L_z_h(double *x, double *par)
{
    double val = OverlapIntegral_xh(BeamProfile::cX, x[0],bs_xyz_dh_bsep_coord, par[0],true);
    if (debug >= 11)
        cout << "Value = " << val << endl;
    return val;
}

void bs_dzdh() {
    cout << endl;
    cout << "********************************************************************************" << endl;
    cout << "STUDYING MOVEMENT OF BEAMSPOT <Z> AS FUNCTION OF BEAM SEPARATION IN Y" << endl;
    cout << "********************************************************************************" << endl;
    cout << endl;

    //hack beam properties (otherwise set by command-line arguments)

    B1->Print();
    B2->Print();

    //Scan as function of h
    const double min_h=0;
    const double max_h=0.1;
    const double step_h=0.03;

    TF1 *f_bs_Lzh = new TF1("f_bs_Lzh", bs_L_z_h, -300, 300,1);
    TGraphErrors *h_bs_L_z_h = new TGraphErrors((max_h-min_h)/step_h);
    h_bs_L_z_h->SetName("h_bs_L_z_h");
    h_bs_L_z_h->SetTitle("<z> as function of vertical beam separation;Beam Separation (mm); <z> (mm)");

    Int_t idx=0;
    for (double h=min_h; h<max_h; h+=step_h) {
        f_bs_Lzh->SetParameter(0, h);

        if (debug >= 5)
            cout << "Drawing h=" << h << "...";

        if (debug >= 10) {
            TCanvas *c = new TCanvas();
            f_bs_Lzh->Draw();
            cout << "Done. Press key to continue." << endl;
            c->WaitPrimitive();
        }
        //find maximum of L(z,h) in z
        double max = f_bs_Lzh->GetMaximumX();
        if (debug >= 5)
            cout << "Maximum <z> = " << max << endl;

        //Fill histogram
        h_bs_L_z_h->SetPoint(idx, h, max);
        h_bs_L_z_h->SetPointError(idx, 0.0, 600. / f_bs_Lzh->GetNpx() / 200); //half of points interval divided by 100 as estimation of the error
        idx++;
    }
    delete f_bs_Lzh;

    //Draw z(h), and fit for dz/dh. Expected to be a straight line. Get its slope
    TCanvas *c_bs_dzdh = new TCanvas("c_bs_dzdh", "Beamspot dz/dh");
    h_bs_L_z_h->Draw("AP");
    h_bs_L_z_h->Fit("pol1");

    //apply a straight line fit
    /*
       TF1 *fit_pol1 = h_bs_L_z_h->GetFunction("pol1");
       double slope=0.0;
       double slopeError=0.0;
       if (fit_pol1) {
       slope=fit_pol1->GetParameter(1);
       slopeError = fit_pol1->GetParError(1);  
       } else {
       cerr << "ERROR: Cannot retrieve fit function: pol1" << endl;
       }
       c_bs_dzdh->Update();
       cout << "RESULT:" << endl;
       cout << "-------------" << endl;
       cout << "dz/dh = " << slope << " +/- " << slopeError << endl;
     */

    //Save results
    TFile *outRootF = TFile::Open("bs_dzdh.root", "RECREATE");
    c_bs_dzdh->Write();
    h_bs_L_z_h->Write();
    outRootF->Close();


}


// Study x/y movement with displacement in y/x of h
//----------------

/** Set the transverse coordinate to do not integrate.
 * Beams will be displaced by the other one.
 */
const BeamProfile::Coord xy_coord = BeamProfile::cY;

const BeamProfile::Coord yx_coord = (xy_coord == BeamProfile::cX) ? BeamProfile::cY : BeamProfile::cX;

/// Beam overlap density as function of x/y and h
double bs_L_xy_h(double *x, double *par)
{  
    double val = OverlapIntegral_xh(xy_coord, x[0], yx_coord, par[0],true);
    if (debug >= 11)
        cout << "Value = " << val << endl;
    return val;
}

void bs_dxydh() {

    cout << endl;
    cout << "********************************************************************************" << endl;
    cout << "STUDYING MOVEMENT OF BEAMSPOT <X> (or <Y>) AS FUNCTION OF BEAM SEPARATION IN Y (X)" << endl;
    cout << "********************************************************************************" << endl;
    cout << endl;

    //hack beam properties (otherwise set by command-line arguments)
    //B1->SetBeamXYZ(0.05,0.05, 60.0);
    //B2->SetBeamXYZ(0.05,0.05, 60.0);

    B1->Print();
    B2->Print();

    //Scan as function of h
    const double min_h=0;
    const double max_h=0.1;
    const double step_h=0.01;
    TF1 *f_bs_Lxyh = new TF1("f_bs_Lxyh", bs_L_xy_h, -5, 5,1);
    TString title = " as function of orthogonal beam separation;(Orthogonal) Beam Separation (mm); <z> (mm)";
    if (xy_coord == BeamProfile::cX)
        title = TString("<x> ") + title;
    else
        title = TString("<y> ") + title;    
    TGraphErrors *h_bs_L_xy_h = new TGraphErrors((max_h-min_h)/step_h);
    h_bs_L_xy_h->SetName("h_bs_L_xy_h");
    h_bs_L_xy_h->SetTitle(title.Data());
    Int_t idx=0;
    for (double h=min_h; h<max_h; h+=step_h) {
        f_bs_Lxyh->SetParameter(0, h);

        if (debug >= 5)
            cout << "Drawing h=" << h << "...";

        if (debug >= 10) {
            TCanvas *c = new TCanvas();
            f_bs_Lxyh->Draw();
            cout << "Done. Press key to continue." << endl;
            c->WaitPrimitive();
        }
        //find maximum of L(x,h) in x (or x<->y)
        double max = f_bs_Lxyh->GetMaximumX();
        if (debug >= 5)
            cout << "Maximum <x> = " << max << endl;

        //Fill histogram
        h_bs_L_xy_h->SetPoint(idx, h, max);
        h_bs_L_xy_h->SetPointError(idx, 0.0, 10. / f_bs_Lxyh->GetNpx() / 200); //half of points interval divided by 100 as estimation of the error
        idx++;
    }
    delete f_bs_Lxyh;

    //Draw z(h), and fit for dz/dh. Expected to be a straight line. Get its slope
    TCanvas *c_bs_dxydh = new TCanvas("c_bs_dxydh", "Beamspot dx/dh");
    h_bs_L_xy_h->Draw("AP");

    /* 
       h_bs_L_xy_h->Fit("pol1");
       TF1 *fit_pol1 = h_bs_L_xy_h->GetFunction("pol1");
       double slope=fit_pol1->GetParameter(1);
       double slopeError = fit_pol1->GetParError(1);  
       c_bs_dxydh->Update();
       cout << "RESULT:" << endl;
       cout << "-------------" << endl;
       cout << "dx/dh = " << slope << " +/- " << slopeError << endl;
     */

    //Save results
    TFile *outRootF = TFile::Open("bs_dxydh.root", "RECREATE");
    c_bs_dxydh->Write();
    h_bs_L_xy_h->Write();
    outRootF->Close();


}


// Study x,y,z beamspot movements with displacement in y/x of h
//----------------
/// Set the transverse coordinate to displace beams
//BeamProfile::Coord bs_xyz_dh_bsep_coord = BeamProfile::cY;
//double bs_L_xyz_h_par_h = 0.0; ///< Set parameter h

double bs_L_xyz_h_par_min_h=0;
double bs_L_xyz_h_par_max_h=0.1;
double bs_L_xyz_h_par_step_h=0.01;

double bs_L_xyz_h_par_offset = 0; ///< External Offset on the other transverse coordinate



int num_gaus = 1; //number of gaussians 
int bs_def = 1; //definition of beam spot 0->min; 1->mean; 2->Sampling; 3->Manual
TString triple = "no"; //triple gaussian or not

//Output file name
TString outputfilename = "default";

//HACK: Put external offset
//double bs_L_xyz_h_par_offset = 0.3; ///< Offset on the other transverse coordinate

/// Beam overlap density as function of (x,y,z;h). Change sign to minimize instead of maximize


//my attempt
/////////////////////////////////////////////


//This is the variable that contains the position at which the convolved distribution is evaluated.
double x_convolution[3];

double bs_L_xyz_h(const double *x)
{
    //  return -OverlapIntegral_xyz_h(bs_xyz_dh_bsep_coord, bs_L_xyz_h_par_h, x, true);
    double val = OverlapIntegral_xyz_h(bs_xyz_dh_bsep_coord, bs_L_xyz_h_par_h, x, bs_L_xyz_h_par_offset,bs_L_xyz_h_par_driftx1,bs_L_xyz_h_par_driftx2, bs_L_xyz_h_par_drifty1, bs_L_xyz_h_par_drifty2, true);

    //Define the parameters of the smearing matrix (3D gaussian)
    TMatrixD Cov_c;
    Cov_c.ResizeTo(3,3);
    Cov_c[0][0] = 0.1;
    Cov_c[1][1] = 0.1;
    Cov_c[2][2] = 0.1;
    Cov_c[0][1] = 0;
    Cov_c[1][0] = 0;
    Cov_c[1][2] = 0;
    Cov_c[2][1] = 0;
    Cov_c[0][2] = 0;
    Cov_c[2][0] = 0;

    double detCov_c = Cov_c.Determinant();
    TVectorD p(3);
    TVectorD p2(3); 
    p2[0]=x_convolution[0];
    p2[1]=x_convolution[1];
    p2[2]=x_convolution[2];
    p[0]=x[0];
    p[1]=x[1];
    p[2]=x[2];

    double rho_c;
    rho_c = 1./( Power(2*Pi(),3./2) * Sqrt(Abs(detCov_c)) ); //normalization
    rho_c = rho_c*Exp(-1./2 * (p2-p)*((Cov_c.Invert()*(p2-p)))); //Gaussian 


    //Multiply value of output of overlap integral by the smearing function
    //  val *= rho_c;

    if (debug >= 100)
        cout << "Value @ (" << x[0] << ","<< x[1] << "," << x[2] << ") = " << val << endl;
    if (bs_def == 0) return -val;
    if (bs_def == 2) return -val;
    if (bs_def == 1) return val;

    if (bs_def > 2)  return val;


}

//Perform the integration in order to smear overlap integral

ROOT::Math::IntegratorMultiDim *ig_conv;
ROOT::Math::Functor *integrandFunctor_conv; ///< cached integrand

double bs_L_xyz_h_conv(const double *x)
{

    x_convolution[0] = x[0];
    x_convolution[1] = x[1];
    x_convolution[2] = x[2];

    if (integrandFunctor_conv)
        delete integrandFunctor_conv;

    ig_conv = new ROOT::Math::IntegratorMultiDim(); //by default: ROOT::Math::IntegrationMultiDim::kADAPTIVE 
    ig_conv->SetRelTolerance(0.001);
    //ig_conv->SetAbsTolerance(0.001);
    integrandFunctor_conv = new ROOT::Math::Functor(&bs_L_xyz_h, 3);
    ig_conv->SetFunction(*integrandFunctor_conv);
    double integralRangeUp1[3];
    double integralRangeLow1[3];

    integralRangeUp1[0] = 1.5;
    integralRangeLow1[0] = -1.5;
    integralRangeUp1[1] = 1.5;
    integralRangeLow1[1] = -1.5;
    integralRangeUp1[2] = 1.5;
    integralRangeLow1[2] = -1.5;

    //   if (bs_def == 0)  return (ig_conv->Integral(integralRangeLow1, integralRangeUp1))*-1;
    //   else return ig_conv->Integral(integralRangeLow1, integralRangeUp1);
    return ig_conv->Integral(integralRangeLow1, integralRangeUp1);
}


//Wrap of bs_L_xyz_h fixing z=0 for TF2 display 
double bs_L_xyz0_h_TF2(double *x, double *par)
{
    double xnew[3];
    xnew[0] = x[0]; xnew[1] = x[1]; 
    xnew[2] = par[0]; //fix z=par[0]
    return bs_L_xyz_h(xnew);
}

//Wrap of bs_L_xyz_h for TF3 display 
double bs_L_xyz0_h_TF3(double *x, double *par)
{
    return bs_L_xyz_h(x);
}


//Wrap of bs_L_xyz_h_conv fixing z=0 for TF2 display 
double bs_L_xyz0_h_TF2_conv(double *x, double *par)
{
    double xnew[3];
    xnew[0] = x[0]; xnew[1] = x[1]; 
    xnew[2] = par[0]; //fix z=par[0]
    return bs_L_xyz_h_conv(xnew);
}

//Analytic convolution assuming double gaussian beam profiles. At present does not take into account crossing angles
double overlap_num(const double *x)
{

    //4 terms for (narrow1*narrow2), (broad1*narrow2), (narrow1*broad2), (broad1*broad2)

    //intermediate refers to distribution before smearing and Convolution is after smearing
    double intermediate1;
    double intermediate2;
    double intermediate3;
    double intermediate4;

    double Convolution1;
    double Convolution2;
    double Convolution3;
    double Convolution4;

    TVectorD t(3);
    t[0] = x[0];
    t[1] = x[1];
    t[2] = x[2];

    TMatrixD a1 = B1->GetCovMatrix();  
    TMatrixD a2 = B1->GetCovMatrix_b();  
    TMatrixD b1 = B2->GetCovMatrix();       
    TMatrixD b2 = B2->GetCovMatrix_b();       

    TMatrixD a1inv(a1);
    TMatrixD a2inv(a2);
    TMatrixD b1inv(b1);
    TMatrixD b2inv(b2);
    a1inv.Invert();
    a2inv.Invert();
    b1inv.Invert();
    b2inv.Invert();

    double W1 = B1->GetWeightN() * B2->GetWeightN();
    double W2 = B1->GetWeightN() * (1-(B2->GetWeightN()));
    double W3 = (1-(B1->GetWeightN())) * B2->GetWeightN();
    double W4 = (1-(B1->GetWeightN())) * (1-(B2->GetWeightN()));


    //Set the parameters of the smearing matrix
    TMatrixD SigmaRes;
    SigmaRes.ResizeTo(3,3);
    for (unsigned int cr=0; cr < 3; cr++){
        for (unsigned int cc=0; cc < 3; cc++){
            SigmaRes[cr][cc] = 0.0;
        }
    }

    SigmaRes[0][0] = 0.01;
    SigmaRes[1][1] = 0.01;
    SigmaRes[2][2] = 0.01;

    TMatrixD SigmaResinv(SigmaRes);
    SigmaResinv.Invert();

    TVectorD mu1(3);
    TVectorD mu2(3);
    mu1[0] =  bs_L_xyz_h_par_h/2;
    mu2[0] = -bs_L_xyz_h_par_h/2;
    mu1[1] = 0;
    mu2[1] = 0;
    mu1[2] = 0;
    mu2[2] = 0;

    TMatrixD scalar;
    scalar.ResizeTo(3,3);
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            scalar[i][j]=0;
        }
    }

    scalar[0][0]=10000;
    scalar[1][1]=10000;
    scalar[2][2]=10000;

    TMatrixD C1inv = (a1inv+b1inv);
    TMatrixD C2inv = (a1inv+b2inv);
    TMatrixD C3inv = (a2inv+b1inv);
    TMatrixD C4inv = (a2inv+b2inv);

    TMatrixD C1(C1inv);
    TMatrixD C2(C2inv);
    TMatrixD C3(C3inv);
    TMatrixD C4(C4inv);

    C1.Invert();
    C2.Invert();
    C3.Invert();
    C4.Invert();

    double alpha1 = exp((-0.5)*(mu1*(a1inv*mu1)+mu2*(b1inv*mu2)-(a1inv*mu1+b1inv*mu2)*(C1*(a1inv*mu1+b1inv*mu2))));
    double alpha2 = exp((-0.5)*(mu1*(a1inv*mu1)+mu2*(b2inv*mu2)-(a1inv*mu1+b2inv*mu2)*(C2*(a1inv*mu1+b2inv*mu2))));
    double alpha3 = exp((-0.5)*(mu1*(a2inv*mu1)+mu2*(b1inv*mu2)-(a2inv*mu1+b1inv*mu2)*(C3*(a2inv*mu1+b1inv*mu2))));
    double alpha4 = exp((-0.5)*(mu1*(a2inv*mu1)+mu2*(b2inv*mu2)-(a2inv*mu1+b2inv*mu2)*(C4*(a2inv*mu1+b2inv*mu2))));

    TVectorD M1 = C1*(a1inv*mu1+b1inv*mu2);
    TVectorD M2 = C2*(a1inv*mu1+b2inv*mu2);
    TVectorD M3 = C3*(a2inv*mu1+b1inv*mu2);
    TVectorD M4 = C4*(a2inv*mu1+b2inv*mu2);

    TMatrixD D1inv = (SigmaResinv+C1inv);
    TMatrixD D2inv = (SigmaResinv+C2inv);
    TMatrixD D3inv = (SigmaResinv+C3inv);
    TMatrixD D4inv = (SigmaResinv+C4inv);

    TMatrixD D1(D1inv);
    TMatrixD D2(D2inv);
    TMatrixD D3(D3inv);
    TMatrixD D4(D4inv);
    D1.Invert();
    D2.Invert();
    D3.Invert();
    D4.Invert();

    TMatrixD E1inv  = (SigmaResinv-SigmaResinv*D1*SigmaResinv);
    TMatrixD E2inv  = (SigmaResinv-SigmaResinv*D2*SigmaResinv);
    TMatrixD E3inv  = (SigmaResinv-SigmaResinv*D3*SigmaResinv);
    TMatrixD E4inv  = (SigmaResinv-SigmaResinv*D4*SigmaResinv);

    TMatrixD E1(E1inv);
    TMatrixD E2(E2inv);
    TMatrixD E3(E3inv);
    TMatrixD E4(E4inv);
    E1.Invert();
    E2.Invert();
    E3.Invert();
    E4.Invert();

    TVectorD Q1 = E1*SigmaResinv*D1*C1inv*M1;
    TVectorD Q2 = E2*SigmaResinv*D2*C2inv*M2;
    TVectorD Q3 = E3*SigmaResinv*D3*C3inv*M3;
    TVectorD Q4 = E4*SigmaResinv*D4*C4inv*M4;

    double gamma1  = exp(-(1/2)*(M1* (C1inv*M1)-M1*(C1inv*(D1*(C1inv*M1)))-Q1*(E1inv*Q1)));   
    double gamma2  = exp(-(1/2)*(M2* (C2inv*M2)-M2*(C2inv*(D2*(C2inv*M2)))-Q2*(E2inv*Q2)));   
    double gamma3  = exp(-(1/2)*(M3* (C3inv*M3)-M3*(C3inv*(D3*(C3inv*M3)))-Q3*(E3inv*Q3)));   
    double gamma4  = exp(-(1/2)*(M4* (C4inv*M4)-M4*(C4inv*(D4*(C4inv*M4)))-Q4*(E4inv*Q4)));   

    double delta1 = (-0.5)*(t-Q1)*(E1inv*(t-Q1));
    double delta2 = (-0.5)*(t-Q2)*(E2inv*(t-Q2));
    double delta3 = (-0.5)*(t-Q3)*(E3inv*(t-Q3));
    double delta4 = (-0.5)*(t-Q4)*(E4inv*(t-Q4));

    double beta1 = gamma1 * exp(delta1);
    double beta2 = gamma2 * exp(delta2);
    double beta3 = gamma3 * exp(delta3);
    double beta4 = gamma4 * exp(delta4);

    Convolution1 = (1./( Power(2*Pi(),3)))*sqrt(D1.Determinant())*pow((sqrt(a1.Determinant())*sqrt(b1.Determinant())*sqrt(SigmaRes.Determinant())),-1)*alpha1*beta1;
    Convolution2 = (1./( Power(2*Pi(),3)))*sqrt(D2.Determinant())*pow((sqrt(a1.Determinant())*sqrt(b2.Determinant())*sqrt(SigmaRes.Determinant())),-1)*alpha2*beta2;
    Convolution3 = (1./( Power(2*Pi(),3)))*sqrt(D3.Determinant())*pow((sqrt(a2.Determinant())*sqrt(b1.Determinant())*sqrt(SigmaRes.Determinant())),-1)*alpha3*beta3;
    Convolution4 = (1./( Power(2*Pi(),3)))*sqrt(D4.Determinant())*pow((sqrt(a2.Determinant())*sqrt(b2.Determinant())*sqrt(SigmaRes.Determinant())),-1)*alpha4*beta4;

    double betai1 = exp((-0.5)*(t-M1)*(C1inv*(t-M1)));
    double betai2 = exp((-0.5)*(t-M2)*(C2inv*(t-M2)));
    double betai3 = exp((-0.5)*(t-M3)*(C3inv*(t-M3)));
    double betai4 = exp((-0.5)*(t-M4)*(C4inv*(t-M4)));

    intermediate1 = (1./(Power(2*Pi(),3)))*pow((sqrt(a1.Determinant())*sqrt(b1.Determinant())),-1)*alpha1*betai1;
    intermediate2 =  (1./(Power(2*Pi(),3)))*pow((sqrt(a1.Determinant())*sqrt(b2.Determinant())),-1)*alpha2*betai2;
    intermediate3 =  (1./(Power(2*Pi(),3)))*pow((sqrt(a2.Determinant())*sqrt(b1.Determinant())),-1)*alpha3*betai3;
    intermediate4 =  (1./(Power(2*Pi(),3)))*pow((sqrt(a2.Determinant())*sqrt(b2.Determinant())),-1)*alpha4*betai4;

    //  return -(W1*intermediate1+W2*intermediate2+W3*intermediate3+W4*intermediate4);
    return -(W1*Convolution1+W2*Convolution2+W3*Convolution3+W4*Convolution4);
}

//Wrap of analytic convolution fixing z=0 for TF2 display 
double overlap_num_TF2(double *x, double *par)
{
    double xnew[3];
    xnew[0] = x[0]; xnew[1] = x[1]; 
    xnew[2] = par[0]; //fix z=par[0]
    return overlap_num(xnew);
}

//Utility for drawing beam profiles
double beam1Profile_TF2(double *x, double*par)
{
    return B1->rho(x[0], x[1], par[0], par[1]);
}

//Utility for drawing beam profiles
double beam2Profile_TF2(double *x, double*par)
{
    return B2->rho(x[0], x[1], par[0], par[1]);
}

/** Study movement of beamspot when beams are displaced.
 * Compared to previous functions, we use a different approach.
 * Integrate only over BeamProfile::cS, then find the maximum
 * in (X,Y,Z) and plot these positions as function of the displacement h.
 *
 * Beam positions are given in the ATLAS coordinate system and 
 * separations in the LHC coordinate system 
 */
void bs_xyz_dh() {

    cout << endl;
    cout << "********************************************************************************" << endl;
    cout << "STUDYING MOVEMENT OF BEAMSPOT <X>,<Y>,<Z> AS FUNCTION OF BEAM SEPARATION IN ";
    if (bs_xyz_dh_bsep_coord == BeamProfile::cX)
        cout << "X";
    else
        cout << "Y";
    cout << endl;
    cout << "********************************************************************************" << endl;
    cout << endl;

    //hack beam properties (otherwise set by command-line arguments or default)
    //B1->SetBeamXY(0.0652+TMath::Sin(0.0001)*73.54,0.0552+TMath::Sin(0.0001)*73.54, 0.0);
    //B2->SetBeamXY(0.0452+TMath::Sin(0.0001)*73.54,0.0552+TMath::Sin(0.0001)*73.54, 0.0);

    B1->Print();
    B2->Print();

    //Scan as function of h
    Int_t nbins = (bs_L_xyz_h_par_max_h - bs_L_xyz_h_par_min_h) / bs_L_xyz_h_par_step_h;
    TGraphErrors *h_bs_L_x_h = new TGraphErrors(nbins);  
    h_bs_L_x_h->GetHistogram()->GetXaxis()->SetTitle("Separation (mm)");
    h_bs_L_x_h->GetHistogram()->GetYaxis()->SetTitle("<X> (mm)");
    h_bs_L_x_h->SetName("x_pos");
    TGraphErrors *h_bs_L_y_h = new TGraphErrors(nbins);  
    h_bs_L_y_h->GetHistogram()->GetXaxis()->SetTitle("Separation (mm)");
    h_bs_L_y_h->GetHistogram()->GetYaxis()->SetTitle("<Y> (mm)");
    h_bs_L_y_h->SetName("y_pos");
    TGraphErrors *h_bs_L_z_h = new TGraphErrors(nbins); 
    h_bs_L_z_h->GetHistogram()->GetXaxis()->SetTitle("Separation (mm)");
    h_bs_L_z_h->GetHistogram()->GetYaxis()->SetTitle("<Z> (mm)");
    h_bs_L_z_h->SetName("z_pos");

    //width plots
    TGraphErrors *h_width_L_x_h = new TGraphErrors(nbins);  
    h_width_L_x_h->GetHistogram()->GetXaxis()->SetTitle("Separation (mm)");
    h_width_L_x_h->GetHistogram()->GetYaxis()->SetTitle("sigmaX (mm)");
    h_width_L_x_h->SetName("x_width");
    TGraphErrors *h_width_L_y_h = new TGraphErrors(nbins);  
    h_width_L_y_h->GetHistogram()->GetXaxis()->SetTitle("Separation (mm)");
    h_width_L_y_h->GetHistogram()->GetYaxis()->SetTitle("sigmaY (mm)");
    h_width_L_y_h->SetName("y_width");
    TGraphErrors *h_width_L_z_h = new TGraphErrors(nbins); 
    h_width_L_z_h->GetHistogram()->GetXaxis()->SetTitle("Separation (mm)");
    h_width_L_z_h->GetHistogram()->GetYaxis()->SetTitle("sigmaZ (mm)");
    h_width_L_z_h->SetName("z_width");

    Int_t idx=0;

    //Declare root files and trees for fitting analysis

    TTree *Vertices = new TTree("Vertices","Test");

    //sampled numbers x,y,z
    double samnumx;
    double samnumy;
    double samnumz;
    double vxx;
    double vyy;
    double vzz;
    double vxy;
    double vxz;
    double vyz;

    //Declare branches required for fitting program

    TBranch *b_x = Vertices->Branch("x",&samnumx);
    TBranch *b_y = Vertices->Branch("y",&samnumy);
    TBranch *b_z = Vertices->Branch("z",&samnumz);
    TBranch *b_vxx = Vertices->Branch("vxx",&vxx);
    TBranch *b_vyy = Vertices->Branch("vyy",&vyy);
    TBranch *b_vzz = Vertices->Branch("vzz",&vzz);
    TBranch *b_vxy = Vertices->Branch("vxy",&vxy);
    TBranch *b_vxz = Vertices->Branch("vxz",&vxz);
    TBranch *b_vyz = Vertices->Branch("vyz",&vyz);

    //using the old minimisation mathod

    ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2");
    if (!min) {
        cerr << "ERROR in creating minimizer." << endl;
        throw;
    }

    min->SetMaxFunctionCalls(1000000);
    min->SetTolerance(0.000001);
    min->SetMaxIterations(100000); //only if using GSL algo

    //   min->SetMaxFunctionCalls(100);
    //   min->SetMaxIterations(100); //only if using GSL algo
    //   min->SetTolerance(0.1);

    if (debug >= 10)
        min->SetPrintLevel(10);
    //smeared distribution (numerical)  
    //      ROOT::Math::Functor f_bs_L_xyz(&bs_L_xyz_h_conv, 3);

    //unsmeared numerical distribution
    ROOT::Math::Functor f_bs_L_xyz(&bs_L_xyz_h, 3);

    //analytical distribution
    //  ROOT::Math::Functor f_bs_L_xyz(&overlap_num, 3);

    double steps[3] = {0.0001, 0.0001, 0.1};
    double start[3] = {0.0, 0.0, 0.0};  

    min->SetFunction(f_bs_L_xyz);
    min->SetVariable(0, "X", start[0], steps[0]);
    min->SetVariable(1, "Y", start[1], steps[1]);
    min->SetVariable(2, "Z", start[2], steps[2]); 


    for (double h=bs_L_xyz_h_par_min_h; h<bs_L_xyz_h_par_max_h; h+=bs_L_xyz_h_par_step_h) {
        if (h > -0.000001 && h < 0.000001) h = 0;

        cout << h << endl;
        bs_L_xyz_h_par_h = h;

        bs_L_xyz_h_par_driftx1 = bs_L_xyz_h_par_driftx1_gradient * idx * bs_L_xyz_h_par_step_h / 0.05;
        bs_L_xyz_h_par_driftx2 = bs_L_xyz_h_par_driftx2_gradient * idx * bs_L_xyz_h_par_step_h / 0.05;
        bs_L_xyz_h_par_drifty1 = bs_L_xyz_h_par_drifty1_gradient * idx * bs_L_xyz_h_par_step_h / 0.05;
        bs_L_xyz_h_par_drifty2 = bs_L_xyz_h_par_drifty2_gradient * idx * bs_L_xyz_h_par_step_h / 0.05;

        //     cout << bs_L_xyz_h_par_driftx1 << endl;
        //     cout << bs_L_xyz_h_par_driftx2 << endl;


        // Choose definition of the beam spot, either the minimum of overlap integral, the mean, the mean from sampling., or manually input

        if (bs_def == 0){

            if (debug >= 5)
                cout << "Maximizing h=" << h << "...";

            min->Minimize();      
            const double *resXYZ = min->X();
            const double *resErrorXYZ = min->Errors();

            if (debug >=5)
                cout << "X= " << resXYZ[0] << ", Y=" << resXYZ[1] << ", Z=" << resXYZ[2] << ", Value=" << -min->MinValue() << endl;    
            h_bs_L_x_h->SetPoint(idx, h, resXYZ[0]);
            h_bs_L_y_h->SetPoint(idx, h, resXYZ[1]);
            h_bs_L_z_h->SetPoint(idx, h, resXYZ[2]);
            //       h_bs_L_x_h->SetPointError(idx, 0.0, resErrorXYZ[0]);    
            //       h_bs_L_y_h->SetPointError(idx, 0.0, resErrorXYZ[1]);    
            //       h_bs_L_z_h->SetPointError(idx, 0.0, resErrorXYZ[2]);    


        }

        if (bs_def == 1){

            double x_val[3]={0,0,0};
            double val = OverlapIntegral_xyz_h(bs_xyz_dh_bsep_coord, bs_L_xyz_h_par_h, x_val, bs_L_xyz_h_par_offset,bs_L_xyz_h_par_driftx1,bs_L_xyz_h_par_driftx2, bs_L_xyz_h_par_drifty1,bs_L_xyz_h_par_drifty2, true);
            double overlap = OverlapIntegral(false); 
            double resXYZ[3]={0};
            double resXYZ2[3]={0};
            resXYZ[0] = OverlapIntegral_mean(BeamProfile::cX,false)/overlap;
            resXYZ[1] = OverlapIntegral_mean(BeamProfile::cY,false)/overlap;
            resXYZ[2] = OverlapIntegral_mean(BeamProfile::cZ,false)/overlap;

            resXYZ2[0] = OverlapIntegral_x2(BeamProfile::cX,false)/overlap;
            resXYZ2[1] = OverlapIntegral_x2(BeamProfile::cY,false)/overlap;
            resXYZ2[2] = OverlapIntegral_x2(BeamProfile::cZ,false)/overlap;

            if (debug >=5)
                cout << "X= " << resXYZ[0] << ", Y=" << resXYZ[1] << ", Z=" << resXYZ[2] << endl;

            h_bs_L_x_h->SetPoint(idx, h, resXYZ[0]);
            h_bs_L_y_h->SetPoint(idx, h, resXYZ[1]);
            h_bs_L_z_h->SetPoint(idx, h, resXYZ[2]);


            h_width_L_x_h->SetPoint(idx, h, sqrt(resXYZ2[0]-pow(resXYZ[0],2)));
            h_width_L_y_h->SetPoint(idx, h, sqrt(resXYZ2[1]-pow(resXYZ[1],2)));
            h_width_L_z_h->SetPoint(idx, h, sqrt(resXYZ2[2]-pow(resXYZ[2],2)));

        }

        if (bs_def == 2){

            //Name the output root files containing the Vertices tree
            ostringstream filename;

            filename << "scanx-" <<  h*100 << ".root";

            //Perform sampling and fill a tree if desired

            TF3 * f_L_xyz0_h_TF3 = new TF3("f_L_xyz0_h_TF3", bs_L_xyz0_h_TF3,-0.2,0.2,-0.2,0.2,-100,100);

            TFile *outRoot = TFile::Open(filename.str().c_str(), "RECREATE");
            Vertices->Reset();

            int nsamples = 10000;
            double samplednumbersx[nsamples];
            double samplednumbersy[nsamples];
            double samplednumbersz[nsamples];
            double meanx = 0;
            double meany = 0;
            double meanz = 0;

            for (int i = 0;i<nsamples;i++){
                //Get a random point from overlap distribution
                f_L_xyz0_h_TF3->GetRandom3(samplednumbersx[i],samplednumbersy[i],samplednumbersz[i]);

                samnumx = samplednumbersx[i];
                samnumy = samplednumbersy[i];
                samnumz = samplednumbersz[i];

                meanx += samplednumbersx[i];
                meany += samplednumbersy[i];
                meanz += samplednumbersz[i];

                //Set vertex errors
                vxx = 0;
                vyy = 0;
                vzz = 0;
                vxy = 0;
                vxz = 0;
                vyz = 0;

                Vertices->Fill();

            }

            meanx = meanx/nsamples;
            meany = meany/nsamples;
            meanz = meanz/nsamples;

            Vertices->Write();
            outRoot->Close();


            cout << "X= " << meanx << ", Y=" << meany << ", Z=" << meanz << endl;

            h_bs_L_x_h->SetPoint(idx, h, meanx);
            h_bs_L_y_h->SetPoint(idx, h, meany);
            h_bs_L_z_h->SetPoint(idx, h, meanz); 

        }    


        if (bs_def == 3){

            //      Set the positions of the beam spot manually
            double resX_user[30]={0.0510,0.05,0.0446,0.042,0.0373,0.0312,0.0282,0.0224,0.0148,0.0093,0.00504,0.00219,0.00161,0.000618,0.000328,0,-0.000398,-0.000678,-0.00141,-0.00283,-0.00457,-0.00941,-0.0151,-0.05,-0.0289,-0.0335,-0.0381,-0.0419,-0.05,-0.0486};
            double resY_user[30]={0};
            double resZ_user[30]={0};

            //Manually input the values by user
            h_bs_L_x_h->SetPoint(idx, h, resX_user[idx]);
            h_bs_L_y_h->SetPoint(idx, h, resY_user[idx]);
            h_bs_L_z_h->SetPoint(idx, h, resZ_user[idx]);
        }

        //display, if requested
        if (debug >= 6) {
            double minAxis=-0.2;
            double maxAxis=0.2;

            TF2 *f_bs_L_xyz0_h = new TF2("f_bs_L_xyz0_h",bs_L_xyz0_h_TF2,minAxis, maxAxis, minAxis, maxAxis, 1);
            f_bs_L_xyz0_h->SetParameter(0, 0.0); //fix z=0

            TF2 *f_bs_L_xyz0_h_analytical = new TF2("f_bs_L_xyz0_h_analytical",overlap_num_TF2,minAxis, maxAxis, minAxis, maxAxis, 1);
            f_bs_L_xyz0_h_analytical->SetParameter(0, 0.0); //fix z=0

            TF2 *f_bs_L_xyz0_h_b1 = new TF2("f_bs_L_xyz0_h_b1",beam1Profile_TF2,minAxis, maxAxis, minAxis, maxAxis, 2);
            f_bs_L_xyz0_h_b1->SetParameter(0, 0.0);  //fix z=0
            f_bs_L_xyz0_h_b1->SetParameter(1, 0.0);  //fix s=0

            TF2 *f_bs_L_xyz0_h_b2 = new TF2("f_bs_L_xyz0_h_b2",beam2Profile_TF2,minAxis, maxAxis, minAxis, maxAxis, 2);
            f_bs_L_xyz0_h_b2->SetParameter(0, 0.0);  //fix z=0
            f_bs_L_xyz0_h_b2->SetParameter(1, 0.0);  //fix s=0


            gStyle->SetPalette(1); //Rainbow palette

            TCanvas *c_bs_L_xyz0_h = new TCanvas();
            c_bs_L_xyz0_h->Divide(2,2);
            c_bs_L_xyz0_h->cd(1);


            TPad* c1_1=(TPad*)(c_bs_L_xyz0_h->GetPad(1));
            c1_1->SetGridx();
            c1_1->SetGridy();
            f_bs_L_xyz0_h_b1->Draw("CONTZ");
            c_bs_L_xyz0_h->cd(2);

            TPad* c1_2=(TPad*)(c_bs_L_xyz0_h->GetPad(2));
            c1_2->SetGridx();
            c1_2->SetGridy();
            f_bs_L_xyz0_h_b2->Draw("CONTZ");
            c_bs_L_xyz0_h->cd(3);

            TPad* c1_3=(TPad*)(c_bs_L_xyz0_h->GetPad(3));
            c1_3->SetGridx();
            c1_3->SetGridy();
            f_bs_L_xyz0_h->Draw("CONTZ");
            // testtest->Draw("BOX");


            //       c_bs_L_xyz0_h->cd(4);
            //       TPad* c1_4=(TPad*)(c_bs_L_xyz0_h->GetPad(4));
            //       c1_4->SetGridx();
            //       c1_4->SetGridy();
            //       f_bs_L_xyz0_h_analytical->Draw("CONTZ");

            cout << "Done. Press key to continue." << endl;

            c_bs_L_xyz0_h->SetGridy();
            c_bs_L_xyz0_h->SetGridx();

            c_bs_L_xyz0_h->WaitPrimitive();
            delete f_bs_L_xyz0_h;
            //      delete f_bs_L_xyz0_h_analytical;
            delete f_bs_L_xyz0_h_b1;
            delete f_bs_L_xyz0_h_b2;
            delete c_bs_L_xyz0_h;
        }

        idx++;
    }
    delete min;

    //Draw x,y,z (h), and make a linear fit for slopes. Expected to be ~ straight line.
    TCanvas *c_bs_xyz_dh = new TCanvas("c_bs_xyz_dh", "Beamspot d(x,y,z)/dh");
    c_bs_xyz_dh->Divide(2,2);

    c_bs_xyz_dh->cd(1);
    h_bs_L_x_h->Draw("AP");

    /*  h_bs_L_x_h->Fit("pol1");

        TF1 *fit_pol1_x = h_bs_L_x_h->GetFunction("pol1");
        double slope_x=0.0;
        double slopeError_x = 0.0;
        if (fit_pol1_x) {
        slope_x = fit_pol1_x->GetParameter(1);
        slopeError_x = fit_pol1_x->GetParError(1);  
        }
     */
    c_bs_xyz_dh->cd(2);
    h_bs_L_y_h->Draw("AP");

    /*  h_bs_L_y_h->Fit("pol1");

        TF1 *fit_pol1_y = h_bs_L_y_h->GetFunction("pol1");
        double slope_y=0.0;
        double slopeError_y = 0.0;
        if (fit_pol1_y) {
        slope_y = fit_pol1_y->GetParameter(1);
        slopeError_y = fit_pol1_y->GetParError(1);  
        }
     */

    c_bs_xyz_dh->cd(3);
    h_bs_L_z_h->Draw("AP");

    /*  h_bs_L_z_h->Fit("pol1");
        TF1 *fit_pol1_z = h_bs_L_z_h->GetFunction("pol1");
        double slope_z=0.0;
        double slopeError_z = 0.0;
        if (fit_pol1_z) {
        slope_z = fit_pol1_z->GetParameter(1);
        slopeError_z = fit_pol1_z->GetParError(1);  
        }
     */

    c_bs_xyz_dh->Update();

    //   cout << "********************************************************************************" << endl;
    //   cout << "RESULTS:" << endl;
    //   cout << "********************************************************************************" << endl;
    //   cout << "dx/dh = " << slope_x << " +/- " << slopeError_x << endl;
    //   cout << "dy/dh = " << slope_y << " +/- " << slopeError_y << endl;
    //   cout << "dz/dh = " << slope_z << " +/- " << slopeError_z << endl;

    cout << "********************************************************************************" << endl;
    cout << "FINISHED" << endl;
    cout << "********************************************************************************" << endl;
    cout << "FINISHED" << endl;


    //Save results


    //  TFile *outRootF = TFile::Open("bs_dxyzdh.root", "RECREATE");
    //  TFile *outRootF = TFile::Open("bs_dxyzdh_set6_offset_x.root", "RECREATE");

    if (outputfilename == "default") TFile *outRootF = TFile::Open("bs_dxyzdh.root", "RECREATE");

    TFile *outRootF = TFile::Open(outputfilename, "RECREATE");
    c_bs_xyz_dh->Write();
    h_bs_L_x_h->Write();
    h_bs_L_y_h->Write();
    h_bs_L_z_h->Write();
    h_width_L_x_h->Write();
    h_width_L_y_h->Write();
    h_width_L_z_h->Write();

    outRootF->Close();

    //  output->Write();

}

// Study Luminosity variation with as function of beam parameters
// Implemented: rho
//----------------
enum L_scan_param {
    lsp_rho,
    lsp_sepX,
    lsp_sepY
};
double minLScan=0;
double maxLScan=0;
double stepLScan=1;

// Choosing whether to scan as function of x, y or rho, based on user inputs
void L_scan(L_scan_param p)
{
    TString pName="";
    //   TString pNameOffset="nooffset";
    //   if (bs_L_xyz_h_par_offset != 0)
    //     pNameOffset = "offset";
    if (p == lsp_rho)
        pName = "rho";
    else if (p == lsp_sepX)
        pName = "sepX";
    else if (p == lsp_sepY)
        pName = "sepY";

    cout << endl;
    cout << "********************************************************************************" << endl;
    cout << "STUDYING LUMINOSITY VARIATIONS AS FUNCTION OF PARAMETER: " << pName << endl;
    cout << "********************************************************************************" << endl;
    cout << endl;

    cout << "Initial beam configuration:" << endl;
    cout << "----------" << endl;

    //print the configuration setting, (covariance matrix, crossing angle etc) for both beams.

    //   B1->Print();
    //   B2->Print();

    cout << "Scanning " << pName << ": " << minLScan << " - " << maxLScan << ", step = " << stepLScan << endl;

    //Making sure that the maximum scan position minus the minimum scan position divided by the step size is an integer.

    Double_t nbinsSep = ((maxLScan - minLScan) / stepLScan) + 1;

    if ( abs(((maxLScan - minLScan) / stepLScan) - floor(((maxLScan - minLScan) / stepLScan))) > 0.000001){
        if ( abs (((maxLScan - minLScan) / stepLScan) - ceil(((maxLScan - minLScan) / stepLScan))) > 0.000001){
            cerr << "Please divide scan space evenly" << endl;
            throw (1);
        }
    }

    if ( abs (((maxLScan - minLScan) / stepLScan) - ceil(((maxLScan - minLScan) / stepLScan))) > 0.000001){
        if ( abs(((maxLScan - minLScan) / stepLScan) - floor(((maxLScan - minLScan) / stepLScan))) > 0.000001){
            cerr << "Please divide scan space evenly" << endl;
            throw (1);
        }
    }

    if (nbinsSep < 0.5 + floor(nbinsSep)) nbinsSep = floor(nbinsSep);
    if (nbinsSep > 0.5 + floor(nbinsSep)) nbinsSep = ceil(nbinsSep);

    cout << "Bins " << nbinsSep << " from " << minLScan - stepLScan/2 << " to " << maxLScan+stepLScan/2 << endl;
    TH1F *L_rel_sep = new TH1F("L_rel_sep", "Luminosity as function of separation", nbinsSep, minLScan-stepLScan/2, maxLScan+stepLScan/2);
    L_rel_sep->SetXTitle(pName);

    for (double pv=minLScan; pv <= maxLScan; pv+=stepLScan) {
        //p states whether the separation is in x, y or rho
        switch (p) {
            case lsp_rho:

                B1->SetBeamXY(B1->GetSigmaX(), B1->GetSigmaY(), pv);
                B1->SetBeamXY_b(B1->GetSigmaXb(), B1->GetSigmaY(), pv);

                B2->SetBeamXY(B2->GetSigmaX(), B2->GetSigmaY(), pv);
                B2->SetBeamXY_b(B2->GetSigmaXb(), B2->GetSigmaY(), pv);
                //B1->Print();
                //B2->Print();
                break;
            case lsp_sepX:
                //whether beam 1 or 2 is moving forward or backward
                //Beam positions in ATLAS coordinate system, separations in LHC coordinate system

                B1->SetBeamOffsetX(pv/2);
                B2->SetBeamOffsetX(-pv/2);
                B1->SetBeamOffsetY(-bs_L_xyz_h_par_offset/2);
                B2->SetBeamOffsetY(+bs_L_xyz_h_par_offset/2);
                break;
            case lsp_sepY:
                B1->SetBeamOffsetY(-pv/2);
                B2->SetBeamOffsetY(+pv/2);
                B1->SetBeamOffsetX(bs_L_xyz_h_par_offset/2);
                B2->SetBeamOffsetX(-bs_L_xyz_h_par_offset/2);
                break;
            default:
                cerr << "UNDEFINED SCAN PARAMETER: " << lsp_rho << endl;
                throw;
        }

        double Lrel = OverlapIntegral(false);

        L_rel_sep->SetBinContent(L_rel_sep->FindBin(pv), Lrel);
        cout << pName << " = " << pv << " -> L = " << Lrel << " (bin = " << L_rel_sep->FindBin(pv) << ")" << endl;

    }

    TGraph * L_rel_sep_graph = new TGraph(L_rel_sep);

    //Fit results
    TCanvas *c_L_rel_sep = new TCanvas();
    L_rel_sep->Draw();

    TF1 *doublegaussian = new TF1("doublegaussian","gaus(0)+gaus(3)");//,minLScan - stepLScan/2,maxLScan+stepLScan/2);
    doublegaussian->SetParameters(L_rel_sep->GetMaximum(), L_rel_sep->GetMean(), L_rel_sep->GetRMS(), L_rel_sep->GetMaximum()/50, L_rel_sep->GetMean(), L_rel_sep->GetRMS()*5); 

    if (num_gaus == 1){  
        //        L_rel_sep->Fit("gaus");
        //    L_rel_sep_graph->Fit("gaus");

        //     cout << endl << "----------" << endl << "SIGMA = " << L_rel_sep->GetFunction("gaus")->GetParameter(2) << endl;
        //     cout << endl << "----------" << endl << "CHI SQUARE = " << L_rel_sep->GetFunction("gaus")->GetChisquare() << endl;
        //     cout << endl << "----------" << endl << "NO. DOF = " << L_rel_sep->GetFunction("gaus")->GetNDF() << endl;
        //     cout << endl << "----------" << endl << "CHI2/DOF = " << (L_rel_sep->GetFunction("gaus")->GetChisquare())/(L_rel_sep->GetFunction("gaus")->GetNDF()) << endl;
    }

    if (num_gaus == 2){
        L_rel_sep->Fit(doublegaussian);
        cout << endl << "----------" << endl << "SIGMA A = " << doublegaussian->GetParameter(2) << endl;
        cout << endl << "----------" << endl << "SIGMA B = " << doublegaussian->GetParameter(5) << endl;
        cout << endl << "----------" << endl << "CHI SQUARE = " << doublegaussian->GetChisquare() << endl;
        cout << endl << "----------" << endl << "NO. DOF = " << doublegaussian->GetNDF() << endl;
        cout << endl << "----------" << endl << "CHI2/DOF = " << (doublegaussian->GetChisquare())/(doublegaussian->GetNDF()) << endl;
    }

    c_L_rel_sep->Update();

    //Save results

    TString outRootFName = "Lsp_";
    outRootFName += pName;
    outRootFName += ".root";
    if (outputfilename == "default") TFile *outRootF = TFile::Open(outRootFName.Data(), "RECREATE");

    TFile *outRootF = TFile::Open(outputfilename, "RECREATE");

    L_rel_sep->Write();
    L_rel_sep_graph->Write();
    c_L_rel_sep->Write();
    outRootF->Close();

    cout << "Integral = " << L_rel_sep_graph->Integral() << endl;

    cout << "End lsp study" << endl;

}


void usage(char **argv) {
    std::cerr << "Usage: " << argv[0] << " [options] action" << endl;

    std::cerr << "Actions: " << endl;
    std::cerr << " bs_xyz_h Study beamspot X,Y,Z movements as function of beam displacement" << endl;
    std::cerr << " lspRho,lspX,lspY Study luminosity variations as function of parameter rho, separation X and Y" << endl;

    std::cerr << "Generic options:" << endl;
    std::cerr << "\t -h Print help" << endl;
    //  std::cerr << "\t -g, --numgaus: Set no. gaussians (default: 1)" << endl;
    std::cerr << "\t -x, --sigmax: Set sigma X a of both beams" << endl;  
    std::cerr << "\t -y, --sigmay: Set sigma Y a of both beams" << endl;  
    std::cerr << "\t -z, --sigmaz: Set sigma Z a of both beams" << endl;  
    std::cerr << "\t -s, --sigma: Set sigma X,Y a of both beams" << endl;  
    std::cerr << "\t --sigmaxa1: Set sigma X a of Beam 1" << endl;  
    std::cerr << "\t --sigmaya1: Set sigma Y a of Beam 1" << endl;  
    std::cerr << "\t --sigmaza1: Set sigma Z a of Beam 1" << endl;  
    std::cerr << "\t --sigmaxa2: Set sigma X a of Beam 2" << endl;  
    std::cerr << "\t --sigmaya2: Set sigma Y a of Beam 2" << endl;  
    std::cerr << "\t --sigmaza2: Set sigma Z a of Beam 2" << endl;  
    std::cerr << "\t --sigmaxb1: Set sigma X b of Beam 1" << endl;  
    std::cerr << "\t --sigmayb1: Set sigma Y b of Beam 1" << endl;  
    std::cerr << "\t --sigmazb1: Set sigma Z b of Beam 1" << endl;  
    std::cerr << "\t --sigmaxb2: Set sigma X b of Beam 2" << endl;  
    std::cerr << "\t --sigmayb2: Set sigma Y b of Beam 2" << endl;  
    std::cerr << "\t --sigmazb2: Set sigma Z b of Beam 2" << endl;  
    std::cerr << "\t -r, --rho: Set XY rho a of covariance matrix for both beams" << endl;  
    std::cerr << "\t --rho1: Set XY rho of covariance matrix for Beam 1" << endl;  
    std::cerr << "\t --rho2: Set XY rho of covariance matrix for Beam 2" << endl;  
    std::cerr << "\t --rhoa1: Set XY rho of covariance matrix a for Beam 1" << endl;  
    std::cerr << "\t --rhoa2: Set XY rho of covariance matrix a for Beam 2" << endl;  
    std::cerr << "\t --rhob1: Set XY rho of covariance matrix b for Beam 1" << endl;  
    std::cerr << "\t --rhob2: Set XY rho of covariance matrix b for Beam 2" << endl;  
    std::cerr << "\t --weighta1: Set weight of 'A' Gaussian for Beam 1" << endl;  
    std::cerr << "\t --weighta2: Set weight of 'A' Gaussian for Beam 2" << endl;  
    std::cerr << "\t -a, --xingyz: Set crossing angle on YZ plane of each beam (rad) (half of beam-beam crossing angle)" << endl;  
    std::cerr << "\t --xingxz: Set crossing angle on XZ plane of each beam (rad) (half of beam-beam crossing angle)" << endl;  
    std::cerr << "\t -o, --beamoffset: Set offset on transverse coordinate (mm)" << endl;  
    std::cerr << "\t --dx1: Set drift gradient (units mm per scan step)" << endl;  
    std::cerr << "\t --dx2: Set drift gradient (units mm per scan step)" << endl;  
    std::cerr << "\t --dy1: Set drift gradient (units mm per scan step)" << endl;  
    std::cerr << "\t --dy2: Set drift gradient (units mm per scan step)" << endl;  

    std::cerr << "\t -f, --filename: Set output file name" << endl;  

    std::cerr << "Specific options:" << endl;

    std::cerr << " bs_xyz_h:" << endl;
    std::cerr << "\t -c (--coord): Set coordinate for beam offset (x or y)" << endl;
    std::cerr << "\t --minh: Minimum beam separation (default: "<< bs_L_xyz_h_par_min_h << " mm)" << endl;
    std::cerr << "\t --maxh: Maximum beam separation (default: "<< bs_L_xyz_h_par_max_h << " mm)" << endl;
    std::cerr << "\t --steph: Step of beam separation (default: "<< bs_L_xyz_h_par_step_h << " mm)" << endl;

    std::cerr << " lspRho, lspX, lspY:" << endl;
    std::cerr << "\t --min: Minimum beam separation (default: "<< minLScan << " mm)" << endl;
    std::cerr << "\t --max: Maximum beam separation (default: "<< maxLScan << " mm)" << endl;
    std::cerr << "\t --step: Step of beam separation (default: "<< stepLScan << " mm)" << endl;
}

int main(int argc, char **argv)
{   
    TApplication *theApp=0;
    if (debug>= 6) 
        theApp = new TApplication("App", &argc, argv);

    //Flags for command-line argument and their default values
    double sigma_x_1 = 0.0652; // sigma_X, Beam 1
    double sigma_x_b_1 = 0.0652; // sigma_X_b, Beam 1
    double sigma_x_2 = 0.0652; // sigma_X, Beam 2
    double sigma_x_b_2 = 0.0652; // sigma_X_b, Beam 2
    double sigma_y_1 = 0.0652; // sigma_Y, Beam 1
    double sigma_y_b_1 = 0.0652; // sigma_Y_b, Beam 1
    double sigma_y_2 = 0.0652; // sigma_Y, Beam 2
    double sigma_y_b_2 = 0.0652; // sigma_Y_b, Beam 2
    double sigma_z_1 = 73.54; // sigma_Z, Beam 1
    double sigma_z_b_1 = 73.54; // sigma_Z_b, Beam 1
    double sigma_z_2 = 73.54; // sigma_Z, Beam 2
    double sigma_z_b_2 = 73.54; // sigma_Z_b, Beam 2
    double rho_xy_1 = 0.0; // rho_XY, Beam 1
    double rho_xy_b_1 = 0.0; // rho_XY_b, Beam 1
    double rho_xy_2 = 0.0; // rho_XY, Beam 2
    double rho_xy_b_2 = 0.0; // rho_XY_b, Beam 2
    double weight_a_1 = 1; // weight Gaussian a, Beam 1
    double weight_a_2 = 1; // weight Gaussian a, Beam 2
    double weight_b_1 = 0; // weight Gaussian b, Beam 1
    double weight_b_2 = 0; // weight Gaussian b, Beam 2
    double xing_yz = 0.0001; // crossing angle, YZ plane
    double xing_xz = 0.0; // crossing angle, XZ plane  

    //triple gaussian
    double sigma_x_c_1 = 0.05; // sigma_X_c, Beam 1
    double sigma_x_c_2 = 0.05; // sigma_X_c, Beam 2
    double sigma_y_c_1 = 0.04; // sigma_Y_c, Beam 1
    double sigma_y_c_2 = 0.06;  // sigma_Y_c, Beam 2
    double sigma_z_c_1 = 85; // sigma_Z_c, Beam 1
    double sigma_z_c_2 = 85; // sigma_Z_c, Beam 2
    double rho_xy_c_1 = 0.0; // rho_XY_c, Beam 1
    double rho_xy_c_2 = 0.0; // rho_XY_c, Beam 2


    //Process command-line arguments
    while (1) {
        //Flags
        static struct option long_options[] = {
            //General options
            {"help", no_argument, 0, 'h'},
            {"sigmax", required_argument, 0, 'x'},
            {"sigmay", required_argument, 0, 'y'},
            {"sigmaz", required_argument, 0, 'z'},
            {"sigma", required_argument, 0, 's'},
            {"sigmaxa1", required_argument, 0, 0},
            {"sigmaya1", required_argument, 0, 0},
            {"sigmaza1", required_argument, 0, 0},
            {"sigmaxa2", required_argument, 0, 0},
            {"sigmaya2", required_argument, 0, 0},
            {"sigmaza2", required_argument, 0, 0},

            {"sigmaxb1", required_argument, 0, 0},
            {"sigmayb1", required_argument, 0, 0},
            {"sigmazb1", required_argument, 0, 0},
            {"sigmaxb2", required_argument, 0, 0},
            {"sigmayb2", required_argument, 0, 0},
            {"sigmazb2", required_argument, 0, 0},

            {"sigmaxc1", required_argument, 0, 0},
            {"sigmayc1", required_argument, 0, 0},
            {"sigmazc1", required_argument, 0, 0},
            {"sigmaxc2", required_argument, 0, 0},
            {"sigmayc2", required_argument, 0, 0},
            {"sigmazc2", required_argument, 0, 0},

            {"rho", required_argument, 0, 'r'},
            {"rhoa1", required_argument, 0, 0},
            {"rhoa2", required_argument, 0, 0},
            {"rhob1", required_argument, 0, 0},
            {"rhob2", required_argument, 0, 0},
            {"rhoc1", required_argument, 0, 0},
            {"rhoc2", required_argument, 0, 0},

            {"weighta1", required_argument, 0, 0},
            {"weighta2", required_argument, 0, 0},
            {"weightb1", required_argument, 0, 0},
            {"weightb2", required_argument, 0, 0},

            {"xingyz", required_argument, 0, 'a'},
            {"xingxz", required_argument, 0, 0},

            {"numgaus", required_argument, 0, 'g'},
            {"beamoffset", required_argument, 0, 'o'},

            {"dx1", required_argument, 0, 0},
            {"dx2", required_argument, 0, 0},
            {"dy1", required_argument, 0, 0},
            {"dy2", required_argument, 0, 0},

            {"filename", required_argument, 0, 'f'},
            {"triple", required_argument, 0, 0},

            // Specific options
            {"coord", required_argument, 0, 'c'},      
            {"minh", required_argument, 0, 0},      
            {"maxh", required_argument, 0, 0},      
            {"steph", required_argument, 0, 0},      
            {"min", required_argument, 0, 0},      
            {"max", required_argument, 0, 0},      
            {"step", required_argument, 0, 0},      
            //end of options
            {0, 0, 0, 0}
        };

        int option_index = 0;
        char c = getopt_long (argc, argv, "hx:y:z:s:r:a:c:g:o:f:", long_options, &option_index);

        if (c == -1) break;

        switch (c) {

            //Process options which have a short format
            case 'h':
                usage(argv);
                return 0;
                break;
            case 'x':
                sigma_x_1 = TString(optarg).Atof();
                sigma_x_2 = TString(optarg).Atof();
                break;
            case 'y':
                sigma_y_1 = TString(optarg).Atof();
                sigma_y_2 = TString(optarg).Atof();
                break;
            case 'z':
                sigma_z_1 = TString(optarg).Atof();
                sigma_z_2 = TString(optarg).Atof();
                break;
            case 's':
                sigma_x_1 = TString(optarg).Atof();
                sigma_x_2 = TString(optarg).Atof();
                sigma_y_1 = TString(optarg).Atof();
                sigma_y_2 = TString(optarg).Atof();
                break;
            case 'r':
                rho_xy_1 = TString(optarg).Atof();
                rho_xy_2 = TString(optarg).Atof();
                break;
            case 'a':
                xing_yz = TString(optarg).Atof();
                break;
            case 'g':
                num_gaus = TString(optarg).Atof();
                break;
            case 'o':
                bs_L_xyz_h_par_offset = TString(optarg).Atof();
                break;
            case 'f':
                outputfilename = TString(optarg);
                break;
            case 'c':
                {
                    TString myc(optarg);
                    myc.ToLower();
                    if (myc == "x")
                        bs_xyz_dh_bsep_coord = BeamProfile::cX;
                    else if (myc == "y")
                        bs_xyz_dh_bsep_coord = BeamProfile::cY;
                    else {
                        cerr <<"Invalid -c option, allowed values are 'x' or 'y', is: " << optarg << endl;
                        return 1;
                    }
                    break;
                }
                //Now process long options too
            case 0:
                if (long_options[option_index].name == "sigmaxa1") {
                    cout << "Setting sxa1 = " << optarg << endl;
                    sigma_x_1 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "sigmaya1") {
                    cout << "Setting sya1 = " << optarg << endl;
                    sigma_y_1 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "sigmaza1") {
                    cout << "Setting sza1 = " << optarg << endl;
                    sigma_z_1 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "sigmaxa2") {
                    cout << "Setting sxa2 = " << optarg << endl;
                    sigma_x_2 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "sigmaya2") {
                    cout << "Setting sya2 = " << optarg << endl;
                    sigma_y_2 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "sigmaza2") {
                    cout << "Setting sza2 = " << optarg << endl;
                    sigma_z_2 = TString(optarg).Atof();

                }else if (long_options[option_index].name == "sigmaxb1") {
                    cout << "Setting sxb1 = " << optarg << endl;
                    sigma_x_b_1 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "sigmayb1") {
                    cout << "Setting syb1 = " << optarg << endl;
                    sigma_y_b_1 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "sigmazb1") {
                    cout << "Setting szb1 = " << optarg << endl;
                    sigma_z_b_1 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "sigmaxb2") {
                    cout << "Setting sxb2 = " << optarg << endl;
                    sigma_x_b_2 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "sigmayb2") {
                    cout << "Setting syb2 = " << optarg << endl;
                    sigma_y_b_2 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "sigmazb2") {
                    cout << "Setting szb2 = " << optarg << endl;
                    sigma_z_b_2 = TString(optarg).Atof();

                }else if (long_options[option_index].name == "sigmaxc1") {
                    cout << "Setting sxc1 = " << optarg << endl;
                    sigma_x_c_1 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "sigmayc1") {
                    cout << "Setting syc1 = " << optarg << endl;
                    sigma_y_c_1 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "sigmazc1") {
                    cout << "Setting szc1 = " << optarg << endl;
                    sigma_z_c_1 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "sigmaxc2") {
                    cout << "Setting sxc2 = " << optarg << endl;
                    sigma_x_c_2 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "sigmayc2") {
                    cout << "Setting syc2 = " << optarg << endl;
                    sigma_y_c_2 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "sigmazc2") {
                    cout << "Setting szc2 = " << optarg << endl;
                    sigma_z_c_2 = TString(optarg).Atof();


                } else if (long_options[option_index].name == "rhoa1") {
                    cout << "Setting rhoxya1 = " << optarg << endl;
                    rho_xy_1 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "rhoa2") {
                    cout << "Setting rhoxya2 = " << optarg << endl;
                    rho_xy_2 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "rhob1") {
                    cout << "Setting rhoxyb1 = " << optarg << endl;
                    rho_xy_b_1 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "rhob2") {
                    cout << "Setting rhoxyb2 = " << optarg << endl;
                    rho_xy_b_2 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "rhoc1") {
                    cout << "Setting rhoxyc1 = " << optarg << endl;
                    rho_xy_c_1 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "rhoc2") {
                    cout << "Setting rhoxyc2 = " << optarg << endl;
                    rho_xy_c_2 = TString(optarg).Atof();

                } else if (long_options[option_index].name == "weighta1") {
                    cout << "Setting weighta1 = " << optarg << endl;
                    weight_a_1 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "weighta2") {
                    cout << "Setting weighta2 = " << optarg << endl;
                    weight_a_2 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "weightb1") {
                    cout << "Setting weightb1 = " << optarg << endl;
                    weight_b_1 = TString(optarg).Atof();
                } else if (long_options[option_index].name == "weightb2") {
                    cout << "Setting weightb2 = " << optarg << endl;
                    weight_b_2 = TString(optarg).Atof();

                } else if (long_options[option_index].name == "xingyz") {
                    cout << "Setting xing_yz = " << optarg << endl;
                    xing_yz = TString(optarg).Atof();
                } else if (long_options[option_index].name == "xingxz") {
                    cout << "Setting xing_xz = " << optarg << endl;
                    xing_xz = TString(optarg).Atof();
                } else if (long_options[option_index].name == "numgaus") {
                    num_gaus = TString(optarg).Atof();
                } else if (long_options[option_index].name == "beamoffset") {
                    bs_L_xyz_h_par_offset = TString(optarg).Atof();
                } else if (long_options[option_index].name == "dx1") {
                    bs_L_xyz_h_par_driftx1_gradient = TString(optarg).Atof();
                } else if (long_options[option_index].name == "dx2") {
                    bs_L_xyz_h_par_driftx2_gradient = TString(optarg).Atof();
                } else if (long_options[option_index].name == "dy1") {
                    bs_L_xyz_h_par_drifty1_gradient = TString(optarg).Atof();
                } else if (long_options[option_index].name == "dy2") {
                    bs_L_xyz_h_par_drifty2_gradient = TString(optarg).Atof();

                } else if (long_options[option_index].name == "filename") {
                    outputfilename = TString(optarg);
                } else if (long_options[option_index].name == "triple") {
                    triple = TString(optarg);

                } else if (long_options[option_index].name == "minh") {
                    bs_L_xyz_h_par_min_h = TString(optarg).Atof();
                } else if (long_options[option_index].name == "maxh") {
                    bs_L_xyz_h_par_max_h = TString(optarg).Atof();
                } else if (long_options[option_index].name == "steph") {
                    bs_L_xyz_h_par_step_h = TString(optarg).Atof();
                } else if (long_options[option_index].name == "min") {
                    minLScan = TString(optarg).Atof();
                } else if (long_options[option_index].name == "max") {
                    maxLScan = TString(optarg).Atof();
                } else if (long_options[option_index].name == "step") {
                    stepLScan = TString(optarg).Atof();
                } else {
                    cerr << "ERROR: Invalid option: " << long_options[option_index].name << endl;
                    throw;
                }
                break;
            case '?':
                cerr << "Unhandled option: " << c << endl;
                usage(argv);
                return 1;
            default:
                cerr << "Unhandled option: " << c << endl;
                usage(argv);
                return 1;
        }

    } //process command-line options

    //get action
    TString action;
    if (argc - optind == 1) {
        action = argv[optind];
    } else {
        //invalid action
        cerr << "Provide one and only one action, please. " << endl;
        usage(argv);
        return 1;
    }

    //Init
    SetAtlasStyle();

    //Create the beam profiles, which allows their properties to be set
    InitOverlapBeams();

    if (action == "lspRho") {
        //reset rho to 0, we'll set them later
        rho_xy_1 = rho_xy_2 = 0;
    }

    cout << "Setting B1 sxa,sxb,sya,syb,sza,szb,rhoxya,rhoxyb = " << sigma_x_1 << ", "<< sigma_x_b_1 << ", " << sigma_y_1 << ", " << sigma_y_b_1 << ", "<< sigma_z_1 << ", " << sigma_z_b_1 << ", " << rho_xy_1 << ", " << rho_xy_b_1 << endl;
    cout << "Setting B2 sxa,sxb,sya,syb,sza,szb,rhoxya,rhoxyb = " << sigma_x_2 << ", "<< sigma_x_b_2 << ", " << sigma_y_2 << ", " << sigma_y_b_2 << ", "<< sigma_z_2 << ", " << sigma_z_b_2 << ", " << rho_xy_2 << ", " << rho_xy_b_2 << endl;


    //Set the beam properies with the input from the user

    B1->SetBeamXY(sigma_x_1, sigma_y_1, rho_xy_1); //sigma_x (a), sigma_y (a), rhoa
    B1->SetBeamXY_b(sigma_x_b_1, sigma_y_b_1, rho_xy_b_1); //sigma_x_b, sigma_y_b, rhob
    B1->SetBeamXY_c(sigma_x_c_1, sigma_y_c_1, rho_xy_c_1); //sigma_x_c, sigma_y_c, rhoc

    B2->SetBeamXY(sigma_x_2, sigma_y_2, rho_xy_2); //sigma_x (a), sigma_y (a), rhoa
    B2->SetBeamXY_b(sigma_x_b_2, sigma_y_b_2, rho_xy_b_2); //sigma_x_b, sigma_y_b, rhob
    B2->SetBeamXY_c(sigma_x_c_2, sigma_y_c_2, rho_xy_c_2); //sigma_x_c, sigma_y_c, rhoc

    B1->SetBeamZ(sigma_z_1); //sigma_z
    B1->SetBeamZ_b(sigma_z_b_1); //sigma_z_b
    B1->SetBeamZ_c(sigma_z_c_1); //sigma_z_c

    B2->SetBeamZ(sigma_z_2); //sigma_z
    B2->SetBeamZ_b(sigma_z_b_2); //sigma_z_b
    B2->SetBeamZ_c(sigma_z_c_2); //sigma_z_c

    B1->SetWeightGausA(weight_a_1); //weight_a
    B2->SetWeightGausA(weight_a_2); //weight_a

    if (triple == "no"){
        cout << "triple=======no" << endl;
        weight_b_1 = 1 - weight_a_1;
        weight_b_2 = 1 - weight_a_2;
    }
    B1->SetWeightGausB(weight_b_1); //weight_b
    B2->SetWeightGausB(weight_b_2); //weight_b
    //  cout << "WEIGHT B -------- = " << weight_b_1 << endl;

    B1->SetXingAngleYZ(xing_yz);
    B2->SetXingAngleYZ(xing_yz);
    B1->SetXingAngleXZ(xing_xz);
    B2->SetXingAngleXZ(xing_xz);
    //B1->SetBetaStar(3500.);
    //B2->SetBetaStar(3500.);

    //Decides which program to run based on user input

    if (action == "bs_xyz_h")
        bs_xyz_dh();
    else if (action == "bs_dzdh") 
        bs_dzdh();
    else if (action == "lspRho")
        L_scan(lsp_rho);
    else if (action == "lspX")
        L_scan(lsp_sepX);
    else if (action == "lspY")
        L_scan(lsp_sepY);
    else {
        cerr << "ERROR: Invalid action: " << action << endl;
        usage(argv);
        return 1;
    }


    //close everything
    DestroyOverlapBeams();
    if (debug>= 6)   
        theApp->Run();
}
