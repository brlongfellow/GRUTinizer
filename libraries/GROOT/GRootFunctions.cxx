#include "GRootFunctions.h"
#include <iostream>
#include "TF1.h"


NamespaceImp(GRootFunctions)


#define PI TMATH::Pi()

Double_t GRootFunctions::PolyBg(Double_t *dim, Double_t *par, Int_t order) {
  Double_t result = 0.0;
  int j=0;
  for(Int_t i=0;i<=order;i++) {
    result += *(par+j) *TMath::Power(dim[0],i);
    j++;
  }
    //result += par[i]*TMath::Power(dim[0]-par[order+1],i);
  return result;
}

Double_t GRootFunctions::LinFit(Double_t *dim, Double_t *par) {
  return PolyBg(dim,par,1);
}

//Double_t GRootFunctions::LinearBG(Double_t *dim,Double_*par) {
//  //  -dim[0]: channels to fit
//  //  -par[0]: offset
//  //  -par[1]: slope
//  //  -par[2]: begin of exclude; should be fixed.
//  //  -par[3]: end of exclude; should be fixed.
// 
//  if(!std::is_nan(par[2]) && !std::is_nan(par[3])) {
//    if(par[2]>par[3])
//      std::swap(par[2],par[3]);
//    if(x[0]>par[2] && x[0]<par[3]) {
//      TF1::RejectPoint();
//      return 0;
//    }
//  }
//  return par[0] + par[1]*x[0];
//}


Double_t GRootFunctions::QuadFit(Double_t *dim, Double_t *par) {
  return PolyBg(dim,par,2);
}

Double_t GRootFunctions::StepFunction(Double_t *dim, Double_t *par) {
  //  -dim[0]: channels to fit
  //  -par[0]: height of peak
  //  -par[1]: centroid of peak
  //  -par[2]: sigma of peak
  //  -par[3]: size of step in step function.

  Double_t x       = dim[0];
  
  Double_t height  = par[0];
  Double_t cent    = par[1];
  Double_t sigma   = par[2];
  //Double_t R       = par[4];
  Double_t step    = par[3];

  //return TMath::Abs(step)*height/100.0*TMath::Erfc((x-cent)/(TMath::Sqrt(2.0)*sigma));
  return height*(step/100.0) *TMath::Erfc((x-cent)/(TMath::Sqrt(2.0)*sigma));

}

Double_t GRootFunctions::StepBG(Double_t *dim, Double_t *par) {
  return StepFunction(dim,par) + PolyBg(dim,(par+4),0);
}

Double_t GRootFunctions::Gaus(Double_t *dim, Double_t *par) {
  // - dim[0]: channels to fit
  // - par[0]: height of peak
  // - par[1]: cent of peak
  // - par[2]: sigma
  // - par[3]: relative height of skewed gaus to gaus

  Double_t x      = dim[0];
  Double_t height = par[0];
  Double_t cent   = par[1];
  Double_t sigma  = par[2];
  Double_t R      = par[3];

  return height*(1.0-R/100.0)*TMath::Gaus(x,cent,sigma);
}

Double_t GRootFunctions::SkewedGaus(Double_t *dim,Double_t *par) {

  // StepFunction(dim,par) + PolyBg
  // - par[0]: height of peak
  // - par[1]: cent of peak
  // - par[2]: sigma
  // - par[3]: relative height of skewed gaus to gaus
  // - par[4]: "skewedness" of the skewed gaussin

  Double_t x      = dim[0]; //channel number used for fitting
  Double_t height = par[0]; //height of photopeak
  Double_t cent   = par[1]; //Peak Centroid of non skew gaus
  Double_t sigma  = par[2]; //standard deviation of gaussian
  Double_t R      = par[3]; //relative height of skewed gaussian
  Double_t beta   = par[4]; //"skewedness" of the skewed gaussian

  double scaling = R*height/100.0;
  //double x_rel = (x - cent)/sigma;

  double fterm = (x-cent)/(sigma*TMath::Sqrt(2.));
  double sterm = sigma /  (beta *TMath::Sqrt(2.));

  return scaling * TMath::Exp((x-cent)/beta) * TMath::Erfc(fterm + sterm); 
}

Double_t GRootFunctions::PhotoPeak(Double_t *dim,Double_t *par) {
  return Gaus(dim,par) + SkewedGaus(dim,par);
}

Double_t GRootFunctions::PhotoPeakBG(Double_t *dim,Double_t *par) {
  // - dim[0]: channels to fit
  // - par[0]: height of peak
  // - par[1]: cent of peak
  // - par[2]: sigma
  // - par[3]: relative height of skewed gaus to gaus
  // - par[4]: "skewedness" of the skewed gaussin
  // - par[5]: size of stepfunction step.

  // - par[6]: base bg height.
  // - par[7]: slope of bg.
  
  double spar[4];
  spar[0] = par[0];
  spar[1] = par[1];
  spar[2] = par[2];
  spar[3] = par[5];  //stepsize;
  return Gaus(dim,par) + SkewedGaus(dim,par) + StepFunction(dim,spar) + PolyBg(dim,par+6,0);
}

Double_t GRootFunctions::PhotoPeakBGExcludeRegion(Double_t *dim,Double_t *par) {
  // - dim[0]: channels to fit
  // - par[0]: height of peak
  // - par[1]: cent of peak
  // - par[2]: sigma
  // - par[3]: relative height of skewed gaus to gaus
  // - par[4]: "skewedness" of the skewed gaussin
  // - par[5]: size of stepfunction step.
  
  // - par[6]: base bg height.
  
  // - par[7]: exclude low;
  // - par[8]: exclude high;

  if(dim[0]>par[7] && dim[0]<par[8]) {
    TF1::RejectPoint();
    return 0;
  }
  double spar[4];
  spar[0] = par[0];
  spar[1] = par[1];
  spar[2] = par[2];
  spar[3] = par[5];  //stepsize;
  return Gaus(dim,par) + SkewedGaus(dim,par) + StepFunction(dim,spar) + PolyBg(dim,par+6,0);
}


// For fitting Ge detector efficiencies.
Double_t GRootFunctions::Efficiency(Double_t *dim, Double_t *par){
  // - dim[0]: energy.
  // - par[0]: zeroth order
  // - par[1]: first order
  // - par[2]: second order
  // - par[3]: inverse energy squared term.
  // - Formula : 10**(0+1*Log(x)+2*Log(x)**2+3/x**2)

  Double_t x  = dim[0];
  Double_t p0 = par[0];
  Double_t p1 = par[1];
  Double_t p2 = par[2];
  Double_t p3 = par[3];

  if(x!=0)
    return pow(10.0,(p0+p1*TMath::Log10(x)+p2*std::pow(TMath::Log10(x),2.0)+p3/(std::pow(x,2.0))));
  else
    return 0;

}


Double_t GRootFunctions::GausExpo(Double_t *x, Double_t *pars) {

  double result;

  // gaus + step*expo conv with a gaus.

  // par[0] = height
  // par[1] = cent
  // par[2] = sigma
  // par[3] = step
  // par[4] = decay parameter

  //result = pars[0]*TMath::Gaus(x[0],pars[1],pars[2])+ pars[0]*(double)(x[0]>pars[1])*TMath::Exp(-pars[3]*x[0]);

  
  double height = pars[0];
  double mu = pars[1];
  double sigma = pars[2];
  double t = pars[3];
  double k = pars[4];

  
  double output = 0;
  const int n = 200;
  const double width = 5*sigma;
  for(double y=0; y<n; y++) {
    double gauss_x = mu - (y-n/2)*width/n;
    double expo_x = x[0] - gauss_x;
    
    double expo = (expo_x>t) ? TMath::Exp(-k*(expo_x-t)) : 0;
    double gauss = TMath::Gaus(x[0], mu, sigma);
    output += height*expo*gauss;
  }
  return output;


  // double p1 = height;
  // double numerator = (-2*k*x[0] + 2*k*mu + sigma*sigma);
  // double denominator = (2.0*k*k);
  // double p2 = TMath::Exp(numerator/denominator);
  // double p3 = TMath::Sqrt(TMath::Pi()/2) * sigma;
  // double p4 = (1-TMath::Erf(((k*(t-x[0]+mu)+sigma*sigma))/(TMath::Sqrt(2)*k*sigma)));

  // result = p1*p2*p3*p4;  
  
  // return result;
}





Double_t GRootFunctions::LanGaus(Double_t *x, Double_t *pars){
   double dy, y, conv, spec, gaus;
   conv = 0;

   for(int i=0; i<10; i++){
    dy = 5*pars[2]/10.0; // truncate the convolution by decreasing number of evaluation points and decreasing range [2.5 sigma still covers 98.8% of gaussian]
    y = x[0]-2.5*pars[2]+dy*i;
    spec = pars[0]+pars[1]*y; // define background SHOULD THIS BE CONVOLUTED ????? *************************************
    //for( int n=0; n<(int)(pars[0]+0.5); n++) // the implementation of landau function should be done using the landau function
      spec +=pars[3]*TMath::Landau(-y,-pars[4],pars[5])/TMath::Landau(0,0,100); // add peaks, dividing by max height of landau
    gaus = TMath::Gaus(-x[0],-y,pars[2])/sqrt(2*TMath::Pi()*pars[2]*pars[2]); // gaus must be normalisd so there is no sigma weighting
    conv += gaus*spec*dy; // now convolve this [integrate the product] with a gaussian centered at x;
  }

  return conv;
}


Double_t GRootFunctions::LanGausHighRes(Double_t *x, Double_t *pars){ // 5x more convolution points with 1.6x larger range
  double dy, y, conv, spec, gaus;
  conv = 0;

  for(int i=0; i<50; i++){
    dy = 8*pars[2]/50.0; // 4 sigma covers 99.99% of gaussian
    y  = x[0]-4*pars[2]+dy*i;

    spec = pars[0]+pars[1]*y;
    //for( int n=0; n<(int)(pars[0]+0.5); n++)
    spec +=pars[3]*TMath::Landau(-y,-pars[4],pars[5])/TMath::Landau(0,0,100);

    gaus = TMath::Gaus(-x[0],-y,pars[2])/sqrt(2*TMath::Pi()*pars[2]*pars[2]);
    conv += gaus*spec*dy;
  }
  return conv;
}


Double_t GRootFunctions::GammaEff(Double_t *x,Double_t *par) {
  // LOG(EFF) = A0 + A1*LOG(E) + A2*LOG(E)^2 + A3/E^2 

  double logE = TMath::Log10(x[0]);
  double temp =  par[0] + par[1]*logE + par[2]*logE*logE +par[3]/(x[0]*x[0]);
  return pow(10,temp);

}


Double_t GRootFunctions::ComptonFormula(Double_t *x,Double_t *par) {

  //par[0] = inital gamma energy; in keV.

  double lower = 1 + (par[0]/511.)*(1-TMath::Cos(TMath::DegToRad()*x[0]));
  return par[0]/lower;

}

Double_t GRootFunctions::AnalyzingPower(Double_t *x,Double_t *par) {

  //par[0] = inital gamma energy; in keV.
  double scattered = ComptonFormula(x,par);

  double sin2 = (TMath::Sin(TMath::DegToRad()*x[0]))*(TMath::Sin(TMath::DegToRad()*x[0]));
  double lower = scattered/par[0] + par[0]/scattered  - sin2;
  return sin2/lower;
}



