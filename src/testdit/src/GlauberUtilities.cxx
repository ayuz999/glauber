//----------------------------------------------------------------------------------------------------
//  Random number generators for
//    - impact parameter (dN/db = b) in 0 < r < 20 fm
//    - radius from several different density profile
//    - phi & theta angles for nucleons
//----------------------------------------------------------------------------------------------------

#include <assert.h>

#include "TClass.h"
#include "TError.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"

#include "GlauberUtilities.h"


namespace GlauberUtilitiesSpace {
  //____________________________________________________________________________________________________
  // Woods-saxon density profile (1D for spherical nuclei)
  Double_t WoodsSaxon(Double_t* x, Double_t* par)
  {
    const Double_t r = x[0] ;
    const Double_t R = par[0] ;
    const Double_t d = par[1] ;

    return r*r/(1.0+TMath::Exp((r-R)/d)) ;
  }

  //____________________________________________________________________________________________________
  // Woods-saxon density profile (2D for deformed nuclei)
  Double_t WoodsSaxon2D(Double_t* x, Double_t* par)
  {
    const Double_t r        = x[0] ;
    const Double_t cosTheta = x[1] ;
    const Double_t R0       = par[0] ;
    const Double_t d        = par[1] ;
    const Double_t beta2    = par[2] ; // Parameter for deformation (beta2=0 --> Standard woods-saxon density)
    const Double_t beta4    = par[3] ; // Parameter for deformation (beta4=0 --> Standard woods-saxon density)
    const Double_t beta6    = 0.0 ;

    const Double_t cosTheta2 = cosTheta * cosTheta ;
    const Double_t cosTheta4 = cosTheta2 * cosTheta2 ;

    // Spherical harmonics (l=2, m=0) Y20(theta) = 1/4 sqrt(5/pi) * (3cos^2(theta) - 1)
    const Double_t Y20  = TMath::Sqrt(5.0/TMath::Pi()) / 4.0 * (3.0 * cosTheta2 - 1.0 ) ;
    const Double_t Y40  = TMath::Sqrt(1.0/TMath::Pi()) * 3.0 / 16.0 * (35.0*cosTheta4 - 30.0*cosTheta2 + 3.0);
    //  const Double_t Y60  = TMath::Sqrt(13.0/TMath::Pi()) / 32.0 * (231.0*cosTheta4*cosTheta2 - 315.0*cosTheta4 + 105.0*cosTheta2 - 5.0) ;
    const Double_t Y60  = 0.0 ;

    const Double_t R = R0 * (1.0 + beta2 * Y20 + beta4 * Y40 + beta6 * Y60 ) ;

    return r*r/(1.0+TMath::Exp((r-R)/d)) ;
  }

  //____________________________________________________________________________________________________
  // Step function
  Double_t StepFunction(Double_t* x, Double_t* par)
  {
    Double_t dx = x[0] ;
    Double_t dy = x[1] ;
    Double_t dz = x[2] ;
    Double_t position = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
    Double_t sigma = par[0] ;

    const Double_t pi = TMath::Pi() ;
    const Double_t r2 = sigma/pi ;
    const Double_t r  = TMath::Sqrt(r2);
    //  return (r - position < 0.0) ? 0.0 : 1.0 ;

    const Double_t V  = 4.0*pi*r2*r/3.0 ;
    return (r - position < 0.0) ? 0.0 : 1.0/V ;
  }

  //____________________________________________________________________________________________________
  // Gaussian function
  Double_t Gaussian(Double_t* x, Double_t* par)
  {
    Double_t dx = x[0] ;
    Double_t dy = x[1] ;
    Double_t dz = x[2] ;
    Double_t r  = TMath::Sqrt(dx*dx + dy*dy + dz*dz) ;
    Double_t sigma2 = par[0]*par[0] ;
    const Double_t norm = 2.0*TMath::Pi()*sigma2 ;

    return 1.0/TMath::Power(norm, 3.0/2.0) * TMath::Exp(-0.5*r*r/sigma2);
  }
}

  GlauberUtilities* GlauberUtilities::mInstance = 0 ;
  const UShort_t GlauberUtilities::mDebugLevel = 20 ;

//____________________________________________________________________________________________________
// Get singleton
GlauberUtilities* GlauberUtilities::instance()
{
  if(!mInstance){
    mInstance = new GlauberUtilities() ;
  }

  return mInstance ;
}

//____________________________________________________________________________________________________
// Default constructor
GlauberUtilities::GlauberUtilities()
{

  /// Define TRandom and seed here, will be used in any classes
  mRandom = new TRandom3() ;
  mRandom->SetSeed(0);
  /// Define impact parameter distribution in 0 < r < 20 fm
  mImpactParameter = new TF1("mImpactParameter", "x", 0, 20.0);
  mDebug = 0 ;
}

//____________________________________________________________________________________________________
// Default destructor
GlauberUtilities::~GlauberUtilities(){}

//____________________________________________________________________________________________________
Double_t GlauberUtilities::GetImpactParameter() const
{
  /// Generate random impact parameter in 0 < r < 20 fm
  /// according to the dN/db = b
  if(!mImpactParameter){
    Error("StGlauberUtilities::GetImpactParameter", "cannot find impact parameter distribution (TF1). abort");
    assert(0);
  }
  const Double_t b = mImpactParameter->GetRandom() ;

  return b ;
}

//____________________________________________________________________________________________________
Double_t GlauberUtilities::GetMaximumImpactParameter() const
{
  return mImpactParameter->GetXmax() ;
}

//____________________________________________________________________________________________________
Double_t GlauberUtilities::GetTheta() const
{
  /// Generate random polar angle in 0 < theta < pi
  /// according to the dN/dtheta = cos(theta)
  const Double_t theta = TMath::ACos(mRandom->Rndm()*2.0-1.0);

  return theta ;
}

//____________________________________________________________________________________________________
Double_t GlauberUtilities::GetPhi() const
{
  /// Generate flat phi angle in -pi < phi < pi
  const Double_t phi = mRandom->Rndm()*TMath::Pi()*2.0 - TMath::Pi() ;

  return phi ;
}

//____________________________________________________________________________________________________
Double_t GlauberUtilities::GetUniform() const
{
  /// Uniform distribution in 0 < x < 1 
  return mRandom->Rndm() ;
}
//____________________________________________________________________________________________________
Double_t GlauberUtilities::GetUniform2() const
{
  /// Uniform distribution in 0 < x < 2
  return mRandom->Uniform(0,2);
}

