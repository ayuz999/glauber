#ifndef __GlauberUtilities_h__
#define __GlauberUtilities_h__

class TF1 ;
class TRandom ;
class TRandom3 ;
#include "Rtypes.h"

//____________________________________________________________________________________________________
// Class GlauberUtilities: Random number, and functions
class GlauberUtilities {
  public:
    static GlauberUtilities* instance();
    virtual ~GlauberUtilities(); /// Default destructor

    /// Get impact parameter
    Double_t GetImpactParameter() const ;

    /// Get maximum impact parameter
    Double_t GetMaximumImpactParameter() const ;

    /// Get theta (polar angle)
    Double_t GetTheta() const ;

    /// Get phi (azimuthal angle)
    Double_t GetPhi() const ;

    /// Get uniform distribution in 0 < x < 1
    Double_t GetUniform() const ;
    Double_t GetUniform2() const ;

    /// Set debug level
    void SetDebug(const UInt_t debug) ;

  private:
    static GlauberUtilities* mInstance ; /// singleton
    GlauberUtilities(); /// Default constructor

    static const UShort_t mDebugLevel ; /// Debug level for StGlauberUtilities
    TF1* mImpactParameter ; /// impact parameter

   // TRandom* mRandom ; /// Random numbers
    TRandom3* mRandom ; /// Random numbers
    UInt_t mDebug ; /// Debug level

};

namespace GlauberUtilitiesSpace {
  // Woods-saxon density profile
  Double_t WoodsSaxon(Double_t* x, Double_t* par) ;

  // Woods-saxon density profile (2D for deformed nuclei)
  Double_t WoodsSaxon2D(Double_t* x, Double_t* par) ;

  // Step function
  Double_t StepFunction(Double_t* x, Double_t* par) ;

  // Gaussian function
  Double_t Gaussian(Double_t* x, Double_t* par) ;
}


#endif

