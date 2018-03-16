
#ifndef __Centrality_h__
#define __Centrality_h__

#include <vector>
#include "TString.h"

//____________________________________________________________________________________________________
// Class StCentrality: Centrality class
//   - Centrality: multiplicity cut, and corresponding centrality range
//   - NBD: npp, k and x parameters
class Centrality {
  public:
    /// Type         description
    /// default      default parameters for a given system
    /// low          low npp, high x (-5% npp)
    /// high         high npp, low x (+5% npp)
    Centrality(const TString system = "AuAu_200GeV", const TString type="default"); /// Default constructor
    virtual ~Centrality(); /// Default destructor

    // Getters

    // mode
    //  0    default
    //  1    -5% total cross section
    //  2    +5% total cross section
    Double_t GetCentrality(const UInt_t multiplicity, const UInt_t mode=0) const ; // Get centrality

    Double_t GetNpp()         const ; /// Get Npp
    Double_t GetK()           const ; /// Get K
    Double_t GetX()           const ; /// Get X
    Double_t GetEfficiency()  const ; /// Get efficiency
    Double_t GetTriggerBias() const ; /// Get trigger bias

    Double_t GetReweighting(const UInt_t multiplicity) const ; /// Re-weighting correction

  private:
    // Functions
    Double_t GetNppError(const Double_t npp) const ; /// +/- 5% error if type is low or high

    // Initialization
    void Init_AuAu54GeV() ; /// Initialization of Au + Au 54 GeV

    // Data members
    enum {
      mNMode = 3
    };

    const TString mType   ; /// Type
    Double_t mNpp         ; /// Average multiplicity in p + p
    Double_t mK           ; /// NBD k parameters
    Double_t mX           ; /// Fraction of hard component
    Double_t mEfficiency  ; /// Efficiency
    Double_t mTriggerBias ; /// Trigger bias
    Double_t mParReweighting[2]; /// Parameters for re-weighting correction 1-exp(-p0*x^{p1})
    std::vector<UInt_t> mMultiplicityCut[mNMode] ; /// Multiplicity cut for each centrality bin
    std::vector<Double_t> mCentralityMin[mNMode] ; /// Minimum centrality
    std::vector<Double_t> mCentralityMax[mNMode] ; /// Maximum centrality

};

inline Double_t Centrality::GetNpp()          const { return mNpp ; }
inline Double_t Centrality::GetK()            const { return mK ; }
inline Double_t Centrality::GetX()            const { return mX ; }
inline Double_t Centrality::GetEfficiency()   const { return mEfficiency ; }
inline Double_t Centrality::GetTriggerBias()  const { return mTriggerBias ; }

#endif
