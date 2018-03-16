
#include <assert.h>

#include "TError.h"
#include "TMath.h"

#include "Centrality.h"


using std::vector ;

//____________________________________________________________________________________________________
// Re-weighting correction
Double_t Reweighting(const Double_t* x, const Double_t* par)
{
  return 1.0 - TMath::Exp(-par[0]*TMath::Power(x[0], par[1]));
}

//____________________________________________________________________________________________________
// Default constructor
Centrality::Centrality(const TString system, const TString type)
  : mType(type), mNpp(0.0), mK(0.0), mX(0.0)
{
  for(UInt_t im=0; im<mNMode; im++){
    mMultiplicityCut[im].clear();
    mCentralityMin[im].clear();
    mCentralityMax[im].clear();
  }

  mParReweighting[0] = 0.0 ;
  mParReweighting[1] = 0.0 ;

  /// Initialize centrality bin
  if( system.CompareTo("auau_54gev", TString::kIgnoreCase) == 0 )      Init_AuAu54GeV() ;
  else{
    Error("Centrality", "can't find system %s. See below for current available systems", system.Data());
    assert(0);
  }

  /// Make sure centrality bin size is equal (and not empty)
  if( mMultiplicityCut[0].empty() ) Warning("Centrality", "Bins for multiplicity cut is empty. Program will be crashed ...");
  if( mCentralityMin[0].empty() )   Warning("Centrality", "Minimum centrality is empty. Program will be crashed ...");
  if( mCentralityMax[0].empty() )   Warning("Centrality", "Maximum centrality is empty. Program will be crashed ...");

  // Cuts
  const Bool_t isBinSizeEqual0 =
    (mMultiplicityCut[0].size() == mMultiplicityCut[1].size())
    && (mMultiplicityCut[1].size() == mMultiplicityCut[2].size())
    && (mMultiplicityCut[2].size() == mMultiplicityCut[0].size())
    ;
  if(!isBinSizeEqual0){
    Error("Centrality", Form("Bin size for arrays are not equal: (default, small, large) = (%3d, %3d, %3d)",
          mMultiplicityCut[0].size(), mMultiplicityCut[1].size(), mMultiplicityCut[2].size())
        );
    assert(isBinSizeEqual0);
  }

  // Among cut, min, max
  const Bool_t isBinSizeEqual1 =
    (mMultiplicityCut[0].size() == mCentralityMin[0].size())
    && (mCentralityMin[0].size() == mCentralityMax[0].size())
    && (mCentralityMax[0].size() == mMultiplicityCut[0].size())
    ;
  if(!isBinSizeEqual1){
    Error("StCentrality", 
        Form("Bin size for arrays are not equal: (mult, centmin, centmax) = (%3d, %3d, %3d)",
          mMultiplicityCut[0].size(), mCentralityMin[0].size(), mCentralityMax[0].size())
        );
    assert(isBinSizeEqual1);
  }


}

//____________________________________________________________________________________________________
// Default destructor
Centrality::~Centrality()
{
  for(UInt_t im=0; im<mNMode; im++){
    mMultiplicityCut[im].clear();
    mCentralityMin[im].clear();
    mCentralityMax[im].clear();
  }
}


//____________________________________________________________________________________________________
Double_t Centrality::GetCentrality(const UInt_t multiplicity, const UInt_t mode) const
{
  /// Get centrality bin from multiplicity
  ///  mode  
  ///   0        default centrality bins
  ///   1        5% lower centrality bins (4.75, 9.5, 19.5, ...)
  ///   2        5% higher centrality bins

  // Check mode
  if( mode >= mNMode ){
    Error("StCentrality::GetCentrality", "Invalid mode, mode=%3d. reutnr -1", mode);
    return -1.0 ;
  }

  // Make sure array has some entries
  if( mMultiplicityCut[mode].empty() ){
    Error("StCentrality::GetCentrality", "multiplicity cuts for centrality have not implemented. return -1");
    return -1.0 ;
  }

  // Find centrality bin
  for(UInt_t icent=0; icent<mMultiplicityCut[mode].size(); icent++){
    if( multiplicity > mMultiplicityCut[mode][icent] ) 
      return (mCentralityMin[mode][icent]+mCentralityMax[mode][icent])/2.0;
  }

  return -1.0 ;
}

//____________________________________________________________________________________________________
Double_t Centrality::GetNppError(const Double_t npp) const
{
  if ( mType.CompareTo("default", TString::kIgnoreCase) == 0 )       return 0.0 ;
  else if ( mType.CompareTo("low", TString::kIgnoreCase) == 0 ){
    return -0.05 * npp ;
  }
  else if ( mType.CompareTo("high", TString::kIgnoreCase) == 0 ){
    return 0.05 * npp ;
  }
  else {
    Warning("StCentrality::GetNppError", "Unknown type, type=%s. Set error = 0", mType.Data());
    return 0.0 ;
  }
}

//____________________________________________________________________________________________________
void Centrality::Init_AuAu54GeV()
{
  /// Initialize Au + Au collisiona at 200 GeV (Run10)
  //----------------------------------------------------------------------------------------------------
  //  Update on Feb/04/2011
  //    npp and x are determined from published spectra paper, PRC79, 034909 (2009)
  //    npp = 2.43
  //    k = 2.0 (probably k dependence is weak, need to study)
  //    x = 0.13 from interpolation between 19.6 and 200 GeV data (PHOBOS, PRC70, 021902, 2004)
  //
  //    - Multiplicity dependent efficiency with d=0.14
  //    - Additional constant efficiency loss 17% (0.17)
  //      (need to be implemented somehow, this is currently done by temporary fix in StNegativeBinomial.cxx)
  //
  //      NOTE: parameters are not final yet (Feb/04/2011)
  //----------------------------------------------------------------------------------------------------
  //  0.00 -  5.00 (%) :  M >  441
  //  5.00 - 10.00 (%) :  M >  375
  // 10.00 - 15.00 (%) :  M >  317
  // 15.00 - 20.00 (%) :  M >  266
  // 20.00 - 25.00 (%) :  M >  221
  // 25.00 - 30.00 (%) :  M >  182
  // 30.00 - 35.00 (%) :  M >  148
  // 35.00 - 40.00 (%) :  M >  118
  // 40.00 - 45.00 (%) :  M >   93
  // 45.00 - 50.00 (%) :  M >   72
  // 50.00 - 55.00 (%) :  M >   55
  // 55.00 - 60.00 (%) :  M >   41
  // 60.00 - 65.00 (%) :  M >   30
  // 65.00 - 70.00 (%) :  M >   21
  // 70.00 - 75.00 (%) :  M >   15
  // 75.00 - 80.00 (%) :  M >   10
  // 80.00 - 85.00 (%) :  M >    6
  // 85.00 - 90.00 (%) :  M >    4
  //

  // Error on x (energy specific ???)
  // NOTE: type 'low' and 'high' are for npp, x is anti-correlated with npp
  Double_t xError = 0.0 ; // absolute error
  if( mType.CompareTo("low", TString::kIgnoreCase) == 0 )       xError = 0.02 ; // npp is low, x is high
  else if( mType.CompareTo("high", TString::kIgnoreCase) == 0 ) xError = -0.02 ; // npp is high, x is low

  const Double_t npp = 2.10 ; // default npp
  mNpp = npp + GetNppError(npp) ; // 5% error on npp
  mK   = 2.00 ;
  mX   = 0.11 + xError ;
  mEfficiency  = 0.14 ;
  mTriggerBias = 1.00 ;
  //  mTriggerBias = 0.835 ; // Run10 200 GeV

  // Final centrality bins
  // Define multiplicity cut from central to peripheral
  mMultiplicityCut[0].push_back( 443 ); mCentralityMin[0].push_back( 0.0 ); mCentralityMax[0].push_back( 5.0 );
  mMultiplicityCut[0].push_back( 376 ); mCentralityMin[0].push_back( 5.0 ); mCentralityMax[0].push_back( 10.0 );
  mMultiplicityCut[0].push_back( 317 ); mCentralityMin[0].push_back( 10.0 ); mCentralityMax[0].push_back( 15.0 );
  mMultiplicityCut[0].push_back( 266 ); mCentralityMin[0].push_back( 15.0 ); mCentralityMax[0].push_back( 20.0 );
  mMultiplicityCut[0].push_back( 221 ); mCentralityMin[0].push_back( 20.0 ); mCentralityMax[0].push_back( 25.0 );
  mMultiplicityCut[0].push_back( 182 ); mCentralityMin[0].push_back( 25.0 ); mCentralityMax[0].push_back( 30.0 );
  mMultiplicityCut[0].push_back( 148 ); mCentralityMin[0].push_back( 30.0 ); mCentralityMax[0].push_back( 35.0 );
  mMultiplicityCut[0].push_back( 119 ); mCentralityMin[0].push_back( 35.0 ); mCentralityMax[0].push_back( 40.0 );
  mMultiplicityCut[0].push_back(  94 ); mCentralityMin[0].push_back( 40.0 ); mCentralityMax[0].push_back( 45.0 );
  mMultiplicityCut[0].push_back(  73 ); mCentralityMin[0].push_back( 45.0 ); mCentralityMax[0].push_back( 50.0 );
  mMultiplicityCut[0].push_back(  56 ); mCentralityMin[0].push_back( 50.0 ); mCentralityMax[0].push_back( 55.0 );
  mMultiplicityCut[0].push_back(  41 ); mCentralityMin[0].push_back( 55.0 ); mCentralityMax[0].push_back( 60.0 );
  mMultiplicityCut[0].push_back(  30 ); mCentralityMin[0].push_back( 60.0 ); mCentralityMax[0].push_back( 65.0 );
  mMultiplicityCut[0].push_back(  21 ); mCentralityMin[0].push_back( 65.0 ); mCentralityMax[0].push_back( 70.0 );
  mMultiplicityCut[0].push_back(  15 ); mCentralityMin[0].push_back( 70.0 ); mCentralityMax[0].push_back( 75.0 );
  mMultiplicityCut[0].push_back(  10 ); mCentralityMin[0].push_back( 75.0 ); mCentralityMax[0].push_back( 80.0 );
  // -5%                                 +5%
  mMultiplicityCut[1].push_back( 447 );  mMultiplicityCut[2].push_back( 440 );
  mMultiplicityCut[1].push_back( 382 );  mMultiplicityCut[2].push_back( 370 );
  mMultiplicityCut[1].push_back( 325 );  mMultiplicityCut[2].push_back( 309 );
  mMultiplicityCut[1].push_back( 275 );  mMultiplicityCut[2].push_back( 256 );
  mMultiplicityCut[1].push_back( 231 );  mMultiplicityCut[2].push_back( 211 );
  mMultiplicityCut[1].push_back( 193 );  mMultiplicityCut[2].push_back( 171 );
  mMultiplicityCut[1].push_back( 159 );  mMultiplicityCut[2].push_back( 137 );
  mMultiplicityCut[1].push_back( 130 );  mMultiplicityCut[2].push_back( 108 );
  mMultiplicityCut[1].push_back( 104 );  mMultiplicityCut[2].push_back(  84 );
  mMultiplicityCut[1].push_back(  83 );  mMultiplicityCut[2].push_back(  64 );
  mMultiplicityCut[1].push_back(  65 );  mMultiplicityCut[2].push_back(  47 );
  mMultiplicityCut[1].push_back(  50 );  mMultiplicityCut[2].push_back(  34 );
  mMultiplicityCut[1].push_back(  37 );  mMultiplicityCut[2].push_back(  24 );
  mMultiplicityCut[1].push_back(  27 );  mMultiplicityCut[2].push_back(  17 );
  mMultiplicityCut[1].push_back(  20 );  mMultiplicityCut[2].push_back(  11 );
  mMultiplicityCut[1].push_back(  14 );  mMultiplicityCut[2].push_back(   7 );
          

	// Set same centrality bins 
	for(UInt_t ic=0; ic<mCentralityMin[0].size(); ic++){
		for(UInt_t im=1; im<mNMode; im++){
			mCentralityMin[im].push_back( mCentralityMin[0][ic] );
			mCentralityMax[im].push_back( mCentralityMax[0][ic] );
		}
	}
}

//____________________________________________________________________________________________________
Double_t Centrality::GetReweighting(const UInt_t multiplicity) const
{
	// Get re-weighting correction
	// parameters depend on the incident energies
	// Check parameters. If all parameters are 0, return 1.0
	if( mParReweighting[0] == 0.0 && mParReweighting[1] == 0.0 ){
		Warning("StCentrality::GetReweighting", "No re-weighting correction implemented. return 1");
		return 1.0 ;
	}
	const Double_t x[] = {multiplicity};
	return Reweighting(x, mParReweighting) ;
}

