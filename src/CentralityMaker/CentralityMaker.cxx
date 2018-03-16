/******************************************************************************
 * $Id: StCentralityMaker.cxx,v 1.4 2012/04/25 05:17:24 hmasui Exp $
 * $Log: StCentralityMaker.cxx,v $
 * Revision 1.4  2012/04/25 05:17:24  hmasui
 * Remove SetModeNBD() and corresponding data member, merged into GetNegativeBinomial() function
 *
******************************************************************************/

#include <assert.h>

#include "TError.h"

#include "NegativeBinomial.h"
#include "Centrality.h"
#include "CentralityMaker.h"

CentralityMaker::CentralityMaker()
  : mNBinomial(0), mCentrality(0) {}

//____________________________________________________________________________________________________
// Initialization inside the constructor
CentralityMaker::CentralityMaker(const Char_t* system)
  : mNBinomial(0), mCentrality(0)
{
  Init(system);
}


//____________________________________________________________________________________________________
const NegativeBinomial* CentralityMaker::GetNegativeBinomial(const UInt_t id) const
{
  if( id >= mNBinomial.size() ){
    Error("StCentralityMaker::GetNegativeBinomial", "Invalid index for negative binomial, id=%3d", id);
    assert(0);
  }
  return mNBinomial[id] ;
}

//____________________________________________________________________________________________________
Centrality* CentralityMaker::GetCentrality(const UInt_t id) const
{
  if( id >= mCentrality.size() ){
    Error("StCentralityMaker::GetCentrality", "Invalid index for centrality, id=%3d", id);
    assert(0);
  }

  return mCentrality[id] ;
}

//____________________________________________________________________________________________________
// Use this function to initialize centrality bin and NBD
void CentralityMaker::Init(const Char_t* system)
{
  /// Initialize centrality and NBD
  const TString systemName(system);
  /// Initialize centrality and NBD parameters
  const TString type[] = {"default", "low", "high"};
  for(UInt_t id=0; id<3; id++){
    mCentrality.push_back( new Centrality(systemName, type[id]) );
  }

  /// Initialize NBD utility
  // Use multiplicity depdent efficiency for all systems
  for(UInt_t id=0; id<3; id++){
    Bool_t isConstEfficiency = kFALSE ;

    mNBinomial.push_back( 
        new NegativeBinomial( mCentrality[id]->GetNpp(), mCentrality[id]->GetK(), mCentrality[id]->GetX(),
        mCentrality[id]->GetEfficiency(), mCentrality[id]->GetTriggerBias(), isConstEfficiency)
        );
  }
}

