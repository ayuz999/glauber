/******************************************************************************
 * $Id: StGlauberConstUtilities.cxx,v 1.2 2012/04/25 04:55:04 hmasui Exp $
 * $Log: StGlauberConstUtilities.cxx,v $
 * Revision 1.2  2012/04/25 04:55:04  hmasui
 * Expand centrality bins
 *
******************************************************************************/

#include <assert.h>
#include "TError.h"
#include "GlauberConstUtilities.h"

//____________________________________________________________________________________________________
namespace GlauberConstUtilities {
  // Impact parameter
  const UInt_t mBBin   = 400 ;
  const Double_t mBMax = 20.0 ;

  // Npart
  const UInt_t mNpartBin = 500 ;
  const Double_t mNpartMax = 500 ;

  // Ncoll
  const UInt_t mNcollBin = 1800 ;
  const Double_t mNcollMax = 1800 ;

  // Multiplicity
  const UInt_t mMultiplicityBin = 2000 ;
  const Double_t mMultiplicityMax = 2000 ;

  // Centrality
  const UInt_t mCentralityBin     = 36 ;
  const Double_t mCentralityMin[] = { 0,  5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75,  80,
     0, 10, 20, 30, 40, 50, 60, 70, 10, 40,  0, 20, 40, 60, 20, 50,  0,  0, 7.5 };
  const Double_t mCentralityMax[] = { 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 100,
    10, 20, 30, 40, 50, 60, 70, 80, 40, 80, 20, 40, 60, 80, 50, 80, 60, 80, 30 };
};

//____________________________________________________________________________________________________
const UInt_t GlauberConstUtilities::GetImpactParameterBin()
{
  return mBBin ;
}

//____________________________________________________________________________________________________
const Double_t GlauberConstUtilities::GetImpactParameterMax()
{
  return mBMax ;
}

//____________________________________________________________________________________________________
const UInt_t GlauberConstUtilities::GetNpartBin()
{
  return mNpartBin ;
}

//____________________________________________________________________________________________________
const Double_t GlauberConstUtilities::GetNpartMax()
{
  return mNpartMax ;
}

//____________________________________________________________________________________________________
const UInt_t GlauberConstUtilities::GetNcollBin()
{
  return mNcollBin ;
}

//____________________________________________________________________________________________________
const Double_t GlauberConstUtilities::GetNcollMax()
{
  return mNcollMax ;
}

//____________________________________________________________________________________________________
const UInt_t GlauberConstUtilities::GetMultiplicityBin()
{
  return mMultiplicityBin ;
}

//____________________________________________________________________________________________________
const Double_t GlauberConstUtilities::GetMultiplicityMax()
{
  return mMultiplicityMax ;
}

//____________________________________________________________________________________________________
const UInt_t GlauberConstUtilities::GetCentralityBin()
{
  return mCentralityBin ;
}

//____________________________________________________________________________________________________
const Double_t GlauberConstUtilities::GetCentralityMin(const UInt_t icent)
{
  if( icent >= mCentralityBin ){
    Error("StGlauberConstUtilities::GetCentralityMin", "Unknown centrality id, icent=%3d", icent);
    assert(0);
  }

  return mCentralityMin[icent] ;
}

//____________________________________________________________________________________________________
const Double_t GlauberConstUtilities::GetCentralityMax(const UInt_t icent)
{
  if( icent >= mCentralityBin ){
    Error("StGlauberConstUtilities::GetCentralityMax", "Unknown centrality id, icent=%3d", icent);
    assert(0);
  }

  return mCentralityMax[icent] ;
}


//____________________________________________________________________________________________________
const Bool_t GlauberConstUtilities::IsCentralityOk(const UInt_t icent, const Double_t centrality)
{
  if(icent>=mCentralityBin){
    Error("StGlauberConstUtilities::IsCentralityOk", "Unknown centrality bin, icent=%3d", icent);
    assert(0);
  }

  return (centrality >= mCentralityMin[icent] && centrality < mCentralityMax[icent] );
}

