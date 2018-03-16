/******************************************************************************
 * $Id: StGlauberHistogramMaker.cxx,v 1.2 2012/04/25 05:03:06 hmasui Exp $
 * $Log: StGlauberHistogramMaker.cxx,v $
 * Revision 1.2  2012/04/25 05:03:06  hmasui
 * Use StCentralityMaker. Added weight for 2D fill. Use STAR logger instead of iostream
 *
******************************************************************************/

#include <assert.h>
#include <fstream>

#include "TError.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"


//#include "CentralityMaker/NegativeBinomial.h"
#include "CentralityMaker/CentralityMaker.h"
#include "CentralityMaker/Centrality.h"
#include "GlauberConstUtilities.h"
#include "GlauberTree/GlauberTree.h"
#include "GlauberHistogramMaker.h"


using std::ofstream ;
using std::endl ;
using std::vector ;

  const TString GlauberHistogramMaker::mXAxisName[]  = { "b", "Npart", "mult", "cent" };
  const TString GlauberHistogramMaker::mXAxisTitle[] = { "impact parameter (fm)", "N_{part}", "multiplicity", "centrality (%)" };
  const Int_t GlauberHistogramMaker::mXAxisBin[]     = { 
    400,  // Impact parameter bin up to 20 fm, 0.2 fm increment
    100,  // Npart bin (5 Npart per bin)
    100, // Multiplicity bin (20 multiplicity per bin)
    GlauberConstUtilities::GetCentralityBin() // Centrality bin
  };

  /// Define maximum (minimum is always 0)
  const Double_t GlauberHistogramMaker::mXAxisMax[] = {
    20.0,     // Maximum impact parameter
    500,      // Maximum Npart
    2000,     // Maximum multiplicity
    GlauberConstUtilities::GetCentralityBin() // Maximum centrality
  };

//____________________________________________________________________________________________________
// Default constructor
GlauberHistogramMaker::GlauberHistogramMaker()
  : mName("hTest"), mTitle("Test"), mYTitle("Test"), mYbin(1), mYmin(0), mYmax(1), mIsUnitWeight(kTRUE)
{
  Init();
}


//____________________________________________________________________________________________________
GlauberHistogramMaker::GlauberHistogramMaker(const TString name, const TString title, const TString ytitle,
    const Int_t ybin, const Double_t ymin, const Double_t ymax, const Bool_t isUnitWeight)
  : mName(name), mTitle(title), mYTitle(ytitle), mYbin(ybin), mYmin(ymin), mYmax(ymax), mIsUnitWeight(isUnitWeight)
{
  Init();
}

//____________________________________________________________________________________________________
// Default destructor
GlauberHistogramMaker::~GlauberHistogramMaker()
{
  Reset() ;
}

//____________________________________________________________________________________________________
void GlauberHistogramMaker::SetTableDirectory(const TString directory)
{
  mTableDirectory = directory ;
}

//____________________________________________________________________________________________________
void GlauberHistogramMaker::Init()
{
  mDebug = 0 ;

  /// Initialize histograms and x-axis
  Reset() ;

  /// Initialize x-axis
  for(UInt_t ix=0; ix<mNXaxis; ix++){
    mXaxis.push_back( -9999.0 );
  }

  /// Initialize histograms
  for(UInt_t ix=0; ix<mNXaxis; ix++){
    // 1D, g + "histogram name" (will be used later to store Profile/Weight)
    // This is the output for final average quantities
    TH1* h1D = (TH1D*) GetTH1D(mName, ix) ;
    mHistogram1D.push_back( (TH1D*) h1D );

    // 2D
    TH2* h2D = (TH2D*) GetTH2D(mName, ix) ;
    mHistogram2D.push_back( (TH2D*) h2D );

    // Profile. Need to be corrected by weight histogram (hWeight)
    TProfile* hProfile = (TProfile*) GetTProfile(mName, ix, "Profile");
    mProfile.push_back( (TProfile*) hProfile );

    /// Initialize weight histogram
    TProfile* hWeight = (TProfile*) GetTProfile(mName, ix, "Weight");
    mWeight.push_back( (TProfile*) hWeight);
  }

  /// Print informations
  for(vector<TProfile*>::iterator iter = mProfile.begin();
      iter != mProfile.end(); iter++){
    TProfile* h = (TProfile*) (*iter);
  }
}

//____________________________________________________________________________________________________
void GlauberHistogramMaker::Reset()
{
  mHistogram1D.clear();
  mHistogram2D.clear();
  mProfile.clear();
  mWeight.clear();
  mXaxis.clear();
}

//____________________________________________________________________________________________________
void GlauberHistogramMaker::SetXaxis(const GlauberTree& tree, const CentralityMaker& centralityMaker, const TString type)
{
  /// Set x-axis (put variable in correct order, see Init())
  mXaxis[0] = tree.GetB()            ; // impact parameter
  mXaxis[1] = tree.GetNpart()        ; // Npart

  // Get multiplicity from tree
  mXaxis[2] = tree.GetMultiplicity() ; // Multiplicity

  /// Set id!=0 if type is "smallTotal" or "largeTotal"
  UInt_t mode = 0 ; // default
  if( type.CompareTo("smalltotal", TString::kIgnoreCase) == 0 )      mode = 1 ; // -5% total Xsection
  else if( type.CompareTo("largetotal", TString::kIgnoreCase) == 0 ) mode = 2 ; // +5% total Xsection

  const Centrality* centrality = centralityMaker.GetCentrality() ;

  mXaxis[3] = centrality->GetCentrality(mXaxis[2], mode) ; // Centrality
}

//____________________________________________________________________________________________________
Bool_t GlauberHistogramMaker::IsXaxisOk() const
{
  UInt_t isXaxisBad = 0 ;

  for(vector<Double_t>::const_iterator iter=mXaxis.begin(); iter != mXaxis.end(); iter++){
    if( (*iter) == -9999. ) isXaxisBad++ ;
  }

  return (isXaxisBad==0) ? kTRUE : kFALSE ;
}

//____________________________________________________________________________________________________
void GlauberHistogramMaker::Fill2D(vector<TH2*> collection, const Double_t y, const Double_t weight)
{
  /// Fill 2D histogram
  if( mXaxis.size() != collection.size() ){
    Error("StGlauberHistogramMaker::Fill2D", Form("Array size is different: (x-axis, 2D) = (%3d, %3d)",
          mXaxis.size(), collection.size()));
    assert(0);
  }

  for(UInt_t ix=0; ix<collection.size(); ix++){
    TH2* h = (TH2D*) collection[ix] ;
    const TString name(h->GetName());

    // Centrality bin is inclusive
    if( name.Contains("cent") ){
      const Double_t centrality = mXaxis[ix] ;
 
      // Loop over all centrality bins
      for(Int_t icent=0; icent<h->GetNbinsX(); icent++){
        // Fill histogram if current centrality is included 
        if(GlauberConstUtilities::IsCentralityOk(icent, centrality)){
          h->Fill(icent+0.5, y, weight);
        }
      }
    }
    else{
      // Others (except for centrality) fill as it is
      h->Fill(mXaxis[ix], y, weight);
    }
  }
}

//____________________________________________________________________________________________________
void GlauberHistogramMaker::FillProfile(vector<TProfile*> collection, const Double_t y)
{
  /// Fill profile histogram
  if( mXaxis.size() != collection.size() ){
    Error("GlauberHistogramMaker::FillProfile", Form("Array size is different: (x-axis, Weight) = (%3d, %3d)",
          mXaxis.size(), collection.size()));
    assert(0);
  }

  for(UInt_t ix=0; ix<collection.size(); ix++){
    TProfile* h = (TProfile*) collection[ix] ;
    const TString name(h->GetName());

    // Centrality bin is inclusive
    if( name.Contains("cent") ){
      const Double_t centrality = mXaxis[ix] ;
 
      // Loop over all centrality bins
      for(Int_t icent=0; icent<h->GetNbinsX(); icent++){
        // Fill histogram if current centrality is included 
        if(GlauberConstUtilities::IsCentralityOk(icent, centrality)){
          h->Fill(icent, y) ;
        }
      }
    }
    else{
      // Others (except for centrality) fill as it is
      h->Fill(mXaxis[ix], y);
    }
  }
}


//____________________________________________________________________________________________________
void GlauberHistogramMaker::Fill(const TString type, const Double_t y, const Double_t weight)
{
  /// Fill histogram
  /// Choose type '2d', 'profile' or 'weight'
  //
  /// (2d)->Fill(x, y, weight);
  /// (profile)->Fill(x, y*weight)
  /// (weight)->Fill(x, weight)
  ///

  if( type.CompareTo("2d", TString::kIgnoreCase) == 0 ){
    Fill2D(mHistogram2D, y, weight) ;
  }
  else if ( type.CompareTo("profile", TString::kIgnoreCase) == 0 ){
    FillProfile(mProfile, y*weight) ;
  }
  else if ( type.CompareTo("weight", TString::kIgnoreCase) == 0 ){
    FillProfile(mWeight, weight) ;
  }
  else{
    Error("GlauberHistogramMaker::Fill", "Unknown option, option=%s", type.Data());
    assert(0);
  }
    
}

//____________________________________________________________________________________________________
void GlauberHistogramMaker::Fill(const Double_t y, const Double_t weight)
{
  /// Fill 'y' value with 'weight'
  /// Re-weighting is only applied for numerator

  /// Check x-axis, abort if all bad
  if(!IsXaxisOk()){
    Error("StGlauberHistogramMaker::Fill", "x-axis has not defined yet for %s", mName.Data());
    assert(0);
  }

  /// Fill
  Fill("2d", y, 1.0) ; // No weight. Don't use this histogram for final result. Just for drawing.
  Fill("profile", y, weight);
  Fill("weight", y, weight) ;
}

//____________________________________________________________________________________________________
void GlauberHistogramMaker::DoWeightCorrection(vector<TH1*> collection1d, vector<TProfile*> collectionp)
{
  /// Profile / Weight
  /// Use TH1 instead of TProfile since TProfile is not accurate to calculate errors in Divide() function

  for(UInt_t ix=0; ix<collectionp.size(); ix++){
    TH1* h1 = collection1d[ix] ;

    /// Profile::Divide() computing error is not accurate, use TH1 instead
    TH1* p1 = (TH1D*) collectionp[ix]->ProjectionX();
    TH1* w1 = (TH1D*) mWeight[ix]->ProjectionX();
    p1->Divide(w1);

    p1->SetXTitle(h1->GetXaxis()->GetTitle());
    p1->SetYTitle(h1->GetXaxis()->GetTitle());

    /// Set y-axis range from 2D
    p1->SetMinimum(mYmin);
    p1->SetMaximum(mYmax);

    /// Copy contenst/errors
    for(Int_t i=0; i<p1->GetNbinsX(); i++){
      collection1d[ix]->SetBinContent(i+1, p1->GetBinContent(i+1));
      collection1d[ix]->SetBinError(i+1, p1->GetBinError(i+1));
    }
    collection1d[ix]->SetEntries(p1->GetEntries());
    collection1d[ix]->Print() ;
    delete p1 ;
    delete w1 ;
  }
}

//____________________________________________________________________________________________________
void GlauberHistogramMaker::WriteTable(vector<TH1*> collection, const TString name)
{
  /// Format will be:
  /// <centrality bin>  <minimum centrality>  <maximum centrality>  <value>  <stat. error>
  const TString tableName(Form("%s/table_%s_vs_centrality.txt", 
        mTableDirectory.Data(), name.Data()));
  ofstream fout(tableName.Data());

  const UInt_t centId = 3 ; // 4-th histogram has centrality in x-axis
  for(UInt_t ic=0; ic<StGlauberConstUtilities::GetCentralityBin(); ic++){
    const Double_t val     = collection[centId]->GetBinContent(ic+1);
    const Double_t error   = collection[centId]->GetBinError(ic+1);

    const Double_t centmin = StGlauberConstUtilities::GetCentralityMin(ic) ;
    const Double_t centmax = StGlauberConstUtilities::GetCentralityMax(ic) ;
    fout << Form("%10d    %1.1f      %1.1f        %1.5f        %1.5f",
        ic, centmin, centmax, val, error) << endl;
  }
  fout << endl << endl ;
  fout << "# <centrality bin>  <minimum centrality>  <maximum centrality>  <value>  <stat. error>" << endl;
  fout << Form("# %s", mTitle.Data()) << endl; // description (type, node name in tree)
}

//____________________________________________________________________________________________________
void GlauberHistogramMaker::WriteGraphs(vector<TH1*> collection)
{
  /// Write TGraphErrors into output ROOT file
  const UInt_t centId = 3; // 4-th histogram has centrality in x-axis

  TGraphErrors* graph = new TGraphErrors() ;

  TString name(collection[centId]->GetName()); // should be g + getHistogramName()
  name.Replace(0, 1, "graph");
  graph->SetName(name);
  graph->SetTitle(collection[centId]->GetTitle());

  // Fill up to 80 %
  const UInt_t ncentrality = 9 ;
  for(UInt_t ic=0; ic<ncentrality; ic++){
    const Double_t centmin = StGlauberConstUtilities::GetCentralityMin(ic) ;
    const Double_t centmax = StGlauberConstUtilities::GetCentralityMax(ic) ;
 
    const Double_t cent = ( centmin + centmax)/2.0 ;
    const Double_t val  = collection[centId]->GetBinContent(ic+1) ;
    const Double_t err  = collection[centId]->GetBinError(ic+1) ;
    graph->SetPoint(ic, cent, val);
    graph->SetPointError(ic, 0.0, err);
  }

  graph->Write() ;

}


//____________________________________________________________________________________________________
void GlauberHistogramMaker::Finish(const TString type)
{
  /// Weight correction
  DoWeightCorrection(mHistogram1D, mProfile) ;

  /// Write down the table
  const TString name(type + "_" + mName);
  WriteTable(mHistogram1D, name) ;

  /// Write graphs
  WriteGraphs(mHistogram1D) ;
}

//____________________________________________________________________________________________________
const TString GlauberHistogramMaker::GetHistogramName(const TString name, const UInt_t ix) const
{
  return Form("%s_%s", name.Data(), mXAxisName[ix].Data());
}

//____________________________________________________________________________________________________
TH1* GlauberHistogramMaker::GetTH1D(const TString name, const UInt_t ix)
{
  TH1* h = new TH1D("g" + GetHistogramName(name, ix), 
      Form("%s vs %s, %s", mYTitle.Data(), mXAxisTitle[ix].Data(), mTitle.Data()), mXAxisBin[ix], 0.0, mXAxisMax[ix]);
  h->SetXTitle(mXAxisTitle[ix].Data());
  h->SetYTitle(mYTitle.Data());

  return h ;
}

//____________________________________________________________________________________________________
TH2* GlauberHistogramMaker::GetTH2D(const TString name, const UInt_t ix)
{
  TH2* h = new TH2D("h" + GetHistogramName(name, ix), 
      Form("%s vs %s, %s (2D)", mYTitle.Data(), mXAxisTitle[ix].Data(), mTitle.Data()), 
      mXAxisBin[ix], 0.0, mXAxisMax[ix], mYbin, mYmin, mYmax);
  h->SetXTitle(mXAxisTitle[ix].Data());
  h->SetYTitle(mYTitle.Data());

  return h ;
}

//____________________________________________________________________________________________________
TProfile* GlauberHistogramMaker::GetTProfile(const TString name, const UInt_t ix, const TString title)
{
  const TString prefix = (title.CompareTo("profile", TString::kIgnoreCase)==0) ? "p" : "w" ;
  TProfile* h = new TProfile(prefix + GetHistogramName(name, ix), 
    Form("%s vs %s, %s (%s)", mYTitle.Data(), mXAxisTitle[ix].Data(), mTitle.Data(), title.Data()), mXAxisBin[ix], 0.0, mXAxisMax[ix]);
  h->SetXTitle(mXAxisTitle[ix].Data());
  h->SetYTitle(mYTitle.Data());

  return h ;
}

//____________________________________________________________________________________________________
const TString GlauberHistogramMaker::GetName() const
{ 
  return mName ;
}

