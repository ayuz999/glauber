#include <assert.h>
#include "TCanvas.h"
#include "TError.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"
#include "TH3.h"
#include "TLine.h"
#include "TStyle.h"
#include "NegativeBinomial.h"
#include "NbdFitMaker.h"
#include <iostream>

using namespace std;
NbdFitMaker::NbdFitMaker()
{
  mNBinomial = 0 ;
  mhRefMult = 0 ;
  mhRefMultSim = 0 ;
  mNData = 0 ;
  mMinimumMultiplicityCut = 100.0 ; // >50 by default. Can be changed by setMinimumMultiplicityCut(const Double_t)

  //mDoCentralityDetermination = kFALSE ;
  mDoCentralityDetermination = kTRUE;

  mCanvas = 0 ;
  mCutOff[0] = 0;
  mCutOff[1] = 0;
  mOneLine = 0 ;
  mChi2Graph = 0;
}

//____________________________________________________________________________________________________
// Default destructor
NbdFitMaker::~NbdFitMaker()
{
  if(mNBinomial) delete mNBinomial ;
  if(mCanvas) delete mCanvas ;

  for(Int_t i=0;i<2;i++){
    if(mCutOff[i]) delete mCutOff[i] ;
  }
  if(mOneLine) delete mOneLine ;

  if(mChi2Graph) delete mChi2Graph ;
}

//____________________________________________________________________________________________________
void NbdFitMaker::DoCentralityDetermination()
{
  mDoCentralityDetermination = kTRUE ;
}

//____________________________________________________________________________________________________
Double_t NbdFitMaker::GetNormalization(const TH1& h1, const TH1& h2, const Double_t min, const Double_t max) const
{
  // Get normalization factor in (min, max)
  Double_t numerator   = 0.0;
  Double_t denominator = 0.0;
  const Int_t mulminbin = h1.GetXaxis()->FindBin(min);
  const Int_t mulmaxbin = h1.GetXaxis()->FindBin(max);

  for(Int_t mul=mulminbin; mul<mulmaxbin; mul++){
    const Double_t n1      = h1.GetBinContent(mul+1);
    const Double_t n1Error = h1.GetBinError(mul+1);

    if( n1 == 0.0 || n1Error == 0.0 ) continue ;

    const Double_t n2 = h2.GetBinContent(mul+1);

    numerator   += n1 * n2 / (n1Error*n1Error) ;
    denominator += n2 * n2 / (n1Error*n1Error) ;
  }
  const Double_t norm = (denominator!=0.0) ? numerator / denominator : 1.0 ;

  return norm ;
}

//____________________________________________________________________________________________________
Double_t NbdFitMaker::CalculateChi2(const TH1& hdata, const TH1& hfunc, const Double_t minimumMultiplicityCut)
{
  /// Calculate chi2 from data and func
  mNData = 0 ;

  Double_t chi2 = 0.0 ;
  for(Int_t ix=0; ix<hdata.GetNbinsX(); ix++){
    // Lower multiplicity cut off
    const Double_t mult      = hdata.GetBinCenter(ix+1);
    if( mult < minimumMultiplicityCut ) continue ;

    // Check data points
    const Double_t dataError = hdata.GetBinError(ix+1);
    if( dataError == 0.0 ) continue ;

    // Calculate chi2
    const Double_t data  = hdata.GetBinContent(ix+1);
    const Double_t func  = hfunc.GetBinContent(ix+1);
    const Double_t delta = (data - func)/dataError ;
    chi2 += delta*delta ;
    mNData++;
  }

  return chi2 ;
}

//____________________________________________________________________________________________________
void NbdFitMaker::CalculateCentrality(const TH1& hdata, const TH1& hmc) const
{
  // - Calculate centrality from the MC multiplicity distribution
  // - Also calculate the re-weighting correction = MC/Data, important for peripheral (typically 60-80%)
  //
  // NOTE: Assume MC multiplicity has been normalized to the real data
  
  const UInt_t ncent = 16 ; // 0-80% (5% increment)
  Int_t centBin[2][3][ncent]; // Final centrality bins

  for(UInt_t i=0; i<2; i++) {
    // For cross check, do centrality determination
    // 1) from peripheral to central
    // 2) from central to peripheral

    // Centrality cut
    Double_t centralityCut[ncent] ;
    Double_t centralityMin[ncent] ;
    Double_t centralityMax[ncent] ;

    for(UInt_t ic=0; ic<ncent; ic++) {
      const Double_t sign         = (i==0) ? -1.0 : 1.0 ;
      const Double_t centStep     = 0.05 ;
      const Double_t centCutStart = (i==0) ? 0.80 : 0.05 ;
      const Double_t centMinStart = (i==0) ? 75.0 : 0.0 ;

      centralityCut[ic] = centCutStart + sign*centStep * ic ;
      centralityMin[ic] = centMinStart + sign*centStep * 100.0 * ic ;
      centralityMax[ic] = (centMinStart + centStep*100.0) + sign*centStep * 100.0 * ic ;

      // Print centrality definitions
    }

    // Start calculation for three different cases, for systematic error study
    // 0: default, 1:-5% total cross section, 2:+5% total cross section
    const Double_t nevent = hmc.Integral() ;
    for(Int_t it=0; it<3; it++){ // vary total cross section by 5 % in it=1,2
      Double_t scale = 1.0 ;
      if( it == 1 ) scale = 0.95 ;
      if( it == 2 ) scale = 1.05 ;

      UInt_t bin = 0 ;
      const Int_t nbin = hmc.GetNbinsX() ;
      Double_t sum = 0.0;
      for(Int_t im=0; im<nbin; im++){
        const Double_t M      = (i==0) ? hmc.GetBinCenter(im+1) : hmc.GetBinCenter(nbin-im) ;
        const Int_t Mint      = (i==0) ? im : nbin-im-1 ;
        const Double_t count  = (i==0) ? hmc.GetBinContent(im+1) : hmc.GetBinContent(nbin-im) ;
        if( count == 0.0 ) continue ;

        sum += count ;
        const Double_t fraction    = sum / nevent ;
        const Double_t fCentBinCut = centralityCut[bin] * scale;
        const Double_t R           = (i==0) ? (1.0 - fraction) : fraction ;
        const Bool_t isCentOk      = (i==0) ? R <= fCentBinCut : R > fCentBinCut ;

        if( isCentOk && bin < ncent ){
          cout << Form("%2.2f - %2.2f (%%) :  M > %4d (im=%3d, M=%1.1f, bin=%4d) (sum, total, fraction>cut) = (%1.3f, %1.3f, %1.3f>%1.3f)",
              TMath::Abs(centralityMin[bin]*scale), TMath::Abs(centralityMax[bin]*scale), Mint, im, M, bin, sum, nevent, R, fCentBinCut) << endl;

          centBin[i][it][bin] = (Double_t)Mint ;
          bin++;
        }
      }// multiplicity loop
    }// different total cross section
  }// from peripheral or central

}

//____________________________________________________________________________________________________
void NbdFitMaker::SetParameters(const Double_t npp, const Double_t k, const Double_t x,
    const Double_t efficiency, const Double_t triggerbias, const Bool_t isConstEfficiency)
{
  if(mNBinomial) delete mNBinomial ;
  mNBinomial = new NegativeBinomial(npp, k, x, efficiency, triggerbias, isConstEfficiency);
}

//____________________________________________________________________________________________________
void NbdFitMaker::SetMinimumMultiplicityCut(const Double_t cut)
{
  mMinimumMultiplicityCut = cut ;
}

//____________________________________________________________________________________________________
void NbdFitMaker::ReadData(const Char_t* data, const Char_t* glauber, const Char_t* dataHistogramName)
{
  // Read real data file
  TFile* inputData = TFile::Open(data);
  if(!inputData || !inputData->IsOpen()){
    Error("NbdFitMaker::readData", "can't open %s", data);
    assert(0);
  }
  
  mhRefMult = 0;
//  if( mNBinomial->useTpc() ){
    /// TPC refMult
    mhRefMult = (TH1D*) inputData->Get(dataHistogramName);
//    mhRefMult = (TH1D*) inputData->Get("hRefMultTpc");
//    mhRefMult = (TH1D*) inputData->Get("hRefMult");
//  }
//  else{
//    /// FTPC refMult
//    mhRefMult = (TH1D*) inputData->Get("hRefMultFTpc");
//  }

  if(!mhRefMult){
    Error("NbdFitMaker::readData", "hRefMult doesn't exist");
    assert(mhRefMult);
  }

  mhRefMult->SetLineColor(1);

  // Define simulated refmult
  const Int_t nbinsx  = mhRefMult->GetNbinsX() ;
  const Double_t xmin = mhRefMult->GetXaxis()->GetXmin() ;
  const Double_t xmax = mhRefMult->GetXaxis()->GetXmax() ;
  mhRefMultSim = new TH1D("hRefMultSim", "", nbinsx, xmin, xmax);
  mhRefMultSim->SetLineColor(2);

  // Sumw2 to calculate error properly
  mhRefMult->Sumw2();
  mhRefMultSim->Sumw2();

  // Read glauber file
  TFile* inputGlauber = TFile::Open(glauber);
  if(!inputGlauber || !inputGlauber->IsOpen()){
    Error("NbdFitMaker::readData", "can't open %s", glauber);
    assert(0);
  }
  
  mhNcoll_Npart = (TH2D*) inputGlauber->Get("hNcoll_Npart");
  if(!mhNcoll_Npart){
    Error("NbdFitMaker::readData", "hNcoll_Npart doesn't exist");
    assert(mhNcoll_Npart);
  }
}

//____________________________________________________________________________________________________
TGraph* NbdFitMaker::Fit(const Int_t nevents, const Char_t* outputFileName)
{
  gStyle->SetOptStat(0);

  /// Fit real data by simulated multiplicity distribution
  /// Make sure the refmult and Ncoll_Npart histograms have benn opened
  if(!mhRefMult){
    Error("NbdFitMaker::Fit", "hRefMult doesn't exist");
    assert(mhRefMult);
  }
  if(!mhNcoll_Npart){
    Error("NbdFitMaker::Fit", "hNcoll_Npart doesn't exist");
    assert(mhNcoll_Npart);
  }
  mhRefMultSim->Reset(); // will clear histogram
  Int_t ievent = 0 ;
  while( ievent < nevents ) {
    Double_t npart, ncoll;
    mhNcoll_Npart->GetRandom2(npart, ncoll);
    //const Bool_t isNpartNcollOk = (npart>=2 && ncoll>=1) ;
    //if ( !isNpartNcollOk ) continue ;
    if(!(npart>=2 && ncoll>=1)) continue;

    const Int_t multiplicity = mNBinomial->GetMultiplicity(npart, static_cast<Int_t>(ncoll));
    mhRefMultSim->Fill(multiplicity);


    ievent++;
  }

  // Normalization
  const Double_t norm = GetNormalization(*mhRefMult, *mhRefMultSim, mMinimumMultiplicityCut, mhRefMult->GetXaxis()->GetXmax());
  mhRefMultSim->Scale(norm);

  // Get chi2
  const Double_t chi2 = CalculateChi2(*mhRefMult, *mhRefMultSim, mMinimumMultiplicityCut);

  //----------------------------------------------------------------------------------------------------
  // Set chi2
  //----------------------------------------------------------------------------------------------------
  if(mChi2Graph) delete mChi2Graph ;
  mChi2Graph = new TGraph();
  mChi2Graph->SetPoint(0, 0, chi2);
  mChi2Graph->SetPoint(1, 1, mNData - 2); // 2 fitting parameters (npp, k)
  mChi2Graph->SetPoint(2, 2, chi2/(mNData-2.0));

  //----------------------------------------------------------------------------------------------------
  // Draw
  //----------------------------------------------------------------------------------------------------
  mhRefMult->SetMinimum(0.1);
  mhRefMult->SetMaximum(mhRefMult->GetMaximum()*10.0);

  mhRefMultSim->SetXTitle("Refmult (MC)");
  mhRefMultSim->SetYTitle("Count");
  mhRefMultSim->SetTitle(Form("npp=%1.2f, k=%1.2f, x=%1.2f",
  mNBinomial->GetNpp(), mNBinomial->GetK(), mNBinomial->GetX()));

  if(mCanvas) delete mCanvas ;
  mCanvas = new TCanvas("c1", "c1", 1200, 500);
  mCanvas->Divide(2, 1);
  mCanvas->cd(1);
  mCanvas->GetPad(1)->SetLogy(1);

  mhRefMult->Draw("h");
  mhRefMultSim->Draw("hsame");

  // Normalization line
  if(mCutOff[0]) delete mCutOff[0] ;
  mCutOff[0] = new TLine( mMinimumMultiplicityCut, mhRefMult->GetMinimum(), mMinimumMultiplicityCut, mhRefMult->GetMaximum());
  mCutOff[0]->SetLineColor(4);
  mCutOff[0]->SetLineStyle(2);
  mCutOff[0]->Draw();

  // Draw ratio of MC to data
  mCanvas->cd(2);
  TH1* hRatio = (TH1D*) mhRefMultSim->Clone();
  hRatio->SetName("hRatio");
  hRatio->Divide(mhRefMult);
  hRatio->SetYTitle("MC/data");

  hRatio->SetMinimum(0);
  hRatio->SetMaximum(2);
  hRatio->Draw();

  if(mOneLine) delete mOneLine ;
  mOneLine = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0, hRatio->GetXaxis()->GetXmax(), 1.0);
  mOneLine->SetLineColor(4);
  mOneLine->SetLineStyle(2);
  mOneLine->Draw();

  if(mCutOff[1]) delete mCutOff[1] ;
  mCutOff[1] = new TLine( mMinimumMultiplicityCut, hRatio->GetMinimum(), mMinimumMultiplicityCut, hRatio->GetMaximum());
  mCutOff[1]->SetLineColor(4);
  mCutOff[1]->SetLineStyle(2);
  mCutOff[1]->Draw();

  mCanvas->cd();
  mCanvas->Update();

  //----------------------------------------------------------------------------------------------------
  // Centrality
  //----------------------------------------------------------------------------------------------------
  if ( mDoCentralityDetermination ) {
    CalculateCentrality(*mhRefMult, *mhRefMultSim) ;
  }

  //----------------------------------------------------------------------------------------------------
  // Write only if outputFileName is given
  //----------------------------------------------------------------------------------------------------
  const TString fileName(outputFileName);
  if(!fileName.IsWhitespace()){
    TFile* output = TFile::Open(outputFileName, "recreate");
    mhRefMult->Write();
    mhRefMultSim->Write();
    hRatio->Write();
    output->Close();
  }


  return mChi2Graph;
}

//____________________________________________________________________________________________________
Int_t NbdFitMaker::Scan(const Int_t nevents,
    const Int_t nppbin, const Double_t nppmin, const Double_t nppmax,
    const Int_t kbin, const Double_t kmin, const Double_t kmax,
    const Int_t xbin, const Double_t xmin, const Double_t xmax,
    //const Double_t x,
    //const Int_t effbin, const Double_t effmin, const Double_t effmax,
    const Double_t efficiency,
    const Double_t triggerbias, const Bool_t isConstEfficiency
){
  /// Loop over all (npp, k, x) to find out minimum chi2
  TH3* hChi2 = new TH3D("hChi2", "#chi^{2}/NDF in (npp, k, x) space",
      nppbin, nppmin, nppmax, kbin, kmin, kmax, xbin, xmin, xmax);
  hChi2->SetXTitle("n_{pp}");
  hChi2->SetYTitle("k");
  hChi2->SetZTitle("x");

  const Double_t nppstep = (nppmax-nppmin)/static_cast<Double_t>(nppbin) ;
  const Double_t kstep   = (kmax-kmin)/static_cast<Double_t>(kbin) ;
  const Double_t xstep   = (xmax-xmin)/static_cast<Double_t>(xbin) ;

  for(Int_t ix=0; ix<nppbin; ix++){
    const Double_t npp = nppmin + nppstep*ix ;

    for(Int_t iy=0; iy<kbin; iy++){
      const Double_t k = kmin + kstep*iy ;
      for(Int_t iz=0; iz<xbin; iz++){
		const Double_t x = xmin + xstep*iz ;
        // Set parameters
        SetParameters(npp, k, x, efficiency, triggerbias, isConstEfficiency);

//add by Lizhu
  
    const Char_t* fileforratio(Form("Ratio_npp%1.3f_k%1.3f_x%1.3f_eff%1.3f.root", npp, k, x, efficiency));
    
        // Fitting,  write ROOT file
    Fit(nevents, fileforratio);
        // Get chi2
    hChi2->SetBinContent(ix+1, iy+1, iz+1, mChi2Graph->GetY()[2]);
      }// x loop
    }// k loop
  }// npp loop

  const Char_t* fileName(Form("chi2_nevents%d_npp%1.3f-%1.3f_k%1.3f-%1.3f_x%1.3f_%1.3f_eff%1.3f.root",
        nevents, nppmin, nppmax, kmin, kmax, xmin, xmax, efficiency));
  TFile* outputFile = TFile::Open(fileName, "recreate");
  hChi2->Write();
  outputFile->Close();
  return 1;
}
