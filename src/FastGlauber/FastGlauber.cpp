#include "TF1.h"
#include "TF3.h"
#include "TH1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TVector3.h"
#include "Nucleon.h"
#include "../CentralityMaker/CentralityMaker.h"
#include "../CentralityMaker/NegativeBinomial.h"
#include "../GlauberTree/GlauberTree.h"
#include "../GlauberUtilities/GlauberUtilities.h"
#include "FastGlauber.h"
#include <string>
#include <iostream>
#include <assert.h>
#include "Rtypes.h"
using namespace std;

FastGlauber::FastGlauber(
        const string outputFileName,
        const int    energy,
        const string type
    ) : mEnergy(energy), mOutputFileName(outputFileName),
  mInelasticNNCrossSection(GetInelasticNNCrossSection(mEnergy, type))
{
  InitAuAu("AuAu");
}

void FastGlauber::SetRepulsionDistance(const Double_t repulsionDistance):mRepulsionDistance(repulsionDistance){}

Int_t FastGlauber::Clear()
{
  mGlauberTree->Clear() ;

  for(UInt_t in=0;in<2;in++){
    for(UInt_t i=0;i<mNucleons[in].size();i++){
      mNucleons[in][i]->Reset() ;
    }
  }

  return 1 ;
}

//____________________________________________________________________________________________________
Int_t FastGlauber::InitTree()
{
  // Initialize Glauber tree in write mode
  mGlauberTree = new GlauberTree(1);
  // Clear data members in tree
  Clear() ;
  // Initialize tree
  mGlauberTree->Open(mOutputFileName) ;
  // Histograms
  for(Int_t in=0;in<4;in++){
    TString nucleus("Projectile");
    if( in % 2 == 1 ) nucleus = "Target";

    TString title(Form("Woods-saxon density profile (%s)", nucleus.Data()));
    if(in>=2) title = Form("Woods-saxon density profile after repulsion (%s)", nucleus.Data());

    mhWoodsSaxon[in] = new TH1D(Form("hWoodsSaxon_%d", in), title, 400, 0, 20);
    mhWoodsSaxon[in]->SetXTitle("r (fm)");
  }

  // Sort output histograms/tree
  mGlauberTree->Sort();

  return 1;
}

//____________________________________________________________________________________________________
Int_t FastGlauber::Init(
        const UInt_t massNumberA,  // Mass number of nucleus for nucleus A
        const Double_t radiusA,    // Radius of nucleus for nucleus A
        const Double_t skinDepthA, // Skin depth of nucleus for nucleus A
        const UInt_t massNumberB,  // Mass number of nucleus for nucleus B
        const Double_t radiusB,    // Radius of nucleus for nucleus B
        const Double_t skinDepthB, // Skin depth of nucleus for nucleus B
        const TString type
){

  mNeventsThrow  = 0 ; // Number of all events
  mNeventsAccept = 0 ; // Number of accepted events (Ncoll>0)

  // Centrlaity maker
  const Char_t* system(Form("%s%s_%dGeV", GetName(massNumberA), GetName(massNumberB), 
        static_cast<Int_t>(mEnergy))) ;
  mCentralityMaker = new CentralityMaker(system);

  /// Repulsion distance is 0fm by default. Use SetRepulsionDistance(const Double_t) method to set finite value
  mRepulsionDistance = 0.0 ;

  mHardCoreSmearing = kFALSE ; /// Hard-core smearing, default is OFF
  mGaussianSmearing = kFALSE ; /// Gaussian smearing, default is OFF

  mCollisionProfile = mkHardCoreProfile  ; /// Default is hard-core collision
  if ( type.CompareTo("gauss", TString::kIgnoreCase) == 0 ){
    //DoGaussianCollision() ;
  }

  /// Hard-core smearing (only use if mHardCoreSmearing is true)
  const Double_t dmaxh = TMath::Sqrt(mInelasticNNCrossSection/TMath::Pi());
  mfHardCore = new TF3("fHardCore", GlauberUtilitiesSpace::StepFunction, -dmaxh, dmaxh, -dmaxh, dmaxh, -dmaxh, dmaxh, 1);
  mfHardCore->SetParameter(0, mInelasticNNCrossSection);

  /// Gaussian smearing (only use if mGaussianSmearing is true)
  const Double_t sigma = 0.79/TMath::Sqrt(3.0) ;
  const Double_t dmaxg = sigma*5.0 ;
  mfGaussian = new TF3("fGaussian", GlauberUtilitiesSpace::Gaussian, -dmaxg, dmaxg, -dmaxg, dmaxg, -dmaxg, dmaxg, 1);
  mfGaussian->SetParameter(0, sigma); // width
      /// Initialize Woods-saxon density profile for deformed nuclei
      mfWoodsSaxon[in] = 0 ; // NULL pointer for spherical woods-saxon

      mfWoodsSaxon2D[in] = new TF2(Form("fWoodsSaxon2D_%d", in), GlauberUtilitiesSpace::WoodsSaxon2D, 0, 20, -1.0, 1.0, 4);
      mfWoodsSaxon2D[in]->SetParName(0, "Radius");
      mfWoodsSaxon2D[in]->SetParName(1, "Skin depth");
      mfWoodsSaxon2D[in]->SetParName(2, "#beta_{2}");
      mfWoodsSaxon2D[in]->SetParName(3, "#beta_{4}");
//      mfWoodsSaxon2D[in]->SetParName(4, "#beta_{6}");

      // Default number of sampling points are too small (probably for y-axis) in TF2
      // this causes the step-like structures if you get (r,theta) from TF2::GetRandom2()
      // --> Increase Npx, Npy by a factor of 2
      mfWoodsSaxon2D[in]->SetNpx(200);
      mfWoodsSaxon2D[in]->SetNpy(200);

      if( in == 0 ){
        mfWoodsSaxon2D[in]->SetParameter(0, radiusA) ;
        mfWoodsSaxon2D[in]->SetParameter(1, skinDepthA) ;
        mfWoodsSaxon2D[in]->SetParameter(2, beta2A) ;
        mfWoodsSaxon2D[in]->SetParameter(3, beta4A) ;
      }
      else{
        mfWoodsSaxon2D[in]->SetParameter(0, radiusB) ;
        mfWoodsSaxon2D[in]->SetParameter(1, skinDepthB) ;
        mfWoodsSaxon2D[in]->SetParameter(2, beta2B) ;
        mfWoodsSaxon2D[in]->SetParameter(3, beta4B) ;
      }

    }
    else{
      /// Initialize Woods-saxon density profile for spherical nuclei
      // Radius and skin depth
      mfWoodsSaxon2D[in] = 0 ; // NULL pointer for deformed woods-saxon

      mfWoodsSaxon[in] = new TF1(Form("fWoodsSaxon_%d", in), GlauberUtilitiesSpace::WoodsSaxon, 0, 20, 2);
      mfWoodsSaxon[in]->SetParName(0, "Radius");
      mfWoodsSaxon[in]->SetParName(1, "Skin depth");
      if( in == 0 ){
        mfWoodsSaxon[in]->SetParameter(0, radiusA) ;
        mfWoodsSaxon[in]->SetParameter(1, skinDepthA) ;
      }
      else{
        mfWoodsSaxon[in]->SetParameter(0, radiusB) ;
        mfWoodsSaxon[in]->SetParameter(1, skinDepthB) ;
      }
    }
    /// Initialize nucleons
    mNucleons[in].clear() ;
    const UInt_t A = (in==0) ? massNumberA : massNumberB ;

    for(UInt_t i=0;i<A;i++){
      mNucleons[in].push_back( new Nucleon() ) ;
    }
  }

  // Initialize output tree
  InitTree() ;

  return 1 ;
}

//____________________________________________________________________________________________________
Int_t FastGlauber::InitAuAu(const TString type)
{
  /// Initialize Au + Au collisions
  const UInt_t A   = 197 ;
  const Double_t R = 6.38 ;
  const Double_t d = 0.535 ;

  // Parameterization with hard-core smearing
  // from PRC79, 064904 (2009) T. Hirano et. al.
//  const Double_t R = 6.42 ;
//  const Double_t d = 0.544 ;
  Double_t RError = 0.0 ; // Absolute error on R
  Double_t dError = 0.0 ; // Absolute error on d


  // Error on R and d is taken from Atomic Data and Nuclear Table 36, 495 (1987)
  // R = 6.38 +/- 0.06, d = 0.535 +/- 0.027
  //   --> take 2 sigma
  if ( type.CompareTo("large", TString::kIgnoreCase) == 0 ){
    RError = 0.12 ;
    dError = -0.054 ;
   // RError = 0.32 ;   //5%
   // dError = -0.0267 ;//5%
  }
  else if ( type.CompareTo("small", TString::kIgnoreCase) == 0 ){
    RError = -0.12 ;
    dError = 0.054 ;
   // RError = -0.32 ;  //5%
   // dError = 0.0267 ; //5%
  }
  else if ( type.CompareTo("gauss", TString::kIgnoreCase) == 0 ){
    DoGaussianCollision() ;
  }
  else if ( type.CompareTo("smallNpp", TString::kIgnoreCase) == 0 ){
    mMode = 1 ;
  }
  else if ( type.CompareTo("largeNpp", TString::kIgnoreCase) == 0 ){
    mMode = 2 ;
  }

  return Init(A, R+RError, d+dError, A, R+RError, d+dError, type);
}

//____________________________________________________________________________________________________
Bool_t FastGlauber::IsCollision(Nucleon* nucleon0, Nucleon* nucleon1) const
{
  const Double_t dx     = nucleon0->GetX() - nucleon1->GetX() ;
  const Double_t dy     = nucleon0->GetY() - nucleon1->GetY() ;
  const Double_t dt     = TMath::Sqrt(dx*dx + dy*dy) ;
  const Double_t cutoff = TMath::Sqrt(mInelasticNNCrossSection/TMath::Pi());

  switch ( mCollisionProfile ){
    /// Hard-core collision profile (default)
    case mkHardCoreProfile:
      {
        return dt <= cutoff ;
      }

    /// Gaussian collision profile
    case mkGaussianProfile:
      {
        const Double_t sigma = cutoff ;
        const Double_t arg   = 0.5*dt/sigma ;
        return ( GlauberUtilities::instance()->GetUniform() <= TMath::Exp(-arg*arg) );
      }

    default:
      {
        assert(0);
      }
  }

  return kFALSE ;
}

//____________________________________________________________________________________________________
Int_t FastGlauber::Make()
{
  /// Clear all data members
  Clear() ;

  GlauberUtilities* glauberUtilities = GlauberUtilities::instance() ;

  /// Impact parameter
  const Double_t impactParameter = glauberUtilities->GetImpactParameter() ;
  mGlauberTree->SetB( impactParameter );

  //----------------------------------------------------------------------------------------------------
  /// 1. Generate nucleon positions (r,theta,phi)
  /// 2. Smearing nucleon positions if switch is ON (either DoHardCoreSmearing() or DoGaussianSmearing())
  /// 3. Distribute (r,theta,phi) with repulsion distance ds if ds != 0
  //----------------------------------------------------------------------------------------------------

  for(UInt_t in=0;in<2;in++){
    const Double_t Theta =  0.0 ;
    const Double_t Phi   =  0.0 ;
    mGlauberTree->SetTheta(in, Theta);
    mGlauberTree->SetPhi(in, Phi);

    /// Determine impact parameter for each nucleus
    const Double_t b = (in==0) ? -impactParameter/2.0 : impactParameter/2.0 ;

    UInt_t nNucleons = 0 ;
    while( nNucleons < mNucleons[in].size() ){
      Double_t r     = 0.0 ;
      Double_t theta = 0.0 ;
      Double_t phi   = 0.0 ;

      // Get (r,theta,phi) for nucleons
      GetRThetaPhi(in, r, theta, phi) ;

      mhWoodsSaxon[in]->Fill(r, 1.0/(r*r)) ;

      if( mRepulsionDistance == 0.0 ){
        /// No repulsion, just store (r,theta,phi)
        mNucleons[in][nNucleons]->Set(r, theta, phi, b, Theta, Phi, kTRUE) ;
        nNucleons++;
      }
      else{
        if(nNucleons==0){
          /// Store first nucleon, and increment total
          mNucleons[in][nNucleons]->Set(r, theta, phi, b, Theta, Phi, kTRUE) ;
          nNucleons++;
        }
        else{
          /// Check all the nucleons stored whether there are overlap each other
          /// defined by mRepulsionDistance. If any nucleons overlap, then discard
          /// current nucleon and try again
          Bool_t isOk = kFALSE ;
          do {
            const Double_t x = r*TMath::Sin(theta)*TMath::Cos(phi) ;
            const Double_t y = r*TMath::Sin(theta)*TMath::Sin(phi) ;
            const Double_t z = r*TMath::Cos(theta);

            Bool_t isOverlap = kFALSE ;
            for(UInt_t i=0; i<nNucleons; i++){
              const Double_t x0 = mNucleons[in][i]->GetX() ;
              const Double_t y0 = mNucleons[in][i]->GetY() ;
              const Double_t z0 = mNucleons[in][i]->GetZ() ;
              const Double_t dx = x - x0 ;
              const Double_t dy = y - y0 ;
              const Double_t dz = z - z0 ;

              // If any nucleons are found to be overlapped to the current one, stop the loop to save time
              if( TMath::Sqrt(dx*dx + dy*dy + dz*dz) <= mRepulsionDistance ){
                isOverlap = kTRUE ;
                break ;
              }
            }// Looping over other nucleons

            if(isOverlap){
              // There are overlap nucleons, try again
              isOk = kFALSE ;

              // Get (r,theta,phi) for nucleons
              GetRThetaPhi(in, r, theta, phi) ;
            }
            else{
              // if no nucleons are overlapped with the current one
              isOk = kTRUE ;
            }

          } while( !isOk ) ; // Find any overlaped nucleons ?

          // Add positions and increment total
          mNucleons[in][nNucleons]->Set(r, theta, phi, b, Theta, Phi, kTRUE) ;
          nNucleons++;
        }// Check overlap from 2nd nucleons
      }// mRepulsionDistance > 0
    }// Looping over all nucleons
  }// Nucleus loop

  //----------------------------------------------------------------------------------------------------
  /// Determine Number of collisions by looping over all pair of nucleons
  //----------------------------------------------------------------------------------------------------
  UInt_t Ncoll = 0 ;
  for(UInt_t in=0; in<mNucleons[0].size(); in++){
    Nucleon* nucleon0 = mNucleons[0][in] ;
    for(UInt_t jn=0; jn<mNucleons[1].size(); jn++){
      Nucleon* nucleon1 = mNucleons[1][jn] ;

      if(IsCollision(nucleon0, nucleon1)){
        nucleon0->IncrementNcoll() ;
        nucleon1->IncrementNcoll() ;

        Ncoll++;
      }
    }
  }

  //----------------------------------------------------------------------------------------------------
  /// Need at least one Ncoll
  //----------------------------------------------------------------------------------------------------
  if( Ncoll == 0 ) return 0 ;

  //----------------------------------------------------------------------------------------------------
  /// Determine Npart
  /// Calculate <x>, <y>, <x^2>, <y^2>, <xy>, ecc_{RP} and ecc_{PP}
  //----------------------------------------------------------------------------------------------------
  Double_t nSum[4] ;
  Double_t sumx[4] ;
  Double_t sumy[4] ;
  Double_t sumx2[4] ;
  Double_t sumy2[4] ;
  Double_t sumxy[4] ;

  for(Int_t i=0;i<4;i++){
    nSum[i] = 0.0 ;
    sumx[i] = 0.0 ;
    sumy[i] = 0.0 ;
    sumx2[i] = 0.0 ;
    sumy2[i] = 0.0 ;
    sumxy[i] = 0.0 ;
  }

  UInt_t Npart = 0 ;
  for(UInt_t in=0;in<2;in++){
    for(UInt_t i=0;i<mNucleons[in].size();i++){
      Nucleon* nucleon = mNucleons[in][i] ;
      const UInt_t ncoll = nucleon->GetNcoll() ;
      if( ncoll > 0 ){
        // Participant weight
        sumx[0]   += nucleon->GetXYZ("x");
        sumy[0]   += nucleon->GetXYZ("y");
        sumx2[0]  += nucleon->GetXYZ("xx");
        sumy2[0]  += nucleon->GetXYZ("yy");
        sumxy[0]  += nucleon->GetXYZ("xy");
        nSum[0]++ ; // Should be identical to Npart

        // Ncoll weight
        sumx[1]   += ncoll * nucleon->GetXYZ("x");
        sumy[1]   += ncoll * nucleon->GetXYZ("y");
        sumx2[1]  += ncoll * nucleon->GetXYZ("xx");
        sumy2[1]  += ncoll * nucleon->GetXYZ("yy");
        sumxy[1]  += ncoll * nucleon->GetXYZ("xy");
        nSum[1]   += ncoll ;

        // Multiplicity weight
        const Double_t mult = mCentralityMaker->GetNegativeBinomial()->GetTwoComponentMultiplicity(1.0, ncoll) ;
        nucleon->SetMultiplicity(mult);
        sumx[2]   += mult * nucleon->GetXYZ("x");
        sumy[2]   += mult * nucleon->GetXYZ("y");
        sumx2[2]  += mult * nucleon->GetXYZ("xx");
        sumy2[2]  += mult * nucleon->GetXYZ("yy");
        sumxy[2]  += mult * nucleon->GetXYZ("xy");
        nSum[2]   += mult ;

        nucleon->IncrementNpart() ;
        Npart++;
      }
      else{
        // Spectator
        sumx[3]   += nucleon->GetXYZ("x");
        sumy[3]   += nucleon->GetXYZ("y");
        sumx2[3]  += nucleon->GetXYZ("xx");
        sumy2[3]  += nucleon->GetXYZ("yy");
        sumxy[3]  += nucleon->GetXYZ("xy");

        nSum[3]++ ; // Should be identical to MassNumber - Npart
      }
    }// Nucleon loop
  }// Nucleus loop

  // Set multiplicity
  mGlauberTree->SetNpart(Npart) ;
  mGlauberTree->SetNcoll(Ncoll) ;

  // Get NegativeBinomial with different (npp, x)
  const NegativeBinomial* nbinomial = mCentralityMaker->GetNegativeBinomial() ;
  mGlauberTree->SetMultiplicity( nbinomial->GetMultiplicity(Npart, Ncoll) );

  // Get average and calculate eccentricity
  for(UInt_t i=0;i<4;i++){
    if(nSum[i]==0.0) continue ;

    sumx[i]  /= nSum[i] ;
    sumy[i]  /= nSum[i] ;
    sumx2[i] /= nSum[i] ;
    sumy2[i] /= nSum[i] ;
    sumxy[i] /= nSum[i] ;

    mGlauberTree->SetSumX(i, sumx[i]);
    mGlauberTree->SetSumY(i, sumy[i]);
    mGlauberTree->SetSumX2(i, sumx2[i]);
    mGlauberTree->SetSumY2(i, sumy2[i]);
    mGlauberTree->SetSumXY(i, sumxy[i]);

    // Eccentricity
    const Double_t sigmaX2  = sumx2[i] - sumx[i]*sumx[i] ;
    const Double_t sigmaY2  = sumy2[i] - sumy[i]*sumy[i] ;
    const Double_t sumSigma = sigmaX2 + sigmaY2 ;
    Double_t eccRP2 = -9999. ;
    Double_t eccPP2 = -9999. ;
    Double_t eccPP3 = -9999. ;
    Double_t eccPP4 = -9999. ;
    Double_t PP2 = -9999. ;
    Double_t PP3 = -9999. ;
    Double_t PP4 = -9999. ;

    if( sumSigma == 0.0 ){
      // Check denominator
      eccRP2 = -9999. ;
      eccPP2 = -9999. ;
      eccPP3 = -9999. ;
      eccPP4 = -9999. ;
      PP2 = -9999. ;
      PP3 = -9999. ;
      PP4 = -9999. ;
    }
    else{
      //----------------------------------------------------------------------------------------------------
      // Reaction plane eccentricity
      //----------------------------------------------------------------------------------------------------
      eccRP2 = (sigmaY2-sigmaX2)/sumSigma ;

      //----------------------------------------------------------------------------------------------------
      // The n-th order participant plane eccentricity
      //----------------------------------------------------------------------------------------------------
      for(Int_t io=0; io<3; io++){
        const Double_t order = io + 2.0 ;
        TGraph* qpart = GetParticipantEccentricity(order, sumx[i], sumy[i], nSum[i], i);
        if( io == 0 ) { PP2 = qpart->GetY()[2] ;  eccPP2 = qpart->GetY()[3] ; }
        if( io == 1 ) { PP3 = qpart->GetY()[2] ;  eccPP3 = qpart->GetY()[3] ; }
        if( io == 2 ) { PP4 = qpart->GetY()[2] ;  eccPP4 = qpart->GetY()[3] ; }
        delete qpart ;
      }
    }

    mGlauberTree->SetEccRP2(i, eccRP2);
    mGlauberTree->SetEccPP2(i, eccPP2);
    mGlauberTree->SetEccPP3(i, eccPP3);
    mGlauberTree->SetEccPP4(i, eccPP4);
    mGlauberTree->SetPP2(i, PP2);
    mGlauberTree->SetPP3(i, PP3);
    mGlauberTree->SetPP4(i, PP4);
  }

  return 1;
}

//____________________________________________________________________________________________________
Int_t FastGlauber::Run(const UInt_t nevents)
{

  while( mNeventsAccept < nevents ){
    if ( Make() == 1 ){
      // Fill tree
      mGlauberTree->Fill() ;
      mNeventsAccept++;
    }
    mNeventsThrow++;
  }
  return 1;
}

//____________________________________________________________________________________________________
Int_t FastGlauber::Finish()
{
  // Store header info.

  // Name of nucleus
  const UInt_t A[2] = {mNucleons[0].size(), mNucleons[0].size()};
  mGlauberTree->SetNameNucleusA(GetName(A[0]));
  mGlauberTree->SetNameNucleusB(GetName(A[1]));

  // Mass numbers
  mGlauberTree->SetMassNumberA(A[0]);
  mGlauberTree->SetMassNumberB(A[1]);

  // Radius
  mGlauberTree->SetRadiusA(mfWoodsSaxon[0]->GetParameter(0));
  mGlauberTree->SetRadiusB(mfWoodsSaxon[1]->GetParameter(0));

  // Skin depth
  mGlauberTree->SetSkinDepthA(mfWoodsSaxon[0]->GetParameter(1));
  mGlauberTree->SetSkinDepthB(mfWoodsSaxon[1]->GetParameter(1));

  // sigmaNN, sqrt(sNN)
  mGlauberTree->SetSigmaNN( mInelasticNNCrossSection * 10.0 ) ; // mb
  mGlauberTree->SetSqrtSNN( mEnergy ) ;

  // Repulsion distance
  mGlauberTree->SetRepulsionD(mRepulsionDistance);

  // Total cross section (mb) skip if no events have been processed
  Double_t totalX      = 0.0 ;
  Double_t totalXError = 0.0 ;
  const Double_t bmax  = GlauberUtilities::instance()->GetMaximumImpactParameter() ;

  if( mNeventsAccept != 0 ){
    const Double_t scaleFactor = TMath::Pi() * bmax * bmax * 10.0 ;
    totalX      = (Double_t)mNeventsAccept / (Double_t)mNeventsThrow * scaleFactor ;
    totalXError = totalX * TMath::Sqrt(1.0/mNeventsAccept + 1.0/mNeventsThrow) ;
  }
  mGlauberTree->SetTotalXsec(totalX) ;
  mGlauberTree->SetTotalXsecError(totalXError) ;

  // Smearing options (0 or 1)
//  UInt_t smearHardCore = (mHardCoreSmearing) ? 1 : 0 ;
//  UInt_t smearGaussian = (mGaussianSmearing) ? 1 : 0 ;
  mGlauberTree->SetSmearHardCore(mHardCoreSmearing);
  mGlauberTree->SetSmearGaussian(mGaussianSmearing);

  // Collision profile options (0 or 1)
  UInt_t collisionHardCore = (mCollisionProfile==mkHardCoreProfile) ? 1 : 0 ;
  UInt_t collisionGaussian = (mCollisionProfile==mkGaussianProfile) ? 1 : 0 ;
  mGlauberTree->SetCollisionHardCore(collisionHardCore);
  mGlauberTree->SetCollisionGaussian(collisionGaussian);

  // Maximum impact parameter
  mGlauberTree->SetBMax(bmax);

  // Number of events
  mGlauberTree->SetNeventsAccept(mNeventsAccept);
  mGlauberTree->SetNeventsThrow(mNeventsThrow);

  // NBD parameters
  const NegativeBinomial* nb = mCentralityMaker->GetNegativeBinomial() ;
  mGlauberTree->SetNpp( nb->GetNpp() );
  mGlauberTree->SetK( nb->GetK() );
  mGlauberTree->SetX( nb->GetX() );
  mGlauberTree->SetEfficiency( nb->GetEfficiency() );
  mGlauberTree->SetIsConstEfficiency( nb->IsConstEfficiency() );

  // version
  mGlauberTree->FillHeader() ;

  mGlauberTree->Close();

  return 1;
}

//____________________________________________________________________________________________________
void FastGlauber::DoHardCoreSmearing()
{
  mHardCoreSmearing = kTRUE ;
  mGaussianSmearing = kFALSE ;
}

//____________________________________________________________________________________________________
void FastGlauber::DoGaussianSmearing()
{
  mGaussianSmearing = kTRUE ;
  /// Turn off Hard-core smearing
  mHardCoreSmearing = kFALSE ;
}

//____________________________________________________________________________________________________
void FastGlauber::DoHardCoreCollision()
{
  /// Hard-core collision
  mCollisionProfile = mkHardCoreProfile ;
}

//____________________________________________________________________________________________________
void FastGlauber::DoGaussianCollision()
{
  /// Gaussian profile collision
  mCollisionProfile = mkGaussianProfile ;
}

//____________________________________________________________________________________________________
const Char_t* FastGlauber::GetName(const UInt_t massNumber) const
{
  switch ( massNumber ){
    case 63: return "Cu";
    case 144: return "Sm";
    case 154: return "Sm";
    case 197: return "Au";
    case 208: return "Pb";
    case 238: return "U";
    default:
      return "";
  }

  return "";
}

//____________________________________________________________________________________________________
Double_t FastGlauber::GetInelasticNNCrossSection(const Double_t energy, const TString type) const
{
  /// Cast to integer
  const UInt_t energyInt = static_cast<UInt_t>(energy) ;
  /// Assign error +/- 1mb if type is either largeXsec (+1mb) or smallXsec (-1mb)
  Double_t error = 0.0 ; // absolute error on cross section
  if( type.CompareTo("smallXsec", TString::kIgnoreCase) == 0 )      error = -0.1 ;
  else if( type.CompareTo("largeXsec", TString::kIgnoreCase) == 0 ) error = 0.1 ;

  Double_t sigmaNN = 0.0 ;
  switch ( energyInt ) {
    case 2760: sigmaNN = 6.4 + error ; break ; // 64 mb
    case 200: sigmaNN = 4.2 + error ; break ; // 42 mb
    case 62:  sigmaNN = 3.6 + error ; break ; // 36 mb
    case 39:  sigmaNN = 3.4 + error ; break ; // 34 mb
    case 27:  sigmaNN = 3.3 + error ; break ; // 33 mb
    case 19:  sigmaNN = 3.2 + error ; break ; // 32 mb
    case 14:  sigmaNN = 3.15 + error; break ; // 31.5 mb (for 14.5 GeV)
    case 11:  sigmaNN = 3.12 + error; break ; // 31.2 mb (for 11.5 GeV)
    case 7:   sigmaNN = 3.08 + error; break ; // 30.8 mb (for 7.7 GeV)
    default:
      assert(0);
  }

  return sigmaNN ;
}

//____________________________________________________________________________________________________
void FastGlauber::GetRThetaPhi(const UInt_t inucleus, Double_t& r, Double_t& theta, Double_t& phi) const
{
    // Spherical nuclei
  r     = mfWoodsSaxon[inucleus]->GetRandom() ;
  theta = GlauberUtilities::instance()->GetTheta() ; // cos(theta) profile in 0 < theta < pi
  phi = GlauberUtilities::instance()->GetPhi() ; // flat phi in -pi < phi < pi

  // Smearing if the switch is ON
  Smearing(r, theta, phi) ;
}

//____________________________________________________________________________________________________
void FastGlauber::Smearing(Double_t& r, Double_t& theta, Double_t& phi) const
{
  // Return original if both smearing options are OFF
  if(!mHardCoreSmearing && !mGaussianSmearing) return;

  // Hard-core smearing
  if(mHardCoreSmearing){
    Double_t dx, dy, dz ; /// Smearing width
    mfHardCore->GetRandom3(dx, dy, dz) ;

    const Double_t x = r * TMath::Sin(theta) * TMath::Cos(phi) + dx ;
    const Double_t y = r * TMath::Sin(theta) * TMath::Sin(phi) + dy ;
    const Double_t z = r * TMath::Cos(theta) + dz ;

    r     = TMath::Sqrt(x*x + y*y + z*z) ;
    theta = TMath::ACos(z/r) ;
    phi   = TMath::ATan2(y, x) ;
    return;
  }

  // Gaussian smearing
  if(mGaussianSmearing){
    Double_t dx, dy, dz ; /// Smearing width
    mfGaussian->GetRandom3(dx, dy, dz) ;

    const Double_t x = r * TMath::Sin(theta) * TMath::Cos(phi) + dx ;
    const Double_t y = r * TMath::Sin(theta) * TMath::Sin(phi) + dy ;
    const Double_t z = r * TMath::Cos(theta) + dz ;

    r     = TMath::Sqrt(x*x + y*y + z*z) ;
    theta = TMath::ACos(z/r) ;
    phi   = TMath::ATan2(y, x) ;
    return;
  }
}

//____________________________________________________________________________________________________
TGraph* FastGlauber::GetParticipantEccentricity(const Double_t order, const Double_t sumx, const Double_t sumy,
    const Double_t sumw, const UInt_t weightId) const
{
  Double_t qx = 0.0 ;
  Double_t qy = 0.0 ;
  Double_t qw = 0.0 ;

  for(UInt_t in=0;in<2;in++){
    for(UInt_t j=0;j<mNucleons[in].size();j++){
      Nucleon* nucleon = mNucleons[in][j] ;
      const UInt_t ncoll = nucleon->GetNcoll() ;

      const Bool_t isOk = (weightId==3 && ncoll==0) // Spectator
        || (weightId!=3 && ncoll > 0) ; // Participants

      if(!isOk) continue ;

      const Double_t x = nucleon->GetX() - sumx ;
      const Double_t y = nucleon->GetY() - sumy ;
      const Double_t z = nucleon->GetZ() ;

      const TVector3 vpart(x, y, z);
      const Double_t rt  = vpart.Perp() ;
      const Double_t phi = vpart.Phi() ;
      Double_t weight = 1.0 ;
      if( weightId == 1 ) weight = ncoll ;
      if( weightId == 2 ) weight = nucleon->GetMultiplicity() ; 
      if( weightId == 3 ) weight = 1.0 ;

      qx += -weight * rt*rt * TMath::Cos(order*phi) ; // -sign added to keep the same definition of participant plane
      qy +=  weight * rt*rt * TMath::Sin(order*phi) ;
      qw +=  weight * rt*rt ;
    }
  }

  qx /= sumw ;
  qy /= sumw ;
  qw /= sumw ;

  TGraph* g = new TGraph() ; // needs to be deleted later
  g->SetPoint(0, 0, qx);
  g->SetPoint(1, 1, qy);
  g->SetPoint(2, 2, TMath::ATan2(qy, qx));
  g->SetPoint(3, 3, TMath::Sqrt(qx*qx + qy*qy)/qw) ;

  return g ;
}

