#ifndef __FastGlauber_h__
#define __FastGlauber_h__

class TF1 ;
class TF2 ;
class TF3 ;
class TGraph ;
class TH1 ;
class TTree ;
class Nucleon ;
class CentralityMaker ;
class GlauberTree ;
class Rtypes;
#include <vector>
#include <string>

using namespace std;

class FastGlauber {
  public:
    FastGlauber();
    // Current available type
    //     default
    //     large              Large R(+2%), small d(-10%)
    //     small              Small R(-2%), large d(+10%)
    //     largeXsec          Large inelastic NN cross section (+1mb)
    //     smallXsec          Small inelastic NN cross section (-1mb)
    //     gauss              Use gaussian collision profile
    FastGlauber(
        const string outputFileName, // Output fileName
        const int    energy,// energy (GeV)
        const string type  // type (see above)
        );

    virtual ~FastGlauber(){}; /// Default destructor

    /// Set repulsion of nucleons (default is 0fm)
    void SetRepulsionDistance(const Double_t repulsionDistance);

    Int_t Make() ; /// Make one event
    Int_t Run(const UInt_t nevents) ; /// Run Make() by nevents
    Int_t Finish() ; /// Finish maker

    /// Hard-core smearing by sigmaNN
    void DoHardCoreSmearing() ; /// Default is OFF

    /// Gaussian smearing by width = 0.79/sqrt(3) from CPC180, 69, 2009
    void DoGaussianSmearing() ; /// Default is OFF

    /// Collision profiles
    void DoHardCoreCollision() ; /// Hard-core collision (default)
    void DoGaussianCollision() ; /// Gaussion profile collision

  private:

    Int_t Clear() ; /// Clear all data members in tree
    Int_t InitTree() ;
    /// Initialization of nucleus and nucleons, and output ROOT file
    Int_t Init(
        const UInt_t massNumberA,  // Mass number of nucleus for nucleus A
        const Double_t radiusA,    // Radius of nucleus for nucleus A
        const Double_t skinDepthA, // Skin depth of nucleus for nucleus A
        const UInt_t massNumberB,  // Mass number of nucleus for nucleus B
        const Double_t radiusB,    // Radius of nucleus for nucleus B
        const Double_t skinDepthB, // Skin depth of nucleus for nucleus B
        const TString type = "default"
        );
    /// Initialization for specific collisions
    Int_t InitAuAu(const TString type) ;
    /// Nucleon-Nucleon collision
    Bool_t IsCollision(Nucleon* nucleon0, Nucleon* nucleon1) const ;
    /// Utilities
    const Char_t* GetName(const UInt_t massNumber) const ; /// Mass number -> nucleus name
    Double_t GetInelasticNNCrossSection(const Double_t energy, const TString type) const ;
    void GetRThetaPhi(const UInt_t inucleus,
        Double_t& r, Double_t& theta, Double_t& phi) const ; /// Arguments are outputs, not inputs

    /// Do smearing (either hard-core or gaussian depending on the switch)
    void Smearing(Double_t& r, Double_t& theta, Double_t& phi) const ; /// NOTE: (r,theta,phi) is output

    /// Participant eccentricity: eps_{part;n} = sqrt(q_{x;n}^2 + q_{y;n}^2)/r_t^2
    /// q_{x;n} = <r_T^2 * cos(n * phi_part)>
    /// q_{y;n} = <r_T^2 * sin(n * phi_part)>
    /// r_T = sqrt(x_{part}^2 + y_{part}^2)
    /// phi_part is the azimth of nucleons with respect to the participant plane
    //  
    //  TGraph::GetY()[0]    <r_T^2 * cos(n*phi)>
    //  TGraph::GetY()[1]    <r_T^2 * sin(n*phi)>
    //  TGraph::GetY()[2]    n * Psi = tan^{-1}(q_{y;n}, q_{x;n})
    //  TGraph::GetY()[3]    eps_{part;n}
    TGraph* GetParticipantEccentricity(const Double_t order, const Double_t sumx, const Double_t sumy,
        const Double_t sumw, const UInt_t weightId) const ;

    //____________________________________________________________________________________________________
    // Data members

    // For collision profile
    enum {
      mkHardCoreProfile = 0,
      mkGaussianProfile = 1
    };


    const Double_t mEnergy ; /// sqrt(sNN)
    const TString mOutputFileName ; /// Output filename
    Double_t mInelasticNNCrossSection ; /// Inelastic NN cross section
    // Switches
    Double_t mRepulsionDistance ; /// Repulsion distance between nucleons (default is 0fm)
    // NOTE: Only one smearing method can be applied
    //       Another option will be switched off when one turn on one of them
    Bool_t mHardCoreSmearing    ; /// true  -> smear nucleons by step function with sigmaNN
                                  /// false -> Use original position (default)
    Bool_t mGaussianSmearing    ; /// true  -> smear nucleons by gaussian distribution with fixed sigma
                                  /// false -> Use original position (default)
    UInt_t mCollisionProfile    ; /// Collision profile flag (mkHardCoreProfile=0, mkGaussianProfile=1)

    UInt_t mMode ; /// 0:default, 1:large npp (small x), 2:small npp (large x)

    CentralityMaker* mCentralityMaker ; /// For multiplicity
    std::vector<Nucleon*> mNucleons[2] ; /// Nucleons
    TF1* mfWoodsSaxon[2]   ; /// Woods-saxon density profile for spherical nuclei
    TF2* mfWoodsSaxon2D[2] ; /// Woods-saxon density profile for deformed nuclei

    UInt_t mNeventsThrow  ; /// Number of all events
    UInt_t mNeventsAccept ; /// Number of accepted events (Ncoll>0)
    TTree* mHeader        ; /// Output ROOT tree (store constant info., R, d, sigma etc)

    // Output data in tree
    GlauberTree* mGlauberTree ; /// MC glauber tree class

    // QA histograms
    TH1* mhWoodsSaxon[4] ; /// Woods-saxon checks. First 2 histograms = all nucleons, Last 2 histograms = for replusive nucleons
                           /// If repulsion distance is 0fm, then First and Last two should be identical

    // Random number generator
    TF1* mfB        ; /// Impact parameter distribtion p(b) ~ b
    TF3* mfHardCore ; /// Hard-core smearing for (x,y,z) of nucleons
    TF3* mfGaussian ; /// Gaussian smearing for (x,y,z) of nucleons

};
#endif
