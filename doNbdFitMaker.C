
//____________________________________________________________________________________________________
void doNbdFitMaker(
	const Char_t* outputFileName = "",
	const Int_t nevents = 1000,
	const Char_t* realData = "refmult3.root",
	const Char_t* glauber  = "ncoll_npart.root",
	const Double_t multiplicityCut = 40,
	const Double_t npp = 1.18,
	const Double_t k = 1.0,
	const Double_t x = 0.13,
	const Double_t efficiency = 1.00,
	const Double_t triggerbias = 1.00,
	const Bool_t isConstEfficiency = kTRUE
	){
  gSystem->Load("St_base");
  gSystem->Load("StUtilities");
  gSystem->Load("StGlauberUtilities");
  gSystem->Load("StCentralityMaker");

  StNbdFitMaker* maker = new StNbdFitMaker();

  // Set parameters
  maker->SetParameters(npp, k, x, efficiency, triggerbias, isConstEfficiency);

  // Set low multiplicity cut off for fitting
  maker->SetMinimumMultiplicityCut(multiplicityCut);

  // Read input files
 //maker->ReadData(realData, glauber, "refmult3.root");

 maker->ReadData(realData, glauber);

  // Fit
  maker->Fit(nevents, outputFileName);
}

//____________________________________________________________________________________________________
void scan(
	const Int_t nevents = 200000,
	const Char_t* realData = "refmult3.root",
	const Char_t* glauber  = "ncoll_npart.root",
	const Double_t multiplicityCut =40,
	const Int_t nppbin = 40, const Double_t nppmin = 1.2, const Double_t nppmax = 1.6,
	const Int_t kbin = 1, const Double_t kmin = 2.0, const Double_t kmax = 2.0,

	const Int_t xbin = 50, const Double_t xmin = 0.08, const Double_t xmax = 0.18,
	const Double_t efficiency = 1.00,
	const Double_t triggerbias = 1.00,
	const Bool_t isConstEfficiency = kTRUE
	){
  gSystem->Load("St_base");
  gSystem->Load("StUtilities");
  gSystem->Load("StGlauberUtilities");
  gSystem->Load("StCentralityMaker");

  StNbdFitMaker* maker = new StNbdFitMaker();

  // Set dummy parameters
  maker->SetParameters(nppmin, kmin, xmin, efficiency, triggerbias, isConstEfficiency);

  // Set low multiplicity cut off for fitting
  maker->SetMinimumMultiplicityCut(multiplicityCut);

  // Read input files
  maker->ReadData(realData, glauber);
  //maker->ReadData(realData, glauber,"refmult3.root");

  // Fit
  maker->Scan(nevents, nppbin, nppmin, nppmax, kbin, kmin, kmax, xbin, xmin, xmax,
	  efficiency, triggerbias, isConstEfficiency);
}

//____________________________________________________________________________________________________
void simulateMultiplicity(
	const Char_t* outputFileName = "test.root",
	const Int_t nevents = 500000,
	const Char_t* glauber  = "ncoll_npart.root",
	const Double_t npp = 1.18,
	const Double_t k = 2.0,
	const Double_t x = 0.20,
	const Double_t efficiency = 1.00,
	const Double_t triggerbias = 1.00,
	const Bool_t isConstEfficiency = kTRUE
	){
  gSystem->Load("St_base");
  gSystem->Load("StUtilities");
  gSystem->Load("StGlauberUtilities");
  gSystem->Load("StCentralityMaker");

  StNbdFitMaker* maker = new StNbdFitMaker();

  // Set parameters
  maker->SetParameters(npp, k, x, efficiency, triggerbias, isConstEfficiency);

  // Read input files
  maker->ReadGlauber(glauber, 3000, 0, 3000);

  // Fit
  maker->Simulate(nevents, outputFileName);
}

