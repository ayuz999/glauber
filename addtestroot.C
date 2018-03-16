{
  TFile *f = new TFile("/star/u/yangzz/yzz/22GeV_cucu/basic_dis/basic_dis.root");
  TH1D *k11=((TH1D *)f->Get("refmult3_dis"))->Clone("hRefMultTpc");

  TFile *f1=new TFile("refmult3.root","recreate");

  f1->cd();
  k11->Write("hRefMultTpc");
  f1->Write();
  f1->Close();
}


