{

  gROOT->Reset();


  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFrameFillColor(10); 
  gStyle->SetPadBorderMode(0);
  gStyle->SetStatColor(0); 
  gStyle->SetOptFit(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillStyle(1001);
  gStyle->SetFrameFillColor(10);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetStatColor(0);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleXOffset(0.);
  gStyle->SetTitleYOffset(0.);
  gStyle->SetPadTickX(3);
  gStyle->SetPadTickY(3);
  gStyle->SetTitleTextColor(1);
  gStyle->SetEndErrorSize(0);


  TCanvas *c1 = new TCanvas("c1","c1",0,0,1000,900);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(0.1);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);
  c1->cd();


  TPad *pad = new TPad("pad", "pad",0.085,0.097,0.9,0.9);
  pad->SetBorderMode(0);
  pad->SetFillColor(kWhite);
  pad->Draw();
  pad->cd();
  pad->Divide(1,2,0.,0.,0.);



  pad->cd(1);
  gPad->SetLogy();
  //gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.01);
  gPad->SetLeftMargin(0.15);

  f1=new TFile("Ratio_npp0.42_k2.00_x0.12_eff0.14.root");
  float norm=1;

  k11=(TH1D *)f1->Get("hRefMultTpc");                                         
  k11->Scale(norm/k11->Integral(60,140));   
  //k11->Rebin(2);
  k11->SetLineColor(1);
  k11->SetLineStyle(1);
  k11->SetLineWidth(2);

  k21=(TH1D *)f1->Get("hRefMultSim");
  k21->Scale(norm/k21->Integral(60,140));
  //k21->Rebin(2);
  k21->SetLineColor(2);
  k21->SetLineStyle(2);
  k21->SetLineWidth(2);


  TH2D *ff1= new TH2D("","",290,0,290,100,0.0000011,0.5);
  ff1->GetXaxis()->SetNdivisions(505);
  ff1->GetXaxis()->SetLabelOffset(0.01);
  ff1->GetXaxis()->SetLabelSize(0.06);
  ff1->GetXaxis()->SetTitleSize(0.07);
  ff1->GetXaxis()->SetTickLength(0.03);
  ff1->GetXaxis()->SetLabelFont(42);


  ff1->GetYaxis()->SetLabelOffset(0.02);
  ff1->GetYaxis()->SetNdivisions(505);
  ff1->GetYaxis()->SetLabelSize(0.08);
  ff1->GetYaxis()->SetTickLength(0.03);
  ff1->GetYaxis()->SetTitleOffset(0.9);
  ff1->GetYaxis()->SetTitleSize(0.09);
  ff1->GetYaxis()->SetTitleFont(42);
  ff1->GetYaxis()->SetLabelFont(42);
  ff1->GetYaxis()->SetTitle("Counts");
  ff1->GetYaxis()->CenterTitle();
  ff1->Draw();

  k11->Draw("lsame");
  k21->Draw("lsame");

  TLine *line1=new TLine(88,0.0000011,88,0.03);
  line1->SetLineColor(1);
  line1->SetLineStyle(2);
  line1->Draw();

  TLatex *tex1=new TLatex(98,0.0001,"0-5%");
  tex1->SetTextSize(0.06);
  tex1->SetTextFont(42);
  tex1->SetTextAngle(90);
  tex1->Draw();

  TLine *line2=new TLine(74,0.0000011,74,0.03);
  line2->SetLineColor(1);
  line2->SetLineStyle(2);
  line2->Draw();

  TLatex *tex2=new TLatex(84,0.0001,"5-10%");
  tex2->SetTextSize(0.06);
  tex2->SetTextFont(42);
  tex2->SetTextAngle(90);
  tex2->Draw();

  TLine *line3=new TLine(51,0.0000011,51,0.04);
  line3->SetLineColor(1);
  line3->SetLineStyle(2);
  line3->Draw();

  TLatex *tex3=new TLatex(60,0.0001,"10-20%");
  tex3->SetTextSize(0.06);
  tex3->SetTextFont(42);
  tex3->SetTextAngle(90);
  tex3->Draw();

  TLine *line4=new TLine(35,0.0000011,35,0.06);
  line4->SetLineColor(1);
  line4->SetLineStyle(2);
  line4->Draw();

  TLatex *tex4=new TLatex(45,0.0001,"20-30%");
  tex4->SetTextSize(0.06);
  tex4->SetTextFont(42);
  tex4->SetTextAngle(90);
  tex4->Draw();

  TLine *line5=new TLine(23,0.0000011,23,0.1);
  line5->SetLineColor(1);
  line5->SetLineStyle(2);
  line5->Draw();

  TLatex *tex5=new TLatex(31,0.0001,"30-40%");
  tex5->SetTextSize(0.06);
  tex5->SetTextFont(42);
  tex5->SetTextAngle(90);
  tex5->Draw();


  TLegend *leg1= new TLegend(0.6,0.5,0.93,0.79,"Centrality Selection");
  leg1->AddEntry(k11,"Data","l");
  leg1->AddEntry(k21,"MC","l");
  leg1->SetTextAlign(22);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.08);
  leg1->SetLineColor(1);
  leg1->SetFillColor(kWhite);
  leg1->SetLineColor(1);
  leg1->SetBorderSize(1);
  leg1->Draw("lt");



  pad->Modified();
  c1->cd();


  pad->cd(2);
  //gPad->SetLogy();
  gPad->SetBottomMargin(0.19);
  gPad->SetRightMargin(0.01);
  gPad->SetLeftMargin(0.15);
  float norm=1;

  k11=(TH1D *)f1->Get("hRatio");
  //k11->Scale(norm/k11->Integral(60,140));
  //k11->Rebin(2);
  k11->SetLineColor(2);
  k11->SetLineStyle(1);
  k11->SetLineWidth(2);


  TH2D *ff1= new TH2D("","",290,0,290,100,0,2.2);
  ff1->GetXaxis()->SetNdivisions(505);
  ff1->GetXaxis()->SetLabelOffset(0.01);
  ff1->GetXaxis()->SetLabelSize(0.08);
  ff1->GetXaxis()->SetTitleSize(0.09);
  ff1->GetXaxis()->SetTickLength(0.03);
  ff1->GetXaxis()->SetLabelFont(42);
  ff1->GetXaxis()->SetTitleOffset(0.9);
  ff1->GetXaxis()->SetTitleFont(42);
  ff1->GetXaxis()->SetTitle("Refmult3");
  ff1->GetXaxis()->CenterTitle();


  ff1->GetYaxis()->SetLabelOffset(0.02);
  ff1->GetYaxis()->SetNdivisions(505);
  ff1->GetYaxis()->SetLabelSize(0.08);
  ff1->GetYaxis()->SetTickLength(0.03);
  ff1->GetYaxis()->SetTitleOffset(0.9);
  ff1->GetYaxis()->SetTitleSize(0.09);
  ff1->GetYaxis()->SetTitleFont(42);
  ff1->GetYaxis()->SetLabelFont(42);
  ff1->GetYaxis()->SetTitle("MC/data");
  ff1->GetYaxis()->CenterTitle();


  ff1->Draw();

  k11->Draw("lsame");
  TLatex *text1=new TLatex(70,1.7,"#splitline{Fixed-Target}{AuAu:4.5GeV}");
  text1->SetTextFont(42);
  text1->SetTextSize(0.06);
  text1->Draw("same");

  TLatex *text2=new TLatex(150,1.7,"npp=0.41, k=2.0, x=0.12");
  text2->SetTextFont(42);
  text2->SetTextSize(0.06);
  text2->Draw("same");

  TLatex *text3=new TLatex(150,1.5,"#chi^{2}/ndf=19.149");
  text3->SetTextFont(42);
  text3->SetTextSize(0.06);
  text3->Draw("same");

  TLatex *text4=new TLatex(150,1.3,"Fit: 60<RefMult3<140");
  text4->SetTextFont(42);
  text4->SetTextSize(0.06);
  text4->Draw("same");
  TLine *line=new TLine(0,1,290,1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw();

  TLine *line1=new TLine(60,0,60,2.2);
  line1->SetLineColor(1);
  line1->SetLineStyle(2);
  line1->SetLineWidth(1);
  line1->Draw();

  TLine *line2=new TLine(140,1,140,2.2);
  line2->SetLineColor(1);
  line2->SetLineStyle(2);
  line2->SetLineWidth(1);
  line2->Draw();

  pad->Modified();
  c1->cd();


}
