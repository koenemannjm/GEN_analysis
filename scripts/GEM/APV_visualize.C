

void APV_visualize(){

  ifstream APV_data("APV_data.txt");

  gStyle->SetOptStat(0);

  double danning = 476.1;
  double histogram = 426.1;
  double sorting = 468.6;

  TH1F *hAPV = new TH1F("hAPV","",128,0,128);

  int ADC;
  int ibin = 0;

  while(APV_data >> ADC){

    hAPV->SetBinContent(ibin++,ADC);
  }

  hAPV->Draw();

  TLine *DLine = new TLine(0,danning,128,danning);
  DLine->SetLineColor(kRed);

  TLine *HLine = new TLine(0,histogram,128,histogram);
  HLine->SetLineColor(kGreen);

  TLine *SLine = new TLine(0,sorting,128,sorting);
  SLine->SetLineColor(kBlue);

  DLine->Draw("same");
  HLine->Draw("same");
  SLine->Draw("same");

  TLegend *legend = new TLegend(0.6,0.7,0.89,0.89);
  legend->AddEntry(DLine,"Danning","l");
  legend->AddEntry(HLine,"Histogramming","l");
  legend->AddEntry(SLine,"Sorting","l");
  legend->SetLineColor(0);
  legend->Draw("same");
}
