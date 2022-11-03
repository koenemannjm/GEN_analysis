
double bg_low;
double bg_high;

double p_low;
double p_high;

double n_low;
double n_high;

double bg_fit(double *x, double *par){

  if(x[0] > bg_low && x[0] < bg_high)
    TF1::RejectPoint();

  return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2));
}

void get_np_spots(TCanvas *c, TString cfg, TFile *file){
  
  TH2D *hdxdy_Wcut = (TH2D*)file->Get("hdxdy_Wcut");

  TH1D *hdx = hdxdy_Wcut->ProjectionY();

  p_low = -3.8;
  p_high = -1.0;

  n_low = -0.8;
  n_high = 1.5;

  double pF_par0 = 9000;
  double pF_par1 = -2.4;
  double pF_par2 = 0.3;
  
  double par[3];
  TF1 *p_xfunc = new TF1("p_xfunc","gaus",p_low,p_high);
  //p_xfunc->SetParameters(pF_par0,pF_par1,pF_par2);
  hdx->Fit(p_xfunc,"qR+");  
  p_xfunc->GetParameters(&par[0]);

  double px_mean = par[1];
  double px_sigma = par[2];

  TF1 *n_xfunc = new TF1("n_xfunc","gaus",n_low,n_high);
  hdx->Fit(n_xfunc,"qR+");  
  n_xfunc->GetParameters(&par[0]);

  double nx_mean = par[1];
  double nx_sigma = par[2];

  double x_cut = -1.0;
  int x_cut_bin = hdxdy_Wcut->GetYaxis()->FindBin(x_cut);

  TH1D *hdy_p = hdxdy_Wcut->ProjectionX("dy_p",0,x_cut_bin);
  TH1D *hdy_n = hdxdy_Wcut->ProjectionX("dy_n",x_cut_bin,-1);


  TF1 *p_yfunc = new TF1("p_yfunc","gaus",-1.5,1.5);
  hdy_p->Fit(p_yfunc,"qR+");  
  p_yfunc->GetParameters(&par[0]);

  double py_mean = par[1];
  double py_sigma = par[2];

  TF1 *n_yfunc = new TF1("n_yfunc","gaus",-1.5,1.5);
  hdy_n->Fit(n_yfunc,"qR+");  
  n_yfunc->GetParameters(&par[0]);

  double ny_mean = par[1];
  double ny_sigma = par[2];

  double nsigma = 1;

  if(cfg == "H2") nsigma = 1;

  TEllipse *p_spot = new TEllipse(py_mean,px_mean,nsigma*py_sigma,nsigma*px_sigma);
  TEllipse *n_spot = new TEllipse(ny_mean,nx_mean,nsigma*ny_sigma,nsigma*nx_sigma);
  p_spot->SetFillStyle(0);
  n_spot->SetFillStyle(0);

  p_spot->SetLineWidth(4);
  n_spot->SetLineWidth(4);

  p_spot->SetLineColor(kRed);
  n_spot->SetLineColor(kRed);

  c->cd();

  gStyle->SetOptStat(0);

  hdxdy_Wcut->SetTitle(cfg + " HCal Elastics;#Deltay (m);#Deltax (m)");
  hdxdy_Wcut->Draw("colz");
  p_spot->Draw("same");
  if(cfg == "He3") n_spot->Draw("same");


  TPaveText *pt = new TPaveText(.65,.8,.88,.88,"ndc");
  pt->AddText("Cuts on good tracks");
  pt->AddText("|W^{2} - 0.88| < 0.4");
  pt->SetFillColor(0);
  pt->Draw("same");

  cout<<cfg<<" data:"<<endl;
  cout<<"proton dx mean = "<<px_mean<<endl;
  cout<<"proton dx sigma = "<<px_sigma<<endl;
  cout<<"neutron dx mean = "<<nx_mean<<endl;
  cout<<"neutron dx sigma = "<<nx_sigma<<endl;

  cout<<"proton dy mean = "<<py_mean<<endl;
  cout<<"proton dy sigma = "<<py_sigma<<endl;
  cout<<"neutron dy mean = "<<ny_mean<<endl;
  cout<<"neutron dy sigma = "<<ny_sigma<<endl;
  cout<<"\n\n";
}

void Elastic_hist(){


  TFile *H2_file = new TFile("elastic_temp_H2.root","read");
  TFile *He3_file = new TFile("elastic_temp_He3.root","read");
  
  TCanvas *c1 = new TCanvas("c1","",800,1000);  
  get_np_spots(c1, "H2", H2_file);

  TCanvas *c2 = new TCanvas("c2","",800,1000);  
  get_np_spots(c2, "He3", He3_file);

  TH2D *hdxdy_H2 = (TH2D*)H2_file->Get("hdxdy_Wcut");
  TH1D *hdx_H2 = hdxdy_H2->ProjectionY("hdx_H2");
  hdx_H2->Scale(1/hdx_H2->GetEntries());
  
  TH2D *hdxdy_He3 = (TH2D*)He3_file->Get("hdxdy_Wcut");
  TH1D *hdx_He3 = hdxdy_He3->ProjectionY("hdx_He3");
  hdx_He3->SetLineColor(kRed);
  hdx_He3->Scale(1/hdx_He3->GetEntries());

  TH1D *hW2_cut_H2 = (TH1D*)H2_file->Get("hW2_cut");
  TH1D *hW2_cut_He3 = (TH1D*)He3_file->Get("hW2_cut");

  TH1D *hW2_all_H2 = (TH1D*)H2_file->Get("hW2_all");
  TH1D *hW2_all_He3 = (TH1D*)He3_file->Get("hW2_all");

  TCanvas *c3 = new TCanvas("c3","",1000,800);  
  hdx_H2->SetTitle("HCal p/n Spots;#Deltax(m);Normalized Entries");  
  hdx_H2->Draw("hist");
  hdx_He3->Draw("same hist");

  TLegend *legend = new TLegend(0.6,0.7,0.89,0.89);
  legend->AddEntry("hdx_H2","H2 Data","l");
  legend->AddEntry("hdx_He3","He3 Data","l");
  legend->SetLineColor(0);
  legend->Draw("same");

  TCanvas *c4 = new TCanvas("c4","",1000,800);  
  hW2_all_H2->SetTitle("H2 Elastic Data;W^{2} (GeV^{2});Entries");  
  hW2_all_H2->Draw();
  hW2_cut_H2->SetLineColor(kRed);
  hW2_cut_H2->Draw("same");

  TLegend *legend2 = new TLegend(0.11,0.75,0.38,0.89);
  legend2->AddEntry("hW2_all","All Events","l");
  legend2->AddEntry("hW2_cut","Events inside proton spot","l");
  legend2->SetLineColor(0);
  legend2->Draw("same");

  TCanvas *c5 = new TCanvas("c5","",1000,800);  
  hW2_all_He3->SetTitle("He3 Elastic Data;W^{2} (GeV^{2});Entries");  
  hW2_all_He3->Draw();
  hW2_cut_He3->SetLineColor(kRed);
  hW2_cut_He3->Draw("same");

  TLegend *legend3 = new TLegend(0.11,0.75,0.38,0.89);
  legend3->AddEntry("hW2_all","All Events","l");
  legend3->AddEntry("hW2_cut","Events inside p/n spot","l");
  legend3->SetLineColor(0);
  legend3->Draw("same");

  TString outputfile = "../plots/GEN_elastics.pdf";
  
  c1->Print(outputfile + "(");
  c2->Print(outputfile);
  c3->Print(outputfile);
  c4->Print(outputfile);
  c5->Print(outputfile + ")");

}
