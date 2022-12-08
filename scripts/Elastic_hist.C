#include "../include/gen-ana.h"
#include "../dflay/src/JSONManager.cxx"

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

void get_np_spots(TCanvas *c, TCanvas *c1, TString cfg, TH2D *hdxdy, TString config){


  TH1D *hdx = hdxdy->ProjectionY("hdx_fits");

  if(config == "GEN2"){
    p_low = -4.0;
    p_high = -1.8;
    
    n_low = -0.6;
    n_high = 0.8;
    
    bg_low = -4.0;
    bg_high = 3.0;
  }

  if(config == "GEN3"){
    p_low = -2.0;
    p_high = -0.8;
    
    n_low = -0.5;
    n_high = 0.5;
    
    bg_low = -4.0;
    bg_high = 3.0;
  }

  double par[11];
  //TF1 *p_xfunc = new TF1("p_xfunc","gaus",p_low,p_high);
  TF1 *p_xfunc = new TF1("p_xfunc","gaus",p_low,p_high);
  TF1 *n_xfunc = new TF1("n_xfunc","gaus",n_low,n_high);  
  TF1 *bg_xfunc = new TF1("bg_xfunc","pol4",bg_low,bg_high);
  TF1 *total_xfunc = new TF1("total_xfunc","gaus(0) + gaus(3) + pol4(6)",bg_low,bg_high);
  

  hdx->Fit(p_xfunc,"qNR");  
  hdx->Fit(n_xfunc,"qNR+");  
  hdx->Fit(bg_xfunc,"qNR+"); 

  p_xfunc->GetParameters(&par[0]);
  n_xfunc->GetParameters(&par[3]);
  bg_xfunc->GetParameters(&par[6]);

  total_xfunc->SetParameters(par);
  hdx->Fit(total_xfunc,"qNR+"); 
  total_xfunc->GetParameters(&par[0]);

  double px_mean = par[1];
  double px_sigma = par[2];
  double nx_mean = par[4];
  double nx_sigma = par[5];

  p_xfunc = new TF1("p_xfunc","gaus",-4,4);
  p_xfunc->SetParameters(&par[0]);

  n_xfunc = new TF1("n_xfunc","gaus",-4,4);
  n_xfunc->SetParameters(&par[3]);
 
  double x_cut = n_low;
  int x_cut_bin = hdxdy->GetYaxis()->FindBin(x_cut);

  TH1D *hdy_p = hdxdy->ProjectionX("dy_p",0,x_cut_bin);
  TH1D *hdy_n = hdxdy->ProjectionX("dy_n",x_cut_bin,-1);
  

  TF1 *p_yfunc = new TF1("p_yfunc","gaus",-1.5,1.5);
  hdy_p->Fit(p_yfunc,"qNR+");  
  p_yfunc->GetParameters(&par[0]);

  double py_mean = par[1];
  double py_sigma = par[2];

  TF1 *n_yfunc = new TF1("n_yfunc","gaus",-1.5,1.5);
  hdy_n->Fit(n_yfunc,"qNR+");  
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

  hdxdy->SetTitle(cfg + " HCal Elastics;#Deltay (m);#Deltax (m)");
  hdxdy->Draw("colz");
  p_spot->Draw("same");
  if(cfg == "He3") n_spot->Draw("same");


  TPaveText *pt = new TPaveText(.65,.8,.88,.88,"ndc");
  pt->AddText("Cuts on good tracks");
  pt->AddText("0.74 < W < 1.14");
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

  c1->cd();
  hdx->SetTitle("He3 HCal p/n Spots;#Deltax(m);Entries");  
  hdx->SetLineColor(kRed);
  hdx->Draw("same hist"); 
  total_xfunc->Draw("same");
  n_xfunc->Draw("same");
  p_xfunc->Draw("same");
 
}


void get_p_spots(TCanvas *c, TString cfg, TH2D *hdxdy,TString config){


  TH1D *hdx = hdxdy->ProjectionY();

  if(config == "GEN2"){
    p_low = -4.0;
    p_high = -1.8;
    
    n_low = -0.6;
    n_high = 0.8;
    
    bg_low = -4.0;
    bg_high = 3.0;
  }

  if(config == "GEN3"){
    p_low = -2.0;
    p_high = -0.8;
    
    n_low = -0.5;
    n_high = 0.5;
    
    bg_low = -4.0;
    bg_high = 3.0;
  }

  double par[11];
  //TF1 *p_xfunc = new TF1("p_xfunc","gaus",p_low,p_high);
  TF1 *p_xfunc = new TF1("p_xfunc","gaus",p_low,p_high);
  TF1 *n_xfunc = new TF1("n_xfunc","gaus",n_low,n_high);  
  TF1 *bg_xfunc = new TF1("bg_xfunc","pol4",bg_low,bg_high);
  TF1 *total_xfunc = new TF1("total_xfunc","gaus(0) + gaus(3) + pol4(6)",-4,3);
  

  hdx->Fit(p_xfunc,"qNR");  
  hdx->Fit(n_xfunc,"qNR+");  
  hdx->Fit(bg_xfunc,"qNR+"); 

  p_xfunc->GetParameters(&par[0]);
  n_xfunc->GetParameters(&par[3]);
  bg_xfunc->GetParameters(&par[6]);

  total_xfunc->SetParameters(par);
  hdx->Fit(total_xfunc,"qR+"); 
  total_xfunc->GetParameters(&par[0]);

  double px_mean = par[1];
  double px_sigma = par[2];
  double nx_mean = par[4];
  double nx_sigma = par[5];

  p_xfunc = new TF1("p_xfunc","gaus",-4,4);
  p_xfunc->SetParameters(&par[0]);

  n_xfunc = new TF1("n_xfunc","gaus",-4,4);
  n_xfunc->SetParameters(&par[3]);
  
  
  double x_cut = -1.0;
  int x_cut_bin = hdxdy->GetYaxis()->FindBin(x_cut);

  TH1D *hdy_p = hdxdy->ProjectionX("dy_p",0,x_cut_bin);
  TH1D *hdy_n = hdxdy->ProjectionX("dy_n",x_cut_bin,-1);
  

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

  
  double nsigma = 2;

  if(cfg == "H2") nsigma = 2;

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

  hdxdy->SetTitle(cfg + " HCal Elastics;#Deltay (m);#Deltax (m)");
  hdxdy->Draw("colz");
  p_spot->Draw("same");
  if(cfg == "He3") n_spot->Draw("same");


  TPaveText *pt = new TPaveText(.65,.8,.88,.88,"ndc");
  pt->AddText("Cuts on good tracks");
  pt->AddText("0.74 < W < 1.14");
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

void Elastic_hist(TString cfg = "GEN2"){


  TFile *H2_file = new TFile("outfiles/QE_test_" + cfg + "_sbs100p_nucleon_p_model2_data.root","read");
  TFile *He3_file = new TFile("outfiles/QE_test_" + cfg + "_sbs100p_nucleon_np_model2_data.root","read");

  TTree *T = (TTree*)H2_file->Get("Tout");

  T->SetBranchStatus("*",0);

  int runnum;   setrootvar::setbranch(T,"runnum","",&runnum);
  bool WCut;   setrootvar::setbranch(T,"WCut","",&WCut);
  bool pCut;   setrootvar::setbranch(T,"pCut","",&pCut);
  bool nCut;   setrootvar::setbranch(T,"nCut","",&nCut);
  bool coinCut;   setrootvar::setbranch(T,"coinCut","",&coinCut);
  double W2;   setrootvar::setbranch(T,"W2","",&W2);
  double dx;   setrootvar::setbranch(T,"dx","",&dx);
  double dy;   setrootvar::setbranch(T,"dy","",&dy);
  double coin_time;   setrootvar::setbranch(T,"coinT_trig","",&coin_time);
  int helicity;   setrootvar::setbranch(T,"helicity","",&helicity);
  int IHWP;   setrootvar::setbranch(T,"IHWP","",&IHWP);

  TH2D *hdxdy_nocut_H2 = new TH2D("h2_dxdy_nocut_H2","",150,-2,2,150,-6,6);
  TH2D *hdxdy_Wcut_H2 = new TH2D("h2_dxdy_Wcut_H2","",150,-2,2,150,-6,6);
  TH2D *hdxdy_coincut_H2 = new TH2D("h2_dxdy_coincut_H2","",150,-2,2,150,-6,6);
  TH1D *hW_all_H2 = new TH1D("hW_all_H2","",200,0,2);
  TH1D *hW_cut_H2 = new TH1D("hW_cut_H2","",200,0,2);

  int nevent = 0;

  while(T->GetEntry(nevent++)){
    double Wrecon = sqrt(max(0., W2));

    if(Wrecon > 0.01) hW_all_H2->Fill(Wrecon);
    if(pCut || nCut) hW_cut_H2->Fill(Wrecon);

    hdxdy_nocut_H2->Fill(dy,dx);
    if(WCut){
      hdxdy_Wcut_H2->Fill(dy,dx);
      if(coinCut) hdxdy_coincut_H2->Fill(dy,dx);
    }
  }
  
  TCanvas *c1 = new TCanvas("c1","",800,1000);  
  get_p_spots(c1, "H2", hdxdy_coincut_H2,cfg);

  T->Delete();

  T = (TTree*)He3_file->Get("Tout");

  T->SetBranchStatus("*",0);

  setrootvar::setbranch(T,"runnum","",&runnum);
  setrootvar::setbranch(T,"WCut","",&WCut);
  setrootvar::setbranch(T,"pCut","",&pCut);
  setrootvar::setbranch(T,"nCut","",&nCut);
  setrootvar::setbranch(T,"coinCut","",&coinCut);
  setrootvar::setbranch(T,"W2","",&W2);
  setrootvar::setbranch(T,"dx","",&dx);
  setrootvar::setbranch(T,"dy","",&dy);
  setrootvar::setbranch(T,"coinT_trig","",&coin_time);
  setrootvar::setbranch(T,"helicity","",&helicity);
  setrootvar::setbranch(T,"IHWP","",&IHWP);
  
  
  TH2D *hdxdy_nocut_He3 = new TH2D("h2_dxdy_nocut_He3","",150,-2,2,150,-6,6);
  TH2D *hdxdy_Wcut_He3 = new TH2D("h2_dxdy_Wcut_He3","",150,-2,2,150,-6,6);
  TH2D *hdxdy_coincut_He3 = new TH2D("h2_dxdy_coincut_He3","",150,-2,2,150,-6,6);
  TH1D *hW_all_He3 = new TH1D("hW_all_He3","",200,0,2);
  TH1D *hW_cut_He3 = new TH1D("hW_cut_He3","",200,0,2);
 
  nevent = 0;

  while(T->GetEntry(nevent++)){
    double Wrecon = sqrt(max(0., W2));

    if(Wrecon > 0.01) hW_all_He3->Fill(Wrecon);
    if(pCut || nCut) hW_cut_He3->Fill(Wrecon);
    
    hdxdy_nocut_He3->Fill(dy,dx);
    if(WCut){
      hdxdy_Wcut_He3->Fill(dy,dx);
      if(coinCut) hdxdy_coincut_He3->Fill(dy,dx);
    }
  }

  TCanvas *c2 = new TCanvas("c2","",800,1000);  
  TCanvas *c3 = new TCanvas("c3","",1000,800);  
 

  TH1D* hcoin = (TH1D*)He3_file->Get("h_coin_time");
  
  TH1D *hdx_H2 = hdxdy_coincut_H2->ProjectionY("");
  hdx_H2->Scale(1/hdx_H2->GetEntries());
  
  TH1D *hdx_He3 = hdxdy_coincut_He3->ProjectionY("hdx_He3");
  hdx_He3->SetLineColor(kRed);
  hdx_He3->Scale(1/hdx_He3->GetEntries());

  TH1D *hdx_nocut_He3 = hdxdy_nocut_He3->ProjectionY("hdx_nocut_He3");
  TH1D *hdx_Wcut_He3 = hdxdy_Wcut_He3->ProjectionY("hdx_Wcut_He3");
  TH1D *hdx_coincut_He3 = hdxdy_coincut_He3->ProjectionY("hdx_coincut_He3");
  hdx_coincut_He3->SetLineColor(kRed);

  
  c3->cd();
  hdx_Wcut_He3->Draw("same");
  get_np_spots(c2, c3, "He3", hdxdy_coincut_He3,cfg);  
  //hdx_coincut_He3->Draw("same hist");

  TLegend *legend = new TLegend(0.5,0.75,0.89,0.89);
  legend->AddEntry("hdx_Wcut_He3","Good Tracks","l");
  legend->AddEntry("hdx_fits","Good Tracks & Coincidence","l");
  legend->SetLineColor(0);
  legend->Draw("same");
  
  TCanvas *c4 = new TCanvas("c4","",1000,800);  
  hdx_H2->SetTitle("HCal p/n Spots;#Deltax(m);Normalized Entries");  
  hdx_H2->Draw("hist");
  hdx_He3->Draw("same hist");

  TLegend *legend2 = new TLegend(0.6,0.7,0.89,0.89);
  legend2->AddEntry("hdx_H2","H2 Data","l");
  legend2->AddEntry("hdx_He3","He3 Data","l");
  legend2->SetLineColor(0);
  legend2->Draw("same");

  TCanvas *c5 = new TCanvas("c5","",1000,800);  
  hW_all_H2->SetTitle("H2 Elastic Data;W (GeV);Entries");  
  hW_all_H2->Draw();
  hW_cut_H2->SetLineColor(kRed);
  hW_cut_H2->Draw("same");

  TLegend *legend3 = new TLegend(0.11,0.75,0.50,0.89);
  legend3->AddEntry("hW_all_H2","All Events","l");
  legend3->AddEntry("hW_cut_H2","Events inside proton spot","l");
  legend3->SetLineColor(0);
  legend3->Draw("same");

  TCanvas *c6 = new TCanvas("c6","",1000,800);  
  hW_all_He3->SetTitle("He3 Elastic Data;W (GeV);Entries");  
  hW_all_He3->Draw();
  hW_cut_He3->SetLineColor(kRed);
  hW_cut_He3->Draw("same");

  TLegend *legend4 = new TLegend(0.11,0.75,0.50,0.89);
  legend4->AddEntry("hW_all_He3","All Events","l");
  legend4->AddEntry("hW_cut_He3","Events inside p/n spot","l");
  legend4->SetLineColor(0);
  legend4->Draw("same");

  TString outputfile = "../plots/" + cfg + "_elastics.pdf";
  
  c1->Print(outputfile + "(");
  c2->Print(outputfile);
  c3->Print(outputfile);
  c4->Print(outputfile);
  c5->Print(outputfile);
  c6->Print(outputfile + ")");
  
}
