#include "../include/gen-ana.h"
#include "../dflay/src/JSONManager.cxx"


TF1 *neutron_yield(TH1D *hdx, TString config){

  double p_low = 0;
  double p_high = 0;
  double n_low = 0;
  double n_high = 0;
  double bg_low = 0;
  double bg_high = 0;
  double par[11];

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
    
    n_low = -0.4;
    n_high = 0.4;
    
    bg_low = -4.0;
    bg_high = 3.0;

    par[0] = 600;
    par[1] = -1.5;
    par[2] = 0.4;

    par[3] = 150;
    par[4] = 0;
    par[5] = 0.4;
  }

  
  //TF1 *p_xfunc = new TF1("p_xfunc","gaus",p_low,p_high);
  TF1 *p_xfunc = new TF1("p_xfunc","gaus",p_low,p_high);
  TF1 *n_xfunc = new TF1("n_xfunc","gaus",n_low,n_high);  
  TF1 *bg_xfunc = new TF1("bg_xfunc","pol4",bg_low,bg_high);
  TF1 *total_xfunc = new TF1("total_xfunc","gaus(0) + gaus(3) + pol4(6)",bg_low,bg_high);

  p_xfunc->SetParameters(&par[0]);
  n_xfunc->SetParameters(&par[3]);
  

  hdx->Fit(p_xfunc,"qNR");  
  hdx->Fit(n_xfunc,"qNR+");  
  hdx->Fit(bg_xfunc,"qNR+"); 

  p_xfunc->GetParameters(&par[0]);
  n_xfunc->GetParameters(&par[3]);
  bg_xfunc->GetParameters(&par[6]);

  total_xfunc->SetParameters(par);
  hdx->Fit(total_xfunc,"qNR+"); 
  total_xfunc->GetParameters(&par[0]);

  par[0] = abs(par[0]);
  par[3] = abs(par[3]);

  p_xfunc = new TF1("p_xfunc","gaus",-4,4);
  p_xfunc->SetParameters(&par[0]);

  n_xfunc = new TF1("n_xfunc","gaus",-4,4);
  n_xfunc->SetParameters(&par[3]);
 
 
  return n_xfunc;

}

double Yield(TF1 *fit, TH1D *hdx){

  double binw = hdx->GetBinWidth(0);

  return fit->Integral(-4,4)/binw;

}

void do_fit_pol4(TCanvas *c, TH1D *h1,double bg_low, double bg_high, double signal_low,double signal_high){

  c->cd();
  double par[8];
  TF1 *gaus_func = new TF1("gaus_func","gaus",signal_low,signal_high);
  TF1 *bg_func = new TF1("bg_func","pol4",bg_low,bg_high);
  TF1 *total_func = new TF1("total_func","gaus(0) + pol4(3)",bg_low,bg_high);

  h1->Fit("gaus_func","qNR");
  h1->Fit("bg_func","qNR+");

  gaus_func->GetParameters(&par[0]);
  bg_func->GetParameters(&par[3]);

  total_func->SetParameters(par);
  h1->Fit(total_func,"qNR+"); 
  total_func->GetParameters(&par[0]);
  //gaus_func->SetParameters(&par[0]);

  total_func->Draw("same");

  cout<<par[1]<<" +/- "<<par[2]<<endl;

}

void do_fit_pol1(TCanvas *c, TH1D *h1,double bg_low, double bg_high, double signal_low,double signal_high){

  c->cd();
  double par[5];
  TF1 *gaus_func = new TF1("gaus_func","gaus",signal_low,signal_high);
  TF1 *bg_func = new TF1("bg_func","pol1",bg_low,bg_high);
  TF1 *total_func = new TF1("total_func","gaus(0) + pol1(3)",bg_low,bg_high);

  h1->Fit("gaus_func","qNR");
  h1->Fit("bg_func","qNR+");

  gaus_func->GetParameters(&par[0]);
  bg_func->GetParameters(&par[3]);

  total_func->SetParameters(par);
  h1->Fit(total_func,"qNR+"); 
  total_func->GetParameters(&par[0]);
  //gaus_func->SetParameters(&par[0]);

  total_func->Draw("same");

  cout<<par[1]<<" +/- "<<par[2]<<endl;

}


void coinc_analysis(TString cfg = "GEN2"){

  TFile *H2_file = new TFile("outfiles/QE_test_" + cfg + "_sbs100p_nucleon_p_model2_data.root","read");
  TFile *He3_file = new TFile("outfiles/QE_test_" + cfg + "_sbs100p_nucleon_np_model2_data.root","read");

  TTree *T = (TTree*)H2_file->Get("Tout");

  T->SetBranchStatus("*",0);

  const int maxClus = 1000;
  int runnum;   setrootvar::setbranch(T,"runnum","",&runnum);
  bool WCut;   setrootvar::setbranch(T,"WCut","",&WCut);
  bool pCut;   setrootvar::setbranch(T,"pCut","",&pCut);
  bool nCut;   setrootvar::setbranch(T,"nCut","",&nCut);
  bool coinCut;   setrootvar::setbranch(T,"coinCut","",&coinCut);
  double W2;   setrootvar::setbranch(T,"W2","",&W2);
  double dx;   setrootvar::setbranch(T,"dx","",&dx);
  double dy;   setrootvar::setbranch(T,"dy","",&dy);
  double coin_time;   setrootvar::setbranch(T,"coinT_trig","",&coin_time);
  double hcal_time;   setrootvar::setbranch(T,"hcal_time","",&hcal_time);
  double bbcal_time;   setrootvar::setbranch(T,"bbcal_time","",&bbcal_time);
  int nhodo_clus;   setrootvar::setbranch(T,"nhodo_clus","",&nhodo_clus);
  double hodo_time[maxClus];   setrootvar::setbranch(T,"hodo_time","",&hodo_time);
  int helicity;   setrootvar::setbranch(T,"helicity","",&helicity);
  int IHWP;   setrootvar::setbranch(T,"IHWP","",&IHWP);
  
  TH1D *h_coinc_H2_trig = new TH1D("h_coinc_H2_trig","",150,400,650);
  TH1D *h_coinc_H2_cal = new TH1D("h_coinc_H2_cal","",150,-100,150);
  TH1D *h_coinc_H2_hodo = new TH1D("h_coinc_H2_hodo","",150,-100,150);

  int nevent = 0;

  while(T->GetEntry(nevent++)){
    double Wrecon = sqrt(max(0., W2));

    if(WCut){
      h_coinc_H2_trig->Fill(coin_time);
      h_coinc_H2_cal->Fill(bbcal_time - hcal_time);
      h_coinc_H2_hodo->Fill(hodo_time[0] - hcal_time);
    }

  }

  
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
  setrootvar::setbranch(T,"hcal_time","",&hcal_time);
  setrootvar::setbranch(T,"bbcal_time","",&bbcal_time);
  setrootvar::setbranch(T,"nhodo_clus","",&nhodo_clus);
  setrootvar::setbranch(T,"hodo_time","",&hodo_time);
  setrootvar::setbranch(T,"helicity","",&helicity);
  setrootvar::setbranch(T,"IHWP","",&IHWP);
  
  TH1D *h_coinc_He3_trig = new TH1D("h_coinc_He3_trig","",150,-50,150);
  TH1D *h_coinc_He3_cal = new TH1D("h_coinc_He3_cal","",150,-100,100);
  TH1D *h_coinc_He3_hodo = new TH1D("h_coinc_He3_hodo","",150,-100,100);

  TH1D *hdx_He3_nocut = new TH1D("hdx_He3_nocut","",150,-6,4);
  TH1D *hdx_He3_trig = new TH1D("hdx_He3_trig","",150,-6,4);
  TH1D *hdx_He3_hodo = new TH1D("hdx_He3_hodo","",150,-6,4);
  TH1D *hdx_He3_Wtest = new TH1D("hdx_He3_Wtest","",150,-6,4);

  nevent = 0;

  while(T->GetEntry(nevent++)){
    double Wrecon = sqrt(max(0., W2));
    
    double coin_hodo = hodo_time[0] - hcal_time;

   if(WCut){

      h_coinc_He3_trig->Fill(bbcal_time - hcal_time);
      h_coinc_He3_cal->Fill(bbcal_time - hcal_time);
      h_coinc_He3_hodo->Fill(hodo_time[0] - hcal_time);
    }

   //if(Wrecon > 0.74 && Wrecon < 1.14) {
   if(WCut) {
     hdx_He3_nocut->Fill(dx);
     if(coin_time > 488 && coin_time < 514) 
       hdx_He3_trig->Fill(dx);
     if(coin_hodo > 3.5 && coin_hodo < 16.3) 
       hdx_He3_hodo->Fill(dx);  
   }

   if(Wrecon > 0.8 && Wrecon < 1.1) 
     if(coin_hodo > 1 && coin_hodo < 17) 
       hdx_He3_Wtest->Fill(dx);  

  }


  TCanvas *c = new TCanvas("c","",800,600);
  h_coinc_He3_trig->Draw();
  h_coinc_He3_trig->SetTitle("He3 Coincidence Time;HCal trigger - BBCal trigger (ns); Entries");
  do_fit_pol4(c,h_coinc_He3_trig,-40,140,30,45);

  TCanvas *c1 = new TCanvas("c1","",800,600);
  h_coinc_He3_cal->Draw();
  do_fit_pol4(c1,h_coinc_He3_cal,-200,200,25,55);

  TCanvas *c2 = new TCanvas("c2","",800,600);
  h_coinc_He3_hodo->Draw();
  h_coinc_He3_hodo->SetTitle("He3 Coincidence Time;Hodo Cluster Time - HCal Cluster Time (ns); Entries"); 
  do_fit_pol4(c2,h_coinc_He3_hodo,-200,200,-5,20);

  
  TCanvas *c3 = new TCanvas("c3","",800,600);
  hdx_He3_hodo->SetLineColor(kRed);
  hdx_He3_trig->SetLineColor(kGreen);
  hdx_He3_nocut->SetTitle("HCal p/n Spots;#Deltax;Entries");
  hdx_He3_nocut->Draw();    
  hdx_He3_hodo->Draw("same");  
  hdx_He3_trig->Draw("same");
  //TF1 *fit_trig = neutron_yield(hdx_He3_trig, cfg);  
  //fit_trig->Draw("same");
  //cout<<Yield(fit_trig,hdx_He3_trig)<<endl;

  //TF1 *fit_hodo = neutron_yield(hdx_He3_hodo, cfg);  
  //fit_hodo->Draw("same");
  //cout<<Yield(fit_hodo,hdx_He3_hodo)<<endl;

  TLegend *legend = new TLegend(0.6,0.7,0.89,0.89);
  legend->AddEntry("hdx_He3_nocut","No Timing Cuts","l");
  legend->AddEntry("hdx_He3_trig","HCal/BBCal Trig Diff < 2#sigma","l");
  legend->AddEntry("hdx_He3_hodo","HCal/Hodo Clus Diff < 2#sigma","l");
  legend->SetLineColor(0);
  legend->Draw("same");

}
