#include "../include/gen-ana.h"
#include "../dflay/src/JSONManager.cxx"

double bg_low;
double bg_high;

double p_low;
double p_high;

double n_low;
double n_high;

double W2min = 0;
double W2max = 1.6;


//Function to get fits for the p/n spots on HCal
void get_np_spots(TCanvas *c, TString cfg, TH2D *hdxdy,TString config, TF1 **fit_n,TF1 **fit_result){

  //First project to x
  TH1D *hdx = hdxdy->ProjectionY();

  //These parameters are determined by looking at the peaks on the plots by eye
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
    
    bg_low = -5.0;
    bg_high = 3.0;
  }

  //total function = proton gaus + neutron gaus + 4th order bkgd
  double par[11];
  TF1 *p_xfunc = new TF1("p_xfunc","gaus",p_low,p_high);
  TF1 *n_xfunc = new TF1("n_xfunc","gaus",n_low,n_high);  
  TF1 *bg_xfunc = new TF1("bg_xfunc","pol4",bg_low,bg_high);
  TF1 *total_xfunc = new TF1("total_xfunc","gaus(0) + gaus(3) + pol4(6)",-4,3);
  
  //Do fit but do not plot the results
  hdx->Fit(p_xfunc,"qNR");  
  hdx->Fit(n_xfunc,"qNR+");  
  hdx->Fit(bg_xfunc,"qNR+"); 

  //Put the fit parameters into the array
  p_xfunc->GetParameters(&par[0]);
  n_xfunc->GetParameters(&par[3]);
  bg_xfunc->GetParameters(&par[6]);

  //Set the parameters in the total function using the results above
  total_xfunc->SetParameters(par);
  hdx->Fit(total_xfunc,"qR+"); 
  
  //Get the fit results
  total_xfunc->GetParameters(&par[0]);

  double px_mean = par[1];
  double px_sigma = par[2];
  double nx_mean = par[4];
  double nx_sigma = par[5];

  //For plotting purposes set the p/n function parameters from the total fit
  p_xfunc = new TF1("p_xfunc","gaus",-4,4);
  p_xfunc->SetParameters(&par[0]);

  n_xfunc = new TF1("n_xfunc","gaus",-4,4);
  n_xfunc->SetParameters(&par[3]);
  
  //Cut between p/n spot for the y-direction fitting
  double x_cut = -1.0;
  int x_cut_bin = hdxdy->GetYaxis()->FindBin(x_cut);

  TH1D *hdy_p = hdxdy->ProjectionX("dy_p",0,x_cut_bin);
  TH1D *hdy_n = hdxdy->ProjectionX("dy_n",x_cut_bin,-1);
  
  //This direction is just a gaus fit
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

  //Right now for He3 use 1 sigma and 2 sigma for H2
  double nsigma = 1;

  if(cfg == "H2") nsigma = 2;

  //Draw ellipses to show the fit results
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
  pt->AddText(Form("%g < W^{2} < %g",W2min,W2max));
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
  cout<<"neutron yield: "<<abs(n_xfunc->Integral(-4,4)/hdx->GetBinWidth(0))<<endl;
  cout<<"\n\n";

  //Save the fit result for use later
  *fit_n = n_xfunc;
 
}



void get_n_yield(TH2D *hdxdy,TString config, TF1 **fit_bg,TF1 **fit_n){

  //First project to x
  TH1D *hdx = hdxdy->ProjectionY();

  //These parameters are determined by looking at the peaks on the plots by eye
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
    
    bg_low = -5.0;
    bg_high = 3.0;
  }

  //total function = proton gaus + neutron gaus + 4th order bkgd
  double par[11];
  TF1 *p_xfunc = new TF1("p_xfunc","gaus",p_low,p_high);
  TF1 *n_xfunc = new TF1("n_xfunc","gaus",n_low,n_high);  
  TF1 *bg_xfunc = new TF1("bg_xfunc","pol4",bg_low,bg_high);
  TF1 *total_xfunc = new TF1("total_xfunc","gaus(0) + gaus(3) + pol4(6)",-5,3);
  
  //Do fit but do not plot the results
  hdx->Fit(p_xfunc,"qNR");  
  hdx->Fit(n_xfunc,"qNR+");  
  hdx->Fit(bg_xfunc,"qNR+"); 

  //Put the fit parameters into the array
  p_xfunc->GetParameters(&par[0]);
  n_xfunc->GetParameters(&par[3]);
  bg_xfunc->GetParameters(&par[6]);

  //Set the parameters in the total function using the results above
  total_xfunc->SetParameters(par);
  hdx->Fit(total_xfunc,"qNR+"); 
  
  //Get the fit results
  total_xfunc->GetParameters(&par[0]);

  //For plotting purposes set the p/n function parameters from the total fit
  p_xfunc = new TF1("p_xfunc","gaus",-4,4);
  p_xfunc->SetParameters(&par[0]);

  n_xfunc = new TF1("n_xfunc","gaus",-4,4);
  n_xfunc->SetParameters(&par[3]);

  bg_xfunc = new TF1("bg_xfunc","pol4",-5,3);
  bg_xfunc->SetParameters(&par[6]);

  //Save the fit result for use later
  *fit_bg = bg_xfunc;
  *fit_n = n_xfunc;
 
}

void Elastic_hist(TString cfg = "GEN2"){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  //Read the He3 run and H2 run files
  TFile *H2_file = new TFile("outfiles/QE_test_" + cfg + "_sbs100p_nucleon_p_model2_data.root","read");
  TFile *He3_file = new TFile("outfiles/QE_test_" + cfg + "_sbs100p_nucleon_np_model2_data.root","read");
  
  JSONManager *jmgr_H2 = new JSONManager("../config/" + cfg + "_H2.cfg");
  vector<double> dx_p_H2; jmgr_H2->GetVectorFromKey<double>("dx_p", dx_p_H2);
  vector<double> dy_p_H2; jmgr_H2->GetVectorFromKey<double>("dy_p", dy_p_H2);
  double Nsigma_cut_dx_p_H2 = jmgr_H2->GetValueFromKey<double>("Nsigma_cut_dx_p");
  double Nsigma_cut_dy_p_H2 = jmgr_H2->GetValueFromKey<double>("Nsigma_cut_dy_p");
  

  JSONManager *jmgr_He3 = new JSONManager("../config/" + cfg + "_He3.cfg");
  vector<double> dx_p_He3; jmgr_He3->GetVectorFromKey<double>("dx_p", dx_p_He3);
  vector<double> dy_p_He3; jmgr_He3->GetVectorFromKey<double>("dy_p", dy_p_He3);
  double Nsigma_cut_dx_p_He3 = jmgr_He3->GetValueFromKey<double>("Nsigma_cut_dx_p");
  double Nsigma_cut_dy_p_He3 = jmgr_He3->GetValueFromKey<double>("Nsigma_cut_dy_p");
  vector<double> dx_n_He3; jmgr_He3->GetVectorFromKey<double>("dx_n", dx_n_He3);
  vector<double> dy_n_He3; jmgr_He3->GetVectorFromKey<double>("dy_n", dy_n_He3);
  double Nsigma_cut_dx_n_He3 = jmgr_He3->GetValueFromKey<double>("Nsigma_cut_dx_n");
  double Nsigma_cut_dy_n_He3 = jmgr_He3->GetValueFromKey<double>("Nsigma_cut_dy_n");
  


   ////////////////////// ~~~~~~~~First analyze the H2 file~~~~~~~~  //////////////////////
  TTree *T = (TTree*)H2_file->Get("Tout");

  //Set the histograms that will be filled
  TH2D *hdxdy_nocut_H2 = new TH2D("hdxdy_nocut_H2","",150,-2,2,150,-6,6);
  TH2D *hdxdy_W2cut_H2 = new TH2D("hdxdy_W2cut_H2","",150,-2,2,150,-6,6);
  TH2D *hdxdy_coin_H2 = new TH2D("hdxdy_coin_H2","",150,-2,2,150,-6,6);
  TH1D *hW2_all_H2 = new TH1D("hW2_all_H2","",200,-0.5,3);
  TH1D *hW2_cut_H2 = new TH1D("hW2_cut_H2","",200,-0.5,3);

  //Fill histograms directly from the tree using Draw funcitons
  T->Draw("W2>>hW2_cut_H2","pCut");
  T->Draw("W2>>hW2_all_H2");
  T->Draw("dx:dy>>hdxdy_nocut_H2");
  T->Draw("dx:dy>>hdxdy_W2cut_H2",Form("W2 > %g && W2 < %g",W2min,W2max));
  //T->Draw("dx:dy>>hdxdy_coin_H2","coinCut && WCut");
  T->Draw("dx:dy>>hdxdy_coin_H2",Form("coinCut &&  W2 > %g && W2 < %g",W2min,W2max));


   ////////////////////// ~~~~~~~~Delete tree and now switch to He3 file~~~~~~~~  //////////////////////
  T->Delete();

  T = (TTree*)He3_file->Get("Tout");
  
  //Set histograms to be filled
  TH2D *hdxdy_nocut_He3 = new TH2D("hdxdy_nocut_He3","",150,-2,2,150,-6,6);
  TH2D *hdxdy_Wcut_He3 = new TH2D("hdxdy_Wcut_He3","",150,-2,2,150,-6,6);
  TH2D *hdxdy_coin_He3 = new TH2D("hdxdy_coin_He3","",150,-2,2,150,-6,6);
  TH2D *hdxdy_coinhodo_He3 = new TH2D("hdxdy_coinhodo_He3","",150,-2,2,150,-6,6);
  TH1D *hW2_all_He3 = new TH1D("hW2_all_He3","",200,-0.5,3);
  TH1D *hW2_cut_He3 = new TH1D("hW2_cut_He3","",200,-0.5,3);

  //Fill histograms same as above for H2
  T->Draw("W2>>hW2_cut_He3","(pCut || nCut) && coinCut");
  T->Draw("W2>>hW2_all_He3");
  T->Draw("dx:dy>>hdxdy_nocut_He3");
  T->Draw("dx:dy>>hdxdy_Wcut_He3",Form("W2 > %g && W2 < %g",W2min,W2max));
  //T->Draw("dx:dy>>hdxdy_coin_He3","coinCut && WCut");
  T->Draw("dx:dy>>hdxdy_coinhodo_He3","WCut && (hodo_time[0] - hcal_time) > 4 && (hodo_time[0] - hcal_time) < 17");
  T->Draw("dx:dy>>hdxdy_coin_He3",Form("!coinCut && W2 > %g && W2 < %g",W2min,W2max));
  
  
  
 
  ////////////////////// ~~~~~~~Plot the results~~~~~~~~  //////////////////////

  //Get some histograms that are projections of other histograms
  TH1D *hdx_H2 = hdxdy_coin_H2->ProjectionY("");
  hdx_H2->Scale(1/hdx_H2->GetEntries());
  
  TH1D *hdx_He3 = hdxdy_coin_He3->ProjectionY("hdx_He3");
  hdx_He3->SetLineColor(kRed);
  hdx_He3->Scale(1/hdx_He3->GetEntries());

  TH1D *hdx_nocut_He3 = hdxdy_nocut_He3->ProjectionY("hdx_nocut_He3");
  TH1D *hdx_Wcut_He3 = hdxdy_Wcut_He3->ProjectionY("hdx_Wcut_He3");
  TH1D *hdx_coin_He3 = hdxdy_coin_He3->ProjectionY("hdx_coin_He3");
  hdx_coin_He3->SetLineColor(kRed);

  
  TF1 *fit_H2;
  TCanvas *c1 = new TCanvas("c1","",800,1000);  
  //get_np_spots(c1, "H2", hdxdy_coin_H2,cfg,&fit_H2);
  hdxdy_coin_H2->Draw("colz");
  hdxdy_coin_H2->SetTitle("H2 HCal Elastics;#Deltay (m);#Deltax (m)");
  
  TPaveText *pt10 = new TPaveText(.55,.8,.88,.88,"ndc");
  pt10->AddText("Cuts on good tracks");
  pt10->AddText(Form("%g < W^{2} < %g",W2min,W2max));
  pt10->SetFillColor(0);
  pt10->Draw("same");

  //Draw ellipses
  TEllipse *p_spot_H2 = new TEllipse(dy_p_H2[0],dx_p_H2[0],Nsigma_cut_dy_p_H2*dy_p_H2[1],Nsigma_cut_dx_p_H2*dx_p_H2[1]);
  p_spot_H2->SetFillStyle(0);
  p_spot_H2->SetLineWidth(4);
  p_spot_H2->SetLineColor(kRed);
  p_spot_H2->Draw("same");
  

  TF1 *fit_He3;
  TF1 *fit_bg;
  get_n_yield(hdxdy_coin_He3,cfg,&fit_bg,&fit_He3);

  int n_yield = abs(fit_He3->Integral(-4,4)/hdx_He3->GetBinWidth(0));
  cout<<"neutron yield: "<<n_yield<<endl;


  TCanvas *c2 = new TCanvas("c2","",800,1000);  
  //get_np_spots(c2, "He3", hdxdy_coin_He3,cfg,&fit_He3);
  hdxdy_coin_He3->Draw("colz");
  hdxdy_coin_He3->SetTitle("He3 HCal Elastics;#Deltay (m);#Deltax (m)");

  TPaveText *pt11 = new TPaveText(.55,.8,.88,.88,"ndc");
  pt11->AddText("Cuts on good tracks");
  pt11->AddText(Form("%g < W^{2} < %g",W2min,W2max));
  pt11->SetFillColor(0);
  pt11->Draw("same");
  
  //Draw ellipses
  TEllipse *p_spot_He3 = new TEllipse(dy_p_He3[0],dx_p_He3[0],Nsigma_cut_dy_p_He3*dy_p_He3[1],Nsigma_cut_dx_p_He3*dx_p_He3[1]);
  TEllipse *n_spot_He3 = new TEllipse(dy_n_He3[0],dx_n_He3[0],Nsigma_cut_dy_n_He3*dy_n_He3[1],Nsigma_cut_dx_n_He3*dx_n_He3[1]);


  p_spot_He3->SetFillStyle(0);
  n_spot_He3->SetFillStyle(0);

  p_spot_He3->SetLineWidth(4);
  n_spot_He3->SetLineWidth(4);

  p_spot_He3->SetLineColor(kRed);
  n_spot_He3->SetLineColor(kRed);
  
  p_spot_He3->Draw("same");
  n_spot_He3->Draw("same");

  TCanvas *c3 = new TCanvas("c3","",1000,800);  
  hdx_Wcut_He3->Draw();
  hdx_coin_He3->SetLineColor(kRed);
  hdx_coin_He3->Draw("same");
  hdx_Wcut_He3->SetTitle("HCal He3 Data;#Deltax;Entries");
  //fit_He3->Draw("same");
  //fit_bg->SetLineStyle(2);
  //fit_bg->Draw("same");
  //hdx_coinhodo_He3->Draw("same hist");

  TLegend *legend = new TLegend(0.43,0.75,0.89,0.89);
  legend->AddEntry("hdx_Wcut_He3",Form("Good Tracks & %g < W^{2} < %g",W2min,W2max),"l");
  legend->AddEntry("hdx_coin_He3","Anti Coincidence Cut","l");
  //legend->AddEntry("hdx_coin_He3",Form("#splitline{Good Tracks & %g < W^{2} < %g}{& Good Coincidence}",W2min,W2max),"l");
  //legend->AddEntry("bg_xfunc","Background Fit","l");
  legend->SetLineColor(0);
  legend->Draw("same");

  TPaveText *pt = new TPaveText(.5,.6,.89,.7,"ndc");
  pt->AddText(Form("%i Neutrons in Gaussian",n_yield));
  pt->SetFillColor(0);
  //pt->Draw("same");
  
  TCanvas *c4 = new TCanvas("c4","",1000,800);  
  hdx_H2->SetTitle("HCal p/n Spots;#Deltax(m);Normalized Entries");  
  hdx_H2->Draw("hist");
  hdx_He3->Draw("same hist");

  TLegend *legend2 = new TLegend(0.55,0.7,0.89,0.89);
  legend2->AddEntry("hdx_H2","H2 Data","l");
  legend2->AddEntry("hdx_He3","He3 Data","l");
  legend2->SetLineColor(0);
  legend2->Draw("same");

  TCanvas *c5 = new TCanvas("c5","",1000,800);  
  hW2_all_H2->SetTitle("H2 Elastic Data;W^{2} (GeV^{2});Entries");  
  hW2_all_H2->Draw();
  hW2_cut_H2->SetLineColor(kRed);
  //hW2_cut_H2->Scale(5);   
  hW2_cut_H2->Draw("same hist");

  TLegend *legend3 = new TLegend(0.55,0.75,0.89,0.89);
  legend3->AddEntry("hW2_all_H2","All Events","l");
  //legend3->AddEntry("hW2_cut_H2","Inside p spot, coincidence, scaled by x5","l");
  legend3->AddEntry("hW2_cut_H2","Inside p spot","l");
  legend3->SetLineColor(0);
  legend3->Draw("same");

  TPaveText *pt2 = new TPaveText(.11,.6,.38,.73,"ndc");
  pt2->AddText(Form("Proton spot #sigma_{x} = %g m",dx_p_H2[1]));
  pt2->AddText(Form("Proton spot #sigma_{y} = %g m",dy_p_H2[1]));
  pt2->SetFillColor(0);
  pt2->Draw("same");

  TCanvas *c6 = new TCanvas("c6","",1000,800);  
  hW2_all_He3->SetTitle("He3 Elastic Data;W^{2} (GeV^{2});Entries");  
  hW2_all_He3->Draw();
  hW2_cut_He3->SetLineColor(kRed);
  hW2_cut_He3->Scale(5); 
  hW2_cut_He3->Draw("hist same");

  TLegend *legend4 = new TLegend(0.11,0.75,0.70,0.89);
  legend4->AddEntry("hW2_all_He3","All Events","l");
  legend4->AddEntry("hW2_cut_He3","Inside p/n spot, scaled by x5","l");
  //legend4->AddEntry("hW2_cut_He3","Inside p/n spot","l");
  legend4->SetLineColor(0);
  legend4->Draw("same");
  
  TPaveText *pt3 = new TPaveText(.11,.6,.38,.73,"ndc");
  pt3->AddText(Form("Proton spot #sigma_{x} = %g m",dx_p_He3[1]));
  pt3->AddText(Form("Proton spot #sigma_{y} = %g m",dy_p_He3[1]));
  pt3->AddText(Form("Neutron spot #sigma_{x} = %g m",dx_n_He3[1]));
  pt3->AddText(Form("Neutron spot #sigma_{y} = %g m",dy_n_He3[1]));
  pt3->SetFillColor(0);
  pt3->Draw("same");

  TString outputfile = "../plots/" + cfg + "_elastics.pdf";
  
  c1->Print(outputfile + "(");
  c2->Print(outputfile);
  c3->Print(outputfile);
  c4->Print(outputfile);
  c5->Print(outputfile);
  c6->Print(outputfile + ")");
  
}
