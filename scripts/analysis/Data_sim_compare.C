//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified February 2, 2024
//
//
//   The purpose of this script is to compare real data and
//   simulated data for the same kinematic point
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

#include "../../include/gen-ana.h"



void Data_sim_compare(TString cfg = "GEN2"){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  distribution_fits *dists = new distribution_fits();

  if(cfg == "GEN2") dists->SetBgShapeOption("pol2");
  else dists->SetBgShapeOption("from data");

  bool use_dy_cut = false;

  TString jmgr_file = "../../config/" + cfg + "_He3.cfg";

  Utilities::KinConf kin_info = Utilities::LoadKinConfig(jmgr_file,1);

  analyzed_tree *T_data = Utilities::LoadAnalyzedRootFiles(kin_info,1,0);
  analyzed_tree *T_sim = Utilities::LoadAnalyzedRootFiles(kin_info,0,0);

  // elastic cut limits
  double W2min = kin_info.W2min;
  double W2max = kin_info.W2max;
  double dymin = kin_info.dymin;
  double dymax = kin_info.dymax;

  vector<double> coin_time_cut = kin_info.coin_time_cut; 
  
  
  /////Set the histograms
  int nbins = 100;
  double xmin = -4;
  double xmax = 2.5;

  if(cfg == "GEN2"){
    xmin = -6;
    xmax = 3;
  }
  
  //dx
  TH1F *hdx_data = new TH1F("hdx_data","",nbins,xmin,xmax);
  TH1F *hdx_sim_p = new TH1F("hdx_sim_p","",nbins,xmin,xmax);
  TH1F *hdx_sim_n = new TH1F("hdx_sim_n","",nbins,xmin,xmax);
  TH1F *hdx_bg_data = new TH1F("hdx_bg_data","",nbins,xmin,xmax);

  TCut W2Cut = Form("W2 > %g && W2 < %g",W2min,W2max);
  TCut dyCut = Form("dy > %g && dy < %g",dymin,dymax);
  //TCut coinCut = Form("abs(coin_time - %g) < %g*2",coin_time_cut[0],coin_time_cut[1]);
  TCut coinCut = "coinCut";
  
  TCut DataCut = W2Cut && coinCut;
  TCut CutSimP = Form("(W2 > %g && W2 < %g && fnucl == 1) * weight",W2min,W2max);
  TCut CutSimN = Form("(W2 > %g && W2 < %g && fnucl == 0) * weight",W2min,W2max);
  
  if(use_dy_cut){
    DataCut = DataCut && dyCut;
    CutSimP = Form("(W2 > %g && W2 < %g && dy > %g && dy < %g && fnucl == 1) * weight",W2min,W2max,dymin,dymax);
    CutSimN = Form("(W2 > %g && W2 < %g && dy > %g && dy < %g && fnucl == 0) * weight",W2min,W2max,dymin,dymax);
  }

  T_data->fChain->Draw("dx>>hdx_data",DataCut);
  T_data->fChain->Draw("dx>>hdx_bg_data",coinCut && W2Cut && !dyCut);
  T_sim->fChain->Draw("dx>>hdx_sim_p",CutSimP);
  T_sim->fChain->Draw("dx>>hdx_sim_n",CutSimN);
  
  dists->SetDataShape(hdx_data);
  dists->SetPShape(hdx_sim_p);
  dists->SetNShape(hdx_sim_n);
  dists->SetBgShape(hdx_bg_data);
  
  dists->He3_fit_dists();
  
  //Copy all the result histograms
  TH1F *hdx_data_plot = dists->GetDataHist();
  TH1F *hdx_sim_p_plot = dists->GetPHist();
  TH1F *hdx_sim_n_plot = dists->GetNHist();
  TH1F *hdx_bg_plot = dists->GetBgHist();
  TH1F *hdx_total_fit_plot = dists->GetTotalHist();

  
  gStyle->SetOptFit(0);
  
  hdx_data_plot->SetTitle("Data/Simulation Comparisons " + cfg + ";#Deltax (m);Entries");

  hdx_data_plot->SetMarkerStyle(kFullCircle);
  hdx_total_fit_plot->SetFillColorAlpha(30,0.5);
  hdx_sim_p_plot->SetFillColorAlpha(kRed,0.3);
  hdx_sim_n_plot->SetFillColorAlpha(kBlue,0.3);
  hdx_bg_plot->SetFillColorAlpha(kMagenta,0.3);

  hdx_total_fit_plot->SetLineStyle(7);
  hdx_sim_p_plot->SetLineStyle(7);
  hdx_sim_n_plot->SetLineStyle(7);
  hdx_bg_plot->SetLineStyle(7);
  
  hdx_total_fit_plot->SetLineColor(30);
  hdx_sim_p_plot->SetLineColor(kRed);
  hdx_sim_n_plot->SetLineColor(kBlue);
  hdx_bg_plot->SetLineColor(kMagenta);
  
  TCanvas *c = new TCanvas("c","",800,600);
  hdx_data_plot->Draw();
  hdx_total_fit_plot->Draw("same hist");
  hdx_sim_p_plot->Draw("same hist");
  hdx_sim_n_plot->Draw("same hist");
  hdx_bg_plot->Draw("same hist");


  TLegend *legend = new TLegend(0.65,0.72,0.89,0.89);
  legend->AddEntry("hdx_data","Data","p");
  legend->AddEntry("hdx_total_fit","MC Fit","lf");
  legend->AddEntry("hdx_sim_p","MC p","lf");
  legend->AddEntry("hdx_sim_n","MC n","lf");
  legend->AddEntry("hdx_bg","Background","lf");
  legend->SetLineColor(0);
  legend->Draw("same");

  TPaveText *pt = new TPaveText(.65,.50,.88,.70,"ndc");
  pt->AddText("Cuts on good tracks");
  pt->AddText("Coincidence Cuts");
  pt->AddText(Form("%g < W^{2} < %g",W2min,W2max));
  if(use_dy_cut) pt->AddText(Form("%g < #Deltay < %g",dymin,dymax));
  pt->AddText(Form("%i Neutrons",(int)hdx_sim_n_plot->GetSumOfWeights()))->SetTextColor(kBlue);
  pt->SetFillColor(0);
  pt->Draw("same");


  TString output = "Data_sim_"+cfg+".pdf";
  
  //c->SaveAs("../../plots/" + output);
  
}
