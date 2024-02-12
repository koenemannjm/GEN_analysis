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


void Data_sim_compare(TString input = "GEN2",TString tgt = "He3"){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  Analysis::bg_option = "pol4";
  bool use_dy_cut = false;

  TString cfg = input;
  TString cfg_sim = input;
  if(input == "GEN4all"){
    cfg = "GEN4";
    cfg_sim = "GEN4";
  }
  else if(input == "GEN4b") cfg_sim = "GEN4";


  TString jmgr_file = "../../config/" + cfg + "_" + tgt + ".cfg";
  JSONManager *jmgr = new JSONManager(jmgr_file);

   // elastic cut limits
  double W2min = jmgr->GetValueFromKey<double>("W2min");
  double W2max = jmgr->GetValueFromKey<double>("W2max");
  //W2min = 0.48;
  //W2max = 1.28;

  double dymin = jmgr->GetValueFromKey<double>("dymin");
  double dymax = jmgr->GetValueFromKey<double>("dymax");
  //dymin = -1.2;
  //dymax = 1.2;

  vector<double> coin_time_cut; jmgr->GetVectorFromKey<double>("coin_time", coin_time_cut);

  TString reaction = "np";
  if(tgt == "H2") reaction = "p";
  TString model = "2";
  if(tgt == "H2") model = "1";

  //Read the He3 run and H2 data files
  TFile *data_file = new TFile("../outfiles/QE_data_" + cfg + "_sbs100p_nucleon_" + reaction + "_model" + model + ".root","read");

  //Read the He3 run and H2 simulation files
  TFile *sim_file = new TFile("../outfiles/QE_sim_" + cfg_sim + "_sbs100p_nucleon_" + reaction + "_model" + model + ".root","read");

  //Get the TTrees
  TTree *T_sim = (TTree*)sim_file->Get("Tout");
  
  TChain *T_data = new TChain("Tout");
  if(input == "GEN4all"){
    T_data->Add("../outfiles/QE_data_GEN4_sbs100p_nucleon_np_model2.root");
    T_data->Add("../outfiles/QE_data_GEN4b_sbs100p_nucleon_np_model2.root");
  }
  else {
    T_data->Add("../outfiles/QE_data_" + cfg + "_sbs100p_nucleon_np_model2.root");
  }

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
  Analysis::hdx_sim_p = new TH1F("hdx_sim_p","",nbins,xmin,xmax);
  Analysis::hdx_sim_n = new TH1F("hdx_sim_n","",nbins,xmin,xmax);
  Analysis::hdx_bg_data = new TH1F("hdx_bg_data","",nbins,xmin,xmax);

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

  T_data->Draw("dx>>hdx_data",DataCut);
  T_data->Draw("dx>>hdx_bg_data",coinCut && W2Cut && !dyCut);
  T_sim->Draw("dx>>hdx_sim_p",CutSimP);
  T_sim->Draw("dx>>hdx_sim_n",CutSimN);
  
  double scale = hdx_data->Integral();

  hdx_data->Scale(1.0/hdx_data->Integral());
  Analysis::hdx_sim_p->Scale(1.0/Analysis::hdx_sim_p->Integral());
  Analysis::hdx_sim_n->Scale(1.0/Analysis::hdx_sim_n->Integral());
  Analysis::hdx_bg_data->Scale(1.0/Analysis::hdx_bg_data->Integral());

  Analysis::He3_sim_fit(hdx_data);

  //Copy all the result histograms
  TH1F *hdx_sim_p = new TH1F(*Analysis::hdx_sim_p);
  TH1F *hdx_sim_n = new TH1F(*Analysis::hdx_sim_n);
  TH1F *hdx_bg = new TH1F(*Analysis::hdx_bg);
  TH1F *hdx_total_fit = new TH1F(*Analysis::hdx_total_fit);


  gStyle->SetOptFit(0);
  
  hdx_data->SetTitle("Data/Simulation Comparisons " + cfg + ";#Deltax (m);Entries");

  hdx_data->Scale(scale);
  hdx_total_fit->Scale(scale);
  hdx_sim_p->Scale(scale);
  hdx_sim_n->Scale(scale);
  hdx_bg->Scale(scale);

  hdx_data->SetMarkerStyle(kFullCircle);
  hdx_total_fit->SetFillColorAlpha(30,0.5);
  hdx_sim_p->SetFillColorAlpha(kRed,0.3);
  hdx_sim_n->SetFillColorAlpha(kBlue,0.3);
  hdx_bg->SetFillColorAlpha(kMagenta,0.3);

  hdx_total_fit->SetLineStyle(7);
  hdx_sim_p->SetLineStyle(7);
  hdx_sim_n->SetLineStyle(7);
  hdx_bg->SetLineStyle(7);
  
  hdx_total_fit->SetLineColor(30);
  hdx_sim_p->SetLineColor(kRed);
  hdx_sim_n->SetLineColor(kBlue);
  hdx_bg->SetLineColor(kMagenta);
  
  TCanvas *c = new TCanvas("c","",800,600);
  hdx_data->Draw();
  hdx_total_fit->Draw("same hist");
  hdx_sim_p->Draw("same hist");
  hdx_sim_n->Draw("same hist");
  hdx_bg->Draw("same hist");

  TString bg_fit_type;

  if(Analysis::bg_option == "pol4")
    bg_fit_type = "BG 4th od fit";
  else if(Analysis::bg_option == "from data")
    bg_fit_type = "BG from data";
  else if(Analysis::bg_option == "pol3")
    bg_fit_type = "BG 3rd od fit";
  else if(Analysis::bg_option == "gaus")
    bg_fit_type = "BG Gaussian";

  TLegend *legend = new TLegend(0.65,0.72,0.89,0.89);
  legend->AddEntry("hdx_data","Data","p");
  legend->AddEntry("hdx_total_fit","MC Fit","lf");
  legend->AddEntry("hdx_sim_p","MC p","lf");
  legend->AddEntry("hdx_sim_n","MC n","lf");
  legend->AddEntry("hdx_bg",bg_fit_type,"lf");
  legend->SetLineColor(0);
  legend->Draw("same");

  TPaveText *pt = new TPaveText(.65,.50,.88,.70,"ndc");
  pt->AddText("Cuts on good tracks");
  pt->AddText("Coincidence Cuts");
  pt->AddText(Form("%g < W^{2} < %g",W2min,W2max));
  if(use_dy_cut) pt->AddText(Form("%g < #Deltay < %g",dymin,dymax));
  pt->AddText(Form("%i Neutrons",(int)hdx_sim_n->GetSumOfWeights()))->SetTextColor(kBlue);
  pt->SetFillColor(0);
  pt->Draw("same");


  TString output = "Data_sim_"+input+".pdf";
  
  c->SaveAs("../../plots/" + output);


  TCanvas *c2 = new TCanvas("c2","",800,600);
  Analysis::hdx_bg_data->Draw("hist");
}
