//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified July 7, 2023
//
//
//   The purpose of this script is to compare real data and
//   simulated data for the same kinematic point
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

#include "../../include/gen-ana.h"
#include "../../dflay/src/JSONManager.cxx"


double GetYmax(TH1D *h1, TH1D *h2){

  double height1 = h1->GetBinContent(h1->GetMaximumBin());
  double height2 = h2->GetBinContent(h2->GetMaximumBin());

  double ymax = height1 > height2 ? height1 : height2;

  return ymax*1.1;
}


void Data_sim_compare(TString cfg = "GEN2"){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  double W2min = 0;
  double W2max = 1.6;
  double HCal_voffset = 0.28;

  //Read the He3 run and H2 data files
  TFile *H2_data_file = new TFile("../outfiles/QE_data_" + cfg + "_sbs100p_nucleon_p_model2.root","read");
  TFile *He3_data_file = new TFile("../outfiles/QE_data_" + cfg + "_sbs100p_nucleon_np_model2.root","read");

  //Read the He3 run and H2 simulation files
  TFile *H2_sim_file = new TFile("../outfiles/QE_sim_" + cfg + "_sbs100p_nucleon_p_model2.root","read");
  TFile *He3_sim_file = new TFile("../outfiles/QE_sim_" + cfg + "_sbs100p_nucleon_np_model2.root","read");

  //Get the TTrees
  TTree *T_H2_data = (TTree*)H2_data_file->Get("Tout");
  TTree *T_He3_data = (TTree*)He3_data_file->Get("Tout");
  TTree *T_H2_sim = (TTree*)H2_sim_file->Get("Tout");
  TTree *T_He3_sim = (TTree*)He3_sim_file->Get("Tout");

  /////Set the histograms

  //W2
  TH1D *hW2_H2_sim = new TH1D("hW2_H2_sim","",200,-0.5,3);
  TH1D *hW2_H2_data = new TH1D("hW2_H2_data","",200,-0.5,3);
  TH1D *hW2_He3_sim = new TH1D("hW2_He3_sim","",200,-0.5,3);
  TH1D *hW2_He3_data = new TH1D("hW2_He3_data","",200,-0.5,3);
  //dx
  TH1D *hdx_H2_data = new TH1D("hdx_H2_data","",200,-6,6);
  TH1D *hdx_H2_sim = new TH1D("hdx_H2_sim","",200,-6,6);
  TH1D *hdx_He3_data = new TH1D("hdx_He3_data","",200,-6,6);
  TH1D *hdx_He3_sim = new TH1D("hdx_He3_sim","",200,-6,6);
  //vertex z
  TH1D *hvz_H2_data = new TH1D("hvz_H2_data","",200,-0.3,0.3);
  TH1D *hvz_H2_sim = new TH1D("hvz_H2_sim","",200,-0.3,0.3);

  
  

  //Fill histograms directly from trees
  T_H2_data->Draw("W2>>hW2_H2_data","pCut && coinCut");
  T_H2_sim->Draw("W2>>hW2_H2_sim","weight");
  T_He3_data->Draw("W2>>hW2_He3_data","(pCut || nCut) && coinCut");
  T_He3_sim->Draw("W2>>hW2_He3_sim","weight");
  
  T_H2_data->Draw("dx>>hdx_H2_data",Form("coinCut &&  W2 > %g && W2 < %g",W2min,W2max));
  T_H2_sim->Draw("dx>>hdx_H2_sim",Form("(W2 > %g && W2 < %g) * weight",W2min,W2max));
  T_He3_data->Draw("dx>>hdx_He3_data",Form("coinCut &&  W2 > %g && W2 < %g",W2min,W2max));
  T_He3_sim->Draw("dx>>hdx_He3_sim",Form("(W2 > %g && W2 < %g) * weight",W2min,W2max));
  
  T_H2_data->Draw("vz>>hvz_H2_data",Form("W2 > %g && W2 < %g",W2min,W2max));
  T_H2_sim->Draw("vz>>hvz_H2_sim",Form("(W2 > %g && W2 < %g) * weight",W2min,W2max));

  //Draw all our canvases and histograms
  TCanvas *c1 = new TCanvas("c1","",1100,400);
  c1->Divide(2,1);
  
  TPaveText *pt2 = new TPaveText(.65,.6,.88,.75,"ndc");
  pt2->AddText("Cuts on good tracks");
  pt2->SetFillColor(0);

  c1->cd(1);
  hW2_H2_data->SetTitle("H2 W^{2} Data/Sim Comparison;W^{2} (GeV^{2});Normalized Entries");
  hW2_H2_data->Scale(1.0/hW2_H2_data->Integral());
  hW2_H2_sim->Scale(1.0/hW2_H2_sim->Integral());
  hW2_H2_data->Draw("hist");
  hW2_H2_sim->Draw("hist same");
  hW2_H2_sim->SetLineColor(kRed);
  
  TLegend *legend = new TLegend(0.6,0.75,0.89,0.89);
  legend->AddEntry("hW2_H2_data","Data","l");
  legend->AddEntry("hW2_H2_sim","Simulation","l");
  legend->SetLineColor(0);
  legend->Draw("same");
  pt2->Draw("same");
  
  c1->cd(2);
  hW2_He3_data->SetTitle("He3 W^{2} Data/Sim Comparison;W^{2} (GeV^{2});Normalized Entries");  
  hW2_He3_data->Scale(1.0/hW2_He3_data->Integral());
  hW2_He3_sim->Scale(1.0/hW2_He3_sim->Integral());
  hW2_He3_data->Draw("hist");
  hW2_He3_sim->Draw("hist same");
  hW2_He3_sim->SetLineColor(kRed);

  TLegend *legend2 = new TLegend(0.6,0.75,0.89,0.89);
  legend2->AddEntry("hW2_He3_data","Data","l");
  legend2->AddEntry("hW2_He3_sim","Simulation","l");
  legend2->SetLineColor(0);
  legend2->Draw("same");
  pt2->Draw("same");  

  TCanvas *c2 = new TCanvas("c2","",1100,400);
  c2->Divide(2,1);

  c2->cd(1);
  hdx_H2_data->SetTitle("H2 HCal #Deltax Data/Sim Comparison;#Deltax (m);Normalized Entries");  
  hdx_H2_data->Scale(1.0/hdx_H2_data->Integral());
  hdx_H2_sim->Scale(1.0/hdx_H2_sim->Integral());
  hdx_H2_data->Draw("hist");
  hdx_H2_sim->Draw("hist same");
  hdx_H2_sim->SetLineColor(kRed);

  hdx_H2_data->GetYaxis()->SetRangeUser(0,GetYmax(hdx_H2_data,hdx_H2_sim));

  TLegend *legend3 = new TLegend(0.6,0.75,0.89,0.89);
  legend3->AddEntry("hdx_H2_data","Data","l");
  legend3->AddEntry("hdx_H2_sim","Simulation","l");
  legend3->SetLineColor(0);
  legend3->Draw("same");

  TPaveText *pt = new TPaveText(.65,.6,.88,.75,"ndc");
  pt->AddText("Cuts on good tracks");
  pt->AddText(Form("%g < W^{2} < %g",W2min,W2max));
  pt->SetFillColor(0);
  pt->Draw("same");


  c2->cd(2);
  hdx_He3_data->SetTitle("He3 HCal #Deltax Data/Sim Comparison;#Deltax (m);Normalized Entries");  
  hdx_He3_data->Scale(1.0/hdx_He3_data->Integral());
  hdx_He3_sim->Scale(1.0/hdx_He3_sim->Integral());
  hdx_He3_data->Draw("hist");
  hdx_He3_sim->Draw("hist same");
  hdx_He3_sim->SetLineColor(kRed);

  hdx_He3_data->GetYaxis()->SetRangeUser(0,GetYmax(hdx_He3_data,hdx_He3_sim));

  TLegend *legend4 = new TLegend(0.6,0.75,0.89,0.89);
  legend4->AddEntry("hdx_He3_data","Data","l");
  legend4->AddEntry("hdx_He3_sim","Simulation","l");
  legend4->SetLineColor(0);
  legend4->Draw("same");
  
  pt->Draw("same");

  TCanvas *c3 = new TCanvas("c3","",800,600);
  hvz_H2_data->SetTitle("H2 HCal #Deltax Data/Sim Comparison;#Deltax (m);Normalized Entries");  
  hvz_H2_data->Scale(1.0/hvz_H2_data->Integral());
  hvz_H2_sim->Scale(1.0/hvz_H2_sim->Integral());
  hvz_H2_data->Draw("hist");
  hvz_H2_sim->Draw("hist same");
  hvz_H2_sim->SetLineColor(kRed);

  hvz_H2_data->GetYaxis()->SetRangeUser(0,GetYmax(hvz_H2_data,hvz_H2_sim));
  

}
