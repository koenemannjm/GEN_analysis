// This macro estimates the single arm (HCAL) trigger rate in GMN experiment using minimum bias simulation data (sets the threshold @ choosen efficiency of elastic data)
//#include "/home/uconn421/G4SBS/install/root_macros/gmn_tree.C"
#include "/work/halla/sbs/jeffas/SBS_GEANT/g4sbs/install/root_macros/gmn_tree.C"
//#include "/w/halla-scshelf2102/sbs/jeffas/SBS_GEANT/g4sbs/install/include/G4SBSRunData.hh"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom.h"
#include "TGraph.h"
#include <iostream>
#include "TSpline.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TLatex.h"
 

void hc_tgrt_beam_gmn(int Q2 = 368)
{
  TString Rootfiles = "/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/simulation/";
  
  TChain *C = new TChain("T");
  C->Add(Rootfiles + Form("gen_%i_beam_bkgd_*job*.root",Q2));
  

  double Q2_double = Q2*1.0/10;
  double_t threshold = 0.061; //GeV

  if(Q2 == 29) threshold = 0.023;
  if(Q2 == 66) threshold = 0.104;
  if(Q2 == 97) threshold = 0.218;

  G4SBSRunData *rd;
  long ngen = 0;
  int nfiles = 0;
  TObjArray *FileList = C->GetListOfFiles();
  TIter next(FileList);  
  TChainElement *chEl = 0;
  set<TString> bad_file_list;

  while( (chEl=(TChainElement*)next() )){
    TFile newfile(chEl->GetTitle());
    newfile.GetObject("run_data",rd);
    if( rd ){
      ngen += rd->fNtries;
      nfiles++;
    } else {
      bad_file_list.insert( chEl->GetTitle());
    }
    //cout << chEl->GetTitle() << endl;
  }

  double Ibeam = 60e-6; //A
  double weight = Ibeam/double(ngen)/1.602e-19;

  
  cout << "Total number of generated events = " << ngen << endl;


  gmn_tree *G = new gmn_tree(C); //defining a pointer to the skeleton class

  TCanvas *c1 = new TCanvas("c1","c1",1200,900);
 
  Double_t Harm_HCalScint_det_esum_sm;
  Long64_t nevent = C->GetEntries();
  TRandom3 num(0);

  TH1D *hesum_true = new TH1D("hesum_true","",500,0.0,0.31); //true histogram
  TH1D *hesum_sm = new TH1D("hesum_sm","",500,0.0,0.31); //smeared histogram


  for (Long64_t i=0;i<nevent;i++) {
    C->GetEntry(i);

    if(i%100000 == 0) cout<<i*1.0/nevent*100<<"%"<<endl;
    //Let's write an if statement containing all the necessary cuts for threshold analysis
    if(G->Harm_HCalScint_det_esum>0.00){
      double_t edep1 = G->Harm_HCalScint_det_esum;
      double_t mean1 = 3000.0*edep1;

      Harm_HCalScint_det_esum_sm = num.Gaus(mean1,sqrt(mean1))/3000.0;

      hesum_sm->Fill(Harm_HCalScint_det_esum_sm, weight);
      hesum_true->Fill( edep1, weight );
    }
  }

  double_t maxIn = hesum_sm->Integral(hesum_sm->FindFixBin(0),hesum_sm->FindFixBin(0.31),"");

  cout << "Total Integral = " << maxIn << endl;
  Double_t error = 0.0;
  cout << "Trigger Rate (Beam) = " << hesum_sm->IntegralAndError(hesum_sm->FindFixBin(threshold),hesum_sm->FindFixBin(0.31),error,"") << " +- " << error << " Hz" << endl;

  c1->cd();

  gStyle->SetOptStat(1111111);
  gPad->SetLogy();
  hesum_sm->SetMarkerStyle(20);
  hesum_sm->SetMarkerColor(1);
  hesum_sm->SetLineColor(1);
  hesum_sm->SetTitle(Form("Q^2=%gGeV^2    Minimum Bias Simulation    Threshold @ 90pct eff.(Elastic)=%gGeV",Q2_double,threshold));
  hesum_sm->GetXaxis()->SetTitle("Energy sum (HCAL) in GeV");
  hesum_sm->GetYaxis()->SetTitle("Rate/bin in Hz");

  hesum_true->SetMarkerStyle(21);
  hesum_true->SetMarkerColor(2);
  hesum_true->SetLineColor(2);
  hesum_sm->Draw();
  hesum_true->Draw("SAME");

  TLegend *l1 = new TLegend(0.47,0.8,0.73,0.61);  //(0.22,0.7,0.51,0.89); 
  l1->AddEntry((TObject*)0, Form("Trriger Rate = %.2e #pm %.2eHz",hesum_sm->IntegralAndError(hesum_sm->FindFixBin(threshold),hesum_sm->FindFixBin(0.31),error,""),error), "");
  //l1->AddEntry((TObject*)0, "Threshold @ 90pct efficiency(Elastic) = 0.103GeV", "");
  //l1->AddEntry((TObject*)0, Form("Threshold for 90pct efficiency = %fGeV",thresh95), "");
  l1->AddEntry((TObject*)0, "Beam Current = 60uA", "");
  l1->AddEntry((TObject*)0, "Target Length = 60cm", "");
  l1->AddEntry(hesum_true,"True Energy","ep");
  l1->AddEntry(hesum_sm,"Smeared Energy","ep");
  //l1->AddEntry(fR,"Gaussian Fit","l");
  l1->Draw();


  TString output = "../plots/";

  c1->Print(output + Form("hcal_rate_%gGeV.pdf",Q2_double));

}

