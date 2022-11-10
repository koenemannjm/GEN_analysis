// This macro finds the appropriate threshold for GMN experiment using elastic kinematics for various Q^2 values. 
#include "/work/halla/sbs/jeffas/SBS_GEANT/g4sbs/install/root_macros/gmn_tree.C"
//#include "/home/uconn421/G4SBS/install/root_macros/gmn_tree.C"
#include "TChain.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "TGraph.h"
#include <iostream>
#include "TSpline.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TF1.h"

void hc_thrsh_vs_effi(int Q2 = 368)
{
  TString Rootfiles = "/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/simulation/";
  
  //TFile *f = new TFile("gmn_4.5GeV2_elas.root","update");
  //TFile *f = new TFile(Rootfiles + Form("gen%i_elastic.root",Q2),"update");
  TChain *T = new TChain("T");

  //TTree *T = (TTree*)f->Get("T");
  //T->Add(Rootfiles + Form("gen_%i_elastic_*job*.root",Q2));
  T->Add(Rootfiles + Form("gen_%i_elastic_*job*.root",Q2));

  gmn_tree *G = new gmn_tree(T); //defining a pointer to the skeleton class

  TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  TCanvas *c2 = new TCanvas("c2","c2",1200,900);
  TCanvas *c3 = new TCanvas("c3","c3",1200,900);

  Double_t Harm_HCalScint_det_esum_sm;

  Long64_t nevent = T->GetEntries();
  TRandom3 num(0);

  double right_bound = 0.58;
  if(Q2 == 97) right_bound = 0.8;

  TH1D *hesum_true = new TH1D("hesum_true","",500,0.0,right_bound); //true histogram
  TH1D *hesum_sm = new TH1D("hesum_sm","",500,0.0,right_bound); //smeared histogram
  TH1D *heres = new TH1D("heres","",500,-1.,1.); //difference between them

  for (Long64_t i=5;i<nevent;i++) {
    T->GetEntry(i);
    //Let's write an if statement containing all the necessary cuts for threshold analysis
    //if(G->Harm_HCalScint_det_esum>0 && G->Earm_BBGEM_Track_ntracks==1){
    
    if(G->ev_fnucl == 0 && G->Harm_HCalScint_det_esum > 0.02 && (*(G->Earm_BBGEM_Track_MID))[0]==0 && (*(G->Earm_BBGEM_Track_NumHits))[0]==5 && G->Earm_BBGEM_Track_ntracks==1 && (*(G->Earm_BBGEM_Track_P))[0]/(G->ev_ep)>0.99){  //used to reproduce the result got by uding det.edum

    
      double_t edep1 = G->Harm_HCalScint_det_esum;
      double_t mean1 = 3000.0*edep1;
      // Earm_BBPSTF1_det_esum_sm = num.PoissonD(mean1)/500.0; //giving some weird discontinuity error below 200.
      // Earm_BBSHTF1_det_esum_sm = num.PoissonD(mean2)/500.0;

      Harm_HCalScint_det_esum_sm = num.Gaus(mean1,sqrt(mean1))/3000.0;

      //filling new branches
      //b_Earm_BBPSTF1_det_esum_sm->Fill();
      //b_Earm_BBSHTF1_det_esum_sm->Fill();
      
      hesum_sm->Fill(Harm_HCalScint_det_esum_sm,G->ev_rate);
      heres->Fill( Harm_HCalScint_det_esum_sm/edep1-1.0, G->ev_rate );
      hesum_true->Fill( edep1, G->ev_rate );
    }
    //cout << "Earm.BBSHTF1.det.esum=" << edep2 << "\t" << "Earm.BBSHTF1.det.esum.sm=" << Earm_BBSHTF1_det_esum_sm << endl;
  }
  hesum_sm->Fit("gaus","","",0.01,0.4);
  TF1 *fR = hesum_sm->GetFunction("gaus");
  fR->SetLineColor(4);
  cout << "Mean=" << fR->GetParameter(1) << endl;
  
  //updating the tree with new branches
  //T->Print(); 
  //T->Write(); 
  //delete f;   

  int nbins = hesum_sm->GetNbinsX();

  double threshold[nbins];
  double efficiency[nbins];
  double_t maxIn = hesum_sm->Integral(hesum_sm->FindFixBin(0),hesum_sm->FindFixBin(0.58),"");
  double thresh95 = 0.58;
  
  for(int i=1; i<=nbins; i++){
    double thresh = hesum_sm->GetBinLowEdge(i);
    double int_above = hesum_sm->Integral(i,nbins);
    double eff = int_above/maxIn;
    //Need to find the threshold at 95% efficiency
    if( eff < 0.90 && thresh < thresh95 ){
      thresh95 = thresh;
    }

    threshold[i-1] = thresh;
    efficiency[i-1] = eff;
  }

  c1->cd(); //cd() -> change to new canvas

  double Q2_double = Q2*1.0/10;

  TGraph *effic_vs_thresh_fine = new TGraph(nbins, threshold, efficiency);
  gStyle->SetOptStat(1111111);
  effic_vs_thresh_fine->SetTitle(Form("Efficiency vs. Threshold for Q^2 = %gGeV^2",Q2_double));
  effic_vs_thresh_fine->GetXaxis()->SetTitle("Threshold in GeV");
  effic_vs_thresh_fine->GetYaxis()->SetTitle("Efficiency");
  effic_vs_thresh_fine->SetMarkerStyle(21);
  effic_vs_thresh_fine->SetMarkerColor(2);
  effic_vs_thresh_fine->SetLineColor(2);
  effic_vs_thresh_fine->Draw("AP");

  cout << "threshold for 90% efficiency = " << thresh95 << endl;
  cout << "Total Integral = " << maxIn << endl;
  cout << "Trigger Rate (Elastic) = " << hesum_sm->Integral(hesum_sm->FindFixBin(thresh95),hesum_sm->FindFixBin(0.58),"") << " Hz" << endl;

  c2->cd();

  hesum_sm->SetMarkerStyle(20);
  hesum_sm->SetMarkerColor(1);
  hesum_sm->SetLineColor(1);
  hesum_sm->SetTitle(Form("Q^2 = %gGeV^2",Q2_double));
  hesum_sm->GetXaxis()->SetTitle("Energy sum (HCAL) in GeV");
  hesum_sm->GetYaxis()->SetTitle("Rate/bin in Hz");
  // hesum_sm->Fit("gaus","","",0.4,0.58);
  //hesum_sm->GetFunction("gaus")->SetLineColor(3);
  // TF1 *fR = hesum_sm->GetFunction("gaus");
  // fR->SetLineColor(3);
  //  cout << "Mean=" << fR->GetParameter(1) << endl;

  hesum_true->SetMarkerStyle(21);
  hesum_true->SetMarkerColor(2);
  hesum_true->SetLineColor(2);
  hesum_sm->Draw();
  hesum_true->Draw("SAME");

  TLegend *l1 = new TLegend(0.47,0.8,0.73,0.61);  //(0.22,0.7,0.51,0.89); 
  l1->AddEntry((TObject*)0, Form("Threshold for 90pct efficiency = %fGeV",thresh95), "");
  //l1->AddEntry((TObject*)0, "Beam Current = 30uA", "");
  //l1->AddEntry((TObject*)0, "Target Length = 15cm", "");
  l1->AddEntry(hesum_true,"True Energy","ep");
  l1->AddEntry(hesum_sm,"Smeared Energy","ep");
  l1->AddEntry(fR,"Gaussian Fit","l");
  l1->Draw();

  c3->cd();
  heres->Draw();
  
  //TSpline didn't work!
  //TSpline3 *effspline = new TSpline3("effspline",threshold,efficiency,nbins);
  //  double thresh95 = effspline.FindX(0.95);
  //cout << "threshold for 95% efficiency for elastics = " << effspline->Eval(0.95) << endl;
  //  TCanvas *c2 = new TCanvas("c2","c2",1000,750);
  //effspline->Draw("AC3");

  TString output = "../plots/";

  c1->Print(output + Form("th_ef_hc_%gGeV.pdf(",Q2_double));
  c2->Print(output + Form("th_ef_hc_%gGeV.pdf",Q2_double));
  c3->Print(output + Form("th_ef_hc_%gGeV.pdf)",Q2_double));
}

