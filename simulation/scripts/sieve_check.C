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

void sieve_check(int Q2 = 368)
{
  TString Rootfiles = "/lustre19/expphy/volatile/halla/sbs/jeffas/GMN_root/simulation/";
  
  //TFile *f = new TFile("gmn_4.5GeV2_elas.root","update");
  //TFile *f = new TFile(Rootfiles + Form("gen%i_elastic.root",Q2),"update");
  TChain *T = new TChain("T");

  //TTree *T = (TTree*)f->Get("T");
  //T->Add(Rootfiles + Form("gen_%i_elastic_*job*.root",Q2));
  T->Add(Rootfiles + Form("gmn_%i_elastic_sieve_*job*.root",Q2));

  gmn_tree *G = new gmn_tree(T); //defining a pointer to the skeleton class

  T->Draw("Earm.BBGEM.Track.X:Earm.BBGEM.Track.Y","Earm.BBGEM.Track.ntracks == 1 && Earm.BBGEM.Track.NumHits == 5");

  Long64_t nevent = T->GetEntries();
 
  /*
  for (Long64_t i=5;i<nevent;i++) {
    T->GetEntry(i);
          

  }
  */

}
