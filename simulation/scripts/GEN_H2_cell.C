//This macro estimates the single arm BB trigger rate for GMN using full trigger logic
//Note: PS rows go from 0-25 whereas SH rows fo from 0-26

#include "/work/halla/sbs/jeffas/SBS_GEANT/g4sbs/install/root_macros/gmn_tree.C"
//#include "/home/deep/G4SBS/build/root_macros/gmn_tree.C"
#include "G4SBSRunData.hh"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TLatex.h"
#include <iostream>

void GEN_H2_cell(int Q2 = 17){

  TString Rootfiles = "/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/simulation/";
  
  const int nfiles = 40;

  TChain *T = new TChain("T");

  TFile *ftemp;
  G4SBSRunData *rdtemp;

  double Luminosity;
  int Nevents_total = 0;

  for(int ifile=1; ifile < nfiles; ifile++){
    
    TString filename = Rootfiles + Form("gen_%i_H2_elastic_50000events_job%i.root",Q2,ifile);
    
    ftemp = new TFile(filename,"READ");
    if( !ftemp->IsZombie() ){
      ftemp->GetObject("run_data",rdtemp);
      if( rdtemp != NULL ){
	int Ngen = rdtemp->fNtries;
	Nevents_total += Ngen;
	Luminosity= rdtemp->fLuminosity;

	T->Add(filename);	

      }
    }
  }

  gmn_tree *G = new gmn_tree(T);

  
  int nevent = 5;
  TRandom3 num(0);
 
  TH1D *hvz = new TH1D("hvz","H2 10 atm Target;Vertex V (m);Rate/bin (Hz)",100,-0.4,0.4);
  
  
  while(T->GetEntry(nevent++)){

    double weight = G->ev_sigma * G->ev_solang * Luminosity / Nevents_total;

    hvz->Fill(G->ev_vz, G->ev_rate * weight);

  }
  
  
  hvz->Draw("hist");
  
} //end of function
