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

void EArm_test(int Q2 = 29){

  TString Rootfiles = "/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/simulation/";
  
  int nfiles = 3;
  double right_bound = 4.0;
  if(Q2 == 66) right_bound = 5.0;
  if(Q2 == 97) right_bound = 5.0;

  TChain *T = new TChain("T");

  TFile *ftemp;
  G4SBSRunData *rdtemp;

  double Luminosity;
  int Nevents_total = 0;

  for(int ifile=1; ifile < nfiles; ifile++){
    
    TString filename = Rootfiles + Form("gen_%i_elastic_100000events_job%i.root",Q2,ifile);
    
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
  

  TH1D *maxsig = new TH1D("maxsig",";PS + SH energy;",500,1.0,3.5);
  
  while(T->GetEntry(nevent++)){
    //  for (Long64_t i=10; i<nevent; i++) {
    //T->GetEntry(i);
        
    if((G->Earm_BBPSTF1_det_esum + G->Earm_BBSHTF1_det_esum)>1.0 ){  


      double weight = G->ev_sigma * G->ev_solang * Luminosity / Nevents_total;
      maxsig->Fill(G->Earm_BBPSTF1_det_esum + G->Earm_BBSHTF1_det_esum, G->ev_rate * weight);
      
    
    } //for the if statement on top 
     
  } //end of event loop

  maxsig->Draw();
  cout<<"Luminosity "<<Luminosity<<endl;
  cout<<"Total Events "<<Nevents_total<<endl;
  cout<<"Rate "<<maxsig->Integral()<<endl;
  
} //end of function
