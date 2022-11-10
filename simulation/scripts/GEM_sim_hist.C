
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

void GEM_sim_hist(int Q2 = 29){

  
  TString Rootfiles = "/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/simulation/";
  
  double right_bound = 4.0;
  if(Q2 == 66) right_bound = 5.0;
  if(Q2 == 97) right_bound = 5.0;

  TChain *T = new TChain("T");

  TFile *ftemp;
  G4SBSRunData *rdtemp;

  double Luminosity;
  int Nevents_total = 0;
  int nfiles = 200;

  for(int ifile=1; ifile <= nfiles; ifile++){
    
    TString filename = Rootfiles + Form("gen_%i_beam_bkgd_15000000events_job%i.root",Q2,ifile);
    
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

  cout<<Nevents_total<<endl;
  gmn_tree *G = new gmn_tree(T);

  double Ibeam = 1e-6; //A
  double weight = Ibeam/double(Nevents_total)/1.602e-19;

  double BBGEM_area_cm2[5] = {40.0*150.0, 40.0*150.0, 40.0*150.0, 40.0*150.0, 60.0*200.0 };
  double BBGEM_LX[5] = {150.,150.,150.,150.,200.};
  double BBGEM_LY[5] = {40.,40.,40.,40.,60.};

  double Xbinwidth = 210.0/100.0;
  double Ybinwidth = 62.0/100.0;

  
  TH1D *hitrate_vs_layer_BBGEM = new TH1D("hitrate_vs_layer_BBGEM","",5,0.5,5.5);
  TH2D *hitrate_vs_X_BBGEM = new TH2D("hitrate_vs_X_BBGEM","Hit rate (Hz/cm^{2}/uA), 15-cm LD2",5,0.5,5.5,100,-1.05,1.05);
  TH2D *hitrate_vs_Y_BBGEM = new TH2D("hitrate_vs_Y_BBGEM","Hit rate (Hz/cm^{2}/uA), 15-cm LD2",5,0.5,5.5,100,-0.31,0.31);

  int nevent = 5;
  
  while( T->GetEntry( nevent++ ) && nevent < 5){
    
    //Add GEMs:
    for( int ihit=0; ihit<G->Earm_BBGEM_hit_nhits; ihit++ ){
      int plane = (*(G->Earm_BBGEM_hit_plane))[ihit];
      hitrate_vs_layer_BBGEM->Fill( plane, weight/BBGEM_area_cm2[plane-1] );
      hitrate_vs_X_BBGEM->Fill( plane, (*(G->Earm_BBGEM_hit_x))[ihit], weight/Xbinwidth/BBGEM_LY[plane-1] );
      hitrate_vs_Y_BBGEM->Fill( plane, (*(G->Earm_BBGEM_hit_y))[ihit], weight/Ybinwidth/BBGEM_LX[plane-1] );
      
      
    }
  }
  

  TString output = "Rootfiles/";

  TFile *f = new TFile(output + Form("GEM_sim_rates_%i_GeV.root",Q2),"recreate");  
  hitrate_vs_layer_BBGEM->Write("hhitrate_vs_layer_BBGEM");
  hitrate_vs_X_BBGEM->Write("hhitrate_vs_X_BBGEM");
  hitrate_vs_Y_BBGEM->Write("hhitrate_vs_Y_BBGEM");
  

  
} //end of function
