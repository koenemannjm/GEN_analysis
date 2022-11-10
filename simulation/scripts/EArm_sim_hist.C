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

void EArm_sim_hist(int Q2 = 29){

  TString Rootfiles = "/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/simulation/";
  
  const int nfiles = 40;
  const int ngen = 5;
  TString gen_name[ngen] = {"elastic","inelastic","wiser_pi0","wiser_pip","wiser_pim"};
  
  double right_bound = 4.0;
  if(Q2 == 66) right_bound = 5.0;
  if(Q2 == 97) right_bound = 5.0;

  TChain *T = new TChain("T");

  TChain *T_gen[ngen];
  gmn_tree *G_gen[ngen];
  TH1D *maxsig_gen[ngen];
  
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

  double Luminosity_gen[ngen];
  int Nevents_total_gen[ngen] = {0};

  for( int igen = 0; igen < ngen; igen++){
    T_gen[igen] = new TChain("T");
    
    for(int ifile=1; ifile < nfiles; ifile++){
      
      TString filename = Rootfiles + Form("gen_%i_" + gen_name[igen] + "_100000events_job%i.root",Q2,ifile);
      
      ftemp = new TFile(filename,"READ");
      if( !ftemp->IsZombie() ){
	ftemp->GetObject("run_data",rdtemp);
	if( rdtemp != NULL ){
	  int Ngen = rdtemp->fNtries;
	  Nevents_total_gen[igen] += Ngen;
	  Luminosity_gen[igen] = rdtemp->fLuminosity;
	  
	  T_gen[igen]->Add(filename);	
	  
	}
      }
    }

    G_gen[igen] = new gmn_tree(T_gen[igen]);
    maxsig_gen[igen] = new TH1D("maxsig_" + gen_name[igen],"",500,0.0,right_bound);
  }
  
  int nevent = 5;
  TRandom3 num(0);
  

  TH1D *maxsig = new TH1D("maxsig","",500,0.0,right_bound);
  
  while(T->GetEntry(nevent++)){
    //  for (Long64_t i=10; i<nevent; i++) {
    //T->GetEntry(i);
        
    if((G->Earm_BBPSTF1_det_esum + G->Earm_BBSHTF1_det_esum)>0.1 && (*(G->Earm_BBGEM_Track_MID))[0]==0 && (*(G->Earm_BBGEM_Track_NumHits))[0]==5 && G->Earm_BBGEM_Track_ntracks==1 && (*(G->Earm_BBGEM_Track_P))[0]/(G->ev_ep)>0.99){  //used to reproduce the result got by uding det.edum

      // *************** Pre-Shower *****************
      
      //create an array to hold the sums(logic groups) and initi. them to 0.
      double PSTSum [26]; 
      for(int j=0; j<26; j++) PSTSum[j] = 0.0;
    
      for (int ihit=0; ihit<G->Earm_BBPSTF1_hit_nhits; ihit++){
	int row = (*(G->Earm_BBPSTF1_hit_row))[ihit];
	//cout << "Event = " << i << " Hit = " << ihit << " Row = " << row << endl;
	    
	double psmean = 191.0*(*(G->Earm_BBPSTF1_hit_sumedep))[ihit];
	double psedep_sm = num.Gaus(psmean,sqrt(psmean))/191.0;
	 
	//1. Define logic groups
	for(int irow=0; irow<25; irow++){
	  if(row == (irow) || row == (irow+1)) PSTSum[irow] += psedep_sm;
	}
      }  //end of PS loop
    
      // ***************** Shower *****************
       
      //create an array to hold the sums(logic groups) and initi. them to 0.
      double SHTSum [27];
      for(int j=0; j<27; j++) SHTSum[j] = 0.0;
       
      for(int ihit=0; ihit<G->Earm_BBSHTF1_hit_nhits; ihit++){
	int row = (*(G->Earm_BBSHTF1_hit_row))[ihit];

	double shmean = 191.0*(*(G->Earm_BBSHTF1_hit_sumedep))[ihit];
	double shedep_sm = num.Gaus(shmean,sqrt(shmean))/191.0;
				  
	//1. Define logic groups
	for(int irow=0; irow<26; irow++){
	  if(row == (irow) || row == (irow+1)) SHTSum[irow] += shedep_sm;
	}		        
      } //end of SH loop

       
      //2. Create a look-up table (use a STL map)
      //Look-up table based on geometry
      /* PS Gr. | SH Gr. */
      /*  -------------  */
      /*  0-8      0-8   */
      
      /*   9       9+10  */
      /*   .        .    */   //in this region 1 PS gr. corrs. to 2 SH grs.
      /*  17      17+18  */
      
      /*  18-25   19-26  */
      /*  -------------  */

      double maxsigva = 0.0, maxsigva0 = 0.0, maxsigva1 = 0.0;
      for(int gr=0; gr<8; gr++){
	double smgrsum = PSTSum[gr] + SHTSum[gr];
	if(smgrsum > maxsigva) maxsigva = smgrsum;
      }
      for(int gr=9; gr<17; gr++){
	double smgrsum = PSTSum[gr] + SHTSum[gr];
	double crsgrsum = PSTSum[gr] + SHTSum[gr+1];
	if(smgrsum > maxsigva0) maxsigva0 = smgrsum;
	if(crsgrsum > maxsigva0) maxsigva0 = crsgrsum;
      }
      for(int gr=18; gr<25; gr++){
	double crsgrsum = PSTSum[gr] + SHTSum[gr+1];
	if(crsgrsum > maxsigva1) maxsigva1 = crsgrsum;
      }
      if(maxsigva0 > maxsigva) maxsigva = maxsigva0;
      if(maxsigva1 > maxsigva) maxsigva = maxsigva1;

      double weight = G->ev_sigma * G->ev_solang * Luminosity / Nevents_total;
      
      //create a simple hist to hold the max
      maxsig->Fill(maxsigva, weight);
    
    } //for the if statement on top 
     
  } //end of event loop

  // c1->cd();
  // gStyle->SetOptStat(1111111);
  // maxsig->GetXaxis()->SetTitle("PS + SH (GeV)");
  // maxsig->Draw();
  
  maxsig->Fit("gaus","","",0.0,right_bound);
  TF1 *fR = maxsig->GetFunction("gaus");
  fR->SetLineColor(4);

  int nbins = maxsig->GetNbinsX();

  double threshold[nbins];
  double efficiency[nbins];
  double_t maxIn = maxsig->Integral(maxsig->FindFixBin(0),maxsig->FindFixBin(3.0),"");
  double thresh95 = 3.0;
  
  for(int i=1; i<=nbins; i++){
    double thresh = maxsig->GetBinLowEdge(i);
    double int_above = maxsig->Integral(i,nbins);
    double eff = int_above/maxIn;
    //Need to find the threshold at 95% efficiency
    if( eff < 0.95 && thresh < thresh95 ){
      thresh95 = thresh;
    }

    threshold[i-1] = thresh;
    efficiency[i-1] = eff;
  }

  
  TGraph *effic_vs_thresh_fine = new TGraph(nbins, threshold, efficiency);
  gStyle->SetOptStat(1111111);
  effic_vs_thresh_fine->SetTitle("Efficiency vs. Threshold");
  effic_vs_thresh_fine->GetXaxis()->SetTitle("Threshold in GeV");
  effic_vs_thresh_fine->GetYaxis()->SetTitle("Efficiency");
  effic_vs_thresh_fine->SetMarkerStyle(21);
  effic_vs_thresh_fine->SetMarkerColor(2);
  effic_vs_thresh_fine->SetLineColor(2);
  effic_vs_thresh_fine->Draw("AP");

 
  double Q2_double = Q2*1.0/10;

  maxsig->SetMarkerStyle(21);
  maxsig->SetMarkerColor(2);
  maxsig->SetLineColor(2);
  maxsig->SetTitle(Form("GEN   Q^2 = %gGeV^2   Threshold Estimation for BigBite Rates   Kine=Elastic",Q2_double));
  maxsig->GetXaxis()->SetTitle("Energy sum (PS+SH) in GeV");
  maxsig->GetYaxis()->SetTitle("Rate/bin in Hz");
  maxsig->Draw();



  for(int igen = 0; igen < ngen; igen++){
    //Long64_t nevent_el = T_gen[igen]->GetEntries();

    int nevent_el = 5;

   
    while(T_gen[igen]->GetEntry(nevent_el++)){
      //for (Long64_t i=0; i<nevent_el; i++) {
      //T_gen[igen]->GetEntry(i);
      if((G_gen[igen]->Earm_BBPSTF1_det_esum + G_gen[igen]->Earm_BBSHTF1_det_esum)>0.1){  //used to reproduce the result got by uding det.edum
	// *************** Pre-Shower *****************
	//create an array to hold the sums(logic groups) and initi. them to 0.
	double PSTSum [26]; 
	for(int j=0; j<26; j++) PSTSum[j] = 0.0;
	//loop through all hits
	for (int ihit=0; ihit<G_gen[igen]->Earm_BBPSTF1_hit_nhits; ihit++){
	  int row = (*(G_gen[igen]->Earm_BBPSTF1_hit_row))[ihit];    
	  double psmean = 191.0*(*(G_gen[igen]->Earm_BBPSTF1_hit_sumedep))[ihit];
	  double psedep_sm = num.Gaus(psmean,sqrt(psmean))/191.0;
	  //1. Define logic groups
	  for(int irow=0; irow<25; irow++){
	    if(row == (irow) || row == (irow+1)) PSTSum[irow] += psedep_sm;
	  }
	}  //end of PS loop
	// ***************** Shower *****************
	//create an array to hold the sums(logic groups) and initi. them to 0.
	double SHTSum [27];
	for(int j=0; j<27; j++) SHTSum[j] = 0.0;
	for(int ihit=0; ihit<G_gen[igen]->Earm_BBSHTF1_hit_nhits; ihit++){
	  int row = (*(G_gen[igen]->Earm_BBSHTF1_hit_row))[ihit];
	  double shmean = 191.0*(*(G_gen[igen]->Earm_BBSHTF1_hit_sumedep))[ihit];
	  double shedep_sm = num.Gaus(shmean,sqrt(shmean))/191.0;
	//1. Define logic groups
	  for(int irow=0; irow<26; irow++){
	    if(row == (irow) || row == (irow+1)) SHTSum[irow] += shedep_sm;
	  }		        
	} //end of SH loop
	//2. Create a look-up table (use a STL map)
	//Look-up table based on geometry
	/* PS Gr. | SH Gr. */
	/*  -------------  */
	/*  0-8      0-8   */
	
	/*   9       9+10  */
	/*   .        .    */   //in this region 1 PS gr. corrs. to 2 SH grs.
	/*  17      17+18  */
	
	/*  18-25   19-26  */
	/*  -------------  */
	double maxsigva = 0.0;   
	for(int gr=0; gr<26; gr++){
	  double smgrsum = PSTSum[gr] + SHTSum[gr];
	  double crsgrsum = PSTSum[gr] + SHTSum[gr+1];
	  if(smgrsum > maxsigva) maxsigva = smgrsum;
	  if(crsgrsum > maxsigva) maxsigva = crsgrsum;
	}

	double weight = G_gen[igen]->ev_sigma * G_gen[igen]->ev_solang * Luminosity_gen[igen] / Nevents_total_gen[igen];
	
	//create a simple hist to hold the max
	maxsig_gen[igen]->Fill(maxsigva, weight);
      } //for the if statement on top 
    } //end of event loop

    maxsig_gen[igen]->SetMarkerStyle(20);
    maxsig_gen[igen]->SetMarkerColor(2);
    maxsig_gen[igen]->SetLineColor(2);
    maxsig_gen[igen]->GetXaxis()->SetTitle("Energy sum (PS+SH) in GeV");
    maxsig_gen[igen]->GetYaxis()->SetTitle("Rate/bin in Hz");
  }



  TString output = "Rootfiles/";

  TFile *f = new TFile(output + Form("EArm_sim_rates_%i_GeV.root",Q2),"recreate");  
  effic_vs_thresh_fine->Write("hrate_vs_th_good_elastics");  
  maxsig->Write("hrate_good_elastics");

  for(int igen = 0; igen < ngen; igen++)
    maxsig_gen[igen]->Write("hrate_" + gen_name[igen]);

  
  
} //end of function
