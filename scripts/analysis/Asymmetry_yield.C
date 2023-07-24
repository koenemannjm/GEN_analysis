
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified July 7, 2023
//
//
//   The purpose of this script is to take analyzed data and
//   to calculate the asymmetry. It requires the output root
//   file from QuasiElastic_ana.C.
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
#include <vector>
#include <iostream>

#include "TCut.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TChain.h"
#include "TVector3.h"
#include "TStopwatch.h"
#include "TTreeFormula.h"
#include "TLorentzVector.h"

#include "../../include/gen-ana.h"
#include "../../dflay/src/JSONManager.cxx"

void Asymmetry_yield(const char *configfilename,std::string filebase="../outfiles/QE_data")
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // Define a clock to get macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);

  // seting up the desired SBS configuration
  TString conf = jmgr->GetValueFromKey_str("GEN_config");
  int sbsmag = jmgr->GetValueFromKey<int>("SBS_magnet_percent");
  SBSconfig sbsconf(conf, sbsmag);
  sbsconf.Print();

  // choosing nucleon type 
  std::string Ntype = jmgr->GetValueFromKey_str("Ntype");

  int model = jmgr->GetValueFromKey<int>("model");

  int IHWP_Flip = jmgr->GetValueFromKey<int>("IHWP_Flip");

  TString inFile = Form("%s_" + sbsconf.GetSBSconf() + "_sbs%dp_nucleon_%s_model%d.root", 
			 filebase.c_str(),   sbsconf.GetSBSmag(), Ntype.c_str(), model);
  TFile *fin = new TFile(inFile.Data(), "read");
  TTree *T = (TTree*)fin->Get("Tout");

  T->SetBranchStatus("*",0);

  // This is the variables we need from the analyzed root file
  int runnum;   setrootvar::setbranch(T,"runnum","",&runnum);
  bool WCut;   setrootvar::setbranch(T,"WCut","",&WCut);
  bool pCut;   setrootvar::setbranch(T,"pCut","",&pCut);
  bool nCut;   setrootvar::setbranch(T,"nCut","",&nCut);
  bool coinCut;   setrootvar::setbranch(T,"coinCut","",&coinCut);
  double W2;   setrootvar::setbranch(T,"W2","",&W2);
  double dx;   setrootvar::setbranch(T,"dx","",&dx);
  double dy;   setrootvar::setbranch(T,"dy","",&dy);
  double coin_time;   setrootvar::setbranch(T,"coinT_trig","",&coin_time);
  double hcal_time;   setrootvar::setbranch(T,"hcal_time","",&hcal_time);
  double hodo_time[1000];   setrootvar::setbranch(T,"hodo_time","",&hodo_time);
  int helicity;   setrootvar::setbranch(T,"helicity","",&helicity);
  int IHWP;   setrootvar::setbranch(T,"IHWP","",&IHWP);

  //Totals for all runs
  int Yp_n_total = 0;
  int Ym_n_total = 0;
  int ncut_n_total = 0;
  int Yp_p_total = 0;
  int Ym_p_total = 0;
  int ncut_p_total = 0;

  //Y+ and Y- neutron helicity yields
  //For individual runs
  int Yp_n = 0;
  int Ym_n = 0;
  int ncut_n = 0;
  
  //Y+ and Y- proton helicity yields
  //For individual runs
  int Yp_p = 0;
  int Ym_p = 0;
  int ncut_p = 0;

  // Setup some variables for helicity per run
  vector<double> runs;
  vector<double> A_n;
  vector<double> A_n_err;
  vector<double> A_p;
  vector<double> A_p_err;
  vector<double> p_nevents_points, n_nevents_points, A_n_nevents, A_p_nevents, A_n_err_nevents, A_p_err_nevents;

  int nevent = 0;
  int currentrunnum = -1;
  int plotrunnum = 10;
  int ndata_per_point = 3000;    //Used to make the plot of A over time
  int ndata_per_point_min = 2000;//Don't plot A if its less than this many events
  int ndata  = 0;
  int past_IHWP = 0;


  //This will loop over all events and record the helicity inside the proton or
  //neutron spots. It will also record the asymmetry over a number of events to 
  //make plots of asymmetry over time
  while(T->GetEntry(nevent++)){

    if(IHWP == 1) helicity *= -1*IHWP_Flip;
    else if(IHWP == -1) helicity *= 1*IHWP_Flip;
    else continue;

    if(past_IHWP == 0) past_IHWP = IHWP;  //initialize the IHWP state

    //If the IHWP state has changed or we have reached the desired number of events then we want to add this point to the plot
    if(IHWP != past_IHWP || ndata == ndata_per_point){      
      
      //This is used to add space between IHWP states on the plot so it looks nicer
      if(IHWP != past_IHWP) plotrunnum += 20;

      //Calculate the Asymmetries for this data subset
      double Asym_n = (Yp_n - Ym_n)*1.0/(Yp_n + Ym_n);
      double Error_n = sqrt((1 - Asym_n*Asym_n)/ncut_n);
      
      double Asym_p = (Yp_p - Ym_p)*1.0/(Yp_p + Ym_p);
      double Error_p = sqrt((1 - Asym_p*Asym_p)/ncut_p);
      
      // Add them to the array if they are above the minimum
      if(ndata > ndata_per_point_min){
	runs.push_back(plotrunnum++);
	A_n.push_back(Asym_n*100);
	A_n_err.push_back(Error_n*100);
	
	A_p.push_back(Asym_p*100);
	A_p_err.push_back(Error_p*100);

      }      

      // Set everything back to 0 for the next data point
      Yp_n = 0;
      Ym_n = 0;
      Yp_p = 0;
      Ym_p = 0;
      ncut_n = 0;
      ncut_p = 0;
      ndata = 0;      

    } //End loop over plotting data points

    //Cuts for good elastic events
    if(W2 > 0 && W2 < 1.6 && coinCut){

      if(nCut){ //Cut on neutron spot
	if(helicity == 1){
	    Yp_n++;
	    ncut_n++;
	    Yp_n_total++;
	    ncut_n_total++;
	    ndata++;
	  }
	  if(helicity == -1){
	    Ym_n++;
	    ncut_n++;
	    Ym_n_total++;
	    ncut_n_total++;
	    ndata++;
	  }
      }
      if(pCut){ //Cut on proton spot
	if(helicity == 1){
	    Yp_p++;
	    ncut_p++;
	    Yp_p_total++;
	    ncut_p_total++;
	  }
	  if(helicity == -1){
	    Ym_p++;
	    ncut_p++;
	    Ym_p_total++;
	    ncut_p_total++;
	  }
      }

      //These variables are used to plot the Asymmetry vs total events in increments of 10000 events
      if(pCut && ncut_p_total % 10000 == 0){
       
	double Asym_p = (Yp_p_total - Ym_p_total)*1.0/(Yp_p_total + Ym_p_total);
	double Error_p = sqrt((1 - Asym_p*Asym_p)/ncut_p_total);

	p_nevents_points.push_back(ncut_p_total);
	A_p_nevents.push_back(Asym_p*100);	
	A_p_err_nevents.push_back(Error_p*100);

      }

      if(nCut && ncut_n_total % 10000 == 0){
	
	double Asym_n = (Yp_n_total - Ym_n_total)*1.0/(Yp_n_total + Ym_n_total);
	double Error_n = sqrt((1 - Asym_n*Asym_n)/ncut_n_total);
      
	n_nevents_points.push_back(ncut_n_total);
	A_n_nevents.push_back(Asym_n*100);	
	A_n_err_nevents.push_back(Error_n*100);
	
      }

    } // End cut over good elastic events

    past_IHWP = IHWP;

  }

  // Draw plots of Asymmetry vs run number
  TGraphErrors *g_A_n = new TGraphErrors(runs.size(),&runs[0],&A_n[0],0,&A_n_err[0]);
  TGraphErrors *g_A_p = new TGraphErrors(runs.size(),&runs[0],&A_p[0],0,&A_p_err[0]);

  TCanvas *c = new TCanvas("c","",800,600);
  g_A_n->Draw("AP");
  g_A_n->SetMarkerStyle(8);
  g_A_n->SetTitle("Asymmetry vs Run Number;Run Number;Asymmetry (%)");
    
  g_A_p->Draw("P same");
  g_A_p->SetMarkerStyle(8);
  g_A_p->SetMarkerColor(kRed);

  //g_A_n->GetYaxis()->SetRangeUser(-15,15);
  g_A_n->GetYaxis()->SetRangeUser(-8,8);
  g_A_n->GetXaxis()->SetLimits(0,160);

  
  TLegend *legend = new TLegend(0.6,0.79,0.89,0.89);
  legend->AddEntry(g_A_n,"(e,e'n) events","p");
  legend->AddEntry(g_A_p,"(e,e'p) events","p");
  legend->SetLineColor(0);
  legend->Draw("same");


  //Draw plots of Asymmetry vs number of events
  TGraphErrors *g_A_n_nevents = new TGraphErrors(n_nevents_points.size(),&n_nevents_points[0],&A_n_nevents[0],0,&A_n_err_nevents[0]);
  TGraphErrors *g_A_p_nevents = new TGraphErrors(p_nevents_points.size(),&p_nevents_points[0],&A_p_nevents[0],0,&A_p_err_nevents[0]);
  
  
  TCanvas *c2 = new TCanvas("c2","",800,600);
  g_A_p_nevents->Draw("AP");
  g_A_p_nevents->SetMarkerStyle(8);
  g_A_p_nevents->SetTitle("Asymmetry vs Total Events;Events Analyzed;Asymmetry (%)");
  g_A_p_nevents->SetMarkerColor(kRed);
    
  g_A_n_nevents->Draw("P same");
  g_A_n_nevents->SetMarkerStyle(8);

  //g_A_n_nevents->GetYaxis()->SetRangeUser(-15,15);
  //g_A_n_nevents->GetYaxis()->SetRangeUser(-8,8);
  //g_A_n_nevents->GetXaxis()->SetLimits(0,250);

  /*
  TLegend *legend = new TLegend(0.6,0.79,0.89,0.89);
  legend->AddEntry(g_A_n,"(e,e'n) events","p");
  legend->AddEntry(g_A_p,"(e,e'p) events","p");
  legend->SetLineColor(0);
  */
  legend->Draw("same");
  

  //Calculate total asymmetry and print it to the screen
  double A_n_total = (Yp_n_total - Ym_n_total)*1.0/(Yp_n_total + Ym_n_total);
  double A_p_total = (Yp_p_total - Ym_p_total)*1.0/(Yp_p_total + Ym_p_total);
  
  double Err_n = sqrt((1 - A_n_total*A_n_total)/ncut_n_total);
  double Err_p = sqrt((1 - A_p_total*A_p_total)/ncut_p_total);

  cout<<"Total Neutron Events "<<ncut_n_total<<endl;
  cout<<"Y+ "<<Yp_n_total<<endl;
  cout<<"Y- "<<Ym_n_total<<endl;
  cout<<"A = "<<A_n_total*100<<"% +/- "<<Err_n*100<<"%"<<endl;
  cout<<endl;
  cout<<"Total Proton Events "<<ncut_p_total<<endl;
  cout<<"Y+ "<<Yp_p_total<<endl;
  cout<<"Y- "<<Ym_p_total<<endl;
  cout<<"A = "<<A_p_total*100<<"% +/- "<<Err_p*100<<"%"<<endl;

  //Save the pdfs
  c->Print("../../plots/" + sbsconf.GetSBSconf() + "_Asym.png");


}
