
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

#include "../include/gen-ana.h"
#include "../dflay/src/JSONManager.cxx"

void Asymmetry_yield(const char *configfilename,std::string filebase="outfiles/QE_test")
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // Define a clock to get macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);

  // seting up the desired SBS configuration
  int conf = jmgr->GetValueFromKey<int>("GEN_config");
  int sbsmag = jmgr->GetValueFromKey<int>("SBS_magnet_percent");
  SBSconfig sbsconf(conf, sbsmag);
  sbsconf.Print();

  // choosing nucleon type 
  std::string Ntype = jmgr->GetValueFromKey_str("Ntype");

  int model = jmgr->GetValueFromKey<int>("model");

  TString inFile = Form("%s_GEN%d_sbs%dp_nucleon_%s_model%d_data.root", 
			 filebase.c_str(), sbsconf.GetSBSconf(),  sbsconf.GetSBSmag(), Ntype.c_str(), model);
  TFile *fin = new TFile(inFile.Data(), "read");
  TTree *T = (TTree*)fin->Get("Tout");

  T->SetBranchStatus("*",0);

  int runnum;   setrootvar::setbranch(T,"runnum","",&runnum);
  bool WCut;   setrootvar::setbranch(T,"WCut","",&WCut);
  bool pCut;   setrootvar::setbranch(T,"pCut","",&pCut);
  bool nCut;   setrootvar::setbranch(T,"nCut","",&nCut);
  bool coinCut;   setrootvar::setbranch(T,"coinCut","",&coinCut);
  double dx;   setrootvar::setbranch(T,"dx","",&dx);
  double dy;   setrootvar::setbranch(T,"dy","",&dy);
  double coin_time;   setrootvar::setbranch(T,"coinT_trig","",&coin_time);
  int helicity;   setrootvar::setbranch(T,"helicity","",&helicity);
  int IHWP;   setrootvar::setbranch(T,"IHWP","",&IHWP);

  //For all runs
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

  vector<double> runs;
  vector<double> A_n;
  vector<double> A_n_err;
  vector<double> A_p;
  vector<double> A_p_err;

  int nevent = 0;
  int currentrunnum = -1;

  while(T->GetEntry(nevent++)){

    if(WCut){
      if(IHWP == 1) helicity *= 1;
      else if(IHWP == -1) helicity *= -1;
      else continue;

      if(nCut){
	if(helicity == 1){
	    Yp_n++;
	    ncut_n++;
	  }
	  if(helicity == -1){
	    Ym_n++;
	    ncut_n++;
	  }
      }
      if(pCut){
	if(helicity == 1){
	    Yp_p++;
	    ncut_p++;
	  }
	  if(helicity == -1){
	    Ym_p++;
	    ncut_p++;
	  }
      }

    }

    if(nevent == 1) currentrunnum = runnum;

    if(runnum != currentrunnum){
      
      Yp_n_total += Yp_n;
      Ym_n_total += Ym_n;
      ncut_n_total += ncut_n;
      
      Yp_p_total += Yp_p;
      Ym_p_total += Ym_p;
      ncut_p_total += ncut_p;
      
      double Asym_n = (Yp_n - Ym_n)*1.0/(Yp_n + Ym_n);
      double Error_n = sqrt((1 - Asym_n*Asym_n)/ncut_n);
      
      double Asym_p = (Yp_p - Ym_p)*1.0/(Yp_p + Ym_p);
      double Error_p = sqrt((1 - Asym_p*Asym_p)/ncut_p);
      
      runs.push_back(currentrunnum);
      A_n.push_back(Asym_n*100);
      A_n_err.push_back(Error_n*100);
      
      A_p.push_back(Asym_p*100);
      A_p_err.push_back(Error_p*100);
      
      Yp_n = 0;
      Ym_n = 0;
      Yp_p = 0;
      Ym_p = 0;
      ncut_n = 0;
      ncut_p = 0;
      
      currentrunnum = runnum;
    }
  }

  TGraphErrors *g_A_n = new TGraphErrors(runs.size(),&runs[0],&A_n[0],0,&A_n_err[0]);
  TGraphErrors *g_A_p = new TGraphErrors(runs.size(),&runs[0],&A_p[0],0,&A_p_err[0]);

  TCanvas *c = new TCanvas("c","",800,600);
  g_A_n->Draw("AP");
  g_A_n->SetMarkerStyle(8);
  g_A_n->SetTitle("Asymmetry vs Run Number;Run Number;Asymmetry (%)");
    
  g_A_p->Draw("P same");
  g_A_p->SetMarkerStyle(8);
  g_A_p->SetMarkerColor(kRed);

  g_A_n->GetYaxis()->SetRangeUser(-15,15);

  
  TLegend *legend = new TLegend(0.6,0.79,0.89,0.89);
  legend->AddEntry(g_A_n,"(e,e'n) events","p");
  legend->AddEntry(g_A_p,"(e,e'p) events","p");
  legend->SetLineColor(0);
  legend->Draw("same");
  

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

  c->Print(Form("../plots/GEN%i_Asym.png",conf));


}
