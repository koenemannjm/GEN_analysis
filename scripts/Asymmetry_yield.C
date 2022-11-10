#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
#include "TEventList.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "../include/gen-ana.h"






void Asymmetry_yield(TString cfg){

  TString configfile = "../config/" + cfg;

  LoadConfig(configfile);
  
  sbstheta *= constant::pi/180.0;
  bbtheta *= constant::pi/180.0;
  
  ifstream runlist("../config/GEN2_run_list.cfg");

  TString run;
  TString Rootfiles = "/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/Rootfiles/GEN2/He3/";
  
  vector<TString> run_list;
  vector<double> helicity_ratio;
  
  
  int MAXNTRACKS=10;
  
  //variables we need are BigBite track px,py,pz and sbs hcal x, y, e
  
  double ntrack;
  
  double epx[MAXNTRACKS];
  double epy[MAXNTRACKS];
  double epz[MAXNTRACKS];
  double ep[MAXNTRACKS];
  
  double vx[MAXNTRACKS];
  double vy[MAXNTRACKS];
  double vz[MAXNTRACKS];
  
  double xhcal,yhcal,ehcal,ps_e;
  double helicity, IHWP;

  TVector3 hcal_origin( -hcaldist*sin(sbstheta), 0, hcaldist*cos(sbstheta) );
  
  TVector3 hcal_zaxis = hcal_origin.Unit();
  TVector3 hcal_xaxis(0,-1,0);
  TVector3 hcal_yaxis = hcal_zaxis.Cross( hcal_xaxis ).Unit();
  
  int Yp_n_total = 0;
  int Ym_n_total = 0;
  int ncut_n_total = 0;
  int Yp_p_total = 0;
  int Ym_p_total = 0;
  int ncut_p_total = 0;
  int nevents_total = 0;

  vector<double> runs;
  vector<double> A_n;
  vector<double> A_n_err;
  vector<double> A_p;
  vector<double> A_p_err;

  TString currentline;
  while( currentline.ReadLine( runlist ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      
      int ntokens = tokens->GetEntries();
      if( ntokens == 1 ){
	TString skey = ( (TObjString*) (*tokens)[0] )->GetString();
	run = skey;

	TChain *T = new TChain("T");
	
	T->Add(Rootfiles + "/good_helicity/e1209016_fullreplay_*" + run + "*");
	
	TEventList *elist = new TEventList("elist","");
	
	T->Draw(">>elist",globalcut);
		
	T->SetBranchStatus("*",0);
	
	std::vector<std::string> trvar = {"n","vz","px","py","pz","p"};
	std::vector<void*> trvar_mem = {&ntrack,&vz,&epx,&epy,&epz,&ep};
	setrootvar::setbranch(T,"bb.tr",trvar,trvar_mem);
	
	std::vector<std::string> hcalvar = {"x","y","e"};
	std::vector<void*> hcalvar_mem = {&xhcal,&yhcal,&ehcal};
	setrootvar::setbranch(T,"sbs.hcal",hcalvar,hcalvar_mem);
	
	std::vector<std::string> helicityvar = {"hel"};
	std::vector<void*> helicityvar_mem = {&helicity};
	setrootvar::setbranch(T,"scalhel",helicityvar,helicityvar_mem);
	
	std::string bbcalvar = "e";
	void* bbcalvar_mem = &ps_e;
	setrootvar::setbranch(T,"bb.ps",bbcalvar,bbcalvar_mem);
	
	string hallvar = "IGL1I00OD16_16"; //1 is IN, 0 is OUT
	void* hallvar_mem = &IHWP;
	setrootvar::setbranch(T,hallvar,"",hallvar_mem);
	
	
	TLorentzVector Pbeam(0,0,Ebeam,Ebeam);
	TLorentzVector Ptarg(0,0,0,0.5*(0.938272+0.939565));

	long nevent=0;
		
	int ntotal = T->GetEntries();
	nevents_total += ntotal;
	cout<<ntotal<<" events found"<<endl;
	

	//Y+ and Y- neutron helicity yields
	int Yp_n = 0;
	int Ym_n = 0;
	int ncut_n = 0;

	//Y+ and Y- proton helicity yields
	int Yp_p = 0;
	int Ym_p = 0;
	int ncut_p = 0;
	

	
	while( T->GetEntry(elist->GetEntry(nevent++)) ){
	  //if(nevent % 100000 == 0) cout<<nevent<<endl;
	  
	  if( ntrack == 1.0 ){
	    TLorentzVector kprime( epx[0], epy[0], epz[0], ep[0] );
	    TLorentzVector q = Pbeam - kprime;
	    
	    TVector3 qdir = q.Vect().Unit();
	    
	    TVector3 vertex(0,0,vz[0]);
	    
	    double sintersect = (hcal_origin-vertex).Dot( hcal_zaxis )/qdir.Dot( hcal_zaxis );

	    TVector3 hcal_intersect = vertex + sintersect * qdir; 
	    
	    double xhcal_expect = hcal_intersect.Dot( hcal_xaxis );
	    double yhcal_expect = hcal_intersect.Dot( hcal_yaxis );
	    
	    double dy = yhcal - yhcal_expect;
	    double dx = xhcal - xhcal_expect;
	    
	    double W2recon = (Ptarg + q).M2();
	    
	    //cut on neutron spot
	    if(pow(dy - ny_mean,2)/ny_sigma*ny_sigma + pow(dx - nx_mean,2)/nx_sigma*nx_sigma < nsigma*nsigma && ntrack == 1 && ps_e > 0.15 && W2recon < 2.0){
	      
	      if(IHWP == 1) helicity *= -1;
	      else if(IHWP == -1) helicity *= 1;
	      else continue;
	      

	      if(helicity == 1){
		Yp_n++;
		ncut_n++;
	      }
	      if(helicity == -1){
		Ym_n++;
		ncut_n++;
	      }
	    }
	    if(pow(dy - py_mean,2)/py_sigma*py_sigma + pow(dx - px_mean,2)/px_sigma*px_sigma < nsigma*nsigma && ntrack == 1 && ps_e > 0.15 && W2recon < 2.0){
	      
	      if(IHWP == 1) helicity *= -1;
	      else if(IHWP == -1) helicity *= 1;
	      else continue;
	    

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
	}
	
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

	runs.push_back(run.Atof());
	A_n.push_back(Asym_n*100);
	A_n_err.push_back(Error_n*100);

	A_p.push_back(Asym_p*100);
	A_p_err.push_back(Error_p*100);

	elist->Delete(); 
	
      }
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

  g_A_n->GetYaxis()->SetRangeUser(-7,7);

  
  TLegend *legend = new TLegend(0.6,0.7,0.89,0.89);
  legend->AddEntry(g_A_n,"(e,e'n) events","p");
  legend->AddEntry(g_A_p,"(e,e'p) events","p");
  legend->SetLineColor(0);
  legend->Draw("same");
  

  double A_n_total = (Yp_n_total - Ym_n_total)*1.0/(Yp_n_total + Ym_n_total);
  double A_p_total = (Yp_p_total - Ym_p_total)*1.0/(Yp_p_total + Ym_p_total);
  
  double Err_n = sqrt((1 - A_n_total*A_n_total)/ncut_n_total);
  double Err_p = sqrt((1 - A_p_total*A_p_total)/ncut_p_total);

  cout<<"\n\n"<<"Total Events Analyzed "<<nevents_total<<"\n\n";

  cout<<"Total Neutron Events "<<ncut_n_total<<endl;
  cout<<"Y+ "<<Yp_n_total<<endl;
  cout<<"Y- "<<Ym_n_total<<endl;
  cout<<"A = "<<A_n_total*100<<"% +/- "<<Err_n*100<<"%"<<endl;
  cout<<endl;
  cout<<"Total Proton Events "<<ncut_p_total<<endl;
  cout<<"Y+ "<<Yp_p_total<<endl;
  cout<<"Y- "<<Ym_p_total<<endl;
  cout<<"A = "<<A_p_total*100<<"% +/- "<<Err_p*100<<"%"<<endl;

  c->Print("../plots/GEN2_Asym.png");

} 
