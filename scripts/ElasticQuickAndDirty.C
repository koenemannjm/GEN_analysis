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

void ElasticQuickAndDirty( TString cfg, double Ebeam=4.291, double bbtheta=29.5, double sbstheta=34.7, double sbsdist=2.8, double hcaldist=17.0, double sbsfieldscale=1.0){

  sbstheta *= TMath::Pi()/180.0;
  bbtheta *= TMath::Pi()/180.0;

  ifstream configfile("GEN2_" + cfg + ".cfg");

  TChain *C = new TChain("T");

  TString currentline;

  cout << endl << "Chaining all the ROOT files.." << endl;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline);
    }
  }

  TCut globalcut = "bb.ps.e>0.15&&abs(bb.tr.vz)<0.27&&sbs.hcal.nclus>0&&bb.tr.n==1";
  
  TEventList *elist = new TEventList("elist","");
  
  C->Draw(">>elist",globalcut);

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

  double xhcal,yhcal,ehcal;

  C->SetBranchStatus("*",0);
  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.vz",1);
  C->SetBranchStatus("bb.tr.px",1);
  C->SetBranchStatus("bb.tr.py",1);
  C->SetBranchStatus("bb.tr.pz",1);
  C->SetBranchStatus("bb.tr.p",1);

  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);

  C->SetBranchAddress("bb.tr.n",&ntrack);
  C->SetBranchAddress("bb.tr.vz",vz);
  C->SetBranchAddress("bb.tr.px",epx);
  C->SetBranchAddress("bb.tr.py",epy);
  C->SetBranchAddress("bb.tr.pz",epz);
  C->SetBranchAddress("bb.tr.p",ep);
  C->SetBranchAddress("sbs.hcal.x",&xhcal);
  C->SetBranchAddress("sbs.hcal.y",&yhcal);
  C->SetBranchAddress("sbs.hcal.e",&ehcal);

  TLorentzVector Pbeam(0,0,Ebeam,Ebeam);
  TLorentzVector Ptarg(0,0,0,0.5*(0.938272+0.939565));

  long nevent=0;

  double W2min = 0.88-0.4;
  double W2max = 0.88+0.4;

  double np_xcut = -1.0;
  double nsigma = 1;

  double px_mean = -2.379;
  double px_sigma = 0.629;
  double nx_mean = 0.259;
  double nx_sigma = 0.782;

  double py_mean = -0.193;
  double py_sigma = 0.624;
  double ny_mean = -0.285;
  double ny_sigma = 0.773;

  
  if(cfg == "H2"){
    nsigma = 1;
    px_mean = -2.43;
    px_sigma = 0.299;
    py_mean = -0.095;
    py_sigma = 0.445;
  }

  TFile *fout = new TFile("elastic_temp_" + cfg + ".root","RECREATE");

  TH2D *hdxdy_all = new TH2D("hdxdy_all","All events;#Deltay (m);#Deltax (m)",125,-2,2,125,-4,6);
  TH2D *hdxdy_Wcut = new TH2D("hdxdy_Wcut","|W^{2}-0.88|<0.4;Deltay (m);#Deltax (m)",125,-2,2,125,-4,6);
  TH2D *hdxdy_Wanticut = new TH2D("hdxdy_Wanticut","|W^{2}-0.88|<0.4;Deltay (m);#Deltax (m)",125,-2,2,125,-4,6);

  TH1D *hW2_all = new TH1D("hW2_all","All events; W^{2} (GeV^{2});", 250, -1, 4 );
  TH1D *hW2_cut = new TH1D("hW2_cut","HCal p/n spot; W^{2} (GeV^{2});", 250, -1, 4 );

  TH1D *hEHCAL_all = new TH1D("hEHCAL_all","All events; HCAL energy sum (GeV);",250,0,0.5);
  TH1D *hEHCAL_Wcut = new TH1D("hEHCAL_Wcut","|W^{2}-0.88|<0.4;HCAL energy sum (GeV);",250,0,0.5);

  

  TVector3 hcal_origin( -hcaldist*sin(sbstheta), 0, hcaldist*cos(sbstheta) );

  TVector3 hcal_zaxis = hcal_origin.Unit();
  TVector3 hcal_xaxis(0,-1,0);
  TVector3 hcal_yaxis = hcal_zaxis.Cross( hcal_xaxis ).Unit();

  int ntotal = C->GetEntries();
  cout<<ntotal<<" events found"<<endl;
  
  while( C->GetEntry(elist->GetEntry(nevent++) ) ){
    if(nevent % 100000 == 0) cout<<nevent<<endl;

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

      hdxdy_all->Fill( dy, dx );
      double W2recon = (Ptarg + q).M2();

      hW2_all->Fill( W2recon );

      if(cfg == "H2"){
	if(pow(dy - py_mean,2)/py_sigma*py_sigma + pow(dx - px_mean,2)/px_sigma*px_sigma < nsigma*nsigma)
	  hW2_cut->Fill(W2recon);
      }
      else{
	if(pow(dy - py_mean,2)/py_sigma*py_sigma + pow(dx - px_mean,2)/px_sigma*px_sigma < nsigma*nsigma || pow(dy - ny_mean,2)/ny_sigma*ny_sigma + pow(dx - nx_mean,2)/nx_sigma*nx_sigma < nsigma*nsigma)
	  hW2_cut->Fill(W2recon);
      }

      hEHCAL_all->Fill( ehcal );

      if( W2recon > W2min && W2recon < W2max ){
	hdxdy_Wcut->Fill( dy, dx );
	hEHCAL_Wcut->Fill( ehcal );
      } else {
	hdxdy_Wanticut->Fill( dy, dx );
      }
    }
  }

  elist->Delete(); 
  fout->Write();

} 
