
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

int nmodules = 8;

bool analyze_alignment(TString rootname, TCut globalcut){

  bool good_align = true;

  //TFile *file = new TFile(rootname,"read");
  //TTree *T = (TTree*)file->Get("T");
  TChain *T = new TChain("T");
  T->Add(rootname);

  TH2F *h2 = new TH2F("h2","",nmodules,0,nmodules,100,-0.001,0.001);
  T->Draw("bb.gem.hit.residu:bb.gem.hit.module>>h2",globalcut,"goff");
  
  for(int imod = 0; imod < nmodules; imod++){
    
    TH1D *hproj = h2->ProjectionY(Form("p%i",imod),imod+1,imod+1);
    
    if(abs(hproj->GetMean()) > 1e-4 || hproj->GetStdDev() > 3e-4)
      good_align = false;
    
  }
  h2->Delete();
  return good_align;
}


void check_alignment(const char *configfilename){

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);

  // parsing trees
  std::string rootfile_dir = jmgr->GetValueFromKey_str("rootfile_dir");
  std::vector<int> runnums; jmgr->GetVectorFromKey<int>("runnums",runnums);
  int nruns = jmgr->GetValueFromKey<int>("Nruns_to_ana"); // # runs to analyze
  
  // setting up global cuts
  std::string gcut = jmgr->GetValueFromKey_str("global_cut");
  TCut globalcut = gcut.c_str();

  if (nruns < 1 || nruns > runnums.size()) nruns = runnums.size();
  for (int i=0; i<nruns; i++) {
    std::string rfname = rootfile_dir + Form("/e1209016_fullreplay_%i_stream0_2_seg0_*.root",runnums[i]);

    bool good_align = analyze_alignment(rfname, globalcut);
  
    if(!good_align) cout<<"Bad Run: "<<runnums[i]<<endl;
  }



}
