//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified February 2, 2024
//
//
//   The purpose of this script is to compare real data and
//   simulated data for the same kinematic point
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

#include "../../include/gen-ana.h"

map<int,double> fGoodHel;

void getDB(TString DB_file){

  // Entries should follow this form:
  //{var name, variable pointer, vairable description, 1/0 (mandatory/not mandatory variable)}
  DBparse::DBRequest request[] = {
    {"Good Helicity", &fGoodHel, "Is the helicity readback good (0/1 = no/yes)", 1}
  };
  
  const int nvar = sizeof(request) / sizeof(request[0]);
  
  DB_load(DB_file,request,nvar);
}

void Asymmetry_checks(TString input = "GEN2",TString tgt = "He3"){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  TString DB_file = "../../DB/Helicity_quality.csv";
  getDB(DB_file);

  TString cfg = input;
  TString cfg_sim = input;
  if(input == "GEN4all"){
    cfg = "GEN4";
    cfg_sim = "GEN4";
  }
  else if(input == "GEN4b") cfg_sim = "GEN4";


  TString jmgr_file = "../../config/" + cfg + "_" + tgt + ".cfg";
  JSONManager *jmgr = new JSONManager(jmgr_file);

  std::vector<int> runnums; jmgr->GetVectorFromKey<int>("runnums",runnums);

  if(input == "GEN4all"){
    TString jmgr_file_4b = "../../config/GEN4b_" + tgt + ".cfg";
    JSONManager *jmgr_4b = new JSONManager(jmgr_file_4b);   
    std::vector<int> runnums_4b; jmgr_4b->GetVectorFromKey<int>("runnums",runnums_4b);

    runnums.insert( runnums.end(), runnums_4b.begin(), runnums_4b.end() );
  }


   // elastic cut limits
  double W2min = jmgr->GetValueFromKey<double>("W2min");
  double W2max = jmgr->GetValueFromKey<double>("W2max");
  //W2min = 0.48;
  //W2max = 1.28;

  double dymin = jmgr->GetValueFromKey<double>("dymin");
  double dymax = jmgr->GetValueFromKey<double>("dymax");
  //dymin = -1.2;
  //dymax = 1.2;

  vector<double> dx_n; jmgr->GetVectorFromKey<double>("dx_n", dx_n);
  double Nsigma_cut_dx_n = jmgr->GetValueFromKey<double>("Nsigma_cut_dx_n");
  double n_dxmin = dx_n[0] - Nsigma_cut_dx_n*dx_n[1];
  double n_dxmax = dx_n[0] + Nsigma_cut_dx_n*dx_n[1];

  int IHWP_Flip = jmgr->GetValueFromKey<int>("IHWP_Flip");

  vector<double> coin_time_cut; jmgr->GetVectorFromKey<double>("coin_time", coin_time_cut);

  TString reaction = "np";
  if(tgt == "H2") reaction = "p";
  TString model = "2";
  if(tgt == "H2") model = "1";

  //Read the He3 run and H2 data files
  TFile *data_file = new TFile("../outfiles/QE_data_" + cfg + "_sbs100p_nucleon_" + reaction + "_model" + model + ".root","read");
  
  TChain *T = new TChain("Tout");
  if(input == "GEN4all"){
    T->Add("../outfiles/QE_data_GEN4_sbs100p_nucleon_np_model2.root");
    T->Add("../outfiles/QE_data_GEN4b_sbs100p_nucleon_np_model2.root");
  }
  else {
    T->Add("../outfiles/QE_data_" + cfg + "_sbs100p_nucleon_np_model2.root");
  }

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
  double coin_time;   setrootvar::setbranch(T,"coin_time","",&coin_time);
  double hcal_time;   setrootvar::setbranch(T,"hcal_time","",&hcal_time);
  double hodo_time[1000];   setrootvar::setbranch(T,"hodo_time","",&hodo_time);
  int helicity;   setrootvar::setbranch(T,"helicity","",&helicity);
  int IHWP;   setrootvar::setbranch(T,"IHWP","",&IHWP);

  /////Set the histograms
  int nbinsdx = 100;
  double xmin = -4;
  double xmax = 2.5;
  double W2min_r = -1;
  double W2max_r = 5;
  double W2bin_size = 0.2;

  if(cfg == "GEN2"){
    xmin = -6;
    xmax = 3;
    W2min_r = -1;
    W2max_r = 5;
  }
  else if(cfg == "GEN3"){
    xmin = -6;
    xmax = 3;
    W2min_r = -1;
    W2max_r = 6;
  }
  else if(cfg == "GEN4"){
    xmin = -6;
    xmax = 3;
    W2min_r = -2;
    W2max_r = 6;
  }

  const int nbinsW2 = (int)(W2max_r - W2min_r) / W2bin_size;
  
  
  //dx
  TH2F *hW2_x = new TH2F("hW2_x","",nbinsW2,W2min_r,W2max_r,nbinsdx,xmin,xmax);
  
  double W2_array[nbinsW2];
  double A_array[nbinsW2];
  double A_err_array[nbinsW2];
  int Yp_array[nbinsW2];
  int Ym_array[nbinsW2];

  for(int ibin=0; ibin < nbinsW2; ibin++){
    Yp_array[ibin] = 0;
    Ym_array[ibin] = 0;
  }

 
  //map<int,bool> good_run;
  //Analysis::check_run_helicity(T,runnums,&good_run);
  
  int nevent = 0;
  int total = 0;
  
  while(T->GetEntry(nevent++)){
   
    if(helicity != -1 && helicity != 1) continue;
    if(W2 < W2min_r || W2 > W2max_r) continue;
    if(!coinCut) continue;
    if(!fGoodHel[runnum]) continue;
    
    hW2_x->Fill(W2,dx);
 
    int W2bin = (int) ((W2 - W2min_r) / W2bin_size);
    helicity *= -1*IHWP*IHWP_Flip;

    if(dx > n_dxmin && dx < n_dxmax){ //Cut around neutron in dx
      if(W2bin == 10) total++;
      if(helicity == 1)
	Yp_array[W2bin]++;
      else if(helicity == -1)
	Ym_array[W2bin]++;
    }
  }
  
  for(int ibin=0; ibin < nbinsW2; ibin++){
    W2_array[ibin] = W2min_r + (ibin + 1) * W2bin_size - W2bin_size/2;
    A_array[ibin] = (Yp_array[ibin] - Ym_array[ibin])*1.0/(Yp_array[ibin] + Ym_array[ibin]);
    A_err_array[ibin] = sqrt((1 - A_array[ibin]*A_array[ibin])/(Yp_array[ibin] + Ym_array[ibin])) * 100;
    A_array[ibin] *= 100; //convert to percentage
  }

  TGraphErrors *gA = new TGraphErrors(nbinsW2,W2_array,A_array,0,A_err_array);
 
  gA->SetTitle(";W^{2} (GeV^{2});A (%)");
  hW2_x->SetTitle(cfg + " Kinematics with Asymmetry;;#Deltax (m)");
  gA->SetMarkerStyle(8);

  gA->GetXaxis()->SetLabelSize(0.1);
  gA->GetYaxis()->SetLabelSize(0.1);
  hW2_x->GetYaxis()->SetLabelSize(0.05);
  gA->GetXaxis()->SetTitleSize(0.10);
  gA->GetYaxis()->SetTitleSize(0.10);
  hW2_x->GetYaxis()->SetTitleSize(0.05);
  gA->GetYaxis()->SetTitleOffset(0.28);
  hW2_x->GetYaxis()->SetTitleOffset(0.35);

  TCanvas *c = new TCanvas("c","",800,600);
  
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
  
  pad1->SetBottomMargin(0.00); // Set bottom margin for pad1
  pad2->SetTopMargin(0.00); // Set top margin for pad2
  pad2->SetBottomMargin(0.30);
  
  pad1->Draw();
  pad2->Draw();

  TLine *lymin = new TLine(W2min_r,n_dxmin,W2max_r,n_dxmin);
  TLine *lymax = new TLine(W2min_r,n_dxmax,W2max_r,n_dxmax);
  lymin->SetLineColor(kRed);
  lymax->SetLineColor(kRed);
  
  pad1->cd();
  pad1->SetGridx();
  hW2_x->Draw("colz");
  lymin->Draw("same");
  lymax->Draw("same");

  TPaveText *pt = new TPaveText(.13,.10,.38,.20,"ndc");
  pt->AddText("Cuts on good tracks");
  pt->AddText("Coincidence Cuts");
  pt->SetFillColor(0);
  pt->Draw("same");

  TLegend *legend = new TLegend(0.13,0.03,0.38,0.10);
  legend->AddEntry(lymin,"#Deltax cut applied below","l");
  legend->SetLineColor(0);
  legend->Draw("same");
  
  pad2->cd();
  pad2->SetGridx();
  gA->Draw("AP");
  gA->GetXaxis()->SetLimits(W2min_r,W2max_r);
  gA->GetYaxis()->SetRangeUser(-3,8);
  
  TLine *l0 = new TLine(W2min_r,0,W2max_r,0);
  l0->Draw("same");

  TString output = "Asymmetry_W2bins_"+input+".pdf";
  
  c->SaveAs("../../plots/" + output);
  
}
