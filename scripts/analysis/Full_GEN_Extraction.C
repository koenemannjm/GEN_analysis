//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified March 4, 2024
//
//
//   The purpose of this script is to calculate the asymmetry
//   given simulated p/n data and real He3 data. This script
//   assumes the data files are already calibrated and analyzed
//   using the script Quasielastic_ana.C
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

#include "../../include/gen-ana.h"


DBparse::DBInfo DBInfo;

// Load database files
void getDB(TString cfg){
  
  cout<<"Attempting to load DB File"<<endl;
  cout<<"---------------------------------------------------------------"<<endl;

   vector<DBparse::DBrequest> request = {
    {"Beam Polarization","Beam Polarization values",1},
    {"Helicity Quality","Helicity readback good? (0/1 = bad/good)",1},
    {"Moller Quality","Moller measurements known? (0/1 = no/yes)",1},
    {"Field Measurement","Magnetic field direction",1},
    {"Asymmetry Correction","All asymmetry correction parameters",1}
  };

  DBInfo.cfg = cfg;
  DBInfo.var_req = request;

  DB_load(DBInfo);

  cout<<"---------------------------------------------------------------"<<endl;

}

// Input types: GEN2, GEN3, GEN4, GEN4b
void Full_GEN_Extraction(TString cfg){

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  getDB(cfg);

  // Load the kinematics and data
  TString jmgr_file = "../../config/" + cfg + "_He3.cfg";
  Utilities::KinConf kin_info = Utilities::LoadKinConfig(jmgr_file,1);

  // Set variables for cuts
  double W2min = kin_info.W2min;
  double W2max = kin_info.W2max;

  double dy_bg_min = kin_info.dymin;
  double dy_bg_max = kin_info.dymax;

  vector<double> dx_n = kin_info.dx_n;
  double Nsigma_dx_n = kin_info.Nsigma_dx_n;
  vector<double> dy_n = kin_info.dy_n;
  double Nsigma_dy_n = kin_info.Nsigma_dy_n;
  double dxmin = dx_n[0] - dx_n[1];
  double dxmax = dx_n[0] + dx_n[1];
  double dymin = dy_n[0] - dy_n[1];
  double dymax = dy_n[0] + dy_n[1];

  int IHWP_Flip = kin_info.IHWP_Flip;

  double coin_min = kin_info.coin_time_cut[0] - kin_info.Nsigma_coin_time*kin_info.coin_time_cut[1];
  double coin_max = kin_info.coin_time_cut[0] + kin_info.Nsigma_coin_time*kin_info.coin_time_cut[1];
  
  // These cuts are for accidental background
  double coin_bg_min = coin_max + 7*kin_info.coin_time_cut[1];
  double coin_bg_max = coin_bg_min + 2*kin_info.Nsigma_coin_time*kin_info.coin_time_cut[1];

  // Set the analyzed tree data
  analyzed_tree *T_data = Utilities::LoadAnalyzedRootFiles(kin_info,1,0);
  analyzed_tree *T_sim = Utilities::LoadAnalyzedRootFiles(kin_info,0,0);

  // Set the distributions
  distribution_fits *dists = new distribution_fits();

  // Set background shape type
  if(cfg == "GEN2") dists->SetBgShapeOption("pol2");
  else dists->SetBgShapeOption("from data");
  

  /////Set the histograms
  int nbins = 100;
  double hxmin = -4;
  double hxmax = 2.5;

  if(cfg == "GEN2"){
    hxmin = -6;
    hxmax = 3;
  }
  
  //dx histograms
  TH1F *hdx_data = new TH1F("hdx_data","",nbins,hxmin,hxmax);
  TH1F *hdx_sim_p = new TH1F("hdx_sim_p","",nbins,hxmin,hxmax);
  TH1F *hdx_sim_n = new TH1F("hdx_sim_n","",nbins,hxmin,hxmax);
  TH1F *hdx_bg_data = new TH1F("hdx_bg_data","",nbins,hxmin,hxmax);
  
  // Cut on the simulated data for QE events from p/n
  // It is also weighted by the "weight" factor from simulation
  TCut CutSimP = Form("(W2 > %g && W2 < %g && dy > %g && dy < %g && fnucl == 1) * weight",W2min,W2max,dymin,dymax);
  TCut CutSimN = Form("(W2 > %g && W2 < %g && dy > %g && dy < %g && fnucl == 0) * weight",W2min,W2max,dymin,dymax);
  
  // Fill the histograms by drawing
  T_sim->fChain->Draw("dx>>hdx_sim_p",CutSimP);
  T_sim->fChain->Draw("dx>>hdx_sim_n",CutSimN);

  int nevent = 0;
  int maxevent = T_data->fChain->GetEntries(); // Used to loop through the tree
  analyzed_info *Asym_total = new analyzed_info(); // Object for total asymmetry

  while(nevent < maxevent){
    T_data->GetEntry(nevent++);  //Get data for one event

    ////// Define all the cuts we will use on the data  ////////////////
    bool good_hel = DBInfo.GoodHel[T_data->runnum] && (T_data->helicity == -1 || T_data->helicity == 1);
    bool good_moller = DBInfo.GoodMoller[T_data->runnum];
    bool good_He3 = T_data->He3Pol > 0.01;
    bool good_W2 = T_data->W2 > W2min && T_data->W2 < W2max;
    bool dy_bg_cut = T_data->dy < dy_bg_min || T_data->dy > dy_bg_max;
    bool good_dy_elas = T_data->dy > dymin && T_data->dy < dymax;
    bool good_dx_elas = T_data->dx > dxmin && T_data->dx < dxmax;
    bool good_coin_time = T_data->coin_time > coin_min && T_data->coin_time < coin_max;
    //bool bad_coin_time = (T_data->coin_time > DBInfo.AccidentalMin && T_data->coin_time < coin_min) || (T_data->coin_time > coin_max && T_data->coin_time < DBInfo.AccidentalMax);
    //bool bad_coin_time_frac = T_data->coin_time > coin_bg_min && T_data->coin_time < coin_bg_max;
    //////////////////////////////////////////////////////////////////////

    int helicity = T_data->helicity;
    helicity *= -1*T_data->IHWP*IHWP_Flip; 

    if(!good_hel) continue;  //Remove events with bad helicity
    if(!good_moller || !good_He3) continue;  //Remove runs with a bad moller measurement or bad He3 measurements
    if(!good_W2) continue;   //Remove events outside of W2 region

    // Use wide dy cut to get bg for fitting
    if(dy_bg_cut && good_coin_time)
      hdx_bg_data->Fill(T_data->dx);  //Fill Bg histogram for fitting later

    if(!good_dy_elas) continue;  //Remove events outside of neutron dy region
    if(!good_coin_time) continue;  //Remove events outside of coincidence time

    hdx_data->Fill(T_data->dx);  // Fill data histogram after all QE cuts

    if(good_dx_elas){ //Cut around neutron in dx
      
      UpdateAverageKinematics(T_data);
      
      // Add raw helicity counts
      Asym_total->IterateRawCount(T_data->runnum, helicity);
      Asym_total->SetPolBeam(T_data->runnum, T_data->datetime, DBInfo.BeamPolTime, DBInfo.BeamPolValue);
      Asym_total->UpdatePolHe3(T_data->runnum, T_data->He3Pol);
      
    }
  }

  // Set the distributions used for fitting
  dists->SetDataShape(hdx_data);
  dists->SetPShape(hdx_sim_p);
  dists->SetNShape(hdx_sim_n);
  dists->SetBgShape(hdx_bg_data);
  
  dists->He3_fit_dists();    // Fit our histograms
  
  ///////// Make plots for the dx results ////////////////////////

  //Copy all the result histograms
  TH1F *hdx_data_plot = dists->GetDataHist();
  TH1F *hdx_sim_p_plot = dists->GetPHist();
  TH1F *hdx_sim_n_plot = dists->GetNHist();
  TH1F *hdx_bg_plot = dists->GetBgHist();
  TH1F *hdx_total_fit_plot = dists->GetTotalHist();
  
  ///// Formatting for all histograms ////////////////////////////
  hdx_data_plot->SetMarkerStyle(kFullCircle);
  hdx_total_fit_plot->SetFillColorAlpha(30,0.5);
  hdx_sim_p_plot->SetFillColorAlpha(kRed,0.3);
  hdx_sim_n_plot->SetFillColorAlpha(kBlue,0.3);
  hdx_bg_plot->SetFillColorAlpha(kMagenta,0.3);

  hdx_total_fit_plot->SetLineStyle(7);
  hdx_sim_p_plot->SetLineStyle(7);
  hdx_sim_n_plot->SetLineStyle(7);
  hdx_bg_plot->SetLineStyle(7);
  
  hdx_total_fit_plot->SetLineColor(30);
  hdx_sim_p_plot->SetLineColor(kRed);
  hdx_sim_n_plot->SetLineColor(kBlue);
  hdx_bg_plot->SetLineColor(kMagenta);

  ///////////////////////////////////////////////////////////////
  
  // Draw histograms
  TCanvas *c1 = new TCanvas("c1","",800,600);
  hdx_data_plot->SetTitle(cfg + " Data/Simulation Comparisons;#Deltax (m);Entries");
  hdx_data_plot->Draw();
  hdx_total_fit_plot->Draw("hist same");
  hdx_sim_p_plot->Draw("hist same");
  hdx_sim_n_plot->Draw("hist same");
  hdx_bg_plot->Draw("hist same");

  c1->Update();

  TLine *lmin = new TLine(dxmin,c1->GetUymin(),dxmin,c1->GetUymax()/3);
  TLine *lmax = new TLine(dxmax,c1->GetUymin(),dxmax,c1->GetUymax()/3);
  lmin->SetLineColor(kRed);
  lmax->SetLineColor(kRed);
  lmin->SetLineWidth(4);
  lmax->SetLineWidth(4);
  lmin->Draw("same");
  lmax->Draw("same");

  TLegend *legend1 = new TLegend(0.59,0.64,0.89,0.89);
  legend1->AddEntry("hdx_data","Data","p");
  legend1->AddEntry("hdx_total_fit","MC Fit","lf");
  legend1->AddEntry("hdx_sim_p","MC p","lf");
  legend1->AddEntry("hdx_sim_n","MC n","lf");
  legend1->AddEntry("hdx_bg","BG from data","lf");
  legend1->AddEntry(lmin,"QE Cut","l");
  legend1->SetLineColor(0);
  legend1->Draw("same");

  TPaveText *pt = new TPaveText(.64,.46,.87,.63,"ndc");
  pt->AddText("Cuts on good tracks");
  pt->AddText("Coincidence Cuts");
  pt->AddText(Form("%g < W^{2} < %g",W2min,W2max));
  pt->AddText(Form("%g < #Deltay < %g",dymin,dymax));
  pt->SetFillColor(0);
  pt->Draw("same");

  //////////////////////////////////////////////////////////////////////

  
  // Set protons and inelastic yields from the histogram fits
  Asym_total->SetNProton(dists->GetPYield(dxmin, dxmax));
  
  // Accidental asymmetry correction, loaded from a different script
  Asym_total->A_acc = DBInfo.AccidentalAsymmetry;
  Asym_total->A_acc_err = DBInfo.AccidentalAsymmetryErr;
  Asym_total->f_acc = DBInfo.AccidentalFraction;
  Asym_total->f_acc_err = DBInfo.AccidentalFractionErr;
  
  // Pion asymmetry correction, loaded from a different script
  Asym_total->A_pion = DBInfo.PionAsymmetry;
  Asym_total->A_pion_err = DBInfo.PionAsymmetryErr;
  Asym_total->f_pion = DBInfo.PionFraction;
  Asym_total->f_pion_err = DBInfo.PionFractionErr;

  // Inelastic asymmetry correction, loaded from a different script
  Asym_total->A_in = DBInfo.InelasticAsymmetry;
  Asym_total->A_in_err = DBInfo.InelasticAsymmetryErr;
  Asym_total->f_in = DBInfo.InelasticFraction;
  Asym_total->f_in_err = DBInfo.InelasticFractionErr;

  // Nitrogen asymmetry correction, loaded from a different script
  Asym_total->f_N2 = DBInfo.NitrogenFraction;
  Asym_total->f_N2_err = DBInfo.NitrogenFractionErr;

  Asym_total->CalcAsymVals();


  //////// Print out the results //////////////////////////////
  cout<<"\n";
  cout<<"-------------------- Result Summary --------------------"<<endl;
  cout<<Form("Q2 = %.2f",Q2_avg)<<endl;
  cout<<"N QE Events = "<<Asym_total->N_raw_p + Asym_total->N_raw_m<<endl;
  cout<<Form("A_phys = %.4f +/- %.4f +/- %.4f (result +/- stat +/- sys)",Asym_total->A_phys,Asym_total->A_phys_stat_err,Asym_total->A_phys_sys_err)<<endl;
  cout<<Form("GE/GM = %.4f +/- %.4f +/- %.4f (result +/- stat +/- sys)",GEGM,GEGM_stat_err,GEGM_sys_err)<<endl;
  
  cout<<"\n";
  cout<<"-------------------- Detailed Results --------------------"<<endl;
  cout<<"------------- Polarizations -------------"<<endl;
  cout<<Form("P He3 Avg = %.4f +/- %.4f",P_He3_avg,P_He3_avg_err)<<endl;
  cout<<Form("P Beam Avg = %.4f +/- %.4f",P_beam_avg,P_beam_avg_err)<<endl;
  cout<<"P Neutron = "<<P_n<<endl;
  cout<<"------------- Asymmetries -------------"<<endl;
  cout<<Form("A Raw = %.4f +/- %.4f",Asym_total->A_raw,Asym_total->A_raw_err)<<endl;
  cout<<Form("A Accidetal = %.4f +/- %.4f",Asym_total->A_acc,Asym_total->A_acc_err)<<endl;
  cout<<Form("A Pion = %.4f +/- %.4f",Asym_total->A_pion,Asym_total->A_pion_err)<<endl;
  cout<<Form("A Proton = %.4f +/- %.4f",Asym_total->A_p,Asym_total->A_p_err)<<endl;
  cout<<Form("A Inelastic = %.4f +/- %.4f",Asym_total->A_in,Asym_total->A_in_err)<<endl;
  cout<<"------------- Fractions -------------"<<endl;
  cout<<Form("F neutron = %.4f",Asym_total->f_n)<<endl;
  cout<<Form("F Accidetal = %.4f +/- %.4f",Asym_total->f_acc,Asym_total->f_acc_err)<<endl;
  cout<<Form("F Nitrogen = %.4f +/- %.4f",Asym_total->f_N2,Asym_total->f_N2_err)<<endl;
  cout<<Form("F Pion = %.4f +/- %.4f",Asym_total->f_pion,Asym_total->f_pion_err)<<endl;
  cout<<Form("F Proton = %.4f +/- %.4f",Asym_total->f_p,Asym_total->f_p_err)<<endl;
  cout<<Form("F Inelastic = %.4f +/- %.4f",Asym_total->f_in,Asym_total->f_in_err)<<endl;


  TString plot_dir = "../../plots/";
  TString plot_name = "Asymmetry_Full_" + cfg + ".pdf";

  TString outputfile = plot_dir + plot_name;

  c1->Print(outputfile);

  TString text_dir = "../outfiles/asym_results/";
  TString text_name = "Asymmetry_" + cfg + ".txt";

  TString outputfile_text = text_dir + text_name;
  
  ofstream textfile;
  textfile.open(outputfile_text);
  
  textfile<<"-------------------- Result Summary --------------------"<<endl;
  textfile<<Form("Q2 = %.2f",Q2_avg)<<endl;
  textfile<<"N QE Events = "<<Asym_total->N_raw_p + Asym_total->N_raw_m<<endl;
  textfile<<Form("A_phys = %.4f +/- %.4f +/- %.4f (result +/- stat +/- sys)",Asym_total->A_phys,Asym_total->A_phys_stat_err,Asym_total->A_phys_sys_err)<<endl;
  textfile<<Form("GE/GM = %.4f +/- %.4f +/- %.4f (result +/- stat +/- sys)",GEGM,GEGM_stat_err,GEGM_sys_err)<<endl;
  
  textfile<<"\n";
  textfile<<"-------------------- Detailed Results --------------------"<<endl;
  textfile<<"------------- Polarizations -------------"<<endl;
  textfile<<Form("P He3 Avg = %.4f +/- %.4f",P_He3_avg,P_He3_avg_err)<<endl;
  textfile<<Form("P Beam Avg = %.4f +/- %.4f",P_beam_avg,P_beam_avg_err)<<endl;
  textfile<<"P Neutron = "<<P_n<<endl;
  textfile<<"------------- Asymmetries -------------"<<endl;
  textfile<<Form("A Raw = %.4f +/- %.4f",Asym_total->A_raw,Asym_total->A_raw_err)<<endl;
  textfile<<Form("A Accidetal = %.4f +/- %.4f",Asym_total->A_acc,Asym_total->A_acc_err)<<endl;
  textfile<<Form("A Pion = %.4f +/- %.4f",Asym_total->A_pion,Asym_total->A_pion_err)<<endl;
  textfile<<Form("A Proton = %.4f +/- %.4f",Asym_total->A_p,Asym_total->A_p_err)<<endl;
  textfile<<Form("A Inelastic = %.4f +/- %.4f",Asym_total->A_in,Asym_total->A_in_err)<<endl;
  textfile<<"------------- Fractions -------------"<<endl;
  textfile<<Form("F neutron = %.4f",Asym_total->f_n)<<endl;
  textfile<<Form("F Accidetal = %.4f +/- %.4f",Asym_total->f_acc,Asym_total->f_acc_err)<<endl;
  textfile<<Form("F Nitrogen = %.4f +/- %.4f",Asym_total->f_N2,Asym_total->f_N2_err)<<endl;
  textfile<<Form("F Pion = %.4f +/- %.4f",Asym_total->f_pion,Asym_total->f_pion_err)<<endl;
  textfile<<Form("F Proton = %.4f +/- %.4f",Asym_total->f_p,Asym_total->f_p_err)<<endl;
  textfile<<Form("F Inelastic = %.4f +/- %.4f",Asym_total->f_in,Asym_total->f_in_err)<<endl;
  
  textfile.close();

  cout<<"\nInfo also written to "<<outputfile_text<<endl;

 
}


