//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified February 29, 2024
//
//
//   The purpose of this script is to calculate the helicity
//   of inelastic data.
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
    {"Asymmetry Correction","All asymmetry correction parameters",1}
  };

  DBInfo.cfg = cfg;
  DBInfo.var_req = request;

  DB_load(DBInfo);

  cout<<"---------------------------------------------------------------"<<endl;

}


void Inelastic_contamination(TString cfg = "GEN2"){

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  getDB(cfg);

  TString jmgr_file = "../../config/" + cfg + "_He3.cfg";
  Utilities::KinConf kin_info = Utilities::LoadKinConfig(jmgr_file,1);

   // elastic cut limits
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
  
  double coin_bg_min = kin_info.coin_time_cut[0] + (1 + 3)*kin_info.Nsigma_coin_time*kin_info.coin_time_cut[1];
  double coin_bg_max = kin_info.coin_time_cut[0] + (3 + 3)*kin_info.Nsigma_coin_time*kin_info.coin_time_cut[1];

  analyzed_tree *T_data = Utilities::LoadAnalyzedRootFiles(kin_info,1,0);
  analyzed_tree *T_sim = Utilities::LoadAnalyzedRootFiles(kin_info,0,0);

  distribution_fits *dists = new distribution_fits();

  if(cfg == "GEN2") dists->SetBgShapeOption("pol2");
  else dists->SetBgShapeOption("from data");
  


  /////Set the histograms
  int nbins = 100;
  double hxmin = -4;
  double hxmax = 2.5;
  double W2min_r = -1;
  double W2max_r = 7;
  double W2bin_size = 0.2;

  if(cfg == "GEN2"){
    hxmin = -6;
    hxmax = 3;
    W2min_r = -1;
    W2max_r = 7;
  }
  else if(cfg == "GEN3"){
    W2min_r = -1;
    W2max_r = 7;
  }
  else if(cfg == "GEN4"){
    W2min_r = -2;
    W2max_r = 7;
  }

  const int nbinsW2 = (int)(W2max_r - W2min_r) / W2bin_size;

  double W2_array[nbinsW2];
  double A_elas_array[nbinsW2];
  double A_err_elas_array[nbinsW2];
  int Yp_elas_array[nbinsW2];
  int Ym_elas_array[nbinsW2];
  double A_bg_array[nbinsW2];
  double A_err_bg_array[nbinsW2];
  int Yp_bg_array[nbinsW2];
  int Ym_bg_array[nbinsW2];

  for(int ibin=0; ibin < nbinsW2; ibin++){
    Yp_elas_array[ibin] = 0;
    Ym_elas_array[ibin] = 0;
    Yp_bg_array[ibin] = 0;
    Ym_bg_array[ibin] = 0;
  }
  
  //dx
  TH1F *hdx_data = new TH1F("hdx_data","",nbins,hxmin,hxmax);
  TH1F *hdx_sim_p = new TH1F("hdx_sim_p","",nbins,hxmin,hxmax);
  TH1F *hdx_sim_n = new TH1F("hdx_sim_n","",nbins,hxmin,hxmax);
  TH1F *hdx_bg_data = new TH1F("hdx_bg_data","",nbins,hxmin,hxmax);

  //dxdy
  TH2F *hdxdy_elas = new TH2F("hdxdy_elas","^{3}He HCal Elastics;#Deltay (m);#Deltax (m)",150,-2,2,150,-6,6);

  //W2
  TH1F *hW2_all = new TH1F("hW2_all",cfg +" W^{2} With Different Cuts;W^{2};",nbinsW2,-1,7); // All W2 in the dx cuts
  TH1F *hW2_cut1 = new TH1F("hW2_cut1","",nbinsW2,-1,7); // W2 in inelastic cuts
  TH1F *hW2_cut2 = new TH1F("hW2_cut2","",nbinsW2,-1,7); // W2 in elastic cuts
  
  //2D histogram
  TH2F *hW2_x_elas = new TH2F("hW2_x_elas","",nbinsW2,W2min_r,W2max_r,nbins,hxmin,hxmax);
  TH2F *hW2_x_bg = new TH2F("hW2_x_bg","",nbinsW2,W2min_r,W2max_r,nbins,hxmin,hxmax);

  //Coincidence histogram
  TH1F *hcoin_time_elas = new TH1F("hcoin_time_elas","Elastic Coincidence Time;Coincidence Time (ns);Entries",100,0,200);
  TH1F *hcoin_time_inel = new TH1F("hcoin_time_inel","Inelastic Coincidence Time;Coincidence Time (ns);Entries",100,0,200);

  TCut CutSimP = Form("(W2 > %g && W2 < %g && dy > %g && dy < %g && fnucl == 1) * weight",W2min,W2max,dymin,dymax);
  TCut CutSimN = Form("(W2 > %g && W2 < %g && dy > %g && dy < %g && fnucl == 0) * weight",W2min,W2max,dymin,dymax);
  
  T_sim->fChain->Draw("dx>>hdx_sim_p",CutSimP);
  T_sim->fChain->Draw("dx>>hdx_sim_n",CutSimN);

  int nevent = 0;
  int maxevent = T_data->fChain->GetEntries();
  int N_p = 0, N_m = 0;
  int N_QE = 0;


  while(nevent < maxevent){
    //while(nevent < 10000000){
    T_data->GetEntry(nevent++);  

    ////// Define all the cuts we will use on the data  ////////////////
    bool good_hel = DBInfo.GoodHel[T_data->runnum] && (T_data->helicity == -1 || T_data->helicity == 1);
    bool good_moller = DBInfo.GoodMoller[T_data->runnum];
    bool good_He3 = T_data->He3Pol > 0.01;
    bool good_W2 = T_data->W2 > W2min && T_data->W2 < W2max;
    bool W2_in_range = T_data->W2 > W2min_r && T_data->W2 < W2max_r;
    bool dy_bg_cut = T_data->dy < dy_bg_min || T_data->dy > dy_bg_max;
    bool good_dy_elas = T_data->dy > dymin && T_data->dy < dymax;
    bool good_dx_elas = T_data->dx > dxmin && T_data->dx < dxmax;
    bool good_coin_time = T_data->coin_time > coin_min && T_data->coin_time < coin_max;
    bool bad_coin_time = T_data->coin_time < coin_min - 3*kin_info.coin_time_cut[1] || T_data->coin_time > coin_max + 3*kin_info.coin_time_cut[1];
    bool bad_coin_time_cont = T_data->coin_time > coin_bg_min && T_data->coin_time < coin_bg_max;
    //////////////////////////////////////////////////////////////////////

    int helicity = T_data->helicity;
    helicity *= -1*T_data->IHWP*IHWP_Flip; 

    if(!good_hel) continue;  //Remove events with bad helicity
    if(!good_moller || !good_He3) continue;  //Remove runs with a bad moller measurement or bad He3 measurements

    if(good_dx_elas && dy_bg_cut) hcoin_time_inel->Fill(T_data->coin_time);
    if(good_dx_elas && good_dy_elas) hcoin_time_elas->Fill(T_data->coin_time);

    if(!good_coin_time) continue; 
    if(!W2_in_range) continue; //Cut bins outside of histogram range

    hW2_x_elas->Fill(T_data->W2,T_data->dx);
    if(dy_bg_cut){
      hW2_x_bg->Fill(T_data->W2,T_data->dx);
      if(good_W2) hdx_bg_data->Fill(T_data->dx);
    }
    int W2bin = (int) ((T_data->W2 - W2min_r) / W2bin_size);

    //This is cuts to make our dxdy histogram
    if(good_W2)
      hdxdy_elas->Fill(T_data->dy,T_data->dx);

    if(good_W2 && good_dy_elas) {
      hdx_data->Fill(T_data->dx);
      if(good_dx_elas)
	N_QE++;
    }

    // Here we will different histograms for different cuts
    if(good_dx_elas){
      hW2_all->Fill(T_data->W2);
      if(helicity == 1){
	Yp_elas_array[W2bin]++;
      }
      else if(helicity == -1){
	Ym_elas_array[W2bin]++;
      }
      
      // We count inelastics outside a wide dy region
      if(dy_bg_cut){	
	hW2_cut1->Fill(T_data->W2);
	if(helicity == 1){
	  N_p++;
	  Yp_bg_array[W2bin]++;
	}
	else if(helicity == -1){
	  N_m++;
	  Ym_bg_array[W2bin]++;
	}
      }
      if(good_dy_elas){
	hW2_cut2->Fill(T_data->W2);
      }
    }
  }
  
  for(int ibin=0; ibin < nbinsW2; ibin++){
    W2_array[ibin] = W2min_r + (ibin + 1) * W2bin_size - W2bin_size/2;
    A_elas_array[ibin] = (Yp_elas_array[ibin] - Ym_elas_array[ibin])*1.0/(Yp_elas_array[ibin] + Ym_elas_array[ibin]);
    A_err_elas_array[ibin] = sqrt((1 - A_elas_array[ibin]*A_elas_array[ibin])/(Yp_elas_array[ibin] + Ym_elas_array[ibin])) * 100;
    A_elas_array[ibin] *= 100; //convert to percentage

    A_bg_array[ibin] = (Yp_bg_array[ibin] - Ym_bg_array[ibin])*1.0/(Yp_bg_array[ibin] + Ym_bg_array[ibin]);
    A_err_bg_array[ibin] = sqrt((1 - A_bg_array[ibin]*A_bg_array[ibin])/(Yp_bg_array[ibin] + Ym_bg_array[ibin])) * 100;
    A_bg_array[ibin] *= 100; //convert to percsentage

  }

  TGraphErrors *gA_elas = new TGraphErrors(nbinsW2,W2_array,A_elas_array,0,A_err_elas_array);
  TGraphErrors *gA_bg = new TGraphErrors(nbinsW2,W2_array,A_bg_array,0,A_err_bg_array);


  /////////////////////////////////////////////////////////////////////////

  hW2_x_elas->SetTitle(cfg + " with Elastic Cuts;#Deltax (m);W^{2} (GeV^{2})");
  hW2_x_bg->SetTitle(cfg + " with Background Cuts;#Deltax (m);W^{2} (GeV^{2})");

  TLine *lymin = new TLine(W2min_r,dxmin,W2max_r,dxmin);
  TLine *lymax = new TLine(W2min_r,dxmax,W2max_r,dxmax);
  lymin->SetLineColor(kRed);
  lymax->SetLineColor(kRed);

  TPaveText *pt1 = new TPaveText(.64,.21,.89,.33,"ndc");
  pt1->AddText("Cuts on good tracks");
  pt1->AddText("Coincidence Cuts");
  pt1->AddText(Form("|#Deltay| < %g",dy_n[1]));
  pt1->SetFillColor(0);

  TPaveText *pt2 = new TPaveText(.64,.21,.89,.33,"ndc");
  pt2->AddText("Cuts on good tracks");
  pt2->AddText("Coincidence Cuts");
  pt2->AddText(Form("|#Deltay| > %g",dy_bg_max));
  pt2->SetFillColor(0);
  
  TLegend *legend3 = new TLegend(0.64,0.14,0.89,0.21);
  legend3->AddEntry(lymin,Form("|#Deltax| < %g",dxmax),"l");
  legend3->SetLineColor(0);

  TCanvas *c1 = new TCanvas("c1","",1600,600);
  c1->Divide(2,1);
  c1->cd(1);
  hW2_x_elas->Draw("colz");
  lymin->Draw("same");
  lymax->Draw("same");
  pt1->Draw("same");
  legend3->Draw("same");

  c1->cd(2);
  hW2_x_bg->Draw("colz");
  lymin->Draw("same");
  lymax->Draw("same");
  pt2->Draw("same");
  legend3->Draw("same");
  
  ////////////////////////////////////////////////////////////////

  gA_elas->SetTitle(";W^{2} (GeV^{2});A (%)");
  gA_bg->SetTitle(";W^{2} (GeV^{2});A (%)");
  
  gA_elas->SetMarkerStyle(8);
  gA_bg->SetMarkerStyle(8);

  gA_elas->SetMarkerColor(kBlack);
  gA_elas->SetLineColor(kBlack);
  gA_bg->SetMarkerColor(kRed);
  gA_bg->SetLineColor(kRed);

  hW2_all->SetLineColor(kBlack);
  hW2_cut1->SetLineColor(kRed);
  hW2_cut2->SetLineColor(kGreen);
  
  gA_elas->GetXaxis()->SetLabelSize(0.1);
  gA_elas->GetYaxis()->SetLabelSize(0.1);
  hW2_all->GetYaxis()->SetLabelSize(0.05);
  gA_elas->GetXaxis()->SetTitleSize(0.10);
  gA_elas->GetYaxis()->SetTitleSize(0.10);
  hW2_all->GetYaxis()->SetTitleSize(0.05);
  gA_elas->GetYaxis()->SetTitleOffset(0.28);
  hW2_all->GetYaxis()->SetTitleOffset(0.35);

  TCanvas *c2 = new TCanvas("c2","",800,600);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
  
  pad1->SetBottomMargin(0.00); // Set bottom margin for pad1
  pad2->SetTopMargin(0.00); // Set top margin for pad2
  pad2->SetBottomMargin(0.30);
  
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  pad1->SetGridx(); 
  hW2_all->Draw();
  hW2_cut1->Draw("same");
  hW2_cut2->Draw("same"); 

  TLegend *legend4 = new TLegend(0.11,0.72,0.45,0.89);
  legend4->AddEntry("hW2_all",Form("|#Deltax| < %g",dxmax),"l");
  legend4->AddEntry("hW2_cut1",Form("|#Deltax| < %g & |#Deltay| > %g",dxmax,dy_bg_max),"l");
  legend4->AddEntry("hW2_cut2",Form("|#DeltaR| < %g",dxmax),"l");
  legend4->SetLineColor(0);
  legend4->Draw("same");
  
  pad2->cd();
  pad2->SetGridx();
  gA_elas->Draw("AP");
  gA_bg->Draw("P");

  gA_elas->GetXaxis()->SetLimits(W2min_r,W2max_r);
  gA_elas->GetYaxis()->SetRangeUser(-3,8);
  
  
  //////////////////////////////////////////////////////////////////////
 
  dists->SetDataShape(hdx_data);
  dists->SetPShape(hdx_sim_p);
  dists->SetNShape(hdx_sim_n);
  dists->SetBgShape(hdx_bg_data);
  
  dists->He3_fit_dists();

  //Copy all the result histograms
  TH1F *hdx_data_plot = dists->GetDataHist();
  TH1F *hdx_sim_p_plot = dists->GetPHist();
  TH1F *hdx_sim_n_plot = dists->GetNHist();
  TH1F *hdx_bg_plot = dists->GetBgHist();
  TH1F *hdx_total_fit_plot = dists->GetTotalHist();

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

  
  TCanvas *c3 = new TCanvas("c3","",800,600);
  TBox *box1 = new TBox(-2,-6,dy_bg_min,6);
  TBox *box2 = new TBox(dy_bg_max,-6,2,6);
  box1->SetFillColorAlpha(kRed, 0.4);
  box2->SetFillColorAlpha(kRed, 0.4);

  TPaveText *pt3 = new TPaveText(.374,.762,.625,.892,"ndc");
  pt3->AddText("Cuts on good tracks");
  pt3->AddText("Coincidence Cuts");
  pt3->AddText(Form("%g < W^{2} < %g",W2min,W2max));
  pt3->AddText(Form("|#Deltay| > %g",dy_bg_max));
  pt3->SetFillColor(0);
  
  TLegend *legend5 = new TLegend(0.374,0.717,0.625,0.764);
  legend5->AddEntry(box1,"Bkgd Cut","f");
  legend5->SetLineColor(0);

  hdxdy_elas->Draw("colz");
  box1->Draw("same");
  box2->Draw("same");
  pt3->Draw("same");
  legend5->Draw("same");

  TCanvas *c4 = new TCanvas("c4","",800,600);
  hdx_data_plot->SetTitle(cfg + " Data/Simulation Comparisons;#Deltax (m);Entries");
  hdx_data_plot->Draw();
  hdx_total_fit_plot->Draw("hist same");
  hdx_sim_p_plot->Draw("hist same");
  hdx_sim_n_plot->Draw("hist same");
  hdx_bg_plot->Draw("hist same");

  c4->Update();

  TLine *lmin = new TLine(dxmin,c4->GetUymin(),dxmin,c4->GetUymax()/2);
  TLine *lmax = new TLine(dxmax,c4->GetUymin(),dxmax,c4->GetUymax()/2);
  lmin->SetLineColor(kRed);
  lmax->SetLineColor(kRed);
  lmin->SetLineWidth(4);
  lmax->SetLineWidth(4);
  lmin->Draw("same");
  lmax->Draw("same");

  TLegend *legend = new TLegend(0.59,0.64,0.89,0.89);
  legend->AddEntry("hdx_data","Data","p");
  legend->AddEntry("hdx_total_fit","MC Fit","lf");
  legend->AddEntry("hdx_sim_p","MC p","lf");
  legend->AddEntry("hdx_sim_n","MC n","lf");
  legend->AddEntry("hdx_bg","Background","lf");
  legend->AddEntry(lmin,"QE Cut","l");
  legend->SetLineColor(0);
  legend->Draw("same");

  TCanvas *c5 = new TCanvas("c5","",800,600);
  hcoin_time_inel->Scale(1.0/hcoin_time_inel->Integral());
  hcoin_time_elas->Scale(1.0/hcoin_time_elas->Integral());
  hcoin_time_inel->SetLineColor(kRed);
  hcoin_time_elas->Draw("hist");
  hcoin_time_inel->Draw("same hist");

  double N_in = dists->GetBgYield(dxmin, dxmax);
  
  N_in -= N_QE*(DBInfo.AccidentalFraction + DBInfo.NitrogenFraction + DBInfo.PionFraction);

  double Delta_in = N_p - N_m;
  double Sigma_in = N_p + N_m;
  double A_in = Delta_in / Sigma_in;
  double f_in = 1.0*N_in / N_QE;
  double Aerr = CalcAsymErr(N_p,N_m);
  double Ferr = CalcFractionErr(N_in, N_QE);

  cout<<"Inelastic Asymmetry = "<<A_in<<"  Aerr = "<<Aerr<<endl;
  cout<<"Inelastic Fraction = "<<f_in<<"  Ferr = "<<Ferr<<endl;


  TString plot_dir = "../../plots/";
  TString plot_name = "Inelastic_cont_" + cfg + ".pdf";

  TString outputfile = plot_dir + plot_name;

  c1->Print(outputfile + "(");
  c2->Print(outputfile);
  c3->Print(outputfile);
  c4->Print(outputfile);
  c5->Print(outputfile + ")");

}


