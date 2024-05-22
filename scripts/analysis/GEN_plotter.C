

#include "../../include/gen-ana.h"

const int npoints = 4;
const int npoints_new = 3;

double Q2[npoints] = {1.16,1.72,2.48,3.42};
double GEGM_old[npoints] = {-0.1247 ,-0.141,-0.209,-0.247};
double GEGM_old_err[npoints] = {0.0150,0.0180,0.0294,0.0394};


double Q2_new[npoints_new] = {2.90,6.62,9.48};
//double GEGM_new[npoints_new] = {-0.2169,-0.6422,-0.3592};
double GEGM_new[npoints_new] = {-0.2169,-0.3101,-0.3592};
double GEGM_stat_err_new[npoints_new] = {0.0213,0.0998,0.3991};
double GEGM_sys_err_new[npoints_new] = {0.0143,0.0251,0.0731};
double GEGM_err_new[npoints_new];

double u = -1.9130427;
double kappa_d = -2.03;
double kappa_u = 1.67;


Double_t PSM_theory(Double_t *x,Double_t *par){
  return (0.4133*x[0] - 0.1038*pow(x[0],2) + 0.01252*pow(x[0],3) + 0.007748*pow(x[0],4))/(1 +  0.9771*x[0] - 0.3506*pow(x[0],2) + 0.09863*pow(x[0],3) + 0.009927*pow(x[0],4) + 0.000002173*pow(x[0],5));
}

double F1F1_to_GEGM(double F1, double F2, double Q2){

  double tau = Q2 / (4 * constant::Mn * constant::Mn);
  double GE = F1 - tau*F2;
  double GM = F1 + F2;
  
  return u*GE / GM;
}


TGraph *make_graph(TString var_name){
  
  TString DB_dir = "../../DB/";
  fstream file; file.open(DB_dir + var_name);

  vector<double> Q2, FF;

  string line;
  int iline = 0;
      
  while (getline(file, line)) {
    // Create a stringstream to parse the line
    stringstream ss(line);
    string cell;
    iline++;
	
    vector<string> val;

    // Split the line into cells using a comma as a delimiter
    while (getline(ss, cell, ',')) {
      val.push_back(cell);  // Put one line into vectros
    }
    if(val[0].substr(0, 1) == "#") continue;
    
    Q2.push_back(stod(val[0]));
    FF.push_back(stod(val[1]));
   
  }

  TGraph *g = new TGraph(Q2.size(),&Q2[0],&FF[0]);
  
  return g;
}


TGraph *RCQM_theory(){

  vector<double> Q2_plot;
  double Q2_max = 12;
  int nQ2 = 20;
  double Q2_bin_size = Q2_max / nQ2;

  for(int iQ2 = 0; iQ2 <= nQ2; iQ2++)
    Q2_plot.push_back(iQ2*Q2_bin_size + 0.12);
  
  TGraph *g_F1u = make_graph("F1u.csv");
  TGraph *g_F1d = make_graph("F1d.csv");
  TGraph *g_F2u = make_graph("F2u.csv");
  TGraph *g_F2d = make_graph("F2d.csv");

  vector<double> GEGM_RCQM;
  
  for(int iQ2 = 0; iQ2 < Q2_plot.size(); iQ2++){
    double Q2 = Q2_plot[iQ2];

    double F1u = g_F1u->Eval(Q2) / (pow(Q2,4));
    double F1d = g_F1d->Eval(Q2) / (pow(Q2,4));
    double F2u = g_F2u->Eval(Q2) * kappa_u / (pow(Q2,4));
    double F2d = g_F2d->Eval(Q2) * kappa_d / (pow(Q2,4));
    
    double F1 = 2.0 / 3.0 * F1d - 1.0 / 3.0 * F1u;
    double F2 = 2.0 / 3.0 * F2d - 1.0 / 3.0 * F2u;

    GEGM_RCQM.push_back(F1F1_to_GEGM(F1, F2, Q2));
    
  }

  TGraph *g = new TGraph(Q2_plot.size(),&Q2_plot[0],&GEGM_RCQM[0]);

  return g;
}

void GEN_plotter(){

  TF1 *PSM_theory_curve = new TF1("PSM_theory_curve", PSM_theory, 0, 12,1);

  TGraph *g_RCQM = RCQM_theory();

  for(int i=0; i < npoints; i++){
    GEGM_old[i] *= u;
    GEGM_old_err[i] *= abs(u);
  }

  for(int i=0; i < npoints_new; i++){
    GEGM_new[i] *= u;
    GEGM_err_new[i] = sqrt( GEGM_stat_err_new[i]*GEGM_stat_err_new[i] + GEGM_sys_err_new[i]*GEGM_sys_err_new[i] );
    GEGM_err_new[i] *= abs(u);
    cout<<GEGM_new[i]<<" "<<GEGM_err_new[i]<<endl;
  }

  TGraphErrors *g = new TGraphErrors(npoints,Q2,GEGM_old,0,GEGM_old_err);
  g->SetMarkerStyle(8);  
  
  TGraphErrors *g_new = new TGraphErrors(npoints_new,Q2_new,GEGM_new,0,GEGM_err_new);
  g_new->SetMarkerStyle(8);  
  g_new->SetMarkerColor(kRed);  
  g_new->SetLineColor(kRed);  
  
  TCanvas *c = new TCanvas("c","",800,600);
  g->SetTitle("GEn-II Results;Q^{2} (GeV^{2});#muG_{E}^{n}/G_{M}^{n}");
  g->Draw("AP");
  g_new->Draw("P");

  PSM_theory_curve->Draw("same");
  g_RCQM->Draw("same");
  
  g->GetXaxis()->SetLimits(0,12);
  g->GetYaxis()->SetRangeUser(-0.5,1.5);

  TLegend *legend1 = new TLegend(0.11,0.64,0.41,0.89);
  legend1->AddEntry(g,"GEn-I Data","p");
  legend1->AddEntry(g_new,"GEn-II Data","p");
  legend1->AddEntry(g_RCQM,"RCQM - Miller","l");
  legend1->AddEntry("PSM_theory_curve","DSE - Roberts","l");
  legend1->SetLineColor(0);
  legend1->Draw("same");

  //g->GetXaxis()->SetLimits(0,4);
  //g->GetYaxis()->SetRangeUser(0,1.2);


}
