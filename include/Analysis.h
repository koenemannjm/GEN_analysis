#ifndef Analysis_H
#define Analysis_H

// These are all assumptions, to be refined later
double P_n = 0.86;    //Set for now, but check later
double P_p = -0.03;    //Set for now, but check later
double P_beam = 0.85; //Set for now, but check later
double P_He3 = 0.45;  //Set for now, but check later
double theta_spin = 0;//Set for now, but check later
double phi_spin = 0;  //Set for now, but check later

// Values for expansion to calculate GE/GM
const int nexp = 6;
double T_avg[nexp] = {0};
double T_1_Q2_avg = 0;
double Q2_avg = 0;
int count_avg = 0;

TVector3 TargetPolDirection;

// Basic function to get the integral from a histogram range
double GetYield(TH1F *h, double xmin, double xmax){
  int binmin = h->FindBin(xmin);
  int binmax = h->FindBin(xmax);
  return h->Integral(binmin, binmax);
}

// Sets the target polarization direction
void SetHe3Pol(){
  TargetPolDirection.SetMagThetaPhi(1.0, theta_spin, phi_spin);
}


void UpdateExpansionCoefficients(analyzed_tree *T);
void AsymExpansion();
void GetGEGM();

// Class to handle the distributions and the fitting procedure
class distribution_fits {
 public:

  TString bg_shape_option;   // sets the type of background shape used
  bool bg_shape_set = false;

  // Histograms for all the fit distributions
  TH1F *hdx_data = NULL;
  TH1F *hdx_sim_p = NULL;
  TH1F *hdx_sim_n = NULL;
  TH1F *hdx_bg_fit;
  TH1F *hdx_bg_data = NULL;
  TH1F *hdx_total_fit;

  // Function sets the background shape from a string. 
  // Options: pol3, pol4, gaus, from data
  void SetBgShapeOption(TString bg_shape){
    bg_shape_option = bg_shape;
    bg_shape_set = true;
  };


  // Functions get the integrals from the fits of neutron/proton/bg distributinos
  int GetNYield(double xmin, double xmax){ return GetYield(hdx_sim_n,xmin,xmax);};
  int GetPYield(double xmin, double xmax){ return GetYield(hdx_sim_p,xmin,xmax);};
  int GetBgYield(double xmin, double xmax){ return GetYield(hdx_bg_fit,xmin,xmax);};

  // Set histogram shapes
  void SetDataShape(TH1F *h){hdx_data = h;};
  void SetPShape(TH1F *h){hdx_sim_p = h;};
  void SetNShape(TH1F *h){hdx_sim_n = h;};
  void SetBgShape(TH1F *h){hdx_bg_data = h;};
  
  // Get the histogram data
  TH1F *GetDataHist(){return hdx_data;};
  TH1F *GetPHist(){return hdx_sim_p;};
  TH1F *GetNHist(){return hdx_sim_n;};
  TH1F *GetBgHist(){return hdx_bg_fit;};
  TH1F *GetTotalHist(){return hdx_total_fit;};

  TString GetBgShape(){return bg_shape_option;};

  double fitbg_pol4(double *x, double *par);
  double fitsim(double *x, double *par);
  double fitbg_pol3(double *x, double *par);
  double fitbg_gaus(double *x, double *par);  
  
  void He3_fit_dists();

  // Constructor
  distribution_fits() { }
  
};

// Class to handle analyzing the asymmetry calculations
class analyzed_info {
 public:

  // raw +/- helicities
  int N_raw_p = 0;
  int N_raw_m = 0;

  // accidental background +/- helicities
  int N_bg_p = 0;
  int N_bg_m = 0;

  // Pion +/- helicities
  int N_pion_p = 0;
  int N_pion_m = 0;

  // Number of protons
  int N_proton = 0;

  // Inelastic +/- helicities
  int N_in_p = 0;
  int N_in_m = 0;
  int N_in_cont = 0;
  
  // FSI +/- helicities
  int N_FSI_p = 0;
  int N_FSI_m = 0;

  // Asymmetry values for each contribution
  double A_raw;
  double A_bg;
  double A_p;
  double A_pion;
  double A_in;
  double A_FSI;
  double A_phys;

  // Signal fractions for each contribution
  double f_bg;
  double f_N2;
  double f_p;
  double f_pion;
  double f_in;
  double f_FSI;
  double f_n;

  // Dilutions for each contribution (outdated, no longer in use)
  double D_bg;
  double D_N2;
  double D_p;
  double D_pion;
  double D_in;

  // Asymmetry errors
  double A_raw_err;
  double A_bg_err;
  double A_p_err;
  double A_phys_stat_err;
  double A_phys_sys_err;
  
  // Dilution errors (outdated, no longer in use)
  double D_bg_err;
  double D_N2_err;
  double D_p_err;
  double D_in_err;

  // Signal fraction errors
  double f_bg_err;
  double f_N2_err;
  double f_p_err;
  double f_pion_err;
  double f_in_err;
  double f_FSI_err;


  void SetNProton(int N){ 
    N_proton = N;
  }
  void IterateRawCount(int helicity){ // Add helicity raw counts
    if(helicity == 1) N_raw_p++;
    else if(helicity == -1) N_raw_m++;
  } 
  void IterateBgCount(int helicity){ // Add helicity accidental backgroun counts
    if(helicity == 1) N_bg_p++;
    else if(helicity == -1) N_bg_m++;
  } 
  void IterateInelasticCount(int helicity){ // Add helicity inelastic ounts
   if(helicity == 1) N_in_p++;
    else if(helicity == -1) N_in_m++;
  } 
  void SetNInelContamination(int N){ // Set helicity number inside cuts
    N_in_cont = N;
  }
  

  void CalcAsymValsOld();
  void CalcAsymVals();

  // Constructor
  analyzed_info() { }

};


#endif
