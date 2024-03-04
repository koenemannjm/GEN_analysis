#include <iostream>

#include "../include/Analysis.h"


// This script calculates the coefficients for GE/GM extraction.
// The formalism here is taken from Seamus/Freddy's thesis. It is a 
// fifth order taylor expansion around GE/GM for the asymmetry
void UpdateExpansionCoefficients(analyzed_tree *T){
 
  SetHe3Pol();  // Sets the He3 polarization value/direction

  TLorentzVector Pe(0,0,T->ebeam,T->ebeam);   // incoming e-
  TLorentzVector Peprime(T->trPx,   // scattered e-
			 T->trPy,
			 T->trPz,
			 T->trP);
  
  // Calculate the scattering plane, which we will use to get angles
  TLorentzVector q = Pe - Peprime;
  TVector3 Pe_vect = Pe.Vect();
  TVector3 Peprime_vect = Peprime.Vect();
  TVector3 q_vect = q.Vect();
  TVector3 normal = Pe_vect.Cross(Peprime_vect);
  normal = normal.Unit();
  q_vect = q_vect.Unit();
  

  double m = constant::Mp;
  double tau = T->Q2 / (4 * m * m);   // Defenition of tau

  // These two give the polarization components 
  double P_x = normal.Dot(q_vect.Cross(TargetPolDirection)); 
  double P_z = q_vect.Dot(TargetPolDirection);
  
  // These variables are the expansion 
  double B = -2 * sqrt( tau * (1 + tau) ) * tan(T->etheta/2) * P_x;
  double C = -2 * tau * sqrt(1 + tau + pow((1 + tau) * tan(T->etheta/2),2) ) * tan(T->etheta/2) * P_z; 
  double D = tau + 2*tau * (1 + tau) * pow(tan(T->etheta/2),2);

  double T_0 = C / D;
  double T_1 = B / D;
  double T_2 = -C / (D*D);
  double T_3 = -B / (D*D);
  double T_4 = C / (D*D*D);
  double T_5 = B / (D*D*D);

  // Now we update the average for each value
  count_avg++;
  T_avg[0] += (T_0 - T_avg[0]) / count_avg;
  T_avg[1] += (T_1 - T_avg[1]) / count_avg;
  T_avg[2] += (T_2 - T_avg[2]) / count_avg;
  T_avg[3] += (T_3 - T_avg[3]) / count_avg;
  T_avg[4] += (T_4 - T_avg[4]) / count_avg;
  T_avg[5] += (T_5 - T_avg[5]) / count_avg;
  T_1_Q2_avg += (T_1*T->Q2 - T_1_Q2_avg) / count_avg;
  Q2_avg = T_1_Q2_avg / T_avg[1];

}


// This gets the GE/GM value from parameterizations
// This can be found in literature or Seamus/Freddy's theses 
double GetGEGM(bool is_neutron, double Q2){

  double m = constant::Mp;
  double tau = Q2 / (4 * m * m);
  double GD = pow(1.0 + Q2/(0.71), -2.0);  // Dipole approximation
  double GE,GM;

  if(is_neutron){ // Neutron parameterizations
    // Seamus Fit
    GE = (1.520*tau + 2.629*tau*tau + 3.055*tau*tau*tau)*GD/(1.0+5.222*tau+0.040*tau*tau+11.438*tau*tau*tau);
    // Kelly
    GM = -1.913*(1.0+2.33*tau)/(1.0 + 14.72*tau + 24.20*tau*tau + 84.1*tau*tau*tau );
  }
  else{ // Proton parameterizations
    // Kelly
    GE = (1.0-0.24*tau)/(1.0 + 10.98*tau + 12.82*tau*tau + 21.97*tau*tau*tau );
    // Kelly
    GM = 2.79*(1.0+0.12*tau)/(1.0 + 10.97*tau + 18.86*tau*tau + 6.55*tau*tau*tau );
  }

  return GE / GM;
}


double distribution_fits::fitbg_pol4( double *x, double *par ){
  double dx = x[0];

  // Let's use 4th order polynomial for the background:
  double bg = 0.0;
  for( int i = 0; i<5; i++ ){
    bg += par[i]*pow(dx,i);
  }

  return bg; 
}

double distribution_fits::fitbg_pol3( double *x, double *par ){
  double dx = x[0];

  // Let's use 3rd order polynomial for the background:
  double bg = 0.0;
  for( int i = 0; i<4; i++ ){
    bg += par[i]*pow(dx,i);
  }

  return bg; 
}

double distribution_fits::fitbg_gaus( double *x, double *par ){
  double dx = x[0];

  // Let's use gaussian for the background:
  double bg = par[0]*exp(-0.5*pow((dx-par[1])/par[2],2));

  return bg; 
}

// Custom fit for the He3 distributions. This assumes we have a proton,
// neutron, and background distribution. These data is then fitted to these
// three distributions
double distribution_fits::fitsim( double *x, double *par){
  double dx = x[0];

  double Norm_overall = par[0]; // Normalization for the spectrum
  double R_pn = par[1];         // Ratio between proton/neutron peak
  double Bg_norm = par[2];      // Normalization of the background

  double bg = 0.0;
  
  // Get the background function:
  if(bg_shape_option == "pol4") bg = fitbg_pol4(x,&par[3]);
  else if(bg_shape_option == "pol3") bg = fitbg_pol3(x,&par[3]);
  else if(bg_shape_option == "gaus") bg = fitbg_gaus(x,&par[3]);
  else if(bg_shape_option == "from data") bg = hdx_bg_data->Interpolate(dx);
  
  // The fit is now the sum of these three distributions
  // fit = N * (p_dist + R * n_dist + N_bg * bg_dist)
  double simu = Norm_overall * (hdx_sim_p->Interpolate(dx) + R_pn * hdx_sim_n->Interpolate(dx) + Bg_norm * bg);
  return simu;   
}


// This function is used to fit the He3 data
void distribution_fits::He3_fit_dists(){

  ////////////First do some checks to make sure that things are set///////
  ////////////////////////////////////////////////////////////////////////
  // Check bg shape is set
  if(!bg_shape_set){
    cout<<"Error: [distribution_fits::He3_sim_fit] bg shape has not been set!!!"<<endl;
    exit(0);
  }
  else {
    if(bg_shape_option != "pol3" && bg_shape_option != "pol4" && bg_shape_option != "gaus" && bg_shape_option != "from data"){
      cout<<"Error: [distribution_fits::He3_sim_fit] bg shape option " + bg_shape_option + " is not supported!!!"<<endl;
      exit(0);
    }
  }

  //Check if all the histograms are set
  if(hdx_data == NULL){
    cout<<"Error: [distribution_fits::He3_sim_fit] hdx_data has not been set!!!"<<endl;
    exit(0);
  }
  if(hdx_sim_p == NULL){
    cout<<"Error: [distribution_fits::He3_sim_fit] hdx_sim_p has not been set!!!"<<endl;
    exit(0);
  }
  if(hdx_sim_n == NULL){
    cout<<"Error: [distribution_fits::He3_sim_fit] hdx_sim_n has not been set!!!"<<endl;
    exit(0);
  }
  if(hdx_bg_data == NULL){
    cout<<"Error: [distribution_fits::He3_sim_fit] hdx_bg_data has not been set!!!"<<endl;
    exit(0);
  }

  ////////////////////////////////////////////////////////////////////////
  
  //Normalize the histograms so they are all of a similar scale
  double scale = hdx_data->Integral();
  hdx_data->Scale(1.0/hdx_data->Integral());
  hdx_sim_p->Scale(1.0/hdx_sim_p->Integral());
  hdx_sim_n->Scale(1.0/hdx_sim_n->Integral());
  hdx_bg_data->Scale(1.0/hdx_bg_data->Integral());

  //Get histogram data so all future histograms match
  int nbins = hdx_data->GetNbinsX();
  double xmin = hdx_data->GetXaxis()->GetBinLowEdge(1);
  double xmax = hdx_data->GetXaxis()->GetBinUpEdge(nbins);
  
  //Set background type and check if the option exists
  int npar = -1;
  if(bg_shape_option == "pol4") npar = 3 + 5;
  else if(bg_shape_option == "pol3") npar = 3 + 4;
  else if(bg_shape_option == "gaus") npar = 3 + 3;
  else if(bg_shape_option == "from data") npar = 3;  //Our only bg par is the scale of the bg

  //Set fit to function fitsim
  TF1 *FitFunc = new TF1( "FitFunc", this, &distribution_fits::fitsim,xmin,xmax,npar);
      
  //Set some arbitrary starting values
  FitFunc->SetNpx(1000);
  double startpar[] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  FitFunc->SetParameters(startpar);
  FitFunc->SetParLimits(0,0.1,100);
  FitFunc->SetParLimits(1,0.1,100);
  FitFunc->SetParLimits(2,0,100);
    
  hdx_data->Fit(FitFunc,"q0","",xmin,xmax);
    
  //Scale the p/n functions
  //This assumes our fit was like N*(p_dist + R*n_dist + N_bg * bg_dist)
  hdx_sim_p->Scale( FitFunc->GetParameter(0) );
  hdx_sim_n->Scale( FitFunc->GetParameter(0)*FitFunc->GetParameter(1) );

  hdx_bg_fit = new TH1F("hdx_bg","",nbins,xmin,xmax);
  hdx_total_fit = new TH1F("hdx_total_fit","",nbins,xmin,xmax);

  //Set the backgrond function using the fit we already found
  TF1 *fit_bg;

  if(bg_shape_option == "pol4") 
    fit_bg = new TF1("fit_bg",this, &distribution_fits::fitbg_pol4,xmin,xmax,npar - 3);
  else if(bg_shape_option == "gaus") 
    fit_bg = new TF1("fit_bg",this, &distribution_fits::fitbg_gaus,xmin,xmax,npar - 3);
  else
    fit_bg = new TF1("fit_bg",this, &distribution_fits::fitbg_pol3,xmin,xmax,npar - 3);

  for(int i=0; i<npar - 3; i++)
    fit_bg->SetParameter(i,FitFunc->GetParameter(i+3));
    
  //Fill the bg hist and total hist
  for(int ibin = 0; ibin < nbins;ibin++){
    //Again normalize bg assuming our fit is N*(p_dist + R*n_dist + N_bg * bg_dist)

    if(bg_shape_option == "from data"){
      hdx_bg_fit->SetBinContent(ibin,FitFunc->GetParameter(0) * FitFunc->GetParameter(2) * hdx_bg_data->GetBinContent(ibin));
    }
    else {
      hdx_bg_fit->SetBinContent(ibin,FitFunc->GetParameter(0) * FitFunc->GetParameter(2) * fit_bg->Eval(hdx_total_fit->GetBinCenter(ibin)));
    }
      
	
    hdx_total_fit->SetBinContent(ibin,hdx_sim_n->GetBinContent(ibin) + hdx_sim_p->GetBinContent(ibin) + hdx_bg_fit->GetBinContent(ibin));
  }


  // Last step is to scale everything back up to the original data scale
  hdx_data->Scale(scale);
  hdx_sim_p->Scale(scale);
  hdx_sim_n->Scale(scale);
  hdx_bg_fit->Scale(scale);
  hdx_total_fit->Scale(scale);

}


// This funciton is out dated. It used Seamus' method of dilutions which I am
// no longer following
void analyzed_info::CalcAsymValsOld(){

  // Make sure all the variables have been loaded
  if(N_raw_p == 0 && N_raw_m == 0){
    cout<<"Error: [analyzed_info::CalcAsymVals] N_raw_p and N_raw_n are both zero!!!"<<endl;
    exit(0);
  }
  if(N_bg_p == 0 && N_bg_m == 0){
    cout<<"Error: [analyzed_info::CalcAsymVals] N_bg_p and N_bg_n are both zero!!!"<<endl;
    exit(0);
  }
  if(N_proton == 0){
    cout<<"Error: [analyzed_info::CalcAsymVals] N_proton has not been set!!!"<<endl;
    exit(0);
  }
 
  double N_N2 = 0.01;  //Currently set these but should be calcuated later
  double N_ot = 0.01;  //Currently set these but should be calcuated later

  // Get total numbers
  double Sigma_raw = N_raw_p + N_raw_m;
  double Sigma_bg =  N_bg_p + N_bg_m;
  double Sigma_p =   N_proton;
  double Sigma_N2 =  N_N2;
  double Sigma_in =  N_in_p + N_in_m;
  double Sigma_in_cont =  N_in_cont;

  // Get differences in helicities
  double Delta_raw = N_raw_p - N_raw_m;
  double Delta_bg = N_bg_p - N_bg_m; 
  double Delta_in = N_in_p - N_in_m; 

  // Calculate the dilutions, formula's from Seamus' thesis
  D_bg = 1 - Sigma_bg / Sigma_raw;
  D_N2 = 1 - Sigma_N2 / (Sigma_raw - Sigma_bg);
  D_p = 1 - Sigma_p / (Sigma_raw - Sigma_N2 - Sigma_bg);
  D_in = 1 - Sigma_in_cont / (Sigma_raw - Sigma_p - Sigma_N2 - Sigma_bg);
  
  // Calculate the asymmetry values
  A_raw = Delta_raw / Sigma_raw;
  A_bg = Delta_bg / Sigma_bg;
  A_in = Delta_in / Sigma_in;

  // Calculate proton asymmetry from the expansion
  A_p = 0;
  for(int i=0; i < nexp; i++)
    A_p += T_avg[i] * pow(GetGEGM(0,Q2_avg),i);
  
  A_p *= (1 - D_p) / (D_bg*D_N2)*P_He3*P_beam*P_p;
  
  // Calculate the total physics asymmetry
  double denom = P_He3 * P_n * P_beam * D_bg * D_N2 * D_p * D_in;
  A_phys = ( A_raw - A_bg - A_p - A_in ) / denom;

  // Calculate the errors
  A_raw_err = 1.0 / sqrt(Sigma_raw);
  A_bg_err = sqrt( Sigma_bg / (2*Sigma_raw*Sigma_raw) + Delta_bg*Delta_bg / (12*Sigma_raw*Sigma_raw) + Delta_bg*Delta_bg / (4*Sigma_raw*Sigma_raw*Sigma_raw) );
  A_p_err = 0;   // Set to 0 for now


  D_bg_err = sqrt( Sigma_bg / (2*Sigma_raw*Sigma_raw) + Sigma_bg*Sigma_bg / (12*Sigma_raw*Sigma_raw) + Sigma_bg*Sigma_bg / (4*Sigma_raw*Sigma_raw*Sigma_raw) );   // Set to 0 for now
  D_N2_err = 0;   // Set to 0 for now
  D_p_err = 0;    // Set to 0 for now
  D_in_err = 0;   // Set to 0 for now
  
  double D_errors = pow(D_bg_err/D_bg,2) + pow(D_N2_err/D_N2,2) + pow(D_p_err/D_p,2) + pow(D_in_err/D_in,2);

  A_phys_stat_err = A_raw_err / denom;
  A_phys_sys_err = sqrt( ( A_p_err*A_p_err + A_bg_err*A_bg_err) / (denom*denom) + A_phys*A_phys*D_errors );
  
}

// This function calculates the asymmetry values. This formalism follows a more
// standard formalism of calculating each asymmetry contribution and scaling 
// it by the fraction inside the quasielastic cuts
void analyzed_info::CalcAsymVals(){

  ////////////////// Make sure all the variables have been loaded /////////////
  if(N_raw_p == 0 && N_raw_m == 0){
    cout<<"Error: [analyzed_info::CalcAsymVals] N_raw_p and N_raw_n are both zero!!!"<<endl;
    exit(0);
  }
  if(N_bg_p == 0 && N_bg_m == 0){
    cout<<"Warning: [analyzed_info::CalcAsymVals] N_bg_p and N_bg_n are both zero!!!"<<endl;
  }
  if(N_pion_p == 0 && N_pion_m == 0){
    cout<<"Warning: [analyzed_info::CalcAsymVals] N_pion_p and N_pion_n are both zero!!!"<<endl;
  }
  if(N_in_p == 0 && N_in_m == 0){
    cout<<"Warning: [analyzed_info::CalcAsymVals] N_in_p and N_in_n are both zero!!!"<<endl;
  }
  if(N_in_cont == 0){
    cout<<"Warning: [analyzed_info::CalcAsymVals] N_in_cont is zero!!!"<<endl;
  }
  if(N_proton == 0){
    cout<<"Warning: [analyzed_info::CalcAsymVals] N_proton is zero!!!"<<endl;
  }
  if(N_FSI_p == 0 && N_FSI_m == 0){
    cout<<"Warning: [analyzed_info::CalcAsymVals] N_FSI_p and N_FSI_n are both zero!!!"<<endl;
  }
  /////////////////////////////////////////////////////////////////////////


  int N_N2 = 0;  //Currently set these but should be calcuated later

  // Get the total values for each contribution
  double Sigma_raw = N_raw_p + N_raw_m;
  double Sigma_bg =  N_bg_p + N_bg_m;
  double Sigma_p =   N_proton;
  double Sigma_N2 =  N_N2;
  double Sigma_pion =  N_pion_p + N_pion_m;
  double Sigma_in =  N_in_p + N_in_m;
  double Sigma_FSI =  N_FSI_p + N_FSI_m;

  // Get the difference in helicity for each contribution
  double Delta_raw = N_raw_p - N_raw_m;
  double Delta_bg = N_bg_p - N_bg_m; 
  double Delta_pion = N_pion_p - N_pion_m; 
  double Delta_in = N_in_p - N_in_m; 
  double Delta_FSI = N_FSI_p - N_FSI_m; 

  // Calculate the fractions
  f_bg = Sigma_bg / Sigma_raw;
  f_N2 = Sigma_N2 / Sigma_raw;
  f_p = Sigma_p / Sigma_raw;
  f_pion = Sigma_pion / Sigma_raw;
  f_in = N_in_cont / Sigma_raw;
  f_FSI = Sigma_FSI / Sigma_raw;
  f_n = 1 - f_bg - f_N2 - f_p - f_pion - f_in - f_FSI;
  
  // Calculate the asymmetry values
  A_raw = Delta_raw / Sigma_raw;
  A_bg = Delta_bg / Sigma_bg;
  A_pion = Delta_pion / Sigma_pion;
  A_in = Delta_in / Sigma_in;
  A_FSI = Delta_FSI / Sigma_FSI;
  
  // Calculate the proton asymmetry from the expansion
  A_p = 0;
  for(int i=0; i < nexp; i++)
    A_p += T_avg[i] * pow(GetGEGM(0,Q2_avg),i);
  
  A_p *= P_He3*P_beam*P_p;
  
  // Calculate the total physics asymmetry
  double denom = P_He3 * P_n * P_beam * f_n;
  A_phys = ( A_raw - f_bg*A_bg - f_p*A_p - f_pion*A_pion - f_in*A_in - f_FSI*A_FSI ) / denom;

  // Calculate the errors
  A_raw_err = 1.0 / sqrt(Sigma_raw);
  A_bg_err = 0;  // Set to 0 for now
  A_p_err = 0;   // Set to 0 for now

  f_bg_err = 0;    // Set to 0 for now
  f_N2_err = 0;    // Set to 0 for now
  f_p_err = 0;     // Set to 0 for now
  f_pion_err = 0;  // Set to 0 for now
  f_in_err = 0;    // Set to 0 for now
  f_FSI_err = 0;   // Set to 0 for now
  
  // Separate out the statistical and systematic errors
  A_phys_stat_err = A_raw_err / denom;
  A_phys_sys_err = 0;
  
}


