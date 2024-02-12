#include <iostream>

#include "../include/Analysis.h"

namespace Analysis {

  double fitbg_pol4( double *x, double *par ){
    double dx = x[0];

    // Let's use (up to) 4th order polynomial for the background:
    double bg = 0.0;
    for( int i = 0; i<5; i++ ){
      bg += par[i]*pow(dx,i);
    }

    return bg; 
  }

  double fitbg_pol3( double *x, double *par ){
    double dx = x[0];

    // Let's use (up to) 3rd order polynomial for the background:
    double bg = 0.0;
    for( int i = 0; i<4; i++ ){
      bg += par[i]*pow(dx,i);
    }

    return bg; 
  }

  double fitbg_gaus( double *x, double *par ){
    double dx = x[0];

    // Let's use gaussian for the background:
    double bg = par[0]*exp(-0.5*pow((dx-par[1])/par[2],2));

    return bg; 
  }

  double fitsim( double *x, double *par ){
    double dx = x[0];

    double Norm_overall = par[0];
    double R_pn = par[1];
    double Bg_norm = par[2];

    double bg = 0.0;

    // Get the background function:
    if(bg_option == "pol4") bg = fitbg_pol4(x,&par[3]);
    else if(bg_option == "pol3") bg = fitbg_pol3(x,&par[3]);
    else if(bg_option == "gaus") bg = fitbg_gaus(x,&par[3]);
    else if(bg_option == "from data") bg = hdx_bg_data->Interpolate(dx);
    else{
      cout<<"[Analysis:fitsim] bg fit option is not valid!!!"<<endl;
      exit(0);
    }
  

    double simu = Norm_overall * (hdx_sim_p->Interpolate(dx) + R_pn * hdx_sim_n->Interpolate(dx) + Bg_norm * bg);

    return simu;   
  }
  
  void He3_fit(TH1F *hdx,TString config, TF1 **fit_bg,TF1 **fit_n,TF1 **fit_p,TF1 **fit_total){
    double bg_low;
    double bg_high;

    double p_low;
    double p_high;

    double n_low;
    double n_high;

    //These parameters are determined by looking at the peaks on the plots by eye
    if(config == "GEN2"){
      p_low = -3.5;
      p_high = -2.2;
    
      n_low = -0.8;
      n_high = 0.5;
    
      bg_low = -7.0;
      bg_high = 3.0;
    }

    if(config == "GEN3"){
      p_low = -2.0;
      p_high = -0.8;
    
      n_low = -0.4;
      n_high = 0.4;
    
      bg_low = -4.5;
      bg_high = 2.5;
    }

    if(config == "GEN4" || config == "GEN4b"){
      p_low = -1.4;
      p_high = -0.9;
    
      n_low = -0.3;
      n_high = 0.1;
    
      bg_low = -4.0;
      bg_high = 1.5;
    }

    //total function = proton gaus + neutron gaus + 4th order bkgd
    double par[11];
    TF1 *p_xfunc = new TF1("p_xfunc","gaus",p_low,p_high);
    TF1 *n_xfunc = new TF1("n_xfunc","gaus",n_low,n_high);  
    TF1 *bg_xfunc = new TF1("bg_xfunc","pol4",bg_low,bg_high);
    TF1 *total_xfunc = new TF1("total_xfunc","gaus(0) + gaus(3) + pol4(6)",bg_low,bg_high);
  

    //Do fit but do not plot the results
    hdx->Fit(p_xfunc,"qNR");  
    hdx->Fit(n_xfunc,"qNR+");  
    hdx->Fit(bg_xfunc,"qNR+"); 

    //Put the fit parameters into the array
    p_xfunc->GetParameters(&par[0]);
    n_xfunc->GetParameters(&par[3]);
    bg_xfunc->GetParameters(&par[6]);

    //Set the parameters in the total function using the results above
    total_xfunc->SetParameters(par);
    if(config == "GEN3") {
      total_xfunc->SetParLimits(4,-0.2,0.2);
      total_xfunc->SetParLimits(5,0.2,0.6);
    }
    if(config == "GEN4") {
      total_xfunc->SetParLimits(4,-0.2,0.1);
      total_xfunc->SetParLimits(5,0.2,0.6);
    }
    hdx->Fit(total_xfunc,"qNR+"); 
    //Get the fit results
    total_xfunc->GetParameters(&par[0]);

    //For plotting purposes set the p/n function parameters from the total fit
    p_xfunc = new TF1("p_xfunc","gaus",-5,5);
    p_xfunc->SetParameters(&par[0]);

    n_xfunc = new TF1("n_xfunc","gaus",-5,5);
    n_xfunc->SetParameters(&par[3]);

    bg_xfunc = new TF1("bg_xfunc","pol4",bg_low,bg_high);
    bg_xfunc->SetParameters(&par[6]);

    cout<<total_xfunc->GetParameter(4)<<" "<<total_xfunc->GetParameter(5)<<endl;

    //Save the fit result for use later
    *fit_bg = bg_xfunc;
    *fit_n = n_xfunc;
    *fit_p = p_xfunc;
    *fit_total = total_xfunc;

  }
  
  void He3_sim_fit(TH1F *hdx_data){

    //Get histogram data so all future histograms match
    int nbins = hdx_data->GetNbinsX();
    double xmin = hdx_data->GetXaxis()->GetBinLowEdge(1);
    double xmax = hdx_data->GetXaxis()->GetBinUpEdge(nbins);

    //Set background type and check if the option exists
    int npar = -1;
    if(bg_option == "pol4") npar = 3 + 5;
    else if(bg_option == "pol3") npar = 3 + 4;
    else if(bg_option == "gaus") npar = 3 + 3;
    else if(bg_option == "from data"){
      if(hdx_bg_data == nullptr){
	cout<<"[Analysis:He3_sim_fit] bg_data histogram is empty!!!"<<endl;
	exit(0);
      }
      npar = 3;  //Our only bg par is the scale of the bg
    }
    else {
      cout<<"[Analysis:He3_sim_fit] bgfit is not a valid option!!!"<<endl;
      exit(0);
    }

    //Set fit to function fitsim
    TF1 *FitFunc = new TF1( "FitFunc", fitsim,xmin,xmax,npar);
    
    //Set some arbitrary starting values, should not be hardcoding this
    FitFunc->SetNpx(1000);
    double startpar[] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    FitFunc->SetParameters(startpar);
    FitFunc->SetParLimits(0,0.1,100);
    FitFunc->SetParLimits(1,0.1,100);
    FitFunc->SetParLimits(2,0,100);
    
    
    hdx_data->Fit(FitFunc,"0","",xmin,xmax);
    
    //Scale the p/n functions
    //This assumes our fit was like N*(p_dist + R*n_dist + N_bg * bg_dist)
    hdx_sim_p->Scale( FitFunc->GetParameter(0) );
    hdx_sim_n->Scale( FitFunc->GetParameter(0)*FitFunc->GetParameter(1) );

    hdx_bg = new TH1F("hdx_bg","",nbins,xmin,xmax);
    hdx_total_fit = new TH1F("hdx_total_fit","",nbins,xmin,xmax);

    //Set the backgrond function using the fit we already found
    TF1 *fit_bg;

    if(bg_option == "pol4") 
      fit_bg = new TF1("fit_bg",fitbg_pol4,xmin,xmax,npar - 3);
    else if(bg_option == "gaus") 
      fit_bg = new TF1("fit_bg",fitbg_gaus,xmin,xmax,npar - 3);
    else
      fit_bg = new TF1("fit_bg",fitbg_pol3,xmin,xmax,npar - 3);

    for(int i=0; i<npar - 3; i++)
      fit_bg->SetParameter(i,FitFunc->GetParameter(i+3));
    
    //Fill the bg hist and total hist
    for(int ibin = 0; ibin < nbins;ibin++){
      //Again normalize bg assuming our fit is N*(p_dist + R*n_dist + N_bg * bg_dist)

      if(bg_option == "from data"){
	hdx_bg->SetBinContent(ibin,FitFunc->GetParameter(0) * FitFunc->GetParameter(2) * hdx_bg_data->GetBinContent(ibin));
      }
      else {
	hdx_bg->SetBinContent(ibin,FitFunc->GetParameter(0) * FitFunc->GetParameter(2) * fit_bg->Eval(hdx_total_fit->GetBinCenter(ibin)));
      }
      
	
      hdx_total_fit->SetBinContent(ibin,hdx_sim_n->GetBinContent(ibin) + hdx_sim_p->GetBinContent(ibin) + hdx_bg->GetBinContent(ibin));
    }

  }
 
}
