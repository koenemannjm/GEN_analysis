#ifndef Analysis_H
#define Analysis_H


namespace Analysis {

  void He3_fit(TH1F *hdx,TString config, TF1 **fit_bg,TF1 **fit_n,TF1 **fit_p);
  void He3_sim_fit(TH1F *hdx_data);

  TH1F *hdx_sim_p;
  TH1F *hdx_sim_n;
  TH1F *hdx_bg;
  TH1F *hdx_bg_data;
  TH1F *hdx_total_fit;

  TString bg_option;
  double fitbg_pol3(double *x, double *par);
  double fitbg_pol4(double *x, double *par);
  double fitbg_gaus(double *x, double *par);
  double fitsim(double *x, double *par);

}

#endif
