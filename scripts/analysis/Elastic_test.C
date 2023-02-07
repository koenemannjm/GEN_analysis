
double bg_low;
double bg_high;

double p_low;
double p_high;

double n_low;
double n_high;

double bg_fit(double *x, double *par){

  if(x[0] > bg_low && x[0] < bg_high)
    TF1::RejectPoint();

  return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2));
}

void get_np_spots(TCanvas *c, TString cfg, TFile *file){
  
  TH2D *hdxdy_Wcut = (TH2D*)file->Get("hdxdy_Wcut");

  TH1D *hdx = hdxdy_Wcut->ProjectionY();

  p_low = -4.0;
  p_high = -1.5;

  n_low = -1.5;
  n_high = 2;
  
  bg_low = -4.0;
  bg_high = 3.0;

  double pF_par0 = 9000;
  double pF_par1 = -2.4;
  double pF_par2 = 0.3;
  
  double par[11];
  //TF1 *p_xfunc = new TF1("p_xfunc","gaus",p_low,p_high);
  TF1 *p_xfunc = new TF1("p_xfunc","gaus",p_low,p_high);
  TF1 *n_xfunc = new TF1("n_xfunc","gaus",n_low,n_high);  
  TF1 *bg_xfunc = new TF1("bg_xfunc","pol4",bg_low,bg_high);
  TF1 *total_xfunc = new TF1("total_xfunc","gaus(0) + gaus(3) + pol4(6)",-4,3);
  

  hdx->Fit(p_xfunc,"NR");  
  hdx->Fit(n_xfunc,"NR+");  
  hdx->Fit(bg_xfunc,"NR+"); 

  p_xfunc->GetParameters(&par[0]);
  n_xfunc->GetParameters(&par[3]);
  bg_xfunc->GetParameters(&par[6]);

  total_xfunc->SetParameters(par);
  hdx->Fit(total_xfunc,"R+"); 
  total_xfunc->GetParameters(&par[0]);


  double px_mean = par[1];
  double px_sigma = par[2];
  double nx_mean = par[4];
  double nx_sigma = par[5];

  cout<<px_mean<<" "<<px_sigma<<endl;
  cout<<nx_mean<<" "<<nx_sigma<<endl;
  
  p_xfunc = new TF1("p_xfunc","gaus",-4,4);
  p_xfunc->SetParameters(&par[0]);

  n_xfunc = new TF1("n_xfunc","gaus",-4,4);
  n_xfunc->SetParameters(&par[3]);

  
  
 
  double x_cut = -0.4;
  int x_cut_bin = hdxdy_Wcut->GetYaxis()->FindBin(x_cut);

  TH1D *hdy_p = hdxdy_Wcut->ProjectionX("dy_p",0,x_cut_bin);
  TH1D *hdy_n = hdxdy_Wcut->ProjectionX("dy_n",x_cut_bin,-1);

  

  /*

  TF1 *p_yfunc = new TF1("p_yfunc","gaus",-1.5,1.5);
  hdy_p->Fit(p_yfunc,"qR+");  
  p_yfunc->GetParameters(&par[0]);

  double py_mean = par[1];
  double py_sigma = par[2];

  TF1 *n_yfunc = new TF1("n_yfunc","gaus",-1.5,1.5);
  hdy_n->Fit(n_yfunc,"qR+");  
  n_yfunc->GetParameters(&par[0]);

  double ny_mean = par[1];
  double ny_sigma = par[2];

  double nsigma = 1;

  if(cfg == "H2") nsigma = 1;

  TEllipse *p_spot = new TEllipse(py_mean,px_mean,nsigma*py_sigma,nsigma*px_sigma);
  TEllipse *n_spot = new TEllipse(ny_mean,nx_mean,nsigma*ny_sigma,nsigma*nx_sigma);
  p_spot->SetFillStyle(0);
  n_spot->SetFillStyle(0);

  p_spot->SetLineWidth(4);
  n_spot->SetLineWidth(4);

  p_spot->SetLineColor(kRed);
  n_spot->SetLineColor(kRed);

  c->cd();

  gStyle->SetOptStat(0);

  hdxdy_Wcut->SetTitle(cfg + " HCal Elastics;#Deltay (m);#Deltax (m)");
  hdxdy_Wcut->Draw("colz");
  p_spot->Draw("same");
  if(cfg == "He3") n_spot->Draw("same");


  TPaveText *pt = new TPaveText(.65,.8,.88,.88,"ndc");
  pt->AddText("Cuts on good tracks");
  pt->AddText("|W^{2} - 0.88| < 0.4");
  pt->SetFillColor(0);
  pt->Draw("same");

  cout<<cfg<<" data:"<<endl;
  cout<<"proton dx mean = "<<px_mean<<endl;
  cout<<"proton dx sigma = "<<px_sigma<<endl;
  cout<<"neutron dx mean = "<<nx_mean<<endl;
  cout<<"neutron dx sigma = "<<nx_sigma<<endl;

  cout<<"proton dy mean = "<<py_mean<<endl;
  cout<<"proton dy sigma = "<<py_sigma<<endl;
  cout<<"neutron dy mean = "<<ny_mean<<endl;
  cout<<"neutron dy sigma = "<<ny_sigma<<endl;
  cout<<"\n\n";

  */
}

void Elastic_test(){


  TFile *He3_file = new TFile("outfiles/QE_test_GEN2_sbs100p_model2_data.root","read");
  
  TCanvas *c1 = new TCanvas("c1","",800,1000);  
  get_np_spots(c1, "He3", He3_file);

  TH2D *hdxdy_He3 = (TH2D*)He3_file->Get("h2_dxdyHCAL");
  TH1D *hdx_He3 = hdxdy_He3->ProjectionY("hdx_He3");
  hdx_He3->SetLineColor(kRed);
  hdx_He3->Scale(1/hdx_He3->GetEntries());

  //hdx_He3->Draw("hist");
}
