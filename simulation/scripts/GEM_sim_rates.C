


void GEM_sim_rates(int Q2 = 29){


  TFile *f = new TFile(Form("Rootfiles/GEN_background_rates_%i_GeV.root",Q2),"read");

  TH1D *hbb_gem_rate = (TH1D*)f->Get("hitrate_vs_layer_BBGEM");
  TH1D *hsbs_gem_rate = (TH1D*)f->Get("hitrate_vs_layer_SBSGEM");

  
  hbb_gem_rate->SetTitle("BigBite GEM rates,  60 cm He3 target;Layer;Hit Rate (Hz/cm^{2}/#mu A)");
  hbb_gem_rate->GetYaxis()->SetTitleOffset(1.4); 

  hsbs_gem_rate->SetTitle("Super BigBite GEM rates,  60 cm He3 target;Layer;Hit Rate (Hz/cm^{2}/#mu A)");

  TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  hbb_gem_rate->Draw(); 

  TCanvas *c2 = new TCanvas("c2","c2",1200,900);
  hsbs_gem_rate->Draw(); 


  TString output = "../plots/";
  double Q2_double = Q2*1.0/10;


  c1->Print(output + Form("GEM_rate_%gGeV.pdf(",Q2_double));
  c2->Print(output + Form("GEM_rate_%gGeV.pdf)",Q2_double));
      
}
