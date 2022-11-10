


void EArm_sim_rates(int Q2 = 29){

  const int ngen = 5;
  TString gen_name[ngen] = {"elastic","inelastic","wiser_pi0","wiser_pip","wiser_pim"};
  TString Kine[ngen] = {"Elastic","Inelastic","Pi0","Pi+","Pi-"};

  double right_bound = 4.0;  
  if(Q2 == 66) right_bound = 5.0;
  if(Q2 == 97) right_bound = 5.0;

  double Q2_double = Q2*1.0/10;

  TFile *f = new TFile(Form("Rootfiles/EArm_sim_rates_%i_GeV.root",Q2),"read");

  TH1D *hrate_good_elastics = (TH1D*)f->Get("hrate_good_elastics");

  TCanvas *c1 = new TCanvas("c1","c1",1200,900);

  hrate_good_elastics->SetTitle(Form("GEN  Q^{2}=%gGeV^{2}  Threshold Estimation for BigBite Rates, Good Elastics",Q2_double));
  hrate_good_elastics->Draw();

  gStyle->SetOptStat(0);

  hrate_good_elastics->Fit("gaus","q","",0.0,right_bound);
  TF1 *fR = hrate_good_elastics->GetFunction("gaus");
  fR->SetLineColor(4);

  double thresh[3] = {fR->GetParameter(1) - 2.0*fR->GetParameter(2), fR->GetParameter(1) - 2.5*fR->GetParameter(2), fR->GetParameter(1) - 3.0*fR->GetParameter(2)};


  TLegend *l1 = new TLegend(0.11,0.8,0.5,0.61);  //(0.22,0.7,0.51,0.89); 
  l1->AddEntry((TObject*)0, Form("Threshold @ 3#sigma = %.3fGeV",thresh[2]), "");
  l1->AddEntry((TObject*)0, Form("Threshold @ 2.5#sigma = %.3fGeV",thresh[1]), "");
  l1->AddEntry((TObject*)0, Form("Threshold @ 2#sigma = %.3fGeV",thresh[0]), "");  
  l1->AddEntry(hrate_good_elastics,"Smeared Energy","ep");
  l1->AddEntry(fR,"Gaussian Fit","l");
  l1->Draw();


  TH1D *hrate_gen[ngen];
  TCanvas *c[ngen];

  for(int igen = 0; igen < ngen; igen++){

    c[igen] = new TCanvas(Form("c_%i",igen),"",1200,900);

    if(gen_name[igen] != "elastic") c[igen]->SetLogy();

    hrate_gen[igen] = (TH1D*)f->Get("hrate_" + gen_name[igen]);

    hrate_gen[igen]->SetTitle(Form("GEN  Q^{2}=%gGeV^{2}  Kine=" + Kine[igen] +" | Trig. Log.",Q2_double));
    hrate_gen[igen]->SetMarkerStyle(20);
    hrate_gen[igen]->SetMarkerColor(2);
    hrate_gen[igen]->SetLineColor(2);
    hrate_gen[igen]->GetXaxis()->SetTitle("Energy sum (PS+SH) in GeV");
    hrate_gen[igen]->GetYaxis()->SetTitle("Rate/bin in Hz");
    hrate_gen[igen]->Draw();


    TPaveText *pt = new TPaveText(0.55,0.68,0.9,0.9,"ndc");
    pt->AddText("Beam Current = 60uA");
    pt->AddText("Target Length = 60cm"); 
    pt->AddText(Form("Thresh@3#sigma = %.2gGeV | Trig rate = %.2eHz",thresh[2],hrate_gen[igen]->Integral(hrate_gen[igen]->FindFixBin(thresh[2]),hrate_gen[igen]->FindFixBin(right_bound))));    
    pt->AddText(Form("Thresh@2.5#sigma = %.2gGeV | Trig rate = %.2eHz",thresh[1],hrate_gen[igen]->Integral(hrate_gen[igen]->FindFixBin(thresh[1]),hrate_gen[igen]->FindFixBin(right_bound))));
    pt->AddText(Form("Thresh@2#sigma = %.2gGeV | Trig rate = %.2eHz",thresh[0],hrate_gen[igen]->Integral(hrate_gen[igen]->FindFixBin(thresh[0]),hrate_gen[igen]->FindFixBin(right_bound))));
    
    pt->SetFillColor(0);
    pt->SetBorderSize(1);
    pt->Draw();

  }

  TString output = "../plots/";

  c1->Print(output + Form("EArm_rate_%gGeV.pdf(",Q2_double));
  
  for( int igen = 0; igen < ngen; igen++){
    if(igen == ngen - 1) c[igen]->Print(output + Form("EArm_rate_%gGeV.pdf)",Q2_double));
    else c[igen]->Print(output + Form("EArm_rate_%gGeV.pdf",Q2_double));
  }
  
}
