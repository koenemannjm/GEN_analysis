


#include "../../include/gen-ana.h"

double getGEGM(double t, double e, double Px, double Pz, double A_phys){

  double A = e / t * A_phys;
  double B = Px*sqrt(2*e*(1 - e) / t);
  double C = A_phys + Pz*sqrt(1 - e*e);

  // We are solving the quadratic equation to get GE/GM
  double R1 = (-B + sqrt(B*B - 4*A*C)) / (2*A);
  double R2 = (-B - sqrt(B*B - 4*A*C)) / (2*A);
  GEGM = R1;  // I think this root is correct

  return GEGM;

}


void test_sim(){


  TString rootdir = "/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/simulation/GEN2/TEST/";

  TString filename = "replayed_GEN2_He3_elastic_job_0.root";

  TFile *file = new TFile(rootdir + filename,"read");
  TTree *T = (TTree*)file->Get("T");
  
  T->SetBranchStatus("*",0);

  //MC variables
  double mc_sigma, mc_omega, mc_fnucl, mc_simc_weight, mc_Q2, mc_Pl, mc_Pt, mc_epx, mc_epz, mc_Apar, mc_Aperp;
  std::vector<std::string> mc = {"mc_sigma","mc_omega","mc_fnucl","simc_Weight","mc_Q2","mc_Pl","mc_Pt","mc_epx","mc_epz","mc_Apar", "mc_Aperp"};
  std::vector<void*> mc_mem = {&mc_sigma,&mc_omega,&mc_fnucl,&mc_simc_weight,&mc_Q2,&mc_Pl,&mc_Pt,&mc_epx,&mc_epz,&mc_Apar,&mc_Aperp};
  setrootvar::setbranch(T,"MC",mc,mc_mem);

  int nevent = 0;
  double m = constant::Mp;
  int count_avg = 0;
  double sim_tau_avg = 0, sim_tan_avg = 0, sim_eps_avg = 0, sim_Px_avg = 0, sim_Pz_avg = 0, sim_Aphys_avg = 0, GEGM_avg = 0;

  while(T->GetEntry(nevent++)){
    
    if(mc_fnucl == 1) continue;
    
    double theta = atan(mc_epx/mc_epz);
    double tan_th = tan(theta / 2);
    double tau = mc_Q2 / (4 * m * m);
    double e =  1.0 / (1 + 2*(1 + tau)*pow(tan_th,2));
    double A_phys = mc_Apar*mc_Pl + mc_Aperp*mc_Pt;
    double GEGM = getGEGM(tau, e, mc_Pt, mc_Pl, A_phys);
    
    count_avg++;
    sim_tau_avg += (tau - sim_tau_avg) / count_avg;
    sim_tan_avg += (tan_th - sim_tan_avg) / count_avg;
    sim_eps_avg += (e - sim_eps_avg) / count_avg;
    sim_Px_avg += (mc_Pt - sim_Px_avg) / count_avg;
    sim_Pz_avg += (mc_Pl - sim_Pz_avg) / count_avg;
    sim_Aphys_avg += (A_phys - sim_Aphys_avg) / count_avg;
    GEGM_avg += (GEGM - GEGM_avg) / count_avg;

    double GEGM_test = getGEGM(sim_tau_avg, sim_eps_avg, sim_Px_avg, sim_Pz_avg, sim_Aphys_avg);
    //cout<<GEGM_avg<<" "<<GEGM_test<<" "<<abs((GEGM_avg - GEGM_test) / GEGM_avg)<<endl;
  }

  double GEGM_kin = getGEGM(sim_tau_avg, sim_eps_avg, sim_Px_avg, sim_Pz_avg, sim_Aphys_avg);
  
  cout<<count_avg<<endl;
  cout<<"Averaged GEGM: "<<GEGM_avg<<endl;
  cout<<"Averaged Kin: "<<GEGM_kin<<endl;
  cout<<"Diff: "<<abs((GEGM_avg - GEGM_kin) / GEGM_avg)<<endl;

}
