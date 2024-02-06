

double l_BB = 1.63;
double l_HCAL = 17;
double p_e[3] = {2.69,2.73,3.21}; //In GeV
double p_n[3] = {2.37,4.51,6.11};
double c = 2.988e8;
double m_n = 0.928;
double m_e = 0.000511;


double beta(double p, double m){
  return p / sqrt(p*p + m*m);
}


double tof_calc(double p, double m, bool electron){

  double beta = p / sqrt(p*p + m*m);
  double length = l_BB;
  if(!electron) length = l_HCAL;
  
  return length / (beta*c) * 1e9; //convert result to ns

}

void TOF(){

  double p_off = 0.2; //200 MeV

  cout<<"GEN2: "<<(tof_calc(p_n[0],m_n,false) - tof_calc(p_n[0] + p_off,m_n,false)) * 1000<<" ps"<<endl;
  cout<<"GEN3: "<<(tof_calc(p_n[1],m_n,false) - tof_calc(p_n[1] + p_off,m_n,false)) * 1000<<" ps"<<endl;
  cout<<"GEN4: "<<(tof_calc(p_n[2],m_n,false) - tof_calc(p_n[2] + p_off,m_n,false)) * 1000<<" ps"<<endl;
  


}
