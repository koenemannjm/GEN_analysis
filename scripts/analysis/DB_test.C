
#include "../../include/gen-ana.h"

map<int,double> fHe3Pol;
map<int,double> fBeamPol;

void getDB(TString DB_file){

  // Entries should follow this form:
  //{var name, variable pointer, vairable description, 1/0 (mandatory/not mandatory variable)}
  DBparse::DBRequest request[] = {
    {"He3 Polarization", &fHe3Pol, "He3 Polarization", 1},
    {"Beam Polarization", &fBeamPol, "Beam Polarization", 1}
  };
  
  const int nvar = sizeof(request) / sizeof(request[0]);
  
  DB_load(DB_file,request,nvar);
}


void DB_test(TString conf = "GEN2"){

  TString DB_file = "../../DB/GEN2.csv";
 
  getDB(DB_file);

  cout<<fHe3Pol[2084]<<endl;
  cout<<fBeamPol[2084]<<endl;

}
