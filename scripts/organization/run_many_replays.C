#include "../include/gen-ana.h"
#include "../dflay/src/JSONManager.cxx"

void run_many_replays(const char *configfilename){

  JSONManager *jmgr = new JSONManager(configfilename);

  std::string rootfile_dir = jmgr->GetValueFromKey_str("rootfile_dir");
  int kine = jmgr->GetValueFromKey<int>("GEN_config");
  std::vector<int> runnums; jmgr->GetVectorFromKey<int>("runnums",runnums);
  int nruns = jmgr->GetValueFromKey<int>("Nruns_to_ana"); // # runs to analyze
  TChain *C = new TChain("T");
  if (nruns < 1 || nruns > runnums.size()) nruns = runnums.size();
  for (int i=0; i<nruns; i++) {
    gSystem->Exec(Form("run_replay.sh GEN%i %i",kine,runnums[i]));

  }
  
}
