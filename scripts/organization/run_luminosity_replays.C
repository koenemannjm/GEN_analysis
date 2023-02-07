


void run_luminosity_replays(){


  ifstream infile("GEM_luminosity_good_events.txt");

  
  string first, second;
  int run = -1;
  int first_event, last_event;
  while(!infile.eof()){ 
    infile>>first;
    infile>>second;
       
    if(first == "run") run = stoi(second);
    else{
      first_event = stoi(first);
      last_event = stoi(second);
      
      gSystem->Exec(Form("run_luminosity_replay.sh %i %i %i",run,first_event,last_event));
    }

  }





}
