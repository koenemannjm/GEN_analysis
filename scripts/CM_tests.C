const int nfiber = 24;
const int nadc = 16;
const int nchan = 128;

const int fN_MPD_TIME_SAMP = 6;
const int fNeventsCommonModeLookBack = 100;
int fNeventsRollingAverage_by_APV[nfiber][nadc];
deque<double> fCommonModeResultContainer_by_APV[nfiber][nadc];
double fCommonModeRollingAverage_by_APV[nfiber][nadc];
double fCommonModeRollingRMS_by_APV[nfiber][nadc];
double CommonModeCorrection[nfiber][nadc][fN_MPD_TIME_SAMP];
int fCorrectCommonModeMinStrips = 20;
double fCorrectCommonMode_Nsigma = 5.0;
int fCommonModeMinStripsInRange = 10;
double sorting_strip_low = 54;
double sorting_strip_high = 54;
int fCommonModeNumIterations = 3;
int fZeroSuppressRMS = 3;
double fRMS_ConversionFactor = sqrt(fN_MPD_TIME_SAMP);

double strip_ADC[nfiber][nadc][nchan][fN_MPD_TIME_SAMP];

double PedMean[nfiber][nadc][nchan] = {0};
double PedRMS[nfiber][nadc][nchan] = {0};
double CMmean[nfiber][nadc] = {0};
double CMrms[nfiber][nadc] = {0};


double CM_sorting_calc[15][6];
double CM_Danning_calc[15][6];

const int nruns = 1;
//int runs[nruns] = {13516, 13420, 13598, 13599, 13600, 13601, 13602,13603,13604};
//TString titles[nruns] = {"Pedestal Run","Cosmic Run","1 uA on LD2, op V","3 uA on LD2, op V","3 uA on LD2, 1500 V","5 uA on LD2, op V","7 uA on LD2, op V","7 uA on LD2, 1500 V","7 uA on LD2, 0 V"};
//double currents[nruns] = {0,0,1,3,3,5,7, 7, 7};

//int runs[nruns] = {13420, 13598, 13599, 13601, 13602, 13791, 13792};
int runs[nruns] = {13602};
//TString titles[nruns] = {"Cosmic Run","1 uA on LD2","3 uA on LD2","5 uA on LD2","7 uA on LD2","20 uA on LD2","34.5 uA on LD2"};
TString titles[nruns] = {"7 uA on LD2"};
double currents[nruns];

int APVMAP[128];


void InitAPVMAP(){
  
  for( UInt_t i=0; i<128; i++ ){
    Int_t strip1 = 32*(i%4) + 8*(i/4) - 31*(i/16);
    Int_t strip2 = strip1 + 1 + strip1 % 4 - 5 * ( ( strip1/4 ) % 2 );
    Int_t strip3 = ( strip2 % 2 == 0 ) ? strip2/2 + 32 : ( (strip2<64) ? (63 - strip2)/2 : 127 + (65-strip2)/2 ); 
    APVMAP[i] = strip3;
  }

}


void LoadPedestals( TString pedfilename ){

  
  std::ifstream pedfile( pedfilename.Data() );

  if( !pedfile.good() ){
    pedfile.close();

    std::cout << "Warning: could not find ped file " << pedfilename << " in working directory, pedestals not loaded" << std::endl;

    return;
  } else {
    std::cout << "Found pedestal file " << pedfilename << endl;
  }

  //temporary storage for pedestals loaded from file:
  //vector<int> Slot, MPD, ADC_ch, APV_ch;
  //vector<double> pedmean, pedrms;

  //let's define a unique index as
  // index = apvchan + 128*adc_ch + 16*128*MPD +
  //map by slot, MPD, and adc_ch
  //std::map<int, std::map<int,std::vector<int> > > Slot;
  
  // std::map<int, std::map<int,std::vector<int> > > APVChan;

  //parse the file: Let's do this a bit more intelligently using TString:
  std::string currentline;

  int crate=0, slot=0, mpd=0, adc_ch=0;
  
  while( std::getline(pedfile, currentline) ){
    //TString currentline;
    if( pedfile.eof() ) break;

    if( currentline[0] != '#' ){    
      
      std::istringstream is(currentline);

      string dummy;
      
      if ( currentline.find("APV") == 0 ){
	is >> dummy >> crate >> slot >> mpd >> adc_ch;
	//std::cout << "crate, slot, mpd, adc_ch = " << crate << ", " << slot << ", " << mpd << ", " << adc_ch << std::endl;
      } else {
      
	int apvchan;
	double mean, rms;
	//for( UInt_t i=0; i<128; i++ ){
	is >> apvchan >> mean >> rms;
	
	PedMean[mpd][adc_ch][apvchan] = mean;
	PedRMS[mpd][adc_ch][apvchan] = rms;


	// std::cout << "mapped value of (apvchan, mean, rms) = ( "
	// 	  << APVChan[crate][slot][index].back() << ", "
	// 	  << PedMean[crate][slot][index].back() << ", "
	// 	  << PedRMS[crate][slot][index].back() << ")" << std::endl;
	
      }
    }
  }

}


void LoadCM( TString CMfilename ){

  
  std::ifstream CMfile( CMfilename.Data() );

  if( !CMfile.good() ){
    CMfile.close();

    std::cout << "Warning: could not find CM file " << CMfilename << " in working directory, pedestals not loaded" << std::endl;

    return;
  } else {
    std::cout << "Found CM file " << CMfilename << endl;
  }


  //parse the file: Let's do this a bit more intelligently using TString:
  std::string currentline;

  int crate=0, slot=0, mpd=0, adc_ch=0;
  
  while( std::getline(CMfile, currentline) ){
    //TString currentline;
    if( CMfile.eof() ) break;

    if( currentline[0] != '#' ){    
      
      std::istringstream is(currentline);

      string dummy;
      
      double mean, rms;
   
      is >> crate >> slot >> mpd >> adc_ch >> mean >> rms;
       
      CMmean[mpd][adc_ch] = mean;
      CMrms[mpd][adc_ch] = rms;
    }
  }

}




double Sorting_CM(TH1F *hAPV, int isamp){

  vector<double> adc;

  for(int istrip = 0; istrip < 128; istrip++){
    adc.push_back(hAPV->GetBinContent(istrip + 129*isamp));
  }

  std::sort( adc.begin(), adc.end() );

  double CM = 0;
  int n_keep = 0;

  for(int ihit = sorting_strip_low; ihit < 128 - sorting_strip_high; ihit++){
    CM += adc[ihit];
    n_keep++;
  }

  return CM/n_keep;

}



double Danning_CM_online(TH1F *hAPV, int fiber, int adc_ch, int isamp){

  double CM_1 = 0;
  double CM_2 = 0;
  int n_keep = 0;

  //cout<<CMmean[fiber][adc_ch]<<" "<<CMrms[fiber][adc_ch]<<endl;
  for(int istrip = 0; istrip < 128; istrip++){
    double adc = hAPV->GetBinContent(istrip + 129*isamp);
    //if(adc > CMmean[fiber][adc_ch] - 5*CMrms[fiber][adc_ch] && adc < CMmean[fiber][adc_ch] + 5*CMrms[fiber][adc_ch]){
      if(adc < CMmean[fiber][adc_ch] + 5*CMrms[fiber][adc_ch]){
      CM_1 += adc;
      n_keep++;
    }
  }

  CM_1 /= n_keep;
  n_keep = 0;


  for(int istrip = 0; istrip < 128; istrip++){
    double adc = hAPV->GetBinContent(istrip + 129*isamp);
    //if(adc > CM_1 - 3*PedRMS[fiber][adc_ch][istrip] && adc < CM_1 + 3*PedRMS[fiber][adc_ch][istrip]){
      if(adc < CM_1 + PedRMS[fiber][adc_ch][APVMAP[istrip]]*3){
      CM_2 += adc;
      n_keep++;
    }
  }



  return CM_2/n_keep;

}
  


//double Danning_CM_offline(TH1F *hAPV, int fiber, int adc_ch, int isamp){
double Danning_CM_offline(double strip_ADC[][fN_MPD_TIME_SAMP], int fiber, int adc_ch, int isamp, int flag){

  double CM_1 = 0;
  double CM_2 = 0;
  int n_keep = 0;
  
  double CM_mean_offline = CMmean[fiber][adc_ch];
  double CM_rms_offline = CMrms[fiber][adc_ch];
  
  if( fNeventsRollingAverage_by_APV[fiber][adc_ch] >= std::min(100, fN_MPD_TIME_SAMP*fNeventsCommonModeLookBack ) ){
    CM_mean_offline = fCommonModeRollingAverage_by_APV[fiber][adc_ch];
    CM_rms_offline = fCommonModeRollingRMS_by_APV[fiber][adc_ch];
  }
  
  //if(fiber == 7 && adc_ch == 8 && flag == 1) cout<<fNeventsRollingAverage_by_APV[fiber][adc_ch]<<" "<<CM_mean_offline<<" "<<CM_rms_offline<<endl;
  double cm_min = CM_mean_offline - fZeroSuppressRMS*CM_rms_offline;
  double cm_max = CM_mean_offline + fZeroSuppressRMS*CM_rms_offline;
  
  double cm_temp = 0.0;
  int nstrips_final;
  
  if(flag == 2) cout<<"\n\n\n";
  for( int iter=0; iter<fCommonModeNumIterations; iter++ ){
    int nstripsinrange=0;
    double sumADCinrange=0.0;
    
    
    if(flag == 2) cout<<"\n\n\n";
    if(flag == 2) cout<<"iteration "<<iter<<endl;
    for( int istrip=0; istrip<nchan; istrip++ ){

      //double ADCtemp =  hAPV->GetBinContent(istrip + 129*isamp);
      double ADCtemp =  strip_ADC[istrip][isamp];
      
      

      //on iterations after the first iteration, reject strips with signals above nsigma * pedrms:
      double rmstemp = PedRMS[fiber][adc_ch][istrip];      
      

      if(flag == 1){
	double strip_sum = 0;
	for(int itsamp=0; itsamp < 6; itsamp++)
	  strip_sum += ADCtemp - CM_Danning_calc[adc_ch][itsamp];
	if(strip_sum/fN_MPD_TIME_SAMP < 3*rmstemp) continue;
      }


      double mintemp = cm_min;
      double maxtemp = cm_max;
      
      
      if( iter > 0 ) {
	maxtemp = cm_temp + fZeroSuppressRMS*rmstemp*fRMS_ConversionFactor; //2.45 = sqrt(6), don't want to calculate sqrt every time
	//mintemp = 0.0;
	mintemp = cm_temp - fZeroSuppressRMS*rmstemp*fRMS_ConversionFactor;
      }
      if(flag == 2) cout<<"ADC "<<ADCtemp<<"    cm_min "<<mintemp<<"    cm_max "<<maxtemp<<endl;
      
  
      
      if( ADCtemp >= mintemp && ADCtemp <= maxtemp ){
	nstripsinrange++;
	sumADCinrange += ADCtemp;
      }
    }
    if(flag == 2) cout<<"nstrips in range  "<<nstripsinrange<<endl;
    if( nstripsinrange >= 1 ){ //require minimum 10 strips in range:
      cm_temp = sumADCinrange/double(nstripsinrange);
      if(flag == 2) cout<<"CM iter "<<iter<<" "<<cm_temp<<endl;
    } else if( iter==0 ){ //not enough strips on FIRST iteration, use mean from sorting-method:
      
      //return Sorting_CM( hAPV, isamp);
    }
    nstrips_final = nstripsinrange;
  } //loop over iterations for "Danning method" CM calculation
  
    //std::cout << std::endl;
  if(flag == 2) cout<<"CM result  "<<cm_temp<<endl;

  if(flag == 1) {
    double occu = nstrips_final*1.0 / 128;
    //if(occu > 1) occu = 1;
    //cm_temp -= 3*8.3*(1 - occu); //Online sigma cut * average pedestal RMS
  }

  return cm_temp;
}










void UpdateRollingCommonModeAverage( int ifiber, int iapv, double CM_sample ){
  //This gets called for each time sample for each APV or each full readout event or whenever BUILD_ALL_SAMPLES is true and CM_ENABLED is false:
  //There are two cases to handle:
  // 1) before the container is full, meaning less than fN_MPD_TIME_SAMP*fNeventsCommonModeLookBack have been added. In this case we increment everything.
  // 2) after the container is full, meaning the earliest sample needs to roll off the average and one new sample has to be added at the end.
  
  UInt_t N = fNeventsRollingAverage_by_APV[ifiber][iapv];
  UInt_t Nmax = fN_MPD_TIME_SAMP*fNeventsCommonModeLookBack;
  
  double sum, sum2;
  
  // For now "pos" and "axis" are unused. Comment out to suppress compiler warning:
  // UInt_t pos = fMPDmap[iapv].pos;
  // UInt_t axis = fMPDmap[iapv].axis;

  // For now these are unused. Comment out to suppress compiler warning:
  // double cm_mean_from_DB = (axis == SBSGEM::kUaxis) ? fCommonModeMeanU[pos] : fCommonModeMeanV[pos];
  // double cm_rms_from_DB = (axis == SBSGEM::kUaxis) ? fCommonModeRMSU[pos] : fCommonModeRMSV[pos];   
  
  //std::cout << "Updating rolling average common mode from full readout events for module " << GetName() << " iapv = " << iapv << std::endl;
  
  if( N < Nmax ){
    //before reaching the size of the look back window, we just add the common-mode samples onto the end of the array

    fCommonModeResultContainer_by_APV[ifiber][iapv][N] = CM_sample;

    if( N == 0 ){ //First sample, initialize all sums/averages:
      fCommonModeRollingAverage_by_APV[ifiber][iapv] = CM_sample;
      fCommonModeRollingRMS_by_APV[ifiber][iapv] = 0.0;
      sum = CM_sample;
      sum2 = pow(CM_sample,2);
    } else { //Second and subsequent samples: increment sums, recalculate
      double oldavg = fCommonModeRollingAverage_by_APV[ifiber][iapv];
      double oldrms = fCommonModeRollingRMS_by_APV[ifiber][iapv];
      
      sum = N*oldavg + CM_sample;
      sum2 = N * (pow(oldrms,2) + pow(oldavg,2)) + pow(CM_sample,2);

      double newavg = sum/double(N+1);
      double newrms = sqrt( sum2/double(N+1) - pow(newavg,2) );

      fCommonModeRollingAverage_by_APV[ifiber][iapv] = newavg;
      fCommonModeRollingRMS_by_APV[ifiber][iapv] = newrms;


    }

    // std::cout << "(N, average, rms)=(" << N << ", " << fCommonModeRollingAverage_by_APV[ifiber][iapv]
    // 	      << ", " << fCommonModeRollingRMS_by_APV[ifiber][iapv] << ")" << std::endl;
    fNeventsRollingAverage_by_APV[ifiber][iapv] = N+1;
  
  } else {
      
    //grab the earliest sample in the rolling average:
    double oldfirstsample = fCommonModeResultContainer_by_APV[ifiber][iapv].front();
    
    //The net result of the following two operations should be to keep the container size the same:
    fCommonModeResultContainer_by_APV[ifiber][iapv].pop_front(); //remove oldest sample
    fCommonModeResultContainer_by_APV[ifiber][iapv].push_back( CM_sample ); //Insert newest sample at the end
    //we only need to update the calculation for the fact that the
    //earliest sample rolled off and a new sample was added: 
    double oldavg = fCommonModeRollingAverage_by_APV[ifiber][iapv];
    double oldsum = oldavg * Nmax;

    double oldrms = fCommonModeRollingRMS_by_APV[ifiber][iapv];
    // RMS^2 = sum^2/N - avg^2 --> sum^2 = N * (RMS^2 + avg^2)
    double oldsum2 = Nmax * ( pow(oldrms,2) + pow(oldavg,2) );

    //double lastsample = fCommonModeResultContainer_by_APV[ifiber][iapv].back();
    double lastsample = CM_sample;
    
    double newsum = oldsum - oldfirstsample + lastsample;
    double newsum2 = oldsum2 - pow(oldfirstsample,2) + pow(lastsample,2);

    double newavg = newsum/double( Nmax );
    double newrms = sqrt( newsum2/double( Nmax ) - pow(newavg,2) );
    
    fCommonModeRollingAverage_by_APV[ifiber][iapv] = newavg;
    fCommonModeRollingRMS_by_APV[ifiber][iapv] = newrms;


  }
}



//void process_run(int run,bool event_display,bool flipped_events,TH2F *hcompare, TH2F *hped, TH2F *hoccu_corr, TH1F *htime_pos, TH1F *htime_neg, TH2F *hadc_corr){
void process_run(int irun, int run,bool event_display,TH2F *hCM_corr, TH2F *hCM_not_corr){



  //if(event_display || flipped_events) gROOT->SetBatch(kTRUE);
  if(event_display) gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  int signal_threshold = 0;
  if(run == 13420) signal_threshold = 7;
  //if(flipped_events) signal_threshold = 7;

  TString Rootfile = Form("../../Rootfiles/hit_0_e1209019_%i.root",run);

  TString outputfile = Form("../plots/CM_event_display_run%i.pdf",run);
  TString outputfile_flip = Form("../plots/Signal_flip_event_display_run%i.pdf",run);

  TChain *t = new TChain("GEMHit");
  
  t->Add(Rootfile);

  int maxch = 60000;

  int evtID;
  int nch;
  int mpd[maxch];
  int adc_ch[maxch];
  int strip[maxch];
  int adc_all[6][maxch];

  t->SetBranchAddress("evtID",&evtID);
  t->SetBranchAddress("nch",&nch);
  t->SetBranchAddress("MPD",mpd); 
  t->SetBranchAddress("ADC_ch",adc_ch); 
  t->SetBranchAddress("strip",strip); 
  t->SetBranchAddress("adc0",adc_all[0]); 
  t->SetBranchAddress("adc1",adc_all[1]); 
  t->SetBranchAddress("adc2",adc_all[2]); 
  t->SetBranchAddress("adc3",adc_all[3]); 
  t->SetBranchAddress("adc4",adc_all[4]); 
  t->SetBranchAddress("adc5",adc_all[5]);  


  TCanvas *c = new TCanvas("c","",1600,1000);
  c->Divide(2,2);

  TCanvas *c4 = new TCanvas("c4","",1600,1000);
  c4->Divide(2,2);

  int fMPD = 7;
  int nsamp = 0;
  int npos = 0;
  int nneg = 0;
  
  TH1F *hAPV_pedsub[15];   //APV raw view after pedestal subtracted
  TH1F *hAPV_zerosup[15];  //APV raw view after old method of CM correction
  TH1F *hAPV_zerosup_CMadded[15];  //APV raw view after online CM added back in
  TH1F *hAPV_zerosup_corr[15];   //APV raw view after new method of CM correction
  TH2F *hCM_online_sorting = new TH2F(titles[irun],"Danning Online - Sorting;APV;ADC",15,0,15,500,-100,100);
  TH2F *hCM_offline_sorting = new TH2F(titles[irun],"Danning Offline - Sorting;APV;ADC",15,0,15,500,-100,100);
  TH2F *htest = new TH2F("htest","",128,0,128,500,-100,100);
  //TH2F *hCM_offline_sorting_corr = new TH2F("hadc_offline_corr","Danning Online with Corrections - Sorting;APV;ADC",15,0,15,500,-100,100);
  //TH2F *hCM_offline_sorting_nocorr = new TH2F("hadc_offline_nocorr","Danning Online without Corrections - Sorting;APV;ADC",15,0,15,500,-100,100);
  TH1F *hflip = new TH1F("hflip","Strips with Flipped Signal;Time Sample;ADC",6,0,6);



  for(int i=0; i < 15; i++){
    hAPV_pedsub[i] = new TH1F(Form("hAPV_pedsub_%i",i),"Raw APV Data;Strip;ADC",768,0,768);
    hAPV_zerosup[i] = new TH1F(Form("hAPV_zerosup_%i",i),"Online Zero Suppression;Strip;ADC",768,0,768);
    hAPV_zerosup_CMadded[i] = new TH1F(Form("hAPV_zerosup_CMadded_%i",i),"Zero Sup with CM Added Back;Strip;ADC",768,0,768);
    hAPV_zerosup_corr[i] = new TH1F(Form("hAPV_zerosup_corr_%i",i),"Zero Sup with New CM;Strip;ADC",768,0,768);
  }

  int ievent = 0;
  int iflip = 0;
  int nflip = 300;
  int idisplay = 0;
  int ndisplay = 10;

  int ntotal = t->GetEntries();
  
  while(t->GetEntry(ievent++) && ievent < 1000){

    
    for(int ifiber=0; ifiber < nfiber; ifiber++){  
      for(int iapv=0; iapv < nadc; iapv++){
	fCommonModeResultContainer_by_APV[ifiber][iapv].resize(fN_MPD_TIME_SAMP*fNeventsCommonModeLookBack);
	for( int isamp=0; isamp<fN_MPD_TIME_SAMP; isamp++ ){
	  CommonModeCorrection[ifiber][iapv][isamp] = 0.0;
	}
      }
    }  
    
    int neg_occu[15] = {0};
    int pos_occu[15] = {0};
    double adc_neg[15] = {0};
    double adc_pos[15] = {0};

    bool full_readout = false;
    bool CM_corrected[16] = {false};

    if(ievent%1000 == 0)
      cout<<ievent<<endl;

    if(ievent < 130)
      full_readout = true;
    
    for(int ich=0; ich < nch; ich++){

      if(mpd[ich] == fMPD){

	for(int isamp=0; isamp < 6; isamp++){
	  double adc_pedsub = adc_all[isamp][ich] - PedMean[mpd[ich]][adc_ch[ich]][strip[ich]];
	  strip_ADC[mpd[ich]][adc_ch[ich]][strip[ich]][isamp] = adc_pedsub;

	  int real_strip = APVMAP[strip[ich]];

	  hAPV_pedsub[adc_ch[ich]]->SetBinContent(real_strip + 129*isamp,adc_pedsub);
	}
      }
    }
    
    for(int i=0; i < 15; i++){
     

      for(int isamp=0; isamp < 6; isamp++){
	hAPV_pedsub[i]->SetBinContent(128*isamp,0);
	
	double CM_Danning = Danning_CM_online(hAPV_pedsub[i],fMPD,i,isamp);
	double CM_Sorting = Sorting_CM(hAPV_pedsub[i],isamp);
	CM_sorting_calc[i][isamp] = CM_Sorting;
	CM_Danning_calc[i][isamp] = CM_Danning;
	
	
	double CM_diff = CM_Danning - CM_Sorting;
	
	nsamp++;
	//cout<<CM_Sorting<<" "<<CMmean[fMPD][i]<<endl;
	//hped->Fill(i,CM_Sorting - CMmean[fMPD][i]);
	//hcompare->Fill(i,CM_diff);
	//hCM_online_sorting->Fill(i,CM_diff);
	//hCM_offline_sorting->Fill(i,Danning_CM_offline(hAPV_pedsub[i],fMPD,i,isamp) - CM_Sorting);	
	if(full_readout) UpdateRollingCommonModeAverage( fMPD, i, Danning_CM_offline(strip_ADC[fMPD][i],fMPD,i,isamp,0));
	//if(full_readout && i == 10) cout<<evtID<<" "<<fCommonModeRollingAverage_by_APV[fMPD][i]<<endl;
      }
    }

    
    if(!full_readout){
      for(int ifiber=0; ifiber < nfiber; ifiber++){
	for(int iapv=0; iapv < nadc; iapv++){
	  for( int isamp=0; isamp<fN_MPD_TIME_SAMP; isamp++ ){
	    
	    double CM_meas = CM_Danning_calc[iapv][isamp];
	    double CM_expect_mean, CM_expect_rms;
	    if(fNeventsRollingAverage_by_APV[ifiber][iapv] >= std::min(100,fNeventsCommonModeLookBack*fN_MPD_TIME_SAMP) ){
	      //If we have a critical mass of events in the rolling CM average for this to be a reliable estimate, use it: 
	      CM_expect_mean = fCommonModeRollingAverage_by_APV[ifiber][iapv];
	      CM_expect_rms = fCommonModeRollingRMS_by_APV[ifiber][iapv]; 
	    } else { //use database value:
	      CM_expect_mean = CMmean[ifiber][iapv];
	      CM_expect_rms = CMrms[ifiber][iapv];
	    }

	    //if(ifiber == fMPD && iapv == 10) cout<<fNeventsRollingAverage_by_APV[ifiber][iapv]<<" "<<CM_expect_mean<<" "<<fCommonModeRollingAverage_by_APV[ifiber][iapv]<<endl;
	    //if(ifiber == fMPD && iapv == 10) cout<<ifiber<<" "<<iapv<<" "<<isamp<<" "<<CM_meas<<" "<<CM_expect_mean<<" "<<fCorrectCommonMode_Nsigma * CM_expect_rms<<" "<<CM_expect_mean - fCorrectCommonMode_Nsigma * CM_expect_rms<<endl;
	    //if(ifiber == fMPD && iapv == 10) cout<<ifiber<<" "<<iapv<<" "<<isamp<<" "<<fNeventsRollingAverage_by_APV[ifiber][iapv]<<" "<<CM_meas<<" "<<CM_expect_mean<<" "<<CMmean[ifiber][iapv]<<" "<<fCommonModeRollingAverage_by_APV[ifiber][iapv]<<endl;
	    //if(ifiber == fMPD && evtID < 142 && iapv == 8) cout<<fCommonModeRollingAverage_by_APV[ifiber][iapv]<<endl;
	    //if( CM_meas < CM_expect_mean - fCorrectCommonMode_Nsigma * CM_expect_rms ){
	    if( 1 == 1 ){
	      //if(ifiber == fMPD && iapv == 10) cout<<ifiber<<" "<<iapv<<" "<<isamp<<" "<<CM_meas<<" "<<CM_expect_mean<<" "<<fCorrectCommonMode_Nsigma * CM_expect_rms<<" "<<CM_expect_mean - fCorrectCommonMode_Nsigma * CM_expect_rms<<endl;
	      // The online common mode appears to have a large negative bias relative to the expectation. 
	      // Try to correct the common-mode. To calculate the correction requires us to loop on all the strips on this APV that passed
	      // online zero suppression.
	      // The simplest approach is just to take a simple average of all the strips within +/- some number of standard deviations of the
	      // *EXPECTED* common-mode mean, but this is a biased approach.
	      //if(ifiber == fMPD && iapv == 8) cout<<"\n"<<iapv<<" "<<isamp<<" "<<CM_meas<<" "<<CM_expect_mean<<" "<<fCorrectCommonMode_Nsigma * CM_expect_rms<<endl;
	      UInt_t NstripsInRange = 0; 
	      
	      //Loop on all strips on this APV and calculate raw ADC values from 
	      for(int istrip=0; istrip<nchan; ++istrip ){
		
		Int_t ADCtemp = strip_ADC[fMPD][iapv][istrip][isamp];

		double strip_sum = 0;
		for(int itsamp=0; itsamp < 6; itsamp++)
		  strip_sum += strip_ADC[fMPD][iapv][istrip][itsamp] - CM_Danning_calc[iapv][itsamp];
		if(strip_sum/fN_MPD_TIME_SAMP < 3*PedRMS[ifiber][iapv][istrip]) continue;
		//if(ADCtemp - CM_Danning_calc[iapv][isamp] < 3*PedRMS[ifiber][iapv][istrip]) continue;
		
		//The pedestal has already been subtracted in this case, but we want to add back the common-mode:
		double pedsubADC = ADCtemp; //this is the one that goes into the common-mode calculation

		//if(ifiber == fMPD && iapv == 1) cout<<ifiber<<" "<<iapv<<" "<<isamp<<" "<<pedsubADC - CM_expect_mean<<" "<<fCorrectCommonMode_Nsigma * CM_expect_rms<<endl;
		if( fabs( pedsubADC - CM_expect_mean ) <= fCorrectCommonMode_Nsigma * CM_expect_rms ) NstripsInRange++;
		
	      }
	      //if(ifiber == fMPD && iapv == 1) cout<<isamp<<" nstrips "<<NstripsInRange<<endl;
	      if( NstripsInRange >= fCommonModeMinStripsInRange ){
		//Correction to be ADDED to ADC value to get corrected value:
		CommonModeCorrection[ifiber][iapv][isamp] = CM_meas -  Danning_CM_offline(strip_ADC[fMPD][iapv],fMPD,iapv,isamp, 1) + 3*8.3*(1 - NstripsInRange*1.0 / 128);   
		//CommonModeCorrection[ifiber][iapv][isamp] = CM_meas -  Danning_CM_offline(strip_ADC[fMPD][iapv],fMPD,iapv,isamp, 1);   
		//CommonModeCorrection[ifiber][iapv][isamp] = correction_temp + 3*8.3;   
		//CommonModeCorrection[ifiber][iapv][isamp] = correction_temp;   
		//CommonModeCorrection[ifiber][iapv][isamp] = correction_temp + exp(-0.5*pow(NstripsInRange/,2));   
		//if(ifiber == fMPD && iapv == 12 && correction_temp < -10*8.3) cout<<isamp<<" nstrips "<<NstripsInRange<<" "<<correction_temp<<endl;
		if(ifiber == fMPD && iapv < 15) htest->Fill(NstripsInRange, CM_Danning_calc[iapv][isamp] - CommonModeCorrection[ifiber][iapv][isamp] - CM_sorting_calc[iapv][isamp]);

		//CommonModeCorrection[ifiber][iapv][isamp] -= 8.3*(CommonModeCorrection[ifiber][iapv][isamp]
		if(ifiber == fMPD) CM_corrected[iapv] = true;
		//if(ifiber == fMPD && evtID == 141) cout<<ifiber<<" "<<iapv<<" "<<isamp<<" "<<CM_meas<<" "<<CM_expect_mean - fCorrectCommonMode_Nsigma * CM_expect_rms<<" "<<NstripsInRange<<endl;
		// Uncorrected ADC value = Raw ADC - pedestal - CM_meas
		// Corrected ADC value = Raw ADC - pedestal - CM_corrected
		// Corrected ADC value = Uncorrected ADC value - CM_corrected + CM_meas = Uncorrected ADC value + CommonModeCorrection.
		// Therefore, when common-mode is corrected this way, 
	      } else {
		CommonModeCorrection[ifiber][iapv][isamp] = 0.0; //To be added to ADC value! 
	      }
	    }
	    if(ifiber == fMPD && iapv < 15){
	      hCM_online_sorting->Fill(iapv,CM_Danning_calc[iapv][isamp] - CM_sorting_calc[iapv][isamp]);
	      hCM_offline_sorting->Fill(iapv,Danning_CM_offline(strip_ADC[fMPD][iapv],fMPD,iapv,isamp,0) - CM_sorting_calc[iapv][isamp]);
	      if(CommonModeCorrection[ifiber][iapv][isamp] != 0.0){
		hCM_corr->Fill(iapv,CM_Danning_calc[iapv][isamp] - CommonModeCorrection[ifiber][iapv][isamp] - CM_sorting_calc[iapv][isamp]);
	      }
	      if(CommonModeCorrection[ifiber][iapv][isamp] == 0.0) hCM_not_corr->Fill(iapv,CM_Danning_calc[iapv][isamp] - CommonModeCorrection[ifiber][iapv][isamp] - CM_sorting_calc[iapv][isamp]);

	      //cout<<iapv<<" "<<isamp<<" "<<CM_meas<<" "<<CM_expect_mean<<" "<<fCorrectCommonMode_Nsigma * CM_expect_rms<<endl;
	    }
	  }
	}
      }
    }


    for(int ich=0; ich < nch; ich++){


      if(mpd[ich] == fMPD){

	double strip_sum = 0;
	double strip_sum_corr = 0;
	//cout<<"\n\n strip "<<strip[ich]<<endl;
	for(int isamp=0; isamp < 6; isamp++){
	  strip_sum += strip_ADC[mpd[ich]][adc_ch[ich]][strip[ich]][isamp] - CM_Danning_calc[adc_ch[ich]][isamp];
	  strip_sum_corr += strip_ADC[mpd[ich]][adc_ch[ich]][strip[ich]][isamp] - CM_Danning_calc[adc_ch[ich]][isamp] + CommonModeCorrection[mpd[ich]][adc_ch[ich]][isamp];
	  //cout<<isamp<<" "<<strip_ADC[mpd[ich]][adc_ch[ich]][strip[ich]][isamp] - CM_Danning_calc[adc_ch[ich]][isamp]<<endl;
	}
	

	for(int isamp=0; isamp < 6; isamp++){
	  double adc_final = strip_ADC[mpd[ich]][adc_ch[ich]][strip[ich]][isamp] - CM_Danning_calc[adc_ch[ich]][isamp];
	  int real_strip = APVMAP[strip[ich]];	  
	  
	  if(strip_sum/fN_MPD_TIME_SAMP < 3*PedRMS[mpd[ich]][adc_ch[ich]][strip[ich]]) adc_final = 0;
	  //cout<<isamp<<" "<<adc_final<<endl;
	  hAPV_zerosup[adc_ch[ich]]->SetBinContent(real_strip + 129*isamp,adc_final);
	  
	  if(adc_final != 0) hAPV_zerosup_CMadded[adc_ch[ich]]->SetBinContent(real_strip + 129*isamp,adc_final + CM_Danning_calc[adc_ch[ich]][isamp]);
	  else
	    hAPV_zerosup_CMadded[adc_ch[ich]]->SetBinContent(real_strip + 129*isamp,0);
	  
	  adc_final += CommonModeCorrection[mpd[ich]][adc_ch[ich]][isamp];

	  if(strip_sum_corr/fN_MPD_TIME_SAMP < 3*PedRMS[mpd[ich]][adc_ch[ich]][strip[ich]]) adc_final = 0;
	  
	  hAPV_zerosup_corr[adc_ch[ich]]->SetBinContent(real_strip + 129*isamp,adc_final);
	}
      }
    }
    

    
      
    if(ievent > 130 && event_display){
    for(int iapv = 0; iapv < 15; iapv++){

      TLegend *leg = new TLegend(0.65,0.7,0.9,0.9);
      TLegend *leg2 = new TLegend(0.65,0.7,0.9,0.9);
      
      if(!CM_corrected[iapv]) continue;
            
      hAPV_pedsub[iapv]->SetMinimum(0);
      hAPV_zerosup[iapv]->SetMinimum(0);
      hAPV_zerosup_corr[iapv]->SetMinimum(0);
      
      //cout<<CMmean[fMPD][iapv]<<" "<<fCommonModeRollingAverage_by_APV[fMPD][iapv]<<" "<<CM_Danning_calc[iapv][0]<<" "<<CM_sorting_calc[iapv][0]<<endl;

      TLine *l_Danning[6];
      TLine *l_Sorting[6];
      TLine *l_Danning_offline[6];
      TLine *l_low;
      TLine *l_high;
        
      c->cd(1);
      hAPV_pedsub[iapv]->Draw();
      cout<<"APV "<<iapv<<endl;
      for(int isamp=0; isamp < 6; isamp++){
	//l_Danning[isamp] = new TLine(128*isamp, Danning_CM_offline(strip_ADC[fMPD][iapv],fMPD,iapv,isamp,0), 128*(isamp+1), Danning_CM_offline(strip_ADC[fMPD][iapv],fMPD,iapv,isamp,0));
	l_Danning[isamp] = new TLine(128*isamp, CM_Danning_calc[iapv][isamp], 128*(isamp+1), CM_Danning_calc[iapv][isamp]);
	l_Danning[isamp]->SetLineColor(kRed);
	l_Danning[isamp]->Draw("same");

	l_Sorting[isamp] = new TLine(128*isamp, CM_sorting_calc[iapv][isamp], 128*(isamp+1), CM_sorting_calc[iapv][isamp]);
	l_Sorting[isamp]->SetLineColor(kGreen);
	l_Sorting[isamp]->Draw("same");

	l_Danning_offline[isamp] = new TLine(128*isamp, Danning_CM_offline(strip_ADC[fMPD][iapv],fMPD,iapv,isamp,1), 128*(isamp+1), Danning_CM_offline(strip_ADC[fMPD][iapv],fMPD,iapv,isamp,1));
	l_Danning_offline[isamp]->SetLineColor(kRed);

	cout<<isamp<<" "<<CM_Danning_calc[iapv][isamp] - CM_sorting_calc[iapv][isamp]<<" "<<CommonModeCorrection[fMPD][iapv][isamp]<<endl;
	Danning_CM_offline(strip_ADC[fMPD][iapv],fMPD,iapv,isamp, 1); 
	cout<<isamp<<" "<<CM_sorting_calc[iapv][isamp]<<" "<<CM_Danning_calc[iapv][isamp]<<" "<<Danning_CM_offline(strip_ADC[fMPD][iapv],fMPD,iapv,isamp,1)<<endl;

      }

      l_low = new TLine(0, fCommonModeRollingAverage_by_APV[fMPD][iapv] - fCorrectCommonMode_Nsigma * fCommonModeRollingRMS_by_APV[fMPD][iapv], 128*6, fCommonModeRollingAverage_by_APV[fMPD][iapv] - fCorrectCommonMode_Nsigma * fCommonModeRollingRMS_by_APV[fMPD][iapv]);
      l_high = new TLine(0, fCommonModeRollingAverage_by_APV[fMPD][iapv] + fCorrectCommonMode_Nsigma * fCommonModeRollingRMS_by_APV[fMPD][iapv], 128*6, fCommonModeRollingAverage_by_APV[fMPD][iapv] + fCorrectCommonMode_Nsigma * fCommonModeRollingRMS_by_APV[fMPD][iapv]);
      l_low->SetLineColor(kBlue);
      l_high->SetLineColor(kBlue);
      
      leg->AddEntry(l_Sorting[0],"Sorting CM","l");
      leg->AddEntry(l_Danning[0],"Online Danning CM","l");
      leg->Draw("same");

      c->cd(2);
      hAPV_zerosup[iapv]->Draw();

      c->cd(3);
      hAPV_zerosup_CMadded[iapv]->Draw();
      l_low->Draw("same");
      l_high->Draw("same");
      
      for(int isamp=0; isamp < 6; isamp++)
	l_Danning_offline[isamp]->Draw("same");

      leg2->AddEntry(l_Danning[0],"Offline Danning CM","l");
      leg2->AddEntry(l_low,"CM Cuts","l");
      leg2->Draw("same");

      c->cd(4);
      hAPV_zerosup_corr[iapv]->Draw();

      c->Update();

      cout<<"Showing event "<<evtID<<endl;
      cout << "press any key to continue (q to quit):" << endl;
      TString reply;
      reply.ReadLine(cin,kFALSE);
      
      if( reply.BeginsWith("q") ) break;

    }
    }
    

    /*
    hAPV[i]->SetMinimum(0);
    
    if(event_display){
      hAPV[i]->Draw();
      
      for(int isamp=0; isamp < 6; isamp++){
	l_Danning[isamp]->Draw("same");
	l_sorting[isamp]->Draw("same");	  
	//l_high[isamp]->Draw("same");	  
      }
      
      if(ievent == 1){      
	if(i == 0){      
	  leg->AddEntry(l_Danning[0],"Danning Method","l");
	  leg->AddEntry(l_sorting[0],"Sorting Method","l");
	}
	
	leg->Draw("same");
      }
    }
    
    
    
    
    if(event_display){
      if(event_signal || run != 13420){
	idisplay++;
	
	c->Update();
	
	if(idisplay == 1)
	  c->Print(outputfile + "(");
	else if(idisplay == ndisplay){
	  c->Print(outputfile + ")");
	  break;
	}
	else
	  c->Print(outputfile);   
      }
      else if(!event_signal && run == 13420 && ievent == ntotal)
	c->Print(outputfile + ")");
    }
  
    */  

    TString reply = "no";;
    for(int iapv = 0; iapv < 15; iapv++){
      
      if(!CM_corrected[iapv]) continue;
         
      if(ievent > 130 && event_display){
	cout<<"Showing event "<<evtID<<endl;
	cout << "press any key to continue (q to quit):" << endl;
	reply.ReadLine(cin,kFALSE);
	if( reply.BeginsWith("q") ) break;
	
      }
    }
    if( reply.BeginsWith("q") ) break;
  }
  
  c4->cd(1);
  hCM_online_sorting->Draw("colz");

  c4->cd(2);
  hCM_offline_sorting->Draw("colz");
  //htest->Draw("colz");
    
  c4->cd(3);
  hCM_corr->Draw("colz");
 
  c4->cd(4);
  hCM_not_corr->Draw("colz");

  if(irun == 0) c4->Print("../plots/CM_bias_correction.pdf(");
  else c4->Print("../plots/CM_bias_correction.pdf");

}

  
  




void CM_tests(bool event_display = false, bool flipped_events = false){

  TString output = "../plots/neg_pulse_study.pdf";
  
  int ped_run = 13594;

  TString Pedfile = Form("../../Rootfiles/daq_ped_bb_gem_run%i.dat",ped_run);
  TString CMfile = Form("../../Rootfiles/CommonModeRange_bb_gem_run%i.txt",ped_run); 
  
  LoadPedestals(Pedfile);
  LoadCM(CMfile);
  InitAPVMAP();

  TH2F *hcompare[nruns];
  TH2F *hoccu_pos[nruns];
  TH2F *hoccu_neg[nruns];
  
  TH1F *htime_pos[nruns];
  TH1F *htime_neg[nruns];
  TH2F *hped[nruns];
  TH2F *hoccu_corr[nruns];
  TH2F *hadc_corr[nruns];

  TH2F *hCM_corr[nruns];
  TH2F *hCM_not_corr[nruns];

  for( int irun = 0; irun < nruns; irun++){
    hcompare[irun] = new TH2F(Form("hcompare_%i",irun),"Common Mode Method Comparison;APV;Danning - Sorting Method CM (ADC)",15,0,15,100,-80,80);
    hped[irun] = new TH2F(Form("hped_%i",irun),"Common Mode Comparison to Pedestal;APV;Sorting - CM from Pedestal (ADC)",15,0,15,100,-80,80);

    hoccu_pos[irun] = new TH2F(Form("hoccu_pos_%i",irun),"APV Positive Signal Occupancy;APV;Occupancy",15,0,15,100,0,1);
    hoccu_neg[irun] = new TH2F(Form("hoccu_neg_%i",irun),"APV Negative Signal Occupancy;APV;Occupancy",15,0,15,100,0,1);
    
    htime_pos[irun] = new TH1F(Form("htime_pos_%i",irun),"APV Positive Signal Peak Time;Peak Time Sample;",6,0,6);
    htime_neg[irun] = new TH1F(Form("htime_neg_%i",irun),"APV Negative Signal Peak Time;Peak Time Sample;",6,0,6);  
    hoccu_corr[irun] = new TH2F(Form("hoccu_corr_%i",irun),"Occupancy Correlations;Positive Occupancy;Negative Occupancy",30,0,0.5,30,0,0.5);

    hadc_corr[irun] = new TH2F(Form("hadc_corr_%i",irun),"ADC Correlations;Positive Sinal Mean ADC Sum;Negative Signal Mean ADC Sum",200,-100,4000,200,-1400,50);

    hCM_corr[irun] = new TH2F(titles[irun],"Danning Online with Corrections - Sorting;APV;ADC",15,0,15,500,-100,100);    
    hCM_not_corr[irun] = new TH2F(titles[irun],"Danning Online without Corrections - Sorting;APV;ADC",15,0,15,500,-100,100);    

    cout<<"Processing run "<<runs[irun]<<endl;
    //process_run(runs[irun],event_display,flipped_events,hcompare[irun],hped[irun],hoccu_corr[irun],htime_pos[irun],htime_neg[irun],hadc_corr[irun]);
    process_run(irun, runs[irun],event_display,hCM_corr[irun], hCM_not_corr[irun]);

    //gPad->Print("../plots/CM_bias_correction.pdf");
  }

    
  double y_CM_corr[nruns][nadc], y_CM_not_corr[nruns][nadc], x_APV[nadc];
  double y_CM_corr_err[nruns][nadc], y_CM_not_corr_err[nruns][nadc];

  for(int iapv=0; iapv < nadc-1; iapv++) x_APV[iapv] = iapv + 1;

  TCanvas *c3 = new TCanvas("c3","",1400,600);
  c3->Divide(2,1);

  TGraphErrors *g_CM_corr[nruns];
  TGraphErrors *g_CM_not_corr[nruns];

  //  g_neg->GetYaxis()->SetRangeUser(0,0.12);
  //g_pos->GetYaxis()->SetRangeUser(0,0.12);

  TLegend *legend = new TLegend(0.60,0.64,0.9,0.9);
  int icolor = 0;
  
  for(int irun = 0; irun < nruns; irun++){
    icolor++;
    if(icolor == 5 || icolor == 10) icolor ++;

    for(int iapv=0; iapv < nadc; iapv++){
      y_CM_corr[irun][iapv] = hCM_corr[irun]->ProjectionY("",x_APV[iapv],x_APV[iapv])->GetMean();
      y_CM_corr_err[irun][iapv] = hCM_corr[irun]->ProjectionY("",x_APV[iapv],x_APV[iapv])->GetRMS();
      
      y_CM_not_corr[irun][iapv] = hCM_not_corr[irun]->ProjectionY("",x_APV[iapv],x_APV[iapv])->GetMean();
      y_CM_not_corr_err[irun][iapv] = hCM_not_corr[irun]->ProjectionY("",x_APV[iapv],x_APV[iapv])->GetRMS();
    }
    
    g_CM_corr[irun] = new TGraphErrors(nadc-1,x_APV,y_CM_corr[irun],0,y_CM_corr_err[irun]);
    g_CM_corr[irun]->SetTitle("CM Corrected Result - Sorting Method;APV;ADC");
    g_CM_corr[irun]->SetMarkerStyle(8);
    g_CM_corr[irun]->SetMarkerColor(icolor);
    g_CM_corr[irun]->SetLineColor(icolor);
    legend->AddEntry(g_CM_corr[irun],titles[irun],"p");
    
    g_CM_not_corr[irun] = new TGraphErrors(nadc-1,x_APV,y_CM_not_corr[irun],0,y_CM_not_corr_err[irun]);
    g_CM_not_corr[irun]->SetTitle("CM Not Corrected Result - Sorting Method;APV;ADC");
    g_CM_not_corr[irun]->SetMarkerStyle(8);
    g_CM_not_corr[irun]->SetMarkerColor(icolor);
    g_CM_not_corr[irun]->SetLineColor(icolor); 
     
    c3->cd(1);
    if(irun == 0)   
      g_CM_corr[irun]->Draw("AP");
    else
      g_CM_corr[irun]->Draw("P");


    c3->cd(2);
    if(irun == 0)
      g_CM_not_corr[irun]->Draw("AP");
    else
      g_CM_not_corr[irun]->Draw("P");

    
    g_CM_corr[irun]->GetYaxis()->SetRangeUser(-30,60);
    g_CM_not_corr[irun]->GetYaxis()->SetRangeUser(-80,20);
    

  }

  c3->cd(1);
  legend->SetTextSize(0.035);
  legend->Draw("same");

  c3->Print("../plots/CM_bias_correction.pdf)");

  
}
