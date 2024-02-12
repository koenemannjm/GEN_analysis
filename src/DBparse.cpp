#include "../include/DBparse.h"

namespace DBparse {

  void DB_load(TString file, DBRequest request[], const int nvar){

    fstream file_csv; file_csv.open(file);
    cout<<"Attempting to load DB File"<<endl;
    cout<<"---------------------------------------------------------------"<<endl;
    //Check if file exists
    if (!file_csv.is_open()) {
      cout<<"!!!!WARNING, ERROR, DBparse::BD_load: File "<<file<<" is not found"<<endl;
      exit(0);
    }
    else {
      cout<<"DB file found: "<<file<<endl;
    }

    string line;
    int iline=0;

    bool var_found[nvar];   // Is the variable found?
    int var_index[nvar];       // Variable index, -1 means it is not found

    while (getline(file_csv, line)) {
      // Create a stringstream to parse the line
      stringstream ss(line);
      string cell;
      vector<string> row;
      iline++;   // Keep track of the number of lines
      // Split the line into cells using a comma as a delimiter
      while (getline(ss, cell, ',')) {
	row.push_back(cell);
      }

      // If we are on line 1 we look for the variable name
      // If the variable does not exist we check if it is mandatory
      if(iline == 1){ 
	for(int ivar=0; ivar < nvar; ivar++){
	  
	  var_found[ivar] = false;
	  var_index[ivar] = -1;   
	  
	  // Now loop over the DB and see if this variable is there
	  for(int irow=0; irow < row.size();irow++){
	      
	    if(irow == 0) continue;   //Ignore the run number cell
	    
	    if(row[irow] == request[ivar].name){
	      var_found[ivar] = true;
	      var_index[ivar] = irow;
	    }

	  }
	  
	  // Check if variable is found and mandatory
	  if(var_found[ivar] && request[ivar].mandatory) cout<<"DB var loaded: "<<request[ivar].name<<endl;
	  if(!var_found[ivar] && request[ivar].mandatory){
	    cout<<"!!!!WARNING, ERROR: Mandatory var not found: "<<request[ivar].name<<endl;
	    exit(0);
	  }
	}
      } // End checks on line one
      else { // Lines except line 1
      
	//Now inser the data into the map
	//This map is from run number to the variable in question
	for(int ivar=0; ivar < nvar; ivar++)
	  if(var_found[ivar])
	    request[ivar].var_map->insert(std::make_pair(stoi(row[0]),stod(row[var_index[ivar]])));
	
      }
    }
  }
  
}
