#ifndef DBparse_H
#define DBparse_H


namespace DBparse {

  struct DBRequest{
    TString           name;      // Variable name, should match the DB file
    map<int,double>*  var_map;   // Data from the DB file
    TString           info;      // Extra info about this variable
    bool              mandatory; // 
  };

  void DB_load(TString file, DBRequest request[],int nvar);
  
}

#endif
