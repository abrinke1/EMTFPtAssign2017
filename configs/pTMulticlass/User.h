
// *** Default user settings *** //
TString OUT_DIR_NAME  = ".";  // Directory for output ROOT file
TString OUT_FILE_NAME = "pTMulticlass";  // Name base for output ROOT file
TString EOS_DIR_NAME  = "root://eoscms.cern.ch//store/user/abrinke1/EMTF/Emulator/ntuples";  // Input directory in eos

namespace pTMulticlass_cfg {
  
  inline void ConfigureUser( const TString USER ) {
    
    std::cout << "\nConfiguring pTMulticlass code for user " << USER << std::endl;
    
    if (USER == "WEI") {
      EOS_DIR_NAME = "root://eoscms.cern.ch//store/user/abrinke1/EMTF/Emulator/ntuples";  // Input directory in eos
      OUT_DIR_NAME = "/afs/cern.ch/work/w/wshi/public/TMVA_2017/TMVA/EMTFPtAssign2017/";
      OUT_FILE_NAME = "pTMulticlass";
    }
    
  } // End function: inline void ConfigureUser()    
  
} // End namespace pTMulticlass_cfg
