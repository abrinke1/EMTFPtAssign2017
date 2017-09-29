
// *** Default track-building settings *** //
int MIN_CSC = 2;  // Minimum number of CSC LCTs per track
int MAX_RPC = 2;  // Maximum number of RPC hits per track

namespace pTMulticlass_cfg {

  inline void ConfigureMode( const int MODE ) {
    
    std::cout << "\nRunning training for mode " << MODE << std::endl;
    
    // 4-station mode
    if (MODE == 15) {
      MIN_CSC = 3; 
      MAX_RPC = 1; // Saves time on training ... ideally should change back to 2? - AWB 02.06.17
    } 
    
    // 3-station modes
    if (MODE == 14 || MODE == 13 || MODE == 11 || MODE == 7) {
      MIN_CSC = 2;
      MAX_RPC = 1;
    }
    
    if (MODE == 14) {
    } 
    if (MODE == 13) {
    } 
    if (MODE == 11) {
    } 
    if (MODE ==  7) {
    } 
    
    // 2-station modes
    if (MODE == 12 || MODE == 10 || MODE == 9 || MODE == 6 || MODE == 5 || MODE == 3) {
      MIN_CSC = 1;
      MAX_RPC = 1;
    }
    
    if (MODE == 12) {
    } 
    if (MODE == 10) {
    } 
    if (MODE ==  9) {
    } 
    if (MODE ==  6) {
    } 
    if (MODE ==  5) {
    } 
    if (MODE ==  3) {
    } 
    
  }  // End function: inline void ConfigureMode( const int MODE )

} // End namespace PtRegression_Apr_2017_cfg
