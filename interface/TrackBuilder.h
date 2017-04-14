
void BuildTracks( std::vector< std::vector<int> >& trks_hits,  // Vector of tracks, with hit indices by station
		  std::vector< std::vector<int> >& trks_modes, // Mode, CSC mode, RPC mode, and sumAbsDPhi of tracks
		  const std::array< std::array< std::vector<int>, 4>, 12> id, // All hit index values, by sector and station
                  const std::array< std::array< std::vector<int>, 4>, 12> ph, // All full-precision integer phi values
                  const std::array< std::array< std::vector<int>, 4>, 12> th, // All full-precision integer theta values
                  const std::array< std::array< std::vector<int>, 4>, 12> dt, // All detector values (0 for none, 1 for CSC, 2 for RPC)
                  const int mode,            // Mode of track we're building
                  const int maxRPC = 0,      // Maximum # of stations with RPC hits
                  const int minCSC = 2,      // Minimum # of stations with CSC hits
		  const int max_dPh=1024,    // Maximum dPhi between any two hits
		  const int max_dTh=8        // Maximum dTheta between any two hits
                  );


void BuiltTrackMode( int& mode, int& mode_CSC, int& mode_RPC, 
		     int& sumAbsDPh, int& sumAbsDTh,
		     const std::vector<int> phs, // Full-precision integer phi by station
		     const std::vector<int> ths, // Full-precision integer theta by station
		     const std::vector<int> dts, // Detector by station (0 for none, 1 for CSC, 2 for RPC)
		     const int maxRPC=0,         // Maximum # of stations with RPC hits
		     const int minCSC=2,         // Minimum # of stations with CSC hits
		     const int max_dPh=1024,     // Maximum dPhi between any two hits
		     const int max_dTh=8         // Maximum dTheta between any two hits
		     );

void SelectTracks( std::vector< std::vector<int> >& s_trks_hits,  // Vector of tracks, with hit indices by station 
		   std::vector< std::vector<int> >& s_trks_modes  // Mode, CSC mode, RPC mode, and sumAbsDPhi of tracks
		   );
