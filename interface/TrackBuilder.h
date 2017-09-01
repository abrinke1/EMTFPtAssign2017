
void BuildTracks( std::vector< std::array<int, 4> >& trks_hits,  // Vector of tracks, with hit indices by station
		  std::vector< std::array<int, 5> >& trks_modes, // Mode, CSC mode, RPC mode, and sumAbsDPhi/Theta of tracks
		  const std::array< std::array< std::array<int, 4>, 4>, 12> id, // All hit index values, by sector and station
                  const std::array< std::array< std::array<int, 4>, 4>, 12> ph, // All full-precision integer phi values
                  const std::array< std::array< std::array<int, 4>, 4>, 12> th, // All full-precision integer theta values
                  const std::array< std::array< std::array<int, 4>, 4>, 12> dt, // All detector values (0 for none, 1 for CSC, 2 for RPC)
                  const int mode,            // Mode of track we're building
                  const int maxRPC = 0,      // Maximum # of stations with RPC hits
                  const int minCSC = 2,      // Minimum # of stations with CSC hits
		  const int max_dPh=1024,    // Maximum dPhi between any two hits
		  const int max_dTh=8        // Maximum dTheta between any two hits
                  );


void BuiltTrackMode( int& mode, int& mode_CSC, int& mode_RPC, 
		     int& sumAbsDPh, int& sumAbsDTh,
		     const std::array<int, 4> phs, // Full-precision integer phi by station
		     const std::array<int, 4> ths, // Full-precision integer theta by station
		     const std::array<int, 4> dts, // Detector by station (0 for none, 1 for CSC, 2 for RPC)
		     const int maxRPC=0,         // Maximum # of stations with RPC hits
		     const int minCSC=2,         // Minimum # of stations with CSC hits
		     const int max_dPh=1024,     // Maximum dPhi between any two hits
		     const int max_dTh=8         // Maximum dTheta between any two hits
		     );

void SelectTracks( std::vector< std::array<int, 4> >& s_trks_hits,  // Vector of tracks, with hit indices by station 
		   std::vector< std::array<int, 5> >& s_trks_modes  // Mode, CSC mode, RPC mode, and sumAbsDPhi/Theta of tracks
		   );
