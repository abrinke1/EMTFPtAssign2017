
#include "../interface/TrackBuilder.h"

void BuildTracks( std::vector< std::array<int, 4> >& trks_hits,  // Vector of tracks, with hit indices by station
		  std::vector< std::array<int, 5> >& trks_modes, // Mode, CSC mode, RPC mode, and sumAbsDPhi/Theta of tracks
		  const std::array< std::array< std::vector<int>, 4>, 12> id, // All hit index values, by sector and station
		  const std::array< std::array< std::vector<int>, 4>, 12> ph, // All full-precision integer phi values
		  const std::array< std::array< std::vector<int>, 4>, 12> th, // All full-precision integer theta values
		  const std::array< std::array< std::vector<int>, 4>, 12> dt, // All detector values (0 for none, 1 for CSC, 2 for RPC)
		  const int mode,            // Mode of track we're building
		  const int maxRPC,          // Maximum # of stations with RPC hits
		  const int minCSC,          // Minimum # of stations with CSC hits
		  const int max_dPh=1024,    // Maximum dPhi between any two hits
		  const int max_dTh=8        // Maximum dTheta between any two hits
		  ) {

  trks_hits.clear();
  trks_modes.clear();

  std::array<int, 4> phs, ths, dts; // Phi, theta, and detector by station
  for (int i = 0; i < 4; i++) {
    phs.at(i) = -99;
    ths.at(i) = -99;
    dts.at(i) =   0;
  }

  // Loop over the sectors
  for (UInt_t iSc = 0; iSc < 12; iSc++) {
    // std::cout << "  * Looking at sector " << iSc << std::endl;

    // Check for consistency
    UInt_t nHitsTot = 0;
    for (UInt_t iSt = 0; iSt < 4; iSt++) {
      UInt_t nHits = id.at(iSc).at(iSt).size();
      assert( ph.at(iSc).at(iSt).size() == nHits && 
	      th.at(iSc).at(iSt).size() == nHits && 
	      dt.at(iSc).at(iSt).size() == nHits );
      nHitsTot += nHits;
    }
    if (nHitsTot == 0) 
      continue;

    // Reset the hits
    for (int i = 0; i < 4; i++) {
      phs.at(i) = -99;
      ths.at(i) = -99;
      dts.at(i) = 0;
    }

    // Create new vectors of tracks for this sector
    std::vector< std::array<int, 4> > s_trks_hits;
    std::vector< std::array<int, 5> > s_trks_modes;

    // Loop over station 1 hits
    for (UInt_t i1 = 0; i1 < max(int(id.at(iSc).at(0).size()), 1); i1++) {
      if (mode >= 8 && id.at(iSc).at(0).size() > 0) {
	phs.at(0) = ph.at(iSc).at(0).at(i1);
	ths.at(0) = th.at(iSc).at(0).at(i1);
	dts.at(0) = dt.at(iSc).at(0).at(i1);
	// std::cout << "\n    - Found a station 1 hit with phi = " << phs.at(0) << ", theta = " << ths.at(0) << std::endl;
      }
      
      // Loop over station 2 hits
      for (UInt_t i2 = 0; i2 < max(int(id.at(iSc).at(1).size()), 1); i2++) {
	if ( (mode % 8) / 4 > 0 && id.at(iSc).at(1).size() > 0) {
	  phs.at(1) = ph.at(iSc).at(1).at(i2);
	  ths.at(1) = th.at(iSc).at(1).at(i2);
	  dts.at(1) = dt.at(iSc).at(1).at(i2);
	  // std::cout << "    - Found a station 2 hit with phi = " << phs.at(1) << ", theta = " << ths.at(1) << std::endl;
	}
	
	// Loop over station 3 hits
	for (UInt_t i3 = 0; i3 < max(int(id.at(iSc).at(2).size()), 1); i3++) {
	  if ( (mode % 4) / 2 > 0 && id.at(iSc).at(2).size() > 0) {
	    phs.at(2) = ph.at(iSc).at(2).at(i3);
	    ths.at(2) = th.at(iSc).at(2).at(i3);
	    dts.at(2) = dt.at(iSc).at(2).at(i3);
	    // std::cout << "    - Found a station 3 hit with phi = " << phs.at(2) << ", theta = " << ths.at(2) << std::endl;
	  }
	  
	  // Loop over station 4 hits
	  for (UInt_t i4 = 0; i4 < max(int(id.at(iSc).at(3).size()), 1); i4++) {
	    if ( (mode % 2) > 0 && id.at(iSc).at(3).size() > 0) {
	      phs.at(3) = ph.at(iSc).at(3).at(i4);
	      ths.at(3) = th.at(iSc).at(3).at(i4);
	      dts.at(3) = dt.at(iSc).at(3).at(i4);
	      // std::cout << "    - Found a station 4 hit with phi = " << phs.at(3) << ", theta = " << ths.at(3) << std::endl;
	    }

	    int _mode, _mode_CSC, _mode_RPC;
	    int _sumAbsDPh, _sumAbsDTh;
	    
	    BuiltTrackMode( _mode, _mode_CSC, _mode_RPC, _sumAbsDPh, _sumAbsDTh,
			    phs, ths, dts, maxRPC, minCSC, max_dPh, max_dTh );
	    // std::cout << "    - BuiltTrackMode = " << _mode << ", CSC = " << _mode_CSC
	    // 	      << ", RPC = " << _mode_RPC << ", sumAbsDPh = " << _sumAbsDPh << std::endl;
	    // std::cout << "    - i1 = " << i1 << ", i2 = " << i2 << ", i3 = " << i3 << ", i4 = " << i4 << std::endl;
	    if (_mode == mode) {
	      std::array<int, 4> trk_hits;
	      std::array<int, 5> trk_modes;
	      
	      trk_hits.at(0) = (dts.at(0) > 0 ? id.at(iSc).at(0).at(i1) : -99);
	      trk_hits.at(1) = (dts.at(1) > 0 ? id.at(iSc).at(1).at(i2) : -99);
	      trk_hits.at(2) = (dts.at(2) > 0 ? id.at(iSc).at(2).at(i3) : -99);
	      trk_hits.at(3) = (dts.at(3) > 0 ? id.at(iSc).at(3).at(i4) : -99);

	      trk_modes.at(0) = (_mode);
	      trk_modes.at(1) = (_mode_CSC);
	      trk_modes.at(2) = (_mode_RPC);
	      trk_modes.at(3) = (_sumAbsDPh);
	      trk_modes.at(4) = (_sumAbsDTh);
	      
	      s_trks_hits.push_back(trk_hits);
	      s_trks_modes.push_back(trk_modes);
	    } // End conditional: if (_mode == mode)

	  } // End loop: for (UInt_t i4 = 0; i4 < nHits; i4++)
	} // End loop: for (UInt_t i3 = 0; i3 < nHits; i3++)
      } // End loop: for (UInt_t i2 = 0; i2 < nHits; i2++)
    } // End loop: for (UInt_t i1 = 0; i1 < nHits; i1+)

    // Select tracks with lowest sumAbsDPh, then sumAbsDTh, for a given CSC mode and RPC mode
    SelectTracks( s_trks_hits, s_trks_modes );
    // std::cout << "  * Selected " << s_trks_hits.size() << " tracks" << std::endl;
 
    trks_hits.insert( trks_hits.end(), s_trks_hits.begin(), s_trks_hits.end() );
    trks_modes.insert( trks_modes.end(), s_trks_modes.begin(), s_trks_modes.end() );

  } // End loop: for (UInt_t iSc = 0; iSc < 12; iSc++)

} // End function: void BuildTracks()



void BuiltTrackMode( int& mode, int& mode_CSC, int& mode_RPC, int& sumAbsDPh, int& sumAbsDTh,
		     const std::array<int, 4> phs, // Full-precision integer phi by station
		     const std::array<int, 4> ths, // Full-precision integer theta by station
		     const std::array<int, 4> dts, // Detector by station (0 for none, 1 for CSC, 2 for RPC)
		     const int maxRPC, const int minCSC,
		     const int max_dPh, const int max_dTh
		     ) {
  
  assert( phs.size() == 4 && ths.size() == 4 && dts.size() == 4 );
  
  mode      = 0;
  mode_CSC  = 0;
  mode_RPC  = 0;
  sumAbsDPh = 0;
  sumAbsDTh = 0;

  int nCSC = 0;
  int nRPC = 0;
  
  // Loop over hits from station 4 to station 2
  for (int iSt = 3; iSt > 0; iSt--) {
    if (dts.at(iSt) <= 0) continue;
    
    Bool_t pass_dPh = true;
    Bool_t pass_dTh = true;
    // Loop over other hits from station 3 to station 1
    for (int jSt = iSt - 1; jSt >= 0; jSt--) {
      if (dts.at(jSt) <= 0) continue;
      
      if ( abs(phs.at(iSt) - phs.at(jSt)) >= max_dPh ) pass_dPh = false;
      if ( abs(ths.at(iSt) - ths.at(jSt)) >  max_dTh ) pass_dTh = false;
      sumAbsDPh += abs(phs.at(iSt) - phs.at(jSt));
      sumAbsDTh += abs(ths.at(iSt) - ths.at(jSt));
    }

    if (pass_dPh && pass_dTh) {

      // Adjust track mode, including CSC/RPC mode
      mode += pow(2, 3 - iSt);
      if (dts.at(iSt) == 1) {
	nCSC += 1;
	mode_CSC += pow(2, 3 - iSt);
      } else if (dts.at(iSt) == 2) {
	nRPC += 1;
	mode_RPC += pow(2, 3 - iSt);
      }

    }
  }
  
  // Adjust modes for hit in station 1
  if (dts.at(0) >  0) mode += 8;
  if (dts.at(0) == 1) {
    nCSC += 1;
    mode_CSC += 8;
  } else if (dts.at(0) == 2) {
    nRPC += 1;
    mode_RPC += 8;
  }

  // std::cout << "nCSC = " << nCSC << ", nRPC = " << nRPC << std::endl;

  if (nCSC < minCSC) mode = 0;
  if (nRPC > maxRPC) mode = 0;
  
} // End function: int BuiltTrackMode()


void SelectTracks( std::vector< std::array<int, 4> >& s_trks_hits,  // Vector of tracks, with hit indices by station
		   std::vector< std::array<int, 5> >& s_trks_modes  // Mode, CSC mode, RPC mode, and sumAbsDPh/Th of tracks
		   ) {

  // Indices of selected tracks
  std::vector<int> trk_idxs;
  for (int iMode = 0; iMode < 16; iMode++) {  // Loop over CSC modes
    for (int jMode = 0; jMode < 16; jMode++) {  // Loop over RPC modes

      // Set default quantities
      int minSumAbsDPh = 999999;
      int minSumAbsDTh = 999999;
      int trk_idx = -99;

      // Loop over all built tracks
      for (int iTrk = 0; iTrk < s_trks_hits.size(); iTrk++) {

	assert(s_trks_hits.at(iTrk).size() == 4 && s_trks_modes.at(iTrk).size() == 5);

	if (s_trks_modes.at(iTrk).at(1) != iMode) continue;  // Check CSC mode
	if (s_trks_modes.at(iTrk).at(2) != jMode) continue;  // Check RPC mode

	if ( s_trks_modes.at(iTrk).at(3)  <  minSumAbsDPh ||
	     (s_trks_modes.at(iTrk).at(3) == minSumAbsDPh && 
	      s_trks_modes.at(iTrk).at(4) <  minSumAbsDTh) ) {
	  minSumAbsDPh = s_trks_modes.at(iTrk).at(3);
	  minSumAbsDTh = s_trks_modes.at(iTrk).at(4);
	  trk_idx = iTrk;
	}
      }
      // Save the best track for this CSC/RPC mode
      if (trk_idx >= 0) 
	trk_idxs.push_back(trk_idx);

    } // End loop over RPC modes
  } // End loop over CSC modes
  
  std::vector< std::array<int, 4> > s_trks_hits_sel;
  std::vector< std::array<int, 5> > s_trks_modes_sel;

  // Save the best tracks from all CSC/RPC modes
  for (int iIdx = 0; iIdx < trk_idxs.size(); iIdx++) {
    for (int iTrk = 0; iTrk < s_trks_hits.size(); iTrk++) {
      if (trk_idxs.at(iIdx) == iTrk) {
	s_trks_hits_sel.push_back( s_trks_hits.at(iTrk) );
	s_trks_modes_sel.push_back( s_trks_modes.at(iTrk) );
      }
    }
  }

  // Replace set of all tracks with selected tracks
  s_trks_hits  = s_trks_hits_sel;
  s_trks_modes = s_trks_modes_sel;

} // End function: void SelectTracks()
	





      


