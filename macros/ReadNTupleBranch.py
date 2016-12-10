#!/usr/bin/env python

#########################################################
###         Simple macro to read EMTF NTuples         ###
###            Andrew Brinkerhoff 10.12.16            ###
###                                                   ###
###  Load branch structure explicitly                 ###
###  NTuple format in interface/PtLutInputBranches.hh ###
#########################################################

from ROOT import *

def main():

    gROOT.ProcessLine("struct GenMuonBranch { int nMuons; float pt[2], eta[2], theta[2], phi[2]; int charge[2]; };")
    
    in_file_name = '/afs/cern.ch/work/a/abrinke1/public/EMTF/Analyzer/ntuples/EMTF_MC_NTuple_SingleMu_RPC_300k.root'    
    in_file = TFile.Open(in_file_name)    
    in_tree = in_file.Get('ntuple/tree')

    # ## Can also use TChain
    # in_tree = TChain('ntuple/tree')
    # in_tree.Add('/afs/cern.ch/work/a/abrinke1/public/EMTF/Analyzer/ntuples/EMTF_MC_NTuple_SingleMu_noRPC_300k.root')

    muon_br  = GenMuonBranch()
    in_tree.SetBranchAddress('muon', AddressOf(muon_br, 'nMuons') )    
    # hit_br   = EMTFHitBranch()
    # in_tree.SetBranchAddress('hit', AddressOf(hit_br, 'nHits') )    
    # track_br = EMTFTrackBranch()
    # in_tree.SetBranchAddress('track', AddressOf(track_br, 'nTracks') )    

    print '\n******* About to enter the event loop *******' 
    for iEvt in range(in_tree.GetEntries()):

        if iEvt > 10: break
        if iEvt % 1 is 0: 
            print '\n*********************'
            print 'Looking at event %d' % iEvt
            print '*********************'

        in_tree.GetEntry(iEvt)

        nMuons = muon_br.nMuons
        # nTrks  = int(track_br.nTracks)
        # nHits  = int(hit_br.nHits)

        print 'nMuons = %d' % nMuons

        ## Print info for GEN level muons
        for iMu in range(nMuons):
            mu_pt  = muon_br.pt[iMu]
            mu_eta = muon_br.eta[iMu]
            mu_phi = muon_br.phi[iMu]
            print 'Muon %d  has pT = %.1f, eta = %.2f, phi = %.2f' % ( iMu+1, mu_pt, mu_eta, mu_phi )


    print '\n******* Leaving the event loop *******'
    
    print '\nExiting main() in ReadNTuple.py\n'

if __name__ == '__main__':
    main()
