#!/usr/bin/env python

from ROOT import *

def main():

    file_name = '/afs/cern.ch/work/a/abrinke1/public/EMTF/Analyzer/ntuples/EMTF_MC_NTuple_SingleMu_RPC_300k.root'
    
    in_file = TFile.Open(file_name)
    
    tree = in_file.Get('ntuple/tree')
    
    # chain = TChain("ntuple/tree")
    # chain.Add("/afs/cern.ch/work/a/abrinke1/public/EMTF/Analyzer/ntuples/EMTF_MC_NTuple_SingleMu_noRPC_300k.root")

    print 'Entering event loop' 
    for iEvt in range(tree.GetEntries()):

        if iEvt > 10: break
        if iEvt % 1 is 0: print 'Event #', iEvt

        tree.GetEntry(iEvt)

        muons = tree.GetBranch('muon')
        trks  = tree.GetBranch('track')
        hits  = tree.GetBranch('hit')

        nMuons = int(muons.GetLeaf('nMuons').GetValue())
        nTrks  = int(trks.GetLeaf('nTracks').GetValue())
        nHits  = int(hits.GetLeaf('nHits').GetValue())

        for iMu in range(nMuons):
            mu_pt     = muons.GetLeaf('pt').GetValue(iMu)
            mu_eta    = muons.GetLeaf('eta').GetValue(iMu)
            mu_charge = muons.GetLeaf('charge').GetValue(iMu)

            print 'GEN muon %d: pT = %f, eta = %f, charge = %d' % ( iMu+1, mu_pt, mu_eta, mu_charge ) 


if __name__ == '__main__':
    main()
