#!/usr/bin/env python

#########################################################
###         Simple macro to read EMTF NTuples         ###
###            Andrew Brinkerhoff 10.12.16            ###
###                                                   ###
###  NTuple format in interface/PtLutInputBranches.hh ###
#########################################################

from ROOT import *

def main():

    in_file_name = '/afs/cern.ch/work/a/abrinke1/public/EMTF/Analyzer/ntuples/EMTF_MC_NTuple_SingleMu_RPC_300k.root'    
    in_file = TFile.Open(in_file_name)    
    in_tree = in_file.Get('ntuple/tree')

    # ## Can also use TChain
    # in_tree = TChain('ntuple/tree')
    # in_tree.Add('/afs/cern.ch/work/a/abrinke1/public/EMTF/Analyzer/ntuples/EMTF_MC_NTuple_SingleMu_noRPC_300k.root')
    
    print '\n******* About to enter the event loop *******' 
    for iEvt in range(in_tree.GetEntries()):

        if iEvt > 10: break
        if iEvt % 1 is 0: 
            print '\n*********************'
            print 'Looking at event %d' % iEvt
            print '*********************'

        in_tree.GetEntry(iEvt)

        muons = in_tree.GetBranch('muon')
        trks  = in_tree.GetBranch('track')
        hits  = in_tree.GetBranch('hit')

        nMuons = int(muons.GetLeaf('nMuons').GetValue())
        nTrks  = int(trks.GetLeaf('nTracks').GetValue())
        nHits  = int(hits.GetLeaf('nHits').GetValue())

        ## Print info for GEN level muons
        for iMu in range(nMuons):
            mu_pt  = muons.GetLeaf('pt').GetValue(iMu)
            mu_eta = muons.GetLeaf('eta').GetValue(iMu)
            mu_phi = muons.GetLeaf('phi').GetValue(iMu)
            print '\nMuon %d  has pT = %.1f, eta = %.2f, phi = %.2f' % ( iMu+1, mu_pt, mu_eta, mu_phi )
            
            
            ## Print info for hits in the same endcap
            for iHit in range(nHits):
                hit_station = int( hits.GetLeaf('station').GetValue(iHit) )
                hit_eta       = hits.GetLeaf('eta').GetValue(iHit)
                hit_phi       = hits.GetLeaf('phi').GetValue(iHit)

                if (hit_station < 0): continue
                if ( (hit_eta > 0) != (mu_eta > 0) ): continue
                print '    * Station %d hit has eta = %.2f, phi = %.2f' % ( hit_station, hit_eta, hit_phi )

            ## Print info for EMTF tracks in the same endcap
            for iTrk in range(nTrks):
                trk_sector = int( trks.GetLeaf('sector').GetValue(iTrk) )
                trk_mode   = int( trks.GetLeaf('mode').GetValue(iTrk) )
                trk_pt     = trks.GetLeaf('pt').GetValue(iTrk)
                trk_eta    = trks.GetLeaf('eta').GetValue(iTrk)
                trk_phi    = trks.GetLeaf('phi').GetValue(iTrk)

                if (trk_sector < 0): continue
                if ( (trk_eta > 0) != (mu_eta > 0) ): continue
                print '  * Sector %d track has pT = %.1f, eta = %.2f, phi = %.2f, mode = %d' % ( trk_sector, trk_pt, trk_eta, trk_phi, trk_mode )

                ## Print info for hits in the EMTF tracks
                for iTrkHit in range(4*iTrk + 4):
                    trk_hit_station = int( trks.GetLeaf('hit_station').GetValue(iTrkHit) )
                    trk_hit_eta       = trks.GetLeaf('hit_eta').GetValue(iTrkHit)
                    trk_hit_phi       = trks.GetLeaf('hit_phi').GetValue(iTrkHit)
                    
                    if (trk_hit_station < 0): continue
                    if ( (trk_hit_eta > 0) != (mu_eta > 0) ): continue
                    print '    * Station %d hit has eta = %.2f, phi = %.2f' % ( trk_hit_station, trk_hit_eta, trk_hit_phi )


    print '\n******* Leaving the event loop *******'
    
    print '\nExiting main() in ReadNTuple.py\n'

if __name__ == '__main__':
    main()
