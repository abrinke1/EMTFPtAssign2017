#!/usr/bin/env python

from ROOT import *
gROOT.ProcessLine("struct GenMuonBranch { int nMuons; float pt[2], eta[2], theta[2], phi[2]; int charge[2]; };")

chain = TChain("ntuple/tree")
chain.Add("/afs/cern.ch/work/a/abrinke1/public/EMTF/Analyzer/ntuples/EMTF_MC_NTuple_SingleMu_noRPC_300k.root")

muon_br = GenMuonBranch()

chain.SetBranchAddress("muon", AddressOf(muon_br, "nMuons") )

for i in xrange(chain.GetEntries()):
    if i > 10: break
    chain.GetEntry()
    print '\nPrinting event %d, with %d muons' % (i+1, muon_br.nMuons)
    for j in range(muon_br.nMuons):
        print 'Muon %d has pT = %f, eta = %f' % ( j+1, muon_br.pt[j], muon_br.eta[j] )

