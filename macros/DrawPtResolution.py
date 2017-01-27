#! /usr/bin/env python

#########################################################
###          Macro to draw pT resolution and          ###
###       resolution score # score ratio plots        ###
###           Andrew Brinkerhoff 05.01.17             ###
#########################################################

import ROOT
from collections import OrderedDict

def main():

    print 'Inside DrawPtResolution.py'
    ROOT.gROOT.SetBatch(ROOT.kTRUE)

###################
## Initialize files
###################

    print 'Accessing files'

    dir_name = '/afs/cern.ch/user/a/abrinke1/TMVA/EMTFPtAssign2017/plots/'
    # file_name = dir_name+'PtResolution_AWB_v1_17_01_23_400_trees_0p002_node.root'
    # file_name = dir_name+'PtResolution_AWB_v1_17_01_24_vars_all.root'
    # file_name = dir_name+'PtResolution_AWB_v1_17_01_24_vars_best.root'
    # file_name = dir_name+'PtResolution_AWB_v1_17_01_24_FRs.root'
    file_name = dir_name+'PtResolution_AWB_v1_17_01_24_bends.root'

    in_file = ROOT.TFile.Open(file_name)

    out_file = ROOT.TFile('plots/DrawPtResolution.root', 'recreate')
    out_dir = 'plots/png/'

##################
## Plot categories
##################

    # colors = [ROOT.kBlack, ROOT.kViolet+2, ROOT.kViolet-2, ROOT.kBlue, ROOT.kSpring, ROOT.kOrange, ROOT.kRed, ROOT.kMagenta]
    colors = [ROOT.kBlack, ROOT.kViolet, ROOT.kBlue, ROOT.kSpring, ROOT.kRed]

    weights = OrderedDict()
    weights['']     = [''           , 0.25] ## Plot central 99.5%
    weights['_wgt'] = [' (weighted)', 1.00] ## Plot central 98.0%

    pt_bins = OrderedDict()
    pt_bins['all']      = [  0, 1000, '1 < p_{T} < 1000 GeV']

    # pt_bins['1_8']      = [  1,    8, '1 < p_{T} < 8 GeV']
    # pt_bins['8_30']     = [  1,   30, '1 < p_{T} < 30 GeV']
    # pt_bins['30_120']   = [ 30,  120, '30 < p_{T} < 120 GeV']
    # pt_bins['120_1000'] = [250, 1000, '120 < p_{T} < 1000 GeV']

    pt_bins['1_4']      = [  1,    4, '1 < p_{T} < 4 GeV']
    pt_bins['4_8']      = [  4,    8, '4 < p_{T} < 8 GeV']
    pt_bins['8_15']     = [  8,   15, '8 < p_{T} < 15 GeV']
    pt_bins['15_30']    = [ 15,   30, '15 < p_{T} < 30 GeV']
    pt_bins['30_60']    = [ 30,   60, '30 < p_{T} < 60 GeV']
    pt_bins['60_120']   = [ 60,  120, '60 < p_{T} < 120 GeV']
    pt_bins['120_250']  = [120,  250, '120 < p_{T} < 250 GeV']
    pt_bins['250_1000'] = [250, 1000, '250 < p_{T} < 1000 GeV']

    eta_bins = OrderedDict()
    eta_bins['all']       = [1.20, 2.40, '1.2 < |#eta| < 2.4']
    # eta_bins['1p2_1p55']  = [1.20, 1.55, '1.2 < |#eta| < 1.55']
    # eta_bins['1p55_1p85'] = [1.55, 1.85, '1.55 < |#eta| < 1.85']
    # eta_bins['1p85_2p1']  = [1.85, 2.10, '1.85 < |#eta| < 2.1']
    # eta_bins['2p1_2p4']   = [2.10, 2.40, '2.1 < |#eta| < 2.4']

    Facts = OrderedDict()

    # Facts['f_0x001f01fd_0x2'] = ['d#phi12/23/34 + comb, #theta, FR1']
    Facts['f_0x001f01ff_0x4_invPt'] = ['log_{2}p_{T}, 1/p_{T} wgt / EMTF vars + comb + St1 ring']

    # # Facts['f_0x0000011d_0x2'      ] = ['1/p_{T} /  EMTF vars']
    # # Facts['f_0x0000011d_0x4'      ] = ['log_{2}p_{T}, EMTF vars']
    # Facts['f_0x001f01fd_0x2'      ] = ['1/p_{T} /  EMTF vars + comb']
    # Facts['f_0x001f01fd_0x4'      ] = ['log_{2}p_{T}, EMTF vars + comb']
    # # Facts['f_0x001fffff_0x2'      ] = ['1/p_{T} /  EMTF vars + comb + more']
    # # Facts['f_0x001fffff_0x4'      ] = ['log_{2}p_{T}, EMTF vars + comb + more']
    # # Facts['f_0x0000011d_0x2_invPt'] = ['1/p_{T} /  1/p_{T} wgt / EMTF vars']
    # # Facts['f_0x0000011d_0x4_invPt'] = ['log_{2}p_{T}, 1/p_{T} wgt / EMTF vars']
    # # Facts['f_0x001f01fd_0x2_invPt'] = ['1/p_{T} /  1/p_{T} wgt / EMTF vars + comb']
    # Facts['f_0x001f01fd_0x4_invPt'] = ['log_{2}p_{T}, 1/p_{T} wgt / EMTF vars + comb']
    # Facts['f_0x001fffff_0x2_invPt'] = ['1/p_{T} /  EMTF 1/p_{T} wgt / vars + comb + more']
    # # Facts['f_0x001fffff_0x4_invPt'] = ['log_{2}p_{T}, 1/p_{T} wgt / EMTF vars + comb + more']

    # Facts['f_0x00000004_0x4_invPt'] = ['d#phi12']
    # Facts['f_0x00000005_0x4_invPt'] = ['d#phi12, #theta']
    # Facts['f_0x0000000d_0x4_invPt'] = ['d#phi12/23, #theta']
    # Facts['f_0x00000085_0x4_invPt'] = ['d#phi12/24, #theta']
    # Facts['f_0x0000001d_0x4_invPt'] = ['d#phi12/23/34, #theta']
    # Facts['f_0x001f00fd_0x4_invPt'] = ['d#phi12/23/34 + comb, #theta']
    # Facts['f_0x001f0ffd_0x4_invPt'] = ['d#phi12/23/34 + comb, #theta, FRs']
    # Facts['f_0x001f0fff_0x4_invPt'] = ['d#phi12/23/34 + comb, #theta, FRs, St1 ring']
    # Facts['f_0x001fffff_0x4_invPt'] = ['d#phi12/23/34 + comb, #theta, FRs, St1 ring, bends']
    # Facts['f_0x8fff0fff_0x4_invPt'] = ['d#phi12/23/34 + comb, #theta, FRs, St1 ring, d#theta vars']

    # Facts['f_0x001f00ff_0x4_invPt'] = ['d#phi12/23/34 + comb, #theta, St1 ring']
    # Facts['f_0x001f01ff_0x4_invPt'] = ['d#phi12/23/34 + comb, #theta, St1 ring, FR1']
    # Facts['f_0x001f03ff_0x4_invPt'] = ['d#phi12/23/34 + comb, #theta, St1 ring, FR1/2']
    # # Facts['f_0x001f05ff_0x4_invPt'] = ['d#phi12/23/34 + comb, #theta, St1 ring, FR1/3']
    # # Facts['f_0x001f09ff_0x4_invPt'] = ['d#phi12/23/34 + comb, #theta, St1 ring, FR1/4']
    # Facts['f_0x001f0fff_0x4_invPt'] = ['d#phi12/23/34 + comb, #theta, St1 ring, all FRs']

    Facts['f_0x001f01ff_0x4_invPt'] = ['d#phi12/23/34 + comb, #theta, St1 ring, FR1']
    Facts['f_0x001f11ff_0x4_invPt'] = ['d#phi12/23/34 + comb, #theta, St1 ring, FR1, bend 1']
    Facts['f_0x001f31ff_0x4_invPt'] = ['d#phi12/23/34 + comb, #theta, St1 ring, FR1, bend 1/2']
    Facts['f_0x001f51ff_0x4_invPt'] = ['d#phi12/23/34 + comb, #theta, St1 ring, FR1, bend 1/3']
    Facts['f_0x001f91ff_0x4_invPt'] = ['d#phi12/23/34 + comb, #theta, St1 ring, FR1, bend 1/4']
    Facts['f_0x001ff1ff_0x4_invPt'] = ['d#phi12/23/34 + comb, #theta, St1 ring, FR1, all bends']


    MVAs = OrderedDict()
    MVAs['EMTF_pt']                = ['EMTF']

    MVAs['BDTG_AWB']               = ['BDT AWB']


#############
## Draw plots
#############

    # print 'Resolution histograms: looping over pT and eta bins'
    # for iWgt in weights.keys():
    #     for iPt in pt_bins.keys():
    #         for iEta in eta_bins.keys():
    #             DrawRes(in_file, out_dir, colors, iWgt, weights[iWgt], 
    #                     iPt, pt_bins[iPt], iEta, eta_bins[iEta], Facts, MVAs)
    #             DrawResRatio(in_file, out_dir, colors, iWgt, weights[iWgt], 
    #                          iPt, pt_bins[iPt], iEta, eta_bins[iEta], Facts, MVAs)

    print 'Score ratio graphs: looping over pT bins'
    # DrawRatioGraph(in_file, out_dir, '')
    for iPt in pt_bins.keys():
        DrawRatioGraph(in_file, out_dir, '_pt_%s' % iPt, Facts, MVAs)

    

def DrawRes(in_file, out_dir, colors, iWgt, weight, 
            iPt, pt_bin, iEta, eta_bin, Facts, MVAs):

    iStr = 'pt_'+iPt+'_eta_'+iEta
    c_res = ROOT.TCanvas('res_'+iStr)
    s_res_e = ROOT.THStack('res_'+iStr+'_EMTF', 'res_'+iStr+'_EMTF')
    s_res_m = ROOT.THStack('res_'+iStr+'_MVA',  'res_'+iStr+'_MVA' )
    
    # c_res.SetLeftMargin(0.08)
    c_res.SetRightMargin(0.04)
    # ## Doesn't actually work, despite helpful documentation here: https://root.cern.ch/phpBB3/viewtopic.php?t=4585
    # ROOT.gStyle.SetTitleOffset(0.5, 't') ## Doesn't have to be 't' - anything except 'x', 'y', or 'z'
    ROOT.gStyle.SetTitleOffset(1.02, 'x')
    ROOT.gStyle.SetTitleOffset(1.18, 'y')
    ROOT.gStyle.SetTitleX(0.53)
    ROOT.gStyle.SetTitleY(0.965)
    ROOT.gStyle.SetTitleW(0.9)
    ROOT.gStyle.SetTitleH(0.053)
    
    firstFact = next(iter(Facts))
    firstMVA = next(iter(MVAs))
    
    nFM = -1
    x_min = 999
    x_max = -999
    for iFact in Facts.keys():
        for iMVA in MVAs.keys():
            if ( iFact == firstFact and iMVA == firstMVA ): 
                ## Don't include factory name in EMTF histogram
                h_res = in_file.Get('h_res_'+iMVA+'_test_'+iStr+iWgt)
                [x_min, x_max] = CalcPercentile( h_res, weight[1], weight[1], x_min, x_max )
                if (abs(x_min) > abs(x_max)):
                    l_res = ROOT.TLegend(0.12, 0.66, 0.49, 0.885)
                else:
                    l_res = ROOT.TLegend(0.58, 0.66, 0.95, 0.885)
                l_res.SetMargin(0.15)  ## Default 0.25 (line length)
            elif ( iFact != firstFact and iMVA == firstMVA ):
                continue ## Don't draw subsequent EMTF histograms - they're identical
            else:
                h_res = in_file.Get('h_res_'+iFact+'_'+iMVA+'_test_'+iStr+iWgt)
                [x_min, x_max] = CalcPercentile( h_res, weight[1], weight[1], x_min, x_max )

            nFM += 1
            h_res.SetLineColor(colors[nFM % len(colors)])
            if (iMVA == firstMVA):
                l_res.AddEntry(h_res, MVAs[iMVA][0] )
            elif len(MVAs) < 3:
                l_res.AddEntry(h_res, Facts[iFact][0] )
            elif len(Facts) < 2:
                l_res.AddEntry(h_res, MVAs[iMVA][0] )
            else:
                l_res.AddEntry(h_res, Facts[iFact][0]+' / '+MVAs[iMVA][0] )
            
            if ( iFact == firstFact and iMVA == firstMVA ): 
                h_res.SetLineWidth(3)
                # t_str = Facts[iFact][0]+' '+iMVA+', '+pt_bin[2]+', '+eta_bin[2]
                t_str = 'p_{T} resolution, %s, %s' % (pt_bin[2], eta_bin[2])
                x_str = h_res.GetXaxis().GetTitle()
                s_res_e.Add(h_res)
            else:
                h_res.SetLineWidth(2)
                s_res_m.Add(h_res)
        ## End loop: for iMVA in MVAs.keys()
    ## End loop: for iFact in Facts.keys()

    s_res_e.SetMaximum( 1.1 * max( s_res_e.GetMaximum("nostack"), s_res_m.GetMaximum("nostack") ) )
    s_res_e.Draw("nostack")
    s_res_e.GetXaxis().SetLimits(x_min, x_max)
    s_res_e.SetTitle( t_str )
    s_res_e.GetXaxis().SetTitle( x_str )
    s_res_e.GetXaxis().SetTitleSize(0.04)
    s_res_e.GetXaxis().SetLabelSize(0.038)
    s_res_e.GetYaxis().SetTitle('Events'+weight[0])
    s_res_e.GetYaxis().SetTitleSize(0.038)
    s_res_m.Draw("histsamenostack")
    l_res.Draw()
    c_res.SaveAs(out_dir+c_res.GetName()+iWgt+'.png')

## End function DrawRes(in_file, out_dir, colors, iWgt, weight, 
##                      iPt, pt_bin, iEta, eta_bin, Facts, MVAs):


def DrawResRatio(in_file, out_dir, colors, iWgt, weight, 
                 iPt, pt_bin, iEta, eta_bin, Facts, MVAs):

    iStr = 'pt_'+iPt+'_eta_'+iEta
    c_res_rat = ROOT.TCanvas('res_rat_'+iStr)
    s_res_rat = ROOT.THStack('res_rat_'+iStr,  'res_rat_'+iStr )
    
    # c_res_rat.SetLeftMargin(0.08)
    c_res_rat.SetRightMargin(0.04)
    # ## Doesn't actually work, despite helpful documentation here: https://root.cern.ch/phpBB3/viewtopic.php?t=4585
    # ROOT.gStyle.SetTitleOffset(0.5, 't') ## Doesn't have to be 't' - anything except 'x', 'y', or 'z'
    ROOT.gStyle.SetTitleOffset(1.02, 'x')
    ROOT.gStyle.SetTitleOffset(1.18, 'y')
    ROOT.gStyle.SetTitleX(0.53)
    ROOT.gStyle.SetTitleY(0.965)
    ROOT.gStyle.SetTitleW(0.9)
    ROOT.gStyle.SetTitleH(0.053)
    
    firstFact = next(iter(Facts))
    firstMVA = next(iter(MVAs))
    
    nFM = 0
    x_min = 999
    x_max = -999
    for iFact in Facts.keys():
        for iMVA in MVAs.keys():
            if ( iFact == firstFact and iMVA == firstMVA ): 
                ## Don't include factory name in EMTF histogram
                h_res_e = in_file.Get('h_res_'+iMVA+'_test_'+iStr+iWgt)
                [x_min, x_max] = CalcPercentile( h_res_e, weight[1], weight[1], x_min, x_max )
                ## t_str = (h_res_e.GetTitle()).replace('EMTF pt (test) ','ratio to EMTF')
                t_str = 'p_{T} resolution ratio to EMTF, %s, %s' % ( pt_bin[2], eta_bin[2] )
                x_str = h_res_e.GetXaxis().GetTitle()
                if (abs(x_min) > abs(x_max)):
                    l_res_rat = ROOT.TLegend(0.12, 0.66, 0.49, 0.885)
                else:
                    l_res_rat = ROOT.TLegend(0.58, 0.66, 0.95, 0.885)
                l_res_rat.SetMargin(0.15)  ## Default 0.25 (line length)
                continue ## Don't draw; just save for ratio
            elif ( iFact != firstFact and iMVA == firstMVA ):
                continue ## Don't access subsequent EMTF histograms - they're identical
            else:
                h_res_m = in_file.Get('h_res_'+iFact+'_'+iMVA+'_test_'+iStr+iWgt)
                [x_min, x_max] = CalcPercentile( h_res_m, weight[1], weight[1], x_min, x_max )
                h_res_m.Divide(h_res_e)

            nFM += 1
            h_res_m.SetLineColor(colors[nFM % len(colors)])
            if (iMVA == firstMVA):
                l_res_rat.AddEntry(h_res_m, MVAs[iMVA][0] )
            elif len(MVAs) < 3:
                l_res_rat.AddEntry(h_res_m, Facts[iFact][0] )
            elif len(Facts) < 2:
                l_res_rat.AddEntry(h_res_m, MVAs[iMVA][0] )
            else:
                l_res_rat.AddEntry(h_res_m, Facts[iFact][0]+' / '+MVAs[iMVA][0] )

            h_res_m.SetLineWidth(2)
            s_res_rat.Add(h_res_m)

        ## End loop: for iMVA in MVAs.keys()
    ## End loop: for iFact in Facts.keys()

    s_res_rat.SetMaximum( min(10, 1.1*s_res_rat.GetMaximum("nostack")) )
    s_res_rat.Draw("nostack")
    s_res_rat.GetXaxis().SetLimits(x_min, x_max)
    s_res_rat.SetTitle( t_str )
    s_res_rat.GetXaxis().SetTitle( x_str )
    s_res_rat.GetXaxis().SetTitleSize(0.04)
    s_res_rat.GetXaxis().SetLabelSize(0.038)
    s_res_rat.GetYaxis().SetTitle('Ratio to EMTF'+weight[0])
    s_res_rat.GetYaxis().SetTitleSize(0.038)
    l_res_rat.Draw()
    c_res_rat.SaveAs(out_dir+c_res_rat.GetName()+iWgt+'.png')

## End function DrawResRatio(in_file, out_dir, colors, iWgt, weight, 
##                           iPt, pt_bin, iEta, eta_bin, Facts, MVAs):


def DrawRatioGraph(in_file, out_dir, pt_str, Facts, MVAs):

    ## Create canvas
    c_rat = ROOT.TCanvas('ratio_graph%s' % pt_str)
    # l_res = ROOT.TLegend(0.66, 0.71, 0.95, 0.885)
    
    ## Get ratio graphs
    h_rat_tr = in_file.Get('h_ratio_train%s' % pt_str)
    h_rat_te = in_file.Get('h_ratio_test%s' % pt_str)
    # l_res.AddEntry(h_res, iMVA )

    # t_str = (h_res.GetTitle()).replace('EMTF (test) ','')
    # x_str = h_res.GetXaxis().GetTitle()

    ## 1 unit per factory / MVA / pT bin
    x_min = min( min(h_rat_tr.GetX()), min(h_rat_te.GetX()) ) - 0.5
    x_max = max( max(h_rat_tr.GetX()), max(h_rat_te.GetX()) ) + 0.5 
    h_rat_tr.GetXaxis().SetRangeUser(x_min, x_max)

    ## Resize canvas
    c_rat.SetGrid()
    c_rat.SetCanvasSize  (800+32*int(x_max - x_min - 1.2), 600)  ## Default 700, 500
    c_rat.SetLeftMargin  (0.11 - 0.0018*(x_max - x_min - 1.2))   ## Default 0.10
    c_rat.SetRightMargin (0.30 - 0.006 *(x_max - x_min - 1.2))   ## Default 0.10
    c_rat.SetBottomMargin(0.28)
    # # ## Doesn't actually work, despite helpful documentation here: https://root.cern.ch/phpBB3/viewtopic.php?t=4585
    # # ROOT.gStyle.SetTitleOffset(0.5, 't') ## Doesn't have to be 't' - anything except 'x', 'y', or 'z'
    # ROOT.gStyle.SetTitleOffset(1.02, 'x')
    # ROOT.gStyle.SetTitleOffset(1.18, 'y')
    # ROOT.gStyle.SetTitleX(0.53)
    # ROOT.gStyle.SetTitleY(0.965)
    # ROOT.gStyle.SetTitleW(0.9)
    # ROOT.gStyle.SetTitleH(0.053)
    

    # t_str = (h_res.GetTitle()).replace('EMTF (test) ','')
    used_facts = []
    used_MVAs = []
    for iBin in range(1, h_rat_tr.GetXaxis().GetNbins()+1):
        for iFact in Facts.keys():
            if iFact in used_facts: continue
            for iMVA in MVAs.keys():
                if iMVA in used_MVAs: continue
                if len(MVAs) < 3:
                    if iFact in h_rat_tr.GetXaxis().GetBinLabel(iBin):
                        h_rat_tr.GetXaxis().SetBinLabel( iBin, Facts[iFact][0] )
                        used_facts.append(iFact)
                elif len(Facts) < 2:
                    if iMVA in h_rat_tr.GetXaxis().GetBinLabel(iBin):
                        h_rat_tr.GetXaxis().SetBinLabel( iBin, MVAs[iMVA][0] )
                        used_MVAs.append(iMVA)
                else:
                    h_rat_tr.GetXaxis().SetBinLabel( iBin, (h_rat_tr.GetXaxis().GetBinLabel(iBin)).replace(iFact, Facts[iFact][0]+' / ') )
                    h_rat_tr.GetXaxis().SetBinLabel( iBin, (h_rat_tr.GetXaxis().GetBinLabel(iBin)).replace(MVA, MVAs[iMVA][0]) )
    
    y_min = min( min(h_rat_tr.GetY()) - max(h_rat_tr.GetEY()), min(h_rat_te.GetY()) - max(h_rat_te.GetEY()) )
    y_max = max( max(h_rat_tr.GetY()) + max(h_rat_tr.GetEY()), max(h_rat_te.GetY()) + max(h_rat_te.GetEY()) )
    h_rat_tr.GetYaxis().SetRangeUser( y_min - 0.1*(y_max - y_min), y_max + 0.1*(y_max - y_min) )

    h_rat_tr.GetYaxis().SetTitleOffset((1.3 - 0.021*(x_max - x_min - 1.2))) ## Default 1.0
    h_rat_tr.GetXaxis().LabelsOption('d')  ## Should be able to change angle with "ChangeLabel", but doesn't work
    h_rat_tr.GetXaxis().SetLabelSize(0.04)
    h_rat_tr.GetXaxis().SetLabelOffset(0.01)
    if (pt_str == ''): h_rat_tr.GetXaxis().SetTitle('MVA  [p_{T} range]')
    else:              h_rat_tr.GetXaxis().SetTitle('MVA')
    h_rat_tr.GetXaxis().SetTitleSize(0.04)
    h_rat_tr.GetXaxis().SetTitleOffset(3.8)
    h_rat_tr.GetXaxis().SetTickLength(0.00001)
    
    h_grid = ROOT.TH2C('h_grid', '', 8, -5., 35., 6 , 0.5, 1.1);   
    h_grid.GetXaxis().SetNdivisions(108);
    h_grid.GetYaxis().SetNdivisions(106);
    h_grid.GetXaxis().SetTickLength(0.1)
    # hgrid->GetYaxis()->SetLabelOffset(999.);
    # hgrid->GetXaxis()->SetLabelOffset(999.);
    
    ROOT.gStyle.SetTitleAlign(23)
    ROOT.gStyle.SetTitleX(0.5*(1 + c_rat.GetLeftMargin() - c_rat.GetRightMargin()))
    h_rat_tr.Draw('APsame')
    h_rat_te.Draw('Psame')
    h_grid.Draw('same') ## Adjust y axis ranges?
    
    # h_rat_tr.SetTitle( t_str )
    # h_rat_tr.GetXaxis().SetTitle( x_str )
    # h_rat_tr.GetXaxis().SetLabelSize(0.038)
    # h_rat_tr.GetYaxis().SetTitle('Events (weighted)')
    # h_rat_tr.GetYaxis().SetTitleSize(0.038)    
    
    # l_res.Draw()
    c_rat.SaveAs(out_dir+c_rat.GetName()+'.png')
    
## End function: DrawRatioGraph(in_file, out_dir)


def CalcPercentile(hist, pct_low, pct_hi, x_min = 999999, x_max = -999999):
    for iBin in range(1, hist.GetNbinsX()+1):
        if ( hist.Integral(1, iBin) / hist.Integral() ) > (pct_low / 100.):
            x_min = min( x_min, hist.GetXaxis().GetBinLowEdge(iBin) )
        if ( hist.Integral(1, iBin) / hist.Integral() ) > 1. - (pct_hi / 100.):
            x_max = max( x_max, hist.GetXaxis().GetBinLowEdge(iBin+1) )
            break
        if iBin == hist.GetNbinsX():
            x_max = max( x_max, hist.GetXaxis().GetBinLowEdge(iBin) + hist.GetXaxis().GetBinWidth(iBin) )

    return [x_min, x_max]

## End function: CalcPercentile(hist, pct_low, pct_hi)


if __name__ == '__main__':
    main()
    
