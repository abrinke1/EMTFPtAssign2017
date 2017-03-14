
void BookScoreHist( const int iFM, const int iMVA,
                    const std::vector<std::tuple<TString, float, float, TString>> pt_bins,
                    const std::vector<std::tuple<TString, float, float, TString>> eta_bins,
                    const TString ft_name, const TString mva_name, 
		    std::vector< std::vector< std::pair<TH2D*, TH2D*> > >& h_score );

void BookResHist( const int iFM, const int iMVA, const int iPt, const int iEta,
                  const std::vector<std::tuple<TString, float, float, TString>> pt_bins, 
                  const std::vector<std::tuple<TString, float, float, TString>> eta_bins, 
                  const TString ft_name, const TString mva_name,
		  std::vector< std::vector< std::vector< std::vector< std::pair<TH1D*, TH1D*> > > > >& h_res );

void LoopOverEvents( TChain* chain, const TString ft_name, const TString trgPt, const int iFM, const TString tr_te,
                     const std::vector<std::tuple<TString, float, float, TString>> pt_bins,
                     const std::vector<std::tuple<TString, float, float, TString>> eta_bins,
                     std::vector<std::tuple<TString, float, float, float>>& MVAs,
                     std::vector< std::vector< std::vector< std::vector< std::pair<TH1D*, TH1D*> > > > >& h_res );

void StoreMedians( const TString ft_name, const int iFM, 
		   const std::vector<std::tuple<TString, float, float, TString>> pt_bins,
                   const std::vector<std::tuple<TString, float, float, TString>> eta_bins,
                   std::vector< std::vector< std::vector< std::vector< std::pair<TH1D*, TH1D*> > > > >& h_res,
		   std::vector<std::tuple<TString, float, float, float>>& MVAs );

void FillScores( const TString ft_name, const int iFM,
		 const std::vector<std::tuple<TString, float, float, float>> MVAs,
                 const std::vector<std::tuple<TString, float, float, TString>> pt_bins,
                 const std::vector<std::tuple<TString, float, float, TString>> eta_bins,
                 const int num_FM, const int num_pt, const int num_MVA,
                 std::vector<float>& EMTF_vec_tr,  std::vector<float>& EMTF_vec_tr_err,
                 std::vector<float>& EMTF_vec_te,  std::vector<float>& EMTF_vec_te_err,
                 std::vector<float>& ratio_vec_tr, std::vector<float>& ratio_vec_tr_err,
                 std::vector<float>& ratio_vec_te, std::vector<float>& ratio_vec_te_err,
                 const std::vector< std::vector< std::vector< std::vector< std::pair<TH1D*, TH1D*> > > > > h_res,
                 std::vector< std::vector< std::pair<TH2D*, TH2D*> > >& h_score );

void FillRatioGraphs( const std::vector<TString> fact_names, const std::vector<TString> MVA_names,
		      const std::vector<std::tuple<TString, float, float, TString>> pt_bins,
                      const std::vector<std::tuple<TString, float, float, TString>> eta_bins,
                      const int num_FM, const int num_pt, const int num_MVA,
                      const std::vector<float> EMTF_vec_tr,  const std::vector<float> EMTF_vec_tr_err,
                      const std::vector<float> EMTF_vec_te,  const std::vector<float> EMTF_vec_te_err,
                      const std::vector<float> ratio_vec_tr, const std::vector<float> ratio_vec_tr_err,
                      const std::vector<float> ratio_vec_te, const std::vector<float> ratio_vec_te_err );

void PrintRatios( const TString ft_name, const int iFM,
		  const std::vector<std::tuple<TString, float, float, float>> MVAs,
                  const std::vector<std::tuple<TString, float, float, TString>> pt_bins,
                  const std::vector<std::tuple<TString, float, float, TString>> eta_bins,
                  const std::vector< std::vector< std::vector< std::vector< std::pair<TH1D*, TH1D*> > > > > h_res );

