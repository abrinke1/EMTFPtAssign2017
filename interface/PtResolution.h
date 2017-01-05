
void BookScoreHist( const std::vector<std::tuple<TString, TString, float, TString, float>> pt_br, const int iBr,
                    const std::vector<std::tuple<TString, float, float, TString>> pt_bins,
                    const std::vector<std::tuple<TString, float, float, TString>> eta_bins,
                    TString& br_name, std::vector<TH2D*>& h_score );

void BookResHist( const std::vector<std::tuple<TString, TString, float, TString, float>> pt_br, const int iBr,
                  const std::vector<std::tuple<TString, float, float, TString>> pt_bins, const int iPt,
                  const std::vector<std::tuple<TString, float, float, TString>> eta_bins, const int iEta,
                  TString& br_name, std::vector< std::vector< std::vector<TH1D*> > >& h_res );

void LoopOverEvents( TChain* chain, const TString opt_str,
                     const std::vector<std::tuple<TString, float, float, TString>> pt_bins,
                     const std::vector<std::tuple<TString, float, float, TString>> eta_bins,
                     std::vector<std::tuple<TString, TString, float, TString, float>>& pt_br,
                     std::vector< std::vector< std::vector<TH1D*> > >& h_res );

void StoreMedians( const std::vector<std::tuple<TString, float, float, TString>> pt_bins,
                   const std::vector<std::tuple<TString, float, float, TString>> eta_bins,
                   const std::vector< std::vector< std::vector<TH1D*> > > h_res,
                   std::vector<std::tuple<TString, TString, float, TString, float>>& pt_br );

void FillScores( const std::vector<std::tuple<TString, TString, float, TString, float>> pt_br,
                 const std::vector<std::tuple<TString, float, float, TString>> pt_bins,
                 const std::vector<std::tuple<TString, float, float, TString>> eta_bins,
                 const int num_pt_bins, const int num_trg,
                 std::vector<float>& EMTF_vec_train,  std::vector<float>& EMTF_vec_train_err,
                 std::vector<float>& EMTF_vec_test,   std::vector<float>& EMTF_vec_test_err,
                 std::vector<float>& ratio_vec_train, std::vector<float>& ratio_vec_train_err,
                 std::vector<float>& ratio_vec_test,  std::vector<float>& ratio_vec_test_err,
                 const std::vector< std::vector< std::vector<TH1D*> > > h_res,
                 std::vector<TH2D*>& h_score );

void FillRatioGraphs( const std::vector<std::tuple<TString, TString, float, TString, float>> pt_br,
                      const std::vector<std::tuple<TString, float, float, TString>> pt_bins,
                      const std::vector<std::tuple<TString, float, float, TString>> eta_bins,
                      const int num_pt_bins, const int num_trg,
                      const std::vector<float> EMTF_vec_train,  const std::vector<float> EMTF_vec_train_err,
                      const std::vector<float> EMTF_vec_test,   const std::vector<float> EMTF_vec_test_err,
                      const std::vector<float> ratio_vec_train, const std::vector<float> ratio_vec_train_err,
                      const std::vector<float> ratio_vec_test,  const std::vector<float> ratio_vec_test_err );

void PrintRatios( const std::vector<std::tuple<TString, TString, float, TString, float>> pt_br,
                  const std::vector<std::tuple<TString, float, float, TString>> pt_bins,
                  const std::vector<std::tuple<TString, float, float, TString>> eta_bins,
                  const std::vector< std::vector< std::vector<TH1D*> > > h_res );

