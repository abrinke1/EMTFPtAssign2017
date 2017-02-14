
void BookPtHist( const int iFM, const int iMVA, const TString ft_name, const TString mva_name,
                 std::vector< std::vector<std::pair<TH2D*, TH2D*> > >& h_pt );

void BookEffHist( const int iFM, const int iMVA, const TString ft_name, const TString mva_name,
		  std::vector< std::vector<std::pair<TH2D*, TH2D*> > >& h_eff );

void BookCountHist( const int iFM, const int iMVA, const TString ft_name, const TString mva_name,
                    std::vector< std::vector<TH1D*> >& h_count );

void BookTurnOnHist( const int iFM, const int iMVA, const TString ft_name, const TString mva_name,
		     std::vector< std::vector<TH1D*> >& h_turn_on );

void BookRateHist( const int iFM, const int iMVA, const int iEff, 
		   const TString ft_name, const TString mva_name, const int eff_cut,
                   std::vector< std::vector< std::vector< std::pair<TH1D*, TH1D*> > > >& h_rate );

void LoopOverEvents( TChain* chain, const TString ft_name, const TString trgPt, const int iFM, const TString tr_te,
                     std::vector< std::tuple<TString, float> >& MVAs,
                     std::vector< std::vector< std::pair<TH2D*, TH2D*> > >& h_pt,
		     std::vector< std::vector<TH1D*> >& h_count );
