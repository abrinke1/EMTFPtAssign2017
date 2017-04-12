
#include "../interface/MacroHelper.h"


/////////////////////////////////////////
// Compute median value of a 1D histogram
/////////////////////////////////////////
Float_t GetMedian(const TH1D* hist) {

  UInt_t  nBins    = hist->GetNbinsX();
  Float_t xMin     = hist->GetBinLowEdge(1);
  Float_t xMax     = hist->GetBinLowEdge(nBins) + hist->GetBinWidth(nBins);
  Float_t integral = hist->Integral();
  Float_t median   = -999;
  // std::cout << "Histogram with " << nBins << " bins from " << xMin << " to " << xMax << " has integral " << integral << std::endl;

  for (UInt_t iBin = 1; iBin <= nBins; iBin++) {
    if (hist->Integral(0, iBin) > 0.5*integral) {
      Float_t int_low = 0.5*integral - hist->Integral(0, iBin-1);
      Float_t int_hi  = hist->Integral(0, iBin) - 0.5*integral;
      assert( fabs( 1.0 - (int_low + int_hi) / hist->GetBinContent(iBin) ) < 0.001 );
      median = hist->GetBinLowEdge(iBin) + (int_low / (int_low + int_hi)) * hist->GetBinWidth(iBin);

      // std::cout << "At iBin = " << iBin << ", int_low = " << int_low << ", int_hi = " << int_hi << ", median = " << median << std::endl;
      // std::cout << "Low edge = " << hist->GetBinLowEdge(iBin) << ", high edge = " << hist->GetBinLowEdge(iBin) + hist->GetBinWidth(iBin)
      // 		<< ", content = " << hist->GetBinContent(iBin) << std::endl;
      break;
    }
  }

  return median;
}


///////////////////////////////////////////////////////////////////
// Reweight log2( calc pT / true pT ) histogram by resolution score
///////////////////////////////////////////////////////////////////
void WeightByResScore(TH1D& hist, const Float_t med_ratio) {
       
  UInt_t  nBins     = hist.GetNbinsX();

  for (UInt_t iBin = 1; iBin <= nBins; iBin++) {
    if (hist.GetBinContent(iBin) == 0) continue;

    Float_t ratio   = pow(2, hist.GetBinCenter(iBin));
    Float_t content = hist.GetBinContent(iBin);
    Float_t error   = hist.GetBinError(iBin) / content;
    
    Float_t med_ratio_1 = 1.0; // Set constant for now - AWB 20.03.17

    if (ratio > med_ratio_1)
      // hist.SetBinContent( iBin, content * (ratio / med_ratio_1) );
      hist.SetBinContent( iBin, content * pow((ratio / med_ratio_1), 3) );
    else
      // hist.SetBinContent( iBin, content * (med_ratio_1 / ratio) );
      hist.SetBinContent( iBin, content * pow((med_ratio_1 / ratio), 3) );
    
    hist.SetBinError( iBin, error * hist.GetBinContent(iBin) );
  }

}


////////////////////////////////////////////////////////////////////////
// Compute the resolution score of a log2( calc pT / true pT ) histogram
////////////////////////////////////////////////////////////////////////
Float_t GetResScore(const TH1D* hist, const Float_t med_ratio) {
  
  UInt_t  nBins     = hist->GetNbinsX();
  Float_t integral  = hist->Integral();
  Float_t score     = 0;

  Float_t med_ratio_1 = 1.0; // Set constant for now - AWB 20.03.17
    
  for (UInt_t iBin = 1; iBin <= nBins; iBin++) {
    Float_t ratio = pow(2, hist->GetBinCenter(iBin));
    if (ratio > med_ratio_1)
      // score += ( (ratio / med_ratio_1) * hist->GetBinContent(iBin) );
      score += ( pow((ratio / med_ratio_1), 3) * hist->GetBinContent(iBin) );
    else
      // score += ( (med_ratio_1 / ratio) * hist->GetBinContent(iBin) );
      score += ( pow((med_ratio_1 / ratio), 3) * hist->GetBinContent(iBin) );
  }

  return (score / integral);
}


///////////////////////////////////////////////////////////////////////////////////////////
// Compute the uncertainty on the resolution score of a log2( calc pT / true pT ) histogram
///////////////////////////////////////////////////////////////////////////////////////////
Float_t GetResScoreErr(const TH1D* hist, const Float_t med_ratio) {

  Float_t med_ratio_1 = 1.0; // Set constant for now - AWB 20.03.17

  TH1D* hist_wgt = (TH1D*) hist->Clone();
  WeightByResScore( (*hist_wgt), med_ratio_1 );

  Double_t error = -999;
  hist_wgt->IntegralAndError(0, hist_wgt->GetNbinsX(), error);

  return ( error / hist->Integral() );
}

