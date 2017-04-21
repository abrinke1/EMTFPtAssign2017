
#include "../interface/MacroHelper.h" // PtAlgo class

void BookPtHist( PtAlgo& algo );

void BookEffHist( PtAlgo& algo );

void BookCountHist( PtAlgo& algo );

void BookTurnOnHist( PtAlgo& algo, const int iPt, const int pt_cut );

void BookRateHist( PtAlgo& algo, const int iEff, const int eff_cut);

void LoopOverEvents( PtAlgo& algo, const TString tr_te );
