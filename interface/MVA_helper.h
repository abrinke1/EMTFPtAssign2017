
// Defines single input variable to MVA in TMVA macro
class MVA_var {

 public:

  // Default constructor
  MVA_var( TString _name = "", TString _descr = "", TString _unit = "", 
	   char _type = 'F', Double_t _def_val = -99 ) {
    name    = _name;
    descr   = _descr;
    unit    = _unit;
    type    = _type;
    def_val = _def_val;
  } // End default constructor MVA_var()

  TString name;
  TString descr;
  TString unit;
  char type;
  Double_t def_val;
  
}; // End class MVA_var

