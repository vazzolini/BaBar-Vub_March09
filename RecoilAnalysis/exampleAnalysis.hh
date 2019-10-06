#ifndef EXAMPLEANALYSIS
#define EXAMPLEANALYSIS

#include "baseClass.hh"

// ----------------------------------------------------------------------
class exampleAnalysis: public baseClass {

public:

  exampleAnalysis(TTree *tree=0,int isMC =0, int newFormat = 2);

//redefine methods that exist in the baseClass as needed -- for example recoil()
//
//virtual void recoil();

//or define completely new methods as needed -- for example:
//void MyMethod();

private:

};

#endif
