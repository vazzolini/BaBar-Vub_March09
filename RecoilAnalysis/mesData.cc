#include "mesData.hh"

ClassImp(mesData)

mesData::mesData() {
  _name = " ";
  _sig =  _sigE = _bg = _bgE = _pu = _puE = 0;
  fSigma = fSigmaE = fMean = fMeanE = 0;
}

mesData::mesData(TString name, double sig, double sigE, double bg, double bgE, double pu, double puE){ 
  _name = name;
  _sig = sig;
  _sigE= sigE;
  _bg = bg; 
  _bgE = bgE;
  _pu = pu;
  _puE = puE;
  fSigma = fSigmaE = fMean = fMeanE = 0;

}

mesData::mesData(double sig) {
  _sig = sig; 
  _name = " ";
  _sigE = _bg = _bgE = _pu = _puE = 0;
  fSigma = fSigmaE = fMean = fMeanE = 0;
}

mesData::mesData(const mesData& pippo) {
  _name = pippo._name;
  _sig = pippo._sig;
  _sigE= pippo._sigE;
  _bg = pippo._bg; 
  _bgE = pippo._bgE;
  _pu = pippo._pu;
  _puE = pippo._puE;
  fSigma  = pippo.fSigma;
  fSigmaE = pippo.fSigmaE;
  fMean   = pippo.fMean;
  fMeanE  = pippo.fMeanE;
}

mesData& mesData::operator = (const mesData& pippo) {
  _name = pippo._name;
  _sig = pippo._sig;
  _sigE= pippo._sigE;
  _bg = pippo._bg; 
  _bgE = pippo._bgE;
  _pu = pippo._pu;
  _puE = pippo._puE;
  fSigma  = pippo.fSigma;
  fSigmaE = pippo.fSigmaE;
  fMean   = pippo.fMean;
  fMeanE  = pippo.fMeanE;
  return *this;
}

mesData::~mesData() {}

TString mesData::theName() { 
  return _name;
}

double mesData::theSig() { 
  return _sig;
}

double mesData::theErrSig() { 
  return _sigE;
}

double mesData::theBg() { 
  return _bg;
}

double mesData::theErrBg() { 
  return _bgE;
}

double mesData::thePu() { 
  return _pu;
}

double mesData::theErrPu() { 
  return _puE;
}
