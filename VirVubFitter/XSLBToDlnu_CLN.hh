//*****************************************************************
//   
//   Creation: Martin Simard & David Cote, Universite de Montreal, September 2006 
//   From: EvtGenModels/EvtHQET2FF.cc
//   Description: form factors for B->Dlnu according to HQET with dispersive FF
//   Use dispersion relation parametrization from 
//   I.Caprini, L.Lellouch, M.Neubert, Nucl. Phys. B 530,153(1998)
//
//   Please also see general comments in the XSLEvtFFWeight.hh file or BAD#809

#ifndef XSLBTODLNU_CLN
#define XSLBTODLNU_CLN

#include "XSLPseudoScalarFF.hh"

#include "TLorentzVector.h"
class XSLKin;

class XSLBToDlnu_CLN : public XSLPseudoScalarFF {
public:
  XSLBToDlnu_CLN(double mB, double mD, double q2, double theta_l,double rho2=1.17,double v1_1=1.0);
  XSLBToDlnu_CLN(TLorentzVector BLab, TLorentzVector LepLab, TLorentzVector DLab,double rho2=1.17,double v1_1=1.0);
  XSLBToDlnu_CLN( XSLKin* DecayKin, double rho2=1.17,double v1_1=1.0);
  virtual ~XSLBToDlnu_CLN();
  
  virtual double FromPHSPToThisModel();
  virtual double FromFLATQ2ToThisModel();

protected:
  virtual double GetFplus(double q2);
  virtual void SetNormalizations(std::string mode);

  double _rho2;
  double _v1_1;
};

#endif



