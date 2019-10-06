//*****************************************************************
//   
//   Creation: Martin Simard & David Cote, Universite de Montreal, September 2006 
//   From: EvtGenModels/EvtHQET2FF.cc
//   Description: form factors for B->Dlnu according to HQET with dispersive FF
//   Use dispersion relation parametrization from 
//   I.Caprini, L.Lellouch, M.Neubert, Nucl. Phys. B 530,153(1998)
//
//   Please also see general comments in the XSLEvtFFWeight.hh file or BAD#809
//#include "BaBar/BaBar.hh"

#include "XSLBToDlnu_CLN.hh"

#include "XSLKin.hh"
#include "TLorentzVector.h"
#include <assert.h>
#include <iostream>

using std::cout;
using std::endl;
using std::string;

XSLBToDlnu_CLN::XSLBToDlnu_CLN(double mB, double mD, double q2, double theta_l, double rho2,double v1_1) : 
  XSLPseudoScalarFF(mB,mD,q2,theta_l) 
{ 
  _rho2=rho2;
  _v1_1=v1_1;
  SetNormalizations("Ignore");  
  Compute();  
}


XSLBToDlnu_CLN::XSLBToDlnu_CLN(TLorentzVector BLab, TLorentzVector LepLab, TLorentzVector DLab,double rho2,double v1_1) : 
  XSLPseudoScalarFF( BLab, LepLab, DLab) 
{ 
  _rho2=rho2;
  _v1_1=v1_1;
  SetNormalizations("Ignore");  
  Compute();  
}

XSLBToDlnu_CLN::XSLBToDlnu_CLN( XSLKin* DecayKin,double rho2,double v1_1 ): XSLPseudoScalarFF(DecayKin)
{ 
  _rho2=rho2;
  _v1_1=v1_1;
  SetNormalizations("Ignore");  
  Compute();  
}

XSLBToDlnu_CLN::~XSLBToDlnu_CLN()
{}


double 
XSLBToDlnu_CLN::GetFplus(double q2)
{
  //Note: _mXu=mD in this case
  double w = (_mB*_mB+_mXu*_mXu-q2)/(2.*_mB*_mXu);
  double z = (sqrt(w+1)-sqrt(2.))/(sqrt(w+1)+sqrt(2.));
  double v1 = _v1_1*(1.- 8.*_rho2*z + (51.*_rho2-10.)*z*z - (252.*_rho2-84.)*z*z*z);
  return v1;
}

void
XSLBToDlnu_CLN::SetNormalizations(string mode)
{
  //Note: The normalization constants were determined empirically from MC samples of 500k events
  //      at rho^2={0.63,0.81,0.99,1.08,1.17,1.26,1.35,1.53,1.71}   
  //      The ISGW2 normalization is obtained from a 3rd order polynomial fit to the points above.
  //      See: http://babar-hn.slac.stanford.edu:5090/hn/aux/cote/NormDlnuFF.eps
  //Code: XslFFReweighting/macros/FitNormalizationCurveDlnu.cxx
  //HFAG: rho^2=1.17+-0.18 (winter 06)
  _FLATQ2Normalization=1.0; 
  _PHSPNormalization=1.0; 
  if(_rho2==1.17 && _v1_1==1.0){ _ISGW2Normalization=1.0/0.4917; } //normalization point from 1M events.
  else{ _ISGW2Normalization=(0.73498 +1.103959*_rho2 -0.568666*_rho2*_rho2 +0.490898*_rho2*_rho2*_rho2); }

  if(_v1_1!=1.0){ cout<<"XSLBToDlnu_CLN: Unknown normalization for rho2="<<_rho2<<" and v1_1="<<_v1_1<<endl; exit(1); }
  return;
}

double XSLBToDlnu_CLN::FromPHSPToThisModel()
{ cout<<"XSLBToDlnu_CLN::FromPHSPToThisModel not implemented!   returning 1.0"<<endl; return 1.0; }

double XSLBToDlnu_CLN::FromFLATQ2ToThisModel()
{ cout<<"XSLBToDlnu_CLN::FromFLATQ2ToThisModel not implemented!   returning 1.0"<<endl; return 1.0; }
