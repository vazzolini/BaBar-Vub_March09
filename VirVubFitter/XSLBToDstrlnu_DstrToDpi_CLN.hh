//*****************************************************************
//   
//   Creation: Alexei Volk, TU-Dresden, 6 March 2008 
//   see comments in the XSLEvtFFWeight.hh file or BAD#809

#ifndef XSLBTODSTRLNU_DSTRTODPI_CLN
#define XSLBTODSTRLNU_DSTRTODPI_CLN

#include "XSLEvtFFWeight.hh"

#include <math.h>
#include <vector>
#include <string>
//#include "CLHEP/Vector/LorentzVector.h"
#include "TLorentzVector.h"
class XSLKin;

class XSLBToDstrlnu_DstrToDpi_CLN : public XSLEvtFFWeight {
public:
  XSLBToDstrlnu_DstrToDpi_CLN(double mB, double mDStar, double q2, double ctl, double ctv, double chi, double R1, double R2, double rho2);
  XSLBToDstrlnu_DstrToDpi_CLN(TLorentzVector BLab, TLorentzVector LepLab, TLorentzVector DStarLab, TLorentzVector DStarDaugLab, 
			    double R1, double R2, double rho2);
  XSLBToDstrlnu_DstrToDpi_CLN( XSLKin* DecayKin, double R1, double R2, double rho2 );
  virtual ~XSLBToDstrlnu_DstrToDpi_CLN();


  //Event Weight from specific Generators
  virtual double FromISGW2ToThisModel();
  virtual double FromPHSPToThisModel();
  virtual double FromFLATQ2ToThisModel();
  virtual double FromSP4ToThisModel();
  virtual double FromSP5ToThisModel(){ return FromSP4ToThisModel(); }
  virtual double FromSP6ToThisModel(){ return FromSP4ToThisModel(); }
  virtual double FromSP7ToThisModel(){ return FromSP8ToThisModel(); }
  virtual double FromSP8ToThisModel();

  //Differential Decay Rate
  virtual double dGammadQ2();
  virtual double dGammadQ2dCtl();
  virtual double dGammadQ2dCtv();
  virtual double dGammadQ2dCtldCtv();
  virtual double dGammadQ2dCtldCtvdChi();
  virtual double Gamma(){return _Gamma;};
  
  
  
protected:

  void InitWithKin( XSLKin* MyDecayKin );
 
  double IFF(double R1, double R2, double rho2,std::string model);
  double totRate(double R1, double R2, double rho2, std::string model);
  
  static double dGammadwLinearQ2(double w, const std::vector<double> &vars);
  static double dGammadwCLN(double w, const std::vector<double> &vars);
  
  //access 
  const double& R1()const {return _R1;}
  const double& R2()const {return _R2;}
  const double& rho2()const {return _rho2;}
 

  //angular terms
  inline double gPP(double ctl,double ctd,double chi)const
  {return (1-ctd*ctd)*(1-ctl)*(1-ctl);}

  inline double gMM(double ctl,double ctd,double chi)const
  {return (1-ctd*ctd)*(1+ctl)*(1+ctl);}

  inline double g00(double ctl,double ctd,double chi)const
  {return 4*ctd*ctd*(1-ctl*ctl);}

  inline double gPM(double ctl,double ctd,double chi)const
  {return -2*(1-ctd*ctd)*(1-ctl*ctl)*cos(2*chi);} //not the same as for LinearQ2 class

  inline double gP0(double ctl,double ctd,double chi)const
  { return -4*sqrt(1-ctd*ctd)*ctd*sqrt(1-ctl*ctl)*(1-ctl)*cos(chi); }

  inline double gM0(double ctl,double ctd,double chi)const 
  { return  4*sqrt(1-ctd*ctd)*ctd*sqrt(1-ctl*ctl)*(1+ctl)*cos(chi); }

  double _R1;
  double _R2;
  double _rho2;
  double _R1_SP4;
  double _R2_SP4;
  double _rho2_SP4;
  double _GammaSP4;
  
  double _R1_SP8;
  double _R2_SP8;
  double _rho2_SP8;
  double _GammaSP8;
  
  double _mDStar;
  double _mB;
  double _q2;
  double _w;
  double _ctl;
  double _ctv;
  double _chi;
  double _Gamma;

};

#endif





