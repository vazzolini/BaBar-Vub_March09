//*****************************************************************
//   
//   Creation: Alexei Volk, TU-Dresden, 6 March 2008
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
//   implemetiation of CLN parametrization based on BAD #1586

//#include "BaBar/BaBar.hh"

#include "XSLBToDstrlnu_DstrToDpi_CLN.hh"

#include "XSLItgSimpsonIntegrator.hh"
#include "XSLItgPtrFunction.hh"

#include "XSLKin.hh"

//#include "CLHEP/Vector/LorentzVector.h"
#include "TLorentzVector.h"
#include <assert.h>

#include <iostream> // Stream declarations
using std::cout;
using std::endl;
using std::string;

XSLBToDstrlnu_DstrToDpi_CLN::XSLBToDstrlnu_DstrToDpi_CLN( double mB, double mDStar, double q2, double ctl, double ctv, double chi, 
								    double R1, double R2, double rho2 ) : 
  XSLEvtFFWeight()
{
  _R1=R1;
  _R2=R2;
  _rho2=rho2;
  _R1_SP4=1.18;
  _R2_SP4=0.72;
  _rho2_SP4=0.92;
  _R1_SP8=1.33;
  _R2_SP8=0.92;
  _rho2_SP8=0.77;
    
  _mB=mB;
  _q2 = q2;
  _mDStar = mDStar;
  _w = (_mB*_mB+_mDStar*_mDStar-_q2)/(2*_mB*_mDStar);

  _ctl=ctl;
  _ctv=ctv;
  _chi=chi;
  
  // needed for calculation of the normalization factor Rn 
  _Gamma = totRate(R1, R2, rho2,"CLN");
  
  // alreday calculate values for SP4/SP8
  // you can calculate these values e.g. using totRate(R1_SP8, R2_SP8, rho2_SP8,"LinearQ2");
  // function
  _GammaSP8 = 0.883;
  _GammaSP4 = 0.857;
  
  
  
}

XSLBToDstrlnu_DstrToDpi_CLN::XSLBToDstrlnu_DstrToDpi_CLN( TLorentzVector BLab, TLorentzVector LepLab, 
								    TLorentzVector DStarLab, TLorentzVector DStarDaugLab,
								    double R1, double R2, double rho2 ) : 
  
  XSLEvtFFWeight()
{
  _R1=R1;
  _R2=R2;
  _rho2=rho2;
  
  _R1_SP4=1.18;
  _R2_SP4=0.72;
  _rho2_SP4=0.92;
  
  _R1_SP8=1.33;
  _R2_SP8=0.92;
  _rho2_SP8=0.77;
  
  XSLKin* MyDecayKin= new XSLKin(BLab,LepLab,DStarLab,DStarDaugLab);
  InitWithKin(MyDecayKin);  
  delete MyDecayKin;
  
  // needed for calculation of the normalization factor Rn 
  _Gamma = totRate(R1, R2, rho2,"CLN");
  
  // alreday calculate values for SP4/SP8
  // you can calculate these values e.g. using totRate(R1_SP8, R2_SP8, rho2_SP8,"LinearQ2");
  // function
  _GammaSP8 = 0.883;
  _GammaSP4 = 0.857;
  
}


XSLBToDstrlnu_DstrToDpi_CLN::XSLBToDstrlnu_DstrToDpi_CLN( XSLKin* DecayKin, double R1, double R2, double rho2 ): XSLEvtFFWeight()
{
  _R1=R1;
  _R2=R2;
  _rho2=rho2;
  
  _R1_SP4=1.18;
  _R2_SP4=0.72;
  _rho2_SP4=0.92;
  
  _R1_SP8=1.33;
  _R2_SP8=0.92;
  _rho2_SP8=0.77;
  

  //We haven't required the  constructor to accept only vector XSLKin objects,
  //but we're still gonna make sure it is the case!
  if(!DecayKin->isVector()){ 
    cout<<"An object of type XSLVectorKin is required by XSLVectorFF!!"<<endl;
    assert(0);
  }
  
  InitWithKin(DecayKin);
  // needed for calculation of the normalization factor Rn 
  _Gamma = totRate(R1, R2, rho2,"CLN");
  
  // alreday calculate values for SP4/SP8
  // you can calculate these values e.g. using totRate(R1_SP8, R2_SP8, rho2_SP8,"LinearQ2");
  // function
  _GammaSP8 = 0.883;
  _GammaSP4 = 0.857;
}

void
XSLBToDstrlnu_DstrToDpi_CLN::InitWithKin( XSLKin* MyDecayKin )
{
  _mB = MyDecayKin->BLab().Mag();
  _q2 = MyDecayKin->q2();
  _mDStar = MyDecayKin->XuLab().Mag();
  _w = MyDecayKin->w();
  
  _ctl=MyDecayKin->ctl();
  _ctv=MyDecayKin->ctv();
  _chi=MyDecayKin->chi();
  
  return;
}

XSLBToDstrlnu_DstrToDpi_CLN::~XSLBToDstrlnu_DstrToDpi_CLN()
{}

double 
XSLBToDstrlnu_DstrToDpi_CLN::FromSP4ToThisModel()
{
  // normalization factor
  double Rn(1); 
  
  // weight
  double w = IFF(_R1, _R2, _rho2,"CLN")/IFF(_R1_SP4, _R2_SP4, _rho2_SP4,"LinearQ2");
  
  // calculate new Gamma for SP4 or use already calculated value _GammaSP4
  double GammaSP4new = totRate(_R1_SP4, _R2_SP4, _rho2_SP4,"LinearQ2");
  
  Rn = GammaSP4new/_Gamma;
  w *= Rn;
  
  
  return w;
}

double 
XSLBToDstrlnu_DstrToDpi_CLN::FromSP8ToThisModel()
{
  // normalization factor
  double Rn(1);
  
  // weight
  double w = IFF(_R1, _R2, _rho2,"CLN")/IFF(_R1_SP8, _R2_SP8, _rho2_SP8,"LinearQ2");
  
  // calculate new Gamma for SP8 or use already calculated value _GammaSP8
  double GammaSP8new = totRate(_R1_SP8, _R2_SP8, _rho2_SP8,"LinearQ2");
  
  Rn = GammaSP8new/_Gamma;
  w *= Rn;
  
  return w;
}


double XSLBToDstrlnu_DstrToDpi_CLN::IFF(double R1, double R2, double rho2, string model)
{
  
  
  
  double w=_w;  double ctl=_ctl;  double ctd=_ctv;  double chi=_chi;
  double rate=0.0;
  
  
  double r=_mDStar/_mB;
  double r2=r*r;
  
  
  double hA1p;
  double R1w,R2w;
  
  double z = (sqrt(w+1.)-sqrt(2.))/(sqrt(w+1.)+sqrt(2.));
  
  
  
  if (model=="CLN") {
    // involve w-dependency for R1 and R2
    R1w = R1 - 0.12*(w - 1) + 0.05*(w - 1)*(w - 1);
    R2w = R2 + 0.11*(w - 1) - 0.06*(w - 1)*(w - 1);
    hA1p = 1. -8.*rho2*z +(53.*rho2 - 15.)*z*z -(231.*rho2 - 91.)*z*z*z;
  }
  else {
    // this is parametrization for LinearQ2 model (like in SP8)
    // be aware: rho2 here is no the same as for CLN parametrization
    hA1p=(1-rho2*(w-1));
    R1w = R1;
    R2w = R2;
  }
  
  double hpfac = (1-sqrt((w-1)/(w+1))*R1w);
  double hmfac=(1+sqrt((w-1)/(w+1))*R1w);
  double hzfac=(w - r) - (w - 1)*R2w ;
  
  double qfac = sqrt(1 - 2*w*r+r2);
  if(qfac<0) qfac=0;
  
  double hp = qfac*hpfac;
  double hm = qfac*hmfac;
  double hz=hzfac;
  
  double hphmTerm = hp*hm*gPM(ctl,ctd,chi);
  double hphzTerm = hp*hz*gP0(ctl,ctd,chi);
  double hmhzTerm = hm*hz*gM0(ctl,ctd,chi);
  double hp2Term  = hp*hp*gPP(ctl,ctd,chi);
  double hm2Term  = hm*hm*gMM(ctl,ctd,chi);
  double hz2Term  = hz*hz*g00(ctl,ctd,chi);
  
  rate=(hp2Term+hm2Term+hz2Term+hphmTerm+hphzTerm+hmhzTerm);
  rate*=hA1p*hA1p;
  
  
  return rate;
}


// need this function for numerical integration
double XSLBToDstrlnu_DstrToDpi_CLN::dGammadwCLN(double w, const std::vector<double> &vars)
{
  double R1 = vars[0];
  double R2 = vars[1];
  double rho2 = vars[2];
  double r = vars[3];
   
  double R1w = R1 - 0.12*(w - 1) + 0.05*(w - 1)*(w - 1);
  double R2w = R2 + 0.11*(w - 1) - 0.06*(w - 1)*(w - 1);
  double z=(sqrt(w+1.)-sqrt(2.))/(sqrt(w+1.)+sqrt(2.));
  double hA1p = 1. -8.*rho2*z +(53.*rho2 - 15.)*z*z -(231.*rho2 - 91.)*z*z*z; 
  
  double hpfac = (1-sqrt((w - 1)/(w + 1))*R1w);
  double hmfac = (1+sqrt((w - 1)/(w + 1))*R1w);
  double hzfac = (w - r) - (w - 1)*R2w ;
  
  double qfac = sqrt(1 - 2*w*r+r*r);
  if(qfac<0) qfac=0;
  
  double hp = qfac*hpfac;
  double hm = qfac*hmfac;
  double hz=hzfac;
  double IFFw = (w+1)*(w+1)*hA1p*hA1p*(hp*hp + hm*hm +hz*hz);
  
  return sqrt(w*w-1)*IFFw;
}

// need this function for numerical integration
double XSLBToDstrlnu_DstrToDpi_CLN::dGammadwLinearQ2(double w, const std::vector<double> &vars)
{
  double R1 = vars[0];
  double R2 = vars[1];
  double rho2 = vars[2];
  double r = vars[3];
  
  double hA1p = (1-rho2*(w - 1));
  
  double hpfac = (1-sqrt((w - 1)/(w + 1))*R1);
  double hmfac = (1+sqrt((w - 1)/(w + 1))*R1);
  double hzfac = (w - r) - (w - 1)*R2 ;
  
  double qfac = sqrt(1 - 2*w*r+r*r);
  if(qfac<0) qfac=0;
  
  double hp = qfac*hpfac;
  double hm = qfac*hmfac;
  double hz=hzfac;
  double IFFw = (w+1)*(w+1)*hA1p*hA1p*(hp*hp + hm*hm +hz*hz);
  
  return sqrt(w*w-1)*IFFw;
}


double XSLBToDstrlnu_DstrToDpi_CLN::totRate(double R1, double R2, double rho2, string model)
{
  
  // double mB = 5.279;
//   double mDStar = 2.007;
  std::vector<double> vars(4);
  // 0: R1
  // 1: R2
  // 2: rho2
  // 3: r
  vars[0] = R1;
  vars[1] = R2;
  vars[2] = rho2;
  vars[3] = _mDStar/_mB;
  
  int maxLoop = 20;
  
  // precision for intgral calculation
  double precision = 2.0e-2;
  
  double lowerlim = 1.;
  double upperlim = 1.504;
  double myintegral = 1;
  
  XSLItgPtrFunction *func;
  
  if (model=="CLN")
    func = new XSLItgPtrFunction(&dGammadwCLN, lowerlim, upperlim, vars);
  else func = new XSLItgPtrFunction(&dGammadwLinearQ2, lowerlim, upperlim, vars);
  
  XSLItgSimpsonIntegrator *integ = new XSLItgSimpsonIntegrator(*func, precision, maxLoop);
  myintegral = integ->evaluate(lowerlim, upperlim);
  delete integ;
  delete func;
  return myintegral;

}


//Event Weight from specific Generators
double XSLBToDstrlnu_DstrToDpi_CLN::FromISGW2ToThisModel()
{ cout<<"XSLBToDstrlnu_DstrToDpi_CLN::FromISGW2ToThisModel not implemented!   returning 1.0"<<endl; return 1.0; }

double XSLBToDstrlnu_DstrToDpi_CLN::FromPHSPToThisModel()
{ cout<<"XSLBToDstrlnu_DstrToDpi_CLN::FromPHSPToThisModel not implemented!   returning 1.0"<<endl; return 1.0; }

double XSLBToDstrlnu_DstrToDpi_CLN::FromFLATQ2ToThisModel()
{ cout<<"XSLBToDstrlnu_DstrToDpi_CLN::FromFLATQ2ToThisModel not implemented!   returning 1.0"<<endl; return 1.0; }

//Differential Decay Rates
double XSLBToDstrlnu_DstrToDpi_CLN::dGammadQ2()
{cout<<"XSLBToDstrlnu_DstrToDpi_CLN::dGammadQ2 not implemented!   returning 0"<<endl; return 0;}

double XSLBToDstrlnu_DstrToDpi_CLN::dGammadQ2dCtl()
{cout<<"XSLBToDstrlnu_DstrToDpi_CLN::dGammadQ2dCtl not implemented!   returning 0"<<endl; return 0;}

double XSLBToDstrlnu_DstrToDpi_CLN::dGammadQ2dCtv()
{cout<<"XSLBToDstrlnu_DstrToDpi_CLN::dGammadQ2dCtv not implemented!   returning 0"<<endl; return 0;}

double
XSLBToDstrlnu_DstrToDpi_CLN::dGammadQ2dCtldCtv()
{cout<<"XSLBToDstrlnu_DstrToDpi_CLN::dGammadQ2dCtldCtv not implemented!   returning 0"<<endl; return 0;}

//double XSLBToDstrlnu_DstrToDpi_CLN::Gamma()
//{cout<<"XSLBToDstrlnu_DstrToDpi_CLN::Gamma not implemented!   returning 0"<<endl; return 0;}

double
XSLBToDstrlnu_DstrToDpi_CLN::dGammadQ2dCtldCtvdChi()
{
  cout<<"XSLBToDstrlnu_DstrToDpi_CLN::dGammadQ2dCtldCtv not implemented!   returning 0"<<endl; return 0;
}












