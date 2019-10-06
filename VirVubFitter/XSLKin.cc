//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   see comments in the XSLKin.hh file
// #include "BaBar/BaBar.hh"

#include "XSLKin.hh"


XSLKin::XSLKin( XSLKin* Kin ) :
  _BLab(Kin->BLab()),
  _LepLab(Kin->LepLab()),
  _XuLab(Kin->XuLab()),
  _XuDaugLab(Kin->XuDLab()),
  _isVector(Kin->isVector())
{ 
  Init(); 
  Compute();
}

XSLKin::XSLKin(TLorentzVector BLab, TLorentzVector LepLab, TLorentzVector XuLab) :
  _BLab(BLab),
  _LepLab(LepLab),
  _XuLab(XuLab)
{ 
  _isVector=false;
  _XuDaugLab=TLorentzVector(0,0,0,0);
  Init(); 
  Compute();
}


XSLKin::XSLKin(TLorentzVector BLab, TLorentzVector LepLab, TLorentzVector XuLab, TLorentzVector XuDaughterLab) :
  _BLab(BLab),
  _LepLab(LepLab),
  _XuLab(XuLab),
  _XuDaugLab(XuDaughterLab)
{ 
  _isVector=true;
  Init(); 
  Compute();
}

void
XSLKin::Init()
{
  _ctl=-666;
  _theta_l=-666;
  _q2=-666;
  _w=-666;
  _ctv=-666;
  _theta_v=-666;
  _chi=-666;
}


void 
XSLKin::Compute()
{
  //Xu meson in the B frame
  TLorentzVector XuB = _XuLab;
  XuB.Boost(-_BLab.BoostVector());

  //q2 computed as (B-Xu)^2 in the B frame
  _q2=_BLab.M2() + _XuLab.M2() - 2*_BLab.M()*XuB.E();
  _w=(_BLab.M2()+ _XuLab.M2()-_q2)/(2*_BLab.M()*_XuLab.M());

  //W boson in the Lab frame
  TLorentzVector WLab = _BLab - _XuLab;
  if(WLab.Mag2()<=0)
    {
      //this non-physical case is identified in the ntuple by q2<0 
      //To compute theta_l and chi however, we'll arbitrarly set the W mass to "almost zero" to avoid nan in the ntuple
      WLab.SetVectM(WLab.Vect(),0.000001);   
    }
  TLorentzVector W_B = WLab;
  W_B.Boost(-_BLab.BoostVector());

  //Lepton in the W frame
  TLorentzVector LepW = _LepLab;
  LepW.Boost(-_BLab.BoostVector());
  LepW.Boost(-W_B.BoostVector());

  _ctl=LepW.Vect().Unit()*W_B.Vect().Unit();
  _theta_l=LepW.Vect().Angle(W_B.Vect());


  if(_isVector){

    //Daugther of the Xu meson in Xu meson's frame
    TLorentzVector VD_Xu(_XuDaugLab);
    VD_Xu.Boost(-_BLab.BoostVector());
    VD_Xu.Boost(-XuB.BoostVector());
    
    //ctv
    _ctv=VD_Xu.Vect().Unit()*XuB.Vect().Unit();
    _theta_v=VD_Xu.Vect().Angle(XuB.Vect());
    
    //chi      
    //Art Snyder's (correct) 
    TVector3 p3LepTrans=LepW.Vect()-_ctl*LepW.Rho()*W_B.Vect().Unit();
    TVector3 p3VDTrans=VD_Xu.Vect()-_ctv*VD_Xu.Rho()*XuB.Vect().Unit();
    _chi=p3VDTrans.Angle(p3LepTrans);
    
    //Leif Wilden's  (wrong...)  --> deprecated!      
    //Hep3Vector LepPlane = LepW.vect().cross(W_B.vect());
    //Hep3Vector HadPlane = VD_Xu.vect().cross(XuB.vect());
    //_chib = LepPlane.angle(HadPlane);

  }//if(_isVector)

  return;
}


double
XSLKin::z(){
  double z=(sqrt(_w+1.)-sqrt(2.))/(sqrt(_w+1.)+sqrt(2.));
  return z;
}



