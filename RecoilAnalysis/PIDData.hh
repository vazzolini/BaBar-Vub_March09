#ifndef PIDDATA
#define PIDDATA

#include "TString.h"
#include "TObject.h"
#include <iostream>

class PIDData: public TObject {

public:

  PIDData() {_pmax=_pmin=_tmin=_tmax=_fmin=_fmax=_e=_s=0.;_pass=_tot=0.;};
  PIDData(double p, double pmax, double th, double tmax, double phi, double fmax, 
	  double e = 0., double s = 0., double pass = 0., double tot = 0.);
  PIDData(const PIDData& orig);

  ~PIDData() {;}

  PIDData operator = (const PIDData &cmp); 
  Bool_t operator >= (PIDData &cmp) const;
  Bool_t operator < (const PIDData &cmp) const;
  Bool_t operator == ( const PIDData& cmp ) const;


  Bool_t isCell(double p, double t, double f);  
  void calcEffAndErr();
  void print();
  void increment(int pass, int tot, double w1 = 1., double w2 = 1.);
  void merge(PIDData other);
  Bool_t isZero();
  Bool_t sameValue(PIDData *, double epsilon = 1.e-3);
  
  double getPctr() const  {return (_pmax + _pmin)/2.;}
  double getTctr() const  {return (_tmax + _tmin)/2.;}
  double getFctr() const  {return (_fmax + _fmin)/2.;}

  double getPmax() const  {return _pmax;}
  double getTmax() const  {return _tmax;}
  double getFmax() const  {return _fmax;}

  double getPmin() const  {return _pmin;}
  double getTmin() const  {return _tmin;}
  double getFmin() const  {return _fmin;}

  double getE() const    {return _e;}
  double getS() const    {return _s;}
  double getPass() const {return _pass;}
  double getTot() const  {return _tot;}

  void setPmax(double p)     { _pmax = p;}
  void setTmax(double th)   { _tmax = th;}
  void setFmax(double phi) { _fmax = phi;}

  void setPmin(double p)   { _pmin = p;}
  void setTmin(double th)  { _tmin = th;}
  void setFmin(double phi) { _fmin = phi;}

  void setE(double e)    { _e = e;}
  void setS(double s)    { _s = s;}
  void setPass(double s) { _pass = s;}
  void setTot(double s)  { _tot = s;}

private:
  double _pmin, _pmax;
  double _tmin, _tmax;
  double _fmin, _fmax;
  double _e, _s;
  double _pass, _tot;

  ClassDef(PIDData,1) //Testing PIDData

};

ostream & operator << (ostream& , const PIDData&);


#endif
