#ifndef TRKDATA
#define TRKDATA

#include "TString.h"
#include "TObject.h"
#include <iostream>

class TRKData: public TObject {

public:

  TRKData();
  TRKData(double p, double pmax, double mult, double mmax, double th, double tmax,
	  double phi, double fmax, double e = 0., double s = 0.);
  TRKData(const TRKData& orig);

  ~TRKData() {;}

  TRKData operator = (const TRKData &cmp); 
  Bool_t operator >= (TRKData &cmp) const;
  Bool_t operator < (const TRKData &cmp) const;
  Bool_t operator == ( const TRKData& cmp ) const;


  Bool_t isCell(double p, double m, double t, double f);  
  void print();
  
  double getPctr() const  {return (_pmax + _pmin)/2.;}
  double getMctr() const  {return (_mmax + _mmin)/2.;}
  double getTctr() const  {return (_tmax + _tmin)/2.;}
  double getFctr() const  {return (_fmax + _fmin)/2.;}

  double getPmax() const  {return _pmax;}
  double getMmax() const  {return _mmax;}
  double getTmax() const  {return _tmax;}
  double getFmax() const  {return _fmax;}

  double getPmin() const  {return _pmin;}
  double getMmin() const  {return _mmin;}
  double getTmin() const  {return _tmin;}
  double getFmin() const  {return _fmin;}

  double getE() const    {return _e;}
  double getS() const    {return _s;}

  void setPmax(double p)     { _pmax = p;}
  void setMmax(double m)     { _mmax = m;}
  void setTmax(double th)   { _tmax = th;}
  void setFmax(double phi) { _fmax = phi;}

  void setPmin(double p)   { _pmin = p;}
  void setMmin(double m)  { _mmin = m;}
  void setTmin(double th)  { _tmin = th;}
  void setFmin(double phi) { _fmin = phi;}

  void setE(double e)    { _e = e;}
  void setS(double s)    { _s = s;}

private:
  double _pmin, _pmax;
  double _mmin, _mmax;
  double _tmin, _tmax;
  double _fmin, _fmax;
  double _e, _s;

  ClassDef(TRKData,1) //Testing TRKData

};

ostream & operator << (ostream& , const TRKData&);


#endif
