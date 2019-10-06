#include "TRKData.hh"

ClassImp(TRKData)


// ----------------------------------------------------------------------
TRKData::TRKData() {
  _pmax=_pmin=_mmin=_mmax=_tmin=_tmax=_fmin=_fmax=_e=_s=0.;
}

// ----------------------------------------------------------------------
TRKData::TRKData(double pmin, double pmax, double mmin, double mmax, 
		 double tmin, double tmax, double fmin, double fmax, 
		 double e, double s) {
  _pmin = pmin;
  _mmin = mmin;
  _tmin = tmin;
  _fmin = fmin; 

  _pmax = pmax;
  _mmax = mmax;  
  _tmax = tmax;
  _fmax = fmax; 

  _e = e;
  _s = s;

}

// ----------------------------------------------------------------------
TRKData::TRKData(const TRKData& orig) {
  _pmin = orig.getPmin();
  _mmin = orig.getMmin();
  _tmin = orig.getTmin();
  _fmin = orig.getFmin();

  _pmax = orig.getPmax();
  _mmax = orig.getMmax();
  _tmax = orig.getTmax();
  _fmax = orig.getFmax();

  _e = orig.getE();
  _s = orig.getS();

}


// ----------------------------------------------------------------------
TRKData TRKData::operator = (const TRKData &cmp) {
  setPmax(cmp.getPmax());
  setMmax(cmp.getMmax());
  setTmax(cmp.getTmax());
  setFmax(cmp.getFmax());
  
  setPmin(cmp.getPmin());
  setMmin(cmp.getMmin());
  setTmin(cmp.getTmin());
  setFmin(cmp.getFmin());
  
  setE(cmp.getE());
  setS(cmp.getS());

  return *this;
}


// ----------------------------------------------------------------------
Bool_t TRKData::operator >= (TRKData &cmp) const {
  Bool_t answer = (getPmax() >= cmp.getPmax() && getMmax() >= cmp.getMmax() 
		   && getTmax() >= cmp.getTmax() && getFmax() >= cmp.getFmax());
  return answer;
}


// ----------------------------------------------------------------------
Bool_t TRKData::operator < (const TRKData &cmp) const {
  Bool_t answer(kFALSE);
  int anInt(0);
  if (getPmax() < cmp.getPmax()) {
    answer = kTRUE;
  } else if (getPmax() == cmp.getPmax()) {
    if (getMmax() < cmp.getMmax()) {
      answer = kTRUE; 
    } else if (getMmax() == cmp.getMmax()) {
      if (getTmax() < cmp.getTmax()) {
	answer = kTRUE;
      } else if (getTmax() == cmp.getTmax()) {
	if (getFmax() < cmp.getFmax()) {
	  answer = kTRUE;
        }
      }
    }
  }
  return answer;
}


// ----------------------------------------------------------------------
Bool_t TRKData::operator == (const TRKData &cmp) const {
  return ( getPmax() == cmp.getPmax() && getMmax() == cmp.getMmax() 
	   && getTmax() == cmp.getTmax() && getFmax() == cmp.getFmax() ); 
}


// ----------------------------------------------------------------------
Bool_t TRKData::isCell(double p, double m, double t, double f) {
  return ((getPmin() <= p) && (p < getPmax()) 
	  && (getMmin() <= m) && (m < getMmax()) 
	  && (getTmin() <= t) && (t < getTmax()) 
	  && (getFmin() <= f) && (f < getFmax()));
}


// ----------------------------------------------------------------------
void TRKData::print() {
  char line[200];
  sprintf(line, "%6.3f%6.3f %5.0f %5.0f %5.3f %5.3f %5.3f %5.3f %9.7f %9.7f", 
	  getPmin(), getPmax(), 
	  getMmin(), getMmax(), 
	  getTmin(), getTmax(), 
	  getFmin(), getFmax(), 
	  getE(),    getS());
  std::cout << line << std::endl;
  
}
 

// ===============================================================================
ostream & operator << (ostream& o, const TRKData& cmp) {
  char line[200];
  sprintf(line, "%6.3f %6.3f %5.0f %5.0f %5.3f %5.3f %5.3f %5.3f %9.7f %9.7f", 
	  cmp.getPmin(),cmp.getPmax(), 
	  cmp.getMmin(), cmp.getMmax(), 
	  cmp.getTmin(), cmp.getTmax(), 
	  cmp.getFmin(), cmp.getFmax(), 
	  cmp.getE(), cmp.getS());
  return o << line;
  
}
