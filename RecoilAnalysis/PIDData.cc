#include "PIDData.hh"

ClassImp(PIDData)


// ----------------------------------------------------------------------
PIDData::PIDData(double pmin, double pmax, double tmin, double tmax, double fmin, double fmax, 
		 double e, double s, double pass, double tot) {
  _pmin = pmin;
  _tmin = tmin;
  _fmin = fmin; 
  _pmax = pmax;
  _tmax = tmax;
  _fmax = fmax; 
  
  _e = e;
  _s = s;

  _pass = pass;
  _tot  = tot;

  // cout << "error: " << _e << " -> "; 
  // always recalculate! if (_e < 1.e-5) 
  calcEffAndErr();
  // cout << _e << endl;

}

// ----------------------------------------------------------------------
PIDData::PIDData(const PIDData& orig) {
  _pmin = orig.getPmin();
  _tmin = orig.getTmin();
  _fmin = orig.getFmin();

  _pmax = orig.getPmax();
  _tmax = orig.getTmax();
  _fmax = orig.getFmax();

  _e = orig.getE();
  _s = orig.getS();

  _pass = orig.getPass();
  _tot  = orig.getTot();

}


// ----------------------------------------------------------------------
PIDData PIDData::operator = (const PIDData &cmp) {
  setPmax(cmp.getPmax());
  setTmax(cmp.getTmax());
  setFmax(cmp.getFmax());
  
  setPmin(cmp.getPmin());
  setTmin(cmp.getTmin());
  setFmin(cmp.getFmin());
  
  setE(cmp.getE());
  setS(cmp.getS());

  setPass(cmp.getPass());
  setTot(cmp.getTot());
  
  return *this;
}


// ----------------------------------------------------------------------
Bool_t PIDData::operator >= (PIDData &cmp) const {
  Bool_t answer = (getPmax() >= cmp.getPmax() && getTmax() >= cmp.getTmax() 
		   && getFmax() >= cmp.getFmax());
  return answer;
}


// ----------------------------------------------------------------------
Bool_t PIDData::operator < (const PIDData &cmp) const {
  Bool_t answer(kFALSE);
  int anInt(0);
  if (getPmax() < cmp.getPmax()) {
    answer = kTRUE;
  } else if (getPmax() == cmp.getPmax()) {
    if (getTmax() < cmp.getTmax()) {
      answer = kTRUE;
    } else if (getTmax() == cmp.getTmax()) {
      if (getFmax() < cmp.getFmax()) {
	answer = kTRUE;
        } else {
          anInt += 1;
        }
    } else {
      anInt += 10;
    }
  } else {
    anInt += 100;
    }
  return answer;
}


// ----------------------------------------------------------------------
Bool_t PIDData::operator == (const PIDData &cmp) const {
  return ( getPmax() == cmp.getPmax() && getTmax() == cmp.getTmax() 
	   && getFmax() == cmp.getFmax() ); 
}


// ----------------------------------------------------------------------
void PIDData::increment(int pass, int total, double w1, double w2) {
  double p = w1*_pass + w2*pass;
  double t = w1*_tot +  w2*total;
  _pass = (int)p;
  _tot  = (int)t;  
  calcEffAndErr();
}
 

// ----------------------------------------------------------------------
void PIDData::calcEffAndErr() {
  
  if (_tot > 0.) {
    // _e = (_pass+1)/(_tot+2);
    _e = _pass/_tot;
    _s = sqrt(((_pass+1)*(_tot-_pass+1))/((_tot+3)*(_tot+2)*(_tot+2)));
  }
  else {
    _e = 0.;
    _s = 0.;
  }
}


// ----------------------------------------------------------------------
Bool_t PIDData::isCell(double p, double t, double f) {
  return ((getPmin() <= p) && (p < getPmax()) 
	  && (getTmin() <= t) && (t < getTmax()) 
	  && (getFmin() <= f) && (f < getFmax()));
}


// ----------------------------------------------------------------------
Bool_t PIDData::isZero() {
  return ((TMath::Abs(getE()) < 1.e-6)
	&& (TMath::Abs(getS()) < 1.e-6)
	&& (TMath::Abs(getPass()) < 1.e-6)
	&& (TMath::Abs(getTot()) < 1.e-6));
}


// ----------------------------------------------------------------------
Bool_t PIDData::sameValue(PIDData *cmp, double eps) {
  return ((TMath::Abs(getE()-cmp->getE()) < eps) 
	  && (TMath::Abs(getS()-cmp->getS()) < eps));
}


// ----------------------------------------------------------------------
void PIDData::print() {
  char line[200];
  sprintf(line, "%6.3f%6.3f %5.0f %5.0f %5.0f %5.0f %9.6f %9.6f %11.0f %11.0f", 
	  getPmin(), getPmax(), 
	  getTmin(), getTmax(), 
	  getFmin(), getFmax(), 
	  getE(),    getS(),
	  getPass(), getTot());
  std::cout << line << std::endl;
  
}


// ----------------------------------------------------------------------
void PIDData::merge(PIDData other) {
  _pmin = (other.getPmin() < getPmin() ? other.getPmin() : getPmin());
  _pmax = (other.getPmax() > getPmax() ? other.getPmax() : getPmax());
  _fmin = (other.getFmin() < getFmin() ? other.getFmin() : getFmin());
  _fmax = (other.getFmax() > getFmax() ? other.getFmax() : getFmax());
  _tmin = (other.getTmin() < getTmin() ? other.getTmin() : getTmin());
  _tmax = (other.getTmax() > getTmax() ? other.getTmax() : getTmax());

  _pass += other.getPass();
  _tot += other.getTot();
  calcEffAndErr();
}
  

// ===============================================================================
ostream & operator << (ostream& o, const PIDData& cmp) {
  char line[200];
  sprintf(line, "%6.3f%6.3f %5.0f %5.0f %5.0f %5.0f %9.6f %9.6f %11.0f %11.0f", 
	  cmp.getPmin(),cmp.getPmax(), 
	  cmp.getTmin(), cmp.getTmax(), 
	  cmp.getFmin(), cmp.getFmax(), 
	  cmp.getE(), cmp.getS(),
	  cmp.getPass(), cmp.getTot());
  return o << line;
  
}
