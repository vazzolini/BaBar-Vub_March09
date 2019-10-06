
#include "VirVubFitter/CMClass.hh"

// ############ CMClass Implementation #################

CMClass::CMClass() : ev("ntp1"),cm2(true),varfit(true),wfermivec() {}  // By Default CM2, varfit(mxhadfit,q2fit)

CMClass::CMClass(bool flag, bool f2,const std::vector<float>& wfv) : cm2(flag), varfit(f2), wfermivec(wfv)
{
  UpdateEv(flag);
}

CMClass::~CMClass() {}

void CMClass::SetCM(bool cm2choice)
{
  cm2=cm2choice;
  UpdateEv(cm2choice);
}

void CMClass::UpdateEv(bool iscm2)
{
  if(iscm2)
    ev="ntp1";
  else
    ev="events";
}

void CMClass::SetVF(bool choice)
{
  varfit=choice;
}
