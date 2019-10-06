#ifndef CMClass_hh
#define CMClass_hh

#include <vector>

//############### COMPUTING MODEL CLASS ###########

class CMClass
{
public:
  CMClass();      //def constr.
  CMClass(bool, bool, const std::vector<float>&);  //constr.
  virtual   ~CMClass();     //destr.

  void SetCM(bool);  //Set Computing Model
  void UpdateEv(bool); //Update Tree Var. dep on CM
  void SetVF(bool);  //Set Varfit

  char* ev;  //Tree Var.
  bool cm2;  //flag to control CM
  bool varfit; //flag to control mxhadfit and q2fit
  std::vector<float> wfermivec; //vector of wfermi's
};

#endif // #ifndef CMClass_hh
