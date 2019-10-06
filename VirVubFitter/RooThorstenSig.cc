 /***************************************************************************** 
  * Project: RooFit                                                           * 
  *                                                                           * 
  * Copyright (c) 2000-2005, Regents of the University of California          * 
  *                          and Stanford University. All rights reserved.    * 
  *                                                                           * 
  * Redistribution and use in source and binary forms,                        * 
  * with or without modification, are permitted according to the terms        * 
  * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             * 
  *****************************************************************************/ 

 // -- CLASS DESCRIPTION [PDF] -- 
 // Thorsten Signal Definition as on BAD 1433 v. 6 p. 17

 #include <iostream> 

 #include "VirVubFitter/RooThorstenSig.hh" 
 #include "RooFitCore/RooAbsReal.hh" 

 ClassImp(RooThorstenSig) 

 RooThorstenSig::RooThorstenSig(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _r,
                        RooAbsReal& _sigma_r1,
                        RooAbsReal& _x_c,
                        RooAbsReal& _sigma_r2,
                        RooAbsReal& _sigma_l,
                        RooAbsReal& _n,
                        RooAbsReal& _alpha) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   r("r","r",this,_r),
   sigma_r1("sigma_r1","sigma_r1",this,_sigma_r1),
   x_c("x_c","x_c",this,_x_c),
   sigma_r2("sigma_r2","sigma_r2",this,_sigma_r2),
   sigma_l("sigma_l","sigma_l",this,_sigma_l),
   n("n","n",this,_n),
   alpha("alpha","alpha",this,_alpha)
 { 
 } 


 RooThorstenSig::RooThorstenSig(const RooThorstenSig& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   r("r",this,other.r),
   sigma_r1("sigma_r1",this,other.sigma_r1),
   x_c("x_c",this,other.x_c),
   sigma_r2("sigma_r2",this,other.sigma_r2),
   sigma_l("sigma_l",this,other.sigma_l),
   n("n",this,other.n),
   alpha("alpha",this,other.alpha)
 { 
 } 



 Double_t RooThorstenSig::evaluate() const 
 { 
   Double_t arg = x - x_c;  //ok
   Double_t e1  = arg/sigma_r1; //ok
   Double_t e2  = arg/sigma_r2; //ok
   
   if (x > x_c) {
     return //Need to be normalized?
       (r/sigma_r1)*( exp(-e1)/( (1+exp(-e1))*(1+exp(-e1)) ) )+(1-r)*exp(-0.5*e2*e2)/(sigma_r2); //Ok
   } else { // x<=x_c
     Double_t el = arg/sigma_l; //ok
     Double_t aa = -2*alpha*sigma_l/(1-exp(alpha)); //ok
     Double_t bb = exp(alpha)*pow((aa+alpha*sigma_l),n)/( (1+exp(alpha))*(1+exp(alpha)) ); //ok
     Double_t cc = 0.;

     // ------------> Compute C <------------------ 
     if ((alpha*sigma_l)<0) {
       cc = (r/(4*sigma_r1)+(1-r)/(sigma_r2))*pow(aa,n)/bb;
     } else {
       cc = (r/sigma_r1+(1-r)*4/(sigma_r2));
     }
     // ------------------------------------------ //
     if (x < x_c-alpha*sigma_l) {
       return (cc*bb)/pow((aa+x_c-x),n);
     } else {
       return cc*( exp(-el)/( (1+exp(-el))*(1+exp(-el)) ) );
     }
   }

   return 1.0 ; 
 } 



