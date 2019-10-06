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
 // Your description goes here... 

 #include <iostream> 

 #include "VirVubFitter/RooGauz.hh" 

 ClassImp(RooGauz) 

 RooGauz::RooGauz(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _x_max,
                        RooAbsReal& _sigma1,
	                RooAbsReal& _sigma2,
	                RooAbsReal& _sigma3) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   x_max("x_max","x_max",this,_x_max),
   sigma1("sigma1","sigma1",this,_sigma1),
   sigma2("sigma2","sigma2",this,_sigma2),
   sigma3("sigma3","sigma3",this,_sigma3)
 { 
 } 


 RooGauz::RooGauz(const RooGauz& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   x_max("x_max",this,other.x_max),
   sigma1("sigma1",this,other.sigma1),
   sigma2("sigma2",this,other.sigma2),
   sigma3("sigma3",this,other.sigma3)
 { 
 } 



 Double_t RooGauz::evaluate() const 
 { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
   Double_t arg= x - x_max;  
   Double_t e=arg/sigma2;
   
   // cout<<"x "<<x<<" x_max "<<x_max<<"sigma1 "<<sigma1<<"eval "<<exp(-0.5*arg*arg/(sigma1*sigma1))+1/sigma1*exp(-e)/((1+exp(-e))*(1+exp(-e)))<<endl;
   // cout<<"x "<<x<<" x_max "<<x_max<<" sigma1 "<<sigma1<<" eval "<<exp(-0.5*arg*arg/(sigma2*sigma2))<<endl;
   if(x<x_max)
     return  exp(-0.5*arg*arg/(sigma1*sigma1))+1/sigma2*exp(-e)/((1+exp(-e))*(1+exp(-e)));
   else
     return  (1+1/(4*sigma2))*exp(-0.5*arg*arg/(sigma3*sigma3));
 } 



