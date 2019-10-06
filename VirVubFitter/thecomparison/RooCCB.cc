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

#include "RooCCB.hh" 

#include <iostream> 

#include "RooFitCore/RooAbsReal.hh" 

 ClassImp(RooCCB) 

 RooCCB::RooCCB(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _mean,
                        RooAbsReal& _sigma,
                        RooAbsReal& _alpha,
                        RooAbsReal& _n,
                        RooAbsReal& _endpoint) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   mean("mean","mean",this,_mean),
   sigma("sigma","sigma",this,_sigma),
   alpha("alpha","alpha",this,_alpha),
   n("n","n",this,_n),
   endpoint("endpoint","endpoint",this,_endpoint)
 { 
 } 


 RooCCB::RooCCB(const RooCCB& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   mean("mean",this,other.mean),
   sigma("sigma",this,other.sigma),
   alpha("alpha",this,other.alpha),
   n("n",this,other.n),
   endpoint("endpoint",this,other.endpoint)
 { 
 } 




Double_t RooCCB::evaluate() const 
{ 
  
  Double_t fcb=0;
  Double_t threshold=2*endpoint-x;
  
  if(x>endpoint)  
    return 0;
  
  else {
    if(x > (mean-alpha*sigma))
      fcb=exp(-0.5*((x-mean)*(x-mean))/(sigma*sigma));
    else
      fcb=exp(-0.5*alpha*alpha)*pow((n*sigma/alpha),n)*1/(pow((mean-sigma*alpha+n*sigma/alpha-x),n));
    if(threshold>mean-alpha*sigma)
      return fcb-exp(-500*(endpoint-x))*exp(-0.5*((threshold-mean)*(threshold-mean))/(sigma*sigma));
    else
      return fcb-exp(-500*(endpoint-x))*exp(-0.5*alpha*alpha)*pow((n*sigma/alpha),n)*1/(pow((mean-sigma*alpha+n*sigma/alpha-threshold),n));
  }
} 



