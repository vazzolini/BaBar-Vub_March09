#include "recoilBuSys.hh"

#include <fstream>
#include <iostream>
#include <string>


recoilBuSys::recoilBuSys( const char * name, int seed ){

  // load input file
  readFile(name);
  TRandom rndm(seed);

  for (int j=0;j<7;j++){
    _bweights[j]=randomized(xv[j]/xt[j],xe[j],rndm,seed) ;
  }
}

float recoilBuSys::randomized(float mean,float sig, TRandom& rndm, int seed){
  if(seed == 0 ) return mean;// no smearing required
  float result=rndm.Gaus(mean,sig);
  if (result<0) result = 0;
  return result;
}

void recoilBuSys::readFile(const char *name){

  char  buffer[200];char  modeName[200];
  sprintf(buffer, "%s", name);
  ifstream is(buffer);
  
  for(int y=0; y<8; y++){
    xv[y] = xt[y] = xe[y] = 0;
  }

  while (is.getline(buffer, 200, '\n')) {

    float BrMode,PdgMode,PdgErrMode;
    if (buffer[0] == '#') {continue;}
    
    sscanf(buffer, "%f %f %f %s", &BrMode,&PdgMode,&PdgErrMode,modeName);
    int k = 0;    
    if (!strcmp(modeName, "pilnu")) k = 0;
    if (!strcmp(modeName, "rholnu")) k = 1;
    if (!strcmp(modeName, "omegalnu")) k = 2;
    if (!strcmp(modeName, "etalnu")) k = 3;
    if (!strcmp(modeName, "etaplnu")) k = 4;
    if (!strcmp(modeName, "a(0,1,..),blnu")) k = 5;
    if (!strcmp(modeName, "f1(0,..)h1ln")) k = 6;
    xt[k] = BrMode ; xv[k] = PdgMode; xe[k] = PdgErrMode;
    std::cout << " mode " << modeName << " " << k << " has a DECAY.DEC BR = " << xt[k] << " and a PDG BR = " << xv[k] << " +/- "  << xe[k] << std::endl;
  }
}



float 
recoilBuSys::weight(int bmode){
  float result;
  int imode(0);
  if(bmode==0)return 1;
  if(TMath::Abs(bmode)==11) imode = 0;
  else if(TMath::Abs(bmode)==13) imode = 1;
  else if(TMath::Abs(bmode)==14) imode = 2;
  else if(TMath::Abs(bmode)==12) imode = 3;
  else if(TMath::Abs(bmode)==15) imode = 4;
  else if(TMath::Abs(bmode)==16 || TMath::Abs(bmode)==17 ||TMath::Abs(bmode)==19 ) imode = 5;
  else if(TMath::Abs(bmode)==18 || TMath::Abs(bmode)==20 ||TMath::Abs(bmode)==21 || TMath::Abs(bmode)==22) imode = 6;
  else return 1;

  result =  _bweights[imode];
  if (result<-0.00001) std::cout <<" missing weight for bmode "<< bmode << std::endl;
  
  return result;
}
