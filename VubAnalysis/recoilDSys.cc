#include "recoilDSys.hh"

#include <fstream.h>
#include <iostream.h>
#include <string.h>
#include <iomanip.h>

recoilDSys::recoilDSys(const char *name, int seed, int imode){
  // load input file
  readFile(name,seed,imode);
  // null all B weights
  for (int j=0;j<6;j++){
    _b0weights[j]=0;
    _bchweights[j]=0;
  }
    
}

recoilDSys::recoilDSys( int seed): _nmodes(0){

  TRandom rndm(seed);
  //  TRandom rndm(123);
  //  gRandom->SetSeed(123);
  //  for(int i=1; i<10*seed+1; i++)
  //    gRandom->Rndm(); // Effectively sets the seed.

  float xv0[6]={10.21,2.07,5.70,0.,0.52,0.23}; // existing measurements (order: all semilep, D,D*,nothing,D1*,D2*)
  float xt0[6]={10.4,2.10,5.60,9.,0.56,0.37};// values from decay.dec
  float xle0[6]={0.17,0.15,1.02,0.,0.15,0.23}; // errors on existing measurements downwards
  float xhe0[6]={0.17,0.15,0.53,0.,0.15,0.23}; // errors on existing measurements upwards

  // fill D*lnu numbers ( in %)
  double totbr0(xv0[0]);
  double totbr0t(xt0[0]);
  for (int j=1;j<6;j++){
    _b0weights[j]=Brandomized(xv0[j]/xt0[j],xle0[j]/xt0[j],xhe0[j]/xt0[j],rndm,seed) ;
    if(j!=3){
      totbr0 -=_b0weights[j]*xt0[j];
      totbr0t-=xt0[j];
    }
  }
  _b0weights[0]=Brandomized(totbr0/totbr0t,xle0[0]/totbr0t,xhe0[0]/totbr0t,rndm,seed) ;// this is the weight on the D** non resonant only

  float xvch[6]={11.04,2.24,6.17,0.,0.56,0.30};
  float xtch[6]={10.4,2.10,5.60,9.,0.56,0.37};
  float xlech[6]={0.18,0.16,1.13,0.,0.16,0.3};
  float xhech[6]={0.18,0.15,0.83,0.,0.16,0.3};

  double totbrch(xvch[0]);
  double totbrcht(xtch[0]);
  for (int j=1;j<6;j++){
    _bchweights[j]=Brandomized(xvch[j]/xtch[j],xlech[j]/xtch[j],xhech[j]/xtch[j],rndm,seed) ;
    if(j!=3){
      totbrch -=_bchweights[j]*xtch[j];    
      totbrcht-=xtch[j];
    }
  }
  _bchweights[0]=Brandomized(totbrch/totbrcht,xlech[0]/totbrcht,xhech[0]/totbrcht,rndm,seed) ;// this is the weight on the D** non resonant only
  cout <<"....... bch and b0 weights 0 " << _bchweights[0]<<" " <<_b0weights[0]<<endl;
  cout <<"....... bch and b0 weights 1 " << _bchweights[1]<<" " <<_b0weights[1]<<endl;
  cout <<"....... bch and b0 weights 2 " << _bchweights[2]<<" " <<_b0weights[2]<<endl;
  cout <<"....... bch and b0 weights 3 " << _bchweights[3]<<" " <<_b0weights[3]<<endl;
  cout <<"....... bch and b0 weights 4 " << _bchweights[4]<<" " <<_b0weights[4]<<endl;
  cout <<"....... bch and b0 weights 5 " << _bchweights[5]<<" " <<_b0weights[5]<<endl;

//    // Write out the weights
//    if(seed){
//      char name[200];
//      sprintf(name,"%s%i%s","bweights",seed,".dat");  
//      ofstream outfile(name);
//      outfile << "b weights from run number " << seed << endl;
//      cout << endl;
//      for(int i=0; i<6; i++)
//        outfile << "ch " << i << " " << _bchweights[i] << endl;
//      for(int i=0; i<6; i++)
//        outfile << "z " << i << "  " << _b0weights[i] << endl;
//      outfile.close();
//    }
}

void recoilDSys::recoilDSys2(int num){
  // Constructor for reading the bweights directly from a file
  // null all B weights
  for (int j=0;j<6;j++){
    _b0weights[j]=0;
    _bchweights[j]=0;
  }
  char  buffer[200];
  sprintf(buffer, "%s%i%s","bweights",num,".dat");
  cout << "Reading " << buffer << endl;
  ifstream is(buffer);
  char chg[5]; int ik; float frac;
  while (is.getline(buffer, 200, '\n')) {
      if (buffer[0] == 'b') {continue;}
      if (buffer[0] == ' ') {continue;}
      sscanf(buffer, "%s %i %f",chg,&ik,&frac);
      if(!strcmp(chg,"ch")) _bchweights[ik] = frac;
      if(!strcmp(chg,"z"))  _b0weights[ik]  = frac;
  }          
  for(int i=0; i<6; i++){
    cout << "b weight charged " << i << "  " << _bchweights[i] << endl;
    cout << "b weight neutral " << i << "  " << _b0weights[i] << endl;
  }
}

void recoilDSys::recoilDSys3(int num){
  // Constructor for reading the dweights directly from a file
  // null all D weights
  for (int j=0;j<MAXNDMODES;j++)
    _weightMode[j] = 0.;
  char  buffer[200];
  sprintf(buffer, "%s%i%s","dweights",num,".dat");
  cout << "Reading " << buffer << endl;
  ifstream is(buffer);
  int ik=0; float frac=1;
  while (is.getline(buffer, 200, '\n')) {
      if (buffer[0] == 'd') {continue;}
      if (buffer[0] == ' ') {continue;}
      sscanf(buffer, "%i %f",&ik,&frac);
      _weightMode[ik] = frac;
  }          
  for(int i=0; i<ik+1; i++){
    cout << "d weight " << i << "  " << _weightMode[i] << endl;
  }
}

float recoilDSys::Brandomized(float mean,float lsig, float hsig, TRandom& rndm, int seed){
  if(seed == 0 ) return mean;// no smearing required
  if(lsig==0. && hsig==0.) return mean; // no smearing possible
  // Version with continuous asymmetric gauss
//   double PI = 3.141592654;
//   TF1 *agaus = new TF1("agausl","sqrt(1/(2*[3]))/[1]*exp(-1*pow(x-[0],2)/2/[1]/[1])*(x<[0]) + sqrt(1/(2*[3]))/[1]*exp(-1*pow(x-[0],2)/2/[2]/[2])*(x>=[0])",-1.,10.);
//   agaus->SetParameters(mean,lsig,hsig,PI);
//   agaus->SetNpx(500);
//   float result = agaus->GetRandom();
//   if(result < 0) result = Brandomized(mean,lsig,hsig,rndm,seed); // if negative weight throw the dice again
//   return result;

  // Version with normal distribution first
   float result;
   float dev=rndm.Gaus(0.,1.);
   if(dev<0)
     result = mean + dev*lsig;
   else
     result = mean + dev*hsig;
   if(result < 0) result = Brandomized(mean,lsig,hsig,rndm,seed); // if negative weight throw the dice again
   return result;
}

float recoilDSys::randomized(float mean,float sig, TRandom& rndm, int seed){
  if(seed == 0 ) return mean;// no smearing required
  float result=rndm.Gaus(mean,sig);
  if(result < 0) result = randomized(mean,sig,rndm,seed); // if negative weight throw the dice again
  return result;
}

void recoilDSys::readFile(const char *name, int seed, int iMode){

  //mode flag needed in order to run the exclusive or the inclusive smearing
  //The code still works for NDMODES > NIMODES 

  char  buffer[200];char  modeName[200];
  sprintf(buffer, "%s", name);
  ifstream is(buffer);
  _nmodes=0;
  TRandom rndm(seed);
//    rndm.SetSeed(111);
//    for(int i=1; i<30*seed+1; i++)
//      rndm.Rndm();

  while (is.getline(buffer, 200, '\n')) {

    float BrMode,PdgMode,PdgErrMode;
    if(iMode ==2) {
      if(_nmodes>=MAXNDMODES){
	cout << " too many D modes in table, more than "<<MAXNDMODES<<endl;
      }
      
      if (buffer[0] == '#') {continue;}
      
      sscanf(buffer, "%s %i %i %i %i %f %f %f", modeName,&_nPiMode[_nmodes],&_nKMode[_nmodes],&_nK0Mode[_nmodes],
	     &_nPi0Mode[_nmodes], &BrMode,&PdgMode,&PdgErrMode);
      
      _weightMode[_nmodes]=randomized(PdgMode/BrMode,PdgErrMode/BrMode,rndm,seed);
      
      cout <<" modes "<<_nmodes << " = " <<_weightMode[_nmodes]<<endl;
      _nmodes++;
    } else {
      if(_nmodes>=MAXNIMODES){
	cout << " too many D modes in table, more than "<<MAXNDMODES<<endl;
      }
      
      if (buffer[0] == '#') {continue;}
    
      sscanf(buffer, "%s %f %f %f", modeName, &BrMode,&PdgMode,&PdgErrMode);
      
      _weightMode[_nmodes]=randomized(PdgMode/BrMode,PdgErrMode/BrMode,rndm,seed);
      
      cout <<" modes "<<_nmodes << " = " <<_weightMode[_nmodes]<<endl;
      _nmodes++;
    }

  }
//    // Write out the weights
//    if(seed){
//      char name[200];
//      sprintf(name,"%s%i%s","dweights",seed,".dat");  
//      ofstream outfile(name);
//      outfile << "d weights from run number " << seed << endl;
//      cout << endl;
//      for(int i=0; i<_nmodes; i++)
//        outfile << i << " " << _weightMode[i] << endl;
//      outfile.close();
//    }
}


float 
recoilDSys::weight(int npi,int nk, int nk0, int npi0, int nlep, int imOde, int & mode){
  if(_nmodes==0) cout <<" not configured for D decays!!!!"<<endl;
  if(imOde == 2) {
    if(nlep !=0)return 1;
    
    for (int j = 0; j<_nmodes; j++)
      if(npi==_nPiMode[j] && nk==_nKMode[j] && nk0==_nK0Mode[j] &&npi0==_nPi0Mode[j] ){
	mode = j;
	return _weightMode[j];
      }
  } else {
    if(((nlep + npi + nk) % 2) == 0) {
      //D0
      if((nk > 0) || (nk0 > 0)) {
	if(nk >0) {
	  mode = 0;
	  return _weightMode[0];
	} else {
	  mode = 1;
	  return _weightMode[1];
	} 
      }
    } else {
      //DC
      if((nk > 0) || (nk0 > 0)) {
	if(nk >0) {
	  mode = 3;
	  return _weightMode[3];
	} else {
	  mode = 2;
	  return _weightMode[2];
	} 
      }
    }
  }
  return 1;
}

float 
recoilDSys::weight(int bmode){
  float result;
  if(bmode==0)return 1;
  int ab=TMath::Abs(bmode);
  if (ab > 10) return 1;
  if(ab == 3 || ab>5) {
    result = bmode > 0 ? _b0weights[0] : _bchweights[0] ;
    if (result<0.00001)cout <<" missing weight for bmode "<< bmode << endl;
    
    return result;
  }
  result = bmode > 0 ? _b0weights[ab] : _bchweights[ab] ;
  if (result<0.00001)cout <<" missing weight for bmode "<< bmode << endl;
  
  return result;
}
