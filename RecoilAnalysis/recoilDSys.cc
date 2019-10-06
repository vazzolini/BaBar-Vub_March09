#include "recoilDSys.hh"

#include <fstream>
#include <iostream>
#include <string>

using namespace std;

ClassImp(recoilDSys)

recoilDSys::recoilDSys(const char *name, int seed, int imode){
  // load input file
  readFile(name,seed,imode);
  // null all B weights
  for (int j=0;j<6;j++){
    _b0weights[j]=0;
    _bchweights[j]=0;

    _b0WeightsAllSemilep=0;
    _bchWeightsAllSemilep=0;
  }
    
}

recoilDSys::recoilDSys(int seed, int rel, int sysbdec): _nmodes(0){

  TRandom rndm(seed);
  //  TRandom rndm(123);
  //  gRandom->SetSeed(123);
  //  for(int i=1; i<10*seed+1; i++)
  //  gRandom->Rndm(); // Effectively sets the seed.
  
  float xv0[6]={10.21,2.07,5.70,0.,0.23,0.52}; // existing measurements (order: all semilep, D,D*,nothing,D2*,D1*)
  float xt0[6]={10.4,2.10,5.60,9.,0.37,0.56};  // values from decay.dec.

  float xle0[6]={0.17,0.15,1.02,0.,0.23,0.15}; // errors on existing measurements downwards
  float xhe0[6]={0.17,0.15,0.53,0.,0.23,0.15}; // errors on existing measurements upwards

  if(rel==18  || rel == 22){
    xv0[0]=10.15; xv0[1]=2.12; xv0[2]=5.10; xv0[4]=0.31; xv0[5]=0.39; // final numbers from Vera Luth 09-03-2009 (order: all semilep, D,D*,nothing,D2*,D1*)
    xt0[0]=10.2;  xt0[1]=2.07; xt0[2]=5.70; xt0[4]=0.23; xt0[5]=0.52;   //SP 8 decay.dec
    xle0[0]=0.38; xle0[1]=0.10; xle0[2]=0.14; xle0[3]=0.; xle0[4]=0.04; xle0[5]=0.03; // errors on existing measurements downwards. for D** (place 0) assumed 100% correlation
    xhe0[0]=0.39; xhe0[1]=0.10; xhe0[2]=0.14; xhe0[3]=0.; xhe0[4]=0.04; xhe0[5]=0.03; // errors on existing measurements upwards. for D** (place 0) assumed 100% correlation
  }

  // fill D*lnu numbers ( in %)
  double totbr0(xv0[0]); //xv0[0]
  double totbr0t(xt0[0]); //xt0[0]
  
  for (int j=1;j<6;j++){
    _b0weights[j]=Brandomized(xv0[j]/xt0[j],xle0[j]/xt0[j],xhe0[j]/xt0[j],rndm,seed) ;
    if(j!=3){
      totbr0 -=_b0weights[j]*xt0[j];
      totbr0t-=xt0[j];
    }
  }

  // We don't want to randomize this component to keep the total B->Xclnu constant
  //_b0weights[0] = Brandomized(totbr0/totbr0t,xle0[0]/totbr0t,xhe0[0]/totbr0t,rndm,seed) ; // this is the weight on the D** broad and non resonant
  
  _b0weights[0] = totbr0/totbr0t;

  float b0allsemilepmeas(0),b0allsemileptrue(0);

  for(Int_t i=1;i<6;i++){
    if(i==3) continue;
    b0allsemilepmeas += _b0weights[i] * xt0[i];
    b0allsemileptrue += xt0[i];
    //cout<<_b0weights[i]<<" "<<xt0[i]<<" "<<xv0[i]<<" "<<_b0weights[i]<<endl;
  }
    
  b0allsemilepmeas +=_b0weights[0] * totbr0t;
  b0allsemileptrue += totbr0t;
  
  //ANTONIO: this is to preserve the total Charmed semileptonic BR
  b0allsemilepmeas = xv0[0];
  
  _b0WeightsAllSemilep =  b0allsemilepmeas / b0allsemileptrue;
  
  //  cout<<"totbr0 "<<totbr0<<" totbr0t "<<totbr0t<<" ratio "<<_b0weights[0]<<endl;
 
  float xvch[6]={11.04,2.24,6.17,0.,0.30,0.56}; // existing measurements DL Pegna  AWG 14May2007 (order: all semilep, D,D*,nothing,D2*,D1*)
  float xtch[6]={10.4,2.10,5.60,9.,0.37,0.56};  // values from decay.dec.

  float xlech[6]={0.18,0.16,1.13,0.,0.3,0.16}; // errors on existing measurements downwards
  float xhech[6]={0.18,0.15,0.83,0.,0.3,0.16}; // errors on existing measurements upwards

  if(rel==18 || rel == 22){
    xvch[0]=10.89; xvch[1]=2.27; xvch[2]=5.47; xvch[3]=0.; xvch[4]=0.32; xvch[5]=0.42; // final numbers: Vera Luth 09-03-2009; order: all semilep, D,D*,nothing,D2*,D1
    xtch[0]=11.04; xtch[1]=2.24; xtch[2]=6.17; xtch[4]=0.30;  xtch[5]=0.56; //SP 8 decay.dec
    xlech[0]=0.37; xlech[1]=0.08; xlech[2]=0.27; xlech[3]=0.; xlech[4]=0.04; xlech[5]=0.03;  // for D** (place 0) assumed 100% correlation
    xhech[0]=0.37; xhech[1]=0.08; xhech[2]=0.27; xhech[3]=0.; xhech[4]=0.04; xhech[5]=0.03;  // for D** (place 0) assumed 100% correlation  
  }
  
  double totbrch(xvch[0]); //xvch[0]
  double totbrcht(xtch[0]); //xtch[0]
  
  for (int j=1;j<6;j++){
    _bchweights[j]=Brandomized(xvch[j]/xtch[j],xlech[j]/xtch[j],xhech[j]/xtch[j],rndm,seed) ;
    if(j!=3){
      totbrch -=_bchweights[j]*xtch[j];    
      totbrcht-=xtch[j];
    }
  }

  // We don't want to randomize this component since we want to keep the total B->Xclnu constant.
  //  _bchweights[0] = Brandomized(totbrch/totbrcht,xlech[0]/totbrcht,xhech[0]/totbrcht,rndm,seed) ;// this is the weight on the D** non resonant only
  _bchweights[0] = totbrch/totbrcht;

  float bchallsemilepmeas(0),bchallsemileptrue(0);

  for(Int_t i=1;i<6;i++){
    if(i==3) continue;
    bchallsemilepmeas += _bchweights[i] * xtch[i];
    bchallsemileptrue += xtch[i];
    //  cout<<_bchweights[i]<<" "<<xtch[i]<<" "<<xvch[i]<<" "<<_bchweights[i]<<endl;
  }
    
  bchallsemilepmeas +=_bchweights[0] * totbrcht;
  bchallsemileptrue += totbrcht;

  // ANTONIO: this is to preserve the total charmed semileptonic BR
  bchallsemilepmeas = xvch[0];
  
  _bchWeightsAllSemilep =  bchallsemilepmeas / bchallsemileptrue;
  
  if(sysbdec != -99) reweightForSys(sysbdec);

  //  cout<<"totbrch "<<totbrch<<" totbrcht "<<totbrcht<<" ratio "<<_bchweights[0]<<endl;

  cout <<"....... bch and b0 weights 0 (B-> other D**) " << _bchweights[0]<<" " <<_b0weights[0]<<endl;
  cout <<"....... bch and b0 weights 1 (B-> D l nu) " << _bchweights[1]<<" " <<_b0weights[1]<<endl;
  cout <<"....... bch and b0 weights 2 (B-> D* l nu) " << _bchweights[2]<<" " <<_b0weights[2]<<endl;
  cout <<"....... bch and b0 weights 3 (not used) " << _bchweights[3]<<" " <<_b0weights[3]<<endl;
  cout <<"....... bch and b0 weights 4 (B-> D*2 l nu) " << _bchweights[4]<<" " <<_b0weights[4]<<endl;
  cout <<"....... bch and b0 weights 5 (B-> D1 l nu) " << _bchweights[5]<<" " <<_b0weights[5]<<endl;
  cout <<"....... bch and b0 weights (B-> Xc l nu) " << _bchWeightsAllSemilep<<" "<<_b0WeightsAllSemilep<<endl;

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

recoilDSys::~recoilDSys(){}

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

  float totBr0Mode(0), totPdg0Mode(0), totBrCMode(0), totPdgCMode(0);
  float *totBrPoint = &totBr0Mode;
  float *totPdgPoint = &totPdg0Mode;
  bool flagDc = false;
  bool notchanged = true;
    
  while (is.getline(buffer, 200, '\n')) {
    
    float BrMode,PdgMode,PdgErrMode;
    
    if( iMode == 2 ) {
      if( _nmodes >= MAXNDMODES )
	cout << " too many D modes in table, more than " << MAXNDMODES << endl;
      
      if ( buffer[0] == '#' ) {
	if( buffer[2] == 'D' && buffer[3] == '+' ) flagDc = true;
	continue;
      }
      
      if( flagDc && notchanged) {
	totBrPoint = &totBrCMode;
	totPdgPoint = &totPdgCMode;
	notchanged = false;
      }

      sscanf( buffer, "%s %i %i %i %i %f %f %f", modeName,&_nPiMode[_nmodes],&_nKMode[_nmodes],&_nK0Mode[_nmodes],
	     &_nPi0Mode[_nmodes], &BrMode,&PdgMode,&PdgErrMode );
      
      _weightMode[_nmodes] = randomized( PdgMode/BrMode,PdgErrMode/BrMode,rndm,seed );
      
      *totBrPoint  += BrMode;
      *totPdgPoint += PdgMode;
      
      cout << " modes " << _nmodes << " = " << _weightMode[_nmodes] << endl;
      _nmodes++;
    } else {
      if( _nmodes >= MAXNIMODES )
	cout << " too many D modes in table, more than " << MAXNDMODES << endl;
     
      if (buffer[0] == '#') continue;
    
      sscanf( buffer, "%s %f %f %f", modeName, &BrMode,&PdgMode,&PdgErrMode );
      
      _weightMode[_nmodes] = randomized( PdgMode/BrMode,PdgErrMode/BrMode,rndm,seed );
      
      cout << " modes " << _nmodes << " = " << _weightMode[_nmodes] << endl;
      _nmodes++;
    }
  }

  _weightMode[_nmodes++] = (1 - totPdg0Mode) / (1 - totBr0Mode); //weight for other D0 channels taking into account unitarity
  _weightMode[_nmodes++] = (1 - totPdgCMode) / (1 - totBrCMode); //weight for other D+ channels taking into account unitarity

//   for(Int_t i = 0; i < _nmodes; i++)
//     cout << _nmodes << " i " << i  << " " <<_weightMode[i] << endl;

  // Fix the weight for normalization

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


float recoilDSys::weight(int npi,int nk, int nk0, int npi0, int nlep, int imOde, int & mode){

  if( _nmodes == 0 ) cout <<" not configured for D decays!!!!"<<endl;
  if( imOde == 2 ) {
    if( nlep !=0 ) {
      mode = ((npi + nk + nlep) % 2) == 0 ? _nmodes-2 : _nmodes-1;
      return _weightMode[mode];
      //      return 1;
    }

    for (int j = 0; j < _nmodes; j++) { 
      if( npi == _nPiMode[j] && nk == _nKMode[j] && nk0 == _nK0Mode[j] && npi0 == _nPi0Mode[j] ){
	mode = j;
	return _weightMode[j];
      }
    } 
  
    mode = ((npi + nk) % 2) == 0 ? _nmodes-2 : _nmodes-1;
    return _weightMode[mode];
    
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
  if(bmode==0) return 1;
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

float recoilDSys::getAllSemilepWeight(bool isBch){

//   float tempbch=0;
//   float tempbneu=0;

 //  for(Int_t i=0;i<6;i++) {
//     tempbch +=_bchweights[i];
//     tempbneu += _b0weights[i];
//   }
  return isBch ? _bchWeightsAllSemilep : _b0WeightsAllSemilep;
}

void recoilDSys::reweightForSys(int bdectype){

  // Use recipe from Vera Luth and David Lopes Pegna  
  
  //
  // convention is the same as in weights array
  // place 0 is for D** broad + NR
  // place 1 is for D lnu
  // place 2 is for D*lnu
  // place 3 is not used
  // place 4,5 is for narrow D**
  //
  // variation array holds sigma/BR_meas values
  //---------------------------------------------------

  
  float variation_bch[6]={0., 0.15/2.30, 0.25/5.95, 0, 0.20, 0.20}; 
  float variation_bneu[6]={0., 0.15/2.13, 0.25/5.53, 0, 0.20, 0.20};

  // ---------------------------------------
  // others array holds the variation to be applied to other weights when 
  // studying the i-th 
  // -------------------------------------
  
  float others[6]={0., 0.018, 0.052, 0., 0.02, 0.02};

  //get sign for variation
  float sign = int(bdectype/10) > 0 ? 1. : -1.;
  
  bdectype = bdectype - 10*int(bdectype/10);

  //reset correct index for D** broad + NR 
  if(bdectype==6) bdectype=0;

  cout<<"Doing B -> Dlv decays systematics! Varying B dec type "<<bdectype<<" with "<< (sign >0 ? "+1" : "-1") <<" sigma variation"<<endl;


  // ------------------------------------------
  // (new weight / old weight) = ( (BR_meas+sigma)/BR_sp8 )/( BR_meas / BR_sp8 )
  //
  // new weight = ( 1+sigma/BR_meas )*old weight
  // ------------------------------------------

  _b0weights[bdectype]  *= 1 + sign* variation_bneu[bdectype];
  _bchweights[bdectype] *= 1 + sign* variation_bch[bdectype];
  
  // Handle narrow D** index: we want to reweight both D** narrow together
  //
  if(bdectype == 4 || bdectype == 5){
    int j = bdectype == 4 ? bdectype + 1 : bdectype - 1;
    _b0weights[j]  *= 1 + sign* variation_bneu[j];
    _bchweights[j] *= 1 + sign* variation_bch[j];
  }

  sign*=-1;
  
  // Reweight other components to have constant B->Xclnu

  for(Int_t i=0; i<6 ;i++){
    if(i==bdectype || i==3) continue; //do not reweight the fraction under evaluation
    if( (bdectype==4 && i==5) || (i==4 && bdectype ==5) ) continue; // do not reweight the other narrow D**
    
    _b0weights[i]  *= 1 + sign * others[bdectype];
    _bchweights[i] *= 1 + sign * others[bdectype];
  }
}
