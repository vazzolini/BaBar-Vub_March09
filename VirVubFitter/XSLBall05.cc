//*****************************************************************
//   
//   Creation: David Cote, Universite de Montreal, 12/11/03
//   see comments in the XSLEvtFFWeight.hh file or BAD #809
//
//   From the LCSR calculation of P. Ball, R. Zwicky, Phys.~Rev.~{\bf D71} 014029 (2005), hep-ph/0412079.
//
//#include "BaBar/BaBar.hh"
 
#include "XSLBall05.hh"
#include "TLorentzVector.h"
#include <iostream>
using std::cout;
using std::endl;
using std::string;

XSLBall05::XSLBall05(double mB,double mXu,double q2,double theta_l,double theta_V,double chi,string mode,string ErrorScale) : 
  XSLVectorFF(mB,mXu,q2,theta_l,theta_V,chi) 
{ Init(mode,ErrorScale); }

XSLBall05::XSLBall05( XSLKin* DecayKin,string mode,string ErrorScale ) : XSLVectorFF(DecayKin) 
{ Init(mode,ErrorScale); }

XSLBall05::XSLBall05(TLorentzVector BLab,TLorentzVector LepLab,TLorentzVector XuLab,TLorentzVector XuDaughterLab,
		     string mode,string ErrorScale ) :  XSLVectorFF(BLab, LepLab, XuLab, XuDaughterLab) 
{ Init(mode,ErrorScale); }

void
XSLBall05::Init(string mode, string ErrorScale)
{ 
  //default values: do nothing
  _A1ScaleFunction=&XSLBall05::ReturnOne;
  _A2ScaleFunction=&XSLBall05::ReturnOne;
  _VScaleFunction=&XSLBall05::ReturnOne;

  //Scale one of the FF by some error function...
  if(ErrorScale=="A1Pos"){ _A1ScaleFunction=&XSLBall05::ErrorScaleFactorPos; } 
  else if(ErrorScale=="A1Neg"){ _A1ScaleFunction=&XSLBall05::ErrorScaleFactorNeg; } 
  else if(ErrorScale=="A2Pos"){ _A2ScaleFunction=&XSLBall05::ErrorScaleFactorPos; } 
  else if(ErrorScale=="A2Neg"){ _A2ScaleFunction=&XSLBall05::ErrorScaleFactorNeg; } 
  else if(ErrorScale=="VPos"){ _VScaleFunction=&XSLBall05::ErrorScaleFactorPos; } 
  else if(ErrorScale=="VNeg"){ _VScaleFunction=&XSLBall05::ErrorScaleFactorNeg; } 
  else if(ErrorScale.rfind("Rand",5)==0){ //This means: "if ErrorScale begins with Rand"
    string KeyStr = ErrorScale.substr(4); //This returns all characters from pos 4 to the end.
    const char* KeyChar= KeyStr.c_str();
    int KeyInt=atoi(KeyChar); 
    if(KeyInt<0 || KeyInt>99){cout<<"XslFFReweighting/XSLBall05: Invalid KeyInt "<<KeyInt<<endl; exit(1); }
    InitLocalGaussianRandomNumberBank();
    _sfA1=double(_rnd[3*KeyInt]);
    _sfA2=double(_rnd[3*KeyInt+1]);
    _sfV=double(_rnd[3*KeyInt+2]);
    _A1ScaleFunction=&XSLBall05::ErrorScaleFactorRandA1; 
    _A2ScaleFunction=&XSLBall05::ErrorScaleFactorRandA2; 
    _VScaleFunction=&XSLBall05::ErrorScaleFactorRandV; 
  }
  else if(ErrorScale!="NoError"){ cout << "Unknown error function in XSLBall05!!  ErrorScale = " << ErrorScale << endl; exit(1); }

  _ErrorScale=ErrorScale; //ErrorScale doesn't need to be a class variable, but SetNormalizations needs its value 
                          //and is required to have only one argument, so this is an easy hack.
  SetNormalizations(mode);  
  SetConstants(mode); 
  Compute(); 
}


void
XSLBall05::SetConstants(string mode){

  //By default, the constants are set to B->rholnu, but they also can be set to B->omegalnu
 
  // hep-ph/0412079
  // from Table 8 
  if(mode=="omegalnu") {
    _r2_A1=0.217;
    _mfit2_A1=37.01;
    _r1_A2=0.006;
    _r2_A2=0.192;
    _mfit2_A2=41.24;
    _r1_V=1.006;
    _r2_V=-0.713;
    _mfit2_V=37.45;
  }
  else { //rholnu
    _r2_A1=0.240;
    _mfit2_A1=37.51;
    _r1_A2=0.009;
    _r2_A2=0.212;
    _mfit2_A2=40.82;
    _r1_V=1.045;
    _r2_V=-0.721;
    _mfit2_V=38.34;
  }
  return;
}

void 
XSLBall05::GetAllFF(double *A1, double *A2, double *V)
{
  *A1=GetA1(_q2);
  *A2=GetA2(_q2);
  *V=GetV(_q2);
}


double 
XSLBall05::GetA1(double q2)
{
 // from Eqn (61)
  double a1 = _r2_A1/(1.-q2/_mfit2_A1);
  double f=(this->*_A1ScaleFunction)(q2);
  return f*a1;
}

double 
XSLBall05::GetA2(double q2)
{
  // from Eqn (60)
  double a2 = _r1_A2/(1.-q2/_mfit2_A2) + _r2_A2/pow(1.-q2/_mfit2_A2,2.);
  double f=(this->*_A2ScaleFunction)(q2);
  return f*a2;

}

double 
XSLBall05::GetV(double q2)
{ 
  const double m1 = 5.32; // B* mass
  // from Eqn (59)
  double v = _r1_V/(1.-q2/m1/m1) + _r2_V/(1.-q2/_mfit2_V);
  double f=(this->*_VScaleFunction)(q2);
  return f*v;
}

void
XSLBall05::SetNormalizations(string mode)
{
  //Note: These normalization constants were determined empirically from MC samples of 1000k events.
  //The quoted error is simply the statistical one, but there can be a bigger error coming from rare high weight events for
  //some region of the phase space, mainly for the Vector-ISGW2 generator. 
  
  _FLATQ2Normalization=1.0; _PHSPNormalization=1.0; _ISGW2Normalization=1.0; 
  
  //WARNING: Weights won't be properly normalized if you're using the PHSP generator with ErrorScale!=NoError...
  //FLATQ2 or ISGW2 Generators are fine!    :-)

  if(_ErrorScale.rfind("Rand",5)==0){ //This means: "if ErrorScale begins with Rand"
    if(mode!="rhoClnu"&&mode!="rho0lnu"&&mode!="omegalnu"&&mode!="omegalnu_o"){ cout<<"Error in XSLBall05 "<<_ErrorScale<<" "<<mode<<endl; exit(1); }
    string KeyStr = _ErrorScale.substr(4); //This returns all characters from pos 4 to the end.
    const char* KeyChar= KeyStr.c_str();
    int KeyInt=atoi(KeyChar); 
    _ISGW2Normalization=RandISGW2Norm(KeyInt);
    _FLATQ2Normalization=RandFLATQ2Norm(KeyInt);
  }
  else if(mode=="rhoClnu") {
    if(_ErrorScale=="NoError"){
      _FLATQ2Normalization=0.09232327364;//determined from 100k events: 0.3% error stat
      _PHSPNormalization=0.1393226273;   //determined from 100k events: 0.3% error stat
      _ISGW2Normalization=0.8447503878;  //determined from 100k events: 0.3% error stat 
                                      //but the effect of rare high weight event seem to introduce much larger uncertainties from ISGW2...
    }
    else if(_ErrorScale=="A1Pos"){
      _FLATQ2Normalization=1.0/14.10562000;
      _ISGW2Normalization=1.0/1.5742245;
    }
    else if(_ErrorScale=="A1Neg"){
      _FLATQ2Normalization=1.0/7.443169500;
      _ISGW2Normalization=1.0/0.8317134375;
    }
    else if(_ErrorScale=="A2Pos"){
      _FLATQ2Normalization=1.0/9.222988000;
      _ISGW2Normalization=1.0/1.03149775;
    }
    else if(_ErrorScale=="A2Neg"){
      _FLATQ2Normalization=1.0/11.777671;
      _ISGW2Normalization=1.0/1.313077375;
    }
    else if(_ErrorScale=="VPos"){
      _FLATQ2Normalization=1.0/10.935409;
      _ISGW2Normalization=1.0/1.219753;
    }
    else if(_ErrorScale=="VNeg"){
      _FLATQ2Normalization=1.0/9.928505;
      _ISGW2Normalization=1.0/1.110209375;
    }
    else{ cout<<"XSLBall05: ErrorScale: "<<_ErrorScale<<" not implemented for mode "<<mode<<endl; exit(1); }
  }
  else if(mode=="rho0lnu") {
    if(_ErrorScale=="NoError"){
      _FLATQ2Normalization=0.09248629758;//determined from 100k events: 0.3% error stat
      _PHSPNormalization=0.138725886;    //determined from 100k events: 0.3% error stat
      _ISGW2Normalization=0.8443531312;  //determined from 100k events: 0.3% error stat
                                         //but the effect of rare high weight event seem to introduce much larger uncertainties from ISGW2...
    }
    else if(_ErrorScale=="A1Pos"){
      _FLATQ2Normalization=1.0/14.08362064;
      _ISGW2Normalization=1.0/1.577346125;
    }
    else if(_ErrorScale=="A1Neg"){
      _FLATQ2Normalization=1.0/7.423538135;
      _ISGW2Normalization=1.0/0.8315596636;
    }
    else if(_ErrorScale=="A2Pos"){
      _FLATQ2Normalization=1.0/9.200457279;
      _ISGW2Normalization=1.0/1.031787894;
    }
    else if(_ErrorScale=="A2Neg"){
      _FLATQ2Normalization=1.0/11.75889457;
      _ISGW2Normalization=1.0/1.315351717;
    }
    else if(_ErrorScale=="VPos"){
      _FLATQ2Normalization=1.0/10.91164425;
      _ISGW2Normalization=1.0/1.220436635;
    }
    else if(_ErrorScale=="VNeg"){
      _FLATQ2Normalization=1.0/9.910329046;
      _ISGW2Normalization=1.0/1.111868493;
    }
    else{ cout<<"XSLBall05: ErrorScale: "<<_ErrorScale<<" not implemented for mode "<<mode<<endl; exit(1); }
  }
  else if(mode=="omegalnu" || mode=="omegalnu_o") {
    if(_ErrorScale=="NoError"){ 
      _FLATQ2Normalization=0.111093; //determined from 500k events
      _PHSPNormalization=1.0;        //wrong value! will need to be determined...
      _ISGW2Normalization=1.0;       //wrong value! will need to be determined...
    }
    else{ cout<<"XSLBall05: ErrorScale: "<<_ErrorScale<<" not implemented for mode "<<mode<<endl; exit(1); }
  }
  else if(mode=="pilnu"||mode=="pi0lnu"||mode=="etalnu"||mode=="eta2lnu"||mode=="eta3lnu"||
	  mode=="etaplnu"||mode=="etaplnuE2PP"||mode=="etaplnuE3PP"||mode=="etaplnuRG") { 
    cout<<"XSLBall05::SetNormalizations  -- This model is not intended to be used with mode: "<<mode<<" !!!!"<<endl; 
    exit(1);
  }
  else if(mode!="Ignore") {
    cout<<"XSLBall05::SetNormalizations  -- Unknown mode: "<<mode<<" !!!"<<endl; 
    exit(1);
  }
  
  return;
}


double 
XSLBall05::ErrorMethodA(double q2)
{
  //The error is 10% @ q2==0, 11.5% @ q2==7, 13% @ q2==14, and the extrapolation continues linearly for higher q2.
  //see: http://babar-hn.slac.stanford.edu:5090/HyperNews/get/semi_lept_decays/380.html
  float err=0.1+0.03*q2/14.0; 
  return err;
}

double
XSLBall05::RandISGW2Norm(int key){
  //code to obtain norm[]: XslFFReweighting/macros/GetNormalizationFactorsRholnu.cxx
  double norm[100];
  norm[0]=1.1722;  norm[1]=0.7402; norm[2]=0.8566; norm[3]=0.4162; norm[4]=0.6569; norm[5]=1.1667; norm[6]=0.5217;
  norm[7]=0.7960;  norm[8]=0.7137; norm[9]=0.7394; norm[10]=0.8743; norm[11]=0.5951; norm[12]=0.6029;
  norm[13]=2.2835; norm[14]=1.7679; norm[15]=0.6999; norm[16]=1.3810; norm[17]=1.0882; norm[18]=1.5713;
  norm[19]=0.6867; norm[20]=0.6248; norm[21]=0.6306; norm[22]=0.9934; norm[23]=0.8365; norm[24]=0.7341;
  norm[25]=1.1055; norm[26]=0.7603; norm[27]=0.7592; norm[28]=0.5813; norm[29]=0.7415; norm[30]=1.1361;
  norm[31]=0.9875; norm[32]=0.7441; norm[33]=0.4540; norm[34]=1.4729; norm[35]=1.0266; norm[36]=0.6954;
  norm[37]=0.9955; norm[38]=0.7800; norm[39]=0.6043; norm[40]=0.4862; norm[41]=1.1597; norm[42]=0.7692;
  norm[43]=0.6639; norm[44]=0.7658; norm[45]=1.1988; norm[46]=0.9129; norm[47]=0.7664; norm[48]=1.4905;
  norm[49]=0.7983; norm[50]=1.1441; norm[51]=0.6035; norm[52]=0.7844; norm[53]=0.7457; norm[54]=1.1584;
  norm[55]=0.8056; norm[56]=0.7955; norm[57]=0.8247; norm[58]=1.2805; norm[59]=0.7747; norm[60]=0.4614;
  norm[61]=1.0331; norm[62]=0.9082; norm[63]=0.5799; norm[64]=0.7368; norm[65]=0.7575; norm[66]=0.5719;
  norm[67]=0.7332; norm[68]=1.0655; norm[69]=0.5575; norm[70]=0.9172; norm[71]=0.7201; norm[72]=0.6211;
  norm[73]=0.7508; norm[74]=0.6220; norm[75]=0.7642; norm[76]=0.4940; norm[77]=0.4783; norm[78]=0.7622;
  norm[79]=1.2974; norm[80]=0.7029; norm[81]=1.5619; norm[82]=0.8976; norm[83]=0.6985; norm[84]=1.1374;
  norm[85]=0.9283; norm[86]=1.2607; norm[87]=0.4480; norm[88]=1.5053; norm[89]=1.2643; norm[90]=0.6561;
  norm[91]=0.6151; norm[92]=0.6486; norm[93]=0.8027; norm[94]=0.7552; norm[95]=1.0590; norm[96]=0.4893;
  norm[97]=1.0628; norm[98]=1.0517; norm[99]=0.7131;
  return norm[key];
}

double
XSLBall05::RandFLATQ2Norm(int key){
  //code to obtain norm[]: XslFFReweighting/macros/GetNormalizationFactorsRholnu.cxx
  double norm[100];
  norm[0]=0.1311;  norm[1]=0.0824; norm[2]=0.0956; norm[3]=0.0467; norm[4]=0.0736; norm[5]=0.1306; norm[6]=0.0584;
  norm[7]=0.0892;  norm[8]=0.0799; norm[9]=0.0824; norm[10]=0.0975; norm[11]=0.0665; norm[12]=0.0674;
  norm[13]=0.2555; norm[14]=0.1977; norm[15]=0.0784; norm[16]=0.1542; norm[17]=0.1215; norm[18]=0.1752;
  norm[19]=0.0765; norm[20]=0.0699; norm[21]=0.0704; norm[22]=0.1111; norm[23]=0.0937; norm[24]=0.0821;
  norm[25]=0.1236; norm[26]=0.0851; norm[27]=0.0849; norm[28]=0.0648; norm[29]=0.0830; norm[30]=0.1269;
  norm[31]=0.1101; norm[32]=0.0830; norm[33]=0.0506; norm[34]=0.1644; norm[35]=0.1152; norm[36]=0.0779;
  norm[37]=0.1112; norm[38]=0.0873; norm[39]=0.0674; norm[40]=0.0544; norm[41]=0.1297; norm[42]=0.0862;
  norm[43]=0.0744; norm[44]=0.0861; norm[45]=0.1341; norm[46]=0.1021; norm[47]=0.0857; norm[48]=0.1664;
  norm[49]=0.0896; norm[50]=0.1283; norm[51]=0.0677; norm[52]=0.0879; norm[53]=0.0836; norm[54]=0.1292;
  norm[55]=0.0901; norm[56]=0.0891; norm[57]=0.0923; norm[58]=0.1429; norm[59]=0.0870; norm[60]=0.0515;
  norm[61]=0.1152; norm[62]=0.1013; norm[63]=0.0648; norm[64]=0.0823; norm[65]=0.0846; norm[66]=0.0642;
  norm[67]=0.0820; norm[68]=0.1191; norm[69]=0.0623; norm[70]=0.1026; norm[71]=0.0804; norm[72]=0.0695;
  norm[73]=0.0838; norm[74]=0.0696; norm[75]=0.0856; norm[76]=0.0551; norm[77]=0.0534; norm[78]=0.0853;
  norm[79]=0.1448; norm[80]=0.0784; norm[81]=0.1743; norm[82]=0.1001; norm[83]=0.0780; norm[84]=0.1271;
  norm[85]=0.1039; norm[86]=0.1407; norm[87]=0.0501; norm[88]=0.1681; norm[89]=0.1412; norm[90]=0.0732;
  norm[91]=0.0688; norm[92]=0.0727; norm[93]=0.0896; norm[94]=0.0843; norm[95]=0.1184; norm[96]=0.0548;
  norm[97]=0.1190; norm[98]=0.1175; norm[99]=0.0797;
  return norm[key];
}


void
XSLBall05::InitLocalGaussianRandomNumberBank(){
  //Pre-generated random numbers according to a Gaussian with mean=0.0 and sigma=1.0.
  _rnd[0]=-0.8766; _rnd[1]=-0.0889; _rnd[2]=-0.4784; _rnd[3]=-0.3410; _rnd[4]=-1.4877; _rnd[5]=1.5632; _rnd[6]=0.0534;
  _rnd[7]=0.7292; _rnd[8]=1.3930; _rnd[9]=2.4427; _rnd[10]=-0.6060; _rnd[11]=-1.5728; _rnd[12]=1.5393;
  _rnd[13]=1.7025; _rnd[14]=-0.1173; _rnd[15]=-0.5785; _rnd[16]=0.6296; _rnd[17]=-0.6463; _rnd[18]=1.6969;
  _rnd[19]=0.0115; _rnd[20]=-0.0118; _rnd[21]=0.7762; _rnd[22]=1.2699; _rnd[23]=-0.4365; _rnd[24]=0.9883;
  _rnd[25]=0.9897; _rnd[26]=-0.1149; _rnd[27]=-0.2632; _rnd[28]=-1.7045; _rnd[29]=0.4871; _rnd[30]=-0.3769;
  _rnd[31]=-0.3406; _rnd[32]=1.1936; _rnd[33]=1.5788; _rnd[34]=1.5038; _rnd[35]=1.7055; _rnd[36]=1.1976;
  _rnd[37]=0.3098; _rnd[38]=0.7802; _rnd[39]=-2.3950; _rnd[40]=-0.0098; _rnd[41]=-1.4390; _rnd[42]=-1.4458;
  _rnd[43]=1.8104; _rnd[44]=-0.4762; _rnd[45]=1.0106; _rnd[46]=0.7681; _rnd[47]=-0.4582; _rnd[48]=-1.1637;
  _rnd[49]=1.3989; _rnd[50]=0.8994; _rnd[51]=-0.9582; _rnd[52]=-0.7671; _rnd[53]=-0.2357; _rnd[54]=-2.2961;
  _rnd[55]=-1.3680; _rnd[56]=0.2611; _rnd[57]=-0.1556; _rnd[58]=-2.0678; _rnd[59]=0.3818; _rnd[60]=1.5554;
  _rnd[61]=1.5776; _rnd[62]=0.7897; _rnd[63]=1.0019; _rnd[64]=0.2584; _rnd[65]=0.8990; _rnd[66]=-0.5202;
  _rnd[67]=-0.5086; _rnd[68]=-0.7622; _rnd[69]=0.4002; _rnd[70]=0.6153; _rnd[71]=-0.6254; _rnd[72]=0.1261;
  _rnd[73]=-1.1731; _rnd[74]=-0.7261; _rnd[75]=-0.5443; _rnd[76]=0.6159; _rnd[77]=0.0355; _rnd[78]=0.2026;
  _rnd[79]=-0.7630; _rnd[80]=-0.9337; _rnd[81]=0.1405; _rnd[82]=-0.7751; _rnd[83]=-0.4037; _rnd[84]=0.9132;
  _rnd[85]=-0.5666; _rnd[86]=1.3065; _rnd[87]=0.7122; _rnd[88]=0.5060; _rnd[89]=-0.3868; _rnd[90]=-0.7529;
  _rnd[91]=0.3926; _rnd[92]=0.2950; _rnd[93]=-0.6148; _rnd[94]=0.4391; _rnd[95]=1.8337; _rnd[96]=0.2545;
  _rnd[97]=-0.0558; _rnd[98]=1.2586; _rnd[99]=1.5735; _rnd[100]=-1.2847; _rnd[101]=1.1892; _rnd[102]=-1.8503;
  _rnd[103]=-0.8854; _rnd[104]=-0.1628; _rnd[105]=-0.3569; _rnd[106]=-0.1930; _rnd[107]=-1.8340; _rnd[108]=0.8509;
  _rnd[109]=0.3058; _rnd[110]=-0.4269; _rnd[111]=-0.4398; _rnd[112]=0.3092; _rnd[113]=0.6063; _rnd[114]=0.3746;
  _rnd[115]=-0.0528; _rnd[116]=-0.7102; _rnd[117]=0.9051; _rnd[118]=-0.2800; _rnd[119]=1.1738; _rnd[120]=1.6646;
  _rnd[121]=-0.6815; _rnd[122]=0.1558; _rnd[123]=-0.8883; _rnd[124]=-0.2160; _rnd[125]=-0.5059; _rnd[126]=0.9111;
  _rnd[127]=1.4178; _rnd[128]=-0.1999; _rnd[129]=0.9021; _rnd[130]=-0.0934; _rnd[131]=-0.9219; _rnd[132]=0.7438;
  _rnd[133]=0.4483; _rnd[134]=-1.8227; _rnd[135]=-1.0230; _rnd[136]=-0.5517; _rnd[137]=-1.0102; _rnd[138]=-0.3696;
  _rnd[139]=-0.7109; _rnd[140]=-0.6301; _rnd[141]=0.0516; _rnd[142]=-0.8965; _rnd[143]=-0.3165; _rnd[144]=-1.5857;
  _rnd[145]=0.5497; _rnd[146]=0.6079; _rnd[147]=0.6492; _rnd[148]=0.6561; _rnd[149]=-1.3814; _rnd[150]=-0.5585;
  _rnd[151]=0.2233; _rnd[152]=-1.3767; _rnd[153]=1.5509; _rnd[154]=0.8318; _rnd[155]=-0.5478; _rnd[156]=0.8800;
  _rnd[157]=1.5035; _rnd[158]=-0.2109; _rnd[159]=1.2490; _rnd[160]=2.1356; _rnd[161]=0.0035; _rnd[162]=-1.3651;
  _rnd[163]=-1.0744; _rnd[164]=0.3947; _rnd[165]=0.3149; _rnd[166]=0.4166; _rnd[167]=0.3200; _rnd[168]=0.3874;
  _rnd[169]=0.1536; _rnd[170]=-0.6694; _rnd[171]=0.1818; _rnd[172]=-0.0582; _rnd[173]=-0.5725; _rnd[174]=-1.2494;
  _rnd[175]=0.6066; _rnd[176]=1.1968; _rnd[177]=1.0620; _rnd[178]=1.6643; _rnd[179]=-0.8702; _rnd[180]=1.3024;
  _rnd[181]=-2.0852; _rnd[182]=-0.0689; _rnd[183]=-1.1936; _rnd[184]=-1.6896; _rnd[185]=-0.0666; _rnd[186]=-0.0987;
  _rnd[187]=1.6288; _rnd[188]=2.6450; _rnd[189]=0.7754; _rnd[190]=-1.3733; _rnd[191]=-0.2194; _rnd[192]=0.4635;
  _rnd[193]=0.1021; _rnd[194]=0.4886; _rnd[195]=0.1805; _rnd[196]=-0.4029; _rnd[197]=0.4594; _rnd[198]=1.5834;
  _rnd[199]=0.2308; _rnd[200]=-1.4014; _rnd[201]=0.4948; _rnd[202]=-0.1693; _rnd[203]=-0.4423; _rnd[204]=-0.3931;
  _rnd[205]=0.7901; _rnd[206]=0.1512; _rnd[207]=1.0429; _rnd[208]=-0.9784; _rnd[209]=0.0961; _rnd[210]=0.1602;
  _rnd[211]=1.1444; _rnd[212]=0.3677; _rnd[213]=0.4116; _rnd[214]=-0.2717; _rnd[215]=0.3613; _rnd[216]=0.9374;
  _rnd[217]=-0.3316; _rnd[218]=-0.0312; _rnd[219]=0.0675; _rnd[220]=-0.8692; _rnd[221]=0.1550; _rnd[222]=0.9548;
  _rnd[223]=-0.3256; _rnd[224]=-0.2098; _rnd[225]=0.3590; _rnd[226]=-0.3206; _rnd[227]=-0.9207; _rnd[228]=1.9441;
  _rnd[229]=0.6376; _rnd[230]=1.8238; _rnd[231]=1.8186; _rnd[232]=-0.2664; _rnd[233]=0.9014; _rnd[234]=0.6086;
  _rnd[235]=0.6210; _rnd[236]=0.0600; _rnd[237]=-1.2028; _rnd[238]=1.0722; _rnd[239]=1.4165; _rnd[240]=0.0628;
  _rnd[241]=-1.2827; _rnd[242]=0.5670; _rnd[243]=-2.0726; _rnd[244]=-1.1358; _rnd[245]=-0.2714; _rnd[246]=-0.4136;
  _rnd[247]=-0.3851; _rnd[248]=0.8253; _rnd[249]=1.0592; _rnd[250]=1.9285; _rnd[251]=2.2076; _rnd[252]=-0.8218;
  _rnd[253]=0.0345; _rnd[254]=-0.0421; _rnd[255]=0.0407; _rnd[256]=0.6169; _rnd[257]=-0.3353; _rnd[258]=-1.4929;
  _rnd[259]=-1.0728; _rnd[260]=-0.1832; _rnd[261]=1.6447; _rnd[262]=-1.5327; _rnd[263]=-0.0247; _rnd[264]=-1.6517;
  _rnd[265]=-0.1443; _rnd[266]=-0.1586; _rnd[267]=-1.1189; _rnd[268]=0.0634; _rnd[269]=-0.0243; _rnd[270]=0.4632;
  _rnd[271]=-0.7017; _rnd[272]=1.0515; _rnd[273]=1.2110; _rnd[274]=0.3934; _rnd[275]=0.3631; _rnd[276]=0.9304;
  _rnd[277]=-0.2865; _rnd[278]=-1.1988; _rnd[279]=-0.1992; _rnd[280]=-0.8670; _rnd[281]=0.5327; _rnd[282]=0.4463;
  _rnd[283]=0.3951; _rnd[284]=0.8119; _rnd[285]=-0.5433; _rnd[286]=0.2233; _rnd[287]=-0.0484; _rnd[288]=1.0256;
  _rnd[289]=-2.4949; _rnd[290]=-1.4643; _rnd[291]=-0.0022; _rnd[292]=1.7251; _rnd[293]=-0.1906; _rnd[294]=0.0031;
  _rnd[295]=2.3231; _rnd[296]=0.9741; _rnd[297]=0.7345; _rnd[298]=0.4403; _rnd[299]=0.2537; _rnd[300]=-0.2449;
  _rnd[301]=-0.2134; _rnd[302]=1.5871; _rnd[303]=0.5122; _rnd[304]=0.2070; _rnd[305]=-0.4720; _rnd[306]=0.5002;
  _rnd[307]=0.0848; _rnd[308]=-0.8666; _rnd[309]=1.6471; _rnd[310]=-0.3969; _rnd[311]=-1.0932; _rnd[312]=0.2794;
  _rnd[313]=-1.4782; _rnd[314]=0.9828; _rnd[315]=0.1977; _rnd[316]=0.0260; _rnd[317]=0.2810; _rnd[318]=1.1698;
  _rnd[319]=-0.6583; _rnd[320]=-0.2043; _rnd[321]=1.3122; _rnd[322]=-0.7695; _rnd[323]=-0.8263; _rnd[324]=0.6468;
  _rnd[325]=-0.9753; _rnd[326]=-0.6040; _rnd[327]=-0.5490; _rnd[328]=0.4976; _rnd[329]=0.7711; _rnd[330]=-0.7659;
  _rnd[331]=-0.2915; _rnd[332]=-0.6804; _rnd[333]=1.2337; _rnd[334]=1.0853; _rnd[335]=0.1058; _rnd[336]=-0.2952;
  _rnd[337]=-1.0236; _rnd[338]=-1.9620; _rnd[339]=1.9197; _rnd[340]=-1.3045; _rnd[341]=0.0488; _rnd[342]=1.1171;
  _rnd[343]=0.1240; _rnd[344]=2.2740; _rnd[345]=0.4135; _rnd[346]=-0.1130; _rnd[347]=-0.2872; _rnd[348]=1.0334;
  _rnd[349]=-0.2200; _rnd[350]=-1.1030; _rnd[351]=-0.4407; _rnd[352]=-0.5210; _rnd[353]=-1.5840; _rnd[354]=-1.2382;
  _rnd[355]=0.4596; _rnd[356]=-0.3103; _rnd[357]=-0.2625; _rnd[358]=0.2009; _rnd[359]=0.0230; _rnd[360]=-0.5180;
  _rnd[361]=-0.1774; _rnd[362]=0.2737; _rnd[363]=-0.6573; _rnd[364]=0.0830; _rnd[365]=-0.6585; _rnd[366]=0.3559;
  _rnd[367]=1.3479; _rnd[368]=-0.9662; _rnd[369]=-0.7459; _rnd[370]=-0.3380; _rnd[371]=-0.0180; _rnd[372]=0.2509;
  _rnd[373]=0.3647; _rnd[374]=-0.3878; _rnd[375]=0.8171; _rnd[376]=-1.6338; _rnd[377]=0.4193; _rnd[378]=-2.0254;
  _rnd[379]=-0.5669; _rnd[380]=0.0938; _rnd[381]=-0.3300; _rnd[382]=-0.2998; _rnd[383]=1.2319; _rnd[384]=-0.7417;
  _rnd[385]=0.6375; _rnd[386]=-1.8402; _rnd[387]=-1.4920; _rnd[388]=-0.1103; _rnd[389]=-0.1049; _rnd[390]=-1.7383;
  _rnd[391]=-1.2308; _rnd[392]=1.2244; _rnd[393]=-0.4083; _rnd[394]=-0.4276; _rnd[395]=-0.4532; _rnd[396]=0.7934;
  _rnd[397]=-2.6660; _rnd[398]=0.6281; _rnd[399]=-0.4752; _rnd[400]=1.2439; _rnd[401]=0.4645; _rnd[402]=-0.7066;
  _rnd[403]=-1.6760; _rnd[404]=-1.5745; _rnd[405]=0.1898; _rnd[406]=-0.0233; _rnd[407]=1.4248; _rnd[408]=-0.8844;
  _rnd[409]=0.4519; _rnd[410]=-1.4373; _rnd[411]=0.3028; _rnd[412]=2.2053; _rnd[413]=-1.2866; _rnd[414]=0.9364;
  _rnd[415]=-0.3570; _rnd[416]=-0.0028; _rnd[417]=-0.2860; _rnd[418]=1.6341; _rnd[419]=1.7346; _rnd[420]=1.3204;
  _rnd[421]=0.5739; _rnd[422]=0.0142; _rnd[423]=-1.5052; _rnd[424]=-1.0151; _rnd[425]=0.4362; _rnd[426]=-0.1947;
  _rnd[427]=0.0272; _rnd[428]=0.7934; _rnd[429]=0.0554; _rnd[430]=-1.1192; _rnd[431]=0.2690; _rnd[432]=-0.9673;
  _rnd[433]=-0.2973; _rnd[434]=-0.7961; _rnd[435]=1.4551; _rnd[436]=-0.4654; _rnd[437]=-0.9932; _rnd[438]=-1.5481;
  _rnd[439]=-0.3311; _rnd[440]=0.0946; _rnd[441]=-1.8256; _rnd[442]=0.0610; _rnd[443]=0.6799; _rnd[444]=-1.8908;
  _rnd[445]=0.9026; _rnd[446]=0.7203; _rnd[447]=-0.6906; _rnd[448]=-0.9062; _rnd[449]=0.9055; _rnd[450]=0.0122;
  _rnd[451]=-2.0666; _rnd[452]=-0.0302; _rnd[453]=0.3732; _rnd[454]=1.4499; _rnd[455]=-1.0339; _rnd[456]=0.2045;
  _rnd[457]=0.0426; _rnd[458]=0.0408; _rnd[459]=0.4840; _rnd[460]=0.3074; _rnd[461]=0.0067; _rnd[462]=0.8228;
  _rnd[463]=-0.7768; _rnd[464]=1.0403; _rnd[465]=0.7617; _rnd[466]=0.5999; _rnd[467]=-0.2402; _rnd[468]=1.8528;
  _rnd[469]=0.7603; _rnd[470]=0.7524; _rnd[471]=0.0517; _rnd[472]=2.7086; _rnd[473]=0.5925; _rnd[474]=0.2453;
  _rnd[475]=0.8861; _rnd[476]=-0.6090; _rnd[477]=-1.2157; _rnd[478]=0.8045; _rnd[479]=-0.2026; _rnd[480]=0.4605;
  _rnd[481]=1.0665; _rnd[482]=2.0682; _rnd[483]=-0.4110; _rnd[484]=0.6199; _rnd[485]=1.9008; _rnd[486]=0.7841;
  _rnd[487]=-0.3811; _rnd[488]=0.0430; _rnd[489]=0.7354; _rnd[490]=-0.2017; _rnd[491]=-1.0250; _rnd[492]=0.5512;
  _rnd[493]=-0.1928; _rnd[494]=0.3036; _rnd[495]=1.0383; _rnd[496]=0.5620; _rnd[497]=-0.6949; _rnd[498]=0.8846;
  _rnd[499]=-0.2742; _rnd[500]=-0.2517; _rnd[501]=-0.8065; _rnd[502]=0.7048; _rnd[503]=0.7590; _rnd[504]=-0.7695;
  _rnd[505]=1.0952; _rnd[506]=1.5478; _rnd[507]=0.3113; _rnd[508]=-0.9879; _rnd[509]=1.5226; _rnd[510]=-0.9905;
  _rnd[511]=-0.1147; _rnd[512]=0.5830; _rnd[513]=0.2602; _rnd[514]=0.7030; _rnd[515]=0.4502; _rnd[516]=0.5244;
  _rnd[517]=0.3131; _rnd[518]=-0.2141; _rnd[519]=-0.3155; _rnd[520]=-2.2485; _rnd[521]=1.1038; _rnd[522]=-0.0025;
  _rnd[523]=-0.3657; _rnd[524]=-0.6867; _rnd[525]=0.5622; _rnd[526]=-0.0343; _rnd[527]=-0.3727; _rnd[528]=0.7694;
  _rnd[529]=-0.5084; _rnd[530]=-3.5204; _rnd[531]=-0.5303; _rnd[532]=0.0502; _rnd[533]=-0.8866; _rnd[534]=-1.5060;
  _rnd[535]=1.3474; _rnd[536]=1.6158; _rnd[537]=-0.1883; _rnd[538]=2.4761; _rnd[539]=-0.6325; _rnd[540]=-0.2165;
  _rnd[541]=0.4222; _rnd[542]=-0.0978; _rnd[543]=0.4319; _rnd[544]=-2.0583; _rnd[545]=-1.7938; _rnd[546]=-1.1924;
  _rnd[547]=1.6228; _rnd[548]=1.3944; _rnd[549]=-1.4302; _rnd[550]=-0.0526; _rnd[551]=-2.4387; _rnd[552]=-0.1687;
  _rnd[553]=1.2948; _rnd[554]=-0.3640; _rnd[555]=-0.2156; _rnd[556]=2.4920; _rnd[557]=-0.1367; _rnd[558]=-0.6627;
  _rnd[559]=0.0261; _rnd[560]=-0.1895; _rnd[561]=-1.2647; _rnd[562]=-1.4514; _rnd[563]=0.0378; _rnd[564]=-0.3925;
  _rnd[565]=0.3708; _rnd[566]=1.0713; _rnd[567]=-0.1759; _rnd[568]=-1.8599; _rnd[569]=1.0671; _rnd[570]=0.2425;
  _rnd[571]=-1.4763; _rnd[572]=0.1535; _rnd[573]=1.5930; _rnd[574]=-0.2828; _rnd[575]=0.8887; _rnd[576]=-1.1020;
  _rnd[577]=-0.9100; _rnd[578]=-1.4716; _rnd[579]=-0.6122; _rnd[580]=1.7982; _rnd[581]=0.1863; _rnd[582]=0.0333;
  _rnd[583]=-0.9980; _rnd[584]=-0.2636; _rnd[585]=0.0683; _rnd[586]=-1.3176; _rnd[587]=-2.1389; _rnd[588]=0.4678;
  _rnd[589]=-0.6368; _rnd[590]=0.4047; _rnd[591]=-1.2433; _rnd[592]=-0.0559; _rnd[593]=0.5002; _rnd[594]=2.0440;
  _rnd[595]=0.5904; _rnd[596]=0.5093; _rnd[597]=-1.1321; _rnd[598]=-0.1634; _rnd[599]=-1.4498; _rnd[600]=0.2668;
  _rnd[601]=-0.3987; _rnd[602]=-0.0431; _rnd[603]=0.9677; _rnd[604]=-0.4586; _rnd[605]=0.5841; _rnd[606]=-0.8937;
  _rnd[607]=1.2570; _rnd[608]=-2.0434; _rnd[609]=-0.0402; _rnd[610]=0.9443; _rnd[611]=-0.7386; _rnd[612]=0.3742;
  _rnd[613]=0.4464; _rnd[614]=0.1608; _rnd[615]=-1.0337; _rnd[616]=0.5304; _rnd[617]=0.6479; _rnd[618]=-1.4796;
  _rnd[619]=-0.3567; _rnd[620]=1.5978; _rnd[621]=0.6693; _rnd[622]=0.1334; _rnd[623]=0.8055; _rnd[624]=-0.0620;
  _rnd[625]=-0.3675; _rnd[626]=-1.3725; _rnd[627]=-0.7873; _rnd[628]=-0.4344; _rnd[629]=-1.8806; _rnd[630]=-2.0714;
  _rnd[631]=-0.4092; _rnd[632]=-0.7789; _rnd[633]=0.2411; _rnd[634]=-0.5105; _rnd[635]=0.0919; _rnd[636]=0.1013;
  _rnd[637]=-0.6059; _rnd[638]=-0.0309; _rnd[639]=-0.7635; _rnd[640]=0.0872; _rnd[641]=0.6013; _rnd[642]=1.7504;
  _rnd[643]=-0.6408; _rnd[644]=-0.6719; _rnd[645]=0.0275; _rnd[646]=-0.0031; _rnd[647]=-2.0722; _rnd[648]=0.6162;
  _rnd[649]=0.4184; _rnd[650]=-0.0374; _rnd[651]=0.1234; _rnd[652]=-0.2135; _rnd[653]=1.2331; _rnd[654]=-1.0843;
  _rnd[655]=0.8806; _rnd[656]=-2.1666; _rnd[657]=-0.5398; _rnd[658]=0.9286; _rnd[659]=0.3014; _rnd[660]=0.4800;
  _rnd[661]=-0.5976; _rnd[662]=-0.4269; _rnd[663]=1.3150; _rnd[664]=0.8482; _rnd[665]=-0.6813; _rnd[666]=-1.1107;
  _rnd[667]=-1.6175; _rnd[668]=1.6722; _rnd[669]=1.8953; _rnd[670]=-0.8533; _rnd[671]=0.2983; _rnd[672]=-0.0386;
  _rnd[673]=-1.0326; _rnd[674]=-0.5834; _rnd[675]=0.6276; _rnd[676]=1.6447; _rnd[677]=0.4408; _rnd[678]=-0.4669;
  _rnd[679]=-0.2625; _rnd[680]=-1.1679; _rnd[681]=-1.0743; _rnd[682]=-1.9064; _rnd[683]=-0.3052; _rnd[684]=-0.2980;
  _rnd[685]=0.8523; _rnd[686]=1.7074; _rnd[687]=0.7351; _rnd[688]=-1.4694; _rnd[689]=-1.3341; _rnd[690]=-2.2380;
  _rnd[691]=2.7290; _rnd[692]=-0.3787; _rnd[693]=-1.1834; _rnd[694]=0.4576; _rnd[695]=-1.2604; _rnd[696]=-1.0961;
  _rnd[697]=-0.1293; _rnd[698]=0.0269; _rnd[699]=1.1675; _rnd[700]=-0.6281; _rnd[701]=0.3886; _rnd[702]=-0.5455;
  _rnd[703]=-1.0303; _rnd[704]=0.7834; _rnd[705]=1.7829; _rnd[706]=-0.8545; _rnd[707]=0.3178; _rnd[708]=1.2923;
  _rnd[709]=1.5357; _rnd[710]=-0.8040; _rnd[711]=-0.5072; _rnd[712]=-0.0356; _rnd[713]=1.2770; _rnd[714]=0.1334;
  _rnd[715]=-1.0085; _rnd[716]=-1.6667; _rnd[717]=-0.0527; _rnd[718]=-3.0432; _rnd[719]=1.5277; _rnd[720]=0.0745;
  _rnd[721]=0.2186; _rnd[722]=-1.3375; _rnd[723]=0.2067; _rnd[724]=0.0352; _rnd[725]=-0.7738; _rnd[726]=-0.9940;
  _rnd[727]=-1.5946; _rnd[728]=-0.3765; _rnd[729]=-1.8789; _rnd[730]=0.6042; _rnd[731]=-0.0030; _rnd[732]=-1.7079;
  _rnd[733]=-0.3378; _rnd[734]=1.3221; _rnd[735]=0.7646; _rnd[736]=0.0123; _rnd[737]=0.9940; _rnd[738]=1.3251;
  _rnd[739]=-0.0951; _rnd[740]=0.8402; _rnd[741]=-1.0654; _rnd[742]=1.4734; _rnd[743]=-0.6086; _rnd[744]=-0.4007;
  _rnd[745]=-0.8824; _rnd[746]=0.2648; _rnd[747]=-1.0430; _rnd[748]=-0.8176; _rnd[749]=-0.8194; _rnd[750]=-0.1402;
  _rnd[751]=-0.1906; _rnd[752]=-0.5319; _rnd[753]=1.0442; _rnd[754]=0.6170; _rnd[755]=-0.9954; _rnd[756]=0.3305;
  _rnd[757]=-0.6467; _rnd[758]=-0.1778; _rnd[759]=0.2733; _rnd[760]=0.3384; _rnd[761]=-1.0562; _rnd[762]=-0.1843;
  _rnd[763]=0.5355; _rnd[764]=0.4702; _rnd[765]=-0.4092; _rnd[766]=-0.0453; _rnd[767]=-1.3264; _rnd[768]=-0.2130;
  _rnd[769]=0.5122; _rnd[770]=0.4723; _rnd[771]=0.4026; _rnd[772]=0.8022; _rnd[773]=-0.0829; _rnd[774]=0.6449;
  _rnd[775]=-0.0648; _rnd[776]=0.9649; _rnd[777]=-0.5386; _rnd[778]=-0.8334; _rnd[779]=-1.5222; _rnd[780]=-0.4840;
  _rnd[781]=-2.0166; _rnd[782]=-0.3542; _rnd[783]=3.0657; _rnd[784]=0.1365; _rnd[785]=-0.4580; _rnd[786]=0.4288;
  _rnd[787]=-0.3985; _rnd[788]=0.4223; _rnd[789]=-0.4883; _rnd[790]=-0.9012; _rnd[791]=0.7559; _rnd[792]=-0.8515;
  _rnd[793]=-1.5583; _rnd[794]=-0.1305; _rnd[795]=0.2308; _rnd[796]=0.1662; _rnd[797]=-0.7269; _rnd[798]=2.0569;
  _rnd[799]=-1.9836; _rnd[800]=-0.6485; _rnd[801]=-1.2795; _rnd[802]=-0.8758; _rnd[803]=0.2691; _rnd[804]=0.7055;
  _rnd[805]=0.7812; _rnd[806]=0.1472; _rnd[807]=1.3969; _rnd[808]=0.3340; _rnd[809]=-0.8044; _rnd[810]=-0.4335;
  _rnd[811]=0.3360; _rnd[812]=-0.5957; _rnd[813]=1.4145; _rnd[814]=0.9581; _rnd[815]=1.4618; _rnd[816]=0.0101;
  _rnd[817]=-0.5267; _rnd[818]=-1.2375; _rnd[819]=0.3619; _rnd[820]=-2.4133; _rnd[821]=1.0097; _rnd[822]=-1.4705;
  _rnd[823]=1.2123; _rnd[824]=-1.1130; _rnd[825]=-0.2740; _rnd[826]=0.1294; _rnd[827]=1.6230; _rnd[828]=-0.2718;
  _rnd[829]=-0.0089; _rnd[830]=-0.7158; _rnd[831]=0.0028; _rnd[832]=-1.8749; _rnd[833]=-1.9890; _rnd[834]=0.5647;
  _rnd[835]=-0.1198; _rnd[836]=0.9557; _rnd[837]=-1.3539; _rnd[838]=-0.9375; _rnd[839]=-0.2098; _rnd[840]=0.2455;
  _rnd[841]=2.8199; _rnd[842]=1.1964; _rnd[843]=1.4195; _rnd[844]=-0.5023; _rnd[845]=-2.4899; _rnd[846]=0.0945;
  _rnd[847]=-2.0580; _rnd[848]=-0.8490; _rnd[849]=-1.3591; _rnd[850]=-0.4318; _rnd[851]=-0.9867; _rnd[852]=-1.7031;
  _rnd[853]=-1.6191; _rnd[854]=2.2114; _rnd[855]=-1.5556; _rnd[856]=0.9839; _rnd[857]=-0.6331; _rnd[858]=2.3881;
  _rnd[859]=-0.1146; _rnd[860]=0.7587; _rnd[861]=1.2853; _rnd[862]=-0.6285; _rnd[863]=2.1710; _rnd[864]=1.9655;
  _rnd[865]=-0.8708; _rnd[866]=0.8074; _rnd[867]=0.0992; _rnd[868]=-0.4684; _rnd[869]=-1.0866; _rnd[870]=-0.4638;
  _rnd[871]=-1.2406; _rnd[872]=-2.0262; _rnd[873]=-0.0800; _rnd[874]=1.0894; _rnd[875]=0.0978; _rnd[876]=1.0198;
  _rnd[877]=-1.5346; _rnd[878]=-0.5992; _rnd[879]=1.6797; _rnd[880]=1.3828; _rnd[881]=0.4969; _rnd[882]=-1.3309;
  _rnd[883]=1.4391; _rnd[884]=-0.1544; _rnd[885]=0.0015; _rnd[886]=-1.7434; _rnd[887]=0.6656; _rnd[888]=-0.2817;
  _rnd[889]=0.9738; _rnd[890]=0.2789; _rnd[891]=2.8550; _rnd[892]=-0.4724; _rnd[893]=0.3643; _rnd[894]=0.2078;
  _rnd[895]=-0.3890; _rnd[896]=-0.5451; _rnd[897]=-1.2858; _rnd[898]=-1.6643; _rnd[899]=0.7655; _rnd[900]=-2.1660;
  _rnd[901]=-0.3191; _rnd[902]=-0.2404; _rnd[903]=0.4279; _rnd[904]=1.1741; _rnd[905]=1.1174; _rnd[906]=1.7610;
  _rnd[907]=1.5078; _rnd[908]=1.9920; _rnd[909]=0.5298; _rnd[910]=1.1952; _rnd[911]=1.7944; _rnd[912]=-0.6017;
  _rnd[913]=1.8417; _rnd[914]=-0.6887; _rnd[915]=0.4933; _rnd[916]=0.7427; _rnd[917]=-1.8009; _rnd[918]=0.4304;
  _rnd[919]=0.3772; _rnd[920]=0.7140; _rnd[921]=-1.2689; _rnd[922]=1.5347; _rnd[923]=-0.9440; _rnd[924]=-0.2566;
  _rnd[925]=-1.6047; _rnd[926]=0.5504; _rnd[927]=0.2689; _rnd[928]=-0.3642; _rnd[929]=1.1004; _rnd[930]=-1.1040;
  _rnd[931]=0.0668; _rnd[932]=-0.5187; _rnd[933]=0.4265; _rnd[934]=-0.6113; _rnd[935]=0.7172; _rnd[936]=0.7532;
  _rnd[937]=0.4954; _rnd[938]=-3.2072; _rnd[939]=-0.5708; _rnd[940]=-1.1468; _rnd[941]=0.7909; _rnd[942]=-0.5778;
  _rnd[943]=-0.6679; _rnd[944]=-0.0025; _rnd[945]=-1.8467; _rnd[946]=-0.8809; _rnd[947]=0.5863; _rnd[948]=0.1608;
  _rnd[949]=0.1744; _rnd[950]=0.0022; _rnd[951]=0.6691; _rnd[952]=-0.6760; _rnd[953]=-0.5031; _rnd[954]=0.9381;
  _rnd[955]=-0.5320; _rnd[956]=-1.0899; _rnd[957]=-0.4087; _rnd[958]=1.8167; _rnd[959]=0.9623; _rnd[960]=0.6563;
  _rnd[961]=-0.9592; _rnd[962]=0.6900; _rnd[963]=-0.1618; _rnd[964]=-0.2720; _rnd[965]=0.2310; _rnd[966]=1.5633;
  _rnd[967]=-1.3548; _rnd[968]=0.4832; _rnd[969]=0.8332; _rnd[970]=-1.5828; _rnd[971]=-0.8817; _rnd[972]=0.2759;
  _rnd[973]=0.5774; _rnd[974]=0.9016; _rnd[975]=0.0722; _rnd[976]=0.3466; _rnd[977]=-1.0456; _rnd[978]=0.0627;
  _rnd[979]=0.2091; _rnd[980]=0.4558; _rnd[981]=-2.4464; _rnd[982]=-0.6837; _rnd[983]=-1.9024; _rnd[984]=-1.0624;
  _rnd[985]=-1.3494; _rnd[986]=0.3033; _rnd[987]=-0.6359; _rnd[988]=1.9533; _rnd[989]=-0.1501; _rnd[990]=-0.5495;
  _rnd[991]=-1.7394; _rnd[992]=3.2566; _rnd[993]=1.5056; _rnd[994]=-0.6082; _rnd[995]=-0.3644; _rnd[996]=-1.8308;
  _rnd[997]=0.5356; _rnd[998]=-1.7947; _rnd[999]=0.4854;

  //Code snipet to generate the above code.
  /*
  int ppp=0;
  while(ppp<1000){
    cout<<Form("_rnd[%i]=%.4f; ",ppp,RandGauss::shoot(0,1));
    if(ppp!=0 && ppp%6==0){ cout<<endl; }
    ppp+=1;
  }
  cout<<endl;
  */
  return;
}
