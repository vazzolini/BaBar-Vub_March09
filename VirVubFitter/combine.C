#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TArrow.h"
#include "TMinuit.h"
#include "TMatrixD.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "TMatrixD.h"

using namespace std;

double pi0lnuq28[6];
double pi0lnuq2816[6];
double pi0lnuq216[6];
                 
double pi0lnuq28scaled[6];
double pi0lnuq2816scaled[6];
double pi0lnuq216scaled[6];
                 
double pilnuq28[6];
double pilnuq2816[6];
double pilnuq216[6];
 
double combpilnuq28[5];
double combpilnuq2816[5];
double combpilnuq216[5];
 
double totalpi0lnu[5];
double totalpilnu[5];
double combtotalpilnu[5];


const double BALLpar[2] = {5.44,1.43};
const double hpqcdpar[2] = {1.29,0.32};
const double fnalpar[2] ={1.83,0.50};

const double BALLfull[2] ={7.74,2.32};
const double hpqcdfull[2] ={5.70,1.71};
const double fnalfull[2] ={6.24,2.12};

double covmat[6][6];
double invcovmat[6][6];
double covmatpie[6][6];

double q21,q21err;
double q22,q22err;
double q23,q23err;

double q21stat,q21errstat;
double q22stat,q22errstat;
double q23stat,q23errstat;

double q21errsyst;
double q22errsyst;
double q23errsyst;

// all results
int use[3]= {1,1,1};

void chi2Hist(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
  Double_t chi2(0.);

  double parfit[6];
  double valfit[6];
  
  valfit[0] = pi0lnuq28[0]         ;   parfit[0] = par[0]/1.84;
  valfit[1] = pi0lnuq2816[0]         ;   parfit[1] = par[1]/1.84;                
  valfit[2] = pi0lnuq216[0]         ;   parfit[2] = par[2]/1.84;         

  valfit[3] = pilnuq28[0]         ;   parfit[3] = par[0];
  valfit[4] = pilnuq2816[0]         ;   parfit[4] = par[1];                
  valfit[5] = pilnuq216[0]         ;   parfit[5] = par[2];         

  for (int i=0; i<6; i++){
    for (int j=0; j<6; j++){
      chi2 += (valfit[i] - parfit[i]) * invcovmat[i][j] *  (valfit[j] - parfit[j]);
    }
  }
  

  f=chi2;
  cout << chi2 << endl;
}


void
BABAR()
{
  // use the 'plain' style for plots (white backgrounds, etc)
  gROOT->SetStyle("Plain");

  TStyle *babarStyle= new TStyle("BABAR","BaBar approved plots style");

  // use plain black on white colors
  babarStyle->SetFrameBorderMode(0);
  babarStyle->SetCanvasBorderMode(0);
  babarStyle->SetPadBorderMode(0);
  babarStyle->SetPadColor(0);
  babarStyle->SetCanvasColor(0);
  //babarStyle->SetTitleColor(0);
  babarStyle->SetStatColor(0);
  babarStyle->SetTitleFillColor(0);

  // set the paper & margin sizes
  babarStyle->SetPaperSize(20,26);
  babarStyle->SetPadTopMargin(0.05);
  babarStyle->SetPadRightMargin(0.05);
  babarStyle->SetPadBottomMargin(0.16);
  babarStyle->SetPadLeftMargin(0.16);

  // use large Times-Roman fonts
  babarStyle->SetTitleFont(132,"xyz");  // set the all 3 axes title font
  babarStyle->SetTitleFont(132," ");    // set the pad title font
  babarStyle->SetTitleSize(0.06,"xyz"); // set the 3 axes title size
  babarStyle->SetTitleSize(0.06," ");   // set the pad title size
  babarStyle->SetLabelFont(132,"xyz");
  babarStyle->SetLabelSize(0.05,"xyz");
  babarStyle->SetTextFont(132);
  babarStyle->SetTextSize(0.08);
  babarStyle->SetStatFont(132);

  // use bold lines and markers
  babarStyle->SetMarkerStyle(8);
  babarStyle->SetHistLineWidth(2);
  babarStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  //..Get rid of X error bars
  //babarStyle->SetErrorX(0.001);

  // do not display any of the standard histogram decorations
  babarStyle->SetOptTitle(0);
  babarStyle->SetOptStat(0);
  babarStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  babarStyle->SetPadTickX(1);
  babarStyle->SetPadTickY(1);

  gROOT->SetStyle("BABAR");
}

// ----------------------------------------------------------------------
void dumpBR() {
  cout << endl;
  cout << "====================================" << endl;
  cout << endl;
  cout << "pi0lnuq28    :"  <<pi0lnuq28[0] << "  " <<pi0lnuq28[1] << "  " <<pi0lnuq28[2] << "  " <<pi0lnuq28[3] << "  " <<pi0lnuq28[4] <<   endl;  
  cout << "pi0lnuq2816  :"  <<pi0lnuq2816[0] << "  " <<pi0lnuq2816[1] << "  " <<pi0lnuq2816[2] << "  " <<pi0lnuq2816[3] << "  " <<pi0lnuq2816[4] <<   endl;  
  cout << "pi0lnuq28    :"  <<pi0lnuq216[0] << "  " <<pi0lnuq216[1] << "  " <<pi0lnuq216[2] << "  " <<pi0lnuq216[3] << "  " <<pi0lnuq216[4] <<   endl;  
  cout << endl;
  cout << "pilnuq28    :"  <<pilnuq28[0] << "  " <<pilnuq28[1] << "  " <<pilnuq28[2] << "  " <<pilnuq28[3] << "  " <<pilnuq28[4] <<   endl;  
  cout << "pilnuq2816  :"  <<pilnuq2816[0] << "  " <<pilnuq2816[1] << "  " <<pilnuq2816[2] << "  " <<pilnuq2816[3] << "  " <<pilnuq2816[4] <<   endl;  
  cout << "pilnuq28    :"  <<pilnuq216[0] << "  " <<pilnuq216[1] << "  " <<pilnuq216[2] << "  " <<pilnuq216[3] << "  " <<pilnuq216[4] <<   endl;  
}

void read(TString  file){
  char  buffer[200];
  sprintf(buffer, "%s", file.Data());
  ifstream is(buffer);
  char BRName[100], BREntry[100];
  float Value;
  int ok(0);
  int index;

  while (is.getline(buffer, 200, '\n')) {
    ok = 0;
    index = -1;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %s %f", BRName, BREntry, &Value);
    // -- 
    if (!strcmp(BREntry, "value")) {
      index = 0;
      Value = Value * 10000.;
    }
    if (!strcmp(BREntry, "stat")) {
      index = 1;
      Value = Value * 10000.;
    }
    if (!strcmp(BREntry, "statmc")) {
      index = 2;
      Value = Value * 10000.;
    }
    if (!strcmp(BREntry, "expected")) {
      index = 5;
      Value = Value * 10000.;
    }
    if (!strcmp(BREntry, "uncorrsyst")) index = 3;
    if (!strcmp(BREntry, "corrsyst")) index = 4;
    if (index<0)  cout << "==> WRONG VALUE " << BREntry <<  " ... prabably I am going to crash" << endl;

    double SLB0(0.1117);
    double SLBP(0.1029);
    if(!strcmp(BREntry, "uncorrsyst") || !strcmp(BREntry, "corrsyst")) {SLB0 = 1; SLBP = 1;}

    if (!strcmp(BRName, "pi0lnuq28")) { pi0lnuq28[index] = Value * SLBP; ok = 1;}
    if (!strcmp(BRName, "pi0lnuq2816")) { pi0lnuq2816[index] = Value * SLBP; ok = 1;}
    if (!strcmp(BRName, "pi0lnuq216")) { pi0lnuq216[index] = Value * SLBP; ok = 1;}

    if (!strcmp(BRName, "pilnuq28")) { pilnuq28[index] = Value * SLB0; ok = 1;}
    if (!strcmp(BRName, "pilnuq2816")) { pilnuq2816[index] = Value * SLB0; ok = 1;}
    if (!strcmp(BRName, "pilnuq216")) { pilnuq216[index] = Value * SLB0; ok = 1;}
  }

  dumpBR();

  for (int i=3;i<5;i++){    
    pi0lnuq28[i] =         pi0lnuq28[i] *       pi0lnuq28[0]/100.;
    pi0lnuq2816[i]  =      pi0lnuq2816[i] *     pi0lnuq2816[0]/100.;
    pi0lnuq216[i] =        pi0lnuq216[i] *      pi0lnuq216[0]/100.;

    pilnuq28[i] =         pilnuq28[i] *        pilnuq28[0]/100.;
    pilnuq2816[i]  =      pilnuq2816[i] *      pilnuq2816[0]/100.;
    pilnuq216[i] =        pilnuq216[i] *       pilnuq216[0]/100.;
 
  }
  
 
  for (int i=0;i<3;i++){    
    pi0lnuq28scaled[i] =         pi0lnuq28[i] *       1.84;
    pi0lnuq2816scaled[i]  =      pi0lnuq2816[i] *     1.84;
    pi0lnuq216scaled[i] =        pi0lnuq216[i] *      1.84;
  } 

}

void avecorr(double &ave, double &erruncorr, double &errcorr, double X1, double S1U, double S1C, double X2, double S2U, double S2C){
  
  double A = S1U*S1U + S1C*S1C;
  double B = S1C*S2C;
  double C = S1C*S2C;
  double D = S2U*S2U + S2C*S2C;

  double det = A*D - B*C;
  
  double invA = D/det;
  double invB = -C/det;
  double invC = -B/det;
  double invD = A/det;
  

  ave = X1*invA + X1*invB + X2*invC +X2*invD;
  ave = ave / (invA + invB + invC +invD);
  
  double errave = sqrt(1 / (invA + invB + invC +invD));
 
//   cout << invA << " " << invB << " " << invC << " " << invD <<  endl;
 
  erruncorr = sqrt(1 / ( 1/(S1U*S1U) + 1/(S2U*S2U)));
  errcorr = sqrt(errave*errave - erruncorr*erruncorr);

//    cout << "X1 = " << X1 << " +/- " << S1U << "(uncorr) +/- " << S1C << "(corr)" << endl;
//    cout << "X2 = " << X2 << " +/- " << S2U << "(uncorr) +/- " << S2C << "(corr)" << endl;
 

  
//   cout << "Xmean = " << ave << " +/- " << erruncorr << "(uncorr.) +/- " <<  errcorr << "(corr.)" << endl;
//   cout << endl;

}


void invertmat(){
  
  TMatrixD cov(6,6);

  cout << endl;
  cout << endl;
  for (int i=0; i<6; i++){
    for (int j=0; j<6; j++){
      cov(i,j) = covmat[i][j];
      cout << covmat[i][j] << "  " ;
    }
    cout << endl;
  }


  Double_t det1;
  TMatrixD invcov = cov;
  invcov.Invert(&det1);

  cout << endl;
  cout << endl;
  for (int i=0; i<6; i++){
    for (int j=0; j<6; j++){
      invcovmat[i][j] = invcov(i,j);
      cout << invcovmat[i][j] << "  " ;
    }
    cout << endl;
  }
  
  cout << endl;
  cout << endl;
}

double setelement(double* vector1, double* vector2, int e1, int e2, int bitword){

  double sigma1(0), sigma2(0);

  covmatpie[e1][e2] = bitword;
  if ( int(bitword/1       )%10 == 1 )   { sigma1 += pow(vector1[1],2); sigma2 += pow(vector2[1],2); }
  if ( int(bitword/10      )%10 == 1 )   { sigma1 += pow(vector1[2],2); sigma2 += pow(vector2[2],2); }
  if ( int(bitword/100     )%10 == 1 )   { sigma1 += pow(vector1[3],2); sigma2 += pow(vector2[3],2); }
  if ( int(bitword/1000    )%10 == 1 )   { sigma1 += pow(vector1[4],2); sigma2 += pow(vector2[4],2); }

  sigma1 = sqrt(sigma1);  
  sigma2 = sqrt(sigma2);

  double scaleforhfag = sqrt(vector1[5]*vector2[5])/sqrt(vector1[0]*vector2[0]);
  cout << e1 << " " << e2 << " " << scaleforhfag << endl;

  return sigma1*sigma2*scaleforhfag;

} 

void buildcov(){

  for (int i=0; i<6; i++){
    for (int j=0; j<6; j++){
      covmat[i][j] = 0;
    }
  }

  covmat[0][0]  = setelement( pi0lnuq28 , pi0lnuq28 ,0,0,  1111);
  covmat[0][1]  = setelement( pi0lnuq28 , pi0lnuq2816 ,0,1, 1100 );
  covmat[0][2]  = setelement( pi0lnuq28 , pi0lnuq216 ,0,2, 1100 );
  covmat[0][3]  = setelement( pi0lnuq28 , pilnuq28 ,0,3,  1000);
  covmat[0][4]  = setelement( pi0lnuq28 , pilnuq2816 ,0,4, 1000 );
  covmat[0][5]  = setelement( pi0lnuq28 , pilnuq216 ,0,5, 1000 );

  covmat[1][1]  = setelement( pi0lnuq2816 , pi0lnuq2816 ,1,1, 1111 );
  covmat[1][2]  = setelement( pi0lnuq2816 , pi0lnuq216 ,1,2, 1100 );
  covmat[1][3]  = setelement( pi0lnuq2816 , pilnuq28 ,1,3,  1000);
  covmat[1][4]  = setelement( pi0lnuq2816 , pilnuq2816 ,1,4, 1000 );
  covmat[1][5]  = setelement( pi0lnuq2816 , pilnuq216 ,1,5, 1000 );

  covmat[2][2]  = setelement( pi0lnuq216 , pi0lnuq216 ,2,2, 1111);
  covmat[2][3]  = setelement( pi0lnuq216 , pilnuq28 ,2,3,  1000);
  covmat[2][4]  = setelement( pi0lnuq216 , pilnuq2816 ,2,4, 1000 );
  covmat[2][5]  = setelement( pi0lnuq216 , pilnuq216 ,2,5, 1000 );

  covmat[3][3]  = setelement( pilnuq28 , pilnuq28 ,3,3,  1111);
  covmat[3][4]  = setelement( pilnuq28 , pilnuq2816 ,3,4, 1000 );
  covmat[3][5]  = setelement( pilnuq28 , pilnuq216 ,3,5, 1000 );
				
  covmat[4][4]  = setelement( pilnuq2816 , pilnuq2816 ,4,4, 1111 );
  covmat[4][5]  = setelement( pilnuq2816 , pilnuq216 ,4,5, 1000 );
				
  covmat[5][5]  = setelement( pilnuq216 , pilnuq216 ,5,5, 1111 );
  //

  for (int i=0; i<6; i++){
    for (int j=0; j<6; j++){
      if(i>j)
	covmat[i][j] = covmat[j][i];
    }
  }

  //  dumpcovmat();
}

void buildcov2(){

  for (int i=0; i<6; i++){
    for (int j=0; j<6; j++){
      covmat[i][j] = 0;
    }
  }

  covmat[0][0]  = setelement( pi0lnuq28 , pi0lnuq28 ,0,0,  1);
  covmat[0][1]  = setelement( pi0lnuq28 , pi0lnuq2816 ,0,1, 0 );
  covmat[0][2]  = setelement( pi0lnuq28 , pi0lnuq216 ,0,2, 0 );
  covmat[0][3]  = setelement( pi0lnuq28 , pilnuq28 ,0,3,  0);
  covmat[0][4]  = setelement( pi0lnuq28 , pilnuq2816 ,0,4, 0 );
  covmat[0][5]  = setelement( pi0lnuq28 , pilnuq216 ,0,5, 0 );

  covmat[1][1]  = setelement( pi0lnuq2816 , pi0lnuq2816 ,1,1, 1 );
  covmat[1][2]  = setelement( pi0lnuq2816 , pi0lnuq216 ,1,2, 0 );
  covmat[1][3]  = setelement( pi0lnuq2816 , pilnuq28 ,1,3,  0);
  covmat[1][4]  = setelement( pi0lnuq2816 , pilnuq2816 ,1,4, 0 );
  covmat[1][5]  = setelement( pi0lnuq2816 , pilnuq216 ,1,5, 0 );

  covmat[2][2]  = setelement( pi0lnuq216 , pi0lnuq216 ,2,2, 1);
  covmat[2][3]  = setelement( pi0lnuq216 , pilnuq28 ,2,3,  0);
  covmat[2][4]  = setelement( pi0lnuq216 , pilnuq2816 ,2,4, 0 );
  covmat[2][5]  = setelement( pi0lnuq216 , pilnuq216 ,2,5, 0 );

  covmat[3][3]  = setelement( pilnuq28 , pilnuq28 ,3,3,  1);
  covmat[3][4]  = setelement( pilnuq28 , pilnuq2816 ,3,4, 0 );
  covmat[3][5]  = setelement( pilnuq28 , pilnuq216 ,3,5, 0 );
				
  covmat[4][4]  = setelement( pilnuq2816 , pilnuq2816 ,4,4, 1 );
  covmat[4][5]  = setelement( pilnuq2816 , pilnuq216 ,4,5, 0 );
				
  covmat[5][5]  = setelement( pilnuq216 , pilnuq216 ,5,5, 1 );
  //

  for (int i=0; i<6; i++){
    for (int j=0; j<6; j++){
      if(i>j)
	covmat[i][j] = covmat[j][i];
    }
  }

  //  dumpcovmat();
}

void combBR(){



  double temperrpi0(0.);
  double temperrpi(0.);
  double foomean, foocor;

  // full error
  temperrpi0 = sqrt(pow(pi0lnuq28[1],2) + pow(pi0lnuq28[2],2) + pow(pi0lnuq28[3],2));
  temperrpi  = sqrt(pow(pilnuq28[1],2) + pow(pilnuq28[2],2) + pow(pilnuq28[3],2));
  avecorr(combpilnuq28[0], combpilnuq28[3], combpilnuq28[4], pi0lnuq28[0]*1.84, temperrpi0*1.84, pi0lnuq28[4]*1.84, pilnuq28[0], temperrpi, pilnuq28[4]);

  // corr syst
  temperrpi0 = sqrt(pow(pi0lnuq28[1],2) + pow(pi0lnuq28[2],2));
  temperrpi  = sqrt(pow(pilnuq28[1],2) + pow(pilnuq28[2],2));
  avecorr(foomean, combpilnuq28[2], foocor, pi0lnuq28[0]*1.84, temperrpi0*1.84, 0., pilnuq28[0], temperrpi, 0.);
  combpilnuq28[3] = sqrt(pow(combpilnuq28[3],2) - pow(combpilnuq28[2],2));

  // MC error and stat error
  avecorr(foomean, combpilnuq28[1], foocor, pi0lnuq28[0]*1.84, pi0lnuq28[1]*1.84, 0., pilnuq28[0], pilnuq28[1], 0.);
  combpilnuq28[2] = sqrt(pow(combpilnuq28[2],2) - pow(combpilnuq28[1],2));


  // full error
  temperrpi0 = sqrt(pow(pi0lnuq2816[1],2) + pow(pi0lnuq2816[2],2) + pow(pi0lnuq2816[3],2));
  temperrpi  = sqrt(pow(pilnuq2816[1],2) + pow(pilnuq2816[2],2) + pow(pilnuq2816[3],2));
  avecorr(combpilnuq2816[0], combpilnuq2816[3], combpilnuq2816[4], pi0lnuq2816[0]*1.84, temperrpi0*1.84, pi0lnuq2816[4]*1.84, pilnuq2816[0], temperrpi, pilnuq2816[4]);

  // corr syst
  temperrpi0 = sqrt(pow(pi0lnuq2816[1],2) + pow(pi0lnuq2816[2],2));
  temperrpi  = sqrt(pow(pilnuq2816[1],2) + pow(pilnuq2816[2],2));
  avecorr(foomean, combpilnuq2816[2], foocor, pi0lnuq2816[0]*1.84, temperrpi0*1.84, 0., pilnuq2816[0], temperrpi, 0.);
  combpilnuq2816[3] = sqrt(pow(combpilnuq2816[3],2) - pow(combpilnuq2816[2],2));

  // MC error and stat error
  avecorr(foomean, combpilnuq2816[1], foocor, pi0lnuq2816[0]*1.84, pi0lnuq2816[1]*1.84, 0., pilnuq2816[0], pilnuq2816[1], 0.);
  combpilnuq2816[2] = sqrt(pow(combpilnuq2816[2],2) - pow(combpilnuq2816[1],2));


  // full error
  temperrpi0 = sqrt(pow(pi0lnuq216[1],2) + pow(pi0lnuq216[2],2) + pow(pi0lnuq216[3],2));
  temperrpi  = sqrt(pow(pilnuq216[1],2) + pow(pilnuq216[2],2) + pow(pilnuq216[3],2));
  avecorr(combpilnuq216[0], combpilnuq216[3], combpilnuq216[4], pi0lnuq216[0]*1.84, temperrpi0*1.84, pi0lnuq216[4]*1.84, pilnuq216[0], temperrpi, pilnuq216[4]);

  // corr syst
  temperrpi0 = sqrt(pow(pi0lnuq216[1],2) + pow(pi0lnuq216[2],2));
  temperrpi  = sqrt(pow(pilnuq216[1],2) + pow(pilnuq216[2],2));
  avecorr(foomean, combpilnuq216[2], foocor, pi0lnuq216[0]*1.84, temperrpi0*1.84, 0., pilnuq216[0], temperrpi, 0.);
  combpilnuq216[3] = sqrt(pow(combpilnuq216[3],2) - pow(combpilnuq216[2],2));

  // MC error and stat error
  avecorr(foomean, combpilnuq216[1], foocor, pi0lnuq216[0]*1.84, pi0lnuq216[1]*1.84, 0., pilnuq216[0], pilnuq216[1], 0.);
  combpilnuq216[2] = sqrt(pow(combpilnuq216[2],2) - pow(combpilnuq216[1],2));



  totalpi0lnu[0] = pi0lnuq28[0] + pi0lnuq2816[0] + pi0lnuq216[0];
  totalpilnu[0] = pilnuq28[0] + pilnuq2816[0] + pilnuq216[0];
  combtotalpilnu[0] = combpilnuq28[0] + combpilnuq2816[0] + combpilnuq216[0];

  totalpi0lnu[1] = sqrt(pow(pi0lnuq28[1],2) + pow(pi0lnuq2816[1],2) + pow(pi0lnuq216[1],2));
  totalpi0lnu[2] = sqrt(pow(pi0lnuq28[2],2) + pow(pi0lnuq2816[2],2) + pow(pi0lnuq216[2],2));
  totalpi0lnu[4] = sqrt(pow(pi0lnuq28[3] + pi0lnuq2816[3] + pi0lnuq216[3],2)+ pow(pi0lnuq28[4] + pi0lnuq2816[4] + pi0lnuq216[4],2));

  totalpilnu[1] = sqrt(pow(pilnuq28[1],2) + pow(pilnuq2816[1],2) + pow(pilnuq216[1],2));
  totalpilnu[2] = sqrt(pow(pilnuq28[2],2) + pow(pilnuq2816[2],2) + pow(pilnuq216[2],2));
  totalpilnu[4] = sqrt(pow(pilnuq28[3] + pilnuq2816[3] + pilnuq216[3],2) + pow(pilnuq28[4] + pilnuq2816[4] + pilnuq216[4],2));

  combtotalpilnu[1] = sqrt(pow(combpilnuq28[1],2) + pow(combpilnuq2816[1],2) + pow(combpilnuq216[1],2));
  combtotalpilnu[2] = sqrt(pow(combpilnuq28[2],2) + pow(combpilnuq2816[2],2) + pow(combpilnuq216[2],2));
  combtotalpilnu[4] = combpilnuq28[3] + combpilnuq2816[3] + combpilnuq216[3] + combpilnuq28[4] + combpilnuq2816[4] + combpilnuq216[4];

  cout << endl;
  cout << endl;
  cout << "total BR(pi0lnu)     =  " << totalpi0lnu[0] << " +/- " << totalpi0lnu[1] << "(stat) +/- " << totalpi0lnu[2] << "(MC stat) +/- " << totalpi0lnu[4] << "(syst)" << endl;
  cout << "total BR(pilnu)      =  " << totalpilnu[0] << " +/- " << totalpilnu[1] << "(stat) +/- " << totalpilnu[2] << "(MC stat) +/- " << totalpilnu[4] << "(syst)" << endl;
  cout << endl;
  cout << endl;
  cout << "SIMPLE COMBINATION METHOD" << endl;
  cout << endl;
  cout << "COMBO BR(pilnu,q2<8)    =  " << combpilnuq28[0] << " +/- " << combpilnuq28[1] << "(stat) +/- " << combpilnuq28[2] << "(MC stat) +/- " << combpilnuq28[3] << "(syst uncorr) +/- " << combpilnuq28[4] << "(syst corr)" << endl;
  cout << "COMBO BR(pilnu,8<q2<16) =  " << combpilnuq2816[0] << " +/- " << combpilnuq2816[1] << "(stat) +/- " << combpilnuq2816[2] << "(MC stat) +/- " << combpilnuq2816[3] << "(syst uncorr) +/- " << combpilnuq2816[4] << "(syst corr)" << endl;
  cout << "COMBO BR(pilnu,16<q2)   =  " << combpilnuq216[0] << " +/- " << combpilnuq216[1] << "(stat) +/- " << combpilnuq216[2] << "(MC stat) +/- " << combpilnuq216[3] << "(syst uncorr) +/- " << combpilnuq216[4] << "(syst corr)" << endl;
  cout << "COMBO total BR(pilnu)   =  " << combtotalpilnu[0] << " +/- " << combtotalpilnu[1] << "(stat) +/- " << combtotalpilnu[2] << "(MC stat) +/- " << combtotalpilnu[4] << "(syst)"<< endl;

  buildcov();
  invertmat();

  TMinuit aMinuit(3);
  int ierflg;
  // initialize minuit
  aMinuit.SetFCN(chi2Hist);
  aMinuit.mnparm(0, "q21",  .1, 0.0001, 0., 5.0, ierflg); 
  aMinuit.mnparm(1, "q22",  .2, 0.0001, 0., 5.0, ierflg); 
  aMinuit.mnparm(2, "q23",  .3, 0.0001, 0., 5.0, ierflg); 

  Double_t arglis[5];
  arglis[0] = 100000; // maxcalls
  arglis[1] = 0.0001;  // tolerance

  aMinuit.mnexcm("MIGRAD", arglis, 2, ierflg);
  aMinuit.mnexcm("MINOS", arglis, 2, ierflg);

  aMinuit.GetParameter(0,q21,q21err);
  aMinuit.GetParameter(1,q22,q22err);
  aMinuit.GetParameter(2,q23,q23err);


  TMatrixD invmatres(3,3);
  
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      invmatres(i,j) =  pow(1/1.84,2)*invcovmat[i][j] + invcovmat[i+3][j+3] + (1/1.84)*(invcovmat[i][j+3]+invcovmat[j][i+3]);
      //      cout << invmatres(i,j) << "  " ;
    }
    //    cout << endl;
    
  }
  Double_t det2;
  TMatrixD matres = invmatres;
  matres.Invert(&det2);
 
  double totalBRerr = 0;
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      totalBRerr += matres(i,j);
      //      cout << totalBRerr << endl;
    }
  }
  totalBRerr = sqrt(totalBRerr);
 
  

  buildcov2();
  invertmat();

  TMinuit aMinuit2(3);
  int ierflg;
  // initialize minuit
  aMinuit2.SetFCN(chi2Hist);
  aMinuit2.mnparm(0, "q21",  .1, 0.0001, 0., 5.0, ierflg); 
  aMinuit2.mnparm(1, "q22",  .2, 0.0001, 0., 5.0, ierflg); 
  aMinuit2.mnparm(2, "q23",  .3, 0.0001, 0., 5.0, ierflg); 

  Double_t arglis[5];
  arglis[0] = 100000; // maxcalls
  arglis[1] = 0.0001;  // tolerance

  aMinuit2.mnexcm("MIGRAD", arglis, 2, ierflg);
  aMinuit2.mnexcm("MINOS", arglis, 2, ierflg);

  aMinuit2.GetParameter(0,q21stat,q21errstat);
  aMinuit2.GetParameter(1,q22stat,q22errstat);
  aMinuit2.GetParameter(2,q23stat,q23errstat);

  q21errsyst = sqrt(pow(q21err,2)-pow(q21errstat,2));
  q22errsyst = sqrt(pow(q22err,2)-pow(q22errstat,2));
  q23errsyst = sqrt(pow(q23err,2)-pow(q23errstat,2));

  cout << endl;
  cout << endl;
  cout << "FULL METHOD" << endl;
  cout << endl;
  cout << "COMBO BR(pilnu,q2<8)       =  " << q21 << " +/- " << q21errstat << "(stat) +/- " << q21errsyst << "(syst corr)" << endl;
  cout << "COMBO BR(pilnu,8<q2<16)    =  " << q22 << " +/- " << q22errstat << "(stat) +/- " << q22errsyst << "(syst corr)" << endl;
  cout << "COMBO BR(pilnu,16<q2)      =  " << q23 << " +/- " << q23errstat << "(stat) +/- " << q23errsyst << "(syst corr)" << endl;
  
  double totalBR    = q21 + q22 + q23;

  TMatrixD invmatres2(3,3);
  
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      invmatres2(i,j) =   pow(1/1.84,2)*invcovmat[i][j] +invcovmat[i+3][j+3] + (1/1.84)*(invcovmat[i][j+3]+invcovmat[j][i+3]);
      //      cout << invmatres2(i,j) << "  " ;
    }
    //    cout << endl;
    
  }
  Double_t det3;
  TMatrixD matres2 = invmatres2;
  matres2.Invert(&det3);
 
  double totalBRerrstat = 0;
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      totalBRerrstat += matres2(i,j);
      //      cout << totalBRerrstat << endl;
    }
  }
  totalBRerrstat = sqrt(totalBRerrstat);

  double totalerrsyst = sqrt(pow(totalBRerr,2)-pow(totalBRerrstat,2));

  cout << endl;
  cout << "COMBO total BR             =  " << totalBR << " +/- " << totalBRerrstat << "(stat) +/- " << totalerrsyst << "(syst corr)" << endl;
  cout << endl;
  
  

  // Vub

  double tb0 = 1.536;
  
  double vubfnalfull          = sqrt(totalBR/(fnalfull[0]*tb0));
  double errstatvubfnalfull   = (sqrt((totalBR+totalBRerrstat)/(fnalfull[0]*tb0))-sqrt((totalBR-totalBRerrstat)/(fnalfull[0]*tb0)))/2;
  double errsystvubfnalfull   = (sqrt((totalBR+totalerrsyst)/(fnalfull[0]*tb0))-sqrt((totalBR-totalerrsyst)/(fnalfull[0]*tb0)))/2;
  double errthvubfnalfull_1   = sqrt(totalBR/((fnalfull[0]-fnalfull[1])*tb0))-vubfnalfull;
  double errthvubfnalfull_2   = sqrt(totalBR/((fnalfull[0]+fnalfull[1])*tb0))-vubfnalfull;

  cout << endl;
  cout << endl;
  cout << "Vub (FNAL full) = " << vubfnalfull << " +/- " << errstatvubfnalfull << "(stat.) " << " +/- " << errsystvubfnalfull << "(syst.) " << " +" << errthvubfnalfull_1 <<errthvubfnalfull_2 << "(th.syst) " << endl;

  double vubBALLfull          = sqrt(totalBR/(BALLfull[0]*tb0));
  double errstatvubBALLfull   = (sqrt((totalBR+totalBRerrstat)/(BALLfull[0]*tb0))-sqrt((totalBR-totalBRerrstat)/(BALLfull[0]*tb0)))/2;
  double errsystvubBALLfull   = (sqrt((totalBR+totalerrsyst)/(BALLfull[0]*tb0))-sqrt((totalBR-totalerrsyst)/(BALLfull[0]*tb0)))/2;
  double errthvubBALLfull_1   = sqrt(totalBR/((BALLfull[0]-BALLfull[1])*tb0))-vubBALLfull;
  double errthvubBALLfull_2   = sqrt(totalBR/((BALLfull[0]+BALLfull[1])*tb0))-vubBALLfull;

  cout << "Vub (BALL full) = " << vubBALLfull << " +/- " << errstatvubBALLfull << "(stat.) " << " +/- " << errsystvubBALLfull << "(syst.) " << " +" << errthvubBALLfull_1 <<errthvubBALLfull_2 << "(th.syst) " << endl;

  double vubhpqcdfull          = sqrt(totalBR/(hpqcdfull[0]*tb0));
  double errstatvubhpqcdfull   = (sqrt((totalBR+totalBRerrstat)/(hpqcdfull[0]*tb0))-sqrt((totalBR-totalBRerrstat)/(hpqcdfull[0]*tb0)))/2;
  double errsystvubhpqcdfull   = (sqrt((totalBR+totalerrsyst)/(hpqcdfull[0]*tb0))-sqrt((totalBR-totalerrsyst)/(hpqcdfull[0]*tb0)))/2;
  double errthvubhpqcdfull_1   = sqrt(totalBR/((hpqcdfull[0]-hpqcdfull[1])*tb0))-vubhpqcdfull;
  double errthvubhpqcdfull_2   = sqrt(totalBR/((hpqcdfull[0]+hpqcdfull[1])*tb0))-vubhpqcdfull;

  cout << "Vub (HPQCD full) = " << vubhpqcdfull << " +/- " << errstatvubhpqcdfull << "(stat.) " << " +/- " << errsystvubhpqcdfull << "(syst.) " << " +" << errthvubhpqcdfull_1 <<errthvubhpqcdfull_2 << "(th.syst) " << endl;



  double vubfnalpar          = sqrt(q23/(fnalpar[0]*tb0));
  double errstatvubfnalpar   = (sqrt((q23+q23errstat)/(fnalpar[0]*tb0))-sqrt((q23-q23errstat)/(fnalpar[0]*tb0)))/2;
  double errsystvubfnalpar   = (sqrt((q23+q23errsyst)/(fnalpar[0]*tb0))-sqrt((q23-q23errsyst)/(fnalpar[0]*tb0)))/2;
  double errthvubfnalpar_1   = sqrt(q23/((fnalpar[0]-fnalpar[1])*tb0))-vubfnalpar;
  double errthvubfnalpar_2   = sqrt(q23/((fnalpar[0]+fnalpar[1])*tb0))-vubfnalpar;

  cout << "Vub (FNAL par) = " << vubfnalpar << " +/- " << errstatvubfnalpar << "(stat.) " << " +/- " << errsystvubfnalpar << "(syst.) " << " +" << errthvubfnalpar_1 <<errthvubfnalpar_2 << "(th.syst) " << endl;

  double vubBALLpar          = sqrt((q21+q22)/(BALLpar[0]*tb0));
  double errstatvubBALLpar   = (sqrt((q21+q22+q21errstat+q22errstat)/(BALLpar[0]*tb0))-sqrt((q21+q22-combpilnuq28[1]-combpilnuq2816[1])/(BALLpar[0]*tb0)))/2;
  double errsystvubBALLpar   = (sqrt((q21+q22+sqrt(pow(q21errsyst,2)+pow(q22errsyst,2)))/(BALLpar[0]*tb0))-
				sqrt((q21+q22-sqrt(pow(q21errsyst,2)+pow(q22errsyst,2)))/(BALLpar[0]*tb0)))/2;
  double errthvubBALLpar_1   = sqrt((q21+q22)/((BALLpar[0]-BALLpar[1])*tb0))-vubBALLpar;
  double errthvubBALLpar_2   = sqrt((q21+q22)/((BALLpar[0]+BALLpar[1])*tb0))-vubBALLpar;

  cout << "Vub (BALL par) = " << vubBALLpar << " +/- " << errstatvubBALLpar << "(stat.) " << " +/- " << errsystvubBALLpar << "(syst.) " << " +" << errthvubBALLpar_1 <<errthvubBALLpar_2 << "(th.syst) " << endl;


  double vubhpqcdpar          = sqrt(q23/(hpqcdpar[0]*tb0));
  double errstatvubhpqcdpar   = (sqrt((q23+q23errstat)/(hpqcdpar[0]*tb0))-sqrt((q23-q23errstat)/(hpqcdpar[0]*tb0)))/2;
  double errsystvubhpqcdpar   = (sqrt((q23+q23errsyst)/(hpqcdpar[0]*tb0))-sqrt((q23-q23errsyst)/(hpqcdpar[0]*tb0)))/2;
  double errthvubhpqcdpar_1   = sqrt(q23/((hpqcdpar[0]-hpqcdpar[1])*tb0))-vubhpqcdpar;
  double errthvubhpqcdpar_2   = sqrt(q23/((hpqcdpar[0]+hpqcdpar[1])*tb0))-vubhpqcdpar;

  cout << "Vub (HPQCD par) = " << vubhpqcdpar << " +/- " << errstatvubhpqcdpar << "(stat.) " << " +/- " << errsystvubhpqcdpar << "(syst.) " << " +" << errthvubhpqcdpar_1 <<errthvubhpqcdpar_2 << "(th.syst) " << endl;
  cout << endl;
  cout << endl;


}


void
pilnuall(char *file)
{
  
  read(file);
  combBR();

  BABAR();

  TLatex* latex = new TLatex;
  latex->SetNDC(true);
  //char hname[100];

  TCanvas* can = new TCanvas("can","can",400,400);

  // draw the box first
  TH1F* hbox = new TH1F("hbox","box",1,0,26.4);
  hbox->SetXTitle("q^{2} (GeV^{2}/c^{4})");
  //hbox->SetLabelSize(0.1);
  TGaxis::SetMaxDigits(3);
  hbox->SetTitleSize(0.05,"X");
  hbox->SetMaximum(11);
  hbox->Draw("AXIS");
  latex->SetTextSize(0.05);
  latex->SetTextAngle(90);
  latex->SetTextAlign(33);
  latex->DrawLatex(0.0,0.95,"dB(B^{0} #rightarrow #pi^{-}l^{+}#nu)/dq^{2} (10^{-6}/GeV^{2}/c^{4})");

  // BF results
  const double brnu[5][4] = { // BABAR v-reco
    { 0.30, 0.05, 0.05, 0.02 },
    { 0.32, 0.05, 0.03, 0.02 },
    { 0.23, 0.05, 0.03, 0.01 },
    { 0.27, 0.05, 0.02, 0.02 },
    { 0.26, 0.03, 0.04, 0.02 }};
  const double brsl[3][4] = { // BABAR sl
    { 0.48, 0.17, 0.06, 0.00 },
    { 0.34, 0.16, 0.07, 0.00 },
    { 0.21, 0.14, 0.06, 0.01 }};
  const double brhdpi[3][4] = { // BABAR Breco
    { pilnuq28[0], sqrt(pow(pilnuq28[1],2)+pow(pilnuq28[2],2)), pilnuq28[3], pilnuq28[4]},
    { pilnuq2816[0], sqrt(pow(pilnuq2816[1],2)+pow(pilnuq2816[2],2)),pilnuq2816[3], pilnuq2816[4] },
    { pilnuq216[0], sqrt(pow(pilnuq216[1],2)+pow(pilnuq216[2],2)),pilnuq216[3], pilnuq216[4] }};
  const double brhdpi0[3][4] = { // BABAR Breco
    { pi0lnuq28[0], sqrt(pow(pi0lnuq28[1],2)+pow(pi0lnuq28[2],2)), pi0lnuq28[3], pi0lnuq28[4] },
    { pi0lnuq2816[0], sqrt(pow(pi0lnuq2816[1],2)+pow(pi0lnuq2816[2],2)), pi0lnuq2816[3], pi0lnuq2816[4] },
    { pi0lnuq216[0], sqrt(pow(pi0lnuq216[1],2)+pow(pi0lnuq216[2],2)), pi0lnuq216[3], pi0lnuq216[4] }};
  const double brhd[3][4] = { // BABAR Breco
    { q21, q21errstat, q21errsyst, 0},
    { q22, q22errstat, q22errsyst, 0},
    { q23, q23errstat, q23errsyst, 0}};

  double x[5],y[5],exl[5],exh[5],eyl[5],eyh[5];
  TGraphAsymmErrors* gbr[5];
  const double bin3[4] = { 0, 8, 16, 26.4 };
  const double bin6[4] = { 0, 8, 16, 26.4 };
  const double bin5[6] = { 0, 5, 10, 15, 20, 25 };
  const double scale = 100.0;
  // BABAR v-reco
  for (int i=0; i<5; i++) {
    double wbin = bin5[i+1]-bin5[i];
    x[i] = bin5[i] + 0.5*wbin;
    y[i] = scale*brnu[i][0]/wbin;
    exl[i] = exh[i] = 0.5*wbin;
    eyl[i] = eyh[i] = scale*
      sqrt(pow(brnu[i][1],2)+pow(brnu[i][2],2)+pow(brnu[i][3],2))/wbin;
  }
  gbr[0] = new TGraphAsymmErrors(5,x,y,exl,exh,eyl,eyh);
  // BABAR sl
  for (int i=0; i<3; i++) {
    double wbin = bin3[i+1]-bin3[i];
    double offset = -0.3;
    x[i] = bin3[i] + 0.5*wbin + offset;
    y[i] = scale*brsl[i][0]/wbin;
    exl[i] = 0.5*wbin + offset;
    exh[i] = 0.5*wbin - offset;
    eyl[i] = eyh[i] = scale*
      sqrt(pow(brsl[i][1],2)+pow(brsl[i][2],2)+pow(brsl[i][3],2))/wbin;
  }
  gbr[1] = new TGraphAsymmErrors(3,x,y,exl,exh,eyl,eyh);
  // BABAR Breco
  for (int i=0; i<3; i++) {
    double wbin = bin6[i+1]-bin6[i];
    double offset = +0.3;
    x[i] = bin6[i] + 0.5*wbin + offset;
    y[i] = scale*brhd[i][0]/wbin;
    exl[i] = 0.5*wbin + offset;
    exh[i] = 0.5*wbin - offset;
    eyl[i] = eyh[i] = scale*
      sqrt(pow(brhd[i][1],2)+pow(brhd[i][2],2)+pow(brhd[i][3],2))/wbin;
  }
  gbr[2] = new TGraphAsymmErrors(3,x,y,exl,exh,eyl,eyh);

  // BABAR Breco pi0
  for (int i=0; i<3; i++) {
    double wbin = bin6[i+1]-bin6[i];
    double offset = +0.3;
    x[i] = bin6[i] + 0.5*wbin + offset;
    y[i] = scale*brhdpi0[i][0]/wbin;
    exl[i] = 0.5*wbin + offset;
    exh[i] = 0.5*wbin - offset;
    eyl[i] = eyh[i] = scale*
      sqrt(pow(brhdpi0[i][1],2)+pow(brhdpi0[i][2],2)+pow(brhdpi0[i][3],2))/wbin;
  }
  gbr[3] = new TGraphAsymmErrors(3,x,y,exl,exh,eyl,eyh);

  // BABAR Breco pi
  for (int i=0; i<3; i++) {
    double wbin = bin6[i+1]-bin6[i];
    double offset = +0.3;
    x[i] = bin6[i] + 0.5*wbin + offset;
    y[i] = scale*brhdpi[i][0]/wbin;
    exl[i] = 0.5*wbin + offset;
    exh[i] = 0.5*wbin - offset;
    eyl[i] = eyh[i] = scale*
      sqrt(pow(brhdpi[i][1],2)+pow(brhdpi[i][2],2)+pow(brhdpi[i][3],2))/wbin;
  }
  gbr[4] = new TGraphAsymmErrors(3,x,y,exl,exh,eyl,eyh);


  const int color[5] = { 4, 2, 1 , 6, 38};
  const int type[5] = { 24,25 , 0 , 0, 0};
  for (int i=0; i<3; i++) {
    gbr[i]->SetMarkerStyle(20);
    gbr[i]->SetMarkerSize(1.5);
    gbr[i]->SetMarkerStyle(type[i]);
    gbr[i]->SetMarkerColor(color[i]);
    gbr[i]->SetLineColor(color[i]);
    gbr[i]->SetLineWidth(2);
    gbr[i]->Draw("P");
  }
  TLegend* legend = new TLegend(0.55,0.72,0.91,0.91);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->AddEntry(gbr[0],"BABAR #nu reco.","p");
  legend->AddEntry(gbr[1],"BABAR sl. tag","p");
  legend->AddEntry(gbr[2],"BABAR had. tag","p");

  // Draw the label 
  double xpos = 0.43;
  double ypos = 0.88;
  double scale = .7; 
  double scale2 = .5; 
  TString str="preliminary";
  TLatex *babar = new TLatex();
  Double_t cheburashkaFactorX=1, cheburashkaFactorY=1, padSizeX=500, padSizeY=500, xpos2, ypos2, xposL;
  babar->SetNDC(kTRUE);
  babar->SetTextFont(32); // Bold-Italic Times
  babar->SetTextAlign(31); // Right-Bottom
  padSizeX = gPad->GetWw()*gPad->GetWNDC(); // Get pad's dimensions
  padSizeY = gPad->GetWh()*gPad->GetHNDC();
  if (padSizeX>padSizeY) cheburashkaFactorX=padSizeY/padSizeX;
  if (padSizeX<padSizeY) cheburashkaFactorY=padSizeX/padSizeY;
  //xpos2=xpos-0.185*scale*cheburashkaFactorX;
  xpos2=xpos-0.188*scale*cheburashkaFactorX;
  ypos2=ypos-0.0620*scale*cheburashkaFactorY;
  xpos2=xpos-0.253*scale*cheburashkaFactorX;
  babar->SetTextSize(0.10*scale); // Beginning to draw "BaBar"
  babar->DrawText(xpos2,ypos2,"B");
  babar->SetTextSize(0.075*scale);
  babar->DrawText(xpos2+0.039*scale*cheburashkaFactorX,ypos2,"A");
  babar->SetTextSize(0.10*scale);
  babar->DrawText(xpos2+0.1015*scale*cheburashkaFactorX,ypos2,"B");
  babar->SetTextSize(0.075*scale);
  babar->DrawText(xpos2+0.1875*scale*cheburashkaFactorX,ypos2,"AR");
  babar->SetTextFont(42); // Helvetica (medium, upright) 
  babar->SetTextSize(0.1*scale2);
  babar->SetTextAlign(33); // Right-Top
  babar->DrawLatex(xpos,ypos2-0.02*scale2*cheburashkaFactorY,str);

  legend->Draw();

  can->Print("pilnuall.eps");

  can->Clear();

  TLatex *text = new TLatex();
  hbox->Draw("AXIS");
  latex->DrawLatex(0.0,0.95,"dB(B^{0} #rightarrow #pi^{-}l^{+}#nu)/dq^{2} (10^{-6}/GeV^{2}/c^{4})");

  gbr[2]->SetMarkerStyle(20);
  gbr[2]->SetMarkerSize(1.5);
  gbr[2]->SetMarkerColor(kBlack);
  gbr[2]->SetLineColor(kBlack);
  gbr[2]->SetLineWidth(2);
  gbr[2]->Draw("P");

  TLegend* legend2 = new TLegend(0.65,0.72,0.91,0.91);
  legend2->SetBorderSize(0);
  legend2->SetFillColor(0);
  legend2->AddEntry(gbr[2],"combo (#pi^{0}l#nu, #pi^{+}l#nu)","p");
  //  legend2->Draw();
  text->SetTextSize(0.07);
  text->DrawLatex(13.,9.5.,"Combination");
  text->SetTextSize(0.045);
  text->DrawLatex(13.,8.5,"B^{+}#rightarrow#pi^{0}l^{+}#nu & B^{0}#rightarrow#pi^{-}l^{+}#nu");
  text->SetTextSize(0.07);

  babar->SetNDC(kTRUE);
  babar->SetTextFont(32); // Bold-Italic Times
  babar->SetTextAlign(31); // Right-Bottom
  babar->SetTextSize(0.10*scale); // Beginning to draw "BaBar"
  babar->DrawText(xpos2,ypos2,"B");
  babar->SetTextSize(0.075*scale);
  babar->DrawText(xpos2+0.039*scale*cheburashkaFactorX,ypos2,"A");
  babar->SetTextSize(0.10*scale);
  babar->DrawText(xpos2+0.1015*scale*cheburashkaFactorX,ypos2,"B");
  babar->SetTextSize(0.075*scale);
  babar->DrawText(xpos2+0.1875*scale*cheburashkaFactorX,ypos2,"AR");
  babar->SetTextFont(42); // Helvetica (medium, upright) 
  babar->SetTextSize(0.1*scale2);
  babar->SetTextAlign(33); // Right-Top
  babar->DrawLatex(xpos,ypos2-0.02*scale2*cheburashkaFactorY,str);

  can->Print("combopilnu.eps");

  can->Clear();

  hbox->Draw("AXIS");
  latex->DrawLatex(0.0,0.95,"dB(B^{+} #rightarrow #pi^{0}l^{+}#nu)/dq^{2} (10^{-6}/GeV^{2}/c^{4})");

  gbr[3]->SetMarkerStyle(20);
  gbr[3]->SetMarkerSize(1.5);
//   gbr[3]->SetMarkerColor(color[3]);
//   gbr[3]->SetLineColor(color[3]);
  gbr[3]->SetLineWidth(2);
  gbr[3]->Draw("P");


  text->DrawLatex(13.,9.,"B^{+}#rightarrow#pi^{0}l^{+}#nu");
  //  legend3->Draw();

  babar->SetNDC(kTRUE);
  babar->SetTextFont(32); // Bold-Italic Times
  babar->SetTextAlign(31); // Right-Bottom
  babar->SetTextSize(0.10*scale); // Beginning to draw "BaBar"
  babar->DrawText(xpos2,ypos2,"B");
  babar->SetTextSize(0.075*scale);
  babar->DrawText(xpos2+0.039*scale*cheburashkaFactorX,ypos2,"A");
  babar->SetTextSize(0.10*scale);
  babar->DrawText(xpos2+0.1015*scale*cheburashkaFactorX,ypos2,"B");
  babar->SetTextSize(0.075*scale);
  babar->DrawText(xpos2+0.1875*scale*cheburashkaFactorX,ypos2,"AR");
  babar->SetTextFont(42); // Helvetica (medium, upright) 
  babar->SetTextSize(0.1*scale2);
  babar->SetTextAlign(33); // Right-Top
  babar->DrawLatex(xpos,ypos2-0.02*scale2*cheburashkaFactorY,str);

  can->Print("pi0lnu.eps");

  can->Clear();

  hbox->Draw("AXIS");
  latex->DrawLatex(0.0,0.95,"dB(B^{0} #rightarrow #pi^{-}l^{+}#nu)/dq^{2} (10^{-6}/GeV^{2}/c^{4})");

  gbr[4]->SetMarkerStyle(20);
  gbr[4]->SetMarkerSize(1.5);
//   gbr[4]->SetMarkerColor(color[4]);
//   gbr[4]->SetLineColor(color[4]);
  gbr[4]->SetLineWidth(2);
  gbr[4]->Draw("P");

  TLegend* legend3 = new TLegend(0.65,0.72,0.91,0.91);
  legend3->SetBorderSize(0);
  legend3->SetFillColor(0);
  legend3->AddEntry(gbr[4],"#pi^{+}lnu","p");
  //  legend3->Draw();
  text->DrawLatex(13.,9.,"B^{0}#rightarrow#pi^{-}l^{+}#nu");

  babar->SetNDC(kTRUE);
  babar->SetTextFont(32); // Bold-Italic Times
  babar->SetTextAlign(31); // Right-Bottom
  babar->SetTextSize(0.10*scale); // Beginning to draw "BaBar"
  babar->DrawText(xpos2,ypos2,"B");
  babar->SetTextSize(0.075*scale);
  babar->DrawText(xpos2+0.039*scale*cheburashkaFactorX,ypos2,"A");
  babar->SetTextSize(0.10*scale);
  babar->DrawText(xpos2+0.1015*scale*cheburashkaFactorX,ypos2,"B");
  babar->SetTextSize(0.075*scale);
  babar->DrawText(xpos2+0.1875*scale*cheburashkaFactorX,ypos2,"AR");
  babar->SetTextFont(42); // Helvetica (medium, upright) 
  babar->SetTextSize(0.1*scale2);
  babar->SetTextAlign(33); // Right-Top
  babar->DrawLatex(xpos,ypos2-0.02*scale2*cheburashkaFactorY,str);

  can->Print("pilnu.eps");


}
