#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <sstream>
#include <stdlib.h>

void fitDRplot(const TString filename, TString outfile="double_ratiointp050", TString fncn = "pol2"){

  std::vector<double> bins;
  std::vector<double> dratio;
  std::vector<double> drerro;
  std::ifstream is(filename);
  if (!is) { 
    std::cout << "Error: Can't open file " << filename << ". Stopping execution!" << std::endl;
    exit(EXIT_FAILURE);
  } else {
    std::cout << " reading values from text file: " << filename << std::endl;
  }

  bins.resize(0);
  dratio.resize(0);
  drerro.resize(0);


  char buffer[200]; 
  float mx_l, mx_h, corr, err_corr;
  char pm[20];
  while(is.getline(buffer,200, '\n')) {
    //    if (buffer.empty()) continue;
    if (buffer[0] == '#') continue;
    if (buffer[0] == '/') continue;
    sscanf(buffer, "%f %f %f %s %f", &mx_l,&mx_h,&corr,&pm,&err_corr);
    bins.push_back(mx_l);
    dratio.push_back(corr);
    drerro.push_back(err_corr);
  }
  
  //get the upper limit of the histogram
  bins.push_back(mx_h);
  
  for (int i(0); i<bins.size(); ++i) {
    std::cout << " bin [" << i << "] = " << bins[i] << "   " << dratio[i] << " +/- " << drerro[i] << std::endl;
  }

  Int_t nBins = bins.size();
  TH1D *h = new TH1D("fit results",outfile,nBins-1,&bins[0]);

  for (int i(0); i<nBins; ++i){
    h->SetBinContent(i,dratio[i-1]);
    h->SetBinError(i,drerro[i-1]);
  }

  gStyle->SetOptFit(0111);
  gStyle->SetOptStat(000000011);
  gStyle->SetStatFontSize(0.009);
  h->SetMaximum(7.);
  h->SetMinimum(-1.);
  h->Draw();
  TF1 *fitfnc = new TF1("fitfnc",fncn,1.55,5);
  h->Fit("fitfnc","VEFMR");

  TF1 *fit = h->GetFunction("fitfnc");
  Int_t npars = fit->GetNpar();
  const int npar = npars;

  Double_t matr[npar][npar];
  gMinuit->mnemat(&matr[0][0],npar);
  if (!matr) {
    std::cout << "Error: NO COVARIANCE MATRIX! " << std::endl;
    return;
  }

  // prepare the output file
  char outt[50];
  sprintf(outt,"%s%s%s",&outfile[0],&fncn[0],".txt");
  FILE *output=fopen(outt,"a");
  fprintf(output,"#mx_l mx_h  corr      err_corr \n");

  //get value at center of bins
  for (Int_t i(0); i<nBins; ++i){
    float bincenter = (bins[i]+bins[i-1])/2.;
    if(bincenter == 0) continue;
    Double_t fitted = fit->Eval(bincenter);
    std::cout << "mx = " << bincenter << " input value " << dratio[i-1] << " fitted value: " << fitted << std::endl;

    //get error for each bin
    float errbin(0.);
    for (Int_t irow=0; irow<npar; irow++){
      float derirow = pow(bincenter,irow);
      for (Int_t icol=0; icol<npar; icol++){
	float dericol = pow(bincenter,icol);
      	float add = derirow*dericol*matr[irow][icol];
	//	float add = fit->GetParError(irow)*fit->GetParError(icol)*matr[irow][icol];
	std::cout << "row " << irow << " col " << icol << 
	  " err.par" << irow << " " << fit->GetParError(irow) <<
	  " err.par" << icol << " " << fit->GetParError(icol) <<
	  " CovMat " << matr[irow][icol] << 
	  " der.par" << irow << " " << derirow << 
	  " der. par" << icol << " " << dericol << 
	  " adding " << add << std::endl;
	errbin+= add; 
      }
    }

    std::cout << " ErrBin = " << sqrt(errbin) << std::endl;
    fprintf(output,"%3.2f\t %3.2f\t  %5.3f +- %5.3f \n",bins[i-1],bins[i],fitted,sqrt(errbin));

  }

  
  fclose(output);
  
  char out[50];
  sprintf(out,"%s%s%s",&outfile[0],&fncn[0],".eps");
  c1->Print(out);

}


