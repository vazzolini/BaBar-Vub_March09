void fitmatch(char *inname, char *histo, double BR = 1, int opt = 0, double mean = 0, double sig = 0.01){
  
  gROOT->SetStyle("Plain");
  c1 = new TCanvas("c1"," ",200,10,800,800);   
  c1->Clear();
  char name[100];
  sprintf(name,"%s%s%s","output_",inname,".root");
  TFile g(name);
  Int_t npar = 7;
  Double_t params[7] = {1,mean,sig,-.01,0.02,.8,1.};
  TF1 *f2 = new TF1("f2",sumgaus,-.01,0.01,npar);
  f2->SetParameters(params);
  f2->SetParLimits(0, 0., 10000.);
  f2->SetParLimits(1, -0.01, 0.3);
  f2->SetParLimits(2, 0., 0.03);
  f2->SetParLimits(3, -0.01, 0.3);
  f2->SetParLimits(4, 0., 0.05);
  f2->SetParLimits(5, 0., 1.);
  f2->SetParLimits(6, 0., 100.);
  if (opt==1) {
    f2->FixParameter(3,0);
    f2->FixParameter(4,1);
    f2->FixParameter(5,1);    
  }
  f2->SetLineColor(kRed);
  gStyle->SetOptStat(0);
  ((TH1D*)gDirectory->Get(histo))->SetMaximum(((TH1D*)gDirectory->Get(histo))->GetMaximum()*1.5);  
  ((TH1D*)gDirectory->Get(histo))->Fit(f2,"lmq","pe");
  ((TH1D*)gDirectory->Get(histo))->SetStats(0);
  double mean1 = f2->GetParameter(1);
  double sigma1 = f2->GetParameter(2);
  double errmean1 = f2->GetParError(1);
  double errsigma1 = f2->GetParError(2);
  double mean2 = f2->GetParameter(3);
  double sigma2 = f2->GetParameter(4);
  double errmean2 = f2->GetParError(3);
  double errsigma2 = f2->GetParError(4);
  double fraction = f2->GetParameter(5);
  double errfraction = f2->GetParError(5);
  double chisqua = f2->GetChisquare();
  double ndof = f2->GetNDF();
  double totarea = f2->GetParameter(0)/((TH1D*)gDirectory->Get(histo))->GetBinWidth(1);
  double errtotarea = f2->GetParError(0)/((TH1D*)gDirectory->Get(histo))->GetBinWidth(1);
  TLatex tl;
  tl.SetTextSize(.03);
  double x = 0.13;
  char myline[100];
  
  if(opt!=1){
    sprintf(myline, "area = %6.2f +/- %5.2f", totarea, errtotarea ); 
    tl.DrawTextNDC(x, 0.85, myline);    
    sprintf(myline, "mean1 = %6.4f +/- %5.4f", mean1, errmean1); 
    tl.DrawTextNDC(x, 0.81, myline);
    sprintf(myline, "sigma1 = %6.4f +/- %5.4f", sigma1, errsigma1); 
    tl.DrawTextNDC(x, 0.77, myline);
    sprintf(myline, "mean2 = %6.4f +/- %5.4f", mean2, errmean2); 
    tl.DrawTextNDC(x, 0.73, myline);
    sprintf(myline, "sigma2 = %6.4f +/- %5.4f", sigma2, errsigma2); 
    tl.DrawTextNDC(x, 0.69, myline);
    sprintf(myline, "area1/area2 = %6.4f +/- %5.4f", fraction, errfraction); 
    tl.DrawTextNDC(x, 0.65, myline);
    sprintf(myline, "chi2 = %7.2f ", chisqua/ndof);
    tl.DrawTextNDC(x, 0.61, myline);
  }else{
    sprintf(myline, "area = %6.2f +/- %5.2f", totarea, errtotarea ); 
    tl.DrawTextNDC(x, 0.85, myline);    
    sprintf(myline, "mean = %6.4f +/- %5.4f", mean1, errmean1); 
    tl.DrawTextNDC(x, 0.81, myline);
    sprintf(myline, "sigma = %6.4f +/- %5.4f", sigma1, errsigma1); 
    tl.DrawTextNDC(x, 0.77, myline);
    sprintf(myline, "chi2 = %7.2f ", chisqua/ndof);
    tl.DrawTextNDC(x, 0.73, myline);    
  }
  sprintf(name,"%s%s%s%s%s","fitted_",histo,"_",inname,".eps");
  c1.SaveAs(name);
  
  if(!strcmp(histo, "resophi")){
    double den = ((TH1D*)gDirectory->Get("entr"))->Integral();
    double num = totarea;
    double eff = num/(den*BR);
    double erreff = sqrt(eff*(1-eff)/(den*BR));
    cout<< endl;
    cout << "Efficiency for " << inname << " is = " << eff << " +/- " << erreff << endl; 
  }

}

Double_t gaussian(Double_t *x, Double_t *par)
// The signal function: a gaussian
{
   Double_t arg = 0;
   if (par[2]) arg = (x[0] - par[1])/par[2];

   Double_t sig = par[0]*TMath::Exp(-0.5*arg*arg);
   return sig;
}


Double_t sumgaus(Double_t *x, Double_t *par)
{

   Double_t arg = 0;

   if (par[2]) arg = (x[0] - par[1])/par[2];

   Double_t gaus1 = TMath::Exp(-0.5*arg*arg)/(sqrt(2*3.1416)*par[2]);

   arg = 0;

   if (par[4]) arg = (x[0] - par[3])/par[4];

   Double_t gaus2 = TMath::Exp(-0.5*arg*arg)/(sqrt(2*3.1416)*par[4]);

   Double_t tot = par[0]*((par[5])*gaus1+(1-par[5])*gaus2) + par[6];
   return tot;
}
