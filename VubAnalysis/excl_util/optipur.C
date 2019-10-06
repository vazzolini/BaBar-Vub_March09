gSystem.Load("libPhysics.so");
gSystem.Load("../RecoilAnalysis/libRecoilAnalysis.so");

void optipur(TString thefile1, char* type){

  char name[200];

  TFile data(thefile1.Data());
  
  double pur[100],purerr[100];
  double dataall[100];
  double dataallerr[100];
  double backall[100];
  double backallerr[100];
  double ssbdataall[100], ssbdataallerr[100];

  recoilAnalysis b;
  mesData themes;

  double thesigma = 0.003;
  double themean = -1111111;
  double theargus = -1111111;
  double thealpha = -1111111;
  double then = 5;
  double resmean, ressigma, resalpha, resn;
  
  double ssbmax=0;

  data.cd();
  for(int i=1; i<101; i++){
    double cut = i/200.; 
    pur[i-1]=cut;
    purerr[i-1]=0.001;
    sprintf(name, "h%d", 11000+i); 
    mesData themes(*(b.vubMes((TH1D*)gDirectory->Get(name), resmean, ressigma, resalpha, resn, 1, 1, themean, thesigma, thealpha, then, theargus)));
    dataall[i-1]=themes.theSig(); dataallerr[i-1]=themes.theErrSig(); 
    if(dataallerr[i-1]/dataall[i-1]>.1) dataallerr[i-1] = sqrt(dataall[i-1])  ;
    backall[i-1]=themes.theBg(); backallerr[i-1]=themes.theErrBg(); 
    ssbdataall[i-1]= dataall[i-1]/sqrt(dataall[i-1]+backall[i-1]);
    ssbdataallerr[i-1]= dataallerr[i-1]/sqrt(dataall[i-1]+backall[i-1]);
    if(ssbdataall[i-1]>ssbmax) ssbmax=ssbdataall[i-1];
  }


  char thetyp[200];
  sprintf (thetyp,type);

  c0.Clear();
  sprintf (name,"%s%s",thetyp,"datassb.eps");
  TGraphErrors *results;
  results = new TGraphErrors(100,pur,ssbdataall,purerr,ssbdataallerr);
  results->SetMinimum(0.);
  results->SetMaximum(ssbmax*1.3);
  results->SetMarkerSize(1.);
  results->Draw("AP");
  c0->Update();
  c0->GetFrame()->SetFillColor(0);
  c0->GetFrame()->SetBorderSize(12);
  results->GetHistogram()->SetTitleOffset(1., "Y");
  results->GetHistogram()->SetXTitle("purity");
  results->GetHistogram()->SetYTitle("stat signif.");
  c0->Modified();
  c0.SaveAs(name);
  c0.Clear();



}


Double_t allfunction(Double_t x)
{
   Double_t xx =x;
   Double_t theomb =  thep4(xx,3.94762e-01+0.074,-8.49734e-02,-1.99211e-01, 1.06066e-01, -1.49296e-02,1.) ;
   Double_t theoa = thep4(xx, 1.89914e+00+0.05,-3.05645e+00,1.81981e+00, -4.73468e-01, 4.53754e-02,1.);
   //  Double_t stat = thep4(xx,10.79, -30.99, 33.47, -15.97, 2.842, 1.2);
   Double_t stat = 0;
   Double_t systBDdecays = thep4(xx, -5.67, 11.96, -8.34, 1.94, 0., 1.);
   Double_t othsyst = thesyst(xx,0.18);
   return sqrt(theomb*theomb+theoa*theoa+stat*stat+systBDdecays*systBDdecays+othsyst*othsyst);
   //return 0;
}

Double_t myfunction(Double_t *x, Double_t *par)
{
  
   Float_t xx =x[0];
   Double_t f = thep4(xx, par[0], par[1], par[2], par[3], par[4], par[5]);
   return f;
}

Double_t SSB(Double_t *x, Double_t *par)
{  
  Float_t xx =x[0];
  Double_t f = thesyst(xx,par[0]);
  return f;
}

Double_t thesyst(Double_t thex, double par0)
{
  Double_t back = -1700. + 4018. * thex - 3179. * thex * thex + 867. * thex * thex * thex;
  Double_t sig = -284. + 326. * thex;   
  Double_t f;
  if((back+sig)==0) {f = 0;}
  else{f = par0*back/(back+sig);}
  return f;
}

Double_t thep4(Double_t thex, double par0, double par1, double par2, double par3, double par4, double scale)
{
   Double_t f = scale * TMath::Abs(par0 + par1 * thex + par2 * thex * thex + par3 * thex * thex * thex +  par4 * thex * thex * thex * thex);
   return f;
}

// ----------------------------------------------------------------------
void SetTitles(TH1 *h, const char *sx, const char *sy, float size, 
               float xoff, float yoff, float lsize, int font) {
  if (h == 0) {
    cout << " Histogram not defined" << endl;
  } else {
    h->SetXTitle(sx);                  h->SetYTitle(sy); 
    h->SetTitleOffset(xoff, "x");      h->SetTitleOffset(yoff, "y");
    h->SetTitleSize(size, "x");        h->SetTitleSize(size, "y");
    h->SetLabelSize(lsize, "x");       h->SetLabelSize(lsize, "y");
    h->SetLabelFont(font, "x");        h->SetLabelFont(font, "y");
    h->GetXaxis()->SetTitleFont(font); h->GetYaxis()->SetTitleFont(font);
    h->SetNdivisions(508, "X");
  }
}
