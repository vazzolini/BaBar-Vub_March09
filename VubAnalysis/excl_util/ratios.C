gSystem.Load("libPhysics.so");
gSystem.Load("../RecoilAnalysis/libRecoilAnalysis.so");

void ratios(TString thefile1, TString thefile2, int type){

  char name[200];

  TFile data(thefile1.Data());
  TFile gene(thefile2.Data());
  
  double datalep[100], datamm2[100], dataall[100];
  double genelep[100], genemm2[100], geneall[100];
  double dataleperr[100], datamm2err[100], dataallerr[100];
  double geneleperr[100], genemm2err[100], geneallerr[100];
  double pur[100],purerr[100];

  double ratiodatamm2[100], ratiodataall[100];
  double ratiogenemm2[100], ratiogeneall[100];
  double ratiodatamm2err[100], ratiodataallerr[100];
  double ratiogenemm2err[100], ratiogeneallerr[100];

  double ratiolep[100], ratiomm2[100], ratioall[100];
  double ratioleperr[100], ratiomm2err[100], ratioallerr[100];

  recoilAnalysis b;
  mesData themes;

  int bias = type * 10000;

  double thesigma = 0.003;
  double themean = -1111111;
  double theargus = -1111111;
  double thealpha = -1111111;
  double then = 5;
  double resmean, ressigma, resalpha, resn;

  data.cd();
  for(int i=1; i<101; i++){
    double cut = i/200.; 
    pur[i-1]=cut;
    purerr[i-1]=0.001;
    sprintf(name, "h%d", 11000+i+bias); 
    mesData themes(*(b.vubMes((TH1D*)gDirectory->Get(name), resmean, ressigma, resalpha, resn, 1, 1, themean, thesigma, thealpha, then, theargus)));
    datalep[i-1]=themes.theSig(); dataleperr[i-1]=themes.theErrSig(); 
    if(dataleperr[i-1]/datalep[i-1]>.1) dataleperr[i-1] = sqrt(datalep[i-1])  ;
    sprintf(name, "h%d", 12000+i+bias); 
    mesData themes(*(b.vubMes((TH1D*)gDirectory->Get(name), resmean, ressigma, resalpha, resn, 1, 1, themean, thesigma, thealpha, then, theargus)));
    datamm2[i-1]=themes.theSig(); datamm2err[i-1]=themes.theErrSig(); 
    if(datamm2err[i-1]/datamm2[i-1]>.1) datamm2err[i-1] = sqrt(datamm2[i-1])  ;
    sprintf(name, "h%d", 13000+i+bias); 
    mesData themes(*(b.vubMes((TH1D*)gDirectory->Get(name), resmean, ressigma, resalpha, resn, 1, 1, themean, thesigma, thealpha, then, theargus)));
    dataall[i-1]=themes.theSig(); dataallerr[i-1]=themes.theErrSig(); 
    if(dataallerr[i-1]/dataall[i-1]>.1) dataallerr[i-1] = sqrt(dataall[i-1])  ;
  }

  gene.cd();
  for(int i=1; i<101; i++){
    sprintf(name, "h%d", 11000+i+bias); 
    mesData themes(*(b.vubMes((TH1D*)gDirectory->Get(name), resmean, ressigma, resalpha, resn, 1, 1, themean, thesigma, thealpha, then, theargus)));
    genelep[i-1]=themes.theSig(); geneleperr[i-1]=themes.theErrSig(); 
    if(geneleperr[i-1]/genelep[i-1]>.1) geneleperr[i-1] = sqrt(genelep[i-1])  ;
    sprintf(name, "h%d", 12000+i+bias); 
    mesData themes(*(b.vubMes((TH1D*)gDirectory->Get(name), resmean, ressigma, resalpha, resn, 1, 1, themean, thesigma, thealpha, then, theargus)));
    genemm2[i-1]=themes.theSig(); genemm2err[i-1]=themes.theErrSig(); 
    if(genemm2err[i-1]/genemm2[i-1]>.1) genemm2err[i-1] = sqrt(genemm2[i-1])  ;
    sprintf(name, "h%d", 13000+i+bias); 
    mesData themes(*(b.vubMes((TH1D*)gDirectory->Get(name), resmean, ressigma, resalpha, resn, 1, 1, themean, thesigma, thealpha, then, theargus)));
    geneall[i-1]=themes.theSig(); geneallerr[i-1]=themes.theErrSig(); 
    if(geneallerr[i-1]/geneall[i-1]>.1) geneallerr[i-1] = sqrt(geneall[i-1])  ;
  }

  for(int i=1; i<101; i++){

    ratiodatamm2[i-1] = datamm2[i-1]/datalep[i-1];
    ratiodataall[i-1] = dataall[i-1]/datalep[i-1];
    ratiodatamm2err[i-1] = datamm2err[i-1]/datalep[i-1];
    ratiodataallerr[i-1] = dataallerr[i-1]/datalep[i-1];

    ratiogenemm2[i-1] = genemm2[i-1]/genelep[i-1];
    ratiogeneall[i-1] = geneall[i-1]/genelep[i-1];
    ratiogenemm2err[i-1] = genemm2err[i-1]/genelep[i-1];
    ratiogeneallerr[i-1] = geneallerr[i-1]/genelep[i-1];

    ratiolep[i-1] = datalep[i-1]/genelep[i-1];
    ratiomm2[i-1] = datamm2[i-1]/genemm2[i-1];
    ratioall[i-1] = dataall[i-1]/geneall[i-1];
    ratioleperr[i-1] = sqrt( (datalep[i-1]*datalep[i-1]/(pow(genelep[i-1],4))) * geneleperr[i-1]*geneleperr[i-1] +  1/(pow(genelep[i-1],2)) * dataleperr[i-1]*dataleperr[i-1]  );
    ratiomm2err[i-1] = sqrt( (datamm2[i-1]*datamm2[i-1]/(pow(genemm2[i-1],4))) * genemm2err[i-1]*genemm2err[i-1] +  1/(pow(genemm2[i-1],2)) * datamm2err[i-1]*datamm2err[i-1]  );
    ratioallerr[i-1] = sqrt( (dataall[i-1]*dataall[i-1]/(pow(geneall[i-1],4))) * geneallerr[i-1]*geneallerr[i-1] +  1/(pow(geneall[i-1],2)) * dataallerr[i-1]*dataallerr[i-1]  );
        
  }

  char thetyp[200];
  sprintf (thetyp,"");
  if(type==1)  sprintf (thetyp,"B0");
  if(type==2)  sprintf (thetyp,"Bch");

  c0.Clear();
  sprintf (name,"%s%s",thetyp,"datadep.eps");
  TGraphErrors *results;
  results = new TGraphErrors(100,pur,datalep,purerr,dataleperr);
  results->SetMinimum(0.);
  results->SetMaximum(datalep[0]*1.3);
  //  results->SetTitle(title);
  results->SetMarkerSize(.5);
  results->Draw("AP");
  c0->Update();
  c0->GetFrame()->SetFillColor(0);
  c0->GetFrame()->SetBorderSize(12);
  results->GetHistogram()->SetTitleOffset(2., "Y");
  results->GetHistogram()->SetXTitle("purity");
  results->GetHistogram()->SetYTitle("yield");
  c0->Modified();
  c0.SaveAs(name);
  c0.Clear();

  c0.Clear();
  sprintf (name,"%s%s",thetyp,"datamm2.eps");
  TGraphErrors *results;
  results = new TGraphErrors(100,pur,datamm2,purerr,datamm2err);
  results->SetMinimum(0.);
  results->SetMaximum(datamm2[0]*1.3);
  //  results->SetTitle(title);
  results->SetMarkerSize(.5);
  results->Draw("AP");
  c0->Update();
  c0->GetFrame()->SetFillColor(0);
  c0->GetFrame()->SetBorderSize(12);
  results->GetHistogram()->SetTitleOffset(2., "Y");
  results->GetHistogram()->SetXTitle("purity");
  results->GetHistogram()->SetYTitle("yield");
  c0->Modified();
  c0.SaveAs(name);
  c0.Clear();

  c0.Clear();
  sprintf (name,"%s%s",thetyp,"dataall.eps");
  TGraphErrors *results;
  results = new TGraphErrors(100,pur,dataall,purerr,dataallerr);
  results->SetMinimum(0.);
  results->SetMaximum(dataall[0]*1.3);
  //  results->SetTitle(title);
  results->SetMarkerSize(.5);
  results->Draw("AP");
  c0->Update();
  c0->GetFrame()->SetFillColor(0);
  c0->GetFrame()->SetBorderSize(12);
  results->GetHistogram()->SetTitleOffset(2., "Y");
  results->GetHistogram()->SetXTitle("purity");
  results->GetHistogram()->SetYTitle("yield");
  c0->Modified();
  c0.SaveAs(name);
  c0.Clear();



  c0.Clear();
  sprintf (name,"%s%s",thetyp,"genedep.eps");
  TGraphErrors *results;
  results = new TGraphErrors(100,pur,genelep,purerr,geneleperr);
  results->SetMinimum(0.);
  results->SetMaximum(genelep[0]*1.3);
  //  results->SetTitle(title);
  results->SetMarkerSize(.5);
  results->Draw("AP");
  c0->Update();
  c0->GetFrame()->SetFillColor(0);
  c0->GetFrame()->SetBorderSize(12);
  results->GetHistogram()->SetTitleOffset(2., "Y");
  results->GetHistogram()->SetXTitle("purity");
  results->GetHistogram()->SetYTitle("yield");
  c0->Modified();
  c0.SaveAs(name);
  c0.Clear();

  c0.Clear();
  sprintf (name,"%s%s",thetyp,"genemm2.eps");
  TGraphErrors *results;
  results = new TGraphErrors(100,pur,genemm2,purerr,genemm2err);
  results->SetMinimum(0.);
  results->SetMaximum(genemm2[0]*1.3);
  //  results->SetTitle(title);
  results->SetMarkerSize(.5);
  results->Draw("AP");
  c0->Update();
  c0->GetFrame()->SetFillColor(0);
  c0->GetFrame()->SetBorderSize(12);
  results->GetHistogram()->SetTitleOffset(2., "Y");
  results->GetHistogram()->SetXTitle("purity");
  results->GetHistogram()->SetYTitle("yield");
  c0->Modified();
  c0.SaveAs(name);
  c0.Clear();

  c0.Clear();
  sprintf (name,"%s%s",thetyp,"geneall.eps");
  TGraphErrors *results;
  results = new TGraphErrors(100,pur,geneall,purerr,geneallerr);
  results->SetMinimum(0.);
  results->SetMaximum(geneall[0]*1.3);
  //  results->SetTitle(title);
  results->SetMarkerSize(.5);
  results->Draw("AP");
  c0->Update();
  c0->GetFrame()->SetFillColor(0);
  c0->GetFrame()->SetBorderSize(12);
  results->GetHistogram()->SetTitleOffset(2., "Y");
  results->GetHistogram()->SetXTitle("purity");
  results->GetHistogram()->SetYTitle("yield");
  c0->Modified();
  c0.SaveAs(name);
  c0.Clear();


  c0.Clear();
  sprintf (name,"%s%s",thetyp,"ratiodatamm2.eps");
  TGraphErrors *results;
  results = new TGraphErrors(100,pur,ratiodatamm2,purerr,ratiodatamm2err);
  results->SetMinimum(0.);
  results->SetMaximum(ratiodatamm2[50]*2.);
  //  results->SetTitle(title);
  results->SetMarkerSize(.5);
  results->Draw("AP");
  c0->Update();
  c0->GetFrame()->SetFillColor(0);
  c0->GetFrame()->SetBorderSize(12);
  results->GetHistogram()->SetTitleOffset(2., "Y");
  results->GetHistogram()->SetXTitle("purity");
  results->GetHistogram()->SetYTitle("ratio");
  c0->Modified();
  c0.SaveAs(name);
  c0.Clear();

  c0.Clear();
  sprintf (name,"%s%s",thetyp,"ratiodataall.eps");
  TGraphErrors *results;
  results = new TGraphErrors(100,pur,ratiodataall,purerr,ratiodataallerr);
  results->SetMinimum(0.);
  results->SetMaximum(ratiodataall[50]*2.);
  //  results->SetTitle(title);
  results->SetMarkerSize(.5);
  results->Draw("AP");
  c0->Update();
  c0->GetFrame()->SetFillColor(0);
  c0->GetFrame()->SetBorderSize(12);
  results->GetHistogram()->SetTitleOffset(2., "Y");
  results->GetHistogram()->SetXTitle("purity");
  results->GetHistogram()->SetYTitle("ratio");
  c0->Modified();
  c0.SaveAs(name);
  c0.Clear();

  c0.Clear();
  sprintf (name,"%s%s",thetyp,"ratiogenemm2.eps");
  TGraphErrors *results;
  results = new TGraphErrors(100,pur,ratiogenemm2,purerr,ratiogenemm2err);
  results->SetMinimum(0.);
  results->SetMaximum(ratiogenemm2[50]*2.);
  //  results->SetTitle(title);
  results->SetMarkerSize(.5);
  results->Draw("AP");
  c0->Update();
  c0->GetFrame()->SetFillColor(0);
  c0->GetFrame()->SetBorderSize(12);
  results->GetHistogram()->SetTitleOffset(2., "Y");
  results->GetHistogram()->SetXTitle("purity");
  results->GetHistogram()->SetYTitle("ratio");
  c0->Modified();
  c0.SaveAs(name);
  c0.Clear();

  c0.Clear();
  sprintf (name,"%s%s",thetyp,"ratiogeneall.eps");
  TGraphErrors *results;
  results = new TGraphErrors(100,pur,ratiogeneall,purerr,ratiogeneallerr);
  results->SetMinimum(0.);
  results->SetMaximum(ratiogeneall[50]*2.);
  //  results->SetTitle(title);
  results->SetMarkerSize(.5);
  results->Draw("AP");
  c0->Update();
  c0->GetFrame()->SetFillColor(0);
  c0->GetFrame()->SetBorderSize(12);
  results->GetHistogram()->SetTitleOffset(2., "Y");
  results->GetHistogram()->SetXTitle("purity");
  results->GetHistogram()->SetYTitle("ratio");
  c0->Modified();
  c0.SaveAs(name);
  c0.Clear();

  c0.Clear();
  sprintf (name,"%s%s",thetyp,"ratiolep.eps");
  TGraphErrors *results;
  results = new TGraphErrors(100,pur,ratiolep,purerr,ratioleperr);
  results->SetMinimum(0.);
  results->SetMaximum(ratiolep[50]*2.);
  //  results->SetTitle(title);
  results->SetMarkerSize(.5);
  results->Draw("AP");
  c0->Update();
  c0->GetFrame()->SetFillColor(0);
  c0->GetFrame()->SetBorderSize(12);
  results->GetHistogram()->SetTitleOffset(2., "Y");
  results->GetHistogram()->SetXTitle("purity");
  results->GetHistogram()->SetYTitle("ratio");
  c0->Modified();
  c0.SaveAs(name);
  c0.Clear();

  c0.Clear();
  sprintf (name,"%s%s",thetyp,"ratiomm2.eps");
  TGraphErrors *results;
  results = new TGraphErrors(100,pur,ratiomm2,purerr,ratiomm2err);
  results->SetMinimum(0.);
  results->SetMaximum(ratiomm2[50]*2.);
  //  results->SetTitle(title);
  results->SetMarkerSize(.5);
  results->Draw("AP");
  c0->Update();
  c0->GetFrame()->SetFillColor(0);
  c0->GetFrame()->SetBorderSize(12);
  results->GetHistogram()->SetTitleOffset(2., "Y");
  results->GetHistogram()->SetXTitle("purity");
  results->GetHistogram()->SetYTitle("ratio");
  c0->Modified();
  c0.SaveAs(name);
  c0.Clear();

  c0.Clear();
  sprintf (name,"%s%s",thetyp,"ratioall.eps");
  TGraphErrors *results;
  results = new TGraphErrors(100,pur,ratioall,purerr,ratioallerr);
  results->SetMinimum(0.);
  results->SetMaximum(ratioall[40]*2.);
  //  results->SetTitle(title);
  results->SetMarkerSize(.5);
  results->Draw("AP");
  c0->Update();
  c0->GetFrame()->SetFillColor(0);
  c0->GetFrame()->SetBorderSize(12);
  results->GetHistogram()->SetTitleOffset(2., "Y");
  results->GetHistogram()->SetXTitle("purity");
  results->GetHistogram()->SetYTitle("ratio");
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
