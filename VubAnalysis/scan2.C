gSystem.Load("libPhysics.so");
gSystem.Load("../shlib/Linux2/libRecoilAnalysis.so");

void scan2(TString thefile, char *thename, double ref, char * type = ""){
   float var, errvar, errvarMC;
   char thevar[300];
   double BRBR[300], errBRBR[300], errBRBRMC[300], errsystBRBR[300];
   double diffBRBR[300], errdiffBRBR[300];
   double optBRBR[300], erroptBRBR[300];
   double vubcomp[300], errvubcomp[300];
   double vcbcomp[300], errvcbcomp[300];
   double othcomp[300], errothcomp[300];
   double data1bin[300], errdata1bin[300];
   double vub1bin[300], errvub1bin[300], errvub1binMC[300];
   double vcb1bin[300], errvcb1bin[300], errvcb1binMC[300];
   double oth1bin[300], erroth1bin[300], erroth1binMC[300];
   double epsu[300], errepsuMC[300];
   double epsmx[300], errepsmxMC[300];
   double epstot[300], errepstotMC[300];
   double nsl[300];
   double nslmc[300];
   double fact[300];
   double pstarfact[300];
   double chisq[300];

   double thescan[300],errthescan[300];

   char buffer[300];
   char name[300];
   char namefile[300];
   char namein[300];
   char axistit[300];
   char title[300], temptit[300];
   int counter = 0;
   int thecounter = 0;
   recoilAnalysis r;
   sprintf(namein,thefile.Data());
   ifstream infile(namein);
   cout << "FILE " << namein << endl;

   TH1D diff("diff","diff",15,1.1,1.9);
   TH1D pulls("pulls","pulls",30,-4.,4.);
   TH1D brbr("brbr","brbr",30,0.020,0.032);
   pulls.SetStats(0);
   brbr.SetStats(0);

   while (infile.getline(buffer, 200, '\n')) {
     if (buffer[0] == '#') {continue;}
     //     if (buffer[0] == '/') {continue;}
     sscanf(buffer, "%s", &namefile);

     ifstream infile2(namefile);
     cout << namefile << endl;

     int thecounter = 0;
     while (infile2.getline(buffer, 200, '\n')&&thecounter<40) {

       if (buffer[0] == '#') {continue;}
       sscanf(buffer, "%s %f %f %f ", &thevar, &var, &errvar, &errvarMC);
       if(!strcmp(thename,thevar)){
 	thescan[counter] = var;
	errthescan[counter] = 0.00001;
	if(!strcmp(thename,"q2Cut")){
	  thescan[counter] = 2 * thescan[counter];
	}	
       }
       cout << thevar << endl;
       if(!strcmp(thevar,"BRBR")){
 	BRBR[counter]=var;
	// 	errBRBRMC[counter]=errvarMC;
 	errBRBRMC[counter]=0;	
 	errBRBR[counter]=errvar;
       }
       if(!strcmp(thevar,"vubcomp")){
 	vubcomp[counter]=var;
 	errvubcomp[counter]=errvar;
       }
       if(!strcmp(thevar,"vcbcomp")){
 	vcbcomp[counter]=var;
 	errvcbcomp[counter]=errvar;
       }
       if(!strcmp(thevar,"othcomp")){
 	othcomp[counter]=var;
 	errothcomp[counter]=errvar;
       }
       if(!strcmp(thevar,"data1bin")){
 	data1bin[counter]=var;
 	errdata1bin[counter]=errvar;
       }
       if(!strcmp(thevar,"vub1bin")){
 	vub1bin[counter]=var;
 	errvub1bin[counter]=errvar;
 	errvub1binMC[counter]=errvarMC;
       }
       if(!strcmp(thevar,"vcb1bin")){
 	vcb1bin[counter]=var;
 	errvcb1bin[counter]=errvar;
 	errvcb1binMC[counter]=errvarMC;
       }
       if(!strcmp(thevar,"oth1bin")){
 	oth1bin[counter]=var;
 	erroth1bin[counter]=errvar;
 	erroth1binMC[counter]=errvarMC;
       }
       if(!strcmp(thevar,"epsu")){
 	epsu[counter]=var;
 	errepsuMC[counter]=errvar;
       }
       if(!strcmp(thevar,"epstot")){
 	epstot[counter]=var;
 	errepstotMC[counter]=errvar;
       }
       if(!strcmp(thevar,"epsmx")){
 	epsmx[counter]=var;
 	errepsmxMC[counter]=errvar;
       }
       if(!strcmp(thevar,"nsl")){
 	nsl[counter]=var;
       }
       if(!strcmp(thevar,"nslmc")){
 	nslmc[counter]=var;
       }
       if(!strcmp(thevar,"fact")){
 	fact[counter]=var;
       }
       if(!strcmp(thevar,"pstarfact")){
 	pstarfact[counter]=var;
       }
       if(!strcmp(thevar,"chisq")){
 	chisq[counter]=var;
       }
       thecounter++;
     }
     counter++;
   }

   // reference error and central value

   double refvalue,errrefvalue,differrrefvalue,errsystrefvalue;
   double refa, refmb, refsyst, refoth;
   double min = 1000.;
   double max = -1000.;
   int iref;
   for (int i=0; i<counter; i++){
     if(TMath::Abs(thescan[i]-ref)<0.0001) {
       refvalue=BRBR[i];
       // referece stat (+MC) error
       errrefvalue=sqrt(errBRBR[i]*errBRBR[i]+errBRBRMC[i]*errBRBRMC[i]);

       // reference syst error calculation
       refmb = thep4(thescan[i],3.94762e-01+0.074,-8.49734e-02,-1.99211e-01, 1.06066e-01, -1.49296e-02,1.) * refvalue ;  // syst mb
       refa =  thep4(thescan[i], 1.89914e+00+0.05,-3.05645e+00,1.81981e+00, -4.73468e-01, 4.53754e-02,1.) * refvalue;    // syst a
       refsyst = thep4(thescan[i], -5.67, 11.96, -8.34, 1.94, 0., 1.) * refvalue;                                        // syst B-D 
       refoth = thesyst(thescan[i],0.18) * refvalue;                                                                     // syst exp (trk,neu,ecc...)
       iref = i; 
     }
     if(thescan[i]>max) max=thescan[i];
     if(thescan[i]<min) min=thescan[i];  
   }


   // difference and its error
  
   double optvalue;
   double diffa, diffmb, diffsyst, diffoth, difftot;

   for (int i=0; i<counter; i++){
     diffBRBR[i] = BRBR[i]-refvalue;
     
     // syst error difference
     diffmb = TMath::Abs(thep4(thescan[i],3.94762e-01+0.074,-8.49734e-02,-1.99211e-01, 1.06066e-01, -1.49296e-02,1.)*refvalue - refmb); // syst mb		     
     diffa =  TMath::Abs(thep4(thescan[i], 1.89914e+00+0.05,-3.05645e+00,1.81981e+00, -4.73468e-01, 4.53754e-02,1.)*refvalue - refa);	// syst a		     
     diffsyst =  TMath::Abs(thep4(thescan[i], -5.67, 11.96, -8.34, 1.94, 0., 1.)*refvalue - refsyst);					// syst B-D 		     
     diffoth = TMath::Abs(thesyst(thescan[i],0.18)*refvalue - refoth);									// syst exp (trk,neu,ecc...)
     // adding in quadrature
     difftot = sqrt(diffmb*diffmb + diffa*diffa + diffsyst*diffsyst + diffoth*diffoth);


     // here uncorrelated stat error!!!!
     int which(1);// 0: original 1: complete uncorrelated
     if(which==0) {
       errdiffBRBR[i] = sqrt(TMath::Abs(errBRBR[i]*errBRBR[i]+errBRBRMC[i]*errBRBRMC[i]+difftot*difftot-errrefvalue*errrefvalue));
     } else {
       if(erroth1bin[i]>errdata1bin[i]) erroth1bin[i]=0;
       double dnd=sqrt(TMath::Abs(pow(errdata1bin[i],2)-pow(errdata1bin[iref],2)));
       double nb=vcb1bin[i]+oth1bin[i];
       double dnb=sqrt(TMath::Abs(pow(errvcb1bin[i],2)-pow(errvcb1bin[iref],2)+pow(erroth1bin[i],2)-pow(erroth1bin[iref],2)));
       errdiffBRBR[i] =sqrt(dnd*dnd+dnb*dnb)*BRBR[i]/(data1bin[i]-nb);
       //RF 12 feb 03 ;remove theo uncertainty from diff
       //       if(!strcmp(thename,"mxCut")){
       //	 errdiffBRBR[i] =sqrt(pow(errdiffBRBR[i],2) + pow(difftot,2));
       //       }
       //epstot[i]/pstarfact[i]/nsl[i]/fact[i];
       cout <<" diff "<<diffBRBR[i]<<" +/- " << errdiffBRBR[i]<<endl;
     }

     cout << BRBR[i]<<" +/- " << errBRBR[i]<<endl;
     cout << "HERE " << endl;
     optBRBR[i] = errBRBR[i]/BRBR[i];
     erroptBRBR[i] = (optBRBR[i]/40.);
     if(TMath::Abs(thescan[i]-ref)<0.0001) {
       optvalue = diffBRBR[i];
     }
   }


  
   for(int i=0;i<counter; i++){
     cout << BRBR[i]<<  thescan[i]<<endl;
     pulls.Fill((BRBR[i]-0.0185)/errBRBR[i]); 
     brbr.Fill(BRBR[i]); 
     double sigmas = 0;
     if(errdiffBRBR[i]) sigmas = diffBRBR[i]/errdiffBRBR[i];
     cout << diffBRBR[i] <<  " +- " << errdiffBRBR[i] <<  "   max diff   "<< sigmas << " sigmas"<< endl;
   }

   f0 = new TF1("f0", "gaus", 0., .1);
   f0->SetLineColor(kRed);
   pulls.Fit(f0,"l","pe");
   double mean = f0->GetParameter(1);
   double sigma = f0->GetParameter(2);
   double errmean = f0->GetParError(1);
   double errsigma = f0->GetParError(2);
   double chisqua = f0->GetChisquare();
   double ndof = f0->GetNDF();
   
   TLatex tl;
   double x = 0.15;
   char myline[100];
   sprintf(myline, "mean = %6.5f +/- %5.5f", mean, errmean); 
   tl.DrawTextNDC(x, 0.8, myline);
   sprintf(myline, "sigma = %6.5f +/- %5.5f", sigma, errsigma); 
   tl.DrawTextNDC(x, 0.7, myline);
   //   sprintf(myline, "chi2 = %7.2f ", chisqua/ndof);
   tl.DrawTextNDC(x, 0.6, myline);
   c0.SaveAs("pulls.eps");
   brbr.Fit(f0,"l","pe");
   mean = f0->GetParameter(1);
   sigma = f0->GetParameter(2);
   errmean = f0->GetParError(1);
   errsigma = f0->GetParError(2);
   chisqua = f0->GetChisquare();
   ndof = f0->GetNDF();
   
   sprintf(myline, "mean = %6.5f +/- %5.5f", mean, errmean); 
   tl.DrawTextNDC(x, 0.8, myline);
   sprintf(myline, "sigma = %6.5f +/- %5.5f", sigma, errsigma); 
   tl.DrawTextNDC(x, 0.7, myline);
   //   sprintf(myline, "chi2 = %7.2f ", chisqua/ndof);
   //   tl.DrawTextNDC(x, 0.6, myline);
   c0.SaveAs("BRBR.eps");
   
   //titles

   if(!strcmp(thename,"mnuSqHigh")){
     sprintf(axistit,"Mv^2 upper cut(GeV^2)");
     sprintf(title,"mm2 scan");
   }
   if(!strcmp(thename,"leptonPCut")){
     sprintf(axistit,"p*(GeV)");
     sprintf(title,"pcms scan");
   }
   if(!strcmp(thename,"prmm2cut")){
				     cout << "here"<< endl;
     sprintf(axistit,"prmm2(GeV^2)");
     sprintf(title,"pr mm2 scan");
   }
   if(!strcmp(thename,"mxCut")){
     sprintf(axistit,"Mx(GeV)");
     sprintf(title,"Mx scan");
   }
   if(!strcmp(thename,"q2Cut")){
     sprintf(axistit,"Q^2(GeV^2)");
     sprintf(title,"Q^2 scan");
   }
   if(!strcmp(thename,"chHigh")){
     sprintf(axistit,"Q cut (|Q|<Qcut)");
     sprintf(title,"Qtot scan");
   }
   if(!strcmp(thename,"maxintpur")){
     sprintf(axistit,"max integrated purity");
     sprintf(title,"int purity scan");
   }
   if(!strcmp(thename,"minintpur")){ 
     sprintf(axistit,"min integrated purity"); 
     sprintf(title,"int purity scan"); 
   } 
   sprintf(temptit,"%s%s",title," optimization");
   
   c0.Clear();
   sprintf (name,"%s%s%s",namein,thename,".eps");
   TGraphErrors *results;
   results = new TGraphErrors(counter,thescan,BRBR,errthescan,errBRBR);
   results->SetMinimum(0.);
   results->SetMaximum(refvalue*2);
   results->SetTitle(title);
   results->SetMarkerSize(1.5);
   results->Draw("AP");
   TLine line(min,refvalue,max,refvalue);
   line.SetLineColor(kRed); line.SetLineWidth(2); line.Draw();
   TLine lineupstat(min,refvalue*1.128,max,refvalue*1.128);
   lineupstat.SetLineColor(kRed); lineupstat.SetLineStyle(2);lineupstat.SetLineWidth(2); lineupstat.Draw();
   TLine linedownstat(min,refvalue*.872,max,refvalue*.872);
   linedownstat.SetLineColor(kRed); linedownstat.SetLineStyle(2);linedownstat.SetLineWidth(2); linedownstat.Draw();
   TLine lineuptot(min,refvalue*1.225,max,refvalue*1.225);
   lineuptot.SetLineColor(kGreen); lineuptot.SetLineStyle(3);lineuptot.SetLineWidth(2); lineuptot.Draw();
   TLine linedowntot(min,refvalue*.795,max,refvalue*.795);
   linedowntot.SetLineColor(kGreen); linedowntot.SetLineStyle(3);linedowntot.SetLineWidth(2); linedowntot.Draw();
   TLegendEntry *legge; 
   leg = new TLegend(0.6,0.7,0.88,0.89);
   c0->Update();
   c0->GetFrame()->SetFillColor(0);
   c0->GetFrame()->SetBorderSize(12);
   results->GetHistogram()->SetTitleOffset(2., "Y");
   results->GetHistogram()->SetXTitle(axistit);
   results->GetHistogram()->SetYTitle("ratio(BR)");
   c0->Modified();
   TLegendEntry *legge; 
   leg = new TLegend(0.6,0.7,0.88,0.89);
   leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.1); 
   leg->SetFillColor(0); 
   legge = leg->AddEntry(results, type, "n"); 
   leg->Draw();
   c0.SaveAs(name);
   c0.Clear();
   sprintf (name,"%s%s%s",namein,thename,"diff.eps");
   TGraphErrors *diffresults;
   diffresults = new TGraphErrors(counter,thescan,diffBRBR,errthescan,errdiffBRBR);
   diffresults->SetMinimum(-.03);
   diffresults->SetMaximum(.03);
   diffresults->SetTitle(title);
   diffresults->SetMarkerSize(1.5);
   diffresults->Draw("AP");
   TLine line2(min,0.,max,0.);
   line2.SetLineColor(kRed); line2.SetLineWidth(2); line2.Draw();
   c0->Update();
   c0->GetFrame()->SetFillColor(0);
   c0->GetFrame()->SetBorderSize(12);
   diffresults->GetHistogram()->SetTitleOffset(2., "Y");
   diffresults->GetHistogram()->SetXTitle(axistit);
   diffresults->GetHistogram()->SetYTitle("D(ratio(BR))");
   c0->Modified();
   TLegendEntry *legge; 
   leg = new TLegend(0.6,0.7,0.88,0.89);
   leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.1); 
   leg->SetFillColor(0); 
   legge = leg->AddEntry(results, type, "n"); 
   leg->Draw();
   c0.SaveAs(name);
   TGraphErrors *optresults;
   c0.Clear();
   
   sprintf (name,"%s%s%s",namein,thename,"opt.eps");
   optresults = new TGraphErrors(counter,thescan,optBRBR,errthescan,erroptBRBR);
   optresults->SetTitle(temptit);
   optresults->Draw("AP");
   //   optresults->Fit("pol4");
   c0->Update();
   c0->GetFrame()->SetFillColor(0);
   c0->GetFrame()->SetBorderSize(30);   
   optresults->GetHistogram()->SetTitleOffset(2., "Y");
   optresults->GetHistogram()->SetXTitle(axistit);
   optresults->GetHistogram()->SetYTitle("relative ratio(BR) error");
   c0->Modified();
   c0.SaveAs(name);
   
   c0.Clear();
   sprintf (name,"%s%s%s",namein,thename,"vubcomp.eps");
   optresults = new TGraphErrors(counter,thescan,vubcomp,errthescan,errvubcomp);
   optresults->Draw("AP");
   c0.SaveAs(name);
   
   c0.Clear();
   sprintf (name,"%s%s%s",namein,thename,"vcbcomp.eps");
   optresults = new TGraphErrors(counter,thescan,vcbcomp,errthescan,errvcbcomp);
   optresults->Draw("AP");
   c0.SaveAs(name);
   
   c0.Clear();
   sprintf (name,"%s%s%s",namein,thename,"othcomp.eps");
   optresults = new TGraphErrors(counter,thescan,othcomp,errthescan,errothcomp);
   optresults->Draw("AP");
   c0.SaveAs(name);

   c0.Clear();
   sprintf (name,"%s%s%s",namein,thename,"vub1bin.eps");
   optresults = new TGraphErrors(counter,thescan,vub1bin,errthescan,errvub1bin);
   optresults->SetMaximum(200.);
   optresults->SetMinimum(0.);
   optresults->Draw("AP");
   //   optresults->Fit("pol1"); 
   c0.SaveAs(name);

   c0.Clear();
   sprintf (name,"%s%s%s",namein,thename,"vcb1bin.eps");
   optresults = new TGraphErrors(counter,thescan,vcb1bin,errthescan,errvcb1bin);
   optresults->Draw("AP");
   optresults->Fit("pol3"); 
   c0.SaveAs(name);

   c0.Clear();
   sprintf (name,"%s%s%s",namein,thename,"oth1bin.eps");
   optresults = new TGraphErrors(counter,thescan,oth1bin,errthescan,erroth1bin);
   optresults->Draw("AP");
   c0.SaveAs(name);

   c0.Clear();
   sprintf (name,"%s%s%s",namein,thename,"epsu.eps");
   optresults = new TGraphErrors(counter,thescan,epsu,errthescan,errepsuMC);
   optresults->Draw("AP");
   c0.SaveAs(name);

   c0.Clear();
   sprintf (name,"%s%s%s",namein,thename,"epstot.eps");
   optresults = new TGraphErrors(counter,thescan,epstot,errthescan,errepstotMC);
   optresults->Draw("AP");
   c0.SaveAs(name);

   c0.Clear();
   sprintf (name,"%s%s%s",namein,thename,"epsmx.eps");
   optresults = new TGraphErrors(counter,thescan,epsmx,errthescan,errepsmxMC);
   optresults->Draw("AP");
   c0.SaveAs(name);

   c0.Clear();
   sprintf (name,"%s%s%s",namein,thename,"nsl.eps");
   optresults = new TGraphErrors(counter,thescan,nsl,errthescan,errthescan);
   optresults->Draw("AP");
   c0.SaveAs(name);

   c0.Clear();
   sprintf (name,"%s%s%s",namein,thename,"nslmc.eps");
   optresults = new TGraphErrors(counter,thescan,nslmc,errthescan,errthescan);
   optresults->Draw("AP");
   c0.SaveAs(name);

   c0.Clear();
   sprintf (name,"%s%s%s",namein,thename,"fact.eps");
   optresults = new TGraphErrors(counter,thescan,fact,errthescan,errthescan);
   optresults->Draw("AP");
   c0.SaveAs(name);

   c0.Clear();
   sprintf (name,"%s%s%s",namein,thename,"pstarfact.eps");
   optresults = new TGraphErrors(counter,thescan,pstarfact,errthescan,errthescan);
   optresults->Draw("AP");
   c0.SaveAs(name);

   c0.Clear();
   sprintf (name,"%s%s%s",namein,thename,"chisq.eps");
   optresults = new TGraphErrors(counter,thescan,chisq,errthescan,errthescan);
   optresults->Draw("AP");
   c0.SaveAs(name);

   
   c1 = new TCanvas("c1"," ",200,10,1300,520); 
   sprintf (name,"%s%s%s",namein,thename,"threeplots.eps");
   c1.Clear();
   c1.Divide(3,1);
   c1.cd(1);
   optresults = new TGraphErrors(counter,thescan,optBRBR,errthescan,erroptBRBR);
   optresults->SetTitle(temptit);
   optresults->SetMarkerColor(kRed);
   optresults->Draw("AP");
   c1->Update();
   c1->GetFrame()->SetFillColor(0);
   c1->GetFrame()->SetBorderSize(10);
   optresults->SetMinimum(0.05);
   //   optresults->SetMaximum(0.15);
   optresults->GetHistogram()->SetTitleOffset(2., "Y");
   optresults->GetHistogram()->SetXTitle(axistit);
   optresults->GetHistogram()->SetYTitle("relative ratio(BR) error");
   optresults->GetHistogram()->SetTitleSize(.07);        
   optresults->GetHistogram()->SetTitleSize(.07, "x");        
   optresults->GetHistogram()->SetTitleSize(.07, "y");
   c1->Modified();
   c1.cd(2);
   optresults = new TGraphErrors(counter,thescan,vub1bin,errthescan,errvub1bin);
   optresults->SetTitle("");
   optresults->SetMarkerColor(kRed);
   optresults->Draw("AP");
   c1->Update();
   c1->GetFrame()->SetFillColor(0);
   c1->GetFrame()->SetBorderSize(10);
   optresults->SetMinimum(0.0);
   optresults->GetHistogram()->SetTitleOffset(2., "Y");
   optresults->GetHistogram()->SetXTitle(axistit);
   optresults->GetHistogram()->SetYTitle("signal events");
   optresults->GetHistogram()->SetTitleSize(.07, "x");        
   optresults->GetHistogram()->SetTitleSize(.07, "y");
   c1->Modified();
   c1.cd(3);
   optresults = new TGraphErrors(counter,thescan,vcb1bin,errthescan,errvcb1bin);
   optresults->SetTitle("");
   optresults->SetMarkerColor(kRed);
   optresults->Draw("AP");
   c1->Update();
   c1->GetFrame()->SetFillColor(0);
   c1->GetFrame()->SetBorderSize(10);
   optresults->SetMinimum(0.0);
   optresults->GetHistogram()->SetTitleOffset(2., "Y");
   optresults->GetHistogram()->SetXTitle(axistit);
   optresults->GetHistogram()->SetYTitle("background events");
   optresults->GetHistogram()->SetTitleSize(.07, "x");        
   optresults->GetHistogram()->SetTitleSize(.07, "y");
   c1->Modified();
   c1.SaveAs(name);
   cout << type << endl;

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
