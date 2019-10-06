void fitMx(){

  // Load RooFit libraries

  gSystem->Load("../../RooFitCore/tmp/libRooFitCore.so");
  gSystem->Load("../../RooFitModels/tmp/libRooFitModels.so");

  using namespace RooFit;

  // Load Data

  TFile *datafile = new TFile("RooRootFiles/signal.root");
  TTree *tree = (TTree*)datafile->Get("ntp2");

  // List of variables saved in the dataset

  // Breco candidate
  RooRealVar* mes = new RooRealVar("mes","mes",5.2,5.3,"GeV/c^{2}");
  RooRealVar* ass_deltapB = new RooRealVar("ass_deltapB","ass_deltapB",-999,999,"GeV/c");
  RooRealVar* brecocharge = new RooRealVar("brecocharge","Breco charge",-1,1,"");
  RooRealVar* brecoflav = new RooRealVar("brecoflav","Breco flavor (norm)",-1,1,"");

  // Best Lepton
  RooRealVar* nle = new RooRealVar("nle","number of leptons",-999,999,"");
  RooRealVar* lcharge = new RooRealVar("lcharge","lepton charge",-999,999,"");
  RooRealVar* isele = new RooRealVar("isele","lepton type",-999,999,""); 
  RooRealVar* pcms = new RooRealVar("pcms","momentum CMS",-999,999,"GeV/c");
  RooRealVar* plab = new RooRealVar("plab","momentum LAB",-999,999,"GeV/c");
  RooRealVar* tlab = new RooRealVar("tlab","theta LAB",-999,999,"");

  // MC truth
  RooRealVar* Gvxbtyp = new RooRealVar("Gvxbtyp","MC Truth",-999,999,"");

  // ETA L NU
  RooRealVar* nrecoEta = new RooRealVar("nrecoEta","Number of #eta",0,100,"");
  RooRealVar* indexbestEta = new RooRealVar("indexbestEta","Index of best #eta",-999,999,"");
  RooRealVar* modebestEta = new RooRealVar("modebestEta","Decay mode of best #eta",-999,999,"");
  RooRealVar* barembestEta = new RooRealVar("barembestEta","M(#eta)",-999,999,"GeV/c^{2}");
  
  // ETA PRIME L NU
  RooRealVar* nrecoEtap = new RooRealVar("nrecoEtap","Number of #eta'",0,100,"");
  RooRealVar* indexbestEtap = new RooRealVar("indexbestEtap","Index of best #eta'",-999,999,"");
  RooRealVar* modebestEtap = new RooRealVar("modebestEtap","Decay mode of best #eta'",-999,999,"");
  RooRealVar* barembestEtap = new RooRealVar("barembestEtap","M(#eta')",-999,999,"GeV/c^{2}");
  RooRealVar* Rho0massdaubestEtap = new RooRealVar("Rho0massdaubestEtap","M(#rho0)",-999,999,"GeV/c^{2}");
  RooRealVar* EtamassdaubestEtap = new RooRealVar("EtamassdaubestEtap","M(#eta)",-999,999,"GeV/c^{2}");

  // Add variables to the dataset
  
  RooArgSet myVars;

  // Breco  
  myVars.add(*mes);
  myVars.add(*ass_deltapB);
  myVars.add(*brecocharge);
  myVars.add(*brecoflav);
  // Best lepton
  myVars.add(*nle);
  myVars.add(*lcharge);
  myVars.add(*isele); 
  myVars.add(*pcms);
  myVars.add(*plab);
  myVars.add(*tlab);
  // MC truth
  myVars.add(*Gvxbtyp);
  // ETA L NU
  myVars.add(*nrecoEta);
  myVars.add(*indexbestEta);
  myVars.add(*modebestEta);
  myVars.add(*barembestEta);
  // ETA PRIME L NU
  myVars.add(*nrecoEtap);
  myVars.add(*indexbestEtap);
  myVars.add(*modebestEtap);
  myVars.add(*barembestEtap);
  myVars.add(*Rho0massdaubestEtap);
  myVars.add(*EtamassdaubestEtap);
  
  // Create ROOT file containing the dataset

  RooDataSet* mydata = new RooDataSet("dataset","dataset",tree,myVars);
  
  TFile f("RooDataSets/dataset.root","RECREATE");
  mydata->Write();
  f.Close();
  
  // Create variables, parameters and pdf's for fit...

  //---------------------------------------------//
  //                 ETA L NU                    //
  //---------------------------------------------//

  RooRealVar EtaMass("barembestEta","M(#eta)",0.45,0.65,"GeV/c^{2}");
  RooRealVar EtaMean("EtaMean","Mean",0.5478,0.546,0.549);
  RooRealVar EtaSigma("EtaSigma","Sigma",0.005,0.001,0.10);
  RooRealVar BkgEtaSlope("BkgEtaSlope","Bkg Slope",-100,100);
  RooRealVar NEtaSig("NEtaSig","NSig",0.,1000000.);
  RooRealVar NEtaBkg("NEtaBkg","NBkg",0.,1000000.);

  RooGaussian EtaMassCore("EtaMassCore","EtaMassCore",EtaMass,EtaMean,EtaSigma);
  RooPolynomial EtaBkg("EtaBkg","EtaBkg",EtaMass,BkgEtaSlope);
  
  RooAddPdf EtaMassPdf("EtaMassPdf","EtaMassPdf",RooArgList(EtaMassCore,EtaBkg),                             RooArgList(NEtaSig,NEtaBkg));
  
  // Reduce the dataset

  TFile myetadataset("RooDataSets/dataset.root");
  RooDataSet *myetadata = (RooDataSet*)myetadataset.Get("dataset");
  RooDataSet *etaslice = myetadata->reduce("(abs(brecocharge)==1)&&!(abs(brecocharge)!=0&&(lcharge + brecoflav)!=0)&&(nle<=1)&&(tlab<2.37)&&(tlab>0.36)&&(plab>0.5)&&((isele==1&&pcms>0.5)||(isele==0&&pcms>0.8))&&Gvxbtyp==-12&&mes>=5.2&&mes<=5.3&&nrecoEta>0&&(ass_deltapB<0.2)&&(indexbestEta>-999)");
  RooDataSet *etaslicegg = etaslice->reduce("modebestEta==1");
  RooDataSet *etaslicepipipi0 = etaslice->reduce("modebestEta==2");
  RooDataSet *etaslicepi0pi0pi0 = etaslice->reduce("modebestEta==3");
  
  // Fit data and plot the results

  RooFitResult* rEta = EtaMassPdf.fitTo(*etaslice,"qhevr");
  rEta->Print("v");

  RooPlot * frVarEta = EtaMass.frame(40);
  etaslice->plotOn(frVarEta);
  EtaMassPdf.plotOn(frVarEta);
  EtaMean.setPlotLabel("mean");
  EtaSigma.setPlotLabel("sigma");
  NEtaSig.setPlotLabel("S");
  NEtaBkg.setPlotLabel("B");
  BkgEtaSlope.setConstant();
  EtaMassPdf.paramOn(frVarEta,Format("NEU"),Layout(0.65,0.9,0.9));
  BkgEtaSlope.setConstant(kFALSE);
  frVarEta->Draw();  
  frVarEta->SetTitle("");
  TPaveText *tpteta = frVarEta->findObject("TPave");
  tpteta->SetTextSize(0.030);
  c1->SaveAs("RooPlots/EtaMass.eps");

  // eta -> gamma gamma

  RooFitResult* rEtagg = EtaMassPdf.fitTo(*etaslicegg,"qhevr");
  rEtagg->Print("v");

  RooPlot * frVarEtagg = EtaMass.frame(40);
  etaslicegg->plotOn(frVarEtagg);
  EtaMassPdf.plotOn(frVarEtagg);
  BkgEtaSlope.setConstant();
  EtaMassPdf.paramOn(frVarEtagg,Format("NELU"),Layout(0.65,0.9,0.9));
  BkgEtaSlope.setConstant(kFALSE);
  frVarEtagg->Draw();
  frVarEtagg->SetTitle("#eta #rightarrow #gamma #gamma");
  TPaveText *tptgg = frVarEtagg->findObject("TPave");
  tptgg->SetTextSize(0.030);
  c1->SaveAs("RooPlots/Eta_ggMass.eps");

  // eta -> pi pi pi0

  RooFitResult* rEtapipipi0 = EtaMassPdf.fitTo(*etaslicepipipi0,"qhevr");
  rEtagg->Print("v");

  RooPlot * frVarEtapipipi0 = EtaMass.frame(40);
  etaslicepipipi0->plotOn(frVarEtapipipi0);
  EtaMassPdf.plotOn(frVarEtapipipi0);
  BkgEtaSlope.setConstant();
  EtaMassPdf.paramOn(frVarEtapipipi0,Format("NELU"),Layout(0.65,0.9,0.9));
  BkgEtaSlope.setConstant(kFALSE);
  frVarEtapipipi0->Draw();
  frVarEtapipipi0->SetTitle("#eta #rightarrow #pi^{+} #pi^{-} #pi^{0}");
  TPaveText *tptpipipi0 = frVarEtapipipi0->findObject("TPave");
  tptpipipi0->SetTextSize(0.030);
  c1->SaveAs("RooPlots/Eta_pipipi0Mass.eps");

  // eta -> pi0 pi0 pi0

  RooFitResult* rEtapi0pi0pi0 = EtaMassPdf.fitTo(*etaslicepi0pi0pi0,"qhevr");
  rEtapi0pi0pi0->Print("v");

  RooPlot * frVarEtapi0pi0pi0 = EtaMass.frame(40);
  etaslicepi0pi0pi0->plotOn(frVarEtapi0pi0pi0);
  EtaMassPdf.plotOn(frVarEtapi0pi0pi0);
  BkgEtaSlope.setConstant();
  EtaMassPdf.paramOn(frVarEtapi0pi0pi0,Format("NELU"),Layout(0.65,0.9,0.9));
  BkgEtaSlope.setConstant(kFALSE);
  frVarEtapi0pi0pi0->Draw();
  frVarEtapi0pi0pi0->SetTitle("#eta #rightarrow #pi^{0} #pi^{0} #pi^{0}");
  TPaveText *tptpi0pi0pi0 = frVarEtapi0pi0pi0->findObject("TPave");
  tptpi0pi0pi0->SetTextSize(0.030);
  c1->SaveAs("RooPlots/Eta_pi0pi0pi0Mass.eps");


  //---------------------------------------------//
  //               ETA PRIME L NU                //
  //---------------------------------------------//

  // eta' mass

  RooRealVar EtapMass("barembestEtap","M(#eta')",0.86,1.06,"GeV/c^{2}");
  RooRealVar EtapMean("EtapMean","Mean",0.9578,0.956,0.959);
  RooRealVar EtapSigma("EtapSigma","Sigma",0.005,0.001,0.10);
  RooRealVar BkgEtapSlope("BkgEtapSlope","Bkg Slope",-100,100);
  RooRealVar NEtapSig("NEtapSig","NSig",0.,1000000.);
  RooRealVar NEtapBkg("NEtapBkg","NBkg",0.,1000000.);

  RooGaussian EtapMassCore("EtapMassCore","EtapMassCore",EtapMass,EtapMean,EtapSigma);
  RooPolynomial EtapBkg("EtapBkg","EtapBkg",EtapMass,BkgEtapSlope);
  RooAddPdf EtapMassPdf("EtapMassPdf","EtapMassPdf",RooArgList(EtapMassCore,EtapBkg),RooArgList(NEtapSig,NEtapBkg));

  // rho0 mass

  RooRealVar Rho0Mass("Rho0massdaubestEtap","M(#rho^{0})",0.2,1.0,"GeV/c^{2}");
  RooRealVar Rho0Mean("Rho0Mean","Mean",0.7758,0.6,1.0);
  RooRealVar Rho0Sigma("Rho0Sigma","Sigma",0.05,0.001,0.5);
  RooRealVar CoreRho0Frac("CoreRho0Frac","Core Fraction",0.8,0.,1.);
  RooRealVar BkgRho0Slope("BkgRho0Slope","Bkg Slope",-100,100);
  RooRealVar NRho0Sig("NRho0Sig","NSig",0.,1000000.);
  RooRealVar NRho0Bkg("NRho0Bkg","NBkg",0.,1000000.);

  RooGaussian Rho0MassCore("Rho0MassCore","Rho0MassCore",Rho0Mass,Rho0Mean,Rho0Sigma);  
  RooPolynomial Rho0Bkg("Rho0Bkg","Rho0Bkg",Rho0Mass,BkgRho0Slope);
  RooAddPdf Rho0MassPdf("Rho0MassPdf","Rho0MassPdf",RooArgList(Rho0MassCore,Rho0Bkg),RooArgList(NRho0Sig,NRho0Bkg));

  // eta mass

  RooRealVar EtadauMass("EtamassdaubestEtap","M(#eta)",0.45,0.65,"GeV/c^{2}");
  RooRealVar EtadauMean("EtadauMean","Mean",0.5478,0.546,0.549);
  RooRealVar EtadauSigma("EtadauSigma","Sigma",0.005,0.001,0.10);
  RooRealVar CoreEtadauFrac("CoreEtadauFrac","Core Fraction",0.8,0.,1.);
  RooRealVar BkgEtadauSlope("BkgEtadauSlope","Bkg Slope",-100,100);
  RooRealVar NEtadauSig("NEtadauSig","NSig",0.,1000000.);
  RooRealVar NEtadauBkg("NEtadauBkg","NBkg",0.,1000000.);

  RooGaussian EtadauMassCore("EtadauMassCore","EtadauMassCore",EtadauMass,EtadauMean,EtadauSigma);
  RooPolynomial EtadauBkg("EtadauBkg","EtadauBkg",EtadauMass,BkgEtadauSlope);
  RooAddPdf EtadauMassPdf("EtadauMassPdf","EtadauMassPdf",RooArgList(EtadauMassCore,EtadauBkg),RooArgList(NEtadauSig,NEtadauBkg));


  // Reduce the dataset

  TFile myetapdataset("RooDataSets/dataset.root");
  RooDataSet *myetapdata = (RooDataSet*)myetapdataset.Get("dataset");
  RooDataSet *etapslice = myetapdata->reduce("(abs(brecocharge)==1)&&!(abs(brecocharge)!=0&&(lcharge + brecoflav)!=0)&&(nle<=1)&&(tlab<2.37)&&(tlab>0.36)&&(plab>0.5)&&(pcms>1.)&&Gvxbtyp==-15&&mes>=5.2&&mes<=5.3&&nrecoEtap>0&&(ass_deltapB<0.2)&&(indexbestEtap>-999)");
  RooDataSet *etapslicerho0g = etapslice->reduce("modebestEtap==1");
  RooDataSet *etapsliceetagg = etapslice->reduce("modebestEtap==2");
  RooDataSet *etapsliceetapipipi0 = etapslice->reduce("modebestEtap==3");
  RooDataSet *etapsliceetapi0pi0pi0 = etapslice->reduce("modebestEtap==4");
  RooDataSet *etapsliceetaall = etapslice->reduce("(modebestEtap==2)||(modebestEtap==3)||(modebestEtap==4)"); 

  // Fit data and plot the results

  RooFitResult* rEtap = EtapMassPdf.fitTo(*etapslice,"qhevr");
  rEtap->Print("v");

  RooPlot * frVarEtap = EtapMass.frame(40);
  etapslice->plotOn(frVarEtap);
  EtapMassPdf.plotOn(frVarEtap);
  EtapMean.setPlotLabel("mean");
  EtapSigma.setPlotLabel("sigma");
  NEtapSig.setPlotLabel("S");
  NEtapBkg.setPlotLabel("B");
  BkgEtapSlope.setConstant();
  EtapMassPdf.paramOn(frVarEtap,Format("NELU"),Layout(0.65,0.9,0.9));
  BkgEtapSlope.setConstant(kFALSE);
  frVarEtap->Draw();  
  frVarEtap->SetTitle("");
  TPaveText *tptetap = frVarEtap->findObject("TPave");
  tptetap->SetTextSize(0.030);
  c1->SaveAs("RooPlots/EtapMass.eps");

  // eta' -> rho0 gamma

  RooFitResult* rEtaprho0g = EtapMassPdf.fitTo(*etapslicerho0g,"qhevr");
  rEtaprho0g->Print("v");

  RooPlot * frVarEtaprho0g = EtapMass.frame(40);
  etapslicerho0g->plotOn(frVarEtaprho0g);
  EtapMassPdf.plotOn(frVarEtaprho0g);
  EtapMean.setPlotLabel("mean");
  EtapSigma.setPlotLabel("sigma");
  BkgEtapSlope.setConstant();
  EtapMassPdf.paramOn(frVarEtaprho0g,Format("NELU"),Layout(0.65,0.9,0.9));
  BkgEtapSlope.setConstant(kFALSE);
  frVarEtaprho0g->Draw();  
  frVarEtaprho0g->SetTitle("#eta' #rightarrow #rho^{0} #gamma");
  TPaveText *tptetaprho0g = frVarEtaprho0g->findObject("TPave");
  tptetaprho0g->SetTextSize(0.030);
  c1->SaveAs("RooPlots/Etap_rho0gMass.eps");

  RooFitResult* rRho0 = Rho0MassPdf.fitTo(*etapslicerho0g,"qhevr");
  rRho0->Print("v");

  RooPlot * frVarRho0 = Rho0Mass.frame(40);
  etapslicerho0g->plotOn(frVarRho0);
  Rho0MassPdf.plotOn(frVarRho0);
  Rho0Mean.setPlotLabel("mean");
  Rho0Sigma.setPlotLabel("sigma");
  NRho0Sig.setPlotLabel("S");
  NRho0Bkg.setPlotLabel("B");
  BkgRho0Slope.setConstant();
  Rho0MassPdf.paramOn(frVarRho0,Format("NELU"),Layout(0.1,0.4,0.9));
  BkgRho0Slope.setConstant(kFALSE);
  frVarRho0->Draw();  
  frVarRho0->SetTitle("#eta' #rightarrow #rho^{0} #gamma");
  TPaveText *tptrho0 = frVarRho0->findObject("TPave");
  tptrho0->SetTextSize(0.030);
  c1->SaveAs("RooPlots/Rho0Mass.eps");

  // eta' -> eta pi pi

  RooFitResult* rEtadau = EtadauMassPdf.fitTo(*etapsliceetaall,"qhevr");
  rEtadau->Print("v");

  RooPlot * frVarEtadau = EtadauMass.frame(40);
  etapsliceetaall->plotOn(frVarEtadau);
  EtadauMassPdf.plotOn(frVarEtadau);
  EtadauMean.setPlotLabel("mean");
  EtadauSigma.setPlotLabel("sigma");
  NEtadauSig.setPlotLabel("S");
  NEtadauBkg.setPlotLabel("B");
  BkgEtadauSlope.setConstant();
  EtadauMassPdf.paramOn(frVarEtadau,Format("NELU"),Layout(0.65,0.9,0.9));
  BkgEtadauSlope.setConstant(kFALSE);
  frVarEtadau->Draw();
  frVarEtadau->SetTitle("#eta' #rightarrow #eta #pi^{+} #pi^{-}");
  TPaveText *tptetadau = frVarEtadau->findObject("TPave");
  tptetadau->SetTextSize(0.030);
  c1->SaveAs("RooPlots/EtadauMass.eps");

  // eta' -> eta pi pi, eta -> gamma gamma

  RooFitResult* rEtapetagg = EtapMassPdf.fitTo(*etapsliceetagg,"qhevr");
  rEtapetagg->Print("v");

  RooPlot * frVarEtapetagg = EtapMass.frame(40);
  etapsliceetagg->plotOn(frVarEtapetagg);
  EtapMassPdf.plotOn(frVarEtapetagg);
  EtapMean.setPlotLabel("mean");
  EtapSigma.setPlotLabel("sigma");
  BkgEtapSlope.setConstant();
  EtapMassPdf.paramOn(frVarEtapetagg,Format("NELU"),Layout(0.65,0.9,0.9));
  BkgEtapSlope.setConstant(kFALSE);
  frVarEtapetagg->Draw();  
  frVarEtapetagg->SetTitle("#eta' #rightarrow #eta #pi^{+} #pi^{-}, #eta #rightarrow #gamma #gamma");
  TPaveText *tptetapetagg = frVarEtapetagg->findObject("TPave");
  tptetapetagg->SetTextSize(0.030);
  c1->SaveAs("RooPlots/Etap_etapipi_ggMass.eps");

  RooFitResult* rEtadaugg = EtadauMassPdf.fitTo(*etapsliceetagg,"qhevr");
  rEtadaugg->Print("v");

  RooPlot * frVarEtadaugg = EtadauMass.frame(40);
  etapsliceetagg->plotOn(frVarEtadaugg);
  EtadauMassPdf.plotOn(frVarEtadaugg);
  EtadauMean.setPlotLabel("mean");
  EtadauSigma.setPlotLabel("sigma");
  NEtadauSig.setPlotLabel("S");
  NEtadauBkg.setPlotLabel("B");
  BkgEtadauSlope.setConstant();
  EtadauMassPdf.paramOn(frVarEtadaugg,Format("NELU"),Layout(0.65,0.9,0.9));
  BkgEtadauSlope.setConstant(kFALSE);
  frVarEtadaugg->Draw();
  frVarEtadaugg->SetTitle("#eta' #rightarrow #eta #pi^{+} #pi^{-}, #eta #rightarrow #gamma #gamma");
  TPaveText *tptetadaugg = frVarEtadaugg->findObject("TPave");
  tptetadaugg->SetTextSize(0.030);
  c1->SaveAs("RooPlots/Etadau_ggMass.eps");

  // eta' -> eta pi pi, eta -> pi pi pi0

  RooFitResult* rEtapetapipipi0 = EtapMassPdf.fitTo(*etapsliceetapipipi0,"qhevr");
  rEtapetapipipi0->Print("v");

  RooPlot * frVarEtapetapipipi0 = EtapMass.frame(40);
  etapsliceetapipipi0->plotOn(frVarEtapetapipipi0);
  EtapMassPdf.plotOn(frVarEtapetapipipi0);
  EtapMean.setPlotLabel("mean");
  EtapSigma.setPlotLabel("sigma");
  BkgEtapSlope.setConstant();
  EtapMassPdf.paramOn(frVarEtapetapipipi0,Format("NELU"),Layout(0.65,0.9,0.9));
  BkgEtapSlope.setConstant(kFALSE);
  frVarEtapetapipipi0->Draw();  
  frVarEtapetapipipi0->SetTitle("#eta' #rightarrow #eta #pi^{+} #pi^{-}, #eta #rightarrow #pi^{+} #pi^{-} #pi^{0}");
  TPaveText *tptetapetapipipi0 = frVarEtapetapipipi0->findObject("TPave");
  tptetapetapipipi0->SetTextSize(0.030);
  c1->SaveAs("RooPlots/Etap_etapipi_pipipi0Mass.eps");

  RooFitResult* rEtadaupipipi0 = EtadauMassPdf.fitTo(*etapsliceetapipipi0,"qhevr");
  rEtadaupipipi0->Print("v");

  RooPlot * frVarEtadaupipipi0 = EtadauMass.frame(40);
  etapsliceetapipipi0->plotOn(frVarEtadaupipipi0);
  EtadauMassPdf.plotOn(frVarEtadaupipipi0);
  EtadauMean.setPlotLabel("mean");
  EtadauSigma.setPlotLabel("sigma");
  NEtadauSig.setPlotLabel("S");
  NEtadauBkg.setPlotLabel("B");
  BkgEtadauSlope.setConstant();
  EtadauMassPdf.paramOn(frVarEtadaupipipi0,Format("NELU"),Layout(0.65,0.9,0.9));
  BkgEtadauSlope.setConstant(kFALSE);
  frVarEtadaupipipi0->Draw();
  frVarEtadaupipipi0->SetTitle("#eta' #rightarrow #eta #pi^{+} #pi^{-}, #eta #rightarrow #pi^{+} #pi^{-} #pi^{0}");
  TPaveText *tptetadaupipipi0 = frVarEtadaupipipi0->findObject("TPave");
  tptetadaupipipi0->SetTextSize(0.030);
  c1->SaveAs("RooPlots/Etadau_pipipi0Mass.eps");

  // eta' -> eta pi pi, eta -> pi0 pi0 pi0

  RooFitResult* rEtapetapi0pi0pi0 = EtapMassPdf.fitTo(*etapsliceetapi0pi0pi0,"qhevr");
  rEtapetapi0pi0pi0->Print("v");

  RooPlot * frVarEtapetapi0pi0pi0 = EtapMass.frame(40);
  etapsliceetapi0pi0pi0->plotOn(frVarEtapetapi0pi0pi0);
  EtapMassPdf.plotOn(frVarEtapetapi0pi0pi0);
  EtapMean.setPlotLabel("mean");
  EtapSigma.setPlotLabel("sigma");
  BkgEtapSlope.setConstant();
  EtapMassPdf.paramOn(frVarEtapetapi0pi0pi0,Format("NELU"),Layout(0.65,0.9,0.9));
  BkgEtapSlope.setConstant(kFALSE);
  frVarEtapetapi0pi0pi0->Draw();  
  frVarEtapetapi0pi0pi0->SetTitle("#eta' #rightarrow #eta #pi^{+} #pi^{-}, #eta #rightarrow #pi^{0} #pi^{0} #pi^{0}");
  TPaveText *tptetapetapi0pi0pi0 = frVarEtapetapi0pi0pi0->findObject("TPave");
  tptetapetapi0pi0pi0->SetTextSize(0.030);
  c1->SaveAs("RooPlots/Etap_etapipi_pi0pi0pi0Mass.eps");

  RooFitResult* rEtadaupi0pi0pi0 = EtadauMassPdf.fitTo(*etapsliceetapi0pi0pi0,"qhevr");
  rEtadaupi0pi0pi0->Print("v");

  RooPlot * frVarEtadaupi0pi0pi0 = EtadauMass.frame(40);
  etapsliceetapi0pi0pi0->plotOn(frVarEtadaupi0pi0pi0);
  EtadauMassPdf.plotOn(frVarEtadaupi0pi0pi0);
  EtadauMean.setPlotLabel("mean");
  EtadauSigma.setPlotLabel("sigma");
  NEtadauSig.setPlotLabel("S");
  NEtadauBkg.setPlotLabel("B");
  BkgEtadauSlope.setConstant();
  EtadauMassPdf.paramOn(frVarEtadaupi0pi0pi0,Format("NELU"),Layout(0.65,0.9,0.9));
  BkgEtadauSlope.setConstant(kFALSE);
  frVarEtadaupi0pi0pi0->Draw();
  frVarEtadaupi0pi0pi0->SetTitle("#eta' #rightarrow #eta #pi^{+} #pi^{-}, #eta #rightarrow #pi^{0} #pi^{0} #pi^{0}");
  TPaveText *tptetadaupi0pi0pi0 = frVarEtadaupi0pi0pi0->findObject("TPave");
  tptetadaupi0pi0pi0->SetTextSize(0.030);
  c1->SaveAs("RooPlots/Etadau_pi0pi0pi0Mass.eps");

}

