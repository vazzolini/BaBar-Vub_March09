void fitMes(){

  ifstream settingfile;
  settingfile.open("RooTxtFiles/fitsettings.txt");
  char type[100];
  settingfile >> type; 
  if(!settingfile.good()) break;
  settingfile.close();

  char title[100]; 
  if(!strcmp(type,"signal")) sprintf(title,"Semileptonic B^{+} events - Signal MC"); 
  if(!strcmp(type,"cocktail")) sprintf(title,"Semileptonic B^{+} events - Cocktail MC");
  if(!strcmp(type,"generic")) sprintf(title,"Semileptonic B^{+} events - Generic MC"); 
  if(!strcmp(type,"data")) sprintf(title,"Semileptonic B^{+} events - Data");

  char name[100]; 
  sprintf(name,"%s%s%s%s","RooPlots/","mes_",type,".eps"); 


  // Load RooFit libraries

  gSystem->Load("../../RooFitCore/tmp/libRooFitCore.so");
  gSystem->Load("../../RooFitModels/tmp/libRooFitModels.so");

  using namespace RooFit;

  // Load Data

  TFile *datafile = new TFile("RooRootFiles/mydata.root");
  TTree *tree = (TTree*)datafile->Get("ntp2");

  gROOT->cd();

  // List of variables saved in the dataset
  // Integer variables are now RooCategories

  // Breco candidate
  RooRealVar* mes = new RooRealVar("mes","m_{ES}",5.2,5.3,"GeV/c^{2}");
  RooRealVar* ass_deltapB = new RooRealVar("ass_deltapB","ass_deltapB",-999,999,"GeV/c");

  RooCategory brecocharge("brecocharge", "Breco charge");
  brecocharge.defineType("Bp",1); brecocharge.defineType("Bm",-1); brecocharge.defineType("Bz",0);

  RooCategory brecoflav("brecoflav","Breco flavor (norm)");
  brecoflav.defineType("b",-1); brecoflav.defineType("b bar",1);

  RooCategory nle("nle", "number of leptons above 1 GeV");
  nle.defineType("nolep",0);
  for(int i=1; i<11; i++){
    char numlep[100];
    sprintf(numlep,"%d",i);
    nle.defineType(numlep,i);
  }

  RooCategory nlept500("nlept500", "number of leptons above 500 MeV");
  nlept500.defineType("nolep",0);
  for(int i=1; i<11; i++){
    char numlep500[100];
    sprintf(numlep500,"%d",i);
    nlept500.defineType(numlep500,i);
  }

  RooCategory lcharge("lcharge", "lepton charge");
  lcharge.defineType("lepton",-1);  lcharge.defineType("anti-lepton",1);

  RooCategory isele("isele", "lepton type");
  isele.defineType("electron",1);  isele.defineType("muon",0);

  RooRealVar* pcms = new RooRealVar("pcms","momentum CMS",-999,999,"GeV/c");
  RooRealVar* plab = new RooRealVar("plab","momentum LAB",-999,999,"GeV/c");
  RooRealVar* tlab = new RooRealVar("tlab","theta LAB",-999,999,"");

  RooCategory Gvxbtyp("Gvxbtyp", "MC truth");
  Gvxbtyp.defineType("noMCTruth",-999);
  for(int i=-30; i<31; i++){
    char mcid[100];
    sprintf(mcid,"%d",i);
    Gvxbtyp.defineType(mcid,i);
  }

  // ETA L NU 
  RooCategory nrecoEta("nrecoEta", "Number of #eta");
  nrecoEta.defineType("norecoEta",0);
  for(int i=1; i<101; i++){
    char numeta[100];
    sprintf(numeta,"%d",i);
    nrecoEta.defineType(numeta,i);
  }

  RooCategory indexbestEta("indexbestEta", "Index of best #eta");
  indexbestEta.defineType("noindexEta",-999);
  for(int i=0; i<100; i++){
    char indexeta[100];
    sprintf(indexeta,"%d",i);
    indexbestEta.defineType(indexeta,i);
  }

  RooCategory modebestEta("modebestEta", "Eta decay mode");
  modebestEta.defineType("nomodeEta",-999); modebestEta.defineType("gammagamma",1);
  modebestEta.defineType("pipipi0",2); modebestEta.defineType("pi0pi0pi0",3);

  RooRealVar* barembestEta = new RooRealVar("barembestEta","M(#eta)",-999,999,"GeV/c^{2}"); 
   
  // ETA PRIME L NU 
  RooCategory nrecoEtap("nrecoEtap", "Number of #eta'");
  nrecoEtap.defineType("norecoEtap",0);
  for(int i=1; i<101; i++){
    char numetap[100];
    sprintf(numetap,"%d",i);
    nrecoEtap.defineType(numetap,i);
  }

  RooCategory indexbestEtap("indexbestEtap", "Index of best #eta'");
  indexbestEtap.defineType("noindexEtap",-999);
  for(int i=0; i<100; i++){
    char indexetap[100];
    sprintf(indexetap,"%d",i);
    indexbestEtap.defineType(indexetap,i);
  }

  RooCategory modebestEtap("modebestEtap", "Etap decay mode");
  modebestEtap.defineType("nomodeEtap",-999);
  modebestEtap.defineType("rho0gamma",1);  modebestEtap.defineType("etapipi_gammagamma",2);
  modebestEtap.defineType("etapipi_pipipi0",3);  modebestEtap.defineType("etapipi_pi0pi0pi0",4);

  RooRealVar* barembestEtap = new RooRealVar("barembestEtap","M(#eta')",-999,999,"GeV/c^{2}"); 
  RooRealVar* Rho0massdaubestEtap = new RooRealVar("Rho0massdaubestEtap","M(#rho0)",-999,999,"GeV/c^{2}"); 
  RooRealVar* EtamassdaubestEtap = new RooRealVar("EtamassdaubestEtap","M(#eta)",-999,999,"GeV/c^{2}"); 

  // Add variables to the dataset

  RooArgSet myVars;

  // Breco  
  myVars.add(*mes);
  myVars.add(*ass_deltapB);
  myVars.add(brecocharge);
  myVars.add(brecoflav);
  // Best lepton
  myVars.add(nle);
  myVars.add(nlept500); 
  myVars.add(lcharge);
  myVars.add(isele); 
  myVars.add(*pcms);
  myVars.add(*plab);
  myVars.add(*tlab);
  // MC truth
  myVars.add(Gvxbtyp);
  // ETA L NU                                                                    
  myVars.add(nrecoEta); 
  myVars.add(indexbestEta);                                                     
  myVars.add(modebestEta); 
  myVars.add(*barembestEta);                                                     
  // ETA PRIME L NU 
  myVars.add(nrecoEtap);                                                        
  myVars.add(indexbestEtap); 
  myVars.add(modebestEtap);                                                     
  myVars.add(*barembestEtap); 
  myVars.add(*Rho0massdaubestEtap);                                              
  myVars.add(*EtamassdaubestEtap); 
  
  // Create ROOT file containing the dataset

  RooDataSet* mydata = new RooDataSet("dataset","dataset",tree,myVars);
  
  TFile f("RooDataSets/dataset.root","RECREATE");
  mydata->Write();
  f.Close();
  
  // Create variables, parameters and pdf's for fit...

  RooRealVar mES("mes","m_{ES}",5.2,5.3,"GeV/c^{2}");

  RooRealVar Rm("mean","mean of gaussian 1",5.28,5.275,5.285) ;
  RooRealVar Rs("sigma","width of gaussian",.003,.002,.004) ;
  RooRealVar Ra("alpha","alpha parameter",1.3,0.,10.) ;
  RooRealVar Rn("n","n parameter",3.46,1., 7.) ;
  RooCBShape cb("cb","Crystal Ball",mES,Rm,Rs,Ra,Rn) ;  

  RooRealVar argpar("argus shape","argus shape parameter",-60.,-95.,-10.) ;
  RooRealVar cutoff("cutoff","argus cutoff",5.29) ;
  RooArgusBG a("a","Argus PDF",mES,cutoff,argpar) ;

  //Signal and Background events
  RooRealVar nsig("S","number of sig events",100.,0.,1500000.);
  RooRealVar nbkg("B","number of bkg events",1.,0.,500000.);

  // Add the components for the extended likelihood fit
  mES.setRange("mesint",5.27,5.29);     
  RooExtendPdf cbe("cbe","cbe",cb,nsig,"mesint");
  RooExtendPdf ae("ae","ae",a,nbkg,"mesint");
  RooAddPdf model("model","a+cb",RooArgList(cbe,ae)) ;

  // Reduce the dataset 
  TFile mydataset("RooDataSets/dataset.root"); 
  RooDataSet *mydata = (RooDataSet*)mydataset.Get("dataset"); 
  RooDataSet *dataslice = mydata->reduce("(abs(brecocharge)==1)&&!(abs(brecocharge)!=0&&(lcharge + brecoflav)!=0)&&(nlept500>0)&&(nle<=1)&&(tlab<2.37)&&(tlab>0.36)&&(plab>0.5)&&((isele==1&&pcms>0.5)||(isele==0&&pcms>0.8))"); 


  // Fit 

  RooFitResult* r = model.fitTo(*dataslice,"rmhe"); r->Print("v");
  RooPlot * xframe = mES.frame(40);
  xframe->SetTitleOffset(0.09,"y");
  dataslice->plotOn(xframe) ; 
  model.plotOn(xframe,Components(RooArgSet(cbe,ae))) ; 

  // Parameters not shown in the fit box: uncomment lines below
  //  nsig.setConstant();
  //  nbkg.setConstant();
  //  Rm.setConstant(); 
  //  Rs.setConstant(); 
  //  Ra.setConstant();
  //  Rn.setConstant();
  //  argpar.setConstant();

  // Show fit box: uncomment line below
  model.paramOn(xframe,Format("NEU"),Layout(0.10,0.60,0.8)); 
 
  // Save plot

  xframe->Draw();  
  xframe->SetTitle(title);
  c1->SaveAs(name);

}

