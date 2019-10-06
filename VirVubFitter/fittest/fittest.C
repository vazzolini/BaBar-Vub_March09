#define fittest_cxx
#include "fittest.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace RooFit;

void fittest::WriteDataSet(Int_t sel, Int_t runperiod)
{
  TFile *file;
  if(sel==0){
    switch(runperiod){
    case 12: file= new TFile("/afs/slac.stanford.edu/u/br/petrella/scra/genDSet-R18-Run12.root","RECREATE");break;
    case 3: file=new TFile("/afs/slac.stanford.edu/u/br/petrella/scra/genDSet-R18-Run3.root","RECREATE"); break;
    case 4: file=new TFile("/afs/slac.stanford.edu/u/br/petrella/scra/genDSet-R18-Run4.root","RECREATE"); break;
    case 5: file= new TFile("/afs/slac.stanford.edu/u/br/petrella/scra/genDSet-R18-Run5.root","RECREATE"); break;
    case 14: file= new TFile("/afs/slac.stanford.edu/u/br/petrella/scra/genDSet-R18-Run14.root","RECREATE");break;
    case 15: file= new TFile("/afs/slac.stanford.edu/u/br/petrella/scra/genDSet-R18-Run15.root","RECREATE");break;
    }
    
    
    this->NewTree(BBbar,runperiod);  //BBbar MC
    this->Loop(BBbar);     //BBbar MC
  } else {
    switch(runperiod){
    case 12: file= new TFile("/afs/slac.stanford.edu/u/br/petrella/scra/onpeakDSet-R18-Run12.root","RECREATE");break;
    case 3: file=new TFile("/afs/slac.stanford.edu/u/br/petrella/scra/onpeakDSet-R18-Run3.root","RECREATE"); break;
    case 4: file=new TFile("/afs/slac.stanford.edu/u/br/petrella/scra/onpeakDSet-R18-Run4.root","RECREATE"); break;
    case 5: file= new TFile("/afs/slac.stanford.edu/u/br/petrella/scra/onpeakDSet-R18-Run5.root","RECREATE"); break;
    case 14: file= new TFile("/afs/slac.stanford.edu/u/br/petrella/scra/onpeakDSet-R18-Run14.root","RECREATE");break;
    case 15: file= new TFile("/afs/slac.stanford.edu/u/br/petrella/scra/onpeakDSet-R18-Run15.root","RECREATE");break;
    }
    
    this->NewTree(Data,runperiod);   //Data
    this->Loop(Data);      //Data
  }
  
  datadata->Write();
  Vmes->Write();
  file->Write();
}


void fittest::Loop(Int_t dsetselect)
{
  Double_t BTYPE=2;
  Double_t MININTPUR = 0.;
  Double_t PRMM2CUT=-3;
  int flavB;
  
  if (fChain == 0) return;

  Vmes    = new RooRealVar("mes","mes (GeV)",5.22,5.291);
  VlepYes = new RooRealVar("lepYes","SL cut",0,1); 
  VlepVub = new RooRealVar("lepVub","lepVub",0,1);
  VlepVcb = new RooRealVar("lepVcb","lepVcb",0,1);
  VlepVubSB = new RooRealVar("lepVubSB","lepVubSB",0,1);
  VflavB   = new RooRealVar("flavB","flavB",0,5);
  VlepYaSe = new RooRealVar("lepYaSe","All Cuts",0,1);
  Vwe      = new RooRealVar("weight","weight",0.,100.);
  Vchop= new RooRealVar("chop","m_X (GeV)",-10000.,10000.);
  Vq2 = new RooRealVar("q2","q2 (GeV)",-10000.,10000.);
  Vmultcat = new RooRealVar("multcat","multcat",0,5);
  Vksele = new RooRealVar("ksele","ksele",0,1);
  Vpplus = new RooRealVar("pplus","pplus",-10000.,10000.);
  Vintpur = new RooRealVar("intpur","intpur",0,1);
  //  Vtrumtch = new RooRealVar("trumtch","trumtch",0,2);
  Visch = new RooRealVar("isBch","Breco charge",0,1);
  Vmm2 = new RooRealVar("mm2","mm2",-20,30);
  Vch = new RooRealVar("ch","Total event charge",0,10);
  VmodeB = new RooRealVar("modeB","modeB",-100000.,100000.);
  VtruemodeB = new RooRealVar("truemodeB","truemodeB",-100000.,100000.);
  Vch1B = new RooRealVar("ch1B","ch1B",0,15);
  Vneu1B = new RooRealVar("ne1B","neu1B",0,15);
    
  RooArgSet args=RooArgSet(*Vmes,*Vchop,*Vq2,*Vwe,*VlepYes,*VflavB,*Vksele,*VlepYaSe);
  args.add(*Vmultcat); args.add(*Vintpur);  args.add(*Visch); args.add(*Vpplus); 
  args.add(*Vmm2); args.add(*Vch); args.add(*VtruemodeB); args.add(*VmodeB); args.add(*Vch1B); args.add(*Vneu1B);
  datadata=new RooDataSet("old","DATA",args,"weight");
  
  int ch;
  bool rCh, isLp, WdeltaCut, flav, ksele, lPYes, isLpSig, lPYesSig, lepVubSB, lepVub, lepVcb, AllCut, ischarged, wrongcharge;
  
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"Nentries = "<<nentries<<endl;
  Int_t nbytes = 0, nb = 0;
  
  Int_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if(jentry%10000==0)
      cout<<jentry<<endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
     
    if (!(TMath::Abs(brecocharge) == BTYPE || BTYPE == 2)) continue;
    ch = TMath::Abs(xcharge + brecocharge);  // total charge   
     
    rCh = (TMath::Abs(brecocharge) == BTYPE || BTYPE == 2);
    isLp = ((nle > 0) && (tlab<2.37) && (tlab>0.36) && (plab>0.5) && rCh && intpur>MININTPUR);   
    WdeltaCut = (wdeltam>PRMM2CUT && brecocharge == 0); //last cut based on wdelta
    //---- charge correlation
    flav =  lcharge + brecoflav;
    ksele = (nkp + nks) > 0;   // fit on the depleted sample?
    lPYes = (pcms > 1. && isLp && !(TMath::Abs(brecocharge)!=0 && (flav)!=0) );
    isLpSig = ((nle == 1) && (tlab<2.37) && (tlab>0.36) && (plab>0.5) && rCh && intpur>MININTPUR);  
    lPYesSig = (pcms > 1. && isLpSig && !(TMath::Abs(brecocharge)!=0 && (flav)!=0) );
    
    lepVub = (vub==1 && pcms > 1. && isLp && !(TMath::Abs(brecocharge)!=0 && (flav)!=0));
    lepVcb = (vcb==1 && pcms > 1. && isLp && !(TMath::Abs(brecocharge)!=0 && (flav)!=0));
    //bool AllCut= (lPYesSig && q2fit>0 && mm2 < 0.5 && mm2 > -100000  &&  ch <0.5 && ch > -0.5  &&  ksele ==0 && mxhadfit>0.&& mxhadfit < 5 && !(WdeltaCut));   
    AllCut= (lPYesSig && q2fit>0 && mm2 < 0.5 && mm2 > -100000  &&  ch <0.5 && ch > -0.5  && mxhadfit>0.&& mxhadfit < 5 && !(WdeltaCut));   

    flavB = 5;
    if(TMath::Abs(brecocharge)==1) flavB = 3;
    if(TMath::Abs(brecocharge)==0 && flav==0)flavB = 4;

    int mult=0;
     
    if(nchg == 1 && nneu > 0) mult = 1;
    if((nchg == 2||nchg==3) && nneu == 0) mult = 2;
    if((nchg == 2||nchg==3) && nneu > 0) mult = 3;
    if(nchg > 3 && nneu == 0) mult = 4;
    if(nchg > 3  && nneu > 0) mult = 5;
     

    ////////////////////////////////////////// my charge reconstruction true. /////////////////
    
//     wrongcharge=2;
    
//     //if(xchargegen!=-999)
//     if(brecoidtrue == brecoid)
//       wrongcharge=0;
//     else
//       wrongcharge=1;
    
    //////////////////////////////////////////////////////////////////////////////////

    if(TMath::Abs(brecocharge)==0) ischarged=false;
    if(TMath::Abs(brecocharge)==1) ischarged=true;
     
    //  mymes=mes/(sqrt(eUps*eUps-pUps*pUps)/2);  //REDUCED MES
    
    // OLD TRUTH-MATCHING
    // int combinatorial = (truemodeB != modeB) ? 1 : 0; //PEAKING + COMBINATORIAL
    // int combinatorial = (truemodeB == -1) ? 1 : 0; //COMBINATORIAL
    // int combinatorial = (truemodeB == modeB) ? 1 : 0; //SIGNAL
    // int combinatorial = ((modeB!=truemodeB)&&truemodeB != -1) ? 1 : 0; //PEAKING;
     
    if(mes>5.22 && lPYes) {
      Vmes->setVal(mes);
      VlepVub->setVal(lepVub);  
      VlepVcb->setVal(lepVcb); 
      VlepVubSB->setVal(lepVubSB); 
      Vksele->setVal(ksele);
      Vpplus->setVal(pplus);
       
      Vq2->setVal(q2);
      Vchop->setVal(mxhad);
      VlepYes->setVal(lPYes);  
      VflavB->setVal(flavB);
      VlepYaSe->setVal(AllCut);
      Vmultcat->setVal(mult);
      Vintpur->setVal(intpur);
      // Vtrumtch->setVal(istmatch);
      Visch->setVal(ischarged);
      Vmm2->setVal(mm2);
      Vch->setVal(ch);
      VtruemodeB->setVal(truemodeB);
      VmodeB->setVal(modeB);
      Vch1B->setVal(ch1B);
      Vneu1B->setVal(neu1B);

      switch(dsetselect){
      case ccbar: datadata->add(args,3.085252); break; //ccbar weight to scale to BBbar Lumi
      case udsudsbar: datadata->add(args,2.893258); break; //uds weight to scale to BBbar Lumi
      case BBbar: datadata->add(args,1.); break; //bbbar
      case Data: datadata->add(args,1.); //data No weights on data
      }
    }
  }
}
 
void fittest::Test(bool isdata,int pdfmodel,int tcmp,bool isallcut,bool isdepl,bool dumpparams,const char* pfile,float intp,int chb,float mx_low, float mx_high, float q2_low, float q2_high, float ppl_low, float ppl_high)
{
  dumppar = dumpparams;
  parfile = TString(pfile);

  TFile *file = new TFile("/nfs/farm/babar/AWG69/gaglio/FIT/DataSets/newrun16.root","READ"); 
  datadata = NULL;
  somecut = new TString("");  

  if(isdata) {
    datadata = (RooDataSet*)file->Get("DATA");
    fname = new TString("data"); tcmp = all; 
  } else {
    switch(tcmp){
    case sigonly: { fname=new TString("sig");
    datadata = (RooDataSet*)file->Get("VubOUT"); } break;
    case bkgonly: { fname=new TString("bkg"); 
    datadata = (RooDataSet*)file->Get("Vcb"); } break;
    case all: { fname=new TString("allMC"); 
    datadata = (RooDataSet*)file->Get("PStar");} break;
    }
  }

  cout << " LOADING dataset " << datadata->GetName() << " from file " <<file->GetName()<<" with total entries " << datadata->numEntries(1) << endl;
  
  switch(chb){
  case 0: { somecut->Append(" && isBch==0"); }  break;
  case 1: { somecut->Append(" && isBch==1"); } break;
  }
  

  char name[200],tmp[10];
  Vmes = NULL;
    
  RooArgSet *varlist = datadata->get();
  Vmes = (RooRealVar*)varlist->find("mes");

  //  if(isdata)
    Vmes->setRange(5.22,5.291);
  //  else
  //  Vmes->setRange(5.22,5.2891);

  x = new RooRealVar(*Vmes);
  RooRealVar bmes(*x);
  int nBins;
  nBins = int((5.291-5.22)/0.0005);
  bmes.setBins(nBins);  

  if(isdata){
    x->setRange("mesint",5.27,5.291); //AGGIUNTI nel passaggio 14-24
    bmes.setRange("bmesint",5.27,5.291);
  } else {
    x->setRange("mesint",5.27,5.291); //AGGIUNTI nel passaggio 14-24
    bmes.setRange("bmesint",5.27,5.291);
  }
  
  //ENRICHED: ksele==0; DEPLETED: ksele!=0
  if( isdepl ){
    sprintf(tmp,"ksele != 0"); //depleted sample cut
    fname->Append("_depl");}
  else {
    sprintf(tmp,"ksele == 0"); //enriched sample cut
    fname->Append("_enr");
  }
  if( isallcut )
    fname->Append("_AC");
  else
    fname->Append("_SL");

  char ip[10];
  sprintf(ip,".intp%3.2f",intp);
  //  fname->Append(ip);

  sprintf(name,"mes > 5.22  && mx > %f && mx < %f && q2 > %f && q2 < %f && pplus > %f && pplus < %f && intpur > %f %s",mx_low,mx_high,q2_low,q2_high,ppl_low,ppl_high,intp,somecut->Data());
  if(isallcut)
    sprintf(name,"mes> 5.22 && mx > %f && mx < %f && q2 > %f && q2 < %f && pplus > %f && pplus < %f && lepYaSe>0 && intpur > %f %s",mx_low,mx_high,q2_low,q2_high,ppl_low,ppl_high,intp,somecut->Data());
  
  delete somecut;
  
  cout << "Using cuts " << name << endl;
  RooDataSet *datred = dynamic_cast<RooDataSet*>(datadata->reduce(name));
  datadata = datred;
  
  int dsetentries = datadata->numEntries(1); //1 takes into account weighted events
  cout << "Number of entries in dataset after applying above selection: " << dsetentries;

  //  nsig->setVal(0.6*datared->sumEntries((std::string(simply)+"&&mes>5.27").c_str()));
  cout << "; for mES > 5.27 " << datadata->sumEntries("mes>5.27") << endl;;
  
  RooDataHist* bdatared = new RooDataHist("bmes","bmes",RooArgSet(bmes),*datadata);
  bdatared->weightError(RooAbsData::SumW2);
  
  BuildModel(isdata,pdfmodel,0,isallcut); //0 extendedmodel
  
  //   RooBinning tbins(5.22,5.2891);
  //   tbins.addUniform(100,5.22,5.2891);
  //   RooPlot* xframe = x->frame() ;
  //   RooPlot* nllframe =Thosigma_r2->frame(0.00150,0.0040,50);
  
  TCanvas *canva=new TCanvas("cannva","canva",0,0,800,600);
  //  canva->Divide(2,2);

  RooFitResult *vvf0;
  FitModel(bdatared,canva->cd(1),0,vvf0,mx_low,mx_high,q2_low,q2_high,ppl_low,ppl_high,tcmp,pdfmodel);

  // ====================== PUMP MY ERRORS SECTION ========================= //

  //   RooDataHist* newdatahist = new RooDataHist("bbmes","bbmes",RooArgSet(bmes),*datadata);

  //   TH1* histo;
  //   histo=newdatahist->createHistogram("histo_newdatahist",bmes);
  //   histo->Sumw2();
  //   //  chisq=sqrt(chisq);
  //   for(int i=0;i<132;i++){
  //     if(i==1)
  //       histo->SetBinError(i,sqrt(histo->GetBinContent(i)));
  //     histo->SetBinError(i,histo->GetBinError(i)*sqrt(4.2306/6.4844));
  //   }
  //   RooDataHist* ultimate= new RooDataHist("cmes","cmes",RooArgList(bmes), histo);  
  //   TH1* histo2;
  //   histo2=ultimate->createHistogram("daultimate",bmes); 
  //   histo2->Sumw2();
  
  //   ultimate->weightError(RooAbsData::SumW2);
  
  //   cout<<"================================================================"<<endl;
  //   cout<<"=======            NLL FIT with modified Errors   (a la VVF) =="<<endl;
  //   cout<<"================================================================"<<endl;

  //   RooFitResult *vvf1;
  //   BuildModel(isdata,0); //Extended Model
  //   FitModel(ultimate,canva->cd(2),0,vvf1);

  //   cout<<"================================================================"<<endl;
  //   cout<<"=======              CHI2 FIT with modified Errors      ========"<<endl;
  //   cout<<"================================================================"<<endl;
  
  //   RooFitResult *fitchi2;
  //   BuildModel(isdata,1); //1 Non Extended model
  //   FitModel(ultimate,canva->cd(3),1,fitchi2);

  //   cout<<"================================================================"<<endl;
  //   cout<<"=======                   NLL FIT  with modified Errors ========"<<endl;
  //   cout<<"================================================================"<<endl;

  //   RooFitResult *fitnll;
  //   BuildModel(isdata,2); //2 non extended model
  //   FitModel(ultimate,canva->cd(4),2,fitnll);
  
  cout << " Covariance Matrix ";
  switch(vvf0->covQual()) {
  case 0: cout << "Not calculated at all" ; break ;
  case 1: cout << "Approximation only, not accurate" ; break ;
  case 2: cout << "Full matrix, but forced positive-definite" ; break ;
  case 3: cout << "Full, accurate covariance matrix" ; break ;
  }
  //  cout << " OUTPUT a la VVF  " << endl;
  vvf0->Print("v");

  //cout<<"OUTPUT a la VVF 2 "<<endl;
  //vvf1->Print("v");
  //cout<<"OUTPUT CHI2 "<<endl;
  //fitchi2->Print("v");
  //cout<<"output NLL"<<endl;
  //fitnll->Print("v");

}

void fittest::FitModel(RooDataHist* rdh, TVirtualPad* pad, Int_t opt, RooFitResult*& r,float mxl, float mxh, float q2l, float q2h, float ppl, float pph,int tcmp,int pdfmodel)
{
  RooPlot *xframe=x->frame();
  char name[50];
  
  RooAbsPdf *sigpdf=NULL;
  RooAbsPdf *bkgpdf=NULL;
  RooAbsPdf *pkgpdf=NULL;

  //Argus
  bkgpdf = (RooAbsPdf*)a;
  
  if(pdfmodel == newmodel)
    sigpdf = (RooAbsPdf*)thosig;
  else
    sigpdf = (RooAbsPdf*)crystalball;


  switch (opt){
  case 0: //FitTo
    {
      r = tmodel->fitTo(*rdh,"rmhe");
      rdh->plotOn(xframe,MarkerSize(0.25)) ;
      
      Int_t floatParms = (r->floatParsFinal()).getSize();

      tmodel->plotOn(xframe,Components(RooArgSet(*bkgpdf)),LineColor(50));
      tmodel->plotOn(xframe,Components(RooArgSet(*sigpdf)),LineColor(8));
      if(pkgpdf!=NULL)
	tmodel->plotOn(xframe,Components(RooArgSet(*pkgpdf)),LineColor(2));
      
      tmodel->plotOn(xframe);
      xframe->SetTitle(fname->Data());
      double  chisq= xframe->chiSquare(floatParms);
      cout << "fitTo xframe chisquare = " <<chisq << endl;
    } break;
  case 1: //Chi2
    {
      cout<<" +_+_+_+_+_+_+_+_+_+_+_+_ CHI2 "<<endl;
      RooChi2Var rcvar("rcvar","rcvar",*tmodel,*rdh,DataError(1));
      RooMinuit minuit(rcvar);
      minuit.migrad();
      minuit.hesse();
      r=minuit.save();
      Int_t floatParms = (r->floatParsFinal()).getSize();
      
      rdh->plotOn(xframe,MarkerSize(0.25),DataError(1)); //ARGHHHH DataError(1)=to plot Symm sqrt(n) errors!
      tmodel->plotOn(xframe,Components(RooArgSet(*a)),LineColor(5));
      tmodel->plotOn(xframe,Components(RooArgSet(*thosig)),LineColor(3));
      tmodel->plotOn(xframe);
      xframe->SetTitle("Chi2");
      double  chisq= xframe->chiSquare(floatParms);
      cout <<"Chi2 xframe chisquare = " <<chisq << endl;
    }  break;
  case 2: // NLL
    {
      cout<<" +_+_+_+_+_+_+_+_+_+_+_+_ NLL "<<endl;
      RooNLLVar nll("nll","nll",*tmodel,*rdh);
      RooMinuit minuitnll(nll);
      minuitnll.migrad();
      minuitnll.hesse();
      r=minuitnll.save();

      rdh->plotOn(xframe,MarkerSize(0.25),DataError(1)); //ARGHHHH DataError(1)=to plot Symm sqrt(n) errors!
      tmodel->plotOn(xframe);  
      Int_t floatParms = (r->floatParsFinal()).getSize();
      xframe->SetTitle("NLL");
      double  chisq= xframe->chiSquare(floatParms);
      cout << "NLL xframe chisquare = " <<chisq << endl;
    }  break;
  }    
  //compute single contributions to chisquare...
  RooHist* pullhisto = xframe->pullHist();
  pullframe = x->frame();
  pullframe->addPlotable(pullhisto);
  pullframe->SetMarkerStyle(24);
  pullframe->SetNdivisions(504, "Y");
  pullframe->SetLabelSize(0.22, "X");  pullframe->SetLabelSize(0.15, "Y");
  pullframe->SetStats(0);
  pullframe->SetTitle("");
  pullframe->SetTitleSize(0.22, "Y");
  pullframe->SetTitleOffset(0.22, "Y");
  pullframe->SetYTitle("Pull");
  //pullhisto->Print("Verbose");
  
  //compute chisquare probability
      
  Double_t probchi = TMath::Prob((xframe->chiSquare(floatParms)) * 
				 (pullhisto->GetN()-floatParms), (pullhisto->GetN()-floatParms));
  cout << "chisq probability = " << probchi << endl;
      
  TPaveText* tbox1 = new TPaveText(0.59, 0.94, 0.79, 0.99, "BRNDC");
  char line1[50];
  sprintf(line1, "#chi^{2} = %5.4f",  xframe->chiSquare(floatParms));
  tbox1->AddText(line1);
  
  cout << " control : number of floating params ----> " << floatParms << endl; 

  xframe->addObject(tbox1);
  pad->cd();
      
  //CB put also pull distributions
  mPad = new TPad("mES plot", "", 0.00, 0.25, 0.99, 0.99);   mPad->Draw(); 
  pPad = new TPad("pull plot", "", 0.00, 0.05, 0.99, 0.245);   pPad->Draw(); 
  
  mPad->cd();    shrinkPad(0.1,0.1,0.1,0.1);
  //      xframe->SetTitle("All MC: Argus+ccb+Thorsten sig");
  xframe->Draw();
      
  pPad->cd();    shrinkPad(0.1,0.1,0.1,0.1);
  pullframe->GetYaxis()->SetTitleSize(0.4);
  pullframe->GetXaxis()->SetTitleSize(0.4);
  pullframe->GetYaxis()->SetLabelSize(0.13);
  pullframe->GetXaxis()->SetLabelSize(0.13);
  pullframe->Draw();

  //==============================================================================//
  //==================== SOME USEFUL OUTPUT  =====================================//
  //==============================================================================//
  
  Double_t sigfra, arfra;
  
  Double_t sigfra = sigpdf->createIntegral(*x,NormSet(*x),Range("mesint"))->getVal();
  Double_t arfra = bkgpdf->createIntegral(*x,NormSet(*x),Range("mesint"))->getVal();
  
  Float_t ccbY, ccbYerr, arY, arYerr, sigY, sigYerr, soverp(0), errsoverp(0);
  Float_t arshpar = argpar->getVal();

  
    
  cout << "============== WHOLE MES RANGE ============" << endl;
  cout << "ARGUS " << argusfrac->getVal() << " +- " << argusfrac->getError() << endl;
  cout << "SIGNAL " << sigpdffraction->getVal() << " +- " << sigpdffraction->getError() << endl;
  cout << "=========== REDUCED MES RANGE: " << x->getMin() << " < mes < " << x->getMax() << " ==============" << endl;
  cout << "Fraction of ARGUS PDF in reduced range: " << arfra << endl;
  cout << "Fraction of Signal PDF in reduced range: " << sigfra << endl;
  cout << "--------------- ENTRIES NOW " << endl;
  cout << " Events in dataset mES > 5.27: " << rdh->sumEntries(0,"bmesint") << " TOTAL (whole mES): " << rdh->numEntries(1) << endl;
  cout << "ARGUS " << arfra*argusfrac->getVal() << " +- " << arfra*argusfrac->getError() << endl;
  cout << "SIGNAL " << sigfra*sigpdffraction->getVal() << " +- " << sigfra*sigpdffraction->getError() << endl;
  
  sprintf(name,"%s.txt",fname->Data());
  FILE *output = fopen(name,"a");
  arY = arfra*argusfrac->getVal();  arYerr = arfra*argusfrac->getError();
  sigY = sigfra*sigpdffraction->getVal(); sigYerr = sigfra*sigpdffraction->getError();
  arshpar = argpar->getVal();
  
  //Compute S/P and it's error
  
  if(arY !=0 && sigY !=0){
    soverp = sigY/arY;
    errsoverp = soverp*sqrt(sigYerr*sigYerr/(sigY*sigY)+arYerr*arYerr/(arY*arY));
  }
  if(q2l!=-999) // mx-q2 
    fprintf(output,"%3.2f\t %3.2f\t %3.2f\t %3.2f\t %5.0f +- %5.0f %5.0f +- %5.0f %5.3f +- %5.3f %f %f\n",mxl,mxh,q2l,q2h,arY,arYerr,sigY,sigYerr,soverp,errsoverp,chisq,arshpar);
  else 
    if(ppl!=-100) // pplus
      fprintf(output,"%3.2f\t %3.2f\t %5.0f +- %5.0f %5.0f +- %5.0f %5.3f +- %5.3f %f %f\n",ppl,pph,arY,arYerr,sigY,sigYerr,soverp,errsoverp,chisq,arshpar);
    else // mx
      fprintf(output,"%3.2f\t %3.2f\t %5.0f +- %5.0f %5.0f +- %5.0f %5.3f +- %5.3f %f %f\n",mxl,mxh,arY,arYerr,sigY,sigYerr,soverp,errsoverp,chisq,arshpar);
  
  //    RooArgList sigparlist=RooArgList(*CB_mean,*CB_sigma,*CB_alpha,*CB_n,*sigpdffraction);
  RooArgList sigparlist = RooArgList(*Thor,*Thosigma_r1,*Thoxc,*Thosigma_r2,*Thosigma_L,*Thon,*Thoalpha,*sigpdffraction);
  RooArgList argusparlist = RooArgList(*cutoff,*argpar,*argusfrac);

  
  fclose(output);
  
  if(dumppar){
    sprintf(name,"%s_params.txt",fname->Data());
    fout=fopen(name,"w");
    fprintf(fout,"# Name\t\t fitted value\t\t\t\t Min\t Max\t IsFixed\t FixedValue\n");
  }
  
  controlParameters(sigparlist,name,tcmp);
  controlParameters(argusparlist,name,tcmp);
  //  controlParameters(ccbparlist,name,tcmp);
  
 //  switch(tcmp){
//   case sigonly: {controlParameters(sigparlist,name,tcmp); } break;
//   case bkgonly: {
//     if(pdfmodel==twopdfs)
//       controlParameters(argusparlist,name,tcmp);
//     else {
//       controlParameters(argusparlist,name,tcmp); controlParameters(ccbparlist,name,tcmp); }
//   } break;
//   case all: { 
//     if(pdfmodel==twopdfs){
//       controlParameters(sigparlist,name,tcmp); controlParameters(argusparlist,name,tcmp);
//     } else {
//       controlParameters(sigparlist,name,tcmp); controlParameters(argusparlist,name,tcmp); controlParameters(ccbparlist,name,tcmp);}
//   } break;
//   }
  
  delete tempsigy; delete tempargy;
  
  if(dumppar){
    fclose(fout);
    cout<< " PDF parameters written in "<<name<<endl;
  }
      
  if(q2l != -999)
    // This is for mx,q2 scans
    sprintf(name,"%s_%3.2f%3.2f%3.2f%3.2f.eps",fname->Data(),mxl,mxh,q2l,q2h);
  else 
    if(ppl != -100) //This is for pplus
      sprintf(name,"%s_%3.2f%3.2f.eps",fname->Data(),ppl,pph);
    else
      // This is for mx scan
      sprintf(name,"%s_%3.2f%3.2f.eps",fname->Data(),mxl,mxh);
      
  pad->SaveAs(name);
}

void fittest::BuildModel(bool isdata,int pdfmodel,int fittype,bool isallcut)
{
  ResetArgus(isdata,fittype,isallcut);
  if(pdfmodel == oldmodel)
    ResetCrystalBall(isdata,fittype,isallcut);
  else
    ResetThosig(isdata,fittype,isallcut);
  
  if(a != NULL) delete a;
  if(thosig != NULL) delete thosig;
  if(crystalball != NULL) delete crystalball;
  if(tmodel != NULL) delete tmodel;

  if(ae != NULL) delete ae;
  if(se != NULL) delete se;

  a = new RooArgusBG("a","Argus PDF",*x,*cutoff,*argpar);
  ae = new RooExtendPdf("ae", "ae", *a, *argusfrac);

  if(pdfmodel == oldmodel)
    crystalball = new RooCBShape("Crystalball sig","Crystalball sig",*x,*CB_mean,*CB_sigma,*CB_alpha,*CB_n);
  else
    thosig = new thosig("Thorsten sig","Thorsten sig",*x,*Thor,*Thosigma_r1,*Thoxc,*Thosigma_r2,*Thosigma_L,*Thon,*Thoalpha);
  
  cout << "PDF MODEL " << pdfmodel << endl;

  switch(pdfmodel) {
  case oldmodel: se = new RooExtendPdf("se","signal extended",*crystalball,*sigpdffraction); break;
  case newmodel: se = new RooExtendPdf("se","signal extended",*thosig, *sigpdffraction); break;
  }

  if(fittype<1) //if 0 is extended, if >0 is non extended
    switch(pdfmodel){
    case oldmodel: tmodel=new RooAddPdf("oldmodel","a+CB",RooArgList(*ae,*se),RooArgList(*argusfrac,*sigpdffraction)); break;
    case newmodel: tmodel=new RooAddPdf("newmodel","a+frank",RooArgList(*ae,*se),RooArgList(*argusfrac,*sigpdffraction)); break;
    }
}
// ==========================================================================
// ==== This function resets the parmeters for ARGUS PDF ====================
// ==========================================================================
void fittest::ResetArgus(bool isdata, int fittype,bool isallcut)
{
  if(argpar!=NULL) delete argpar;
  if(cutoff!=NULL) delete cutoff;
  if(argusfrac!=NULL) delete argusfrac;

  // Read Parameters
  argpar = readParameters(parfile.Data(),"ar"); argpar->SetTitle("Argus Shape Parameter");
  cutoff = readParameters(parfile.Data(),"cutoff"); cutoff->SetTitle("Argus Cutoff"); 
  argusfrac = readParameters(parfile.Data(),"argusfrac"); argusfrac->SetTitle("fraction of argus");

  tempargy = new RooRealVar(*argusfrac);

  argusfrac->setMax(datadata->numEntries(kTRUE));
  argusfrac->setVal(datadata->sumEntries("mes<5.27"));

  if(fittype>1){
    argusfrac->setMax(1);
    argusfrac->setMin(0);
    argusfrac->setVal(0.5);
  }
//    if( !isdata )
//      { cutoff->setVal(5.2891); cutoff->setConstant(); }
  
}

void fittest::ResetThosig(bool isdata, int fittype,bool isallcut)
{
  if(Thosigma_L != NULL) delete Thosigma_L;
  if(Thosigma_r1 != NULL) delete Thosigma_r1;
  if(Thosigma_r2 != NULL) delete Thosigma_r2;
  if(Thor != NULL) delete Thor;
  if(Thon != NULL) delete Thon;
  if(Thoalpha != NULL) delete Thoalpha;
  if(Thoxc != NULL) delete Thoxc;
  if(sigpdffraction != NULL) delete sigpdffraction;

  Thosigma_L = readParameters(parfile.Data(),"sigma_L"); Thosigma_L->SetTitle("sigma_L");
  Thosigma_r1 = readParameters(parfile.Data(),"sigma_r1"); Thosigma_r1->SetTitle("sigma_r1");
  Thosigma_r2 = readParameters(parfile.Data(),"sigma_r2"); Thosigma_r2->SetTitle("sigma_r2");
  Thor = readParameters(parfile.Data(),"Thosig_r"); Thor->SetTitle("Thosig r");
  Thon = readParameters(parfile.Data(),"Thosig_n"); Thon->SetTitle("Thosig n");
  Thoalpha = readParameters(parfile.Data(),"Thosig_alpha"); Thoalpha->SetTitle("Thosig alpha");
  Thoxc = readParameters(parfile.Data(),"Thosig_xc"); Thoxc->SetTitle("Thosig xc");
  sigpdffraction = readParameters(parfile.Data(),"sigpdffraction"); sigpdffraction->SetTitle("fraction of thorsten");

  tempsigy = new RooRealVar(*sigpdffraction);

  sigpdffraction->setMax(datadata->numEntries(kTRUE));
  //  sigpdffraction->setMin(sigpdffraction->getMin()*scalefactor);
  sigpdffraction->setVal(0.6*datadata->sumEntries("mes>5.27"));
  
  if(fittype>1){
    sigpdffraction->setMax(1);
    sigpdffraction->setMin(0);
    sigpdffraction->setVal(0.5);
  }
  //PARAMETERS THAT I GOT FITTING ALL CUTS MC ENRICHED
  //   Thoxc->setVal(5.27968); Thoxc->setConstant();
  //   Thon->setVal(1.19604); Thon->setConstant();
  //   Thosigma_r2->setVal(0.00236417); Thosigma_r2->setConstant();
  //   Thor->setVal(0.707672); Thor->setConstant();
  //   Thoalpha->setVal(4.39726); Thoalpha->setConstant();
  //   Thosigma_L->setVal(0.00192949); 
  //   Thosigma_r1->setVal(0.00177180); 

  //PARAMETERS THAT I GOT FITTING ALL CUTS MC DEPLETED
  // Thoxc->setVal(5.27968); Thoxc->setConstant();
  //   Thon->setVal(1.34627); Thon->setConstant();
  //   Thosigma_r2->setVal(0.00236417); Thosigma_r2->setConstant();
  //   Thor->setVal(0.755292); Thor->setConstant();
  //   Thoalpha->setVal(4.40201); Thoalpha->setConstant();
  //   Thosigma_L->setVal(0.00192597); 
  //   Thosigma_r1->setVal(0.00169737); 

  //CB for data depleted: refit sigma_L and sigma_r1 on entire sample 
  //CB and allcuts, take other parameters from entire MC depl all cuts
  // Thoxc->setVal(5.27952); //Thoxc->setConstant();
  //   Thon->setVal(1.39619); //Thon->setConstant();
  //   Thosigma_r2->setVal(0.002496); //Thosigma_r2->setConstant();
  //   Thor->setVal(0.633032); //Thor->setConstant();
  //   Thoalpha->setVal(4.24602); //Thoalpha->setConstant();
  //   Thosigma_L->setVal(0.0015321); 
  //   Thosigma_r1->setVal(0.0022168); 
  
  if(datadata->sumEntries("mes>5.27") < 1000.) {
    Thosigma_L->setConstant();
    Thosigma_r1->setConstant();
  }
  

  // Thon->setVal(1.24698); Thon->setConstant();
  // Thosigma_L->setVal(0.001977); Thosigma_L->setConstant();
  // Thosigma_r1->setVal(0.0016393); Thosigma_r1->setConstant();
  // Thosigma_r2->setVal(0.00226044); Thosigma_r2->setConstant();
  // Thor->setVal(0.936655); Thor->setConstant();
  // Thoalpha->setVal(4.19962); Thoalpha->setConstant();
  
  //    WARNING FOR SCAN ONLY REMOVE DATA DEPLETED!!!
  //    Thoxc->setVal(5.27950); Thoxc->setConstant();
  //    Thosigma_L->setVal(0.00157); Thosigma_L->setConstant();
  //    Thosigma_r1->setVal(0.00188); Thosigma_r1->setConstant();
  //    cutoff->setVal(5.28931); cutoff->setConstant();
  //    Thoxc->setVal(5.27971); Thoxc->setConstant(); //FOR ALLCUTS ENRICHED MC 5.27975
  //    Thosigma_L->setVal(0.00201350); Thosigma_L->setConstant();  //FOR ALLCUTS ENRICHED MC 0.00201657
  //    Thosigma_r1->setVal(0.00156159); Thosigma_r1->setConstant();  //FOR ALLCUTS ENRICHED MC 0.00150645
  //    Thoxc->setVal(5.27970); Thoxc->setConstant(); //FOR ALLCUTS DEPLTEED MC
  //    Thosigma_L->setVal(0.00201657); Thosigma_L->setConstant();  //FOR ALLCUTS DEPLETED MC
  //    Thosigma_r1->setVal(0.00150645); Thosigma_r1->setConstant();  //FOR ALLCUTS DEPLETED MC
}
void fittest::ResetCrystalBall(bool isdata,int fittype, bool isallcut)
{
  if(CB_mean!=NULL) delete CB_mean;
  if(CB_sigma!=NULL) delete CB_sigma;
  if(CB_alpha!=NULL) delete CB_alpha;
  if(CB_n!=NULL) delete CB_n;
  if(sigpdffraction!=NULL) delete sigpdffraction;

  CB_mean=readParameters(parfile.Data(),"mean_CB"); CB_mean->SetTitle("Mean of Cry.Ball");
  CB_sigma=readParameters(parfile.Data(),"sigma_CB"); CB_sigma->SetTitle("Sigma of Cry.Ball");
  CB_alpha=readParameters(parfile.Data(),"alpha_CB"); CB_alpha->SetTitle("alpha of Cry.Ball");
  CB_n=readParameters(parfile.Data(),"n_CB"); CB_n->SetTitle("n of Cry.Ball");
  sigpdffraction=readParameters(parfile.Data(),"sigpdffraction"); sigpdffraction->SetTitle("fraction of Cry.Ball");

  tempcrystalbally=new RooRealVar(*sigpdffraction);

  sigpdffraction->setMax(datadata->numEntries(kTRUE));
  sigpdffraction->setVal(datadata->sumEntries("mes > 5.27"));
      
  if(fittype>1){
    sigpdffraction->setMax(1);
    sigpdffraction->setMin(0);
    sigpdffraction->setVal(0.5);
  }
}

//========================================================================
//=====  Function to control PDF's parameters after a fit  ===============
//========================================================================
void fittest::controlParameters(const RooArgList& var,const char* filename, int tcmp)
{
  RooRealVar *current;

  for(int i=0;i<var.getSize();i++)
    {
      current=(RooRealVar*)var[i];
      float passo=(current->getMax()-current->getMin())/10.;
      bool insideup = current->getVal() < (current->getMax()-passo);
      bool insidedown = current->getVal() > (current->getMin()+passo);
      if(!insideup){
	cout<<" !!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!! "<<current->getTitle()<<" is near ( < 10%) to its upper limit: "<<current->getVal();
	cout<<" Range= ["<<current->getMin()<<","<<current->getMax()<<"]"<<endl;
      }
      if(!insidedown){
	cout<<" !!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!! "<<current->getTitle()<<" is near ( < 10%) to its lower limit: "<<current->getVal();
	cout<<" Range= ["<<current->getMin()<<","<<current->getMax()<<"]"<<endl;
      }
      if(dumppar){
	if(strcmp(current->GetName(),"sigpdffraction")==0) current=tempsigy;
	if(strcmp(current->GetName(),"argusfrac")==0) current=tempargy;
	//	if(strcmp(current->GetName(),"ccbfrac")==0) current=tempccby;
	
	if(tcmp==all){
	  fprintf(fout,"%8s\t %f\t +- %f\t %f\t %f\t %d\t%f\n",current->GetName(),
		  current->getVal(),current->getError(),current->getMin(),current->getMax(),
		  current->isConstant(),current->getVal()); 
	} else if (tcmp == sigonly) {
	  if(strcmp(current->GetName(),"sigpdffraction")==0) continue;
	  fprintf(fout,"%8s\t %f\t +- %f\t %f\t %f\t %d\t%f\n",current->GetName(),
		  current->getVal(),current->getError(),current->getMin(),current->getMax(),
		  1,current->getVal()); 	
	} else if (tcmp == bkgonly) {
	  //	  if(strcmp(current->GetName(),"ccbfrac")==0 || strcmp(current->GetName(),"argusfrac")==0) continue;
	  fprintf(fout,"%8s\t %f\t +- %f\t %f\t %f\t %d\t%f\n",current->GetName(),
		  current->getVal(),current->getError(),current->getMin(),current->getMax(),
		  0,current->getVal()); 	
	  else {
	    fprintf(fout,"%8s\t %f\t +- %f\t %f\t %f\t %d\t%f\n",current->GetName(),
		    current->getVal(),current->getError(),current->getMin(),current->getMax(),
		    1,current->getVal()); 	
	  }
	}
      }
    }
}
//========================================================================
//=====  Function to read PDF's parameters from a txt file ===============
//========================================================================
RooRealVar* fittest::readParameters(const char* filename,const char* varname)
{
  ifstream fin;
  fin.open(filename);
  char buffer[100];
  char name[20];
  float start,min,max,fixval;
  int isfix;

  RooRealVar *n=NULL;

  while(fin.getline(buffer,100,'\n')){
    if(buffer[0]=='#') continue;
    sscanf(buffer,"%s %f %*s %*f %f %f %d %f %*s",name,&start,&min,&max,&isfix,&fixval);
    //    cout<<"INIZIALIZING "<<name<<" "<<start<<" "<<min<<" "<<max<<" "<<isfix<<" "<<fixval<<endl;
    if(strcmp(varname,name)==0)
      if(isfix==1){
	n=new RooRealVar(name,"",fixval,min,max);
	n->setVal(fixval);
	n->setConstant();
      } 
      else
	n=new RooRealVar(name,"",start,min,max);
  }
  
  if (n==NULL)
    cout<<"WARNING "<<varname<<" NOT FOUND in "<<filename<<endl;
  return n;
}

//Functions for graphics management
// ----------------------------------------------------------------------
void fittest::shrinkPad(double b, double l, double r, double t) {
  gPad->SetBottomMargin(b); 
  gPad->SetLeftMargin(l);
  gPad->SetRightMargin(r);
  gPad->SetTopMargin(t);
}
