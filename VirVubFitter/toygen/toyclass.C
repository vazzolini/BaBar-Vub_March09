#define toyclass_cxx
#include "toyclass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "RooThorstenSig.hh"
#include "RooFitModels/RooArgusBG.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooExtendPdf.hh"
#include "RooFitCore/RooAddPdf.hh"
#include <cassert>

toyclass::toyclass(const TString& filename) {
  //  cout << "toyclas "<<filename.Data()<<endl;

  TFile *file = new TFile(filename,"READ");
  if(!file) {
    cout << filename.Data() << " NOT FOUND !!! " << endl;
    assert(file->IsOpen());
  }  
  TTree *t = (TTree*)file->Get("PdfTree");
  if(t == NULL)
    cout << "ERROR! PdfTree Tree not found in file " << filename.Data() << endl;
  else {
    cout << t->GetName() << " tree successfully found in " << filename.Data() << endl;
    Init(t);
  }
  RetrieveInfo();
}

void toyclass::Loop(Double_t sfd, Double_t sfg, Double_t sfs, Int_t start, Int_t stop)
{
  TString name = Form("PdfGenFile_%d-%d.root",start,stop);
  
  TFile *fout = new TFile(name,"RECREATE");
  
  Double_t scalefactor;

  RooRealVar mes = RooRealVar("mes","mes",5.22,5.2979);
  if (fChain == 0) return;

  RooRealVar r = RooRealVar("r","r",0);
  RooRealVar s_r1 =  RooRealVar("sigma_r1","sigma_r1",0);
  RooRealVar xc =  RooRealVar("xc","xc",0);
  RooRealVar s_r2 = RooRealVar("sigma_r2","sigma_r2",0);
  RooRealVar s_l =  RooRealVar("sigma_l","sigma_l",0);
  RooRealVar n =  RooRealVar("n","n",0);
  RooRealVar alpha = RooRealVar("alpha","alpha",0);
  
  RooRealVar pArgPar = RooRealVar("ar", "argus shape parameter", -100., -10.);
  RooRealVar pCutOff = RooRealVar("cutoff", "argus cutoff", 5.288, 5.292);
  
  RooRealVar nsig = RooRealVar("S","number of sig events", 0., 500000000.);
  RooRealVar nbkg = RooRealVar("B","number of bkg events", 0., 500000000.);

  RooDataSet *rds = NULL;

  Long64_t nentries = fChain->GetEntriesFast();

  Int_t nbytes = 0, nb = 0;

  stop = stop > nentries ? nentries : stop;
  start = start < 0 ? 0 : start;

  for (Long64_t jentry = start; jentry < stop+1; jentry++) {
    //    if(jentry % 10 == 0) 
    cout << " Processing mesFit " << jentry << " / " << nentries << " ( " << Form("%2.1f",(Double_t)jentry/nentries*100) << "% )" << endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
    
    switch(identifier) {
    case 0: { scalefactor = sfd; name = "Data";} break;
    case 1: { scalefactor = sfg; name = "PStar";} break;
    case 2: { scalefactor = sfs; name = "VubIN"; } break;
    case 3: { scalefactor = sfs; name = "VubOUT"; } break;
    case 4: { scalefactor = sfg; name = "Vcb"; } break;
    case 5: { scalefactor = sfg; name = "VcbOth"; } break;
    default: scalefactor = 0; break;
    }
    
    r.setVal(ThoSigR);
    s_r1.setVal(sigma_r1);
    xc.setVal(ThoSigXc);
    s_r2.setVal(sigma_r2);
    s_l.setVal(sigma_l);
    n.setVal(ThoSigN);
    alpha.setVal(ThoSigAlpha);

    pArgPar.setVal(ar);
    pCutOff.setVal(cutoff);

    nsig.setVal(S);
    nbkg.setVal(B);
    
    RooThorstenSig *thosig = new RooThorstenSig("thosig","thosig",mes,r,s_r1,xc,s_r2,s_l,n,alpha);
    RooArgusBG *ar = new RooArgusBG("argus","argus",mes,pCutOff,pArgPar);
    
    RooExtendPdf* ae = new RooExtendPdf("ae","ae", *ar, nbkg);
    RooExtendPdf* se = new RooExtendPdf("se","se", *thosig, nsig);

    RooAddPdf *model = new RooAddPdf("model","se+ae",RooArgList(*ae,*se),RooArgList(nbkg,nsig));
    
    Int_t togenerate = scalefactor*(S+B);
    cout << "Dataset " << name.Data() << " has " << S+B << " events. Scalefactor " << scalefactor
	 << " Now generating " << togenerate << " events " << endl;

//      if( togenerate == 0 )
//        continue;

    togenerate = 0;
    
    rds = model->generate(mes,togenerate);
    rds->SetTitle(name);
    name = Form("dataset_%d",mesfitIndex);
    rds->SetName(name);
    rds->Write();
    
    delete model;
    delete se; delete ae;
    delete ar; delete thosig;
  }
  fout->Close();
}

void toyclass::RetrieveInfo() {

  TH1D *histo;
  fChain->Draw("luminosity >> htemp","identifier == 0");
  histo = (TH1D*)gDirectory->Get("htemp");
  LUMI_DATA_OLD = histo->GetMean();

  if(histo->GetRMS() > 0.0001) 
    cout << "POSSIBLE ERRORS!! LUMI_DATA_OLD has RMS! It is supposed to have no RMS! " << endl;
  histo->Reset();
  delete histo; 

  fChain->Draw("luminosity >> htemp1","identifier == 1");
  histo = (TH1D*)gDirectory->Get("htemp1");
  LUMI_GENERIC_OLD = histo->GetMean();
  
  if(histo->GetRMS() > 0.0001) 
    cout << "POSSIBLE ERRORS!! LUMI_GENERIC_OLD has RMS: " << histo->GetRMS() << " It is supposed to have no RMS! " << endl;

  delete histo;

  fChain->Draw("luminosity >> htemp2","identifier == 2");
  histo = (TH1D*)gDirectory->Get("htemp2");
  LUMI_SIGNAL_OLD = histo->GetMean();

  if(histo->GetRMS() > 0.0001) 
    cout << "POSSIBLE ERRORS!! LUMI_SIGNAL_OLD has RMS: " << histo->GetRMS() << " ! It is supposed to have no RMS! " << endl;

}
