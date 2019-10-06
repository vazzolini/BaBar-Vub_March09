#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <math.h>
#include <stdlib.h>

#include <TStyle.h>
#include <TCanvas.h>
#include <TVirtualPad.h>  // access to gPad

#include "b2uQA.hh"
#include "mesFit.hh"

// Driverfile for the class in 
// ~ursl/macros/AnaClasses/b2uQA.{hh,cc}
//


// ----------------------------------------------------------------------
void runMerged(int mode = 0, int redo = 1, int Sys = 0, TString cutfile = "BLAH", int chains = 0, int reweight = 0, int areanorm = 0) {

  char line[1000];
  sprintf(line, "summary-%d.txt", mode);
  ofstream OUT(line);

  b2uQA *ad;
  if (redo == 1) {
    if (chains == 0) {
      ad = new b2uQA("csx-data.root", Sys, cutfile, 0);
      ad->makeAll(mode);
    } else {
      ad = new b2uQA("csx-data-*", Sys, cutfile, 0);
      ad->makeAll(mode);      
    }
  } else {
    sprintf(line, "qa-%d-csx-data.root", mode);
    ad = new b2uQA(line, Sys, cutfile, 0);
    ad->readHist();
  }

  b2uQA *a1; 
  if (redo == 1) {
    if (chains == 0) {
      a1 = new b2uQA("csx-genb-new.root", Sys, cutfile, reweight);
      a1->makeAll(mode);
    } else {
      a1 = new b2uQA("csx-genb*", Sys, cutfile, reweight);
      a1->makeAll(mode);
    }
  } else {
    sprintf(line, "qa-%d-csx-genb-new.root", mode);
    a1 = new b2uQA(line, Sys, cutfile, reweight);
    a1->readHist();
  }    


  b2uQA *a2=0; 
  if (redo == 1) {
    if (chains == 0) {
      a2 = new b2uQA("csx-cocktail-new.root", Sys, cutfile, reweight);
      a2->makeAll(mode);
    } else {
      a2 = new b2uQA("csx-b*cock*", Sys, cutfile, reweight);
      a2->makeAll(mode);
    }
  } else {
    sprintf(line, "qa-%d-csx-cocktail-new.root", mode);
    a2 = new b2uQA(line, Sys, cutfile, reweight);
    a2->readHist();
  }    

  if (mode%2 == 0) {
    ad->fillStrings("b#rightarrow u enh.", "Data", "Gen.  MC", "Cock. MC");
  } else {
    ad->fillStrings("b#rightarrow u dep.", "Data", "Gen.  MC", "Cock. MC");
  }


  // -- Normalization
  // ----------------
  TCanvas c0("c0");

  TH1D *hdMes = ad->fHmes;
  mesFit md(hdMes, "cb");
  double nhd = md.getSig();
  cout << "From mES fit, Data: " << nhd << endl;


  TH1D *h1Mes = a1->fHmes;
  mesFit m1(h1Mes, "cb");
  double nh1 = m1.getSig();
  cout << "From mES fit, MC 1: " << nh1 << endl;


  TH1D *h2Mes;
  double nh2(1.);
  if (a2 != 0) {
    h2Mes = a2->fHmes;
    mesFit m2(h2Mes, "cb");
    nh2 = m2.getSig();
  }
  cout << "From mES fit, MC 2: " << nh2 << endl;


  // ---------------------------------
  TH1D *hd = ad->fH1000->signalHist();
  //  nhd = hd->GetSumOfWeights();
  if (areanorm==0) {
    hd->Scale(1./nhd);
  } else {
    hd->Scale(1./(hd->Integral()));
  }
  //  hd->SetName("hd");

  TH1D *h1 = a1->fH1000->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  if (areanorm==0) {
    h1->Scale(1./nh1);
  } else {
    h1->Scale(1./(h1->Integral()));
  }
  //  h1->SetName("h1");

  TH1D *h2 = 0;
  if (a2 != 0) {
    h2 = a2->fH1000->signalHist();
    if (areanorm==0) {
      h2->Scale(1./nh2);
    } else {
      h2->Scale(1./(h2->Integral()));
    }
    //  h2->SetName("h2");
  }

  setTitles(hd, "p^{*} [GeV/c]", "");
  ad->show(hd, h1, h2);
  //   OUT << "pcms: " << ad->fNrProb << " " << chi2Prob(ad->fCChi2, ad->fCDof) 
  //       << " " << ad->fCfit << "+/-" << ad->fCfitE << " " << ad->fCChi2 << "/" << ad->fCDof << endl;
  OUT << "pcms: " 
      << ad->fNrProbA << " " << ad->fNrProbB << " " << ad->fKstA << " " << ad->fKstB
      << endl;

  
  sprintf(line, "h1000-%d.eps", mode);
  TCanvas *c6 = (TCanvas*)gROOT->FindObject("c6");
  if (c6) {
    c6->SaveAs(line);
  } else {
    cout << "No canvas found, no output ... " << endl;
  }




  // ---------------------------------
  hd = ad->fH900->signalHist();
  //  nhd = hd->GetSumOfWeights();
  if (areanorm==0) {
    hd->Scale(1./nhd);
  } else {
    hd->Scale(1./(hd->Integral()));
  }
  //  hd->SetName("hd");

  h1 = a1->fH900->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  if (areanorm==0) {
    h1->Scale(1./nh1);
  } else {
    h1->Scale(1./(h1->Integral()));
  }
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH900->signalHist();
    if (areanorm==0) {
      h2->Scale(1./nh2);
    } else {
      h2->Scale(1./(h2->Integral()));
    }
    //  h2->SetName("h2");
  }

  setTitles(hd, "pcms [GeV/c]", "");
  ad->show(hd, h1, h2);
  OUT << "pcmsloose: " 
      << ad->fNrProbA << " " << ad->fNrProbB << " " << ad->fKstA << " " << ad->fKstB
      << endl;
  sprintf(line, "h900-%d.eps", mode);
  c6 = (TCanvas*)gROOT->FindObject("c6");
  if (c6) c6->SaveAs(line);


  if (mode==0) {
    // ---------------------------------
    hd = ad->fH800->signalHist();
    //  nhd = hd->GetSumOfWeights();
    if (areanorm==0) {
      hd->Scale(1./nhd);
    } else {
      hd->Scale(1./(hd->Integral()));
    }
    //  hd->SetName("hd");

    h1 = a1->fH800->signalHist();
    //  nh1 = h1->GetSumOfWeights();
    if (areanorm==0) {
      h1->Scale(1./nh1);
    } else {
      h1->Scale(1./(h1->Integral()));
    }
    //  h1->SetName("h1");
    
    if (a2 != 0) {
      h2 = a2->fH800->signalHist();
      if (areanorm==0) {
	h2->Scale(1./nh2);
      } else {
	h2->Scale(1./(h2->Integral()));
      }
      //  h2->SetName("h2");
    }
    
    setTitles(hd, "pcms [GeV/c]", "");
    ad->show(hd, h1, h2);
    OUT << "pcmsvloose: " 
	<< ad->fNrProbA << " " << ad->fNrProbB << " " << ad->fKstA << " " << ad->fKstB
	<< endl;
    sprintf(line, "h800.eps");
    c6 = (TCanvas*)gROOT->FindObject("c6");
    if (c6) c6->SaveAs(line);
  }



  // ---------------------------------
  hd = ad->fH2000->signalHist();
  //  nhd = hd->GetSumOfWeights();
  if (areanorm==0) {
    hd->Scale(1./nhd);
  } else {
    hd->Scale(1./(hd->Integral()));
  }
  //  hd->SetName("hd");

  h1 = a1->fH2000->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  if (areanorm==0) {
    h1->Scale(1./nh1);
  } else {
    h1->Scale(1./(h1->Integral()));
  }
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH2000->signalHist();
    if (areanorm==0) {
      h2->Scale(1./nh2);
    } else {
      h2->Scale(1./(h2->Integral()));
    }
    //  h2->SetName("h2");
  }

  setTitles(hd, "m_{X} [GeV/c^{2}]", "");
  ad->show(hd, h1, h2);
  //   OUT << "mX:   " << ad->fNrProb << " " << chi2Prob(ad->fCChi2, ad->fCDof) 
  //       << " " << ad->fCfit << "+/-" << ad->fCfitE << " " << ad->fCChi2 << "/" << ad->fCDof << endl;
  OUT << "mX: " 
      << ad->fNrProbA << " " << ad->fNrProbB << " " << ad->fKstA << " " << ad->fKstB
      << endl;
  sprintf(line, "h2000-%d.eps", mode);
  c6 = (TCanvas*)gROOT->FindObject("c6");
  if (c6) c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH2100->signalHist();
  //  nhd = hd->GetSumOfWeights();
  if (areanorm==0) {
    hd->Scale(1./nhd);
  } else {
    hd->Scale(1./(hd->Integral()));
  }
  //  hd->SetName("hd");

  h1 = a1->fH2100->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  if (areanorm==0) {
    h1->Scale(1./nh1);
  } else {
    h1->Scale(1./(h1->Integral()));
  }
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH2100->signalHist();
    if (areanorm==0) {
      h2->Scale(1./nh2);
    } else {
      h2->Scale(1./(h2->Integral()));
    }
    //  h2->SetName("h2");
  }

  setTitles(hd, "m_{X}^{fit} [GeV/c^{2}]", "");
  ad->show(hd, h1, h2);
  //   OUT << "mXfit: " << ad->fNrProb << " " << chi2Prob(ad->fCChi2, ad->fCDof) 
  //       << " " << ad->fCfit << "+/-" << ad->fCfitE << " " << ad->fCChi2 << "/" << ad->fCDof << endl;
  OUT << "mXfit: " 
      << ad->fNrProbA << " " << ad->fNrProbB << " " << ad->fKstA << " " << ad->fKstB
      << endl;
  sprintf(line, "h2100-%d.eps", mode);
  c6 = (TCanvas*)gROOT->FindObject("c6");
  if (c6) c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH2200->signalHist();
  //  nhd = hd->GetSumOfWeights();
  if (areanorm==0) {
    hd->Scale(1./nhd);
  } else {
    hd->Scale(1./(hd->Integral()));
  }
  //  hd->SetName("hd");

  h1 = a1->fH2200->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  if (areanorm==0) {
    h1->Scale(1./nh1);
  } else {
    h1->Scale(1./(h1->Integral()));
  }
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH2200->signalHist();
    if (areanorm==0) {
      h2->Scale(1./nh2);
    } else {
      h2->Scale(1./(h2->Integral()));
    }
    //  h2->SetName("h2");
  }

  setTitles(hd, "m_{X}^{trk} [GeV/c^{2}]", "");
  ad->show(hd, h1, h2);
  //   OUT << "mXtrk: " << ad->fNrProb << " " << chi2Prob(ad->fCChi2, ad->fCDof) 
  //       << " " << ad->fCfit << "+/-" << ad->fCfitE << " " << ad->fCChi2 << "/" << ad->fCDof << endl;
  OUT << "mXtrk: " 
      << ad->fNrProbA << " " << ad->fNrProbB << " " << ad->fKstA << " " << ad->fKstB
      << endl;
  sprintf(line, "h2200-%d.eps", mode);
  c6 = (TCanvas*)gROOT->FindObject("c6");
  if (c6) c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH2300->signalHist();
  //  nhd = hd->GetSumOfWeights();
  if (areanorm==0) {
    hd->Scale(1./nhd);
  } else {
    hd->Scale(1./(hd->Integral()));
  }
  //  hd->SetName("hd");

  h1 = a1->fH2300->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  if (areanorm==0) {
    h1->Scale(1./nh1);
  } else {
    h1->Scale(1./(h1->Integral()));
  }
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH2300->signalHist();
    if (areanorm==0) {
      h2->Scale(1./nh2);
    } else {
      h2->Scale(1./(h2->Integral()));
    }
    //  h2->SetName("h2");
  }

  setTitles(hd, "m_{X}^{nut} [GeV/c^{2}]", "");
  ad->show(hd, h1, h2);
  //   OUT << "mXnut: " << ad->fNrProb << " " << chi2Prob(ad->fCChi2, ad->fCDof) 
  //       << " " << ad->fCfit << "+/-" << ad->fCfitE << " " << ad->fCChi2 << "/" << ad->fCDof << endl;
  OUT << "mXnut: " 
      << ad->fNrProbA << " " << ad->fNrProbB << " " << ad->fKstA << " " << ad->fKstB
      << endl;
  sprintf(line, "h2300-%d.eps", mode);
  c6 = (TCanvas*)gROOT->FindObject("c6");
  if (c6) c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH3000->signalHist();
  //  nhd = hd->GetSumOfWeights();
  if (areanorm==0) {
    hd->Scale(1./nhd);
  } else {
    hd->Scale(1./(hd->Integral()));
  }
  //  hd->SetName("hd");

  h1 = a1->fH3000->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  if (areanorm==0) {
    h1->Scale(1./nh1);
  } else {
    h1->Scale(1./(h1->Integral()));
  }
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH3000->signalHist();
    if (areanorm==0) {
      h2->Scale(1./nh2);
    } else {
      h2->Scale(1./(h2->Integral()));
    }
    //  h2->SetName("h2");
  }
  setTitles(hd, "m_{miss}^{2} [GeV/c^{2}]", "");
  ad->show(hd, h1, h2);
  //   OUT << "mm2: " << ad->fNrProb << " " << chi2Prob(ad->fCChi2, ad->fCDof) 
  //       << " " << ad->fCfit << "+/-" << ad->fCfitE << " " << ad->fCChi2 << "/" << ad->fCDof << endl;
  OUT << "mm2: " 
      << ad->fNrProbA << " " << ad->fNrProbB << " " << ad->fKstA << " " << ad->fKstB
      << endl;
  sprintf(line, "h3000-%d.eps", mode);
  c6 = (TCanvas*)gROOT->FindObject("c6");
  if (c6) c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH3100->signalHist();
  //  nhd = hd->GetSumOfWeights();
  if (areanorm==0) {
    hd->Scale(1./nhd);
  } else {
    hd->Scale(1./(hd->Integral()));
  }
  //  hd->SetName("hd");

  h1 = a1->fH3100->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  if (areanorm==0) {
    h1->Scale(1./nh1);
  } else {
    h1->Scale(1./(h1->Integral()));
  }
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH3100->signalHist();
    if (areanorm==0) {
      h2->Scale(1./nh2);
    } else {
      h2->Scale(1./(h2->Integral()));
    }
    //  h2->SetName("h2");
  }
  setTitles(hd, "p_{miss} [GeV/c]", "");
  ad->show(hd, h1, h2);
  //   OUT << "pmiss: " << ad->fNrProb << " " << chi2Prob(ad->fCChi2, ad->fCDof) 
  //       << " " << ad->fCfit << "+/-" << ad->fCfitE << " " << ad->fCChi2 << "/" << ad->fCDof << endl;
  OUT << "pmiss: " 
      << ad->fNrProbA << " " << ad->fNrProbB << " " << ad->fKstA << " " << ad->fKstB
      << endl;
  sprintf(line, "h3100-%d.eps", mode);
  c6 = (TCanvas*)gROOT->FindObject("c6");
  if (c6) c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH3200->signalHist();
  //  nhd = hd->GetSumOfWeights();
  if (areanorm==0) {
    hd->Scale(1./nhd);
  } else {
    hd->Scale(1./(hd->Integral()));
  }
  //  hd->SetName("hd");

  h1 = a1->fH3200->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  if (areanorm==0) {
    h1->Scale(1./nh1);
  } else {
    h1->Scale(1./(h1->Integral()));
  }
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH3200->signalHist();
    if (areanorm==0) {
      h2->Scale(1./nh2);
    } else {
      h2->Scale(1./(h2->Integral()));
    }
    //  h2->SetName("h2");
  }

  setTitles(hd, "#theta_{miss} [rad]", "");
  ad->show(hd, h1, h2);
  //   OUT << "tmiss: " << ad->fNrProb << " " << chi2Prob(ad->fCChi2, ad->fCDof) 
  //       << " " << ad->fCfit << "+/-" << ad->fCfitE << " " << ad->fCChi2 << "/" << ad->fCDof << endl;
  OUT << "tmiss: " 
      << ad->fNrProbA << " " << ad->fNrProbB << " " << ad->fKstA << " " << ad->fKstB
      << endl;
  sprintf(line, "h3200-%d.eps", mode);
  c6 = (TCanvas*)gROOT->FindObject("c6");
  if (c6) c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH3300->signalHist();
  //  nhd = hd->GetSumOfWeights();
  if (areanorm==0) {
    hd->Scale(1./nhd);
  } else {
    hd->Scale(1./(hd->Integral()));
  }
  //  hd->SetName("hd");

  h1 = a1->fH3300->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  if (areanorm==0) {
    h1->Scale(1./nh1);
  } else {
    h1->Scale(1./(h1->Integral()));
  }
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH3300->signalHist();
    if (areanorm==0) {
      h2->Scale(1./nh2);
    } else {
      h2->Scale(1./(h2->Integral()));
    }
    //  h2->SetName("h2");
  }
  setTitles(hd, "Q^{2} [GeV^{2}/c^{4}]", "");
  ad->show(hd, h1, h2);
  //   OUT << "Q2: " << ad->fNrProb << " " << chi2Prob(ad->fCChi2, ad->fCDof) 
  //       << " " << ad->fCfit << "+/-" << ad->fCfitE << " " << ad->fCChi2 << "/" << ad->fCDof << endl;
  OUT << "Q2: " 
      << ad->fNrProbA << " " << ad->fNrProbB << " " << ad->fKstA << " " << ad->fKstB
      << endl;
  sprintf(line, "h3300-%d.eps", mode);
  c6 = (TCanvas*)gROOT->FindObject("c6");
  if (c6) c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH3400->signalHist();
  //  nhd = hd->GetSumOfWeights();
  if (areanorm==0) {
    hd->Scale(1./nhd);
  } else {
    hd->Scale(1./(hd->Integral()));
  }
  //  hd->SetName("hd");

  h1 = a1->fH3400->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  if (areanorm==0) {
    h1->Scale(1./nh1);
  } else {
    h1->Scale(1./(h1->Integral()));
  }
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH3400->signalHist();
    if (areanorm==0) {
      h2->Scale(1./nh2);
    } else {
      h2->Scale(1./(h2->Integral()));
    }
    //  h2->SetName("h2");
  }
  setTitles(hd, "E_{miss} - p_{miss} [GeV]", "");
  ad->show(hd, h1, h2);
  //   OUT << "emiss-pmiss: " << ad->fNrProb << " " << chi2Prob(ad->fCChi2, ad->fCDof) 
  //       << " " << ad->fCfit << "+/-" << ad->fCfitE << " " << ad->fCChi2 << "/" << ad->fCDof << endl;
  OUT << "emiss-pmiss: " 
      << ad->fNrProbA << " " << ad->fNrProbB << " " << ad->fKstA << " " << ad->fKstB
      << endl;
  sprintf(line, "h3400-%d.eps", mode);
  c6 = (TCanvas*)gROOT->FindObject("c6");
  if (c6) c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH3500->signalHist();
  //  nhd = hd->GetSumOfWeights();
  if (areanorm==0) {
    hd->Scale(1./nhd);
  } else {
    hd->Scale(1./(hd->Integral()));
  }
  //  hd->SetName("hd");

  h1 = a1->fH3500->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  if (areanorm==0) {
    h1->Scale(1./nh1);
  } else {
    h1->Scale(1./(h1->Integral()));
  }
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH3500->signalHist();
    if (areanorm==0) {
      h2->Scale(1./nh2);
    } else {
      h2->Scale(1./(h2->Integral()));
    }
    //  h2->SetName("h2");
  }
  setTitles(hd, "cos(#theta_{miss})", "");
  ad->show(hd, h1, h2);
  //   OUT << "cos(tmiss): " << ad->fNrProb << " " << chi2Prob(ad->fCChi2, ad->fCDof) 
  //       << " " << ad->fCfit << "+/-" << ad->fCfitE << " " << ad->fCChi2 << "/" << ad->fCDof << endl;
  OUT << "cos(tmiss): " 
      << ad->fNrProbA << " " << ad->fNrProbB << " " << ad->fKstA << " " << ad->fKstB
      << endl;
  sprintf(line, "h3500-%d.eps", mode);
  c6 = (TCanvas*)gROOT->FindObject("c6");
  if (c6) c6->SaveAs(line);



  // ---------------------------------
  hd = ad->fH4000->signalHist();
  //  nhd = hd->GetSumOfWeights();
  if (areanorm==0) {
    hd->Scale(1./nhd);
  } else {
    hd->Scale(1./(hd->Integral()));
  }
  //  hd->SetName("hd");

  h1 = a1->fH4000->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  if (areanorm==0) {
    h1->Scale(1./nh1);
  } else {
    h1->Scale(1./(h1->Integral()));
  }
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH4000->signalHist();
    if (areanorm==0) {
      h2->Scale(1./nh2);
    } else {
      h2->Scale(1./(h2->Integral()));
    }
    //  h2->SetName("h2");
  }
  setTitles(hd, "N_{TRK}", "");
  ad->show(hd, h1, h2);
  //   OUT << "Ntrk: " << ad->fNrProb << " " << chi2Prob(ad->fCChi2, ad->fCDof) 
  //       << " " << ad->fCfit << "+/-" << ad->fCfitE << " " << ad->fCChi2 << "/" << ad->fCDof << endl;
  OUT << "Ntrk: " 
      << ad->fNrProbA << " " << ad->fNrProbB << " " << ad->fKstA << " " << ad->fKstB
      << endl;
  sprintf(line, "h4000-%d.eps", mode);
  c6 = (TCanvas*)gROOT->FindObject("c6");
  if (c6) c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH4100->signalHist();
  //  nhd = hd->GetSumOfWeights();
  if (areanorm==0) {
    hd->Scale(1./nhd);
  } else {
    hd->Scale(1./(hd->Integral()));
  }
  //  hd->SetName("hd");

  h1 = a1->fH4100->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  if (areanorm==0) {
    h1->Scale(1./nh1);
  } else {
    h1->Scale(1./(h1->Integral()));
  }
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH4100->signalHist();
    if (areanorm==0) {
      h2->Scale(1./nh2);
    } else {
      h2->Scale(1./(h2->Integral()));
    }
    //  h2->SetName("h2");
  }
  setTitles(hd, "N_{Neutrals}", "");
  ad->show(hd, h1, h2);
  //   OUT << "Nnut: " << ad->fNrProb << " " << chi2Prob(ad->fCChi2, ad->fCDof) 
  //       << " " << ad->fCfit << "+/-" << ad->fCfitE << " " << ad->fCChi2 << "/" << ad->fCDof << endl;
  OUT << "Nnut: " 
      << ad->fNrProbA << " " << ad->fNrProbB << " " << ad->fKstA << " " << ad->fKstB
      << endl;
  sprintf(line, "h4100-%d.eps", mode);
  c6 = (TCanvas*)gROOT->FindObject("c6");
  if (c6) c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH4200->signalHist();
  //  nhd = hd->GetSumOfWeights();
  if (areanorm==0) {
    hd->Scale(1./nhd);
  } else {
    hd->Scale(1./(hd->Integral()));
  }
  //  hd->SetName("hd");

  h1 = a1->fH4200->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  if (areanorm==0) {
    h1->Scale(1./nh1);
  } else {
    h1->Scale(1./(h1->Integral()));
  }
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH4200->signalHist();
    if (areanorm==0) {
      h2->Scale(1./nh2);
    } else {
      h2->Scale(1./(h2->Integral()));
    }
    //  h2->SetName("h2");
  }
  setTitles(hd, "Q_{Total}", "");
  ad->show(hd, h1, h2);
  //   OUT << "qtot: " << ad->fNrProb << " " << chi2Prob(ad->fCChi2, ad->fCDof) 
  //       << " " << ad->fCfit << "+/-" << ad->fCfitE << " " << ad->fCChi2 << "/" << ad->fCDof << endl;
  OUT << "qtot: " 
      << ad->fNrProbA << " " << ad->fNrProbB << " " << ad->fKstA << " " << ad->fKstB
      << endl;
  sprintf(line, "h4200-%d.eps", mode);
  c6 = (TCanvas*)gROOT->FindObject("c6");
  if (c6) c6->SaveAs(line);


  delete ad;
  delete a1;
  if (a2 != 0) {
    delete a2;
  }


//   TFile f1("histograms.root", "RECREATE");
//   hd->SetDirectory(gDirectory);
//   h1->SetDirectory(gDirectory);
//   h2->SetDirectory(gDirectory);
//   f1.Write();
//   f1.Close();
  
}



// ----------------------------------------------------------------------
int main(int argc, char *argv[]) {

  gROOT->SetBatch(kTRUE);
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111111);  // Show overflow, underflow + SumOfWeights 
  gStyle->SetStatStyle(0);      // for a completely transparent stat box
  gStyle->SetOptFit(111110); 
  gStyle->SetOptFile(1); 
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(.8);
  gStyle->SetMarkerColor(1);
  gStyle->SetNdivisions(505, "x");
  gStyle->SetTickLength(-0.02, "x");
  gStyle->SetTickLength(-0.02, "y");
  gStyle->SetLabelOffset(0.02, "x");
  gStyle->SetLabelOffset(0.02, "y");
  gStyle->SetPadLeftMargin(0.15); 
  gStyle->SetPadBottomMargin(015); 
  gStyle->SetTitleBorderSize(0);
  gStyle->SetStatFont(132); 
  gStyle->SetTextFont(132); 
  gStyle->SetLabelFont(132, "X"); 
  gStyle->SetLabelFont(132, "Y"); 
  gStyle->SetTitleFont(132); 

  gROOT->ForceStyle();

  int all = 0;
  int redo = 1;
  int Sys = 0;
  TString cutFile("BLAH");
  int chains = 0;
  int reweight = 0;
  int areanorm = 0;

  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if(strcmp(argv[i],"-a")  == 0)  {all = 1; }
    if(strcmp(argv[i],"-r")  == 0)  {redo = atoi(argv[++i]); }
    if(strcmp(argv[i],"-s")  == 0)  {Sys = atoi(argv[++i]); }
    if(strcmp(argv[i],"-C")  == 0)  {cutFile    = TString(argv[++i]); }
    if(strcmp(argv[i],"-c")  == 0)  {chains = 1; }
    if(strcmp(argv[i],"-rw")  == 0) {reweight = 1; }
    if(strcmp(argv[i],"-an")  == 0) {areanorm = 1; }
  }


  //   // -- testing with enhanced only
  //   runMerged(0, redo);
  //   return 0;



  // -- Enhanced all 
  runMerged(0, redo, Sys, cutFile, chains, reweight, areanorm);
  // -- Depleted all 
  runMerged(1, redo, Sys, cutFile, chains, reweight, areanorm);

  if (all == 0) {
    return 0;
  }

  // -- Electrons
  runMerged(20, redo, Sys, cutFile, chains, reweight, areanorm);
  runMerged(21, redo, Sys, cutFile, chains, reweight, areanorm);

  // -- Muons
  runMerged(30, redo, Sys, cutFile, chains, reweight, areanorm);
  runMerged(31, redo, Sys, cutFile, chains, reweight, areanorm);

  // -- B0
  runMerged(200, redo, Sys, cutFile, chains, reweight, areanorm);
  runMerged(201, redo, Sys, cutFile, chains, reweight, areanorm);

  // -- B+
  runMerged(300, redo, Sys, cutFile, chains, reweight, areanorm);
  runMerged(301, redo, Sys, cutFile, chains, reweight, areanorm);

  return 0;
}
