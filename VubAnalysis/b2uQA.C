
// Driverfile for the class in 
// b2uQA.{hh,cc}
//
// Do on shire or on noric:
//
//  root [0] gSystem.Load(".../shlib/SunOS58/libVubAnalysis.so");
//  root [0] .L .../b2uQA.C  
//  root [2] runMerged() 
//
//



// ----------------------------------------------------------------------
void makeAll(int all = 1, int redo = 1) {
  // -- Enhanced all 
  runMerged(0, redo);
  // -- Depleted all 
  runMerged(1, redo);

  if (all == 0) return;

  // -- Electrons
  runMerged(20, redo);
  runMerged(21, redo);

  // -- Muons
  runMerged(30, redo);
  runMerged(31, redo);

  // -- B0
  runMerged(200, redo);
  runMerged(201, redo);

  // -- B+
  runMerged(300, redo);
  runMerged(301, redo);
}


// ----------------------------------------------------------------------
void runMerged(int mode = 0, int redo = 1) {

  char line[1000];
  sprintf(line, "summary-%d.txt", mode);
  ofstream OUT(line);

  b2uQA *ad;
  if (redo == 1) {
    ad = new b2uQA("csx-data*");
    ad->makeAll(mode);
  } else {
    sprintf(line, "qa-%d-csx-data.root", mode);
    ad = new b2uQA(line);
    ad->readHist();
  }

  b2uQA *a1; 
  if (redo == 1) {
    a1 = new b2uQA("csx-genb*");
    a1->makeAll(mode);
  } else {
    sprintf(line, "qa-%d-csx-genb-new.root", mode);
    a1 = new b2uQA(line);
    a1->readHist();
  }    

  b2uQA *a2=0; 
  if (redo == 1) {
    a2 = new b2uQA("csx-b*");
    a2->makeAll(mode);
  } else {
    sprintf(line, "qa-%d-csx-cocktail-new.root", mode);
    a2 = new b2uQA(line);
    a2->readHist();
  }    

  if (mode%2 == 0) {
    ad->fillStrings("b#rightarrow u enh.", "Data", "Gen.  MC", "Cock. MC");
  } else {
    ad->fillStrings("b#rightarrow u dep.", "Data", "Gen.  MC", "Cock. MC");
  }


  // -- Normalization
  // ----------------
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
  hd->Scale(1./nhd);
  //  hd->SetName("hd");

  TH1D *h1 = a1->fH1000->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  h1->Scale(1./nh1);
  //  h1->SetName("h1");

  TH1D *h2 = 0;
  if (a2 != 0) {
    h2 = a2->fH1000->signalHist();
    h2->Scale(1./nh2);
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
  c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH2000->signalHist();
  //  nhd = hd->GetSumOfWeights();
  hd->Scale(1./nhd);
  //  hd->SetName("hd");

  h1 = a1->fH2000->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  h1->Scale(1./nh1);
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH2000->signalHist();
    h2->Scale(1./nh2);
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
  c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH2100->signalHist();
  //  nhd = hd->GetSumOfWeights();
  hd->Scale(1./nhd);
  //  hd->SetName("hd");

  h1 = a1->fH2100->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  h1->Scale(1./nh1);
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH2100->signalHist();
    h2->Scale(1./nh2);
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
  c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH2200->signalHist();
  //  nhd = hd->GetSumOfWeights();
  hd->Scale(1./nhd);
  //  hd->SetName("hd");

  h1 = a1->fH2200->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  h1->Scale(1./nh1);
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH2200->signalHist();
    h2->Scale(1./nh2);
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
  c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH2300->signalHist();
  //  nhd = hd->GetSumOfWeights();
  hd->Scale(1./nhd);
  //  hd->SetName("hd");

  h1 = a1->fH2300->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  h1->Scale(1./nh1);
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH2300->signalHist();
    h2->Scale(1./nh2);
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
  c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH3000->signalHist();
  //  nhd = hd->GetSumOfWeights();
  hd->Scale(1./nhd);
  //  hd->SetName("hd");

  h1 = a1->fH3000->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  h1->Scale(1./nh1);
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH3000->signalHist();
    h2->Scale(1./nh2);
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
  c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH3100->signalHist();
  //  nhd = hd->GetSumOfWeights();
  hd->Scale(1./nhd);
  //  hd->SetName("hd");

  h1 = a1->fH3100->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  h1->Scale(1./nh1);
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH3100->signalHist();
    h2->Scale(1./nh2);
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
  c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH3200->signalHist();
  //  nhd = hd->GetSumOfWeights();
  hd->Scale(1./nhd);
  //  hd->SetName("hd");

  h1 = a1->fH3200->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  h1->Scale(1./nh1);
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH3200->signalHist();
    h2->Scale(1./nh2);
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
  c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH3300->signalHist();
  //  nhd = hd->GetSumOfWeights();
  hd->Scale(1./nhd);
  //  hd->SetName("hd");

  h1 = a1->fH3300->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  h1->Scale(1./nh1);
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH3300->signalHist();
    h2->Scale(1./nh2);
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
  c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH3400->signalHist();
  //  nhd = hd->GetSumOfWeights();
  hd->Scale(1./nhd);
  //  hd->SetName("hd");

  h1 = a1->fH3400->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  h1->Scale(1./nh1);
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH3400->signalHist();
    h2->Scale(1./nh2);
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
  c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH3500->signalHist();
  //  nhd = hd->GetSumOfWeights();
  hd->Scale(1./nhd);
  //  hd->SetName("hd");

  h1 = a1->fH3500->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  h1->Scale(1./nh1);
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH3500->signalHist();
    h2->Scale(1./nh2);
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
  c6->SaveAs(line);



  // ---------------------------------
  hd = ad->fH4000->signalHist();
  //  nhd = hd->GetSumOfWeights();
  hd->Scale(1./nhd);
  //  hd->SetName("hd");

  h1 = a1->fH4000->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  h1->Scale(1./nh1);
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH4000->signalHist();
    h2->Scale(1./nh2);
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
  c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH4100->signalHist();
  //  nhd = hd->GetSumOfWeights();
  hd->Scale(1./nhd);
  //  hd->SetName("hd");

  h1 = a1->fH4100->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  h1->Scale(1./nh1);
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH4100->signalHist();
    h2->Scale(1./nh2);
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
  c6->SaveAs(line);


  // ---------------------------------
  hd = ad->fH4200->signalHist();
  //  nhd = hd->GetSumOfWeights();
  hd->Scale(1./nhd);
  //  hd->SetName("hd");

  h1 = a1->fH4200->signalHist();
  //  nh1 = h1->GetSumOfWeights();
  h1->Scale(1./nh1);
  //  h1->SetName("h1");

  if (a2 != 0) {
    h2 = a2->fH4200->signalHist();
    h2->Scale(1./nh2);
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
  c6->SaveAs(line);


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



