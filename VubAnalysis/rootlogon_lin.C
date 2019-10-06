{

  //  gROOT->Reset();

  
 gSystem.Load("libPhysics.so");
 gSystem.Load("../shlib/$BFARCH/libRecoilAnalysis.so");
 gSystem.Load("../shlib/$BFARCH/libVubAnalysis.so");

  //gSystem.Load("/u/ec/ursl/macros/lib/libEmcUtil-30207.so");
  //gSystem.Load("/u/ec/ursl/macros/lib/libFsxUtil-30207.so");

  gROOT->SetStyle("Plain");

  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111111);  // Show overflow, underflow + SumOfWeights 
  //  gStyle->SetStatStyle(0); // for a completely transparent stat box
  gStyle->SetOptFit(111110); 
  gStyle->SetOptFile(1); 

  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(.8);
  gStyle->SetMarkerColor(1);

  gStyle->SetTitleBorderSize(0);  // no border around histogram title (font size can't be calculated anyways ...)

  gROOT->ForceStyle();

  // --- Cleanup if this is not the first call to rootlogon.C
  TCanvas *c = 0;
  c = (TCanvas*)gROOT->FindObject("c0"); if (c) c->Delete(); c = 0;
  p = (TPad*)gROOT->FindObject("p0"); if (p) p->Delete(); p = 0;
  // --- Create a new canvas.
  TCanvas c0("c0","--c0--",472,0,800,900);
  //  TCanvas c0("c0","--c0--",356,0,656,700);
  c0->ToggleEventStatus();

  TLatex tl;
  tl.SetNDC(kTRUE);


}



