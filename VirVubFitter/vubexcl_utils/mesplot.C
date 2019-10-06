void mesplot(){

  //  TFile g1("/nfs/farm/babar/AWGsemilep01/gallof/fit/fiteta/data/testdatafitresult.root");

  char name[100];
  int print = 0;

  c1 = new TCanvas("c1"," ",200,10,800,640);
  c1->Clear();

  nsldata->SetTitle("");
  nsldata->SetXTitle("m_{ES}[GeV/c^{2}]");
  nsldata->SetYTitle("Events/(2.5 MeV/c^{2})");
  nsldata->SetTitleOffset(1.);
  nsldata->SetTitleSize(.07);

  recoilAnalysis k;
  //  double thesigma = 0.003;
  double thesigma = -1111111;
  double themean = -1111111;
  double theargus = -1111111;
  double thealpha = -1111111;
  double then = -1111111;
  double resmean, ressigma, resalpha, resn;
  k.vubMes(&nsldata, resmean, ressigma, resalpha, resn, print, 1, themean, thesigma, thealpha, 5., theargus);
  
  double xpos = 0.35;
  double ypos = 0.83;
  double scale = 1; 
  double scale2 = .5; 
  TString str="preliminary";
  
  // Draw the label 
  TLatex *babar = new TLatex();
  Double_t cheburashkaFactorX=1, cheburashkaFactorY=1, padSizeX=500, padSizeY=500, xpos2, ypos2, xposL;
  babar->SetNDC(kTRUE);
  babar->SetTextFont(32); // Bold-Italic Times
  babar->SetTextAlign(31); // Right-Bottom
  padSizeX = gPad->GetWw()*gPad->GetWNDC(); // Get pad's dimensions
  padSizeY = gPad->GetWh()*gPad->GetHNDC();
  if (padSizeX>padSizeY) cheburashkaFactorX=padSizeY/padSizeX;
  if (padSizeX<padSizeY) cheburashkaFactorY=padSizeX/padSizeY;
  //xpos2=xpos-0.185*scale*cheburashkaFactorX;
  xpos2=xpos-0.188*scale*cheburashkaFactorX;
  ypos2=ypos-0.0620*scale*cheburashkaFactorY;
  xposL=xpos-0.253*scale*cheburashkaFactorX;
  babar->SetTextSize(0.10*scale); // Beginning to draw "BaBar"
  babar->DrawText(xpos2,ypos2,"B");
  babar->SetTextSize(0.075*scale);
  babar->DrawText(xpos2+0.039*scale*cheburashkaFactorX,ypos2,"A");
  babar->SetTextSize(0.10*scale);
  babar->DrawText(xpos2+0.1015*scale*cheburashkaFactorX,ypos2,"B");
  babar->SetTextSize(0.075*scale);
  babar->DrawText(xpos2+0.1875*scale*cheburashkaFactorX,ypos2,"AR");
  babar->SetTextFont(42); // Helvetica (medium, upright) 
  babar->SetTextSize(0.1*scale2);
  babar->SetTextAlign(33); // Right-Top
  babar->DrawLatex(xpos,ypos2-0.02*scale2*cheburashkaFactorY,str);

  c1->SaveAs("nsldata.eps");

  c1->Clear();

  allcutsdata->SetTitle("");
  allcutsdata->SetXTitle("m_{ES}[GeV/c^{2}]");
  allcutsdata->SetTitleOffset(1.);
  allcutsdata->SetTitleSize(.07);

  recoilAnalysis k;
  //  double thesigma = 0.003;
  double thesigma = -1111111;
  double themean = -1111111;
  double theargus = -1111111;
  double thealpha = -1111111;
  double then = -1111111;
  double resmean2, ressigma2, resalpha2, resn2;
  k.vubMes(&allcutsdata, resmean2, ressigma2, resalpha2, resn2, print, 1, resmean, ressigma, resalpha, resn, theargus);
  
  double xpos = 0.35;
  double ypos = 0.83;
  double scale = 1; 
  double scale2 = .5; 
  TString str="preliminary";
  
  // Draw the label 
  TLatex *babar = new TLatex();
  Double_t cheburashkaFactorX=1, cheburashkaFactorY=1, padSizeX=500, padSizeY=500, xpos2, ypos2, xposL;
  babar->SetNDC(kTRUE);
  babar->SetTextFont(32); // Bold-Italic Times
  babar->SetTextAlign(31); // Right-Bottom
  padSizeX = gPad->GetWw()*gPad->GetWNDC(); // Get pad's dimensions
  padSizeY = gPad->GetWh()*gPad->GetHNDC();
  if (padSizeX>padSizeY) cheburashkaFactorX=padSizeY/padSizeX;
  if (padSizeX<padSizeY) cheburashkaFactorY=padSizeX/padSizeY;
  //xpos2=xpos-0.185*scale*cheburashkaFactorX;
  xpos2=xpos-0.188*scale*cheburashkaFactorX;
  ypos2=ypos-0.0620*scale*cheburashkaFactorY;
  xposL=xpos-0.253*scale*cheburashkaFactorX;
  babar->SetTextSize(0.10*scale); // Beginning to draw "BaBar"
  babar->DrawText(xpos2,ypos2,"B");
  babar->SetTextSize(0.075*scale);
  babar->DrawText(xpos2+0.039*scale*cheburashkaFactorX,ypos2,"A");
  babar->SetTextSize(0.10*scale);
  babar->DrawText(xpos2+0.1015*scale*cheburashkaFactorX,ypos2,"B");
  babar->SetTextSize(0.075*scale);
  babar->DrawText(xpos2+0.1875*scale*cheburashkaFactorX,ypos2,"AR");
  babar->SetTextFont(42); // Helvetica (medium, upright) 
  babar->SetTextSize(0.1*scale2);
  babar->SetTextAlign(33); // Right-Top
  babar->DrawLatex(xpos,ypos2-0.02*scale2*cheburashkaFactorY,str);

  c1->SaveAs("allcutsdata.eps");

}
