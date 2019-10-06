void BABARSmartLabel(Double_t xpos=0.9, Double_t ypos=0.9, Double_t scale=1.0, TString str="null", Double_t scale2=0.5, TString align="R") {
  // Making -1 a placeholder for function's default value 
  if (xpos == -1) xpos = 0.9;
  if (ypos == -1) ypos = 0.9;
  if (scale == -1) scale = 1; 
  if (str == "-1") str = "null"; 
  if (scale2 == -1) scale2 = 0.5; 
  if (align == "-1") align = "R"; 
  // A few predefined labels to go to the second line of text
  if (str == "~1") str = "preliminary"; 
  if (str == "~2") str = "very preliminary"; 
  if (str == "~2000") str = "year 2000 preliminary"; 
  if (str == "~2001") str = "year 2000 preliminary"; 
  if (str == "~25") str = "25 fb^{-1}";
  if (str == "~B->fcp") str = "#font[12]{#Gamma#font[132]{(}B^{#font[12]{0}}_{#font[132]{phys}}#font[132]{(}t#font[132]{)} #rightarrow f_{CP}#font[132]{)} = #left|A_{f_{CP}}#right|^{#font[132]{2}} e^{-#Gamma t} #left[#frac{1+|#lambda_{f_{CP}}|^{#font[132]{2}}}{#font[132]{2}} + #frac{#font[132]{1}-|#lambda_{f_{CP}}|^{#font[132]{2}}}{#font[132]{2}} #font[132]{cos}#font[132]{(}#DeltaMt#font[132]{)} - #font[132]{Im }#lambda_{f_{CP}}#font[132]{sin}#font[132]{(}#DeltaMt#font[132]{)}#right]}";

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
  if (align == "L")  {
    babar->SetTextAlign(13); // Left-Top
    babar->DrawLatex(xposL,ypos2-0.02*scale2*cheburashkaFactorY,str);
  }
  else {
    babar->SetTextAlign(33); // Right-Top
    babar->DrawLatex(xpos,ypos2-0.02*scale2*cheburashkaFactorY,str);
  }
  delete babar;
}

void diffBR(double valq21, double errvalq21, double valq22, double errvalq22, double valq23, double errvalq23, char *title){

  TCanvas c1("c1","--c1--",472,0,800,600);
  

  TH1D diffBR("diffBR","differential BR",3,0.,24.);
 
  diffBR.SetBinContent(1,valq21);
  diffBR.SetBinError(1,errvalq21);
  diffBR.SetBinContent(2,valq22);
  diffBR.SetBinError(2,errvalq22);
  diffBR.SetBinContent(3,valq23);
  diffBR.SetBinError(3,errvalq23);
  diffBR.SetLabelOffset(.01,"Y");
  diffBR.SetTitleOffset(.95,"Y");
  diffBR.SetTitleSize(0.055,"Y");
  diffBR.SetTitleOffset(1.05,"X");
  diffBR.SetTitleSize(0.055,"X");
  diffBR.SetTitle("");
  diffBR.SetYTitle(title);
  diffBR.SetXTitle("q^{2}[GeV^{2}]");
  diffBR.SetBarWidth(.5);
  //  diffBR.SetNdivisions(114,"X");

  diffBR.SetStats(0);
  diffBR.SetMarkerSize(1.4);

  diffBR.Draw("e1p");

//   TGaxis tempaxis1(0.,0.01,30.,0.01,0.,30.,103);
//   tempaxis1.Draw();

  BABARSmartLabel(0.86,0.87,1,"preliminary");
  
  c1.SaveAs("test.eps");

}
