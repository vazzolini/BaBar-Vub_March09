

void makefinalplots(){

  //  TFile g("/u/br/dorazioa/scra/newfit/pidefault/pilnudatafitresult.root");
  c1 = new TCanvas("c1"," ",200,10,800,800);
  c1->Clear();
  char namevar[100],namepar[100],nameuni[100],nameaxis[100];
  char name[100], name2[100];
  char var[100];
  sprintf(namevar, "m_{miss}^{2}");
  sprintf(nameuni, "[GeV^{2}/c^{4}]");
  //  sprintf(namepar, "%s%s",namevar,"(#pi^{+})");
  sprintf(namepar, "%s%s",namevar,"(#eta)");
  sprintf(nameaxis, "%s%s",namepar,nameuni);
  sprintf(name, "%s%s",var,"dataexcl");
  sprintf(var,"mm2");
  sprintf(name, "%s%s",var,"dataexcl");
  double max(0.);
  for (int i=6; i<27; i++){
    double tempmax = ((TH1D*)gDirectory->Get(name))->GetBinContent(i);
    if(tempmax>max) max = tempmax;
  }
  ((TH1D*)gDirectory->Get(name))->SetTitle("");
  ((TH1D*)gDirectory->Get(name))->SetMaximum(max*1.2);
  ((TH1D*)gDirectory->Get(name))->SetMarkerStyle(8);
  ((TH1D*)gDirectory->Get(name))->SetMarkerSize(1.);
  ((TH1D*)gDirectory->Get(name))->SetLabelSize(.035,"X");
  ((TH1D*)gDirectory->Get(name))->SetLabelSize(.035,"Y");
  ((TH1D*)gDirectory->Get(name))->SetMarkerColor(kBlack);
  ((TH1D*)gDirectory->Get(name))->SetXTitle(nameaxis);
  ((TH1D*)gDirectory->Get(name))->SetMinimum(0.);
  ((TH1D*)gDirectory->Get(name))->SetAxisRange(-1.625,3.874);
  ((TH1D*)gDirectory->Get(name))->SetStats(0);
  ((TH1D*)gDirectory->Get(name))->Draw();
  sprintf(name, "%s%s",var,"allmcexcl");
  ((TH1D*)gDirectory->Get(name))->SetLineWidth(3.);
  ((TH1D*)gDirectory->Get(name))->SetLineColor(kBlack);
  ((TH1D*)gDirectory->Get(name))->SetFillColor(kWhite);
  ((TH1D*)gDirectory->Get(name))->SetStats(0);
  ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
  sprintf(name, "%s%s",var,"allbkgexcl");
  ((TH1D*)gDirectory->Get(name))->SetFillColor(18);
  ((TH1D*)gDirectory->Get(name))->SetStats(0);
  ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
  sprintf(name, "%s%s",var,"nonvubexcl");
  ((TH1D*)gDirectory->Get(name))->SetFillColor(16);
  ((TH1D*)gDirectory->Get(name))->SetStats(0);
  ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
  sprintf(name, "%s%s",var,"scaleothexcl");
  ((TH1D*)gDirectory->Get(name))->SetFillColor(13);
  ((TH1D*)gDirectory->Get(name))->SetStats(0);
  ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
  sprintf(name, "%s%s",var,"dataexcl");
  ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
  TLegendEntry *legge;
  TLegend *leg;
  leg = new TLegend(0.13,0.45,0.35,0.85);
  leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  sprintf(name, "%s%s",var,"allmcexcl");
  legge = leg->AddEntry(((TH1D*)gDirectory->Get(name)), "signal", "f");
  sprintf(name, "%s%s",var,"allbkgexcl");
  legge = leg->AddEntry(((TH1D*)gDirectory->Get(name)), "b#rightarrowul#nu", "f");
  sprintf(name, "%s%s",var,"nonvubexcl");
  legge = leg->AddEntry(((TH1D*)gDirectory->Get(name)), "b#rightarrowcl#nu", "f");
  sprintf(name, "%s%s",var,"scaleothexcl");
  legge = leg->AddEntry(((TH1D*)gDirectory->Get(name)), "other", "f");
  sprintf(name, "%s%s",var,"dataexcl");
  legge = leg->AddEntry(((TH1D*)gDirectory->Get(name)), "data", "p");
  leg->Draw();

  TGaxis tempaxis1(3.874,0.,3.874,int(max*1.2),0.,int(max*1.2),510,"+U");
  tempaxis1.Draw();

  double xpos = 0.75;
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

  cout << xpos << " " << ypos2-0.02 << " " << scale << " " << cheburashkaFactorY << " " << str << endl;
  
  c1->SaveAs("mm2.eps");
  
  sprintf(namevar, "m");
  sprintf(nameuni, "[GeV/c^{2}]");
  sprintf(namepar, "%s%s",namevar,"(#pi^{+})");
  sprintf(nameaxis, "%s%s",namepar,nameuni);
  sprintf(name, "%s%s",var,"dataexcl");
  sprintf(var,"mass");
  sprintf(name, "%s%s",var,"dataexcl");
  double max(0.);
  for (int i=1; i<11; i++){
    double tempmax = ((TH1D*)gDirectory->Get(name))->GetBinContent(i);
    if(tempmax>max) max = tempmax;
  }
  ((TH1D*)gDirectory->Get(name))->SetTitle("");
  ((TH1D*)gDirectory->Get(name))->SetMaximum(max*1.2);
  ((TH1D*)gDirectory->Get(name))->SetMarkerStyle(8);
  ((TH1D*)gDirectory->Get(name))->SetMarkerSize(1.);
  ((TH1D*)gDirectory->Get(name))->SetMarkerColor(kBlack);
  ((TH1D*)gDirectory->Get(name))->SetLabelSize(.035,"X");
  ((TH1D*)gDirectory->Get(name))->SetLabelSize(.035,"Y");
  ((TH1D*)gDirectory->Get(name))->SetXTitle(nameaxis);
  ((TH1D*)gDirectory->Get(name))->SetMinimum(0.);
  ((TH1D*)gDirectory->Get(name))->SetStats(0);
  ((TH1D*)gDirectory->Get(name))->Draw();
  sprintf(name, "%s%s",var,"allmcexcl");
  ((TH1D*)gDirectory->Get(name))->SetLineWidth(3.);
  ((TH1D*)gDirectory->Get(name))->SetLineColor(kBlack);
  ((TH1D*)gDirectory->Get(name))->SetFillColor(kWhite);
  ((TH1D*)gDirectory->Get(name))->SetStats(0);
  ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
  sprintf(name, "%s%s",var,"allbkgexcl");
  ((TH1D*)gDirectory->Get(name))->SetFillColor(18);
  ((TH1D*)gDirectory->Get(name))->SetStats(0);
  ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
  sprintf(name, "%s%s",var,"nonvubexcl");
  ((TH1D*)gDirectory->Get(name))->SetFillColor(16);
  ((TH1D*)gDirectory->Get(name))->SetStats(0);
  ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
  sprintf(name, "%s%s",var,"scaleothexcl");
  ((TH1D*)gDirectory->Get(name))->SetFillColor(13);
  ((TH1D*)gDirectory->Get(name))->SetStats(0);
  ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
  sprintf(name, "%s%s",var,"dataexcl");
  ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
  TLegendEntry *legge;
  TLegend *leg;
  leg = new TLegend(0.13,0.45,0.35,0.85);
  leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  sprintf(name, "%s%s",var,"allmcexcl");
  legge = leg->AddEntry(((TH1D*)gDirectory->Get(name)), "signal", "f");
  sprintf(name, "%s%s",var,"allbkgexcl");
  legge = leg->AddEntry(((TH1D*)gDirectory->Get(name)), "b#rightarrowul#nu", "f");
  sprintf(name, "%s%s",var,"nonvubexcl");
  legge = leg->AddEntry(((TH1D*)gDirectory->Get(name)), "b#rightarrowcl#nu", "f");
  sprintf(name, "%s%s",var,"scaleothexcl");
  legge = leg->AddEntry(((TH1D*)gDirectory->Get(name)), "other", "f");
  sprintf(name, "%s%s",var,"dataexcl");
  legge = leg->AddEntry(((TH1D*)gDirectory->Get(name)), "data", "p");
  leg->Draw();

  xpos = 0.85;

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

  
  c1->SaveAs("mass.eps");
}

