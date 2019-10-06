void domesplots(char* chaininput, char* modename){

  
  TCanvas c1("c1","--c1--",472,0,700,900);
  char name[100],title[100],outputname[100];
  TH1D mesbch("mesbch", "", 40, 5.2, 5.3);  mesbch.Sumw2();
  sprintf(title,"%s%s","    B^{+} mes for ", modename);
  mesbch.SetTitle(title);
  TH1D mesb0("mesb0", "", 40, 5.2, 5.3);  mesb0.Sumw2();
  sprintf(title,"%s%s","    B^{0} mes for ", modename);
  mesb0.SetTitle(title);
//   sprintf(name,"mesb0os");  sprintf(title, "mes b0 os");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
//   sprintf(name,"mesb0ss");  sprintf(title, "mes b0 ss");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  
  TChain *pippo = getchain(chaininput);
  
  pippo.Project("mesbch","mes","pcms>1.&&brecocharge!=0&&pur>.06");
  pippo.Project("mesb0","mes","pcms>1.&&brecocharge==0&&pur>0.06");
  
  recoilAnalysis b;
  double resmean, ressigma, resalpha, resn;
  double thesigma = -1111111;
  double themean = 5.2795;
  double theargus = -1111111;
  double thealpha = -1111111;
  double then = 5;
  b.vubMes((TH1D*)gDirectory->Get("mesbch"), resmean, ressigma, resalpha, resn, 1, 1, themean, thesigma, thealpha, then, theargus);
  sprintf(outputname,"%s%s",chaininput,"_bch.eps");
  c1.SaveAs(outputname);

  b.vubMes((TH1D*)gDirectory->Get("mesb0"), resmean, ressigma, resalpha, resn, 1, 1, themean, thesigma, thealpha, then, theargus);
  sprintf(outputname,"%s%s",chaininput,"_b0.eps");
  c1.SaveAs(outputname);
  
}

// ----------------------------------------------------------------------
TChain * getchain(char *thechain) {

  TChain *chain = new TChain("events");
  cout << "Chaining ... " << thechain << endl;
  char pName[2000]; 
  char buffer[200];
  sprintf(buffer, "%s", thechain);
  ifstream is(buffer);  
  cout << "files " << buffer <<  endl;
  while(is.getline(buffer, 200, '\n')){
    //    if (buffer[0] == '#') continue; 
    sscanf(buffer, "%s", pName); 
    cout << "   Add: " << buffer << endl; 
    chain->Add(pName); 
  }
  is.close();
  return chain;

}
