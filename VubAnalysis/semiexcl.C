#include "../RecoilAnalysis/mesData.hh"
#include <TLatex.h>

recoilAnalysis a;

TLegend *leg;
TLegendEntry *entry;



// ----------------------------------------------------------------------
void loadMergedFiles(TString inputfile) {
  double dl    = 20.7; // data lumi
  double b0cl  = (3551000./1000000.)*198.; // b0cocktail lumi
  double bpcl  = (3551000./1000000.)*198.; // bpcocktail lumi
  double vubcl = (3551000./1000000.)*198.; // brecovubmix lumi

//   a.loadMc("b0-cocktail.root", b0cl);    //0
//   a.loadMc("bp-cocktail.root", bpcl);    //1
//   a.loadMc("b0-brecovubmix.root", vubcl); //2
//   a.loadMc("bp-brecovubmix.root", vubcl); //3
//   a.loadMc("b0-genbnu.root", 15.4);      //4
//   a.loadMc("bp-genbch.root", 20.9);      //5

//   a.loadDa("b0-run1.root", dl);   //0
//   a.loadDa("bp-run1.root", dl);   //1
//   a.loadDa("dstar-run1.root", dl); //2
//   a.loadDa("dstar0-run1.root", dl); //3
//   a.loadDa("dc-run1.root", dl);    //4 
//   a.loadDa("d0-run1.root", dl);    //5
  a.loadDa(inputfile.Data(), dl);   //0
}

// ----------------------------------------------------------------------
void loadFiles() {
  a.loadMc("dstar-b0cocktail.root");   //0
  a.loadMc("dc-b0cocktail.root");      //1
  a.loadMc("dstar0-bpcocktail.root");  //2
  a.loadMc("d0-bpcocktail.root");      //3
  a.loadMc("dstar-brecovubmix.root");  //4
  a.loadMc("dc-brecovubmix.root");     //5
  a.loadMc("dstar0-brecovubmix.root"); //6
  a.loadMc("d0-brecovubmix.root");     //7

}

// ----------------------------------------------------------------------
void plotMes(TString inputfile, TString outputdir, int file = 2) {

  char name[1000], command[1000], bmode[100], dmode3[100], dmode2[100],  dmode[100], thename[100], thename2[100], thename3[100];
  int min,  max;
  
  // dstar:  file 2, modes 13000 ... 14000
  if (file == 2) {
    min = 13000.; 
    max = 13400.;
    sprintf(dmode2, "DstarX"); 
    sprintf(dmode3, "%s%s","B0->",dmode2); 
  }

  // dc:  file 4, modes 12000 ... 13000
  if (file == 4) {
    min = 12000.; 
    max = 12500.;
    sprintf(dmode2, "DcX"); 
    sprintf(dmode3, "%s%s","B0->",dmode2); 
  }

  // dstar0:  file 3, modes 14000 ... 15000
  if (file == 3) {
    min = 14000.; 
    max = 15400.;
    sprintf(dmode2, "Dstar0X"); 
    sprintf(dmode3, "%s%s","B+->",dmode2); 
  }

  // d0:  file 5, modes 11000 ... 12000
  if (file == 5) {
    min = 11000.; 
    max = 11400.;
    sprintf(dmode2, "D0X"); 
    sprintf(dmode3, "%s%s","B+->",dmode2); 
  }
  
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  if (gFile == 0) {
    loadMergedFiles(inputfile);
  }
  a.fDa[0]->cd();

  zone(1);
  TH1D *h;
  mesData *mes;  

  TH1D messum("messum","sum of all modes",40,5.2,5.3);
  TH1D messum10("messum10","sum of all modes",40,5.2,5.3);
  TH1D messum20("messum20","sum of all modes",40,5.2,5.3);

  TH1D messumlep("messumlep","sum of all modes lepton cuts",40,5.2,5.3);
  TH1D messumall("messumall","sum of all modes all cuts",40,5.2,5.3);
  TH1D messumallmx("messumallmx","sum of all modes all cuts + mx cut",40,5.2,5.3);

  sprintf(name, "%s%s%d%s",outputdir.Data(),"/sxModes-",file,".html");
  ofstream OUT(name);
  OUT << "<!DOCTYPE HTML PUBLIC \"-//IETF//DTD HTML//EN\">" << endl;
  OUT << "<html>" << endl;
  OUT << "<head>" << endl;
  OUT << "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">" << endl;
  OUT << "<title>All B modes: </title>" << endl;
  OUT << "</head>" << endl;
  OUT << "<body bgcolor=white>" << endl;
  OUT << "<br>" << endl;
  OUT << "<br>" << endl;
  OUT << "<center><b><font size=+4> SemiExclusive Reco: " << dmode3 << " </font></b></center>" << endl;
  OUT << "<br>" << endl;
  OUT << "<br>" << endl;
  OUT << "<br>" << endl;
  OUT << "<br>" << endl;
  OUT << "<caption></caption>" << endl;
  OUT << "<table border=1  align=\"center\">" << endl;
  OUT << "<tr>" << endl;
  OUT << "<td><b><font color=red> <font size=+2>B mode</font></font></b></td>" 
      << "<td><b><font color=red> <font size=+2>D mode </font></font></b></td>"
      << "<td><b><font color=blue> <font size=+2>Yield </font></font></b></td>"
      << "<td><b><font color=blue> <font size=+2>Purity %</font></font></b></td>"
      << "<td><b><font color=blue> <font size=+2>S/sqrt(S+B)%</font></font></b></td>" 
      << endl;
  OUT << "</tr>" << endl;

  char signal[50], purity[50], ssb[50];
  double intpurity;
  int themode[600],thebmode[600], thedmode[600],  position[600];
  double thesig[600], theback[600], thepur[600], intpur[600];
  double thesigma = -1111111;
  double themean = -1111111;
  double theargus = -1111111;
  double thealpha = -1111111;
  double then = -1111111;
  double resmean, ressigma, resalpha, resn;
  int indexmode = 0;

  for (int i = 0 ; i<600; i++){
 
    thesig[i] = 0.;
    theback[i] = 100000.;
    themode[i] = intpur[i] = -1111;    
    thepur[i] = 0.;
  }

  for (int mode = min; mode < max; ++mode) {
    sprintf(name, "breco/h%d", mode); 
    h = (TH1D*)gFile->Get(name);
    if (!h) continue;
    //if (h->GetSumOfWeights() < 1.) continue;
    if(file == 4){
      if(mode%100>53) continue;
    }else{
      if(mode%100>52) continue;
    }
    if(mode%100==0) continue;
    if(mode>14360&& mode<15000) continue;

    mes = a.vubMes(h, resmean, ressigma, resalpha, resn, 1, 0, 5.28, 0.003, thealpha, 5., theargus);
    messum.Add(h,1);
    if (mes.thePu()>.1) messum10.Add(h,1);
    if (mes.thePu()>.2) messum20.Add(h,1);
 
    sprintf(name, "breco/ae%d", mode); 
    h = (TH1D*)gFile->Get(name);
    messumlep.Add(h,1);
    
    sprintf(name, "breco/ac%d", mode); 
    h = (TH1D*)gFile->Get(name);
    messumall.Add(h,1);
    
    sprintf(name, "breco/acm%d", mode); 
    h = (TH1D*)gFile->Get(name);
    messumallmx.Add(h,1);
    

    sprintf(name, "%s%s%d%s",outputdir.Data(),"/mes-",mode,".eps"); 
    c0->SaveAs(name);
    
    sprintf(signal, "%6.1f +/- %6.1f", mes.theSig(), mes.theErrSig());
    sprintf(purity, "%5.1f +/- %5.1f", 100.*mes.thePu(), 100.*mes.theErrPu());
    sprintf(ssb, "0.");

    Dmode(mode,dmode);
    Bmode(mode,bmode);

    sprintf(command, "gzip -f %s", name);
    cout << command << endl;
    gSystem->Exec(command);
    sprintf(name, "mes-%d.eps.gz", mode); 
    OUT << "<tr><td><a href=\"" << name << "\">" << bmode << "</td>" 
	<< "<td>" << dmode << "</td>"
	<< "<td>" << signal << "</td>"
	<< "<td>" << purity << "</td>"
	<< "<td>" << ssb << "</td>"
	<< "</tr>" << endl;
    thebmode[indexmode] =  mode%100;
    thedmode[indexmode] = mode/100 - 100;
    thesig[indexmode] = mes.theSig();
    theback[indexmode] = mes.theBg();  
    if(mes.theBg()==0) theback[indexmode] = 100000.;
    themode[indexmode] = mode;
    thepur[indexmode] = mes.thePu();
    indexmode ++;
  }

  for (int i=0; i<indexmode;i++){
     sprintf(name, "%s%d","mess", i); 
     h = new TH1D(name, name,40,5.2,5.3);  h->Sumw2();  
     sprintf(name, "%s%d","messlep", i); 
     h = new TH1D(name, name,40,5.2,5.3);  h->Sumw2();  
     sprintf(name, "%s%d","messall", i); 
     h = new TH1D(name, name,40,5.2,5.3);  h->Sumw2();  
     sprintf(name, "%s%d","messallmx", i); 
     h = new TH1D(name, name,40,5.2,5.3);  h->Sumw2();  
  }

  sprintf(name, "%s","mesall"); 
  h = new TH1D(name, name,40,5.2,5.3);  h->Sumw2();  

  sprintf(name, "%s","mesalllep"); 
  h = new TH1D(name, name,40,5.2,5.3);  h->Sumw2();  

  sprintf(name, "%s","mesallall"); 
  h = new TH1D(name, name,40,5.2,5.3);  h->Sumw2();  

  sprintf(name, "%s","mesallallmx"); 
  h = new TH1D(name, name,40,5.2,5.3);  h->Sumw2();  

  double maxpurity(100.);
  int counter(0);
  int preindex=-1;
 
  for (int i=0; i<indexmode; i++){
    double maxtemppurity(-1.);

      cout << preindex <<  "   " << thepur[i] << "   " << maxpurity<<endl; 
    for (int j=0; j<indexmode; j++){

      if(thepur[j]>maxtemppurity && thepur[j]<maxpurity) {
	preindex = j;
	maxtemppurity = thepur[j];
	position[i] = j;
      }
      if(thepur[j]==0&&maxpurity==0&&j>preindex){
	maxtemppurity = thepur[j];
	position[i] = j;
	preindex = j;
	break;
      }
    }

    maxpurity = maxtemppurity;
   }  
  

  for(int y=0;y<3;y++){
    if(y==0) mes = a.vubMes(&messum, resmean, ressigma, resalpha, resn, 1, 1, 5.28, 0.003, thealpha, 5., theargus);
    if(y==1) mes = a.vubMes(&messum10, resmean, ressigma, resalpha, resn, 1, 1, 5.28, 0.003, thealpha, 5., theargus);
    if(y==2) mes = a.vubMes(&messum20, resmean, ressigma, resalpha, resn, 1, 1, 5.28, 0.003, thealpha, 5., theargus);
    sprintf(name,"%s%s%d%s%s",outputdir.Data(),"/messum",y*10,dmode2,".eps", mode); 
    c0->SaveAs(name);
    
    sprintf(signal, "%6.1f +/- %6.1f", mes.theSig(), mes.theErrSig());
    sprintf(purity, "%5.1f +/- %5.1f", 100.*mes.thePu(), 100.*mes.theErrPu());
    sprintf(ssb, "0.");
    
    sprintf(ssb, "%5.1f", 100.*mes.theSig()/TMath::Abs(mes.theBg()*mes.theSig()));
    sprintf(command, "gzip -f %s", name);
    cout << command << endl;
    gSystem->Exec(command);
    sprintf(name, "%s%d%s%s","messum",y*10,dmode2,".eps.gz"); 
    OUT << "<tr><td><a href=\"" << name << "\">" << "Sum of all modes (ipur > " << y*10 << "%)" << "</td>" 
	<< "<td>" << " " << "</td>"
	<< "<td>" << signal << "</td>"
	<< "<td>" << purity << "</td>"
	<< "<td>" << ssb << "</td>"
	<< "</tr>" << endl;
  }
  
  OUT << "</table>" << endl;
  OUT << "</body>" << endl;
  OUT << "</html>" << endl;


  mes = a.vubMes(&messum, resmean, ressigma, resalpha, resn, 1, 1, 5.28, 0.003, thealpha, 5., theargus);    
    
  TH2D ssbvsyield("ssbvsyield","stat sign vs yield",200,0.,mes.theSig()*1.2,200,0, mes.theSig()/sqrt(mes.theSig()+mes.theBg())*1.4);

  TH2D purvsyield("purvsyield","stat sign vs yield",200,0.,mes.theSig()*1.2,200,0.,1.);

  double maxy = mes.theSig()*1.2;


  mes = a.vubMes(&messumlep, resmean, ressigma, resalpha, resn, 1, 1, 5.28, 0.003, thealpha, 5., theargus);    

  double maxylep = mes.theSig()*1.2;

  TH2D ssbvsyieldlep("ssbvsyieldlep","stat sign vs yield lep cuts",200,0.,mes.theSig()*1.2,200,0, mes.theSig()/sqrt(mes.theSig()+mes.theBg())*1.4);

  TH2D purvsyieldlep("purvsyieldlep","stat sign vs yield all cuts",200,0.,mes.theSig()*1.2,200,0.,1.);


  mes = a.vubMes(&messumall, resmean, ressigma, resalpha, resn, 1, 1, 5.28, 0.003, thealpha, 5., theargus);    

  TH2D ssbvsyieldall("ssbvsyieldall","stat sign vs yield all cuts",200,0.,mes.theSig()*1.2,200,0, mes.theSig()/sqrt(mes.theSig()+mes.theBg())*1.4);

  TH2D purvsyieldall("purvsyieldall","stat sign vs yield all cuts",200,0.,mes.theSig()*1.2,200,0.,1.);

  TH2D yvsy("yvsy","yield no cut vs yield all cuts",200,0.,maxy,200,0.,mes.theSig()*1.2);

  TH2D yvsylep("yvsylep","yield lep cut vs yield all cuts",200,0.,maxylep,200,0.,mes.theSig()*1.2);

  TH2D ratio("ratio","yield all cuts/yields no cut vs yield no cuts",200,0.,maxy,200,0.,.035);

  TH2D ratiolep("ratiolep","yield all cuts/yields lep cut vs yields lep",200,0.,maxylep,200,0.,0.35);


  mes = a.vubMes(&messumallmx, resmean, ressigma, resalpha, resn, 1, 1, 5.28, 0.003, thealpha, 5., theargus);    

  TH2D ssbvsyieldallmx("ssbvsyieldallmx","stat sign vs yield all cuts+mx",200,0.,mes.theSig()*1.2,200,0, mes.theSig()/sqrt(mes.theSig()+mes.theBg())*1.4);

  TH2D purvsyieldallmx("purvsyieldallmx","stat sign vs yield all cuts+mx",200,0.,mes.theSig()*1.2,200,0.,1.);


  sprintf(name, "%s%s%d%s",outputdir.Data(),"/table-",file,".dat");
  ofstream table(name);  
  
  for (int i=0; i<indexmode; i++){      

    sprintf(name, "breco/h%d", themode[position[i]]); 
    h = (TH1D*)gFile->Get(name);
    ((TH1D*)gDirectory->Get("mesall"))->Add(((TH1D*)gDirectory->Get(name)),1);
    
    sprintf(name, "%s%d","mess", i);  
    ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get("mesall")),1);
    
    mes = a.vubMes(((TH1D*)gDirectory->Get(name)), resmean, ressigma, resalpha, resn, 1, 1, 5.28, 0.003, thealpha, 5., theargus);  
    
    sprintf(name, "%s%s%s%s%d%s",outputdir.Data(),"/",name,"-",file,".eps"); 
    c0->SaveAs(name);
    sprintf(command, "gzip -f %s", name);
    cout << command << endl;
    gSystem->Exec(command);

    intpur[position[i]] = mes.thePu();
    
    double tempy = mes.theSig();
    
    ssbvsyield.Fill(mes.theSig(),mes.theSig()/sqrt(mes.theSig()+mes.theBg()));
    purvsyield.Fill(mes.theSig(),mes.thePu());
     
    
    sprintf(name, "breco/ae%d", themode[position[i]]); 
    h = (TH1D*)gFile->Get(name);
    ((TH1D*)gDirectory->Get("mesalllep"))->Add(((TH1D*)gDirectory->Get(name)),1);

    sprintf(name, "%s%d","messlep", i);  
    ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get("mesalllep")),1);
    
    mes = a.vubMes(((TH1D*)gDirectory->Get(name)), resmean, ressigma, resalpha, resn, 1, 1, 5.28, 0.003, thealpha, 5., theargus);  
    
    sprintf(name, "%s%s%s%s%d%s",outputdir.Data(),"/",name,"-",file,".eps"); 
    c0->SaveAs(name);
    sprintf(command, "gzip -f %s", name);
    cout << command << endl;
    gSystem->Exec(command);
    
    ssbvsyieldlep.Fill(mes.theSig(),mes.theSig()/sqrt(mes.theSig()+mes.theBg()));
    purvsyieldlep.Fill(mes.theSig(),mes.thePu());
    
    double tempylep = mes.theSig();


    sprintf(name, "breco/ac%d", themode[position[i]]); 
    h = (TH1D*)gFile->Get(name);
    ((TH1D*)gDirectory->Get("mesallall"))->Add(((TH1D*)gDirectory->Get(name)),1);

    sprintf(name, "%s%d","messall", i);  
    ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get("mesallall")),1);
    
    mes = a.vubMes(((TH1D*)gDirectory->Get(name)), resmean, ressigma, resalpha, resn, 1, 1, 5.28, 0.003, thealpha, 5., theargus);  
    
    sprintf(name, "%s%s%s%s%d%s",outputdir.Data(),"/",name,"-",file,".eps"); 
    c0->SaveAs(name);
    sprintf(command, "gzip -f %s", name);
    cout << command << endl;
    gSystem->Exec(command);
         
    ssbvsyieldall.Fill(mes.theSig(),mes.theSig()/sqrt(mes.theSig()+mes.theBg()));
    purvsyieldall.Fill(mes.theSig(),mes.thePu());

    yvsy.Fill(tempy, mes.theSig());

    yvsylep.Fill(tempylep, mes.theSig());
    
    ratio.Fill(tempy,mes.theSig()/tempy);

    ratiolep.Fill(tempylep,mes.theSig()/tempylep);
    
 
    sprintf(name, "breco/acm%d", themode[position[i]]); 
    h = (TH1D*)gFile->Get(name);
    ((TH1D*)gDirectory->Get("mesallallmx"))->Add(((TH1D*)gDirectory->Get(name)),1);

    sprintf(name, "%s%d","messallmx", i);  
    ((TH1D*)gDirectory->Get(name))->Add(((TH1D*)gDirectory->Get("mesallallmx")),1);
    
    mes = a.vubMes(((TH1D*)gDirectory->Get(name)), resmean, ressigma, resalpha, resn, 1, 1, 5.28, 0.003, thealpha, 5., theargus);  
 
    ssbvsyieldallmx.Fill(mes.theSig(),mes.theSig()/sqrt(mes.theSig()+mes.theBg()));
    purvsyieldallmx.Fill(mes.theSig(),mes.thePu());

  }
    
  ssbvsyield.Draw();
  sprintf(name, "%s%s%s%d%s",outputdir.Data(),"/","ssbvsyield-",file,".eps");
  c0.SaveAs(name);

  purvsyield.Draw();
  sprintf(name, "%s%s%s%d%s",outputdir.Data(),"/","puvsyieldr-",file,".eps");
  c0.SaveAs(name);
  

  ssbvsyieldlep.Draw();
  sprintf(name, "%s%s%s%d%s",outputdir.Data(),"/","ssbvsyieldlep-",file,".eps");   
  c0.SaveAs(name);

  purvsyieldlep.Draw();
  sprintf(name, "%s%s%s%d%s",outputdir.Data(),"/","puvsyieldlep-",file,".eps");
  c0.SaveAs(name);

  
  ssbvsyieldall.Draw();
  sprintf(name, "%s%s%s%d%s",outputdir.Data(),"/","ssbvsyieldall-",file,".eps");   
  c0.SaveAs(name);

  purvsyieldall.Draw();
  sprintf(name, "%s%s%s%d%s",outputdir.Data(),"/","puvsyieldall-",file,".eps");
  c0.SaveAs(name);

  
  yvsy.Draw();
  sprintf(name, "%s%s%s%d%s",outputdir.Data(),"/","yvsy-",file,".eps");   
  c0.SaveAs(name);

  yvsylep.Draw();
  sprintf(name, "%s%s%s%d%s",outputdir.Data(),"/","yvsylep-",file,".eps");   
  c0.SaveAs(name);

  ratio.Draw();
  sprintf(name, "%s%s%s%d%s",outputdir.Data(),"/","ratio-",file,".eps");   
  c0.SaveAs(name);

  ratiolep.Draw();
  sprintf(name, "%s%s%s%d%s",outputdir.Data(),"/","ratiolep-",file,".eps");   
  c0.SaveAs(name);


  ssbvsyieldallmx.Draw();
  sprintf(name, "%s%s%s%d%s",outputdir.Data(),"/","ssbvsyieldallmx-",file,".eps");   
  c0.SaveAs(name);

  purvsyieldallmx.Draw();
  sprintf(name, "%s%s%s%d%s",outputdir.Data(),"/","puvsyieldallmx-",file,".eps");
  c0.SaveAs(name);

  
  for (int i=0; i<indexmode; i++){      
    cout << " i " << i << "position " << position[i] << "   " <<themode[position[i]]<<endl;
  }
  
  for (int i=0; i<indexmode; i++) {
    
    double sblock;

    if(intpur[i]>.8) {sblock = 1;}
    else if(intpur[i]>.5) {sblock = 2;}
    else if(thepur[i]>.1) {sblock = 3;}
    else {sblock = 4;}
    
    table << thebmode[i] << "   " <<  thedmode[i] << "   " << thesig[i] << "   " <<  theback[i] << "   " <<  intpur[i] << "   " <<  sblock <<endl;
   
  }

  cout << indexmode << endl;

  sprintf(name, "%s%s%d%s",outputdir.Data(),"/sx-oredered-Modes-",file,".html");
  ofstream OUT(name);

  OUT << "<!DOCTYPE HTML PUBLIC \"-//IETF//DTD HTML//EN\">" << endl;
  OUT << "<html>" << endl;
  OUT << "<head>" << endl;
  OUT << "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">" << endl;
  OUT << "<title>All B modes ordered by purity: </title>" << endl;
  OUT << "</head>" << endl;
  OUT << "<body bgcolor=white>" << endl;
  OUT << "<br>" << endl;
  OUT << "<br>" << endl;
  OUT << "<center><b><font size=+4> SemiExclusive Reco: " << dmode3 << " </font></b></center>" << endl;
  OUT << "<br>" << endl;
  OUT << "<br>" << endl;
  OUT << "<br>" << endl;
  OUT << "<br>" << endl;
  OUT << "<caption></caption>" << endl;
  OUT << "<table border=1  align=\"center\">" << endl;
  OUT << "<tr>" << endl;
  OUT << "<td><b><font color=red> <font size=+2>B mode</font></font></b></td>" 
      << "<td><b><font color=red> <font size=+2>D mode </font></font></b></td>"
      << "<td><b><font color=blue> <font size=+2>Yield </font></font></b></td>"
      << "<td><b><font color=blue> <font size=+2>Purity %</font></font></b></td>"
      << "<td><b><font color=blue> <font size=+2>Int.purity%</font></font></b></td>" 
       << "<td><b><font color=blue> <font size=+2>Int.plot</font></font></b></td>" 
       << "<td><b><font color=blue> <font size=+2>Int.plot leptcuts</font></font></b></td>" 
       << "<td><b><font color=blue> <font size=+2>Int.plot allcuts</font></font></b></td>" 
      << endl;
  OUT << "</tr>" << endl;
   
  for (int i=0; i<indexmode; i++) {
  
    
    Dmode(themode[position[i]],dmode);
    Bmode(themode[position[i]],bmode);

    sprintf(name, "mes-%d.eps.gz",themode[position[i]]); 
    sprintf(thename, "%s%d%s%d%s","mess",i,"-",file,".eps.gz"); 
    sprintf(thename2, "%s%d%s%d%s","messlep",i,"-",file,".eps.gz"); 
    sprintf(thename3, "%s%d%s%d%s","messall",i,"-",file,".eps.gz"); 
    OUT << "<tr><td><a href=\"" << name << "\">" << bmode << "</td>" 
	<< "<td>" << dmode << "</td>"
	<< "<td>" << thesig[position[i]] << "</td>"
	<< "<td>" << thepur[position[i]] << "</td>"
	<< "<td>" << intpur[position[i]] << "</td>" 
	<< "<td>" << "<a href=\"" << thename << "\">plot</td>"
	<< "<td>" << "<a href=\"" << thename2 << "\">plot</td>"
	<< "<td>" << "<a href=\"" << thename3 << "\">plot</td>"
	<< "</tr>" << endl;
    
  }

  OUT << "</table>" << endl;
  OUT << "</body>" << endl;
  OUT << "</html>" << endl;


}



void  Bmode(int theid, char filename[100]){
  
  char filename[100];

  int bmode = theid%100;  
  
  if(bmode == 1)  sprintf(filename,"%s%s","B->D","pi");
  if(bmode == 2)  sprintf(filename,"%s%s","B->D","k");
  if(bmode == 3)  sprintf(filename,"%s%s","B->D","pipi0_<1.5GeV");
  if(bmode == 4)  sprintf(filename,"%s%s","B->D","kpi0_<1.5GeV");
  if(bmode == 5)  sprintf(filename,"%s%s","B->D","piks");
  if(bmode == 6)  sprintf(filename,"%s%s","B->D","kks");
  if(bmode == 7)  sprintf(filename,"%s%s","B->D","pi2pi0_<1.5GeV");
  if(bmode == 8)  sprintf(filename,"%s%s","B->D","k2pi0_<1.5GeV");
  if(bmode == 9)  sprintf(filename,"%s%s","B->D","3pi_<1.5GeV");
  if(bmode == 10)  sprintf(filename,"%s%s","B->D","k2pi_<1.5GeV");
  if(bmode == 11)  sprintf(filename,"%s%s","B->D","2kpi_Ds");
  if(bmode == 12)  sprintf(filename,"%s%s","B->D","omegah");
  if(bmode == 13)  sprintf(filename,"%s%s","B->D","k2pipi0_<2.2GeV");
  if(bmode == 14)  sprintf(filename,"%s%s","B->D","2kpipi0_Ds*");
  if(bmode == 15)  sprintf(filename,"%s%s","B->D","pipi0ks");
  if(bmode == 16)  sprintf(filename,"%s%s","B->D","kpi0ks_<1.8GeV");
  if(bmode == 17)  sprintf(filename,"%s%s","B->D","k2pi0ks_1.8-2.2GeV");
  if(bmode == 18)  sprintf(filename,"%s%s","B->D","2ksX");
  if(bmode == 19)  sprintf(filename,"%s%s","B->D","3pi2pi0_<2.2GeV");
  if(bmode == 20)  sprintf(filename,"%s%s","B->D","k2pi2pi0_<2.2GeV");
  if(bmode == 21)  sprintf(filename,"%s%s","B->D","2kpi2pi0_Ds*");
  if(bmode == 22)  sprintf(filename,"%s%s","B->D","5pi_<2.3GeV");
  if(bmode == 23)  sprintf(filename,"%s%s","B->D","k4p_<2.7GeV");
  if(bmode == 24)  sprintf(filename,"%s%s","B->D","2K3pi_<2.7GeV");
  if(bmode == 25)  sprintf(filename,"%s%s","B->D","5pipi0_<2.2GeV");
  if(bmode == 26)  sprintf(filename,"%s%s","B->D","k4pipi0_<2.2GeV");
  if(bmode == 27)  sprintf(filename,"%s%s","B->D","2k3pipi0_<2.5GeV");
  if(bmode == 28)  sprintf(filename,"%s%s","B->D","3piks_D*");
  if(bmode == 29)  sprintf(filename,"%s%s","B->D","3pikspi0_D*");
  if(bmode == 30)  sprintf(filename,"%s%s","B->D","k2piks_D*");
  if(bmode == 31)  sprintf(filename,"%s%s","B->D","D*_Dpi0");
  if(bmode == 32)  sprintf(filename,"%s%s","B->D","pipi0_>1.5GeV");
  if(bmode == 33)  sprintf(filename,"%s%s","B->D","kpi0_>1.5GeV");
  if(bmode == 34)  sprintf(filename,"%s%s","B->D","pi2pi0_1.5-2GeV");
  if(bmode == 35)  sprintf(filename,"%s%s","B->D","k2pi0_>1.5GeV");
  if(bmode == 36)  sprintf(filename,"%s%s","B->D","3pi_1.5-2GeV");
  if(bmode == 37)  sprintf(filename,"%s%s","B->D","k2pi_>1.5GeV");
  if(bmode == 38)  sprintf(filename,"%s%s","B->D","2kpi_K*");
  if(bmode == 39)  sprintf(filename,"%s%s","B->D","2kpi_other");
  if(bmode == 40)  sprintf(filename,"%s%s","B->D","3pipi0_<1.6GeV");
  if(bmode == 41)  sprintf(filename,"%s%s","B->D","3pipi0_1.6-2.2GeV");
  if(bmode == 42)  sprintf(filename,"%s%s","B->D","k2pipi0_>2.2GeV");
  if(bmode == 43)  sprintf(filename,"%s%s","B->D","2kpipi0_other");
  if(bmode == 44)  sprintf(filename,"%s%s","B->D","kpi0ks_>1.8GeV");
  if(bmode == 45)  sprintf(filename,"%s%s","B->D","3pi2pi0_>2.2GeV");
  if(bmode == 46)  sprintf(filename,"%s%s","B->D","k2pi2pi0_>2.2GeV");
  if(bmode == 47)  sprintf(filename,"%s%s","B->D","2kpi2pi0_other");
  if(bmode == 48)  sprintf(filename,"%s%s","B->D","5pi_>2.3GeV");
  if(bmode == 49)  sprintf(filename,"%s%s","B->D","k4p_>2.7GeV");
  if(bmode == 50)  sprintf(filename,"%s%s","B->D","2K3pi_>2.7GeV");
  if(bmode == 51)  sprintf(filename,"%s%s","B->D","5pipi0_>2.2GeV");
  if(bmode == 52)  sprintf(filename,"%s%s","B->D","3piks_noD*");
  if(bmode == 53)  sprintf(filename,"%s%s","B->D","3pikspi0_noD*");

  return;

}

void Dmode(int theid, char filename[100]){
  
  char filename[100];
  char seed[100];

  int dmode = theid/100;
  
  if(dmode == 110)   sprintf(filename,"D0->kpi");
  if(dmode == 111)   sprintf(filename,"D0->kpipi0");
  if(dmode == 112)   sprintf(filename,"D0->k3pi");
  if(dmode == 113)   sprintf(filename,"D0->kspipi");
  if(dmode == 130)   sprintf(filename,"D*,D0->kpi");
  if(dmode == 131)   sprintf(filename,"D*,D0->kpipi0");
  if(dmode == 132)   sprintf(filename,"D*,D0->k3pi");
  if(dmode == 133)   sprintf(filename,"D*,D0->kspipi");
  if(dmode == 120)   sprintf(filename,"Dc->kspi");
  if(dmode == 121)   sprintf(filename,"Dc->kpipi");
  if(dmode == 122)   sprintf(filename,"Dc->kspipi0");
  if(dmode == 123)   sprintf(filename,"Dc->kpipipi0");
  if(dmode == 124)   sprintf(filename,"Dc->kspipipi");
  if(dmode == 140)   sprintf(filename,"D*0->D0pi0,D0->kpi");
  if(dmode == 141)   sprintf(filename,"D*0->D0pi0,D0->kpipi0");
  if(dmode == 142)   sprintf(filename,"D*0->D0pi0,D0->k3pi");
  if(dmode == 143)   sprintf(filename,"D*0->D0pi0,D0->kspipi");
  if(dmode == 150)   sprintf(filename,"D*0->D0gamma,D0->kpi");
  if(dmode == 151)   sprintf(filename,"D*0->D0gamma,D0->kpipi0");
  if(dmode == 152)   sprintf(filename,"D*0->D0gamma,D0->k3pi");
  if(dmode == 153)   sprintf(filename,"D*0->D0gamma,D0->kspipi");

  return;

}

