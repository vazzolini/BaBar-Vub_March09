#include <fstream.h>

#include "theanal.hh"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "../RecoilAnalysis/mesData.hh" 
#include <TVector2.h>
theanal::theanal(char *var, char *dir, double mi, double ma, int b)
{
  char name[100];
  sprintf(name, "%s%s%s%s",dir,"/",var,".root");
  fHistFile = new TFile(name, "RECREATE");
   themin = mi;
   themax = ma;
   thebins = b;
   thevar = var ;
   TFile dal("dalitz.root");
   mydalitz = (*(TH2D*)gDirectory->Get("truedalitzsigforcut"));
   //    Dvar = new recoilDSys("ddecay.table",0,0);
   //    Bsem = new recoilDSys(0);
   
}
 
void theanal::Loop(int nevents)
{

  fHistFile->cd();
  
  char name[100];
  char le[100];
  
  //id = Mx category == Mx bin
  int id(0);
  int idflav;
  int group;        
  
  if (fChain == 0) return;
  
  Int_t nentries = Int_t(fChain->GetEntries());
  
  if( nentries > nevents) nentries = nevents;
  
  cout <<  nentries << endl;

  Int_t nbytes = 0, nb = 0;
  for (Int_t jentry=0; jentry<nentries;jentry++) {
    LoadTree(jentry); //in case of a TChain, ientry is the entry number in the current file
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%10000 == 0) cout << "Event " << jentry << endl;

    int flav =  lcharge + brecoflav; // charge correlation
    bool ksele = nkp + nks;          // fit on the depleted sample?
    int ch = xcharge + brecocharge;  // total charge
    int cutflav = 0;              // flavor category (3 = charged B, 4 = neutral B OS, 5 = neutral B SS)  

    double themm2=mm2;
    if(FITCATEGORY == -11 && nneu !=1) themm2 = mm2bestPi0;
    if(FITCATEGORY == -11 && nneu ==1) themm2 = mm2gamma;
    if(FITCATEGORY == 11) themm2 = mm2bestPi;
    if(FITCATEGORY == -12) themm2 = mm2bestEta;
    if(FITCATEGORY == -13) themm2 = mm2rho0;
    if(FITCATEGORY == -14) themm2 = mm2omega;
    if(FITCATEGORY == 13) themm2 = mm2rho;
    if(FITCATEGORY == -15) themm2 = mm2bestEtap;
    if(FITCATEGORY == -19) themm2 = mm2a0;
    if(FITCATEGORY == 19) themm2 = mm2a0p;
    
    bool isLp = ((nlept500> 0) && (tlab<2.37) && (tlab>0.36) && (plab>0.5));
    bool PCMS =  isLp && (pcms>LEPTONPCUT);
    if(isele == 1) PCMS = isLp && (pcms > ELECTRONPCUT);
    if(isele == 0) PCMS = isLp && (pcms > MUONPCUT);

    bool NLE = nle <= 1;
    bool MM2 = (themm2 < MNUSQHIGH && themm2 > MNUSQLOW);
    if(FITCATEGORY == -11 && nneu ==1)  MM2 = (themm2 < 2.);
    bool CH = (ch >CHLOW  && ch < CHHIGH);
    //bool KSELE = ksele >0;
    bool KSELE = 0;
    int FLAV = !(TMath::Abs(brecocharge)!=0 && (lcharge + brecoflav)!=0);
    bool BCH = (TMath::Abs(brecocharge) == BTYPE);     
    bool IPUR = intpur>MININTPUR && intpur<MAXINTPUR;
    bool NCHG = nchg-1 >= NCHGLOWEXCL && nchg-1 <= NCHGHIGHEXCL;
    bool MX = mxhadfit < MXCUTHIGHEXCL && mxhadfit > MXCUTLOWEXCL;
    bool OTHER(1);
    bool DECAYMODE(1);
    if(FITCATEGORY == -11) OTHER = nrecoPi0>0;
    if(FITCATEGORY == 11) {
      OTHER = nrecoPi>0;
      OTHER = OTHER && (Eneualt<MAXENEU ) && (TMath::Abs(baremassjpsiPi[indexbestPi]-3.1)>JPSIWIN);
      OTHER = OTHER && ((chbestPi*lcharge)<0);
    }
    if(FITCATEGORY == -12) OTHER = nrecoEta>0;
    if(FITCATEGORY == -15) OTHER = nrecoEtap>0;
    int CATEGORY=0;
    if (Gvxbtyp == FITCATEGORY) CATEGORY = 1;
    if (vub && Gvxbtyp != FITCATEGORY) CATEGORY = 2;
    if (vcb) CATEGORY = 3;
    if (vub + vcb == 0) CATEGORY = 4;
    if( CATEGORY==0 ) cout << "no category: some mismatch..." <<  endl;
    bool MOM1(1), MOM2(1); 
    bool COMB = 1;
    bool PIZ = 1;
    bool DALITZ = 1;
    bool WdeltaCut = (wdeltam>PRMM2CUT && brecocharge == 0); //last cut based on wdelta
    //PIZ =  nrecoPi0 >= NPI0LOWEXCL && nrecoPi0 <= NPI0HIGHEXCL;
    bool DAUETA = 1;
    bool DAURHO = 1;
    bool DAUGAMMOM = 1;
    if(FITCATEGORY == -11){
      MOM1 = Estar1dauPi0[indexbestPi0]>MOM1MIN;
      MOM2 = Estar2dauPi0[indexbestPi0]>MOM2MIN;        
      MX = (barembestPi0 < MXCUTHIGHEXCL && barembestPi0 > MXCUTLOWEXCL);
    }
    if(FITCATEGORY == -12){
      MX = (barembestEta < MXCUTHIGHEXCL && barembestEta > MXCUTLOWEXCL);
      if(modeEta[indexbestEta]==1) MX = (barembestEta < MXCUTHIGHEXCL1 && barembestEta > MXCUTLOWEXCL1);
      if(modeEta[indexbestEta]==2) MX = (barembestEta < MXCUTHIGHEXCL2 && barembestEta > MXCUTLOWEXCL2);
      if(modeEta[indexbestEta]==3) MX = (barembestEta < MXCUTHIGHEXCL3 && barembestEta > MXCUTLOWEXCL3);
      if(DORIGHTNCHG){
	NCHG = 0;
	if(modeEta[indexbestEta] == 2) { NCHG = nchg-1 == 2; }
	else { NCHG = nchg-1 == 0; }
      }
    }
    if(FITCATEGORY == -13){
      MOM1 = mom1pirho0>MOM1MIN;
      MOM2 = mom2pirho0>MOM2MIN;        
      KSELE = 0;
      MX = (mrho0 < MXCUTHIGHEXCL && mrho0 > MXCUTLOWEXCL);
      COMB =  nrecorho0 >= NCOMBLOWEXCL && nrecorho0 <= NCOMBHIGHEXCL;
    }      
    if(FITCATEGORY == 13){
      MOM1 = mompirho>MOM1MIN;
      MOM2 = mompi0rho>MOM2MIN;        
      MX = (mrho < MXCUTHIGHEXCL && mrho > MXCUTLOWEXCL);
      COMB =  nrecorho >= NCOMBLOWEXCL && nrecorho <= NCOMBHIGHEXCL;
    }      
    if(FITCATEGORY == -14){
      MOM1 = mom1piome>MOM1MIN;
      MOM2 = mom2piome>MOM2MIN;        
      MX = (momega < MXCUTHIGHEXCL && momega > MXCUTLOWEXCL);
      COMB =  nrecoomega >= NCOMBLOWEXCL && nrecoomega <= NCOMBHIGHEXCL;
      int  xbindal = int(dalitzpi1pi0ome/0.016);
      int  ybindal = int(dalitzpi1pi2ome/0.016);
      if(xbindal<0) xbindal = 1;
      if(ybindal<0) ybindal = 1;      
      if(mydalitz.GetBinContent(xbindal,ybindal)<DALITZCUT) DALITZ = 0;
    }      
    if(FITCATEGORY == -15){
      MX = (barembestEtap < MXCUTHIGHEXCL && barembestEtap > MXCUTLOWEXCL);      
      if(modeEtap[indexbestEtap]==1) MX = (barembestEtap < MXCUTHIGHEXCL1 && barembestEtap > MXCUTLOWEXCL1);
      if(modeEtap[indexbestEtap]==2) MX = (barembestEtap < MXCUTHIGHEXCL2 && barembestEtap > MXCUTLOWEXCL2);
      if(modeEtap[indexbestEtap]==3) MX = (barembestEtap < MXCUTHIGHEXCL3 && barembestEtap > MXCUTLOWEXCL3);
      if(modeEtap[indexbestEtap]==4) MX = (barembestEtap < MXCUTHIGHEXCL4 && barembestEtap > MXCUTLOWEXCL4);
      DAUGAMMOM = GammamomdauEtap[indexbestEtap]>DAUGAMMAMOM || GammamomdauEtap[indexbestEtap]<-10.;
      if(modeEtap[indexbestEtap] == 1) {
	DAURHO = TMath::Abs(Rho0massdauEtap[indexbestEtap]-0.775) < DAURHOMASS;
      }else{
	DAUETA = TMath::Abs(EtamassdauEtap[indexbestEtap]-0.54775)<DAUETAMASS;
	if(modeEtap[indexbestEtap]==2) DAUETA= TMath::Abs(EtamassdauEtap[indexbestEtap]-0.54775)<DAUETAMASS2;
        if(modeEtap[indexbestEtap]==3) DAUETA= TMath::Abs(EtamassdauEtap[indexbestEtap]-0.54775)<DAUETAMASS3;
        if(modeEtap[indexbestEtap]==4) DAUETA= TMath::Abs(EtamassdauEtap[indexbestEtap]-0.54775)<DAUETAMASS4;
      }
      if(DORIGHTNCHG){
	NCHG = 0;
	if(modeEtap[indexbestEtap] == 3) { NCHG = nchg-1 == 4; }	
	else { NCHG = nchg-1 == 2; }
      }
    }
    if(FITCATEGORY == -19){
      MX = (ma0 < MXCUTHIGHEXCL && ma0 > MXCUTLOWEXCL);
      DAUETA = TMath::Abs(a0massetadau-0.54775)<DAUETAMASS;
      if(DORIGHTNCHG){
	NCHG = 0;
	if(modea0 == 2) { NCHG = nchg-1 == 2; }
	else { NCHG = nchg-1 == 0; }
      }
    }
    if(FITCATEGORY == 19){
      MX = (ma0p < MXCUTHIGHEXCL && ma0p > MXCUTLOWEXCL);
      DAUETA = TMath::Abs(a0pmassetadau-0.54775)<DAUETAMASS;
      if(DORIGHTNCHG){
	NCHG = 0;
	if(modea0p == 2) { NCHG = nchg-1 == 3; }
	else { NCHG = nchg-1 == 1; }
      }
    }

    // mm2 cut in order to remove Vub crossfeed
    bool MM2CROSS = 1;
    bool MM2PI0 = (mm2bestPi0 >= MNUSQPI0HIGH || mm2bestPi0 <= MNUSQPI0LOW);
    bool MM2ETA = (mm2bestEta >= MNUSQETAHIGH || mm2bestEta <= MNUSQETALOW);
    bool MM2RHO = (mm2rho >= MNUSQRHOHIGH || mm2rho <= MNUSQRHOLOW);
    bool MM2RHO0 = (mm2rho0 >= MNUSQRHO0HIGH || mm2rho0 <= MNUSQRHO0LOW);
    bool MM2OMEGA = (mm2omega >= MNUSQOMEGAHIGH || mm2omega <= MNUSQOMEGALOW);
    MM2CROSS = MM2PI0*MM2ETA;
    

    if (BTYPE == 2) BCH = 1;
    double myvar = -1;
    if (strcmp(thevar,"mes") == 0) {myvar = mes; sprintf(xname,"m_{es} [GeV]");}
    if (strcmp(thevar,"mm2") == 0) {myvar = mm2; MM2=1;  sprintf(xname,"m_{miss}^{2} [GeV^{2}]");}
    if (strcmp(thevar,"prmm2") == 0) {myvar = wdeltam; WdeltaCut=0;  sprintf(xname,"m_{missPR}^{2} [GeV^{2}]");}
    if (strcmp(thevar,"nchg") == 0) { myvar = nchg-1; NCHG = 1; CH = 1;  sprintf(xname,"N_{chg}");}
    if (strcmp(thevar,"nneu") == 0) { myvar = nneu; sprintf(xname,"N_{neu}");}
    if (strcmp(thevar,"eneu") == 0) { myvar = Eneualt; OTHER = 1; sprintf(xname,"E_{neutral} [GeV]");} 
    if (strcmp(thevar,"nrecopi0") == 0) { myvar = nrecoPi0; PIZ = 1; sprintf(xname,"N_{#pi0}");}
    if (strcmp(thevar,"pcms") == 0) {myvar = pcms; PCMS = 1; sprintf(xname,"P_{cms} [GeV]");}
    if (strcmp(thevar,"qtot") == 0) {myvar = ch; CH = 1;NCHG = 1; sprintf(xname,"#Sigma_{i} chg_{i}");}
    if (strcmp(thevar,"mm2pi") == 0) {myvar = mm2bestPi; MM2 =1; sprintf(xname,"m_{miss}^{2}(#pi^{+}) [GeV^{2}]");}
    if (strcmp(thevar,"mpi0") == 0) {myvar = barembestPi0; MX =1;sprintf(xname,"m(#pi0) [GeV]");}
    if (strcmp(thevar,"resompi0") == 0) {myvar = barembestPi0-mxhadgen; MX =1;sprintf(xname,"m_{#pi0}-m^{true}_{#pi0} [GeV]");}
    if (strcmp(thevar,"mm2pi0") == 0) {myvar = mm2bestPi0; MM2 =1; MM2PI0=1; sprintf(xname,"m_{miss}^{2}(#pi^{0}) [GeV^{2}]");}
    if (strcmp(thevar,"mm2gamma") == 0) {myvar = mm2gamma; MM2 =1;sprintf(xname,"mm2_{#gamma} [GeV^{2}]");}
    if (strcmp(thevar,"mom1phpi0") == 0) {myvar = Estar1dauPi0[indexbestPi0]; MOM1 =1;sprintf(xname,"p^{*}_{#gamma1}(#pi^{0}) [GeV]");}
    if (strcmp(thevar,"mom2phpi0") == 0) {myvar = Estar2dauPi0[indexbestPi0]; MOM2 =1;sprintf(xname,"mom2#gamma_{#pi0} [GeV]");}
    if (strcmp(thevar,"truemom1phpi0") == 0) {myvar = truemom1phpi0; MOM1 =1;sprintf(xname,"truemom1#gamma_{#pi0} [GeV]");}
    if (strcmp(thevar,"truemom2phpi0") == 0) {myvar = truemom2phpi0; MOM2 =1;sprintf(xname,"truemom2#gamma_{#pi0} [GeV]");}
    if (strcmp(thevar,"nrecorho0") == 0) { myvar = nrecorho0; COMB = 1; sprintf(xname,"N_{#rho0}");}
    if (strcmp(thevar,"mrho0") == 0) {myvar = mrho0; MX =1;sprintf(xname,"m_{#rho0} [GeV]");}
    if (strcmp(thevar,"resomrho0") == 0) {myvar = mrho0-mxhadgen; MX =1;sprintf(xname,"m_{#rho0}-m^{true}_{#rho0} [GeV]");}
    if (strcmp(thevar,"mm2rho0") == 0) {myvar = mm2rho0; MM2 =1; MM2RHO0=1; sprintf(xname,"mm2_{#rho0} [GeV^{2}]");}
    if (strcmp(thevar,"mom1pirho0") == 0) {myvar = mom1pirho0; MOM1 =1;sprintf(xname,"mom1#pi_{#rho0} [GeV]");}
    if (strcmp(thevar,"mom2pirho0") == 0) {myvar = mom2pirho0; MOM2 =1;sprintf(xname,"mom2#pi_{#rho0} [GeV]");}
    if (strcmp(thevar,"truemom1pirho0") == 0) {myvar = truemom1pirho0; MOM1 =1;sprintf(xname,"truemom1#pi_{#rho0} [GeV]");}
    if (strcmp(thevar,"truemom2pirho0") == 0) {myvar = truemom2pirho0; MOM2 =1;sprintf(xname,"truemom2#pi_{#rho0} [GeV]");}
    if (strcmp(thevar,"nrecorho") == 0) { myvar = nrecorho; COMB = 1; sprintf(xname,"N_{#rho}");}
    if (strcmp(thevar,"mrho") == 0) {myvar = mrho; MX =1;sprintf(xname,"m_{#rho} [GeV]");}
    if (strcmp(thevar,"resomrho") == 0) {myvar = mrho-mxhadgen; MX =1;sprintf(xname,"m_{#rho}-m^{true}_{#rho} [GeV]");}
    if (strcmp(thevar,"mm2rho") == 0) {myvar = mm2rho; MM2 =1; MM2RHO=1; sprintf(xname,"mm2_{#rho} [GeV^{2}]");}
    if (strcmp(thevar,"mompirho") == 0) {myvar = mompirho; MOM1 =1;sprintf(xname,"mom#pi_{#rho} [GeV]");}
    if (strcmp(thevar,"mompi0rho") == 0) {myvar = mompi0rho; MOM2 =1;sprintf(xname,"mom#pi0_{#rho} [GeV]");}
    if (strcmp(thevar,"truemompirho") == 0) {myvar = truemompirho; MOM1 =1;sprintf(xname,"truemom#pi_{#rho} [GeV]");}
    if (strcmp(thevar,"truemompi0rho") == 0) {myvar = truemompi0rho; MOM2 =1;sprintf(xname,"truemom#pi0_{#rho} [GeV]");}
    if (strcmp(thevar,"nrecoomega") == 0) { myvar = nrecoomega; COMB = 1; sprintf(xname,"N_{#omega}");}
    if (strcmp(thevar,"momega") == 0) {myvar = momega; MX =1;sprintf(xname,"m_{#omega} [GeV]");}
    if (strcmp(thevar,"resomomega") == 0) {myvar = momega-mxhadgen; MX =1;sprintf(xname,"m_{#omega}-m^{true}_{#omega} [GeV]");}
    if (strcmp(thevar,"mm2omega") == 0) {myvar = mm2omega; MM2 =1; MM2OMEGA=1; sprintf(xname,"mm2_{#omega} [GeV^{2}]");}
    if (strcmp(thevar,"mom1piome") == 0) {myvar = mom1piome; MOM1 =1;sprintf(xname,"mom1#pi_{#omega} [GeV]");}
    if (strcmp(thevar,"mom2piome") == 0) {myvar = mom2piome; MOM2 =1;sprintf(xname,"mom2#pi_{#omega} [GeV]");}
    if (strcmp(thevar,"truemom1piome") == 0) {myvar = truemom1piome; MOM1 =1;sprintf(xname,"truemom1#pi_{#omega} [GeV]");}
    if (strcmp(thevar,"truemom2piome") == 0) {myvar = truemom2piome; MOM2 =1;sprintf(xname,"truemom2#pi_{#omega} [GeV]");}
    if (strcmp(thevar,"costhome") == 0) {myvar = costhome;sprintf(xname,"cos(#theta) ");}
    if (strcmp(thevar,"meta") == 0) {myvar = barembestEta; MX =1;sprintf(xname,"m_{#eta} [GeV]");}
    if (strcmp(thevar,"metagg") == 0) {myvar = barembestEta; MX =1; DECAYMODE=(modeEta[indexbestEta]==1); sprintf(xname,"m_{#eta} [GeV]");}
    if (strcmp(thevar,"metapipipi0") == 0) {myvar = barembestEta; MX =1; DECAYMODE=(modeEta[indexbestEta]==2); sprintf(xname,"m_{#eta} [GeV]");}
    if (strcmp(thevar,"metapi0pi0pi0") == 0) {myvar = barembestEta; MX =1; DECAYMODE=(modeEta[indexbestEta]==3); sprintf(xname,"m_{#eta} [GeV]");}
    if (strcmp(thevar,"mm2eta") == 0) {myvar = mm2bestEta; MM2 =1; sprintf(xname,"m^{2}_{miss}(#eta) [GeV^{2}]");}
    if (strcmp(thevar,"mm2etagg") == 0) {myvar = mm2bestEta; MM2 =1; DECAYMODE=(modeEta[indexbestEta]==1); sprintf(xname,"m^{2}_{miss}(#eta) [GeV^{2}]");}
    if (strcmp(thevar,"mm2etapipipi0") == 0) {myvar = mm2bestEta; MM2 =1; DECAYMODE=(modeEta[indexbestEta]==2); sprintf(xname,"m^{2}_{miss}(#eta) [GeV^{2}]");}
    if (strcmp(thevar,"mm2etapi0pi0pi0") == 0) {myvar = mm2bestEta; MM2 =1; DECAYMODE=(modeEta[indexbestEta]==3); sprintf(xname,"m^{2}_{miss}(#eta) [GeV^{2}]");}
    if (strcmp(thevar,"mm2etapi0hyp") == 0) {myvar = mm2bestPi0; MM2CROSS=1; MM2PI0 = 1; sprintf(xname,"m_{miss}^{2}(#eta with #pi^{0} hypothesis) [GeV^{2}]");}
    if (strcmp(thevar,"mm2etapi0hypgg") == 0) {myvar = mm2bestPi0; MM2CROSS=1; MM2PI0 = 1; DECAYMODE=(modeEta[indexbestEta]==1); sprintf(xname,"m_{miss}^{2}(#eta, #pi^{0} hypothesis) [GeV^{2}]");}
    if (strcmp(thevar,"mm2etapi0hyppipipi0") == 0) {myvar = mm2bestPi0; MM2CROSS=1; MM2PI0 = 1; DECAYMODE=(modeEta[indexbestEta]==2); sprintf(xname,"m_{miss}^{2}(#eta with #pi^{0} hypothesis) [GeV^{2}]");}
    if (strcmp(thevar,"mm2etapi0hyppi0pi0pi0") == 0) {myvar = mm2bestPi0; MM2CROSS=1; MM2PI0 = 1; DECAYMODE=(modeEta[indexbestEta]==3); sprintf(xname,"m_{miss}^{2}(#eta with #pi^{0} hypothesis) [GeV^{2}]");}
    if (strcmp(thevar,"modeeta") == 0) {myvar = modeEta[indexbestEta];sprintf(xname,"#eta decay mode");}
    if (strcmp(thevar,"resometa") == 0) {myvar = barembestEta-mxhadgen; MX =1;sprintf(xname,"m_{#eta}-m^{true}_{#eta} [GeV]");}
    if (strcmp(thevar,"metap") == 0) {myvar = barembestEtap; MX =1;sprintf(xname,"m_{#eta'} [GeV]");}
    if (strcmp(thevar,"metaprho0g") == 0) {myvar = barembestEtap; MX =1; DECAYMODE=(modeEtap[indexbestEtap]==1); sprintf(xname,"m_{#eta'} [GeV]");}
    if (strcmp(thevar,"metapetagg") == 0) {myvar = barembestEtap; MX =1; DECAYMODE=(modeEtap[indexbestEtap]==2); sprintf(xname,"m_{#eta'} [GeV]");}
    if (strcmp(thevar,"metapetapipipi0") == 0) {myvar = barembestEtap; MX =1; DECAYMODE=(modeEtap[indexbestEtap]==3); sprintf(xname,"m_{#eta'} [GeV]");}
    if (strcmp(thevar,"metapetapi0pi0pi0") == 0) {myvar = barembestEtap; MX =1; DECAYMODE=(modeEtap[indexbestEtap]==4); sprintf(xname,"m_{#eta'} [GeV]");}
    if (strcmp(thevar,"mm2etap") == 0) {myvar = mm2bestEtap; MM2 =1;sprintf(xname,"m^{2}_{miss}(#eta') [GeV^{2}]");}
    if (strcmp(thevar,"mm2etaprho0g") == 0) {myvar = mm2bestEtap; MM2 =1; DECAYMODE=(modeEtap[indexbestEtap]==1); sprintf(xname,"m^{2}_{miss}(#eta') [GeV^{2}]");}
    if (strcmp(thevar,"mm2etapetagg") == 0) {myvar = mm2bestEtap; MM2 =1; DECAYMODE=(modeEtap[indexbestEtap]==2); sprintf(xname,"m^{2}_{miss}(#eta') [GeV^{2}]");}
    if (strcmp(thevar,"mm2etapetapipipi0") == 0) {myvar = mm2bestEtap; MM2 =1; DECAYMODE=(modeEtap[indexbestEtap]==3); sprintf(xname,"m^{2}_{miss}(#eta') [GeV^{2}]");}
    if (strcmp(thevar,"mm2etapetapi0pi0pi0") == 0) {myvar = mm2bestEtap; MM2 =1; DECAYMODE=(modeEtap[indexbestEtap]==4); sprintf(xname,"m^{2}_{miss}(#eta') [GeV^{2}]");}
    if (strcmp(thevar,"modeetap") == 0) {myvar = modeEtap[indexbestEtap];sprintf(xname,"#eta' decay mode");}
    if (strcmp(thevar,"resometap") == 0) {myvar = barembestEtap-mxhadgen; MX =1;sprintf(xname,"m_{#eta'}-m^{true}_{#eta'} [GeV]");}
    if (strcmp(thevar,"momrho0ph") == 0) {myvar = GammamomdauEtap[indexbestEtap];DAUGAMMOM = 1; sprintf(xname,"mom#gamma_{#eta'}[GeV]");}
    if (strcmp(thevar,"etapmassetadau") == 0) {myvar = EtamassdauEtap[indexbestEtap];DAUETA=1; sprintf(xname,"m_{#eta}#eta' dau[GeV]");}
    if (strcmp(thevar,"etapmassrho0dau") == 0) {myvar = Rho0massdauEtap[indexbestEtap];DAURHO=1; sprintf(xname,"m_{#rho^{0}}#eta' dau[GeV]");}
    if (strcmp(thevar,"ma0p") == 0) {myvar = ma0p; MX =1;sprintf(xname,"m_{a_{0}^{#pm}}} [GeV]");}
    if (strcmp(thevar,"mm2a0p") == 0) {myvar = mm2a0p; MM2 =1;sprintf(xname,"mm2_{a_{0}^{#pm}} [GeV^{2}]");}
    if (strcmp(thevar,"resoma0p") == 0) {myvar = ma0p-mxhadgen; MX =1;sprintf(xname,"m_{a_{0}^{#pm}-m^{true}_{a_{0}^{#pm'}}  [GeV]");}
    if (strcmp(thevar,"a0massetadau") == 0) {myvar = a0massetadau;DAUETA=1;sprintf(xname,"m_{#eta}a_{0} dau[GeV]");}
    if (strcmp(thevar,"modea0") == 0) {myvar = modea0;sprintf(xname,"mode_{a_{0}^{0}}");}
    if (strcmp(thevar,"ma0") == 0) {myvar = ma0; MX =1;sprintf(xname,"m_{a_{0}'} [GeV]");}
    if (strcmp(thevar,"mm2a0") == 0) {myvar = mm2a0; MM2 =1;sprintf(xname,"mm2_{a_{0}'} [GeV^{2}]");}
    if (strcmp(thevar,"resoma0") == 0) {myvar = ma0-mxhadgen; MX =1;sprintf(xname,"m_{a_{0}'}-m^{true}_{a_{0}'} [GeV]");}
    if (strcmp(thevar,"a0pmassetadau") == 0) {myvar = a0pmassetadau;DAUETA=1;sprintf(xname,"m_{#eta}a_{0}^{+-} dau[GeV]");}
    if (strcmp(thevar,"modea0p") == 0) {myvar = modea0p;sprintf(xname,"mode_{a_{0}^{+-}}");}
    
    if(PCMS && NLE && FLAV && BCH && IPUR && OTHER && DECAYMODE) {
      
      sprintf(le, "h");
      
      // assign flavor category
      cutflav = 5;
      if(TMath::Abs(brecocharge)) cutflav = 3;
      if(TMath::Abs(brecocharge)==0 && flav==0) cutflav = 4;	  
      
      id = hist(myvar);
      idflav = hist(myvar) + cutflav * 10000;    
      group = 200 + CATEGORY * 1000 + 100000;
      sprintf(name, "%s%d",le,group+idflav);
      if (id!=-1) ((TH1D*)gDirectory->Get(name))->Fill(mes);  
      if((!KSELE) && MM2 && CH && NCHG && !WdeltaCut && MX && MOM1 && MOM2 && COMB && PIZ && DALITZ && MM2CROSS && DAUETA && DAURHO && DAUGAMMOM){  	
	group = 200 + CATEGORY * 1000;
	sprintf(name, "%s%d",le,group+idflav);
	if (id!=-1) ((TH1D*)gDirectory->Get(name))->Fill(mes); 
      }
    }
    
    if(PCMS && NLE && !flav && BCH && IPUR && OTHER && DECAYMODE){  
      sprintf(le, "h");
      group = 90206 + CATEGORY * 1000 + 100000;
      sprintf(name, "%s%d",le,group);
      if (id!=-1) ((TH1D*)gDirectory->Get(name))->Fill(mes);  
      if((!KSELE) && MM2 && CH && NCHG && !WdeltaCut && MX && MOM1 && MOM2 && COMB && PIZ && DALITZ && MM2PI0 && MM2ETA && MM2RHO && MM2RHO0 && MM2OMEGA && DAUETA && DAUGAMMOM){  
	group = 90206 + CATEGORY * 1000;
	sprintf(name, "%s%d",le,group);
	if (id!=-1) ((TH1D*)gDirectory->Get(name))->Fill(mes);  
      }
    }
    
  }
   


}
// ----------------------------------------------------------------------
double theanal::getBsysweight(int decType,int thevub) {
  double theweight;
  theweight = 1.;  
  theweight = Bsem->weight(decType); 
  if(thevub) theweight = 1.;
  return theweight;
}
// ----------------------------------------------------------------------
double theanal::getDsysweight(int decDpi,int decDk,int decDks,int decDpiz,int decDlep,int decImode,int thevub) {
  double theweight;
  theweight = 1.;  
  int bla(1);
  theweight = Dvar->weight(decDpi,decDk,decDks,decDpiz,decDlep,decImode,bla);    
  if(thevub) theweight = 1.;
  return theweight;
}

double theanal::FermiWeight(double kp, double deltamb, double deltaa){

  Double_t BMASS, BQMASS, A0;
  BMASS   = 5.2792;
  BQMASS  = 4.800;
  A0      = 1.29;
  double mb1 = BQMASS + deltamb;
  
  double a1 = A0 + deltaa ;
  double w81;

  Int_t nIntegral = 10000;
  const double kmin = -5.;
  const double kmax =  1.;
  double Nold = 0;
  double Nnew1(0);

  // -- Compute reweighting ratios
  for (int i = 0; i < nIntegral; ++i) {
    double kplus = (i+0.5)/((double)nIntegral)*(kmax-kmin)+kmin;
    Nold += (kmax-kmin)/((double)nIntegral)*fermi(kplus, BQMASS, A0);
    Nnew1 += (kmax-kmin)/((double)nIntegral)*fermi(kplus, mb1, a1);
  }
  double kpold = kp; 
  double kpnew1 = 5.279 - mb1 - (5.279 - BQMASS - kp);
  
  double fold = fermi(kpold, BQMASS, A0)/Nold; 
  double fnew1 = fermi(kpnew1, mb1, a1)/Nnew1;
  if(fold<=0)cout<<"?? "<<fold<<" "<<kpold<<" "<<BQMASS<<" "<<A0<<endl;
  w81 = fnew1/fold; 
  return w81;
}

  // -- Shape function
double theanal::fermi(double kp, double m, double a) {
  double BMASS   = 5.2792;
  double x = kp/(BMASS - m);
  if ((kp>-m) && (x <= 1.)) {
    return TMath::Power((1-x), a) * TMath::Exp((1+a)*x); 
  } 
  return 0.;
}


void theanal::Bookhist()
{
  fHistFile->cd();
  TH1 *h;
  char name[100], title[100], number[100];
  int lo; 
  sprintf(name, "h100000");  sprintf(title, "%s%s", thevar, " signal events  after all cuts: enriched");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
  sprintf(name, "h200000");  sprintf(title, "%s%s", thevar, " vub events after all cuts: enriched");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
  sprintf(name, "h300000");  sprintf(title, "%s%s", thevar, " vcb events  after all cuts: enriched");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
  sprintf(name, "h400000");  sprintf(title, "%s%s", thevar, " other events after all cuts: enriched");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
  sprintf(name, "h91206");  sprintf(title, "mes signal after all cuts: enriched");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "h92206");  sprintf(title, "mes vub after all cuts: enriched");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "h93206");  sprintf(title, "mes vcb after all cuts: enriched");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "h94206");  sprintf(title, "mes other after all cuts: enriched");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "h1100000");  sprintf(title, "%s%s", thevar, " signal events after lepton cuts: enriched");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
  sprintf(name, "h1200000");  sprintf(title, "%s%s", thevar, " vub events after lepton cuts: enriched");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
  sprintf(name, "h1300000");  sprintf(title, "%s%s", thevar, " vcb events after lepton cuts: enriched");  h = new TH1D(name, title, thebins, themin, themax  );  h->Sumw2();
  sprintf(name, "h1400000");  sprintf(title, "%s%s", thevar, " other events after lepton cuts: enriched");  h = new TH1D(name, title, thebins, themin, themax );  h->Sumw2();
  sprintf(name, "h191206");  sprintf(title, "mes signal after lepton cuts: enriched");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "h192206");  sprintf(title, "mes vub after lepton cuts: enriched");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "h193206");  sprintf(title, "mes vcb after lepton cuts: enriched");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
  sprintf(name, "h194206");  sprintf(title, "mes other after lepton cuts: enriched");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();


  sprintf(name, "h7100000");  sprintf(title, "%s%s", thevar, " all component");  h = new TH1D(name, title, thebins, themin, themax  );  
  sprintf(name, "h7200000");  sprintf(title, "%s%s", thevar, " all comp - signal");  h = new TH1D(name, title, thebins, themin, themax );  
  sprintf(name, "h7300000");  sprintf(title, "%s%s", thevar, " all comp - signal - vub");  h = new TH1D(name, title, thebins, themin, themax  ); 
  sprintf(name, "h7400000");  sprintf(title, "%s%s", thevar, " other");  h = new TH1D(name, title, thebins, themin, themax );

  sprintf(name, "h8100000");  sprintf(title, "%s%s", thevar, " all vub component");  h = new TH1D(name, title, thebins, themin, themax  );  
  sprintf(name, "h8200000");  sprintf(title, "%s%s", thevar, " vub");  h = new TH1D(name, title, thebins, themin, themax );

  sprintf(name, "h17100000");  sprintf(title, "%s%s", thevar, " all component");  h = new TH1D(name, title, thebins, themin, themax  );  
  sprintf(name, "h17200000");  sprintf(title, "%s%s", thevar, " all comp - signal");  h = new TH1D(name, title, thebins, themin, themax );  
  sprintf(name, "h17300000");  sprintf(title, "%s%s", thevar, " all comp - signal - vub");  h = new TH1D(name, title, thebins, themin, themax  ); 
  sprintf(name, "h17400000");  sprintf(title, "%s%s", thevar, " other");  h = new TH1D(name, title, thebins, themin, themax );

  sprintf(name, "h18100000");  sprintf(title, "%s%s", thevar, " all vub component");  h = new TH1D(name, title, thebins, themin, themax  );  
  sprintf(name, "h18200000");  sprintf(title, "%s%s", thevar, " vub");  h = new TH1D(name, title, thebins, themin, themax );

  for (lo=1;lo<thebins+1;lo++) {
    
    sprintf(number, "%d",  200+lo);     
    sprintf(name,"%s%s" , "h1",number);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);         h->Sumw2();
    sprintf(name,"%s%s" , "h31",number);  sprintf(title, "mes data bch");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h41",number);  sprintf(title, "mes data b0 os");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h51",number);  sprintf(title, "mes data b0 ss");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h2",number);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);         h->Sumw2();
    sprintf(name,"%s%s" , "h32",number);  sprintf(title, "mes data bch");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h42",number);  sprintf(title, "mes data b0 os");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h52",number);  sprintf(title, "mes data b0 ss");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h3",number);  sprintf(title, "mes MC");  h = new TH1D(name, title, 40, 5.2, 5.3);         h->Sumw2();
    sprintf(name,"%s%s" , "h33",number);  sprintf(title, "mes MC bch");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h43",number);  sprintf(title, "mes MC b0 os");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h53",number);  sprintf(title, "mes MC b0 ss");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h4",number);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);         h->Sumw2();
    sprintf(name,"%s%s" , "h34",number);  sprintf(title, "mes data bch");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h44",number);  sprintf(title, "mes data b0 os");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h54",number);  sprintf(title, "mes data b0 ss");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();

    sprintf(name,"%s%s" , "h11",number);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);         h->Sumw2();
    sprintf(name,"%s%s" , "h131",number);  sprintf(title, "mes data bch");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h141",number);  sprintf(title, "mes data b0 os");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h151",number);  sprintf(title, "mes data b0 ss");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h12",number);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);         h->Sumw2();
    sprintf(name,"%s%s" , "h132",number);  sprintf(title, "mes data bch");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h142",number);  sprintf(title, "mes data b0 os");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h152",number);  sprintf(title, "mes data b0 ss");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h13",number);  sprintf(title, "mes MC");  h = new TH1D(name, title, 40, 5.2, 5.3);         h->Sumw2();
    sprintf(name,"%s%s" , "h133",number);  sprintf(title, "mes MC bch");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h143",number);  sprintf(title, "mes MC b0 os");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h153",number);  sprintf(title, "mes MC b0 ss");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h14",number);  sprintf(title, "mes data");  h = new TH1D(name, title, 40, 5.2, 5.3);         h->Sumw2();
    sprintf(name,"%s%s" , "h134",number);  sprintf(title, "mes data bch");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h144",number);  sprintf(title, "mes data b0 os");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();
    sprintf(name,"%s%s" , "h154",number);  sprintf(title, "mes data b0 ss");  h = new TH1D(name, title, 40, 5.2, 5.3);  h->Sumw2();

  }
   sprintf(name, "h8888");  sprintf(title, "test gauss");  h = new TH1D(name, title, 100, -.05, .05);    h->Sumw2(); 
}

int theanal::hist(double mx){

  // categories
  int bin = int((mx-themin)/((themax-themin)/thebins)+1);
  if (mx>themax || mx<themin) bin = -1;
  if( mx==themax) bin = thebins;
  if (mx==themin) bin = 1;
  return bin;
}

void theanal::Fitmes(int cat, int cut){

  fHistFile->cd();
  // fit to mes distribution and fill of the Mx plots...

  //  gROOT->cd();

  int group = 200 + cat * 1000;
   
  char name[100];
  char namebch[100];
  char nameb0os[100];
  char nameb0ss[100];
  char le[100];
  double sigs[1000];
  double errsigs[1000];
  double tempbin;
  double tempbinchb;
  double tempbinb0os;
  double tempbinb0ss;
  double temperr;
  double temperrchb;
  double temperrb0os;
  double temperrb0ss;
  double chid = 0.181;
  int is;
  int thegroup = group;      
  double foomean;
  double foosigma;
  double fooalpha;
  double foon;
  double usedmean;
  double usedsigma;
  double usedalpha;
  double usedn;
  int i;
  int addcut = 0; 
  if(cut) addcut = 1;
  int const theb = thebins + 1;
  sprintf(le, "h");

  //extracting the signal fit parameters 
  sprintf(name, "%s%d",le,group+90006+addcut*100000);  
  cout << name << endl;
  cout << ((TH1D*)gDirectory->Get(name))->Integral() << endl;
  double dummy1,dummy2;
  sighisto(dummy1,dummy2,(TH1D*)gDirectory->Get(name),usedmean,usedsigma,usedalpha,usedn,1,-11111111.,-1111111.,-1111111.,-1111111.,-1111111.);
  cout << "mes result for: " << name << " is, MEAN " << usedmean << " SIGMA " << usedsigma << " ALPHA " << usedalpha << " N " << usedn << endl;
  
  group = thegroup;
  for (i = 1;  i < theb; i++) { 
    is = i;
    sprintf(name, "%s%d",le,group+is+addcut*10000);  
    
    sprintf(namebch, "%s%d",le,group+is+30000+addcut*100000);  
    sprintf(nameb0os, "%s%d",le,group+is+40000+addcut*100000);  
    sprintf(nameb0ss, "%s%d",le,group+is+50000+addcut*100000);  
    
    // mixing correction      
    
    for(int k=1;k<40;k++){
      tempbinchb = ((TH1D*)gDirectory->Get(namebch))->GetBinContent(k);
      tempbinb0os = ((TH1D*)gDirectory->Get(nameb0os))->GetBinContent(k);
      tempbinb0ss = ((TH1D*)gDirectory->Get(nameb0ss))->GetBinContent(k);
      temperrchb = ((TH1D*)gDirectory->Get(namebch))->GetBinError(k);
      temperrb0os = ((TH1D*)gDirectory->Get(nameb0os))->GetBinError(k);
      temperrb0ss = ((TH1D*)gDirectory->Get(nameb0ss))->GetBinError(k);
      tempbin = tempbinchb + ((1-chid)/(1-2*chid)) * tempbinb0os - (chid/(1-2*chid)) * tempbinb0ss;
      temperr = sqrt(temperrchb*temperrchb + ((1-chid)/(1-2*chid)) * ((1-chid)/(1-2*chid)) * tempbinb0os + (chid/(1-2*chid)) * (chid/(1-2*chid))* tempbinb0ss);
      ((TH1D*)gDirectory->Get(name))->SetBinContent(k, tempbin);
      ((TH1D*)gDirectory->Get(name))->SetBinError(k, temperr); 
    }	
    
    // put the result of the fit in each Mx bin
    
    sighisto(sigs[is-1],errsigs[is-1],(TH1D*)gDirectory->Get(name),foomean,foosigma,fooalpha,foon,1,usedmean,usedsigma,usedalpha,usedn,-1111111);
    if(!strcmp(thevar,"mes")){      
      sigs[is-1] =((TH1D*)gDirectory->Get(name))->Integral();
      errsigs[is-1] = sqrt(((TH1D*)gDirectory->Get(name))->Integral());     
    }
    int title = cat * 100000 + addcut * 1000000;
    sprintf(name, "%s%d",le,title);  
    ((TH1D*)gDirectory->Get(name))->SetBinContent(is, sigs[is-1]);
    ((TH1D*)gDirectory->Get(name))->SetBinError(is, errsigs[is-1]);
  }
}

void 
theanal::sighisto(double&signal, double& signalErr, TH1D *histo, double &resmean, double &ressigma, double &resalpha, double &resn, int fixpar, double mean, double sigma, double alpha,  double n, double argus){

  recoilAnalysis b;
  mesData themes;
  double thesigma = -1111111;
  double themean = -1111111;
  double theargus = -1111111;
  double thealpha = -1111111;
  double then = -1111111;
  if (fixpar){
    themean = mean;
    thesigma = sigma;
    thealpha = alpha;
    then = n;
    mesData tempthemes(*(b.vubMes(histo, resmean, ressigma, resalpha, resn, 1, 1, themean, thesigma, thealpha, then, theargus)));
    themes = tempthemes;
  }else{
    mesData tempthemes(*(b.vubMes(histo, resmean, ressigma, resalpha, resn, 1,1, themean, thesigma, thealpha, then, theargus)));
    themes = tempthemes;
  }
  signal=themes.theSig(); signalErr=themes.theErrSig(); 

}

void theanal::overlap(int typ, double relnorm, TString dir){
   
//    //plot in mX with all components

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(11111111);
  fHistFile->cd();
  char addit[100];
  sprintf(addit,"h");
  if(typ==1) sprintf(addit,"h1");
  char name[100], name2[100];
  int const theb = thebins + 1;
  sprintf(name, "%s%s",addit,"7400000");
  sprintf(name2, "%s%s",addit,"400000");
  int y;
  for(y=1;y<theb;y++){
    ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*relnorm);    
  }
  
  sprintf(name, "%s%s",addit,"7300000");
  sprintf(name2, "%s%s",addit,"300000");
  for(y=1;y<theb;y++){
    ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*relnorm);
  }
  sprintf(name2, "%s%s",addit,"400000");
  for(y=1;y<theb;y++){
    double it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
    ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*relnorm);
  }
  
  sprintf(name, "%s%s",addit,"7200000");
  sprintf(name2, "%s%s",addit,"200000");
  for(y=1;y<theb;y++){
    ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y));
  }
  sprintf(name2, "%s%s",addit,"300000");
  for(y=1;y<theb;y++){
    double it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
    ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*relnorm);
  }
  sprintf(name2, "%s%s",addit,"400000");
  for(y=1;y<theb;y++){
    double it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
    ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*relnorm);
  }

  sprintf(name, "%s%s",addit,"7100000");
  sprintf(name2, "%s%s",addit,"100000");
  for(y=1;y<theb;y++){
    ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y));
  }
  sprintf(name2, "%s%s",addit,"200000");
  for(y=1;y<theb;y++){
    double it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
    ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y));
  }
  sprintf(name2, "%s%s",addit,"300000");
  for(y=1;y<theb;y++){
    double it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
    ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*relnorm);
  }
  sprintf(name2, "%s%s",addit,"400000");
  for(y=1;y<theb;y++){
    double it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
    ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y)*relnorm);
  }

   c1 = new TCanvas("c1"," ",200,10,800,800);   
   c1->Clear();

   sprintf(name, "%s%s",addit,"7100000");
   ((TH1D*)gDirectory->Get(name))->SetLineWidth(3.);
   ((TH1D*)gDirectory->Get(name))->SetLineColor(kBlack);
   ((TH1D*)gDirectory->Get(name))->SetFillColor(kWhite);
   ((TH1D*)gDirectory->Get(name))->SetStats(0);
   themax = ((TH1D*)gDirectory->Get(name))->GetMaximum() * 1.2; 
   ((TH1D*)gDirectory->Get(name))->SetMaximum(themax); 
   ((TH1D*)gDirectory->Get(name))->SetTitle("");
   ((TH1D*)gDirectory->Get(name))->SetXTitle(xname);
   ((TH1D*)gDirectory->Get(name))->Draw();
   sprintf(name, "%s%s",addit,"7200000");
   ((TH1D*)gDirectory->Get(name))->SetFillColor(38);
   ((TH1D*)gDirectory->Get(name))->SetStats(0);
   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
   sprintf(name, "%s%s",addit,"7300000");
   ((TH1D*)gDirectory->Get(name))->SetFillColor(kYellow);
   ((TH1D*)gDirectory->Get(name))->SetStats(0);
   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
   sprintf(name, "%s%s",addit,"7400000");
   ((TH1D*)gDirectory->Get(name))->SetFillColor(13);
   ((TH1D*)gDirectory->Get(name))->SetStats(0);
   ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
   TLegendEntry *legge19; 
   TLegend *leg19;
   //leg19 = new TLegend(0.13,0.45,0.35,0.85);
   leg19 = new TLegend(0.58,0.48,0.88,0.88);
    leg19->SetFillStyle(0); leg19->SetBorderSize(0.); leg19->SetTextSize(0.035); 
   leg19->SetFillColor(0); 
   sprintf(name, "%s%s",addit,"7100000");   
   legge19 = leg19->AddEntry(((TH1D*)gDirectory->Get(name)), "signal", "f"); 
   sprintf(name, "%s%s",addit,"7200000");   
   legge19 = leg19->AddEntry(((TH1D*)gDirectory->Get(name)), "other b #rightarrow u l #nu", "f"); 
   sprintf(name, "%s%s",addit,"7300000");   
   legge19 = leg19->AddEntry(((TH1D*)gDirectory->Get(name)), "b #rightarrow c l #nu", "f"); 
   sprintf(name, "%s%s",addit,"7400000");   
   legge19 = leg19->AddEntry(((TH1D*)gDirectory->Get(name)), "other", "f"); 
   leg19->Draw();
   sprintf(name, "%s%s%d%s", dir.Data(), thevar,typ,"-allcomp.eps");

   c1->SaveAs(name);  

  sprintf(name, "%s%s",addit,"8100000");
  sprintf(name2, "%s%s",addit,"100000");
  for(y=1;y<theb;y++){
    ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y));
  }
  sprintf(name2, "%s%s",addit,"200000");
  for(y=1;y<theb;y++){
    double it = ((TH1D*)gDirectory->Get(name))->GetBinContent(y);
    ((TH1D*)gDirectory->Get(name))->SetBinContent(y,it + ((TH1D*)gDirectory->Get(name2))->GetBinContent(y));
  }
  

  sprintf(name, "%s%s",addit,"8200000");
  sprintf(name2, "%s%s",addit,"200000");
  for(int y=1;y<theb;y++){
    ((TH1D*)gDirectory->Get(name))->SetBinContent(y,((TH1D*)gDirectory->Get(name2))->GetBinContent(y));
  }
  
  c1 = new TCanvas("c1"," ",200,10,800,800);   
  c1->Clear();
  
  sprintf(name, "%s%s",addit,"8100000");
  ((TH1D*)gDirectory->Get(name))->SetLineWidth(3.);
  ((TH1D*)gDirectory->Get(name))->SetLineColor(kBlack);
  ((TH1D*)gDirectory->Get(name))->SetFillColor(kWhite);
  ((TH1D*)gDirectory->Get(name))->SetStats(0);
  themax = ((TH1D*)gDirectory->Get(name))->GetMaximum() * 1.2; 
  ((TH1D*)gDirectory->Get(name))->SetMaximum(themax); 
  ((TH1D*)gDirectory->Get(name))->SetTitle("");
  ((TH1D*)gDirectory->Get(name))->SetXTitle(xname);
  ((TH1D*)gDirectory->Get(name))->Draw();
  ((TH1D*)gDirectory->Get(name))->DrawCopy("same");  
  sprintf(name, "%s%s",addit,"8200000");
  ((TH1D*)gDirectory->Get(name))->SetFillColor(38);
  ((TH1D*)gDirectory->Get(name))->SetStats(0);
  ((TH1D*)gDirectory->Get(name))->DrawCopy("same");
  TLegendEntry *legge18; 
  TLegend *leg18;
  //  leg18 = new TLegend(0.13,0.45,0.35,0.85);
  leg18 = new TLegend(0.58,0.68,0.88,0.88);
  leg18->SetFillStyle(0); leg18->SetBorderSize(0.); leg18->SetTextSize(0.035); 
  leg18->SetFillColor(0); 
  sprintf(name, "%s%s",addit,"8100000");   
  legge18 = leg18->AddEntry(((TH1D*)gDirectory->Get(name)), "signal", "f"); 
  sprintf(name, "%s%s",addit,"8200000");   
  legge18 = leg18->AddEntry(((TH1D*)gDirectory->Get(name)), "other b #rightarrow u l #nu", "f"); 
  leg18->Draw();

  sprintf(name, "%s%s%d%s", dir.Data(), thevar,typ,"-vubcomp.eps");  
  c1->SaveAs(name);  

  c1 = new TCanvas("c1"," ",200,10,800,800);   
  c1->Clear();
  
  sprintf(name, "%s%s",addit,"100000");
  ((TH1D*)gDirectory->Get(name))->SetTitle("");
  ((TH1D*)gDirectory->Get(name))->SetXTitle(xname);
  ((TH1D*)gDirectory->Get(name))->SetMarkerStyle(20);
  ((TH1D*)gDirectory->Get(name))->SetMarkerSize(1.);	
  ((TH1D*)gDirectory->Get(name))->SetStats(1);
  themax = ((TH1D*)gDirectory->Get(name))->GetMaximum() * 1.2;  
  ((TH1D*)gDirectory->Get(name))->SetMaximum(themax);  
  ((TH1D*)gDirectory->Get(name))->Draw();

  sprintf(name, "%s%s%d%s", dir.Data(), thevar,typ,"-signal.eps");  
  c1->SaveAs(name);  
 
  
  c1 = new TCanvas("c1"," ",200,10,800,800);   
  c1->Clear();
  
  sprintf(name, "%s%s",addit,"200000");
  ((TH1D*)gDirectory->Get(name))->SetTitle("");
  ((TH1D*)gDirectory->Get(name))->SetXTitle(xname);
  ((TH1D*)gDirectory->Get(name))->SetMarkerStyle(20);
  ((TH1D*)gDirectory->Get(name))->SetMarkerSize(1.);	
  ((TH1D*)gDirectory->Get(name))->SetStats(11);
  ((TH1D*)gDirectory->Get(name))->Draw();

  sprintf(name, "%s%s%d%s", dir.Data(), thevar,typ,"-vub.eps");  
  c1->SaveAs(name);  
  
  c1 = new TCanvas("c1"," ",200,10,800,800);   
  c1->Clear();
  
  sprintf(name, "%s%s",addit,"300000");
  ((TH1D*)gDirectory->Get(name))->SetTitle("");
  ((TH1D*)gDirectory->Get(name))->SetXTitle(xname);
  ((TH1D*)gDirectory->Get(name))->SetMarkerStyle(20);
  ((TH1D*)gDirectory->Get(name))->SetMarkerSize(1.);	
  ((TH1D*)gDirectory->Get(name))->SetStats(1);
  ((TH1D*)gDirectory->Get(name))->Draw();

  sprintf(name, "%s%s%d%s", dir.Data(), thevar,typ,"-vcb.eps");  
  c1->SaveAs(name);  
  
  sprintf(name, "%s%s",addit,"400000");
  ((TH1D*)gDirectory->Get(name))->SetTitle("");
  ((TH1D*)gDirectory->Get(name))->SetXTitle(xname);
  ((TH1D*)gDirectory->Get(name))->SetMarkerStyle(20);
  ((TH1D*)gDirectory->Get(name))->SetMarkerSize(1.);	
  ((TH1D*)gDirectory->Get(name))->SetStats(1);
  ((TH1D*)gDirectory->Get(name))->Draw();

  sprintf(name, "%s%s%d%s", dir.Data(), thevar,typ,"-other.eps");  
  c1->SaveAs(name);  
  
}

double theanal::chisq(TH1 *h1, TH1 *h2){
  double chisq = 0; 
  double tempchisq = 0;
  double inter;
  double temperr1;
  double temperr2;
  inter = h1->Integral();
  double inter2;
  inter2 = h2->Integral();
  int nbins = h1->GetNbinsX();
  for(int i=1;i<nbins;i++) {
    tempchisq = 0;
    temperr1 = 1;
    temperr2 = 1;
    if (h1->GetBinError(i)>.9) temperr1 = h1->GetBinError(i);
    if (h2->GetBinError(i)>.9) temperr2 = h2->GetBinError(i);
    if (((h1->GetBinError(i))*(h1->GetBinError(i))+(h2->GetBinError(i))*(h2->GetBinError(i)))>0) tempchisq = ((h1->GetBinContent(i)) - (h2->GetBinContent(i))*inter/inter2) / sqrt((temperr1*temperr1)+(temperr2*temperr2)*(inter/inter2)*(inter/inter2));  
    chisq = chisq + tempchisq*tempchisq;
  }
  chisq = chisq / nbins;
  return chisq;
}


// ----------------------------------------------------------------------
Double_t theanal::smeargauss(double invalue, double mean, double sigma){    

  // smears using a gaussian distribution

  double returnvalue=0.;
  double xsmear = gRandom->Gaus(SHIFTNEUT, SIGMANEUT); 
  returnvalue = invalue + xsmear;
  //cout << sigma << mean << " " << invalue << " " << returnvalue-mean-invalue << endl;
  ((TH1D*)gDirectory->Get("h8888"))->Fill(returnvalue - invalue);
  return returnvalue;
}

void theanal::Init(TTree *tree)
{
//   Set branch addresses
   if (tree == 0) return;
   fChain    = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("mes",&mes);
   fChain->SetBranchAddress("de",&de);
   fChain->SetBranchAddress("pur",&pur);
   fChain->SetBranchAddress("intpur",&intpur);
   fChain->SetBranchAddress("mode",&mode);
   fChain->SetBranchAddress("brecoflav",&brecoflav);
   fChain->SetBranchAddress("brecocharge",&brecocharge);
   fChain->SetBranchAddress("mxhad",&mxhad);
   fChain->SetBranchAddress("mxhadgen",&mxhadgen);
   fChain->SetBranchAddress("q2",&q2);
   fChain->SetBranchAddress("q2fit",&q2fit);
   fChain->SetBranchAddress("q2gen",&q2gen);
   fChain->SetBranchAddress("pcmsgen",&pcmsgen);
   fChain->SetBranchAddress("vub",&vub);
   fChain->SetBranchAddress("vcb",&vcb);
   fChain->SetBranchAddress("other",&other);
   fChain->SetBranchAddress("xcharge",&xcharge);
   fChain->SetBranchAddress("mxhad",&mxhad);
   fChain->SetBranchAddress("mxhadfit",&mxhadfit);
   fChain->SetBranchAddress("mxhadchg",&mxhadchg);
   fChain->SetBranchAddress("lcharge",&lcharge);
   fChain->SetBranchAddress("plab",&plab);
   fChain->SetBranchAddress("tlab",&tlab);
   fChain->SetBranchAddress("flab",&flab);
   fChain->SetBranchAddress("pcms",&pcms);
   fChain->SetBranchAddress("tcms",&tcms);
   fChain->SetBranchAddress("fcms",&fcms);
   fChain->SetBranchAddress("ecms",&ecms);
   fChain->SetBranchAddress("nle",&nle);
   fChain->SetBranchAddress("nlept500",&nlept500);
   fChain->SetBranchAddress("nel",&nel);
   fChain->SetBranchAddress("nmu",&nmu);
   fChain->SetBranchAddress("isele",&isele);
   fChain->SetBranchAddress("nchg",&nchg);
   fChain->SetBranchAddress("nneu",&nneu);
   fChain->SetBranchAddress("nkp",&nkp);
   fChain->SetBranchAddress("nks",&nks);
   fChain->SetBranchAddress("mm2",&mm2);
   fChain->SetBranchAddress("wdeltam",&wdeltam);
   fChain->SetBranchAddress("Gvxbtyp",&Gvxbtyp);
   fChain->SetBranchAddress("nrecoPi",&nrecoPi);
   fChain->SetBranchAddress("indexbestPi",&indexbestPi);
   fChain->SetBranchAddress("mm2bestPi",&mm2bestPi);
   fChain->SetBranchAddress("Eneualt",&Eneualt);
   fChain->SetBranchAddress("baremassjpsiPi",baremassjpsiPi);
   fChain->SetBranchAddress("chbestPi",&chbestPi);
   fChain->SetBranchAddress("nrecoPi0",&nrecoPi0);
   fChain->SetBranchAddress("indexbestPi0",&indexbestPi0);
   fChain->SetBranchAddress("barembestPi0",&barembestPi0);
   fChain->SetBranchAddress("mm2bestPi0",&mm2bestPi0);
   fChain->SetBranchAddress("mm2gamma",&mm2gamma);
   fChain->SetBranchAddress("truemom1phpi0",&truemom1phpi0);
   fChain->SetBranchAddress("truemom2phpi0",&truemom2phpi0);
   fChain->SetBranchAddress("truemomlab1phpi0",&truemomlab1phpi0);
   fChain->SetBranchAddress("truemomlab2phpi0",&truemomlab2phpi0);
   fChain->SetBranchAddress("trueth1phpi0",&trueth1phpi0);
   fChain->SetBranchAddress("trueth2phpi0",&trueth2phpi0);
   fChain->SetBranchAddress("Estar1dauPi0",Estar1dauPi0);
   fChain->SetBranchAddress("Estar2dauPi0",Estar2dauPi0);
   fChain->SetBranchAddress("truemrho",&truemrho);
   fChain->SetBranchAddress("nrecorho",&nrecorho);
   fChain->SetBranchAddress("mrho",&mrho);
   fChain->SetBranchAddress("mm2rho",&mm2rho);
   fChain->SetBranchAddress("truemompirho",&truemompirho);
   fChain->SetBranchAddress("truemompi0rho",&truemompi0rho);
   fChain->SetBranchAddress("mompirho",&mompirho);
   fChain->SetBranchAddress("mompi0rho",&mompi0rho);
   fChain->SetBranchAddress("truemrho0",&truemrho0);
   fChain->SetBranchAddress("nrecorho0",&nrecorho0);
   fChain->SetBranchAddress("mrho0",&mrho0);
   fChain->SetBranchAddress("mm2rho0",&mm2rho0);
   fChain->SetBranchAddress("truemom1pirho0",&truemom1pirho0);
   fChain->SetBranchAddress("truemom2pirho0",&truemom2pirho0);
   fChain->SetBranchAddress("mom1pirho0",&mom1pirho0);
   fChain->SetBranchAddress("mom2pirho0",&mom2pirho0);
   fChain->SetBranchAddress("truemomega",&truemomega);
   fChain->SetBranchAddress("nrecoomega",&nrecoomega);
   fChain->SetBranchAddress("momega",&momega);
   fChain->SetBranchAddress("mm2omega",&mm2omega);
   fChain->SetBranchAddress("truemom1piome",&truemom1piome);
   fChain->SetBranchAddress("truemom2piome",&truemom2piome);
   fChain->SetBranchAddress("truemompi0ome",&truemompi0ome);
   fChain->SetBranchAddress("truedalitzpi1pi2ome",&truedalitzpi1pi2ome);
   fChain->SetBranchAddress("truedalitzpi1pi0ome",&truedalitzpi1pi0ome);
   fChain->SetBranchAddress("truecosthome",&truecosthome);
   fChain->SetBranchAddress("mom1piome",&mom1piome);
   fChain->SetBranchAddress("mom2piome",&mom2piome);
   fChain->SetBranchAddress("mompi0ome",&mompi0ome);
   fChain->SetBranchAddress("dalitzpi1pi2ome",&dalitzpi1pi2ome);
   fChain->SetBranchAddress("dalitzpi1pi0ome",&dalitzpi1pi0ome);
   fChain->SetBranchAddress("costhome",&costhome);
   fChain->SetBranchAddress("nrecoEta",&nrecoEta);
   fChain->SetBranchAddress("indexbestEta",&indexbestEta);
   fChain->SetBranchAddress("barembestEta",&barembestEta);
   fChain->SetBranchAddress("mm2bestEta",&mm2bestEta);
   fChain->SetBranchAddress("modeEta",modeEta);
   fChain->SetBranchAddress("nrecoEtap",&nrecoEtap);
   fChain->SetBranchAddress("indexbestEtap",&indexbestEtap);
   fChain->SetBranchAddress("barembestEtap",&barembestEtap);
   fChain->SetBranchAddress("mm2bestEtap",&mm2bestEtap);
   fChain->SetBranchAddress("modeEtap",modeEtap);
   fChain->SetBranchAddress("EtamassdauEtap",EtamassdauEtap);
   fChain->SetBranchAddress("Rho0massdauEtap",Rho0massdauEtap);
   fChain->SetBranchAddress("GammamomdauEtap",GammamomdauEtap);
   fChain->SetBranchAddress("ma0",&ma0);
   fChain->SetBranchAddress("mm2a0",&mm2a0);
   fChain->SetBranchAddress("a0massetadau",&a0massetadau);
   fChain->SetBranchAddress("modea0",&modea0);
   fChain->SetBranchAddress("ma0p",&ma0p);
   fChain->SetBranchAddress("mm2a0p",&mm2a0p);
   fChain->SetBranchAddress("a0pmassetadau",&a0pmassetadau);
   fChain->SetBranchAddress("modea0p",&modea0p);
   Notify();
}

Bool_t theanal::Notify()
{
//   called when loading a new file
//   get branch pointers
   b_mes = fChain->GetBranch("mes");
   b_de = fChain->GetBranch("de");
   b_pur = fChain->GetBranch("pur");
   b_intpur = fChain->GetBranch("intpur");
   b_mode = fChain->GetBranch("mode");
   b_brecoflav = fChain->GetBranch("brecoflav");
   b_brecocharge = fChain->GetBranch("brecocharge");
   b_mxhadgen = fChain->GetBranch("mxhadgen");
   b_q2 = fChain->GetBranch("q2");
   b_q2fit = fChain->GetBranch("q2fit");
   b_q2gen = fChain->GetBranch("q2gen");
   b_pcmsgen = fChain->GetBranch("pcmsgen");
   b_vub = fChain->GetBranch("vub");
   b_vcb = fChain->GetBranch("vcb");
   b_other = fChain->GetBranch("other");
   b_xcharge = fChain->GetBranch("xcharge");
   b_mxhadfit = fChain->GetBranch("mxhadfit");
   b_mxhadchg = fChain->GetBranch("mxhadchg");
   b_lcharge = fChain->GetBranch("lcharge");
   b_plab = fChain->GetBranch("plab");
   b_tlab = fChain->GetBranch("tlab");
   b_flab = fChain->GetBranch("flab");
   b_pcms = fChain->GetBranch("pcms");
   b_tcms = fChain->GetBranch("tcms");
   b_fcms = fChain->GetBranch("fcms");
   b_ecms = fChain->GetBranch("ecms");
   b_nle = fChain->GetBranch("nle");
   b_nlept500 = fChain->GetBranch("nlept500");
   b_nel = fChain->GetBranch("nel");
   b_nmu = fChain->GetBranch("nmu");
   b_isele = fChain->GetBranch("isele");
   b_nchg = fChain->GetBranch("nchg");
   b_nneu = fChain->GetBranch("nneu");
   b_nkp = fChain->GetBranch("nkp");
   b_nks = fChain->GetBranch("nks");
   b_mm2 = fChain->GetBranch("mm2");
   b_wdeltam = fChain->GetBranch("wdeltam");
   b_vxbtyp = fChain->GetBranch("Gvxbtyp");
   b_nrecopi = fChain->GetBranch("nrecoPi");
   b_eneualt = fChain->GetBranch("Eneualt");
   b_baremassjpsipi = fChain->GetBranch("baremassjpsiPi");
   b_mm2bestpi = fChain->GetBranch("mm2bestPi");
   b_indexbestpi = fChain->GetBranch("indexbestPi");
   b_chbestpi = fChain->GetBranch("chbestPi");
   b_nrecopi0 = fChain->GetBranch("nrecoPi0");
   b_indexbestpi0 = fChain->GetBranch("indexbestPi0");
   b_barembestpi0 = fChain->GetBranch("barembestPi0");
   b_mm2bestpi0 = fChain->GetBranch("mm2bestPi0");
   b_mm2gamma = fChain->GetBranch("mm2gamma");
   b_truemom1phpi0 = fChain->GetBranch("truemom1phpi0");
   b_truemom2phpi0 = fChain->GetBranch("truemom2phpi0");
   b_truemomlab1phpi0 = fChain->GetBranch("truemomlab1phpi0");
   b_truemomlab2phpi0 = fChain->GetBranch("truemomlab2phpi0");
   b_trueth1phpi0 = fChain->GetBranch("trueth1phpi0");
   b_trueth2phpi0 = fChain->GetBranch("trueth2phpi0");
   b_Estar1daupi0 = fChain->GetBranch("Estar2dauPi0");
   b_Estar2daupi0 = fChain->GetBranch("Estar2dauPi0");
   b_truemrho = fChain->GetBranch("truemrho");
   b_nrecorho = fChain->GetBranch("nrecorho");
   b_mrho = fChain->GetBranch("mrho");
   b_mm2rho = fChain->GetBranch("mm2rho");
   b_truemompirho = fChain->GetBranch("truemompirho");
   b_truemompi0rho = fChain->GetBranch("truemompi0rho");
   b_mompirho = fChain->GetBranch("mompirho");
   b_mompi0rho = fChain->GetBranch("mompi0rho");
   b_truemrho0 = fChain->GetBranch("truemrho0");
   b_nrecorho0 = fChain->GetBranch("nrecorho0");
   b_mrho0 = fChain->GetBranch("mrho0");
   b_mm2rho0 = fChain->GetBranch("mm2rho0");
   b_truemom1pirho0 = fChain->GetBranch("truemom1pirho0");
   b_truemom2pirho0 = fChain->GetBranch("truemom2pirho0");
   b_mom1pirho0 = fChain->GetBranch("mom1pirho0");
   b_mom2pirho0 = fChain->GetBranch("mom2pirho0");
   b_truemomega = fChain->GetBranch("truemomega");
   b_nrecoomega = fChain->GetBranch("nrecoomega");
   b_momega = fChain->GetBranch("momega");
   b_mm2omega = fChain->GetBranch("mm2omega");
   b_truemom1piome = fChain->GetBranch("truemom1piome");
   b_truemom2piome = fChain->GetBranch("truemom2piome");
   b_truemompi0ome = fChain->GetBranch("truemompi0ome");
   b_truedalitzpi1pi2ome = fChain->GetBranch("truedalitzpi1pi2ome");
   b_truedalitzpi1pi0ome = fChain->GetBranch("truedalitzpi1pi0ome");
   b_truecosthome = fChain->GetBranch("truecosthome");
   b_mom1piome = fChain->GetBranch("mom1piome");
   b_mom2piome = fChain->GetBranch("mom2piome");
   b_mompi0ome = fChain->GetBranch("mompi0ome");
   b_dalitzpi1pi2ome = fChain->GetBranch("dalitzpi1pi2ome");
   b_dalitzpi1pi0ome = fChain->GetBranch("dalitzpi1pi0ome");
   b_costhome = fChain->GetBranch("costhome");
   b_nrecoeta = fChain->GetBranch("nrecoEta");
   b_indexbesteta = fChain->GetBranch("indexbestEta");
   b_barembesteta = fChain->GetBranch("barembestEta");
   b_mm2besteta = fChain->GetBranch("mm2bestEta");
   b_modeeta = fChain->GetBranch("modeEta");
   b_nrecoetap = fChain->GetBranch("nrecoEtap");
   b_indexbestetap = fChain->GetBranch("indexbestEtap");
   b_barembestetap = fChain->GetBranch("barembestEtap");
   b_mm2bestetap = fChain->GetBranch("mm2bestEtap");
   b_modeetap = fChain->GetBranch("modeEtap");
   b_etamassdauetap = fChain->GetBranch("EtamassdauEtap");
   b_rho0massdauetap = fChain->GetBranch("Rho0massdauEtap");
   b_gammamomdauetap = fChain->GetBranch("GammamomdauEtap");
   b_ma0 = fChain->GetBranch("ma0");
   b_mm2a0 = fChain->GetBranch("mm2a0");
   b_a0massetadau = fChain->GetBranch("a0massetadau");
   b_modea0 = fChain->GetBranch("modea0");
   b_ma0p = fChain->GetBranch("ma0p");
   b_mm2a0p = fChain->GetBranch("mm2a0p");
   b_a0pmassetadau = fChain->GetBranch("a0pmassetadau");
   b_modea0p = fChain->GetBranch("modea0p");
   return kTRUE;
}
// ----------------------------------------------------------------------
void theanal::readCuts(TString filename) {
  char  buffer[200];
  sprintf(buffer, "%s", filename.Data());
  ifstream is(buffer);
  char CutName[100];
  float CutValue;
  int ok(0);
  TOTALSTAT = 270. ;
  TOTALSTATMODEL = 270.;
  BRRATIOGENVALUE = 0.0017 ;
  BRRATIOVALUETAIL = 0.0017 ;
  PSTARFACT = 1.17;
  FITCATEGORY = 11;
  Q2CUT = 0;
  MXCUT = 1.6;
  USECB = 1;
  FIXMEANVALUE = 1;
  FIXSIGMA = 0;
  FIXARGUS1 = 0;
  FIXARGUS2 = 0;
  FIXCB1 = 0;
  FIXCB2 = 0;
  LEPTONPCUT = 1.;
  ELECTRONPCUT = 1.;
  MUONPCUT = 1.;
  DELTAM = 0.1;
  PRMM2CUT = 1000000.;
  MNUSQLOW = -1000.;
  MNUSQHIGH = 0.5;
  CHLOW = -.5;
  CHHIGH = .5;
  DEPL = 1.;
  SSBAR = 0;
  BTYPE = 2;
  LEPTTYPE = 2;
  FITTOTSHAPE = 0;
  MIXCORR = 2;
  FITMC = 1;
  FITOPT = 0;
  TOYMC = 0.;
  MULTIFIT = 0;
  BLINDING = 1;
  RANDOMSEED = 990717;
  BLINDSIZE = 0.51;
  ISSMEARALL = 0;              
  SMEARALLMEANVALUE = 0;
  SMEARALLSIGMA   = 1;      
  ISSMEARBKG      = 0;        
  SMEARBKGMEANVALUE = 0;
  SMEARBKGSIGMA    = 1; 
  DOBRECOWEIGHT = 0;
  DOBDECWEIGHT = 0;
  DODDECWEIGHT = 0;
  DOTRKWEIGHT = 0;
  DONEUWEIGHT = 0;
  MAXINTPUR = 1000.;
  MININTPUR = 0.;
  RUN = 0;
  CUTNNPI0 = -1000;
  DOFERMI = 0;
  DOTHEO = 0;
  FERMIAPP = 0;
  DELTAMB = 0;
  DELTAA = 0;


  MXCUTLOWEXCL=-1000.;
  MXCUTHIGHEXCL=1000.;
  MXCUTLOWEXCL1=-1000.;
  MXCUTHIGHEXCL1=1000.;
  MXCUTLOWEXCL2=-1000.;
  MXCUTHIGHEXCL2=1000.;
  MXCUTLOWEXCL3=-1000.;
  MXCUTHIGHEXCL3=1000.;
  MXCUTLOWEXCL4=-1000.;
  MXCUTHIGHEXCL4=1000.;
  NCHGLOWEXCL=0;
  NCHGHIGHEXCL=1000;
  MOM1MIN=-1000.;
  MOM2MIN=-1000.;
  NCOMBLOWEXCL=0;
  NCOMBHIGHEXCL=1000;
  NPI0LOWEXCL=0;
  NPI0HIGHEXCL=1000;
  DALITZCUT=1.;
  DAUETAMASS=1000.;
  DAUETAMASS2=1000.;
  DAUETAMASS3=1000.;
  DAUETAMASS4=1000.;
  DAURHOMASS=1000.;
  DAUGAMMAMOM = -1000.;
  DORIGHTNCHG = 0;

  MNUSQPI0LOW =  0.;
  MNUSQPI0HIGH = 0.;
  MNUSQETALOW =  0.;
  MNUSQETAHIGH = 0.;
  MNUSQETALOW =  0.;
  MNUSQETAHIGH = 0.;
  MNUSQRHOLOW =  0.;
  MNUSQRHOHIGH = 0.;
  MNUSQRHO0LOW =  0.;
  MNUSQRHO0HIGH = 0.;
  MNUSQOMEGALOW =  0.;
  MNUSQOMEGAHIGH = 0.;
  MAXENEU = 1000.;
  JPSIWIN = 1000.;

  while (is.getline(buffer, 200, '\n')) {
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);
    // -- 
    if (!strcmp(CutName, "totalStat")) { TOTALSTAT  = CutValue; ok = 1;}
    if (!strcmp(CutName, "totalStatModel")) { TOTALSTATMODEL = CutValue; ok = 1;}
    if (!strcmp(CutName, "BRRatioGenValue")) { BRRATIOGENVALUE  = CutValue; ok = 1;}
    if (!strcmp(CutName, "BRRatioValueTail")) { BRRATIOVALUETAIL  = CutValue; ok = 1;}
    if (!strcmp(CutName, "pstarfact")) { PSTARFACT  = CutValue; ok = 1;}
    if (!strcmp(CutName, "fitcategory")) { FITCATEGORY= int(CutValue); ok = 1;}
    if (!strcmp(CutName, "mxCut")) { MXCUT  = CutValue; ok = 1;}
    if (!strcmp(CutName, "q2Cut")) { Q2CUT = CutValue; ok = 1;}
    if (!strcmp(CutName, "useCB")) { USECB= int(CutValue); ok = 1;}
    if (!strcmp(CutName, "fixMeanValue")) { FIXMEANVALUE= int(CutValue); ok = 1;}
    if (!strcmp(CutName, "fixSigma")) { FIXSIGMA= int(CutValue); ok = 1;}
    if (!strcmp(CutName, "fixArgus1")) { FIXARGUS1= int(CutValue); ok = 1;}
    if (!strcmp(CutName, "fixArgus2")) { FIXARGUS2 = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "fixCB1")) { FIXCB1 = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "fixCB2")) { FIXCB2 = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "leptonPCut")) { LEPTONPCUT = CutValue; ok = 1;}
    if (!strcmp(CutName, "electronPCut")) { ELECTRONPCUT = CutValue; ok = 1;}
    if (!strcmp(CutName, "muonPCut")) { MUONPCUT = CutValue; ok = 1;}
    if (!strcmp(CutName, "deltam")) { DELTAM= CutValue; ok = 1;}
    if (!strcmp(CutName, "prmm2Cut")) { PRMM2CUT = CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqLow")) { MNUSQLOW= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqHigh")) { MNUSQHIGH= CutValue; ok = 1;}
    if (!strcmp(CutName, "chLow")) { CHLOW= CutValue; ok = 1;}
    if (!strcmp(CutName, "chHigh")) { CHHIGH= CutValue; ok = 1;}
    if (!strcmp(CutName, "depl")) { DEPL = CutValue; ok = 1;}
    if (!strcmp(CutName, "ssbar")) { SSBAR = CutValue; ok = 1;}
    if (!strcmp(CutName, "Btype")) { BTYPE = CutValue; ok = 1;}
    if (!strcmp(CutName, "lepttype")) { LEPTTYPE = CutValue; ok = 1;}
    if (!strcmp(CutName, "fittotshape")) { FITTOTSHAPE= int(CutValue); ok = 1;}
    if (!strcmp(CutName, "mixcorr")) { MIXCORR = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "fitMC")) { FITMC= int(CutValue); ok = 1;}
    if (!strcmp(CutName, "fitOption")) { FITOPT= int(CutValue); ok = 1;}
    if (!strcmp(CutName, "toyMC")) { TOYMC= CutValue; ok = 1;}
    if (!strcmp(CutName, "multifit")) { MULTIFIT = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "blinding")) { BLINDING = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "randomseed")) { RANDOMSEED = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "blindsize")) { BLINDSIZE = CutValue; ok = 1; }
    if (!strcmp(CutName, "issmearAll")) {ISSMEARALL = int(CutValue); ok = 1;}       
    if (!strcmp(CutName, "smearAllMeanValue")) {SMEARALLMEANVALUE = CutValue; ok = 1;} 
    if (!strcmp(CutName, "smearAllSigma")) {SMEARALLSIGMA = CutValue; ok = 1;}     
    if (!strcmp(CutName, "issmearBkg")) {ISSMEARBKG = int(CutValue); ok = 1;}       
    if (!strcmp(CutName, "smearBkgMeanValue")) {SMEARBKGMEANVALUE = CutValue; ok = 1;} 
    if (!strcmp(CutName, "smearBkgSigma")) {SMEARBKGSIGMA = CutValue; ok = 1;}    
    if (!strcmp(CutName, "dobrecoreweight")) { DOBRECOWEIGHT = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "doBdecreweight")) { DOBDECWEIGHT = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "doDdecreweight")) { DODDECWEIGHT = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "dotrkreweight")) { DOTRKWEIGHT = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "doneureweight")) { DONEUWEIGHT = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "maxintpur")) {MAXINTPUR = CutValue; ok = 1;}
    if (!strcmp(CutName, "minintpur")) {MININTPUR = CutValue; ok = 1;}
    if (!strcmp(CutName, "run")) {RUN = CutValue; ok = 1;}
    if (!strcmp(CutName, "nnpi0")) {CUTNNPI0 = int(CutValue); ok = 1;}
    if (!strcmp(CutName, "fermiapp")) {FERMIAPP = CutValue; ok = 1;}
    if (!strcmp(CutName, "deltamb")) {DELTAMB = CutValue; ok = 1;}
    if (!strcmp(CutName, "deltaa")) {DELTAA = CutValue; ok = 1;}

    if (!strcmp(CutName, "mxCutLowExcl")) { MXCUTLOWEXCL  = CutValue; ok = 1;}
    if (!strcmp(CutName, "mxCutHighExcl")) { MXCUTHIGHEXCL  = CutValue; ok = 1;}
    if (!strcmp(CutName, "mxCutLowExcl1")) { MXCUTLOWEXCL1  = CutValue; ok = 1;}
    if (!strcmp(CutName, "mxCutHighExcl1")) { MXCUTHIGHEXCL1  = CutValue; ok = 1;}
    if (!strcmp(CutName, "mxCutLowExcl2")) { MXCUTLOWEXCL2  = CutValue; ok = 1;}
    if (!strcmp(CutName, "mxCutHighExcl2")) { MXCUTHIGHEXCL2  = CutValue; ok = 1;}
    if (!strcmp(CutName, "mxCutLowExcl3")) { MXCUTLOWEXCL3  = CutValue; ok = 1;}
    if (!strcmp(CutName, "mxCutHighExcl3")) { MXCUTHIGHEXCL3  = CutValue; ok = 1;}
    if (!strcmp(CutName, "mxCutLowExcl4")) { MXCUTLOWEXCL4  = CutValue; ok = 1;}
    if (!strcmp(CutName, "mxCutHighExcl4")) { MXCUTHIGHEXCL4  = CutValue; ok = 1;}
    if (!strcmp(CutName, "nchgLowExcl")) { NCHGLOWEXCL = CutValue; ok = 1;}
    if (!strcmp(CutName, "nchgHighExcl")) { NCHGHIGHEXCL = CutValue; ok = 1;}
    if (!strcmp(CutName, "ncombLowExcl")) { NCOMBLOWEXCL = CutValue; ok = 1;}
    if (!strcmp(CutName, "ncombHighExcl")) { NCOMBHIGHEXCL = CutValue; ok = 1;}
    if (!strcmp(CutName, "npi0LowExcl")) { NPI0LOWEXCL = CutValue; ok = 1;}
    if (!strcmp(CutName, "npi0HighExcl")) { NPI0HIGHEXCL = CutValue; ok = 1;}
    if (!strcmp(CutName, "mom1min")) { MOM1MIN = CutValue; ok = 1;}
    if (!strcmp(CutName, "mom2min")) { MOM2MIN = CutValue; ok = 1;}
    if (!strcmp(CutName, "dalitzcut")) { DALITZCUT = CutValue; ok = 1;}
    if (!strcmp(CutName, "dauetamass")) { DAUETAMASS = CutValue; ok = 1;}    
    if (!strcmp(CutName, "dauetamass2")) { DAUETAMASS2 = CutValue; ok = 1;}
    if (!strcmp(CutName, "dauetamass3")) { DAUETAMASS3 = CutValue; ok = 1;}
    if (!strcmp(CutName, "dauetamass4")) { DAUETAMASS4 = CutValue; ok = 1;}
    if (!strcmp(CutName, "daurhomass")) { DAURHOMASS = CutValue; ok = 1;}    
    if (!strcmp(CutName, "daugammamom")) { DAUGAMMAMOM = CutValue; ok = 1;}    
    if (!strcmp(CutName, "dorightnchg")) { DORIGHTNCHG = CutValue; ok = 1;}    

    if (!strcmp(CutName, "mnuSqpi0Low")) { MNUSQPI0LOW= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqpi0High")) { MNUSQPI0HIGH= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqetaLow")) { MNUSQETALOW= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqetaHigh")) { MNUSQETAHIGH= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqetaLow")) { MNUSQETALOW= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqetaHigh")) { MNUSQETAHIGH= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqrhoLow")) { MNUSQRHOLOW= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqrhoHigh")) { MNUSQRHOHIGH= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqrho0Low")) { MNUSQRHO0LOW= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqrho0High")) { MNUSQRHO0HIGH= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqomegaLow")) { MNUSQOMEGALOW= CutValue; ok = 1;}
    if (!strcmp(CutName, "mnuSqomegaHigh")) { MNUSQOMEGAHIGH= CutValue; ok = 1;}    
    if (!strcmp(CutName, "maxeneu")) {MAXENEU  = CutValue; ok = 1;}
    if (!strcmp(CutName, "jpsiwindow")) {JPSIWIN  = CutValue; ok = 1;}

    if (ok == 0)  cout << "==> exclfitNtp::readCuts() Error: Don't know about variable " << CutName << endl;
  }
  dumpCuts();

}

// ----------------------------------------------------------------------
void theanal::dumpCuts() {
  cout << "====================================" << endl;
  cout << " fitcategory         :"  << FITCATEGORY   <<  endl;  
  cout << " q2Cut               :"  << Q2CUT   <<  endl;  
  cout << " leptonPCut          :"  << LEPTONPCUT   <<  endl;  
  cout << " electronPCut        :"  << ELECTRONPCUT   <<  endl;
  cout << " muonPCut            :"  << MUONPCUT   <<  endl;
  cout << " deltam              :"  << DELTAM  <<  endl;  
  cout << " prmm2Cut            :"  << PRMM2CUT   <<  endl;  
  cout << " mnuSqLow            :"  << MNUSQLOW   <<  endl;  
  cout << " mnuSqHigh           :"  << MNUSQHIGH   <<  endl;  
  cout << " chLow               :"  << CHLOW   <<  endl;  
  cout << " chHigh              :"  << CHHIGH   <<  endl;  
  cout << " Btype               :"  << BTYPE   <<  endl;  
  cout << " lepttype            :"  << LEPTTYPE   <<  endl;  
  cout << " maxintpur           :"  << MAXINTPUR    <<  endl;
  cout << " minintpur           :"  << MININTPUR    <<  endl;
  cout << " run                 :"  << RUN    <<  endl;
  cout << " mxCutLowExcl        :"  << MXCUTLOWEXCL   <<  endl;  
  cout << " mxCutHighExcl       :"  << MXCUTHIGHEXCL   <<  endl;  
  cout << " mxCutLowExcl1        :"  << MXCUTLOWEXCL1   <<  endl;
  cout << " mxCutHighExcl1       :"  << MXCUTHIGHEXCL1   <<  endl;
  cout << " mxCutLowExcl2        :"  << MXCUTLOWEXCL2   <<  endl;
  cout << " mxCutHighExcl2       :"  << MXCUTHIGHEXCL2   <<  endl;
  cout << " mxCutLowExcl3        :"  << MXCUTLOWEXCL3   <<  endl;
  cout << " mxCutHighExcl3       :"  << MXCUTHIGHEXCL3   <<  endl;
  cout << " mxCutLowExcl4        :"  << MXCUTLOWEXCL4   <<  endl;
  cout << " mxCutHighExcl4       :"  << MXCUTHIGHEXCL4   <<  endl;
  cout << " nchgLowExcl         :"  << NCHGLOWEXCL   <<  endl;  
  cout << " nchgHighExcl        :"  << NCHGHIGHEXCL   <<  endl;  
  cout << " ncombLowExcl        :"  << NCOMBLOWEXCL   <<  endl;  
  cout << " ncombHighExcl       :"  << NCOMBHIGHEXCL   <<  endl;  
  cout << " npi0LowExcl         :"  << NPI0LOWEXCL   <<  endl;  
  cout << " npi0HighExcl        :"  << NPI0HIGHEXCL   <<  endl;  
  cout << " mom1min             :"  << MOM1MIN   <<  endl;  
  cout << " mom2min             :"  << MOM2MIN   <<  endl;  
  cout << " dalitzcut           :"  << DALITZCUT   <<  endl;  
  cout << " dauetamass          :"  << DAUETAMASS   <<  endl;  
  cout << " dauetamass2          :"  << DAUETAMASS2   <<  endl;
  cout << " dauetamass3          :"  << DAUETAMASS3   <<  endl;
  cout << " dauetamass4          :"  << DAUETAMASS4   <<  endl;
  cout << " daurhomass          :"  << DAURHOMASS   <<  endl;  
  cout << " daugammamom         :"  << DAUGAMMAMOM   <<  endl;  
  cout << " dorightnchg         :"  << DORIGHTNCHG   <<  endl;  
  cout << " mnuSqpi0Low         :"  << MNUSQPI0LOW   <<  endl;  
  cout << " mnuSqpi0High        :"  << MNUSQPI0HIGH   <<  endl;  
  cout << " mnuSqetaLow         :"  << MNUSQETALOW   <<  endl;  
  cout << " mnuSqetaHigh        :"  << MNUSQETAHIGH   <<  endl;  
  cout << " mnuSqetaLow         :"  << MNUSQETALOW   <<  endl;  
  cout << " mnuSqetaHigh        :"  << MNUSQETAHIGH   <<  endl;  
  cout << " mnuSqrhoLow         :"  << MNUSQRHOLOW   <<  endl;  
  cout << " mnuSqrhoHigh        :"  << MNUSQRHOHIGH   <<  endl;  
  cout << " mnuSqrho0Low        :"  << MNUSQRHO0LOW   <<  endl;  
  cout << " mnuSqrho0High       :"  << MNUSQRHO0HIGH   <<  endl;  
  cout << " mnuSqomegaLow       :"  << MNUSQOMEGALOW   <<  endl;  
  cout << " mnuSqomegaHigh      :"  << MNUSQOMEGAHIGH   <<  endl;  
  cout << " maxeneu             :"  << MAXENEU   <<  endl;  
  cout << " jpsiwindow          :"  << JPSIWIN   <<  endl;  
  cout << "====================================" << endl;
  cout << "" << endl;
  cout << "" << endl;
}


theanal::~theanal()
{
  fHistFile->cd();
  fHistFile->Write();
  fHistFile->Close();
  delete fHistFile;
  delete Dvar;
  delete Bsem;
   if (!fChain) return;
   //   delete fChain->GetCurrentFile();
}

Int_t theanal::GetEntry(Int_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Int_t theanal::LoadTree(Int_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Int_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->IsA() != TChain::Class()) return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void theanal::Show(Int_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t theanal::Cut(Int_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
