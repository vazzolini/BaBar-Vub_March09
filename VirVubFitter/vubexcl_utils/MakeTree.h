//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 13 05:49:56 2006 by ROOT version 4.01/02
// from TChain ntp1/
//////////////////////////////////////////////////////////

#ifndef MakeTree_h
#define MakeTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class MakeTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Int_t           run;
   Int_t           upper;
   Int_t           lower;
   Int_t           nchg;
   Int_t           nneu;
   Int_t           nB;
   Int_t           IdB;
   Int_t           brecoflav;
   Int_t           brecocharge;
   Int_t           xcharge;
   Int_t           mode;
   Double_t        mes;
   Double_t        de;
   Double_t        pB;
   Double_t        eB;
   Double_t        thetaB;
   Double_t        phiB;
   Double_t        pur;
   Double_t        intpur;
   Int_t           nle;
   Int_t           nel;
   Int_t           nmu;
   Int_t           nkp;
   Int_t           nks;
   Int_t           nlept500;
   Int_t           nelec500;
   Int_t           nmu500;
   Int_t           brems;
   Double_t        plab;
   Double_t        tlab;
   Double_t        flab;
   Double_t        elab;
   Double_t        pcms;
   Double_t        ecms;
   Double_t        tcms;
   Double_t        fcms;
   Int_t           lcharge;
   Int_t           isele;
   Double_t        wdeltam;
   Double_t        Eneualt;
   Int_t           indexbestPi;
   Int_t           chbestPi;
   Double_t        mbestPi;
   Double_t        barembestPi;
   Double_t        mfitbestPi;
   Double_t        baremfitbestPi;
   Double_t        mm2bestPi;
   Double_t        q2bestPi;
   Double_t        q2fitbestPi;
   Double_t        mm2fitbestPi;
   Int_t           nrecoPi;
   Float_t         mm2Pi[100];   //[nrecoPi]
   Float_t         q2Pi[100];   //[nrecoPi]
   Float_t         pPi[100];   //[nrecoPi]
   Float_t         thPi[100];   //[nrecoPi]
   Float_t         phiPi[100];   //[nrecoPi]
   Float_t         nominalmPi[100];   //[nrecoPi]
   Float_t         baremPi[100];   //[nrecoPi]
   Float_t         ElabPi[100];   //[nrecoPi]
   Int_t           chPi[100];   //[nrecoPi]
   Float_t         mm2fitPi[100];   //[nrecoPi]
   Float_t         q2fitPi[100];   //[nrecoPi]
   Float_t         massfitPi[100];   //[nrecoPi]
   Float_t         pfitPi[100];   //[nrecoPi]
   Float_t         thetafitPi[100];   //[nrecoPi]
   Float_t         phifitPi[100];   //[nrecoPi]
   Float_t         nominalmfitPi[100];   //[nrecoPi]
   Float_t         baremfitPi[100];   //[nrecoPi]
   Float_t         ElabfitPi[100];   //[nrecoPi]
   Float_t         chisqPi[100];   //[nrecoPi]
   Float_t         globchisqPi[100];   //[nrecoPi]
   Float_t         probchisqPi[100];   //[nrecoPi]
   Int_t           ndofPi[100];   //[nrecoPi]
   Float_t         pstarfitleptPi[100];   //[nrecoPi]
   Float_t         pfitleptPi[100];   //[nrecoPi]
   Float_t         thetafitleptPi[100];   //[nrecoPi]
   Float_t         phifitleptPi[100];   //[nrecoPi]
   Float_t         pfitBPi[100];   //[nrecoPi]
   Float_t         massBfitPi[100];   //[nrecoPi]
   Float_t         thetafitBPi[100];   //[nrecoPi]
   Float_t         phifitBPi[100];   //[nrecoPi]
   Int_t           lepmapPi[100];   //[nrecoPi]
   Float_t         massjpsiPi[100];   //[nrecoPi]
   Float_t         baremassjpsiPi[100];   //[nrecoPi]
   Int_t           indexbestPi0;
   Int_t           chbestPi0;
   Double_t        mbestPi0;
   Double_t        barembestPi0;
   Double_t        mfitbestPi0;
   Double_t        baremfitbestPi0;
   Double_t        mm2bestPi0;
   Double_t        q2bestPi0;
   Double_t        q2fitbestPi0;
   Double_t        mm2fitbestPi0;
   Int_t           nrecoPi0;
   Float_t         mm2Pi0[100];   //[nrecoPi0]
   Float_t         q2Pi0[100];   //[nrecoPi0]
   Float_t         pPi0[100];   //[nrecoPi0]
   Float_t         thPi0[100];   //[nrecoPi0]
   Float_t         phiPi0[100];   //[nrecoPi0]
   Float_t         nominalmPi0[100];   //[nrecoPi0]
   Float_t         baremPi0[100];   //[nrecoPi0]
   Float_t         ElabPi0[100];   //[nrecoPi0]
   Int_t           chPi0[100];   //[nrecoPi0]
   Float_t         mm2fitPi0[100];   //[nrecoPi0]
   Float_t         q2fitPi0[100];   //[nrecoPi0]
   Float_t         massfitPi0[100];   //[nrecoPi0]
   Float_t         pfitPi0[100];   //[nrecoPi0]
   Float_t         thetafitPi0[100];   //[nrecoPi0]
   Float_t         phifitPi0[100];   //[nrecoPi0]
   Float_t         nominalmfitPi0[100];   //[nrecoPi0]
   Float_t         baremfitPi0[100];   //[nrecoPi0]
   Float_t         ElabfitPi0[100];   //[nrecoPi0]
   Float_t         chisqPi0[100];   //[nrecoPi0]
   Float_t         globchisqPi0[100];   //[nrecoPi0]
   Float_t         probchisqPi0[100];   //[nrecoPi0]
   Int_t           ndofPi0[100];   //[nrecoPi0]
   Float_t         pstarfitleptPi0[100];   //[nrecoPi0]
   Float_t         pfitleptPi0[100];   //[nrecoPi0]
   Float_t         thetafitleptPi0[100];   //[nrecoPi0]
   Float_t         phifitleptPi0[100];   //[nrecoPi0]
   Float_t         pfitBPi0[100];   //[nrecoPi0]
   Float_t         massBfitPi0[100];   //[nrecoPi0]
   Float_t         thetafitBPi0[100];   //[nrecoPi0]
   Float_t         phifitBPi0[100];   //[nrecoPi0]
   Int_t           ndauPi0[100];   //[nrecoPi0]
   Float_t         Estar1dauPi0[100];   //[nrecoPi0]
   Float_t         Estar2dauPi0[100];   //[nrecoPi0]
   Float_t         Elab1dauPi0[100];   //[nrecoPi0]
   Float_t         Elab2dauPi0[100];   //[nrecoPi0]
   Float_t         Estar1fitdauPi0[100];   //[nrecoPi0]
   Float_t         Estar2fitdauPi0[100];   //[nrecoPi0]
   Float_t         Elab1fitdauPi0[100];   //[nrecoPi0]
   Float_t         Elab2fitdauPi0[100];   //[nrecoPi0]
   Int_t           indexbestEta;
   Int_t           chbestEta;
   Double_t        mbestEta;
   Double_t        barembestEta;
   Double_t        mfitbestEta;
   Double_t        baremfitbestEta;
   Double_t        mm2bestEta;
   Double_t        q2bestEta;
   Double_t        q2fitbestEta;
   Double_t        mm2fitbestEta;
   Int_t           nrecoEta;
   Float_t         mm2Eta[100];   //[nrecoEta]
   Float_t         q2Eta[100];   //[nrecoEta]
   Float_t         pEta[100];   //[nrecoEta]
   Float_t         thEta[100];   //[nrecoEta]
   Float_t         phiEta[100];   //[nrecoEta]
   Float_t         nominalmEta[100];   //[nrecoEta]
   Float_t         baremEta[100];   //[nrecoEta]
   Float_t         ElabEta[100];   //[nrecoEta]
   Int_t           chEta[100];   //[nrecoEta]
   Float_t         mm2fitEta[100];   //[nrecoEta]
   Float_t         q2fitEta[100];   //[nrecoEta]
   Float_t         massfitEta[100];   //[nrecoEta]
   Float_t         pfitEta[100];   //[nrecoEta]
   Float_t         thetafitEta[100];   //[nrecoEta]
   Float_t         phifitEta[100];   //[nrecoEta]
   Float_t         nominalmfitEta[100];   //[nrecoEta]
   Float_t         baremfitEta[100];   //[nrecoEta]
   Float_t         ElabfitEta[100];   //[nrecoEta]
   Float_t         chisqEta[100];   //[nrecoEta]
   Float_t         globchisqEta[100];   //[nrecoEta]
   Float_t         probchisqEta[100];   //[nrecoEta]
   Int_t           ndofEta[100];   //[nrecoEta]
   Float_t         pstarfitleptEta[100];   //[nrecoEta]
   Float_t         pfitleptEta[100];   //[nrecoEta]
   Float_t         thetafitleptEta[100];   //[nrecoEta]
   Float_t         phifitleptEta[100];   //[nrecoEta]
   Float_t         pfitBEta[100];   //[nrecoEta]
   Float_t         massBfitEta[100];   //[nrecoEta]
   Float_t         thetafitBEta[100];   //[nrecoEta]
   Float_t         phifitBEta[100];   //[nrecoEta]
   Int_t           ndauEta[100];   //[nrecoEta]
   Int_t           modeEta[100];   //[nrecoEta]
   Float_t         Estar1dauEta[100];   //[nrecoEta]
   Float_t         Estar2dauEta[100];   //[nrecoEta]
   Float_t         Estar3dauEta[100];   //[nrecoEta]
   Float_t         Elab1dauEta[100];   //[nrecoEta]
   Float_t         Elab2dauEta[100];   //[nrecoEta]
   Float_t         Elab3dauEta[100];   //[nrecoEta]
   Float_t         Estar1fitdauEta[100];   //[nrecoEta]
   Float_t         Estar2fitdauEta[100];   //[nrecoEta]
   Float_t         Estar3fitdauEta[100];   //[nrecoEta]
   Float_t         Elab1fitdauEta[100];   //[nrecoEta]
   Float_t         Elab2fitdauEta[100];   //[nrecoEta]
   Float_t         Elab3fitdauEta[100];   //[nrecoEta]
   Int_t           indexbestEtap;
   Int_t           chbestEtap;
   Double_t        mbestEtap;
   Double_t        barembestEtap;
   Double_t        mfitbestEtap;
   Double_t        baremfitbestEtap;
   Double_t        mm2bestEtap;
   Double_t        q2bestEtap;
   Double_t        q2fitbestEtap;
   Double_t        mm2fitbestEtap;
   Int_t           nrecoEtap;
   Float_t         mm2Etap[100];   //[nrecoEtap]
   Float_t         q2Etap[100];   //[nrecoEtap]
   Float_t         pEtap[100];   //[nrecoEtap]
   Float_t         thEtap[100];   //[nrecoEtap]
   Float_t         phiEtap[100];   //[nrecoEtap]
   Float_t         nominalmEtap[100];   //[nrecoEtap]
   Float_t         baremEtap[100];   //[nrecoEtap]
   Float_t         ElabEtap[100];   //[nrecoEtap]
   Int_t           chEtap[100];   //[nrecoEtap]
   Float_t         mm2fitEtap[100];   //[nrecoEtap]
   Float_t         q2fitEtap[100];   //[nrecoEtap]
   Float_t         massfitEtap[100];   //[nrecoEtap]
   Float_t         pfitEtap[100];   //[nrecoEtap]
   Float_t         thetafitEtap[100];   //[nrecoEtap]
   Float_t         phifitEtap[100];   //[nrecoEtap]
   Float_t         nominalmfitEtap[100];   //[nrecoEtap]
   Float_t         baremfitEtap[100];   //[nrecoEtap]
   Float_t         ElabfitEtap[100];   //[nrecoEtap]
   Float_t         chisqEtap[100];   //[nrecoEtap]
   Float_t         globchisqEtap[100];   //[nrecoEtap]
   Float_t         probchisqEtap[100];   //[nrecoEtap]
   Int_t           ndofEtap[100];   //[nrecoEtap]
   Float_t         pstarfitleptEtap[100];   //[nrecoEtap]
   Float_t         pfitleptEtap[100];   //[nrecoEtap]
   Float_t         thetafitleptEtap[100];   //[nrecoEtap]
   Float_t         phifitleptEtap[100];   //[nrecoEtap]
   Float_t         pfitBEtap[100];   //[nrecoEtap]
   Float_t         massBfitEtap[100];   //[nrecoEtap]
   Float_t         thetafitBEtap[100];   //[nrecoEtap]
   Float_t         phifitBEtap[100];   //[nrecoEtap]
   Int_t           ndauEtap[100];   //[nrecoEtap]
   Int_t           modeEtap[100];   //[nrecoEtap]
   Float_t         EtamassdauEtap[100];   //[nrecoEtap]
   Float_t         Rho0massdauEtap[100];   //[nrecoEtap]
   Float_t         GammamomdauEtap[100];   //[nrecoEtap]
   Float_t         Estar1dauEtap[100];   //[nrecoEtap]
   Float_t         Estar2dauEtap[100];   //[nrecoEtap]
   Float_t         Estar3dauEtap[100];   //[nrecoEtap]
   Float_t         Elab1dauEtap[100];   //[nrecoEtap]
   Float_t         Elab2dauEtap[100];   //[nrecoEtap]
   Float_t         Elab3dauEtap[100];   //[nrecoEtap]
   Float_t         Estar1fitdauEtap[100];   //[nrecoEtap]
   Float_t         Estar2fitdauEtap[100];   //[nrecoEtap]
   Float_t         Estar3fitdauEtap[100];   //[nrecoEtap]
   Float_t         Elab1fitdauEtap[100];   //[nrecoEtap]
   Float_t         Elab2fitdauEtap[100];   //[nrecoEtap]
   Float_t         Elab3fitdauEtap[100];   //[nrecoEtap]
   Int_t           isassocB;
   Int_t           isassocB_GHIT;
   Double_t        ass_deltapB;
   Int_t           ch1B;
   Int_t           ch2B;
   Int_t           chunm;
   Int_t           neu1B;
   Int_t           neu2B;
   Int_t           neuunm;
   Int_t           brecoqual;
   Double_t        brecoqualangle;
   Int_t           chgdaugen;
   Int_t           neudaugen;
   Int_t           vub;
   Int_t           vcb;
   Int_t           other;
   Int_t           nvubexcl;
   Int_t           nvubnres;
   Double_t        mxhadgen;
   Double_t        mxhadgenwoph;
   Double_t        pcmsgen;
   Double_t        ecmsgen;
   Double_t        tcmsgen;
   Double_t        fcmsgen;
   Double_t        pxhadgen;
   Double_t        exhadgen;
   Double_t        fxhadgen;
   Double_t        txhadgen;
   Double_t        q2Gen;
   Double_t        ctvgen;
   Double_t        ctlgen;
   Double_t        chigen;
   Int_t           Gvxbtyp;
   Int_t           GSem;
   Int_t           GfDpi;
   Int_t           GfDpiz;
   Int_t           GfDk;
   Int_t           GfDks;
   Int_t           GfDkl;
   Int_t           GfDlep;
   Int_t           GfDgam;
   Int_t           GfDnu;
   Int_t           GfD0Ds;
   Int_t           GfDDs;
   Int_t           GfDkspiopio;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_upper;   //!
   TBranch        *b_lower;   //!
   TBranch        *b_nchg;   //!
   TBranch        *b_nneu;   //!
   TBranch        *b_nB;   //!
   TBranch        *b_IdB;   //!
   TBranch        *b_brecoflav;   //!
   TBranch        *b_brecocharge;   //!
   TBranch        *b_xcharge;   //!
   TBranch        *b_mode;   //!
   TBranch        *b_mes;   //!
   TBranch        *b_de;   //!
   TBranch        *b_pB;   //!
   TBranch        *b_eB;   //!
   TBranch        *b_thetaB;   //!
   TBranch        *b_phiB;   //!
   TBranch        *b_pur;   //!
   TBranch        *b_intpur;   //!
   TBranch        *b_nle;   //!
   TBranch        *b_nel;   //!
   TBranch        *b_nmu;   //!
   TBranch        *b_nkp;   //!
   TBranch        *b_nks;   //!
   TBranch        *b_nlept500;   //!
   TBranch        *b_nelec500;   //!
   TBranch        *b_nmu500;   //!
   TBranch        *b_brems;   //!
   TBranch        *b_plab;   //!
   TBranch        *b_tlab;   //!
   TBranch        *b_flab;   //!
   TBranch        *b_elab;   //!
   TBranch        *b_pcms;   //!
   TBranch        *b_ecms;   //!
   TBranch        *b_tcms;   //!
   TBranch        *b_fcms;   //!
   TBranch        *b_lcharge;   //!
   TBranch        *b_isele;   //!
   TBranch        *b_wdeltam;   //!
   TBranch        *b_Eneualt;   //!
   TBranch        *b_indexbestPi;   //!
   TBranch        *b_chbestPi;   //!
   TBranch        *b_mbestPi;   //!
   TBranch        *b_barembestPi;   //!
   TBranch        *b_mfitbestPi;   //!
   TBranch        *b_baremfitbestPi;   //!
   TBranch        *b_mm2bestPi;   //!
   TBranch        *b_q2bestPi;   //!
   TBranch        *b_q2fitbestPi;   //!
   TBranch        *b_mm2fitbestPi;   //!
   TBranch        *b_nrecoPi;   //!
   TBranch        *b_mm2Pi;   //!
   TBranch        *b_q2Pi;   //!
   TBranch        *b_pPi;   //!
   TBranch        *b_thPi;   //!
   TBranch        *b_phiPi;   //!
   TBranch        *b_nominalmPi;   //!
   TBranch        *b_baremPi;   //!
   TBranch        *b_ElabPi;   //!
   TBranch        *b_chPi;   //!
   TBranch        *b_mm2fitPi;   //!
   TBranch        *b_q2fitPi;   //!
   TBranch        *b_massfitPi;   //!
   TBranch        *b_pfitPi;   //!
   TBranch        *b_thetafitPi;   //!
   TBranch        *b_phifitPi;   //!
   TBranch        *b_nominalmfitPi;   //!
   TBranch        *b_baremfitPi;   //!
   TBranch        *b_ElabfitPi;   //!
   TBranch        *b_chisqPi;   //!
   TBranch        *b_globchisqPi;   //!
   TBranch        *b_probchisqPi;   //!
   TBranch        *b_ndofPi;   //!
   TBranch        *b_pstarfitleptPi;   //!
   TBranch        *b_pfitleptPi;   //!
   TBranch        *b_thetafitleptPi;   //!
   TBranch        *b_phifitleptPi;   //!
   TBranch        *b_pfitBPi;   //!
   TBranch        *b_massBfitPi;   //!
   TBranch        *b_thetafitBPi;   //!
   TBranch        *b_phifitBPi;   //!
   TBranch        *b_lepmapPi;   //!
   TBranch        *b_massjpsiPi;   //!
   TBranch        *b_baremassjpsiPi;   //!
   TBranch        *b_indexbestPi0;   //!
   TBranch        *b_chbestPi0;   //!
   TBranch        *b_mbestPi0;   //!
   TBranch        *b_barembestPi0;   //!
   TBranch        *b_mfitbestPi0;   //!
   TBranch        *b_baremfitbestPi0;   //!
   TBranch        *b_mm2bestPi0;   //!
   TBranch        *b_q2bestPi0;   //!
   TBranch        *b_q2fitbestPi0;   //!
   TBranch        *b_mm2fitbestPi0;   //!
   TBranch        *b_nrecoPi0;   //!
   TBranch        *b_mm2Pi0;   //!
   TBranch        *b_q2Pi0;   //!
   TBranch        *b_pPi0;   //!
   TBranch        *b_thPi0;   //!
   TBranch        *b_phiPi0;   //!
   TBranch        *b_nominalmPi0;   //!
   TBranch        *b_baremPi0;   //!
   TBranch        *b_ElabPi0;   //!
   TBranch        *b_chPi0;   //!
   TBranch        *b_mm2fitPi0;   //!
   TBranch        *b_q2fitPi0;   //!
   TBranch        *b_massfitPi0;   //!
   TBranch        *b_pfitPi0;   //!
   TBranch        *b_thetafitPi0;   //!
   TBranch        *b_phifitPi0;   //!
   TBranch        *b_nominalmfitPi0;   //!
   TBranch        *b_baremfitPi0;   //!
   TBranch        *b_ElabfitPi0;   //!
   TBranch        *b_chisqPi0;   //!
   TBranch        *b_globchisqPi0;   //!
   TBranch        *b_probchisqPi0;   //!
   TBranch        *b_ndofPi0;   //!
   TBranch        *b_pstarfitleptPi0;   //!
   TBranch        *b_pfitleptPi0;   //!
   TBranch        *b_thetafitleptPi0;   //!
   TBranch        *b_phifitleptPi0;   //!
   TBranch        *b_pfitBPi0;   //!
   TBranch        *b_massBfitPi0;   //!
   TBranch        *b_thetafitBPi0;   //!
   TBranch        *b_phifitBPi0;   //!
   TBranch        *b_ndauPi0;   //!
   TBranch        *b_Estar1dauPi0;   //!
   TBranch        *b_Estar2dauPi0;   //!
   TBranch        *b_Elab1dauPi0;   //!
   TBranch        *b_Elab2dauPi0;   //!
   TBranch        *b_Estar1fitdauPi0;   //!
   TBranch        *b_Estar2fitdauPi0;   //!
   TBranch        *b_Elab1fitdauPi0;   //!
   TBranch        *b_Elab2fitdauPi0;   //!
   TBranch        *b_indexbestEta;   //!
   TBranch        *b_chbestEta;   //!
   TBranch        *b_mbestEta;   //!
   TBranch        *b_barembestEta;   //!
   TBranch        *b_mfitbestEta;   //!
   TBranch        *b_baremfitbestEta;   //!
   TBranch        *b_mm2bestEta;   //!
   TBranch        *b_q2bestEta;   //!
   TBranch        *b_q2fitbestEta;   //!
   TBranch        *b_mm2fitbestEta;   //!
   TBranch        *b_nrecoEta;   //!
   TBranch        *b_mm2Eta;   //!
   TBranch        *b_q2Eta;   //!
   TBranch        *b_pEta;   //!
   TBranch        *b_thEta;   //!
   TBranch        *b_phiEta;   //!
   TBranch        *b_nominalmEta;   //!
   TBranch        *b_baremEta;   //!
   TBranch        *b_ElabEta;   //!
   TBranch        *b_chEta;   //!
   TBranch        *b_mm2fitEta;   //!
   TBranch        *b_q2fitEta;   //!
   TBranch        *b_massfitEta;   //!
   TBranch        *b_pfitEta;   //!
   TBranch        *b_thetafitEta;   //!
   TBranch        *b_phifitEta;   //!
   TBranch        *b_nominalmfitEta;   //!
   TBranch        *b_baremfitEta;   //!
   TBranch        *b_ElabfitEta;   //!
   TBranch        *b_chisqEta;   //!
   TBranch        *b_globchisqEta;   //!
   TBranch        *b_probchisqEta;   //!
   TBranch        *b_ndofEta;   //!
   TBranch        *b_pstarfitleptEta;   //!
   TBranch        *b_pfitleptEta;   //!
   TBranch        *b_thetafitleptEta;   //!
   TBranch        *b_phifitleptEta;   //!
   TBranch        *b_pfitBEta;   //!
   TBranch        *b_massBfitEta;   //!
   TBranch        *b_thetafitBEta;   //!
   TBranch        *b_phifitBEta;   //!
   TBranch        *b_ndauEta;   //!
   TBranch        *b_modeEta;   //!
   TBranch        *b_Estar1dauEta;   //!
   TBranch        *b_Estar2dauEta;   //!
   TBranch        *b_Estar3dauEta;   //!
   TBranch        *b_Elab1dauEta;   //!
   TBranch        *b_Elab2dauEta;   //!
   TBranch        *b_Elab3dauEta;   //!
   TBranch        *b_Estar1fitdauEta;   //!
   TBranch        *b_Estar2fitdauEta;   //!
   TBranch        *b_Estar3fitdauEta;   //!
   TBranch        *b_Elab1fitdauEta;   //!
   TBranch        *b_Elab2fitdauEta;   //!
   TBranch        *b_Elab3fitdauEta;   //!
   TBranch        *b_indexbestEtap;   //!
   TBranch        *b_chbestEtap;   //!
   TBranch        *b_mbestEtap;   //!
   TBranch        *b_barembestEtap;   //!
   TBranch        *b_mfitbestEtap;   //!
   TBranch        *b_baremfitbestEtap;   //!
   TBranch        *b_mm2bestEtap;   //!
   TBranch        *b_q2bestEtap;   //!
   TBranch        *b_q2fitbestEtap;   //!
   TBranch        *b_mm2fitbestEtap;   //!
   TBranch        *b_nrecoEtap;   //!
   TBranch        *b_mm2Etap;   //!
   TBranch        *b_q2Etap;   //!
   TBranch        *b_pEtap;   //!
   TBranch        *b_thEtap;   //!
   TBranch        *b_phiEtap;   //!
   TBranch        *b_nominalmEtap;   //!
   TBranch        *b_baremEtap;   //!
   TBranch        *b_ElabEtap;   //!
   TBranch        *b_chEtap;   //!
   TBranch        *b_mm2fitEtap;   //!
   TBranch        *b_q2fitEtap;   //!
   TBranch        *b_massfitEtap;   //!
   TBranch        *b_pfitEtap;   //!
   TBranch        *b_thetafitEtap;   //!
   TBranch        *b_phifitEtap;   //!
   TBranch        *b_nominalmfitEtap;   //!
   TBranch        *b_baremfitEtap;   //!
   TBranch        *b_ElabfitEtap;   //!
   TBranch        *b_chisqEtap;   //!
   TBranch        *b_globchisqEtap;   //!
   TBranch        *b_probchisqEtap;   //!
   TBranch        *b_ndofEtap;   //!
   TBranch        *b_pstarfitleptEtap;   //!
   TBranch        *b_pfitleptEtap;   //!
   TBranch        *b_thetafitleptEtap;   //!
   TBranch        *b_phifitleptEtap;   //!
   TBranch        *b_pfitBEtap;   //!
   TBranch        *b_massBfitEtap;   //!
   TBranch        *b_thetafitBEtap;   //!
   TBranch        *b_phifitBEtap;   //!
   TBranch        *b_ndauEtap;   //!
   TBranch        *b_modeEtap;   //!
   TBranch        *b_EtamassdauEtap;   //!
   TBranch        *b_Rho0massdauEtap;   //!
   TBranch        *b_GammamomdauEtap;   //!
   TBranch        *b_Estar1dauEtap;   //!
   TBranch        *b_Estar2dauEtap;   //!
   TBranch        *b_Estar3dauEtap;   //!
   TBranch        *b_Elab1dauEtap;   //!
   TBranch        *b_Elab2dauEtap;   //!
   TBranch        *b_Elab3dauEtap;   //!
   TBranch        *b_Estar1fitdauEtap;   //!
   TBranch        *b_Estar2fitdauEtap;   //!
   TBranch        *b_Estar3fitdauEtap;   //!
   TBranch        *b_Elab1fitdauEtap;   //!
   TBranch        *b_Elab2fitdauEtap;   //!
   TBranch        *b_Elab3fitdauEtap;   //!
   TBranch        *b_isassocB;   //!
   TBranch        *b_isassocB_GHIT;   //!
   TBranch        *b_ass_deltapB;   //!
   TBranch        *b_ch1B;   //!
   TBranch        *b_ch2B;   //!
   TBranch        *b_chunm;   //!
   TBranch        *b_neu1B;   //!
   TBranch        *b_neu2B;   //!
   TBranch        *b_neuunm;   //!
   TBranch        *b_brecoqual;   //!
   TBranch        *b_brecoqualangle;   //!
   TBranch        *b_chgdaugen;   //!
   TBranch        *b_neudaugen;   //!
   TBranch        *b_vub;   //!
   TBranch        *b_vcb;   //!
   TBranch        *b_other;   //!
   TBranch        *b_nvubexcl;   //!
   TBranch        *b_nvubnres;   //!
   TBranch        *b_mxhadgen;   //!
   TBranch        *b_mxhadgenwoph;   //!
   TBranch        *b_pcmsgen;   //!
   TBranch        *b_ecmsgen;   //!
   TBranch        *b_tcmsgen;   //!
   TBranch        *b_fcmsgen;   //!
   TBranch        *b_pxhadgen;   //!
   TBranch        *b_exhadgen;   //!
   TBranch        *b_fxhadgen;   //!
   TBranch        *b_txhadgen;   //!
   TBranch        *b_q2Gen;   //!
   TBranch        *b_ctvgen;   //!
   TBranch        *b_ctlgen;   //!
   TBranch        *b_chigen;   //!
   TBranch        *b_Gvxbtyp;   //!
   TBranch        *b_GSem;   //!
   TBranch        *b_GfDpi;   //!
   TBranch        *b_GfDpiz;   //!
   TBranch        *b_GfDk;   //!
   TBranch        *b_GfDks;   //!
   TBranch        *b_GfDkl;   //!
   TBranch        *b_GfDlep;   //!
   TBranch        *b_GfDgam;   //!
   TBranch        *b_GfDnu;   //!
   TBranch        *b_GfD0Ds;   //!
   TBranch        *b_GfDDs;   //!
   TBranch        *b_GfDkspiopio;   //!

   MakeTree(TTree *tree=0);
   virtual ~MakeTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MakeTree_cxx
MakeTree::MakeTree(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f) {
         f = new TFile("Memory Directory");
         f->cd("Rint:/");
      }
      tree = (TTree*)gDirectory->Get("ntp1");

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("ntp1","");
      chain->Add("$sl_awg_disc/francesco/test_prod/root/signal/SP-3617-BSemiExcl-R14-1.root/ntp1");
      chain->Add("$sl_awg_disc/francesco/test_prod/root/signal/SP-3617-BSemiExcl-R14-2.root/ntp1");
      chain->Add("$sl_awg_disc/francesco/test_prod/root/signal/SP-3617-BSemiExcl-R14-3.root/ntp1");
      chain->Add("$sl_awg_disc/francesco/test_prod/root/signal/SP-3617-BSemiExcl-R14-4.root/ntp1");
      chain->Add("$sl_awg_disc/francesco/test_prod/root/signal/SP-3617-BSemiExcl-R14-5.root/ntp1");
      chain->Add("$sl_awg_disc/francesco/test_prod/root/signal/SP-3617-BSemiExcl-R14-6.root/ntp1");
      chain->Add("$sl_awg_disc/francesco/test_prod/root/signal/SP-3618-BSemiExcl-R14-1.root/ntp1");
      chain->Add("$sl_awg_disc/francesco/test_prod/root/signal/SP-3618-BSemiExcl-R14-2.root/ntp1");
      chain->Add("$sl_awg_disc/francesco/test_prod/root/signal/SP-3618-BSemiExcl-R14-3.root/ntp1");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

MakeTree::~MakeTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MakeTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MakeTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->IsA() != TChain::Class()) return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MakeTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses of the tree
   // will be set. It is normaly not necessary to make changes to the
   // generated code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running with PROOF.

   // Set branch addresses
   if (tree == 0) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run",&run);
   fChain->SetBranchAddress("upper",&upper);
   fChain->SetBranchAddress("lower",&lower);
   fChain->SetBranchAddress("nchg",&nchg);
   fChain->SetBranchAddress("nneu",&nneu);
   fChain->SetBranchAddress("nB",&nB);
   fChain->SetBranchAddress("IdB",&IdB);
   fChain->SetBranchAddress("brecoflav",&brecoflav);
   fChain->SetBranchAddress("brecocharge",&brecocharge);
   fChain->SetBranchAddress("xcharge",&xcharge);
   fChain->SetBranchAddress("mode",&mode);
   fChain->SetBranchAddress("mes",&mes);
   fChain->SetBranchAddress("de",&de);
   fChain->SetBranchAddress("pB",&pB);
   fChain->SetBranchAddress("eB",&eB);
   fChain->SetBranchAddress("thetaB",&thetaB);
   fChain->SetBranchAddress("phiB",&phiB);
   fChain->SetBranchAddress("pur",&pur);
   fChain->SetBranchAddress("intpur",&intpur);
   fChain->SetBranchAddress("nle",&nle);
   fChain->SetBranchAddress("nel",&nel);
   fChain->SetBranchAddress("nmu",&nmu);
   fChain->SetBranchAddress("nkp",&nkp);
   fChain->SetBranchAddress("nks",&nks);
   fChain->SetBranchAddress("nlept500",&nlept500);
   fChain->SetBranchAddress("nelec500",&nelec500);
   fChain->SetBranchAddress("nmu500",&nmu500);
   fChain->SetBranchAddress("brems",&brems);
   fChain->SetBranchAddress("plab",&plab);
   fChain->SetBranchAddress("tlab",&tlab);
   fChain->SetBranchAddress("flab",&flab);
   fChain->SetBranchAddress("elab",&elab);
   fChain->SetBranchAddress("pcms",&pcms);
   fChain->SetBranchAddress("ecms",&ecms);
   fChain->SetBranchAddress("tcms",&tcms);
   fChain->SetBranchAddress("fcms",&fcms);
   fChain->SetBranchAddress("lcharge",&lcharge);
   fChain->SetBranchAddress("isele",&isele);
   fChain->SetBranchAddress("wdeltam",&wdeltam);
   fChain->SetBranchAddress("Eneualt",&Eneualt);
   fChain->SetBranchAddress("indexbestPi",&indexbestPi);
   fChain->SetBranchAddress("chbestPi",&chbestPi);
   fChain->SetBranchAddress("mbestPi",&mbestPi);
   fChain->SetBranchAddress("barembestPi",&barembestPi);
   fChain->SetBranchAddress("mfitbestPi",&mfitbestPi);
   fChain->SetBranchAddress("baremfitbestPi",&baremfitbestPi);
   fChain->SetBranchAddress("mm2bestPi",&mm2bestPi);
   fChain->SetBranchAddress("q2bestPi",&q2bestPi);
   fChain->SetBranchAddress("q2fitbestPi",&q2fitbestPi);
   fChain->SetBranchAddress("mm2fitbestPi",&mm2fitbestPi);
   fChain->SetBranchAddress("nrecoPi",&nrecoPi);
   fChain->SetBranchAddress("mm2Pi",mm2Pi);
   fChain->SetBranchAddress("q2Pi",q2Pi);
   fChain->SetBranchAddress("pPi",pPi);
   fChain->SetBranchAddress("thPi",thPi);
   fChain->SetBranchAddress("phiPi",phiPi);
   fChain->SetBranchAddress("nominalmPi",nominalmPi);
   fChain->SetBranchAddress("baremPi",baremPi);
   fChain->SetBranchAddress("ElabPi",ElabPi);
   fChain->SetBranchAddress("chPi",chPi);
   fChain->SetBranchAddress("mm2fitPi",mm2fitPi);
   fChain->SetBranchAddress("q2fitPi",q2fitPi);
   fChain->SetBranchAddress("massfitPi",massfitPi);
   fChain->SetBranchAddress("pfitPi",pfitPi);
   fChain->SetBranchAddress("thetafitPi",thetafitPi);
   fChain->SetBranchAddress("phifitPi",phifitPi);
   fChain->SetBranchAddress("nominalmfitPi",nominalmfitPi);
   fChain->SetBranchAddress("baremfitPi",baremfitPi);
   fChain->SetBranchAddress("ElabfitPi",ElabfitPi);
   fChain->SetBranchAddress("chisqPi",chisqPi);
   fChain->SetBranchAddress("globchisqPi",globchisqPi);
   fChain->SetBranchAddress("probchisqPi",probchisqPi);
   fChain->SetBranchAddress("ndofPi",ndofPi);
   fChain->SetBranchAddress("pstarfitleptPi",pstarfitleptPi);
   fChain->SetBranchAddress("pfitleptPi",pfitleptPi);
   fChain->SetBranchAddress("thetafitleptPi",thetafitleptPi);
   fChain->SetBranchAddress("phifitleptPi",phifitleptPi);
   fChain->SetBranchAddress("pfitBPi",pfitBPi);
   fChain->SetBranchAddress("massBfitPi",massBfitPi);
   fChain->SetBranchAddress("thetafitBPi",thetafitBPi);
   fChain->SetBranchAddress("phifitBPi",phifitBPi);
   fChain->SetBranchAddress("lepmapPi",lepmapPi);
   fChain->SetBranchAddress("massjpsiPi",massjpsiPi);
   fChain->SetBranchAddress("baremassjpsiPi",baremassjpsiPi);
   fChain->SetBranchAddress("indexbestPi0",&indexbestPi0);
   fChain->SetBranchAddress("chbestPi0",&chbestPi0);
   fChain->SetBranchAddress("mbestPi0",&mbestPi0);
   fChain->SetBranchAddress("barembestPi0",&barembestPi0);
   fChain->SetBranchAddress("mfitbestPi0",&mfitbestPi0);
   fChain->SetBranchAddress("baremfitbestPi0",&baremfitbestPi0);
   fChain->SetBranchAddress("mm2bestPi0",&mm2bestPi0);
   fChain->SetBranchAddress("q2bestPi0",&q2bestPi0);
   fChain->SetBranchAddress("q2fitbestPi0",&q2fitbestPi0);
   fChain->SetBranchAddress("mm2fitbestPi0",&mm2fitbestPi0);
   fChain->SetBranchAddress("nrecoPi0",&nrecoPi0);
   fChain->SetBranchAddress("mm2Pi0",mm2Pi0);
   fChain->SetBranchAddress("q2Pi0",q2Pi0);
   fChain->SetBranchAddress("pPi0",pPi0);
   fChain->SetBranchAddress("thPi0",thPi0);
   fChain->SetBranchAddress("phiPi0",phiPi0);
   fChain->SetBranchAddress("nominalmPi0",nominalmPi0);
   fChain->SetBranchAddress("baremPi0",baremPi0);
   fChain->SetBranchAddress("ElabPi0",ElabPi0);
   fChain->SetBranchAddress("chPi0",chPi0);
   fChain->SetBranchAddress("mm2fitPi0",mm2fitPi0);
   fChain->SetBranchAddress("q2fitPi0",q2fitPi0);
   fChain->SetBranchAddress("massfitPi0",massfitPi0);
   fChain->SetBranchAddress("pfitPi0",pfitPi0);
   fChain->SetBranchAddress("thetafitPi0",thetafitPi0);
   fChain->SetBranchAddress("phifitPi0",phifitPi0);
   fChain->SetBranchAddress("nominalmfitPi0",nominalmfitPi0);
   fChain->SetBranchAddress("baremfitPi0",baremfitPi0);
   fChain->SetBranchAddress("ElabfitPi0",ElabfitPi0);
   fChain->SetBranchAddress("chisqPi0",chisqPi0);
   fChain->SetBranchAddress("globchisqPi0",globchisqPi0);
   fChain->SetBranchAddress("probchisqPi0",probchisqPi0);
   fChain->SetBranchAddress("ndofPi0",ndofPi0);
   fChain->SetBranchAddress("pstarfitleptPi0",pstarfitleptPi0);
   fChain->SetBranchAddress("pfitleptPi0",pfitleptPi0);
   fChain->SetBranchAddress("thetafitleptPi0",thetafitleptPi0);
   fChain->SetBranchAddress("phifitleptPi0",phifitleptPi0);
   fChain->SetBranchAddress("pfitBPi0",pfitBPi0);
   fChain->SetBranchAddress("massBfitPi0",massBfitPi0);
   fChain->SetBranchAddress("thetafitBPi0",thetafitBPi0);
   fChain->SetBranchAddress("phifitBPi0",phifitBPi0);
   fChain->SetBranchAddress("ndauPi0",ndauPi0);
   fChain->SetBranchAddress("Estar1dauPi0",Estar1dauPi0);
   fChain->SetBranchAddress("Estar2dauPi0",Estar2dauPi0);
   fChain->SetBranchAddress("Elab1dauPi0",Elab1dauPi0);
   fChain->SetBranchAddress("Elab2dauPi0",Elab2dauPi0);
   fChain->SetBranchAddress("Estar1fitdauPi0",Estar1fitdauPi0);
   fChain->SetBranchAddress("Estar2fitdauPi0",Estar2fitdauPi0);
   fChain->SetBranchAddress("Elab1fitdauPi0",Elab1fitdauPi0);
   fChain->SetBranchAddress("Elab2fitdauPi0",Elab2fitdauPi0);
   fChain->SetBranchAddress("indexbestEta",&indexbestEta);
   fChain->SetBranchAddress("chbestEta",&chbestEta);
   fChain->SetBranchAddress("mbestEta",&mbestEta);
   fChain->SetBranchAddress("barembestEta",&barembestEta);
   fChain->SetBranchAddress("mfitbestEta",&mfitbestEta);
   fChain->SetBranchAddress("baremfitbestEta",&baremfitbestEta);
   fChain->SetBranchAddress("mm2bestEta",&mm2bestEta);
   fChain->SetBranchAddress("q2bestEta",&q2bestEta);
   fChain->SetBranchAddress("q2fitbestEta",&q2fitbestEta);
   fChain->SetBranchAddress("mm2fitbestEta",&mm2fitbestEta);
   fChain->SetBranchAddress("nrecoEta",&nrecoEta);
   fChain->SetBranchAddress("mm2Eta",mm2Eta);
   fChain->SetBranchAddress("q2Eta",q2Eta);
   fChain->SetBranchAddress("pEta",pEta);
   fChain->SetBranchAddress("thEta",thEta);
   fChain->SetBranchAddress("phiEta",phiEta);
   fChain->SetBranchAddress("nominalmEta",nominalmEta);
   fChain->SetBranchAddress("baremEta",baremEta);
   fChain->SetBranchAddress("ElabEta",ElabEta);
   fChain->SetBranchAddress("chEta",chEta);
   fChain->SetBranchAddress("mm2fitEta",mm2fitEta);
   fChain->SetBranchAddress("q2fitEta",q2fitEta);
   fChain->SetBranchAddress("massfitEta",massfitEta);
   fChain->SetBranchAddress("pfitEta",pfitEta);
   fChain->SetBranchAddress("thetafitEta",thetafitEta);
   fChain->SetBranchAddress("phifitEta",phifitEta);
   fChain->SetBranchAddress("nominalmfitEta",nominalmfitEta);
   fChain->SetBranchAddress("baremfitEta",baremfitEta);
   fChain->SetBranchAddress("ElabfitEta",ElabfitEta);
   fChain->SetBranchAddress("chisqEta",chisqEta);
   fChain->SetBranchAddress("globchisqEta",globchisqEta);
   fChain->SetBranchAddress("probchisqEta",probchisqEta);
   fChain->SetBranchAddress("ndofEta",ndofEta);
   fChain->SetBranchAddress("pstarfitleptEta",pstarfitleptEta);
   fChain->SetBranchAddress("pfitleptEta",pfitleptEta);
   fChain->SetBranchAddress("thetafitleptEta",thetafitleptEta);
   fChain->SetBranchAddress("phifitleptEta",phifitleptEta);
   fChain->SetBranchAddress("pfitBEta",pfitBEta);
   fChain->SetBranchAddress("massBfitEta",massBfitEta);
   fChain->SetBranchAddress("thetafitBEta",thetafitBEta);
   fChain->SetBranchAddress("phifitBEta",phifitBEta);
   fChain->SetBranchAddress("ndauEta",ndauEta);
   fChain->SetBranchAddress("modeEta",modeEta);
   fChain->SetBranchAddress("Estar1dauEta",Estar1dauEta);
   fChain->SetBranchAddress("Estar2dauEta",Estar2dauEta);
   fChain->SetBranchAddress("Estar3dauEta",Estar3dauEta);
   fChain->SetBranchAddress("Elab1dauEta",Elab1dauEta);
   fChain->SetBranchAddress("Elab2dauEta",Elab2dauEta);
   fChain->SetBranchAddress("Elab3dauEta",Elab3dauEta);
   fChain->SetBranchAddress("Estar1fitdauEta",Estar1fitdauEta);
   fChain->SetBranchAddress("Estar2fitdauEta",Estar2fitdauEta);
   fChain->SetBranchAddress("Estar3fitdauEta",Estar3fitdauEta);
   fChain->SetBranchAddress("Elab1fitdauEta",Elab1fitdauEta);
   fChain->SetBranchAddress("Elab2fitdauEta",Elab2fitdauEta);
   fChain->SetBranchAddress("Elab3fitdauEta",Elab3fitdauEta);
   fChain->SetBranchAddress("indexbestEtap",&indexbestEtap);
   fChain->SetBranchAddress("chbestEtap",&chbestEtap);
   fChain->SetBranchAddress("mbestEtap",&mbestEtap);
   fChain->SetBranchAddress("barembestEtap",&barembestEtap);
   fChain->SetBranchAddress("mfitbestEtap",&mfitbestEtap);
   fChain->SetBranchAddress("baremfitbestEtap",&baremfitbestEtap);
   fChain->SetBranchAddress("mm2bestEtap",&mm2bestEtap);
   fChain->SetBranchAddress("q2bestEtap",&q2bestEtap);
   fChain->SetBranchAddress("q2fitbestEtap",&q2fitbestEtap);
   fChain->SetBranchAddress("mm2fitbestEtap",&mm2fitbestEtap);
   fChain->SetBranchAddress("nrecoEtap",&nrecoEtap);
   fChain->SetBranchAddress("mm2Etap",mm2Etap);
   fChain->SetBranchAddress("q2Etap",q2Etap);
   fChain->SetBranchAddress("pEtap",pEtap);
   fChain->SetBranchAddress("thEtap",thEtap);
   fChain->SetBranchAddress("phiEtap",phiEtap);
   fChain->SetBranchAddress("nominalmEtap",nominalmEtap);
   fChain->SetBranchAddress("baremEtap",baremEtap);
   fChain->SetBranchAddress("ElabEtap",ElabEtap);
   fChain->SetBranchAddress("chEtap",chEtap);
   fChain->SetBranchAddress("mm2fitEtap",mm2fitEtap);
   fChain->SetBranchAddress("q2fitEtap",q2fitEtap);
   fChain->SetBranchAddress("massfitEtap",massfitEtap);
   fChain->SetBranchAddress("pfitEtap",pfitEtap);
   fChain->SetBranchAddress("thetafitEtap",thetafitEtap);
   fChain->SetBranchAddress("phifitEtap",phifitEtap);
   fChain->SetBranchAddress("nominalmfitEtap",nominalmfitEtap);
   fChain->SetBranchAddress("baremfitEtap",baremfitEtap);
   fChain->SetBranchAddress("ElabfitEtap",ElabfitEtap);
   fChain->SetBranchAddress("chisqEtap",chisqEtap);
   fChain->SetBranchAddress("globchisqEtap",globchisqEtap);
   fChain->SetBranchAddress("probchisqEtap",probchisqEtap);
   fChain->SetBranchAddress("ndofEtap",ndofEtap);
   fChain->SetBranchAddress("pstarfitleptEtap",pstarfitleptEtap);
   fChain->SetBranchAddress("pfitleptEtap",pfitleptEtap);
   fChain->SetBranchAddress("thetafitleptEtap",thetafitleptEtap);
   fChain->SetBranchAddress("phifitleptEtap",phifitleptEtap);
   fChain->SetBranchAddress("pfitBEtap",pfitBEtap);
   fChain->SetBranchAddress("massBfitEtap",massBfitEtap);
   fChain->SetBranchAddress("thetafitBEtap",thetafitBEtap);
   fChain->SetBranchAddress("phifitBEtap",phifitBEtap);
   fChain->SetBranchAddress("ndauEtap",ndauEtap);
   fChain->SetBranchAddress("modeEtap",modeEtap);
   fChain->SetBranchAddress("EtamassdauEtap",EtamassdauEtap);
   fChain->SetBranchAddress("Rho0massdauEtap",Rho0massdauEtap);
   fChain->SetBranchAddress("GammamomdauEtap",GammamomdauEtap);
   fChain->SetBranchAddress("Estar1dauEtap",Estar1dauEtap);
   fChain->SetBranchAddress("Estar2dauEtap",Estar2dauEtap);
   fChain->SetBranchAddress("Estar3dauEtap",Estar3dauEtap);
   fChain->SetBranchAddress("Elab1dauEtap",Elab1dauEtap);
   fChain->SetBranchAddress("Elab2dauEtap",Elab2dauEtap);
   fChain->SetBranchAddress("Elab3dauEtap",Elab3dauEtap);
   fChain->SetBranchAddress("Estar1fitdauEtap",Estar1fitdauEtap);
   fChain->SetBranchAddress("Estar2fitdauEtap",Estar2fitdauEtap);
   fChain->SetBranchAddress("Estar3fitdauEtap",Estar3fitdauEtap);
   fChain->SetBranchAddress("Elab1fitdauEtap",Elab1fitdauEtap);
   fChain->SetBranchAddress("Elab2fitdauEtap",Elab2fitdauEtap);
   fChain->SetBranchAddress("Elab3fitdauEtap",Elab3fitdauEtap);
   fChain->SetBranchAddress("isassocB",&isassocB);
   fChain->SetBranchAddress("isassocB_GHIT",&isassocB_GHIT);
   fChain->SetBranchAddress("ass_deltapB",&ass_deltapB);
   fChain->SetBranchAddress("ch1B",&ch1B);
   fChain->SetBranchAddress("ch2B",&ch2B);
   fChain->SetBranchAddress("chunm",&chunm);
   fChain->SetBranchAddress("neu1B",&neu1B);
   fChain->SetBranchAddress("neu2B",&neu2B);
   fChain->SetBranchAddress("neuunm",&neuunm);
   fChain->SetBranchAddress("brecoqual",&brecoqual);
   fChain->SetBranchAddress("brecoqualangle",&brecoqualangle);
   fChain->SetBranchAddress("chgdaugen",&chgdaugen);
   fChain->SetBranchAddress("neudaugen",&neudaugen);
   fChain->SetBranchAddress("vub",&vub);
   fChain->SetBranchAddress("vcb",&vcb);
   fChain->SetBranchAddress("other",&other);
   fChain->SetBranchAddress("nvubexcl",&nvubexcl);
   fChain->SetBranchAddress("nvubnres",&nvubnres);
   fChain->SetBranchAddress("mxhadgen",&mxhadgen);
   fChain->SetBranchAddress("mxhadgenwoph",&mxhadgenwoph);
   fChain->SetBranchAddress("pcmsgen",&pcmsgen);
   fChain->SetBranchAddress("ecmsgen",&ecmsgen);
   fChain->SetBranchAddress("tcmsgen",&tcmsgen);
   fChain->SetBranchAddress("fcmsgen",&fcmsgen);
   fChain->SetBranchAddress("pxhadgen",&pxhadgen);
   fChain->SetBranchAddress("exhadgen",&exhadgen);
   fChain->SetBranchAddress("fxhadgen",&fxhadgen);
   fChain->SetBranchAddress("txhadgen",&txhadgen);
   fChain->SetBranchAddress("q2Gen",&q2Gen);
   fChain->SetBranchAddress("ctvgen",&ctvgen);
   fChain->SetBranchAddress("ctlgen",&ctlgen);
   fChain->SetBranchAddress("chigen",&chigen);
   fChain->SetBranchAddress("Gvxbtyp",&Gvxbtyp);
   fChain->SetBranchAddress("GSem",&GSem);
   fChain->SetBranchAddress("GfDpi",&GfDpi);
   fChain->SetBranchAddress("GfDpiz",&GfDpiz);
   fChain->SetBranchAddress("GfDk",&GfDk);
   fChain->SetBranchAddress("GfDks",&GfDks);
   fChain->SetBranchAddress("GfDkl",&GfDkl);
   fChain->SetBranchAddress("GfDlep",&GfDlep);
   fChain->SetBranchAddress("GfDgam",&GfDgam);
   fChain->SetBranchAddress("GfDnu",&GfDnu);
   fChain->SetBranchAddress("GfD0Ds",&GfD0Ds);
   fChain->SetBranchAddress("GfDDs",&GfDDs);
   fChain->SetBranchAddress("GfDkspiopio",&GfDkspiopio);
   Notify();
}

Bool_t MakeTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. Typically here the branch pointers
   // will be retrieved. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed.

   // Get branch pointers
   b_run = fChain->GetBranch("run");
   b_upper = fChain->GetBranch("upper");
   b_lower = fChain->GetBranch("lower");
   b_nchg = fChain->GetBranch("nchg");
   b_nneu = fChain->GetBranch("nneu");
   b_nB = fChain->GetBranch("nB");
   b_IdB = fChain->GetBranch("IdB");
   b_brecoflav = fChain->GetBranch("brecoflav");
   b_brecocharge = fChain->GetBranch("brecocharge");
   b_xcharge = fChain->GetBranch("xcharge");
   b_mode = fChain->GetBranch("mode");
   b_mes = fChain->GetBranch("mes");
   b_de = fChain->GetBranch("de");
   b_pB = fChain->GetBranch("pB");
   b_eB = fChain->GetBranch("eB");
   b_thetaB = fChain->GetBranch("thetaB");
   b_phiB = fChain->GetBranch("phiB");
   b_pur = fChain->GetBranch("pur");
   b_intpur = fChain->GetBranch("intpur");
   b_nle = fChain->GetBranch("nle");
   b_nel = fChain->GetBranch("nel");
   b_nmu = fChain->GetBranch("nmu");
   b_nkp = fChain->GetBranch("nkp");
   b_nks = fChain->GetBranch("nks");
   b_nlept500 = fChain->GetBranch("nlept500");
   b_nelec500 = fChain->GetBranch("nelec500");
   b_nmu500 = fChain->GetBranch("nmu500");
   b_brems = fChain->GetBranch("brems");
   b_plab = fChain->GetBranch("plab");
   b_tlab = fChain->GetBranch("tlab");
   b_flab = fChain->GetBranch("flab");
   b_elab = fChain->GetBranch("elab");
   b_pcms = fChain->GetBranch("pcms");
   b_ecms = fChain->GetBranch("ecms");
   b_tcms = fChain->GetBranch("tcms");
   b_fcms = fChain->GetBranch("fcms");
   b_lcharge = fChain->GetBranch("lcharge");
   b_isele = fChain->GetBranch("isele");
   b_wdeltam = fChain->GetBranch("wdeltam");
   b_Eneualt = fChain->GetBranch("Eneualt");
   b_indexbestPi = fChain->GetBranch("indexbestPi");
   b_chbestPi = fChain->GetBranch("chbestPi");
   b_mbestPi = fChain->GetBranch("mbestPi");
   b_barembestPi = fChain->GetBranch("barembestPi");
   b_mfitbestPi = fChain->GetBranch("mfitbestPi");
   b_baremfitbestPi = fChain->GetBranch("baremfitbestPi");
   b_mm2bestPi = fChain->GetBranch("mm2bestPi");
   b_q2bestPi = fChain->GetBranch("q2bestPi");
   b_q2fitbestPi = fChain->GetBranch("q2fitbestPi");
   b_mm2fitbestPi = fChain->GetBranch("mm2fitbestPi");
   b_nrecoPi = fChain->GetBranch("nrecoPi");
   b_mm2Pi = fChain->GetBranch("mm2Pi");
   b_q2Pi = fChain->GetBranch("q2Pi");
   b_pPi = fChain->GetBranch("pPi");
   b_thPi = fChain->GetBranch("thPi");
   b_phiPi = fChain->GetBranch("phiPi");
   b_nominalmPi = fChain->GetBranch("nominalmPi");
   b_baremPi = fChain->GetBranch("baremPi");
   b_ElabPi = fChain->GetBranch("ElabPi");
   b_chPi = fChain->GetBranch("chPi");
   b_mm2fitPi = fChain->GetBranch("mm2fitPi");
   b_q2fitPi = fChain->GetBranch("q2fitPi");
   b_massfitPi = fChain->GetBranch("massfitPi");
   b_pfitPi = fChain->GetBranch("pfitPi");
   b_thetafitPi = fChain->GetBranch("thetafitPi");
   b_phifitPi = fChain->GetBranch("phifitPi");
   b_nominalmfitPi = fChain->GetBranch("nominalmfitPi");
   b_baremfitPi = fChain->GetBranch("baremfitPi");
   b_ElabfitPi = fChain->GetBranch("ElabfitPi");
   b_chisqPi = fChain->GetBranch("chisqPi");
   b_globchisqPi = fChain->GetBranch("globchisqPi");
   b_probchisqPi = fChain->GetBranch("probchisqPi");
   b_ndofPi = fChain->GetBranch("ndofPi");
   b_pstarfitleptPi = fChain->GetBranch("pstarfitleptPi");
   b_pfitleptPi = fChain->GetBranch("pfitleptPi");
   b_thetafitleptPi = fChain->GetBranch("thetafitleptPi");
   b_phifitleptPi = fChain->GetBranch("phifitleptPi");
   b_pfitBPi = fChain->GetBranch("pfitBPi");
   b_massBfitPi = fChain->GetBranch("massBfitPi");
   b_thetafitBPi = fChain->GetBranch("thetafitBPi");
   b_phifitBPi = fChain->GetBranch("phifitBPi");
   b_lepmapPi = fChain->GetBranch("lepmapPi");
   b_massjpsiPi = fChain->GetBranch("massjpsiPi");
   b_baremassjpsiPi = fChain->GetBranch("baremassjpsiPi");
   b_indexbestPi0 = fChain->GetBranch("indexbestPi0");
   b_chbestPi0 = fChain->GetBranch("chbestPi0");
   b_mbestPi0 = fChain->GetBranch("mbestPi0");
   b_barembestPi0 = fChain->GetBranch("barembestPi0");
   b_mfitbestPi0 = fChain->GetBranch("mfitbestPi0");
   b_baremfitbestPi0 = fChain->GetBranch("baremfitbestPi0");
   b_mm2bestPi0 = fChain->GetBranch("mm2bestPi0");
   b_q2bestPi0 = fChain->GetBranch("q2bestPi0");
   b_q2fitbestPi0 = fChain->GetBranch("q2fitbestPi0");
   b_mm2fitbestPi0 = fChain->GetBranch("mm2fitbestPi0");
   b_nrecoPi0 = fChain->GetBranch("nrecoPi0");
   b_mm2Pi0 = fChain->GetBranch("mm2Pi0");
   b_q2Pi0 = fChain->GetBranch("q2Pi0");
   b_pPi0 = fChain->GetBranch("pPi0");
   b_thPi0 = fChain->GetBranch("thPi0");
   b_phiPi0 = fChain->GetBranch("phiPi0");
   b_nominalmPi0 = fChain->GetBranch("nominalmPi0");
   b_baremPi0 = fChain->GetBranch("baremPi0");
   b_ElabPi0 = fChain->GetBranch("ElabPi0");
   b_chPi0 = fChain->GetBranch("chPi0");
   b_mm2fitPi0 = fChain->GetBranch("mm2fitPi0");
   b_q2fitPi0 = fChain->GetBranch("q2fitPi0");
   b_massfitPi0 = fChain->GetBranch("massfitPi0");
   b_pfitPi0 = fChain->GetBranch("pfitPi0");
   b_thetafitPi0 = fChain->GetBranch("thetafitPi0");
   b_phifitPi0 = fChain->GetBranch("phifitPi0");
   b_nominalmfitPi0 = fChain->GetBranch("nominalmfitPi0");
   b_baremfitPi0 = fChain->GetBranch("baremfitPi0");
   b_ElabfitPi0 = fChain->GetBranch("ElabfitPi0");
   b_chisqPi0 = fChain->GetBranch("chisqPi0");
   b_globchisqPi0 = fChain->GetBranch("globchisqPi0");
   b_probchisqPi0 = fChain->GetBranch("probchisqPi0");
   b_ndofPi0 = fChain->GetBranch("ndofPi0");
   b_pstarfitleptPi0 = fChain->GetBranch("pstarfitleptPi0");
   b_pfitleptPi0 = fChain->GetBranch("pfitleptPi0");
   b_thetafitleptPi0 = fChain->GetBranch("thetafitleptPi0");
   b_phifitleptPi0 = fChain->GetBranch("phifitleptPi0");
   b_pfitBPi0 = fChain->GetBranch("pfitBPi0");
   b_massBfitPi0 = fChain->GetBranch("massBfitPi0");
   b_thetafitBPi0 = fChain->GetBranch("thetafitBPi0");
   b_phifitBPi0 = fChain->GetBranch("phifitBPi0");
   b_ndauPi0 = fChain->GetBranch("ndauPi0");
   b_Estar1dauPi0 = fChain->GetBranch("Estar1dauPi0");
   b_Estar2dauPi0 = fChain->GetBranch("Estar2dauPi0");
   b_Elab1dauPi0 = fChain->GetBranch("Elab1dauPi0");
   b_Elab2dauPi0 = fChain->GetBranch("Elab2dauPi0");
   b_Estar1fitdauPi0 = fChain->GetBranch("Estar1fitdauPi0");
   b_Estar2fitdauPi0 = fChain->GetBranch("Estar2fitdauPi0");
   b_Elab1fitdauPi0 = fChain->GetBranch("Elab1fitdauPi0");
   b_Elab2fitdauPi0 = fChain->GetBranch("Elab2fitdauPi0");
   b_indexbestEta = fChain->GetBranch("indexbestEta");
   b_chbestEta = fChain->GetBranch("chbestEta");
   b_mbestEta = fChain->GetBranch("mbestEta");
   b_barembestEta = fChain->GetBranch("barembestEta");
   b_mfitbestEta = fChain->GetBranch("mfitbestEta");
   b_baremfitbestEta = fChain->GetBranch("baremfitbestEta");
   b_mm2bestEta = fChain->GetBranch("mm2bestEta");
   b_q2bestEta = fChain->GetBranch("q2bestEta");
   b_q2fitbestEta = fChain->GetBranch("q2fitbestEta");
   b_mm2fitbestEta = fChain->GetBranch("mm2fitbestEta");
   b_nrecoEta = fChain->GetBranch("nrecoEta");
   b_mm2Eta = fChain->GetBranch("mm2Eta");
   b_q2Eta = fChain->GetBranch("q2Eta");
   b_pEta = fChain->GetBranch("pEta");
   b_thEta = fChain->GetBranch("thEta");
   b_phiEta = fChain->GetBranch("phiEta");
   b_nominalmEta = fChain->GetBranch("nominalmEta");
   b_baremEta = fChain->GetBranch("baremEta");
   b_ElabEta = fChain->GetBranch("ElabEta");
   b_chEta = fChain->GetBranch("chEta");
   b_mm2fitEta = fChain->GetBranch("mm2fitEta");
   b_q2fitEta = fChain->GetBranch("q2fitEta");
   b_massfitEta = fChain->GetBranch("massfitEta");
   b_pfitEta = fChain->GetBranch("pfitEta");
   b_thetafitEta = fChain->GetBranch("thetafitEta");
   b_phifitEta = fChain->GetBranch("phifitEta");
   b_nominalmfitEta = fChain->GetBranch("nominalmfitEta");
   b_baremfitEta = fChain->GetBranch("baremfitEta");
   b_ElabfitEta = fChain->GetBranch("ElabfitEta");
   b_chisqEta = fChain->GetBranch("chisqEta");
   b_globchisqEta = fChain->GetBranch("globchisqEta");
   b_probchisqEta = fChain->GetBranch("probchisqEta");
   b_ndofEta = fChain->GetBranch("ndofEta");
   b_pstarfitleptEta = fChain->GetBranch("pstarfitleptEta");
   b_pfitleptEta = fChain->GetBranch("pfitleptEta");
   b_thetafitleptEta = fChain->GetBranch("thetafitleptEta");
   b_phifitleptEta = fChain->GetBranch("phifitleptEta");
   b_pfitBEta = fChain->GetBranch("pfitBEta");
   b_massBfitEta = fChain->GetBranch("massBfitEta");
   b_thetafitBEta = fChain->GetBranch("thetafitBEta");
   b_phifitBEta = fChain->GetBranch("phifitBEta");
   b_ndauEta = fChain->GetBranch("ndauEta");
   b_modeEta = fChain->GetBranch("modeEta");
   b_Estar1dauEta = fChain->GetBranch("Estar1dauEta");
   b_Estar2dauEta = fChain->GetBranch("Estar2dauEta");
   b_Estar3dauEta = fChain->GetBranch("Estar3dauEta");
   b_Elab1dauEta = fChain->GetBranch("Elab1dauEta");
   b_Elab2dauEta = fChain->GetBranch("Elab2dauEta");
   b_Elab3dauEta = fChain->GetBranch("Elab3dauEta");
   b_Estar1fitdauEta = fChain->GetBranch("Estar1fitdauEta");
   b_Estar2fitdauEta = fChain->GetBranch("Estar2fitdauEta");
   b_Estar3fitdauEta = fChain->GetBranch("Estar3fitdauEta");
   b_Elab1fitdauEta = fChain->GetBranch("Elab1fitdauEta");
   b_Elab2fitdauEta = fChain->GetBranch("Elab2fitdauEta");
   b_Elab3fitdauEta = fChain->GetBranch("Elab3fitdauEta");
   b_indexbestEtap = fChain->GetBranch("indexbestEtap");
   b_chbestEtap = fChain->GetBranch("chbestEtap");
   b_mbestEtap = fChain->GetBranch("mbestEtap");
   b_barembestEtap = fChain->GetBranch("barembestEtap");
   b_mfitbestEtap = fChain->GetBranch("mfitbestEtap");
   b_baremfitbestEtap = fChain->GetBranch("baremfitbestEtap");
   b_mm2bestEtap = fChain->GetBranch("mm2bestEtap");
   b_q2bestEtap = fChain->GetBranch("q2bestEtap");
   b_q2fitbestEtap = fChain->GetBranch("q2fitbestEtap");
   b_mm2fitbestEtap = fChain->GetBranch("mm2fitbestEtap");
   b_nrecoEtap = fChain->GetBranch("nrecoEtap");
   b_mm2Etap = fChain->GetBranch("mm2Etap");
   b_q2Etap = fChain->GetBranch("q2Etap");
   b_pEtap = fChain->GetBranch("pEtap");
   b_thEtap = fChain->GetBranch("thEtap");
   b_phiEtap = fChain->GetBranch("phiEtap");
   b_nominalmEtap = fChain->GetBranch("nominalmEtap");
   b_baremEtap = fChain->GetBranch("baremEtap");
   b_ElabEtap = fChain->GetBranch("ElabEtap");
   b_chEtap = fChain->GetBranch("chEtap");
   b_mm2fitEtap = fChain->GetBranch("mm2fitEtap");
   b_q2fitEtap = fChain->GetBranch("q2fitEtap");
   b_massfitEtap = fChain->GetBranch("massfitEtap");
   b_pfitEtap = fChain->GetBranch("pfitEtap");
   b_thetafitEtap = fChain->GetBranch("thetafitEtap");
   b_phifitEtap = fChain->GetBranch("phifitEtap");
   b_nominalmfitEtap = fChain->GetBranch("nominalmfitEtap");
   b_baremfitEtap = fChain->GetBranch("baremfitEtap");
   b_ElabfitEtap = fChain->GetBranch("ElabfitEtap");
   b_chisqEtap = fChain->GetBranch("chisqEtap");
   b_globchisqEtap = fChain->GetBranch("globchisqEtap");
   b_probchisqEtap = fChain->GetBranch("probchisqEtap");
   b_ndofEtap = fChain->GetBranch("ndofEtap");
   b_pstarfitleptEtap = fChain->GetBranch("pstarfitleptEtap");
   b_pfitleptEtap = fChain->GetBranch("pfitleptEtap");
   b_thetafitleptEtap = fChain->GetBranch("thetafitleptEtap");
   b_phifitleptEtap = fChain->GetBranch("phifitleptEtap");
   b_pfitBEtap = fChain->GetBranch("pfitBEtap");
   b_massBfitEtap = fChain->GetBranch("massBfitEtap");
   b_thetafitBEtap = fChain->GetBranch("thetafitBEtap");
   b_phifitBEtap = fChain->GetBranch("phifitBEtap");
   b_ndauEtap = fChain->GetBranch("ndauEtap");
   b_modeEtap = fChain->GetBranch("modeEtap");
   b_EtamassdauEtap = fChain->GetBranch("EtamassdauEtap");
   b_Rho0massdauEtap = fChain->GetBranch("Rho0massdauEtap");
   b_GammamomdauEtap = fChain->GetBranch("GammamomdauEtap");
   b_Estar1dauEtap = fChain->GetBranch("Estar1dauEtap");
   b_Estar2dauEtap = fChain->GetBranch("Estar2dauEtap");
   b_Estar3dauEtap = fChain->GetBranch("Estar3dauEtap");
   b_Elab1dauEtap = fChain->GetBranch("Elab1dauEtap");
   b_Elab2dauEtap = fChain->GetBranch("Elab2dauEtap");
   b_Elab3dauEtap = fChain->GetBranch("Elab3dauEtap");
   b_Estar1fitdauEtap = fChain->GetBranch("Estar1fitdauEtap");
   b_Estar2fitdauEtap = fChain->GetBranch("Estar2fitdauEtap");
   b_Estar3fitdauEtap = fChain->GetBranch("Estar3fitdauEtap");
   b_Elab1fitdauEtap = fChain->GetBranch("Elab1fitdauEtap");
   b_Elab2fitdauEtap = fChain->GetBranch("Elab2fitdauEtap");
   b_Elab3fitdauEtap = fChain->GetBranch("Elab3fitdauEtap");
   b_isassocB = fChain->GetBranch("isassocB");
   b_isassocB_GHIT = fChain->GetBranch("isassocB_GHIT");
   b_ass_deltapB = fChain->GetBranch("ass_deltapB");
   b_ch1B = fChain->GetBranch("ch1B");
   b_ch2B = fChain->GetBranch("ch2B");
   b_chunm = fChain->GetBranch("chunm");
   b_neu1B = fChain->GetBranch("neu1B");
   b_neu2B = fChain->GetBranch("neu2B");
   b_neuunm = fChain->GetBranch("neuunm");
   b_brecoqual = fChain->GetBranch("brecoqual");
   b_brecoqualangle = fChain->GetBranch("brecoqualangle");
   b_chgdaugen = fChain->GetBranch("chgdaugen");
   b_neudaugen = fChain->GetBranch("neudaugen");
   b_vub = fChain->GetBranch("vub");
   b_vcb = fChain->GetBranch("vcb");
   b_other = fChain->GetBranch("other");
   b_nvubexcl = fChain->GetBranch("nvubexcl");
   b_nvubnres = fChain->GetBranch("nvubnres");
   b_mxhadgen = fChain->GetBranch("mxhadgen");
   b_mxhadgenwoph = fChain->GetBranch("mxhadgenwoph");
   b_pcmsgen = fChain->GetBranch("pcmsgen");
   b_ecmsgen = fChain->GetBranch("ecmsgen");
   b_tcmsgen = fChain->GetBranch("tcmsgen");
   b_fcmsgen = fChain->GetBranch("fcmsgen");
   b_pxhadgen = fChain->GetBranch("pxhadgen");
   b_exhadgen = fChain->GetBranch("exhadgen");
   b_fxhadgen = fChain->GetBranch("fxhadgen");
   b_txhadgen = fChain->GetBranch("txhadgen");
   b_q2Gen = fChain->GetBranch("q2Gen");
   b_ctvgen = fChain->GetBranch("ctvgen");
   b_ctlgen = fChain->GetBranch("ctlgen");
   b_chigen = fChain->GetBranch("chigen");
   b_Gvxbtyp = fChain->GetBranch("Gvxbtyp");
   b_GSem = fChain->GetBranch("GSem");
   b_GfDpi = fChain->GetBranch("GfDpi");
   b_GfDpiz = fChain->GetBranch("GfDpiz");
   b_GfDk = fChain->GetBranch("GfDk");
   b_GfDks = fChain->GetBranch("GfDks");
   b_GfDkl = fChain->GetBranch("GfDkl");
   b_GfDlep = fChain->GetBranch("GfDlep");
   b_GfDgam = fChain->GetBranch("GfDgam");
   b_GfDnu = fChain->GetBranch("GfDnu");
   b_GfD0Ds = fChain->GetBranch("GfD0Ds");
   b_GfDDs = fChain->GetBranch("GfDDs");
   b_GfDkspiopio = fChain->GetBranch("GfDkspiopio");

   return kTRUE;
}

void MakeTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MakeTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MakeTree_cxx
