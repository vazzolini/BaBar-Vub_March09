//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 27 04:20:29 2006 by ROOT version 4.04/02b
// from TChain ntp1/
//////////////////////////////////////////////////////////

#ifndef MesPlots_h
#define MesPlots_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class MesPlots {
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
   Int_t           indexbestRho0;
   Int_t           chbestRho0;
   Double_t        mbestRho0;
   Double_t        barembestRho0;
   Double_t        mfitbestRho0;
   Double_t        baremfitbestRho0;
   Double_t        mm2bestRho0;
   Double_t        q2bestRho0;
   Double_t        q2fitbestRho0;
   Double_t        mm2fitbestRho0;
   Int_t           nrecoRho0;
   Float_t         mm2Rho0[100];   //[nrecoRho0]
   Float_t         q2Rho0[100];   //[nrecoRho0]
   Float_t         pRho0[100];   //[nrecoRho0]
   Float_t         thRho0[100];   //[nrecoRho0]
   Float_t         phiRho0[100];   //[nrecoRho0]
   Float_t         nominalmRho0[100];   //[nrecoRho0]
   Float_t         baremRho0[100];   //[nrecoRho0]
   Float_t         ElabRho0[100];   //[nrecoRho0]
   Int_t           chRho0[100];   //[nrecoRho0]
   Float_t         mm2fitRho0[100];   //[nrecoRho0]
   Float_t         q2fitRho0[100];   //[nrecoRho0]
   Float_t         massfitRho0[100];   //[nrecoRho0]
   Float_t         pfitRho0[100];   //[nrecoRho0]
   Float_t         thetafitRho0[100];   //[nrecoRho0]
   Float_t         phifitRho0[100];   //[nrecoRho0]
   Float_t         nominalmfitRho0[100];   //[nrecoRho0]
   Float_t         baremfitRho0[100];   //[nrecoRho0]
   Float_t         ElabfitRho0[100];   //[nrecoRho0]
   Float_t         chisqRho0[100];   //[nrecoRho0]
   Float_t         globchisqRho0[100];   //[nrecoRho0]
   Float_t         probchisqRho0[100];   //[nrecoRho0]
   Int_t           ndofRho0[100];   //[nrecoRho0]
   Float_t         pstarfitleptRho0[100];   //[nrecoRho0]
   Float_t         pfitleptRho0[100];   //[nrecoRho0]
   Float_t         thetafitleptRho0[100];   //[nrecoRho0]
   Float_t         phifitleptRho0[100];   //[nrecoRho0]
   Float_t         pfitBRho0[100];   //[nrecoRho0]
   Float_t         massBfitRho0[100];   //[nrecoRho0]
   Float_t         thetafitBRho0[100];   //[nrecoRho0]
   Float_t         phifitBRho0[100];   //[nrecoRho0]
   Int_t           ndauRho0[100];   //[nrecoRho0]
   Float_t         Estar1dauRho0[100];   //[nrecoRho0]
   Float_t         Estar2dauRho0[100];   //[nrecoRho0]
   Float_t         Elab1dauRho0[100];   //[nrecoRho0]
   Float_t         Elab2dauRho0[100];   //[nrecoRho0]
   Float_t         Estar1fitdauRho0[100];   //[nrecoRho0]
   Float_t         Estar2fitdauRho0[100];   //[nrecoRho0]
   Float_t         Elab1fitdauRho0[100];   //[nrecoRho0]
   Float_t         Elab2fitdauRho0[100];   //[nrecoRho0]
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
   TBranch        *b_indexbestRho0;   //!
   TBranch        *b_chbestRho0;   //!
   TBranch        *b_mbestRho0;   //!
   TBranch        *b_barembestRho0;   //!
   TBranch        *b_mfitbestRho0;   //!
   TBranch        *b_baremfitbestRho0;   //!
   TBranch        *b_mm2bestRho0;   //!
   TBranch        *b_q2bestRho0;   //!
   TBranch        *b_q2fitbestRho0;   //!
   TBranch        *b_mm2fitbestRho0;   //!
   TBranch        *b_nrecoRho0;   //!
   TBranch        *b_mm2Rho0;   //!
   TBranch        *b_q2Rho0;   //!
   TBranch        *b_pRho0;   //!
   TBranch        *b_thRho0;   //!
   TBranch        *b_phiRho0;   //!
   TBranch        *b_nominalmRho0;   //!
   TBranch        *b_baremRho0;   //!
   TBranch        *b_ElabRho0;   //!
   TBranch        *b_chRho0;   //!
   TBranch        *b_mm2fitRho0;   //!
   TBranch        *b_q2fitRho0;   //!
   TBranch        *b_massfitRho0;   //!
   TBranch        *b_pfitRho0;   //!
   TBranch        *b_thetafitRho0;   //!
   TBranch        *b_phifitRho0;   //!
   TBranch        *b_nominalmfitRho0;   //!
   TBranch        *b_baremfitRho0;   //!
   TBranch        *b_ElabfitRho0;   //!
   TBranch        *b_chisqRho0;   //!
   TBranch        *b_globchisqRho0;   //!
   TBranch        *b_probchisqRho0;   //!
   TBranch        *b_ndofRho0;   //!
   TBranch        *b_pstarfitleptRho0;   //!
   TBranch        *b_pfitleptRho0;   //!
   TBranch        *b_thetafitleptRho0;   //!
   TBranch        *b_phifitleptRho0;   //!
   TBranch        *b_pfitBRho0;   //!
   TBranch        *b_massBfitRho0;   //!
   TBranch        *b_thetafitBRho0;   //!
   TBranch        *b_phifitBRho0;   //!
   TBranch        *b_ndauRho0;   //!
   TBranch        *b_Estar1dauRho0;   //!
   TBranch        *b_Estar2dauRho0;   //!
   TBranch        *b_Elab1dauRho0;   //!
   TBranch        *b_Elab2dauRho0;   //!
   TBranch        *b_Estar1fitdauRho0;   //!
   TBranch        *b_Estar2fitdauRho0;   //!
   TBranch        *b_Elab1fitdauRho0;   //!
   TBranch        *b_Elab2fitdauRho0;   //!
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

   MesPlots(TTree *tree=0);
   virtual ~MesPlots();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MesPlots_cxx
MesPlots::MesPlots(TTree *tree)
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
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-1.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-10.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-100.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-101.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-102.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-103.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-104.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-105.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-106.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-107.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-108.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-109.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-11.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-110.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-111.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-112.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-113.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-114.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-115.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-116.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-12.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-13.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-14.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-15.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-16.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-17.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-18.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-19.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-2.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-20.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-21.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-22.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-23.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-24.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-25.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-26.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-27.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-28.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-29.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-3.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-30.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-31.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-32.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-33.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-34.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-35.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-36.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-37.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-38.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-39.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-4.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-40.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-41.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-42.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-43.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-44.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-45.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-46.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-47.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-48.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-49.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-5.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-50.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-51.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-52.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-53.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-54.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-55.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-56.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-57.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-58.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-59.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-6.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-60.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-61.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-62.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-63.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-64.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-65.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-66.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-67.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-68.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-69.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-7.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-70.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-71.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-72.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-73.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-74.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-75.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-76.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-77.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-78.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-79.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-8.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-80.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-81.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-82.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-83.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-84.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-85.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-86.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-87.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-88.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-89.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-9.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-90.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-91.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-92.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-93.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-94.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-95.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-96.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-97.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-98.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v03-99.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v04-1.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v04-2.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v04-3.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run1-OnPeak-R18b-v04-4.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-1.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-10.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-100.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-101.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-102.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-103.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-104.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-105.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-106.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-107.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-108.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-109.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-11.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-110.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-111.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-112.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-113.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-114.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-115.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-116.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-117.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-118.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-119.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-12.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-120.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-121.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-122.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-123.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-124.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-125.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-126.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-127.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-128.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-129.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-13.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-130.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-131.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-132.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-133.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-134.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-135.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-136.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-137.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-138.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-139.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-14.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-140.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-141.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-142.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-143.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-144.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-145.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-146.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-147.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-148.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-149.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-15.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-150.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-151.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-152.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-153.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-154.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-155.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-156.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-157.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-158.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-159.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-16.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-160.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-161.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-162.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-163.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-164.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-165.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-166.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-167.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-168.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-169.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-17.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-170.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-171.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-172.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-173.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-174.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-175.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-176.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-177.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-178.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-179.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-18.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-180.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-181.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-182.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-183.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-184.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-185.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-186.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-187.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-188.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-189.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-19.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-190.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-191.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-192.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-193.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-194.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-195.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-196.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-197.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-198.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-199.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-2.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-20.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-200.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-201.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-202.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-203.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-204.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-205.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-206.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-207.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-208.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-209.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-21.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-210.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-211.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-212.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-213.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-214.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-215.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-216.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-217.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-218.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-219.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-22.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-220.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-221.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-222.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-223.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-224.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-225.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-226.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-227.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-228.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-229.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-23.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-230.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-231.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-232.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-233.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-234.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-235.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-236.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-237.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-238.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-239.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-24.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-240.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-241.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-242.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-243.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-244.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-245.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-246.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-247.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-248.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-249.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-25.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-250.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-251.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-252.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-253.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-254.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-255.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-256.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-257.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-258.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-259.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-26.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-260.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-261.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-262.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-263.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-264.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-265.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-266.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-267.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-268.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-269.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-27.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-270.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-271.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-272.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-273.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-274.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-275.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-276.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-277.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-278.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-279.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-28.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-280.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-281.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-282.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-283.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-284.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-285.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-286.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-287.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-288.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-289.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-29.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-290.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-291.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-292.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-293.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-294.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-295.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-296.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-297.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-298.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-299.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-3.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-30.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-300.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-301.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-302.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-303.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-304.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-305.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-306.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-307.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-308.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-309.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-31.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-310.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-311.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-312.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-313.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-314.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-315.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-316.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-317.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-318.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-319.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-32.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-320.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-321.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-322.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-323.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-324.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-325.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-326.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-327.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-328.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-329.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-33.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-330.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-331.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-332.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-333.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-334.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-335.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-336.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-337.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-338.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-339.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-34.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-340.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-341.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-342.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-343.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-344.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-345.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-346.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-347.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-348.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-349.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-35.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-350.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-351.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-352.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-353.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-354.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-355.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-356.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-357.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-36.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-37.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-38.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-39.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-4.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-40.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-41.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-42.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-43.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-44.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-45.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-46.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-47.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-48.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-49.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-5.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-50.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-51.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-52.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-53.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-54.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-55.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-56.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-57.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-58.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-59.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-6.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-60.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-61.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-62.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-63.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-64.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-65.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-66.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-67.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-68.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-69.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-7.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-70.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-71.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-72.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-73.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-74.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-75.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-76.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-77.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-78.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-79.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-8.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-80.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-81.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-82.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-83.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-84.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-85.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-86.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-87.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-88.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-89.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-9.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-90.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-91.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-92.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-93.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-94.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-95.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-96.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-97.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-98.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v03-99.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v04-1.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run2-OnPeak-R18b-v04-2.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-1.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-100.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-101.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-102.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-103.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-104.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-105.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-106.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-107.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-108.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-109.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-110.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-111.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-112.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-113.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-114.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-115.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-116.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-117.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-118.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-119.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-120.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-121.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-122.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-123.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-124.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-125.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-126.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-127.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-128.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-129.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-130.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-131.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-132.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-133.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-134.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-135.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-136.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-137.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-138.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-139.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-140.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-141.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-142.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-143.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-144.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-145.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-146.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-147.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-148.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-149.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-150.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-151.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-152.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-153.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-154.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-155.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-156.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-157.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-158.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-159.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-160.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-161.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-162.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-163.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-164.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-165.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-166.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-167.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-168.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-169.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-17.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-170.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-171.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-172.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-173.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-174.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-175.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-176.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-177.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-178.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-179.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-18.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-180.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-181.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-182.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-183.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-184.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-185.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-186.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-187.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-188.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-189.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-19.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-190.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-191.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-192.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-193.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-2.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-20.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-21.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-22.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-23.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-24.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-25.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-26.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-27.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-28.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-29.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-3.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-30.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-31.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-32.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-33.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-34.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-35.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-36.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-37.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-38.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-39.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-4.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-40.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-41.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-42.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-43.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-44.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-45.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-46.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-47.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-48.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-49.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-5.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-50.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-51.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-52.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-53.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-54.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-55.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-56.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-57.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-58.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-59.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-6.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-60.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-61.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-62.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-63.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-64.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-65.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-66.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-67.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-68.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-69.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-7.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-70.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-71.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-72.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-73.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-74.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-75.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-76.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-77.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-78.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-79.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-80.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-81.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-82.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-83.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-84.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-85.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-86.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-95.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-96.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-97.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-98.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v03-99.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v04-1.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v04-10.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v04-11.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v04-12.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v04-13.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v04-14.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v04-15.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v04-16.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v04-2.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v04-3.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v04-4.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v04-5.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v04-6.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v04-7.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v04-8.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run3-OnPeak-R18b-v04-9.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-1.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-10.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-100.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-101.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-102.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-103.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-104.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-105.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-106.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-107.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-108.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-109.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-11.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-110.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-111.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-112.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-113.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-114.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-115.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-116.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-117.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-118.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-119.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-12.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-120.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-121.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-122.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-123.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-124.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-125.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-126.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-127.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-128.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-129.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-13.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-130.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-131.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-132.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-133.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-134.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-135.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-136.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-137.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-138.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-139.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-14.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-140.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-141.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-142.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-143.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-144.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-145.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-146.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-147.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-148.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-149.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-15.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-150.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-151.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-152.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-153.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-154.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-155.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-156.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-157.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-158.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-159.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-16.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-160.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-161.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-162.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-163.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-164.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-165.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-166.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-167.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-168.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-169.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-17.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-170.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-171.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-172.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-173.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-174.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-175.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-176.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-177.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-178.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-179.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-18.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-180.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-181.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-182.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-183.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-184.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-185.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-186.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-187.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-188.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-189.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-19.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-190.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-191.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-192.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-193.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-194.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-195.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-196.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-197.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-198.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-199.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-2.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-20.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-200.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-201.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-202.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-203.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-204.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-205.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-206.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-207.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-208.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-209.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-21.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-210.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-211.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-212.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-213.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-214.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-215.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-216.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-217.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-218.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-219.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-22.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-220.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-221.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-222.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-223.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-224.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-225.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-226.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-227.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-228.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-229.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-23.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-230.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-231.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-232.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-233.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-234.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-235.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-236.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-237.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-238.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-239.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-24.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-240.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-241.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-242.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-243.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-244.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-245.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-246.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-247.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-248.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-249.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-25.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-250.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-251.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-252.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-253.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-254.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-255.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-256.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-257.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-258.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-259.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-26.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-260.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-261.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-262.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-263.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-264.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-265.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-266.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-267.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-268.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-269.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-27.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-270.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-271.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-272.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-273.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-274.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-275.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-276.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-277.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-278.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-279.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-28.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-280.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-281.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-282.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-283.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-284.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-285.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-286.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-287.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-288.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-289.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-29.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-290.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-291.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-292.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-293.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-294.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-295.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-296.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-297.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-298.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-299.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-3.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-30.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-300.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-301.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-302.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-303.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-304.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-305.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-306.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-307.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-308.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-309.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-31.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-310.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-311.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-312.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-313.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-314.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-315.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-316.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-317.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-318.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-319.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-32.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-320.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-321.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-322.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-323.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-324.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-325.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-326.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-327.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-328.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-329.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-33.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-330.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-331.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-332.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-333.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-334.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-335.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-336.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-337.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-338.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-339.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-34.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-340.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-341.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-342.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-343.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-344.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-345.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-346.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-347.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-348.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-349.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-35.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-350.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-351.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-352.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-353.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-354.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-355.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-356.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-357.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-358.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-359.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-36.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-360.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-361.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-362.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-363.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-364.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-365.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-366.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-367.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-368.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-369.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-37.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-370.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-371.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-372.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-373.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-374.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-375.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-376.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-377.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-378.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-379.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-38.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-380.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-381.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-382.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-383.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-384.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-385.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-386.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-387.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-388.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-389.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-39.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-390.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-391.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-392.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-393.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-394.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-395.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-396.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-397.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-398.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-399.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-4.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-40.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-400.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-401.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-402.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-403.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-404.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-405.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-406.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-407.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-408.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-409.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-41.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-410.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-411.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-412.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-413.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-414.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-415.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-416.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-417.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-418.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-419.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-42.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-420.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-421.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-422.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-423.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-424.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-425.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-426.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-427.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-428.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-429.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-43.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-430.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-431.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-432.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-433.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-434.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-435.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-436.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-437.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-438.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-439.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-44.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-440.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-441.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-442.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-443.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-444.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-445.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-446.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-447.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-448.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-449.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-45.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-450.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-451.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-452.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-453.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-454.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-455.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-456.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-457.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-458.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-459.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-46.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-460.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-461.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-462.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-463.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-464.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-465.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-466.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-467.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-468.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-469.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-47.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-470.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-471.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-472.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-473.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-474.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-475.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-476.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-477.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-478.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-479.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-48.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-480.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-481.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-482.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-483.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-484.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-485.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-486.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-487.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-488.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-489.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-49.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-490.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-491.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-492.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-493.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-494.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-495.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-496.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-497.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-498.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-499.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-5.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-50.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-500.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-501.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-502.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-503.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-504.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-505.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-506.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-507.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-508.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-509.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-51.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-510.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-511.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-512.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-513.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-514.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-515.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-516.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-517.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-518.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-519.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-52.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-520.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-521.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-522.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-523.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-524.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-525.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-526.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-527.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-528.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-529.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-53.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-530.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-531.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-532.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-533.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-534.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-535.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-536.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-537.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-538.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-539.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-54.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-540.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-541.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-542.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-543.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-544.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-545.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-546.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-547.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-548.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-549.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-55.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-550.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-551.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-552.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-553.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-554.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-555.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-556.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-557.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-558.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-559.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-56.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-560.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-561.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-562.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-563.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-564.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-565.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-566.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-567.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-568.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-569.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-57.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-570.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-571.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-572.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-573.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-574.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-575.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-576.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-577.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-578.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-579.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-58.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-580.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-581.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-582.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-583.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-584.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-585.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-586.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-587.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-588.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-589.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-59.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-590.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-591.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-592.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-593.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-594.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-595.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-596.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-597.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-598.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-599.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-6.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-60.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-600.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-601.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-602.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-603.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-604.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-605.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-606.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-607.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-608.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-609.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-61.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-610.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-611.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-612.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-613.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-614.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-615.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-616.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-617.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-618.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-619.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-62.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-620.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-621.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-63.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-64.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-65.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-66.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-67.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-68.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-69.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-7.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-70.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-71.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-72.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-73.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-74.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-75.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-76.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-77.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-78.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-79.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-8.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-80.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-81.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-82.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-83.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-84.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-85.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-86.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-87.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-88.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-89.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-9.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-90.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-91.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-92.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-93.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-94.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-95.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-96.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-97.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-98.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v03-99.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v04-1.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run4-OnPeak-R18b-v04-2.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-1.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-10.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-100.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-101.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-102.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-103.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-104.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-105.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-106.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-107.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-108.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-109.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-11.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-110.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-111.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-112.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-113.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-114.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-115.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-116.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-117.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-118.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-119.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-12.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-120.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-121.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-122.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-123.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-124.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-125.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-126.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-127.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-128.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-129.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-13.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-130.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-131.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-132.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-133.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-134.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-135.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-136.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-137.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-138.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-139.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-14.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-140.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-141.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-142.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-143.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-144.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-145.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-146.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-147.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-148.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-149.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-15.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-150.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-151.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-152.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-153.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-154.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-155.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-156.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-157.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-158.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-159.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-16.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-160.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-161.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-162.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-163.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-164.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-165.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-166.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-167.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-168.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-169.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-17.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-170.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-171.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-172.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-173.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-174.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-175.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-176.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-177.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-178.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-179.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-18.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-180.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-181.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-182.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-183.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-184.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-185.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-186.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-187.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-188.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-189.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-19.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-190.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-191.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-192.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-193.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-194.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-195.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-196.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-197.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-198.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-199.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-2.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-20.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-200.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-201.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-202.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-203.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-204.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-205.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-206.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-207.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-208.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-209.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-21.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-210.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-211.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-212.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-213.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-214.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-215.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-216.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-217.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-218.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-219.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-22.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-220.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-221.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-222.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-223.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-224.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-225.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-226.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-227.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-228.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-229.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-23.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-230.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-231.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-232.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-233.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-234.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-235.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-236.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-237.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-238.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-239.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-24.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-240.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-241.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-242.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-243.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-244.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-245.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-246.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-247.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-248.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-249.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-25.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-250.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-251.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-252.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-253.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-254.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-255.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-256.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-257.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-258.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-259.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-26.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-260.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-261.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-262.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-263.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-264.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-265.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-266.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-267.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-268.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-269.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-27.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-270.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-271.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-272.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-273.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-274.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-275.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-276.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-277.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-278.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-279.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-28.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-280.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-281.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-282.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-283.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-284.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-285.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-286.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-287.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-288.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-289.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-29.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-290.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-291.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-292.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-293.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-294.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-295.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-296.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-297.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-298.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-299.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-3.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-30.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-300.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-301.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-302.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-303.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-304.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-305.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-306.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-307.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-308.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-309.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-31.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-310.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-311.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-312.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-313.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-314.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-315.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-316.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-317.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-318.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-319.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-32.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-320.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-321.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-322.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-323.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-324.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-325.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-326.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-327.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-328.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-329.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-33.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-330.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-331.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-332.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-333.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-334.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-335.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-336.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-337.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-338.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-339.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-34.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-340.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-341.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-342.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-343.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-344.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-345.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-346.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-347.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-348.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-349.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-35.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-350.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-351.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-352.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-353.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-354.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-355.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-356.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-357.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-358.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-359.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-36.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-360.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-361.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-362.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-363.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-364.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-365.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-366.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-367.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-368.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-369.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-37.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-370.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-371.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-372.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-373.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-374.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-375.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-376.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-377.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-378.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-379.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-38.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-380.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-381.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-382.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-383.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-384.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-385.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-386.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-387.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-388.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-389.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-39.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-390.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-391.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-392.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-393.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-394.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-395.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-396.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-397.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-398.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-399.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-4.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-40.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-400.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-401.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-402.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-403.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-404.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-405.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-406.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-407.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-408.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-409.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-41.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-410.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-411.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-412.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-413.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-414.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-415.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-416.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-417.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-418.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-419.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-42.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-420.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-421.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-422.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-423.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-424.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-425.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-426.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-427.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-428.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-429.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-43.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-430.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-431.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-432.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-433.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-434.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-435.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-436.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-437.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-438.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-439.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-44.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-440.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-441.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-442.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-443.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-444.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-445.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-446.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-447.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-448.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-449.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-45.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-450.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-451.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-452.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-453.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-454.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-455.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-456.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-46.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-47.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-48.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-49.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-5.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-50.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-51.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-52.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-53.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-54.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-55.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-56.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-57.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-58.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-59.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-6.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-60.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-61.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-62.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-63.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-64.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-65.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-66.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-67.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-68.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-69.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-7.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-70.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-71.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-72.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-73.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-74.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-75.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-76.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-77.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-78.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-79.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-8.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-80.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-81.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-82.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-83.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-84.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-85.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-86.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-87.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-88.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-89.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-9.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-90.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-91.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-92.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-93.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-94.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-95.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-96.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-97.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-98.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v03-99.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-1.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-10.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-100.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-101.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-102.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-103.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-104.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-105.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-106.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-107.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-108.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-109.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-11.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-110.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-111.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-112.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-113.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-114.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-115.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-116.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-117.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-118.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-119.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-12.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-120.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-121.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-122.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-123.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-124.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-125.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-126.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-127.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-128.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-129.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-13.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-130.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-131.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-132.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-133.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-134.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-135.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-136.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-137.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-138.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-139.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-14.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-140.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-15.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-16.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-17.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-18.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-19.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-2.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-20.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-21.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-22.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-23.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-24.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-25.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-26.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-27.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-28.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-29.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-3.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-30.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-31.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-32.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-33.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-34.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-35.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-36.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-37.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-38.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-39.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-4.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-40.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-41.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-42.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-43.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-44.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-45.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-46.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-47.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-48.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-49.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-5.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-50.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-51.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-52.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-53.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-54.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-55.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-56.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-57.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-58.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-59.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-6.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-60.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-61.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-62.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-63.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-64.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-65.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-66.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-67.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-68.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-69.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-7.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-70.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-71.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-72.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-73.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-74.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-75.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-76.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-77.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-78.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-79.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-8.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-80.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-81.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-82.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-83.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-84.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-85.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-86.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-87.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-88.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-89.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-9.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-90.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-91.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-92.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-93.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-94.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-95.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-96.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-97.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-98.root/ntp1");
      chain->Add("$sl_awg_disc/dorazioa/dataSP8/BSemiExcl-Run5-OnPeak-R18b-v04-99.root/ntp1");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

MesPlots::~MesPlots()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MesPlots::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MesPlots::LoadTree(Long64_t entry)
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

void MesPlots::Init(TTree *tree)
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
   fChain->SetBranchAddress("indexbestRho0",&indexbestRho0);
   fChain->SetBranchAddress("chbestRho0",&chbestRho0);
   fChain->SetBranchAddress("mbestRho0",&mbestRho0);
   fChain->SetBranchAddress("barembestRho0",&barembestRho0);
   fChain->SetBranchAddress("mfitbestRho0",&mfitbestRho0);
   fChain->SetBranchAddress("baremfitbestRho0",&baremfitbestRho0);
   fChain->SetBranchAddress("mm2bestRho0",&mm2bestRho0);
   fChain->SetBranchAddress("q2bestRho0",&q2bestRho0);
   fChain->SetBranchAddress("q2fitbestRho0",&q2fitbestRho0);
   fChain->SetBranchAddress("mm2fitbestRho0",&mm2fitbestRho0);
   fChain->SetBranchAddress("nrecoRho0",&nrecoRho0);
   fChain->SetBranchAddress("mm2Rho0",mm2Rho0);
   fChain->SetBranchAddress("q2Rho0",q2Rho0);
   fChain->SetBranchAddress("pRho0",pRho0);
   fChain->SetBranchAddress("thRho0",thRho0);
   fChain->SetBranchAddress("phiRho0",phiRho0);
   fChain->SetBranchAddress("nominalmRho0",nominalmRho0);
   fChain->SetBranchAddress("baremRho0",baremRho0);
   fChain->SetBranchAddress("ElabRho0",ElabRho0);
   fChain->SetBranchAddress("chRho0",chRho0);
   fChain->SetBranchAddress("mm2fitRho0",mm2fitRho0);
   fChain->SetBranchAddress("q2fitRho0",q2fitRho0);
   fChain->SetBranchAddress("massfitRho0",massfitRho0);
   fChain->SetBranchAddress("pfitRho0",pfitRho0);
   fChain->SetBranchAddress("thetafitRho0",thetafitRho0);
   fChain->SetBranchAddress("phifitRho0",phifitRho0);
   fChain->SetBranchAddress("nominalmfitRho0",nominalmfitRho0);
   fChain->SetBranchAddress("baremfitRho0",baremfitRho0);
   fChain->SetBranchAddress("ElabfitRho0",ElabfitRho0);
   fChain->SetBranchAddress("chisqRho0",chisqRho0);
   fChain->SetBranchAddress("globchisqRho0",globchisqRho0);
   fChain->SetBranchAddress("probchisqRho0",probchisqRho0);
   fChain->SetBranchAddress("ndofRho0",ndofRho0);
   fChain->SetBranchAddress("pstarfitleptRho0",pstarfitleptRho0);
   fChain->SetBranchAddress("pfitleptRho0",pfitleptRho0);
   fChain->SetBranchAddress("thetafitleptRho0",thetafitleptRho0);
   fChain->SetBranchAddress("phifitleptRho0",phifitleptRho0);
   fChain->SetBranchAddress("pfitBRho0",pfitBRho0);
   fChain->SetBranchAddress("massBfitRho0",massBfitRho0);
   fChain->SetBranchAddress("thetafitBRho0",thetafitBRho0);
   fChain->SetBranchAddress("phifitBRho0",phifitBRho0);
   fChain->SetBranchAddress("ndauRho0",ndauRho0);
   fChain->SetBranchAddress("Estar1dauRho0",Estar1dauRho0);
   fChain->SetBranchAddress("Estar2dauRho0",Estar2dauRho0);
   fChain->SetBranchAddress("Elab1dauRho0",Elab1dauRho0);
   fChain->SetBranchAddress("Elab2dauRho0",Elab2dauRho0);
   fChain->SetBranchAddress("Estar1fitdauRho0",Estar1fitdauRho0);
   fChain->SetBranchAddress("Estar2fitdauRho0",Estar2fitdauRho0);
   fChain->SetBranchAddress("Elab1fitdauRho0",Elab1fitdauRho0);
   fChain->SetBranchAddress("Elab2fitdauRho0",Elab2fitdauRho0);
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
   Notify();
}

Bool_t MesPlots::Notify()
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
   b_indexbestRho0 = fChain->GetBranch("indexbestRho0");
   b_chbestRho0 = fChain->GetBranch("chbestRho0");
   b_mbestRho0 = fChain->GetBranch("mbestRho0");
   b_barembestRho0 = fChain->GetBranch("barembestRho0");
   b_mfitbestRho0 = fChain->GetBranch("mfitbestRho0");
   b_baremfitbestRho0 = fChain->GetBranch("baremfitbestRho0");
   b_mm2bestRho0 = fChain->GetBranch("mm2bestRho0");
   b_q2bestRho0 = fChain->GetBranch("q2bestRho0");
   b_q2fitbestRho0 = fChain->GetBranch("q2fitbestRho0");
   b_mm2fitbestRho0 = fChain->GetBranch("mm2fitbestRho0");
   b_nrecoRho0 = fChain->GetBranch("nrecoRho0");
   b_mm2Rho0 = fChain->GetBranch("mm2Rho0");
   b_q2Rho0 = fChain->GetBranch("q2Rho0");
   b_pRho0 = fChain->GetBranch("pRho0");
   b_thRho0 = fChain->GetBranch("thRho0");
   b_phiRho0 = fChain->GetBranch("phiRho0");
   b_nominalmRho0 = fChain->GetBranch("nominalmRho0");
   b_baremRho0 = fChain->GetBranch("baremRho0");
   b_ElabRho0 = fChain->GetBranch("ElabRho0");
   b_chRho0 = fChain->GetBranch("chRho0");
   b_mm2fitRho0 = fChain->GetBranch("mm2fitRho0");
   b_q2fitRho0 = fChain->GetBranch("q2fitRho0");
   b_massfitRho0 = fChain->GetBranch("massfitRho0");
   b_pfitRho0 = fChain->GetBranch("pfitRho0");
   b_thetafitRho0 = fChain->GetBranch("thetafitRho0");
   b_phifitRho0 = fChain->GetBranch("phifitRho0");
   b_nominalmfitRho0 = fChain->GetBranch("nominalmfitRho0");
   b_baremfitRho0 = fChain->GetBranch("baremfitRho0");
   b_ElabfitRho0 = fChain->GetBranch("ElabfitRho0");
   b_chisqRho0 = fChain->GetBranch("chisqRho0");
   b_globchisqRho0 = fChain->GetBranch("globchisqRho0");
   b_probchisqRho0 = fChain->GetBranch("probchisqRho0");
   b_ndofRho0 = fChain->GetBranch("ndofRho0");
   b_pstarfitleptRho0 = fChain->GetBranch("pstarfitleptRho0");
   b_pfitleptRho0 = fChain->GetBranch("pfitleptRho0");
   b_thetafitleptRho0 = fChain->GetBranch("thetafitleptRho0");
   b_phifitleptRho0 = fChain->GetBranch("phifitleptRho0");
   b_pfitBRho0 = fChain->GetBranch("pfitBRho0");
   b_massBfitRho0 = fChain->GetBranch("massBfitRho0");
   b_thetafitBRho0 = fChain->GetBranch("thetafitBRho0");
   b_phifitBRho0 = fChain->GetBranch("phifitBRho0");
   b_ndauRho0 = fChain->GetBranch("ndauRho0");
   b_Estar1dauRho0 = fChain->GetBranch("Estar1dauRho0");
   b_Estar2dauRho0 = fChain->GetBranch("Estar2dauRho0");
   b_Elab1dauRho0 = fChain->GetBranch("Elab1dauRho0");
   b_Elab2dauRho0 = fChain->GetBranch("Elab2dauRho0");
   b_Estar1fitdauRho0 = fChain->GetBranch("Estar1fitdauRho0");
   b_Estar2fitdauRho0 = fChain->GetBranch("Estar2fitdauRho0");
   b_Elab1fitdauRho0 = fChain->GetBranch("Elab1fitdauRho0");
   b_Elab2fitdauRho0 = fChain->GetBranch("Elab2fitdauRho0");
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

   return kTRUE;
}

void MesPlots::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MesPlots::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MesPlots_cxx
