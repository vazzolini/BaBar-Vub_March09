
#include "VirVubFitter/b2uClass.hh"

#include "TH2.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLegendEntry.h"
#include "TLegend.h"
#include <TRandom.h>
#include <TVector2.h>
#include <TLatex.h>
#include <TMinuit.h>
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooAddPdf.hh"
#include "RooFitCore/RooNLLVar.hh"
#include "RooFitCore/RooMinuit.hh"
#include "RooFitCore/RooFitResult.hh"
#include "RooFitCore/RooPlot.hh"
#include "RooFitCore/RooArgSet.hh"
#include "RooFitCore/RooBinning.rdl"
#include "RooFitModels/RooCBShape.hh"
#include "RooFitModels/RooArgusBG.hh"
#include "RooFitModels/RooGaussian.hh"

// ----------------------------------------------------------------------
b2uClass::b2uClass() {
}

// ----------------------------------------------------------------------
b2uClass::b2uClass(TTree *tree, TString filename, int Sys, int q2F, int comb, int un, int unfB, double hiunfB, int me, int mu, int usemxfit, bool iscm2,bool varfit,const vector<float>& wFermivec, int newbin):
  VirClass(tree, filename, Sys, q2F, comb, un, unfB, hiunfB, me, mu, iscm2, varfit, wFermivec, newbin)
{
  VlepmxYes = new RooRealVar("lepmxYes","lepmxYes",0,1);
  VmxYes =    new RooRealVar("mxYes"   ,"mxYes"   ,0,1);
  VlepYaSe2 = new RooRealVar("lepYaSe2","lepYaSe2",0,1);
  VlepYaSe3 = new RooRealVar("lepYaSe3","lepYaSe3",0,1);
  datamcvub_mx = new RooDataSet("","",RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*Vwe,*VlepmxYes,*VflavB,*VlepYaSe,*Vmultcat),"weight");
  datavubin_mx = new RooDataSet("","",RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*Vwe,*VlepmxYes,*VflavB,*VlepYaSe,*Vmultcat),"weight");
  datamcvub_r =  new RooDataSet("","",RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*Vwe,*VmxYes,*VflavB,*VlepYaSe,*Vmultcat),"weight");
  datavubin_r =  new RooDataSet("","",RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*Vwe,*VmxYes,*VflavB,*VlepYaSe,*Vmultcat),"weight");
  datamcvubtrue = new RooDataSet("","",RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*Vwe,*VlepYes,*VflavB,*VlepYaSe2,*Vmultcat),"weight");
  datamcvubfalse = new RooDataSet("","",RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*Vwe,*VlepYes,*VflavB,*VlepYaSe3,*Vmultcat),"weight");
  datamcvcbNC = new RooDataSet("","",RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*Vwe,*VlepYes,*VflavB,*VlepYaSe,*Vmultcat),"weight");
  datamcvubNC = new RooDataSet("","",RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*Vwe,*VlepYes,*VflavB,*VlepYaSe,*Vmultcat),"weight");

  fusemxfit = usemxfit;

}

// ----------------------------------------------------------------------
b2uClass::~b2uClass() {
  b2uresultDumping();
}

void b2uClass::b2uresultDumping() {

  double dGam, errdGam;
  double mymb, mya;
  char wfout[200];
  sprintf(wfout, "%s", fWeightFile.Data());
  if (!strcmp(wfout, "sysWd/belle3d_0.inc")) { mymb = 4.62; mya = 2.27;}
  if (!strcmp(wfout, "sysWd/belle3d_1.inc")) { mymb = 4.679992; mya = 2.868765;}
  if (!strcmp(wfout, "sysWd/belle3d_2.inc")) { mymb = 4.639684; mya = 3.888155;}
  if (!strcmp(wfout, "sysWd/belle3d_3.inc")) { mymb = 4.639684; mya = 1.446391;}
  if (!strcmp(wfout, "sysWd/belle3d_4.inc")) { mymb = 4.601767; mya = 3.204753;}
  if (!strcmp(wfout, "sysWd/belle3d_5.inc")) { mymb = 4.601767; mya = 1.232287;}
  if (!strcmp(wfout, "sysWd/belle3d_6.inc")) { mymb = 4.565861; mya = 2.448922;}
  if (!strcmp(wfout, "sysWd/belle3d_7.inc")) { mymb = 4.565861; mya = 1.202151;}
  if (!strcmp(wfout, "sysWd/belle3d_8.inc")) { mymb = 4.531675; mya = 1.533402;}
  if (!strcmp(wfout, "sysWd/vir_cleo3d_0.inc"))  { mymb = 4.735; mya = 1.600;}
  if (!strcmp(wfout, "sysWd/vir_cleo3d_1.inc"))  { mymb = 4.845; mya = 2.550;}
  if (!strcmp(wfout, "sysWd/vir_cleo3d_2.inc"))  { mymb = 4.480; mya = 0.574;}
  if (!strcmp(wfout, "sysWd/vir_cleo3d_3.inc"))  { mymb = 4.785; mya = 1.150;}
  if (!strcmp(wfout, "sysWd/vir_cleo3d_4.inc"))  { mymb = 4.785; mya = 3.590;}
  if (!strcmp(wfout, "sysWd/vir_cleo3d_5.inc"))  { mymb = 4.735; mya = 0.896;}
  if (!strcmp(wfout, "sysWd/vir_cleo3d_6.inc"))  { mymb = 4.690; mya = 0.680;}
  if (!strcmp(wfout, "sysWd/vir_cleo3d_7.inc"))  { mymb = 4.690; mya = 2.050;}
  if (!strcmp(wfout, "sysWd/vir_cleo3d_8.inc"))  { mymb = 4.580; mya = 0.550;}
  if (!strcmp(wfout, "sysWd/vir_cleo3d_9.inc"))  { mymb = 4.580; mya = 1.130;}
  if (!strcmp(wfout, "sysWd/vir_cleo3d_10.inc")) { mymb = 4.530; mya = 0.560;}

  double factor = (1+(N2/N1)-(N3/N1));
  double factorerr = (1/(N1*N1)) * sqrt(N1*N1*(N2err*N2err+N3err*N3err)+(N2*N2+N3*N3)*N1err*N1err);

  dGam = (.108*S*factor/(tot*fact*epsselbr*calcpstarfact));
  errdGam = dGam*sqrt((errS/S)*(errS/S)+(.0018/.108)*(.0018/.108));

  char name[200];
  sprintf(name,"%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),"b2uresults.dat");  
  ofstream outfile(name);
  outfile << "############################################################" << endl;
  outfile << "###        Output of b2uFit results for mma notebook     ###" << endl;
  outfile << "############################################################" << endl;
  outfile << "##" << endl;
  outfile << "## Weightfile used:" << endl;
  outfile << wfout << endl;
  outfile << "mb = " << mymb << endl;
  outfile << "a = " << mya << endl;
  outfile << "##" << endl;
  outfile << "## b->u yield - Nu" << endl;
  outfile << "Nu = " << S << " err " << errS << endl;
  outfile << "## N(SL)" << endl;
  outfile << "NSL = " << tot*fact << " err " << errtot*fact << endl;
  outfile << "## selection efficiencies " << endl;
  outfile << "## first, the BRBR/Nsl method " << endl;
  outfile << "factor = " << factor << " +/- " << factorerr << endl;
  outfile << "eps1 = " << epsumx << " err " << errepsumx << endl;
  outfile << "eps_mx = " << epschop << " err " << errepschop << endl;
  outfile << "## our eff in this method is the product of eps1 and eps_mx" << endl;
  outfile << "epsu = " <<  epsselbr << " err " << errepsselbr << endl;
  outfile << "## ratio of the lepton efficiencies, epssl" << endl;
  outfile << "epssl = " << (1/calcpstarfact) << " err " << errcalcpstarfact << endl;
  outfile << "## ratio of the reco efficiencies, epsbr - always 1" << endl;
  outfile << "epsbr = 1.00 err 0.03" << endl;
  outfile << "## Now, the Nbreco way:" << endl;
  outfile << "##epssel = " << epssel << " err " << errepssel << endl;
  outfile << "##epsbrbias = 0.969 err 0.257" << endl;
  outfile << "##Nbreco = 323657 err 1963 err 5076" << endl;
  outfile << "############################################################" << endl;
  outfile << "############################################################" << endl;
  outfile << "## conversion of the mX cut to LLR variable, c" << endl;
  outfile << "mXcut = " << MXBIN << endl;
  outfile << "c = " << (MXBIN*MXBIN/(5.279*5.279)) << endl;
  outfile << "## calc of dGam for the ratio of BR method" << endl;
  outfile << "dGam(" << (MXBIN*MXBIN/(5.279*5.279)) << ") = " << dGam << " err " << errdGam << endl;
  outfile << "############################################################" << endl;
  outfile << "## calc of dGam for the Nbreco method" << endl;
  outfile << "## dGam(" << (MXBIN*MXBIN/(5.279*5.279)) << ") = " << (S*0.969/(epssel*323657)) << endl;
  outfile.close();

}

void b2uClass::Loop(int isdata, int icat, int isMC, int nres, int comb)

  //---- flags legenda:
  // isdata  0 = MC to get the shapes and the efficiencies, 1 = data 
  // icat    0 = vub component, 1 = vcb and other component, 
  // nevents = number of events
  // isMC    1 = fit on the MC, 0 = fit on data

{

  fHistFile->cd();

  int flavB, type; 

  int countFB = 0;
  if (fChain == 0) return;
  
  Int_t nentries = Int_t(fChain->GetEntries());
  int nevents = 1111111111;
  if(isdata == 0 && icat == 0) {
    nevents = nVUB;
  } else if (isdata == 0 && icat == 1) {
    nevents = nVCB;
  } else {
    if (isdata == 1 && isMC == 1) {
      if(icat == 1){
	nevents = nVCBDATA;
      } else if (icat == 0) {
	nevents = nVUBDATA;
      }
    } else if (isMC ==0){
      nevents = nDATA;
    }
  }

  cout<<"nevents:"<<nevents<<endl;

  if( nentries > nevents) nentries = nevents;

  cout <<  nentries << " Entries"<<endl;
   
  Int_t nbytes = 0, nb = 0;

  
  f1all = new TF1("gaussall","gaus",-5*SMEARALLSIGMA,5*SMEARALLSIGMA);
  f1bkg = new TF1("gaussbck","gaus",-5*SMEARBKGSIGMA,5*SMEARBKGSIGMA);
  TRandom rand(0);
  double w, wfermi, wnre;
  for (Int_t jentry=0; jentry<nentries;jentry++) {     
    Int_t ientry = LoadTree(jentry); //---- in case of a TChain, ientry is the entry number in the current file
    nb = fChain->GetEntry(jentry);  
    nbytes += nb;

    if (!(TMath::Abs(brecocharge) == BTYPE || BTYPE == 2)) continue;

    //Calculate ev-by-ev reweightings
    w = wfermi = wnre = 1;           //Weight initialization
    

    if( !isdata ){
      //Vcb and oth weights
      w *= getTrackingWeight();                        //trk weighting
      w *= getNeutralWeight();                         //neu weighting
      w *= getBrecoWeight(intpur);                     //breco weighting  
      if(icat != 0) {
	w *= getBsysweight(Gvxbtyp,vub);                                //Bdec weighting
	w *= getDsysweight(GfDpi,GfDk,GfDks,GfDpiz,GfDlep,dImode,vub);  //Ddec weighting
      }
      wfermi = w;
      //Star vub non-resonantspecific part

      if(vub)  {

	//CB IMPLEMENT NEW REWEIGHTING FROM D. FORTIN
	//NEW REWEIGHTING IS APPLIED ALSO ON WIDE RESONANT EVENTS 
	//ONLY rho, pi, eta, eta', omega ARE CONSIDERED AS RESONANT...
	//rescale BRs according to the table presented in 
	//http://www.slac.stanford.edu/BFROOT/www/Physics/Analysis/AWG/InclusiveSL/common/dominique_hybrid.html
	if(fprlRew==10){
	  if(Gvxbtyp == 11 ){            // B0bar ->pi+ l- nubar + CC
	    //CB am I applying this also for generic BBbar events? 
	    //wfermi *= 1.33/1.8;
	    wfermi *= fweights[0];
	  }else if(Gvxbtyp == -11 ){     // B- -> pi0 l- nubar + CC
	    //wfermi *= 0.72/0.9;
	    wfermi *= fweights[1];
	  }else if(Gvxbtyp == -12){      // B- ->eta0 l- nubar + CC
	    //wfermi *= 0.84/0.3;
	    wfermi *= fweights[2];
	  }else if(Gvxbtyp == 13) {      // B0bar -> rho+ l- nubar + CC
	    //wfermi *= 2.69/2.6;
	    wfermi *= fweights[3];
	  }else if(Gvxbtyp == -13) {     // B- -> rho0 l- nubar + CC
	    //wfermi *= 1.45/1.3;
	    wfermi *= fweights[4];
	  }else if(Gvxbtyp == -14) {     // B- -> omega0 l- nubar + CC
	    //wfermi *= 1.34/1.3;
	    wfermi *= fweights[5];
	  }else if(Gvxbtyp == -15) {     // B- -> eta_prime l- nubar + CC
	    //wfermi *= 0.84/0.6;
	    wfermi *= fweights[6];
	  }else if(Gvxbtyp == 7)   {     // B0bar -> Xu+ l- nubar + CC
	    if(nres == 0) continue;
	    //
	    //	    wfermi *= 17.48/13.65; actually I am taking non-res from a different datafile...
	    //                             so I have to care about nonres/res normalization 
            //compute normalization factor between resonant and old hybrid...
            //RIC says the ratio between resonant (narrow+wide) and non-resonant events in root files is 0.91
            //I assume is the same for B+- and B0 (?)
            //so I have to scale 0.91 with (BR_res^NEW / BR_res^OLD) / (BR_res^NEW / BR_nres^NEW)=BR_nres^NEW /BR_res^OLD=1.71/0.735

	    // this is the old mysterious value
	    //wfermi *= 2.33;
	    // this is the value for CLEO SF
	    //wfermi *= 2.16;
	    // this is the value for BELLE SF
	    //wfermi *= 2.80;
	    wfermi *= fweights[7];

	    int ys = b2unewrHistmx(mxhadgen)+rHistq2(q2Gen)*8 + rHistel(ecmsgen)*8*8;
	    //
	    //	    cout << " B0 ys mx q2 el bins:   " << ys << "  " 
	    //		 << b2unewrHistmx(mxhadgen) << "  "
	    //		 << rHistq2(q2Gen) << "  "
	    //		 << rHistel(ecmsgen) << endl;
	    //
	    wfermi *=newMatrixW[ys];
	    if(FERMIAPP) {
	      wfermi *=FermiWeight(fkplus,DELTAMB,DELTAA);

	    }
	  }else if(Gvxbtyp == -7) {      // B- -> Xu0 l- nubar + CC
	    if(nres == 0) continue;
	    //
	    //	    wfermi *= 17.90/13.70; actually I am taking non-res from a different datafile...
	    //                             so I have to care about nonres/res normalization 
	    //so I have to scale 0.91 with (BR_res^NEW / BR_res^OLD) / (BR_res^NEW / BR_nres^NEW)=BR_nres^NEW /BR_res^OLD=1.71/0.735
	    //wfermi *= 2.33;

	    // this is the old mysterious value
	    //wfermi *= 2.33;
            // this is the value for CLEO SF
	    //wfermi *= 2.20;
	    // this is the value for BELLE SF
	    //wfermi *= 2.90;
	    wfermi *= fweights[8];

	    int ys = (b2unewrHistmx(mxhadgen)+rHistq2(q2Gen)*8 + rHistel(ecmsgen)*8*8)+512;
	    //
	    wfermi *= newMatrixW[ys];
	    if(FERMIAPP) {
	      wfermi *=FermiWeight(fkplus,DELTAMB,DELTAA);
	    }
	  }else{
	    //reject all other vubish decays
	    continue;
	  }
	}else{
	  //old implementation
	  if(TMath::Abs(Gvxbtyp)==7 ){
	    if(nres == 0) continue;
	    if(FERMIAPP) {
	      wnre *= FermiWeight(fkplus,DELTAMB,DELTAA);
	      wfermi *= wnre;
	    }
            if(fprlRew==1) {
              wfermi *= getTrueMxWeight(mxhadgen,Gvxbtyp);
              //            cout<<"Weight :: "<<getTrueMxWeight(mxhadgen,Gvxbtyp)<<" "<<mxhadgen<<" "<<Gvxbtyp<<endl;
            } else if(fprlRew==3) {
              int ys = rHistmx(mxhadgen)+ rHistq2(q2Gen)*14 + rHistel(ecmsgen)*14*8;
              //int ys = rHistmx(mxhadgen);
              wfermi *= 1.7755*MatrixW[ys];
              if(MULFAC != 0) {
                wfermi *= 1/MULFAC;
                //      cout<<"MatrixW[ys]:: "<<MatrixW[ys]<<" "<<ys<<"  "<<MULFAC<<"            nel LOOP"<<endl;
              }
            }else {
                           cout<<"not reweighting Vub events"<<endl;
            }
	  }
 	}
	if(TMath::Abs(Gvxbtyp) == 7){
	  //	  cout  << mxhadgen << "  " <<q2Gen<< "  "<<ecmsgen<< "  "<<wfermi<<endl;
	  ((TH1D*)gDirectory->Get("inclwei"))->Fill(wfermi, 1.);
	} else {
	  ((TH1D*)gDirectory->Get("exclwei"))->Fill(wfermi, 1.);
	}
      }
    }

    //---- tag the event with the lepton type

    if (fusemxfit==0) {
      if (FITQ2 ==0){
	if(TMath::IsNaN(mxhad)) {
	  cout << "MEZZEGA: NAN in MXHAD!!" << endl;
	  mxhad = -999.;
	}   
      }else{
	if(TMath::IsNaN(q2fit)) {
	  cout << "MEZZEGA: NAN in Q2FIT!!" << endl;
	  q2fit = -999.;
	}   
      }
    } else {
      if (FITQ2 ==0){
	if(TMath::IsNaN(mxhadfit)) {
	  cout << "MEZZEGA: NAN in MXHADFIT!!" << endl;
	  mxhadfit = -999.;
	}   
      }else{
	if(TMath::IsNaN(q2fit)) {
	  cout << "MEZZEGA: NAN in Q2FIT!!" << endl;
	  q2fit = -999.;
	}   
      }
    }

    bool rCh = (TMath::Abs(brecocharge) == BTYPE || BTYPE == 2);
    //Semileptonic events selection
    bool isLp = ((nle > 0) && rCh);  
    if(LEPTTYPE == 0)
      isLp = (nel > 0);
    if(LEPTTYPE == 1) 
      isLp = (nmu > 0);
    //Signal Events Selection    
    bool isLpSig = ((nle == 1) && rCh);  
    if(LEPTTYPE == 0)
      isLpSig = (nel == 1);
    if(LEPTTYPE == 1) 
      isLpSig = (nmu == 1);
    

    //---- charge correlation
    int flav =  lcharge + brecoflav;
    bool lPmxYes, mxYes;
    bool lPYes = (pcms > LEPTONPCUT && isLp && !(TMath::Abs(brecocharge)!=0 && (flav)!=0) );
    if (fusemxfit==0) {
      lPmxYes = (pcms > LEPTONPCUT && isLp && !(TMath::Abs(brecocharge)!=0 && (flav)!=0) && mxhad>0. && mxhad<MXBIN);
    } else {
      lPmxYes = (pcms > LEPTONPCUT && isLp && !(TMath::Abs(brecocharge)!=0 && (flav)!=0) && mxhadfit>0. && mxhadfit<MXBIN);
    }
    if (fusemxfit==0) {
      mxYes = (mxhad<MXBIN);
    } else {
      mxYes = (mxhadfit<MXBIN);
    }
    bool lPYesSig = (pcms > LEPTONPCUT && isLpSig && !(TMath::Abs(brecocharge)!=0 && (flav)!=0) );

    //---- flavor category 
    //   3 = charged B,  4 = neutral B OS,  5 = neutral B SS
    
    flavB = 5;
    if(TMath::Abs(brecocharge)==1)flavB = 3;
    if(TMath::Abs(brecocharge)==0 && flav==0)flavB = 4;
    
    int type = 0;
    int area = 0;

    //dataset component type
    //at this point: isdata is = 0,icat !=0  
    if(isMC) {
      if(vcb && icat!=0) type = 1;              //vcb
      if(vub && icat ==0) type = 2;             //vub
      if((vcb + vub) == 0 && icat!=0) type = 3; //other 
      if(Sun && type == 2){
	if(comb){	
	  if(mxhadgen < MXBIN && q2Gen > Q2BIN){ 
	    area = 1;                             //vubin  2D
	  }else{ 
	    area = 2;                             //vubout 2D
	  }
	} else {
	  if(mxhadgen < MXBIN){ 
	    area = 1;                             //vubin  1D
	  }else{ 
	    area = 2;                             //vubout 1D
	  }
	}
      }
    }
 
    int ch = TMath::Abs(xcharge + brecocharge);  // total charge

    // int ibch = 1;
    // AllCut  # survivors lepYes to mass nu cut  

    bool ksele = (nkp + nks) > 0;   // fit on the depleted sample?

    //Multiplicity category (X l nu side)
    int mult=0;
    if(nchg == 1 && nneu > 0) mult = 1;
    if((nchg == 2||nchg==3) && nneu == 0) mult = 2;
    if((nchg == 2||nchg==3) && nneu > 0) mult = 3;
    if(nchg > 3 && nneu == 0) mult = 4;
    if(nchg > 3  && nneu > 0) mult = 5;

    //Generated multiplicity category (X l nu side)
    int multgen=0;
    if(chgdaugen == 1 && neudaugen > 0) multgen = 1;
    if((chgdaugen == 2||chgdaugen==3) && neudaugen == 0) multgen = 2;
    if((chgdaugen == 2||chgdaugen==3) && neudaugen > 0) multgen = 3;
    if(chgdaugen > 3 && neudaugen == 0) multgen = 4;
    if(chgdaugen > 3  && neudaugen > 0) multgen = 5;

    //bool WdeltaCut = (mm2pr>PRMM2CUT && brecocharge == 0); //last cut based on mm2pr
    bool WdeltaCut = ((mm1pr>PRMM1CUT) || (mm3pr>PRMM3CUT) || (mm2pr>PRMM2CUT));   // new PR cuts EJH
    bool AllCut, AllCut2, AllCut3;
    if(MU){
      if (fusemxfit==0) {
	AllCut= (lPYesSig && q2fit>Q2CUT && mm2 < MNUSQHIGH && mm2 > MNUSQLOW &&  ch < CHHIGH && ch > CHLOW &&  ksele == DEPL && mxhad>0. && mxhad < MXCUT && !(WdeltaCut) && mult==MU && pmiss>PMISSCUTLO && (cos(tmiss))<COSTMISSCUTHI && (cos(tmiss))>COSTMISSCUTLO && pcmstrklo>PCMSTRKLOCUT);
	AllCut2= (lPYesSig && q2fit>Q2CUT && mm2 < MNUSQHIGH && mm2 > MNUSQLOW &&  ch < CHHIGH && ch > CHLOW &&  ksele == DEPL && mxhad>0. && mxhad < MXCUT && !(WdeltaCut) && mult==MU && pmiss>PMISSCUTLO && (cos(tmiss))<COSTMISSCUTHI && (cos(tmiss))>COSTMISSCUTLO && pcmstrklo>PCMSTRKLOCUT && mxhadgen<MXBIN);
	AllCut3= (lPYesSig && q2fit>Q2CUT && mm2 < MNUSQHIGH && mm2 > MNUSQLOW &&  ch < CHHIGH && ch > CHLOW &&  ksele == DEPL && mxhad>0. && mxhad < MXCUT && !(WdeltaCut) && mult==MU && pmiss>PMISSCUTLO && (cos(tmiss))<COSTMISSCUTHI && (cos(tmiss))>COSTMISSCUTLO && pcmstrklo>PCMSTRKLOCUT && mxhadgen>MXBIN);
      } else {
	AllCut= (lPYesSig && q2fit>Q2CUT && mm2 < MNUSQHIGH && mm2 > MNUSQLOW &&  ch < CHHIGH && ch > CHLOW &&  ksele == DEPL && mxhadfit>0. && mxhadfit < MXCUT && !(WdeltaCut) && mult==MU && pmiss>PMISSCUTLO && (cos(tmiss))<COSTMISSCUTHI && (cos(tmiss))>COSTMISSCUTLO && pcmstrklo>PCMSTRKLOCUT);
	AllCut2= (lPYesSig && q2fit>Q2CUT && mm2 < MNUSQHIGH && mm2 > MNUSQLOW &&  ch < CHHIGH && ch > CHLOW &&  ksele == DEPL && mxhadfit>0. && mxhadfit < MXCUT && !(WdeltaCut) && mult==MU && pmiss>PMISSCUTLO && (cos(tmiss))<COSTMISSCUTHI && (cos(tmiss))>COSTMISSCUTLO && pcmstrklo>PCMSTRKLOCUT && mxhadgen<MXBIN);
	AllCut3= (lPYesSig && q2fit>Q2CUT && mm2 < MNUSQHIGH && mm2 > MNUSQLOW &&  ch < CHHIGH && ch > CHLOW &&  ksele == DEPL && mxhadfit>0. && mxhadfit < MXCUT && !(WdeltaCut) && mult==MU && pmiss>PMISSCUTLO && (cos(tmiss))<COSTMISSCUTHI && (cos(tmiss))>COSTMISSCUTLO && pcmstrklo>PCMSTRKLOCUT && mxhadgen>MXBIN);
      }
    } else{
      if (fusemxfit==0) {
	AllCut= (lPYesSig && q2fit>Q2CUT && mm2 < MNUSQHIGH && mm2 > MNUSQLOW &&  ch < CHHIGH && ch > CHLOW &&  ksele == DEPL && mxhad>0.&& mxhad < MXCUT && !(WdeltaCut) && pmiss>PMISSCUTLO && (cos(tmiss))<COSTMISSCUTHI && (cos(tmiss))>COSTMISSCUTLO && pcmstrklo>PCMSTRKLOCUT);
	AllCut2= (lPYesSig && q2fit>Q2CUT && mm2 < MNUSQHIGH && mm2 > MNUSQLOW &&  ch < CHHIGH && ch > CHLOW &&  ksele == DEPL && mxhad>0.&& mxhad < MXCUT && !(WdeltaCut) && pmiss>PMISSCUTLO && (cos(tmiss))<COSTMISSCUTHI && (cos(tmiss))>COSTMISSCUTLO && pcmstrklo>PCMSTRKLOCUT && mxhadgen<MXBIN);
	AllCut3= (lPYesSig && q2fit>Q2CUT && mm2 < MNUSQHIGH && mm2 > MNUSQLOW &&  ch < CHHIGH && ch > CHLOW &&  ksele == DEPL && mxhad>0.&& mxhad < MXCUT && !(WdeltaCut) && pmiss>PMISSCUTLO && (cos(tmiss))<COSTMISSCUTHI && (cos(tmiss))>COSTMISSCUTLO && pcmstrklo>PCMSTRKLOCUT && mxhadgen>MXBIN);
      } else {
	AllCut= (lPYesSig && q2fit>Q2CUT && mm2 < MNUSQHIGH && mm2 > MNUSQLOW &&  ch < CHHIGH && ch > CHLOW &&  ksele == DEPL && mxhadfit>0.&& mxhadfit < MXCUT && !(WdeltaCut) && pmiss>PMISSCUTLO && (cos(tmiss))<COSTMISSCUTHI && (cos(tmiss))>COSTMISSCUTLO && pcmstrklo>PCMSTRKLOCUT);
	AllCut2= (lPYesSig && q2fit>Q2CUT && mm2 < MNUSQHIGH && mm2 > MNUSQLOW &&  ch < CHHIGH && ch > CHLOW &&  ksele == DEPL && mxhadfit>0.&& mxhadfit < MXCUT && !(WdeltaCut) && pmiss>PMISSCUTLO && (cos(tmiss))<COSTMISSCUTHI && (cos(tmiss))>COSTMISSCUTLO && pcmstrklo>PCMSTRKLOCUT && mxhadgen<MXBIN);
	AllCut3= (lPYesSig && q2fit>Q2CUT && mm2 < MNUSQHIGH && mm2 > MNUSQLOW &&  ch < CHHIGH && ch > CHLOW &&  ksele == DEPL && mxhadfit>0.&& mxhadfit < MXCUT && !(WdeltaCut) && pmiss>PMISSCUTLO && (cos(tmiss))<COSTMISSCUTHI && (cos(tmiss))>COSTMISSCUTLO && pcmstrklo>PCMSTRKLOCUT && mxhadgen>MXBIN);
      }
    }

    if(mes>5.2) {      
      Vmes->setVal(mes);
      if (fusemxfit==0) {
	Vmx->setVal(mxhad);
      } else {
	Vmx->setVal(mxhadfit);
      }
      Vq2->setVal(q2fit);
      if (FITQ2 ==0){
	if (fusemxfit==0) {
	  Vchop->setVal(mxhad);
	} else {
	  Vchop->setVal(mxhadfit);
	}
      }else{
	Vchop->setVal(q2fit);
      }
      VlepYes->setVal(lPYes);  
      VlepmxYes->setVal(lPmxYes);  
      VflavB->setVal(flavB);
      VlepYaSe->setVal(AllCut);
      VlepYaSe2->setVal(AllCut2);
      VlepYaSe3->setVal(AllCut3);
      Vmultcat->setVal(mult);
    }
    VmxYes->setVal(mxYes);
    if(mes>5.2 && isdata!=1 && mxYes) {
      if(!Sun){
	if(type == 2)  datamcvub_r->add(RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*VmxYes,*VflavB, *VlepYaSe , *Vmultcat),wfermi);
      }else{
	if(area == 1) datavubin_r->add(RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*VmxYes,*VflavB, *VlepYaSe, *Vmultcat),wfermi);
      }

    }

    if(mes>5.2) {      
      if(isdata != 1)  {
	  if(type == 2)  datamcvubNC->add(RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*VlepYes,*VflavB, *VlepYaSe, *Vmultcat),wfermi);
	  if(type == 1)  datamcvcbNC->add(RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*VlepYes,*VflavB, *VlepYaSe, *Vmultcat),w);
      }
    }

    if(mes>5.2 && lPYes) {      
      if(isdata == 1)  {
	datadata->add(RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*VlepYes,*VflavB, *VlepYaSe, *Vmultcat),1.);
      } else {
	if(!Sun){
	  if(type == 2)  datamcvub->add(RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*VlepYes,*VflavB, *VlepYaSe , *Vmultcat),wfermi);
	  if(type == 2)  datamcvubtrue->add(RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*VlepYes,*VflavB, *VlepYaSe2 , *Vmultcat),wfermi);
	  if(type == 2)  datamcvubfalse->add(RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*VlepYes,*VflavB, *VlepYaSe3 , *Vmultcat),wfermi);
	  if(type == 2&&lPmxYes)  datamcvub_mx->add(RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*VlepmxYes,*VflavB, *VlepYaSe , *Vmultcat),wfermi);
	  //    if(type == 1)  datamcvcb->add(RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*VlepYes,*VflavB, *VlepYaSe, *Vmultcat),w);
	  if(type == 3)  datamcoth->add(RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*VlepYes,*VflavB, *VlepYaSe, *Vmultcat),w);
	}else{
	  ///FILL RESONANT SIGNAL ONLY IF RES OPTION IS SELECTED
	  if(type == 1 || type == 3) datavcboth->add(RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*VlepYes,*VflavB, *VlepYaSe, *Vmultcat),w);
	  if(area == 1) datavubin->add(RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*VlepYes,*VflavB, *VlepYaSe, *Vmultcat),wfermi);
	  if(area == 1&&lPmxYes) datavubin_mx->add(RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*VlepmxYes,*VflavB, *VlepYaSe, *Vmultcat),wfermi);
	  if(area == 2) datavubout->add(RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*VlepYes,*VflavB, *VlepYaSe, *Vmultcat),wfermi);	
	  //   if(nres == 0 && TMath::Abs(Gvxbtyp)==7 ){
	  // if(area == 1) datavubin->add(RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*VlepYes,*VflavB, *VlepYaSe, *Vmultcat),wfermi);
	  // if(area == 2) datavubout->add(RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*VlepYes,*VflavB, *VlepYaSe, *Vmultcat),wfermi);	
	  //}
	}
	if(type == 1)   datamcvcb->add(RooArgSet(*Vmes,*Vchop,*Vmx,*Vq2,*VlepYes,*VflavB, *VlepYaSe, *Vmultcat),w);
      }
    }
    if(UNFBINNING&&type==2){
      Vallmes->setVal(mes);
      if (fusemxfit==0) {
	Vchop->setVal(mxhad);
      } else {
	Vchop->setVal(mxhadfit);
      }
      VlepYaSe->setVal(AllCut);
      VlepYaSe2->setVal(AllCut2);
      VlepYaSe3->setVal(AllCut3);
      Vmxgenwoph->setVal(mxhadgenwoph);
      Vbrecoqualy->setVal(brecoqual);
      Vmultcat->setVal(mult);
      Vmultcatgen->setVal(multgen);
      Vvxbtyp->setVal(Gvxbtyp);
      unfmcvub->add(RooArgSet(*Vallmes,*Vchop,*VlepYaSe,*Vmxgenwoph,*Vbrecoqualy, *Vmultcat, *Vmultcatgen,*Vvxbtyp),wfermi);
    }
  }
}

void b2uClass::b2ueff(int su){

  vubmcmx = 0; errvubmcmx = 0; vubmcr = 0; errvubmcr = 0; 
  epsumx = 0; errepsumx = 0; epssel = 0; errepssel = 0; epsselbr = 0; errepsselbr = 0;

  double Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar;
  char simply[100], line[200];
  double chid = 0.181;
  //TVector2 signalsig;
  TVector2 signal;
  TCanvas *c2  = new TCanvas("c2"," ",200,10,1200,1000); 
  c2->Clear();
  c2->Divide(3, 3);

  c2->cd(1); 
  //----------vub lepton mc bch----------- 
  sprintf(simply,"(flavB==3) && lepmxYes");
  if(su == 1){
    signal = sighistounb(datavubin_mx, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, mesvubMC[0], mesvubMC[1], mesvubMC[2], 5.,-60., 6000, 1000, simply, USECB);
    cout<<"entrati:::"<<endl;
  } else{
    signal = sighistounb(datamcvub_mx, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, mesvubMC[0], mesvubMC[1], mesvubMC[2], 5.,-60., 6000, 1000, simply, USECB);
  }
  cout<<"Fitting b2ueff events on vub mc3!!! :: "<<  Nresmean<<" "<< Nressigma<<" "<< Nresalpha<<" "<< Nresn<<" "<< Nresargpar <<endl;
  double  vub3 = signal.X();                    double  errvub3 = signal.Y();       
  if(BINNED){}else{xframe->Draw();};

  c2->cd(2); 
  //----------vub lepton mc bos----------- 
  sprintf(simply,"(flavB==4) && lepmxYes");
  if(su == 1){
    signal = sighistounb(datavubin_mx, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, mesvubMC[0], mesvubMC[1], mesvubMC[2], 5.,-60., 6000, 1000, simply, USECB);
  } else{
    signal = sighistounb(datamcvub_mx, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, mesvubMC[0], mesvubMC[1], mesvubMC[2], 5.,-60., 6000, 1000, simply, USECB);
  }
  cout<<"Fitting b2ueff events on vub mc4!!! :: "<<  Nresmean<<" "<< Nressigma<<" "<< Nresalpha<<" "<< Nresn<<" "<< Nresargpar <<endl;
  double  vub4 = signal.X();                    double errvub4 = signal.Y();       
  if(BINNED){}else{xframe->Draw();};

  c2->cd(3);
  //----------vub lepton mc bss----------- 
  sprintf(simply,"(flavB==5) && lepmxYes");
  if(su == 1){
    signal = sighistounb(datavubin_mx, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, mesvubMC[0], mesvubMC[1], mesvubMC[2], 5.,-60., 6000, 1000, simply, USECB);
  } else{
    signal = sighistounb(datamcvub_mx, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, mesvubMC[0], mesvubMC[1], mesvubMC[2], 5.,-60., 6000, 1000, simply, USECB);
  }
  cout<<"Fitting b2ueff events on vub mc5!!! :: "<<  Nresmean<<" "<< Nressigma<<" "<< Nresalpha<<" "<< Nresn<<" "<< Nresargpar <<endl;
  double  vub5 = signal.X();                    double  errvub5 = signal.Y();       
  if(BINNED){}else{xframe->Draw();};
  
  if(MIXCORR==0){
    vubmcmx = vub3 + vub4 +vub5;
    errvubmcmx = sqrt(errvub3 * errvub3 + errvub4 * errvub4 + errvub5 * errvub5);
  }else{      
    vubmcmx = vub3 + ((1-chid)/(1-2*chid)) * vub4 - (chid/(1-2*chid)) * vub5;
    errvubmcmx = sqrt(errvub3 * errvub3 +
		    ((1-chid)/(1-2*chid)) * ((1-chid)/(1-2*chid)) * errvub4 * errvub4 +
		    (chid/(1-2*chid)) * (chid/(1-2*chid)) * errvub5 * errvub5);
    cout<<"control:::::"<<"vub3  "<<vub3<<"  errvub3  "<<errvub3<<"vub4  "<<vub4<<"  errvub4  "<<errvub4<<"vub5  "<<vub5<<"  errvub5  "<<errvub5<<endl;    
    cout<<"control:::::"<<"           vubmcmx  "<<vubmcmx<<"  errvubmcmx  "<<errvubmcmx<<endl;    
  }


  /*

  c2->cd(4); 
  //----------vub lepton mc bch----------- 
  sprintf(simply,"(flavB==3) && mxYes");
  if(su == 1){
    signal = sighistounb(datavubin_r, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, mesvubMC[0], mesvubMC[1], mesvubMC[2], 5.,-60., 6000, 1000, simply, USECB);
    cout<<"entrati:::"<<endl;
  } else{
    signal = sighistounb(datamcvub_r, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, mesvubMC[0], mesvubMC[1], mesvubMC[2], 5.,-60., 6000, 1000, simply, USECB);
  }
  cout<<"Fitting semileptonic events on vub mc3!!! :: "<<  Nresmean<<" "<< Nressigma<<" "<< Nresalpha<<" "<< Nresn<<" "<< Nresargpar <<endl;
  double  vubr3 = signal.X();                    double  errvubr3 = signal.Y();       
  if(BINNED){}else{xframe->Draw();};

  c2->cd(5); 
  //----------vub lepton mc bos----------- 
  sprintf(simply,"(flavB==4) && mxYes");
  if(su == 1){
    signal = sighistounb(datavubin_r, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, mesvubMC[0], mesvubMC[1], mesvubMC[2], 5.,-60., 6000, 1000, simply, USECB);
  } else{
    signal = sighistounb(datamcvub_r, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, mesvubMC[0], mesvubMC[1], mesvubMC[2], 5.,-60., 6000, 1000, simply, USECB);
  }
  cout<<"Fitting semileptonic events on vub mc4!!! :: "<<  Nresmean<<" "<< Nressigma<<" "<< Nresalpha<<" "<< Nresn<<" "<< Nresargpar <<endl;
  double  vubr4 = signal.X();                    double errvubr4 = signal.Y();       
  if(BINNED){}else{xframe->Draw();};

  c2->cd(6);
  //----------vub lepton mc bss----------- 
  sprintf(simply,"(flavB==5) && mxYes");
  if(su == 1){
    signal = sighistounb(datavubin_r, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, mesvubMC[0], mesvubMC[1], mesvubMC[2], 5.,-60., 6000, 1000, simply, USECB);
  } else{
    signal = sighistounb(datamcvub_r, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, mesvubMC[0], mesvubMC[1], mesvubMC[2], 5.,-60., 6000, 1000, simply, USECB);
  }
  cout<<"Fitting semileptonic events on vub mc5!!! :: "<<  Nresmean<<" "<< Nressigma<<" "<< Nresalpha<<" "<< Nresn<<" "<< Nresargpar <<endl;
  double  vubr5 = signal.X();                    double  errvubr5 = signal.Y();       
  if(BINNED){}else{xframe->Draw();};
  
  if(MIXCORR==0){
    vubmcr = vubr3 + vubr4 +vubr5;
    errvubmcr = sqrt(errvubr3 * errvubr3 + errvubr4 * errvubr4 + errvubr5 * errvubr5);
  }else{      
    vubmcr = vubr3 + ((1-chid)/(1-2*chid)) * vubr4 - (chid/(1-2*chid)) * vubr5;
    errvubmcr = sqrt(errvubr3 * errvubr3 +
		    ((1-chid)/(1-2*chid)) * ((1-chid)/(1-2*chid)) * errvubr4 * errvubr4 +
		    (chid/(1-2*chid)) * (chid/(1-2*chid)) * errvubr5 * errvubr5);
    cout<<"control:::::"<<"vubr3  "<<vubr3<<"  errvubr3  "<<errvubr3<<"vubr4  "<<vubr4<<"  errvubr4  "<<errvubr4<<"vubr5  "<<vubr5<<"  errvubr5  "<<errvubr5<<endl;    
    cout<<"control:::::"<<"           vubmcr  "<<vubmcr<<"  errvubmcr  "<<errvubmcr<<endl;    
  }


  */
  vubmcr = 1; errvubmcr=1;

  epsumx = vubmcallforeff/vubmcmx;
  errepsumx = sqrt(epsumx*(1 - epsumx)/vubmcmx) ;

  epssel = vubmcSB/vubmcr;
  errepssel = sqrt(epssel*(1 - epssel)/vubmcr) ;

  if(vubmcSB > vubmcallforeff){
    epsselbr = vubmcSB/vubmcmx;
    errepsselbr = sqrt(epsselbr*(1 - epsselbr)/vubmcmx);
  } else{
      epsselbr = epsumx * epschop;
      errepsselbr = sqrt(epsumx*epsumx*errepschop*errepschop + epschop*epschop*errepsumx*errepsumx);
  }


  c2->cd(7);
  sprintf(simply,"%s%f"," lepYaSe && chop>0. && chop<",MXBIN);
  signal= sighistounb(datamcvub, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, mesvubcuts[0],mesvubcuts[1],mesvubcuts[2], 5.,-60., 2500, 200, simply, USECB);
  N1 = signal.X(); N1err = signal.Y();  if(BINNED){}else{xframe->Draw();};
  cout<< "....N1  " << N1 << " +/- " << N1err << endl;

  c2->cd(8);
  sprintf(simply,"%s%f"," lepYaSe2 && chop>",MXBIN);
  signal= sighistounb(datamcvubtrue, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, mesvubcuts[0],mesvubcuts[1],mesvubcuts[2], 5.,-60., 2500, 200, simply, USECB);
  N2 = signal.X(); N2err = signal.Y();  if(BINNED){}else{xframe->Draw();};
  cout<< "....N2  " << N2 << " +/- " << N2err << endl;

  c2->cd(9);
  sprintf(simply,"%s%f"," lepYaSe3 && chop>0. && chop<",MXBIN);
  signal= sighistounb(datamcvubfalse, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, mesvubcuts[0],mesvubcuts[1],mesvubcuts[2], 5.,-60., 2500, 200, simply, USECB);
  N3 = signal.X();  N3err = signal.Y(); if(BINNED){}else{xframe->Draw();};
  cout<< "....N3  " << N3 << " +/- " << N3err << endl;



  sprintf(line,"%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),"mES_b2ueff.eps");
  c2->Print(line);
}

//------------------------------------------------------------------------
void b2uClass::b2umxBinning1d(double width){

  cout << "Binning with width " << width << endl;
  int ini=0;
  double recenter = 0;
  for (int k=0; k<100; k++) {

    ini++;
    mxB1[k]=recenter+k*width;

    if ((MXBIN<mxB1[k])&&(MXBIN>mxB1[k-1])) {
      recenter = (MXBIN-mxB1[k-1]);
      mxB1[k-1] = MXBIN;
      mxB1[k] = mxB1[k-1] + width;
    }

    if (mxB1[k]>5.0) {
      mxB1[k] = 5.0;
      break;
    }

  }

  for (int k=0; k<ini; k++) {
    cout << " 1d Mx bin [" << k << "] = " << mxB1[k] << endl;
  }
  nB = ini;
}

//------------------------------------------------------------------------
void b2uClass::b2uGetWeights(TString b2uwFile){
  ifstream fin;
  fin.open(b2uwFile);
  char buff[250];
  float tmp;
  int i=0;
  while(fin.getline(buff,250,'\n'))
    {
      if(buff[0]!='#')
	{
	  sscanf(buff, "%f",&tmp);
	  fweights[i] = tmp;
	  cout << "setting fweight[" << i << "] = " << tmp << endl;
	  i++;
	}
    }

}

// ----------------------------------------------------------------------
int b2uClass::b2unewrHistmx(double cmx){

  //double   histRbA[8] = {0.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.5};	  
  double   histRbA[8] = {0.0, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5} ; 
  // Mx categories
  for(int i =0; i <8; i++)   if(cmx<histRbA[i]) return i-1;
  
  return 8;
  
}

// ----------------------------------------------------------------------
void b2uClass::b2upstar(int su){

  cout << "b2upstar :: " << vubmc << " " << vcbmc << endl;

  double Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar;
  char simply[100], line[200];
  double chid = 0.181;
  //TVector2 signalsig;
  TVector2 signal;
  TCanvas *c2  = new TCanvas("c2"," ",200,10,1200,1000); 
  c2->Clear();
  c2->Divide(3, 2);

  c2->cd(1); 
  sprintf(simply,"(flavB==3)");
  signal = sighistounb(datamcvubNC, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, 5.27936, 0.0028, 1.28, 5.,-60., -6000, -1000, simply, USECB);
  //signal = sighistounb(datamcvubNC, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, mesvubMC[0], mesvubMC[1], mesvubMC[2], 5.,-60., 6000, 1000, simply, USECB);
  cout<<"b2uFitting NoCuts events on vub mc3!!! :: "<<  Nresmean<<" "<< Nressigma<<" "<< Nresalpha<<" "<< Nresn<<" "<< Nresargpar <<endl;
  double  vub3 = signal.X();                    double  errvub3 = signal.Y();       
  if(BINNED){}else{xframe->Draw();};

  c2->cd(2); 
  sprintf(simply,"(flavB==4)");
  //signal = sighistounb(datamcvubNC, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, mesvubMC[0], mesvubMC[1], mesvubMC[2], 5.,-60., 6000, 1000, simply, USECB);
  signal = sighistounb(datamcvubNC, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, 5.27936, 0.0028, 1.28, 5.,-60., 6000, 1000, simply, USECB);
  cout<<"b2uFitting NoCuts events on vub mc4!!! :: "<<  Nresmean<<" "<< Nressigma<<" "<< Nresalpha<<" "<< Nresn<<" "<< Nresargpar <<endl;
  double  vub4 = signal.X();                    double errvub4 = signal.Y();       
  if(BINNED){}else{xframe->Draw();};

  c2->cd(3);
  sprintf(simply,"(flavB==5)");
  //signal = sighistounb(datamcvubNC, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, mesvubMC[0], mesvubMC[1], mesvubMC[2], 5.,-60., 6000, 1000, simply, USECB);
  signal = sighistounb(datamcvubNC, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, 5.27936, 0.0028, 1.28, 5.,-60., 6000, 1000, simply, USECB);
  cout<<"b2uFitting NoCuts events on vub mc5!!! :: "<<  Nresmean<<" "<< Nressigma<<" "<< Nresalpha<<" "<< Nresn<<" "<< Nresargpar <<endl;
  double  vub5 = signal.X();                    double  errvub5 = signal.Y();       
  if(BINNED){}else{xframe->Draw();};
  
  if(MIXCORR==0){
    vubmcNC = vub3 + vub4 +vub5;
    errvubmcNC = sqrt(errvub3 * errvub3 + errvub4 * errvub4 + errvub5 * errvub5);
  }else{      
    vubmcNC = vub3 + ((1-chid)/(1-2*chid)) * vub4 - (chid/(1-2*chid)) * vub5;
    errvubmcNC = sqrt(errvub3 * errvub3 +
		    ((1-chid)/(1-2*chid)) * ((1-chid)/(1-2*chid)) * errvub4 * errvub4 +
		    (chid/(1-2*chid)) * (chid/(1-2*chid)) * errvub5 * errvub5);
  }



  c2->cd(4);
  sprintf(simply,"(flavB==3)");
  //signal = sighistounb(datamcvcbNC, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, mesvcbMC[0], mesvcbMC[1], mesvcbMC[2], 5.,-60., 100000, 25000, simply, USECB);
  signal = sighistounb(datamcvcbNC, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, 5.27948, 0.00288, 1.22, 5.,-60., 100000, 25000, simply, USECB);
  cout<<"b2uFitting NoCuts events on vcb mc3!!! :: "<<  Nresmean<<" "<< Nressigma<<" "<< Nresalpha<<" "<< Nresn<<" "<< Nresargpar <<endl;
  double  vcb3 = signal.X();               double  errvcb3 = signal.Y();       
  if(BINNED){}else{xframe->Draw();};

  c2->cd(5);
  sprintf(simply,"(flavB==4)");
  //signal = sighistounb(datamcvcbNC, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, mesvcbMC[0], mesvcbMC[1], mesvcbMC[2], 5.,-60., 100000, 25000, simply, USECB);
  signal = sighistounb(datamcvcbNC, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, 5.27948, 0.00288, 1.22, 5.,-60., 100000, 25000, simply, USECB);
  cout<<"b2uFitting NoCuts events on vcb mc4!!! :: "<<  Nresmean<<" "<< Nressigma<<" "<< Nresalpha<<" "<< Nresn<<" "<< Nresargpar <<endl;
  double  vcb4 = signal.X();                   double   errvcb4 = signal.Y();       
  if(BINNED){}else{xframe->Draw();};

  c2->cd(6);
  sprintf(simply,"(flavB==5)");
  //signal = sighistounb(datamcvcbNC, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, mesvcbMC[0], mesvcbMC[1], mesvcbMC[2], 5.,-60., 100000, 25000, simply, USECB);
  signal = sighistounb(datamcvcbNC, Vmes, Nresmean, Nressigma, Nresalpha, Nresn, Nresargpar, 5.27948, 0.00288, 1.22, 5.,-60., 100000, 25000, simply, USECB);
  cout<<"b2uFitting NoCuts events on vcb mc5!!! :: "<<  Nresmean<<" "<< Nressigma<<" "<< Nresalpha<<" "<< Nresn<<" "<< Nresargpar <<endl;
  double  vcb5 = signal.X();                    double  errvcb5 = signal.Y();       
  if(BINNED){}else{xframe->Draw();};

  if(MIXCORR==0){
    vcbmcNC = vcb3 + vcb4 +vcb5;
    errvcbmcNC = sqrt(errvcb3 * errvcb3 + errvcb4 * errvcb4 + errvcb5 * errvcb5);
  }else{      
    vcbmcNC = vcb3 + ((1-chid)/(1-2*chid)) * vcb4 - (chid/(1-2*chid)) * vcb5;
    errvcbmcNC = sqrt(errvcb3 * errvcb3 +
		    ((1-chid)/(1-2*chid)) * ((1-chid)/(1-2*chid)) * errvcb4 * errvcb4 +
		    (chid/(1-2*chid)) * (chid/(1-2*chid)) * errvcb5 * errvcb5);
  }

  cout << "b2upstar :: " << vubmcNC << " " << vcbmcNC << endl;

  calcpstarfact  = (vubmc/vubmcNC)/(vcbmc/vcbmcNC);
  //errcalcpstarfact = TMath::Sqrt((errvubmc/vubmc)*(errvubmc/vubmc) +  (errvcbmc/vcbmc)*(errvcbmc/vcbmc)); // approximately
  errcalcpstarfact = TMath::Sqrt((errvubmc/vubmc)*(errvubmc/vubmc) +  (errvcbmc/vcbmc)*(errvcbmc/vcbmc) + (errvubmcNC/vubmcNC)*(errvubmcNC/vubmcNC) +  (errvcbmcNC/vcbmcNC)*(errvcbmcNC/vcbmcNC) ); // approximately

  sprintf(line,"%s%s%s",DIRNAME.Data(),PREFIXOUT.Data(),"mES_SLeff.eps");
  c2->Print(line);

}
