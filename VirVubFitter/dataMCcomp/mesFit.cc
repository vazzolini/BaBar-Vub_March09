#include "mesFit.hh"

#include "RooFitCore/RooAbsPdf.hh"
#include "RooFitCore/RooRealVar.hh"
#include "RooFitCore/RooExtendPdf.hh"
#include "RooFitCore/RooAddPdf.hh"

#include "RooFitModels/RooArgusBG.hh"
#include "RooFitModels/RooCBShape.hh"
#include "RooFitModels/RooGaussian.hh"

#include "RooFitCore/RooGlobalFunc.hh"
#include "RooFitCore/RooBinning.hh"
#include "RooFitCore/RooFitResult.hh"

#include "RooThorstenSig.hh"
#include "RooFitCore/RooHist.hh"
#include "RooCCB.hh"

#include "TCanvas.h"
#include "TPaveText.h"
#include "TPad.h"

//#include "TMesCor.hh"
#include <iostream>
#include <string>
#include <sstream>

ClassImp(mesFit);

mesFit::mesFit() :  mesFitModel(ArgusAndThosig), mes(NULL), mesparsetting("mesparsetting_R22.dat") {}
// -------------------------------------------------------------------------------------------------
mesFit::mesFit(ModelType model, RooRealVar *m, TFile *fout) :  mesFitModel(model), mes(m), mesparsetting("mesparsetting_R22.dat"), outfile(fout) {
  
}
// -------------------------------------------------------------------------------------------------

mesFit::~mesFit(){
  
  
}

TVector2* mesFit::fitModel(RooDataSet* dataset, int fixpar, vector<Double_t>& results, const TString& title, bool isMC, bool isAllCut){
  
  //Signal and Background events
  RooRealVar* nsig = new RooRealVar("S","number of sig events", 1., 0., 800000.);
  RooRealVar* nbpk = new RooRealVar("P","number of bpk events", 1., 0., 800000.);
  RooRealVar* nbkg = new RooRealVar("B","number of bkg events", 1., 0., 800000.);
  
  RooAbsPdf* model = 0;
  RooExtendPdf *sigpointer,*arguspointer;
  
  if(results.size() == 0) readMesParamFile();

  nsig->setMax(dataset->numEntries(kTRUE));
  nbpk->setMax(dataset->numEntries(kTRUE));
  nbkg->setMax(dataset->numEntries(kTRUE));
  
  nsig->setVal(0.6*dataset->sumEntries("mes > 5.27"));
  nbkg->setVal(dataset->sumEntries("mes > 5.27")); 
  
  switch(mesFitModel){
    
  case GaussFit: {
    RooAbsPdf* a = createArgus(*mes);
    
    // Gaussian for signal component
    RooRealVar* Rm = new RooRealVar("mean","mean of gaussian 1",5.28,5.275,5.285);
    RooRealVar* Rs = new RooRealVar("sigma","width of gaussians",.003,.002,.004);
    RooGaussian* g = new RooGaussian("g", "Gaussian", *mes, *Rm, *Rs);
    
    
    RooExtendPdf* ae = new RooExtendPdf("ae","ae", *a, *nbkg);
    RooExtendPdf* ge = new RooExtendPdf("ge","ge", *g, *nsig);
    //------------------------

    /* --- build final model --- */
    model = new RooAddPdf("gmodel","a+g",RooArgList(*ge,*ae),RooArgList(*nsig,*nbkg));
    sigpointer = ge; arguspointer = ae;
  } break;
  
  case ArgusAndCB: {
    RooAbsPdf* a = createArgus(*mes);
    RooAbsPdf* cb = createCB(*mes);

    RooExtendPdf* ae = new RooExtendPdf("ae","ae", *a, *nbkg);
    RooExtendPdf* cbe = new RooExtendPdf("cbe","cbe", *cb, *nsig);
    
    model = new RooAddPdf("model","argus+CB",RooArgList(*cbe,*ae),RooArgList(*nsig,*nbkg));

    sigpointer = cbe; arguspointer = ae;

  } break;
    
  case ArgusAndThosig: { //ARGUS And Thorsten Signal Function AKA Frankenstein
    
    /* --- Build Argus background PDF --- */
    RooAbsPdf* a = createArgus(*mes);
    
    /* --- Build Thorsten signal (AKA Frankenstein) PDF --- */
    RooAbsPdf* pSignal = createThorsten(*mes);

    /* --- build extended pdfs --- */
    //-----------------------ANTONIO 23-APR-2007
    // RooExtendPdf* ae  = new RooExtendPdf("ae", "ae", *a, *nbkg, "mesint");
    // RooExtendPdf* se = new RooExtendPdf("se","signal extended", *pSignal, *nsig, "mesint");
    
    RooExtendPdf* ae  = new RooExtendPdf("ae", "ae", *a, *nbkg);
    RooExtendPdf* se = new RooExtendPdf("se","signal extended", *pSignal, *nsig);
    //------------------------
    
    /* --- build final model --- */
    model = new RooAddPdf("model","a+tho",RooArgList(*se,*ae),RooArgList(*nsig,*nbkg));

    sigpointer = se; arguspointer = ae;
     
  } break;
    
  default:
    cout << "W A R N I N G ! No known mes fit model detected! This is an error!" << endl;
    break;
  }
  
  vector<Double_t> inputPar;
  if(results.size() == 0) {
    if(!isMC && !isAllCut) inputPar = mesNsl;
    else if(!isMC &&  isAllCut) inputPar = mesdatacuts;
    else if( isMC && !isAllCut) inputPar = mespstarMC;
    else if( isMC &&  isAllCut) inputPar = mespstarcuts;
  } else
    inputPar = results;

  int nBins = int((5.291-5.22)/0.0005);
  mes->setBins(nBins);
  mes->setRange("mesint", 5.27, 5.291);

  if (mesFitModel == GaussFit) {

    // get hook to parameters
    RooRealVar* ar     = getPointer(model, "ar");
    RooRealVar* cutoff = getPointer(model, "cutoff");

    // get hook to parameters
    RooRealVar* Rm     = getPointer(model, "mean");
    RooRealVar* Rs     = getPointer(model, "sigma");

    cutoff->setVal(5.29); cutoff->setConstant();
    
    // preset parameters and fix some if requested
    if (fixpar == 0) { // preset parameters only
      Rm->setVal(inputPar[iMean]); 
      Rs->setVal(inputPar[iSigma]);
      ar->setVal(inputPar[iArgus]);
    } else if (fixpar == 1) { // set and fix all parameters (if reasonable value is given)
      if (inputPar[iMean]  > 0.)   { Rm->setVal(inputPar[iMean]);  Rm->setConstant(); }
      if (inputPar[iSigma] > 0.)   { Rs->setVal(inputPar[iSigma]); Rs->setConstant(); }
      if (inputPar[iArgus] > -55.) { ar->setVal(inputPar[iArgus]); }
      
      //CB we fit the argus shape with USECB=1 ar.setConstant();
    }
    std::cout<< "Parameters before fitting:: " << Rm->getVal() << " " << Rs->getVal() << " " << ar->getVal() 
	     << " " << nsig->getVal() << " " << nbkg->getVal() << std::endl;
  }
  
  /* ========= Set Parameter for Argus and CB Model  ============ */
  
  if(mesFitModel == ArgusAndCB){

    // get hook to parameters
    RooRealVar* ar     = getPointer(model, "ar");
    RooRealVar* cutoff = getPointer(model, "cutoff");

    // get hook to parameters
    RooRealVar* Rm = getPointer(model, "mean");
    RooRealVar* Rs = getPointer(model, "sigma");
    RooRealVar* Ra = getPointer(model, "alpha");
    RooRealVar* Rn = getPointer(model, "n");


    cout<<"CHECKING PARAMETERS WHAT DO WE USE? " << Rm->GetName() << "= " << Rm->getVal() << "; "
	<< Rs->GetName() << "= " << Rs->getVal() << "; " << Ra->GetName() << "= " << Ra->getVal() << "; "
	<< Rn->GetName() << "= " << Rn->getVal() << endl;

    
    // DO WE WANT TO FIX THIS ON DATA TOO???//
    cutoff->setVal(5.2895); 
    cutoff->setConstant();
    
    // preset parameters and fix some if requested by flag fixpar
    if (fixpar == 0) {   // preset paramters and fix necessary parameters
      Rm->setVal(inputPar[iMean]); 
      Rs->setVal(inputPar[iSigma]);  
      Ra->setVal(inputPar[iAlpha]);
      Rn->setVal(inputPar[iN]);  //   Rn->setConstant();
      ar->setVal(inputPar[iArgus]);
    } else if(fixpar == 1) {   // set and fix all parameters (if reasonable value is given)
      if (inputPar[iMean]  > 0.)   { Rm->setVal(inputPar[iMean]); Rm->setConstant();  }
      if (inputPar[iSigma] > 0.)   { Rs->setVal(inputPar[iSigma]); Rs->setConstant();  }
      if (inputPar[iN] > 0.)       { Rn->setVal(inputPar[iN]);      Rn->setConstant();  }
      if (inputPar[iAlpha] > 0.)   { Ra->setVal(inputPar[iAlpha]);  Ra->setConstant();  }
      if (inputPar[iArgus] > -55.) { ar->setVal(inputPar[iArgus]); }
      //CB we fit the argus shape with USECB=1 ar.setConstant();
    }
    std::cout<< "Parameters before fitting:: " << Rm->getVal() << " " << Rs->getVal() << " " << Ra->getVal() << " " << Rn->getVal() << " " << ar->getVal() 
	     << " " << nsig->getVal() << " " << nbkg->getVal() << std::endl;
  }

  /* ========= Set Parameter for Argus and Thosig Model (AKA Frankenstein)  ============ */
  
  if(mesFitModel == ArgusAndThosig) {
    // get hook to parameters
    RooRealVar* ar     = getPointer(model, "ar");
    RooRealVar* cutoff = getPointer(model, "cutoff");
    
    // get hook to parameters
    RooRealVar* Thor        = getPointer(model,"ThoSigR"); 
    RooRealVar* Thosigma_r1 = getPointer(model,"sigma_r1"); 
    RooRealVar* Thoxc       = getPointer(model,"ThoSigXc"); 
    RooRealVar* Thosigma_r2 = getPointer(model,"sigma_r2"); 
    RooRealVar* Thosigma_l  = getPointer(model,"sigma_l"); 
    RooRealVar* Thon        = getPointer(model,"ThoSigN"); 
    RooRealVar* Thoalpha    = getPointer(model,"ThoSigAlpha"); 
    
    
    //     if (data != datadata) { // fix endpoint on MC
    //       cutoff->setVal(5.2895); cutoff->setConstant();
    //     }
    
    cutoff->setVal(5.2895); 
    cutoff->setConstant();
    
    // preset parameters and fix some if requested by flag fixpar
    if (fixpar == 0) {   // preset paramters and fix necessary parameters
      if (inputPar[iThoSigR] > 0.)     { Thor->setVal(inputPar[iThoSigR]);    Thor->setConstant(); }
      if (inputPar[iSigma_r1] > 0.)    { Thosigma_r1->setVal(inputPar[iSigma_r1]); }
      if (inputPar[iThoSigXc] > 0.)    { Thoxc->setVal(inputPar[iThoSigXc]); }          
      if (inputPar[iSigma_r2] > 0.)    { Thosigma_r2->setVal(inputPar[iSigma_r2]); Thosigma_r2->setConstant(); }
      if (inputPar[iSigma_l]  > 0.)    { Thosigma_l->setVal(inputPar[iSigma_l]); }
      if (inputPar[iThoSigN] > 0.)     { Thon->setVal(inputPar[iThoSigN]);         Thon->setConstant(); }
      if (inputPar[iThoSigAlpha] > 0.) { Thoalpha->setVal(inputPar[iThoSigAlpha]); Thoalpha->setConstant(); }
      ar->setVal(inputPar[iArgus]);
      
    } else if(fixpar == 1) {   // set and fix all parameters (if reasonable value is given)
      if (inputPar[iThoSigR] > 0.)     { Thor->setVal(inputPar[iThoSigR]);         Thor->setConstant(); }
      if (inputPar[iSigma_r1] > 0.)    { Thosigma_r1->setVal(inputPar[iSigma_r1]); Thosigma_r1->setConstant(); }
      if (inputPar[iThoSigXc] > 0.)    { Thoxc->setVal(inputPar[iThoSigXc]);       Thoxc->setConstant(); }
      if (inputPar[iSigma_r2] > 0.)    { Thosigma_r2->setVal(inputPar[iSigma_r2]); Thosigma_r2->setConstant(); }
      if (inputPar[iSigma_l]  > 0.)    { Thosigma_l->setVal(inputPar[iSigma_l]);   Thosigma_l->setConstant(); }
      if (inputPar[iThoSigN] > 0.)     { Thon->setVal(inputPar[iThoSigN]);         Thon->setConstant(); }
      if (inputPar[iThoSigAlpha] > 0.) { Thoalpha->setVal(inputPar[iThoSigAlpha]); Thoalpha->setConstant(); }
      if (inputPar[iArgus] > -55.) { ar->setVal(inputPar[iArgus]); }
      //CB we fit the argus shape with USECB=1 ar.setConstant();
    }
    
    std::cout<< "Parameters before fitting:: " << Thor->GetName()<<": "<<Thor->getVal() << " " 
	     << Thosigma_r1->GetName()<<": "<<Thosigma_r1->getVal() << " " 
	     << Thoxc->GetName()<<": "<<Thoxc->getVal() << " " 
	     << Thosigma_r2->GetName()<<": "<<Thosigma_r2->getVal()<< " " 
	     << Thosigma_l->GetName()<<": "<<Thosigma_l->getVal() << " " 
	     << Thon->GetName()<<": "<<Thon->getVal() << " " 
	     << Thoalpha->GetName()<<": "<<Thoalpha->getVal()<< " " 
	     << ar->GetName()<<": "<<ar->getVal() << " " 
	     << cutoff->GetName()<<": "<<cutoff->getVal() << std::endl;
  }
    
  cout << "========= mesFit: processing dataset " << dataset->GetName() << " with " << dataset->GetTitle() <<" cuts" 
       << ". Initialization values for nsig = " << nsig->getVal() 
       << "; nbkg = " << nbkg ->getVal() << endl; 
  
  // DO the actual fit

  switch(mesFitModel){
  case GaussFit: std::cout << "Fitting Gaussian instead of Crystal Ball" << std::endl; break;
  case ArgusAndCB:  std::cout << "Fitting Argus and Crystal Ball" << std::endl; break;
  case ArgusAndThosig: std::cout << "Fitting Argus and Frankenstein" << std::endl; break;
  }
  
  RooFitResult* r = 0; 
  r = model->fitTo(*dataset, RooFit::Save(kTRUE),RooFit::Minos(kFALSE),RooFit::Hesse(kTRUE),RooFit::Extended(kTRUE));

  Double_t sig_fraction =  sigpointer->createIntegral(*mes,RooArgSet(*mes),"mesint")->getVal();
  
  RooArgSet* comps = model->getComponents();

  std::cout << std::endl;
  std::cout << "Full tree for " << model->ClassName() << "::" << model->GetName();
  std::cout << " Components: "; comps->Print(); std::cout << std::endl;

  TIterator* iter= comps->createIterator();
  RooAbsArg* pdf = 0;
  while ((pdf=(RooAbsArg*)iter->Next())) {
    
    if (std::string(pdf->GetName()) == std::string(model->GetName())) continue;
    
    RooArgSet* params = pdf->getParameters(RooArgSet());
    
    std::cout << " Component " << pdf->ClassName() << "::" << pdf->GetName();
    std::cout << " Variables: "; params->Print("v"); std::cout << std::endl;
    
    delete params;
  }
  delete iter; delete comps;
  
  //determine nDOF
  Int_t floatParms = (r->floatParsFinal()).getSize();
  std::cout << " control : number of floating params ----> " << floatParms << std::endl; 

  // make plot
  xframe = mes->frame();
  xframe->SetName(dataset->GetName());
  xframe->SetTitle(dataset->GetName());
  
  // plot data on it
  if (dataset->IsA() == RooDataHist::Class()) 
    dataset->plotOn(xframe, RooFit::MarkerSize(0.25), RooFit::DataError(RooAbsData::SumW2));
  else {
    RooBinning tbins(mes->getBins(), mes->getMin(), mes->getMax(), "plotbinning");
    dataset->plotOn(xframe, RooFit::Binning(tbins), RooFit::MarkerSize(0.5), RooFit::DataError(RooAbsData::SumW2));
  }

  model->plotOn(xframe, RooFit::Components("ae"), RooFit::LineWidth(2), RooFit::LineColor(kGreen)) ;
  
  // plot total model
  model->plotOn(xframe, RooFit::LineWidth(2));
  model->paramOn(xframe, dataset, "Fit Results",1,"ne",0.15,0.55,0.8);

  // plot data on it (for display purposes?)
  if (dataset->IsA() == RooDataHist::Class()) {
    dataset->plotOn(xframe, RooFit::MarkerSize(0.25), RooFit::DataError(RooAbsData::SumW2));
  } else {
    RooBinning tbins(mes->getBins(), mes->getMin(), mes->getMax(), "plotbinning");
    dataset->plotOn(xframe, RooFit::Binning(tbins), RooFit::MarkerSize(0.25), RooFit::DataError(RooAbsData::SumW2));
  }
  
  std::cout << "xframe chisquare = " << xframe->chiSquare(floatParms) << std::endl;

  TPaveText*  tbox1 = new TPaveText(0.59, 0.94, 0.79, 0.99, "BRNDC");
  char line1[50];
  sprintf(line1, "#chi^{2} = %5.4f",  xframe->chiSquare(floatParms));
  tbox1->AddText(line1);
  xframe->addObject(tbox1);
  
  RooHist* pullhisto = xframe->pullHist();// ((RooHist*)(xframe->getObject(0)))->GetName(),((RooCurve*)(xframe->getObject(2)))->GetName() );
  RooPlot* pullframe = mes->frame();
  // pullframe->addObject(pullhisto); 
  pullframe->addPlotable(pullhisto);//,"",kFALSE,kTRUE);
  // pullframe->updateYAxis(pullhisto->getYAxisMin(),pullhisto->getYAxisMax(),pullhisto->getYAxisLabel());
  // pullframe->updateFitRangeNorm(pullhisto,kFALSE) ; //?? kTRUE ??
  pullframe->SetMarkerStyle(24);
  pullframe->SetNdivisions(504, "Y");
  pullframe->SetLabelSize(0.13, "Y");
  pullframe->SetLabelSize(0.13, "X");  
  pullframe->SetStats(0);
  pullframe->SetTitle("");
  pullframe->SetTitleSize(0.20, "Y");
  pullframe->SetTitleOffset(0.20, "Y");
  pullframe->SetYTitle("Pull");
  //  pullhisto->Print("Verbose");

  //compute chisquare probability
  Double_t probchi = TMath::Prob((xframe->chiSquare(floatParms)) * 
				 (pullhisto->GetN()-floatParms), (pullhisto->GetN()-floatParms));
  std::cout << "chisq probability = " << probchi << std::endl;

  TPaveText*  tbox2 = new TPaveText(0.79, 0.94, 0.99, 0.99, "BRNDC");
  char line2[50];
  sprintf(line2, "Prob(#chi^{2}) = %4.2f", probchi);
  tbox2->AddText(line2);
  xframe->addObject(tbox2);
 
  //CB put also pull distributions
  //  if(mPad != NULL) {delete mPad; mPad = 0;}
  //  if(pPad != NULL) {delete pPad; pPad = 0;}

  TCanvas *canvas = new TCanvas(title,dataset->GetTitle(),0,0,800,600);

  TPad* mPad = new TPad("mES plot",  "", 0.02, 0.30, 0.98, 0.98);  mPad->Draw(); 
  TPad* pPad = new TPad("pull plot", "", 0.02, 0.02, 0.98, 0.28);  pPad->Draw();
  
  mPad->cd();
  xframe->Draw(); 
  mPad->Modified();
  
  pPad->cd();
  pullframe->Draw(); 
  
  pPad->Modified();
  pPad->Update();

  canvas->Write();
  
  r->Print();
  cout << "Fit covQaul " << r->covQual() << " fit Status " << r->status()<< endl;
  //  if( r->covQual() < 3 ) fprintf(ff,"%s %s covQual %d\n",dataset->GetName(),dataset->GetTitle(),r->covQual());
  
  // delete pPad; delete mPad; delete pullframe; 
  // delete tbox1;  delete tbox2; delete xframe;
  

   RooArgSet *fit_parameters = model->getVariables();
   TIterator *iter2 = fit_parameters->createIterator();
   
   //   RooRealVar* current;
   //   while( (current=(RooRealVar*)iter2->Next() ) ){
   //     cout << current->GetName() << " " << current->getVal() << " " << current->getError() << endl;
   //   }
     

   //  for(int i=0;i<fit_parameters.getSize();i++){
   //    current = (RooRealVar*)fit_parameters[i];
   //  cout << current->getVal() <<" +- " << current->getError() <<endl;
   //  }
   
   double p0 = sig_fraction*nsig->getVal(); double Dp0 = sig_fraction*nsig->getError();
   
   if( results.size() ==0 ) {
     //get Best parameters.
     switch(mesFitModel){
     case GaussFit: {
       results.resize(5);
       results[iMean]  = getVal(model, "mean"); //pos 0
       results[iSigma] = getVal(model, "sigma"); //pos 1
     
       results[iArgus]  = getVal(model, "ar"); //pos 4
       results[iCutOff] = getVal(model, "cutoff"); //pos5
     } break;
     case ArgusAndCB:  {
       results.resize(5);
       results[iMean]     = getVal(model, "mean"); //pos 0
       results[iSigma]    = getVal(model, "sigma"); //pos 1
       results[iAlpha]    = getVal(model, "alpha"); //pos 2
       results[iN]        = getVal(model, "n"); // pos3

       results[iArgus]    = getVal(model, "ar"); //pos 4
       results[iCutOff]      = getVal(model, "cutoff");//pos5
     } break;
     case ArgusAndThosig: {
       results.resize(14);
       results[iThoSigR]     = getVal(model, "ThoSigR"); //pos 7
       results[iSigma_r1]    = getVal(model, "sigma_r1"); //pos 8
       results[iThoSigXc]    = getVal(model, "ThoSigXc");//pos 9
       results[iSigma_r2]    = getVal(model, "sigma_r2"); //pos 10
       results[iSigma_l]     = getVal(model, "sigma_l"); //pos 11
       results[iThoSigN]     = getVal(model, "ThoSigN"); //pos 12
       results[iThoSigAlpha] = getVal(model, "ThoSigAlpha"); //pos 13
     
       results[iArgus] =      getVal(model, "ar"); //pos4
       results[iCutOff]=      getVal(model, "cutoff"); //pos5
     } break;
     }
   }
   
   return new TVector2(p0,Dp0);
}

// -------------------------------------------------------------------------------------------------
RooAbsPdf* mesFit::createArgus(RooRealVar& mes)
{
  // return function
  RooAbsPdf* pFunction = 0;
  
  // argus parameters (same for all functions)
  RooRealVar* pArgPar = new RooRealVar("ar", "argus shape parameter", -60., -100., -10.);
  RooRealVar* pCutOff = new RooRealVar("cutoff", "argus cutoff", 5.2891, 5.288, 5.292);

  /*  if (endpointCorrection) {

      std::vector<Double_t> vWeight, vEndpoint;
      Double_t totalWeight = 0.;

      Mes Correction 
      for (size_t i = 0; i < mesCor->mesPeriodMax(); ++i) {
      Double_t endPoint = mesCor->mesEndpoint(i) - 5.2891;
      Double_t weight = mesCor->mesWeight(i);

      if (weight == 0.) continue;
      totalWeight += weight;

      vWeight.push_back(weight);
      vEndpoint.push_back(endPoint);
      }
    
      for (size_t i = 0; i < vWeight.size(); ++i) vWeight[i] /= totalWeight;
    
      pFunction = new RooSumArgusBG("a", "Argus PDF", mes, *pCutOff, *pArgPar, vWeight, vEndpoint);
      pCutOff->setConstant();
    
      } else {
  */
  //    if (_debug > 1) std::cout << "Building single Argus" << std::endl;
  pFunction = new RooArgusBG("a", "Argus PDF", mes, *pCutOff, *pArgPar);
  
  //}
  
  //  if (_debug > 2) pFunction->printCompactTree();
  
  return pFunction;
}

// -------------------------------------------------------------------------------------------------
RooAbsPdf* mesFit::createCB(RooRealVar& mes)
{
  // return function
  RooAbsPdf* pFunction = 0;

  // modified crystal ball parameters (same for all functions)
  RooRealVar* Rm = new RooRealVar("mean","cb: mean of gaussian 1", 5.28, 5.275, 5.285);
  RooRealVar* Rs = new RooRealVar("sigma","cb: width of gaussians", 0.003, 0.002, 0.004);
  RooRealVar* Ra = new RooRealVar("alpha","cb: alpha parameter", 1.3, 0., 10.);
  RooRealVar* Rn = new RooRealVar("n","cb: n parameter", 3.46, 1., 7.);

  //  if (_debug > 1) std::cout << "Building Crystal Ball Function" << std::endl;
  pFunction = new RooCBShape("cb", "Crystal Ball", mes, *Rm, *Rs, *Ra, *Rn);

  //  if (_debug > 2) pFunction->printCompactTree();

  return pFunction;
}
// -------------------------------------------------------------------------------------------------
RooAbsPdf* mesFit::createThorsten(RooRealVar& mes)
{
  // modified gauz parameters
  RooRealVar* Thor =        new RooRealVar("ThoSigR","thosig r",0.9,0.2,1);
  RooRealVar* Thosigma_r1 = new RooRealVar("sigma_r1","sigma_r1",0.00235,0.00050,0.00290);
  RooRealVar* Thoxc =       new RooRealVar("ThoSigXc","Thosig xc",5.27966,5.270,5.2810);
  RooRealVar* Thosigma_r2 = new RooRealVar("sigma_r2","sigma_r2",0.00304,0.0015,0.0035);
  RooRealVar* Thosigma_l  = new RooRealVar("sigma_l","sigma_l",0.0017,0.00120,0.00280);
  RooRealVar* Thon =        new RooRealVar("ThoSigN","thosig n",1.2,0.1,3);
  RooRealVar* Thoalpha =    new RooRealVar("ThoSigAlpha","thosig alpha",4.2,3.4,4.5);

  //  if (_debug > 1) std::cout << "Building Thorstens Modified Triple Function" << std::endl;
  
  RooThorstenSig* pFunction = new RooThorstenSig("ThoSig","Thorstens Triple Function", mes, 
						 *Thor, *Thosigma_r1, *Thoxc, *Thosigma_r2, *Thosigma_l, *Thon, *Thoalpha);

  //  if (_debug > 2) pFunction->printCompactTree();

  return pFunction;
}

// -------------------------------------------------------------------------------------------------
RooRealVar* getPointer(RooAbsPdf* pdf, const char* name)
{
  RooArgSet* list = pdf->getParameters(RooArgSet());
  RooRealVar* var = dynamic_cast<RooRealVar*>(list->find(name));
  delete list;

  return var;
}
// -------------------------------------------------------------------------------------------------
//! helper access method used by vubMesUnb
Double_t getVal(RooAbsPdf* pdf, const char* name)
{
  RooRealVar* var = getPointer(pdf, name);
  return (var != 0 ? var->getVal() : 0.);
}
// -------------------------------------------------------------------------------------------------
//! helper access method used by vubMesUnb
Double_t getError(RooAbsPdf* pdf, const char* name)
{
  RooRealVar* var = getPointer(pdf, name);
  return (var != 0 ? var->getError() : 0.);
}
// -------------------------------------------------------------------------------------------------

void mesFit::setMesParamFile(const TString& filename){
  mesparsetting = filename;
}


void mesFit::readMesParamFile()
{
  std::ifstream is(mesparsetting.Data());
  if (!is) { 
    std::cout << "Error: Can't open file " << mesparsetting.Data() << ". Stopping execution!" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  int pos(0), sign(0);
  std::vector<double> mesSysNsl, mesSysdatacuts, mesSysvubcuts, mesSysvcbcuts, mesSysothcuts;
  std::vector<double> mesSysvubMC, mesSysvcbMC, mesSysvubMCchop, mesSysvuboutMC, mesSyspstarMC, mesSyspstarcuts;

  std::string buffer;
  while(getline(is, buffer, '\n')){
    if (buffer[0] == '#') continue;
    
    std::istringstream bufferstream(buffer.c_str());
    
    int ok = 0;
    
    string tag; 

    //    float input; // xxx change this to double after all modifications
    double input; // Done As prescription
    std::vector<double> inputs;

    bufferstream >> tag;
    while (bufferstream >> input) { inputs.push_back(input); }


    //    if (_debug > 2) std::cout << "Found tag: " << tag << " values: " << inputs << std::endl; 
    
    // hack for compatibility
    if (inputs.size() == 4) inputs.push_back(-60.);

    // look up the right place to store
    if (tag == "Nsl")          { mesNsl = inputs; ok = 1; }
    if (tag == "datacuts")     { mesdatacuts = inputs; ok = 1; }
    if (tag == "vubcuts")      { mesvubcuts = inputs; ok = 1; }
    if (tag == "vcbcuts")      { mesvcbcuts = inputs; ok = 1; }
    if (tag == "othcuts")      { mesothcuts = inputs; ok = 1; }
    if (tag == "NslMC")        { mesNslMC = inputs; ok = 1; }
    if (tag == "vubMC")        { mesvubMC = inputs; ok = 1; }
    if (tag == "vcbMC")        { mesvcbMC = inputs; ok = 1; }
    if (tag == "vubMCchop")    { mesvubMCchop = inputs; ok = 1; }
    if (tag == "vubMCall")     { mesvubMCall = inputs; ok = 1; }
    if (tag == "vcbMCall")     { mesvcbMCall = inputs; ok = 1; }
    if (tag == "vubMClepteff") { mesvubMClepteff = inputs; ok = 1; }
    if (tag == "vubMCalleff")  { mesvubMCalleff = inputs; ok = 1; }
    if (tag == "vuboutMC" )    { mesvuboutMC = inputs; ok = 1; }
    if (tag == "pstarMC" )     { mespstarMC = inputs; ok = 1; }
    if (tag == "pstarcuts" )   { mespstarcuts = inputs; ok = 1; }

    if (tag == "Nslerror")       { mesSysNsl = inputs; ok = 1; }
    if (tag == "datacutserror")  { mesSysdatacuts = inputs; ok = 1; }
    if (tag == "vubcutserror")   { mesSysvubcuts = inputs; ok = 1; }
    if (tag == "vcbcutserror")   { mesSysvcbcuts = inputs; ok = 1; }
    if (tag == "othcutserror")   { mesSysothcuts = inputs; ok = 1; }
    if (tag == "vubMCerror")     { mesSysvubMC = inputs; ok = 1; }
    if (tag == "vcbMCerror")     { mesSysvcbMC = inputs; ok = 1; }
    if (tag == "vubMCchoperror") { mesSysvubMCchop = inputs; ok = 1; }
    if (tag == "vuboutMCerror")  { mesSysvuboutMC = inputs; ok = 1; }
    if (tag == "pstarMCerror")   { mesSyspstarMC = inputs; ok = 1; }
    if (tag == "pstarcutserror") { mesSyspstarcuts = inputs; ok = 1; }

    if (ok == 0)  std::cout << "==> fitNtp::readmesParam() Error: Don't know about variable " << tag << std::endl;
  }
  
  is.close();

  //  for(UInt_t i=0; i< mesNsl.size(); i++)
  //    cout << "mesNsl["<<i<<"]= "<<mesNsl[i]<< endl;
}

