//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: TMesCor.cc,v 1.1 2006/05/18 17:44:56 menges Exp $
//
// Description:
//	Class TMesCor. Wrapper class for mes corrections.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Wolfgang Menges
//
// Copyright Information:
//      Copyright (C) 2006      Queen Mary, University of London
//
//------------------------------------------------------------------------

//-----------------------
// This Class's Header --
//-----------------------
#include "VirVubFitter/TMesCor.hh"

//---------------
// C++ Headers --
//---------------
#include <iostream>

ClassImp(TMesCor)

//----------------
// Constructors --
//----------------

TMesCor::TMesCor(const Int_t runFlag)
{
  initMesCorrections(runFlag);
}

//--------------
// Destructor --
//--------------

TMesCor::~TMesCor()
{
}

//-------------
// Methods   --
//-------------

Int_t TMesCor::mesPeriod(const Int_t runNumber)
{
  if (runNumber <  mesRunLimits[0].first) { 
    std::cout << "Warning: run period (" << runNumber << ") out of range for mes corrections" << std::endl;
    return -1;
  }

  Int_t index = -1;
  for (size_t i = 0;  i < mesRunLimits.size(); ++i) {
    if (runNumber > mesRunLimits[i].second) continue;
    index = i;
    break;
  }

  if (index == -1) {
    std::cout << "Warning: run period (" << runNumber << ") out of range for mes corrections" << std::endl;
    return -1;
  }

  if (runNumber >= mesRunLimits[index].first && runNumber <= mesRunLimits[index].second) return index;

  std::cout << "Warning: run periods for mes corrections is inconsistend: " << index << std::endl;

  return -2;
}

Double_t TMesCor::mesWeight(const Int_t runNumber)
{
  Int_t index = mesPeriod(runNumber);
  
  if (index < 0) {
    return 0.;
  }

  return mesCorTable[index][TMesCor::corWeight];
}

Double_t TMesCor::mesCorrections(const Int_t runNumber)
{
  Int_t index = mesPeriod(runNumber);
  
  if (index < 0) {
    return -5.;
  }

  return mesCorTable[index][TMesCor::corMes];
}

Double_t TMesCor::mesEndpoint(const Int_t runNumber)
{
  Int_t index = mesPeriod(runNumber);
  
  if (index < 0) {
    return -5.;
  }

  return mesCorTable[index][TMesCor::corEndpoint];
}

void TMesCor::initMesCorrections(const Int_t runFlag)
{
  // 
  mesRunLimits.resize(0);

  mesRunLimits.push_back(std::pair<double, double>(9923,  11341));
  mesRunLimits.push_back(std::pair<double, double>(11342, 12339));
  mesRunLimits.push_back(std::pair<double, double>(12340, 13340));
  mesRunLimits.push_back(std::pair<double, double>(13341, 14036));
  mesRunLimits.push_back(std::pair<double, double>(14037, 14462));
  mesRunLimits.push_back(std::pair<double, double>(14463, 15090));
  mesRunLimits.push_back(std::pair<double, double>(15091, 15827));
  mesRunLimits.push_back(std::pair<double, double>(15828, 16827));
  mesRunLimits.push_back(std::pair<double, double>(16828, 17106));
  mesRunLimits.push_back(std::pair<double, double>(17107, 18148));
  mesRunLimits.push_back(std::pair<double, double>(18149, 19194));
  mesRunLimits.push_back(std::pair<double, double>(19195, 20294));
  mesRunLimits.push_back(std::pair<double, double>(20295, 21251));
  mesRunLimits.push_back(std::pair<double, double>(21252, 22252));
  mesRunLimits.push_back(std::pair<double, double>(22253, 22538));
  mesRunLimits.push_back(std::pair<double, double>(22539, 23540));
  mesRunLimits.push_back(std::pair<double, double>(23541, 23987));
  mesRunLimits.push_back(std::pair<double, double>(23988, 25305));
  mesRunLimits.push_back(std::pair<double, double>(25306, 25758));
  mesRunLimits.push_back(std::pair<double, double>(25759, 26225));
  mesRunLimits.push_back(std::pair<double, double>(26226, 27361));
  mesRunLimits.push_back(std::pair<double, double>(27362, 28467));
  mesRunLimits.push_back(std::pair<double, double>(28468, 29435));
  mesRunLimits.push_back(std::pair<double, double>(29436, 32954));
  mesRunLimits.push_back(std::pair<double, double>(32955, 33958));
  mesRunLimits.push_back(std::pair<double, double>(33959, 34959));
  mesRunLimits.push_back(std::pair<double, double>(34960, 35645));
  mesRunLimits.push_back(std::pair<double, double>(35646, 36645));
  mesRunLimits.push_back(std::pair<double, double>(36646, 37279));
  mesRunLimits.push_back(std::pair<double, double>(37280, 38280));
  mesRunLimits.push_back(std::pair<double, double>(38281, 39320));
  mesRunLimits.push_back(std::pair<double, double>(39321, 40002));
  mesRunLimits.push_back(std::pair<double, double>(40003, 41229));
  mesRunLimits.push_back(std::pair<double, double>(41230, 42252));
  mesRunLimits.push_back(std::pair<double, double>(42253, 42676));
  mesRunLimits.push_back(std::pair<double, double>(42677, 43677));
  mesRunLimits.push_back(std::pair<double, double>(43678, 44678));
  mesRunLimits.push_back(std::pair<double, double>(44679, 45645));
  mesRunLimits.push_back(std::pair<double, double>(45646, 46646));
  mesRunLimits.push_back(std::pair<double, double>(46647, 47647));
  mesRunLimits.push_back(std::pair<double, double>(47648, 48652));
  mesRunLimits.push_back(std::pair<double, double>(48653, 49653));
  mesRunLimits.push_back(std::pair<double, double>(49654, 50635));

  mesCorTable.resize(0);

  std::vector<Double_t> dummy; 
  dummy.resize(3);

  dummy[0] = -0.00031; dummy[1] = 5.290625; dummy[2] =   468.80; mesCorTable.push_back(dummy);
  dummy[0] = +0.00015; dummy[1] = 5.291125; dummy[2] =  1343.90; mesCorTable.push_back(dummy);
  dummy[0] = -0.00037; dummy[1] = 5.291375; dummy[2] =  4393.10; mesCorTable.push_back(dummy);
  dummy[0] = -0.00039; dummy[1] = 5.289125; dummy[2] =  2089.90; mesCorTable.push_back(dummy);
  dummy[0] = -0.00131; dummy[1] = 5.289125; dummy[2] =  1763.60; mesCorTable.push_back(dummy);
  dummy[0] = -0.00123; dummy[1] = 5.289625; dummy[2] =  1860.80; mesCorTable.push_back(dummy);
  dummy[0] = -0.00110; dummy[1] = 5.289875; dummy[2] =  2724.60; mesCorTable.push_back(dummy);
  dummy[0] = -0.00110; dummy[1] = 5.290125; dummy[2] =  3731.60; mesCorTable.push_back(dummy);
  dummy[0] = -0.00117; dummy[1] = 5.289875; dummy[2] =  1114.60; mesCorTable.push_back(dummy);
  dummy[0] =  0.0    ; dummy[1] = 0.      ; dummy[2] =     0.00; mesCorTable.push_back(dummy);
  dummy[0] = -0.00069; dummy[1] = 5.289625; dummy[2] =  1977.60; mesCorTable.push_back(dummy);
  dummy[0] = -0.00091; dummy[1] = 5.289875; dummy[2] =  3658.60; mesCorTable.push_back(dummy);
  dummy[0] = -0.00097; dummy[1] = 5.289625; dummy[2] =  5500.40; mesCorTable.push_back(dummy);
  dummy[0] = -0.00074; dummy[1] = 5.290125; dummy[2] =  5965.00; mesCorTable.push_back(dummy);
  dummy[0] = -0.00093; dummy[1] = 5.290125; dummy[2] =   958.30; mesCorTable.push_back(dummy);
  dummy[0] = -0.00101; dummy[1] = 5.289875; dummy[2] =  7005.30; mesCorTable.push_back(dummy);
  dummy[0] = -0.00106; dummy[1] = 5.289625; dummy[2] =  2726.60; mesCorTable.push_back(dummy);
  dummy[0] = -0.00065; dummy[1] = 5.289625; dummy[2] =  7202.40; mesCorTable.push_back(dummy);
  dummy[0] = -0.00066; dummy[1] = 5.289875; dummy[2] =  1984.60; mesCorTable.push_back(dummy);
  dummy[0] = -0.00071; dummy[1] = 5.290375; dummy[2] =  1724.90; mesCorTable.push_back(dummy);
  dummy[0] = -0.00095; dummy[1] = 5.290625; dummy[2] =  8204.20; mesCorTable.push_back(dummy);
  dummy[0] = -0.00094; dummy[1] = 5.290625; dummy[2] =  7521.50; mesCorTable.push_back(dummy);
  dummy[0] = -0.00081; dummy[1] = 5.290125; dummy[2] =  5837.40; mesCorTable.push_back(dummy);
  dummy[0] =  0.0    ; dummy[1] = 0.      ; dummy[2] =     0.00; mesCorTable.push_back(dummy);
  dummy[0] = -0.00088; dummy[1] = 5.290125; dummy[2] =  3130.80; mesCorTable.push_back(dummy);
  dummy[0] = -0.00094; dummy[1] = 5.290375; dummy[2] =  4218.10; mesCorTable.push_back(dummy);
  dummy[0] = -0.00069; dummy[1] = 5.290375; dummy[2] =  3205.90; mesCorTable.push_back(dummy);
  dummy[0] = -0.00071; dummy[1] = 5.290875; dummy[2] =  6262.20; mesCorTable.push_back(dummy);
  dummy[0] = -0.00090; dummy[1] = 5.289875; dummy[2] =  2443.60; mesCorTable.push_back(dummy);
  dummy[0] = -0.00070; dummy[1] = 5.289625; dummy[2] =  5613.60; mesCorTable.push_back(dummy);
  dummy[0] = -0.00057; dummy[1] = 5.289375; dummy[2] =  6212.10; mesCorTable.push_back(dummy);
  dummy[0] =  0.0    ; dummy[1] = 0.      ; dummy[2] =     0.00; mesCorTable.push_back(dummy);
  dummy[0] = -0.00068; dummy[1] = 5.289375; dummy[2] =  5220.20; mesCorTable.push_back(dummy);
  dummy[0] = -0.00055; dummy[1] = 5.289375; dummy[2] =  6726.70; mesCorTable.push_back(dummy);
  dummy[0] = -0.00050; dummy[1] = 5.289375; dummy[2] =  1372.00; mesCorTable.push_back(dummy);
  dummy[0] = -0.00041; dummy[1] = 5.289375; dummy[2] =  8217.00; mesCorTable.push_back(dummy);
  dummy[0] = -0.00023; dummy[1] = 5.289375; dummy[2] =  9910.50; mesCorTable.push_back(dummy);
  dummy[0] = -0.00024; dummy[1] = 5.289875; dummy[2] =  5678.60; mesCorTable.push_back(dummy);
  dummy[0] = -0.00022; dummy[1] = 5.289625; dummy[2] =  8774.40; mesCorTable.push_back(dummy);
  dummy[0] = -0.00022; dummy[1] = 5.289125; dummy[2] = 14625.40; mesCorTable.push_back(dummy);
  dummy[0] = -0.00032; dummy[1] = 5.289125; dummy[2] = 15121.60; mesCorTable.push_back(dummy);
  dummy[0] = -0.00058; dummy[1] = 5.289375; dummy[2] = 14346.10; mesCorTable.push_back(dummy);
  dummy[0] = -0.00079; dummy[1] = 5.289875; dummy[2] = 10118.30; mesCorTable.push_back(dummy);

  return;
}
