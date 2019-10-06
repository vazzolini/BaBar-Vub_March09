#include "exampleAnalysis.hh"

exampleAnalysis::exampleAnalysis(TTree *tree, int isMC, int newFormat):baseClass(tree, isMC, newFormat) {
  std::cout << "exampleAnlaysis has been called" << std::endl;
}

// ----------------------------------------------------------------------
//void exampleAnalysis::recoil() {
//
//  put your own recoil code here, for example
//
//}

// ----------------------------------------------------------------------
//void exampleAnalysis::MyMethod() {
//
//  define your methods as needed
//  
//}
