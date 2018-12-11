#include <iostream>
#include "fiTQun.h"
#include "fq_NLL_scanner.h"
#include "TMath.h"
#include "TFitter.h"
#include <TMinuit.h>
#include "fiTQun_shared.h"
#include "mcguts.h"
#include "apscndryC.h"
#include <fstream>

using namespace fiTQun_parameters;

extern"C" {

  void mcguts(int *ntracks);
  void trginfo_(float *);
  void vcrdvccm_();

}

fq_NLL_scanner::fq_NLL_scanner(int anpmt) 
  : fiTQun(anpmt)
{
  std::cout << "Starting fq_NLL_scanner" << std::endl;
  //  static_fq_NLL_scanner = this;
  //GetTruefq_NLL_scannerParams(tParams);
}

fq_NLL_scanner::~fq_NLL_scanner() {
  std::cout << "fq_NLL_scanner routine ended." << std::endl;
}


int fq_NLL_scanner::ScanNLL(int uPID, double *snglTrkParams){
  fiTQun_shared::nring = 1;

  int PCflg;


  Resetmu();

  //      0   1   2   3     4      5    6   7   8   9   10   11     12    13  
  // X = {X1, Y1, Z1, T1, theta1, phi1, p1, X2, Y2, Z2, T2, theta2, phi2, p2}

  double X[] = {0., 0., 0., 0., 0., 0., 500};
  int PCflg1 = OneRing(uPID,snglTrkParams,0); 
  double nll = nglogL();
  std::cout << "ScanNLL: nll is " << nll << std::endl;

}

