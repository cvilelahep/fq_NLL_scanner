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
#include <vector>
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


int fq_NLL_scanner::ScanNLL(){
  fiTQun_shared::nring = 1;

  int PCflg;


  // Find true parameters of primary lepton?
  
  bool printStuff = true;
  
  int temp_track = -1;
#ifndef NOSKLIBRARIES
  /* load mc_particle structure */
  int nmctrks;
  mcguts(&nmctrks);

  float trgofst;
  vcrdvccm_();
  trginfo_(&trgofst);

  if (printStuff) std::cout << "This event has " << nmctrks << " true tracks" << std::endl;
  for (int itrk=0; itrk<nmctrks; itrk++) {
    int parentid = mc_particle.parent[itrk];
    int inOutOrSecondary = mc_particle.origin[itrk];
    int pdgCode = mc_particle.type[itrk];
    float mass = mc_particle.mass[itrk];
    float momentum = mc_particle.momentum[itrk];
    float energy = mc_particle.energy[itrk];
    float beta = mc_particle.beta[itrk];
    float dir[3];
    float inivtx[3];
    float finvtx[3];
    for (int ix=0; ix<3; ix++) {
      dir[ix] = mc_particle.dir[itrk][ix];
      /* NUANCE outgoing particles get mislabelled as incoming, so if the
         initial_vertex is not set, use the the_vertex instead */
      inivtx[ix] = mc_particle.initial_vertex[itrk][ix]; // so this is broken for Nuance input files
      finvtx[ix] = mc_particle.final_vertex[itrk][ix];
    } float time = mc_particle.time[itrk];
    if (printStuff) {
      std::cout << "track number " << itrk << " has PDG code, momentum, parent = " << pdgCode << ", " << momentum << ", " << parentid << \
        std::endl;
      std::cout << "direction = (" << dir[0] << "," << dir[1] << "," << dir[2] << ")" << std::endl;
      std::cout << "inivtx = (" << inivtx[0] << "," << inivtx[1] << "," << inivtx[2] << ")" << std::endl;
      std::cout << "finvtx = (" << finvtx[0] << "," << finvtx[1] << "," << finvtx[2] << ")" << std::endl;
      std::cout << "time, trgofst, mass, beta = " << time << ", " << trgofst << ", " << mass << ", " << beta << std::endl;
    }
  }
#endif

  //      0   1   2   3     4      5    6   7   8   9   10   11     12    13  
  // X = {X1, Y1, Z1, T1, theta1, phi1, p1, X2, Y2, Z2, T2, theta2, phi2, p2}
  int itrk = 0;

  std::cout << "TRUE MOMENTUM " << mc_particle.momentum[0] << std::endl;
  
  double Xtrue[] =  {mc_particle.initial_vertex[itrk][0], mc_particle.initial_vertex[itrk][1], mc_particle.initial_vertex[itrk][2],
                     mc_particle.time[itrk],
                     mc_particle.dir[itrk][2], atan2(mc_particle.dir[itrk][0], mc_particle.dir[itrk][1]), // Not sure if direction is consistent with rest of fiTQun... CHECK.
                     mc_particle.momentum[itrk]};
  

  //  double X[] = {0., 0., 0., 0., 0., 0., 500};
  double  snglTrkParams[7];
  for (int i = 0; i < 7; i++) snglTrkParams[i] = Xtrue[i];

  // Scan over momentum range:
  std::vector<double> xx;
  std::vector<double> yy;
  
  for (double momFrac = 0.5; momFrac < 2; momFrac+=0.1){
    Resetmu();

    double thisMom = Xtrue[6]*momFrac;
    xx.push_back(thisMom);

    snglTrkParams[6] = thisMom;
    // MUON IS HARD-CODED, NEED TO FIX IF RUNNING ON OTHER PARTICLES! - CV
    int PCflg1 = OneRing(imu,snglTrkParams,0); 
    double nll = nglogL();
    yy.push_back(nll);
    std::cout << "ScanNLL: nll is " << nll << " for p = " << thisMom << " true momentum " << Xtrue[6] << " mom frac " << momFrac <<  std::endl;

    
  }


  TFile * fOut = new TFile("TestOut.root", "RECREATE");
  TGraph * gOut = new TGraph(xx.size(), &(xx[0]), &(yy[0]));
  gOut->SetName("nllScan");
  gOut->Write();
  fOut->Close();

}

