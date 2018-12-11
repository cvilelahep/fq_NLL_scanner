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

  // Find true parameters of primary lepton?
  
  void PDK_MuGamma::GetTruePDK_MuGammaParams(double* tParams){ //NEEDS A LOT OF CHANGES!!!!!!!!!! not sure how to do them.

    bool printStuff = true;
  
    int temp_track = -1;

    for (int ipar=0; ipar<fiTQun_shared::nPDK_MuGammaParams; ipar++) {
      tParams[ipar] = 0.;
    }

    bool printstuff = true;

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

      if (pdgCode == 321 && parentid == 0)
	temp_track = itrk;

      if (pdgCode == -13 && parentid == temp_track) {
	tParams[0] = inivtx[0];
	tParams[1] = inivtx[1];
	tParams[2] = inivtx[2];
	tParams[3] = trgofst + time;
	tParams[4] = acos(dir[2]);
	tParams[5] = atan2(dir[1],dir[0]);
	tParams[6] = momentum;
      }

      if (pdgCode == 22 && parentid == 0) {
	tParams[7] = inivtx[0];
	tParams[8] = inivtx[1];
	tParams[9] = inivtx[2];
	if (printStuff) std::cout << "setting tParams[3] to " << trgofst + time << std::endl;
	tParams[10] = trgofst + time;
	tParams[11] = acos(dir[2]);
	tParams[12] = atan2(dir[1],dir[0]);
	tParams[13] = momentum;
      }

    }
#endif



  double X[] = {0., 0., 0., 0., 0., 0., 500};
  int PCflg1 = OneRing(uPID,snglTrkParams,0); 
  double nll = nglogL();
  std::cout << "ScanNLL: nll is " << nll << std::endl;

}

