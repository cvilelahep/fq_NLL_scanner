#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>

#include "fiTQun.h"
#include "fiTQun_shared.h"
#include "spliTChan.h"
#include "PDK_MuGamma.h"

#include "fq_NLL_scanner.h"

#include "fQROOTOut.h"

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>

#ifdef NOSKLIBRARIES
#include "WCSimWrap.h"
#else
#include "skheadC.h"
#include "apmringC.h"
#include "apmueC.h"
#include "fillntC.h"
#include "vcworkC.h"
#include "mcguts.h"
#endif

#include "spliTChanOutC.h"
#include "fitqunoutC.h"
#include "TRuntimeParameters.hxx"

#include <time.h>

using namespace fiTQun_parameters;
using namespace std;

extern"C" {
  
#ifndef NOSKLIBRARIES
  void fortinit_(char*, char*, int, int);
  int fortread_();
  void skclosef_(int*);
  
  void trginfo_(float *);
  void vcrdvccm_();
  //  void mcguts(int *ntracks);
  
  void aprstbnk_(int *);
  void filldstvars_();
  void readfqzbsbank_(int*);
  void fillfqzbsbank_(int*);
#endif
}

const double pi=3.1415926535898;

int main(int argc, char* argv[]){
  
  bool printStuff = false;
  
  int idum;
  
  int ilist[500];
  int nlist=0;
  bool flglst=false;
  
  float trgofst;
  
  time_t timer;
  time_t timeevt;
  time_t tcurrevt;
  
  double tmpnglnL;
  int PCdumflg;

  TRuntimeParameters::Get().ParseOptions(argc,argv);
  
  // fiTQun parameter file
  std::string fqParOverrideFileName = TRuntimeParameters::Get().GetFQParOverrideFileName();
  if (fqParOverrideFileName.size() >= 1) {
    TRuntimeParameters::Get().ReadParamOverrideFile(fqParOverrideFileName);
  }
  
  // spliTChan parameter file
  std::string scParOverrideFileName = TRuntimeParameters::Get().GetSCParOverrideFileName();
  if (scParOverrideFileName.size() >= 1) {
    TRuntimeParameters::Get().ReadParamOverrideFile(scParOverrideFileName);
  }
  
  std::string inFileName = TRuntimeParameters::Get().GetInFileName();
  std::string outZBSFileName = TRuntimeParameters::Get().GetOutZBSFileName();
  std::string outROOTFileName = TRuntimeParameters::Get().GetOutROOTFileName();
  std::string eventListFileName = TRuntimeParameters::Get().GetEventListFileName();
  
  int flgRingType[nPID];
  flgRingType[ig] = TRuntimeParameters::Get().GetParameterI("fiTQun.DoGamma1RFit");
  flgRingType[ie] = TRuntimeParameters::Get().GetParameterI("fiTQun.DoElectron1RFit");
  flgRingType[imu] = TRuntimeParameters::Get().GetParameterI("fiTQun.DoMuon1RFit");
  flgRingType[ipip] = TRuntimeParameters::Get().GetParameterI("fiTQun.DoPionPlus1RFit");
  flgRingType[ikp] = TRuntimeParameters::Get().GetParameterI("fiTQun.DoKaonPlus1RFit");
  flgRingType[ip] = TRuntimeParameters::Get().GetParameterI("fiTQun.DoProton1RFit");
  flgRingType[icg] = TRuntimeParameters::Get().GetParameterI("fiTQun.DoConeGenerator1RFit");

  bool flg1Rfit = TRuntimeParameters::Get().GetParameterI("fiTQun.DoSingleRingFits");

  if (eventListFileName.size() >= 1) {
    std::cout << "Event list file: " << eventListFileName << std::endl;
    ifstream flist(eventListFileName.c_str());
    while( !flist.eof() ) {
      flist >> ilist[nlist];
      if ((nlist!=0) && !(ilist[nlist]>ilist[nlist-1])) break;
      std::cout << ilist[nlist] << std::endl;
      nlist++;
    }
    flist.close();
    flglst=true;
    std::cout << "Number of events on list: " << nlist << std::endl;
    
  }
  std::cout << "Use event list = " << flglst << std::endl;
  
  int nevents = TRuntimeParameters::Get().GetReadCount();
  std::cout << "Read count = " << nevents << std::endl;
  
  int SkipToEvent = TRuntimeParameters::Get().GetSkipCount();
  int StopatEvent;
  if (SkipToEvent>0) {
    StopatEvent=SkipToEvent+nevents;
  }
  else {
    StopatEvent=nevents;
  }
    
  // Get file name lengths for input into fortinit.
  char* cinFileName = (char*)inFileName.c_str();
  int inFileLength;
  for (int i=0; i<256; i++) {
    if (cinFileName[i]=='\0') {
      inFileLength=i;
      break;
    }
  }
  
  char* coutZBSFileName = (char*)outZBSFileName.c_str();
  int outZBSFileLength;
  for (int i=0; i<256; i++) {
    if (coutZBSFileName[i]=='\0') {
      outZBSFileLength=i;
      break;
    }
  }
  
  int writeHistos = TRuntimeParameters::Get().GetWriteHistos();

  std::cout << "Input file: " << inFileName.c_str() << ";" << inFileLength << std::endl;
  
#ifdef NOSKLIBRARIES
  // setup geometry info, and read first event
  WCSimWrap::Get( (char *)inFileName.c_str() );
  
  // initialize the fit with number of PMTs to use
  fq_NLL_scanner * thefit = new fq_NLL_scanner( WCSimWrap::Get()->NPMT() );
  
  spliTChan * theSubEvents = new spliTChan(writeHistos, WCSimWrap::Get()->NPMT());
  
#else
  
  fortinit_((char*)inFileName.c_str(),(char*)outZBSFileName.c_str(),inFileLength,outZBSFileLength);
  
  std::cout << "Output file: " << outZBSFileName.c_str() << ";" << outZBSFileLength << std::endl;
  
  fq_NLL_scanner *thefit = new fq_NLL_scanner();
  
  // If memory is not dynamically allocated, program crashes for some reason
  spliTChan *theSubEvents = new spliTChan(writeHistos);
  // Runtime parameters loaded within constructor
  
#endif
  
  fiTQun_shared::Get()->SetWriteHistos(writeHistos);
  int useSubEvents = TRuntimeParameters::Get().GetUseSubEvents();
  int copyTrees = TRuntimeParameters::Get().GetCopyTrees();
  std::string copyTreesList = TRuntimeParameters::Get().GetCopyTreesList();
  int doTOFind    =  TRuntimeParameters::Get().GetParameterI("spliTChan.UseTOFMethod");
  std::cout<<"Use TOF Method? "<<doTOFind<<std::endl;
 
  int fQuiet = TRuntimeParameters::Get().GetParameterI("fiTQun.fQuiet");
  fiTQun_shared::Get()->Cantyoubequiet(fQuiet);

  // Which fits should be run?
  
  bool flgSubEvt = TRuntimeParameters::Get().GetParameterI("fiTQun.DoSubEvent");
  
  // Options for all fits
  bool fFitCprof = TRuntimeParameters::Get().GetParameterI("fiTQun.UseFitCProfile");
  if (fFitCprof) std::cout << "Using fit Cherenkov profile!" << std::endl;
  bool fFittprof = TRuntimeParameters::Get().GetParameterI("fiTQun.UseFittProfile");
  if (fFittprof) std::cout << "Using fit time profile!" << std::endl;
  
  int iDetGeom=TRuntimeParameters::Get().GetParameterI("fiTQun.DetGeomType");
  
  thefit->ReadSharedParams(fiTQun_shared::Get()->GetScatflg(),flgRingType,fFitCprof,fFittprof,iDetGeom);
#ifdef NOSKLIBRARIES
  TString strDetName = "WCSim";
#else
  TString strDetName = Form("SK%d",fiTQun_shared::Get()->GetSKVer());
#endif
  
  int iForceFitDefWnd = TRuntimeParameters::Get().GetParameterI("fiTQun.ForceFitDefaultWindow");
  
  bool flgRePrefitTwnd = TRuntimeParameters::Get().GetParameterI("fiTQun.RePrefitTimeWindow");
  int flgInGatefit = TRuntimeParameters::Get().GetParameterI("fiTQun.DoInGateDcyeFit");
    
  if( TRuntimeParameters::Get().HasParameter("fiTQun.CorrTQ") )
    thefit->SetTQCorrflg(TRuntimeParameters::Get().GetParameterI("fiTQun.CorrTQ"));
  
  int flgVtx4Pk = TRuntimeParameters::Get().GetParameterI("fiTQun.PeakSearchVtx");
  
  thefit->SetRingDirTfmFit(TRuntimeParameters::Get().GetParameterI("fiTQun.DoRingDirTfmFit"));
  
  // Read in vertex pre-fit parameters:
  if( TRuntimeParameters::Get().HasParameter("fiTQun.vtxPreFitNfitr") ){
    fiTQun_shared::npfitr = TRuntimeParameters::Get().GetParameterI("fiTQun.vtxPreFitNfitr");
    char sgmarbuff[50];
    if (fiTQun_shared::npfitr > 10) {
      std::cerr << "Error fiTQun.vtxPreFitNfitr must be <= 10." << std::endl;
      exit(-1);
    }
    for (int iNpfitr = 0; iNpfitr < fiTQun_shared::npfitr; iNpfitr++){
      sprintf(sgmarbuff, "fiTQun.vtxPreFitSigmar%i", iNpfitr);
      if( TRuntimeParameters::Get().HasParameter(sgmarbuff) ){
	fiTQun_shared::sgmar[iNpfitr] = TRuntimeParameters::Get().GetParameterD(sgmarbuff);
      } else {
	std::cerr << "Error parameter " << sgmarbuff << " not found. Need " << fiTQun_shared::npfitr << " values of sgmar" << std::endl;
	exit(-1);
      }
    }
  }

  if( TRuntimeParameters::Get().HasParameter("fiTQun.vtxPreFitGrdszVtx") )
    fiTQun_shared::grdszVtx = TRuntimeParameters::Get().GetParameterD("fiTQun.vtxPreFitGrdszVtx");

  if( TRuntimeParameters::Get().HasParameter("fiTQun.vtxPreFitGrdsztime") )
    fiTQun_shared::grdsztime = TRuntimeParameters::Get().GetParameterD("fiTQun.vtxPreFitGrdsztime");
    
  // Options for single ring fits
  int flgtotq[nPID];
  flgtotq[ig] = TRuntimeParameters::Get().GetParameterI("fiTQun.ChargeConstraintGamma1RFit");
  flgtotq[ie] = TRuntimeParameters::Get().GetParameterI("fiTQun.ChargeConstraintElectron1RFit");
  flgtotq[imu] = TRuntimeParameters::Get().GetParameterI("fiTQun.ChargeConstraintMuon1RFit");
  flgtotq[ipip] = TRuntimeParameters::Get().GetParameterI("fiTQun.ChargeConstraintPionPlus1RFit");
  flgtotq[ikp] = TRuntimeParameters::Get().GetParameterI("fiTQun.ChargeConstraintKaonPlus1RFit");
  flgtotq[ip] = TRuntimeParameters::Get().GetParameterI("fiTQun.ChargeConstraintProton1RFit");
  flgtotq[icg] = TRuntimeParameters::Get().GetParameterI("fiTQun.ChargeConstraintConeGenerator1RFit");
  

  thefit->SetQEEffCorr(TRuntimeParameters::Get().GetParameterD(Form("fiTQun.QEEffCorr%s",strDetName.Data())));
    
  
  int flgQTdist = TRuntimeParameters::Get().GetParameterI("fiTQun.OutputQTdist");
  double d_vtx_for_QT = TRuntimeParameters::Get().GetParameterD("fiTQun.dVtxQT");
 
  //Decay E search related
  int doDecayESearch = TRuntimeParameters::Get().GetParameterI("fiTQun.DecayESearch");
  Int_t tmptwl = TRuntimeParameters::Get().GetParameterI("fiTQun.DecayESearchWindowStart");
  Int_t tmptwr = TRuntimeParameters::Get().GetParameterI("fiTQun.DecayESearchWindowEnd");

  TRuntimeParameters::Get().PrintListOfParameters();

  double cqtot = fiTQun_shared::Get()->GetTotqCoef();
  
  // Check if it's an atmospheric MC APFIT file and extract the run/subrun
  string  sinFileName = inFileName.c_str();
  stringstream ssinFileName;
  ssinFileName << inFileName.c_str();
  bool isATMPD = 0;
  bool isMC = 0;
  int mcrun=-1, mcsubrun=-1;
  if (sinFileName.find("atmpd")!=string::npos) isATMPD = 1;
  if (sinFileName.find("mc")!=string::npos) isMC = 1;
#ifndef NOSKLIBRARIES
  if (isATMPD && isMC && sinFileName.find("apfit")!=string::npos) {
    
    vector<string> parsed;
    while(ssinFileName.good()) {
      string substr;
      getline(ssinFileName, substr, '.' );
      parsed.push_back(substr);
    }
    mcrun = atoi(parsed[2].c_str());
    mcsubrun = atoi(parsed[3].c_str());
    
    cout << "Loading atmospheric neutrino MC run " << mcrun << ", subrun " << mcsubrun << endl;
  }
#endif
  
  fiTQun_shared::Get()->Setofilenm(coutZBSFileName);
  
  const int ntestiter=1000;
  
  // ---Variable declarations
  // Single track fiTQun parameters
  double snglTrkParams[fiTQun_shared::nSnglTrkParams];
  
  enum peakFindingMethods {fqpeak_fqpftof, fqpeak_fq1retof, muechk_notof, muechk_fqpftof, muechk_fq1retof, muechk_aptof};
  
  int nevt=0;
  
  const int nEvisRang = 6;
  Double_t EvisRang[nEvisRang+1] = {0,300,700,1330,2500,5000,15000};
  TH2D *hTQdist[nEvisRang][3] = {};
  
  fQROOTOut *theROOTOut = NULL;
  if (outROOTFileName.size() >= 1) {
    theROOTOut = new fQROOTOut((char*)outROOTFileName.c_str());
      if(copyTrees) theROOTOut->CopyInput( (char *)inFileName.c_str(), copyTreesList );
    theROOTOut->GetTree()->Branch("nevt",&nevt,"nevt/I");
    
#ifndef NOSKLIBRARIES
    
    theROOTOut->GetTree()->Branch("nrun",&apdstnt_.nrun,"nrun/I");
    theROOTOut->GetTree()->Branch("nsub",&apdstnt_.nsub,"nsub/I");
    theROOTOut->GetTree()->Branch("nev",&apdstnt_.nev,"nev/I");
    theROOTOut->GetTree()->Branch("date",apdstnt_.date,"date[3]/I");
    theROOTOut->GetTree()->Branch("time",apdstnt_.time,"time[4]/I");
    
    theROOTOut->GetTree()->Branch("nhit",&apdstnt_.nhit,"nhit/I");
    theROOTOut->GetTree()->Branch("potot",&apdstnt_.potot,"potot/F");
    
    theROOTOut->GetTree()->Branch("nhitac",&apdstnt_.nhitac,"nhitac/I");
    theROOTOut->GetTree()->Branch("wall",&apdstnt_.wall,"wall/F");
    theROOTOut->GetTree()->Branch("evis",&apdstnt_.evis,"evis/F");
    theROOTOut->GetTree()->Branch("nring",&apdstnt_.nring,"nring/I");
    theROOTOut->GetTree()->Branch("ip",apdstnt_.ip,"ip[nring]/I");
    theROOTOut->GetTree()->Branch("pos",apdstnt_.pos,"pos[3]/F");
    theROOTOut->GetTree()->Branch("dir",apdstnt_.dir,"dir[nring][3]/F");
    theROOTOut->GetTree()->Branch("amome",apdstnt_.amome,"amome[nring]/F");
    theROOTOut->GetTree()->Branch("amomm",apdstnt_.amomm,"amomm[nring]/F");
    
    theROOTOut->GetTree()->Branch("nmue",&apntmue_.nmue,"nmue/I");
    theROOTOut->GetTree()->Branch("etype",apntmue_.etype,"etype[nmue]/I");
    theROOTOut->GetTree()->Branch("etime",apntmue_.etime,"etime[nmue]/F");
    theROOTOut->GetTree()->Branch("egood",apntmue_.egood,"egood[nmue]/F");
    theROOTOut->GetTree()->Branch("ehit",apntmue_.ehit,"ehit[nmue]/F");
    
#endif
    
    if (flgQTdist) {
      if (flgQTdist>0) thefit->SetTreeQTdist(theROOTOut->GetTree());
      for (int iErang=0; iErang<nEvisRang; iErang++) {
        for (int iFV=0; iFV<3; iFV++) {
          hTQdist[iErang][iFV] = new TH2D(Form("hTQdist_Erang%d_FV%d",iErang,iFV),"",180,0.,180.,80,-20.,20.);
        }
      }
    }
    
  }
  
  timer = time(NULL);
  
  int nlstctr=0;

  while (1) {
    
    for (int idx=0; idx<20; idx++) {
      fitquninfo_.fqproctime[idx]=0.;
    }
    
    tcurrevt = time(NULL);
#ifndef NOSKLIBRARIES
    
    // Read next ZBS event
    if (nevt>0) {// first event is read in fortinit
      int iret = fortread_();
      if (iret>0) break;
    }
    
    std::cout << std::endl;
    
    std::cout << "Processing ZBS event # " << nevt << std::endl;
    
#else
    
    WCSimWrap * wc = WCSimWrap::Get();
    int iret = wc->LoadEntry( nevt );
    int curevnum = wc->SubEvt(0)->GetHeader()->GetEvtNum();
    
    std::cout <<"================================================="<< std::endl;
    std::cout <<" Load Entry "<<nevt<<" of "<<wc->NEvt()<<" EvtNum "<<curevnum<<" retval="<<iret<<std::endl;
    
    if (iret<0) break;
#endif
    
    if (flglst) {
      if (nlstctr>=nlist) break;
#ifndef NOSKLIBRARIES
      if ( (isATMPD && ilist[nlstctr]==skhead_.nevsk) || (ilist[nlstctr]==nevt) ) {
        nlstctr++;
      } else {
        nevt++;
        continue;
      }
#else
      if ( (isATMPD && ilist[nlstctr]==curevnum) ||
          (ilist[nlstctr]==nevt) ) {
        nlstctr++;
      } else {
        nevt++;
        continue;
      }
#endif
    }
    
    if (SkipToEvent>0) {
      if (nevt<SkipToEvent) {
        nevt++;
        continue;
      }
    }
    if (StopatEvent>0) {
      if (nevt>=StopatEvent) break;
    }
    
#ifndef NOSKLIBRARIES
    /**************Read information from zbs banks************************/
    int inputbnk=0;
    aprstbnk_(&inputbnk);
    filldstvars_();
    
    readfqzbsbank_(&idum);
    
    /* load mc_particle structure */
    int nmctrks;
    //    mcguts(&nmctrks);
    
    vcrdvccm_();
    trginfo_(&trgofst);
    
    /**************End reading zbs banks**********************************/
#endif
    
    thefit->InitEvent(nevt);
    
    /**************Subevent Section***************************************/
    while(1) {
      
      if (!flgSubEvt) break;
      
      thefit->ClrInfo(0);
   
      thefit->SetTimeWindow(0);

      // fiTQun pre-fit vertex on default time window
      double vtx_prefit[4] = {0};
      thefit->DoVtxPrefit(vtx_prefit);
      
      // Do fiTQun 1R e-hypothesis fit for better vertex
      if (flgRingType[ie] && flgVtx4Pk!=0) {
        if (flgtotq[ie]) fiTQun_shared::Get()->SetTotqCoef(cqtot);
        else fiTQun_shared::Get()->SetTotqCoef(0.);
        tmpnglnL=thefit->Do1RFit(ie,snglTrkParams,PCdumflg,1);
        
        // TOF corrected fiTQun/spliTChan peak search with FQ 1Re vertex
        thefit->DoPeakSearch(snglTrkParams);
        theSubEvents->CopyPeakInfo(fqpeak_fq1retof,fitquntwnd_.fqtwnd_npeak[0],fitquntwnd_.fqtwnd_peakt0[0], fitquntwnd_.fqtwnd_peakiness[0]);
      }
      
      // TOF corrected fiTQun/spliTChan peak search with FQ prefit vertex
      thefit->DoPeakSearch(vtx_prefit);
      theSubEvents->CopyPeakInfo(fqpeak_fqpftof,fitquntwnd_.fqtwnd_npeak[0],fitquntwnd_.fqtwnd_peakt0[0], fitquntwnd_.fqtwnd_peakiness[0]);
      
      // Search for clusters with simple scanning algorithm.
      // This fills the Cluster_* variables, which is used for shifting the time
      // window below with "FillTQbanks".
      theSubEvents->Cluster_Search(0);
      
      // Search for decay-e peaks based on original APFIT::Muechk algorithm.
      theSubEvents->MuechkWrap(0,muechk_notof); // No TOF correction
      
      // Use default SW trigger window for prefit and 1Re fit below
      // to mimic behaviour of original Muechk
      
      // TOF corrected Muechk with FQ prefit vertex
      theSubEvents->MuechkWrap(vtx_prefit,muechk_fqpftof);
      
      // TOF corrected Muechk with FQ 1Re vertex
      if (flgRingType[ie] && flgVtx4Pk!=0)
        theSubEvents->MuechkWrap(snglTrkParams, muechk_fq1retof);
      
      // Count the number of peaks found by Muechk in each cluster
      theSubEvents->Cluster_AssociateMuechk();
      cout<<"find clusters using peaks..?";
      if (doTOFind) cout<<"YES!"<<endl;
      if (doTOFind){
        theSubEvents->MakeClusters(fitquntwnd_.fqtwnd_npeak[0],
                                   fitquntwnd_.fqtwnd_peakt0[0],vtx_prefit,nevt);
      }
      
#ifndef NOSKLIBRARIES
      if (fiTQun_shared::Get()->GetSKVer()<=3) theSubEvents->SetATMClusters();
#endif
      
      for (int icluster=0; icluster<theSubEvents->Cluster_nCand; icluster++) {
        fitqunse_.cluster_npeaks[icluster][fqpeak_fqpftof]=0;
      }
      
      // Associate peaks to clusters
      for (int ipeak=0; ipeak<fitquntwnd_.fqtwnd_npeak[0]; ipeak++) {
        vtx_prefit[3]=fitquntwnd_.fqtwnd_peakt0[0][ipeak];
        for (int icluster=0; icluster<theSubEvents->Cluster_nCand; icluster++) {
          
          if (thefit->IsPeakInCluster(vtx_prefit,theSubEvents->Cluster_tStart[icluster],
                                      theSubEvents->Cluster_tEnd[icluster])){//is the peak caused by hits in this cluster?
            
            fitqunse_.muechk_icluster[fqpeak_fqpftof][ipeak] = icluster;
            
            fitqunse_.cluster_ipeak[icluster][fqpeak_fqpftof][fitqunse_.cluster_npeaks[icluster][fqpeak_fqpftof]] = ipeak;
            fitqunse_.cluster_timeofpeak[icluster][fqpeak_fqpftof][fitqunse_.cluster_npeaks[icluster][fqpeak_fqpftof]]=fitquntwnd_.fqtwnd_peakt0[0][ipeak];
            fitqunse_.cluster_npeaks[icluster][fqpeak_fqpftof]++;
            
            break;
          }
        }
      }
      
      // Tag bad clusters as those without a corresponding FQ peak
      if (flgRingType[ie] && flgVtx4Pk==2) theSubEvents->RemoveNoPeakClusters(fqpeak_fq1retof);
      else theSubEvents->RemoveNoPeakClusters(fqpeak_fqpftof);
      
      fiTQun_shared::Get()->ClearPeaks(0);
      
      break;
    }
    /**********End Subevent Section***************************************/
    
    fitquninfo_.fqproctime[1] = time(NULL)-tcurrevt;
    
    fitquntwnd_.trgoff=trgofst;
    
    // Sometimes cluster is not found, or does not exist,
    // but run fiTQun on original TQREAL window anyways
    int ncluster_min = 0;
    if (iForceFitDefWnd) ncluster_min = 1;
    int ncluster = max(ncluster_min,fitqunse_.cluster_ncand);
    
    int icluster_good=0;
    for (int icluster=0; icluster<ncluster; icluster++) {
      
      bool noGoodClusters = 0;
      if (!fitqunse_.cluster_goodflag[icluster]) {
        if (!theSubEvents->fQuiet) cout << "Skipping bad cluster " << icluster << endl;
        
        if (icluster==ncluster-1 && icluster_good==0)
          noGoodClusters = 1;
        else
          continue;
      }
      //      if (!(fitqunse_.cluster_npeaks[icluster][0]>1)) continue;//comment this out!!!!!!!!!!!!!!!!!!!!!!!!!
      bool flgCluster=true;
      // Modify time window and offset
      if (useSubEvents && fitqunse_.cluster_ncand && noGoodClusters==0) {
      	
      	// Shift fiTQun time window and offset
      	thefit->SetTimeWindow(icluster_good,fitqunse_.cluster_tstart[icluster],
                              fitqunse_.cluster_tend[icluster]);
      } 
      else {
        flgCluster=false;
      	thefit->SetTimeWindow(icluster_good);
      }
      
      if (flgSubEvt) {
        if (flgCluster) {
          fitquntwnd_.fqtwnd_iclstr[icluster_good]=icluster;
          thefit->DoIngateFit(flgInGatefit,flgRePrefitTwnd,fitqunse_.cluster_npeaks[icluster][fqpeak_fqpftof],fitqunse_.cluster_timeofpeak[icluster][fqpeak_fqpftof]);
        }
        else {
          fitquntwnd_.fqtwnd_iclstr[icluster_good]=-1;
          thefit->DoIngateFit(flgInGatefit,flgRePrefitTwnd,0,NULL);
        }
      }
      else {
        thefit->MaskIngates();
      }
      
      thefit->ScanNLL(imu,snglTrkParams);
      break;
    }
    nevt++;
  }
}

