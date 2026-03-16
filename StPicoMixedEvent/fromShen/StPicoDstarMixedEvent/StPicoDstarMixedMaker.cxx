/* **************************************************
 *
 *  Authors: Kaifeng Shen
 *           Yuanjing Ji
 *           Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************/


#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

#include "TFile.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"

#include "StEvent/StDcaGeometry.h"
#include "StPhysicalHelixD.hh"
#include "phys_constants.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoDstarMixedMaker.h"
#include "StAnaCuts.h"
#include "StMemStat.h"
//#include "globalconstant.h"
//#include "calmean.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "StLorentzVector.hh"
//#include "../StRefMultCorr/StRefMultCorr.h"
//#include "../StRefMultCorr/CentralityMaker.h"
// #include "TVector3.h"
#ifndef DEBUG
#define DEBUG 1
#endif

// Mixed Event
//mMagBins=2; mCenBins=9; mVzBins=20; mEveBins=24; mMaxEventsInBuffer=100;
float bfield=0;
const Int_t mMagBins = 2; //1;
//const Int_t mCenBins = 16;
const Int_t mCenBins = 17;
const Int_t mVzBins = 20; //10;
//const Int_t mEveBins = 24; //12;
const Int_t mMaxEventsInBuffer = 50; //50;
const Int_t mMaxElectrons = 30;
//const Float_t mPhiVCutMRange = 0.2;
//Float_t current_EQx[mMaxElectrons],current_EQy[mMaxElectrons];
Int_t current_nE,current_nEPlus,current_nEMinus;
TLorentzVector current_ePlus[mMaxElectrons];
TLorentzVector current_eMinus[mMaxElectrons];
// Float_t current_p_E0_ePlus[mMaxElectrons];
// Float_t current_p_E0_eMinus[mMaxElectrons];
// Float_t current_p_E_ePlus[mMaxElectrons];
// Float_t current_p_E_eMinus[mMaxElectrons];
// Float_t current_bemcdz_ePlus[mMaxElectrons];
// Float_t current_bemcdz_eMinus[mMaxElectrons];
// Float_t current_bemcdphi_ePlus[mMaxElectrons];
// Float_t current_bemcdphi_eMinus[mMaxElectrons];
// Float_t current_bemcTowerPhi_ePlus[mMaxElectrons];
// Float_t current_bemcTowerPhi_eMinus[mMaxElectrons];
// Float_t current_bemcTowerEta_ePlus[mMaxElectrons];
// Float_t current_bemcTowerEta_eMinus[mMaxElectrons];
// Float_t current_nSMDphi_ePlus[mMaxElectrons];
// Float_t current_nSMDphi_eMinus[mMaxElectrons];
// Float_t current_nSMDeta_ePlus[mMaxElectrons];
// Float_t current_nSMDeta_eMinus[mMaxElectrons];
Float_t current_MomMag_ePlus[mMaxElectrons];
Float_t current_MomMag_eMinus[mMaxElectrons];
Float_t current_MomPerp_ePlus[mMaxElectrons];
Float_t current_MomPerp_eMinus[mMaxElectrons];
Float_t current_nSigE_eMinus[mMaxElectrons];
Float_t current_nSigE_ePlus[mMaxElectrons];
Bool_t current_phePlusTag[mMaxElectrons];
Bool_t current_pheMinusTag[mMaxElectrons];
Int_t current_bemcId_ePlus[mMaxElectrons];
Int_t current_bemcId_eMinus[mMaxElectrons];

StThreeVectorF current_gMom_ePlus[mMaxElectrons];
StThreeVectorF current_gMom_eMinus[mMaxElectrons];
StThreeVectorF current_origin_ePlus[mMaxElectrons];
StThreeVectorF current_origin_eMinus[mMaxElectrons];



Int_t magBufferPointer, cenBufferPointer, vzBufferPointer;
Int_t nEventsInBuffer[mMagBins][mCenBins][mVzBins];
Bool_t bufferFullFlag[mMagBins][mCenBins][mVzBins];
Int_t buffer_nEPlus[mMagBins][mCenBins][mVzBins][mMaxEventsInBuffer];
Int_t buffer_nEMinus[mMagBins][mCenBins][mVzBins][mMaxEventsInBuffer];
TLorentzVector buffer_ePlus[mMagBins][mCenBins][mVzBins][mMaxEventsInBuffer][mMaxElectrons];
TLorentzVector buffer_eMinus[mMagBins][mCenBins][mVzBins][mMaxEventsInBuffer][mMaxElectrons];
TRandom3* myRandom;
TRandom3* myRandom_phi;
StRefMultCorr* refmultCorrUtil;

Int_t bins[4] = {40,10,40,17};
Double_t xmin[4] = {0,0,-1.0*TMath::Pi(),-1};
Double_t xmax[4] = {4,1,TMath::Pi(),16};
Double_t fill_4[5] = {0,0,0,0,0};//mass,pT,delta phi, centrality, weight

Int_t bins_phi[4] = {10,40,40,17};
Double_t xmin_phi[4] = {0,-1.0*TMath::Pi(),-1.0*TMath::Pi(),-1};
Double_t xmax_phi[4] = {1,1.0*TMath::Pi(),1.0*TMath::Pi(),16};
Double_t fill_4_phi[5] = {0,0,0,0,0};//pT, phi1, phi2, centrality, weight


//nSigmaE corrected
const Int_t nCenBins = 8;
Int_t mCenBinLow[nCenBins]      = {0, 2, 4, 6, 8, 10, 12,14};
Int_t mCenBinHi[nCenBins]       = {1, 3, 5, 7, 9, 11, 13,15};
Int_t CentralityLow[nCenBins] = {70,60, 50, 40, 30, 20, 10, 0};//%
Int_t CentralityHi[nCenBins]  = {80,70, 60, 50, 40, 30, 20, 10};//20220801 eid

const Int_t sys_index = 2;//0->Ru;1->Zr

TF1 * f_PEPlusnSigmaE_eta_2D[sys_index][nCenBins];
TF1 * f_PEMinusnSigmaE_eta_2D[sys_index][nCenBins];

TProfile * pPEPlusnSigmaE_eta_2D[sys_index][nCenBins];
TProfile * pPEMinusnSigmaE_eta_2D[sys_index][nCenBins];

TProfile * pPEPlusnSigmaE_phi_2D[sys_index][nCenBins];
TProfile * pPEMinusnSigmaE_phi_2D[sys_index][nCenBins];


TString sys_name[sys_index] = {"Ru","Zr"};


TFile * file_nSigmaE;

//phi V cut
TF1 *phiVcut;
const Float_t mPhiVCutMRange = 0.2;






ClassImp(StPicoDstarMixedMaker)
  StPicoDstarMixedMaker::StPicoDstarMixedMaker(char const * name, TString const inputFilesList, TString const outFileBaseName, StPicoDstMaker* picoDstMaker):
    StMaker(name), mPicoDstMaker(picoDstMaker),
    mInputFilesList(inputFilesList), mOutFileBaseName(outFileBaseName),mRunbyRunQA(false)
{}

Int_t StPicoDstarMixedMaker::Init()
{
  mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");
  // -------------- USER VARIABLES -------------------------
  //mFile = new TFile(mOutFileBaseName+".QA.root", "RECREATE");
  mFile_analysis = new TFile(mOutFileBaseName+".analysis.root","RECREATE");
  //initialize trees
  initHists();
 //cout<<"end Init"<<endl;
  return kStOK;
}
//-----------------------------------------------------------------------------
StPicoDstarMixedMaker::~StPicoDstarMixedMaker()
{}
//-----------------------------------------------------------------------------
void StPicoDstarMixedMaker::initHists(){
  //int totalNum = 614;//54.4 GeV
  int totalNum = 40;//200 GeV 2019 auau

  char name_RunID[100];

    ifstream readnum;
    readnum.open(mRunNumList);
    int tmp=0;
    if (DEBUG) cout<<"start initial run number..."<<endl;
    for (int i=0;i<totalNum;i++){
      readnum>>tmp;
      runnum.insert(pair<int,int>(tmp,i));

    /*   sprintf(name_RunID,"hinvBetavsP_%d",tmp);
       hinvBetavsP_RunID[i] = new TH2F(name_RunID,"#frac{1}{#beta} vs P;P(GeV/c);#frac{1}{#beta}",300,0,3,200,0.5,1.5);
       memset(name_RunID,0,sizeof(name_RunID));*/

      if (DEBUG) cout <<"run number : " << tmp <<" id :" <<runnum[tmp] <<endl;
    }
    readnum.close();

  myRandom = new TRandom3();
  myRandom_phi = new TRandom3();
  memset(nEventsInBuffer,0,sizeof(nEventsInBuffer));
  memset(bufferFullFlag,0,sizeof(bufferFullFlag));
  memset(buffer_nEPlus,0,sizeof(buffer_nEPlus));
  memset(buffer_nEMinus,0,sizeof(buffer_nEMinus));


  //phi V
  phiVcut=new TF1("phiVcut","0.84326*exp((-49.4819)*x)+(-0.996609)*x+(0.19801)",0.,1.0); 

  // event level QA
  hevt = new TH1D("hevt","hevt",totalNum,0,totalNum);
  hevt_Cent = new TH1D("hevt_Cent","hevt_Cent",16,-0.5,15.5);
  hevt_Cent->Sumw2();
  
  hpassevtcut = new TH1D("hpassevtcut","pass event cut",6  , -0.5 , 5.5 );
 

    hVz = new TH1F("hVz","Vz;Vz(cm);Counts",240,-60,60);
    hVz_total = new TH1F("hVz_total","Vz without event cut;Vz(cm);Counts",400,-100,100);
    hVpdVz = new TH1F("hVpdVz","VpdVz;VpdVz(cm);Counts",240,-60,60);
    hVr = new TH1F("hVr","Vr;Vr(cm);Counts",100,0,1);
    hVzVpdVz = new TH1F("hVzVpdVz","Vz-VpdVz(cm)",120,-8,8);

    hrefmult = new TH1F("hrefmult","hrefmult",700,0,700);

   //tracl level QA
    hnHitsFit = new TH1F("hnHitsFit","nHitsFit;nHitsFit",90,0,90);
    // hpDca = new TH1F("hpDca","pDca",50,0,5);
    hgDca = new TH1F("hgDca","gDca",50,0,5);
    hinvBetavsP = new TH2F("hinvBetavsP","#frac{1}{#beta} vs p;p(GeV/c);#frac{1}{#beta}",600,0,3,1200,0.5,5.5);
    // hinvBetavsY;
    hdEdx = new TH2F("hdEdx","dEdx vs p*charge;p*charge(GeV/c);dEdx",200,-2,2,400,0,25);
    //hNsigEvsinvBeta = new TH3F("hNsigEvsinvBeta","nSigmaE vs 1/#beta;nSigmaE;1/#beta;p",200,-10,10,100,0.8,1.2,100,0.15,2.5);
    hpt = new TH1F("hpt","hpt;p_{T}(GeV/c)",240,0,12);
    hGpt = new TH1F("hGpt","hGpt;global p_{T}(GeV/c)",240,0,12);
    hEta = new TH1F("hEta","Eta;#eta",60,-1.5,1.5);
    hPhi = new TH1F("hPhi","Phi;#phi",80,-4,4);


    //invariant mass electron
    h_electron_count = new TH1F("h_electron_count","electron counts;electron(1) positron(2)",2,0.5,2.5);
    hMeeCount = new TH1F("hMee","hMee;Mee(GeV/c^{2})",400000,0,4);
    hMeeCount_like1 = new TH1F("hMee_like1","hMee like sign electron;Mee(GeV/c^{2})",400000,0,4);
    hMeeCount_like2 = new TH1F("hMee_like2","hMee like sign positron;Mee(GeV/c^{2})",400000,0,4);

    hMeeCount_wo_PE = new TH1F("hMee_wo_PE","hMee;Mee(GeV/c^{2})",400000,0,4);
    hMeeCount_like1_wo_PE = new TH1F("hMee_like1_wo_PE","hMee like sign electron;Mee(GeV/c^{2})",400000,0,4);
    hMeeCount_like2_wo_PE = new TH1F("hMee_like2_wo_PE","hMee like sign positron;Mee(GeV/c^{2})",400000,0,4);

    hMeeCount_helix = new TH1F("hMee_helix","hMee;Mee(GeV/c^{2})",400000,0,4);
    hMeeCount_like1_helix = new TH1F("hMee_like1_helix","hMee like sign electron;Mee(GeV/c^{2})",400000,0,4);
    hMeeCount_like2_helix = new TH1F("hMee_like2_helix","hMee like sign positron;Mee(GeV/c^{2})",400000,0,4);

    hMeeCount_helix_PE = new TH1F("hMee_helix_PE","hMee;Mee(GeV/c^{2})",400000,0,4);
    hMeeCount_like1_helix_PE = new TH1F("hMee_like1_helix_PE","hMee like sign electron;Mee(GeV/c^{2})",400000,0,4);
    hMeeCount_like2_helix_PE = new TH1F("hMee_like2_helix_PE","hMee like sign positron;Mee(GeV/c^{2})",400000,0,4);


    hMee_Pt_Cent = new TH3F("hMee_Pt_Cen","hMee vs Pt vs cent;P_{t};Cent;Mee",100,0,10,17,-1,16,400,0,4);
    hMeelike1_Pt_Cent = new TH3F("hMeelike1_Pt_Cen","hMee like sign electron vs Pt vs cent;P_{t};Cent;Mee",100,0,10,17,-1,16,400,0,4);
    hMeelike2_Pt_Cent = new TH3F("hMeelike2_Pt_Cen","hMee like sign positron vs Pt vs cent;P_{t};Cent;Mee",100,0,10,17,-1,16,400,0,4);
    hMeeMix_Pt_Cent = new TH3F("hMeeMix_Pt_Cen","hMeeMix vs Pt vs cent;P_{t};Cent;Mee",100,0,10,17,-1,16,400,0,4);
    hMeeMixlike1_Pt_Cent = new TH3F("hMeeMixlike1_Pt_Cen","hMeeMix like sign electron vs Pt vs cent;P_{t};Cent;Mee",100,0,10,17,-1,16,400,0,4);
    hMeeMixlike2_Pt_Cent = new TH3F("hMeeMixlike2_Pt_Cen","hMeeMix like sign positron vs Pt vs cent;P_{t};Cent;Mee",100,0,10,17,-1,16,400,0,4);//for low pT

    hMee_Pt_Cent->Sumw2();
    hMeelike1_Pt_Cent->Sumw2();
    hMeelike2_Pt_Cent->Sumw2();
    hMeeMix_Pt_Cent->Sumw2();
    hMeeMixlike1_Pt_Cent->Sumw2();
    hMeeMixlike2_Pt_Cent->Sumw2();
   
   // hMee_Pt_80_100 = new TH2F("hMee_Pt_80_100","hMee vs Pt vs 80;P_{t};Cent;Mee",1000,0,10,400,0,4);
   // hMeelike1_Pt_80_100 = new TH2F("hMeelike1_Pt_80_100","hMee like sign electron vs Pt vs 80;P_{t};Cent;Mee",1000,0,10,400,0,4);
   // hMeelike2_Pt_80_100 = new TH2F("hMeelike2_Pt_80_100","hMee like sign positron vs Pt vs 80;P_{t};Cent;Mee",1000,0,10,400,0,4);
   // hMeeMix_Pt_80_100 = new TH2F("hMeeMix_Pt_80_100","hMeeMix vs Pt vs 80;P_{t};Cent;Mee",1000,0,10,400,0,4);
   // hMeeMixlike1_Pt_80_100 = new TH2F("hMeeMixlike1_Pt_80_100","hMeeMix like sign electron vs Pt vs 80;P_{t};Cent;Mee",1000,0,10,400,0,4);
   // hMeeMixlike2_Pt_80_100 = new TH2F("hMeeMixlike2_Pt_80_100","hMeeMix like sign positron vs Pt vs 80;P_{t};Cent;Mee",1000,0,10,400,0,4);//for low pT

    //THnSparseF hMee_Pt_phi_Cent("hMee_Pt_phi_Cen","hMee vs Pt vs phi vs cent",4,bins;xmin,xmax);
    hMee_Pt_phi_Cent = new THnF("hMee_Pt_phi_Cen","hMee vs Pt vs phi vs cent;Mass;P_{t};#Delta#phi;Cent",4,bins,xmin,xmax);
    //hMee_Pt_phi1_phi2_Cent = new THnF("hMee_Pt_phi1_phi2_Cen","hMee vs Pt vs phi1 vs phi2 vs cent;P_{t};#phi_{1};#phi_{2};Cent",4,bins_phi,xmin_phi,xmax_phi);
    hMeelike1_Pt_phi_Cent = new THnF("hMeelike1_Pt_phi_Cen","hMee like sign electron vs Pt vs phi vs cent;Mass;P_{t};#Delta#phi;Cent",4,bins,xmin,xmax);
    //hMeelike1_Pt_phi1_phi2_Cent = new THnF("hMeelike1_Pt_phi1_phi2_Cen","hMee like sign electron vs Pt vs phi1 vs phi2 vs cent;P_{t};#phi_{1};#phi_{2};Cent",4,bins_phi,xmin_phi,xmax_phi);
    hMeelike2_Pt_phi_Cent = new THnF("hMeelike2_Pt_phi_Cen","hMee like sign positive vs Pt vs phi vs cent;Mass;P_{t};#Delta#phi;Cent",4,bins,xmin,xmax);
    hMeeMix_Pt_phi_Cent = new THnF("hMeeMix_Pt_phi_Cen","hMeeMix vs Pt vs phi vs cent;Mass;P_{t};#Delta#phi;Cent",4,bins,xmin,xmax);
    hMeeMixlike1_Pt_phi_Cent = new THnF("hMeeMixlike1_Pt_phi_Cen","hMeeMix like sign electron vs Pt vs phi vs cent;Mass;P_{t};#Delta#phi;Cent",4,bins,xmin,xmax);
    hMeeMixlike2_Pt_phi_Cent = new THnF("hMeeMixlike2_Pt_phi_Cen","hMeeMix like sign positron vs Pt vs phi vs cent;Mass;P_{t};#Delta#phi;Cent",4,bins,xmin,xmax);//for low pT, delta phi
    //hMee_Pt_phi_Cent->Sumw2();
    //hMeelike1_Pt_phi_Cent->Sumw2();
    //hMeelike2_Pt_phi_Cent->Sumw2();
    //hMeeMix_Pt_phi_Cent->Sumw2();
    //hMeeMixlike1_Pt_phi_Cent->Sumw2();
    //hMeeMixlike2_Pt_phi_Cent->Sumw2();

    h_nSigmaElectron_tof_p = new TH2F("h_nSigmaElectron_tof_p"," ;p (GeV/c);n#sigma_{e}",40,0,4,200,-10,10);
    h_nSigmaElectron_tof_pT = new TH2F("h_nSigmaElectron_tof_pT"," ;p_{T} (GeV/c);n#sigma_{e}",40,0,4,200,-10,10);
   
    

	h_Cen_Vz = new TH2F("h_Cen_Vz"," ;Centrality;Vz (cm)",16,0,16,600,-35,25);
	h_Cen_reweight = new TH2F("h_Cen_reweight"," ;Centrality;reweight(Cent)",16,0,16,1000,-5,5);
	//h_Cen_reweight->Sumw2();
	

	
}
//-----------------------------------------------------------------------------
Int_t StPicoDstarMixedMaker::Finish()
{
    cout<<"finish"<<endl;
    mFile_analysis->cd();
  //write the hists
 //event QA
//    hVxVyVz->Write();
    hVz->Write();
    hVz_total->Write();
    hVpdVz->Write();
    hVr->Write();
    hVzVpdVz->Write();
   // hnEvsEtavsVz->Write();
   // hnEvsPhivsVz->Write();
   // hnTofMulvsRef->Write();
    //hnTofMatvsRef->Write();
   // hnTofHitvsRef->Write();

    hevt->Write();
    hevt_Cent->Write();

    hpassevtcut->Write();
    hrefmult->Write();
   // hNsigEvsinvBeta->Write();
   //tracl level QA
    hnHitsFit->Write();
    // hpDca->Write();
    hgDca->Write();
    hinvBetavsP->Write();
    // hinvBetavsY->Write();
    hdEdx->Write();
    hpt->Write();
    hGpt->Write();
    hEta->Write();
    hPhi->Write();
   
	h_Cen_Vz->Write();
	h_Cen_reweight->Write();
	

  hMeeCount->Write();
  hMeeCount_like1->Write();
  hMeeCount_like2->Write();

   hMeeCount_wo_PE->Write();
   hMeeCount_like1_wo_PE->Write();
   hMeeCount_like2_wo_PE->Write();

   hMeeCount_helix->Write();
   hMeeCount_like1_helix->Write();
   hMeeCount_like2_helix->Write();

   hMeeCount_helix_PE->Write();
   hMeeCount_like1_helix_PE->Write();
   hMeeCount_like2_helix_PE->Write();



  hMee_Pt_Cent->Write();
  hMeelike1_Pt_Cent->Write();
  hMeelike2_Pt_Cent->Write();
  hMeeMix_Pt_Cent->Write();
  hMeeMixlike1_Pt_Cent->Write();
  hMeeMixlike2_Pt_Cent->Write();
  hMee_Pt_phi_Cent->Write();
  hMeelike1_Pt_phi_Cent->Write();
  hMeelike2_Pt_phi_Cent->Write();
  hMeeMix_Pt_phi_Cent->Write();
  hMeeMixlike1_Pt_phi_Cent->Write();
  hMeeMixlike2_Pt_phi_Cent->Write();

  //hMee_Pt_phi1_phi2_Cent->Write();
 // hMeelike1_Pt_phi1_phi2_Cent->Write();

  h_nSigmaElectron_tof_p->Write();
  h_nSigmaElectron_tof_pT->Write();


  mFile_analysis->Close();


  return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoDstarMixedMaker::Make()
{
  // StMemStat mem;
  //cout<<"star make()"<<endl;
  ParticleInfo_Electron particleinfo;
  vector<ParticleInfo_Electron> electroninfo;
  vector<ParticleInfo_Electron> positroninfo;
  iEvent_tof=iEvent_tof+1;


  if (!mPicoDstMaker)
  {
    LOG_WARN << " StPicoDstarMixedMaker - No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }
  StPicoDst const* picoDst = mPicoDstMaker->picoDst();
  if (!picoDst)
  {
    LOG_WARN << "StPicoDstarMixedMaker - No PicoDst! Skip! " << endm;
    return kStWarn;
  }
  // -------------- USER ANALYSIS -------------------------
  StPicoEvent const * picoEvent = picoDst->event();
  //trigger
  hpassevtcut->Fill(0);

  if (!isGoodTrigger(picoEvent)) return 0;

  hpassevtcut->Fill(1);
  mRunId = picoEvent->runId();
  //cout<<"run number(good trigger): "<<mRunId<<endl;
  hevt->Fill(runnum[mRunId]);


  


  TVector3 pVtx = picoEvent->primaryVertex();
  bool vzcut = fabs(pVtx.z()) < 60;
  bool verrcut = !(fabs(pVtx.x()) < anaCuts::Verror && fabs(pVtx.y()) < anaCuts::Verror &&
        fabs(pVtx.z()) < anaCuts::Verror);
  bool vrcut =  sqrt(TMath::Power(pVtx.x(), 2) + TMath::Power(pVtx.y(), 2)) <=  anaCuts::Vr ;
  // bool vpdvzcut = fabs(pVtx.z() - picoEvent->vzVpd()) < 3;
  bool vpdvzcut = true;
  //if (vzcut) hpassevtcut->Fill(2);
  //if (vzcut &&  vrcut) hpassevtcut->Fill(3);
  // if (vzcut && vrcut  &&  vpdvzcut ) hpassevtcut->Fill(4);
  //if (vzcut && vrcut  &&  vpdvzcut && verrcut ) hpassevtcut->Fill(4);
  bool refusepileup = picoEvent->refMult()<picoEvent->btofTrayMultiplicity()*0.36+45;
  bool refusebadtof = picoEvent->refMult()>picoEvent->btofTrayMultiplicity()*0.28-115;
  bool passCentralityCut = refusepileup && refusebadtof  && verrcut && vrcut && fabs(pVtx.z()) < 10;
  if (passCentralityCut) hrefmult->Fill(picoEvent->refMult());

if(!isBadrun(mRunId)){
    //hevt_bad_cut->Fill(runnum[mRunId]);
	hpassevtcut->Fill(2);
    // StThreeVectorF pVtx = picoEvent->primaryVertex();
   // cout<<"good run"<<endl;
    TVector3 pVtx = picoEvent->primaryVertex();
    mVx = pVtx.x();
    mVy = pVtx.y();
    mVz = pVtx.z();
    //cout<<"vz: "<<mVz<<endl;
    hVz_total->Fill(mVz);
 if (isGoodEvent(picoEvent)){
    //cout<<"star event"<<endl;
	hpassevtcut->Fill(3);

    // StThreeVectorF pVtx = picoEvent->primaryVertex();
    TVector3 pVtx = picoEvent->primaryVertex();
    mVx = pVtx.x();
    mVy = pVtx.y();
    mVz = pVtx.z();

    
    bfield = picoEvent->bField();

    mVpdVz = picoEvent->vzVpd();

    hVz->Fill(mVz);
    hVpdVz->Fill(mVpdVz);
//    hVxVyVz->Fill(mVx,mVy,mVz);
    hVr->Fill(sqrt(mVy*mVy+mVx*mVx));
    hVzVpdVz->Fill(mVpdVz-mVz);

    mWeight=1;
    double refmultcor;
    //int mCentrality9 = getCentralityBin(mVz,picoEvent->runId(),picoEvent->refMult(),mWeight,refmultcor);
    //int mCentrality9 = getCentralityBin(mVz,picoEvent->runId(),picoEvent->refMult(),mWeight,refmultcor);
    //cout<<"ustc1"<<endl;
/*    refmultCorrUtil = CentralityMaker::instance()->getRefMultCorr();
   // StRefMultCorr * refmultCorrUtil = new StRefMultCorr("refmult","Isobar");
    //StRefMultCorr * refmultCorrUtil = new StRefMultCorr("refmult","Isobar");
    //cout<<"ustc2"<<endl;
    refmultCorrUtil->init(mRunId);
    //refmultCorrUtil->initEvent(picoEvent->refMult(),mVz,picoEvent->ZDCx());
    Bool_t isBadRun_Cen = refmultCorrUtil->isBadRun(mRunId);
    //if(isBadRun_Cen) return kStOk;
    if(isBadRun_Cen){
    //cout<<"ustc3"<<endl;
    }
    Bool_t isPileUpEvt_Cen = !refmultCorrUtil->passnTofMatchRefmultCut(picoEvent->refMult()*1.0,picoEvent->nBTOFMatch()*1.0);
    if(isPileUpEvt_Cen) return kStOk;
    //cout<<"ustc3"<<endl;
    refmultCorrUtil->initEvent(picoEvent->refMult(),mVz,picoEvent->ZDCx());
    Double_t refmultcor_16 = refmultCorrUtil->getRefMultCorr();
    //Double_t reweight = refmultCorrUtil->getWeight();// add to hevt_cent and invariant mass spectra
    //if(refmultCorrUtil->getCentralityBin16()<0 || refmultCorrUtil->getCentralityBin9()<0) return kStOk;
    //cout<<"ustc4"<<endl;
    Int_t mCentrality = refmultCorrUtil->getCentralityBin16();
    Int_t mCentrality9 = refmultCorrUtil->getCentralityBin9();
    //if(mCentrality<0 || mCentrality9<0) return kStOk;

    Double_t reweight = refmultCorrUtil->getWeight();// add to hevt_cent and invariant mass spectra
    //cout<<"reweight: "<<reweight<<endl;*/

    Int_t mCentrality = 1;
    Double_t reweight = 1.0;


    if(2==mMagBins)
      {
      if(bfield<0.) magBufferPointer = 0;
      if(bfield>0.) magBufferPointer = 1;
      }
      else{magBufferPointer =0;}
   if(magBufferPointer!=0 && magBufferPointer!=1)return kStOk;
   cenBufferPointer = mCentrality + 1;
   if(cenBufferPointer<0 || cenBufferPointer>=mCenBins)return kStOk;
   vzBufferPointer = (Int_t)((mVz - anaCuts::vz_low)/(anaCuts::vz_up - anaCuts::vz_low)*mVzBins);
   if(vzBufferPointer<0 || vzBufferPointer>=mVzBins)return kStOk;
   //cout<<"magBufferPointer: "<<magBufferPointer<<" "<<"cenBufferPointer: "<<cenBufferPointer<<" "<<"vzBufferPointer: "<<vzBufferPointer<<endl;
   current_nE=0;
   current_nEPlus=0;
   current_nEMinus=0;

   //current_nE_PE=0;
   //current_nEPlus_PE=0;
   //current_nEMinus_PE=0;

    hevt_Cent->Fill(mCentrality,reweight);
	h_Cen_Vz->Fill(mCentrality,mVz);
	h_Cen_reweight->Fill(mCentrality,reweight);

	hpassevtcut->Fill(4);

    //hnTofMulvsRef->Fill(picoEvent->btofTrayMultiplicity(),picoEvent->refMult());
   // hnTofMatvsRef->Fill(picoEvent->nBTOFMatch(),picoEvent->refMult());
    double ntofhits = 0;
//    int ntrack_tof_hits =0;

    //cout<<"star reconstruction"<<endl;

    electroninfo.clear();
    positroninfo.clear();
    int nTracks = picoDst->numberOfTracks();
	
	
	

   // cout<<"nTracks: "<<nTracks<<endl;
    for (int itrack=0;itrack<nTracks;itrack++){
      StPicoTrack* trk = picoDst->track(itrack);
      //cout<<"step0"<<endl;
      hgDca->Fill(trk->gDCA(mVx,mVy,mVz));
      //cout<<"step0.5"<<endl;
      bool isprimary = trk->isPrimary();
      bool goodtrack = isGoodTrack(trk,trk->gDCA(mVx,mVy,mVz));
      if (!goodtrack) continue;
      if (!isprimary) continue;
    // cout<<"primary track"<<endl;
    //  int index2tof = trk->bTofPidTraitsIndex();
    //  StPicoBTofPidTraits const* const tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);
    //  int tofid = tofPid->btofCellId();
    //  if(tofid==8994 || tofid==8998 || tofid==8999) continue;
      //cout<<"step01"<<endl;
      //StThreeVectorF mom = trk->pMom();
      TVector3 mom = trk->pMom();
	  TVector3 porigin = trk->origin();
      if(mom.Perp()<0.2) continue;
     
      double beta = getTofBeta(trk);
      bool tofmatch = (beta!=std::numeric_limits<float>::quiet_NaN()) && beta>0;
	  
	  if (tofmatch) {
      
      int index2tof = trk->bTofPidTraitsIndex();
      StPicoBTofPidTraits const* const tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);
      int tofid = tofPid->btofCellId();
	  
//	  if(isBadcell(tofid)) continue;


        hinvBetavsP->Fill(mom.Mag(),1./beta);
        //hNsigEvsinvBeta->Fill(trk->nSigmaElectron(),1./beta,mom.Mag());
      }
	  
	  hpt->Fill(mom.Perp());
      hGpt->Fill(trk->gMom().Perp());
      hPhi->Fill(mom.Phi());
      hEta->Fill(mom.Eta());
      // hpDca->Fill(trk->pDca(mVx,mVy,mVz));
      hnHitsFit->Fill(trk->nHitsFit());
	  
	  
      //choose inclusive electron
      // bool isTPCElectron =  trk->nSigmaElectron()<2 && trk->nSigmaElectron()>0.75;
//      if(tofmatch){
//      int index2tof = trk->bTofPidTraitsIndex();
//      StPicoBTofPidTraits const* const tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);
//      int tofid = tofPid->btofCellId();
//      }
    //  cout<<"step02"<<endl;
      bool isTPCElectron=0;


     //Double_t nSigmaE_corr  = recorrection_nSigmaE(trk->nSigmaElectron(),mom.Eta(),mom.Phi(),mCentrality,trk->charge(),0);//0->Ru,1->Zr
     //Double_t nSigmaE_corr  = 0;//0->Ru,1->Zr
     // cout<<"step03"<<endl;


     

       if (mom.Mag()>0.8) isTPCElectron =  trk->nSigmaElectron()<(2) && trk->nSigmaElectron()>(-0.75);
         else isTPCElectron = trk->nSigmaElectron()<(2) && trk->nSigmaElectron()>(3*mom.Mag()-3.15);//20220806
		 

       

       //double mass_square = 0;
       //mass_square = mom.Mag()*mom.Mag()*(1.0/beta*1.0/beta-1);
       //h_mass->Fill(mass_square);
	
	   
	   
      //cout<<"step04"<<endl;
	   

       bool isTOFElectron = tofmatch?fabs(1./beta-1.)<0.025:false;
	  
	   
	   if(isTOFElectron){
       h_nSigmaElectron_tof_p->Fill(mom.Mag(),trk->nSigmaElectron());
       h_nSigmaElectron_tof_pT->Fill(mom.Perp(),trk->nSigmaElectron());
	   }

       bool isElectron = 0;

       isElectron = isTOFElectron && isTPCElectron;//20211219
	   
     // cout<<"store electron"<<endl;
      if (isElectron ) {
       

        //cout<<"electrons"<<endl;

        if(trk->charge()<0)
         {
            particleinfo.charge = trk->charge();
//            cout<<"charge: "<<trk->charge()<<endl;
            particleinfo.pt = mom.Perp();
            particleinfo.eta = mom.Eta();
            particleinfo.phi = mom.Phi();
            particleinfo.p = mom.Mag();
            particleinfo.nSigmaPi = trk->nSigmaPion();
            particleinfo.beta = beta;
            particleinfo.energy = sqrt(pow(M_electron,2.0)+pow(mom.Mag(),2.0));
            particleinfo.p1 = mom.x();
            particleinfo.p2 = mom.y();
            particleinfo.p3 = mom.z();
            electroninfo.push_back(particleinfo);
//            cout<<"debug01"<<endl;

             current_eMinus[current_nEMinus].SetPx(mom.x());
             current_eMinus[current_nEMinus].SetPy(mom.y());
             current_eMinus[current_nEMinus].SetPz(mom.z());
             current_eMinus[current_nEMinus].SetE(sqrt(pow(M_electron,2.0)+pow(mom.Mag(),2.0)));

             // current_p_E0_eMinus[current_nEMinus] = E0;
             // current_p_E_eMinus[current_nEMinus] = E_bemc;
             // current_bemcdz_eMinus[current_nEMinus] = bemcZ;
             // current_bemcdphi_eMinus[current_nEMinus] = bemcPhi;
             // current_bemcTowerPhi_eMinus[current_nEMinus] = bemcTowerPhi;
             // current_bemcTowerEta_eMinus[current_nEMinus] = bemcTowerEta;
             // current_nSMDphi_eMinus[current_nEMinus] = bemcSmdN_phi;
             // current_nSMDeta_eMinus[current_nEMinus] = bemcSmdN_eta;
             
             current_MomMag_eMinus[current_nEMinus] = mom.Mag();
             current_MomPerp_eMinus[current_nEMinus] = mom.Perp();
             current_nSigE_eMinus[current_nEMinus] = trk->nSigmaElectron();
			 
			 current_gMom_eMinus[current_nEMinus].setX(trk->gMom().x()); 
             current_gMom_eMinus[current_nEMinus].setY(trk->gMom().y());
             current_gMom_eMinus[current_nEMinus].setZ(trk->gMom().z());  					
             current_origin_eMinus[current_nEMinus].setX(trk->origin().x());
	         current_origin_eMinus[current_nEMinus].setY(trk->origin().y());
			 current_origin_eMinus[current_nEMinus].setZ(trk->origin().z());
			 
             current_nEMinus++;
            h_electron_count->Fill(1);
         }

        if(trk->charge()>0)
         {
            particleinfo.charge = trk->charge();
//            cout<<"charge: "<<trk->charge()<<endl;
            particleinfo.pt = mom.Perp();
            particleinfo.eta = mom.Eta();
            particleinfo.phi = mom.Phi();
            particleinfo.p = mom.Mag();
            particleinfo.nSigmaPi = trk->nSigmaPion();
            particleinfo.beta = beta;
            particleinfo.energy = sqrt(pow(M_electron,2.0)+pow(mom.Mag(),2.0));
            particleinfo.p1 = mom.x();
            particleinfo.p2 = mom.y();
            particleinfo.p3 = mom.z();
            positroninfo.push_back(particleinfo);

             current_ePlus[current_nEPlus].SetPx(mom.x());
             current_ePlus[current_nEPlus].SetPy(mom.y());
             current_ePlus[current_nEPlus].SetPz(mom.z());
             current_ePlus[current_nEPlus].SetE(sqrt(pow(M_electron,2.0)+pow(mom.Mag(),2.0)));
             // current_p_E0_ePlus[current_nEPlus] = E0;
             // current_p_E_ePlus[current_nEPlus] = E_bemc;
             // current_bemcdz_ePlus[current_nEPlus] = bemcZ;
             // current_bemcdphi_ePlus[current_nEPlus] = bemcPhi;
             // current_bemcTowerPhi_ePlus[current_nEPlus] = bemcTowerPhi;
             // current_bemcTowerEta_ePlus[current_nEPlus] = bemcTowerEta;
             // current_nSMDphi_ePlus[current_nEPlus] = bemcSmdN_phi;
             // current_nSMDeta_ePlus[current_nEPlus] = bemcSmdN_eta;
             
             current_MomMag_ePlus[current_nEPlus] = mom.Mag();
             current_MomPerp_ePlus[current_nEPlus] = mom.Perp();
             current_nSigE_ePlus[current_nEPlus] = trk->nSigmaElectron();
			 
			 
			 current_gMom_ePlus[current_nEPlus].setX(trk->gMom().x()); 
             current_gMom_ePlus[current_nEPlus].setY(trk->gMom().y());
             current_gMom_ePlus[current_nEPlus].setZ(trk->gMom().z());  					
             current_origin_ePlus[current_nEPlus].setX(trk->origin().x());
	         current_origin_ePlus[current_nEPlus].setY(trk->origin().y());
			 current_origin_ePlus[current_nEPlus].setZ(trk->origin().z());
			 
             current_nEPlus++;
            h_electron_count->Fill(2);
//            cout<<"debug02"<<endl;
         }
         current_nE++;

      }

     
      hdEdx->Fill(mom.Mag()*trk->charge(),trk->dEdx());
    }
   // hnTofHitvsRef->Fill(ntofhits,picoEvent->refMult());
  // cout<<"make pair"<<endl;
   
      int x=0;
      int y=0;
      int num_electron = electroninfo.size();
      int num_positron = positroninfo.size();
   
   

      TVector3 momentum_particle;
      TLorentzVector eepair(0,0,0,0);//e1+e2
      TLorentzVector eepair1(0,0,0,0);//e1-e2
      TLorentzVector particle1_4V(0,0,0,0);
      TLorentzVector particle2_4V(0,0,0,0);
      Float_t diffphi = 0;
      for(x=0;x<num_electron;x++)
         {
			 
             particle1_4V.SetPx(electroninfo[x].p1);
             particle1_4V.SetPy(electroninfo[x].p2);
             particle1_4V.SetPz(electroninfo[x].p3);
             particle1_4V.SetE(electroninfo[x].energy);
               for(y=x+1;y<num_electron;y++)
                  {
					  
                    particle2_4V.SetPx(electroninfo[y].p1);
                    particle2_4V.SetPy(electroninfo[y].p2);
                    particle2_4V.SetPz(electroninfo[y].p3);
                    particle2_4V.SetE(electroninfo[y].energy);
                    eepair = particle1_4V + particle2_4V;
                    eepair1 = particle1_4V - particle2_4V;

                    if(gRandom->Rndm()>0.5){
                     eepair1 = particle2_4V - particle1_4V;
                   }

                    diffphi = eepair1.DeltaPhi(eepair);

                    fill_4[0] = eepair.M();
                    fill_4[1] = eepair.Perp();
                    fill_4[2] = diffphi;
                    fill_4[3] = mCentrality;
                    fill_4[4] = reweight;//weight

                    fill_4_phi[0] = eepair.Perp();
                    fill_4_phi[1] = particle1_4V.Phi();
                    fill_4_phi[2] = particle2_4V.Phi();
                    fill_4_phi[3] = mCentrality;
                    fill_4_phi[4] = 1;

                    //cout<<"step06.4"<<endl;
                    hMeeCount_like1->Fill(eepair.M());
                    //cout<<"step06.5"<<endl;
                    
						
                    //cout<<"step06.6"<<endl;
					   Double_t angleV = phiVAngle(particle1_4V,particle2_4V,-1,-1);
                    //cout<<"step06.65"<<endl;
					   Double_t angleVcut = phiVcut->Eval(eepair.M());
                    //cout<<"step06.7"<<endl;

                  if((eepair.M() < mPhiVCutMRange && angleV >  angleVcut) || eepair.M() >= mPhiVCutMRange){
                      hMeeCount_like1_wo_PE->Fill(eepair.M());
                    }  

                  StPhysicalHelixD const p1Line(current_gMom_eMinus[x],current_origin_eMinus[x],bfield*kilogauss,-1);
                  StPhysicalHelixD const p2Line(current_gMom_eMinus[y],current_origin_eMinus[y],bfield*kilogauss,-1);
                  pair<double,double> const ss = p1Line.pathLengths(p2Line);
                  StThreeVectorF const p1AtDca2p2 = p1Line.at(ss.first);
                  StThreeVectorF const p2AtDca2p1 = p2Line.at(ss.second);
                  float DcaDaughters = (p1AtDca2p2 -p2AtDca2p1).mag();
                  
                  StThreeVectorF const p1MomAtDca = p1Line.momentumAt(ss.first,bfield * kilogauss);
                  StThreeVectorF const p2MomAtDca = p2Line.momentumAt(ss.second,bfield * kilogauss);

                  StLorentzVector<float> p2FourMom(p2MomAtDca,p2MomAtDca.massHypothesis(M_electron));
                  StLorentzVector<float> p1FourMom(p1MomAtDca,p1MomAtDca.massHypothesis(M_electron));
                  StLorentzVector<float> pair_helix = p2FourMom + p1FourMom;

                  bool passEtopocut = fabs(DcaDaughters) < anaCuts::EEdcaDaughter;
                  if(!passEtopocut) continue;

                  hMeeCount_like1_helix->Fill(pair_helix.m());
                     
                  if(eepair.M() < mPhiVCutMRange && angleV < angleVcut){
                     hMeeCount_like1_helix_PE->Fill(pair_helix.m());
                  }
						
                  if(fabs(eepair.Rapidity())<=1){
                           if( (angleV>angleVcut && eepair.M()<mPhiVCutMRange) || eepair.M()>=mPhiVCutMRange ){	
                                 hMeelike1_Pt_Cent->Fill(eepair.Perp(),mCentrality,eepair.M(),reweight);
                                 hMeelike1_Pt_phi_Cent->Fill(fill_4);
                           }
					   	
                       //if(eepair.M()>0.4 && eepair.M()<2.6){hMeelike1_Pt_phi1_phi2_Cent->Fill(fill_4_phi);}
                    }
//                    cout<<"debug03"<<endl;
                  }
         }

      for(x=0;x<num_positron;x++)
         {

			 
             particle1_4V.SetPx(positroninfo[x].p1);
             particle1_4V.SetPy(positroninfo[x].p2);
             particle1_4V.SetPz(positroninfo[x].p3);
             particle1_4V.SetE(positroninfo[x].energy);
               for(y=x+1;y<num_positron;y++)
                  {

					  
                    particle2_4V.SetPx(positroninfo[y].p1);
                    particle2_4V.SetPy(positroninfo[y].p2);
                    particle2_4V.SetPz(positroninfo[y].p3);
                    particle2_4V.SetE(positroninfo[y].energy);
                    eepair = particle1_4V + particle2_4V;
                    eepair1 = particle1_4V - particle2_4V;
                    if(gRandom->Rndm()>0.5){
                       eepair1 = particle2_4V - particle1_4V;
                    }

                    diffphi = eepair1.DeltaPhi(eepair);
                    //diffphi = eepair.Phi() - eepair1.Phi();
                    fill_4[0] = eepair.M();
                    fill_4[1] = eepair.Perp();
                    fill_4[2] = diffphi;
                    fill_4[3] = mCentrality;
                    fill_4[4] = reweight;//weight
                    hMeeCount_like2->Fill(eepair.M());

                    Double_t angleV = phiVAngle(particle1_4V,particle2_4V,1,1);
                    Double_t angleVcut = phiVcut->Eval(eepair.M());	

                    if((eepair.M() < mPhiVCutMRange && angleV > angleVcut) || eepair.M() >= mPhiVCutMRange){
                      hMeeCount_like2_wo_PE->Fill(eepair.M());
                    }

                    StPhysicalHelixD const p1Line(current_gMom_ePlus[x],current_origin_ePlus[x],bfield*kilogauss,1);
                    StPhysicalHelixD const p2Line(current_gMom_ePlus[y],current_origin_ePlus[y],bfield*kilogauss,1);
                    pair<double,double> const ss = p1Line.pathLengths(p2Line);
                    StThreeVectorF const p1AtDca2p2 = p1Line.at(ss.first);
                    StThreeVectorF const p2AtDca2p1 = p2Line.at(ss.second);
                    float DcaDaughters = (p1AtDca2p2 -p2AtDca2p1).mag();
                    
                    StThreeVectorF const p1MomAtDca = p1Line.momentumAt(ss.first,bfield * kilogauss);
                    StThreeVectorF const p2MomAtDca = p2Line.momentumAt(ss.second,bfield * kilogauss);
  
                    StLorentzVector<float> p2FourMom(p2MomAtDca,p2MomAtDca.massHypothesis(M_electron));
                    StLorentzVector<float> p1FourMom(p1MomAtDca,p1MomAtDca.massHypothesis(M_electron));
                    StLorentzVector<float> pair_helix = p2FourMom + p1FourMom;

                    bool passEtopocut = fabs(DcaDaughters) < anaCuts::EEdcaDaughter;
                    if(!passEtopocut) continue;

                    hMeeCount_like2_helix->Fill(pair_helix.m());
                    
                    if(eepair.M() < mPhiVCutMRange && angleV < angleVcut){
                      hMeeCount_like2_helix_PE->Fill(pair_helix.m());
                    }
                    
                    if(fabs(eepair.Rapidity())<=1){
						
                        if( (angleV>angleVcut && eepair.M()<mPhiVCutMRange) || eepair.M()>=mPhiVCutMRange ){	
                           hMeelike2_Pt_Cent->Fill(eepair.Perp(),mCentrality,eepair.M(),reweight);
                           hMeelike2_Pt_phi_Cent->Fill(fill_4);
                           }
                    }
                  }
         }

      for(x=0;x<num_positron;x++)
         {

			 
             particle1_4V.SetPx(positroninfo[x].p1);
             particle1_4V.SetPy(positroninfo[x].p2);
             particle1_4V.SetPz(positroninfo[x].p3);
             particle1_4V.SetE(positroninfo[x].energy);
               for(y=0;y<num_electron;y++)
                  {
	  
					  
                    particle2_4V.SetPx(electroninfo[y].p1);
                    particle2_4V.SetPy(electroninfo[y].p2);
                    particle2_4V.SetPz(electroninfo[y].p3);
                    particle2_4V.SetE(electroninfo[y].energy);
                    eepair = particle1_4V + particle2_4V;
                    eepair1 = particle1_4V - particle2_4V;
                    if(gRandom->Rndm()>0.5){
                       eepair1 = particle2_4V - particle1_4V;
                    }

                    diffphi = eepair1.DeltaPhi(eepair);
                    //diffphi = eepair.Phi() - eepair1.Phi();
                    fill_4[0] = eepair.M();
                    fill_4[1] = eepair.Perp();
                    fill_4[2] = diffphi;
                    fill_4[3] = mCentrality;
                    fill_4[4] = reweight;//weight

                    fill_4_phi[0] = eepair.Perp();
                    fill_4_phi[1] = particle1_4V.Phi();
                    fill_4_phi[2] = particle2_4V.Phi();
                    fill_4_phi[3] = mCentrality;
                    fill_4_phi[4] = 1;

                    hMeeCount->Fill(eepair.M());


                     Double_t angleV = phiVAngle(particle1_4V,particle2_4V,1,-1);
                     Double_t angleVcut = phiVcut->Eval(eepair.M());

                     if((eepair.M() < mPhiVCutMRange && angleV >  angleVcut) || eepair.M() >= mPhiVCutMRange){
                         hMeeCount_wo_PE->Fill(eepair.M());
                        } 
                     
                    StPhysicalHelixD const p1Line(current_gMom_ePlus[x],current_origin_ePlus[x],bfield*kilogauss,1);
                    StPhysicalHelixD const p2Line(current_gMom_eMinus[y],current_origin_eMinus[y],bfield*kilogauss,-1);
                    pair<double,double> const ss = p1Line.pathLengths(p2Line);
                    StThreeVectorF const p1AtDca2p2 = p1Line.at(ss.first);
                    StThreeVectorF const p2AtDca2p1 = p2Line.at(ss.second);
                    float DcaDaughters = (p1AtDca2p2 -p2AtDca2p1).mag();
                    
                    StThreeVectorF const p1MomAtDca = p1Line.momentumAt(ss.first,bfield * kilogauss);
                    StThreeVectorF const p2MomAtDca = p2Line.momentumAt(ss.second,bfield * kilogauss);
  
                    StLorentzVector<float> p2FourMom(p2MomAtDca,p2MomAtDca.massHypothesis(M_electron));
                    StLorentzVector<float> p1FourMom(p1MomAtDca,p1MomAtDca.massHypothesis(M_electron));
                    StLorentzVector<float> pair_helix = p2FourMom + p1FourMom;

                    bool passEtopocut = fabs(DcaDaughters) < anaCuts::EEdcaDaughter;
                    if(!passEtopocut) continue;

                    hMeeCount_helix->Fill(pair_helix.m());

                    if(eepair.M() < mPhiVCutMRange && angleV <  angleVcut){
                      hMeeCount_helix_PE->Fill(pair_helix.m());
                    }
                  
                  if(fabs(eepair.Rapidity())<=1){
                        if( (angleV>angleVcut && eepair.M()<mPhiVCutMRange) || eepair.M()>=mPhiVCutMRange ){	
                           hMee_Pt_Cent->Fill(eepair.Perp(),mCentrality,eepair.M(),reweight);
                           hMee_Pt_phi_Cent->Fill(fill_4);
                        }
                      // if(eepair.M()>0.4 && eepair.M()<2.6){hMee_Pt_phi1_phi2_Cent->Fill(fill_4_phi);}
                       }
                  }
         }
      //cout<<"num_positron: "<<num_positron<<endl;
      //cout<<"num_electron: "<<num_electron+num_positron<<endl;


//mixed event
    //cout<<"mixed"<<endl;
    for(Int_t iBufferEvent=0;iBufferEvent<nEventsInBuffer[magBufferPointer][cenBufferPointer][vzBufferPointer];iBufferEvent++)
       {
           for(x=0;x<current_nEPlus;x++)
              {
                  for(y=0;y<buffer_nEMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][iBufferEvent];y++)
                        {
                             eepair = current_ePlus[x] + buffer_eMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][iBufferEvent][y];
                             eepair1 = current_ePlus[x] - buffer_eMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][iBufferEvent][y];
                             if(gRandom->Rndm()>0.5){
                                //eepair1 = particle2_4V - particle1_4V;
                                eepair1 = buffer_eMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][iBufferEvent][y] - current_ePlus[x];
                             }

                             diffphi = eepair1.DeltaPhi(eepair);
                             //diffphi = eepair.Phi() - eepair1.Phi();
                             fill_4[0] = eepair.M();
                             fill_4[1] = eepair.Perp();
                             fill_4[2] = diffphi;
                             fill_4[3] = mCentrality;
                             fill_4[4] = reweight;//weight
                             if(fabs(eepair.Rapidity())<=1){
								 
							   Double_t angleV = phiVAngle(particle1_4V,particle2_4V,1,-1);
							   Double_t angleVcut = phiVcut->Eval(eepair.M());	 
								
							   	
							   if( (angleV>angleVcut && eepair.M()<mPhiVCutMRange) || eepair.M()>=mPhiVCutMRange ){	
                               hMeeMix_Pt_Cent->Fill(eepair.Perp(),mCentrality,eepair.M(),reweight);
                               hMeeMix_Pt_phi_Cent->Fill(fill_4);
							   }
                             }
                        }
              }

          for(x=0;x<current_nEMinus;x++)
              {
                  for(y=0;y<buffer_nEPlus[magBufferPointer][cenBufferPointer][vzBufferPointer][iBufferEvent];y++)
                        {
                             eepair = current_eMinus[x] + buffer_ePlus[magBufferPointer][cenBufferPointer][vzBufferPointer][iBufferEvent][y];
                             eepair1 = current_eMinus[x] - buffer_ePlus[magBufferPointer][cenBufferPointer][vzBufferPointer][iBufferEvent][y];
                             if(gRandom->Rndm()>0.5){
                                //eepair1 = particle2_4V - particle1_4V;
                                eepair1 = buffer_ePlus[magBufferPointer][cenBufferPointer][vzBufferPointer][iBufferEvent][y] - current_eMinus[x];
                             }

                             diffphi = eepair1.DeltaPhi(eepair);
                             //diffphi = eepair.Phi() - eepair1.Phi();
                             fill_4[0] = eepair.M();
                             fill_4[1] = eepair.Perp();
                             fill_4[2] = diffphi;
                             fill_4[3] = mCentrality;
                             fill_4[4] = reweight;//weight
                             if(fabs(eepair.Rapidity())<=1){
							   
                               Double_t angleV = phiVAngle(particle1_4V,particle2_4V,-1,1);
							   Double_t angleVcut = phiVcut->Eval(eepair.M());							   
							   if( (angleV>angleVcut && eepair.M()<mPhiVCutMRange) || eepair.M()>=mPhiVCutMRange ){
                               hMeeMix_Pt_Cent->Fill(eepair.Perp(),mCentrality,eepair.M(),reweight);
                               hMeeMix_Pt_phi_Cent->Fill(fill_4);
							   }
                             }
                        }
              }

          for(x=0;x<current_nEPlus;x++)
              {
                  for(y=0;y<buffer_nEPlus[magBufferPointer][cenBufferPointer][vzBufferPointer][iBufferEvent];y++)
                        {
                             eepair = current_ePlus[x] + buffer_ePlus[magBufferPointer][cenBufferPointer][vzBufferPointer][iBufferEvent][y];
                             eepair1 = current_ePlus[x] - buffer_ePlus[magBufferPointer][cenBufferPointer][vzBufferPointer][iBufferEvent][y];
                             if(gRandom->Rndm()>0.5){
                                //eepair1 = particle2_4V - particle1_4V;
                                eepair1 = buffer_ePlus[magBufferPointer][cenBufferPointer][vzBufferPointer][iBufferEvent][y] - current_ePlus[x];
                             }

                             diffphi = eepair1.DeltaPhi(eepair);
                             //diffphi = eepair.Phi() - eepair1.Phi();
                             fill_4[0] = eepair.M();
                             fill_4[1] = eepair.Perp();
                             fill_4[2] = diffphi;
                             fill_4[3] = mCentrality;
                             fill_4[4] = reweight;//weight
                             if(fabs(eepair.Rapidity())<=1){
								 
							    Double_t angleV = phiVAngle(particle1_4V,particle2_4V,1,1);
							   Double_t angleVcut = phiVcut->Eval(eepair.M());
                               if( (angleV>angleVcut && eepair.M()<mPhiVCutMRange) || eepair.M()>=mPhiVCutMRange ){							   
                               hMeeMixlike2_Pt_Cent->Fill(eepair.Perp(),mCentrality,eepair.M(),reweight);
                               hMeeMixlike2_Pt_phi_Cent->Fill(fill_4);
							   }
                             }
                        }
              }

          for(x=0;x<current_nEMinus;x++)
              {
                  for(y=0;y<buffer_nEMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][iBufferEvent];y++)
                        {
                             eepair = current_eMinus[x] + buffer_eMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][iBufferEvent][y];
                             eepair1 = current_eMinus[x] - buffer_eMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][iBufferEvent][y];
                             if(gRandom->Rndm()>0.5){
                                //eepair1 = particle2_4V - particle1_4V;
                                eepair1 = buffer_eMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][iBufferEvent][y] - current_eMinus[x];
                             }

                             diffphi = eepair1.DeltaPhi(eepair);
                             //diffphi = eepair.Phi() - eepair1.Phi();
                             fill_4[0] = eepair.M();
                             fill_4[1] = eepair.Perp();
                             fill_4[2] = diffphi;
                             fill_4[3] = mCentrality;
                             fill_4[4] = reweight;//weight
                             if(fabs(eepair.Rapidity())<=1){
								 
							   Double_t angleV = phiVAngle(particle1_4V,particle2_4V,-1,-1);
							   Double_t angleVcut = phiVcut->Eval(eepair.M());
							   if( (angleV>angleVcut && eepair.M()<mPhiVCutMRange) || eepair.M()>=mPhiVCutMRange ){ 
                               hMeeMixlike1_Pt_Cent->Fill(eepair.Perp(),mCentrality,eepair.M(),reweight);
                               hMeeMixlike1_Pt_phi_Cent->Fill(fill_4);
							   }
                             }
                        }
              }
       }

    // cout<<"end mixed even pair"<<endl;
       copyCurrentToBuffer();
  //   cout<<"end mixed Event"<<endl;
    } //Good Event
 }//bad run cut
 

//  cout<<"end make"<<endl;
  return kStOK;
}
bool StPicoDstarMixedMaker::isGoodTrigger(StPicoEvent const* const picoEvent) const
{
  for (auto trg : anaCuts::triggers)
  {
    if (picoEvent->isTrigger(trg)) return true;
  }

  return false;
}
/*bool StPicoDstarMixedMaker::isGoodTrack(StPicoTrack const* trk, float dca) const
{
  // StThreeVectorF const vtx = mPicoDstMaker->picoDst()->event()->primaryVertex();
  if(trk->pMom().Perp()>1){
    return  fabs(trk->nHitsFit()) >= anaCuts::NHitsFit &&
    fabs(trk->pMom().Eta())<anaCuts::Eta &&
    fabs(trk->nHitsFit()*1.0/trk->nHitsMax())>=anaCuts::NHitsFit2Poss &&
    fabs(trk->nHitsDedx())>=anaCuts::NHitsDedx && fabs(dca)<=anaCuts::Dca_hpT;
    // fabs(trk->nHitsDedx())>=anaCuts::NHitsDedx &&
    // fabs( trk->gDCA(vtx.x() , vtx.y(), vtx.z() )) <= anaCuts::Dca;
    }else{
    return  fabs(trk->nHitsFit()) >= anaCuts::NHitsFit &&
    trk->pMom().Perp() >=0.2 &&
    fabs(trk->pMom().Eta())<anaCuts::Eta &&
    fabs(trk->nHitsFit()*1.0/trk->nHitsMax())>=anaCuts::NHitsFit2Poss &&
    fabs(trk->nHitsDedx())>=anaCuts::NHitsDedx && fabs(dca)<=anaCuts::Dca_lpT;
    // fabs(trk->nHitsDedx())>=anaCuts::NHitsDedx &&
    // fabs( trk->gDCA(vtx.x() , vtx.y(), vtx.z() )) <= anaCuts::Dca;
   }
}*/

//////////
bool StPicoDstarMixedMaker::isGoodTrack(StPicoTrack const* trk, float dca) const
{
  // StThreeVectorF const vtx = mPicoDstMaker->picoDst()->event()->primaryVertex();
    return  fabs(trk->nHitsFit()) >= anaCuts::NHitsFit &&
    trk->pMom().Perp() > anaCuts::pPt &&
    fabs(trk->pMom().Eta())<anaCuts::Eta &&
    fabs(trk->nHitsFit()*1.0/trk->nHitsMax())>=anaCuts::NHitsFit2Poss &&
    fabs(trk->nHitsDedx())>=anaCuts::NHitsDedx && fabs(dca)<=anaCuts::Dca;
    // fabs(trk->nHitsDedx())>=anaCuts::NHitsDedx &&
    // fabs( trk->gDCA(vtx.x() , vtx.y(), vtx.z() )) <= anaCuts::Dca;
}
bool StPicoDstarMixedMaker::isGoodQaTrack(StPicoTrack const* const trk) const
{
  // StThreeVectorF vtx = mPicoDstMaker->picoDst()->event()->primaryVertex();
  return trk->gPt() > anaCuts::qaGPt && fabs(trk->nHitsFit()) >= anaCuts::qaNHitsFit &&
    fabs(trk->gMom().Eta())<anaCuts::qaEta &&
    // fabs(trk->nHitsDedx())>=anaCuts::qaNHitsDedx && fabs(trk->gDCA(vtx.x(),vtx.y(),vtx.z()))<=anaCuts::qaDca;
    fabs(trk->nHitsDedx())>=anaCuts::qaNHitsDedx ;
}
bool StPicoDstarMixedMaker::isGoodQaEvent(StPicoEvent const* const picoEvent) const
{
  TVector3 pVtx = picoEvent->primaryVertex();
 // cout<<"is good Qa event Vz: "<<pVtx.z()<<" "<<anaCuts::vz<<endl;
 // StThreeVectorF pVtx = picoEvent->primaryVertex();
  return fabs(pVtx.z()) < anaCuts::qavz &&
    // fabs(pVtx.z() - picoEvent->vzVpd()) < anaCuts::qavzVpdVz &&
    !(fabs(pVtx.x()) < anaCuts::qaVerror && fabs(pVtx.y()) < anaCuts::qaVerror &&
        fabs(pVtx.z()) < anaCuts::qaVerror) &&
    sqrt(TMath::Power(pVtx.x(), 2) + TMath::Power(pVtx.y(), 2)) <=  anaCuts::qaVr;
}
bool StPicoDstarMixedMaker::isGoodEvent(StPicoEvent const* const picoEvent) const
{
  TVector3 pVtx = picoEvent->primaryVertex();
 // StThreeVectorF pVtx = picoEvent->primaryVertex();
 // cout<<"is good event Vz: "<<pVtx.z()<<" "<<anaCuts::vz<<endl;
 // cout<<"is good event Vx: "<<pVtx.x()<<" "<<anaCuts::Verror<<endl;
 // cout<<"is good event Vz: "<<pVtx.y()<<" "<<anaCuts::Vr<<endl;
 // cout<<"is good event VzVpd: "<<picoEvent->vzVpd()<<" "<<anaCuts::vzVpdVz<<endl;
  //return fabs(pVtx.z()) < anaCuts::vz &&
  return pVtx.z() > anaCuts::vz_low &&
         pVtx.z() < anaCuts::vz_up &&
     fabs(pVtx.z() - picoEvent->vzVpd()) < anaCuts::vzVpdVz &&
    !(fabs(pVtx.x()) < anaCuts::Verror && fabs(pVtx.y()) < anaCuts::Verror &&
        fabs(pVtx.z()) < anaCuts::Verror) &&
    sqrt(TMath::Power(pVtx.x(), 2) + TMath::Power(pVtx.y(), 2)) <=  anaCuts::Vr;
}
float StPicoDstarMixedMaker::getTofBeta(StPicoTrack const* const trk) const
{
  int index2tof = trk->bTofPidTraitsIndex();
  float beta = std::numeric_limits<float>::quiet_NaN();
  if (index2tof >= 0)
  {
    StPicoBTofPidTraits const* const tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);
    if (tofPid)
    {
      beta = tofPid->btofBeta();
      if (beta < 1e-4)
      {
        //StThreeVectorF const vtx = mPicoDstMaker->picoDst()->event()->primaryVertex();
        TVector3 const vtx3 = mPicoDstMaker->picoDst()->event()->primaryVertex();
        StThreeVectorF vtx(vtx3.x(),vtx3.y(),vtx3.z());
        TVector3 const btofHitPos3 = tofPid->btofHitPos();
      //   StThreeVectorF const btofHitPos = tofPid->btofHitPos();
       StThreeVectorF btofHitPos(btofHitPos3.x(),btofHitPos3.y(),btofHitPos3.z());
        StPicoPhysicalHelix helix = trk->helix(mPicoDstMaker->picoDst()->event()->bField());
       // StPhysicalHelixD helix = trk->helix(mPicoDstMaker->picoDst()->event()->bField());
        float L = tofPathLength(&vtx, &btofHitPos, helix.curvature());
        float tof = tofPid->btof();
        if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
        else beta = std::numeric_limits<float>::quiet_NaN();
      }
    }
  }
  return beta;
}

int StPicoDstarMixedMaker::getCentralityBin(float vz,int runId,float mult,float &weight,double &refmultcor) const
{
//  cout<<"ok1.1"<<endl;
  float mult_corr = mult;
  //cout<<"ok1.5"<<endl;
  //if trigger 580001, correct the refmult
  mult+=gRandom->Rndm();
  if (runId<18156031) {
    float fvz = 0;
    for (int i=0;i<anaCuts::nparVz_mult;i++){
      fvz +=anaCuts::parVz_mult[i]*std::pow(vz,i);
    }
    mult=mult*1.0*anaCuts::parVz_mult[0]/fvz;
  }
  refmultcor=mult;
  weight = reweight(mult);
//  cout<<"ok1.6"<<endl;
  for (int cent=0;cent<anaCuts::nCent;cent++)
  {
    if (mult<anaCuts::Refmult_cent[cent]) return cent-1;
  }
  return anaCuts::nCent-1;
}
float StPicoDstarMixedMaker::reweight(float x) const
{
  x+=gRandom->Rndm();
  // float p[7] = {3.9,-204.4,1.85,24.3,-0.01746,6405,3.7e-5};
  float p[5] = {0.811,238.9,24.31,-25,6.325e-5};
  // return 1;
  if (x>70) return 1;
  // else return p[0] + p[1]/(p[2]*x + p[3]) + p[4]*(p[2]*x + p[3]) + p[5]/(p[2]*x + p[3])/(p[2]*x + p[3]) + p[6]*(p[2]*x + p[3])*(p[2]*x + p[3]);
  else return p[0] + p[1]/(p[2]*x + p[3]) + p[4]*(p[2]*x + p[3]);
}
void StPicoDstarMixedMaker::copyCurrentToBuffer()
{

         //cout<<"magBufferPointer: "<<magBufferPointer<<" "<<"cenBufferPointer: "<<cenBufferPointer<<" "<<"vzBufferPointer: "<<vzBufferPointer<<endl;
	if(nEventsInBuffer[magBufferPointer][cenBufferPointer][vzBufferPointer]>=mMaxEventsInBuffer) bufferFullFlag[magBufferPointer][cenBufferPointer][vzBufferPointer] = kTRUE;
	//cout<<"nEventsInBuffer: "<<nEventsInBuffer[magBufferPointer][cenBufferPointer][vzBufferPointer]<<endl;
        Int_t eventPointer = -1;
	if(bufferFullFlag[magBufferPointer][cenBufferPointer][vzBufferPointer]){
	//	cout<<"Random star"<<endl;
                eventPointer = (Int_t)myRandom->Uniform(0,mMaxEventsInBuffer-1.e-6);
        //        cout<<"Random: "<<eventPointer<<endl;
	}else{
		eventPointer = nEventsInBuffer[magBufferPointer][cenBufferPointer][vzBufferPointer];
	}
      // cout<<"eventPointer: "<<eventPointer<<"current_nEPlus: "<<current_nEPlus<<"current_nEMinus: "<<current_nEMinus<<endl;
	buffer_nEPlus[magBufferPointer][cenBufferPointer][vzBufferPointer][eventPointer] = current_nEPlus;
	for(Int_t i=0;i<current_nEPlus;i++){
		buffer_ePlus[magBufferPointer][cenBufferPointer][vzBufferPointer][eventPointer][i] = current_ePlus[i];
	}

	buffer_nEMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][eventPointer] = current_nEMinus;
	for(Int_t i=0;i<current_nEMinus;i++){
		buffer_eMinus[magBufferPointer][cenBufferPointer][vzBufferPointer][eventPointer][i] = current_eMinus[i];
	}

	if(nEventsInBuffer[magBufferPointer][cenBufferPointer][vzBufferPointer]<mMaxEventsInBuffer){
		nEventsInBuffer[magBufferPointer][cenBufferPointer][vzBufferPointer]++;
	}
}

Double_t StPicoDstarMixedMaker::shiftPhi(Double_t phi) const
{
    while (phi < 0)
    {
      phi = phi + 2*TMath::Pi();
    }

    while (phi > 2*TMath::Pi())
    {
      phi = phi - 2*TMath::Pi();
    }

    return phi;
}

Double_t StPicoDstarMixedMaker::fold2Pi(Double_t phi) const
{

   if(phi < 0 || phi > 2*TMath::Pi())
   {
     phi = shiftPhi(phi);
   }
   if(phi > TMath::Pi())
   {
     phi = TMath::Pi() * 2 - phi;
   }

   return phi;
}

double StPicoDstarMixedMaker::recorrection_nSigmaE(double nSigmaE,double trk_eta,double trk_phi,int mCentrality,int trk_charge,int sys) const
{

	double nSigmaE_eta = 0;
	double nSigmaE_phi = 0;
	double nSigmaE_phi_mean = 0;
	int cen_index = 0;
	int phi_bin = -1;

	for(int n = 0;n<nCenBins;n++){
		if(mCentrality >=mCenBinLow[n] && mCentrality <=mCenBinHi[n]){
			cen_index = n;
			break;
		}
	}

	if(mCentrality == -1){
		cen_index = 0;
	}

	if(trk_charge>0){

		nSigmaE_eta = f_PEPlusnSigmaE_eta_2D[sys][cen_index]->Eval(trk_eta);

		phi_bin=pPEPlusnSigmaE_phi_2D[sys][cen_index]->FindBin(trk_phi);
		nSigmaE_phi = pPEPlusnSigmaE_phi_2D[sys][cen_index]->GetBinContent(phi_bin);

		nSigmaE_phi_mean = 0;
		for(int m = 1;m <= pPEPlusnSigmaE_phi_2D[sys][cen_index]->GetNbinsX();m++){
			nSigmaE_phi_mean = nSigmaE_phi_mean+pPEPlusnSigmaE_phi_2D[sys][cen_index]->GetBinContent(m);
		}
		nSigmaE_phi_mean = nSigmaE_phi_mean*1.0/pPEPlusnSigmaE_phi_2D[sys][cen_index]->GetNbinsX();

	}else{

		nSigmaE_eta = f_PEMinusnSigmaE_eta_2D[sys][cen_index]->Eval(trk_eta);

		phi_bin=pPEMinusnSigmaE_phi_2D[sys][cen_index]->FindBin(trk_phi);
		nSigmaE_phi = pPEMinusnSigmaE_phi_2D[sys][cen_index]->GetBinContent(phi_bin);

		nSigmaE_phi_mean = 0;
		for(int m = 1;m <= pPEMinusnSigmaE_phi_2D[sys][cen_index]->GetNbinsX();m++){
			nSigmaE_phi_mean = nSigmaE_phi_mean+pPEMinusnSigmaE_phi_2D[sys][cen_index]->GetBinContent(m);
		}
		nSigmaE_phi_mean = nSigmaE_phi_mean*1.0/pPEMinusnSigmaE_phi_2D[sys][cen_index]->GetNbinsX();
	}


   // cout<<"nSigmaE: "<<nSigmaE<<endl;
//	cout<<"nSigmaE_eta: "<<nSigmaE_eta<<" trk_eta: "<<trk_eta<<endl;
//	cout<<"nSigmaE_phi: "<<nSigmaE_phi<<" trk_phi: "<<trk_phi<<" phi_bin: "<<phi_bin<<endl;
//	cout<<"nSigmaE_phi_mean: "<<nSigmaE_phi_mean<<endl;


	return (nSigmaE - (nSigmaE_eta + nSigmaE_phi - nSigmaE_phi_mean));

}

Double_t StPicoDstarMixedMaker::phiVAngle(TLorentzVector e1, TLorentzVector e2, Int_t q1, Int_t q2) const
{
	Double_t pt1 = e1.Pt();
	Double_t eta1 = e1.Eta();
	Double_t phi1 = e1.Phi();

	Double_t pt2 = e2.Pt();
	Double_t eta2 = e2.Eta();
	Double_t phi2 = e2.Phi();

	TVector3 e1Mom,e2Mom;
	if(q1>0&&q2<0){
		e2Mom.SetPtEtaPhi(pt1,eta1,phi1);//e+
		e1Mom.SetPtEtaPhi(pt2,eta2,phi2);//e-
	}else if(q1<0&&q2>0){
		e2Mom.SetPtEtaPhi(pt2,eta2,phi2);//e+
		e1Mom.SetPtEtaPhi(pt1,eta1,phi1);//e-
	}else if(q1==q2&&TMath::Abs(q1)==1){
		Double_t ran = myRandom_phi->Uniform(-1,1);
		if(ran>0){
			e2Mom.SetPtEtaPhi(pt1,eta1,phi1);
			e1Mom.SetPtEtaPhi(pt2,eta2,phi2);
		}
		else{
			e2Mom.SetPtEtaPhi(pt2,eta2,phi2);
			e1Mom.SetPtEtaPhi(pt1,eta1,phi1);
		}
	}else return -1;
	Double_t mN = 0.;
	if(bfield<0.) mN = -1.;
	if(bfield>0.) mN = 1.;

	TVector3 pu=e1Mom+e2Mom;
	TVector3 pv=e1Mom.Cross(e2Mom);
	TVector3 pw=pu.Cross(pv);
	TVector3 pnz(0.,0.,mN);
	TVector3 pwc=pu.Cross(pnz);

	Double_t angleV = pw.Angle(pwc);

	return angleV;
}
