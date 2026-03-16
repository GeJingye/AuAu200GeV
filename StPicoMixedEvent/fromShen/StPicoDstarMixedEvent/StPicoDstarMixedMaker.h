#ifndef StPicoDstarMixedMaker_h
#define StPicoDstarMixedMaker_h

/* **************************************************
 *  A Maker to read a StPicoEvent and StPicoDstarEvent
 *  simultaneously and do analysis.
 *
 *  Authors: Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */

#include "TChain.h"
#include "StMaker.h"
#include "StThreeVectorF.hh"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THn.h"
#include "TF1.h"
class TString;
class TFile;
class TNtuple;
class StPicoTrack;
class StPicoDstMaker;
class StPicoEvent;
class StRefMultCorr;
class CentralityMaker;

struct ParticleInfo_Electron
{
    Int_t charge;
    Float_t pt;
    Float_t eta;
    Float_t phi;
    Float_t p;
    Float_t nSigmaPi;
    Float_t beta;
   // Float_t betaElectron;
   // Float_t dBetaElectron;
    Float_t energy;
    Float_t p1;
    Float_t p2;
    Float_t p3;
};

//vector<ParticleInfo> ElectronInfo;

class StPicoDstarMixedMaker : public StMaker
{
  public:
    StPicoDstarMixedMaker(char const * name, TString const inputFilesList,
        TString const outBaseName, StPicoDstMaker* picoDstMaker);
    virtual ~StPicoDstarMixedMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();
  private:

    StPicoDstarMixedMaker() {}
    void initHists();
    bool isGoodTrigger(StPicoEvent const*) const;
    bool isGoodQaEvent(StPicoEvent const* const picoEvent) const;
    bool isGoodEvent(StPicoEvent const* const picoEvent) const;
    bool isGoodQaTrack(StPicoTrack const* const trk) const;
    bool isGoodTrack(StPicoTrack const* trk,float dca) const;
    float getTofBeta(StPicoTrack const* const trk) const;
    StPicoDstMaker* mPicoDstMaker;
    TString mInputFilesList;
    TString mOutFileBaseName;
    bool isBadrun(Int_t runId);
	bool isBadcell(Int_t cellId);

    int getCentralityBin(float vz,int runId,float mult,float &weight,double &refmultcor) const;
    float reweight(float x) const;
    void copyCurrentToBuffer();
    Double_t shiftPhi(Double_t phi) const;
    Double_t fold2Pi(Double_t phi) const;
    double recorrection_nSigmaE(double nSigmaE,double trk_eta,double trk_phi,int mCentrality,int trk_charge,int sys) const;
    Double_t phiVAngle(TLorentzVector e1, TLorentzVector e2, Int_t q1, Int_t q2) const;
	// -------------- USER variables -------------------------
    // add your member variables here.
    // Remember that ntuples size can be really big, use histograms where appropriate

  public:
    void setRunNumList(string list){
      mRunNumList = list;
    }
    void setRunbyRunQA(bool b){
      mRunbyRunQA = b;
    }
    void setTotalEvents(int events){
      nEvent_tof = events;
    }

    void getBadruns(string inputFileName);
	void getBadcells(string inputFileName);
  private:
    TFile* mFile;

    TFile* mFile_RunID;//check ToF problems
    TH2F* hinvBetavsP_RunID[78];//check ToF problems

    TFile* mFile_analysis;
    TFile* mFile_tof;//tof calibration

    std::map<int,int> runnum;
    string mRunNumList;
    vector<int> mBadRun;
	vector<int> mBadCell;

    //Event level
    int  mRunId;
    float mVpdVz;
    float mRefmult;
    float mVpdHitEast;
    float mVpdHitWest;
    //primaryVertex vertex
    float mVx;
    float mVy;
    float mVz;
    float mVr;
    //event QA
    TH3F* hVxVyVz;
    TH1F* hVz;
    TH1F* hVz_total;
    TH1F* hVpdVz;
    TH1F* hVr;
    TH1F* hVzVpdVz;
    TH2F* hnEvsPhivsVz;
    TH2F* hnEvsEtavsVz;
    TH2F* hnTofMulvsRef;
    TH2F* hnTofMatvsRef;
    TH2F* hnTofHitvsRef;

    TH1D* hevt;
    TH1D* hevt_Cent;
    TH1D* hevtcut;
    TH1D* hevt_bad_cut;
    TH1D* hpassevtcut;
    TH1F* hrefmult;

    TH2F* h_Cen_Vz;
    TH2F* h_Cen_reweight;

   //tracl level QA
    TH3F* hNsigEvsinvBeta;
    TH1F* hnHitsFit;
    TH1F* hpDca;
    TH1F* hgDca;
    TH2F* hinvBetavsP;
    TH2F* hinvBetavsY;
    TH2F* hdEdx;
    TH1F* hpt;
    TH1F* hGpt;
    TH1F* hEta;
    TH1F* hPhi;
    TH3F* hBadTofId;

    TH1F* h_p_E0;
    TH1F* h_p_E0_PE;
    TH1F* h_p_E0_electron;
    TH1F* h_p_E;
    TH1F* h_p_E_PE;
    TH1F* h_p_E_electron;
    TH2F* h_bemcdz_pT;
    TH2F* h_bemcdz_pT_PE;
    TH2F* h_bemcdz_pT_electron;
    TH2F* h_bemcdphi_pT;
    TH2F* h_bemcdphi_pT_PE;
    TH2F* h_bemcdphi_pT_electron;
    TH2F* h_bemcTowerPhi_pT;
    TH2F* h_bemcTowerPhi_pT_PE;
    TH2F* h_bemcTowerPhi_pT_electron;
    TH2F* h_bemcTowerEta_pT;
    TH2F* h_bemcTowerEta_pT_PE;
    TH2F* h_bemcTowerEta_pT_electron;
    TH2F* h_nSMDphi_pT;
    TH2F* h_nSMDphi_pT_PE;
    TH2F* h_nSMDphi_pT_electron;
    TH2F* h_nSMDeta_pT;
    TH2F* h_nSMDeta_pT_PE;
    TH2F* h_nSMDeta_pT_electron;
    TH2F * h_pOverE0_p;
    TH2F * h_pOverE_p;
    TH2F * h_pOverE0_p_PE;
    TH2F * h_pOverE0_p_electron;
    TH2F * h_pOverE_p_PE;
    TH2F * h_pOverE_p_electron;
    TH2F * h_E0Overp_pT;
    TH2F * h_E0Overp_pT_bemc;
    TH2F * h_E0Overp_pT_proton;
    TH2F * h_EOverp_pT;
    TH2F * h_E0Overp_pT_PE;
    TH2F * h_E0Overp_pT_electron;
    TH2F * h_EOverp_pT_PE;
    TH2F * h_EOverp_pT_electron;

    TH2F * h_pvsE0_pion;
    TH2F * h_pvsE0_electron;
    TH2F * h_pvsE0;
    TH2F * h_pvsE_pion;
    TH2F * h_pvsE_electron;
    TH2F * h_pvsE;
    TH3F * h_pvsE0vsE_pion;
    TH3F * h_pvsE0vsE_electron;
    TH3F * h_pvsE0vsE;
    TH2F * h_E_E0_ratio;
    TH2F * h_E_E0_ratio_pion;
    TH2F * h_E_E0_ratio_electron;

    //invariant mass electron
    TH1F* hMeeCount;//unlike sign
    TH1F* hMeeCount_like1;//like sign electron
    TH1F* hMeeCount_like2;//like sign positron

    TH1F* hMeeCount_helix;//unlike sign, using global track and helix
    TH1F* hMeeCount_like1_helix;//like sign electron
    TH1F* hMeeCount_like2_helix;//like sign positron

    TH1F* hMeeCount_wo_PE;//unlike sign, using phiV cut for exclude PE
    TH1F* hMeeCount_like1_wo_PE;//like sign electron
    TH1F* hMeeCount_like2_wo_PE;//like sign positron

    TH1F* hMeeCount_helix_PE;//unlike sign, using global track and helix, using phiV cut to select PE
    TH1F* hMeeCount_like1_helix_PE;//like sign electron
    TH1F* hMeeCount_like2_helix_PE;//like sign positron



    TH3F* hMee_Pt_Cent;
    TH3F* hMeelike1_Pt_Cent;
    TH3F* hMeelike2_Pt_Cent;
    TH1F* h_electron_count;
	TH1F* h_electron_count_after_PE;
    TH2F* h_nSigmaElectron_p;//eID
    TH2F* h_nSigmaElectron_pT;//eID
    TH2F* h_nSigmaElectron_tof_p;//eID
    TH2F* h_nSigmaElectron_tof_pT;//eID
    TH2F* h_nSigmaElectron_bemc_p;//eID
    TH2F* h_nSigmaElectron_bemc_pT;//eID
    TH2F* h_nSigmaElectron_tof_bemc_p;//eID
    TH2F* h_nSigmaElectron_tof_bemc_pT;//eID
    TH1F* h_isBemcMatch;
    TH1F* h_mass;
    TH2F* h_pure_pion_nSigE_pT;
    TH2F* h_pure_pion_nSigE_p;
    TH2F* h_pure_proton_nSigE_pT;
    TH2F* h_pure_proton_nSigE_p;
    TH2F* h_pure_kaon_nSigE_pT;
    TH2F* h_pure_kaon_nSigE_p;
    TH2F* h_pure_Merged_pion_nSigE_pT;
    TH2F* h_pure_Merged_pion_nSigE_p;
    TH2F* h_pure_electron_nSigE_pT;
    TH2F* h_pure_electron_nSigE_p;
	
	TH2F *hULAngleVvsM;
    TH2F *hLPosAngleVvsM;
    TH2F *hLNegAngleVvsM;
	
	TH2F *hULAngleVvsM_pos;
	TH2F *hULAngleVvsM_neg;
	
	TH2F *hLPPAngleVvsM_pos;
	TH2F *hLPPAngleVvsM_neg;
	
	TH2F *hLMMAngleVvsM_pos;
	TH2F *hLMMAngleVvsM_neg;


    int mCentrality=0;
    float mWeight;
    float M_electron=0.000511;//GeV
    float M_pion=0.139570;//GeV
//    TLorentzVector eepair;
//    TLorentzVector particle1_4V;
//    TLorentzVector particle2_4V;
    int magBufferPoint=0;
    int cenBufferPoint=0;
    int vzBufferPoint=0;
    TH3F* hMeeMix_Pt_Cent;
    TH3F* hMeeMixlike1_Pt_Cent;
    TH3F* hMeeMixlike2_Pt_Cent;

    TH2F* h_trackid_vs_pT;
    TH2F* h_trackid_vs_phi;

    THnF* hMee_Pt_phi_Cent;
    THnF* hMeelike1_Pt_phi_Cent;
    THnF* hMeelike2_Pt_phi_Cent;
    THnF* hMeeMix_Pt_phi_Cent;
    THnF* hMeeMixlike1_Pt_phi_Cent;
    THnF* hMeeMixlike2_Pt_phi_Cent;
    THnF* hMee_Pt_phi1_phi2_Cent;
    THnF* hMeelike1_Pt_phi1_phi2_Cent;
    //MRPC ToF
    TH1F* ModuleId_1;//1/Beta 0.8-0.9, P 0.4-3 GeV
    TH1F* TofId_1;
    TH2F* htofhits_Z_phi_1;
    TH1F* TofId_pionband;
    TH1F* TofId_zj;
    TH2F* h_Z_Phi_diff_cell;
    TH2F* h_Z_Phi_pionband;
    TH2F* h_Z_Phi_zj;
    TH2F* h_Z_Phi_zj_cal;// calculated by self
    TH2F* htofhitpos_1;//STAR btofHitPos
    TH3F* htofhit_xyz_1;
    TH2F* h_Z_Phi_cell8994_1;
    TH2F* h_Z_Phi_cell8998_1;
    TH2F* tofhit_Z_Phi_cell8998_1;
    TH2F* h_Z_Phi_cell8999_1;
    TH2F* h_Z_Phi_cell8498_1;

    TH2F* h_nSigmaPion_P;//tof problem cell calibration
    TH2F* h_nSigmaPion_P_01;//P < 1;1.2 < 1/beta <1.8
    TH1F* h_8994_t_calPion;
    TH1F* h_8994_t_tofPion;
    TH1F* h_8994_t_cal_tof;//cal time minus tof time
    TH1F* h_11917_t_calPion;
    TH1F* h_11917_t_tofPion;
    TH1F* h_11917_t_cal_tof;//cal time minus tof time
    TH1F* h_9001_t_calPion;
    TH1F* h_9001_t_tofPion;
    TH1F* h_9001_t_cal_tof;//cal time minus tof time
    float t_cal_Pion=0;
    float t_cal_Pion01=0;
    float t_tof_Pion=0;
    TH2F* h_cellid_time;
    TH2F* h_cellid_time_19086073;//run number is 19086073
    TH1F* h_tof_delta_cal;//calculation
    TH1F* h_tof_delta_tof;//btof
    TH1F* h_tof_delta_track;//track number
    float t_cal[23040];
    float t_tof[23040];
    int i_pion_track[23040];
    int iEvent_tof=0;
    int i_init=0;
    int i_tofcell=0;
    int nEvent_tof=0;
    int ustc=0;
    int ustc01=0;
    int ustc02=0;
    int ustc03=0;
    int ustc04=0;

    TH2F* hZphi_cell12510;
    TH2F* hZphi_cell12511;
    TH2F* hZphi_cell12512;
    TH2F* hZphi_cell12513;
    TH2F* hZphi_cell12514;
    TH2F* hZphi_cell12515;

    float R_rho;//Pt radius
    float R_tof = 200;//cm
    float Z_hit;//cm
    float phi_hit;//
    float tofhit_x;//STAR btofHitPos
    float tofhit_y;
    float tofhit_z;

    TH1F* TrayId_1;
    TH1F* ModuleId_2;//1/Beta 0.82-0.9, P 0.4-3 GeV
    TH1F* TofId_2;
    TH1F* ModuleId_3;//1/Beta 0.82-0.88, P 0.4-3 GeV
    TH1F* TofId_3;
    TH1F* ModuleId_4;//1/Beta 0.84-0.88, P 0.4-3 GeV
    TH1F* TofId_4;
    TH1F* ModuleId_5;//1/Beta 0.84-0.9, P 0.4-3 GeV
    TH1F* TofId_5;

    TH2F* hULMvsPt;
	TH3F* hULMvsPtCen;


    //Run by run QA
    bool mRunbyRunQA;
    TProfile* pVpdVz;
    TProfile* pVzVpdVz;
    TProfile* pRefmult;
    TProfile* pVpdHitEast;
    TProfile* pVpdHitWest;
    TProfile* pVx;
    TProfile* pVy;
    TProfile* pVz;
    TProfile* pVr;
    //Run by run track level
    TProfile* pTof;
    TProfile* pDedx;
    TProfile* pRMSDedx;
    TProfile* pgDCA;
    TProfile* pgPt;
    TProfile* pgPhi;
    TProfile* pgEta;
    TProfile* pNFits;
    TProfile* ppPt;
    TProfile* ppEta;
    TProfile* ppPhi;
    ClassDef(StPicoDstarMixedMaker, 1)
};

inline void StPicoDstarMixedMaker::getBadruns(string inputFileName){
    ifstream fin(inputFileName.c_str());
    if(!fin){
      cout <<"no Bad runs list" << endl;
      return;
    }
    cout << "  " << inputFileName << flush;

    Int_t runId = 0 ;
    while( fin >> runId ) {
      mBadRun.push_back(runId);
    }
    cout << "get Bad runs list [OK]" << endl;
}
inline  bool StPicoDstarMixedMaker::isBadrun(Int_t runId){
    vector<Int_t>::iterator iter = std::find(mBadRun.begin(), mBadRun.end(), runId);
    return ( iter != mBadRun.end() ) ;
}


inline void StPicoDstarMixedMaker::getBadcells(string inputFileName){
    ifstream fin(inputFileName.c_str());
    if(!fin){
      cout <<"no Bad cells list" << endl;
      return;
    }
    cout << "  " << inputFileName << flush;

    Int_t cellId = 0 ;
    while( fin >> cellId ) {
      mBadCell.push_back(cellId);
    }
    cout << "get Bad cells list [OK]" << endl;
}
inline  bool StPicoDstarMixedMaker::isBadcell(Int_t cellId){
    vector<Int_t>::iterator iter = std::find(mBadCell.begin(), mBadCell.end(), cellId);
    return ( iter != mBadCell.end() ) ;
}

#endif
