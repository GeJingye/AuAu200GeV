// ****************************************************
// *                                                  *
// *  Authors: Yuanjing Ji                            *
// *           Guannan Xie <guannanxie@lbl.gov>       *
// *           Mustafa Mustafa <mmustafa@lbl.gov>     *
// *                                                  *
// ****************************************************
 /* ***********************************************************************
  * Extract pure electron to recalibrating nsigama_e vs pT,eta,phi vs Cen *
  * ***********************************************************************
  */
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"

#include "StEvent/StDcaGeometry.h"
#include "StPhysicalHelixD.hh"
#include "phys_constants.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoETofPidTraits.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"
#include "StPicoEvent/StPicoMtdPidTraits.h"
#include "StPicoDstarMixedMaker.h"
#include "StAnaCuts.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"

#ifndef DEBUG
#define DEBUG 1
#endif

/*
命名规则
m = member → 类的变量成员
s = static → 类的静态成员
g = global → 全局变量（极少用）
k = const/constexpr → 编译期常量（如 kStOK）
h = histogram →  ROOT 直方图指针
*/

/*
ClassImp() 是 ROOT 框架中的一个宏，用于向 ROOT 系统注册类，使得该类：
支持 ROOT 的 I / O 系统（如 TFile::Write()、TFile::Read()）
支持 反射机制（如 TClass::GetClass("ClassName")）
支持 交互式命令行（CLING）识别
支持 TTree::Branch()自动序列化
*/
ClassImp(StPicoDstarMixedMaker)

// 成员初始化函数的类外实现
StPicoDstarMixedMaker::StPicoDstarMixedMaker(Char_t const *name, TString const inputFilesList, TString const outFileBaseName, StPicoDstMaker* picoDstMaker):
    StMaker(name),// char*实际上是指向存储字符串的字符数组第一个元素的地址
	mPicoDstMaker(picoDstMaker),
    mInputFilesList(inputFilesList),
	mOutFileBaseName(outFileBaseName)
{}
StPicoDstarMixedMaker::~StPicoDstarMixedMaker() {}

Int_t StPicoDstarMixedMaker::Init()
{
  mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");
  mFile = new TFile(mOutFileBaseName+".root", "RECREATE");

  // initialize histograms and trees
  initHists();
  return kStOK;
}

void StPicoDstarMixedMaker::initHists()
{
	ifstream readnum;
	readnum.open(mRunNumList);// mRunNumList由readPico.cxx通过setRunNumList函数赋值，即runnumber.list
	Int_t runum = 0;
	Int_t totalNum = 0;
	//if (DEBUG) cout << "start initial run number..." << endl;
	while (readnum >> runum) {
		mrunnum[runum] = totalNum;
		totalNum++;
		//if (DEBUG) cout << "run number : " << runum << " id :" << mrunnum[runum] << endl;
	}
	readnum.close();

	mRefMultCorrUtil = CentralityMaker::instance()->getRefMult6Corr();
	mRefMultCorrUtil->setVerbose(kFALSE);

	// Event histograms
	fphiVcut = new TF1("fphiVcut", "0.84326*exp(-49.4819*x)-0.996609*x+0.19801", 0.0, 1.0);
	h_RunNum__gTri = new TH1F("h_RunNum__gTri", "h_RunNum__gTri", totalNum, -0.5, totalNum-0.5);
	h_cen9 = new TH1F("h_cen9", "h_cen9", 9, 0., 9.);
	h_cen9__reWeight = new TH1F("h_cen9__reWeight", "h_cen9__reWeight", 9, 0., 9.);
	h_cen9->GetXaxis()->SetBinLabel(9, "0~5%");
	h_cen9->GetXaxis()->SetBinLabel(8, "5~10%");
	h_cen9->GetXaxis()->SetBinLabel(7, "10~20%");
	h_cen9->GetXaxis()->SetBinLabel(6, "20~30%");
	h_cen9->GetXaxis()->SetBinLabel(5, "30~40%");
	h_cen9->GetXaxis()->SetBinLabel(4, "40~50%");
	h_cen9->GetXaxis()->SetBinLabel(3, "50~60%");
	h_cen9->GetXaxis()->SetBinLabel(2, "60~70%");
	h_cen9->GetXaxis()->SetBinLabel(1, "70~80%");
	h_RefMult = new TH1F("h_RefMult", "h_RefMult", 250, 0, 250);// 参考多重数
	h_nTofMat_RefMul = new TH2F("h_nTofMat_RefMul", "RefMul VS nTofmatch;RefMul;nTofMatch", 250, 0, 250, 250, 0, 250);// 与TOF匹配的径迹数vs参考多重数关系
	h_RefMult6 = new TH1F("h_RefMult6", "h_RefMult6", 250, 0, 250);// 监控事件多重数分布->计算中心度
	h_TotnMIP = new TH1F("h_TotnMIP", "h_TotnMIP", 1000, 0, 1000);// 最小电离粒子（MIP）总数
	h_passEvtcut = new TH1F("h_passEvtcut", "pass event cut", 8, -0.5, 7.5);
	h_passEvtcut->GetXaxis()->SetBinLabel(1, "All");
	h_passEvtcut->GetXaxis()->SetBinLabel(2, "Good Runs");
	h_passEvtcut->GetXaxis()->SetBinLabel(3, "passVz");
	h_passEvtcut->GetXaxis()->SetBinLabel(4, "passVerr");
	h_passEvtcut->GetXaxis()->SetBinLabel(5, "passVr");
	h_passEvtcut->GetXaxis()->SetBinLabel(6, "notPileUp");
	h_passEvtcut->GetXaxis()->SetBinLabel(7, "0-80%");
	h_passEvtcut->GetXaxis()->SetBinLabel(8, "pass VpdVz-Vz");

	h_passTrkcut = new TH1F("h_passTrkcut", "tracks in different conditions", 3, -0.5, 2.5);
	h_passTrkcut->GetXaxis()->SetBinLabel(1, "All");
	h_passTrkcut->GetXaxis()->SetBinLabel(2, "Primary Tracks");
	h_passTrkcut->GetXaxis()->SetBinLabel(3, "Good Tracks");

	//顶点位置
	h_Vx_Vy_Vz = new TH3F("h_Vx_Vy_Vz", "Vz vs Vy vs Vx;Vx(cm);Vy(cm);Vz(cm)", 250, -5, 5, 250, -5, 5, 400, -200, 200);
	h_Vx_Vy = new TH2F("h_Vx_Vy", "Vy vs Vx;Vx(cm);Vy(cm)", 1400, -7, 7, 1400, -7, 7);
	h_Vr = new TH1F("h_Vr", "Vr;Vr(cm);Counts", 400, 0, 4);
	h_Vz = new TH1F("h_Vz", "Vz;Vz(cm);Counts", 400, -200, 200);
	h_VpdVz = new TH1F("h_VpdVz", "VpdVz;VpdVz(cm);Counts", 400, -200, 200);
	h_VpdVzmVz = new TH1F("h_VpdVzmVz", "VpdVz-Vz;VpdVz-Vz(cm);Counts", 2000, -200, 200);
	h_VpdVz_Vz = new TH2F("h_VpdVz_Vz", "Vz vs VpdVz;Vz(cm);VpdVz(cm)", 400, -200, 200, 400, -200, 200);

	// 径迹信息&主径迹信息
	h_nHitsFit = new TH1F("h_nHitsFit", "nHitsFit;nHitsFit", 160, -80., 80.);
	h_nHitsPoss = new TH1F("h_nHitsPoss", "nHitsPoss;nHitsPoss", 160, -80., 80.);
	h_nHitsDEdx = new TH1F("h_nHitsDEdx", "nHitsDedx;nHitsDedx", 160, -80., 80.);
	h_nHitsFit_Pt_Eta = new TH3F("h_nHitsFit_Pt_Eta", "nHitsFit vs p_{T} vs #eta;p_{T} (GeV/c);#eta;nHitsFit", 240, 0., 12., 40, -2., 2., 80, 0., 80.);
	h_nHitsDEdx_Pt_Eta = new TH3F("h_nHitsDEdx_Pt_Eta", "nHitsDedx vs p_{T} vs #eta;p_{T} (GeV/c);#eta;nHitsDedx", 240, 0., 12., 40, -2., 2., 80, 0., 80.);

	h_pDca = new TH1F("h_pDca", "pDca;DCA;counts", 50, 0., 5.);// p代表primary，pDCA指主径迹与重建顶点的最小距离
	h_ppT = new TH1F("h_ppT", "primary p_{T};p_{T} (GeV/c);counts", 1000, 0., 10.);
	h_pP = new TH1F("h_pP", "primary p;p (GeV/c);counts", 1000, 0., 10.);
	h_pP_ppT = new TH2F("h_pP_ppT", "primary p vs primary p_{T};p (GeV/c);p_{T} (GeV/c)", 1000, 0., 10., 1000, 0, 10);
	h_gPt = new TH1F("h_gPt", "global p_{T} of track;global p_{T} (GeV/c);counts", 1000, 0., 10.);
	h_pEta = new TH1F("h_pEta", "primary #eta;#eta;counts", 500, -2.5, 2.5);
	h_pPhi = new TH1F("h_pPhi", "primary #phi;#phi;counts", 80, -4.0, 4.0);
	h_ppTc_pEta = new TH2F("h_ppTc_pEta", "p_{T}*q vs #eta;p_{T}*q (GeV/c);#eta", 400, -10, 10, 500, -2.5, 2.5);
	h_ppTc_pPhi = new TH2F("h_ppTc_pPhi", "p_{T}*q vs #phi;p_{T}*q (GeV/c);#phi", 400, -10, 10, 800, -4, 4);
	h_pDca_Eta_NHitsFit = new TH3F("h_pDca_Eta_NHitsFit", "NHitsFit vs #eta vs DCA;DCA;#eta;NHitsFit", 50, 0., 5., 40, -2., 2., 90, 0, 90);
	h_pDca_Pt_Eta = new TH3F("h_pDca_Pt_Eta", "#eta vs p_{T} vs DCA;DCA;p_{T};#eta", 50, 0., 5.,1000 ,0., 10., 40, -2., 2.);
	//GoodTrack径迹信息
	h_nSigmaE_P__total = new TH2F("h_nSigmaE_P__total", "n#sigma_{e} vs p;p (GeV/c)n#sigma_{e}", 500, 0, 5, 2000, -10, 10);
	h_nSigmaPi_P__total = new TH2F("h_nSigmaPi_P__total", "n#sigma_{pi} vs p;p (GeV/c);n#sigma_{pi}", 500, 0, 5, 2000, -10, 10);
	h_nSigmaK_P__total = new TH2F("h_nSigmaK_P__total", "n#sigma_{K} vs p;p (GeV/c);n#sigma_{K}", 500, 0, 5, 2000, -10, 10);
	h_nSigmaP_P__total = new TH2F("h_nSigmaP_P__total", "n#sigma_{P} vs p;p (GeV/c);n#sigma_{P}", 500, 0, 5, 2000, -10, 10);

	//Recalibrate
	h_Pt_Cen_nSigmaE__PureE = new TH3F("h_Pt_Cen_nSigmaE__PureE", "p_{T} vs Cen vs n#sigma_{e};p_{T};Cen;n#sigma_{e}", 50, 0, 5, 8, 0, 8, 200, -10, 10);
	h_Eta_Cen_nSigmaE__PureE = new TH3F("h_Eta_Cen_nSigmaE__PureE", "#eta vs Cen vs n#sigma_{e};#eta;Cen;n#sigma_{e}", 40, -2, 2, 8, 0, 8, 200, -10, 10);
	h_Phi_Cen_nSigmaE__PureE = new TH3F("h_Phi_Cen_nSigmaE__PureE", "#phi vs Cen vs n#sigma_{e};#phi;Cen;n#sigma_{e}", 64, -3.2, 3.2, 8, 0, 8, 200, -10, 10);
	h_Pt_Eta_Phi__PureE = new TH3F("h_Pt_Eta_Phi__PureE", "p_{T} vs #eta vs #phi;p_{T};#eta;#phi", 50, 0, 5, 40, -2, 2, 64, -3.2, 3.2);

	h_Pt_Cen_nSigmaPi__PurePi = new TH3F("h_Pt_Cen_nSigmaPi__PurePi", "p_{T} vs Cen vs n#sigma_{pi};p_{T};Cen;n#sigma_{pi}", 50, 0, 5, 8, 0, 8, 200, -10, 10);
	h_Eta_Cen_nSigmaPi__PurePi = new TH3F("h_Eta_Cen_nSigmaPi__PurePi", "#eta vs Cen vs n#sigma_{pi};#eta;Cen;n#sigma_{pi}", 40, -2, 2, 8, 0, 8, 200, -10, 10);
	h_Phi_Cen_nSigmaPi__PurePi = new TH3F("h_Phi_Cen_nSigmaPi__PurePi", "#phi vs Cen vs n#sigma_{pi};#phi;Cen;n#sigma_{pi}", 64, -3.2, 3.2, 8, 0, 8, 200, -10, 10);
	h_Pt_Eta_Phi__PurePi = new TH3F("h_Pt_Eta_Phi__PurePi", "p_{T} vs #eta vs #phi;p_{T};#eta;#phi", 50, 0, 5, 40, -2, 2, 64, -3.2, 3.2);

	h_nSigmaElectron_P__1 = new TH2F("h_nSigmaElectron_P__1", "n#sigma_{e} vs p (p_{T}>0.2,|#eta|<0.9);p (GeV/c);n#sigma_{e}", 500, 0, 5, 2000, -10, 10);
	h_nSigmaElectron_P__2 = new TH2F("h_nSigmaElectron_P__2", "n#sigma_{e} vs p (p_{T}>0.2,0.9<|#eta|<1.8);p (GeV/c);n#sigma_{e}", 500, 0, 5, 2000, -10, 10);
	h_nSigmaElectron_P__3 = new TH2F("h_nSigmaElectron_P__3", "n#sigma_{e} vs p (p_{T}<0.2,|#eta|<1.8);p (GeV/c);n#sigma_{e}", 500, 0, 5, 2000, -10, 10);
	h_nSigmaElectron_P__TOFMatch = new TH2F("h_nSigmaElectron_P__TOFMatch", "n#sigma_{e} vs p (p_{T}>0.2,|#eta|<0.9);p (GeV/c);n#sigma_{e}", 500, 0, 5, 2000, -10, 10);
	h_nSigmaElectron_P__PID_2 = new TH2F("h_nSigmaElectron_P__PID_2",		"n#sigma_{e} vs p (p_{T}>0.2,0.9<|#eta|<1.8);p (GeV/c);n#sigma_{e}", 500, 0, 5, 2000, -10, 10);
	h_nSigmaElectron_P__PID_3 = new TH2F("h_nSigmaElectron_P__PID_3",		"n#sigma_{e} vs p (p_{T}<0.2,|#eta|<1.8);p (GeV/c);n#sigma_{e}", 500, 0, 5, 2000, -10, 10);
	h_nSigmaElectron_P__EIDcut_1 = new TH2F("h_nSigmaElectron_P__EIDcut_1", "n#sigma_{e} vs p (p_{T}>0.2,|#eta|<0.9);p (GeV/c);n#sigma_{e}", 500, 0, 5, 2000, -10, 10);
	h_nSigmaElectron_P__EIDcut_2 = new TH2F("h_nSigmaElectron_P__EIDcut_2", "n#sigma_{e} vs p (p_{T}>0.2,0.9<|#eta|<1.8);p (GeV/c);n#sigma_{e}", 500, 0, 5, 2000, -10, 10);
	h_nSigmaElectron_Eta__EIDcut_2 = new TH2F("h_nSigmaElectron_Eta__EIDcut_2", "n#sigma_{e} vs #eta (p_{T}>0.2,0.9<|#eta|<1.8);#eta;n#sigma_{e}", 400, -2, 2, 2000, -10, 10);
	h_nSigmaElectron_P__EIDcut_3 = new TH2F("h_nSigmaElectron_P__EIDcut_3", "n#sigma_{e} vs p (p_{T}<0.2,|#eta|<1.8);p (GeV/c);n#sigma_{e}", 500, 0, 5, 2000, -10, 10);
	h_nSigmaElectron_Eta__EIDcut_3_lowP = new TH2F("h_nSigmaElectron_Eta__EIDcut_3_lowP", "n#sigma_{e} vs #eta (p_{T}<0.2,low p,|#eta|<1.8);#eta;n#sigma_{e}", 400, -2.0, 2.0, 2000, -10, 10);
	h_nSigmaElectron_Eta__EIDcut_3_highP = new TH2F("h_nSigmaElectron_Eta__EIDcut_3_highP", "n#sigma_{e} vs #eta (p_{T}<0.2,high p,|#eta|<1.8);#eta;n#sigma_{e}", 400, -2.0, 2.0, 2000, -10, 10);
	h_nSigmaElectron_P__EIDcut_total = new TH2F("h_nSigmaElectron_P__EIDcut_total", "n#sigma_{e} vs p (p_{T}>0.06,0.0<|#eta|<1.8);p (GeV/c);n#sigma_{e}", 500, 0, 5, 2000, -10, 10);

	h_nSigmaPion_P__2   = new TH2F("h_nSigmaPion_P__2", "n#sigma_{pi} vs p (p_{T}>0.2,0.9<|#eta|<1.8);p (GeV/c);n#sigma_{pi}", 500, 0, 5, 2000, -10, 10);
	h_nSigmaKaon_P__2   = new TH2F("h_nSigmaKaon_P__2",  "n#sigma_{K} vs p (p_{T}>0.2,0.9<|#eta|<1.8);p (GeV/c);n#sigma_{K}", 500, 0, 5, 2000, -10, 10);
	h_nSigmaProton_P__2 = new TH2F("h_nSigmaProton_P__2","n#sigma_{P} vs p (p_{T}>0.2,0.9<|#eta|<1.8);p (GeV/c);n#sigma_{P}", 500, 0, 5, 2000, -10, 10);
	h_nSigmaPion_P__3   = new TH2F("h_nSigmaPion_P__3",  "n#sigma_{pi} vs p (p_{T}<0.2,|#eta|<1.8);p (GeV/c);n#sigma_{pi}", 500, 0, 5, 2000, -10, 10);
	h_nSigmaKaon_P__3   = new TH2F("h_nSigmaKaon_P__3",  "n#sigma_{K} vs p (p_{T}<0.2,|#eta|<1.8);p (GeV/c);n#sigma_{K}", 500, 0, 5, 2000, -10, 10);

	// TOF匹配的径迹信息
	h_invBeta_P__TOFMatch = new TH2F("h_invBeta_P__TOFMatch", "1/#beta vs p;p (GeV/c);1/#beta", 500, 0, 5, 5000, 0, 5);
	h_pT__TOFMatch  = new TH1F("h_pT__TOFMatch", "p_{T} of TOF matched tracks;p_{T} (GeV/c);counts", 1000, 0., 10.);
	h_Eta__TOFMatch = new TH1F("h_Eta__TOFMatch", "#eta of TOF matched tracks;#eta;counts", 500, -2.5, 2.5);
	h_Phi__TOFMatch = new TH1F("h_Phi__TOFMatch", "#phi of TOF matched tracks;#phi;counts", 80, -4, 4);
	h_nSigmaEcorr_P__TOFMatch = new TH2F("h_nSigmaEcorr_P__TOFMatch", "n#sigma_{e} vs p (TOF acceptance);p (GeV/c);n#sigma_{e}", 500, 0, 5, 2000, -10, 10);
	// 经TPC和TOF判选后的电子信息
	h_Rapidity__unlikeSame = new TH1F("h_Rapidity__unlikeSame", "y distribution of e^{+}e^{-};y;counts", 800, -4, 4);
	h_Mee_PhiV__unlikeSame = new TH2F("h_Mee_PhiV__unlikeSame", "Mee vs #phi_{V};Mee(GeV/c^{2});#phi_{V}", 800, 0, 4, 100, 0, 1);
	h_Mee__unlikeSame = new TH1F("h_Mee__unlikeSame", "Mee without #phi_{V} cut;Mee(GeV/c^{2})", 800, 0, 4);
	h_Mee__unlikeSame__w_PhiV_Cut = new TH1F("h_Mee__unlikeSame__w_PhiV_Cut", "Mee with #phi_{V} cut;Mee(GeV/c^{2})", 800, 0, 4);

	h_Mee_Pt_Cen__unlikeSame = new TH3F("h_Mee_Pt_Cen__unlikeSame", "Mee vs p_{T} vs Cen;Mee (GeV/c^{2});p_{T} (GeV/c);Cen", 800, 0, 4, 100, 0, 5, 16, 0, 16);
	h_Mee_Pt_Cen__likemm = new TH3F("h_Mee_Pt_Cen__likemm", "Mee vs p_{T} vs Cen;Mee (GeV/c^{2});p_{T} (GeV/c);Cen", 800, 0, 4, 100, 0, 5, 16, 0, 16);
	h_Mee_Pt_Cen__likepp = new TH3F("h_Mee_Pt_Cen__likepp", "Mee vs p_{T} vs Cen;Mee (GeV/c^{2});p_{T} (GeV/c);Cen", 800, 0, 4, 100, 0, 5, 16, 0, 16);
}// 

Int_t StPicoDstarMixedMaker::Make()
{
  ParticleInfo particleinfo;
  vector<ParticleInfo> electroninfo;
  vector<ParticleInfo> positroninfo;

  // StMemStat mem;
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


  StPicoEvent const *picoEvent = picoDst->event();
  if (!isGoodTrigger(picoEvent)) return kStOK;

  // 清空正负电子信息缓存
  electroninfo.clear();
  positroninfo.clear();

  mRunId = picoEvent->runId();
  h_RunNum__gTri->Fill(mrunnum[mRunId]);// 填充不同run号的事例数；

  TVector3 pVtx = picoEvent->primaryVertex();// 获得TPC重建的该事例顶点的三维坐标
  mVx = pVtx.x();
  mVy = pVtx.y();
  mVz = pVtx.z();
  mVr = sqrt(mVx*mVx+mVy*mVy);
  mVpdVz = picoEvent->vzVpd();// 获得VPD测量的该事例顶点的z坐标
  
  // event and track level QA
  h_passEvtcut->Fill(0);// 原始事例数+1

  if (!isBadrun(mRunId))// bad run list 
  {
	h_passEvtcut->Fill(1);// 通过bad run cut的事例数+1
	h_Vx_Vy_Vz->Fill(mVx, mVy, mVz);
	h_Vz->Fill(mVz);
	h_Vx_Vy->Fill(mVx, mVy);
	h_Vr->Fill(mVr);

	if (fabs(mVpdVz + 999.0) > 1e-2 && fabs(mVpdVz) < 1000.0)// STAR约定：当 VPD 无法重建顶点时，会把 mVpdVz 设为 -999 cm 作为占位标志。
	{
		h_VpdVz->Fill(mVpdVz);
		h_VpdVz_Vz->Fill(mVz, mVpdVz);
		h_VpdVzmVz->Fill(mVpdVz - mVz);
	}

	// 获取多重数
	Int_t Refmult = picoEvent->refMult();//给出的是在线/原始 TPC 多重数（reference multiplicity），未经任何修正
	mRefmult6 = getRefmult6(picoDst, picoEvent);// 用户离线重算的多重数，满足严格几何与质量要求
	mTotnMIP = getTotnMIP(picoDst);// 获取EPD监测的该事例的最小电离粒子总数目

	h_RefMult->Fill(Refmult);
	h_nTofMat_RefMul->Fill(Refmult, picoEvent->nBTOFMatch());
    h_RefMult6->Fill(mRefmult6);
    h_TotnMIP->Fill(mTotnMIP);

	// 通过离线多重数计算该事例的中心度
	mRefMultCorrUtil->init(mRunId);
	mRefMultCorrUtil->initEvent(mRefmult6, mVz);
	Double_t reWeight = mRefMultCorrUtil->getWeight();
	//mCen9 = mRefMultCorrUtil->getCentralityBin9();
    mCen16 = mRefMultCorrUtil->getCentralityBin16();

	// 不同条件cut后的事例数统计
    Bool_t vzcut = mVz < anaCuts::Vz_up && mVz > anaCuts::Vz_low;
    Bool_t verrcut = !(fabs(mVx) < anaCuts::Verror && fabs(mVy) < anaCuts::Verror && fabs(mVz) < anaCuts::Verror);// Vx,Vy,Vz<1.0e-5 cm, why? too small that better than resolution. 
    Bool_t vrcut =  mVr < anaCuts::Vr;
    Bool_t notPileUp = !mRefMultCorrUtil->isPileUpEvent(mRefmult6, picoEvent->nBTOFMatch(), mVz, mTotnMIP);
    Bool_t cen0280cut = mCen16 > -1;
	Bool_t vzvpdvzcut = fabs(mVz - mVpdVz) < anaCuts::vzVpdVz;
	
    if (vzcut) h_passEvtcut->Fill(2);
    if (vzcut && vrcut) h_passEvtcut->Fill(3);
    if (vzcut && vrcut && verrcut ) h_passEvtcut->Fill(4);
    if (vzcut && vrcut && verrcut && notPileUp ) h_passEvtcut->Fill(5);
    if (vzcut && vrcut && verrcut && notPileUp && cen0280cut) h_passEvtcut->Fill(6);
	if (vzcut && vrcut && verrcut && notPileUp && cen0280cut && vzvpdvzcut) h_passEvtcut->Fill(7);

    if (isGoodEvent(picoEvent) && notPileUp)
    {
      mBfield  = picoEvent->bField();// 获取磁场
	  h_cen9->Fill(mCen16);// 填充中心度
	  h_cen9__reWeight->Fill(mCen16, reWeight);
      // ******************以下分析均基于0~80%中心度***********************
      if(!cen0280cut) return kStOK;

	  // ************************开始径迹分析*****************************
      Int_t nTracks = picoDst->numberOfTracks();
	  for (Int_t itrack = 0; itrack < nTracks; itrack++)
	  {
		  h_passTrkcut->Fill(0);
		  StPicoTrack* trk = picoDst->track(itrack);
		  // ******************以下分析均基于主径迹**************************
		  Bool_t isPrimaryTrack = trk->isPrimary();
		  if (!isPrimaryTrack) continue;
		  h_passTrkcut->Fill(1);
		  TVector3 mom = trk->pMom();// 提取当前径迹在“主顶点处”的三维动量矢量？
		  Float_t mgDCAs = trk->gDCA(picoEvent->primaryVertex()).Mag();// 返回该径迹与主顶点的最短距离

		  // 获得nSigmaElectron等图像
		  Double_t nSigmaE = trk->nSigmaElectron();
		  Double_t nSigmaPi = trk->nSigmaPion();
		  Double_t nSigmaK = trk->nSigmaKaon();
		  Double_t nSigmaP = trk->nSigmaProton();
		  Double_t beta = getTofBeta(trk);
		  Double_t m2 = pow(mom.Mag()*1.0 / beta, 2)*(1 - beta * beta);
		  Bool_t isTOFMatch = (beta != std::numeric_limits<Float_t>::quiet_NaN()) && beta > 0;
		  // 填充直方图
		  if (fabs(mgDCAs + 999.0) > 1e-2) h_pDca->Fill(mgDCAs);// STAR约定：当 DCA 无法计算（例如径迹缺 hit、拟合失败）时，把 mgDCAs 设为 - 999 cm 作为无效标志。
		  h_ppT->Fill(mom.Perp());
		  h_pP->Fill(mom.Mag());
		  h_pP_ppT->Fill(mom.Mag(), mom.Perp());
		  h_gPt->Fill(trk->gPt());
		  h_pEta->Fill(mom.Eta());
		  h_pPhi->Fill(mom.Phi());
		  h_ppTc_pEta->Fill(mom.Perp()*trk->charge(), mom.Eta());// mom.Perp()指横动量
		  h_ppTc_pPhi->Fill(mom.Perp()*trk->charge(), mom.Phi());

		  h_nHitsFit->Fill(trk->nHitsFit()*trk->charge());
		  h_nHitsPoss->Fill(trk->nHitsPoss()*trk->charge());
		  h_nHitsDEdx->Fill(trk->nHitsDedx()*trk->charge());
		  h_nHitsFit_Pt_Eta->Fill(mom.Perp(), mom.Eta(), trk->nHitsFit());
		  h_nHitsDEdx_Pt_Eta->Fill(mom.Perp(), mom.Eta(), trk->nHitsDedx());
		  h_pDca_Eta_NHitsFit->Fill(mgDCAs, mom.Eta(), trk->nHitsFit());
		  h_pDca_Pt_Eta->Fill(mgDCAs, mom.Perp(), mom.Eta());
		  // ******************以下分析均基于满足goodTrackCuts的主径迹**************************
		  Bool_t goodtrack = isGoodTrack(trk, picoEvent);
		  if (!goodtrack) continue;
		  h_passTrkcut->Fill(2);

		  h_nSigmaE_P__total->Fill(mom.Mag(),nSigmaE);
		  h_nSigmaPi_P__total->Fill(mom.Mag(),nSigmaPi);
		  h_nSigmaK_P__total->Fill(mom.Mag(),nSigmaK);
		  h_nSigmaP_P__total->Fill(mom.Mag(),nSigmaP);

		  // Particle IDentification
		  Bool_t inRegion__1 = mom.Perp() > 0.2 && fabs(mom.Eta()) < anaCuts::Eta;
		  Bool_t isTOFElectron__1 = kFALSE;
		  Bool_t isTPCElectron__1 = kFALSE;

		  Bool_t inRegion__2 = mom.Perp() > 0.2 && fabs(mom.Eta()) > anaCuts::Eta && fabs(mom.Eta()) < 1.8;
		  Bool_t isTPCPion__2 = kTRUE;
		  Bool_t isTPCKaon__2 = kTRUE;
		  Bool_t isTPCProton__2 = kTRUE;
		  Bool_t isTPCElectron__2 = kFALSE;

		  Bool_t inRegion__3 = mom.Perp() < 0.2 && fabs(mom.Eta()) < 1.8;
		  Bool_t isTPCPion__3 = kTRUE;
		  Bool_t isTPCKaon__3 = kTRUE;
		  Bool_t isTPCElectron__3 = kFALSE;
		  Bool_t isLowPElectron__3 = kFALSE;
		  Bool_t isLowEtaElectron__3 = kFALSE;

		  Bool_t isElectronRegion1 = kFALSE;
		  Bool_t isElectronRegion2 = kFALSE;
		  Bool_t isElectronRegion3 = kFALSE;

		  // TPC cut + TOF cut
		  // region1
		  if (isTOFMatch)
		  {
			  Int_t index2tof = trk->bTofPidTraitsIndex();
			  StPicoBTofPidTraits const* const tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);
			  Int_t tofid = tofPid->btofCellId();

			  if (isGoodBTofCell(tofid) && fabs(tofPid->btofYLocal()) < 1.8)
			  {
				  h_pT__TOFMatch->Fill(mom.Perp());
				  h_Eta__TOFMatch->Fill(mom.Eta());
				  h_Phi__TOFMatch->Fill(mom.Phi());
				  if (mom.Perp() > 0.2 && fabs(mom.Eta()) < anaCuts::Eta) h_invBeta_P__TOFMatch->Fill(mom.Mag(), 1. / beta);
				  if (fabs(1.0 / beta - 1) < 0.025)
				  {
					  isTOFElectron__1 = kTRUE;
				  }
			  }
		  }
		  if (mom.Mag() <= 1.0) isTPCElectron__1 = nSigmaE<3 && nSigmaE>(2.8*mom.Mag() - 4.0);
		  if (mom.Mag() > 1.0)  isTPCElectron__1 = nSigmaE<3 && nSigmaE>-1.2;
		  if (inRegion__1)
		  {
			  h_nSigmaElectron_P__1->Fill(mom.Mag(), nSigmaE);
			  if (isTOFElectron__1)
			  {
				  h_nSigmaElectron_P__TOFMatch->Fill(mom.Mag(), nSigmaE);
				  if (isTPCElectron__1)
				  {
					  h_nSigmaElectron_P__EIDcut_1->Fill(mom.Mag(), nSigmaE);
					  isElectronRegion1 = kTRUE;
				  }
			  }
		  }
		  // group 2 P24ia
		  //isTPCKaon__2 = mom.Mag() < 1.0 ? nSigmaK > -2.*mom.Mag() && nSigmaK<4.5 : nSigmaK>-2.0&&nSigmaK < 4.5;
		  //isTPCProton__2 = mom.Mag() < 1.5 ? nSigmaP > -1.33*mom.Mag() && nSigmaP<-1.33*mom.Mag() + 6 : nSigmaP>-2.&&nSigmaP < 4.;
		  //isTPCElectron__2 = mom.Mag() < 1.5 ? nSigmaE > 2.*mom.Mag() - 3.&&nSigmaE<3. : nSigmaE>0.&&nSigmaK < 3.;
		  // region2 P24iy
		  isTPCPion__2	   = nSigmaPi > -3.5 && nSigmaK < 3.;
		  isTPCKaon__2     = nSigmaK > -2.5 && nSigmaK < 4.0;
		  isTPCProton__2   = nSigmaP > -2.2 && nSigmaP < 4.2;
		  isTPCElectron__2 = nSigmaE > -1.5 && nSigmaE < 2.0;
		  if (inRegion__2)
		  {
			  h_nSigmaElectron_P__2->Fill(mom.Mag(), nSigmaE);
			  h_nSigmaPion_P__2->Fill(mom.Mag(), nSigmaPi);
			  h_nSigmaKaon_P__2->Fill(mom.Mag(), nSigmaK);
			  h_nSigmaProton_P__2->Fill(mom.Mag(), nSigmaP);
			  if (!isTPCKaon__2 && !isTPCProton__2)
			  {
				  h_nSigmaElectron_P__PID_2->Fill(mom.Mag(), nSigmaE);
				  if (isTPCElectron__2)
				  {
					  h_nSigmaElectron_P__EIDcut_2->Fill(mom.Mag(), nSigmaE);
					  h_nSigmaElectron_Eta__EIDcut_2->Fill(mom.Eta(), nSigmaE);
					  isElectronRegion2 = kTRUE;
				  }
			  }
		  }
		  //// group 3 P24ia
		  //isTPCPion__3 = mom.Mag() < 0.2 ? nSigmaPi > -20.*mom.Mag() - 1.0&&nSigmaPi<5. : nSigmaPi>-5.&&nSigmaPi < 5.;
		  //isTPCKaon__3 = nSigmaK > 13.57*mom.Mag() - 6.0 && nSigmaK < 8.0;
		  //if (mom.Mag() > 0.0 && mom.Mag() < 0.14)  isTPCElectron__3 = nSigmaE > -4.0	&& nSigmaE < 1.0;
		  //if (mom.Mag() > 0.14 && mom.Mag() < 0.18) isTPCElectron__3 = nSigmaE > -1.0	&& nSigmaE < 1.0;
		  //if (mom.Mag() > 0.18 && mom.Mag() < 0.7)  isTPCElectron__3 = nSigmaE > -1.0	&& nSigmaE < 4.0;
		  //// region3 P24iy
		  if (mom.Mag() > 0.0 && mom.Mag() <= 0.2) isTPCPion__3 = nSigmaPi > -20.0 * mom.Mag() - 1.0 && nSigmaPi < -5.71 * mom.Mag() + 5.0;
		  if (mom.Mag() > 0.2 && mom.Mag() < 0.7)  isTPCPion__3 = nSigmaPi > -5.0 && nSigmaPi < -5.71 * mom.Mag() + 5.0;
		  if (mom.Mag() > 0.0 && mom.Mag() < 0.7)  isTPCKaon__3 = nSigmaK > 13.57 * mom.Mag() - 6.0 && nSigmaK < 5.0;
		  if (mom.Mag() > 0.0 && mom.Mag() < 0.14)  isTPCElectron__3 = nSigmaE > -4.0	&& nSigmaE < 1.0;
		  if (mom.Mag() > 0.14 && mom.Mag() < 0.18) isTPCElectron__3 = nSigmaE > -2.0	&& nSigmaE < 1.0;
		  if (mom.Mag() > 0.18 && mom.Mag() < 0.7)  isTPCElectron__3 = nSigmaE > -2.0	&& nSigmaE < 3.0;
		  if (inRegion__3)
		  {
			  h_nSigmaElectron_P__3->Fill(mom.Mag(), nSigmaE);
			  h_nSigmaPion_P__3->Fill(mom.Mag(), nSigmaPi);
			  h_nSigmaKaon_P__3->Fill(mom.Mag(), nSigmaK);
			  if (!isTPCPion__3 && !isTPCKaon__3)
			  {
				  h_nSigmaElectron_P__PID_3->Fill(mom.Mag(), nSigmaE);
				  if (isTPCElectron__3)
				  {
					  isLowPElectron__3 = nSigmaE < -75 * mom.Mag() + 12.5;
					  isLowEtaElectron__3 = fabs(mom.Eta()) < 0.1;
					  h_nSigmaElectron_P__EIDcut_3->Fill(mom.Mag(), nSigmaE);
					  isElectronRegion3 = kTRUE;
					  if (isLowPElectron__3)
					  {
						  h_nSigmaElectron_Eta__EIDcut_3_lowP->Fill(mom.Eta(), nSigmaE);
					  }
					  if (!isLowPElectron__3)
					  {
						  h_nSigmaElectron_Eta__EIDcut_3_highP->Fill(mom.Eta(), nSigmaE);

					  }
				  }
			  }
		  }
		  
		  // Pure pion
		  if ((inRegion__1 && nSigmaPi>-3.5&&nSigmaPi<3.5 && m2>0.0&&m2<0.04) || (inRegion__2&&isTPCPion__2) || (inRegion__3&&isTPCPion__3))
		  {
			h_Pt_Cen_nSigmaPi__PurePi->Fill(mom.Perp(), mCen16, nSigmaPi, reWeight);
			h_Eta_Cen_nSigmaPi__PurePi->Fill(mom.Eta(), mCen16, nSigmaPi, reWeight);
			h_Phi_Cen_nSigmaPi__PurePi->Fill(mom.Phi(), mCen16, nSigmaPi, reWeight);
			h_Pt_Eta_Phi__PurePi->Fill(mom.Perp(), mom.Eta(), mom.Phi());
		  }
		  // Pure kaon
		  //if ((inRegion__1 && nSigmaK>-3.5&&nSigmaK<3.5 && m2>0.2&&m2<0.3) || (inRegion__2&&isTPCKaon__2) || (inRegion__3&&isTPCKaon__3))
		  //{
			 // h_Pt_Cen_nSigmaK__PureK->Fill(mom.Perp(), mCen16, nSigmaK, reWeight);
			 // h_Eta_Cen_nSigmaK__PureK->Fill(mom.Eta(), mCen16, nSigmaK, reWeight);
			 // h_Phi_Cen_nSigmaK__PureK->Fill(mom.Phi(), mCen16, nSigmaK, reWeight);
			 // h_Pt_Eta_Phi__PureK->Fill(mom.Perp(), mom.Eta(), mom.Phi());
		  //}
		  //// Pure proton
		  //if ((inRegion__1 && nSigmaP>-3.5&&nSigmaP<3.5 && m2>0.80&&m2<0.96) || (inRegion__2&&isTPCProton__2) || (inRegion__3&&isTPCProton__3))
		  //{
			 // h_Pt_Cen_nSigmaP__PureP->Fill(mom.Perp(), mCen16, nSigmaP, reWeight);
			 // h_Eta_Cen_nSigmaP__PureP->Fill(mom.Eta(), mCen16, nSigmaP, reWeight);
			 // h_Phi_Cen_nSigmaP__PureP->Fill(mom.Phi(), mCen16, nSigmaP, reWeight);
			 // h_Pt_Eta_Phi__PureP->Fill(mom.Perp(), mom.Eta(), mom.Phi());
		  //}
		  // Pure electron
		  if (isElectronRegion1 || isElectronRegion2 || (isElectronRegion3 && !isLowPElectron__3))
		  {
			h_nSigmaElectron_P__EIDcut_total->Fill(mom.Mag(), nSigmaE);

			if (trk->charge() < 0)// electron
			{
				particleinfo.charge = trk->charge();
				particleinfo.pt = mom.Perp();
				particleinfo.eta = mom.Eta();
				particleinfo.phi = mom.Phi();
				particleinfo.p = mom.Mag();
				particleinfo.nSigmaE = nSigmaE;
				particleinfo.beta = beta;
				particleinfo.energy = sqrt(pow(M_electron, 2.0) + pow(mom.Mag(), 2.0));
				particleinfo.p1 = mom.X();
				particleinfo.p2 = mom.Y();
				particleinfo.p3 = mom.Z();
				particleinfo.isPhotonicE = kFALSE;
				particleinfo.isPureE = kFALSE;
				electroninfo.push_back(particleinfo);
			}
			if (trk->charge() > 0) // positron
			{
				particleinfo.charge = trk->charge();
				particleinfo.pt = mom.Perp();
				particleinfo.eta = mom.Eta();
				particleinfo.phi = mom.Phi();
				particleinfo.p = mom.Mag();
				particleinfo.nSigmaE = nSigmaE;
				particleinfo.beta = beta;
				particleinfo.energy = sqrt(pow(M_electron, 2.0) + pow(mom.Mag(), 2.0));
				particleinfo.p1 = mom.X();
				particleinfo.p2 = mom.Y();
				particleinfo.p3 = mom.Z();
				particleinfo.isPhotonicE = kFALSE;
				particleinfo.isPureE = kFALSE;
				positroninfo.push_back(particleinfo);
			}
          }// 填充单径迹的正/负电子信息
      }// 填充单事例所有径迹的正负电子信息

      Int_t x=0;
      Int_t y=0;
      Int_t num_electron = electroninfo.size();
      Int_t num_positron = positroninfo.size();
      TLorentzVector eepair(0,0,0,0);
      TLorentzVector particle1_4V(0,0,0,0);
      TLorentzVector particle2_4V(0,0,0,0);
      
	  // pure electron sample selection
	  for (x = 0; x < num_positron; x++)
	  {
		  particle1_4V.SetPx(positroninfo[x].p1);
		  particle1_4V.SetPy(positroninfo[x].p2);
		  particle1_4V.SetPz(positroninfo[x].p3);
		  particle1_4V.SetE(positroninfo[x].energy);
		  for (y = 0; y < num_electron; y++)
		  {
			  particle2_4V.SetPx(electroninfo[y].p1);
			  particle2_4V.SetPy(electroninfo[y].p2);
			  particle2_4V.SetPz(electroninfo[y].p3);
			  particle2_4V.SetE(electroninfo[y].energy);

			  eepair = particle1_4V + particle2_4V;
			  if (eepair.M() < 0.015)
			  {
				  positroninfo[x].isPureE = kTRUE;
				  electroninfo[y].isPureE = kTRUE;
			  }
		  }
	  }
	  for (x = 0; x < num_positron; x++)//填充nsigma_e vs pT,eta,phi vs Cens
	  {
		  if (positroninfo[x].isPureE)
		  {
			  h_Pt_Cen_nSigmaE__PureE->Fill(positroninfo[x].pt, mCen16, positroninfo[x].nSigmaE, reWeight);
			  h_Eta_Cen_nSigmaE__PureE->Fill(positroninfo[x].eta, mCen16, positroninfo[x].nSigmaE, reWeight);
			  h_Phi_Cen_nSigmaE__PureE->Fill(positroninfo[x].phi, mCen16, positroninfo[x].nSigmaE, reWeight);
			  h_Pt_Eta_Phi__PureE->Fill(positroninfo[x].pt, positroninfo[x].eta, positroninfo[x].phi);
		  }
	  }
	  for (x = 0; x < num_electron; x++)//填充nsigma_e vs pT,eta,phi vs Cen
	  {
		  if (electroninfo[x].isPureE)
		  {
			  h_Pt_Cen_nSigmaE__PureE->Fill(positroninfo[x].pt, mCen16, positroninfo[x].nSigmaE, reWeight);
			  h_Eta_Cen_nSigmaE__PureE->Fill(positroninfo[x].eta, mCen16, positroninfo[x].nSigmaE, reWeight);
			  h_Phi_Cen_nSigmaE__PureE->Fill(positroninfo[x].phi, mCen16, positroninfo[x].nSigmaE, reWeight);
			  h_Pt_Eta_Phi__PureE->Fill(positroninfo[x].pt, positroninfo[x].eta, positroninfo[x].phi);
		  }
	  }

	  // 通过随机组合++, --, +-电子对重建信号
	  // +-
	  for (x = 0; x < num_positron; x++)
	  {
		  particle1_4V.SetPx(positroninfo[x].p1);
		  particle1_4V.SetPy(positroninfo[x].p2);
		  particle1_4V.SetPz(positroninfo[x].p3);
		  particle1_4V.SetE(positroninfo[x].energy);
		  for (y = 0; y < num_electron; y++)
		  {
			  if (positroninfo[x].isPureE && electroninfo[y].isPureE)
			  {
				  particle2_4V.SetPx(electroninfo[y].p1);
				  particle2_4V.SetPy(electroninfo[y].p2);
				  particle2_4V.SetPz(electroninfo[y].p3);
				  particle2_4V.SetE(electroninfo[y].energy);
				  eepair = particle1_4V + particle2_4V;

				  h_Mee_Pt_Cen__unlikeSame->Fill(eepair.M(), eepair.Perp(), mCen16, reWeight);
			  }
		  }
	  }// end: +-
	  // --
	  for (x = 0; x < num_electron; x++)
	  {
		  particle1_4V.SetPx(electroninfo[x].p1);
		  particle1_4V.SetPy(electroninfo[x].p2);
		  particle1_4V.SetPz(electroninfo[x].p3);
		  particle1_4V.SetE(electroninfo[x].energy);
		  for (y = x + 1; y < num_electron; y++)//从x+1开始，避免自组合和重复组合
		  {
			  if (electroninfo[x].isPureE && electroninfo[y].isPureE)
			  {
				  particle2_4V.SetPx(electroninfo[y].p1);
				  particle2_4V.SetPy(electroninfo[y].p2);
				  particle2_4V.SetPz(electroninfo[y].p3);
				  particle2_4V.SetE(electroninfo[y].energy);
				  eepair = particle1_4V + particle2_4V;

				  h_Mee_Pt_Cen__likemm->Fill(eepair.M(), eepair.Perp(), mCen16, reWeight);
			  }
		  }
	  }// end: --
	  // ++
	  for (x = 0; x < num_positron; x++)
	  {
		  particle1_4V.SetPx(positroninfo[x].p1);
		  particle1_4V.SetPy(positroninfo[x].p2);
		  particle1_4V.SetPz(positroninfo[x].p3);
		  particle1_4V.SetE(positroninfo[x].energy);
		  for (y = x + 1; y < num_positron; y++)
		  {
			  if (positroninfo[x].isPureE && positroninfo[y].isPureE)
			  {
				  particle2_4V.SetPx(positroninfo[y].p1);
				  particle2_4V.SetPy(positroninfo[y].p2);
				  particle2_4V.SetPz(positroninfo[y].p3);
				  particle2_4V.SetE(positroninfo[y].energy);
				  eepair = particle1_4V + particle2_4V;

				  h_Mee_Pt_Cen__likepp->Fill(eepair.M(), eepair.Perp(), mCen16, reWeight);
			  }
		  }
	  }// end: ++
    }// Good Event
  }// Good Run

  return kStOK;
}// end Make
Int_t StPicoDstarMixedMaker::Finish()
{
	mFile->cd();// 将当前 ROOT 的输出目录切换到 mFile 所指向的 ROOT 文件

	h_RunNum__gTri->Write();
	h_passEvtcut->Write();
	h_passTrkcut->Write();
	h_cen9->Write();
	h_cen9__reWeight->Write();
	h_RefMult6->Write();
	h_TotnMIP->Write();

	// write the hists  
	// event QA
	h_Vx_Vy_Vz->Write();
	h_Vz->Write();
	h_VpdVz->Write();
	h_Vr->Write();
	h_VpdVz_Vz->Write();
	h_VpdVzmVz->Write();
	h_Vx_Vy->Write();
	h_nTofMat_RefMul->Write();
	h_RefMult->Write();
	h_pDca_Eta_NHitsFit->Write();
	h_pDca_Pt_Eta->Write();
	// track level QA
	h_nHitsFit->Write();
	h_nHitsPoss->Write();
	h_nHitsDEdx->Write();
	h_nHitsFit_Pt_Eta->Write();
	h_nHitsDEdx_Pt_Eta->Write();
	h_invBeta_P__TOFMatch->Write();
	h_pDca->Write();
	h_ppT->Write();
	h_pP->Write();
	h_pP_ppT->Write();
	h_gPt->Write();
	h_pEta->Write();
	h_pPhi->Write();
	h_ppTc_pEta->Write();
	h_ppTc_pPhi->Write();
	// good track
	h_nSigmaE_P__total->Write();
	h_nSigmaPi_P__total->Write();
	h_nSigmaK_P__total-> Write();
	h_nSigmaP_P__total->Write();

	h_pT__TOFMatch->Write();
	h_Eta__TOFMatch->Write();
	h_Phi__TOFMatch->Write();

	//h_Pt_Cen_nSigmaPi__PurePi->Write();
	//h_Eta_Cen_nSigmaPi__PurePi->Write();
	//h_Phi_Cen_nSigmaPi__PurePi->Write();
	//h_Pt_Eta_Phi__PurePi->Write();
	h_Pt_Cen_nSigmaE__PureE->Write();
	h_Eta_Cen_nSigmaE__PureE->Write();
	h_Phi_Cen_nSigmaE__PureE->Write();
	h_Pt_Eta_Phi__PureE->Write();

	h_nSigmaElectron_P__1->Write();
	h_nSigmaElectron_P__2->Write();
	h_nSigmaElectron_P__3->Write();
	h_nSigmaElectron_P__TOFMatch->Write();
	h_nSigmaElectron_P__PID_2->Write();
	h_nSigmaElectron_P__PID_3->Write();
	h_nSigmaElectron_P__EIDcut_1->Write();
	h_nSigmaElectron_P__EIDcut_2->Write();
	h_nSigmaElectron_Eta__EIDcut_2->Write();
	h_nSigmaElectron_P__EIDcut_3->Write();
	h_nSigmaElectron_Eta__EIDcut_3_lowP->Write();
	h_nSigmaElectron_Eta__EIDcut_3_highP->Write();
	h_nSigmaElectron_P__EIDcut_total->Write();
	h_nSigmaPion_P__2->Write();
	h_nSigmaKaon_P__2->Write();
	h_nSigmaPion_P__3->Write();
	h_nSigmaKaon_P__3->Write();
	h_nSigmaProton_P__2->Write();


	h_Mee_Pt_Cen__unlikeSame->Write();
	h_Mee_Pt_Cen__likemm->Write();
	h_Mee_Pt_Cen__likepp->Write();

	mFile->Close();
	return kStOK;
}

Bool_t StPicoDstarMixedMaker::isGoodTrigger(StPicoEvent const* const picoEvent) const
{
  for(auto trg : anaCuts::trigNumber)
  {
    if (picoEvent->isTrigger(trg)) return kTRUE;
  }
  return kFALSE;
}

Bool_t StPicoDstarMixedMaker::isGoodEvent(StPicoEvent const* const picoEvent) const
{
  TVector3 pVtx = picoEvent->primaryVertex();
  return pVtx.z() < anaCuts::Vz_up &&
		 pVtx.z() > anaCuts::Vz_low &&
		 fabs(pVtx.x()) > anaCuts::Verror &&
		 fabs(pVtx.y()) > anaCuts::Verror &&
		 fabs(pVtx.z()) > anaCuts::Verror &&
		 sqrt(pVtx.x()*pVtx.x() + pVtx.y()*pVtx.y()) < anaCuts::Vr &&
		 fabs(pVtx.z() - picoEvent->vzVpd()) < anaCuts::vzVpdVz;
}

Bool_t StPicoDstarMixedMaker::isGoodTrack(StPicoTrack const* trk, StPicoEvent const* const picoEvent) const
{
	TVector3 mom = trk->pMom();
	//trk->gPt() > anaCuts::GPt &&
    //fabs(trk->gMom().Eta()) < anaCuts::Eta &&
	return (
	(mom.Perp() >= 0.2 &&
	fabs(trk->nHitsFit()) >= anaCuts::NHitsFit_highPt &&
	fabs(trk->nHitsDedx()) >= anaCuts::NHitsDedx_highPt &&
	trk->gDCA(picoEvent->primaryVertex()).Mag() <= anaCuts::Dca_highPt)
	||
	(mom.Perp() < 0.2 &&
	fabs(trk->nHitsFit()) >= anaCuts::NHitsFit_lowPt &&
	fabs(trk->nHitsDedx()) >= anaCuts::NHitsDedx_lowPt &&
	trk->gDCA(picoEvent->primaryVertex()).Mag() <= anaCuts::Dca_lowPt)
	) &&
	fabs(trk->nHitsFit()*1.0 / trk->nHitsMax()) >= anaCuts::NHitsFitRatio;
}

Float_t StPicoDstarMixedMaker::getTofBeta(StPicoTrack const* const trk) const
{
  Int_t index2tof = trk->bTofPidTraitsIndex();
  Float_t beta = std::numeric_limits<Float_t>::quiet_NaN();
  if (index2tof >= 0)
  {
    StPicoBTofPidTraits const* const tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);
    if (tofPid)
    {
      beta = tofPid->btofBeta();
      if (beta < 1e-4)
      {
        TVector3 const vtx3 = mPicoDstMaker->picoDst()->event()->primaryVertex();
        StThreeVectorF vtx(vtx3.x(),vtx3.y(),vtx3.z());
        TVector3 const btofHitPos3 = tofPid->btofHitPos();
        StThreeVectorF btofHitPos(btofHitPos3.x(),btofHitPos3.y(),btofHitPos3.z());
        StPicoPhysicalHelix helix = trk->helix(mPicoDstMaker->picoDst()->event()->bField());
        Float_t L = tofPathLength(&vtx, &btofHitPos, helix.curvature());
        Float_t tof = tofPid->btof();
        if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
        else beta = std::numeric_limits<Float_t>::quiet_NaN();
      }
    }
  } 
  return beta;
}

Int_t StPicoDstarMixedMaker::getRefmult6(StPicoDst const* const picoDst, StPicoEvent const* const picoEvent) const
{
  Int_t refMult6 = 0;
  for ( UInt_t iTrk=0; iTrk<picoDst->numberOfTracks(); iTrk++ )// UInt_t是无符号整型
  {
    StPicoTrack* track = picoDst->track( iTrk );
    if ( !track ) continue;
    if ( !track->isPrimary() ) continue; // 只考虑主要径迹
	// 只考虑通过track QA cuts的径迹
    if ( TMath::Abs( track->pMom().Eta() ) > 1.5 ) continue; // if ( TMath::Abs( track->pMom().Eta() ) > anaCuts::Eta ) continue;
    if ( track->pPt() >= 2.0 ) continue; //因为统计量太少？
	if ( track->pPt() <= 0.2) continue; // if (track->pPt() <= anaCuts::GPt) continue;//
    if ( track->gDCA(picoEvent->primaryVertex()).Mag() >= 3. ) continue;
    if ( track->nHitsFit() <= 15 ) continue;
    if ( (Double_t)track->nHitsFit() / track->nHitsPoss() <= 0.52 ) continue;
	refMult6++;
  }
  // cout << "RefMult6: " << refMult6 << endl;
  return refMult6;
}

Double_t StPicoDstarMixedMaker::getTotnMIP(StPicoDst const* const picoDst) const
{
  Double_t totnMIP = 0.;
  Int_t nEpdHits = picoDst->numberOfEpdHits();
  for(Int_t iHit=0; iHit<nEpdHits; iHit++)
  {
    StPicoEpdHit *epdHit = picoDst->epdHit(iHit);
    if( !epdHit ) continue;
	if (epdHit->nMIP() < 0.3) continue;
    else if (epdHit->nMIP() <= 6.) totnMIP += epdHit->nMIP();
	else totnMIP += 6.;
  }
  // cout << "Total nMIP: " << totnMIP << endl;
  return totnMIP;
}