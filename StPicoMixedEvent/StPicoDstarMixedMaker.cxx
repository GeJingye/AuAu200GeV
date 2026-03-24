// ****************************************************
// *                                                  *
// *  Authors: Yuanjing Ji                            *
// *           Guannan Xie <guannanxie@lbl.gov>       *
// *           Mustafa Mustafa <mmustafa@lbl.gov>     *
// *                                                  *
// ****************************************************
 /* ****************************************************************************************
  * read PicoDst document about AuAu 200GeV collision for produciton within iTPC acceptance*
  * ****************************************************************************************
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
#include "../StRefMultCorr/StRefMultCorr.h"
#include "../StRefMultCorr/CentralityMaker.h"
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

const Int_t kMagBins = 2;// 磁场的取向（+/-）
const Int_t kCenBins = 16;// 中心度分成的bin数(0,80,5),80%~100%被忽略了吗？
const Int_t kVzBins = 20;// Vz分成的bin数
Int_t magBufferIndex, cenBufferIndex, vzBufferIndex;// 索引到每个bin
const Int_t kMaxEventsInBuffer = 100;// 每个buffer中事件数最大值
const Int_t kMaxElectrons = 50;// 每个event中正、负电子数的上界

Int_t current_nPositron;// 当前事件中通过筛选的正电子（e+）数
Int_t current_nElectron;// 当前事件中通过筛选的电子（e-）数
TLorentzVector current_positron[kMaxElectrons];// 当前事件中通过筛选的正电子的四动量
TLorentzVector current_electron[kMaxElectrons];// 当前事件中通过筛选的电子的四动量

Int_t nEventsInBuffer[kMagBins][kCenBins][kVzBins];// 每个buffer每个bin中事件数
Bool_t bufferFullFlag[kMagBins][kCenBins][kVzBins];// 每个buffer每个bin是否满员
Int_t buffer_nEPlus[kMagBins][kCenBins][kVzBins][kMaxEventsInBuffer];// 每个buffer每个bin每个事件中正电子数
Int_t buffer_nEMinus[kMagBins][kCenBins][kVzBins][kMaxEventsInBuffer];// 每个buffer每个bin每个事件中电子数
TLorentzVector buffer_ePlus[kMagBins][kCenBins][kVzBins][kMaxEventsInBuffer][kMaxElectrons];// 每个buffer每个bin每个事件每个正电子的四动量
TLorentzVector buffer_eMinus[kMagBins][kCenBins][kVzBins][kMaxEventsInBuffer][kMaxElectrons];// 每个buffer每个bin每个事件每个电子的四动量

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
	readnum.open(mRunNumList);
	Int_t runum = 0;
	Int_t totalNum = 0;
	while (readnum >> runum) {
		mrunnum[runum] = totalNum;
		totalNum++;
	}
	readnum.close();


	// 把一段连续的内存 每个字节 都设为0,例如：Int_t buf[1024]; memset(buf, 0, sizeof(buf));
	memset(nEventsInBuffer, 0, sizeof(nEventsInBuffer));
	memset(bufferFullFlag, 0, sizeof(bufferFullFlag));
	memset(buffer_nEPlus, 0, sizeof(buffer_nEPlus));
	memset(buffer_nEMinus, 0, sizeof(buffer_nEMinus));

	mRefMultCorrUtil = CentralityMaker::instance()->getRefMultCorr();
	mRefMultCorrUtil->setVerbose(kFALSE);

	// Event histograms
	fphiVcut = new TF1("fphiVcut", "0.84326*exp(-49.4819*x)-0.996609*x+0.19801", 0.0, 1.0);
	h_RunNum = new TH1F("h_RunNum", "h_RunNum", totalNum, -0.5, totalNum-0.5);
	h_cen = new TH1F("h_cen", "h_cen", 17, -1.5, 15.5);
	h_cen__reWeight = new TH1F("h_cen__reWeight", "h_cen__reWeight", 17, -1.5, 15.5);
	for (Int_t i = 1; i < h_cen->GetNbinsX() + 1; i++)
	{
		if (i == 1) h_cen->GetXaxis()->SetBinLabel(i, "80-100%");
		else h_cen->GetXaxis()->SetBinLabel(i, Form("%d-%d%%", 85 - 5 * i, 90 - 5 * i));
	}
	h_RefMult = new TH1F("h_RefMult", "h_RefMult", 250, 0, 250);// 参考多重数
	h_nTofMat_RefMul = new TH2F("h_nTofMat_RefMul", "RefMul VS nTofmatch;RefMul;nTofMatch", 250, 0, 250, 250, 0, 250);// 与TOF匹配的径迹数vs参考多重数关系
	h_mRefMult = new TH1F("h_mRefMult", "h_mRefMult", 250, 0, 250);// 监控事件多重数分布->计算中心度
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
	h_nHitsFit_Pt_Eta = new TH3F("h_nHitsFit_Pt_Eta", "nHitsFit vs p_{T} vs #eta;p_{T} (GeV/c);#eta;nHitsFit", 100, 0., 10., 40, -2., 2., 80, 0., 80.);
	h_nHitsDEdx_Pt_Eta = new TH3F("h_nHitsDEdx_Pt_Eta", "nHitsDedx vs p_{T} vs #eta;p_{T} (GeV/c);#eta;nHitsDedx", 100, 0., 10., 40, -2., 2., 80, 0., 80.);

	h_pDca = new TH1F("h_pDca", "pDca;DCA;counts", 50, 0., 5.);// p代表primary，pDCA指主径迹与重建顶点的最小距离
	h_ppT = new TH1F("h_ppT", "primary p_{T};p_{T} (GeV/c);counts", 1000, 0., 10.);
	h_pP = new TH1F("h_pP", "primary p;p (GeV/c);counts", 1000, 0., 10.);
	h_pP_ppT = new TH2F("h_pP_ppT", "primary p vs primary p_{T};p (GeV/c);p_{T} (GeV/c)", 1000, 0., 10., 1000, 0., 10.);
	h_gPt = new TH1F("h_gPt", "global p_{T} of track;global p_{T} (GeV/c);counts", 1000, 0., 10.);
	h_pEta = new TH1F("h_pEta", "primary #eta;#eta;counts", 500, -2.5, 2.5);
	h_pPhi = new TH1F("h_pPhi", "primary #phi;#phi;counts", 80, -4.0, 4.0);
	h_ppTc_pEta = new TH2F("h_ppTc_pEta", "p_{T}*q vs #eta;p_{T}*q (GeV/c);#eta", 400, -10, 10, 500, -2.5, 2.5);
	h_ppTc_pPhi = new TH2F("h_ppTc_pPhi", "p_{T}*q vs #phi;p_{T}*q (GeV/c);#phi", 400, -10, 10, 800, -4, 4);
	h_pDca_Eta_NHitsFit = new TH3F("h_pDca_Eta_NHitsFit", "NHitsFit vs #eta vs DCA;DCA;#eta;NHitsFit", 50, 0., 5., 40, -2., 2., 90, 0, 90);
	h_pDca_Pt_Eta = new TH3F("h_pDca_Pt_Eta", "#eta vs p_{T} vs DCA;DCA;p_{T};#eta", 50, 0., 5.,1000 ,0., 10., 40, -2., 2.);
	// GoodTrack径迹信息
	h_nSigmaElectron_P = new TH2F("h_nSigmaElectron_P", "n#sigma_{e} vs p;p (GeV/c);n#sigma_{e}", 500, 0, 5, 4000, -20, 20);
	h_nSigmaEcorr_P = new TH2F("h_nSigmaEcorr_P", "corrected n#sigma_{e} vs p;p (GeV/c);corrected n#sigma_{e}", 500, 0, 5, 4000, -20, 20);
	h_nSigmaPion_P = new TH2F("h_nSigmaPion_P", "n#sigma_{pi} vs p;p (GeV/c);n#sigma_{pi}", 500, 0, 5, 4000, -20, 20);
	h_nSigmaPicorr_P = new TH2F("h_nSigmaPicorr_P", "corrected n#sigma_{pi} vs p;p (GeV/c);corrected n#sigma_{pi}", 500, 0, 5, 4000, -20, 20);
	h_nSigmaKaon_P = new TH2F("h_nSigmaKaon_P", "n#sigma_{K} vs p;p (GeV/c);n#sigma_{K}", 500, 0, 5, 4000, -20, 20);
	h_nSigmaProton_P = new TH2F("h_nSigmaProton_P", "n#sigma_{P} vs p;p (GeV/c);n#sigma_{P}", 500, 0, 5, 4000, -20, 20);
	h_dEdx_Pc = new TH2F("h_dEdx_Pc", "dE/dx vs p*q;p*q(GeV/c);#frac{dE}{dx} (GeV cm^{2}/g)", 400, -2, 2, 400, 0, 25);
	h_m2    = new TH1F("h_m2", "m^{2};m^{2};counts", 2000, -0.5, 1.5);
	h_m2_Pc = new TH2F("h_m2_Pc", "m^{2} vs p*q;p*q (GeV/c);m^{2} (GeV/c^{2})^{2}", 1000, -5, 5, 2000, -0.5, 1.5);
	// nSigmaE
	h_Pt_Cen_nSigmaE__PureE = new TH3F("h_Pt_Cen_nSigmaE__PureE", "p_{T} vs Cen vs n#sigma_{e};p_{T};Cen;n#sigma_{e}", 50, 0, 5, 8, 0, 8, 200, -10, 10);
	h_Eta_Cen_nSigmaE__PureE = new TH3F("h_Eta_Cen_nSigmaE__PureE", "#eta vs Cen vs n#sigma_{e};#eta;Cen;n#sigma_{e}", 40, -2, 2, 8, 0, 8, 200, -10, 10);
	h_Phi_Cen_nSigmaE__PureE = new TH3F("h_Phi_Cen_nSigmaE__PureE", "#phi vs Cen vs n#sigma_{e};#phi;Cen;n#sigma_{e}", 64, -3.2, 3.2, 8, 0, 8, 200, -10, 10);

	h_Pt_Cen_nSigmaE  = new TH3F("h_Pt_Cen_nSigmaE",  "n#sigma_{e} vs p_{T} vs Cen;p_{T} (GeV/c);Cen;n#sigma_{e}",500, 0, 5, 16 ,0, 16, 2000, -10, 10);
	h_Eta_Cen_nSigmaE = new TH3F("h_Eta_Cen_nSigmaE", "n#sigma_{e} vs #eta vs Cen;#eta;Cen;n#sigma_{e}",         400, -2.0, 2.0, 16, 0, 16, 2000, -10, 10);
	h_Phi_Cen_nSigmaE = new TH3F("h_Phi_Cen_nSigmaE", "n#sigma_{e} vs #phi vs Cen;#phi;Cen;n#sigma_{e}",         640, -3.2, 3.2, 16, 0, 16, 2000, -10, 10);
	h_Pt_Cen_nSigmaEcorr  = new TH3F("h_Pt_Cen_nSigmaEcorr",  "corrected n#sigma_{e} vs p_{T} vs Cen;p_{T} (GeV/c);corrected n#sigma_{e}", 500, 0, 5, 16, 0, 16, 2000, -10, 10);
	h_Eta_Cen_nSigmaEcorr = new TH3F("h_Eta_Cen_nSigmaEcorr", "corrected n#sigma_{e} vs #eta vs Cen;#eta;Cen;corrected n#sigma_{e}",        400, -2.0, 2.0, 16, 0, 16, 2000, -10, 10);
	h_Phi_Cen_nSigmaEcorr = new TH3F("h_Phi_Cen_nSigmaEcorr", "corrected n#sigma_{e} vs #phi vs Cen;#phi;Cen;corrected n#sigma_{e}",        640, -3.2, 3.2, 16, 0, 16, 2000, -10, 10);
	// nSigmaPi
	h_Pt_Cen_nSigmaPion = new TH3F("h_Pt_Cen_nSigmaPion", "n#sigma_{pi} vs p_{T} vs Cen;p_{T} (GeV/c);Cen;n#sigma_{pi}", 500, 0, 5, 16, 0, 16, 2000, -10, 10);
	h_Eta_Cen_nSigmaPion = new TH3F("h_Eta_Cen_nSigmaPion", "n#sigma_{pi} vs #eta vs Cen;#eta;Cen;n#sigma_{pi}", 400, -2.0, 2.0, 16, 0, 16, 2000, -10, 10);
	h_Phi_Cen_nSigmaPion = new TH3F("h_Phi_Cen_nSigmaPion", "n#sigma_{pi} vs #phi vs Cen;#phi;Cen;n#sigma_{pi}", 640, -3.2, 3.2, 16, 0, 16, 2000, -10, 10);
	h_Pt_Cen_nSigmaPicorr = new TH3F("h_Pt_Cen_nSigmaPicorr", "corrected n#sigma_{pi} vs p_{T} vs Cen;p_{T} (GeV/c);Cen;corrected n#sigma_{pi}", 500, 0, 5, 16, 0, 16, 2000, -10, 10);
	h_Eta_Cen_nSigmaPicorr = new TH3F("h_Eta_Cen_nSigmaPicorr", "corrected n#sigma_{pi} vs #eta vs Cen;#eta;Cen;corrected n#sigma_{pi}", 400, -2.0, 2.0, 16, 0, 16, 2000, -10, 10);
	h_Phi_Cen_nSigmaPicorr = new TH3F("h_Phi_Cen_nSigmaPicorr", "corrected n#sigma_{pi} vs #phi vs Cen;#phi;Cen;corrected n#sigma_{pi}", 640, -3.2, 3.2, 16, 0, 16, 2000, -10, 10);

	// group 1 2 3
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
	h_nSigmaKaon_P__2   = new TH2F("h_nSigmaKaon_P__2",  "n#sigma_{K} vs p (p_{T}>0.2,0.9<|#eta|<1.8);p (GeV/c);n#sigma_{K}", 500, 0, 5, 2000, -10, 10);
	h_nSigmaProton_P__2 = new TH2F("h_nSigmaProton_P__2","n#sigma_{P} vs p (p_{T}>0.2,0.9<|#eta|<1.8);p (GeV/c);n#sigma_{P}", 500, 0, 5, 2000, -10, 10);
	h_nSigmaPion_P__3   = new TH2F("h_nSigmaPion_P__3",  "n#sigma_{pi} vs p (p_{T}<0.2,|#eta|<1.8);p (GeV/c);n#sigma_{pi}", 500, 0, 5, 2000, -10, 10);
	h_nSigmaKaon_P__3   = new TH2F("h_nSigmaKaon_P__3",  "n#sigma_{K} vs p (p_{T}<0.2,|#eta|<1.8);p (GeV/c);n#sigma_{K}", 500, 0, 5, 2000, -10, 10);
	h_pT_Eta_Phi__EIDcut_3_lowP_e = new TH3F("h_pT_Eta_Phi__EIDcut_3_lowP_e", "p_{T} vs #eta vs #phi;p_{T} (GeV/c);#eta;#phi", 50, 0., 0.5, 50, -2.5, 2.5, 64,-3.2,3.2);
	h_pT_Eta_Phi__EIDcut_3_lowP_p = new TH3F("h_pT_Eta_Phi__EIDcut_3_lowP_p", "p_{T} vs #eta vs #phi;p_{T} (GeV/c);#eta;#phi", 50, 0., 0.5, 50, -2.5, 2.5, 64,-3.2,3.2);

	// TOF匹配的径迹信息
	h_invBeta_P__TOFMatch = new TH2F("h_invBeta_P__TOFMatch", "1/#beta vs p;p (GeV/c);1/#beta", 500, 0, 5, 5000, 0, 5);
	h_pT__TOFMatch  = new TH1F("h_pT__TOFMatch", "p_{T} of TOF matched tracks;p_{T} (GeV/c);counts", 1000, 0., 10.);
	h_Eta__TOFMatch = new TH1F("h_Eta__TOFMatch", "#eta of TOF matched tracks;#eta;counts", 500, -2.5, 2.5);
	h_Phi__TOFMatch = new TH1F("h_Phi__TOFMatch", "#phi of TOF matched tracks;#phi;counts", 64, -3.2, 3.2);
	h_nSigmaEcorr_P__TOFMatch = new TH2F("h_nSigmaEcorr_P__TOFMatch", "n#sigma_{e} vs p (TOF acceptance);p (GeV/c);n#sigma_{e}", 500, 0, 5, 2000, -10, 10);
	// 经TPC和TOF判选后的电子信息
	h_eNumber_Cen = new TH2F("h_eNumber_Cen", ";Num.;Cen", 50, 0., 50., 16, 0., 16.);
	h_pT__electrons = new TH1F("h_pT__electrons", "p_{T} of electrons;p_{T} (GeV/c)", 100, 0., 5.);
	h_eta__electrons = new TH1F("h_eta__electrons", "#eta of electrons;#eta", 50, -2.5, 2.5);
	h_phi__electrons = new TH1F("h_phi__electrons", "#phi of electrons;#phi", 64,-3.2,3.2);
	h_pT__electrons_w_PhiV_Cut = new TH1F("h_pT__electrons_w_PhiV_Cut", "p_{T} of electrons with #phi_{V} cut;p_{T} (GeV/c)", 100, 0., 5.);//"w"表示with, "wo"表示without
	h_eta__electrons_w_PhiV_Cut = new TH1F("h_eta__electrons_w_PhiV_Cut", "#eta of electrons with #phi_{V} cut;#eta", 50, -2.5, 2.5);
	h_phi__electrons_w_PhiV_Cut = new TH1F("h_phi__electrons_w_PhiV_Cut", "#phi of electrons with #phi_{V} cut;#phi", 64,-3.2,3.2);
	h_pT__positrons = new TH1F("h_pT__positrons", "p_{T} of positrons;p_{T} (GeV/c)", 100, 0, 5.);
	h_eta__positrons = new TH1F("h_eta__positrons", "#eta of positrons;#eta", 50, -2.5, 2.5);
	h_phi__positrons = new TH1F("h_phi__positrons", "#phi of positrons;#phi", 64,-3.2,3.2);
	h_pT__positrons_w_PhiV_Cut = new TH1F("h_pT__positrons_w_PhiV_Cut", "p_{T} of positrons with #phi_{V} cut;p_{T} (GeV/c)", 100, 0., 5.);//"w"表示with, "wo"表示without
	h_eta__positrons_w_PhiV_Cut = new TH1F("h_eta__positrons_w_PhiV_Cut", "#eta of positrons with #phi_{V} cut;#eta", 50, -2.5, 2.5);
	h_phi__positrons_w_PhiV_Cut = new TH1F("h_phi__positrons_w_PhiV_Cut", "#phi of positrons with #phi_{V} cut;#phi", 64,-3.2,3.2);

	h_Rapidity__unlikeSame = new TH1F("h_Rapidity__unlikeSame", "y distribution of e^{+}e^{-};y;counts", 800, -4, 4);
	h_Mee_PhiV__unlikeSame = new TH2F("h_Mee_PhiV__unlikeSame", "Mee vs #phi_{V};Mee(GeV/c^{2});#phi_{V}", 800, 0, 4, 100, 0, 1);
	h_Mee__unlikeSame = new TH1F("h_Mee__unlikeSame", "Mee without #phi_{V} cut;Mee(GeV/c^{2})", 800, 0, 4);
	h_Mee__unlikeSame__w_PhiV_Cut = new TH1F("h_Mee__unlikeSame__w_PhiV_Cut", "Mee with #phi_{V} cut;Mee(GeV/c^{2})", 800, 0, 4);

	h_Mee_Pt_Cen__unlikeSame = new TH3F("h_Mee_Pt_Cen__unlikeSame", "Mee vs p_{T} vs Cen;Mee (GeV/c^{2});p_{T} (GeV/c);Cen", 800, 0, 4, 100, 0, 5, 16, 0, 16);
	h_Mee_Pt_Cen__likemm = new TH3F("h_Mee_Pt_Cen__likemm", "Mee vs p_{T} vs Cen;Mee (GeV/c^{2});p_{T} (GeV/c);Cen", 800, 0, 4, 100, 0, 5, 16, 0, 16);
	h_Mee_Pt_Cen__likepp = new TH3F("h_Mee_Pt_Cen__likepp", "Mee vs p_{T} vs Cen;Mee (GeV/c^{2});p_{T} (GeV/c);Cen", 800, 0, 4, 100, 0, 5, 16, 0, 16);
	h_Mee_Pt_Cen__unlikeMixed = new TH3F("h_Mee_Pt_Cen__unlikeMixed", "Mee vs p_{T} vs Cen;Mee(GeV/c^{2});p_{T} (GeV/c);Cen", 800, 0, 4, 100, 0, 5, 16, 0, 16);
	h_Mee_Pt_Cen__likemmMixed = new TH3F("h_Mee_Pt_Cen__likemmMixed", "Mee vs p_{T} vs Cen;Mee(GeV/c^{2});p_{T} (GeV/c);Cen", 800, 0, 4, 100, 0, 5, 16, 0, 16);
	h_Mee_Pt_Cen__likeppMixed = new TH3F("h_Mee_Pt_Cen__likeppMixed", "Mee vs p_{T} vs Cen;Mee(GeV/c^{2});p_{T} (GeV/c);Cen", 800, 0, 4, 100, 0, 5, 16, 0, 16);

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
  h_RunNum->Fill(mrunnum[mRunId]);// 填充不同run号的事例数；

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

	if (fabs(mVpdVz + 999.0)>1e-2 && fabs(mVpdVz) < 1000.0)// STAR约定：当 VPD 无法重建顶点时，会把 mVpdVz 设为 -999 cm 作为占位标志。
	{
		h_VpdVz->Fill(mVpdVz);
		h_VpdVz_Vz->Fill(mVz, mVpdVz);
		h_VpdVzmVz->Fill(mVpdVz - mVz);
	}

	// 获取多重数
	Refmult = picoEvent->refMult();//给出的是在线/原始 TPC 多重数（reference multiplicity），未经任何修正
	mRefmult = getRefmult(picoDst, picoEvent);// 用户离线重算的多重数，满足严格几何与质量要求

	h_RefMult->Fill(Refmult);
	h_nTofMat_RefMul->Fill(Refmult, picoEvent->nBTOFMatch());
    h_mRefMult->Fill(mRefmult);

	// 通过离线多重数计算该事例的中心度
	mRefMultCorrUtil->init(mRunId);
	mRefMultCorrUtil->initEvent(mRefmult, mVz);
	//Double_t reWeight = 1;
	Double_t reWeight = mRefMultCorrUtil->getWeight();
	//mCen9 = 1;
	mCen9 = mRefMultCorrUtil->getCentralityBin9();
	//mCen16 = 1;
	mCen16 = mRefMultCorrUtil->getCentralityBin16();
	// 不同条件cut后的事例数统计
    Bool_t vzcut = mVz < anaCuts::Vz_up && mVz > anaCuts::Vz_low;
    Bool_t verrcut = !(fabs(mVx) < anaCuts::Verror && fabs(mVy) < anaCuts::Verror && fabs(mVz) < anaCuts::Verror);// Vx,Vy,Vz<1.0e-5 cm, why? too small that better than resolution. 
    Bool_t vrcut =  mVr < anaCuts::Vr;
	Bool_t notPileUp = picoEvent->refMult()<picoEvent->btofTrayMultiplicity()*0.36+45;// mRefMultCorrUtil->isPileUpEvent(mRefmult, picoEvent->nBTOFMatch(), mVz, 1);
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
      h_cen->Fill(mCen16);// 填充中心度
	  h_cen__reWeight->Fill(mCen16, reWeight);
      // ******************以下分析均基于0~80%中心度***********************
      if(!cen0280cut) return kStOK;
	  // 通过mag, cen, vz的值确定magBufferIndex, cenBufferIndex, vzBufferIndex的索引
      if(kMagBins==2)
      {
        if(mBfield<0.) magBufferIndex = 0;
        if(mBfield>0.) magBufferIndex = 1;
      }
      else magBufferIndex = 0;
	  if (magBufferIndex != 0 && magBufferIndex != 1) return kStOK;// 界外判断
      cenBufferIndex = mCen16;
	  if (cenBufferIndex < 0 || cenBufferIndex >= kCenBins) return kStOK;// 界外判断
      vzBufferIndex = static_cast<Int_t>((mVz-anaCuts::Vz_low)/(anaCuts::Vz_up - anaCuts::Vz_low))*kVzBins;
	  if (vzBufferIndex < 0 || vzBufferIndex >= kVzBins) return kStOK;// 界外判断

	  current_nPositron = 0;
      current_nElectron = 0;

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
		  h_pP_ppT->Fill(mom.Mag(), mom.Perp());
		  h_ppT->Fill(mom.Perp());
		  h_pP->Fill(mom.Mag());
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

		  h_nSigmaElectron_P->Fill(mom.Mag(), nSigmaE);
		  h_nSigmaPion_P->Fill(mom.Mag(), nSigmaPi);
		  h_nSigmaKaon_P->Fill(mom.Mag(), nSigmaK);
		  h_nSigmaProton_P->Fill(mom.Mag(), nSigmaP);
		  // ******************以下分析均基于满足goodTrackCuts的主径迹**************************
		  /*    Track with conditions in StAnaCuts.h, such as:
				NHitsFit > 40;
				NHitsFitRatio > 0.52;
				NHitsDedx > 30;
				Dca < 1 cm;
				GPt > 0.2 GeV/c;
			   |Eta| < 1.8;		  */
		  Bool_t goodtrack = isGoodTrack(trk, picoEvent);
		  if (!goodtrack) continue;
		  h_passTrkcut->Fill(2);
		  h_dEdx_Pc->Fill(mom.Mag()*trk->charge(), trk->dEdx());// mom.Mag()指动量模
		  h_m2->Fill(m2);
		  h_m2_Pc->Fill(mom.Mag()*trk->charge(), m2);

		  // Recalibrate nSigmaElectron nSigmaPion
		  Double_t nSigmaE_corrfactor = getNSigmaECorr(mom);
		  Double_t nSigmaE_corr = nSigmaE - nSigmaE_corrfactor;
		  h_nSigmaEcorr_P->Fill(mom.Mag(), nSigmaE_corr);

		  Double_t nSigmaPi_corrfactor = getNSigmaPiKPCorr(trk->charge() > 0 ? 1 : -1, mom);
		  Double_t nSigmaPi_corr = nSigmaPi - nSigmaPi_corrfactor;
		  h_nSigmaPicorr_P->Fill(mom.Mag(), nSigmaPi_corr);
		  //Double_t temp = nSigmaE_corr; nSigmaE_corr = nSigmaE; nSigmaE = temp;//将nsigmaE和nsigmaEcorr置换

	/*	  if (nSigmaPi > -3.5 && nSigmaPi < 3.5)
		  {
			  h_Pt_Cen_nSigmaPion->Fill(mom.Perp(), mCen16, nSigmaPi, reWeight);
			  h_Eta_Cen_nSigmaPion->Fill(mom.Eta(), mCen16, nSigmaPi, reWeight);
			  h_Phi_Cen_nSigmaPion->Fill(mom.Phi(), mCen16, nSigmaPi, reWeight);
			  h_Pt_Cen_nSigmaPicorr->Fill(mom.Perp(), mCen16, nSigmaPi_corr, reWeight);
			  h_Eta_Cen_nSigmaPicorr->Fill(mom.Eta(), mCen16, nSigmaPi_corr, reWeight);
			  h_Phi_Cen_nSigmaPicorr->Fill(mom.Phi(), mCen16, nSigmaPi_corr, reWeight);
		  }*/


		  // Electrons IDentification
		  Bool_t isTOFElectron__1 = kFALSE;
		  Bool_t isTPCElectron__1 = kFALSE;

		  Bool_t isTPCKaon__2 = kTRUE;
		  Bool_t isTPCProton__2 = kTRUE;
		  Bool_t isTPCElectron__2 = kFALSE;

		  Bool_t isTPCPion__3 = kTRUE;
		  Bool_t isTPCKaon__3 = kTRUE;
		  Bool_t isTPCElectron__3 = kFALSE;
		  Bool_t isLowPElectron__3 = kFALSE;
		  Bool_t isSomePhiElectron__3 = kFALSE;
		  Bool_t isLowEtaElectron__3 = kFALSE;
		  Bool_t isValidElectron__lowP_3 = kTRUE;

		  Bool_t isElectronRegion1 = kFALSE;
		  Bool_t isElectronRegion2 = kFALSE;
		  Bool_t isElectronRegion3 = kFALSE;

		  // TPC cut + TOF cut
		  // groupn 1
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
				  if (fabs(1.0 / beta - 1) < 0.02)
				  {
					  isTOFElectron__1 = kTRUE;
					  h_nSigmaEcorr_P__TOFMatch->Fill(mom.Mag(), nSigmaE_corr);
				  }
			  }
		  }
		  // group 1
		  if (mom.Mag() <= 1.0) isTPCElectron__1 = nSigmaE<3.0 && nSigmaE>(2.5*mom.Mag() - 3.5);
		  if (mom.Mag() > 1.0)  isTPCElectron__1 = nSigmaE<3.0 && nSigmaE>-1.0;

		  // group 2
		  isTPCKaon__2     = (mom.Mag() <= 1.0) ? (nSigmaK > -2.0*mom.Mag() && nSigmaK < 3.0)				   : (nSigmaK > -2.0 && nSigmaK < 5.0);
		  isTPCProton__2   = (mom.Mag() <= 1.5) ? (nSigmaP > -1.33*mom.Mag() && nSigmaP < -1.33*mom.Mag()+6.0) : (nSigmaP > -2.0 && nSigmaP < 4.0);
		  isTPCElectron__2 = (mom.Mag() <= 1.5) ? (nSigmaE > 1.5*mom.Mag()-2.25 && nSigmaE < 3.0)			   : (nSigmaE > 0.0 && nSigmaE < 3.0);
		  // group 3
		  if (mom.Mag() > 0.0 && mom.Mag() <= 0.2) isTPCPion__3 = nSigmaPi > -20.0 * mom.Mag() - 1.0 && nSigmaPi < -5.71 * mom.Mag() + 5.5;
		  if (mom.Mag() > 0.2 && mom.Mag() < 0.7)  isTPCPion__3 = nSigmaPi > -5.0 && nSigmaPi < -5.71 * mom.Mag() + 5.5;
		  if (mom.Mag() > 0.0 && mom.Mag() < 0.7)  isTPCKaon__3 = nSigmaK > 13.57 * mom.Mag() - 6.0 && nSigmaK < 5.0;
		  if (mom.Mag() > 0.0 && mom.Mag() < 0.14)  isTPCElectron__3 = nSigmaE > -4.0	&& nSigmaE < 1.0;
		  if (mom.Mag() > 0.14 && mom.Mag() < 0.18) isTPCElectron__3 = nSigmaE > -2.0	&& nSigmaE < 1.0;
		  if (mom.Mag() > 0.18 && mom.Mag() < 0.7)  isTPCElectron__3 = nSigmaE > -2.0	&& nSigmaE < 3.0;

		  // Fill Histogram
		  // group 1
		  if (mom.Perp() > 0.2 && fabs(mom.Eta()) < anaCuts::Eta)
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
		  // group 2
		  if (mom.Perp() > 0.2 && fabs(mom.Eta()) > anaCuts::Eta && fabs(mom.Eta()) < 1.8)
		  { 
			  h_nSigmaElectron_P__2->Fill(mom.Mag(), nSigmaE);
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
		  // group 3
		  if (mom.Perp() < 0.2 && fabs(mom.Eta()) < 1.8)
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
					  isSomePhiElectron__3 = mom.Phi() < -1.1 && mom.Phi() > -1.6;
					  h_nSigmaElectron_P__EIDcut_3->Fill(mom.Mag(), nSigmaE);
					  isElectronRegion3 = kTRUE;
					  if (isLowPElectron__3)
					  {
						  h_nSigmaElectron_Eta__EIDcut_3_lowP->Fill(mom.Eta(), nSigmaE);
						  isValidElectron__lowP_3 = (mom.Phi() < 3.2) && (fabs(mom.Eta()) < 0.8 && fabs(mom.Eta()) > 0.2) && (mom.Perp() > 0.07 && mom.Perp() < 0.14);
						  if (trk->charge() > 0)
						  {
							  h_pT_Eta_Phi__EIDcut_3_lowP_p->Fill(mom.Perp(), mom.Eta(), mom.Phi());
						  }
						  if (trk->charge() < 0)
						  {
							  h_pT_Eta_Phi__EIDcut_3_lowP_e->Fill(mom.Perp(), mom.Eta(), mom.Phi());
						  }

					  }
					  if (!isLowPElectron__3)h_nSigmaElectron_Eta__EIDcut_3_highP->Fill(mom.Eta(), nSigmaE);
				  }
			  }
		  }

		  //if (isElectronRegion1)//
		  if (isElectronRegion1 || isElectronRegion2 || (isElectronRegion3 && isValidElectron__lowP_3))// && !isLowPElectron__3;isLowEtaElectron__3
		  {
			h_Pt_Cen_nSigmaE->Fill(mom.Perp(), mCen16, nSigmaE, reWeight);
			h_Eta_Cen_nSigmaE->Fill(mom.Eta(), mCen16, nSigmaE, reWeight);
			h_Phi_Cen_nSigmaE->Fill(mom.Phi(), mCen16, nSigmaE, reWeight);
			//h_Pt_Cen_nSigmaEcorr->Fill(mom.Perp(), mCen16, nSigmaE_corr, reWeight);
			//h_Eta_Cen_nSigmaEcorr->Fill(mom.Eta(), mCen16, nSigmaE_corr, reWeight);
			//h_Phi_Cen_nSigmaEcorr->Fill(mom.Phi(), mCen16, nSigmaE_corr, reWeight);
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

				//current_electron[current_nElectron].SetPx(mom.X());
				//current_electron[current_nElectron].SetPy(mom.Y());
				//current_electron[current_nElectron].SetPz(mom.Z());
				//current_electron[current_nElectron].SetE(sqrt(pow(M_electron, 2.0) + pow(mom.Mag(), 2.0)));
				//current_nElectron++;
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

				//current_positron[current_nPositron].SetPx(mom.X());
				//current_positron[current_nPositron].SetPy(mom.Y());
				//current_positron[current_nPositron].SetPz(mom.Z());
				//current_positron[current_nPositron].SetE(sqrt(pow(M_electron, 2.0) + pow(mom.Mag(), 2.0)));
				//current_nPositron++;
			}
          }// 填充单径迹的正/负电子信息
      }// 填充单事例所有径迹的正负电子信息

      Int_t x=0;
      Int_t y=0;
      Int_t num_electron = electroninfo.size();
      Int_t num_positron = positroninfo.size();
	  h_eNumber_Cen->Fill(num_electron + num_positron, mCen16, reWeight);
      TLorentzVector eepair(0,0,0,0);
      TLorentzVector particle1_4V(0,0,0,0);
      TLorentzVector particle2_4V(0,0,0,0);
      
	  // \phi_v cut && pure electron sample selection
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
			  //pure e sample
			  if (eepair.M() < 0.015)
			  {
				  positroninfo[x].isPureE = kTRUE;
				  electroninfo[y].isPureE = kTRUE;
			  }
			  //phiV cut
			  Double_t angleV = getPhiVAngle(particle1_4V, particle2_4V, 1, -1);// 注意参数1、-1的选取要求
			  Double_t angleVcut = fphiVcut->Eval(eepair.M()); // 根据fphiVcut关于pair-M的函数取值
			  if (fabs(eepair.Rapidity()) <= 1) h_Mee_PhiV__unlikeSame->Fill(eepair.M(), angleV);
			  if (eepair.M() < anaCuts::PhiVCutMRange && angleV < angleVcut)
			  {
				  positroninfo[x].isPhotonicE = kTRUE;
				  electroninfo[y].isPhotonicE = kTRUE;
			  }
		  }
	  }
	  for (x = 0; x < num_positron; x++)//填充current_positron，只填充isPhotonicE=False的正电子
	  {
		  if (positroninfo[x].isPureE)
		  {
			  //h_Pt_Cen_nSigmaE__PureE->Fill(positroninfo[x].pt, mCen16, positroninfo[x].nSigmaE, reWeight);
			  //h_Eta_Cen_nSigmaE__PureE->Fill(positroninfo[x].eta, mCen16, positroninfo[x].nSigmaE, reWeight);
			  //h_Phi_Cen_nSigmaE__PureE->Fill(positroninfo[x].phi, mCen16, positroninfo[x].nSigmaE, reWeight);
		  }
		  if (!positroninfo[x].isPhotonicE)
		  {
			  current_positron[current_nPositron].SetPx(positroninfo[x].p1);
			  current_positron[current_nPositron].SetPy(positroninfo[x].p2);
			  current_positron[current_nPositron].SetPz(positroninfo[x].p3);
			  current_positron[current_nPositron].SetE(sqrt(pow(M_electron, 2.0) + pow(positroninfo[x].p, 2.0)));
			  current_nPositron++;
		  }
	  }
	  for (x = 0; x < num_electron; x++)//填充current_positron，只填充isPhotonicE=False的电子
	  {
		  if (electroninfo[x].isPureE)
		  {
			  //h_Pt_Cen_nSigmaE__PureE->Fill(positroninfo[x].pt, mCen16, positroninfo[x].nSigmaE, reWeight);
			  //h_Eta_Cen_nSigmaE__PureE->Fill(positroninfo[x].eta, mCen16, positroninfo[x].nSigmaE, reWeight);
			  //h_Phi_Cen_nSigmaE__PureE->Fill(positroninfo[x].phi, mCen16, positroninfo[x].nSigmaE, reWeight);
		  }
		  if (!electroninfo[x].isPhotonicE)
		  {
			  current_electron[current_nElectron].SetPx(electroninfo[x].p1);
			  current_electron[current_nElectron].SetPy(electroninfo[x].p2);
			  current_electron[current_nElectron].SetPz(electroninfo[x].p3);
			  current_electron[current_nElectron].SetE(sqrt(pow(M_electron, 2.0) + pow(electroninfo[x].p, 2.0)));
			  current_nElectron++;
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
			  particle2_4V.SetPx(electroninfo[y].p1);
			  particle2_4V.SetPy(electroninfo[y].p2);
			  particle2_4V.SetPz(electroninfo[y].p3);
			  particle2_4V.SetE(electroninfo[y].energy);
			  eepair = particle1_4V + particle2_4V;
			  //Double_t angleV = getPhiVAngle(particle1_4V, particle2_4V, 1, -1);// 注意参数1、-1的选取要求
			  //Double_t angleVcut = fphiVcut->Eval(eepair.M()); // 根据fphiVcut关于pair-M的函数取值
			  //if (eepair.M() > anaCuts::PhiVCutMRange || angleV > angleVcut)
			  if (fabs(eepair.Rapidity()) <= 1) h_Mee__unlikeSame->Fill(eepair.M());
			  if (!positroninfo[x].isPhotonicE && !electroninfo[y].isPhotonicE)
			  {
				  h_Rapidity__unlikeSame->Fill(eepair.Rapidity());
				  if (fabs(eepair.Rapidity()) <= 1)// 为什么需要在中心快度区？
				  {
					  h_Mee__unlikeSame__w_PhiV_Cut->Fill(eepair.M());
					  h_Mee_Pt_Cen__unlikeSame->Fill(eepair.M(), eepair.Perp(), mCen16, reWeight);
				  }
			  }
		  }
	  }// end: +-

	  // --
      for(x=0;x<num_electron;x++)
      {
		particle1_4V.SetPx(electroninfo[x].p1);
		particle1_4V.SetPy(electroninfo[x].p2);
		particle1_4V.SetPz(electroninfo[x].p3);
		particle1_4V.SetE(electroninfo[x].energy);
		for (y = x + 1; y < num_electron; y++)//从x+1开始，避免自组合和重复组合
		{
			particle2_4V.SetPx(electroninfo[y].p1);
			particle2_4V.SetPy(electroninfo[y].p2);
			particle2_4V.SetPz(electroninfo[y].p3);
			particle2_4V.SetE(electroninfo[y].energy);
			eepair = particle1_4V + particle2_4V;
			//Double_t angleV = getPhiVAngle(particle1_4V, particle2_4V, 1, -1);// 注意参数1、-1的选取要求
			//Double_t angleVcut = fphiVcut->Eval(eepair.M()); // 根据fphiVcut关于pair-M的函数取值
			//if (eepair.M() > anaCuts::PhiVCutMRange || angleV > angleVcut)
			if (!electroninfo[x].isPhotonicE && !electroninfo[y].isPhotonicE)
			{
				if (fabs(eepair.Rapidity()) <= 1)// 判断重建的粒子是否在中心快度区，为什么需要在中心快度区？
				{
					h_Mee_Pt_Cen__likemm->Fill(eepair.M(), eepair.Perp(), mCen16, reWeight);
				}
			}
		}
      }// end: for(x=0;x<num_electron;x++)
	  // ++
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
			//Double_t angleV = getPhiVAngle(particle1_4V, particle2_4V, 1, -1);// 注意参数1、-1的选取要求
			//Double_t angleVcut = fphiVcut->Eval(eepair.M()); // 根据fphiVcut关于pair-M的函数取值
			//if (eepair.M() > anaCuts::PhiVCutMRange || angleV > angleVcut)
			if (!positroninfo[x].isPhotonicE && !positroninfo[y].isPhotonicE)
			{
				if (fabs(eepair.Rapidity()) <= 1)
				{
					h_Mee_Pt_Cen__likepp->Fill(eepair.M(), eepair.Perp(), mCen16, reWeight);
				}
			}
        }
      }// end: for(x=0;x<num_positron;x++)

	  // 读取经过\phi_v cut前、后的正、负电子的横动量，赝快度，方位角
	  for (x = 0; x < num_positron; x++)
	  {
		  h_pT__positrons->Fill(positroninfo[x].pt);
		  h_eta__positrons->Fill(positroninfo[x].eta);
		  h_phi__positrons->Fill(positroninfo[x].phi);
		  if (!positroninfo[x].isPhotonicE)
		  {
			  h_pT__positrons_w_PhiV_Cut->Fill(positroninfo[x].pt);
			  h_eta__positrons_w_PhiV_Cut->Fill(positroninfo[x].eta);
			  h_phi__positrons_w_PhiV_Cut->Fill(positroninfo[x].phi);
		  }
	  }// end: for(x=0;x<num_positron;x++)
	  for (y = 0; y < num_electron; y++)
	  {
		  h_pT__electrons->Fill(electroninfo[y].pt);
		  h_eta__electrons->Fill(electroninfo[y].eta);
		  h_phi__electrons->Fill(electroninfo[y].phi);
		  if (!electroninfo[y].isPhotonicE)
		  {
			  h_pT__electrons_w_PhiV_Cut->Fill(electroninfo[y].pt);
			  h_eta__electrons_w_PhiV_Cut->Fill(electroninfo[y].eta);
			  h_phi__electrons_w_PhiV_Cut->Fill(electroninfo[y].phi);
		  }
	  }//

      // Event Mixing
      for(Int_t iBufferEvent=0;iBufferEvent<nEventsInBuffer[magBufferIndex][cenBufferIndex][vzBufferIndex];iBufferEvent++)
      {
        for(x=0;x<current_nPositron;x++)// unlike-sign +-
        {
			for (y = 0; y < buffer_nEMinus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent]; y++)
			{
				eepair = current_positron[x] + buffer_eMinus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent][y];
				//Double_t angleV = getPhiVAngle(current_positron[x], buffer_eMinus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent][y], 1, -1);// 注意参数1、-1的选取要求
				//Double_t angleVcut = fphiVcut->Eval(eepair.M()); // 根据fphiVcut关于pair-M的函数取值
				//if (eepair.M() > anaCuts::PhiVCutMRange || angleV > angleVcut)
				//{
					if (fabs(eepair.Rapidity()) <= 1)
					{
						h_Mee_Pt_Cen__unlikeMixed->Fill(eepair.M(), eepair.Perp(), mCen16, reWeight);
					}
				//}
			}
        }        
        for(x=0;x<current_nElectron;x++)// unlike-sign -+
        {
			for (y = 0; y < buffer_nEPlus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent]; y++)
			{
				eepair = current_electron[x] + buffer_ePlus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent][y];
				//Double_t angleV = getPhiVAngle(current_electron[x], buffer_ePlus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent][y], 1, -1);// 注意参数1、-1的选取要求
				//Double_t angleVcut = fphiVcut->Eval(eepair.M()); // 根据fphiVcut关于pair-M的函数取值
				//if (eepair.M() > anaCuts::PhiVCutMRange || angleV > angleVcut)
				//{
				if (fabs(eepair.Rapidity()) <= 1)
				{
					h_Mee_Pt_Cen__unlikeMixed->Fill(eepair.M(), eepair.Perp(), mCen16, reWeight);
				}
				//}
			}
        }          
        for(x=0;x<current_nPositron;x++)// like-sign ++
        {
			for (y = 0; y < buffer_nEPlus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent]; y++)
			{
				eepair = current_positron[x] + buffer_ePlus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent][y];
				//Double_t angleV = getPhiVAngle(current_positron[x], buffer_ePlus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent][y], 1, -1);// 注意参数1、-1的选取要求
				//Double_t angleVcut = fphiVcut->Eval(eepair.M()); // 根据fphiVcut关于pair-M的函数取值
				//if (eepair.M() > anaCuts::PhiVCutMRange || angleV > angleVcut)
				//{
					if (fabs(eepair.Rapidity()) <= 1)
					{
						h_Mee_Pt_Cen__likeppMixed->Fill(eepair.M(), eepair.Perp(), mCen16, reWeight);
					}
				//}
			}
        }
        for(x=0;x<current_nElectron;x++)// like-sign --
        {
			for (y = 0; y < buffer_nEMinus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent]; y++)
			{
				eepair = current_electron[x] + buffer_eMinus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent][y];
				//Double_t angleV = getPhiVAngle(current_electron[x], buffer_eMinus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent][y], 1, -1);// 注意参数1、-1的选取要求
				//Double_t angleVcut = fphiVcut->Eval(eepair.M()); // 根据fphiVcut关于pair-M的函数取值
				//if (eepair.M() > anaCuts::PhiVCutMRange || angleV > angleVcut)
				//{
					if (fabs(eepair.Rapidity()) <= 1)
					{
						h_Mee_Pt_Cen__likemmMixed->Fill(eepair.M(), eepair.Perp(), mCen16, reWeight);
					}
				//}
			}
        }
      }// End Mixed Event

      copyCurrentToBuffer();
    }// Good Event
  }// Good Run

  return kStOK;
}// end Make
Int_t StPicoDstarMixedMaker::Finish()
{
	mFile->cd();// 将当前 ROOT 的输出目录切换到 mFile 所指向的 ROOT 文件

	h_RunNum->Write();
	h_passEvtcut->Write();
	h_passTrkcut->Write();
	h_cen->Write();
	h_cen__reWeight->Write();
	h_RefMult->Write();

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
	h_mRefMult->Write();
	h_pDca_Eta_NHitsFit->Write();
	h_pDca_Pt_Eta->Write();
	// track level QA
	h_nHitsFit->Write();
	h_nHitsPoss->Write();
	h_nHitsDEdx->Write();
	h_nHitsFit_Pt_Eta->Write();
	h_nHitsDEdx_Pt_Eta->Write();

	h_invBeta_P__TOFMatch->Write();
	h_dEdx_Pc->Write();
	h_m2->Write();
	h_m2_Pc->Write();
	h_pDca->Write();
	h_ppT->Write();
	h_pP->Write();
	h_pP_ppT->Write();
	h_gPt->Write();
	h_pEta->Write();
	h_pPhi->Write();
	h_ppTc_pEta->Write();
	h_ppTc_pPhi->Write();

	h_pT__TOFMatch->Write();
	h_Eta__TOFMatch->Write();
	h_Phi__TOFMatch->Write();

	h_nSigmaElectron_P->Write();
	h_nSigmaEcorr_P->Write();
	h_nSigmaPion_P->Write();
	h_nSigmaPicorr_P->Write();
	h_nSigmaKaon_P->Write();
	h_nSigmaProton_P->Write();

	//h_Pt_Cen_nSigmaPion->Write();
	//h_Eta_Cen_nSigmaPion->Write();
	//h_Phi_Cen_nSigmaPion->Write();
	//h_Pt_Cen_nSigmaPicorr->Write();
	//h_Eta_Cen_nSigmaPicorr->Write();
	//h_Phi_Cen_nSigmaPicorr->Write();

	//h_Pt_Cen_nSigmaE__PureE->Write();
	//h_Eta_Cen_nSigmaE__PureE->Write();
	//h_Phi_Cen_nSigmaE__PureE->Write();
	h_Pt_Cen_nSigmaE->Write();
	h_Eta_Cen_nSigmaE->Write();
	h_Phi_Cen_nSigmaE->Write();
	//h_Pt_Cen_nSigmaEcorr->Write();
	//h_Eta_Cen_nSigmaEcorr->Write();
	//h_Phi_Cen_nSigmaEcorr->Write();

	//group 1 2 3
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
	h_nSigmaKaon_P__2->Write();
	h_nSigmaPion_P__3->Write();
	h_nSigmaKaon_P__3->Write();
	h_nSigmaProton_P__2->Write();
	h_nSigmaEcorr_P__TOFMatch->Write();
	h_pT_Eta_Phi__EIDcut_3_lowP_e->Write();
	h_pT_Eta_Phi__EIDcut_3_lowP_p->Write();

	// phiV cut
	h_eNumber_Cen->Write();
	h_pT__positrons->Write();
	h_eta__positrons->Write();
	h_phi__positrons->Write();
	h_pT__positrons_w_PhiV_Cut->Write();
	h_eta__positrons_w_PhiV_Cut->Write();
	h_phi__positrons_w_PhiV_Cut->Write();
	h_pT__electrons->Write();
	h_eta__electrons->Write();
	h_phi__electrons->Write();
	h_pT__electrons_w_PhiV_Cut->Write();
	h_eta__electrons_w_PhiV_Cut->Write();
	h_phi__electrons_w_PhiV_Cut->Write();

	h_Rapidity__unlikeSame->Write();
	h_Mee_PhiV__unlikeSame->Write();
	h_Mee__unlikeSame->Write();
	h_Mee__unlikeSame__w_PhiV_Cut->Write();

	h_Mee_Pt_Cen__unlikeSame->Write();
	h_Mee_Pt_Cen__likemm->Write();
	h_Mee_Pt_Cen__likepp->Write();

	h_Mee_Pt_Cen__unlikeMixed->Write();
	h_Mee_Pt_Cen__likemmMixed->Write();
	h_Mee_Pt_Cen__likeppMixed->Write();

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

Double_t StPicoDstarMixedMaker::getPhiVAngle(TLorentzVector e1, TLorentzVector e2, Int_t q1, Int_t q2) const
{
	Double_t pt1 = e1.Pt();
	Double_t eta1 = e1.Eta();
	Double_t phi1 = e1.Phi();

	Double_t pt2 = e2.Pt();
	Double_t eta2 = e2.Eta();
	Double_t phi2 = e2.Phi();

	TVector3 e1Mom,e2Mom;
	if(q1>0&&q2<0){
		e2Mom.SetPtEtaPhi(pt1,eta1,phi1);// e+
		e1Mom.SetPtEtaPhi(pt2,eta2,phi2);// e-
	}else if(q1<0&&q2>0){
		e2Mom.SetPtEtaPhi(pt2,eta2,phi2);// e+
		e1Mom.SetPtEtaPhi(pt1,eta1,phi1);// e-
	}else if(q1==q2&&TMath::Abs(q1)==1){
		Double_t ran = gRandom->Uniform(-1,1);
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
	if(mBfield<0.) mN = -1.;
	if(mBfield>0.) mN = 1.;

	TVector3 pu=e1Mom+e2Mom;
	TVector3 pv=e1Mom.Cross(e2Mom);
	TVector3 pw=pu.Cross(pv);
	TVector3 pnz(0.,0.,mN);
	TVector3 pwc=pu.Cross(pnz);

	Double_t angleV = pw.Angle(pwc);

	return angleV;
}

Int_t StPicoDstarMixedMaker::getRefmult(StPicoDst const* const picoDst, StPicoEvent const* const picoEvent) const
{
	Int_t refMult = 0;
	for (UInt_t iTrk = 0; iTrk < picoDst->numberOfTracks(); iTrk++)// UInt_t是无符号整型
	{
		StPicoTrack* track = picoDst->track(iTrk);
		if (!track) continue;
		if (!track->isPrimary()) continue; // 只考虑主要径迹
		if (TMath::Abs(track->pMom().Eta()) > 0.5) continue;
		refMult++;
	}
	return refMult;
}

Double_t StPicoDstarMixedMaker::getNSigmaECorr(TVector3 mom) const
{
	return 0;
}
Double_t StPicoDstarMixedMaker::getNSigmaPiKPCorr(Int_t num_variable, TVector3 mom) const
{
	return 0;
}

void StPicoDstarMixedMaker::copyCurrentToBuffer()
{
	//判断是否已满
	if (nEventsInBuffer[magBufferIndex][cenBufferIndex][vzBufferIndex] >= kMaxEventsInBuffer)
		bufferFullFlag[magBufferIndex][cenBufferIndex][vzBufferIndex] = kTRUE;

	Int_t eventIndex = -1;
	if (bufferFullFlag[magBufferIndex][cenBufferIndex][vzBufferIndex])
		eventIndex = gRandom->Integer(kMaxEventsInBuffer);
	else 
		eventIndex = nEventsInBuffer[magBufferIndex][cenBufferIndex][vzBufferIndex];
	buffer_nEPlus[magBufferIndex][cenBufferIndex][vzBufferIndex][eventIndex] = current_nPositron;
	for (Int_t i = 0; i < current_nPositron; i++)
		buffer_ePlus[magBufferIndex][cenBufferIndex][vzBufferIndex][eventIndex][i] = current_positron[i];
	buffer_nEMinus[magBufferIndex][cenBufferIndex][vzBufferIndex][eventIndex] = current_nElectron;
	for (Int_t i = 0; i < current_nElectron; i++)
		buffer_eMinus[magBufferIndex][cenBufferIndex][vzBufferIndex][eventIndex][i] = current_electron[i];
	if (nEventsInBuffer[magBufferIndex][cenBufferIndex][vzBufferIndex] < kMaxEventsInBuffer)
		nEventsInBuffer[magBufferIndex][cenBufferIndex][vzBufferIndex]++;
}