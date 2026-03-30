//---------------------从OO.root中提取直方图，并进行设置更改和元素添加----------------------
#include "someFunction.h"
void DrawNSigmaEPiKP(TString inFilename = "roots/6_20260330_Au200_3groups_MoreMaxProFile_P23ie.root", Int_t number = 6) //
{
	// 从root文件中导入待拟合的直方图
	TFile *inFile = new TFile(inFilename);
	if (!inFile)
	{
		cout << "The input file is not found! Exiting..." << endl;
		return;
	}
	TH1F *h_passEvtcut = (TH1F *)inFile->Get("h_passEvtcut");
	TH1F *h_passTrkcut = (TH1F *)inFile->Get("h_passTrkcut");
	TH2F *h_VpdVz_Vz = (TH2F *)inFile->Get("h_VpdVz_Vz");
	TH2F *h_Vx_Vy = (TH2F *)inFile->Get("h_Vx_Vy");
	TH1F *h_VpdVzmVz = (TH1F *)inFile->Get("h_VpdVzmVz");
	TH1F *h_Vr = (TH1F *)inFile->Get("h_Vr");

	TH3F *h_nHitsFit_Pt_Eta = (TH3F *)inFile->Get("h_nHitsFit_Pt_Eta");
	TH2F *h_nHitsFit_Pt = (TH2F *)h_nHitsFit_Pt_Eta->Project3D("zx");
	TH2F *h_nHitsFit_Eta = (TH2F *)h_nHitsFit_Pt_Eta->Project3D("zy");

	TH3F *h_nHitsDEdx_Pt_Eta = (TH3F *)inFile->Get("h_nHitsDEdx_Pt_Eta");
	TH2F *h_nHitsDEdx_Pt = (TH2F *)h_nHitsDEdx_Pt_Eta->Project3D("zx");
	TH2F *h_nHitsDEdx_Eta = (TH2F *)h_nHitsDEdx_Pt_Eta->Project3D("zy");

	TH3F *h_pDca_Pt_Eta = (TH3F *)inFile->Get("h_pDca_Pt_Eta");
	TH2F *h_pDca_Pt = (TH2F *)h_pDca_Pt_Eta->Project3D("xy");
	TH2F *h_pDca_Eta = (TH2F *)h_pDca_Pt_Eta->Project3D("xz");

	TH1F *h_ppT = (TH1F *)inFile->Get("h_ppT");
	TH1F *h_pEta = (TH1F *)inFile->Get("h_pEta");
	TH1F *h_pPhi = (TH1F *)inFile->Get("h_pPhi");
	TH2F *h_pP_ppT = (TH2F *)inFile->Get("h_pP_ppT");
	TH2F *h_ppTc_pEta = (TH2F *)inFile->Get("h_ppTc_pEta");

	TH2F *h_nSigmaElectron_P = (TH2F *)inFile->Get("h_nSigmaElectron_P");
	TH2F *h_nSigmaEcorr_P = (TH2F *)inFile->Get("h_nSigmaEcorr_P");
	TH2F *h_nSigmaPion_P = (TH2F *)inFile->Get("h_nSigmaPion_P");
	TH2F *h_nSigmaPicorr_P = (TH2F *)inFile->Get("h_nSigmaPicorr_P");
	// total
	TH2F *h_nSigmaElectron_P__EIDcut_total = (TH2F *)inFile->Get("h_nSigmaElectron_P__EIDcut_total");
	// group 1
	TH1F *h_pT__TOFMatch = (TH1F *)inFile->Get("h_pT__TOFMatch");
	TH1F *h_Eta__TOFMatch = (TH1F *)inFile->Get("h_Eta__TOFMatch");
	TH1F *h_Phi__TOFMatch = (TH1F *)inFile->Get("h_Phi__TOFMatch");
	TH2F *h_nSigmaElectron_P__1 = (TH2F *)inFile->Get("h_nSigmaElectron_P__1");
	TH2F *h_invBeta_P__TOFMatch = (TH2F *)inFile->Get("h_invBeta_P__TOFMatch");
	TH2F *h_nSigmaElectron_P__TOFMatch = (TH2F *)inFile->Get("h_nSigmaElectron_P__TOFMatch");
	TH2F *h_nSigmaElectron_P__EIDcut_1 = (TH2F *)inFile->Get("h_nSigmaElectron_P__EIDcut_1");
	// group 2
	TH2F *h_nSigmaElectron_P__2 = (TH2F *)inFile->Get("h_nSigmaElectron_P__2");
	TH2F *h_nSigmaKaon_P__2 = (TH2F *)inFile->Get("h_nSigmaKaon_P__2");
	TH2F *h_nSigmaProton_P__2 = (TH2F *)inFile->Get("h_nSigmaProton_P__2");
	TH2F *h_nSigmaElectron_P__PID_2 = (TH2F *)inFile->Get("h_nSigmaElectron_P__PID_2");
	TH2F *h_nSigmaElectron_P__EIDcut_2 = (TH2F *)inFile->Get("h_nSigmaElectron_P__EIDcut_2");
	TH2F *h_nSigmaElectron_Eta__EIDcut_2 = (TH2F *)inFile->Get("h_nSigmaElectron_Eta__EIDcut_2");
	// group 3
	TH2F *h_nSigmaElectron_P__3 = (TH2F *)inFile->Get("h_nSigmaElectron_P__3");
	TH2F *h_nSigmaPion_P__3 = (TH2F *)inFile->Get("h_nSigmaPion_P__3");
	TH2F *h_nSigmaKaon_P__3 = (TH2F *)inFile->Get("h_nSigmaKaon_P__3");
	TH2F *h_nSigmaElectron_P__PID_3 = (TH2F *)inFile->Get("h_nSigmaElectron_P__PID_3");
	TH2F *h_nSigmaElectron_P__EIDcut_3 = (TH2F *)inFile->Get("h_nSigmaElectron_P__EIDcut_3");
	TH2F *h_nSigmaElectron_Eta__EIDcut_3_lowP = (TH2F *)inFile->Get("h_nSigmaElectron_Eta__EIDcut_3_lowP");
	TH2F *h_nSigmaElectron_Eta__EIDcut_3_highP = (TH2F *)inFile->Get("h_nSigmaElectron_Eta__EIDcut_3_highP");
	TH2F *h_nSigmaElectron_P__EIDcutcut_total = (TH2F *)inFile->Get("h_nSigmaElectron_P__EIDcutcut_total");
	TH3F *h_pT_Eta_Phi__EIDcut_3_lowP_e = (TH3F *)inFile->Get("h_pT_Eta_Phi__EIDcut_3_lowP_e");
	TH1F *h_pT__3_lowP_e = (TH1F *)h_pT_Eta_Phi__EIDcut_3_lowP_e->Project3D("x");
	TH1F *h_eta__3_lowP_e = (TH1F *)h_pT_Eta_Phi__EIDcut_3_lowP_e->Project3D("y");
	TH1F *h_phi__3_lowP_e = (TH1F *)h_pT_Eta_Phi__EIDcut_3_lowP_e->Project3D("z");
	TH3F *h_pT_Eta_Phi__EIDcut_3_lowP_p = (TH3F *)inFile->Get("h_pT_Eta_Phi__EIDcut_3_lowP_p");
	TH1F *h_pT__3_lowP_p = (TH1F *)h_pT_Eta_Phi__EIDcut_3_lowP_p->Project3D("x");
	TH1F *h_eta__3_lowP_p = (TH1F *)h_pT_Eta_Phi__EIDcut_3_lowP_p->Project3D("y");
	TH1F *h_phi__3_lowP_p = (TH1F *)h_pT_Eta_Phi__EIDcut_3_lowP_p->Project3D("z");
	TH1F *h_pT_ratio_EperP_3 = (TH1F *)h_pT__3_lowP_e->Clone("h_pT_ratio_EperP_3");
	h_pT_ratio_EperP_3->SetTitle("e^{-}/e^{+} p_{T} ratio;p_{T} [GeV/c];Ratio");
	h_pT_ratio_EperP_3->Divide(h_pT__3_lowP_p);
	TH1F *h_eta_ratio_EperP_3 = (TH1F *)h_eta__3_lowP_e->Clone("h_eta_ratio_EperP_3");
	h_eta_ratio_EperP_3->SetTitle("e^{-}/e^{+} #eta ratio;#eta [GeV/c];Ratio");
	h_eta_ratio_EperP_3->Divide(h_eta__3_lowP_p);
	TH1F *h_phi_ratio_EperP_3 = (TH1F *)h_phi__3_lowP_e->Clone("h_phi_ratio_EperP_3");
	h_phi_ratio_EperP_3->SetTitle("e^{-}/e^{+} #phi ratio;#phi [GeV/c];Ratio");
	h_phi_ratio_EperP_3->Divide(h_phi__3_lowP_p);
	//phiV cut
	TH2F* h_Mee_PhiV__ranComb = (TH2F*)inFile->Get("h_Mee_PhiV__unlikeSame");//h_Mee_PhiV__ranComb
	TH1F* h_Mee__ranComb = (TH1F*)inFile->Get("h_Mee__unlikeSame");//h_Mee__ranComb
	TH1F* h_Mee__ranComb__w_PhiV_Cut = (TH1F*)inFile->Get("h_Mee__unlikeSame__w_PhiV_Cut");//h_Mee__ranComb__w_PhiV_Cut
	auto fphiVcut = new TF1("fphiVcut", "0.84326*exp(-49.4819*x)-0.996609*x+0.19801", 0, 1.0);
	fphiVcut->SetMaximum(1.0);
	fphiVcut->SetNpx(1000);
	// sum
	TH2F *h_eNumber_Cen = (TH2F *)inFile->Get("h_eNumber_Cen");
	TH1F *h_pT__electrons = (TH1F *)inFile->Get("h_pT__electrons");
	TH1F *h_eta__electrons = (TH1F *)inFile->Get("h_eta__electrons");
	TH1F *h_phi__electrons = (TH1F *)inFile->Get("h_phi__electrons");
	TH1F *h_pT__positrons = (TH1F *)inFile->Get("h_pT__positrons");
	TH1F *h_eta__positrons = (TH1F *)inFile->Get("h_eta__positrons");
	TH1F *h_phi__positrons = (TH1F *)inFile->Get("h_phi__positrons");
	TH1F *h_pT_ratio_EperP = (TH1F *)h_pT__electrons->Clone("h_pT_ratio_EperP");
	h_pT_ratio_EperP->SetTitle("e^{-}/e^{+} p_{T} ratio;p_{T} [GeV/c];Ratio");
	h_pT_ratio_EperP->Divide(h_pT__positrons);
	TH1F *h_eta_ratio_EperP = (TH1F *)h_eta__electrons->Clone("h_eta_ratio_EperP");
	h_eta_ratio_EperP->SetTitle("e^{-}/e^{+} #eta ratio;#eta [GeV/c];Ratio");
	h_eta_ratio_EperP->Divide(h_eta__positrons);
	TH1F *h_phi_ratio_EperP = (TH1F *)h_phi__electrons->Clone("h_phi_ratio_EperP");
	h_phi_ratio_EperP->SetTitle("e^{-}/e^{+} #phi ratio;#phi [GeV/c];Ratio");
	h_phi_ratio_EperP->Divide(h_phi__positrons);

	if (0) // Draw clear plot
	{
		TCanvas *c_sum = new TCanvas("c_sum", "c_sum", 1400, 600);
		c_sum->Divide(2);

		//c_sum->cd(1);
		//gPad->SetLogz(1);
		//gPad->SetLeftMargin(0.12);
		//gPad->SetRightMargin(0.12);
		//gStyle->SetOptStat(0);
		//h_nSigmaElectron_P->GetYaxis()->SetRangeUser(-10, 10);
		//h_nSigmaElectron_P->Draw("col z");

		//c_sum->cd(2);
		//gPad->SetLogz(1);
		//gPad->SetLeftMargin(0.12);
		//gPad->SetRightMargin(0.12);
		//gStyle->SetOptStat(0);
		//h_nSigmaElectron_P__EIDcut_total->Draw("col z");

		//c_sum->cd(1);
		//gPad->SetLeftMargin(0.12);
		//gPad->SetRightMargin(0.12);
		//gStyle->SetOptStat(0);
		////h_passEvtcut->SetMaximum(5.5e8);
		//h_passEvtcut->DrawClone("TEXT");
		//h_passEvtcut->DrawClone("same");

		//c_sum->cd(1);
		//gPad->SetLeftMargin(0.12);
		//gPad->SetRightMargin(0.12);
		//gStyle->SetOptStat(0);
		//h_passTrkcut->DrawClone("TEXT0");
		//h_passTrkcut->DrawClone("same");

		//c_sum->cd(1);
		//gPad->SetLogz(1);
		//gPad->SetLeftMargin(0.12);
		//gPad->SetRightMargin(0.12);
		//gStyle->SetOptStat(0);
		//h_Vx_Vy->DrawClone("col z");

		//c_sum->cd(2);
		//gPad->SetLogz(1);
		//gPad->SetLeftMargin(0.12);
		//gPad->SetRightMargin(0.12);
		//gStyle->SetOptStat(0);
		//h_VpdVz_Vz->GetYaxis()->SetTitleOffset(1.5);
		//h_VpdVz_Vz->DrawClone("col z");

		//c_sum->cd(1);
		//gPad->SetLogz(1);
		//gPad->SetLeftMargin(0.12);
		//gPad->SetRightMargin(0.12);
		//gStyle->SetOptStat(0);
		////h_nHitsFit_Pt->SetTitle("nHitsFit vs p_{T}");
		//h_nHitsFit_Pt->GetXaxis()->SetRangeUser(0.0, 3.0);
		//h_nHitsFit_Pt->DrawClone("col z");

		//c_sum->cd(2);
		//gPad->SetLogz(1);
		//gPad->SetLeftMargin(0.12);
		//gPad->SetRightMargin(0.12);
		//gStyle->SetOptStat(0);
		////h_nHitsDEdx_Pt->SetTitle("nHitsDEdx vs p_{T}");
		//h_nHitsDEdx_Pt->GetXaxis()->SetRangeUser(0.0, 3.0);
		//h_nHitsDEdx_Pt->DrawClone("col z");

		//c_sum->cd(1);
		//gPad->SetLogz(1);
		//gPad->SetLeftMargin(0.12);
		//gPad->SetRightMargin(0.12);
		//gStyle->SetOptStat(0);
		//h_nSigmaElectron_Eta__EIDcut_3_lowP->DrawClone("col z");

		//c_sum->cd(2);
		//gPad->SetLogz(1);
		//gPad->SetLeftMargin(0.12);
		//gPad->SetRightMargin(0.12);
		//gStyle->SetOptStat(0);
		//h_nSigmaElectron_Eta__EIDcut_3_highP->DrawClone("col z");

		c_sum->cd(1);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_VpdVzmVz->GetXaxis()->SetRangeUser(-10, 10);
		h_VpdVzmVz->DrawClone("col z");

		c_sum->cd(2);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_Vr->GetXaxis()->SetRangeUser(0, 3.0);
		h_Vr->DrawClone("");

		c_sum->SaveAs(Form("roots/%d_VzVpdVz.png", number));
	}

	if (0) // nsigmaElectron__EID, h_pP_ppT, h_passEvtcut, h_passTrkcut
	{
		// 误差条，设置Marker形状颜色，设置线条颜色，设置图例，设置坐标轴标题，设置对数Y坐标）
		TCanvas *c_temp = new TCanvas("c_temp", "c_temp", 200, 100, 1300, 1200);
		c_temp->Divide(4, 4);

		c_temp->cd(1);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);

		c_temp->cd(2);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);

		c_temp->cd(3);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaElectron_Eta__EIDcut_3_lowP->DrawClone("col z");

		c_temp->cd(4);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaElectron_Eta__EIDcut_3_highP->DrawClone("col z");

		c_temp->cd(5);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_Vx_Vy->DrawClone("col z");

		c_temp->cd(6);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_VpdVz_Vz->GetYaxis()->SetTitleOffset(1.5);
		h_VpdVz_Vz->DrawClone("col z");

		c_temp->cd(7);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_passEvtcut->SetMaximum(5.5e8);
		h_passEvtcut->DrawClone("TEXT");
		h_passEvtcut->DrawClone("same");

		c_temp->cd(8);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_passTrkcut->DrawClone("TEXT0");
		h_passTrkcut->DrawClone("same");

		c_temp->cd(9);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		//h_nHitsFit_Pt->SetTitle("nHitsFit vs p_{T}");
		h_nHitsFit_Pt->GetXaxis()->SetRangeUser(0.0, 3.0);
		h_nHitsFit_Pt->DrawClone("col z");

		c_temp->cd(10);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		//h_nHitsFit_Eta->SetTitle("nHitsFit vs #eta");
		h_nHitsFit_Eta->DrawClone("col z");

		c_temp->cd(11);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		//h_nHitsDEdx_Pt->SetTitle("nHitsDEdx vs p_{T}");
		h_nHitsDEdx_Pt->GetXaxis()->SetRangeUser(0.0, 3.0);
		h_nHitsDEdx_Pt->DrawClone("col z");

		c_temp->cd(12);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		//h_nHitsDEdx_Eta->SetTitle("nHitsDEdx vs #eta");
		h_nHitsDEdx_Eta->DrawClone("col z");

		c_temp->cd(13);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);

		c_temp->cd(14);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		//h_pDca_Pt->SetTitle("DCA vs p_{T}");
		h_pDca_Pt->GetXaxis()->SetRangeUser(0.0, 3.0);
		h_pDca_Pt->DrawClone("col z");

		c_temp->cd(15);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		//h_pDca_Eta->SetTitle("DCA vs #eta");
		h_pDca_Eta->DrawClone("col z");

		c_temp->cd(16);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);

		c_temp->SaveAs(Form("roots/%d_Track_and_Event_Check.png", number));
	}
	if (0) // track QA and TOF track check
	{
		TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
		c1->Divide(3, 2);
		c1->cd(1);
		gPad->SetLogy(1);
		h_ppT->GetXaxis()->SetRangeUser(0., 3.0);
		// h_ppT->GetYaxis()->SetRangeUser(1., 1e8);
		h_ppT->Draw();

		c1->cd(2);
		gPad->SetLogy(1);
		// h_pEta->GetYaxis()->SetRangeUser(1., 1e8);
		h_pEta->Draw();

		c1->cd(3);
		gPad->SetLogy(1);
		h_pPhi->GetXaxis()->SetRangeUser(-TMath::Pi() * 4 / 3, TMath::Pi() * 4 / 3);
		// h_pPhi->GetYaxis()->SetRangeUser(1., 1e8);
		h_pPhi->GetXaxis()->SetNdivisions(-108);						   // 负号 = 用户自定义，202 = π 刻度
		h_pPhi->GetXaxis()->ChangeLabel(1, -1, -1, -1, -1, -1, "-4#pi/3"); // 让 ROOT 自动写 kπ
		h_pPhi->GetXaxis()->ChangeLabel(2, -1, -1, -1, -1, -1, "-#pi");
		h_pPhi->GetXaxis()->ChangeLabel(3, -1, -1, -1, -1, -1, "-2#pi/3");
		h_pPhi->GetXaxis()->ChangeLabel(4, -1, -1, -1, -1, -1, "-#pi/3");
		h_pPhi->GetXaxis()->ChangeLabel(5, -1, -1, -1, -1, -1, "0");
		h_pPhi->GetXaxis()->ChangeLabel(-4, -1, -1, -1, -1, -1, "#pi/3");
		h_pPhi->GetXaxis()->ChangeLabel(-3, -1, -1, -1, -1, -1, "2#pi/3");
		h_pPhi->GetXaxis()->ChangeLabel(-2, -1, -1, -1, -1, -1, "#pi");
		h_pPhi->GetXaxis()->ChangeLabel(-1, -1, -1, -1, -1, -1, "4#pi/3"); // 让 ROOT 自动写 kπ
		h_pPhi->Draw();

		c1->cd(4);
		gPad->SetLogy(1);
		h_pT__TOFMatch->GetXaxis()->SetRangeUser(0., 3.0);
		// h_pT__TOFMatch->GetYaxis()->SetRangeUser(1., 1e8);
		h_pT__TOFMatch->Draw();

		c1->cd(5);
		gPad->SetLogy(1);
		// h_Eta__TOFMatch->GetYaxis()->SetRangeUser(1., 1e8);
		h_Eta__TOFMatch->Draw();

		c1->cd(6);
		gPad->SetLogy(1);
		h_Phi__TOFMatch->GetXaxis()->SetRangeUser(-TMath::Pi() * 4 / 3, TMath::Pi() * 4 / 3);
		// h_Phi__TOFMatch->GetYaxis()->SetRangeUser(1., 1e8);
		h_Phi__TOFMatch->GetXaxis()->SetNdivisions(-108);							// 负号 = 用户自定义，202 = π 刻度
		h_Phi__TOFMatch->GetXaxis()->ChangeLabel(1, -1, -1, -1, -1, -1, "-4#pi/3"); // 让 ROOT 自动写 kπ
		h_Phi__TOFMatch->GetXaxis()->ChangeLabel(2, -1, -1, -1, -1, -1, "-#pi");
		h_Phi__TOFMatch->GetXaxis()->ChangeLabel(3, -1, -1, -1, -1, -1, "-2#pi/3");
		h_Phi__TOFMatch->GetXaxis()->ChangeLabel(4, -1, -1, -1, -1, -1, "-#pi/3");
		h_Phi__TOFMatch->GetXaxis()->ChangeLabel(5, -1, -1, -1, -1, -1, "0");
		h_Phi__TOFMatch->GetXaxis()->ChangeLabel(-4, -1, -1, -1, -1, -1, "#pi/3");
		h_Phi__TOFMatch->GetXaxis()->ChangeLabel(-3, -1, -1, -1, -1, -1, "2#pi/3");
		h_Phi__TOFMatch->GetXaxis()->ChangeLabel(-2, -1, -1, -1, -1, -1, "#pi");
		h_Phi__TOFMatch->GetXaxis()->ChangeLabel(-1, -1, -1, -1, -1, -1, "4#pi/3"); // 让 ROOT 自动写 kπ
		h_Phi__TOFMatch->Draw();

		c1->SaveAs(Form("roots/%d_TrackQA_and_TrackTOFMatch.png", number));
	}
	if (0) // e+ e- track QA
	{
		h_pT__positrons->SetLineColor(kRed);
		h_eta__positrons->SetLineColor(kRed);
		h_phi__positrons->SetLineColor(kRed);
		h_pT__3_lowP_p->SetLineColor(kRed);
		h_eta__3_lowP_p->SetLineColor(kRed);
		h_phi__3_lowP_p->SetLineColor(kRed);

		TCanvas *c_eTrack = new TCanvas("c_eTrack", "c_eTrack", 600, 900);
		c_eTrack->Divide(3, 4);

		c_eTrack->cd(1);
		gPad->SetLogy(1);
		TLegend *leg1 = new TLegend(0.65, 0.75, 0.88, 0.88);
		leg1->AddEntry(h_pT__electrons, "e^{-} (electrons)", "l");
		leg1->AddEntry(h_pT__positrons, "e^{+} (positrons)", "l");
		leg1->SetBorderSize(0);
		leg1->SetFillStyle(0);
		// h_pT__electrons->GetXaxis()->SetRangeUser(0., 0.5);
		h_pT__electrons->Draw();
		h_pT__positrons->Draw("same");
		leg1->Draw("same");

		c_eTrack->cd(2);
		gPad->SetLogy(1);
		// h_eta__electrons->GetXaxis()->SetRangeUser(-2.0, 2.0);
		h_eta__electrons->Draw();
		h_eta__positrons->Draw("same");

		c_eTrack->cd(3);
		gPad->SetLogy(1);
		// h_phi__electrons->GetXaxis()->SetRangeUser(-3.2, 3.2);
		h_phi__electrons->Draw();
		h_phi__positrons->Draw("same");

		c_eTrack->cd(4);
		gPad->SetLogy(0);
		gStyle->SetOptStat(0);
		h_pT_ratio_EperP->SetMaximum(2);
		h_pT_ratio_EperP->SetMinimum(0.5);
		h_pT_ratio_EperP->Draw();

		c_eTrack->cd(5);
		gPad->SetLogy(0);
		gStyle->SetOptStat(0);
		h_eta_ratio_EperP->SetMaximum(2);
		h_eta_ratio_EperP->SetMinimum(0.5);
		h_eta_ratio_EperP->Draw();

		c_eTrack->cd(6);
		gPad->SetLogy(0);
		gStyle->SetOptStat(0);
		h_phi_ratio_EperP->SetMaximum(2);
		h_phi_ratio_EperP->SetMinimum(0.5);
		h_phi_ratio_EperP->Draw();

		c_eTrack->cd(7);
		gPad->SetLogy(1);
		TLegend *leg2 = new TLegend(0.65, 0.75, 0.88, 0.88);
		leg2->AddEntry(h_pT__3_lowP_e, "e^{-} (electrons)", "l");
		leg2->AddEntry(h_pT__3_lowP_p, "e^{+} (positrons)", "l");
		leg2->SetBorderSize(0);
		leg2->SetFillStyle(0);
		// h_pT__3_lowP_e->GetXaxis()->SetRangeUser(0., 0.5);
		h_pT__3_lowP_e->Draw();
		h_pT__3_lowP_p->Draw("same");
		leg2->Draw("same");

		c_eTrack->cd(8);
		gPad->SetLogy(1);
		// h_eta__3_lowP_e->GetXaxis()->SetRangeUser(-2.0, 2.0);
		h_eta__3_lowP_e->Draw();
		h_eta__3_lowP_p->Draw("same");

		c_eTrack->cd(9);
		gPad->SetLogy(1);
		// h_phi__3_lowP_e->GetXaxis()->SetRangeUser(-3.2, 3.2);
		h_phi__3_lowP_e->Draw();
		h_phi__3_lowP_p->Draw("same");

		c_eTrack->cd(10);
		gPad->SetLogy(0);
		gStyle->SetOptStat(0);
		h_pT_ratio_EperP_3->SetMaximum(2);
		h_pT_ratio_EperP_3->SetMinimum(0.5);
		h_pT_ratio_EperP_3->Draw();

		c_eTrack->cd(11);
		gPad->SetLogy(0);
		gStyle->SetOptStat(0);
		h_eta_ratio_EperP_3->SetMaximum(2);
		h_eta_ratio_EperP_3->SetMinimum(0.5);
		h_eta_ratio_EperP_3->Draw();

		c_eTrack->cd(12);
		gPad->SetLogy(0);
		gStyle->SetOptStat(0);
		h_phi_ratio_EperP_3->SetMaximum(2);
		h_phi_ratio_EperP_3->SetMinimum(0.5);
		h_phi_ratio_EperP_3->Draw();

		c_eTrack->SaveAs(Form("roots/%d_ep_TrackQA.png", number));
		// delete c_eTrack;
	}
	if (0) // EID in group1(pT>0.2, |eta|<1)
	{
		// 设置直方图格式
		// 去除误差条，设置Marker形状颜色，设置线条颜色，设置图例，设置坐标轴标题，设置对数Y坐标）
		//h_nSigmaElectron_P__1->SetTitle("n#sigma_{e} in p_{T}>0.2&|#eta|<0.9;p (GeV/c);n#sigma_{e}");
		h_invBeta_P__TOFMatch->GetYaxis()->SetRangeUser(0.0, 3.5);
		//h_nSigmaElectron_P__TOFMatch->SetTitle("n#sigma_{e} with TOF cut;p (GeV/c);n#sigma_{e}");
		// 画图
		TCanvas *c4 = new TCanvas("c4", "c4", 1200, 800);
		c4->Divide(3, 2);
		c4->cd(1); // nSigmaE
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaElectron_P__1->Draw("col z");

		c4->cd(2);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_invBeta_P__TOFMatch->Draw("col z");
		TLine *line1_beta1 = new TLine(0.0, 1 - 0.02, 5.0, 1 - 0.02);
		line1_beta1->SetLineColor(kRed);
		line1_beta1->Draw("same");
		TLine *line1_beta2 = new TLine(0.0, 1 + 0.02, 5.0, 1 + 0.02);
		line1_beta2->SetLineColor(kRed);
		line1_beta2->Draw("same");

		c4->cd(3);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaElectron_P__TOFMatch->Draw("col z");
		TLine *line1_TOF1 = new TLine(0.0, 3.0, 5.0, 3.0);
		line1_TOF1->SetLineColor(kRed);
		line1_TOF1->Draw("same");
		TLine *line1_TOF2 = new TLine(1.0, -1.0, 5.0, -1.0);
		line1_TOF2->SetLineColor(kRed);
		line1_TOF2->Draw("same");
		TLine *line1_TOF3 = new TLine(0.0, -3.5, 1.0, -1.0);
		line1_TOF3->SetLineColor(kRed);
		line1_TOF3->Draw("same");

		c4->cd(4);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaElectron_P__EIDcut_1->Draw("col z");

		c4->SaveAs(Form("roots/%d_group1.png", number));
	}
	if (0) // EID in group2(pT>0.2, 1.0<|eta|<1.8)
	{
		// 设置直方图格式
		//h_nSigmaElectron_P__2->SetTitle("n#sigma_{e} in p_{T}>0.2&0.9<|#eta|<1.8;p (GeV/c);n#sigma_{e}");
		//h_nSigmaKaon_P__2->SetTitle("nSigmaKaon_P;P(GeV/c);nSigmaKaon");

		// 画图
		TCanvas *c5 = new TCanvas("c5", "c5", 1200, 800);
		c5->Divide(3, 2);
		c5->cd(1); // nSigmaE
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaElectron_P__2->GetYaxis()->SetRangeUser(-10, 10);
		h_nSigmaElectron_P__2->Draw("col z");

		c5->cd(2);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaKaon_P__2->Draw("col z");
		TLine *line2_K1 = new TLine(0.0, 0.0, 1.0, -2.0);
		line2_K1->SetLineColor(kRed);
		line2_K1->Draw("same");
		TLine *line2_K2 = new TLine(1.0, -2.0, 5.0, -2.0);
		line2_K2->SetLineColor(kRed);
		line2_K2->Draw("same");
		TLine *line2_K3 = new TLine(0.0, 3.0, 1.0, 3.0);
		line2_K3->SetLineColor(kRed);
		line2_K3->Draw("same");
		TLine *line2_K4 = new TLine(1.0, 4.5, 5.0, 4.5);
		line2_K4->SetLineColor(kRed);
		line2_K4->Draw("same");
		TLine *line2_K5 = new TLine(1.0, 3.0, 1.0, 4.5);
		line2_K5->SetLineColor(kRed);
		line2_K5->Draw("same");

		c5->cd(3);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaProton_P__2->Draw("col z");

		TLine *line2_P1 = new TLine(0.0, 6.0, 1.5, 4.0);
		line2_P1->SetLineColor(kRed);
		line2_P1->Draw("same");
		TLine *line2_P2 = new TLine(1.5, 4.0, 5.0, 4.0);
		line2_P2->SetLineColor(kRed);
		line2_P2->Draw("same");
		TLine *line2_P3 = new TLine(0.0, 0.0, 1.5, -2.0);
		line2_P3->SetLineColor(kRed);
		line2_P3->Draw("same");
		TLine *line2_P4 = new TLine(1.5, -2.0, 5.0, -2.0);
		line2_P4->SetLineColor(kRed);
		line2_P4->Draw("same");

		c5->cd(4);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaElectron_P__PID_2->Draw("col z");
		TLine *line2_E1 = new TLine(0.0, 3.0, 5.0, 3.0);
		line2_E1->SetLineColor(kRed);
		line2_E1->Draw("same");
		TLine *line2_E2 = new TLine(0.0, -2.25, 1.5, 0.0);
		line2_E2->SetLineColor(kRed);
		line2_E2->Draw("same");
		TLine *line2_E3 = new TLine(1.5, 0.0, 5.0, 0.0);
		line2_E3->SetLineColor(kRed);
		line2_E3->Draw("same");

		c5->cd(5);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaElectron_P__EIDcut_2->Draw("col z");

		c5->cd(6);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaElectron_Eta__EIDcut_2->DrawClone("col z");

		c5->SaveAs(Form("roots/%d_group2.png", number));
	}
	if (0) // EID in group3(pT>0.2, 1.0<|eta|<1.8)
	{
		// 设置直方图格式

		// 画图
		TCanvas *c6 = new TCanvas("c6", "c6", 1200, 800);
		c6->Divide(3, 2);
		c6->cd(1); // nSigmaE
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		//h_nSigmaElectron_P__3->SetTitle("n#sigma_{e} vs p (p_{T}<0.2,|#eta|<1.8);p (GeV/c);n#sigma_{e}");
		h_nSigmaElectron_P__3->GetXaxis()->SetRangeUser(0, 0.7);
		h_nSigmaElectron_P__3->GetYaxis()->SetRangeUser(-10, 10);
		h_nSigmaElectron_P__3->Draw("col z");

		c6->cd(2);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaPion_P__3->GetXaxis()->SetRangeUser(0, 0.7);
		h_nSigmaPion_P__3->Draw("col z");

		TLine *line3_Pi1 = new TLine(0.0, 5.5, 0.7, 1.5);
		line3_Pi1->SetLineColor(kRed);
		line3_Pi1->Draw("same");
		TLine *line3_Pi2 = new TLine(0.0, -1.0, 0.2, -5.0);
		line3_Pi2->SetLineColor(kRed);
		line3_Pi2->Draw("same");
		TLine *line3_Pi3 = new TLine(0.2, -5.0, 0.7, -5.0);
		line3_Pi3->SetLineColor(kRed);
		line3_Pi3->Draw("same");

		c6->cd(3);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaKaon_P__3->GetXaxis()->SetRangeUser(0, 0.7);
		h_nSigmaKaon_P__3->Draw("col z");
		TLine *line3_K1 = new TLine(0.0, 5.0, 0.7, 5.0);
		line3_K1->SetLineColor(kRed);
		line3_K1->Draw("same");
		TLine *line3_K2 = new TLine(0.0, -6.0, 0.7, 3.5);
		line3_K2->SetLineColor(kRed);
		line3_K2->Draw("same");

		c6->cd(4);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaElectron_P__PID_3->GetXaxis()->SetRangeUser(0, 0.7);
		h_nSigmaElectron_P__PID_3->Draw("col z");
		TLine *line3_E1 = new TLine(0.0, 1.0, 0.18, 1.0);
		line3_E1->SetLineColor(kRed);
		line3_E1->Draw("same");
		TLine *line3_E2 = new TLine(0.18, 1.0, 0.18, 3.0);
		line3_E2->SetLineColor(kRed);
		line3_E2->Draw("same");
		TLine *line3_E3 = new TLine(0.18, 3.0, 0.7, 3.0);
		line3_E3->SetLineColor(kRed);
		line3_E3->Draw("same");
		TLine *line3_E4 = new TLine(0.0, -4.0, 0.14, -4.0);
		line3_E4->SetLineColor(kRed);
		line3_E4->Draw("same");
		TLine *line3_E5 = new TLine(0.14, -4.0, 0.14, -2.0);
		line3_E5->SetLineColor(kRed);
		line3_E5->Draw("same");
		TLine *line3_E6 = new TLine(0.14, -2.0, 0.7, -2.0);
		line3_E6->SetLineColor(kRed);
		line3_E6->Draw("same");

		c6->cd(5);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaElectron_P__EIDcut_3->GetXaxis()->SetRangeUser(0, 0.7);
		h_nSigmaElectron_P__EIDcut_3->Draw("col z");

		c6->SaveAs(Form("roots/%d_group3.png", number));
	}
	if (1)//check PhiV cut
	{
		h_Mee__ranComb->SetLineColor(kBlack);
		h_Mee__ranComb->GetYaxis()->SetTitleOffset(1.5);
		h_Mee__ranComb__w_PhiV_Cut->SetLineColor(kBlue);
		h_Mee__ranComb->SetTitle("M_{ee} w/wo #phi_{V} cut;M_{ee} (GeV/c^{2});Counts");

		TCanvas* c_temp = new TCanvas("c_temp", "c_temp", 800, 800);
		c_temp->Divide(2, 2);
		c_temp->cd(1);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_Mee_PhiV__ranComb->GetXaxis()->SetRangeUser(0.0, 0.4);
		h_Mee_PhiV__ranComb->GetYaxis()->SetRangeUser(0.0, 1.0);
		fphiVcut->GetXaxis()->SetRangeUser(0, 0.2);
		h_Mee_PhiV__ranComb->DrawClone("col z");
		fphiVcut->DrawClone("same");
		TPaveText *pt2 = new TPaveText(0.15, 0.79, 0.85, 0.86, "NDC NB");
		pt2->SetFillColorAlpha(0, 0);   // 透明底
		pt2->SetBorderSize(0);
		pt2->SetTextFont(42);
		pt2->SetTextSize(0.032);
		pt2->SetTextAlign(12);
		pt2->AddText("#phi_{V}(M_{ee})=0.8433*Exp(-49.4819*M)-0.9966*M+0.1980");
		pt2->DrawClone("same");

		c_temp->cd(2);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_Mee__ranComb->GetXaxis()->SetRangeUser(0, 0.2);
		h_Mee__ranComb__w_PhiV_Cut->GetXaxis()->SetRangeUser(0, 0.2);
		h_Mee__ranComb->DrawClone("");
		h_Mee__ranComb__w_PhiV_Cut->DrawClone("same");
		auto legend = new TLegend(0.53, 0.58, 0.78, 0.68);
		legend->SetFillColor(0); legend->SetBorderSize(0);
		legend->AddEntry(h_Mee__ranComb, "\t without #phi_{V} cut", "lp");
		legend->AddEntry(h_Mee__ranComb__w_PhiV_Cut, "\t with #phi_{V} cut", "lp");
		legend->SetMargin(0.20);
		gStyle->SetLegendTextSize(0.04);
		legend->DrawClone("same");

		c_temp->cd(3);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_Mee__ranComb->GetXaxis()->SetRangeUser(0, 1);
		h_Mee__ranComb->DrawClone("");

		c_temp->cd(4);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);

		c_temp->SaveAs(Form("roots/%d_PhiV_Check.png", number));
	}
}
