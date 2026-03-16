#ifndef StAnaCuts_H
#define StAnaCuts_H

/* **************************************************
 *
 *  Authors: Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */
 /* ***********************************************************************
  * Extract pure electron to recalibrating nsigama_e vs pT,eta,phi vs Cen *
  * ***********************************************************************
  */
#include "Rtypes.h"
#include <string>
#include <array>

namespace anaCuts
{
	const std::array<UInt_t, 4> trigNumber = {860001, 860002, 860011, 860012};// 200GeV RFF+FF OO 2021 MBinimumBias

	// event cuts 
	Float_t const Vz_up = 30;// < cm.
	Float_t const Vz_low = -30;// < cm.
	Float_t const Verror = 1.0e-5; // 
	Float_t const Vr = 2; // cm
	Float_t const vzVpdVz = 3; // cm
	// tracks cuts
	Int_t   const NHitsFit_highPt = 40;//40
	Int_t   const NHitsDedx_highPt = 30;//30
	Int_t   const NHitsFit_lowPt = 20;//13
	Int_t   const NHitsDedx_lowPt = 14;//10
	Float_t const NHitsFitRatio = 0.52;
	Float_t const Dca_highPt = 1;//1
	Float_t const Dca_lowPt = 3;//1
	Float_t const Pt = 0.2;
	Float_t const Eta = 0.9;
    Float_t const PhiVCutMRange = 0.2;
}
#endif
