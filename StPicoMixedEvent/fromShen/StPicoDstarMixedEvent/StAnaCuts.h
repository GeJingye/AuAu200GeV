#ifndef StAnaCuts_H
#define StAnaCuts_H

/* **************************************************
 *
 *  Authors: Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */

#include "Rtypes.h"
#include <string>
#include <array>

namespace anaCuts
{
   std::array<unsigned int, 1> const triggers = {
   700001
   };    
   //cut before QA
   float const qavz = 40.0;// < cm.
   float const qaVerror = 1.0e-5; //
   float const qaVr = 2.0; //cm
   float const qavzVpdVz = 3; //cm

   // QA tracks cuts
   float const qaGPt = 0.20;
   int const qaNHitsFit = 20;
   int const qaNHitsDedx = 12;
   float const qaDca = 3;// < cm
   float const qaEta = 1; 
   float const qaTofPion=4;
   float const qaTpcPion=4;

   //cut 
   float const vz = 60.0;// < cm.
   float const vz_up = 60.0;// < cm.
   float const vz_low = -60.0;// < cm.
   float const Verror = 1.0e-5; //
   float const Vr = 2.0; //cm
   float const vzVpdVz = 3; //cm
   //float const vzVpdVz = 1e3; //cm
   //float EEdcaDaughter = 1; //cm   20200720 cut
   float EEdcaDaughter = 1; //cm   
   // QA tracks cuts
   // float const GPt = 0.20;
   float const GPt = 0.15;
   float const pPt = 0.2;
   int const NHitsFit = 15;
   int const NHitsDedx = 10;
   int const NHitsFit2Poss = 0.52;
   float const Dca = 3;// < cm
   //float const Dca = 1.0;// < cm
   float const Eta = 1.0; 
   float const TofPion=4;
   float const TpcPion=4;
   
   int const nparVz_mult = 7;
   float parVz_mult[nparVz_mult]={435.9,-0.02413,-0.003707,0.0002204,1.487e-5,-2.95e-07,-1.866e-8};
   int const nCent = 9 ;
   float Refmult_cent[nCent] = {7,16,31,54,89,138,205,299,361}; //refmult > par[i],  70-80%, 60-70%, ... ,0-5%
}
#endif
