//insert all headers
#include "TTree.h"
#include "TFile.h"
#include <math.h>
#include <string>
#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TClonesArray.h"
#include "TChain.h"

vector<double> twim_selection(TClonesArray* twim_tclone, TClonesArray* tof_tclone){
	Double_t charge_1 = 0.;
	Double_t charge_2 = 0.;
	vector<double> charge_vec = {charge_1,charge_2};
	Long64_t nr_twim = twim_tclone->GetEntries();
	Long64_t nr_tof = tof_tclone->GetEntries();
	if (nr_twim == 2 && nr_tof > 1){
	R3BSofTwimHitData** softwimhitdata  = new R3BSofTwimHitData*[nr_twim];
	R3BSofTofWHitData** softofwhitdata = new R3BSofTofWHitData*[nr_tof];

	Double_t tof[2] = { 0., 0. };
	int padid[2] = { 0, 0 };
	for (Int_t ihit = 0; ihit < nr_tof; ihit++){
	
			softofwhitdata[ihit] = (R3BSofTofWHitData*)tof_tclone->At(ihit);
			if (padid[0] == 0){
				padid[0] = softofwhitdata[ihit]->GetPaddle();
				tof[0] = softofwhitdata[ihit]->GetTof();
				}
			else {
				if (softofwhitdata[ihit]->GetPaddle() > padid[0] && softofwhitdata[ihit]->GetPaddle() - padid[0] > 1){
					padid[1] = softofwhitdata[ihit]->GetPaddle();
					tof[1] = softofwhitdata[ihit]->GetTof(); // right
					}	
				else if (softofwhitdata[ihit]->GetPaddle() - padid[0] < -1){
					tof[1] = tof[0]; // right
					padid[1] = padid[0];
					tof[0] = softofwhitdata[ihit]->GetTof(); // new left
					padid[0] = softofwhitdata[ihit]->GetPaddle();
					}
				}
			}

		//insert logic from TWIM online...
		Float_t zr = 0., zl = 0.;
		for (Int_t ihit = 0.; ihit < nr_twim; ihit++){
			softwimhitdata[ihit] = (R3BSofTwimHitData*)twim_tclone->At(ihit);
			if (softwimhitdata[ihit]->GetSecID() == 0){
				zl = softwimhitdata[ihit]->GetZcharge();
				}
			else if (softwimhitdata[ihit]->GetSecID() == 2){
				zr = softwimhitdata[ihit]->GetZcharge();
				}
			}
		//end of logic from TWIM online...
		if (zr > 0. && zl > 0. && tof[0] > 0. && tof[1] > 0.){	
			Double_t charge_1 = zl+52.55-(147.39-6.568844*tof[0]+0.1109114*tof[0]*tof[0])-8.82;
			Double_t charge_2 = zr+52.55-(147.39-6.568844*tof[1]+0.1109114*tof[1]*tof[1])-8.35-0.9;
			charge_vec[0] = charge_1;
			charge_vec[1] = charge_2;

			}
		delete [] softwimhitdata;
		delete [] softofwhitdata;
		}
		return charge_vec;

}

