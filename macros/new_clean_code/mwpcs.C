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
//Assignment of the TWIM sections:
//0: Messel Down
//1: Messel Up
//2: Wixhausen Up
//3: Wixhausen Down

vector<vector<double> > mwpcs(TClonesArray* twim_tclone, TClonesArray* mw2_tclone, TClonesArray* mw3_tclone,TClonesArray* tof_tclone, TClonesArray* start_tclone){
	Int_t nr_twim = twim_tclone->GetEntries();
	Int_t nr_mw2 = mw2_tclone->GetEntries();
	Int_t nr_mw3 = mw3_tclone->GetEntries();
	Int_t nr_tof = tof_tclone->GetEntries();
	Int_t nr_start = start_tclone->GetEntries();
	vector<double> low_y {0,0,0,0};
	vector<double> high_y {0,0,0,0};
	vector<vector<double> > mwpc_vector;
	mwpc_vector.push_back(low_y);
	mwpc_vector.push_back(high_y);
	if (nr_twim == 2 && nr_mw2 == 2 && nr_mw3 == 2 && nr_tof == 2  && nr_start > 1){
		R3BSofTwimHitData** softwimhitdata  = new R3BSofTwimHitData*[nr_twim];
		R3BSofMwpcHitData** sofmwpc2hitdata = new R3BSofMwpcHitData*[nr_mw2];
		R3BSofMwpcHitData** sofmwpc3hitdata = new R3BSofMwpcHitData*[nr_mw3];
		R3BSofTofWHitData** softofwhitdata = new R3BSofTofWHitData*[nr_tof];
		R3BSofSciTcalData** sofscitcaldata = new R3BSofSciTcalData*[nr_start];
		vector<double> v_mw2_y;
		vector<double> v_mw3_y;
		vector<double> v_mw2_x;
		vector<double> v_mw3_x;
		vector<double> v_tof;
		vector<double> v_start_times;
		for (Int_t i = 0; i < 2; i++){
			softwimhitdata[i] = (R3BSofTwimHitData*)twim_tclone->At(i);
			sofmwpc2hitdata[i] = (R3BSofMwpcHitData*)mw2_tclone->At(i);
			sofmwpc3hitdata[i] = (R3BSofMwpcHitData*)mw3_tclone->At(i);
			softofwhitdata[i] = (R3BSofTofWHitData*)tof_tclone->At(i);
			sofscitcaldata[i] = (R3BSofSciTcalData*)start_tclone->At(i);
			v_mw2_y.push_back(sofmwpc2hitdata[i]->GetY());
			v_mw3_y.push_back(sofmwpc3hitdata[i]->GetY());
			v_mw2_x.push_back(sofmwpc2hitdata[i]->GetX());
			v_mw3_x.push_back(sofmwpc3hitdata[i]->GetX());
			v_tof.push_back(softofwhitdata[i]->GetX());
			v_start_times.push_back(sofscitcaldata[i]->GetRawTimeNs());
			}
		Double_t start_time = 0.5*(v_start_times[0]+v_start_times[1]);
		v_tof[0] =  v_tof[0] - start_time;
		v_tof[1] = v_tof[1] - start_time;
		Double_t charge_1 = softwimhitdata[0]->GetZcharge();
		Double_t charge_2 = softwimhitdata[1]->GetZcharge();
		Int_t fSection1 = softwimhitdata[0]->GetSecID();
		Int_t fSection2 = softwimhitdata[1]->GetSecID();
		vector<double> v_mw3_y_not_sorted = v_mw3_y;
		sort(v_mw2_y.begin(),v_mw2_y.end());
		sort(v_mw3_y.begin(),v_mw3_y.end());
		if ((charge_1+charge_2 < 100 && ((fSection1 == 1 && fSection2 == 3) || (fSection1 == 3 && fSection2 == 1))) ||       (charge_1+charge_2 < 100 && ((fSection1 == 0 && fSection2 == 2) || (fSection1 == 2 && fSection2 == 0)))){
		if ((fSection1 == 0 || fSection1 == 3) && (fSection2 == 1 || fSection2 == 2)){ //charge1 down, charge2 up
			if (v_mw3_y_not_sorted[0] == v_mw3_y[0]){
			 mwpc_vector[0][0] = charge_1;
			 mwpc_vector[0][1] = v_mw2_x[0];
			 mwpc_vector[0][2] = v_mw3_x[0];
			 mwpc_vector[1][0] = charge_2;
			 mwpc_vector[1][1] = v_mw2_x[1];
			 mwpc_vector[1][2] = v_mw3_x[1];
				}
			else{
			 mwpc_vector[0][0] = charge_1;
			 mwpc_vector[0][1] = v_mw2_x[1];
			 mwpc_vector[0][2] = v_mw3_x[1];
			 mwpc_vector[1][0] = charge_2;
			 mwpc_vector[1][1] = v_mw2_x[0];
			 mwpc_vector[1][2] = v_mw3_x[0];
				}
			if ((v_tof[0] > v_tof[1] && mwpc_vector[0][2] > mwpc_vector[1][2]) || (v_tof[0] < v_tof[1] && mwpc_vector[0][2] < mwpc_vector[1][2])){
				mwpc_vector[0][3] = (softofwhitdata[0]->GetTof());
				mwpc_vector[1][3] = (softofwhitdata[1]->GetTof());
				}
			else {
				mwpc_vector[0][3] = (softofwhitdata[1]->GetTof());
				mwpc_vector[1][3] = (softofwhitdata[0]->GetTof());
				}
			}
		else if ((fSection1 == 1 || fSection1 == 2) && (fSection2 == 0 || fSection2 == 3)){ //charge1 up, charge2 down
			if (v_mw3_y_not_sorted[0] == v_mw3_y[0]){
			 mwpc_vector[0][0] = charge_2;
			 mwpc_vector[0][1] = v_mw2_x[0];
			 mwpc_vector[0][2] = v_mw3_x[0];
			 mwpc_vector[1][0] = charge_1;
			 mwpc_vector[1][1] = v_mw2_x[1];
			 mwpc_vector[1][2] = v_mw3_x[1];
				}
			else {
			 mwpc_vector[0][0] = charge_2;
			 mwpc_vector[0][1] = v_mw2_x[1];
			 mwpc_vector[0][2] = v_mw3_x[1];
			 mwpc_vector[1][0] = charge_1;
			 mwpc_vector[1][1] = v_mw2_x[0];
			 mwpc_vector[1][2] = v_mw3_x[0];
				}

			if ((v_tof[0] > v_tof[1] && mwpc_vector[0][2] > mwpc_vector[1][2]) || (v_tof[0] < v_tof[1] && mwpc_vector[0][2] < mwpc_vector[1][2])){
				mwpc_vector[0][3] = (softofwhitdata[0]->GetTof());
				mwpc_vector[1][3] = (softofwhitdata[1]->GetTof());
				}
			else {
				mwpc_vector[0][3] = (softofwhitdata[1]->GetTof());
				mwpc_vector[1][3] = (softofwhitdata[0]->GetTof());
				}
			}
		}
		}
	return mwpc_vector;
	}
