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


bool sorty( const vector<double>& v1,const vector<double>& v2 ) {
                                    return v1[1] > v2[1];
                                    }


vector<vector<double> > mwpcs(TClonesArray* twim_tclone, TClonesArray* mw2_tclone, TClonesArray* mw3_tclone,TClonesArray* tof_tclone, TClonesArray* start_tclone){
	Int_t nr_twim = twim_tclone->GetEntries();
	Int_t nr_mw2 = mw2_tclone->GetEntries();
	Int_t nr_mw3 = mw3_tclone->GetEntries();
	Int_t nr_tof = tof_tclone->GetEntries();
	Int_t nr_start = start_tclone->GetEntries();
	vector<double> low_y {0,0,0,0,0,0};
	vector<double> high_y {0,0,0,0,0,0};
	vector<vector<double> > mwpc_vector;
	mwpc_vector.push_back(low_y);
	mwpc_vector.push_back(high_y);
	if (nr_twim == 2 && nr_mw2 == 2 && nr_mw3 == 2 && nr_tof == 2  && nr_start > 1){
		R3BSofTwimHitData** softwimhitdata  = new R3BSofTwimHitData*[nr_twim];
		R3BSofMwpcHitData** sofmwpc2hitdata = new R3BSofMwpcHitData*[nr_mw2];
		R3BSofMwpcHitData** sofmwpc3hitdata = new R3BSofMwpcHitData*[nr_mw3];
		R3BSofTofWHitData** softofwhitdata = new R3BSofTofWHitData*[nr_tof];
		R3BSofSciTcalData** sofscitcaldata = new R3BSofSciTcalData*[nr_start];
		vector<vector<double> > v_mw2;
		vector<vector<double> > v_tofw;
		vector<double> v_mw3_x;
		vector<double> v_mw3_y;


		vector<double> v_tof;
		vector<double> v_start_times;
		for (Int_t i = 0; i < 2; i++){
			softwimhitdata[i] = (R3BSofTwimHitData*)twim_tclone->At(i);
			sofmwpc2hitdata[i] = (R3BSofMwpcHitData*)mw2_tclone->At(i);
			sofmwpc3hitdata[i] = (R3BSofMwpcHitData*)mw3_tclone->At(i);
			softofwhitdata[i] = (R3BSofTofWHitData*)tof_tclone->At(i);
			sofscitcaldata[i] = (R3BSofSciTcalData*)start_tclone->At(i);
			
			vector<double> v_mw2_temp;
			vector<double> v_tofw_temp;
			v_mw2_temp.push_back(sofmwpc2hitdata[i]->GetX());
			v_mw2_temp.push_back(sofmwpc2hitdata[i]->GetY());
			v_tofw_temp.push_back(softofwhitdata[i]->GetX());
			v_tofw_temp.push_back(softofwhitdata[i]->GetY());
			v_tofw_temp.push_back(softofwhitdata[i]->GetTof());
			v_mw2.push_back(v_mw2_temp);
			v_tofw.push_back(v_tofw_temp);
			v_mw2_temp.clear();
			v_tofw_temp.clear();
			v_mw3_x.push_back(sofmwpc3hitdata[i]->GetX());
			v_mw3_y.push_back(sofmwpc3hitdata[i]->GetY());
			//ok what this is doing I don't know....
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
		//TODO: sort v_mw2 and v_tof according to the y values
		sort(v_mw2.begin(),v_mw2.end(),sorty);
        sort(v_tofw.begin(),v_tofw.end(),sorty);	
		//TODO: sort v_mx3_x according the entries
		sort(v_mw3_x.begin(),v_mw3_x.end());
		//TODO: sort v_mw3_y according the entries		
		sort(v_mw3_y.begin(),v_mw3_y.end());

		if (!isnan(v_mw2[0][1]) && !isnan(v_mw2[1][1]) && !isnan(v_mw3_y[0]) && !isnan(v_mw3_y[1]) && v_mw2[0][1] != -1000 && v_mw2[1][1] != -1000 && v_mw3_y[0] != -1000 && v_mw3_y[1] != -1000 && !isnan(v_mw2[0][0]) && !isnan(v_mw2[1][0]) && !isnan(v_mw3_x[0]) && !isnan(v_mw3_x[1]) && v_mw2[0][0] != -1000 &&   v_mw2[1][0] != -1000 && v_mw3_x[0] != -1000 && v_mw3_x[1] != -1000){
		if ((charge_1+charge_2 < 100 && ((fSection1 == 1 && fSection2 == 3) || (fSection1 == 3 && fSection2 == 1))) ||       (charge_1+charge_2 < 100 && ((fSection1 == 0 && fSection2 == 2) || (fSection1 == 2 && fSection2 == 0)))){

		mwpc_vector[0][1] = v_mw2[0][0]; //part1, x val, mw2
		mwpc_vector[0][2] = v_mw2[0][1]; //part1, y val, mw2
		mwpc_vector[1][1] = v_mw2[1][0]; //part2, x val, mw2
		mwpc_vector[1][2] = v_mw2[1][1]; //part2, y val, mw2
        mwpc_vector[0][5] = v_tofw[0][2];//part1, time TOF
        mwpc_vector[1][5] = v_tofw[1][2];//part2, time TOF
		if (v_tofw[0][0] < v_tofw[1][0]){
			mwpc_vector[0][3] = v_mw3_x[0]; //part1, x val, mw3
			mwpc_vector[0][4] = v_mw3_y[0]; //part1, y val, mw3
			mwpc_vector[1][3] = v_mw3_x[1]; //part2, x val, mw3
			mwpc_vector[1][4] = v_mw3_y[1]; //part2, y val, mw3
			}
		else if (v_tofw[1][0] < v_tofw[0][0]){
			mwpc_vector[0][3] = v_mw3_x[1];
			mwpc_vector[0][4] = v_mw3_y[0];
			mwpc_vector[1][3] = v_mw3_x[0];
			mwpc_vector[1][4] = v_mw3_y[1];
			}
		if ((fSection1 == 0 || fSection1 == 3) && (fSection2 == 1 || fSection2 == 2)){ //charge1 down, charge2 up
			mwpc_vector[0][0] = charge_1;
			mwpc_vector[1][0] = charge_2;
			}
		else if ((fSection1 == 1 || fSection1 == 2) && (fSection2 == 0 || fSection2 == 3)){ //charge1 up, charge2 down
			mwpc_vector[0][0] = charge_2;
			mwpc_vector[1][0] = charge_1;
			}
			}
		}
	}
	return mwpc_vector;
	}
