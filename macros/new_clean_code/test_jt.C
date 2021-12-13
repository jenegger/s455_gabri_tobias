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
#include "analysis_histos.C"
#include "califa_p2p_selection.C"
#include "twim_selection.C"
#include "t_pattern.C"
#include "mwpcs.C"


#include <stdlib.h>
using namespace std;

void test_jt(const string& input_str){
my_histos();
//Input file
//char fname[500];
vector<double> v_false = {-1,-1,-1,-1};
vector<double> v_3_false = {-1,-1,-1};
//getting the mapped file.....
fstream fin;
fin.open("../parameters/crystal_mapping_s455_03.csv", ios::in); //FIXME:change to right .csv file

//new functionality, using map
map<int,int> map_id_califa_side;
string line, word, temp;
while(fin.peek()!=EOF) {
	getline(fin, line);
	stringstream s(line);
	vector<string> temp_vec;
	while (getline(s, word, ',')) {		temp_vec.push_back(word);
	}
	map_id_califa_side[stoi(temp_vec[10])] = stoi(temp_vec[0]);
	temp_vec.clear();
}
//end of getting mapped califa file....

Long64_t entries_califa = 0;
Long64_t entries_califa_hit = 0;
Long64_t entries_tof = 0;
Long64_t entries_twim = 0;
Long64_t entries_mw2 = 0;
Long64_t entries_mw3 = 0;
Int_t entries_wr = 0;

string fname = string("/u/land/tobias_jenegger/calibrated_s455_data/with_mw_reco/s455_03_273_") + input_str + string("_calibrated.root");
const char* char_fname= fname.c_str();

TChain* chain = new TChain("evt");
chain->Reset();
chain->Add(char_fname);
Long64_t nevents = chain->GetEntries();
cout << "total number of entries:\t" << nevents << endl;

//TClonesArray

TClonesArray* CalifaHitData = new TClonesArray("R3BCalifaHitData",3);
R3BCalifaHitData** califahitdata;
TBranch *branchCalifaHitData = chain->GetBranch("CalifaHitData");
branchCalifaHitData->SetAddress(&CalifaHitData);

R3BEventHeader* DataCA = new R3BEventHeader();
TBranch* branchData = chain->GetBranch("EventHeader.");
branchData->SetAddress(&DataCA);

TClonesArray* SofTwimHitData = new TClonesArray("R3BSofTwimHitData",2);
R3BSofTwimHitData** softwimhitdata;
TBranch *branchSofTwimHitData = chain->GetBranch("TwimHitData");
branchSofTwimHitData->SetAddress(&SofTwimHitData);

TClonesArray* SofToFWHitData = new TClonesArray("R3BSofTofWHitData",2);
R3BSofTofWHitData** softofwhitdata;
TBranch *branchSofToFWHitData = chain->GetBranch("TofWHitData");
branchSofToFWHitData->SetAddress(&SofToFWHitData);

TClonesArray* SofMwpc2HitData = new TClonesArray("R3BSofMwpcHitData",5);
R3BSofMwpcHitData** sofmwpc2hitdata;
TBranch *branchSofMwpc2HitData = chain->GetBranch("Mwpc2HitData");
branchSofMwpc2HitData->SetAddress(&SofMwpc2HitData);

TClonesArray* SofMwpc3HitData = new TClonesArray("R3BSofMwpcHitData",5);
R3BSofMwpcHitData** sofmwpc3hitdata;
TBranch *branchSofMwpc3HitData = chain->GetBranch("Mwpc3HitData");
branchSofMwpc3HitData->SetAddress(&SofMwpc3HitData);

TClonesArray* SofSciTcalData = new TClonesArray("R3BSofSciTcalData",2);
R3BSofSciTcalData** sofscitcaldata;
TBranch *branchSofSciTcalData = chain->GetBranch("SofSciTcalData");
branchSofSciTcalData->SetAddress(&SofSciTcalData);
// END of TClonesArray

const Double_t PI = 3.14159265358979323846;
uint64_t wr_Master_ts;
Long64_t events_califa_no_master = 0;
Long64_t events_califa_with_master = 0;
Long64_t events_total_califa = 0;
Long64_t events_with_cut = 0;
Long64_t events_one_proton = 0;
Long64_t events_two_proton = 0;
Long64_t events_50_41 = 0;
Long64_t events_50_42 = 0;

for(Long64_t i=0;i < nevents;i++){
    Long64_t evtnr = i;
    if (i%100000==0)
        cout<<"Processing event for charge analysis "<<i<<endl;
    chain->GetEvent(i);

    entries_twim = SofTwimHitData->GetEntries();
    entries_califa_hit = CalifaHitData->GetEntries();
    entries_tof = SofToFWHitData->GetEntries();
    entries_mw2 = SofMwpc2HitData->GetEntries();
    entries_mw3 = SofMwpc3HitData->GetEntries();
    //analysis part--------------------------
	int trigger_pattern = t_pattern(DataCA);
	uint64_t wr_master_ts = wrt_ts(DataCA);
//	if (wr_master_ts > 0){
//		int64_t t_diff_events = wr_master_ts - old_wrts;
//		}	

	vector<double> charge_check = twim_calibrated(SofTwimHitData);
	h1_trigger->Fill(trigger_pattern);
	if (trigger_pattern == 1){  // Sofstart
		h1_charge1_plus_charge2_tpat1->Fill(charge_check[0]+charge_check[1]);
		}
	if (trigger_pattern == 2){  // Fission (S+F)
		h1_charge1_plus_charge2_tpat2->Fill(charge_check[0]+charge_check[1]);
		}
	if (trigger_pattern == 3){  // CalifOR (S+C&&)
		h1_charge1_plus_charge2_tpat3->Fill(charge_check[0]+charge_check[1]);
		}
	if (trigger_pattern == 4){  // CalifAND (S+F+C&&)
		h1_charge1_plus_charge2_tpat4->Fill(charge_check[0]+charge_check[1]);
		}
	if (trigger_pattern == 7){ //!S
		h1_charge1_plus_charge2_tpat1->Fill(charge_check[0]+charge_check[1]);
		}
	if (trigger_pattern == 8){ //!S+F
		h1_charge1_plus_charge2_tpat2->Fill(charge_check[0]+charge_check[1]);
		}
	if (trigger_pattern == 9){ // !S + C&&
		h1_charge1_plus_charge2_tpat3->Fill(charge_check[0]+charge_check[1]);
		}
	if (trigger_pattern == 10){ // !S + F + C&&
		h1_charge1_plus_charge2_tpat4->Fill(charge_check[0]+charge_check[1]);
		}
	//Z_sum = 91 ---------------------------
	if (trigger_pattern >= 1 && entries_califa_hit >= 2 && entries_twim == 2 && entries_tof > 1 && (charge_check[0]+charge_check[1]) < 91.5 && (charge_check[0]+charge_check[1]) > 90.5 ){
		//this is just quick test, I will remove it soon
		//Int_t nr_mw3 = SofMwpc3HitData->GetEntries();
		//if (nr_mw3 == 2){
		//R3BSofMwpcHitData** sofmwpc3hitdata = new R3BSofMwpcHitData*[nr_mw3];
		//sofmwpc3hitdata[0] = (R3BSofMwpcHitData*)SofMwpc3HitData->At(0);
		//sofmwpc3hitdata[1] = (R3BSofMwpcHitData*)SofMwpc3HitData->At(1);
		//if (isnan(sofmwpc3hitdata[0]->GetY()) || isnan(sofmwpc3hitdata[1]->GetY())){
		//	cout << "This is eventnr. wit nan in mw3:\t" << evtnr << endl;
		//}
		//}
		//---------------------------------------------
			

		Califa_p2p_selection select_30_30(CalifaHitData,30,30,10);
		Califa_p2p_selection select_30_15(CalifaHitData,30,15,10);
		Califa_p2p_selection select_100_30(CalifaHitData,100,30,10);
		vector<double> v_sel_30_30 = select_30_30.best_comb_hits();
		int first_index_30_30 = v_sel_30_30[0];
		int second_index_30_30 = v_sel_30_30[1];
		vector<double> v_sel_30_15 = select_30_15.best_comb_hits();
		int first_index_30_15 = v_sel_30_15[0];
		int second_index_30_15 = v_sel_30_15[1];
		vector<double> v_sel_100_30 = select_100_30.best_comb_hits();
		int first_index_100_30 = v_sel_100_30[0];
		int second_index_100_30 = v_sel_100_30[1];
		//Find out multiplicity of protons (> 30MeV in lab system);
		int mult_proton_91 = select_30_30.multiplicity_energy_cut(100);
		h1_proton_mult_91->Fill(mult_proton_91);

		vector<vector<double> > array_v_sel_30_30 = select_30_30.get_califa_array();
		vector<vector<double> > array_v_sel_30_15 = select_30_15.get_califa_array();
		vector<vector<double> > array_v_sel_100_30 = select_100_30.get_califa_array();
		if (v_sel_30_30 != v_3_false){	
		if ((array_v_sel_30_30[first_index_30_30][0]+array_v_sel_30_30[second_index_30_30][0]) > 400 && (array_v_sel_30_30[first_index_30_30][0]+array_v_sel_30_30[second_index_30_30][0]) < 700){	
			vector<uint64_t>wrts_califa_vec = select_30_30.get_califa_time();
			int64_t time_diff_proton1 = wr_master_ts - wrts_califa_vec[first_index_30_30];
			int64_t time_diff_proton2 = wr_master_ts - wrts_califa_vec[second_index_30_30];
		h1_theta_sum_z_sum_905_915_30_30->Fill(array_v_sel_30_30[first_index_30_30][1]+array_v_sel_30_30[second_index_30_30][1]);
		h1_wr_diff_905_915_30_30->Fill(time_diff_proton1);
		h1_wr_diff_905_915_30_30->Fill(time_diff_proton2);
			}
			}
		if (v_sel_30_15 != v_3_false){

		//if ((array_v_sel_30_15[first_index_30_15][0]+array_v_sel_30_15[second_index_30_15][0]) > 400 && (array_v_sel_30_15[first_index_30_15][0]+array_v_sel_30_15[second_index_30_15][0]) < 700){
			vector<uint64_t>wrts_califa_vec = select_30_15.get_califa_time();
			int64_t time_diff_proton1 = wr_master_ts - wrts_califa_vec[first_index_30_15];
			int64_t time_diff_proton2 = wr_master_ts - wrts_califa_vec[second_index_30_15];
		h1_theta_sum_z_sum_905_915_30_15->Fill(array_v_sel_30_15[first_index_30_15][1]+array_v_sel_30_15[second_index_30_15][1]);
		h1_wr_diff_905_915_30_15->Fill(time_diff_proton1);
		h1_wr_diff_905_915_30_15->Fill(time_diff_proton2);
			if (time_diff_proton1 > -1200 && time_diff_proton1 < -200 && time_diff_proton2 > -1200 && time_diff_proton2 < -200){
				h1_theta_sum_z_sum_905_915_30_15_wr->Fill(array_v_sel_30_15[first_index_30_15][1]+array_v_sel_30_15[second_index_30_15][1]);
				}
			h2_theta_vs_phi_91->Fill(array_v_sel_30_15[second_index_30_15][1],array_v_sel_30_15[second_index_30_15][2]);
			h2_theta_vs_phi_91->Fill(array_v_sel_30_15[first_index_30_15][1],array_v_sel_30_15[first_index_30_15][2]);
			h2_thetasum_vs_tpat_91->Fill(array_v_sel_30_15[first_index_30_15][1]+array_v_sel_30_15[second_index_30_15][1],trigger_pattern);
			//}
			}
		if (v_sel_100_30 != v_3_false){
		if ((array_v_sel_100_30[first_index_100_30][0]+array_v_sel_100_30[second_index_100_30][0]) > 400 && (array_v_sel_100_30[first_index_100_30][0]+array_v_sel_100_30[second_index_100_30][0]) < 700){
			vector<uint64_t>wrts_califa_vec = select_100_30.get_califa_time();
			int64_t time_diff_proton1 = wr_master_ts - wrts_califa_vec[first_index_100_30];
			int64_t time_diff_proton2 = wr_master_ts - wrts_califa_vec[second_index_100_30];
		h1_theta_sum_z_sum_905_915_100_30->Fill(array_v_sel_100_30[first_index_100_30][1]+array_v_sel_100_30[second_index_100_30][1]);
		if ((array_v_sel_100_30[first_index_100_30][2] >50 && array_v_sel_100_30[first_index_100_30][2] < 120) || ( array_v_sel_100_30[second_index_100_30][2] > 50 && array_v_sel_100_30[second_index_100_30][2] < 120 )){

			//cout << "this is event with high_e hit in gamma range:\t" << evtnr << endl;
				}
		h1_wr_diff_905_915_100_30->Fill(time_diff_proton1);
		h1_wr_diff_905_915_100_30->Fill(time_diff_proton2);
			}
			}


		

		}
	//Z_sum = 92 --------------------------------
	if (trigger_pattern >= 1 && entries_califa_hit >= 2 && entries_twim == 2 && entries_tof > 1 && (charge_check[0]+charge_check[1]) < 92.5 && (charge_check[0]+charge_check[1]) > 91.5 ){
		Califa_p2p_selection select_30_30(CalifaHitData,30,30,10);
		Califa_p2p_selection select_30_15(CalifaHitData,30,15,10);
		Califa_p2p_selection select_100_30(CalifaHitData,100,30,10);
		vector<double> v_sel_30_30 = select_30_30.best_comb_hits();
		int first_index_30_30 = v_sel_30_30[0];
		int second_index_30_30 = v_sel_30_30[1];
		vector<double> v_sel_30_15 = select_30_15.best_comb_hits();
		int first_index_30_15 = v_sel_30_15[0];
		int second_index_30_15 = v_sel_30_15[1];
		vector<double> v_sel_100_30 = select_100_30.best_comb_hits();
		int first_index_100_30 = v_sel_100_30[0];
		int second_index_100_30 = v_sel_100_30[1];

		//Find out multiplicity of protons (> 30MeV in lab system);
		int mult_proton_92 = select_30_30.multiplicity_energy_cut(100);
		h1_proton_mult_92->Fill(mult_proton_92);

		vector<vector<double> > array_v_sel_30_30 = select_30_30.get_califa_array();
		vector<vector<double> > array_v_sel_30_15 = select_30_15.get_califa_array();
		vector<vector<double> > array_v_sel_100_30 = select_100_30.get_califa_array();
		if (v_sel_30_30 != v_3_false){	
		if ((array_v_sel_30_30[v_sel_30_30[0]][0]+array_v_sel_30_30[v_sel_30_30[1]][0]) > 400 && (array_v_sel_30_30[v_sel_30_30[0]][0]+array_v_sel_30_30[v_sel_30_30[1]][0]) < 700){	
			vector<uint64_t>wrts_califa_vec = select_30_30.get_califa_time();
			int64_t time_diff_proton1 = wr_master_ts - wrts_califa_vec[first_index_30_30];
			int64_t time_diff_proton2 = wr_master_ts - wrts_califa_vec[second_index_30_30];
		h1_theta_sum_z_sum_915_925_30_30->Fill(array_v_sel_30_30[v_sel_30_30[0]][1]+array_v_sel_30_30[v_sel_30_30[1]][1]);
		h1_wr_diff_915_925_30_30->Fill(time_diff_proton1);
		h1_wr_diff_915_925_30_30->Fill(time_diff_proton2);
			}
			}
		if (v_sel_30_15 != v_3_false){
		if ((array_v_sel_30_15[v_sel_30_15[0]][0]+array_v_sel_30_15[v_sel_30_15[1]][0]) > 400 && (array_v_sel_30_15[v_sel_30_15[0]][0]+array_v_sel_30_15[v_sel_30_15[1]][0]) < 700){
			vector<uint64_t>wrts_califa_vec = select_30_15.get_califa_time();
			int64_t time_diff_proton1 = wr_master_ts - wrts_califa_vec[first_index_30_15];
			int64_t time_diff_proton2 = wr_master_ts - wrts_califa_vec[second_index_30_15];
		h1_theta_sum_z_sum_915_925_30_15->Fill(array_v_sel_30_15[v_sel_30_15[0]][1]+array_v_sel_30_15[v_sel_30_15[1]][1]);
		h1_wr_diff_915_925_30_15->Fill(time_diff_proton1);
		h1_wr_diff_915_925_30_15->Fill(time_diff_proton2);
			}

		h2_theta_vs_phi_92->Fill(array_v_sel_30_15[v_sel_30_15[0]][1],array_v_sel_30_15[v_sel_30_15[0]][2]);
		h2_theta_vs_phi_92->Fill(array_v_sel_30_15[v_sel_30_15[1]][1],array_v_sel_30_15[v_sel_30_15[1]][2]);
		h2_thetasum_vs_tpat_92->Fill(array_v_sel_30_15[v_sel_30_15[0]][1]+array_v_sel_30_15[v_sel_30_15[1]][1],trigger_pattern);
			}
		if (v_sel_100_30 != v_3_false){
		if ((array_v_sel_100_30[v_sel_100_30[0]][0]+array_v_sel_100_30[v_sel_100_30[1]][0]) > 400 && (array_v_sel_100_30[v_sel_100_30[0]][0]+array_v_sel_100_30[v_sel_100_30[1]][0]) < 700){
			vector<uint64_t>wrts_califa_vec = select_100_30.get_califa_time();
			int64_t time_diff_proton1 = wr_master_ts - wrts_califa_vec[first_index_100_30];
			int64_t time_diff_proton2 = wr_master_ts - wrts_califa_vec[second_index_100_30];
		h1_theta_sum_z_sum_915_925_100_30->Fill(array_v_sel_100_30[v_sel_100_30[0]][1]+array_v_sel_100_30[v_sel_100_30[1]][1]);
		h1_wr_diff_915_925_100_30->Fill(time_diff_proton1);
		h1_wr_diff_915_925_100_30->Fill(time_diff_proton2);
			if (time_diff_proton1 > -1200 && time_diff_proton1 < -200 && time_diff_proton2 > -1200 && time_diff_proton2 < -200){
				h1_theta_sum_z_sum_915_925_100_30_wr->Fill(array_v_sel_100_30[v_sel_100_30[0]][1]+array_v_sel_100_30[v_sel_100_30[1]][1]);
				}
			}
			}


		

		}
		
	//Z_sum = 90---------------------
	if (trigger_pattern >= 1 && entries_califa_hit >= 2 && entries_twim == 2 && entries_tof > 1 && (charge_check[0]+charge_check[1]) < 90.5 && (charge_check[0]+charge_check[1]) > 89.5 ){
		Califa_p2p_selection select_30_30(CalifaHitData,30,30,10);
		//Find out multiplicity of protons (> 30MeV in lab system);
                int mult_proton_90 = select_30_30.multiplicity_energy_cut(100);
		h1_proton_mult_90->Fill(mult_proton_90);
		}

	if(trigger_pattern >= 1 && entries_califa_hit >= 1){
		Califa_p2p_selection select_100_30(CalifaHitData,100,30,10);

		//check now how many protons do we see ... (meaning how many particles with E > 100 MeV
		if ( select_100_30.multiplicity_energy_cut(100) ==1){ // one proton like particle in CALIFA
			h2_charge1_charge_2_diff_vs_z_sum_one_CALIFA->Fill(charge_check[0]-charge_check[1],charge_check[0]+charge_check[1]);
			}
		if ( select_100_30.multiplicity_energy_cut(100) ==2){ // two proton like particles in CALIFA
			h2_charge1_charge_2_diff_vs_z_sum_two_CALIFA->Fill(charge_check[0]-charge_check[1],charge_check[0]+charge_check[1]);
			}
		if ( select_100_30.multiplicity_energy_cut(100) ==3){ // three proton like particles in CALIFA
			h2_charge1_charge_2_diff_vs_z_sum_three_CALIFA->Fill(charge_check[0]-charge_check[1],charge_check[0]+charge_check[1]);
			}
		}
	//part to cut on Z_50 and Z_41 and look at the CALIFA spectrum...
	if ((charge_check[0] <41.5 && charge_check[0]>40.5 && charge_check[1]<50.5 && charge_check[1]>49.5) || 
	    (charge_check[1] <41.5 && charge_check[1]>40.5 && charge_check[0]<50.5 && charge_check[0]>49.5))  {
		events_50_41++;
		Double_t tin_nr = 0;
		vector<vector<double> > v_mwpc = mwpcs(SofTwimHitData,SofMwpc2HitData,SofMwpc3HitData,SofToFWHitData,SofSciTcalData);
		if (v_mwpc[0][0] <50.5 && v_mwpc[0][0] > 49.5){
			tin_nr = 1;
			cout << "this is interesting event:\t" << evtnr << endl;
			//cout << "this is the x component in MW2, particle1:\t" << v_mwpc[0][1] << endl;
			//cout << "this is the x component in MW2, particle2:\t" << v_mwpc[1][1] << endl;
			h2_mw2_x_vs_mw3_x_sn->Fill(v_mwpc[0][1],v_mwpc[0][3]);	
			h2_mw23_diffx_vs_tof_sn->Fill(v_mwpc[0][3]-v_mwpc[0][1],v_mwpc[0][5]);
			h2_mw2_x_vs_mw3_x_41->Fill(v_mwpc[1][1],v_mwpc[1][3]);
			if (v_mwpc[0][1] > 12 && v_mwpc[0][3] < 7){
				cout << "critical eventnr:\t" << evtnr << endl;
				}
			}
		else if (v_mwpc[1][0] <50.5 && v_mwpc[1][0] > 49.5){
			tin_nr = 2;
			h2_mw2_x_vs_mw3_x_sn->Fill(v_mwpc[1][1],v_mwpc[1][3]);
			h2_mw23_diffx_vs_tof_sn->Fill(v_mwpc[1][3]-v_mwpc[1][1],v_mwpc[1][5]);
			h2_mw2_x_vs_mw3_x_41->Fill(v_mwpc[0][1],v_mwpc[0][3]);
			if (v_mwpc[1][1] > 12 && v_mwpc[1][3] < 7){
				cout << "critical eventnr:\t" << evtnr << endl;
				}
			}
		if (entries_califa_hit >= 1){
			Califa_p2p_selection select_100_30(CalifaHitData,100,30,10); //the parameter values are here arbitrary...
			vector<vector<double> > califa_array_doppler = select_100_30.get_califa_array();
			vector<double> califa_array_high_gamma = select_100_30.get_highest_califa_gamma(20);  //vector with highest energy gamma
			vector<uint64_t>califa_times_vec = select_100_30.get_califa_time();
			h1_mult_with_M->Fill(entries_califa_hit);
			if (califa_array_high_gamma != v_false){
			h1_califa_energy_z50_z41_high_gamma->Fill(califa_array_high_gamma[3]);
			h1_califa_energy_z50_z41_no_doppler_high_gamma->Fill(califa_array_high_gamma[0]);
			h2_califa_energy_z50_z41_no_doppler_high_gamma_vs_theta->Fill(califa_array_high_gamma[0],califa_array_high_gamma[1]);
			h2_califa_energy_z50_z41_no_doppler_high_gamma_vs_phi->Fill(califa_array_high_gamma[0],califa_array_high_gamma[2]);
			}
			Int_t multiplicity_gammas = 0;
			for (int n = 0; n< califa_array_doppler.size();n++){
				if (califa_array_doppler[n][0] < 20){
					//fill histo for energies in low regime...
					uint64_t time_califa = califa_times_vec[n];
					int64_t time_diff_wrm_wrc= wr_master_ts - time_califa;
					h1_wr_diff_z_50_z41->Fill(time_diff_wrm_wrc);
					if (califa_array_doppler[n][2] > -90 && califa_array_doppler[n][2] < 90){
						h1_wr_diff_z_50_z41_messel->Fill(time_diff_wrm_wrc);
						}
					if (califa_array_doppler[n][2] < -90 || califa_array_doppler[n][2] > 90){
						h1_wr_diff_z_50_z41_wix->Fill(time_diff_wrm_wrc);
						}

					h1_califa_energy_z50_z41->Fill(califa_array_doppler[n][3]);
					h1_califa_energy_z50_z41_no_doppler->Fill(califa_array_doppler[n][0]);
					if(select_100_30.multiplicity_energy_cut(200) >=2){
						h1_califa_energy_z50_z41_200->Fill(califa_array_doppler[n][3]);
						}
					multiplicity_gammas++;
					}
				}
			//now select only events with multiplicity < 4 && -1200 < wrm -wr_califa < -200 ns
			if (multiplicity_gammas < 4){
				vector<double> califa_low_energies_wr_constraint;
				for (int n = 0; n< califa_array_doppler.size();n++){
					if (califa_array_doppler[n][0] < 20){
						uint64_t time_califa = califa_times_vec[n];
						int64_t wr_time_diff = wr_master_ts - time_califa;
						if (wr_time_diff > -1200 && wr_time_diff < -200){
							h1_califa_energy_z50_z41_mult_wr->Fill(califa_array_doppler[n][3]);
							h1_califa_energy_z50_addup_mult_wr->Fill(califa_array_doppler[n][3]);
							califa_low_energies_wr_constraint.push_back(califa_array_doppler[n][3]);
							//check regions when plotting mw2 mw3/tof:
							//down right
							if (v_mwpc[0][0] <50.5 && v_mwpc[0][0] > 49.5 && v_mwpc[0][1] > 12 && v_mwpc[0][3] < 7){
								h1_califa_energy_z50_z41_mult_wr_rd->Fill(califa_array_doppler[n][3]);
								}
							else if (v_mwpc[1][0] <50.5 && v_mwpc[1][0] > 49.5 && v_mwpc[1][1] > 12 && v_mwpc[1][3] < 7){
								h1_califa_energy_z50_z41_mult_wr_rd->Fill(califa_array_doppler[n][3]);
								}

							//down left
							if (v_mwpc[0][0] <50.5 && v_mwpc[0][0] > 49.5 && v_mwpc[0][1] < 3 && v_mwpc[0][3] < 7){
								h1_califa_energy_z50_z41_mult_wr_ld->Fill(califa_array_doppler[n][3]);	
								}
							else if (v_mwpc[1][0] <50.5 && v_mwpc[1][0] > 49.5 && v_mwpc[1][1] < 3 && v_mwpc[1][3] < 7){
								h1_califa_energy_z50_z41_mult_wr_ld->Fill(califa_array_doppler[n][3]);
								}
							//up right
							if (v_mwpc[0][0] <50.5 && v_mwpc[0][0] > 49.5 && v_mwpc[0][1] > 12 && v_mwpc[0][3] < 50){
								h1_califa_energy_z50_z41_mult_wr_ru->Fill(califa_array_doppler[n][3]);
								}
							else if (v_mwpc[1][0] <50.5 && v_mwpc[1][0] > 49.5 && v_mwpc[1][1] > 12 && v_mwpc[1][3] < 50){
								h1_califa_energy_z50_z41_mult_wr_ru->Fill(califa_array_doppler[n][3]);
								}
							//up left
							if (v_mwpc[0][0] <50.5 && v_mwpc[0][0] > 49.5 && v_mwpc[0][1] < 3 && v_mwpc[0][3] > 50){
								h1_califa_energy_z50_z41_mult_wr_lu->Fill(califa_array_doppler[n][3]);
								}
							else if (v_mwpc[1][0] <50.5 && v_mwpc[1][0] > 49.5 && v_mwpc[1][1] < 3 && v_mwpc[1][3] < 50){
								h1_califa_energy_z50_z41_mult_wr_lu->Fill(califa_array_doppler[n][3]);
								}
							//this here is a check if I see different isotopes of thin
							if(tin_nr == 1){
							if ((v_mwpc[0][3]-v_mwpc[0][1])*(-0.0095976)+37.5814 < v_mwpc[0][5]){
								//fill histogram for values above...
								h1_califa_energy_tin_above->Fill(califa_array_doppler[n][3]);
								}
							else {
								//fill histogram for values below...
								h1_califa_energy_tin_below->Fill(califa_array_doppler[n][3]);
								}

							}
							else if (tin_nr == 2){

							if ((v_mwpc[1][3]-v_mwpc[1][1])*(-0.0095976)+37.5814 < v_mwpc[1][5]){
								//fill histogram for values above...
								h1_califa_energy_tin_above->Fill(califa_array_doppler[n][3]);
								}
							else {
								//fill histogram for values below...
								h1_califa_energy_tin_above->Fill(califa_array_doppler[n][3]);
								}
							}
							//--------------------------------------------------------
							}
						}
					}
				if (califa_low_energies_wr_constraint.size()){
					Double_t max_gamma = *max_element(califa_low_energies_wr_constraint.begin(),califa_low_energies_wr_constraint.end());
					h1_califa_energy_z50_z41_mult_wr_high_e->Fill(max_gamma);
					}
				}
			//end of selection multiplicity && wr difference-----------------------------------
			h1_multiplicity_z50_z41->Fill(multiplicity_gammas);
			}
		}
	//part to cut on Z_50 and Z_42 and look at the CALIFA spectrum
	if ((charge_check[0] <42.5 && charge_check[0]>41.5 && charge_check[1]<50.5 && charge_check[1]>49.5) || 
	    (charge_check[1] <42.5 && charge_check[1]>41.5 && charge_check[0]<50.5 && charge_check[0]>49.5))  {
		events_50_42++;
		Double_t tin_nr = 0;
		vector<vector<double> > v_mwpc = mwpcs(SofTwimHitData,SofMwpc2HitData,SofMwpc3HitData,SofToFWHitData,SofSciTcalData);
		if (v_mwpc[0][0] <50.5 && v_mwpc[0][0] > 49.5){
			tin_nr = 1;
			h2_mw2_x_vs_mw3_x_sn->Fill(v_mwpc[0][1],v_mwpc[0][3]);	
			h2_mw23_diffx_vs_tof_sn->Fill(v_mwpc[0][3]-v_mwpc[0][1],v_mwpc[0][5]);
			h2_mw2_x_vs_mw3_x_42->Fill(v_mwpc[1][1],v_mwpc[1][3]);
			if (v_mwpc[0][1] > 12 && v_mwpc[0][3] < 7){
				cout << "critical eventnr:\t" << evtnr << endl;
				}
			}
		else if (v_mwpc[1][0] <50.5 && v_mwpc[1][0] > 49.5){
			tin_nr = 2;
			h2_mw2_x_vs_mw3_x_sn->Fill(v_mwpc[1][1],v_mwpc[1][3]);
			h2_mw23_diffx_vs_tof_sn->Fill(v_mwpc[1][3]-v_mwpc[1][1],v_mwpc[1][5]);
			h2_mw2_x_vs_mw3_x_42->Fill(v_mwpc[0][1],v_mwpc[0][3]);
			if (v_mwpc[1][1] > 12 && v_mwpc[1][3] < 7){
				cout << "critical eventnr:\t" << evtnr << endl;
				}
			}
		
		if (entries_califa_hit >= 1){
			Califa_p2p_selection select_100_30(CalifaHitData,100,30,10); //the parameter values are here arbitrary...
			vector<vector<double> > califa_array_doppler = select_100_30.get_califa_array();
			vector<double> califa_array_high_gamma = select_100_30.get_highest_califa_gamma(20); //vector with highest energy gamma
			vector<uint64_t>califa_times_vec = select_100_30.get_califa_time();
			if (califa_array_high_gamma != v_false){
			h1_califa_energy_z50_z42_high_gamma->Fill(califa_array_high_gamma[3]);
			h1_califa_energy_z50_z42_no_doppler_high_gamma->Fill(califa_array_high_gamma[0]);
			h2_califa_energy_z50_z42_no_doppler_high_gamma_vs_theta->Fill(califa_array_high_gamma[0],califa_array_high_gamma[1]);
			h2_califa_energy_z50_z42_no_doppler_high_gamma_vs_phi->Fill(califa_array_high_gamma[0],califa_array_high_gamma[2]);
			}
			Int_t multiplicity_gammas = 0;
			for (int n = 0; n< califa_array_doppler.size();n++){
				if (califa_array_doppler[n][0] < 20){
					//fill histo for energies in low regime...
					uint64_t time_califa = califa_times_vec[n];
					int64_t time_diff_wrm_wrc= wr_master_ts - time_califa;
					h1_wr_diff_z_50_z42->Fill(time_diff_wrm_wrc);
					h1_califa_energy_z50_z42->Fill(califa_array_doppler[n][3]);
					h1_califa_energy_z50_z42_no_doppler->Fill(califa_array_doppler[n][0]);
					multiplicity_gammas++;
					}
				}
			//now select only events with multiplicity <= 5 && -1200 < wrm -wr_califa < -200 ns
			if (multiplicity_gammas < 4){
				vector<double> califa_low_energies_wr_constraint;
				for (int n = 0; n< califa_array_doppler.size();n++){
					if (califa_array_doppler[n][0] < 20){
						uint64_t time_califa = califa_times_vec[n];
						int64_t wr_time_diff = wr_master_ts - time_califa;
						if (wr_time_diff > -1200 && wr_time_diff < -200){
							h1_califa_energy_z50_z42_mult_wr->Fill(califa_array_doppler[n][3]);
							h1_califa_energy_z50_addup_mult_wr->Fill(califa_array_doppler[n][3]);
							califa_low_energies_wr_constraint.push_back(califa_array_doppler[n][3]);
							//this here is a check if I see different isotopes of thin
							if(tin_nr == 1){
							if ((v_mwpc[0][3]-v_mwpc[0][1])*(-0.0095976)+37.5814 < v_mwpc[0][5]){
								//fill histogram for values above...
								h1_califa_energy_tin_above->Fill(califa_array_doppler[n][3]);
								}
							else {
								//fill histogram for values below...
								h1_califa_energy_tin_below->Fill(califa_array_doppler[n][3]);
								}

							}
							else if (tin_nr == 2){

							if ((v_mwpc[1][3]-v_mwpc[1][1])*(-0.0095976)+37.5814 < v_mwpc[1][5]){
								//fill histogram for values above...
								h1_califa_energy_tin_above->Fill(califa_array_doppler[n][3]);
								}
							else {
								//fill histogram for values below...
								h1_califa_energy_tin_above->Fill(califa_array_doppler[n][3]);
								}
							}
							//--------------------------------------------------------
							}
						}
					}
				if (califa_low_energies_wr_constraint.size()){
					Double_t max_gamma = *max_element(califa_low_energies_wr_constraint.begin(),califa_low_energies_wr_constraint.end());
					h1_califa_energy_z50_z42_mult_wr_high_e->Fill(max_gamma);
					}
				}
			//end of selection multiplicity && wr difference-----------------------------------
			h1_multiplicity_z50_z42->Fill(multiplicity_gammas);

			}
		}

    //---------------------------------------

    if (entries_twim == 2){
    int trigger_pattern = t_pattern(DataCA);
    if (trigger_pattern >= 1){
    
    vector<double> charge_vector = twim_calibrated(SofTwimHitData);
    if (charge_vector[0] > 0 && charge_vector[1] > 0){
    h2_charge1_charge_2->Fill(charge_vector[0],charge_vector[1]);
    h2_charge1_charge_2_diff_vs_z_sum->Fill(charge_vector[0]-charge_vector[1],charge_vector[0]+charge_vector[1]);
    }
    }
    }
}
cout << "Events with fission products Z=50 && Z=41:\t" << events_50_41 << endl;
cout << "Events with fission products Z=50 && Z=42:\t" << events_50_42 << endl;

//char f_out_name[500];
string f_out_name = string("/u/land/tobias_jenegger/analysis_s455_data/s455_03_273_") + input_str + string("_analysis.root");
//sprintf(f_out_name,"../../file_output/twim_analysis_main0273_10_test.root");
write_to_file(f_out_name);
}
