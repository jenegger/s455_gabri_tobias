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


#include <stdlib.h>
using namespace std;

void test_jt(){
my_histos();
//califa_p2p_selection();
//Input file
char fname[500];
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
Int_t entries_wr = 0;

//sprintf(fname,"/u/land/antiagg/s455/sofia/twim/calibration/4th_total_twim_cal/273_TWIM_corrbytof_hitdata.root");
sprintf(fname,"../../file_src/s455_03_273_10_sq_w_cl_antia.root");

TChain* chain = new TChain("evt");
chain->Reset();
chain->Add(fname);
Long64_t nevents = chain->GetEntries();
cout << "total number of entries:\t" << nevents << endl;

//TClonesArray

TClonesArray* CalifaMappedData = new TClonesArray("R3BCalifaMappedData",3);
R3BCalifaMappedData** califamappeddata;
TBranch *branchCalifaMappedData = chain->GetBranch("CalifaMappedData");
branchCalifaMappedData->SetAddress(&CalifaMappedData);

TClonesArray* WRMasterData = new TClonesArray("R3BWRMasterData",1);
R3BWRMasterData** wrmasterdata;
TBranch *branchWRMasterData = chain->GetBranch("WRMasterData");
branchWRMasterData->SetAddress(&WRMasterData);

TClonesArray* CalifaHitData = new TClonesArray("R3BCalifaHitData",3);
R3BCalifaHitData** califahitdata;
TBranch *branchCalifaHitData = chain->GetBranch("CalifaHitData");
branchCalifaHitData->SetAddress(&CalifaHitData);

R3BEventHeader* DataCA = new R3BEventHeader();
TBranch* branchData = chain->GetBranch("EventHeader.");
branchData->SetAddress(&DataCA);

TClonesArray* SofToFWMappedData = new TClonesArray("R3BSofTofWMappedData",2);
R3BSofTofWMappedData** softofwmappeddata;
TBranch *branchSofToFWMappedData = chain->GetBranch("SofTofWMappedData");
branchSofToFWMappedData->SetAddress(&SofToFWMappedData);

TClonesArray* SofTwimHitData = new TClonesArray("R3BSofTwimHitData",2);
R3BSofTwimHitData** softwimhitdata;
TBranch *branchSofTwimHitData = chain->GetBranch("TwimHitData");
branchSofTwimHitData->SetAddress(&SofTwimHitData);

TClonesArray* SofToFWHitData = new TClonesArray("R3BSofTofWHitData",2);
R3BSofTofWHitData** softofwhitdata;
TBranch *branchSofToFWHitData = chain->GetBranch("TofWHitData");
branchSofToFWHitData->SetAddress(&SofToFWHitData);
// END of TClonesArray

const Double_t PI = 3.14159265358979323846;
uint64_t wr_Master_ts;
Long64_t events_califa_no_master = 0;
Long64_t events_califa_with_master = 0;
Long64_t events_total_califa = 0;
Long64_t events_with_cut = 0;
Long64_t events_one_proton = 0;
Long64_t events_two_proton = 0;

for(Long64_t i=0;i < 0.1*nevents;i++){
    Long64_t evtnr = i;
    if (i%100000==0)
        cout<<"Processing event for charge analysis "<<i<<endl;
    chain->GetEvent(i);

    entries_twim = SofTwimHitData->GetEntries();
    entries_califa_hit = CalifaHitData->GetEntries();
    
    if (entries_twim == 2){
    int trigger_pattern = t_pattern(DataCA);
    if (trigger_pattern >= 1){
    
    vector<double> charge_vector = twim_selection(SofTwimHitData,SofToFWHitData);
//    cout << "this is first entry:\t" << charge_vector[0] << endl;
//    cout << "this is second entry:\t" << charge_vector[1] << endl;
    if (charge_vector[0] > 0 && charge_vector[1] > 0){
//    cout << "before filling..." << endl;
    h2_charge1_charge_2->Fill(charge_vector[0],charge_vector[1]);
    }
   // if (entries_califa_hit > 2){
   //     Califa_p2p_selection my_one(CalifaHitData,40,100, 10);
   //     vector<double> testv = my_one.best_comb_hits();
   //     vector<double> tests = my_one.hight_energy_comb();
   //     int numbers_for_me = my_one.multiplicity_energy_cut(2);
   // }
    }
    }
}
cout << "number of events in histo:\t" << h2_charge1_charge_2->GetEntries() << endl;

char f_out_name[500];
sprintf(f_out_name,"../../file_output/twim_analysis_main0273_test.root");
write_to_file(f_out_name);
}
