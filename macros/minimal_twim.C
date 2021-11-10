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


#include <stdlib.h>
using namespace std;

void minimal_twim(){

//Input file
char fname[500];
char hist_name[500];
Long64_t entries_twim = 0;
sprintf(fname,"/u/land/tobias_jenegger/s455_gabri_tobias/file_src/s455_03_273_10_sq_w_cl_antia_newest_parameters.root");


//histogram

TH2D* h2_charge1_charge_2;
sprintf(hist_name, "Charge 1 vs Charge 2 (TWIM Music)");
h2_charge1_charge_2 = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
h2_charge1_charge_2->GetXaxis()->SetTitle("Charge 1");
h2_charge1_charge_2->GetYaxis()->SetTitle("Charge 2");
h2_charge1_charge_2->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2->GetYaxis()->SetTitleSize(0.045);

TChain* chain = new TChain("evt");
chain->Reset();
chain->Add(fname);
Long64_t nevents = chain->GetEntries();
cout << "total number of entries:\t" << nevents << endl;

//TClonesArray
TClonesArray* SofTwimHitData = new TClonesArray("R3BSofTwimHitData",2);
R3BSofTwimHitData** softwimhitdata;
TBranch *branchSofTwimHitData = chain->GetBranch("TwimHitData");
branchSofTwimHitData->SetAddress(&SofTwimHitData);


for(Long64_t i=0;i< nevents;i++){
    Long64_t evtnr = i;
    if (i%100000==0)
        cout<<"Processing event for charge analysis "<<i<<endl;
    chain->GetEvent(i);

    entries_twim = SofTwimHitData->GetEntries();
    
    if (entries_twim == 2){
	//check twim data
	softwimhitdata = new R3BSofTwimHitData*[2];
	softwimhitdata[0] = (R3BSofTwimHitData*)SofTwimHitData->At(0);
	softwimhitdata[1] = (R3BSofTwimHitData*)SofTwimHitData->At(1);
	Double_t charge_1 = softwimhitdata[0]->GetZcharge();
	Double_t charge_2 = softwimhitdata[1]->GetZcharge();
	Int_t fSection1 = softwimhitdata[0]->GetSecID();
	Int_t fSection2 = softwimhitdata[1]->GetSecID();
	if ((fSection1 == 1 && fSection2 == 3) || (fSection1 == 3 && fSection2 == 1)){	
		h2_charge1_charge_2->Fill(charge_1,charge_2);
		}
	}
}

char f_out_name[500];
sprintf(f_out_name,"../file_output/twim_analysis_main0273_10_with_antia_par_and_unpacker.root");
TFile * f = new TFile(f_out_name,"RECREATE");
TList *l = new TList();
l->Add(h2_charge1_charge_2);
l->Write("histlist", TObject::kSingleKey);
}
