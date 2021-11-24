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
#include <vector>

TH2D* h2_charge1_charge_2;

vector<TH2D*> list_of_histos;
char hist_name[500];
void my_histos(){
sprintf(hist_name, "Charge 1 vs Charge 2 (TWIM Music)");
h2_charge1_charge_2 = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
h2_charge1_charge_2->GetXaxis()->SetTitle("Charge 1");
h2_charge1_charge_2->GetYaxis()->SetTitle("Charge 2");
h2_charge1_charge_2->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2->GetYaxis()->SetTitleSize(0.045);
list_of_histos.push_back(h2_charge1_charge_2);

}
void write_to_file(char* output_file){
	TFile *f = new TFile(output_file,"RECREATE");
	TList *l = new TList();
	for (Int_t i = 0; i < list_of_histos.size(); i++){
		l->Add(list_of_histos[i]);
		}
	l->Write("histlist", TObject::kSingleKey);
	}
