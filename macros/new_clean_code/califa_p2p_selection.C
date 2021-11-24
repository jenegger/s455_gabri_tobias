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
#include<algorithm>
using namespace std;

bool sortcol( const vector<double>& v1,const vector<double>& v2 ) {
	return v1[0] > v2[0];
	}
bool sort_min( const vector<double>& vec1,const vector<double>& vec2 ) {
	return vec1[2] < vec2[2];
}
const Double_t PI = 3.14159265358979323846;
//void califa_p2p_selection(){
class Califa_p2p_selection {
	public:

		Califa_p2p_selection(TClonesArray* califa_tclone, Double_t e_cut, Double_t delta_phi_cut, Double_t delta_theta_cut);
		~Califa_p2p_selection(){delete [] califahitdata;}
		vector<double> best_comb_hits();
		vector<double> hight_energy_comb();
		int  multiplicity_energy_cut(double energy_cut);
		vector<vector<double> >get_califa_array();
		vector<double>get_highest_califa_gamma(double energy_cut_gamma);
		vector<uint64_t> get_califa_time();

	private:
		vector<vector<double> > v_phi_diff;
		Double_t d_phi;
		Double_t d_delta;
		Double_t d_e;
		Int_t entries_califa_hit;
		TClonesArray* my_tclone;
		R3BCalifaHitData** califahitdata;
		vector<vector<double> > v_califa_energy_theta_phi; //energy-theta-phi vector for all hits with energy > e_cut


	};

Califa_p2p_selection::Califa_p2p_selection(TClonesArray* califa_tclone, Double_t e_cut,Double_t delta_phi_cut, Double_t delta_theta_cut){
d_phi = delta_phi_cut;
d_delta = delta_theta_cut;
d_e = e_cut;
entries_califa_hit = califa_tclone->GetEntriesFast();
califahitdata = new R3BCalifaHitData*[entries_califa_hit];
my_tclone = califa_tclone;
	for (Int_t m = 0; m < entries_califa_hit; m++){
		califahitdata[m] = (R3BCalifaHitData*)my_tclone->At(m);
		vector<double> v_temp_califa_ext(4); //vector with energy-theta-phi values and index
		if (((califahitdata[m]->GetEnergy())/1000) > d_e){
			v_temp_califa_ext[0] = (califahitdata[m]->GetEnergy())/1000.;
			v_temp_califa_ext[1] = ((califahitdata[m]->GetTheta()/PI)*180.);
			v_temp_califa_ext[2] = ((califahitdata[m]->GetPhi())/PI)*180;
			v_temp_califa_ext[3] = m;
			v_califa_energy_theta_phi.push_back(v_temp_califa_ext);
			v_temp_califa_ext.clear();
			}

	}




	if (v_califa_energy_theta_phi.size() >= 2){
		sort(v_califa_energy_theta_phi.begin(),v_califa_energy_theta_phi.end(),sortcol);
		for(Int_t v1 = 0; v1 < 2; v1++){
			for (Int_t v2 = v1+1; v2 < v_califa_energy_theta_phi.size();v2++){
				vector<double> v_temp_phi_diff(5);
				v_temp_phi_diff[0] = v1;
				v_temp_phi_diff[1] = v2;
				if (v_califa_energy_theta_phi[v1][2] > v_califa_energy_theta_phi[v2][2]){
					Double_t angular_value1 = 360 + v_califa_energy_theta_phi[v2][2] - v_califa_energy_theta_phi[v1][2];
					Double_t angular_value2 = v_califa_energy_theta_phi[v1][2] - v_califa_energy_theta_phi[v2][2];
						if (angular_value1 < angular_value2){
							 v_temp_phi_diff[2] = 180 - angular_value1;
							}
						else {
							v_temp_phi_diff[2] = 180 - angular_value2;
							}
					}
				if ( v_califa_energy_theta_phi[v2][2] > v_califa_energy_theta_phi[v1][2]){
					Double_t angular_value1 = 360 + v_califa_energy_theta_phi[v1][2] - v_califa_energy_theta_phi[v2][2];
					Double_t angular_value2 = v_califa_energy_theta_phi[v2][2] - v_califa_energy_theta_phi[v1][2];	
						if (angular_value1 < angular_value2){
							v_temp_phi_diff[2] = 180 - angular_value1;
							}
						else {
							v_temp_phi_diff[2] = 180 - angular_value2;
							}

					}
				v_temp_phi_diff[3] = v_califa_energy_theta_phi[v1][3];
				v_temp_phi_diff[4] = v_califa_energy_theta_phi[v2][3];
				this->v_phi_diff.push_back(v_temp_phi_diff);
				v_temp_phi_diff.clear();
				}	
				}
			}
}
vector<double> Califa_p2p_selection::best_comb_hits(){
	vector<double> indices_delta_phi(3);
	if (v_phi_diff.size()){

	if (v_phi_diff[0][2] < d_phi){
		indices_delta_phi[0] = v_phi_diff[0][3]; //index of hit
		indices_delta_phi[1] = v_phi_diff[0][4]; //index of hit
		indices_delta_phi[2] = v_phi_diff[0][2]; //delta_phi
		}	
	if ( v_phi_diff.size() > 1 &&  v_phi_diff[0][2] > d_phi) {
		vector<vector<double> > best_comb_v_phi_diff = v_phi_diff; 
		sort(best_comb_v_phi_diff.begin(),best_comb_v_phi_diff.end(),sort_min);
		if (best_comb_v_phi_diff[0][2] < d_phi){
			indices_delta_phi[0] = best_comb_v_phi_diff[0][3];
			indices_delta_phi[1] = best_comb_v_phi_diff[0][4];
			indices_delta_phi[2] = best_comb_v_phi_diff[0][2];
			}
		else{
			indices_delta_phi={-1,-1,-1};
			}
		}
	if ( v_phi_diff.size() == 1 && v_phi_diff[0][2] > d_phi){
		indices_delta_phi={-1,-1,-1};
		}
	}
	else	{
		indices_delta_phi={-1,-1,-1};
		}
	return indices_delta_phi;
	}

vector<double> Califa_p2p_selection::hight_energy_comb(){
	vector<double> indices_delta_phi_hight(3);
	if (v_phi_diff.size()){
		indices_delta_phi_hight[0] = v_phi_diff[0][3];
		indices_delta_phi_hight[1] = v_phi_diff[0][4];
		indices_delta_phi_hight[2] = v_phi_diff[0][2];
		}
	else {
		indices_delta_phi_hight = {-1,-1,-1};
		}
	return indices_delta_phi_hight;
	}
int Califa_p2p_selection::multiplicity_energy_cut(double energy_cut){
	int multiplicity_value = 0;
        for (Int_t m = 0; m < entries_califa_hit; m++){
	        califahitdata[m] = (R3BCalifaHitData*)my_tclone->At(m);
		if (((califahitdata[m]->GetEnergy())/1000) > energy_cut){	
			multiplicity_value+=1;
		}
	}
	return multiplicity_value;
}

vector<vector<double> > Califa_p2p_selection::get_califa_array(){
	vector<vector<double> >califa_array;
	Double_t gamma = (1+(637./938));
	Double_t beta = sqrt(1-(1/(gamma*gamma)));
	for (Int_t m = 0; m < entries_califa_hit; m++){
		vector<double> m_array(4);
		califahitdata[m] = (R3BCalifaHitData*)my_tclone->At(m);		
		m_array[0] = ((califahitdata[m]->GetEnergy())/1000);		//Energy in MeV
		m_array[1] = ((califahitdata[m]->GetTheta()/PI)*180.);		//Theta in degrees
		m_array[2] = ((califahitdata[m]->GetPhi())/PI)*180;		//Phi in degrees
		m_array[3] = m_array[0]*gamma*(1-beta*cos(m_array[1]*PI/180.));	//Doppler corrected energy in MeV
		califa_array.push_back(m_array);
		}
	return califa_array;
	}
vector<uint64_t> Califa_p2p_selection::get_califa_time(){
	vector <uint64_t> time_array_califa;
	for (Int_t m = 0; m < entries_califa_hit; m++){
		califahitdata[m] = (R3BCalifaHitData*)my_tclone->At(m);
		uint64_t time_califa = (califahitdata[m]->GetTime());
		time_array_califa.push_back(time_califa);
		}
	return time_array_califa;
	}

vector<double> Califa_p2p_selection::get_highest_califa_gamma(double energy_cut_gamma){
	vector<vector<double> >califa_array;
	Double_t gamma = (1+(637./938));
	Double_t beta = sqrt(1-(1/(gamma*gamma)));
	vector<double> high_gamma_hit={-1,-1,-1,-1};
	for (Int_t m = 0; m < entries_califa_hit; m++){
		vector<double> m_array(4);
		califahitdata[m] = (R3BCalifaHitData*)my_tclone->At(m);		
		m_array[0] = ((califahitdata[m]->GetEnergy())/1000);		//Energy in MeV
		m_array[1] = ((califahitdata[m]->GetTheta()/PI)*180.);		//Theta in degrees
		m_array[2] = ((califahitdata[m]->GetPhi())/PI)*180;		//Phi in degrees
		m_array[3] = m_array[0]*gamma*(1-beta*cos(m_array[1]*PI/180.));	//Doppler corrected energy in MeV
		if (m_array[0] < energy_cut_gamma){
		califa_array.push_back(m_array);
		}
	}
	if (califa_array.size() == 1){
		for (Int_t i = 0; i < 4; i++){
			high_gamma_hit[i] = califa_array[0][i];
			}
		}
	if (califa_array.size() > 1){
		sort(califa_array.begin(),califa_array.end(),sortcol);
		for (Int_t i = 0; i < 4; i++){
			high_gamma_hit[i] = califa_array[0][i];
			}
		}
	
	return high_gamma_hit;
	}

void califa_p2p_selection(){

}
