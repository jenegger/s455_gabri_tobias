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

bool sortcol( const vector<double>& v1,const vector<double>& v2 ) {
	                                    return v1[0] > v2[0];
					}
bool sort_min( const vector<double>& vec1,const vector<double>& vec2 ) {
	                                            return vec1[2] < vec2[2];
					}

double array_for_cuts_theta[6][3] = {};
double array_for_cuts_phi[6][3] = {};
double array_for_cuts_energy[6][3] = {};



char hist_name[500];
//similar code to wr_time_diff.C, here with restriction that I expect at least one hit and at most 2two hits in tofw
void corr_tof_check(){


TH1D* h1_wr_ts_califa_diff_WRM_messel_rand;
sprintf(hist_name, "WR Messel -  WR Master (one random hit for event)");
h1_wr_ts_califa_diff_WRM_messel_rand = new TH1D(hist_name,hist_name,500,-5000,5000);
h1_wr_ts_califa_diff_WRM_messel_rand->GetXaxis()->SetTitle(" wr time difference [ns]");
h1_wr_ts_califa_diff_WRM_messel_rand->GetYaxis()->SetTitle("Counts");
h1_wr_ts_califa_diff_WRM_messel_rand->GetXaxis()->CenterTitle(true);
h1_wr_ts_califa_diff_WRM_messel_rand->GetYaxis()->CenterTitle(true);
h1_wr_ts_califa_diff_WRM_messel_rand->GetYaxis()->SetLabelSize(0.045);
h1_wr_ts_califa_diff_WRM_messel_rand->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_wr_ts_califa_diff_WRM_wix_rand;
sprintf(hist_name, "WR Wixhausen - WR Master (one random hit for event)");
h1_wr_ts_califa_diff_WRM_wix_rand = new TH1D(hist_name,hist_name,500,-5000,5000);
h1_wr_ts_califa_diff_WRM_wix_rand->GetXaxis()->SetTitle("wr time difference [ns]");
h1_wr_ts_califa_diff_WRM_wix_rand->GetYaxis()->SetTitle("Counts");
h1_wr_ts_califa_diff_WRM_wix_rand->GetXaxis()->CenterTitle(true);
h1_wr_ts_califa_diff_WRM_wix_rand->GetYaxis()->CenterTitle(true);
h1_wr_ts_califa_diff_WRM_wix_rand->GetYaxis()->SetLabelSize(0.045);
h1_wr_ts_califa_diff_WRM_wix_rand->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_wr_ts_califa_diff_WRM_both_rand;
sprintf(hist_name, "WR CALIFA - WR Master (one random hit for event wixhausen or messel)");
h1_wr_ts_califa_diff_WRM_both_rand = new TH1D(hist_name,hist_name,500,-5000,5000);
h1_wr_ts_califa_diff_WRM_both_rand->GetXaxis()->SetTitle("wr time difference [ns]");
h1_wr_ts_califa_diff_WRM_both_rand->GetYaxis()->SetTitle("Counts");
h1_wr_ts_califa_diff_WRM_both_rand->GetXaxis()->CenterTitle(true);
h1_wr_ts_califa_diff_WRM_both_rand->GetYaxis()->CenterTitle(true);
h1_wr_ts_califa_diff_WRM_both_rand->GetYaxis()->SetLabelSize(0.045);
h1_wr_ts_califa_diff_WRM_both_rand->GetYaxis()->SetTitleSize(0.045);

//Energy vs ts difference for events with wr master trigger
TH1D* h1_wr_ts_califa_vs_E_M;
sprintf(hist_name, "Hit Energy (events with master trigger)");
h1_wr_ts_califa_vs_E_M = new TH1D(hist_name,hist_name,1000,0,10);
h1_wr_ts_califa_vs_E_M->GetXaxis()->SetTitle("Energy [MeV]");
h1_wr_ts_califa_vs_E_M->GetYaxis()->SetTitle("Counts");
h1_wr_ts_califa_vs_E_M->GetXaxis()->CenterTitle(true);
h1_wr_ts_califa_vs_E_M->GetYaxis()->CenterTitle(true);
h1_wr_ts_califa_vs_E_M->GetYaxis()->SetLabelSize(0.045);
h1_wr_ts_califa_vs_E_M->GetYaxis()->SetTitleSize(0.045);

//Energy vs ts difference for events without wr master trigger
TH1D* h1_wr_ts_califa_vs_E_noM;
sprintf(hist_name, "Hit Energy (events NO master trigger)");
h1_wr_ts_califa_vs_E_noM = new TH1D(hist_name,hist_name,1000,0,10);
h1_wr_ts_califa_vs_E_noM->GetXaxis()->SetTitle("Energy [MeV]");
h1_wr_ts_califa_vs_E_noM->GetYaxis()->SetTitle("Counts");
h1_wr_ts_califa_vs_E_noM->GetXaxis()->CenterTitle(true);
h1_wr_ts_califa_vs_E_noM->GetYaxis()->CenterTitle(true);
h1_wr_ts_califa_vs_E_noM->GetYaxis()->SetLabelSize(0.045);
h1_wr_ts_califa_vs_E_noM->GetYaxis()->SetTitleSize(0.045);


//Multiplicity for events with wr master

TH1D* h1_mult_with_M;
sprintf(hist_name, "Multiplicity for events with master trigger");
h1_mult_with_M = new TH1D(hist_name,hist_name,2500,0,2500);
h1_mult_with_M->GetXaxis()->SetTitle("Multiplicity");
h1_mult_with_M->GetYaxis()->SetTitle("Counts");
h1_mult_with_M->GetXaxis()->CenterTitle(true);
h1_mult_with_M->GetYaxis()->CenterTitle(true);
h1_mult_with_M->GetYaxis()->SetLabelSize(0.045);
h1_mult_with_M->GetYaxis()->SetTitleSize(0.045);


//Multiplicity for events without wr master

TH1D* h1_mult_no_M;
sprintf(hist_name, "Multiplicity for events NO master trigger");
h1_mult_no_M = new TH1D(hist_name,hist_name,2500,0,2500);
h1_mult_no_M->GetXaxis()->SetTitle("Multiplicity");
h1_mult_no_M->GetYaxis()->SetTitle("Counts");
h1_mult_no_M->GetXaxis()->CenterTitle(true);
h1_mult_no_M->GetYaxis()->CenterTitle(true);
h1_mult_no_M->GetYaxis()->SetLabelSize(0.045);
h1_mult_no_M->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_mult_e_high;
sprintf(hist_name, "Multiplicity of hits with E_lab > 20 MeV (with master trigger)");
h1_mult_e_high = new TH1D(hist_name,hist_name,2500,0,2500);
h1_mult_e_high->GetXaxis()->SetTitle("Multiplicity");
h1_mult_e_high->GetYaxis()->SetTitle("Counts");
h1_mult_e_high->GetXaxis()->CenterTitle(true);
h1_mult_e_high->GetYaxis()->CenterTitle(true);
h1_mult_e_high->GetYaxis()->SetLabelSize(0.045);
h1_mult_e_high->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_mult_e_high_noM;
sprintf(hist_name, "Multiplicity of hits with E_lab > 20 MeV (events without master trigger)");
h1_mult_e_high_noM = new TH1D(hist_name,hist_name,2500,0,2500);
h1_mult_e_high_noM->GetXaxis()->SetTitle("Multiplicity");
h1_mult_e_high_noM->GetYaxis()->SetTitle("Counts");
h1_mult_e_high_noM->GetXaxis()->CenterTitle(true);
h1_mult_e_high_noM->GetYaxis()->CenterTitle(true);
h1_mult_e_high_noM->GetYaxis()->SetLabelSize(0.045);
h1_mult_e_high_noM->GetYaxis()->SetTitleSize(0.045);


//check Energy versus polar angle distribution
//
TH2D* h2_energy_vs_theta;
sprintf(hist_name, "E_lab vs. polar angle (with master trigger)");
h2_energy_vs_theta = new TH2D(hist_name,hist_name,1000,0,1000,25,22.15,84.65);
h2_energy_vs_theta->GetXaxis()->SetTitle("Energy in laboratory system [MeV]");
h2_energy_vs_theta->GetYaxis()->SetTitle("Polar angle Theta [degrees]");
h2_energy_vs_theta->GetXaxis()->CenterTitle(true);
h2_energy_vs_theta->GetYaxis()->CenterTitle(true);
h2_energy_vs_theta->GetYaxis()->SetLabelSize(0.045);
h2_energy_vs_theta->GetYaxis()->SetTitleSize(0.045);

//check angular distribution in phi and theta

TH2D* h2_phi_vs_theta;
sprintf(hist_name, "Phi vs Theta (with master trigger)");
h2_phi_vs_theta = new TH2D(hist_name,hist_name,25,22.15,84.65,60,-180,180);
h2_phi_vs_theta->GetXaxis()->SetTitle("Polar angle  Theta[degrees]");
h2_phi_vs_theta->GetYaxis()->SetTitle("Arzimuthal angle Phi [degrees]");
h2_phi_vs_theta->GetXaxis()->CenterTitle(true);
h2_phi_vs_theta->GetYaxis()->CenterTitle(true);
h2_phi_vs_theta->GetYaxis()->SetLabelSize(0.045);
h2_phi_vs_theta->GetYaxis()->SetTitleSize(0.045);

//energy versus multiplicity distribution with and without master trigger
TH2D* h2_energy_vs_mult_no_M;
sprintf(hist_name, "Summed up E_lab(mapped level) vs Multiplicity, Uncorrelated");
h2_energy_vs_mult_no_M = new TH2D(hist_name,hist_name,1000,0,1000,400,0,400);
h2_energy_vs_mult_no_M->GetXaxis()->SetTitle("Summed Energy in laboratory system [MeV]");
h2_energy_vs_mult_no_M->GetYaxis()->SetTitle("Multiplicity");
h2_energy_vs_mult_no_M->GetXaxis()->CenterTitle(true);
h2_energy_vs_mult_no_M->GetYaxis()->CenterTitle(true);
h2_energy_vs_mult_no_M->GetYaxis()->SetLabelSize(0.045);
h2_energy_vs_mult_no_M->GetYaxis()->SetTitleSize(0.045);


TH2D* h2_energy_vs_mult_with_M;
sprintf(hist_name, "Summed up E_lab(mapped level) vs Multiplicity, Correlated");
h2_energy_vs_mult_with_M = new TH2D(hist_name,hist_name,1000,0,1000,400,0,400);
h2_energy_vs_mult_with_M->GetXaxis()->SetTitle("Summed Energy in laboratory system [MeV]");
h2_energy_vs_mult_with_M->GetYaxis()->SetTitle("Multiplicity");
h2_energy_vs_mult_with_M->GetXaxis()->CenterTitle(true);
h2_energy_vs_mult_with_M->GetYaxis()->CenterTitle(true);
h2_energy_vs_mult_with_M->GetYaxis()->SetLabelSize(0.045);
h2_energy_vs_mult_with_M->GetYaxis()->SetTitleSize(0.045);


//charge analysis

TH2D* h2_charge1_charge_2;
sprintf(hist_name, "Charge 1 vs Charge 2 (TWIM Music)");
h2_charge1_charge_2 = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
h2_charge1_charge_2->GetXaxis()->SetTitle("Charge 1");
h2_charge1_charge_2->GetYaxis()->SetTitle("Charge 2");
h2_charge1_charge_2->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_charge1_plus_charge2;
sprintf(hist_name, "Charge 1 + Charge 2  (TWIM Music)");
h1_charge1_plus_charge2 = new TH1D(hist_name,hist_name,500,0,100);
h1_charge1_plus_charge2->GetXaxis()->SetTitle("Charge Sum");
h1_charge1_plus_charge2->GetYaxis()->SetTitle("Counts");
h1_charge1_plus_charge2->GetXaxis()->CenterTitle(true);
h1_charge1_plus_charge2->GetYaxis()->CenterTitle(true);
h1_charge1_plus_charge2->GetYaxis()->SetLabelSize(0.045);
h1_charge1_plus_charge2->GetYaxis()->SetTitleSize(0.045);



TH2D* h2_charge1_charge_2_diff_vs_z_sum;
sprintf(hist_name, "Charge 1 - Charge 2 vs Charge_sum(TWIM Music)");
h2_charge1_charge_2_diff_vs_z_sum = new TH2D(hist_name,hist_name,1400,-70,70,1000,0,100);
h2_charge1_charge_2_diff_vs_z_sum->GetXaxis()->SetTitle("Charge 1 - Charge 2");
h2_charge1_charge_2_diff_vs_z_sum->GetYaxis()->SetTitle("Charge 1 + Charge 2");
h2_charge1_charge_2_diff_vs_z_sum->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_diff_vs_z_sum->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_daughter_nuclei_charge;
sprintf(hist_name, "Charge of Daughter Nuclei (TWIM Music)");
h1_daughter_nuclei_charge = new TH1D(hist_name,hist_name,500,0,100);
h1_daughter_nuclei_charge->GetXaxis()->SetTitle("Charge Z");
h1_daughter_nuclei_charge->GetYaxis()->SetTitle("Counts");
h1_daughter_nuclei_charge->GetXaxis()->CenterTitle(true);
h1_daughter_nuclei_charge->GetYaxis()->CenterTitle(true);
h1_daughter_nuclei_charge->GetYaxis()->SetLabelSize(0.045);
h1_daughter_nuclei_charge->GetYaxis()->SetTitleSize(0.045);



TH2D* h2_charge1_charge_2_after;
sprintf(hist_name, "Charge 1 vs Charge 2 (TWIM Music) after the minimal CALIFA cuts");
h2_charge1_charge_2_after = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
h2_charge1_charge_2_after->GetXaxis()->SetTitle("Charge 1");
h2_charge1_charge_2_after->GetYaxis()->SetTitle("Charge 2");
h2_charge1_charge_2_after->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_after->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_after->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_after->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge_sum_after;
sprintf(hist_name, "Charge Sum (TWIM Music) after the minimal CALIFA cuts");
h1_charge_sum_after = new TH1D(hist_name,hist_name,500,0,100);
h1_charge_sum_after->GetXaxis()->SetTitle("Charge 1 + Charge 2");
h1_charge_sum_after->GetYaxis()->SetTitle("Counts");
h1_charge_sum_after->GetXaxis()->CenterTitle(true);
h1_charge_sum_after->GetYaxis()->CenterTitle(true);
h1_charge_sum_after->GetYaxis()->SetLabelSize(0.045);
h1_charge_sum_after->GetYaxis()->SetTitleSize(0.045);


//make charge analysis for different trigger patterns

TH2D* h2_charge1_charge_2_tpat1;
sprintf(hist_name, "Charge 1 vs Charge 2 (TWIM Music) for TPat = SofStart");
h2_charge1_charge_2_tpat1 = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
h2_charge1_charge_2_tpat1->GetXaxis()->SetTitle("Charge 1");
h2_charge1_charge_2_tpat1->GetYaxis()->SetTitle("Charge 2");
h2_charge1_charge_2_tpat1->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_tpat1->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_tpat1->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_tpat1->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_charge1_charge_2_tpat2;
sprintf(hist_name, "Charge 1 vs Charge 2 (TWIM Music) for TPat =  Fission");
h2_charge1_charge_2_tpat2 = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
h2_charge1_charge_2_tpat2->GetXaxis()->SetTitle("Charge 1");
h2_charge1_charge_2_tpat2->GetYaxis()->SetTitle("Charge 2");
h2_charge1_charge_2_tpat2->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_tpat2->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_tpat2->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_tpat2->GetYaxis()->SetTitleSize(0.045);


TH2D* h2_charge1_charge_2_tpat3;
sprintf(hist_name, "Charge 1 vs Charge 2 (TWIM Music) for TPat =  CalifOR (= S+C&&)");
h2_charge1_charge_2_tpat3 = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
h2_charge1_charge_2_tpat3->GetXaxis()->SetTitle("Charge 1");
h2_charge1_charge_2_tpat3->GetYaxis()->SetTitle("Charge 2");
h2_charge1_charge_2_tpat3->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_tpat3->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_tpat3->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_tpat3->GetYaxis()->SetTitleSize(0.045);


TH2D* h2_charge1_charge_2_tpat4;
sprintf(hist_name, "Charge 1 vs Charge 2 (TWIM Music) for TPat =  CalifAND (= S+F+C&&)");
h2_charge1_charge_2_tpat4 = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
h2_charge1_charge_2_tpat4->GetXaxis()->SetTitle("Charge 1");
h2_charge1_charge_2_tpat4->GetYaxis()->SetTitle("Charge 2");
h2_charge1_charge_2_tpat4->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_tpat4->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_tpat4->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_tpat4->GetYaxis()->SetTitleSize(0.045);


TH2D* h2_charge1_charge_2_one_CALIFA;
sprintf(hist_name, "Charge 1 vs Charge 2 (TWIM Music) for one hit > 200 MeV in CALIFA");
h2_charge1_charge_2_one_CALIFA = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
h2_charge1_charge_2_one_CALIFA->GetXaxis()->SetTitle("Charge 1");
h2_charge1_charge_2_one_CALIFA->GetYaxis()->SetTitle("Charge 2");
h2_charge1_charge_2_one_CALIFA->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_one_CALIFA->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_one_CALIFA->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_one_CALIFA->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_charge1_charge_2_two_CALIFA;
sprintf(hist_name, "Charge 1 vs Charge 2 (TWIM Music) for two hits > 200 MeV in CALIFA");
h2_charge1_charge_2_two_CALIFA = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
h2_charge1_charge_2_two_CALIFA->GetXaxis()->SetTitle("Charge 1");
h2_charge1_charge_2_two_CALIFA->GetYaxis()->SetTitle("Charge 2");
h2_charge1_charge_2_two_CALIFA->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_two_CALIFA->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_two_CALIFA->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_two_CALIFA->GetYaxis()->SetTitleSize(0.045);


TH2D* h2_charge1_charge_2_three_CALIFA;
sprintf(hist_name, "Charge 1 vs Charge 2 (TWIM Music) for three hits > 200 MeV in CALIFA");
h2_charge1_charge_2_three_CALIFA = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
h2_charge1_charge_2_three_CALIFA->GetXaxis()->SetTitle("Charge 1");
h2_charge1_charge_2_three_CALIFA->GetYaxis()->SetTitle("Charge 2");
h2_charge1_charge_2_three_CALIFA->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_three_CALIFA->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_three_CALIFA->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_three_CALIFA->GetYaxis()->SetTitleSize(0.045);


TH2D* h2_charge1_charge_2_diff_vs_z_sum_tpat1;
sprintf(hist_name, "Charge 1 - Charge 2 vs Charge_sum(TWIM Music) for TPat = SofStart");
h2_charge1_charge_2_diff_vs_z_sum_tpat1 = new TH2D(hist_name,hist_name,1400,-70,70,1000,0,100);
h2_charge1_charge_2_diff_vs_z_sum_tpat1->GetXaxis()->SetTitle("Charge 1 - Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_tpat1->GetYaxis()->SetTitle("Charge 1 + Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_tpat1->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_tpat1->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_tpat1->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_diff_vs_z_sum_tpat1->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_charge1_charge_2_diff_vs_z_sum_tpat2;
sprintf(hist_name, "Charge 1 - Charge 2 vs Charge_sum(TWIM Music) for TPat =  Fission");
h2_charge1_charge_2_diff_vs_z_sum_tpat2 = new TH2D(hist_name,hist_name,1400,-70,70,1000,0,100);
h2_charge1_charge_2_diff_vs_z_sum_tpat2->GetXaxis()->SetTitle("Charge 1 - Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_tpat2->GetYaxis()->SetTitle("Charge 1 + Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_tpat2->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_tpat2->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_tpat2->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_diff_vs_z_sum_tpat2->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_charge1_charge_2_diff_vs_z_sum_tpat3;
sprintf(hist_name, "Charge 1 - Charge 2 vs Charge_sum(TWIM Music) for TPat =  CalifOR (= S+C&&)");
h2_charge1_charge_2_diff_vs_z_sum_tpat3 = new TH2D(hist_name,hist_name,1400,-70,70,1000,0,100);
h2_charge1_charge_2_diff_vs_z_sum_tpat3->GetXaxis()->SetTitle("Charge 1 - Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_tpat3->GetYaxis()->SetTitle("Charge 1 + Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_tpat3->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_tpat3->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_tpat3->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_diff_vs_z_sum_tpat3->GetYaxis()->SetTitleSize(0.045);


TH2D* h2_charge1_charge_2_diff_vs_z_sum_tpat4;
sprintf(hist_name, "Charge 1 - Charge 2 vs Charge_sum(TWIM Music) for TPat =  CalifAND (= S+F+C&&)");
h2_charge1_charge_2_diff_vs_z_sum_tpat4 = new TH2D(hist_name,hist_name,1400,-70,70,1000,0,100);
h2_charge1_charge_2_diff_vs_z_sum_tpat4->GetXaxis()->SetTitle("Charge 1 - Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_tpat4->GetYaxis()->SetTitle("Charge 1 + Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_tpat4->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_tpat4->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_tpat4->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_diff_vs_z_sum_tpat4->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_charge1_charge_2_diff_vs_z_sum_one_CALIFA;
sprintf(hist_name, "Charge 1 - Charge 2 vs Charge_sum(TWIM Music) for one hit > 200MeV in CALIFA");
h2_charge1_charge_2_diff_vs_z_sum_one_CALIFA = new TH2D(hist_name,hist_name,1400,-70,70,1000,0,100);
h2_charge1_charge_2_diff_vs_z_sum_one_CALIFA->GetXaxis()->SetTitle("Charge 1 - Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_one_CALIFA->GetYaxis()->SetTitle("Charge 1 + Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_one_CALIFA->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_one_CALIFA->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_one_CALIFA->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_diff_vs_z_sum_one_CALIFA->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_charge1_charge_2_diff_vs_z_sum_two_CALIFA;
sprintf(hist_name, "Charge 1 - Charge 2 vs Charge_sum(TWIM Music) for two hits > 200MeV in CALIFA");
h2_charge1_charge_2_diff_vs_z_sum_two_CALIFA = new TH2D(hist_name,hist_name,1400,-70,70,1000,0,100);
h2_charge1_charge_2_diff_vs_z_sum_two_CALIFA->GetXaxis()->SetTitle("Charge 1 - Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_two_CALIFA->GetYaxis()->SetTitle("Charge 1 + Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_two_CALIFA->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_two_CALIFA->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_two_CALIFA->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_diff_vs_z_sum_two_CALIFA->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_charge1_charge_2_diff_vs_z_sum_three_CALIFA;
sprintf(hist_name, "Charge 1 - Charge 2 vs Charge_sum(TWIM Music) for three hits > 200MeV in CALIFA");
h2_charge1_charge_2_diff_vs_z_sum_three_CALIFA = new TH2D(hist_name,hist_name,1400,-70,70,1000,0,100);
h2_charge1_charge_2_diff_vs_z_sum_three_CALIFA->GetXaxis()->SetTitle("Charge 1 - Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_three_CALIFA->GetYaxis()->SetTitle("Charge 1 + Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_three_CALIFA->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_three_CALIFA->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_three_CALIFA->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_diff_vs_z_sum_three_CALIFA->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_charge1_plus_charge2_tpat1;
sprintf(hist_name, "Charge 1 + Charge 2  (TWIM Music) for TPat = SofStart");
h1_charge1_plus_charge2_tpat1 = new TH1D(hist_name,hist_name,500,0,100);
h1_charge1_plus_charge2_tpat1->GetXaxis()->SetTitle("Charge Sum");
h1_charge1_plus_charge2_tpat1->GetYaxis()->SetTitle("Counts");
h1_charge1_plus_charge2_tpat1->GetXaxis()->CenterTitle(true);
h1_charge1_plus_charge2_tpat1->GetYaxis()->CenterTitle(true);
h1_charge1_plus_charge2_tpat1->GetYaxis()->SetLabelSize(0.045);
h1_charge1_plus_charge2_tpat1->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_charge1_plus_charge2_tpat2;
sprintf(hist_name, "Charge 1 + Charge 2  (TWIM Music) for TPat = Fission");
h1_charge1_plus_charge2_tpat2 = new TH1D(hist_name,hist_name,500,0,100);
h1_charge1_plus_charge2_tpat2->GetXaxis()->SetTitle("Charge Sum");
h1_charge1_plus_charge2_tpat2->GetYaxis()->SetTitle("Counts");
h1_charge1_plus_charge2_tpat2->GetXaxis()->CenterTitle(true);
h1_charge1_plus_charge2_tpat2->GetYaxis()->CenterTitle(true);
h1_charge1_plus_charge2_tpat2->GetYaxis()->SetLabelSize(0.045);
h1_charge1_plus_charge2_tpat2->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge1_plus_charge2_tpat3;
sprintf(hist_name, "Charge 1 + Charge 2  (TWIM Music) for TPat =  CalifOR (= S+C&&)");
h1_charge1_plus_charge2_tpat3 = new TH1D(hist_name,hist_name,500,0,100);
h1_charge1_plus_charge2_tpat3->GetXaxis()->SetTitle("Charge Sum");
h1_charge1_plus_charge2_tpat3->GetYaxis()->SetTitle("Counts");
h1_charge1_plus_charge2_tpat3->GetXaxis()->CenterTitle(true);
h1_charge1_plus_charge2_tpat3->GetYaxis()->CenterTitle(true);
h1_charge1_plus_charge2_tpat3->GetYaxis()->SetLabelSize(0.045);
h1_charge1_plus_charge2_tpat3->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge1_plus_charge2_tpat4;
sprintf(hist_name, "Charge 1 + Charge 2  (TWIM Music) for TPat =  CalifAND (= S+F+C&&)");
h1_charge1_plus_charge2_tpat4 = new TH1D(hist_name,hist_name,500,0,100);
h1_charge1_plus_charge2_tpat4->GetXaxis()->SetTitle("Charge Sum");
h1_charge1_plus_charge2_tpat4->GetYaxis()->SetTitle("Counts");
h1_charge1_plus_charge2_tpat4->GetXaxis()->CenterTitle(true);
h1_charge1_plus_charge2_tpat4->GetYaxis()->CenterTitle(true);
h1_charge1_plus_charge2_tpat4->GetYaxis()->SetLabelSize(0.045);
h1_charge1_plus_charge2_tpat4->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_charge1_plus_charge2_one_CALIFA;
sprintf(hist_name, "Charge 1 + Charge 2  (TWIM Music) for one hit > 200 MeV in CALIFA");
h1_charge1_plus_charge2_one_CALIFA = new TH1D(hist_name,hist_name,500,0,100);
h1_charge1_plus_charge2_one_CALIFA->GetXaxis()->SetTitle("Charge Sum");
h1_charge1_plus_charge2_one_CALIFA->GetYaxis()->SetTitle("Counts");
h1_charge1_plus_charge2_one_CALIFA->GetXaxis()->CenterTitle(true);
h1_charge1_plus_charge2_one_CALIFA->GetYaxis()->CenterTitle(true);
h1_charge1_plus_charge2_one_CALIFA->GetYaxis()->SetLabelSize(0.045);
h1_charge1_plus_charge2_one_CALIFA->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge1_plus_charge2_two_CALIFA;
sprintf(hist_name, "Charge 1 + Charge 2  (TWIM Music) for two hits > 200 MeV in CALIFA");
h1_charge1_plus_charge2_two_CALIFA = new TH1D(hist_name,hist_name,500,0,100);
h1_charge1_plus_charge2_two_CALIFA->GetXaxis()->SetTitle("Charge Sum");
h1_charge1_plus_charge2_two_CALIFA->GetYaxis()->SetTitle("Counts");
h1_charge1_plus_charge2_two_CALIFA->GetXaxis()->CenterTitle(true);
h1_charge1_plus_charge2_two_CALIFA->GetYaxis()->CenterTitle(true);
h1_charge1_plus_charge2_two_CALIFA->GetYaxis()->SetLabelSize(0.045);
h1_charge1_plus_charge2_two_CALIFA->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_charge1_plus_charge2_three_CALIFA;
sprintf(hist_name, "Charge 1 + Charge 2  (TWIM Music) for three hits > 200 MeV in CALIFA");
h1_charge1_plus_charge2_three_CALIFA = new TH1D(hist_name,hist_name,500,0,100);
h1_charge1_plus_charge2_three_CALIFA->GetXaxis()->SetTitle("Charge Sum");
h1_charge1_plus_charge2_three_CALIFA->GetYaxis()->SetTitle("Counts");
h1_charge1_plus_charge2_three_CALIFA->GetXaxis()->CenterTitle(true);
h1_charge1_plus_charge2_three_CALIFA->GetYaxis()->CenterTitle(true);
h1_charge1_plus_charge2_three_CALIFA->GetYaxis()->SetLabelSize(0.045);
h1_charge1_plus_charge2_three_CALIFA->GetYaxis()->SetTitleSize(0.045);

//check 3d plot for energy distribution for p2p reactions
Int_t p2p_check = 20;
Int_t p2p_check_one = 20;
Int_t p2p_check_two = 20;
TH2D* h2_energy_angular_dist[p2p_check];
for (Int_t i = 0; i < p2p_check ; i++){
sprintf(hist_name,"Energy distribution over theta_phi for p2p eventnr %i",i);
h2_energy_angular_dist[i] = new TH2D(hist_name,hist_name,52,22.15,152.15,60,-180,180);
h2_energy_angular_dist[i]->GetXaxis()->SetTitle("Theta [degr]");
h2_energy_angular_dist[i]->GetYaxis()->SetTitle("Phi [degr]");
h2_energy_angular_dist[i]->GetXaxis()->CenterTitle(true);
h2_energy_angular_dist[i]->GetYaxis()->CenterTitle(true);
h2_energy_angular_dist[i]->GetXaxis()->SetLabelSize(0.045);
h2_energy_angular_dist[i]->GetYaxis()->SetLabelSize(0.045);
}


TH2D* h2_energy_angular_dist_one[p2p_check_one];
for (Int_t i = 0; i < p2p_check ; i++){
	sprintf(hist_name,"Energy distribution over theta_phi for p2p eventnr  %i with one hit > 100 MeV",i);
	h2_energy_angular_dist_one[i] = new TH2D(hist_name,hist_name,52,22.15,152.15,60,-180,180);
	h2_energy_angular_dist_one[i]->GetXaxis()->SetTitle("Theta [degr]");
	h2_energy_angular_dist_one[i]->GetYaxis()->SetTitle("Phi [degr]");
	h2_energy_angular_dist_one[i]->GetXaxis()->CenterTitle(true);
	h2_energy_angular_dist_one[i]->GetYaxis()->CenterTitle(true);
	h2_energy_angular_dist_one[i]->GetXaxis()->SetLabelSize(0.045);
	h2_energy_angular_dist_one[i]->GetYaxis()->SetLabelSize(0.045);
}

TH2D* h2_energy_angular_dist_two[p2p_check_two];
for (Int_t i = 0; i < p2p_check ; i++){
	sprintf(hist_name,"Energy distribution over theta_phi for p2p eventnr %i with two hits > 100 MeV",i);
	h2_energy_angular_dist_two[i] = new TH2D(hist_name,hist_name,52,22.15,152.15,60,-180,180);
	h2_energy_angular_dist_two[i]->GetXaxis()->SetTitle("Theta [degr]");
	h2_energy_angular_dist_two[i]->GetYaxis()->SetTitle("Phi [degr]");
	h2_energy_angular_dist_two[i]->GetXaxis()->CenterTitle(true);
	h2_energy_angular_dist_two[i]->GetYaxis()->CenterTitle(true);
	h2_energy_angular_dist_two[i]->GetXaxis()->SetLabelSize(0.045);
	h2_energy_angular_dist_two[i]->GetYaxis()->SetLabelSize(0.045);
}

TH2D* h2_charge1_vs_charge2_energy_cut[6];
for (Int_t i = 0; i < 6;i++){
	sprintf(hist_name, "Charge1 vs Charge 2 for Energy cut: %i < E_sum < %i ",i*50+50, 900-i*50);
	h2_charge1_vs_charge2_energy_cut[i] = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
	h2_charge1_vs_charge2_energy_cut[i]->GetXaxis()->SetTitle("Charge 1 ");
	h2_charge1_vs_charge2_energy_cut[i]->GetYaxis()->SetTitle("Charge 2 ");
	h2_charge1_vs_charge2_energy_cut[i]->GetXaxis()->CenterTitle(true);
	h2_charge1_vs_charge2_energy_cut[i]->GetYaxis()->CenterTitle(true);
	h2_charge1_vs_charge2_energy_cut[i]->GetYaxis()->SetLabelSize(0.045);
	h2_charge1_vs_charge2_energy_cut[i]->GetYaxis()->SetTitleSize(0.045);
	}


TH2D* h2_charge1_vs_charge2_theta_cut[6];
for (Int_t i = 0; i < 6;i++){
	sprintf(hist_name, "Charge1 vs Charge 2 for theta cut: %i < Theta1+Theta2 < %i ",70+i, 90-i);
	h2_charge1_vs_charge2_theta_cut[i] = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
	h2_charge1_vs_charge2_theta_cut[i]->GetXaxis()->SetTitle("Charge 1 ");
	h2_charge1_vs_charge2_theta_cut[i]->GetYaxis()->SetTitle("Charge 2");
	h2_charge1_vs_charge2_theta_cut[i]->GetXaxis()->CenterTitle(true);
	h2_charge1_vs_charge2_theta_cut[i]->GetYaxis()->CenterTitle(true);
	h2_charge1_vs_charge2_theta_cut[i]->GetYaxis()->SetLabelSize(0.045);
	h2_charge1_vs_charge2_theta_cut[i]->GetYaxis()->SetTitleSize(0.045);
	}

TH2D* h2_charge1_vs_charge2_phi_cut[6];
for (Int_t i = 0; i < 6;i++){
	sprintf(hist_name, "Charge1 vs Charge 2 for phi cut: %i < delta_phi < %i ",150+i*5, 210-i*5);
	h2_charge1_vs_charge2_phi_cut[i] = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
	h2_charge1_vs_charge2_phi_cut[i]->GetXaxis()->SetTitle("Charge 1 ");
	h2_charge1_vs_charge2_phi_cut[i]->GetYaxis()->SetTitle("Charge 2");
	h2_charge1_vs_charge2_phi_cut[i]->GetXaxis()->CenterTitle(true);
	h2_charge1_vs_charge2_phi_cut[i]->GetYaxis()->CenterTitle(true);
	h2_charge1_vs_charge2_phi_cut[i]->GetYaxis()->SetLabelSize(0.045);
	h2_charge1_vs_charge2_phi_cut[i]->GetYaxis()->SetTitleSize(0.045);
	}

TH1D* h1_mult_100mev;
sprintf(hist_name, "Multiplicity of events with one hit > 100 MeV ");
h1_mult_100mev = new TH1D(hist_name,hist_name,2500,0,2500);
h1_mult_100mev->GetXaxis()->SetTitle("Multiplicity");
h1_mult_100mev->GetYaxis()->SetTitle("Counts");
h1_mult_100mev->GetXaxis()->CenterTitle(true);
h1_mult_100mev->GetYaxis()->CenterTitle(true);
h1_mult_100mev->GetYaxis()->SetLabelSize(0.045);
h1_mult_100mev->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_mult_100mev_two;
sprintf(hist_name, "Multiplicity of events with two hits > 100 MeV ");
h1_mult_100mev_two = new TH1D(hist_name,hist_name,2500,0,2500);
h1_mult_100mev_two->GetXaxis()->SetTitle("Multiplicity");
h1_mult_100mev_two->GetYaxis()->SetTitle("Counts");
h1_mult_100mev_two->GetXaxis()->CenterTitle(true);
h1_mult_100mev_two->GetYaxis()->CenterTitle(true);
h1_mult_100mev_two->GetYaxis()->SetLabelSize(0.045);
h1_mult_100mev_two->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_theta_sum_p2p;
sprintf(hist_name, "Theta1 +  Theta 2 for the two hits > 100 MeV ");
h1_theta_sum_p2p = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_p2p->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_p2p->GetYaxis()->SetTitle("Counts");
h1_theta_sum_p2p->GetXaxis()->CenterTitle(true);
h1_theta_sum_p2p->GetYaxis()->CenterTitle(true);
h1_theta_sum_p2p->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_p2p->GetYaxis()->SetTitleSize(0.045);


//analysis of theta_sum with Z_sum cuts

TH1D* h1_theta_sum_z_sum_915_925;
sprintf(hist_name, "Theta1 +  Theta 2 for 91.5 < Z_sum < 92.5");
h1_theta_sum_z_sum_915_925 = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_915_925->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_915_925->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_915_925->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_915_925->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_915_925->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_915_925->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_theta_sum_z_sum_915_925_phi_cut;
sprintf(hist_name, "Theta1 +  Theta 2 for 91.5 < Z_sum < 92.5 and delta_phi = 180+-15");
h1_theta_sum_z_sum_915_925_phi_cut = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_915_925_phi_cut->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_915_925_phi_cut->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_915_925_phi_cut->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_915_925_phi_cut->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_915_925_phi_cut->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_915_925_phi_cut->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_theta_sum_z_sum_915_925_phi_ener_cut;
sprintf(hist_name, "Theta1 +  Theta 2 for 91.5 < Z_sum < 92.5 and delta_phi = 180+-15 and  200 < E_sum < 700");
h1_theta_sum_z_sum_915_925_phi_ener_cut = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_915_925_phi_ener_cut->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_915_925_phi_ener_cut->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_915_925_phi_ener_cut->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_915_925_phi_ener_cut->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_915_925_phi_ener_cut->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_915_925_phi_ener_cut->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_theta_sum_z_sum_902_914;
sprintf(hist_name, "Theta1 +  Theta 2 for 90.2 < Z_sum < 91.4");
h1_theta_sum_z_sum_902_914 = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_902_914->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_902_914->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_902_914->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_902_914->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_902_914->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_902_914->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_theta_sum_z_sum_902_914_phi_cut;
sprintf(hist_name, "Theta1 +  Theta 2 for 90.2 < Z_sum < 91.4 and delta_phi = 180+-15");
h1_theta_sum_z_sum_902_914_phi_cut = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_902_914_phi_cut->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_902_914_phi_cut->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_902_914_phi_cut->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_902_914_phi_cut->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_902_914_phi_cut->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_902_914_phi_cut->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_theta_sum_z_sum_902_914_phi_ener_cut;
sprintf(hist_name, "Theta1 +  Theta 2 for 90.2 < Z_sum < 91.4 and delta_phi = 180+-15 and 200 < E_sum < 700");
h1_theta_sum_z_sum_902_914_phi_ener_cut = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_902_914_phi_ener_cut->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_902_914_phi_ener_cut->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_902_914_phi_ener_cut->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_902_914_phi_ener_cut->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_902_914_phi_ener_cut->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_902_914_phi_ener_cut->GetYaxis()->SetTitleSize(0.045);



TH1D* h1_theta_sum_z_sum_891_901;
sprintf(hist_name, "Theta1 +  Theta 2 for 89.1 < Z_sum < 90.1");
h1_theta_sum_z_sum_891_901 = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_891_901->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_891_901->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_891_901->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_891_901->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_891_901->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_891_901->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_theta_sum_z_sum_891_901_phi_cut;
sprintf(hist_name, "Theta1 +  Theta 2 for 89.1 < Z_sum < 90.1 and delta_phi = 180+-15");
h1_theta_sum_z_sum_891_901_phi_cut = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_891_901_phi_cut->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_891_901_phi_cut->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_891_901_phi_cut->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_891_901_phi_cut->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_891_901_phi_cut->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_891_901_phi_cut->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_theta_sum_z_sum_891_901_phi_ener_cut;
sprintf(hist_name, "Theta1 +  Theta 2 for 89.1 < Z_sum < 90.1 and delta_phi = 180+-15 and 200 < E_sum < 700");
h1_theta_sum_z_sum_891_901_phi_ener_cut = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_891_901_phi_ener_cut->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_891_901_phi_ener_cut->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_891_901_phi_ener_cut->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_891_901_phi_ener_cut->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_891_901_phi_ener_cut->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_891_901_phi_ener_cut->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_theta_sum_z_sum_879_89;
sprintf(hist_name, "Theta1 +  Theta 2 for 87.9 < Z_sum < 89");
h1_theta_sum_z_sum_879_89 = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_879_89->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_879_89->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_879_89->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_879_89->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_879_89->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_879_89->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_theta_sum_z_sum_879_89_phi_cut;
sprintf(hist_name, "Theta1 +  Theta 2 for 87.9 < Z_sum < 89 and delta_phi = 180+-15");
h1_theta_sum_z_sum_879_89_phi_cut = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_879_89_phi_cut->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_879_89_phi_cut->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_879_89_phi_cut->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_879_89_phi_cut->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_879_89_phi_cut->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_879_89_phi_cut->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_theta_sum_z_sum_879_89_phi_ener_cut;
sprintf(hist_name, "Theta1 +  Theta 2 for 87.9 < Z_sum < 89 and delta_phi = 180+-15 and 200 < E_sum < 700");
h1_theta_sum_z_sum_879_89_phi_ener_cut = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_879_89_phi_ener_cut->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_879_89_phi_ener_cut->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_879_89_phi_ener_cut->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_879_89_phi_ener_cut->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_879_89_phi_ener_cut->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_879_89_phi_ener_cut->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_theta_sum_z_sum_868_879;
sprintf(hist_name, "Theta1 +  Theta 2 for 86.8 < Z_sum < 87.9");
h1_theta_sum_z_sum_868_879 = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_868_879->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_868_879->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_868_879->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_868_879->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_868_879->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_868_879->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_theta_sum_z_sum_868_879_phi_cut;
sprintf(hist_name, "Theta1 +  Theta 2 for 86.8 < Z_sum < 87.9 and delta_phi = 180+-15");
h1_theta_sum_z_sum_868_879_phi_cut = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_868_879_phi_cut->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_868_879_phi_cut->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_868_879_phi_cut->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_868_879_phi_cut->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_868_879_phi_cut->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_868_879_phi_cut->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_theta_sum_z_sum_868_879_phi_ener_cut;
sprintf(hist_name, "Theta1 +  Theta 2 for 86.8 < Z_sum < 87.9 and delta_phi = 180+-15 and 200 < E_sum < 700");
h1_theta_sum_z_sum_868_879_phi_ener_cut = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_868_879_phi_ener_cut->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_868_879_phi_ener_cut->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_868_879_phi_ener_cut->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_868_879_phi_ener_cut->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_868_879_phi_ener_cut->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_868_879_phi_ener_cut->GetYaxis()->SetTitleSize(0.045);

//-------------------------------------

//Theta1+Theta2 analysis with cut on z1/2

TH1D* h1_theta_sum_z_1_2_cut1;
sprintf(hist_name, "Theta1 +  Theta 2 for 42.7 < Z_1 < 43.9 and 46.3 < Z_2 < 47.6");
h1_theta_sum_z_1_2_cut1 = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_1_2_cut1->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_1_2_cut1->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_1_2_cut1->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_1_2_cut1->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_1_2_cut1->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_1_2_cut1->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_theta_sum_z_1_2_cut2;
sprintf(hist_name, "Theta1 +  Theta 2 for 43.9 < Z_1 < 45.3 and 45.3 < Z_2 < 46.3");
h1_theta_sum_z_1_2_cut2 = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_1_2_cut2->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_1_2_cut2->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_1_2_cut2->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_1_2_cut2->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_1_2_cut2->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_1_2_cut2->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_theta_sum_z_1_2_cut3;
sprintf(hist_name, "Theta1 +  Theta 2 for 42.7 < Z_1 < 43.9 and 45.3 < Z_2 < 46.3");
h1_theta_sum_z_1_2_cut3 = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_1_2_cut3->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_1_2_cut3->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_1_2_cut3->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_1_2_cut3->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_1_2_cut3->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_1_2_cut3->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_theta_sum_z_1_2_cut4;
sprintf(hist_name, "Theta1 +  Theta 2 for 43.9 < Z_1 < 45.3 and 46.3 < Z_2 < 47.6");
h1_theta_sum_z_1_2_cut4 = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_1_2_cut4->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_1_2_cut4->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_1_2_cut4->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_1_2_cut4->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_1_2_cut4->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_1_2_cut4->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_theta_sum_z_1_2_cut5;
sprintf(hist_name, "Theta1 +  Theta 2 for 42.7 < Z_1 < 43.9 and 47.8 < Z_2 < 48.7");
h1_theta_sum_z_1_2_cut5 = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_1_2_cut5->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_1_2_cut5->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_1_2_cut5->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_1_2_cut5->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_1_2_cut5->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_1_2_cut5->GetYaxis()->SetTitleSize(0.045);
//-----------------------------------------------------

////Looking at events with only one hit in TWIM:
TH1D* h1_one_hit_twim_charge;
sprintf(hist_name, "Charge in TWIM when just having one hit");
h1_one_hit_twim_charge = new TH1D(hist_name,hist_name,500,0,100);
h1_one_hit_twim_charge->GetXaxis()->SetTitle("Charge of single Particle hit");
h1_one_hit_twim_charge->GetYaxis()->SetTitle("Counts");
h1_one_hit_twim_charge->GetXaxis()->CenterTitle(true);
h1_one_hit_twim_charge->GetYaxis()->CenterTitle(true);
h1_one_hit_twim_charge->GetYaxis()->SetLabelSize(0.045);
h1_one_hit_twim_charge->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_one_hit_twim_charge_tpat1;
sprintf(hist_name, "Charge in TWIM when just having one hit and tpat == 1");
h1_one_hit_twim_charge_tpat1 = new TH1D(hist_name,hist_name,500,0,100);
h1_one_hit_twim_charge_tpat1->GetXaxis()->SetTitle("Charge of single Particle hit");
h1_one_hit_twim_charge_tpat1->GetYaxis()->SetTitle("Counts");
h1_one_hit_twim_charge_tpat1->GetXaxis()->CenterTitle(true);
h1_one_hit_twim_charge_tpat1->GetYaxis()->CenterTitle(true);
h1_one_hit_twim_charge_tpat1->GetYaxis()->SetLabelSize(0.045);
h1_one_hit_twim_charge_tpat1->GetYaxis()->SetTitleSize(0.045);


//---------------------------------------





TH2D* h2_theta1_vs_theta2;
sprintf(hist_name, "Theta 1 vs Theta 2 for the two hits > 100 MeV");
h2_theta1_vs_theta2 = new TH2D(hist_name,hist_name,25,22.15,84.65,25,22.15,84.65);
h2_theta1_vs_theta2->GetXaxis()->SetTitle("Theta1 [degr]");
h2_theta1_vs_theta2->GetYaxis()->SetTitle("Theta2 [degr]");
h2_theta1_vs_theta2->GetXaxis()->CenterTitle(true);
h2_theta1_vs_theta2->GetYaxis()->CenterTitle(true);
h2_theta1_vs_theta2->GetYaxis()->SetLabelSize(0.045);
h2_theta1_vs_theta2->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_theta1_vs_theta2_cluster;
sprintf(hist_name, "Theta 1 vs Theta 2 for the two hits with highest energy");
h2_theta1_vs_theta2_cluster = new TH2D(hist_name,hist_name,70,20,90,70,20,90);
h2_theta1_vs_theta2_cluster->GetXaxis()->SetTitle("Theta1 [degr]");
h2_theta1_vs_theta2_cluster->GetYaxis()->SetTitle("Theta2 [degr]");
h2_theta1_vs_theta2_cluster->GetXaxis()->CenterTitle(true);
h2_theta1_vs_theta2_cluster->GetYaxis()->CenterTitle(true);
h2_theta1_vs_theta2_cluster->GetYaxis()->SetLabelSize(0.045);
h2_theta1_vs_theta2_cluster->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_theta1_vs_theta2_cluster_best_dphi;
sprintf(hist_name, "Theta 1 vs Theta 2 for the two hits with highest energy (> 30MeV) and best dphi");
h2_theta1_vs_theta2_cluster_best_dphi = new TH2D(hist_name,hist_name,70,20,90,70,20,90);
h2_theta1_vs_theta2_cluster_best_dphi->GetXaxis()->SetTitle("Theta1 [degr]");
h2_theta1_vs_theta2_cluster_best_dphi->GetYaxis()->SetTitle("Theta2 [degr]");
h2_theta1_vs_theta2_cluster_best_dphi->GetXaxis()->CenterTitle(true);
h2_theta1_vs_theta2_cluster_best_dphi->GetYaxis()->CenterTitle(true);
h2_theta1_vs_theta2_cluster_best_dphi->GetYaxis()->SetLabelSize(0.045);
h2_theta1_vs_theta2_cluster_best_dphi->GetYaxis()->SetTitleSize(0.045);


TH2D* h2_e1_vs_e2_cluster;
sprintf(hist_name, "Energy 1 vs Energy 2 for the two hits with highest energy");
h2_e1_vs_e2_cluster = new TH2D(hist_name,hist_name,500,0,1000,500,0,1000);
h2_e1_vs_e2_cluster->GetXaxis()->SetTitle("Energy1 [MeV]");
h2_e1_vs_e2_cluster->GetYaxis()->SetTitle("Energy2 [MeV]");
h2_e1_vs_e2_cluster->GetXaxis()->CenterTitle(true);
h2_e1_vs_e2_cluster->GetYaxis()->CenterTitle(true);
h2_e1_vs_e2_cluster->GetYaxis()->SetLabelSize(0.045);
h2_e1_vs_e2_cluster->GetYaxis()->SetTitleSize(0.045);


TH2D* h2_e1_vs_e2_cluster_best_dphi;
sprintf(hist_name, "Energy 1 vs Energy 2 for the two hits with highest energy (> 30MeV) and best dphi");
h2_e1_vs_e2_cluster_best_dphi = new TH2D(hist_name,hist_name,500,0,1000,500,0,1000);
h2_e1_vs_e2_cluster_best_dphi->GetXaxis()->SetTitle("Energy1 [MeV]");
h2_e1_vs_e2_cluster_best_dphi->GetYaxis()->SetTitle("Energy2 [MeV]");
h2_e1_vs_e2_cluster_best_dphi->GetXaxis()->CenterTitle(true);
h2_e1_vs_e2_cluster_best_dphi->GetYaxis()->CenterTitle(true);
h2_e1_vs_e2_cluster_best_dphi->GetYaxis()->SetLabelSize(0.045);
h2_e1_vs_e2_cluster_best_dphi->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_e1_vs_e2_cluster_two_highest;
sprintf(hist_name, "Energy 1 vs Energy 2 for the two hits with highest energy (> 30MeV) and dphi = 180 +-15");
h2_e1_vs_e2_cluster_two_highest = new TH2D(hist_name,hist_name,500,0,1000,500,0,1000);
h2_e1_vs_e2_cluster_two_highest->GetXaxis()->SetTitle("Energy1 [MeV]");
h2_e1_vs_e2_cluster_two_highest->GetYaxis()->SetTitle("Energy2 [MeV]");
h2_e1_vs_e2_cluster_two_highest->GetXaxis()->CenterTitle(true);
h2_e1_vs_e2_cluster_two_highest->GetYaxis()->CenterTitle(true);
h2_e1_vs_e2_cluster_two_highest->GetYaxis()->SetLabelSize(0.045);
h2_e1_vs_e2_cluster_two_highest->GetYaxis()->SetTitleSize(0.045);




TH1D* h1_theta_sum_p2p_best_dphi;
sprintf(hist_name, "Theta1 +  Theta 2 for two hits (> 30 MeV) and best dphi");
h1_theta_sum_p2p_best_dphi = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_p2p_best_dphi->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_p2p_best_dphi->GetYaxis()->SetTitle("Counts");
h1_theta_sum_p2p_best_dphi->GetXaxis()->CenterTitle(true);
h1_theta_sum_p2p_best_dphi->GetYaxis()->CenterTitle(true);
h1_theta_sum_p2p_best_dphi->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_p2p_best_dphi->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_opening_angle;
sprintf(hist_name, "Opening Angle for two hits (> 30 MeV) and best dphi");
h1_opening_angle = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_opening_angle->GetXaxis()->SetTitle("Opening Angle [degr]");
h1_opening_angle->GetYaxis()->SetTitle("Counts");
h1_opening_angle->GetXaxis()->CenterTitle(true);
h1_opening_angle->GetYaxis()->CenterTitle(true);
h1_opening_angle->GetYaxis()->SetLabelSize(0.045);
h1_opening_angle->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_opening_angle_first_two;
sprintf(hist_name, "Opening Angle for the first two high energetic hits (> 30 MeV) and best dphi");
h1_opening_angle_first_two = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_opening_angle_first_two->GetXaxis()->SetTitle("Opening Angle [degr]");
h1_opening_angle_first_two->GetYaxis()->SetTitle("Counts");
h1_opening_angle_first_two->GetXaxis()->CenterTitle(true);
h1_opening_angle_first_two->GetYaxis()->CenterTitle(true);
h1_opening_angle_first_two->GetYaxis()->SetLabelSize(0.045);
h1_opening_angle_first_two->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_opening_angle_one_other;
sprintf(hist_name, "Opening Angle for the first and other (not second) energetic hits (> 30 MeV) and best dphi");
h1_opening_angle_one_other = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_opening_angle_one_other->GetXaxis()->SetTitle("Opening Angle [degr]");
h1_opening_angle_one_other->GetYaxis()->SetTitle("Counts");
h1_opening_angle_one_other->GetXaxis()->CenterTitle(true);
h1_opening_angle_one_other->GetYaxis()->CenterTitle(true);
h1_opening_angle_one_other->GetYaxis()->SetLabelSize(0.045);
h1_opening_angle_one_other->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_opening_angle_two_other;
sprintf(hist_name, "Opening Angle for the second and other (not first) energetic hits (> 30 MeV) and best dphi");
h1_opening_angle_two_other = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_opening_angle_two_other->GetXaxis()->SetTitle("Opening Angle [degr]");
h1_opening_angle_two_other->GetYaxis()->SetTitle("Counts");
h1_opening_angle_two_other->GetXaxis()->CenterTitle(true);
h1_opening_angle_two_other->GetYaxis()->CenterTitle(true);
h1_opening_angle_two_other->GetYaxis()->SetLabelSize(0.045);
h1_opening_angle_two_other->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_opening_angle_proton;
sprintf(hist_name, "Opening Angle for two hits (> 30 MeV) and best dphi only in proton range");
h1_opening_angle_proton = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_opening_angle_proton->GetXaxis()->SetTitle("Opening Angle [degr]");
h1_opening_angle_proton->GetYaxis()->SetTitle("Counts");
h1_opening_angle_proton->GetXaxis()->CenterTitle(true);
h1_opening_angle_proton->GetYaxis()->CenterTitle(true);
h1_opening_angle_proton->GetYaxis()->SetLabelSize(0.045);
h1_opening_angle_proton->GetYaxis()->SetTitleSize(0.045);

TH2D* h2_opening_angle_vs_esum;
sprintf(hist_name, "Opening Angle for two hits (> 30 MeV) and best dphi vs E1+E2");
h2_opening_angle_vs_esum = new TH2D(hist_name,hist_name,52,22.15,152.15,100,0,1000);
h2_opening_angle_vs_esum->GetXaxis()->SetTitle("Opening Angle [degr]");
h2_opening_angle_vs_esum->GetYaxis()->SetTitle("E1+E2 [MeV]");
h2_opening_angle_vs_esum->GetXaxis()->CenterTitle(true);
h2_opening_angle_vs_esum->GetYaxis()->CenterTitle(true);
h2_opening_angle_vs_esum->GetYaxis()->SetLabelSize(0.045);
h2_opening_angle_vs_esum->GetYaxis()->SetTitleSize(0.045);


TH2D* h2_e_vs_phi;
sprintf(hist_name, "E1 or E2 vs angular distribution");
h2_e_vs_phi = new TH2D(hist_name,hist_name,52,22.15,152.15,60,-180,180);
h2_e_vs_phi->GetXaxis()->SetTitle("Theta [degr]");
h2_e_vs_phi->GetYaxis()->SetTitle("Phi [degr]");
h2_e_vs_phi->GetXaxis()->CenterTitle(true);
h2_e_vs_phi->GetYaxis()->CenterTitle(true);
h2_e_vs_phi->GetYaxis()->SetLabelSize(0.045);
h2_e_vs_phi->GetYaxis()->SetTitleSize(0.045);


TH1D* h1_angle_check;
sprintf(hist_name, "angle check ");
h1_angle_check = new TH1D(hist_name,hist_name,3500,0,3.5);
h1_angle_check->GetXaxis()->SetTitle("Polar angle");
h1_angle_check->GetYaxis()->SetTitle("Counts");
h1_angle_check->GetXaxis()->CenterTitle(true);
h1_angle_check->GetYaxis()->CenterTitle(true);
h1_angle_check->GetYaxis()->SetLabelSize(0.045);
h1_angle_check->GetYaxis()->SetTitleSize(0.045);

TH1D* h1_angle_check_arzi;
sprintf(hist_name, "angle check arzimuthal ");
h1_angle_check_arzi = new TH1D(hist_name,hist_name,7000,-3.5,3.5);
h1_angle_check_arzi->GetXaxis()->SetTitle("Arzimuthal angle");
h1_angle_check_arzi->GetYaxis()->SetTitle("Counts");
h1_angle_check_arzi->GetXaxis()->CenterTitle(true);
h1_angle_check_arzi->GetYaxis()->CenterTitle(true);
h1_angle_check_arzi->GetYaxis()->SetLabelSize(0.045);
h1_angle_check_arzi->GetYaxis()->SetTitleSize(0.045);

//Histo with trigger info:
TH1F* h1_trigger = new TH1F("h1_trigger", "Trigger information: Tpat", 17, -0.5, 16.5);
h1_trigger->GetXaxis()->SetTitle("Trigger number (tpat)");
h1_trigger->GetYaxis()->SetTitle("Counts");
h1_trigger->GetXaxis()->CenterTitle(true);
h1_trigger->GetYaxis()->CenterTitle(true);
h1_trigger->GetXaxis()->SetLabelSize(0.04);
h1_trigger->GetXaxis()->SetTitleSize(0.04);
h1_trigger->GetYaxis()->SetTitleOffset(1.1);
h1_trigger->GetXaxis()->SetTitleOffset(1.1);
h1_trigger->GetYaxis()->SetLabelSize(0.04);
h1_trigger->GetYaxis()->SetTitleSize(0.04);
h1_trigger->SetFillColor(kBlue + 2);

fstream fin;
fin.open("../parameters/crystal_mapping_s455_03.csv", ios::in); //FIXME:change to right .csv file

//new functionality, using map
map<int,int> map_id_califa_side;
string line, word, temp;
while(fin.peek()!=EOF) {
	getline(fin, line);
	stringstream s(line);
	vector<string> temp_vec;
	while (getline(s, word, ',')) {
		temp_vec.push_back(word);
	}
	map_id_califa_side[stoi(temp_vec[10])] = stoi(temp_vec[0]);
	temp_vec.clear();
}


char fname[500];
Long64_t entries_califa = 0;
Long64_t entries_califa_hit = 0;
Long64_t entries_tof = 0;
Long64_t entries_twim = 0;
Int_t entries_wr = 0;
//Input file
sprintf(fname,"../file_src/s455_03_273_10_sq_w_cl_antia.root");
TChain* chain = new TChain("evt");
chain->Reset();
chain->Add(fname);
Long64_t nevents = chain->GetEntries();
cout << "total number of entries:\t" << nevents << endl;

//TClonesarrays

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
//
const Double_t PI = 3.14159265358979323846;
uint64_t wr_Master_ts;
Long64_t events_califa_no_master = 0;
Long64_t events_califa_with_master = 0;
Long64_t events_total_califa = 0;
Long64_t events_with_cut = 0;
Long64_t events_one_proton = 0;
Long64_t events_two_proton = 0;
vector<vector<double> > v_califa_energy_theta;
vector<vector<double> > v_califa_energy_theta_phi;
for(Long64_t i=0;i< nevents;i++){
    Long64_t evtnr = i;
    Int_t trigger_pattern = 0;
    if (i%100000==0)
        cout<<"Processing event for charge analysis "<<i<<endl;
    chain->GetEvent(i);

    Int_t tpatbin = 0;
    if (DataCA->GetTpat() > 0){
            for (Int_t i = 0; i < 16; i++){
        	    tpatbin = (DataCA->GetTpat() & (1 << i));
        	    if (tpatbin != 0){
			    trigger_pattern = i + 1;
        		    h1_trigger->Fill(trigger_pattern);
        	    }
            }
    }
    entries_califa = CalifaMappedData->GetEntries();
    entries_califa_hit = CalifaHitData->GetEntries();
    entries_wr = WRMasterData->GetEntries();
    entries_tof = SofToFWHitData->GetEntries();
    entries_twim = SofTwimHitData->GetEntries();
	if (entries_twim == 1){
		softwimhitdata = new R3BSofTwimHitData*[1];
		softwimhitdata[0] = (R3BSofTwimHitData*)SofTwimHitData->At(0);
		Double_t charge_one_particle = softwimhitdata[0]->GetZcharge();
		h1_one_hit_twim_charge->Fill(charge_one_particle);
		if ((DataCA->GetTpat() & (1 << 0))==1){
			h1_one_hit_twim_charge_tpat1->Fill(charge_one_particle);
			}
		delete [] softwimhitdata;
		}
	//before I had: if (entries_califa >= 1 && entries_tof >= 2 && entries_tof <= 4 && entries_twim == 2){
	if (entries_califa >= 1 && entries_twim == 2 && entries_tof==2){
		//check twim data
		softwimhitdata = new R3BSofTwimHitData*[2];
		softofwhitdata = new R3BSofTofWHitData*[2];

		//insert logic from TOFW online...
		Double_t tof[2] = { 0., 0. };
		int padid[2] = { 0, 0 };
		for (Int_t ihit = 0; ihit < entries_tof; ihit++){
			softofwhitdata[ihit] = (R3BSofTofWHitData*)SofToFWHitData->At(ihit);
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
		for (Int_t ihit = 0.; ihit < entries_twim; ihit++){
			softwimhitdata[ihit] = (R3BSofTwimHitData*)SofTwimHitData->At(ihit);
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
		h2_charge1_charge_2->Fill(charge_1,charge_2);
		h2_charge1_charge_2_diff_vs_z_sum->Fill(charge_1-charge_2,charge_1+charge_2);
		h1_charge1_plus_charge2->Fill(charge_1+charge_2);
		h1_daughter_nuclei_charge->Fill(charge_1);
		h1_daughter_nuclei_charge->Fill(charge_2);

		//Check Z-pattern for different tpats...
		if (trigger_pattern == 1 || trigger_pattern == 7){
			h2_charge1_charge_2_tpat1->Fill(charge_1,charge_2);
			h2_charge1_charge_2_diff_vs_z_sum_tpat1->Fill(charge_1-charge_2,charge_1+charge_2);
			h1_charge1_plus_charge2_tpat1->Fill(charge_1+charge_2);
			}
		if (trigger_pattern == 2 || trigger_pattern == 8){
			h2_charge1_charge_2_tpat2->Fill(charge_1,charge_2);
			h2_charge1_charge_2_diff_vs_z_sum_tpat2->Fill(charge_1-charge_2,charge_1+charge_2);
			h1_charge1_plus_charge2_tpat2->Fill(charge_1+charge_2);
			}
		if (trigger_pattern == 3 || trigger_pattern == 9){
			h2_charge1_charge_2_tpat3->Fill(charge_1,charge_2);
			h2_charge1_charge_2_diff_vs_z_sum_tpat3->Fill(charge_1-charge_2,charge_1+charge_2);
			h1_charge1_plus_charge2_tpat3->Fill(charge_1+charge_2);
			}
		if (trigger_pattern == 4 || trigger_pattern == 10){
			h2_charge1_charge_2_tpat4->Fill(charge_1,charge_2);
			h2_charge1_charge_2_diff_vs_z_sum_tpat4->Fill(charge_1-charge_2,charge_1+charge_2);
			h1_charge1_plus_charge2_tpat4->Fill(charge_1+charge_2);
			}

		//cut  on charge of twim
		//if (30 < charge_1 && charge_1 < 60 && 30 < charge_2 && charge_2 < 60 && (charge_1+charge_2) > 89 && (charge_1+charge_2) < 93 ){
		if ((charge_1+charge_2) < 100){
	
		events_total_califa +=1;
		califamappeddata = new R3BCalifaMappedData*[entries_califa];
		wrmasterdata = new R3BWRMasterData*[entries_wr];
		vector<uint64_t> califa_messel_ts;
		vector<uint64_t> califa_wixhausen_ts;
		vector<uint64_t> califa_both_sides_ts;
		int e_high = 0;
		double e_sum_with_M = 0;
		double e_sum_no_M = 0;
		if (entries_wr == 1){
			events_califa_with_master+=1;
			h1_mult_with_M->Fill(entries_califa);
			if (entries_califa_hit){
				califahitdata = new R3BCalifaHitData*[entries_califa_hit];
			Int_t califa_100mev_hits = 0;
			Int_t califa_200mev_hits = 0;
			for (Int_t m = 0; m < entries_califa_hit; m++){
				vector<double> v_temp_califa(2);
				califahitdata[m] = (R3BCalifaHitData*)CalifaHitData->At(m);
				h2_energy_vs_theta->Fill((califahitdata[m]->GetEnergy())/1000,((califahitdata[m]->GetTheta())/PI)*180);
				h2_phi_vs_theta->Fill(((califahitdata[m]->GetTheta())/PI)*180,((califahitdata[m]->GetPhi())/PI)*180);
				h1_angle_check->Fill(califahitdata[m]->GetTheta());
				h1_angle_check_arzi->Fill(califahitdata[m]->GetPhi());
				v_temp_califa[0] = (califahitdata[m]->GetEnergy())/1000.;
				//randomize angles to remove artifacts...
				if (califahitdata[m]->GetTheta() > 0.45){
					Int_t rand_angl = rand() % 4500 + 1;
					v_temp_califa[1] = ((califahitdata[m]->GetTheta() - 0.0225 + 0.00001*rand_angl)/PI)*180.;
					}
				else {
					v_temp_califa[1] = ((califahitdata[m]->GetTheta())/PI)*180.;
					}
				////v_temp_califa[1] =  ((califahitdata[m]->GetTheta())/PI)*180.;
				v_califa_energy_theta.push_back(v_temp_califa);
				v_temp_califa.clear();
				if (((califahitdata[m]->GetEnergy())/1000) > 100){
					califa_100mev_hits++;
				}
				if(((califahitdata[m]->GetEnergy())/1000) > 200){
					califa_200mev_hits++;
					}
				//section to fill 2d vector for califa hits > 30 mev and choose the best ones with diff_phi = 180+-30-------------------
				vector<double> v_temp_califa_ext(3);
				if (((califahitdata[m]->GetEnergy())/1000) > 30){
					v_temp_califa_ext[0] = (califahitdata[m]->GetEnergy())/1000.;
					
					if (califahitdata[m]->GetTheta() > 0.45){
						Int_t rand_angl = rand() % 4500 + 1;
						v_temp_califa_ext[1] = ((califahitdata[m]->GetTheta() - 0.0225 + 0.00001*rand_angl)/PI)*180.;
						}
					else {
						v_temp_califa_ext[1] = ((califahitdata[m]->GetTheta())/PI)*180.;
						}
					
					v_temp_califa_ext[2] = ((califahitdata[m]->GetPhi())/PI)*180;
					v_califa_energy_theta_phi.push_back(v_temp_califa_ext);
					v_temp_califa_ext.clear();
					}

				//-----------------------------------------------------------------------------------------------------------------------
				//fill 3D plot
				if(events_with_cut < 20){
					for (Int_t j = 0; j < round((califahitdata[m]->GetEnergy()));j++){
				h2_energy_angular_dist[events_with_cut]->Fill(((califahitdata[m]->GetTheta())/PI)*180,((califahitdata[m]->GetPhi())/PI)*180);
					}
				}
				}
			vector<double> v_special;
			//here I select the best califa hits with diff_phi = 180+-30---------
			//------Algorithm begin------------" << endl;
			if (v_califa_energy_theta_phi.size() > 2){
				sort(v_califa_energy_theta_phi.begin(),v_califa_energy_theta_phi.end(),sortcol);
				vector<vector<double> > v_phi_diff;
				for(Int_t v1 = 0; v1 < 2; v1++){
					for (Int_t v2 = v1+1; v2 < v_califa_energy_theta_phi.size();v2++){
						vector<double> v_temp_phi_diff(3);
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


						v_phi_diff.push_back(v_temp_phi_diff);
						v_temp_phi_diff.clear();
						}
					}
				Double_t phi_diff_one_two = v_phi_diff[0][2];
				//-----------Algorithm end -------------------------------
				
				//select the two hits with highest energy if phi_diff = 180+-30
				if ( phi_diff_one_two < 30){
					v_special.push_back(v_califa_energy_theta_phi[0][0]+v_califa_energy_theta_phi[1][0]);
					v_special.push_back(v_califa_energy_theta_phi[0][1]+v_califa_energy_theta_phi[1][1]);
					v_special.push_back(phi_diff_one_two);
					int rand_num;
					rand_num = rand() % 2;
					h2_charge1_charge_2_after->Fill(charge_1,charge_2);
					h1_charge_sum_after->Fill(charge_1+charge_2);
					//if ((v_califa_energy_theta_phi[0][0]+v_califa_energy_theta_phi[1][0]) < 650 && (v_califa_energy_theta_phi[0][0]+v_califa_energy_theta_phi[1][0]) > 300 && (v_califa_energy_theta_phi[0][1]+v_califa_energy_theta_phi[1][1]) < 90 && (v_califa_energy_theta_phi[0][1]+v_califa_energy_theta_phi[1][1]) > 70 ){
					//	h2_charge1_charge_2_after->Fill(charge_1,charge_2);	
					//	h1_charge_sum_after->Fill(charge_1+charge_2);
					//	}
					h2_theta1_vs_theta2_cluster_best_dphi->Fill(v_califa_energy_theta_phi[rand_num][1],v_califa_energy_theta_phi[1-rand_num][1]);
					h2_e1_vs_e2_cluster_best_dphi->Fill(v_califa_energy_theta_phi[rand_num][0],v_califa_energy_theta_phi[1-rand_num][0]);
					h1_theta_sum_p2p_best_dphi->Fill(v_califa_energy_theta_phi[0][1]+v_califa_energy_theta_phi[1][1]);
					h2_e1_vs_e2_cluster_two_highest->Fill(v_califa_energy_theta_phi[rand_num][0],v_califa_energy_theta_phi[1-rand_num][0]);
					h2_e_vs_phi->Fill(v_califa_energy_theta_phi[0][1],v_califa_energy_theta_phi[0][2]);
					h2_e_vs_phi->Fill(v_califa_energy_theta_phi[1][1],v_califa_energy_theta_phi[1][2]);
					TVector3 sel_vector_1 (v_califa_energy_theta_phi[0][0]*sin(v_califa_energy_theta_phi[0][1]*PI/180.)*cos(v_califa_energy_theta_phi[0][2]*PI/180.),v_califa_energy_theta_phi[0][0]*sin(v_califa_energy_theta_phi[0][1]*PI/180.)*sin(v_califa_energy_theta_phi[0][2]*PI/180.),v_califa_energy_theta_phi[0][0]*cos(v_califa_energy_theta_phi[0][1]*PI/180.));

					TVector3 sel_vector_2 (v_califa_energy_theta_phi[1][0]*sin(v_califa_energy_theta_phi[1][1]*PI/180.)*cos(v_califa_energy_theta_phi[1][2]*PI/180.),v_califa_energy_theta_phi[1][0]*   sin(v_califa_energy_theta_phi[1][1]*PI/180.)*sin(v_califa_energy_theta_phi[1][2]*PI/180.),v_califa_energy_theta_phi[1][0]*cos(v_califa_energy_theta_phi[1][1]*PI/180.));
					Double_t angle_two_vecs = sel_vector_1.Angle(sel_vector_2);
					h1_opening_angle->Fill(angle_two_vecs*180/PI);
					h1_opening_angle_first_two->Fill(angle_two_vecs*180/PI);
					if (angle_two_vecs*180/PI < 60){
					cout << "event to study:\t" << evtnr << endl;
					}
					//now I select only hits in the PROTON RANGE region
					if( v_califa_energy_theta_phi[0][1] < 43){
						if (v_califa_energy_theta_phi[1][1] < 43){
							h1_opening_angle_proton->Fill(angle_two_vecs*180/PI);
							}
						else {
							if (((v_califa_energy_theta_phi[1][2] > -22.5 && v_califa_energy_theta_phi[1][2] < 22.5) || v_califa_energy_theta_phi[1][2] > 157 || v_califa_energy_theta_phi[1][2] < -   157)){
								h1_opening_angle_proton->Fill(angle_two_vecs*180/PI);
								}
								}
								}
					else{
						if(v_califa_energy_theta_phi[1][1] < 43){
							if(((v_califa_energy_theta_phi[0][2] > -22.5 && v_califa_energy_theta_phi[0][2] < 22.5) || v_califa_energy_theta_phi[0][2] > 157 || v_califa_energy_theta_phi[0][2] < -   157)){
								h1_opening_angle_proton->Fill(angle_two_vecs*180/PI);
							}
							}
						else if(((v_califa_energy_theta_phi[0][2] > -22.5 && v_califa_energy_theta_phi[0][2] < 22.5) || v_califa_energy_theta_phi[0][2] > 157 || v_califa_energy_theta_phi[0][2] < -   157) && ((v_califa_energy_theta_phi[1][2] > -22.5 && v_califa_energy_theta_phi[1][2] < 22.5) || v_califa_energy_theta_phi[1][2] > 157 || v_califa_energy_theta_phi[1][2] < -   157)){
							h1_opening_angle_proton->Fill(angle_two_vecs*180/PI);
							}
					}

					//----END proton range-------------------------------------
						

					if (angle_two_vecs*180/PI < 30){
						cout << "this is the eventnumber with small opening angle:" << evtnr << endl;
						cout << "magnitude of first vector:\t" << sel_vector_1.Mag() << endl;;
						cout << "magnitude of second vector:\t" << sel_vector_2.Mag() << endl;
						cout << "Energy one of v_califa_energy_theta_phi:\t" << v_califa_energy_theta_phi[0][0] << endl;
						cout << "Energy two of v_califa_energy_theta_phi:\t" << v_califa_energy_theta_phi[1][0] << endl;
						cout << "and the according opening angle:\t" << angle_two_vecs*180/PI << endl;
						}
					h2_opening_angle_vs_esum->Fill(angle_two_vecs*180/PI,v_califa_energy_theta_phi[0][0]+v_califa_energy_theta_phi[1][0]);

					}
				sort(v_phi_diff.begin(),v_phi_diff.end(),sort_min);

				//the hits with highest energy do not have phi_diff = 180+-, therefore use result from algorithm above ....
				if (phi_diff_one_two > 30 &&  v_phi_diff[0][2] < 30){
					v_special.push_back(v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][0]+v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][0]);
					v_special.push_back(v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][1]+v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][1]);
					v_special.push_back(v_phi_diff[0][2]);
					h2_charge1_charge_2_after->Fill(charge_1,charge_2);
					h1_charge_sum_after->Fill(charge_1+charge_2);
					//if ((v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][0]+v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][0]) < 650 && (v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][0]+v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][0]) > 300 &&                  (v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][1]+v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][1]) < 90 && (v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][1]+v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][1]) > 70 ){
					//h2_charge1_charge_2_after->Fill(charge_1,charge_2);	
					//h1_charge_sum_after->Fill(charge_1+charge_2);
					//}
					TVector3 sel_vector_1 (v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][0]*sin(v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][1]*PI/180.)*cos(v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][2]*PI/180.),v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][0]*sin(v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][1]*PI/180.)*sin(v_califa_energy_theta_phi[(int)   v_phi_diff[0][0]][2]*PI/180.),v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][0]*cos(v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][1]*PI/180.));

					TVector3 sel_vector_2 (v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][0]*sin(v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][1]*PI/180.)*cos(v_califa_energy_theta_phi[(int)   v_phi_diff[0][1]][2]*PI/180.),v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][0]*sin(v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][1]*PI/180.)*sin(v_califa_energy_theta_phi[(int)   v_phi_diff[0][1]][2]*PI/180.),v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][0]*cos(v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][1]*PI/180.));
					
					h2_e_vs_phi->Fill(v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][1],v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][2]);
					h2_e_vs_phi->Fill(v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][1],v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][2]);
					
					Double_t angle_two_vecs = sel_vector_1.Angle(sel_vector_2);
					h1_opening_angle->Fill(angle_two_vecs*180/PI);
					
					if (v_phi_diff[0][0] == 0 ){
						h1_opening_angle_one_other->Fill(angle_two_vecs*180/PI);
						}
					if (v_phi_diff[0][0] == 1){
						h1_opening_angle_two_other->Fill(angle_two_vecs*180/PI);
						}
					
					//now I select only hits in the PROTON RANGE region
					if (v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][1] < 43){
						if (v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][1] < 43){
							h1_opening_angle_proton->Fill(angle_two_vecs*180/PI);
							}
						else{
							if(((v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][2] > -22.5 && v_califa_energy_theta_phi[(int)                    v_phi_diff[0][1]][2] < 22.5) || v_califa_energy_theta_phi[(int)       v_phi_diff[0][1]][2] > 157 || v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][2] < -157)){
								h1_opening_angle_proton->Fill(angle_two_vecs*180/PI);
								}
							}
						}
					else {
						if (v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][1] < 43){
							if (((v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][2] > -22.5 && v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][2] < 22.5) || v_califa_energy_theta_phi[(int)       v_phi_diff[0][0]][2] > 157 || v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][2] < -157)){
								h1_opening_angle_proton->Fill(angle_two_vecs*180/PI);
								}

							}
						else if (((v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][2] > -22.5 && v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][2] < 22.5) || v_califa_energy_theta_phi[(int)       v_phi_diff[0][0]][2] > 157 || v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][2] < -157) && ((v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][2] > -22.5 && v_califa_energy_theta_phi[(int)                    v_phi_diff[0][1]][2] < 22.5) || v_califa_energy_theta_phi[(int)       v_phi_diff[0][1]][2] > 157 || v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][2] < -157)){
							h1_opening_angle_proton->Fill(angle_two_vecs*180/PI);
							}
						}




					//----END of proton region-------------------------

					if (angle_two_vecs*180/PI < 30){
						cout << "this is the eventnumber with small opening angle:" << evtnr << endl;
						cout << "and the according opening angle:\t" << angle_two_vecs*180/PI << endl;
						}
					h2_opening_angle_vs_esum->Fill(angle_two_vecs*180/PI,v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][0]+v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][0]);

					//fill histo
					int rand_num;
					rand_num = rand() % 2;
					h2_theta1_vs_theta2_cluster_best_dphi->Fill(v_califa_energy_theta_phi[(int) v_phi_diff[0][rand_num]][1],v_califa_energy_theta_phi[(int) v_phi_diff[0][1-rand_num]][1]);
					h2_e1_vs_e2_cluster_best_dphi->Fill(v_califa_energy_theta_phi[(int) v_phi_diff[0][rand_num]][0],v_califa_energy_theta_phi[(int) v_phi_diff[0][1-rand_num]][0]);
					h1_theta_sum_p2p_best_dphi->Fill(v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][1]+v_califa_energy_theta_phi[(int) v_phi_diff[0][1]][1]);
					//cout << "Chosen energies: \t" << v_califa_energy_theta_phi[(int) v_phi_diff[0][0]][0] << " , " << 
					v_califa_energy_theta_phi.clear();
					}	
				else{
					v_califa_energy_theta_phi.clear();
					//continue;
				}
				}
			//if we have exactly two hits with E > 30 MeV
			else if (v_califa_energy_theta_phi.size() == 2){
				sort(v_califa_energy_theta_phi.begin(),v_califa_energy_theta_phi.end(),sortcol);
				
				Double_t abs_phi_diff;

				if (v_califa_energy_theta_phi[0][2] > v_califa_energy_theta_phi[1][2]){
					Double_t angular_value1 = 360 + v_califa_energy_theta_phi[1][2] - v_califa_energy_theta_phi[0][2];
					Double_t angular_value2 = v_califa_energy_theta_phi[0][2] - v_califa_energy_theta_phi[1][2];
						if (angular_value1 < angular_value2){
							 abs_phi_diff = 180 - angular_value1;
							}
						else {
							abs_phi_diff = 180 - angular_value2;
							}
					}
				if ( v_califa_energy_theta_phi[1][2] > v_califa_energy_theta_phi[0][2]){
					Double_t angular_value1 = 360 + v_califa_energy_theta_phi[0][2] - v_califa_energy_theta_phi[1][2];
					Double_t angular_value2 = v_califa_energy_theta_phi[1][2] - v_califa_energy_theta_phi[0][2];	
						if (angular_value1 < angular_value2){
							abs_phi_diff = 180 - angular_value1;
							}
						else {
							abs_phi_diff = 180 - angular_value2;
							}

					}

				if (abs_phi_diff < 30){
					
					v_special.push_back(v_califa_energy_theta_phi[0][0]+v_califa_energy_theta_phi[1][0]);
					v_special.push_back(v_califa_energy_theta_phi[0][1]+v_califa_energy_theta_phi[1][1]);
					v_special.push_back(abs_phi_diff);
					//fill histo
					int rand_num;
					rand_num = rand() % 2;
					h2_charge1_charge_2_after->Fill(charge_1,charge_2);
					h1_charge_sum_after->Fill(charge_1+charge_2);
					//if ((v_califa_energy_theta_phi[0][0]+v_califa_energy_theta_phi[1][0]) < 650 && (v_califa_energy_theta_phi[0][0]+v_califa_energy_theta_phi[1][0]) > 300 && (v_califa_energy_theta_phi[0][1]+v_califa_energy_theta_phi[1][1]) < 90 && (v_califa_energy_theta_phi[0][1]+v_califa_energy_theta_phi[1][1]) > 70 ){
					//	h2_charge1_charge_2_after->Fill(charge_1,charge_2);	
					//	h1_charge_sum_after->Fill(charge_1+charge_2);
					//	}
					h2_theta1_vs_theta2_cluster_best_dphi->Fill(v_califa_energy_theta_phi[rand_num][1],v_califa_energy_theta_phi[1-rand_num][1]);
					h2_e1_vs_e2_cluster_best_dphi->Fill(v_califa_energy_theta_phi[rand_num][0],v_califa_energy_theta_phi[1-rand_num][0]);
					h1_theta_sum_p2p_best_dphi->Fill(v_califa_energy_theta_phi[0][1]+v_califa_energy_theta_phi[1][1]);
					h2_e_vs_phi->Fill(v_califa_energy_theta_phi[0][1],v_califa_energy_theta_phi[0][2]);
					h2_e_vs_phi->Fill(v_califa_energy_theta_phi[1][1],v_califa_energy_theta_phi[1][2]);

					
					TVector3 sel_vector_1 (v_califa_energy_theta_phi[0][0]*sin(v_califa_energy_theta_phi[0][1]*PI/180.)*cos(v_califa_energy_theta_phi[0][2]*PI/180.),v_califa_energy_theta_phi[0][0]*sin(v_califa_energy_theta_phi[0][1]*PI/180.)*sin(v_califa_energy_theta_phi[0][2]*PI/180.),v_califa_energy_theta_phi[0][0]*cos(v_califa_energy_theta_phi[0][1]*PI/180.));

					TVector3 sel_vector_2 (v_califa_energy_theta_phi[1][0]*sin(v_califa_energy_theta_phi[1][1]*PI/180.)*cos(v_califa_energy_theta_phi[1][2]*PI/180.),v_califa_energy_theta_phi[1][0]*   sin(v_califa_energy_theta_phi[1][1]*PI/180.)*sin(v_califa_energy_theta_phi[1][2]*PI/180.),v_califa_energy_theta_phi[1][0]*cos(v_califa_energy_theta_phi[1][1]*PI/180.));

					Double_t angle_two_vecs = sel_vector_1.Angle(sel_vector_2);
					h1_opening_angle->Fill(angle_two_vecs*180/PI);
					h1_opening_angle_first_two->Fill(angle_two_vecs*180/PI);
					
			


					//now I select only hits in the PROTON RANGE region
					if( v_califa_energy_theta_phi[0][1] < 43){
						if (v_califa_energy_theta_phi[1][1] < 43){
							h1_opening_angle_proton->Fill(angle_two_vecs*180/PI);
							}
						else {
							if (((v_califa_energy_theta_phi[1][2] > -22.5 && v_califa_energy_theta_phi[1][2] < 22.5) || v_califa_energy_theta_phi[1][2] > 157 || v_califa_energy_theta_phi[1][2] < -   157)){
								h1_opening_angle_proton->Fill(angle_two_vecs*180/PI);
								}
								}
								}
					else{
						if(v_califa_energy_theta_phi[1][1] < 43){
							if(((v_califa_energy_theta_phi[0][2] > -22.5 && v_califa_energy_theta_phi[0][2] < 22.5) || v_califa_energy_theta_phi[0][2] > 157 || v_califa_energy_theta_phi[0][2] < -   157)){
								h1_opening_angle_proton->Fill(angle_two_vecs*180/PI);
							}
							}
						else if(((v_califa_energy_theta_phi[0][2] > -22.5 && v_califa_energy_theta_phi[0][2] < 22.5) || v_califa_energy_theta_phi[0][2] > 157 || v_califa_energy_theta_phi[0][2] < -   157) && ((v_califa_energy_theta_phi[1][2] > -22.5 && v_califa_energy_theta_phi[1][2] < 22.5) || v_califa_energy_theta_phi[1][2] > 157 || v_califa_energy_theta_phi[1][2] < -   157)){
							h1_opening_angle_proton->Fill(angle_two_vecs*180/PI);
							}
					}

					//----END proton range-------------------------------------



					if (angle_two_vecs*180/PI < 30){
						cout << "--Begin, just two entries-----------" << endl;
						cout << "this is the eventnumber with small opening angle:" << evtnr << endl;
						cout << "magnitude of first vector:\t" << sel_vector_1.Mag() << endl;;
						cout << "magnitude of second vector:\t" << sel_vector_2.Mag() << endl;
						cout << "Energy one of v_califa_energy_theta_phi:\t" << v_califa_energy_theta_phi[0][0] << endl;
						cout << "Energy two of v_califa_energy_theta_phi:\t" << v_califa_energy_theta_phi[1][0] << endl;
						cout << "and the according opening angle:\t" << angle_two_vecs*180/PI << endl;
						cout << "---END-------------------------------" << endl;
						}

					h2_opening_angle_vs_esum->Fill(angle_two_vecs*180/PI,v_califa_energy_theta_phi[0][0]+v_califa_energy_theta_phi[1][0]);
					v_califa_energy_theta_phi.clear();
					}
				else{
					v_califa_energy_theta_phi.clear();
					//continue;
					}
				}
			else{
				v_califa_energy_theta_phi.clear();
				//continue;
				}
			if (v_special.size() == 3){
			for (Int_t q = 0; q < 6; q++){
			if (v_special[0] > (50+50*q) && v_special[0] < (900-50*q)){
				h2_charge1_vs_charge2_energy_cut[q]->Fill(charge_1,charge_2);
				if ((charge_1+charge_2) > 92){  //old settings as slide (charge_1+charge_2) > 93.2
					array_for_cuts_energy[q][2] +=1;
					}
				if ((charge_1+charge_2) < 92 && (charge_1+charge_2) > 90.8){ //old settings as slide (charge_1+charge_2) < 93.2 && (charge_1+charge_2) > 92
					array_for_cuts_energy[q][1] +=1;
					}
				if ((charge_1+charge_2) < 90.8){ //old settings as slide (charge_1+charge_2) < 92
					array_for_cuts_energy[q][0] +=1;
					}

				}

			if (v_special[1] > (70+q*1) && v_special[1] < (90-q*1)){
				h2_charge1_vs_charge2_theta_cut[q]->Fill(charge_1,charge_2);

				if ((charge_1+charge_2) > 92){ //old settings as slide (charge_1+charge_2) > 93.
					array_for_cuts_theta[q][2] +=1;
					}
				if ((charge_1+charge_2) < 92 && (charge_1+charge_2) > 90.8){  //old settings as slide(charge_1+charge_2) < 93.2 && (charge_1+charge_2) > 92
					array_for_cuts_theta[q][1] +=1;
					}
				if ((charge_1+charge_2) < 90.8){ //old settings as slide (charge_1+charge_2) < 92
					array_for_cuts_theta[q][0] +=1;
					}
			}
			if (v_special[2] < (30-q*5)){
				h2_charge1_vs_charge2_phi_cut[q]->Fill(charge_1,charge_2);

				if ((charge_1+charge_2) > 92){  // old settings as slide (charge_1+charge_2) > 93.2
					array_for_cuts_phi[q][2] +=1;
					}
				if ((charge_1+charge_2) < 92 && (charge_1+charge_2) > 90.8){ //old settings as slide (charge_1+charge_2) < 93.2 && (charge_1+charge_2) > 92
					array_for_cuts_phi[q][1] +=1;
					}
				if ((charge_1+charge_2) < 90.8){ //old settings as slide (charge_1+charge_2) < 92
					array_for_cuts_phi[q][0] +=1;
					}
					}
			}
			//systematically check theta1+theta2 for different Z_sums -------------
			Double_t charge_sum = charge_1+charge_2;
			if (charge_sum > 91.5 && charge_sum < 92.5){
				h1_theta_sum_z_sum_915_925->Fill(v_special[1]);	
				//this should be the p2p channel, furter steps and cuts...
				if ( v_special[2] < 15){
				h1_theta_sum_z_sum_915_925_phi_cut->Fill(v_special[1]);
				if (v_special[0] > 200 && v_special[0] < 700){
					h1_theta_sum_z_sum_915_925_phi_ener_cut->Fill(v_special[1]);	
					}
				}


				}
			if (charge_sum > 90.2 && charge_sum < 91.4){
				h1_theta_sum_z_sum_902_914->Fill(v_special[1]);
				if (v_special[2] < 15){
					h1_theta_sum_z_sum_902_914_phi_cut->Fill(v_special[1]);
					if (v_special[0] > 200 && v_special[0] < 700){
						h1_theta_sum_z_sum_902_914_phi_ener_cut->Fill(v_special[1]);
						}
					}
				}
			if (charge_sum > 89.1 && charge_sum < 90.1){
				h1_theta_sum_z_sum_891_901->Fill(v_special[1]);
				if (v_special[2] < 15){
					h1_theta_sum_z_sum_891_901_phi_cut->Fill(v_special[1]);
					if (v_special[0] > 200 && v_special[0] < 700){
						h1_theta_sum_z_sum_891_901_phi_ener_cut->Fill(v_special[1]);
						}
					}
				}
			if (charge_sum > 87.9 && charge_sum < 89){
				h1_theta_sum_z_sum_879_89->Fill(v_special[1]);	
				if (v_special[2] < 15){
					h1_theta_sum_z_sum_879_89_phi_cut->Fill(v_special[1]);
					if (v_special[0] > 200 && v_special[0] < 700){
						h1_theta_sum_z_sum_879_89_phi_ener_cut->Fill(v_special[1]);
						}
					}
				}
			if (charge_sum > 86.8 && charge_sum < 87.9){
				h1_theta_sum_z_sum_868_879->Fill(v_special[1]);
				if (v_special[2] < 15){
					h1_theta_sum_z_sum_868_879_phi_cut->Fill(v_special[1]);
					if (v_special[0] > 200 && v_special[0] < 700){
						h1_theta_sum_z_sum_868_879_phi_ener_cut->Fill(v_special[1]);
						}
					}
				}
			//---------------------------------------------------------------------

			//and now systematically check theta1+theta2 for differnt combinations of charge1 and charge2------
			if( (charge_1 > 42.7 && charge_1 < 43.9 && charge_2 > 46.3 && charge_2 < 47.6) || (charge_2 > 42.7 && charge_2 < 43.9 && charge_1 > 46.3 && charge_1 < 47.6)){
				//cut1
				h1_theta_sum_z_1_2_cut1->Fill(v_special[1]);

				}
			if ((charge_1 > 43.9 && charge_1 < 45.3 && charge_2 > 45.3 && charge_2 < 46.3) || (charge_2 > 43.9 && charge_2 < 45.3 && charge_1 > 45.3 && charge_1 < 46.3)){
				//cut2
				h1_theta_sum_z_1_2_cut2->Fill(v_special[1]);

				}
			if ((charge_1 > 42.7 && charge_1 < 43.9 && charge_2 > 45.3 && charge_2 < 46.3) || (charge_2 > 42.7 && charge_2 < 43.9 && charge_1 > 45.3 && charge_1 < 46.3)){
				//cut3
				h1_theta_sum_z_1_2_cut3->Fill(v_special[1]);

				}
			if ((charge_1 > 43.9 && charge_1 < 45.3 && charge_2 > 46.3 && charge_2 < 47.6) || (charge_2 > 43.9 && charge_2 < 45.3 && charge_1 > 46.3 && charge_1 < 47.6)){
				//cut4
				h1_theta_sum_z_1_2_cut4->Fill(v_special[1]);

				}
			if ((charge_1 > 42.7 && charge_1 < 43.9 && charge_2 > 47.8 && charge_2 < 48.7) || (charge_2 > 42.7 && charge_2 < 43.9 && charge_1 > 47.8 && charge_1 < 48.7)){
				//cut5
				h1_theta_sum_z_1_2_cut5->Fill(v_special[1]);

				}
			//if ((charge_1 > 42.7 && charge_1 < 43.9 && charge_2 > 45.3 && charge_2 < 46.3) || (charge_2 > 42.7 && charge_2 < 43.9 && charge_1 > 45.3 && charge_1 < 46.3)){
			//	//cut6
			//	h1_theta_sum_z_1_2_cut6->Fill(v_special[1]);

			//	}

			//-------------------------------------------------------------------------------------------------

			}
			v_special.clear();


			//-------------------------------------------------------------------
			

			if (v_califa_energy_theta.size() >= 2){
				sort(v_califa_energy_theta.begin(),v_califa_energy_theta.end(),sortcol);
				//fill Cluster histogramm
				h2_theta1_vs_theta2_cluster->Fill(v_califa_energy_theta[0][1],v_califa_energy_theta[1][1]);
				h2_e1_vs_e2_cluster->Fill(v_califa_energy_theta[0][0],v_califa_energy_theta[1][0]);


				}
			v_califa_energy_theta.clear();
			v_califa_energy_theta_phi.clear();
			if (califa_200mev_hits == 1){
				h2_charge1_charge_2_one_CALIFA->Fill(charge_1,charge_2);
				h2_charge1_charge_2_diff_vs_z_sum_one_CALIFA->Fill(charge_1-charge_2,charge_1+charge_2);
				h1_charge1_plus_charge2_one_CALIFA->Fill(charge_1+charge_2);
				}
			if (califa_200mev_hits == 2){
				h2_charge1_charge_2_two_CALIFA->Fill(charge_1,charge_2);
				h2_charge1_charge_2_diff_vs_z_sum_two_CALIFA->Fill(charge_1-charge_2,charge_1+charge_2);
				h1_charge1_plus_charge2_two_CALIFA->Fill(charge_1+charge_2);
				}
			if (califa_200mev_hits == 3){
				h2_charge1_charge_2_three_CALIFA->Fill(charge_1,charge_2);
				h2_charge1_charge_2_diff_vs_z_sum_three_CALIFA->Fill(charge_1-charge_2,charge_1+charge_2);
				h1_charge1_plus_charge2_three_CALIFA->Fill(charge_1+charge_2);
				}
			if (califa_100mev_hits ==1){
				h1_mult_100mev->Fill(entries_califa_hit);
				}
			if (califa_100mev_hits == 2){
				h1_mult_100mev_two->Fill(entries_califa_hit);
				Double_t sum_degree = 0;
				Int_t theta_count = 0;
				Double_t theta1_high;
				Double_t theta2_high;
				for (Int_t m = 0; m < entries_califa_hit; m++){
					if(((califahitdata[m]->GetEnergy())/1000) > 100){
						sum_degree += ((califahitdata[m]->GetTheta())/ PI)*180;
						theta_count++;
						if (theta_count == 1){
							theta1_high = ((califahitdata[m]->GetTheta())/ PI)*180;
							}
						if (theta_count == 2){
							theta2_high = ((califahitdata[m]->GetTheta())/ PI)*180;
							}


						}
					}
				h1_theta_sum_p2p->Fill(sum_degree);
				h2_theta1_vs_theta2->Fill(theta1_high,theta2_high);
				}
			if (califa_100mev_hits ==1 && events_one_proton < 20){
				for (Int_t m = 0; m < entries_califa_hit; m++){
					for (Int_t j = 0; j < round((califahitdata[m]->GetEnergy()));j++){
						h2_energy_angular_dist_one[events_one_proton]->Fill(((califahitdata[m]->GetTheta())/PI)*180,               ((califahitdata[m]->GetPhi())/PI)*180);
						}
					}
				events_one_proton++;
				}
			if (califa_100mev_hits == 2 && events_two_proton < 20){
				for (Int_t m = 0; m < entries_califa_hit; m++){
					for (Int_t j = 0; j < round((califahitdata[m]->GetEnergy()));j++){
						h2_energy_angular_dist_two[events_two_proton]->Fill(((califahitdata[m]->GetTheta())/ PI)*180,               ((califahitdata[m]->GetPhi())/PI)*180);
						}

					}
				events_two_proton++;
				}

			}
			}
		else if(entries_wr == 0){
			events_califa_no_master+=1;
			h1_mult_no_M->Fill(entries_califa);
			}
		for (Int_t m = 0; m < entries_califa; m++){
			califamappeddata[m] = (R3BCalifaMappedData*)CalifaMappedData->At(m);
			Int_t crystal_ID = califamappeddata[m]->GetCrystalId();
			if (crystal_ID == 1 || crystal_ID == 2){
				cout << "these are pulser events"<< endl;
				}
					if(map_id_califa_side[crystal_ID] == 1 && crystal_ID != 1 && crystal_ID != 2){ //-->Wixhausen
						califa_wixhausen_ts.push_back(califamappeddata[m]->GetWrts());
						califa_both_sides_ts.push_back(califamappeddata[m]->GetWrts());
						if (entries_wr == 1){
							h1_wr_ts_califa_vs_E_M->Fill(double((califamappeddata[m]->GetEnergy())/1000.));
							e_sum_with_M += double((califamappeddata[m]->GetEnergy())/1000.);
							if(double((califamappeddata[m]->GetEnergy())/1000.) > 20){
								e_high+=1;
								}
							}
						else if (entries_wr == 0){
							h1_wr_ts_califa_vs_E_noM->Fill(double((califamappeddata[m]->GetEnergy())/1000.));
							e_sum_no_M += double((califamappeddata[m]->GetEnergy())/1000.);
							if(double((califamappeddata[m]->GetEnergy())/1000.) > 20){
								e_high+=1;
								}
							}
						}
					else if ( map_id_califa_side[crystal_ID] ==2 && crystal_ID != 1 && crystal_ID != 2){ 	// --->MESSEL
						califa_messel_ts.push_back(califamappeddata[m]->GetWrts());
						califa_both_sides_ts.push_back(califamappeddata[m]->GetWrts());
						if (entries_wr == 1){
							h1_wr_ts_califa_vs_E_M->Fill((califamappeddata[m]->GetEnergy())/1000);
							e_sum_with_M += double((califamappeddata[m]->GetEnergy())/1000.);
							if ( (califamappeddata[m]->GetEnergy())/1000 > 20){
							e_high+=1;
							}
							}
						else if (entries_wr == 0){
							h1_wr_ts_califa_vs_E_noM->Fill((califamappeddata[m]->GetEnergy())/1000);
							e_sum_no_M += double((califamappeddata[m]->GetEnergy())/1000.);
							if(double((califamappeddata[m]->GetEnergy())/1000.) > 20){
								e_high+=1;
								}
							}
						}

				else
					continue;

			}
		if (entries_wr == 1){
			h1_mult_e_high->Fill(e_high);
			}
		if (entries_wr == 0){
			h1_mult_e_high_noM->Fill(e_high);
			}
		h2_energy_vs_mult_no_M->Fill(e_sum_no_M,entries_califa);
		h2_energy_vs_mult_with_M->Fill(e_sum_with_M,entries_califa);
		//choose random hit in event for both messel and wixhausen:
		int rand_messel;
		int rand_wixh;
		int rand_both = rand() % califa_both_sides_ts.size();
		if (califa_messel_ts.size()){
			rand_messel = rand() % califa_messel_ts.size();
			}
		if (califa_wixhausen_ts.size()){
			rand_wixh = rand() % califa_wixhausen_ts.size();
			}
		if (califa_wixhausen_ts.size()){
		for(Int_t j = 1; j < califa_wixhausen_ts.size();++j) {
			int64_t time_diff_wixhausen = califa_wixhausen_ts[j] - califa_wixhausen_ts[0];
			if (abs(time_diff_wixhausen) > 4000){
				//cout << "this is abs. difference:\t" <<  abs(time_diff_wixhausen) <<endl;
				//cout << "and corresponding eventnr:\t" << evtnr << endl;
				}
			}
			}
		if (entries_wr == 1){
			wrmasterdata[0] = (R3BWRMasterData*)WRMasterData->At(0);
			wr_Master_ts = wrmasterdata[0]->GetTimeStamp();
			if (califa_messel_ts.size()){
				int64_t time_diff_messel_WRM_rand = califa_messel_ts[rand_messel] - wr_Master_ts;
				h1_wr_ts_califa_diff_WRM_messel_rand->Fill(time_diff_messel_WRM_rand);
				}
			if (califa_wixhausen_ts.size()){
				int64_t time_diff_wixhausen_WRM_rand = califa_wixhausen_ts[rand_wixh] - wr_Master_ts;
				h1_wr_ts_califa_diff_WRM_wix_rand->Fill(time_diff_wixhausen_WRM_rand);
				}
			if(califa_both_sides_ts.size()){
				int64_t time_diff_both_WRM_rand = califa_both_sides_ts[rand_both] - wr_Master_ts;
				h1_wr_ts_califa_diff_WRM_both_rand->Fill(time_diff_both_WRM_rand);
				}
			}
		if (entries_wr > 1){
			cout << "more than two wr_master timestamps in one event, can that be??" << endl;
			}
		califa_messel_ts.clear();
		califa_wixhausen_ts.clear();
		califa_both_sides_ts.clear();
		events_with_cut++;
		}
		}
	}
	
	}

char f_out_name[500];
sprintf(f_out_name,"../file_output/output_corr_tof_273_10.root");
TFile * f = new TFile(f_out_name,"RECREATE");
TList *l = new TList();
gStyle->SetOptStat(1111111);
l->Add(h1_wr_ts_califa_vs_E_M);
l->Add(h1_wr_ts_califa_vs_E_noM);
l->Add(h1_mult_with_M);
l->Add(h1_mult_no_M);
l->Add(h1_wr_ts_califa_diff_WRM_messel_rand);
l->Add(h1_wr_ts_califa_diff_WRM_wix_rand);
l->Add(h1_wr_ts_califa_diff_WRM_both_rand);
l->Add(h1_mult_e_high);
l->Add(h1_mult_e_high_noM);
l->Add(h2_energy_vs_theta);
l->Add(h2_phi_vs_theta);
l->Add(h2_energy_vs_mult_no_M);
l->Add(h2_energy_vs_mult_with_M);
l->Add(h2_charge1_charge_2);
l->Add(h1_charge1_plus_charge2);
l->Add(h1_daughter_nuclei_charge);
l->Add(h1_trigger);
for(Int_t i = 0; i < 20; i++){
l->Add(h2_energy_angular_dist[i]);
l->Add(h2_energy_angular_dist_one[i]);
l->Add(h2_energy_angular_dist_two[i]);
}
l->Add(h1_mult_100mev);
l->Add(h1_mult_100mev_two);
l->Add(h1_theta_sum_p2p);
l->Add(h2_theta1_vs_theta2);
l->Add(h1_angle_check);
l->Add(h1_angle_check_arzi);
l->Add(h2_theta1_vs_theta2_cluster);
l->Add(h2_e1_vs_e2_cluster);
l->Add(h2_theta1_vs_theta2_cluster_best_dphi);
l->Add(h2_e1_vs_e2_cluster_best_dphi);
l->Add(h1_theta_sum_p2p_best_dphi);
l->Add(h2_e1_vs_e2_cluster_two_highest);
l->Add(h1_opening_angle);
l->Add(h1_opening_angle_proton);
l->Add(h2_opening_angle_vs_esum);
l->Add(h2_e_vs_phi);
l->Add(h1_opening_angle_first_two);
l->Add(h1_opening_angle_one_other);
l->Add(h1_opening_angle_two_other);
l->Add(h2_charge1_charge_2_after);
l->Add(h1_charge_sum_after);
l->Add(h1_theta_sum_z_sum_915_925);
l->Add(h1_theta_sum_z_sum_915_925_phi_cut);
l->Add(h1_theta_sum_z_sum_915_925_phi_ener_cut);
l->Add(h1_theta_sum_z_sum_902_914);
l->Add(h1_theta_sum_z_sum_902_914_phi_cut);
l->Add(h1_theta_sum_z_sum_902_914_phi_ener_cut);
l->Add(h1_theta_sum_z_sum_891_901);
l->Add(h1_theta_sum_z_sum_891_901_phi_cut);
l->Add(h1_theta_sum_z_sum_891_901_phi_ener_cut);
l->Add(h1_theta_sum_z_sum_879_89);
l->Add(h1_theta_sum_z_sum_879_89_phi_cut);
l->Add(h1_theta_sum_z_sum_879_89_phi_ener_cut);
l->Add(h1_theta_sum_z_sum_868_879);
l->Add(h1_theta_sum_z_sum_868_879_phi_cut);
l->Add(h1_theta_sum_z_sum_868_879_phi_ener_cut);
l->Add(h1_theta_sum_z_1_2_cut1);
l->Add(h1_theta_sum_z_1_2_cut2);
l->Add(h1_theta_sum_z_1_2_cut3);
l->Add(h1_theta_sum_z_1_2_cut4);
l->Add(h1_theta_sum_z_1_2_cut5);
//l->Add(h1_theta_sum_z_1_2_cut6);
l->Add(h1_one_hit_twim_charge);
l->Add(h1_one_hit_twim_charge_tpat1);
l->Add(h2_charge1_charge_2_diff_vs_z_sum);
l->Add(h2_charge1_charge_2_tpat1);
l->Add(h2_charge1_charge_2_tpat2);
l->Add(h2_charge1_charge_2_tpat3);
l->Add(h2_charge1_charge_2_tpat4);
l->Add(h2_charge1_charge_2_diff_vs_z_sum_tpat1);
l->Add(h2_charge1_charge_2_diff_vs_z_sum_tpat2);
l->Add(h2_charge1_charge_2_diff_vs_z_sum_tpat3);
l->Add(h2_charge1_charge_2_diff_vs_z_sum_tpat4);
l->Add(h1_charge1_plus_charge2_tpat1);
l->Add(h1_charge1_plus_charge2_tpat2);
l->Add(h1_charge1_plus_charge2_tpat3);
l->Add(h1_charge1_plus_charge2_tpat4);
l->Add(h2_charge1_charge_2_one_CALIFA);
l->Add(h2_charge1_charge_2_two_CALIFA);
l->Add(h2_charge1_charge_2_three_CALIFA);
l->Add(h2_charge1_charge_2_diff_vs_z_sum_one_CALIFA);
l->Add(h2_charge1_charge_2_diff_vs_z_sum_two_CALIFA);
l->Add(h2_charge1_charge_2_diff_vs_z_sum_three_CALIFA);
l->Add(h1_charge1_plus_charge2_one_CALIFA);
l->Add(h1_charge1_plus_charge2_two_CALIFA);
l->Add(h1_charge1_plus_charge2_three_CALIFA);



//char pic_name[500];
//
//for (Int_t i=0; i< 6; i++){
//l->Add(h2_charge1_vs_charge2_energy_cut[i]);
//sprintf(pic_name,"step_energy%i.pdf",i);
//h2_charge1_vs_charge2_energy_cut[i]->SaveAs(pic_name,"pdf");
//}
//for (Int_t i=0; i< 6; i++){
//l->Add(h2_charge1_vs_charge2_theta_cut[i]);
//sprintf(pic_name,"step_theta%i.pdf",i);
//h2_charge1_vs_charge2_theta_cut[i]->SaveAs(pic_name,"pdf");
//}
//for (Int_t i=0; i< 6; i++){
//l->Add(h2_charge1_vs_charge2_phi_cut[i]);
//sprintf(pic_name,"step_phi%i.pdf",i);
//h2_charge1_vs_charge2_phi_cut[i]->SaveAs(pic_name,"pdf");
//}
l->Write("histlist", TObject::kSingleKey);



cout << "----SUMMARY---------------------" << endl;
cout << "Califa events with master trigger:\t" << events_califa_with_master <<  " " << double(events_califa_with_master)/double(events_total_califa) << " % " << endl;
cout << "Califa events no master trigger:\t" << events_califa_no_master << " " << double(events_califa_no_master)/double(events_total_califa) << " % " << endl;
cout << "Total number of Califa events:\t" << events_total_califa << endl;
cout << "--------------------------------" << endl;

//cout << "Array values for Energy parameter change: --------------" << endl;
//cout<< "under the line\t" << "line\t" << "over the line\t" << endl;
//Double_t x[5] = {1,2,3,4,5};
//Double_t y_energy[5];
//Double_t y_theta[5];
//Double_t y_phi[5];
//Double_t y_rel_energy[5];
//Double_t y_rel_theta[5];
//Double_t y_rel_phi[5];
//
//for (Int_t i = 0; i< 6; i++){
//if (i+1 < 6){
//y_energy[i] = array_for_cuts_energy[i+1][0]-array_for_cuts_energy[i][0];
//y_rel_energy[i] = (array_for_cuts_energy[i+1][0]/array_for_cuts_energy[i][0])/(array_for_cuts_energy[i+1][1]/array_for_cuts_energy[i][1] );
//}
//cout << array_for_cuts_energy[i][0] << "\t" << array_for_cuts_energy[i][1] << "\t" << array_for_cuts_energy[i][2] << endl; 
//}
//cout << "Array values for theta parameter change: --------------" << endl;
//cout<< "under the line\t" << "line\t" << "over the line\t" << endl;
//for (Int_t i = 0; i< 6; i++){
//if (i+1 < 6){
//y_theta[i] = array_for_cuts_theta[i+1][0]-array_for_cuts_theta[i][0];
//y_rel_theta[i] = (array_for_cuts_theta[i+1][0]/array_for_cuts_theta[i][0])/(array_for_cuts_theta[i+1][1]/array_for_cuts_theta[i][1]);
//}
//cout << array_for_cuts_theta[i][0] << "\t" << array_for_cuts_theta[i][1] << "\t" << array_for_cuts_theta[i][2] << endl; 
//}
//
//cout << "Array values for phi parameter change: --------------" << endl;
//cout<< "under the line\t" << "line\t" << "over the line\t" << endl;
//for (Int_t i = 0; i< 6; i++){
//if (i+1 < 6){
//y_phi[i] = array_for_cuts_phi[i+1][0]-array_for_cuts_phi[i][0];
//y_rel_phi[i] = (array_for_cuts_phi[i+1][0]/array_for_cuts_phi[i][0])/(array_for_cuts_phi[i+1][1]/array_for_cuts_phi[i][1]);
//}
//cout << array_for_cuts_phi[i][0] << "\t" << array_for_cuts_phi[i][1] << "\t" << array_for_cuts_phi[i][2] << endl; 
//}
//auto c3 = new TCanvas("c3","c3",2400, 1600);
//TGraph* gr_energy = new TGraph(5,x,y_energy);
//gr_energy->SetTitle("Stepwise energy cuts from 50 < E_sum < 900 to 300 < E_sum < 650");
//gr_energy->SetLineColor(kBlue);
//gr_energy->SetLineWidth(3);
//gr_energy->SetMarkerStyle(21);
//
//TGraph* gr_theta = new TGraph(5,x,y_theta);
//gr_theta->SetTitle("Stepwise theta cus from 70 < Theta_sum < 90 to 75 < Theta_sum < 85");
//gr_theta->SetLineColor(kRed);
//gr_theta->SetLineWidth(3);
//gr_theta->SetMarkerStyle(22);
//
//TGraph* gr_phi = new TGraph(5,x,y_phi);
//gr_phi->SetTitle("Stepwise phi cuts from delta_phi = 180+-30 to 180+-5");
//gr_phi->SetLineColor(kGreen);
//gr_phi->SetLineWidth(3);
//gr_phi->SetMarkerStyle(23);
//
//auto mg = new TMultiGraph();
//mg->SetTitle("Parameter cut analysis in region Z_sum < 90.8; Stepnumber; Decrease of events after cut-step");
//mg->Add(gr_energy,"PL");
//mg->Add(gr_theta,"PL");
////gr_phi->Draw("ALP");
//mg->Add(gr_phi,"PL");
//mg->Draw("A pmc plc");
//
//c3->BuildLegend();
//c3->Print("step_approach.png");
//
//TList *l1 = new TList();
//l1->Add(mg);
//l1->Write("cut_steps", TObject::kSingleKey);
//
//
//auto c4 = new TCanvas("c3","c3",2400, 1600);
//TGraph* gr_rel_energy = new TGraph(5,x,y_rel_energy);
//gr_rel_energy->SetTitle("Stepwise energy cuts from 50 < E_sum < 900 to 300 < E_sum < 650 relative decr.");
//gr_rel_energy->SetLineColor(kBlue);
//gr_rel_energy->SetLineWidth(3);
//gr_rel_energy->SetMarkerStyle(21);
//
//TGraph* gr_rel_theta = new TGraph(5,x,y_rel_theta);
//gr_rel_theta->SetTitle("Stepwise theta cus from 70 < Theta_sum < 90 to 75 < Theta_sum < 85 relative decr.");
//gr_rel_theta->SetLineColor(kRed);
//gr_rel_theta->SetLineWidth(3);
//gr_rel_theta->SetMarkerStyle(22);
//
//TGraph* gr_rel_phi = new TGraph(5,x,y_rel_phi);
//gr_rel_phi->SetTitle("Stepwise phi cuts from delta_phi = 180+-30 to 180+-5 relative decr.");
//gr_rel_phi->SetLineColor(kGreen);
//gr_rel_phi->SetLineWidth(3);
//gr_rel_phi->SetMarkerStyle(23);
//
//auto mg2 = new TMultiGraph();
//mg2->SetTitle("Relative decrease of events (lower triangle/correlation line) between step-cuts; Stepnumber; relative decrease lower triangle vs correlation line");
//mg2->Add(gr_rel_energy,"PL");
//mg2->Add(gr_rel_theta,"PL");
//mg2->Add(gr_rel_phi,"PL2");
//mg2->Draw("A pmc plc");
//
//c4->BuildLegend();
//c4->Print("rel_step_approach.png");
//
//TList* l2 = new TList();
//l2->Add(mg2);
//l2->Write("rel_step_numbers",TObject::kSingleKey);

}
