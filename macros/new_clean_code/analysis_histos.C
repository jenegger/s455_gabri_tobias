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

vector<TH2D*> list_of_2D_histos;
vector<TH1D*> list_of_1D_histos;
//Declaration of 1D histos
TH1D* h1_trigger;
TH1D* h1_wr_ts_califa_diff_WRM_messel_rand;
TH1D* h1_wr_ts_califa_diff_WRM_wix_rand;
TH1D* h1_wr_ts_califa_diff_WRM_both_rand;
TH1D* h1_wr_ts_califa_vs_E_M;
TH1D* h1_wr_ts_califa_vs_E_noM;
TH1D* h1_mult_with_M;
TH1D* h1_mult_no_M;
TH1D* h1_mult_e_high;
TH1D* h1_mult_e_high_noM;
TH1D* h1_charge1_plus_charge2;
TH1D* h1_daughter_nuclei_charge;
TH1D* h1_charge_sum_after;
TH1D* h1_charge1_plus_charge2_tpat1;
TH1D* h1_charge1_plus_charge2_tpat2;
TH1D* h1_charge1_plus_charge2_tpat3;
TH1D* h1_charge1_plus_charge2_tpat4;
TH1D* h1_charge1_plus_charge2_one_CALIFA;
TH1D* h1_charge1_plus_charge2_two_CALIFA;
TH1D* h1_charge1_plus_charge2_three_CALIFA;
TH1D* h1_mult_100mev;
TH1D* h1_mult_100mev_two;
TH1D* h1_theta_sum_p2p;
TH1D* h1_theta_sum_z_sum_915_925;
TH1D* h1_theta_sum_z_sum_915_925_phi_cut;
TH1D* h1_theta_sum_z_sum_915_925_phi_ener_cut;
TH1D* h1_theta_sum_z_sum_902_914;
TH1D* h1_theta_sum_z_sum_902_914_phi_cut;
TH1D* h1_theta_sum_z_sum_891_901;
TH1D* h1_theta_sum_z_sum_902_914_phi_ener_cut;
TH1D* h1_theta_sum_z_sum_879_89;
TH1D* h1_theta_sum_z_sum_891_901_phi_cut;
TH1D* h1_theta_sum_z_sum_879_89_phi_ener_cut;
TH1D* h1_theta_sum_z_sum_891_901_phi_ener_cut;
TH1D* h1_theta_sum_z_sum_868_879;
TH1D* h1_theta_sum_z_sum_879_89_phi_cut;
TH1D* h1_theta_sum_z_sum_868_879_phi_cut;
TH1D* h1_theta_sum_z_sum_868_879_phi_ener_cut;
TH1D* h1_theta_sum_z_1_2_cut1;
TH1D* h1_theta_sum_z_1_2_cut2;
TH1D* h1_theta_sum_z_1_2_cut3;
TH1D* h1_theta_sum_z_1_2_cut4;
TH1D* h1_theta_sum_z_1_2_cut5;
TH1D* h1_one_hit_twim_charge;
TH1D* h1_one_hit_twim_charge_tpat1;
TH1D* h1_theta_sum_p2p_best_dphi;
TH1D* h1_opening_angle;
TH1D* h1_opening_angle_first_two;
TH1D* h1_opening_angle_one_other;
TH1D* h1_opening_angle_two_other;
TH1D* h1_opening_angle_proton;
TH1D* h1_angle_check;
TH1D* h1_angle_check_arzi;
//End of Declaration of 1D histos

//Declaration of 2D histos
TH2D* h2_energy_vs_theta;
TH2D* h2_phi_vs_theta;
TH2D* h2_energy_vs_mult_no_M;
TH2D* h2_energy_vs_mult_with_M;
TH2D* h2_charge1_charge_2;
TH2D* h2_charge1_charge_2_diff_vs_z_sum;
TH2D* h2_charge1_charge_2_after;
TH2D* h2_charge1_charge_2_tpat1;
TH2D* h2_charge1_charge_2_tpat2;
TH2D* h2_charge1_charge_2_tpat3;
TH2D* h2_charge1_charge_2_tpat4;
TH2D* h2_charge1_charge_2_one_CALIFA;
TH2D* h2_charge1_charge_2_two_CALIFA;
TH2D* h2_charge1_charge_2_three_CALIFA;
TH2D* h2_charge1_charge_2_diff_vs_z_sum_tpat1;
TH2D* h2_charge1_charge_2_diff_vs_z_sum_two_CALIFA;
TH2D* h2_charge1_charge_2_diff_vs_z_sum_one_CALIFA;
TH2D* h2_charge1_charge_2_diff_vs_z_sum_tpat4;
TH2D* h2_charge1_charge_2_diff_vs_z_sum_tpat3;
TH2D* h2_charge1_charge_2_diff_vs_z_sum_tpat2;
TH2D* h2_charge1_charge_2_diff_vs_z_sum_three_CALIFA;
TH2D* h2_energy_angular_dist[20];
TH2D* h2_energy_angular_dist_one[20];
TH2D* h2_energy_angular_dist_two[20];
TH2D* h2_charge1_vs_charge2_energy_cut[6];
TH2D* h2_charge1_vs_charge2_theta_cut[6];
TH2D* h2_charge1_vs_charge2_phi_cut[6];
TH2D* h2_theta1_vs_theta2;
TH2D* h2_theta1_vs_theta2_cluster;
TH2D* h2_theta1_vs_theta2_cluster_best_dphi;
TH2D* h2_e1_vs_e2_cluster;
TH2D* h2_e1_vs_e2_cluster_best_dphi;
TH2D* h2_e1_vs_e2_cluster_two_highest;
TH2D* h2_opening_angle_vs_esum;
TH2D* h2_e_vs_phi;
//End of Declaration of 2D histos



char hist_name[500];
void my_histos(){
//1D_histos-----------------------------------
sprintf(hist_name, "WR Messel -  WR Master (one random hit for event)");
h1_wr_ts_califa_diff_WRM_messel_rand = new TH1D(hist_name,hist_name,500,-5000,5000);
h1_wr_ts_califa_diff_WRM_messel_rand->GetXaxis()->SetTitle(" wr time difference [ns]");
h1_wr_ts_califa_diff_WRM_messel_rand->GetYaxis()->SetTitle("Counts");
h1_wr_ts_califa_diff_WRM_messel_rand->GetXaxis()->CenterTitle(true);
h1_wr_ts_califa_diff_WRM_messel_rand->GetYaxis()->CenterTitle(true);
h1_wr_ts_califa_diff_WRM_messel_rand->GetYaxis()->SetLabelSize(0.045);
h1_wr_ts_califa_diff_WRM_messel_rand->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "WR Wixhausen - WR Master (one random hit for event)");
h1_wr_ts_califa_diff_WRM_wix_rand = new TH1D(hist_name,hist_name,500,-5000,5000);
h1_wr_ts_califa_diff_WRM_wix_rand->GetXaxis()->SetTitle("wr time difference [ns]");
h1_wr_ts_califa_diff_WRM_wix_rand->GetYaxis()->SetTitle("Counts");
h1_wr_ts_califa_diff_WRM_wix_rand->GetXaxis()->CenterTitle(true);
h1_wr_ts_califa_diff_WRM_wix_rand->GetYaxis()->CenterTitle(true);
h1_wr_ts_califa_diff_WRM_wix_rand->GetYaxis()->SetLabelSize(0.045);
h1_wr_ts_califa_diff_WRM_wix_rand->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "WR CALIFA - WR Master (one random hit for event wixhausen or messel)");
h1_wr_ts_califa_diff_WRM_both_rand = new TH1D(hist_name,hist_name,500,-5000,5000);
h1_wr_ts_califa_diff_WRM_both_rand->GetXaxis()->SetTitle("wr time difference [ns]");
h1_wr_ts_califa_diff_WRM_both_rand->GetYaxis()->SetTitle("Counts");
h1_wr_ts_califa_diff_WRM_both_rand->GetXaxis()->CenterTitle(true);
h1_wr_ts_califa_diff_WRM_both_rand->GetYaxis()->CenterTitle(true);
h1_wr_ts_califa_diff_WRM_both_rand->GetYaxis()->SetLabelSize(0.045);
h1_wr_ts_califa_diff_WRM_both_rand->GetYaxis()->SetTitleSize(0.045);

//Energy vs ts difference for events with wr master trigger
sprintf(hist_name, "Hit Energy (events with master trigger)");
h1_wr_ts_califa_vs_E_M = new TH1D(hist_name,hist_name,1000,0,10);
h1_wr_ts_califa_vs_E_M->GetXaxis()->SetTitle("Energy [MeV]");
h1_wr_ts_califa_vs_E_M->GetYaxis()->SetTitle("Counts");
h1_wr_ts_califa_vs_E_M->GetXaxis()->CenterTitle(true);
h1_wr_ts_califa_vs_E_M->GetYaxis()->CenterTitle(true);
h1_wr_ts_califa_vs_E_M->GetYaxis()->SetLabelSize(0.045);
h1_wr_ts_califa_vs_E_M->GetYaxis()->SetTitleSize(0.045);

//Energy vs ts difference for events without wr master trigger
sprintf(hist_name, "Hit Energy (events NO master trigger)");
h1_wr_ts_califa_vs_E_noM = new TH1D(hist_name,hist_name,1000,0,10);
h1_wr_ts_califa_vs_E_noM->GetXaxis()->SetTitle("Energy [MeV]");
h1_wr_ts_califa_vs_E_noM->GetYaxis()->SetTitle("Counts");
h1_wr_ts_califa_vs_E_noM->GetXaxis()->CenterTitle(true);
h1_wr_ts_califa_vs_E_noM->GetYaxis()->CenterTitle(true);
h1_wr_ts_califa_vs_E_noM->GetYaxis()->SetLabelSize(0.045);
h1_wr_ts_califa_vs_E_noM->GetYaxis()->SetTitleSize(0.045);


//Multiplicity for events with wr master

sprintf(hist_name, "Multiplicity for events with master trigger");
h1_mult_with_M = new TH1D(hist_name,hist_name,2500,0,2500);
h1_mult_with_M->GetXaxis()->SetTitle("Multiplicity");
h1_mult_with_M->GetYaxis()->SetTitle("Counts");
h1_mult_with_M->GetXaxis()->CenterTitle(true);
h1_mult_with_M->GetYaxis()->CenterTitle(true);
h1_mult_with_M->GetYaxis()->SetLabelSize(0.045);
h1_mult_with_M->GetYaxis()->SetTitleSize(0.045);


//Multiplicity for events without wr master

sprintf(hist_name, "Multiplicity for events NO master trigger");
h1_mult_no_M = new TH1D(hist_name,hist_name,2500,0,2500);
h1_mult_no_M->GetXaxis()->SetTitle("Multiplicity");
h1_mult_no_M->GetYaxis()->SetTitle("Counts");
h1_mult_no_M->GetXaxis()->CenterTitle(true);
h1_mult_no_M->GetYaxis()->CenterTitle(true);
h1_mult_no_M->GetYaxis()->SetLabelSize(0.045);
h1_mult_no_M->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Multiplicity of hits with E_lab > 20 MeV (with master trigger)");
h1_mult_e_high = new TH1D(hist_name,hist_name,2500,0,2500);
h1_mult_e_high->GetXaxis()->SetTitle("Multiplicity");
h1_mult_e_high->GetYaxis()->SetTitle("Counts");
h1_mult_e_high->GetXaxis()->CenterTitle(true);
h1_mult_e_high->GetYaxis()->CenterTitle(true);
h1_mult_e_high->GetYaxis()->SetLabelSize(0.045);
h1_mult_e_high->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Multiplicity of hits with E_lab > 20 MeV (events without master trigger)");
h1_mult_e_high_noM = new TH1D(hist_name,hist_name,2500,0,2500);
h1_mult_e_high_noM->GetXaxis()->SetTitle("Multiplicity");
h1_mult_e_high_noM->GetYaxis()->SetTitle("Counts");
h1_mult_e_high_noM->GetXaxis()->CenterTitle(true);
h1_mult_e_high_noM->GetYaxis()->CenterTitle(true);
h1_mult_e_high_noM->GetYaxis()->SetLabelSize(0.045);
h1_mult_e_high_noM->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 + Charge 2  (TWIM Music)");
h1_charge1_plus_charge2 = new TH1D(hist_name,hist_name,500,0,100);
h1_charge1_plus_charge2->GetXaxis()->SetTitle("Charge Sum");
h1_charge1_plus_charge2->GetYaxis()->SetTitle("Counts");
h1_charge1_plus_charge2->GetXaxis()->CenterTitle(true);
h1_charge1_plus_charge2->GetYaxis()->CenterTitle(true);
h1_charge1_plus_charge2->GetYaxis()->SetLabelSize(0.045);
h1_charge1_plus_charge2->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge of Daughter Nuclei (TWIM Music)");
h1_daughter_nuclei_charge = new TH1D(hist_name,hist_name,500,0,100);
h1_daughter_nuclei_charge->GetXaxis()->SetTitle("Charge Z");
h1_daughter_nuclei_charge->GetYaxis()->SetTitle("Counts");
h1_daughter_nuclei_charge->GetXaxis()->CenterTitle(true);
h1_daughter_nuclei_charge->GetYaxis()->CenterTitle(true);
h1_daughter_nuclei_charge->GetYaxis()->SetLabelSize(0.045);
h1_daughter_nuclei_charge->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge Sum (TWIM Music) after the minimal CALIFA cuts");
h1_charge_sum_after = new TH1D(hist_name,hist_name,500,0,100);
h1_charge_sum_after->GetXaxis()->SetTitle("Charge 1 + Charge 2");
h1_charge_sum_after->GetYaxis()->SetTitle("Counts");
h1_charge_sum_after->GetXaxis()->CenterTitle(true);
h1_charge_sum_after->GetYaxis()->CenterTitle(true);
h1_charge_sum_after->GetYaxis()->SetLabelSize(0.045);
h1_charge_sum_after->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 + Charge 2  (TWIM Music) for TPat = SofStart");
h1_charge1_plus_charge2_tpat1 = new TH1D(hist_name,hist_name,500,0,100);
h1_charge1_plus_charge2_tpat1->GetXaxis()->SetTitle("Charge Sum");
h1_charge1_plus_charge2_tpat1->GetYaxis()->SetTitle("Counts");
h1_charge1_plus_charge2_tpat1->GetXaxis()->CenterTitle(true);
h1_charge1_plus_charge2_tpat1->GetYaxis()->CenterTitle(true);
h1_charge1_plus_charge2_tpat1->GetYaxis()->SetLabelSize(0.045);
h1_charge1_plus_charge2_tpat1->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Charge 1 + Charge 2  (TWIM Music) for TPat = Fission");
h1_charge1_plus_charge2_tpat2 = new TH1D(hist_name,hist_name,500,0,100);
h1_charge1_plus_charge2_tpat2->GetXaxis()->SetTitle("Charge Sum");
h1_charge1_plus_charge2_tpat2->GetYaxis()->SetTitle("Counts");
h1_charge1_plus_charge2_tpat2->GetXaxis()->CenterTitle(true);
h1_charge1_plus_charge2_tpat2->GetYaxis()->CenterTitle(true);
h1_charge1_plus_charge2_tpat2->GetYaxis()->SetLabelSize(0.045);
h1_charge1_plus_charge2_tpat2->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 + Charge 2  (TWIM Music) for TPat =  CalifOR (= S+C&&)");
h1_charge1_plus_charge2_tpat3 = new TH1D(hist_name,hist_name,500,0,100);
h1_charge1_plus_charge2_tpat3->GetXaxis()->SetTitle("Charge Sum");
h1_charge1_plus_charge2_tpat3->GetYaxis()->SetTitle("Counts");
h1_charge1_plus_charge2_tpat3->GetXaxis()->CenterTitle(true);
h1_charge1_plus_charge2_tpat3->GetYaxis()->CenterTitle(true);
h1_charge1_plus_charge2_tpat3->GetYaxis()->SetLabelSize(0.045);
h1_charge1_plus_charge2_tpat3->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 + Charge 2  (TWIM Music) for TPat =  CalifAND (= S+F+C&&)");
h1_charge1_plus_charge2_tpat4 = new TH1D(hist_name,hist_name,500,0,100);
h1_charge1_plus_charge2_tpat4->GetXaxis()->SetTitle("Charge Sum");
h1_charge1_plus_charge2_tpat4->GetYaxis()->SetTitle("Counts");
h1_charge1_plus_charge2_tpat4->GetXaxis()->CenterTitle(true);
h1_charge1_plus_charge2_tpat4->GetYaxis()->CenterTitle(true);
h1_charge1_plus_charge2_tpat4->GetYaxis()->SetLabelSize(0.045);
h1_charge1_plus_charge2_tpat4->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Charge 1 + Charge 2  (TWIM Music) for one hit > 200 MeV in CALIFA");
h1_charge1_plus_charge2_one_CALIFA = new TH1D(hist_name,hist_name,500,0,100);
h1_charge1_plus_charge2_one_CALIFA->GetXaxis()->SetTitle("Charge Sum");
h1_charge1_plus_charge2_one_CALIFA->GetYaxis()->SetTitle("Counts");
h1_charge1_plus_charge2_one_CALIFA->GetXaxis()->CenterTitle(true);
h1_charge1_plus_charge2_one_CALIFA->GetYaxis()->CenterTitle(true);
h1_charge1_plus_charge2_one_CALIFA->GetYaxis()->SetLabelSize(0.045);
h1_charge1_plus_charge2_one_CALIFA->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 + Charge 2  (TWIM Music) for two hits > 200 MeV in CALIFA");
h1_charge1_plus_charge2_two_CALIFA = new TH1D(hist_name,hist_name,500,0,100);
h1_charge1_plus_charge2_two_CALIFA->GetXaxis()->SetTitle("Charge Sum");
h1_charge1_plus_charge2_two_CALIFA->GetYaxis()->SetTitle("Counts");
h1_charge1_plus_charge2_two_CALIFA->GetXaxis()->CenterTitle(true);
h1_charge1_plus_charge2_two_CALIFA->GetYaxis()->CenterTitle(true);
h1_charge1_plus_charge2_two_CALIFA->GetYaxis()->SetLabelSize(0.045);
h1_charge1_plus_charge2_two_CALIFA->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 + Charge 2  (TWIM Music) for three hits > 200 MeV in CALIFA");
h1_charge1_plus_charge2_three_CALIFA = new TH1D(hist_name,hist_name,500,0,100);
h1_charge1_plus_charge2_three_CALIFA->GetXaxis()->SetTitle("Charge Sum");
h1_charge1_plus_charge2_three_CALIFA->GetYaxis()->SetTitle("Counts");
h1_charge1_plus_charge2_three_CALIFA->GetXaxis()->CenterTitle(true);
h1_charge1_plus_charge2_three_CALIFA->GetYaxis()->CenterTitle(true);
h1_charge1_plus_charge2_three_CALIFA->GetYaxis()->SetLabelSize(0.045);
h1_charge1_plus_charge2_three_CALIFA->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Multiplicity of events with one hit > 100 MeV ");
h1_mult_100mev = new TH1D(hist_name,hist_name,2500,0,2500);
h1_mult_100mev->GetXaxis()->SetTitle("Multiplicity");
h1_mult_100mev->GetYaxis()->SetTitle("Counts");
h1_mult_100mev->GetXaxis()->CenterTitle(true);
h1_mult_100mev->GetYaxis()->CenterTitle(true);
h1_mult_100mev->GetYaxis()->SetLabelSize(0.045);
h1_mult_100mev->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Multiplicity of events with two hits > 100 MeV ");
h1_mult_100mev_two = new TH1D(hist_name,hist_name,2500,0,2500);
h1_mult_100mev_two->GetXaxis()->SetTitle("Multiplicity");
h1_mult_100mev_two->GetYaxis()->SetTitle("Counts");
h1_mult_100mev_two->GetXaxis()->CenterTitle(true);
h1_mult_100mev_two->GetYaxis()->CenterTitle(true);
h1_mult_100mev_two->GetYaxis()->SetLabelSize(0.045);
h1_mult_100mev_two->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Theta1 +  Theta 2 for the two hits > 100 MeV ");
h1_theta_sum_p2p = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_p2p->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_p2p->GetYaxis()->SetTitle("Counts");
h1_theta_sum_p2p->GetXaxis()->CenterTitle(true);
h1_theta_sum_p2p->GetYaxis()->CenterTitle(true);
h1_theta_sum_p2p->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_p2p->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta1 +  Theta 2 for 91.5 < Z_sum < 92.5");
h1_theta_sum_z_sum_915_925 = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_915_925->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_915_925->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_915_925->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_915_925->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_915_925->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_915_925->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Theta1 +  Theta 2 for 91.5 < Z_sum < 92.5 and delta_phi = 180+-15");
h1_theta_sum_z_sum_915_925_phi_cut = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_915_925_phi_cut->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_915_925_phi_cut->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_915_925_phi_cut->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_915_925_phi_cut->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_915_925_phi_cut->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_915_925_phi_cut->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Theta1 +  Theta 2 for 91.5 < Z_sum < 92.5 and delta_phi = 180+-15 and  200 < E_sum < 700");
h1_theta_sum_z_sum_915_925_phi_ener_cut = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_915_925_phi_ener_cut->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_915_925_phi_ener_cut->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_915_925_phi_ener_cut->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_915_925_phi_ener_cut->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_915_925_phi_ener_cut->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_915_925_phi_ener_cut->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Theta1 +  Theta 2 for 90.2 < Z_sum < 91.4");
h1_theta_sum_z_sum_902_914 = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_902_914->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_902_914->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_902_914->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_902_914->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_902_914->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_902_914->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta1 +  Theta 2 for 90.2 < Z_sum < 91.4 and delta_phi = 180+-15");
h1_theta_sum_z_sum_902_914_phi_cut = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_902_914_phi_cut->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_902_914_phi_cut->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_902_914_phi_cut->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_902_914_phi_cut->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_902_914_phi_cut->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_902_914_phi_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta1 +  Theta 2 for 90.2 < Z_sum < 91.4 and delta_phi = 180+-15 and 200 < E_sum < 700");
h1_theta_sum_z_sum_902_914_phi_ener_cut = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_902_914_phi_ener_cut->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_902_914_phi_ener_cut->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_902_914_phi_ener_cut->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_902_914_phi_ener_cut->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_902_914_phi_ener_cut->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_902_914_phi_ener_cut->GetYaxis()->SetTitleSize(0.045);



sprintf(hist_name, "Theta1 +  Theta 2 for 89.1 < Z_sum < 90.1");
h1_theta_sum_z_sum_891_901 = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_891_901->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_891_901->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_891_901->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_891_901->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_891_901->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_891_901->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta1 +  Theta 2 for 89.1 < Z_sum < 90.1 and delta_phi = 180+-15");
h1_theta_sum_z_sum_891_901_phi_cut = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_891_901_phi_cut->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_891_901_phi_cut->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_891_901_phi_cut->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_891_901_phi_cut->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_891_901_phi_cut->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_891_901_phi_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta1 +  Theta 2 for 89.1 < Z_sum < 90.1 and delta_phi = 180+-15 and 200 < E_sum < 700");
h1_theta_sum_z_sum_891_901_phi_ener_cut = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_891_901_phi_ener_cut->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_891_901_phi_ener_cut->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_891_901_phi_ener_cut->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_891_901_phi_ener_cut->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_891_901_phi_ener_cut->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_891_901_phi_ener_cut->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Theta1 +  Theta 2 for 87.9 < Z_sum < 89");
h1_theta_sum_z_sum_879_89 = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_879_89->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_879_89->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_879_89->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_879_89->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_879_89->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_879_89->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta1 +  Theta 2 for 87.9 < Z_sum < 89 and delta_phi = 180+-15");
h1_theta_sum_z_sum_879_89_phi_cut = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_879_89_phi_cut->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_879_89_phi_cut->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_879_89_phi_cut->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_879_89_phi_cut->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_879_89_phi_cut->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_879_89_phi_cut->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta1 +  Theta 2 for 87.9 < Z_sum < 89 and delta_phi = 180+-15 and 200 < E_sum < 700");
h1_theta_sum_z_sum_879_89_phi_ener_cut = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_879_89_phi_ener_cut->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_879_89_phi_ener_cut->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_879_89_phi_ener_cut->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_879_89_phi_ener_cut->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_879_89_phi_ener_cut->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_879_89_phi_ener_cut->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Theta1 +  Theta 2 for 86.8 < Z_sum < 87.9");
h1_theta_sum_z_sum_868_879 = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_868_879->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_868_879->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_868_879->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_868_879->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_868_879->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_868_879->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta1 +  Theta 2 for 86.8 < Z_sum < 87.9 and delta_phi = 180+-15");
h1_theta_sum_z_sum_868_879_phi_cut = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_sum_868_879_phi_cut->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_sum_868_879_phi_cut->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_sum_868_879_phi_cut->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_sum_868_879_phi_cut->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_sum_868_879_phi_cut->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_sum_868_879_phi_cut->GetYaxis()->SetTitleSize(0.045);

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

sprintf(hist_name, "Theta1 +  Theta 2 for 42.7 < Z_1 < 43.9 and 46.3 < Z_2 < 47.6");
h1_theta_sum_z_1_2_cut1 = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_1_2_cut1->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_1_2_cut1->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_1_2_cut1->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_1_2_cut1->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_1_2_cut1->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_1_2_cut1->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Theta1 +  Theta 2 for 43.9 < Z_1 < 45.3 and 45.3 < Z_2 < 46.3");
h1_theta_sum_z_1_2_cut2 = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_1_2_cut2->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_1_2_cut2->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_1_2_cut2->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_1_2_cut2->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_1_2_cut2->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_1_2_cut2->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta1 +  Theta 2 for 42.7 < Z_1 < 43.9 and 45.3 < Z_2 < 46.3");
h1_theta_sum_z_1_2_cut3 = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_1_2_cut3->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_1_2_cut3->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_1_2_cut3->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_1_2_cut3->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_1_2_cut3->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_1_2_cut3->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta1 +  Theta 2 for 43.9 < Z_1 < 45.3 and 46.3 < Z_2 < 47.6");
h1_theta_sum_z_1_2_cut4 = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_1_2_cut4->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_1_2_cut4->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_1_2_cut4->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_1_2_cut4->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_1_2_cut4->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_1_2_cut4->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta1 +  Theta 2 for 42.7 < Z_1 < 43.9 and 47.8 < Z_2 < 48.7");
h1_theta_sum_z_1_2_cut5 = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_z_1_2_cut5->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_z_1_2_cut5->GetYaxis()->SetTitle("Counts");
h1_theta_sum_z_1_2_cut5->GetXaxis()->CenterTitle(true);
h1_theta_sum_z_1_2_cut5->GetYaxis()->CenterTitle(true);
h1_theta_sum_z_1_2_cut5->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_z_1_2_cut5->GetYaxis()->SetTitleSize(0.045);

////Looking at events with only one hit in TWIM:
sprintf(hist_name, "Charge in TWIM when just having one hit");
h1_one_hit_twim_charge = new TH1D(hist_name,hist_name,500,0,100);
h1_one_hit_twim_charge->GetXaxis()->SetTitle("Charge of single Particle hit");
h1_one_hit_twim_charge->GetYaxis()->SetTitle("Counts");
h1_one_hit_twim_charge->GetXaxis()->CenterTitle(true);
h1_one_hit_twim_charge->GetYaxis()->CenterTitle(true);
h1_one_hit_twim_charge->GetYaxis()->SetLabelSize(0.045);
h1_one_hit_twim_charge->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Charge in TWIM when just having one hit and tpat == 1");
h1_one_hit_twim_charge_tpat1 = new TH1D(hist_name,hist_name,500,0,100);
h1_one_hit_twim_charge_tpat1->GetXaxis()->SetTitle("Charge of single Particle hit");
h1_one_hit_twim_charge_tpat1->GetYaxis()->SetTitle("Counts");
h1_one_hit_twim_charge_tpat1->GetXaxis()->CenterTitle(true);
h1_one_hit_twim_charge_tpat1->GetYaxis()->CenterTitle(true);
h1_one_hit_twim_charge_tpat1->GetYaxis()->SetLabelSize(0.045);
h1_one_hit_twim_charge_tpat1->GetYaxis()->SetTitleSize(0.045);


//---------------------------------------

sprintf(hist_name, "Theta1 +  Theta 2 for two hits (> 30 MeV) and best dphi");
h1_theta_sum_p2p_best_dphi = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_theta_sum_p2p_best_dphi->GetXaxis()->SetTitle("Theta1 + Theta2 [degr]");
h1_theta_sum_p2p_best_dphi->GetYaxis()->SetTitle("Counts");
h1_theta_sum_p2p_best_dphi->GetXaxis()->CenterTitle(true);
h1_theta_sum_p2p_best_dphi->GetYaxis()->CenterTitle(true);
h1_theta_sum_p2p_best_dphi->GetYaxis()->SetLabelSize(0.045);
h1_theta_sum_p2p_best_dphi->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Opening Angle for two hits (> 30 MeV) and best dphi");
h1_opening_angle = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_opening_angle->GetXaxis()->SetTitle("Opening Angle [degr]");
h1_opening_angle->GetYaxis()->SetTitle("Counts");
h1_opening_angle->GetXaxis()->CenterTitle(true);
h1_opening_angle->GetYaxis()->CenterTitle(true);
h1_opening_angle->GetYaxis()->SetLabelSize(0.045);
h1_opening_angle->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Opening Angle for the first two high energetic hits (> 30 MeV) and best dphi");
h1_opening_angle_first_two = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_opening_angle_first_two->GetXaxis()->SetTitle("Opening Angle [degr]");
h1_opening_angle_first_two->GetYaxis()->SetTitle("Counts");
h1_opening_angle_first_two->GetXaxis()->CenterTitle(true);
h1_opening_angle_first_two->GetYaxis()->CenterTitle(true);
h1_opening_angle_first_two->GetYaxis()->SetLabelSize(0.045);
h1_opening_angle_first_two->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Opening Angle for the first and other (not second) energetic hits (> 30 MeV) and best dphi");
h1_opening_angle_one_other = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_opening_angle_one_other->GetXaxis()->SetTitle("Opening Angle [degr]");
h1_opening_angle_one_other->GetYaxis()->SetTitle("Counts");
h1_opening_angle_one_other->GetXaxis()->CenterTitle(true);
h1_opening_angle_one_other->GetYaxis()->CenterTitle(true);
h1_opening_angle_one_other->GetYaxis()->SetLabelSize(0.045);
h1_opening_angle_one_other->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Opening Angle for the second and other (not first) energetic hits (> 30 MeV) and best dphi");
h1_opening_angle_two_other = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_opening_angle_two_other->GetXaxis()->SetTitle("Opening Angle [degr]");
h1_opening_angle_two_other->GetYaxis()->SetTitle("Counts");
h1_opening_angle_two_other->GetXaxis()->CenterTitle(true);
h1_opening_angle_two_other->GetYaxis()->CenterTitle(true);
h1_opening_angle_two_other->GetYaxis()->SetLabelSize(0.045);
h1_opening_angle_two_other->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Opening Angle for two hits (> 30 MeV) and best dphi only in proton range");
h1_opening_angle_proton = new TH1D(hist_name,hist_name,52,22.15,152.15);
h1_opening_angle_proton->GetXaxis()->SetTitle("Opening Angle [degr]");
h1_opening_angle_proton->GetYaxis()->SetTitle("Counts");
h1_opening_angle_proton->GetXaxis()->CenterTitle(true);
h1_opening_angle_proton->GetYaxis()->CenterTitle(true);
h1_opening_angle_proton->GetYaxis()->SetLabelSize(0.045);
h1_opening_angle_proton->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "angle check ");
h1_angle_check = new TH1D(hist_name,hist_name,3500,0,3.5);
h1_angle_check->GetXaxis()->SetTitle("Polar angle");
h1_angle_check->GetYaxis()->SetTitle("Counts");
h1_angle_check->GetXaxis()->CenterTitle(true);
h1_angle_check->GetYaxis()->CenterTitle(true);
h1_angle_check->GetYaxis()->SetLabelSize(0.045);
h1_angle_check->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "angle check arzimuthal ");
h1_angle_check_arzi = new TH1D(hist_name,hist_name,7000,-3.5,3.5);
h1_angle_check_arzi->GetXaxis()->SetTitle("Arzimuthal angle");
h1_angle_check_arzi->GetYaxis()->SetTitle("Counts");
h1_angle_check_arzi->GetXaxis()->CenterTitle(true);
h1_angle_check_arzi->GetYaxis()->CenterTitle(true);
h1_angle_check_arzi->GetYaxis()->SetLabelSize(0.045);
h1_angle_check_arzi->GetYaxis()->SetTitleSize(0.045);

//Histo with trigger info:
sprintf(hist_name,"Trigger information: Tpat");
h1_trigger = new TH1D(hist_name, hist_name, 17, -0.5, 16.5);
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

//end of 1D_histos-----------------


//2D histos------------------------
sprintf(hist_name, "E_lab vs. polar angle (with master trigger)");
h2_energy_vs_theta = new TH2D(hist_name,hist_name,1000,0,1000,25,22.15,84.65);
h2_energy_vs_theta->GetXaxis()->SetTitle("Energy in laboratory system [MeV]");
h2_energy_vs_theta->GetYaxis()->SetTitle("Polar angle Theta [degrees]");
h2_energy_vs_theta->GetXaxis()->CenterTitle(true);
h2_energy_vs_theta->GetYaxis()->CenterTitle(true);
h2_energy_vs_theta->GetYaxis()->SetLabelSize(0.045);
h2_energy_vs_theta->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Phi vs Theta (with master trigger)");
h2_phi_vs_theta = new TH2D(hist_name,hist_name,25,22.15,84.65,60,-180,180);
h2_phi_vs_theta->GetXaxis()->SetTitle("Polar angle  Theta[degrees]");
h2_phi_vs_theta->GetYaxis()->SetTitle("Arzimuthal angle Phi [degrees]");
h2_phi_vs_theta->GetXaxis()->CenterTitle(true);
h2_phi_vs_theta->GetYaxis()->CenterTitle(true);
h2_phi_vs_theta->GetYaxis()->SetLabelSize(0.045);
h2_phi_vs_theta->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Summed up E_lab(mapped level) vs Multiplicity, Uncorrelated");
h2_energy_vs_mult_no_M = new TH2D(hist_name,hist_name,1000,0,1000,400,0,400);
h2_energy_vs_mult_no_M->GetXaxis()->SetTitle("Summed Energy in laboratory system [MeV]");
h2_energy_vs_mult_no_M->GetYaxis()->SetTitle("Multiplicity");
h2_energy_vs_mult_no_M->GetXaxis()->CenterTitle(true);
h2_energy_vs_mult_no_M->GetYaxis()->CenterTitle(true);
h2_energy_vs_mult_no_M->GetYaxis()->SetLabelSize(0.045);
h2_energy_vs_mult_no_M->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Summed up E_lab(mapped level) vs Multiplicity, Correlated");
h2_energy_vs_mult_with_M = new TH2D(hist_name,hist_name,1000,0,1000,400,0,400);
h2_energy_vs_mult_with_M->GetXaxis()->SetTitle("Summed Energy in laboratory system [MeV]");
h2_energy_vs_mult_with_M->GetYaxis()->SetTitle("Multiplicity");
h2_energy_vs_mult_with_M->GetXaxis()->CenterTitle(true);
h2_energy_vs_mult_with_M->GetYaxis()->CenterTitle(true);
h2_energy_vs_mult_with_M->GetYaxis()->SetLabelSize(0.045);
h2_energy_vs_mult_with_M->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 vs Charge 2 (TWIM Music)");
h2_charge1_charge_2 = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
h2_charge1_charge_2->GetXaxis()->SetTitle("Charge 1");
h2_charge1_charge_2->GetYaxis()->SetTitle("Charge 2");
h2_charge1_charge_2->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 - Charge 2 vs Charge_sum(TWIM Music)");
h2_charge1_charge_2_diff_vs_z_sum = new TH2D(hist_name,hist_name,1400,-70,70,1000,0,100);
h2_charge1_charge_2_diff_vs_z_sum->GetXaxis()->SetTitle("Charge 1 - Charge 2");
h2_charge1_charge_2_diff_vs_z_sum->GetYaxis()->SetTitle("Charge 1 + Charge 2");
h2_charge1_charge_2_diff_vs_z_sum->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_diff_vs_z_sum->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 vs Charge 2 (TWIM Music) after the minimal CALIFA cuts");
h2_charge1_charge_2_after = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
h2_charge1_charge_2_after->GetXaxis()->SetTitle("Charge 1");
h2_charge1_charge_2_after->GetYaxis()->SetTitle("Charge 2");
h2_charge1_charge_2_after->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_after->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_after->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_after->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 vs Charge 2 (TWIM Music) for TPat = SofStart");
h2_charge1_charge_2_tpat1 = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
h2_charge1_charge_2_tpat1->GetXaxis()->SetTitle("Charge 1");
h2_charge1_charge_2_tpat1->GetYaxis()->SetTitle("Charge 2");
h2_charge1_charge_2_tpat1->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_tpat1->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_tpat1->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_tpat1->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 vs Charge 2 (TWIM Music) for TPat =  Fission");
h2_charge1_charge_2_tpat2 = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
h2_charge1_charge_2_tpat2->GetXaxis()->SetTitle("Charge 1");
h2_charge1_charge_2_tpat2->GetYaxis()->SetTitle("Charge 2");
h2_charge1_charge_2_tpat2->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_tpat2->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_tpat2->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_tpat2->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 vs Charge 2 (TWIM Music) for TPat =  CalifOR (= S+C&&)");
h2_charge1_charge_2_tpat3 = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
h2_charge1_charge_2_tpat3->GetXaxis()->SetTitle("Charge 1");
h2_charge1_charge_2_tpat3->GetYaxis()->SetTitle("Charge 2");
h2_charge1_charge_2_tpat3->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_tpat3->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_tpat3->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_tpat3->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 vs Charge 2 (TWIM Music) for TPat =  CalifAND (= S+F+C&&)");
h2_charge1_charge_2_tpat4 = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
h2_charge1_charge_2_tpat4->GetXaxis()->SetTitle("Charge 1");
h2_charge1_charge_2_tpat4->GetYaxis()->SetTitle("Charge 2");
h2_charge1_charge_2_tpat4->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_tpat4->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_tpat4->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_tpat4->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 vs Charge 2 (TWIM Music) for one hit > 200 MeV in CALIFA");
h2_charge1_charge_2_one_CALIFA = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
h2_charge1_charge_2_one_CALIFA->GetXaxis()->SetTitle("Charge 1");
h2_charge1_charge_2_one_CALIFA->GetYaxis()->SetTitle("Charge 2");
h2_charge1_charge_2_one_CALIFA->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_one_CALIFA->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_one_CALIFA->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_one_CALIFA->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 vs Charge 2 (TWIM Music) for two hits > 200 MeV in CALIFA");
h2_charge1_charge_2_two_CALIFA = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
h2_charge1_charge_2_two_CALIFA->GetXaxis()->SetTitle("Charge 1");
h2_charge1_charge_2_two_CALIFA->GetYaxis()->SetTitle("Charge 2");
h2_charge1_charge_2_two_CALIFA->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_two_CALIFA->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_two_CALIFA->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_two_CALIFA->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 vs Charge 2 (TWIM Music) for three hits > 200 MeV in CALIFA");
h2_charge1_charge_2_three_CALIFA = new TH2D(hist_name,hist_name,500,0,100,500,0,100);
h2_charge1_charge_2_three_CALIFA->GetXaxis()->SetTitle("Charge 1");
h2_charge1_charge_2_three_CALIFA->GetYaxis()->SetTitle("Charge 2");
h2_charge1_charge_2_three_CALIFA->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_three_CALIFA->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_three_CALIFA->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_three_CALIFA->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 - Charge 2 vs Charge_sum(TWIM Music) for TPat = SofStart");
h2_charge1_charge_2_diff_vs_z_sum_tpat1 = new TH2D(hist_name,hist_name,1400,-70,70,1000,0,100);
h2_charge1_charge_2_diff_vs_z_sum_tpat1->GetXaxis()->SetTitle("Charge 1 - Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_tpat1->GetYaxis()->SetTitle("Charge 1 + Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_tpat1->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_tpat1->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_tpat1->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_diff_vs_z_sum_tpat1->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 - Charge 2 vs Charge_sum(TWIM Music) for TPat =  Fission");
h2_charge1_charge_2_diff_vs_z_sum_tpat2 = new TH2D(hist_name,hist_name,1400,-70,70,1000,0,100);
h2_charge1_charge_2_diff_vs_z_sum_tpat2->GetXaxis()->SetTitle("Charge 1 - Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_tpat2->GetYaxis()->SetTitle("Charge 1 + Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_tpat2->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_tpat2->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_tpat2->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_diff_vs_z_sum_tpat2->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 - Charge 2 vs Charge_sum(TWIM Music) for TPat =  CalifOR (= S+C&&)");
h2_charge1_charge_2_diff_vs_z_sum_tpat3 = new TH2D(hist_name,hist_name,1400,-70,70,1000,0,100);
h2_charge1_charge_2_diff_vs_z_sum_tpat3->GetXaxis()->SetTitle("Charge 1 - Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_tpat3->GetYaxis()->SetTitle("Charge 1 + Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_tpat3->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_tpat3->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_tpat3->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_diff_vs_z_sum_tpat3->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Charge 1 - Charge 2 vs Charge_sum(TWIM Music) for TPat =  CalifAND (= S+F+C&&)");
h2_charge1_charge_2_diff_vs_z_sum_tpat4 = new TH2D(hist_name,hist_name,1400,-70,70,1000,0,100);
h2_charge1_charge_2_diff_vs_z_sum_tpat4->GetXaxis()->SetTitle("Charge 1 - Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_tpat4->GetYaxis()->SetTitle("Charge 1 + Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_tpat4->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_tpat4->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_tpat4->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_diff_vs_z_sum_tpat4->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 - Charge 2 vs Charge_sum(TWIM Music) for one hit > 200MeV in CALIFA");
h2_charge1_charge_2_diff_vs_z_sum_one_CALIFA = new TH2D(hist_name,hist_name,1400,-70,70,1000,0,100);
h2_charge1_charge_2_diff_vs_z_sum_one_CALIFA->GetXaxis()->SetTitle("Charge 1 - Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_one_CALIFA->GetYaxis()->SetTitle("Charge 1 + Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_one_CALIFA->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_one_CALIFA->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_one_CALIFA->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_diff_vs_z_sum_one_CALIFA->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 - Charge 2 vs Charge_sum(TWIM Music) for two hits > 200MeV in CALIFA");
h2_charge1_charge_2_diff_vs_z_sum_two_CALIFA = new TH2D(hist_name,hist_name,1400,-70,70,1000,0,100);
h2_charge1_charge_2_diff_vs_z_sum_two_CALIFA->GetXaxis()->SetTitle("Charge 1 - Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_two_CALIFA->GetYaxis()->SetTitle("Charge 1 + Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_two_CALIFA->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_two_CALIFA->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_two_CALIFA->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_diff_vs_z_sum_two_CALIFA->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Charge 1 - Charge 2 vs Charge_sum(TWIM Music) for three hits > 200MeV in CALIFA");
h2_charge1_charge_2_diff_vs_z_sum_three_CALIFA = new TH2D(hist_name,hist_name,1400,-70,70,1000,0,100);
h2_charge1_charge_2_diff_vs_z_sum_three_CALIFA->GetXaxis()->SetTitle("Charge 1 - Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_three_CALIFA->GetYaxis()->SetTitle("Charge 1 + Charge 2");
h2_charge1_charge_2_diff_vs_z_sum_three_CALIFA->GetXaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_three_CALIFA->GetYaxis()->CenterTitle(true);
h2_charge1_charge_2_diff_vs_z_sum_three_CALIFA->GetYaxis()->SetLabelSize(0.045);
h2_charge1_charge_2_diff_vs_z_sum_three_CALIFA->GetYaxis()->SetTitleSize(0.045);



for (Int_t i = 0; i < 20 ; i++){
sprintf(hist_name,"Energy distribution over theta_phi for p2p eventnr %i",i);
h2_energy_angular_dist[i] = new TH2D(hist_name,hist_name,52,22.15,152.15,60,-180,180);
h2_energy_angular_dist[i]->GetXaxis()->SetTitle("Theta [degr]");
h2_energy_angular_dist[i]->GetYaxis()->SetTitle("Phi [degr]");
h2_energy_angular_dist[i]->GetXaxis()->CenterTitle(true);
h2_energy_angular_dist[i]->GetYaxis()->CenterTitle(true);
h2_energy_angular_dist[i]->GetXaxis()->SetLabelSize(0.045);
h2_energy_angular_dist[i]->GetYaxis()->SetLabelSize(0.045);
}


for (Int_t i = 0; i < 20 ; i++){
	sprintf(hist_name,"Energy distribution over theta_phi for p2p eventnr  %i with one hit > 100 MeV",i);
	h2_energy_angular_dist_one[i] = new TH2D(hist_name,hist_name,52,22.15,152.15,60,-180,180);
	h2_energy_angular_dist_one[i]->GetXaxis()->SetTitle("Theta [degr]");
	h2_energy_angular_dist_one[i]->GetYaxis()->SetTitle("Phi [degr]");
	h2_energy_angular_dist_one[i]->GetXaxis()->CenterTitle(true);
	h2_energy_angular_dist_one[i]->GetYaxis()->CenterTitle(true);
	h2_energy_angular_dist_one[i]->GetXaxis()->SetLabelSize(0.045);
	h2_energy_angular_dist_one[i]->GetYaxis()->SetLabelSize(0.045);
}

for (Int_t i = 0; i < 20 ; i++){
	sprintf(hist_name,"Energy distribution over theta_phi for p2p eventnr %i with two hits > 100 MeV",i);
	h2_energy_angular_dist_two[i] = new TH2D(hist_name,hist_name,52,22.15,152.15,60,-180,180);
	h2_energy_angular_dist_two[i]->GetXaxis()->SetTitle("Theta [degr]");
	h2_energy_angular_dist_two[i]->GetYaxis()->SetTitle("Phi [degr]");
	h2_energy_angular_dist_two[i]->GetXaxis()->CenterTitle(true);
	h2_energy_angular_dist_two[i]->GetYaxis()->CenterTitle(true);
	h2_energy_angular_dist_two[i]->GetXaxis()->SetLabelSize(0.045);
	h2_energy_angular_dist_two[i]->GetYaxis()->SetLabelSize(0.045);
}

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

sprintf(hist_name, "Theta 1 vs Theta 2 for the two hits > 100 MeV");
h2_theta1_vs_theta2 = new TH2D(hist_name,hist_name,25,22.15,84.65,25,22.15,84.65);
h2_theta1_vs_theta2->GetXaxis()->SetTitle("Theta1 [degr]");
h2_theta1_vs_theta2->GetYaxis()->SetTitle("Theta2 [degr]");
h2_theta1_vs_theta2->GetXaxis()->CenterTitle(true);
h2_theta1_vs_theta2->GetYaxis()->CenterTitle(true);
h2_theta1_vs_theta2->GetYaxis()->SetLabelSize(0.045);
h2_theta1_vs_theta2->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta 1 vs Theta 2 for the two hits with highest energy");
h2_theta1_vs_theta2_cluster = new TH2D(hist_name,hist_name,70,20,90,70,20,90);
h2_theta1_vs_theta2_cluster->GetXaxis()->SetTitle("Theta1 [degr]");
h2_theta1_vs_theta2_cluster->GetYaxis()->SetTitle("Theta2 [degr]");
h2_theta1_vs_theta2_cluster->GetXaxis()->CenterTitle(true);
h2_theta1_vs_theta2_cluster->GetYaxis()->CenterTitle(true);
h2_theta1_vs_theta2_cluster->GetYaxis()->SetLabelSize(0.045);
h2_theta1_vs_theta2_cluster->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "Theta 1 vs Theta 2 for the two hits with highest energy (> 30MeV) and best dphi");
h2_theta1_vs_theta2_cluster_best_dphi = new TH2D(hist_name,hist_name,70,20,90,70,20,90);
h2_theta1_vs_theta2_cluster_best_dphi->GetXaxis()->SetTitle("Theta1 [degr]");
h2_theta1_vs_theta2_cluster_best_dphi->GetYaxis()->SetTitle("Theta2 [degr]");
h2_theta1_vs_theta2_cluster_best_dphi->GetXaxis()->CenterTitle(true);
h2_theta1_vs_theta2_cluster_best_dphi->GetYaxis()->CenterTitle(true);
h2_theta1_vs_theta2_cluster_best_dphi->GetYaxis()->SetLabelSize(0.045);
h2_theta1_vs_theta2_cluster_best_dphi->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Energy 1 vs Energy 2 for the two hits with highest energy");
h2_e1_vs_e2_cluster = new TH2D(hist_name,hist_name,500,0,1000,500,0,1000);
h2_e1_vs_e2_cluster->GetXaxis()->SetTitle("Energy1 [MeV]");
h2_e1_vs_e2_cluster->GetYaxis()->SetTitle("Energy2 [MeV]");
h2_e1_vs_e2_cluster->GetXaxis()->CenterTitle(true);
h2_e1_vs_e2_cluster->GetYaxis()->CenterTitle(true);
h2_e1_vs_e2_cluster->GetYaxis()->SetLabelSize(0.045);
h2_e1_vs_e2_cluster->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Energy 1 vs Energy 2 for the two hits with highest energy (> 30MeV) and best dphi");
h2_e1_vs_e2_cluster_best_dphi = new TH2D(hist_name,hist_name,500,0,1000,500,0,1000);
h2_e1_vs_e2_cluster_best_dphi->GetXaxis()->SetTitle("Energy1 [MeV]");
h2_e1_vs_e2_cluster_best_dphi->GetYaxis()->SetTitle("Energy2 [MeV]");
h2_e1_vs_e2_cluster_best_dphi->GetXaxis()->CenterTitle(true);
h2_e1_vs_e2_cluster_best_dphi->GetYaxis()->CenterTitle(true);
h2_e1_vs_e2_cluster_best_dphi->GetYaxis()->SetLabelSize(0.045);
h2_e1_vs_e2_cluster_best_dphi->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Energy 1 vs Energy 2 for the two hits with highest energy (> 30MeV) and dphi = 180 +-15");
h2_e1_vs_e2_cluster_two_highest = new TH2D(hist_name,hist_name,500,0,1000,500,0,1000);
h2_e1_vs_e2_cluster_two_highest->GetXaxis()->SetTitle("Energy1 [MeV]");
h2_e1_vs_e2_cluster_two_highest->GetYaxis()->SetTitle("Energy2 [MeV]");
h2_e1_vs_e2_cluster_two_highest->GetXaxis()->CenterTitle(true);
h2_e1_vs_e2_cluster_two_highest->GetYaxis()->CenterTitle(true);
h2_e1_vs_e2_cluster_two_highest->GetYaxis()->SetLabelSize(0.045);
h2_e1_vs_e2_cluster_two_highest->GetYaxis()->SetTitleSize(0.045);


sprintf(hist_name, "Opening Angle for two hits (> 30 MeV) and best dphi vs E1+E2");
h2_opening_angle_vs_esum = new TH2D(hist_name,hist_name,52,22.15,152.15,100,0,1000);
h2_opening_angle_vs_esum->GetXaxis()->SetTitle("Opening Angle [degr]");
h2_opening_angle_vs_esum->GetYaxis()->SetTitle("E1+E2 [MeV]");
h2_opening_angle_vs_esum->GetXaxis()->CenterTitle(true);
h2_opening_angle_vs_esum->GetYaxis()->CenterTitle(true);
h2_opening_angle_vs_esum->GetYaxis()->SetLabelSize(0.045);
h2_opening_angle_vs_esum->GetYaxis()->SetTitleSize(0.045);

sprintf(hist_name, "E1 or E2 vs angular distribution");
h2_e_vs_phi = new TH2D(hist_name,hist_name,52,22.15,152.15,60,-180,180);
h2_e_vs_phi->GetXaxis()->SetTitle("Theta [degr]");
h2_e_vs_phi->GetYaxis()->SetTitle("Phi [degr]");
h2_e_vs_phi->GetXaxis()->CenterTitle(true);
h2_e_vs_phi->GetYaxis()->CenterTitle(true);
h2_e_vs_phi->GetYaxis()->SetLabelSize(0.045);
h2_e_vs_phi->GetYaxis()->SetTitleSize(0.045);
//end of 2D_histos-----------------
//}
//void write_to_file(char* output_file){
//Add 1D  histos to vector -> write to file
list_of_1D_histos.push_back(h1_trigger);
list_of_1D_histos.push_back(h1_wr_ts_califa_diff_WRM_messel_rand);
list_of_1D_histos.push_back(h1_wr_ts_califa_diff_WRM_wix_rand);
list_of_1D_histos.push_back(h1_wr_ts_califa_diff_WRM_both_rand);
list_of_1D_histos.push_back(h1_wr_ts_califa_vs_E_M);
list_of_1D_histos.push_back(h1_wr_ts_califa_vs_E_noM);
list_of_1D_histos.push_back(h1_mult_with_M);
list_of_1D_histos.push_back(h1_mult_no_M);
list_of_1D_histos.push_back(h1_mult_e_high);
list_of_1D_histos.push_back(h1_mult_e_high_noM);
list_of_1D_histos.push_back(h1_charge1_plus_charge2);
list_of_1D_histos.push_back(h1_daughter_nuclei_charge);
list_of_1D_histos.push_back(h1_charge_sum_after);
list_of_1D_histos.push_back(h1_charge1_plus_charge2_tpat1);
list_of_1D_histos.push_back(h1_charge1_plus_charge2_tpat2);
list_of_1D_histos.push_back(h1_charge1_plus_charge2_tpat3);
list_of_1D_histos.push_back(h1_charge1_plus_charge2_tpat4);
list_of_1D_histos.push_back(h1_charge1_plus_charge2_one_CALIFA);
list_of_1D_histos.push_back(h1_charge1_plus_charge2_two_CALIFA);
list_of_1D_histos.push_back(h1_charge1_plus_charge2_three_CALIFA);
list_of_1D_histos.push_back(h1_mult_100mev);
list_of_1D_histos.push_back(h1_mult_100mev_two);
list_of_1D_histos.push_back(h1_theta_sum_p2p);
list_of_1D_histos.push_back(h1_theta_sum_z_sum_915_925);
list_of_1D_histos.push_back(h1_theta_sum_z_sum_915_925_phi_cut);
list_of_1D_histos.push_back(h1_theta_sum_z_sum_915_925_phi_ener_cut);
list_of_1D_histos.push_back(h1_theta_sum_z_sum_902_914);
list_of_1D_histos.push_back(h1_theta_sum_z_sum_902_914_phi_cut);
list_of_1D_histos.push_back(h1_theta_sum_z_sum_891_901);
list_of_1D_histos.push_back(h1_theta_sum_z_sum_902_914_phi_ener_cut);
list_of_1D_histos.push_back(h1_theta_sum_z_sum_879_89);
list_of_1D_histos.push_back(h1_theta_sum_z_sum_891_901_phi_cut);
list_of_1D_histos.push_back(h1_theta_sum_z_sum_879_89_phi_ener_cut);
list_of_1D_histos.push_back(h1_theta_sum_z_sum_891_901_phi_ener_cut);
list_of_1D_histos.push_back(h1_theta_sum_z_sum_868_879);
list_of_1D_histos.push_back(h1_theta_sum_z_sum_879_89_phi_cut);
list_of_1D_histos.push_back(h1_theta_sum_z_sum_868_879_phi_cut);
list_of_1D_histos.push_back(h1_theta_sum_z_sum_868_879_phi_ener_cut);
list_of_1D_histos.push_back(h1_theta_sum_z_1_2_cut1);
list_of_1D_histos.push_back(h1_theta_sum_z_1_2_cut2);
list_of_1D_histos.push_back(h1_theta_sum_z_1_2_cut3);
list_of_1D_histos.push_back(h1_theta_sum_z_1_2_cut4);
list_of_1D_histos.push_back(h1_theta_sum_z_1_2_cut5);
list_of_1D_histos.push_back(h1_one_hit_twim_charge);
list_of_1D_histos.push_back(h1_one_hit_twim_charge_tpat1);
list_of_1D_histos.push_back(h1_theta_sum_p2p_best_dphi);
list_of_1D_histos.push_back(h1_opening_angle);
list_of_1D_histos.push_back(h1_opening_angle_first_two);
list_of_1D_histos.push_back(h1_opening_angle_one_other);
list_of_1D_histos.push_back(h1_opening_angle_two_other);
list_of_1D_histos.push_back(h1_opening_angle_proton);
list_of_1D_histos.push_back(h1_angle_check);
list_of_1D_histos.push_back(h1_angle_check_arzi);



//Add 2D  histos to vector -> write to file

list_of_2D_histos.push_back(h2_energy_vs_theta);
list_of_2D_histos.push_back(h2_phi_vs_theta);
list_of_2D_histos.push_back(h2_energy_vs_mult_no_M);
list_of_2D_histos.push_back(h2_energy_vs_mult_with_M);
list_of_2D_histos.push_back(h2_charge1_charge_2);
list_of_2D_histos.push_back(h2_charge1_charge_2_diff_vs_z_sum);
list_of_2D_histos.push_back(h2_charge1_charge_2_after);
list_of_2D_histos.push_back(h2_charge1_charge_2_tpat1);
list_of_2D_histos.push_back(h2_charge1_charge_2_tpat2);
list_of_2D_histos.push_back(h2_charge1_charge_2_tpat3);
list_of_2D_histos.push_back(h2_charge1_charge_2_tpat4);
list_of_2D_histos.push_back(h2_charge1_charge_2_one_CALIFA);
list_of_2D_histos.push_back(h2_charge1_charge_2_two_CALIFA);
list_of_2D_histos.push_back(h2_charge1_charge_2_three_CALIFA);
list_of_2D_histos.push_back(h2_charge1_charge_2_diff_vs_z_sum_tpat1);
list_of_2D_histos.push_back(h2_charge1_charge_2_diff_vs_z_sum_two_CALIFA);
list_of_2D_histos.push_back(h2_charge1_charge_2_diff_vs_z_sum_one_CALIFA);
list_of_2D_histos.push_back(h2_charge1_charge_2_diff_vs_z_sum_tpat4);
list_of_2D_histos.push_back(h2_charge1_charge_2_diff_vs_z_sum_tpat3);
list_of_2D_histos.push_back(h2_charge1_charge_2_diff_vs_z_sum_tpat2);
list_of_2D_histos.push_back(h2_charge1_charge_2_diff_vs_z_sum_three_CALIFA);
//list_of_2D_histos.push_back(h2_energy_angular_dist[20]);
//list_of_2D_histos.push_back(h2_energy_angular_dist_one[20]);
//list_of_2D_histos.push_back(h2_energy_angular_dist_two[20]);
//list_of_2D_histos.push_back(h2_charge1_vs_charge2_energy_cut[6]);
//list_of_2D_histos.push_back(h2_charge1_vs_charge2_theta_cut[6]);
//list_of_2D_histos.push_back(h2_charge1_vs_charge2_phi_cut[6]);
list_of_2D_histos.push_back(h2_theta1_vs_theta2);
list_of_2D_histos.push_back(h2_theta1_vs_theta2_cluster);
list_of_2D_histos.push_back(h2_theta1_vs_theta2_cluster_best_dphi);
list_of_2D_histos.push_back(h2_e1_vs_e2_cluster);
list_of_2D_histos.push_back(h2_e1_vs_e2_cluster_best_dphi);
list_of_2D_histos.push_back(h2_e1_vs_e2_cluster_two_highest);
list_of_2D_histos.push_back(h2_opening_angle_vs_esum);
list_of_2D_histos.push_back(h2_e_vs_phi);

}
void write_to_file(char* output_file){
	cout << "number of events in histo:\t" << h2_charge1_charge_2->GetEntries() << endl;
	h2_charge1_charge_2->Draw();
	TFile *f = new TFile(output_file,"RECREATE");
	TList *l = new TList();
	for (Int_t i = 0; i < list_of_1D_histos.size(); i++){
		l->Add(list_of_1D_histos[i]);
		}

	for (Int_t i = 0; i < list_of_2D_histos.size(); i++){
		l->Add(list_of_2D_histos[i]);
		}
	l->Write("histlist", TObject::kSingleKey);
	}
