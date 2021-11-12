// ------------------------------------------------------------------------------------
// -----  Macro to correct the energy of the TWIM by TOF and calibrate charge  -----
// -----               Created 25/10/2021 by Antia Grana Gonzalez              -----
// ------------------------------------------------------------------------------------

// ROOT headers
#include "TClonesArray.h"
#include "TDecompSVD.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TRandom.h"
#include "TVectorD.h"
R3BTGeoPar* targetPar; // Variables globales
R3BSofTwimCalPar* Cal_Par;
R3BSofGladFieldPar* fGladPar;
const Double_t c = 29.9792458;

void SetParContainers()
{
  // Containers ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  FairRuntimeDb* rtdb = FairRuntimeDb::instance();
  Cal_Par = (R3BSofTwimCalPar*)(rtdb->getContainer("twimCalPar"));
  //R3BTGeoPar* fMw1GeoPar = (R3BTGeoPar*)rtdb->getContainer("Mwpc1GeoPar");
  //R3BTGeoPar* fMw2GeoPar = (R3BTGeoPar*)rtdb->getContainer("Mwpc2GeoPar");
  //R3BTGeoPar* fMw3GeoPar = (R3BTGeoPar*)rtdb->getContainer("Mwpc3GeoPar");
  //R3BTGeoPar* ftofwPar = (R3BTGeoPar*)rtdb->getContainer("TofwGeoPar");
  targetPar = (R3BTGeoPar*)rtdb->getContainer("TargetGeoPar");
  fGladPar = (R3BSofGladFieldPar*)rtdb->getContainer("GladFieldPar");
  //Bool_t kParameterMerged = kFALSE;
  FairParAsciiFileIo* parIo2 = new FairParAsciiFileIo(); // Root file
  parIo2->open("/u/land/antiagg/R3BRoot/sofia/macros/s455Up2p/parameters/CalibParam.par", "in");
  rtdb->setFirstInput(parIo2);
  //rtdb->print();
  UInt_t runId = 1;
  rtdb->initContainers(runId);
  return;

}


void TWIM_calib_TOF_corr()
{
  // Graphs

  // E vs TOF
  TH2F* E_LD_vs_tof0 = new TH2F("E_LD_vs_tof0","E_LD_vs_tof0", 1500, 0, 50, 1500, 0, 90000);
  E_LD_vs_tof0->GetXaxis()->SetTitle("Tof0 [ns]");
  E_LD_vs_tof0->GetYaxis()->SetTitle("Energy [KeV]");

  TH2F* E_RU_vs_tof1 = new TH2F("E_RU_vs_tof1","E_RU_vs_tof1", 1500, 0, 50, 1500, 0, 90000);
  E_RU_vs_tof1->GetXaxis()->SetTitle("Tof1 [ns]");
  E_RU_vs_tof1->GetYaxis()->SetTitle("Energy [KeV]");

  TH2F* E_LU_vs_tof0 = new TH2F("E_LU_vs_tof0","E_LU_vs_tof0", 1500, 0, 50, 1500, 0, 90000);
  E_LU_vs_tof0->GetXaxis()->SetTitle("Tof0 [ns]");
  E_LU_vs_tof0->GetYaxis()->SetTitle("Energy [KeV]");

  TH2F* E_RD_vs_tof1 = new TH2F("E_RD_vs_tof1","E_RD_vs_tof1", 1500, 0, 50, 1500, 0, 90000);
  E_RD_vs_tof1->GetXaxis()->SetTitle("Tof1 [ns]");
  E_RD_vs_tof1->GetYaxis()->SetTitle("Energy [KeV]");

  // E vs TOF CORR
  TH2F* E_LD_vs_tof0_corr = new TH2F("E_LD_vs_tof0_corr","E_LD_vs_tof0_corr", 1500, 0, 50, 1500, 0, 90000);
  E_LD_vs_tof0_corr->GetXaxis()->SetTitle("Tof0 [ns]");
  E_LD_vs_tof0_corr->GetYaxis()->SetTitle("Energy [KeV]");

  TH2F* E_RU_vs_tof1_corr = new TH2F("E_RU_vs_tof1_corr","E_RU_vs_tof1_corr", 1500, 0, 50, 1500, 0, 90000);
  E_RU_vs_tof1_corr->GetXaxis()->SetTitle("Tof1 [ns]");
  E_RU_vs_tof1_corr->GetYaxis()->SetTitle("Energy [KeV]");

  TH2F* E_LU_vs_tof0_corr = new TH2F("E_LU_vs_tof0_corr","E_LU_vs_tof0_corr", 1500, 0, 50, 1500, 0, 90000);
  E_LU_vs_tof0_corr->GetXaxis()->SetTitle("Tof0 [ns]");
  E_LU_vs_tof0_corr->GetYaxis()->SetTitle("Energy [KeV]");

  TH2F* E_RD_vs_tof1_corr = new TH2F("E_RD_vs_tof1_corr","E_RD_vs_tof1_corr", 1500, 0, 50, 1500, 0, 90000);
  E_RD_vs_tof1_corr->GetXaxis()->SetTitle("Tof1 [ns]");
  E_RD_vs_tof1_corr->GetYaxis()->SetTitle("Energy [KeV]");

  // E CORR
  TH1F* Z_LD_pre = new TH1F("Z_LD_pre","Z_LD_pre", 1000, 0, 300);
  Z_LD_pre->GetXaxis()->SetTitle("Energy [KeV]");
  Z_LD_pre->GetYaxis()->SetTitle("Counts");

  TH1F* Z_RU_pre = new TH1F("Z_RU_pre","Z_RU_pre", 1000, 0, 300);
  Z_RU_pre->GetXaxis()->SetTitle("Energy [KeV]");
  Z_RU_pre->GetYaxis()->SetTitle("Counts");

  TH1F* Z_LU_pre = new TH1F("Z_LU_pre","Z_LU_pre", 1000, 0, 300);
  Z_LU_pre->GetXaxis()->SetTitle("Energy [KeV]");
  Z_LU_pre->GetYaxis()->SetTitle("Counts");

  TH1F* Z_RD_pre = new TH1F("Z_RD_pre","Z_RD_pre", 1000, 0, 300);
  Z_RD_pre->GetXaxis()->SetTitle("Energy [KeV]");
  Z_RD_pre->GetYaxis()->SetTitle("Counts");

  // CALIBRATION Z CHANNELS TO Z
  TH2F* Z_LD_calibrating = new TH2F("Z_LD_calibrating","Z_LD_calibrating", 1000, 0, 300, 1000, 30, 80);
  Z_LD_calibrating->GetXaxis()->SetTitle("Z [channels]");
  Z_LD_calibrating->GetYaxis()->SetTitle("Z");

  TH2F* Z_RU_calibrating = new TH2F("Z_RU_calibrating","Z_RU_calibrating", 1000, 0, 300, 1000, 30, 80);
  Z_RU_calibrating->GetXaxis()->SetTitle("Z [channels]");
  Z_RU_calibrating->GetYaxis()->SetTitle("Z");

  TH2F* Z_LU_calibrating = new TH2F("Z_LU_calibrating","Z_LU_calibrating", 1000, 0, 300, 1000, 30, 80);
  Z_LU_calibrating->GetXaxis()->SetTitle("Z [channels]");
  Z_LU_calibrating->GetYaxis()->SetTitle("Z");

  TH2F* Z_RD_calibrating = new TH2F("Z_RD_calibrating","Z_RD_calibrating", 1000, 0, 300, 1000, 30, 80);
  Z_RD_calibrating->GetXaxis()->SetTitle("Z [channels]");
  Z_RD_calibrating->GetYaxis()->SetTitle("Z");

  // Z RAW
  TH1F* Z_LD_raw_graph = new TH1F("Z_LD_raw_graph","Z_LD_raw_graph", 1000, 20, 80);
  Z_LD_raw_graph->GetXaxis()->SetTitle("Z LD raw");
  Z_LD_raw_graph->GetYaxis()->SetTitle("Counts");

  TH1F* Z_RU_raw_graph = new TH1F("Z_RU_raw_graph","Z_RU_raw_graph", 1000, 20, 80);
  Z_RU_raw_graph->GetXaxis()->SetTitle("Z RU raw");
  Z_RU_raw_graph->GetYaxis()->SetTitle("Counts");

  TH1F* Z_LU_raw_graph = new TH1F("Z_LU_raw_graph","Z_LU_raw_graph", 1000, 20, 80);
  Z_LU_raw_graph->GetXaxis()->SetTitle("Z LU raw");
  Z_LU_raw_graph->GetYaxis()->SetTitle("Counts");

  TH1F* Z_RD_raw_graph = new TH1F("Z_RD_raw_graph","Z_RD_raw_graph", 1000, 20, 80);
  Z_RD_raw_graph->GetXaxis()->SetTitle("Z RD raw");
  Z_RD_raw_graph->GetYaxis()->SetTitle("Counts");

  // Z RAW PART 2
  TH2F* Z_LD_vs_RU_raw = new TH2F("Z_LD_vs_RU_raw","Z_LD_vs_RU_raw", 1000, 20, 80, 1000, 20, 80);
  Z_LD_vs_RU_raw->GetXaxis()->SetTitle("Z RU raw");
  Z_LD_vs_RU_raw->GetYaxis()->SetTitle("Z LD raw");

  TH1F* Zsum_LD_RU_raw = new TH1F("Zsum_LD_RU_raw","Zsum_LD_RU_raw", 1000, 85, 120);
  Zsum_LD_RU_raw->GetXaxis()->SetTitle("Z_LD+Z_RU raw");
  Zsum_LD_RU_raw->GetYaxis()->SetTitle("Counts");

  TH2F* Z_LU_vs_RD_raw = new TH2F("Z_LU_vs_RD_raw","Z_LU_vs_RD_raw", 1000, 20, 80, 1000, 20, 80);
  Z_LU_vs_RD_raw->GetXaxis()->SetTitle("Z RD raw");
  Z_LU_vs_RD_raw->GetYaxis()->SetTitle("Z LU raw");

  TH1F* Zsum_LU_RD_raw = new TH1F("Zsum_LU_RD_raw","Zsum_LU_RD_raw", 1000, 85, 120);
  Zsum_LU_RD_raw->GetXaxis()->SetTitle("Z_LU+Z_RD raw");
  Zsum_LU_RD_raw->GetYaxis()->SetTitle("Counts");

  // Z RAW FISSION FRAGMENTS
  TH2F* Z_ff_raw = new TH2F("Z_ff_raw","Z_ff_raw", 1000, 20, 80, 1000, 20, 80);
  Z_ff_raw->GetXaxis()->SetTitle("Z ff1 raw");
  Z_ff_raw->GetYaxis()->SetTitle("Z ff2 raw");

  // Z  RAW FISSIONING SYSTEM
  TH1F* Z_f_raw = new TH1F("Z_f_raw","Z_f_raw", 1000, 85, 105);
  Z_f_raw->GetXaxis()->SetTitle("Z_f raw");
  Z_f_raw->GetYaxis()->SetTitle("Counts");

  // Z GRAPHS
  TH1F* Z_LD_graph = new TH1F("Z_LD_graph","Z_LD_graph", 1000, 20, 80);
  Z_LD_graph->GetXaxis()->SetTitle("Z LD");
  Z_LD_graph->GetYaxis()->SetTitle("Counts");

  TH1F* Z_RU_graph = new TH1F("Z_RU_graph","Z_RU_graph", 1000, 20, 80);
  Z_RU_graph->GetXaxis()->SetTitle("Z RU");
  Z_RU_graph->GetYaxis()->SetTitle("Counts");

  TH1F* Z_LU_graph = new TH1F("Z_LU_graph","Z_LU_graph", 1000, 20, 80);
  Z_LU_graph->GetXaxis()->SetTitle("Z LU");
  Z_LU_graph->GetYaxis()->SetTitle("Counts");

  TH1F* Z_RD_graph = new TH1F("Z_RD_graph","Z_RD_graph", 1000, 20, 80);
  Z_RD_graph->GetXaxis()->SetTitle("Z RD");
  Z_RD_graph->GetYaxis()->SetTitle("Counts");

  // Z GRAPHS PART 2
  TH2F* Z_LD_vs_RU = new TH2F("Z_LD_vs_RU","Z_LD_vs_RU", 1000, 20, 80, 1000, 20, 80);
  Z_LD_vs_RU->GetXaxis()->SetTitle("Z RU");
  Z_LD_vs_RU->GetYaxis()->SetTitle("Z LD");

  TH1F* Zsum_LD_RU = new TH1F("Zsum_LD_RU","Zsum_LD_RU", 1000, 85, 120);
  Zsum_LD_RU->GetXaxis()->SetTitle("Z_LD+Z_RU");
  Zsum_LD_RU->GetYaxis()->SetTitle("Counts");

  TH2F* Z_LU_vs_RD = new TH2F("Z_LU_vs_RD","Z_LU_vs_RD", 1000, 20, 80, 1000, 20, 80);
  Z_LU_vs_RD->GetXaxis()->SetTitle("Z RD");
  Z_LU_vs_RD->GetYaxis()->SetTitle("Z LU");

  TH1F* Zsum_LU_RD = new TH1F("Zsum_LU_RD","Zsum_LU_RD", 1000, 85, 120);
  Zsum_LU_RD->GetXaxis()->SetTitle("Z_LU+Z_RD");
  Zsum_LU_RD->GetYaxis()->SetTitle("Counts");

  // Z FISSION FRAGMENTS
  TH2F* Z_ff = new TH2F("Z_ff","Z_ff", 1000, 20, 80, 1000, 20, 80);
  Z_ff->GetXaxis()->SetTitle("Z ff1");
  Z_ff->GetYaxis()->SetTitle("Z ff2");

  // Z FISSIONING SYSTEM
  TH1F* Z_f = new TH1F("Z_f","Z_f", 1000, 85, 105);
  Z_f->GetXaxis()->SetTitle("Z_f");
  Z_f->GetYaxis()->SetTitle("Counts");


  // INPUT DATA /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  char title1[400] = "273_mw_tof_hitdata.root";
  //char title1[400] = "cal_data_main0273_all_2hitsmwpc_FULL_STAT_68.root";
  TFile* file1 = TFile::Open(title1);
  TTree* tree = (TTree*)file1->Get("evt");
  Float_t nevents = tree->GetEntries();
  std::cout << "Events: " << nevents << std::endl;

  TClonesArray* TwimCalDataCA = new TClonesArray("R3BSofTwimCalData", 10);
  TBranch* branchTwimCalData = tree->GetBranch("TwimCalData");
  branchTwimCalData->SetAddress(&TwimCalDataCA);

  TClonesArray* fTofWHitDataCA = new TClonesArray("R3BSofTofWHitData", 10);
  TBranch* branchTofWHitData = tree->GetBranch("TofWHitData");
  branchTofWHitData->SetAddress(&fTofWHitDataCA);

  Int_t NumSec = 4;
  Int_t NumAnodes = 16; // Porque hay 16 anodos normales + 1 anodo de referencia
  Int_t NumAnodesAngleFit = 0;
  Int_t NumParams = 2;

  SetParContainers();

  // LOOP IN THE EVENTS--------------------------------------------------------------
  for (Int_t y = 0; y < nevents; y++)
  {
    if (y % 1000000 == 0)
        std::cout << "Event:" << y << "\n";
        //cout << "Event = " << y << endl;

    tree->GetEntry(y);

    Int_t nHitTwim = TwimCalDataCA->GetEntries();
    Int_t nHitTofW = fTofWHitDataCA->GetEntries();

    R3BSofTwimCalData** CalDat = new R3BSofTwimCalData*[nHitTwim];
    R3BSofTofWHitData** HitTofW = new R3BSofTofWHitData*[nHitTofW];

    Float_t Esum_LD = 0.;
    Float_t Esum_LU = 0.;
    Float_t Esum_RU = 0.;
    Float_t Esum_RD = 0.;
    Float_t Esum_mean_LD = 0.;
    Float_t Esum_mean_LU = 0.;
    Float_t Esum_mean_RU = 0.;
    Float_t Esum_mean_RD = 0.;
    Float_t E_anodo_LD = 0.;
    Float_t E_anodo_LU = 0.;
    Float_t dt_LD = 0.;
    Float_t dt_mean_LD = 0;
    Float_t VD_mean = 0;
    Double_t E_LD_corr = 0;
    Double_t E_RU_corr = 0;
    Double_t E_LU_corr = 0;
    Double_t E_RD_corr = 0;
    Double_t Z_LD_raw = 0;
    Double_t Z_RU_raw = 0;
    Double_t Z_LU_raw = 0;
    Double_t Z_RD_raw = 0;
    Double_t Z_LD = 0;
    Double_t Z_RU = 0;
    Double_t Z_LU = 0;
    Double_t Z_RD = 0;
    Float_t p0_LD_origen = 15300;
    Float_t p0_LD = -12062.1;
    Float_t p1_LD = 748.278;
    Float_t p2_LD = 0;
    Float_t p0_RU_origen = 15300;
    Float_t p0_RU = -12411;
    Float_t p1_RU = 762.646;
    Float_t p2_RU = 0;
    Float_t p0_LU_origen = 16500;
    Float_t p0_LU = -13609.4;
    Float_t p1_LU = 823.864;
    Float_t p2_LU = 0;
    Float_t p0_RD_origen = 16500;
    Float_t p0_RD = -12702;
    Float_t p1_RD = 779.28;
    Float_t p2_RD = 0;
    Double_t p0_LD_calZ = 4.47155;//4.6023;
    Double_t p1_LD_calZ = 0.360702;//0.359404;
    Double_t p0_RU_calZ = 4.43048;//4.54966;
    Double_t p1_RU_calZ = 0.361024;//0.359826;
    Double_t p0_LU_calZ = 3.66836;//3.71056;
    Double_t p1_LU_calZ = 0.35354;//0.353021;
    Double_t p0_RD_calZ = 1.90584;//1.9143;
    Double_t p1_RD_calZ = 0.367022;//0.366784;


    Int_t secId = 0, anodeId = 0;
    Double_t energyperanode[NumSec][NumAnodes];
    Double_t dt[NumSec][NumAnodes];
    Double_t good_dt[NumAnodes];
    Int_t nbdet = 0;
    Int_t StatusAnodes[NumSec][NumAnodes];

    Int_t pdid[2];
    Double_t tof[2];
    Double_t zf[2];

    for (Int_t j = 0; j < 2; j++)
    {
        pdid[j] = 0;
        tof[j] = 0.;
        zf[j] = 0.;
    }
    for (Int_t i = 0; i < nHitTofW; i++)
    {

        HitTofW[i] = (R3BSofTofWHitData*)(fTofWHitDataCA->At(i));
        if (i == 0)
        {
            pdid[0] = HitTofW[i]->GetPaddle();
            tof[0] = HitTofW[i]->GetTof();
        }
        else
        {
            if (HitTofW[i]->GetPaddle() > pdid[0])
            {
                pdid[1] = HitTofW[i]->GetPaddle();
                tof[1] = HitTofW[i]->GetTof();
            }
            else
            {
                pdid[1] = pdid[0];  //[0]->izquierda, [1]->derecha
                tof[1] = tof[0];
                pdid[0] = HitTofW[i]->GetPaddle();
                tof[0] = HitTofW[i]->GetTof();
            }
        }
    }


// TWIM

    for (Int_t s = 0; s < NumSec; s++)
    {
        for (Int_t i = 0; i < NumAnodes; i++)
        {
            StatusAnodes[s][i] = Cal_Par->GetInUse(s + 1, i + 1);
        }
    }
    for (Int_t j = 0; j < NumAnodes; j++)
    {
        good_dt[j] = 0.;
        for (Int_t i = 0; i < NumSec; i++)
        {
            energyperanode[i][j] = 0.;
            dt[i][j] = 0.;
        }
    }

    for (Int_t i = 0; i < nHitTwim; i++)
    {
        CalDat[i] = (R3BSofTwimCalData*)(TwimCalDataCA->At(i));
        secId = CalDat[i]->GetSecID();
        anodeId = CalDat[i]->GetAnodeID();
        if (energyperanode[secId][anodeId] == 0) // para que la multiplicidad=1
        {
            energyperanode[secId][anodeId] = CalDat[i]->GetEnergy(); //
            dt[secId][anodeId] = CalDat[i]->GetDTime();
            //cout << "hit = " << i << ", secId = " << secId << ", anodeId = " << anodeId << ", Energy = " << CalDat[i]->GetEnergy() << std::endl;
        }
    }

    // calculate truncated dE from 16 anodes, Twim-MUSIC
    for (Int_t i = 0; i < NumSec; i++)
    {
        Double_t nba = 0, theta = -5000., Esum = 0.;
        Int_t NumAnodesAngleFit = 0;
        Double_t nba_LD = 0., nba_LU = 0., nba_RD = 0., nba_RU = 0. ;
        for (Int_t j = 0; j < NumAnodes; j++)
        {
            if (energyperanode[i][j] > 1000 && energyperanode[i][j] < 65535 && StatusAnodes[i][j] == 1)
            {
                if (i==0){
                  Esum_LD += energyperanode[i][j];
                  nba_LD++;
                  if (dt[i][j] > 0.)
                    dt_LD += dt[i][j];
                }
                if (i==1){
                  Esum_LU += energyperanode[i][j];
                  nba_LU++;
                }
                if (i==2){
                  Esum_RU += energyperanode[i][j];
                  nba_RU++;
                }
                if (i==3){
                  Esum_RD += energyperanode[i][j];
                  nba_RD++;
                }
            }
        }

        if (i==0 && nba_LD > 8 && (Esum_LD / nba_LD) > 0.){
          Esum_mean_LD = Esum_LD / nba_LD;
          dt_mean_LD = dt_LD / nba_LD;
          VD_mean = 110/dt_mean_LD;
          for (Int_t i=0; i<nba_LD ; i++){
            Double_t x = VD_mean*dt[0][i];
          }
        }
        if (i==1 && nba_LU > 8 && (Esum_LU / nba_LU) > 0.){
          Esum_mean_LU = Esum_LU / nba_LU;
        }
        if (i==2 && nba_RU > 8 && (Esum_RU / nba_RU) > 0.){
          Esum_mean_RU = Esum_RU / nba_RU;
        }
        if (i==3 && nba_RD > 8 && (Esum_RD / nba_RD) > 0.){
          Esum_mean_RD = Esum_RD / nba_RD;
        }
    }

    if (CalDat)
        delete[] CalDat;
    if (HitTofW)
        delete[] HitTofW;

    if(Esum_mean_LD>0. && Esum_mean_RU>0.){
      Double_t E_LD = Esum_mean_LD;
      Double_t E_RU = Esum_mean_RU;
      //Esum_LD_RU->Fill(E_LD+E_RU);
      E_LD_vs_tof0->Fill(tof[0],E_LD);
      E_RU_vs_tof1->Fill(tof[1],E_RU);
      E_LD_corr = E_LD - (p0_LD - p0_LD_origen + p1_LD * tof[0]);// + p2_LD * tof[0] * tof[0]);
      E_RU_corr = E_RU - (p0_RU - p0_RU_origen + p1_RU * tof[1]);// + p2_RU * tof[1] * tof[1]);
      E_LD_vs_tof0_corr->Fill(tof[0],E_LD_corr);
      E_RU_vs_tof1_corr->Fill(tof[1],E_RU_corr);
      // Q = TMath::Sqrt(E)
      Z_LD_pre->Fill(TMath::Sqrt(E_LD_corr));
      Z_RU_pre->Fill(TMath::Sqrt(E_RU_corr));

      Z_LD_raw = (p0_LD_calZ + p1_LD_calZ*TMath::Sqrt(E_LD_corr));
      Z_RU_raw = (p0_RU_calZ + p1_RU_calZ*TMath::Sqrt(E_RU_corr));

      Z_LD_raw_graph->Fill(Z_LD_raw);
      Z_RU_raw_graph->Fill(Z_RU_raw);
      Z_LD_vs_RU_raw->Fill(Z_RU_raw,Z_LD_raw);
      Zsum_LD_RU_raw->Fill(Z_LD_raw+Z_RU_raw);
      Z_ff_raw->Fill(Z_RU_raw,Z_LD_raw);
      Z_f_raw->Fill(Z_LD_raw+Z_RU_raw);

      Z_LD = Z_LD_raw-5.05;//-5.0685;
      Z_RU = Z_RU_raw-5.05;//-5.0685;

      Z_LD_graph->Fill(Z_LD);
      Z_RU_graph->Fill(Z_RU);
      Z_LD_vs_RU->Fill(Z_RU,Z_LD);
      Zsum_LD_RU->Fill(Z_LD+Z_RU);

      Z_ff->Fill(Z_RU,Z_LD);
      Z_f->Fill(Z_LD+Z_RU);

    }

      if (Esum_mean_LU>0. && Esum_mean_RD>0.){
      Double_t E_LU = Esum_mean_LU;
      Double_t E_RD = Esum_mean_RD;
      //Esum_LU_RD->Fill(E_LU+E_RD);
      E_LU_vs_tof0->Fill(tof[0],E_LU);
      E_RD_vs_tof1->Fill(tof[1],E_RD);
      E_LU_corr = E_LU - (p0_LU - p0_LU_origen + p1_LU * tof[0]);// + p2_LU * tof[0] * tof[0]);
      E_RD_corr = E_RD - (p0_RD - p0_RD_origen + p1_RD * tof[1]);// + p2_RD * tof[1] * tof[1]);
      E_LU_vs_tof0_corr->Fill(tof[0],E_LU_corr);
      E_RD_vs_tof1_corr->Fill(tof[1],E_RD_corr);
      // Q = TMath::Sqrt(E)
      Z_LU_pre->Fill(TMath::Sqrt(E_LU_corr));
      Z_RD_pre->Fill(TMath::Sqrt(E_RD_corr));

      Z_LU_raw = (p0_LU_calZ + p1_LU_calZ*TMath::Sqrt(E_LU_corr));
      Z_RD_raw = (p0_RD_calZ + p1_RD_calZ*TMath::Sqrt(E_RD_corr));

      Z_LU_raw_graph->Fill(Z_LU_raw);
      Z_RD_raw_graph->Fill(Z_RD_raw);
      Z_LU_vs_RD_raw->Fill(Z_RD_raw,Z_LU_raw);
      Zsum_LU_RD_raw->Fill(Z_LU_raw+Z_RD_raw);
      Z_ff_raw->Fill(Z_RD_raw,Z_LU_raw);
      Z_f_raw->Fill(Z_LU_raw+Z_RD_raw);

      Z_LU = Z_LU_raw-4.05;//-4.5-0.4418;
      Z_RD = Z_RD_raw-5.04;//-4.5-0.4418;

      Z_LU_graph->Fill(Z_LU);
      Z_RD_graph->Fill(Z_RD);
      Z_LU_vs_RD->Fill(Z_RD,Z_LU);
      Zsum_LU_RD->Fill(Z_LU+Z_RD);

      Z_ff->Fill(Z_RD,Z_LU);
      Z_f->Fill(Z_LU+Z_RD);

    }
}


Double_t Sigma_LD = 3.0;
Double_t Thr_LD = 0.0001;
Double_t NumPeaks_LD = 22;
int nfound_LD = 0;
TSpectrum* spectrum_LD = new TSpectrum(NumPeaks_LD);
nfound_LD = spectrum_LD->Search(Z_LD_pre, Sigma_LD, "", Thr_LD);
Double_t* ChannelPeaks_LD = (Double_t*)spectrum_LD->GetPositionX();
Int_t idx_LD[nfound_LD];
TMath::Sort(nfound_LD, ChannelPeaks_LD, idx_LD, kTRUE);
Double_t Zfound_LD[nfound_LD + 1];
//cout << "nfound_LD = " << nfound_LD << endl;
for (Int_t j = 0; j < nfound_LD; j++)
{
    Zfound_LD[j] = ChannelPeaks_LD[idx_LD[nfound_LD - j - 1]]; //n canal del pico
    Z_LD_calibrating->Fill(Zfound_LD[j],j+40);
}
Zfound_LD[nfound_LD] = 0.;

Double_t Sigma_LU = 3.0;
Double_t Thr_LU = 0.0001;
Double_t NumPeaks_LU = 22;
int nfound_LU = 0;
TSpectrum* spectrum_LU = new TSpectrum(NumPeaks_LU);
nfound_LU = spectrum_LU->Search(Z_LU_pre, Sigma_LU, "", Thr_LU);
Double_t* ChannelPeaks_LU = (Double_t*)spectrum_LU->GetPositionX();
Int_t idx_LU[nfound_LU];
TMath::Sort(nfound_LU, ChannelPeaks_LU, idx_LU, kTRUE);
Double_t Zfound_LU[nfound_LU + 1];
//cout << "nfound_LU = " << nfound_LU << endl;
for (Int_t j = 0; j < nfound_LU; j++)
{
    Zfound_LU[j] = ChannelPeaks_LU[idx_LU[nfound_LU - j - 1]]; //n canal del pico
    Z_LU_calibrating->Fill(Zfound_LU[j],j+40);
}
Zfound_LU[nfound_LU] = 0.;

Double_t Sigma_RD = 3.0;
Double_t Thr_RD = 0.0001;
Double_t NumPeaks_RD = 22;
int nfound_RD = 0;
TSpectrum* spectrum_RD = new TSpectrum(NumPeaks_RD);
nfound_RD = spectrum_RD->Search(Z_RD_pre, Sigma_RD, "", Thr_RD);
Double_t* ChannelPeaks_RD = (Double_t*)spectrum_RD->GetPositionX();
Int_t idx_RD[nfound_RD];
TMath::Sort(nfound_RD, ChannelPeaks_RD, idx_RD, kTRUE); //los ordena de mayor a menor
Double_t Zfound_RD[nfound_RD + 1];
//cout << "nfound_RD = " << nfound_RD << endl;
for (Int_t j = 0; j < nfound_RD; j++)
{
    Zfound_RD[j] = ChannelPeaks_RD[idx_RD[nfound_RD - j - 1]]; // esto es para que queden ordenados de menor a mayor
    Z_RD_calibrating->Fill(Zfound_RD[j],j+40);
}
Zfound_RD[nfound_RD] = 0.;

Double_t Sigma_RU = 3.0;
Double_t Thr_RU = 0.0001;
Double_t NumPeaks_RU = 22;
int nfound_RU = 0;
TSpectrum* spectrum_RU = new TSpectrum(NumPeaks_RU);
nfound_RU = spectrum_RU->Search(Z_RU_pre, Sigma_RU, "", Thr_RU);
Double_t* ChannelPeaks_RU = (Double_t*)spectrum_RU->GetPositionX();
Int_t idx_RU[nfound_RU];
TMath::Sort(nfound_RU, ChannelPeaks_RU, idx_RU, kTRUE); //los ordena de mayor a menor
Double_t Zfound_RU[nfound_RU + 1];
//cout << "nfound_RU = " << nfound_RU << endl;
for (Int_t j = 0; j < nfound_RU; j++)
{
    Zfound_RU[j] = ChannelPeaks_RU[idx_RU[nfound_RU - j - 1]]; // esto es para que queden ordenados de menor a mayor
    Z_RU_calibrating->Fill(Zfound_RU[j],j+40);
    //cout << "j = " << j << ", nfound_RU - j - 1 = " << nfound_RU - j - 1 << ", Zfound_RU[j] = " << Zfound_RU[j] << endl;
}
Zfound_RU[nfound_RU] = 0.;

TString path = "./figures_TOF_CORR/";
//TString path = "./figures_AUX/";

// E vs TOF
TCanvas* E_LD_vs_tof0_canvas = new TCanvas("E_LD_vs_tof0_canvas", "E_LD_vs_tof0_canvas", 0, 0, 800, 400);
E_LD_vs_tof0_canvas->cd();
E_LD_vs_tof0->Draw("ZCOL");
E_LD_vs_tof0_canvas->SaveAs(path + "E_LD_vs_tof0.C");
TCanvas* E_RU_vs_tof1_canvas = new TCanvas("E_RU_vs_tof1_canvas", "E_RU_vs_tof1_canvas", 0, 0, 800, 400);
E_RU_vs_tof1_canvas->cd();
E_RU_vs_tof1->Draw("ZCOL");
E_RU_vs_tof1_canvas->SaveAs(path + "E_RU_vs_tof1.C");
TCanvas* E_LU_vs_tof0_canvas = new TCanvas("E_LU_vs_tof0_canvas", "E_LU_vs_tof0_canvas", 0, 0, 800, 400);
E_LU_vs_tof0_canvas->cd();
E_LU_vs_tof0->Draw("ZCOL");
E_LU_vs_tof0_canvas->SaveAs(path + "E_LU_vs_tof0.C");
TCanvas* E_RD_vs_tof1_canvas = new TCanvas("E_RD_vs_tof1_canvas", "E_RD_vs_tof1_canvas", 0, 0, 800, 400);
E_RD_vs_tof1_canvas->cd();
E_RD_vs_tof1->Draw("ZCOL");
E_RD_vs_tof1_canvas->SaveAs(path + "E_RD_vs_tof1.C");

// E vs TOF CORRECTED
TCanvas* E_LD_vs_tof0_corr_canvas = new TCanvas("E_LD_vs_tof0_corr_canvas", "E_LD_vs_tof0_corr_canvas", 0, 0, 800, 400);
E_LD_vs_tof0_corr_canvas->cd();
E_LD_vs_tof0_corr->Draw("ZCOL");
E_LD_vs_tof0_corr_canvas->SaveAs(path + "E_LD_vs_tof0_corr.C");
TCanvas* E_RU_vs_tof1_corr_canvas = new TCanvas("E_RU_vs_tof1_corr_canvas", "E_RU_vs_tof1_corr_canvas", 0, 0, 800, 400);
E_RU_vs_tof1_corr_canvas->cd();
E_RU_vs_tof1_corr->Draw("ZCOL");
E_RU_vs_tof1_corr_canvas->SaveAs(path + "E_RU_vs_tof1_corr.C");
TCanvas* E_LU_vs_tof0_corr_canvas = new TCanvas("E_LU_vs_tof0_corr_canvas", "E_LU_vs_tof0_corr_canvas", 0, 0, 800, 400);
E_LU_vs_tof0_corr_canvas->cd();
E_LU_vs_tof0_corr->Draw("ZCOL");
E_LU_vs_tof0_corr_canvas->SaveAs(path + "E_LU_vs_tof0_corr.C");
TCanvas* E_RD_vs_tof1_corr_canvas = new TCanvas("E_RD_vs_tof1_corr_canvas", "E_RD_vs_tof1_corr_canvas", 0, 0, 800, 400);
E_RD_vs_tof1_corr_canvas->cd();
E_RD_vs_tof1_corr->Draw("ZCOL");
E_RD_vs_tof1_corr_canvas->SaveAs(path + "E_RD_vs_tof1_corr.C");

// E CORRECTED
TCanvas* Z_LD_pre_canvas = new TCanvas("Z_LD_pre_canvas", "Z_LD_pre_canvas", 0, 0, 800, 400);
Z_LD_pre_canvas->cd();
Z_LD_pre->Draw("ZCOL");
Z_LD_pre_canvas->SaveAs(path + "Z_LD_pre.C");
TCanvas* Z_RU_pre_canvas = new TCanvas("Z_RU_pre_canvas", "Z_RU_pre_canvas", 0, 0, 800, 400);
Z_RU_pre_canvas->cd();
Z_RU_pre->Draw("ZCOL");
Z_RU_pre_canvas->SaveAs(path + "Z_RU_pre.C");
TCanvas* Z_LU_pre_canvas = new TCanvas("Z_LU_pre_canvas", "Z_LU_pre_canvas", 0, 0, 800, 400);
Z_LU_pre_canvas->cd();
Z_LU_pre->Draw("ZCOL");
Z_LU_pre_canvas->SaveAs(path + "Z_LU_pre.C");
TCanvas* Z_RD_pre_canvas = new TCanvas("Z_RD_pre_canvas", "Z_RD_pre_canvas", 0, 0, 800, 400);
Z_RD_pre_canvas->cd();
Z_RD_pre->Draw("ZCOL");
Z_RD_pre_canvas->SaveAs(path + "Z_RD_pre.C");

// CALIBRATION Z CHANNELS TO Z
TCanvas* Z_LD_calibrating_canvas = new TCanvas("Z_LD_calibrating_canvas", "Z_LD_calibrating_canvas", 0, 0, 800, 400);
Z_LD_calibrating_canvas->cd();
Z_LD_calibrating->Draw("");
TF1* pol1_LD = new TF1("pol1_LD","pol1", 95, 160); //152
Z_LD_calibrating->Fit(pol1_LD,"R");
Z_LD_calibrating_canvas->SaveAs(path + "Z_LD_calibrating.C");

TCanvas* Z_RU_calibrating_canvas = new TCanvas("Z_RU_calibrating_canvas", "Z_RU_calibrating_canvas", 0, 0, 800, 400);
Z_RU_calibrating_canvas->cd();
Z_RU_calibrating->Draw("");
TF1* pol1_RU = new TF1("pol1_RU","pol1", 95, 160); //155
Z_RU_calibrating->Fit(pol1_RU,"R");
Z_RU_calibrating_canvas->SaveAs(path + "Z_RU_calibrating.C");

TCanvas* Z_LU_calibrating_canvas = new TCanvas("Z_LU_calibrating_canvas", "Z_LU_calibrating_canvas", 0, 0, 800, 400);
Z_LU_calibrating_canvas->cd();
Z_LU_calibrating->Draw("");
TF1* pol1_LU = new TF1("pol1_LU","pol1", 95, 160); //154
Z_LU_calibrating->Fit(pol1_LU,"R");
Z_LU_calibrating_canvas->SaveAs(path + "Z_LU_calibrating.C");

TCanvas* Z_RD_calibrating_canvas = new TCanvas("Z_RD_calibrating_canvas", "Z_RD_calibrating_canvas", 0, 0, 800, 400);
Z_RD_calibrating_canvas->cd();
Z_RD_calibrating->Draw("");
TF1* pol1_RD = new TF1("pol1_RD","pol1", 95, 160); //156
Z_RD_calibrating->Fit(pol1_RD,"R");
Z_RD_calibrating_canvas->SaveAs(path + "Z_RD_calibrating.C");

// RAW
TCanvas* Z_LD_raw_canvas = new TCanvas("Z_LD_raw_canvas", "Z_LD_raw_canvas", 0, 0, 800, 400);
Z_LD_raw_canvas->cd();
Z_LD_raw_graph->Draw("ZCOL");
Z_LD_raw_canvas->SaveAs(path + "Z_LD_raw_graph.C");

TCanvas* Z_RU_raw_canvas = new TCanvas("Z_RU_raw_canvas", "Z_RU_raw_canvas", 0, 0, 800, 400);
Z_RU_raw_canvas->cd();
Z_RU_raw_graph->Draw("ZCOL");
Z_RU_raw_canvas->SaveAs(path + "Z_RU_raw_graph.C");

TCanvas* Z_LU_raw_canvas = new TCanvas("Z_LU_raw_canvas", "Z_LU_raw_canvas", 0, 0, 800, 400);
Z_LU_raw_canvas->cd();
Z_LU_raw_graph->Draw("ZCOL");
Z_LU_raw_canvas->SaveAs(path + "Z_LU_raw_graph.C");

TCanvas* Z_RD_raw_canvas = new TCanvas("Z_RD_raw_canvas", "Z_RD_raw_canvas", 0, 0, 800, 400);
Z_RD_raw_canvas->cd();
Z_RD_raw_graph->Draw("ZCOL");
Z_RD_raw_canvas->SaveAs(path + "Z_RD_raw_graph.C");

// Z RAW PART 2

TCanvas* Z_LD_vs_RU_raw_canvas = new TCanvas("Z_LD_vs_RU_raw_canvas", "Z_LD_vs_RU_raw_canvas", 0, 0, 800, 400);
Z_LD_vs_RU_raw_canvas->cd();
Z_LD_vs_RU_raw->Draw("ZCOL");
Z_LD_vs_RU_raw_canvas->SaveAs(path + "Z_LD_vs_RU_raw_canvas.C");

TCanvas* Zsum_LD_RU_raw_canvas = new TCanvas("Zsum_LD_RU_raw_canvas", "Zsum_LD_RU_raw_canvas", 0, 0, 800, 400);
Zsum_LD_RU_raw_canvas->cd();
Zsum_LD_RU_raw->Draw("ZCOL");
Zsum_LD_RU_raw_canvas->SaveAs(path + "Zsum_LD_RU_raw_canvas.C");

TCanvas* Z_LU_vs_RD_raw_canvas = new TCanvas("Z_LU_vs_RD_raw_canvas", "Z_LU_vs_RD_raw_canvas", 0, 0, 800, 400);
Z_LU_vs_RD_raw_canvas->cd();
Z_LU_vs_RD_raw->Draw("ZCOL");
Z_LU_vs_RD_raw_canvas->SaveAs(path + "Z_LU_vs_RD_raw_canvas.C");

TCanvas* Zsum_LU_RD_raw_canvas = new TCanvas("Zsum_LU_RD_raw_canvas", "Zsum_LU_RD_raw_canvas", 0, 0, 800, 400);
Zsum_LU_RD_raw_canvas->cd();
Zsum_LU_RD_raw->Draw("ZCOL");
Zsum_LU_RD_raw_canvas->SaveAs(path + "Zsum_LU_RD_raw_canvas.C");

// Z RAW FISSION FRAGMENTS
TCanvas* Z_ff_raw_canvas = new TCanvas("Z_ff_raw_canvas", "Z_ff_raw_canvas", 0, 0, 800, 400);
Z_ff_raw_canvas->cd();
Z_ff_raw->Draw("ZCOL");
Z_ff_raw_canvas->SaveAs(path + "Z_ff_raw.C");

// Z RAW FISSIONING SYSTEM
TCanvas* Z_f_raw_canvas = new TCanvas("Z_f_raw_canvas", "Z_f_raw_canvas", 0, 0, 800, 400);
Z_f_raw_canvas->cd();
Z_f_raw->Draw("ZCOL");
Z_f_raw_canvas->SaveAs(path + "Z_f_raw.C");

// COMPARISON
TCanvas* Z_LD_RU_comparison_canvas = new TCanvas("Z_LD_RU_comparison_canvas", "Z_LD_RU_comparison_canvas", 0, 0, 800, 400);
Z_LD_RU_comparison_canvas->cd();
Z_LD_raw_graph->Draw("ZCOL");
Z_RU_raw_graph->Draw("SAME");
Z_LD_RU_comparison_canvas->SaveAs(path + "Z_LD_RU_comparison.C");

TCanvas* Z_LU_RD_comparison_canvas = new TCanvas("Z_LU_RD_comparison_canvas", "Z_LU_RD_comparison_canvas", 0, 0, 800, 400);
Z_LU_RD_comparison_canvas->cd();
Z_LU_raw_graph->Draw("ZCOL");
Z_RD_raw_graph->Draw("SAME");
Z_LU_RD_comparison_canvas->SaveAs(path + "Z_LU_RD_comparison.C");

TCanvas* Z_LD_LU_comparison_canvas = new TCanvas("Z_LD_LU_comparison_canvas", "Z_LD_LU_comparison_canvas", 0, 0, 800, 400);
Z_LD_LU_comparison_canvas->cd();
Z_LD_raw_graph->Draw("ZCOL");
Z_LU_raw_graph->Draw("SAME");
Z_LD_LU_comparison_canvas->SaveAs(path + "Z_LD_LU_comparison.C");

TCanvas* Z_RD_RU_comparison_canvas = new TCanvas("Z_RD_RU_comparison_canvas", "Z_RD_RU_comparison_canvas", 0, 0, 800, 400);
Z_RD_RU_comparison_canvas->cd();
Z_RD_raw_graph->Draw("ZCOL");
Z_RU_raw_graph->Draw("SAME");
Z_RD_RU_comparison_canvas->SaveAs(path + "Z_RD_RU_comparison.C");

// Z GRAPHS PART 2
TCanvas* Z_LD_vs_RU_canvas = new TCanvas("Z_LD_vs_RU_canvas", "Z_LD_vs_RU_canvas", 0, 0, 800, 400);
Z_LD_vs_RU_canvas->cd();
Z_LD_vs_RU->Draw("ZCOL");
Z_LD_vs_RU_canvas->SaveAs(path + "Z_LD_vs_RU.C");

TCanvas* Zsum_LD_RU_canvas = new TCanvas("Zsum_LD_RU_canvas", "Zsum_LD_RU_canvas", 0, 0, 800, 400);
Zsum_LD_RU_canvas->cd();
Zsum_LD_RU->Draw("ZCOL");
Zsum_LD_RU_canvas->SaveAs(path + "Zsum_LD_RU.C");

TCanvas* Z_LU_vs_RD_canvas = new TCanvas("Z_LU_vs_RD_canvas", "Z_LU_vs_RD_canvas", 0, 0, 800, 400);
Z_LU_vs_RD_canvas->cd();
Z_LU_vs_RD->Draw("ZCOL");
Z_LU_vs_RD_canvas->SaveAs(path + "Z_LU_vs_RD.C");

TCanvas* Zsum_LU_RD_canvas = new TCanvas("Zsum_LU_RD_canvas", "Zsum_LU_RD_canvas", 0, 0, 800, 400);
Zsum_LU_RD_canvas->cd();
Zsum_LU_RD->Draw("ZCOL");
Zsum_LU_RD_canvas->SaveAs(path + "Zsum_LU_RD.C");

// Z FISSION FRAGMENTS
TCanvas* Z_ff_canvas = new TCanvas("Z_ff_canvas", "Z_ff_canvas", 0, 0, 800, 400);
Z_ff_canvas->cd();
Z_ff->Draw("ZCOL");
Z_ff_canvas->SaveAs(path + "Z_ff.C");

// Z FISSIONING SYSTEM
TCanvas* Z_f_canvas = new TCanvas("Z_f_canvas", "Z_f_canvas", 0, 0, 800, 400);
Z_f_canvas->cd();
Z_f->Draw("ZCOL");
Z_f_canvas->SaveAs(path + "Z_f.C");

}
