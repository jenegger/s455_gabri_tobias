/*
  TPat Matrix:

  // ------------  On spill triggers

  // ----------- When AMS is alive.
  TRIG_LMU_OUT( 1) = BEAM_GATE_AUX and not AMS_DEAD_AUX and in_SOF_START;
	TRIG_LMU_OUT( 2) = BEAM_GATE_AUX and not AMS_DEAD_AUX and in_SOF_START and in_SOF_FISSION;
	TRIG_LMU_OUT( 3) = BEAM_GATE_AUX and not AMS_DEAD_AUX and in_SOF_START and in_CALIFA_AND;
	TRIG_LMU_OUT( 4) = BEAM_GATE_AUX and not AMS_DEAD_AUX and in_SOF_START and in_SOF_FISSION and in_CALIFA_AND;
	TRIG_LMU_OUT( 5) = BEAM_GATE_AUX and not AMS_DEAD_AUX and in_SOF_START and in_CALIFA_OR;
	TRIG_LMU_OUT( 6) = BEAM_GATE_AUX and not AMS_DEAD_AUX and in_SOF_START and in_SOF_FISSION and in_CALIFA_OR;

	// ------------- When AMS is busy but not close to finishing readout.
	TRIG_LMU_OUT( 7) = BEAM_GATE_AUX and NONAMS_BONUS_AUX and in_SOF_START;
	TRIG_LMU_OUT( 8) = BEAM_GATE_AUX and NONAMS_BONUS_AUX and in_SOF_START and in_SOF_FISSION;
	TRIG_LMU_OUT( 9) = BEAM_GATE_AUX and NONAMS_BONUS_AUX and in_SOF_START and in_CALIFA_AND;
	TRIG_LMU_OUT(10) = BEAM_GATE_AUX and NONAMS_BONUS_AUX and in_SOF_START and in_SOF_FISSION and in_CALIFA_AND;
	TRIG_LMU_OUT(11) = BEAM_GATE_AUX and NONAMS_BONUS_AUX and in_SOF_START and in_CALIFA_OR;
	TRIG_LMU_OUT(12) = BEAM_GATE_AUX and NONAMS_BONUS_AUX and in_SOF_START and in_SOF_FISSION and in_CALIFA_OR;

	// ------------ Off spill triggers.
	TRIG_LMU_OUT(13) = not BEAM_GATE_AUX and in_CALIFA_OR;
	TRIG_LMU_OUT(14) = not BEAM_GATE_AUX and in_NEULAND;
*/
void p2p_study()
{

  TStopwatch timer;
  timer.Start();

  // Put here the root file
  TString fileList = "/home/gabri/Analysis/Califa_tests/files/rootfiles/calibrated_main0273_HT.root";

  TFile *eventFile;
  TTree* eventTree;


  eventFile = TFile::Open(fileList);
  eventTree = (TTree*)eventFile->Get("evt");
  eventTree->SetBranchStatus("*",0);

  eventTree->SetBranchStatus("CalifaHitData*",1);
  eventTree->SetBranchStatus("CalifaCrystalCalData*",1);
  eventTree->SetBranchStatus("TwimHitData*",1);
  eventTree->SetBranchStatus("EventHeader.*",1);

  TClonesArray *calCA = new TClonesArray("R3BCalifaCrystalCalData",5);
  TBranch  *calBranch = eventTree->GetBranch("CalifaCrystalCalData");
  calBranch->SetAddress(&calCA);


  TClonesArray *hitCA = new TClonesArray("R3BCalifaHitData",5);
  TBranch  *hitBranch = eventTree->GetBranch("CalifaHitData");
  hitBranch->SetAddress(&hitCA);

  TClonesArray *twimCA = new TClonesArray("R3BSofTwimHitData",5);
  TBranch  *twimBranch = eventTree->GetBranch("TwimHitData");
  twimBranch->SetAddress(&twimCA);


  TBranch  *tpatBranch = eventTree->GetBranch("EventHeader.");
  R3BEventHeader *fEventHeader;
  tpatBranch->SetAddress(&fEventHeader);



  // ------------ Histograms for Correlations --------------
  TCanvas *myFirstCanvas = new TCanvas("myFirstCanvas","myFirstCanvas",2000,2000);
  myFirstCanvas->Divide(3,1);

  TCanvas *mySecondCanvas = new TCanvas("mySecondCanvas","mySecondCanvas",2000,2000);

  TH2F *phiVsPhiHisto = new TH2F("phiVsPhiHisto","P2P : Phi Vs Phi ",300,-190,190,300,-190,190);
  TH2F *energyCorrHisto = new TH2F("energyCorrHisto","P2P : Energy Vs Energy ",400,0,700,400,0,700);
  TH1F *openingAngleHisto = new TH1F("openingAngleHisto","P2P : Opening Angle",300,0,100);

  TH1F *multHisto = new TH1F("multHisto","Crystal Multiplicity",300,0,100);

  Int_t nEvents = eventTree->GetEntries();




  Int_t nCalifaHits,nCalifaCalHits,nTwimHits,fValidHits,fCrystalId,tpatbin;
  Int_t fTPat,fUncorrelatedEvents=0;
  Int_t goodEvent;

  Float_t fEnergy,fCharge,fTheta2,fTheta1,fPhi1,fPhi2,fEnergy1,fEnergy2;


  for (Int_t j = 0; j<nEvents; j++) {

    hitCA->Clear();
    twimCA->Clear();
    calCA->Clear();

    eventTree->GetEvent(j);

    nCalifaHits = hitCA->GetEntries();
    nCalifaCalHits=calCA->GetEntries();
    nTwimHits = twimCA->GetEntries();

    fEnergy=0.0;
    fCharge=0.0;
    goodEvent=0;

    if(!(j%1000000))
     cout<<"Reading event "<<j<<" out of "<<nEvents<<" ("<<100.0*Float_t(j)/Float_t(nEvents)<<" % ) "<<endl;



     if (fEventHeader->GetTpat() > 0)
     {
         for (Int_t i = 0; i < 16; i++)
         {
             tpatbin = (fEventHeader->GetTpat() & (1 << i));
             if (tpatbin != 0){
                 fTPat=i+1;
                }
             }
          }


     else {
      fTPat=0;
      fUncorrelatedEvents++;
     }


     for(Int_t s=0;s<nTwimHits;s++)
       fCharge = fCharge + ((R3BSofTwimHitData*)twimCA->At(s))->GetZcharge();


     for(Int_t z=0;z<nCalifaCalHits;z++){

         if (0.001*((R3BCalifaCrystalCalData*)calCA->At(z))->GetEnergy() >= 100)
         goodEvent++;


     }




     if(nCalifaHits >= 2 && fCharge > 85 && fCharge < 98 && nTwimHits ==2 &&  goodEvent>=2 /*&& (fTPat == 4 || fTPat == 10)*/){

         fPhi1 = TMath::RadToDeg()*((R3BCalifaHitData*)hitCA->At(0))->GetPhi();
         fPhi2 = TMath::RadToDeg()*((R3BCalifaHitData*)hitCA->At(1))->GetPhi();

         fTheta1 = TMath::RadToDeg()*((R3BCalifaHitData*)hitCA->At(0))->GetTheta();
         fTheta2 = TMath::RadToDeg()*((R3BCalifaHitData*)hitCA->At(1))->GetTheta();

         fEnergy1 = 0.001*((R3BCalifaHitData*)hitCA->At(0))->GetEnergy();
         fEnergy2 = 0.001*((R3BCalifaHitData*)hitCA->At(1))->GetEnergy();

         multHisto->Fill(nCalifaCalHits);

         if(( (TMath::Abs(fPhi1) + TMath::Abs(fPhi2)) > 160 && (TMath::Abs(fPhi1) + TMath::Abs(fPhi2) < 200))){

           openingAngleHisto->Fill(fTheta1+fTheta2);
           energyCorrHisto->Fill(fEnergy1,fEnergy2);
           phiVsPhiHisto->Fill(fPhi1,fPhi2);

        }
      }
    }


   myFirstCanvas->cd(1);
   openingAngleHisto->Draw();

   myFirstCanvas->cd(2);
   phiVsPhiHisto->Draw("colz");

   myFirstCanvas->cd(3);
   energyCorrHisto->Draw("colz");

   mySecondCanvas->cd();
   multHisto->Draw();

  cout<<"Uncorrelated Events: "<<fUncorrelatedEvents<<"( "<<100*Float_t(fUncorrelatedEvents)/Float_t(nEvents)<<" )"<<endl;



}
