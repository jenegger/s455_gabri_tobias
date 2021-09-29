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
void p2p_maybe()
{

  TStopwatch timer;
  timer.Start();

  // Put here the root file
  TString fileList = "calibrated_main0273.root";

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




  TCanvas *myCanvas = new TCanvas("myCanvas","myCanvas",2000,2000);

  myCanvas->Divide(3,1);


  TH2F *energyMultHisto = new TH2F("energyMultHisto","P2P Like: Energy Vs Multiplicity",600,0,700,300,0,300);

  TH2F *phiVsPhiHisto = new TH2F("phiVsPhiHisto","P2P : Phi Vs Phi ",300,-190,190,300,-190,190);
  TH2F *thetaVsThetaHisto = new TH2F("thetaVsThetaHisto","P2P : Theta Vs Theta ",300,0,90,300,0,90);
  TH1F *openingAngleHisto = new TH1F("openingAngleHisto","P2P : Opening Angle",300,0,100);
  Int_t nEvents = eventTree->GetEntries();


  Int_t nCalifaHits,nCalifaCalHits,nTwimHits,fValidHits,fCrystalId,tpatbin,fTPat,fUncorrelatedEvents=0,fUncorrelatedEventsFission=0;
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

     else{

     fTPat=0;
     fUncorrelatedEvents++;
 }


     for(Int_t s=0;s<nTwimHits;s++)
       fCharge = fCharge + ((R3BSofTwimHitData*)twimCA->At(s))->GetZcharge();


     for(Int_t z=0;z<nCalifaCalHits;z++){

         if (0.001*((R3BCalifaCrystalCalData*)calCA->At(z))->GetEnergy() >= 100)
         goodEvent++;


     }




     if(nCalifaHits == 2 && fCharge > 80 && fCharge < 100 && nTwimHits ==2 && goodEvent>=2){

         fPhi1 = TMath::RadToDeg()*((R3BCalifaHitData*)hitCA->At(0))->GetPhi();
         fPhi2 = TMath::RadToDeg()*((R3BCalifaHitData*)hitCA->At(1))->GetPhi();

         fTheta1 = TMath::RadToDeg()*((R3BCalifaHitData*)hitCA->At(0))->GetTheta();
         fTheta2 = TMath::RadToDeg()*((R3BCalifaHitData*)hitCA->At(1))->GetTheta();

         fEnergy1 = 0.001*((R3BCalifaHitData*)hitCA->At(0))->GetEnergy();
         fEnergy2 = 0.001*((R3BCalifaHitData*)hitCA->At(1))->GetEnergy();

         openingAngleHisto->Fill(fTheta1+fTheta2);
         thetaVsThetaHisto->Fill(fTheta1,fTheta2);
         phiVsPhiHisto->Fill(fPhi1,fPhi2);

     }






}

  myCanvas->cd(1);
  openingAngleHisto->Draw();

  myCanvas->cd(2);
  phiVsPhiHisto->Draw("colz");

  myCanvas->cd(3);
  thetaVsThetaHisto->Draw("colz");

  cout<<"Uncorrelated Events: "<<fUncorrelatedEvents<<"( "<<100*Float_t(fUncorrelatedEvents)/Float_t(nEvents)<<" )"<<endl;
  cout<<"Uncorrelated Events (FISSION): "<<fUncorrelatedEventsFission<<"( "<<100*Float_t(fUncorrelatedEventsFission)/Float_t(nEvents)<<" )"<<endl;






}
