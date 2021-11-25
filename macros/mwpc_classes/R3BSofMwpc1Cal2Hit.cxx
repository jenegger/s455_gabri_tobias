// ----------------------------------------------------------------------
// -----		         R3BSofMwpc1Cal2Hit                           -----
// -----  Created 13/10/21  by Antía Graña González               -----
// -----  by modifying J.L. Rodriguez-Sanchez classes for Mwpc    -----
// ----------------------------------------------------------------------

// ROOT headers
#include "TClonesArray.h"
#include "TMath.h"

// Fair headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include <iomanip>

// MWPC headers
#include "R3BSofMwpc1Cal2Hit.h"
#include "R3BSofMwpcCalData.h"
#include "R3BSofMwpcHitData.h"

// R3BSofMwpc1Cal2Hit: Default Constructor --------------------------
R3BSofMwpc1Cal2Hit::R3BSofMwpc1Cal2Hit()
    : FairTask("R3B Hit-MWPC1 Task", 1)
    , fMwpcCalDataCA(NULL)
    , fMwpcHitDataCA(NULL)
    , fwx(3.125)   // in mm
    , fwy(5.000)   // in mm
    , fSize(200.0) // in mm
    , fOnline(kFALSE)
{
}

// R3BSofMwpc1Cal2Hit: Standard Constructor --------------------------
R3BSofMwpc1Cal2Hit::R3BSofMwpc1Cal2Hit(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMwpcCalDataCA(NULL)
    , fMwpcHitDataCA(NULL)
    , fwx(3.125)   // in mm
    , fwy(5.000)   // in mm
    , fSize(200.0) // in mm
    , fOnline(kFALSE)
{
}

// Virtual R3BSofMwpc1Cal2Hit: Destructor
R3BSofMwpc1Cal2Hit::~R3BSofMwpc1Cal2Hit()
{
    LOG(INFO) << "R3BSofMwpc1Cal2Hit: Delete instance";
    if (fMwpcCalDataCA)
        delete fMwpcCalDataCA;
    if (fMwpcHitDataCA)
        delete fMwpcHitDataCA;
}

// -----   Public method Init   --------------------------------------------
InitStatus R3BSofMwpc1Cal2Hit::Init()
{
    LOG(INFO) << "R3BSofMwpc1Cal2Hit::Init()";

    // INPUT DATA
    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        return kFATAL;
    }

    fMwpcCalDataCA = (TClonesArray*)rootManager->GetObject("Mwpc1CalData");
    if (!fMwpcCalDataCA)
    {
        return kFATAL;
    }

    // OUTPUT DATA
    // Hit data
    fMwpcHitDataCA = new TClonesArray("R3BSofMwpcHitData", 10);
    rootManager->Register("Mwpc1HitData", "MWPC1 Hit", fMwpcHitDataCA, !fOnline);

    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BSofMwpc1Cal2Hit::ReInit() { return kSUCCESS; }

bool sortPairsmwpc1(const pair<Double_t, Int_t> &x, const pair<Double_t, Int_t> &y)
{
  return x.first > y.first;
}

Double_t R3BSofMwpc1Cal2Hit::GetPositionX(Double_t qmax, Int_t padmax, Double_t qleft, Double_t qright)
{
    Double_t a3 = TMath::Pi() * fwx / (TMath::ACosH(0.5 * (TMath::Sqrt(qmax / qleft) + TMath::Sqrt(qmax / qright))));
    // Double_t a2 = gRandom->Uniform(-fwx / 2,fwx / 2);
    Double_t a2 = (a3 / TMath::Pi()) * TMath::ATanH((TMath::Sqrt(qmax / qleft) - TMath::Sqrt(qmax / qright)) /
                                                    (2 * TMath::SinH(TMath::Pi() * fwx / a3)));

    return (-1. * padmax * fwx + (fSize / 2) - (fwx / 2) - a2); // Left is positive and right negative
}

/* ----   Protected method to obtain the position Y ---- */
Double_t R3BSofMwpc1Cal2Hit::GetPositionY(Double_t qmax, Int_t padmax, Double_t qdown, Double_t qup)
{
    Double_t a3 = TMath::Pi() * fwy / (TMath::ACosH(0.5 * (TMath::Sqrt(qmax / qdown) + TMath::Sqrt(qmax / qup))));
    // Double_t a2 = gRandom->Uniform(-fwy / 2, fwy / 2);
    Double_t a2 = (a3 / TMath::Pi()) * TMath::ATanH((TMath::Sqrt(qmax / qdown) - TMath::Sqrt(qmax / qup)) /
                                                    (2 * TMath::SinH(TMath::Pi() * fwy / a3)));

    return (padmax * fwy - (fSize / 2) + (fwy / 2) + a2);
}


// -----   Public method Execution   --------------------------------------------
void R3BSofMwpc1Cal2Hit::Exec(Option_t* option)
{
    // Reset entries in output arrays, local arrays
    Reset();

    // Reading the Input -- Cal Data --
    Int_t nHits = fMwpcCalDataCA->GetEntries();
    R3BSofMwpcCalData** calData = new R3BSofMwpcCalData*[nHits];

    Int_t planeId = 0;
    Int_t padId = 0;
    Double_t q = 0;
    Int_t planex1 = 0, planex2 = 0;
    Int_t padmx1 = -1, padmy1 = -1, padmx2 = -1, padmy2 = -1;
    Int_t padmx_p1 = -1, padmx_p2 = -1;
    Int_t qmx_p1 = 0., qmx_p2 = 0.;
    Double_t qmx1 = 0., qmy1 = 0., qleft1 = 0., qright1 = 0., qdown1 = 0., qup1 = 0.;
    Double_t qmx2 = 0., qmy2 = 0., qleft2 = 0., qright2 = 0., qdown2 = 0., qup2 = 0.;
    Double_t x1 = -1000, y1 = -1000, x2 = -1000, y2 = -1000;
    Int_t nx_p1 = 0, nx_p2 = 0, ny = 0;

    vector<pair<Double_t, Int_t>> QpadX_p1, QpadX_p2, QpadY;
    QpadX_p1.clear();
    QpadX_p2.clear();
    QpadY.clear();
    //Double_t fx[Mw1PadsX], fy[Mw1PadsY];
    //cout << "NUEVO EVENTO" << endl;
    fx_p1[Mw1PadsX] = {0.};
    fx_p2[Mw1PadsX] = {0.};
    fy[Mw1PadsY] = {0.};
    for (Int_t i = 0; i < nHits; i++){
      calData[i] = (R3BSofMwpcCalData*)(fMwpcCalDataCA->At(i));
      planeId = calData[i]->GetPlane();
      padId = calData[i]->GetPad() - 1;
      q = calData[i]->GetQ();
      //cout << "i = " << i << ", q = " << q << ", padId = " << padId << ", planeId = " << planeId << endl;
      pair<Double_t, Int_t> hit_pair = make_pair(q,padId);
      if (planeId == 1)
      {
        fx_p1[padId] = q;
        QpadX_p1.push_back(hit_pair);
        nx_p1 = nx_p1 + 1;
      }
      if (planeId == 2)
      {
        fx_p2[padId] = q;
        QpadX_p2.push_back(hit_pair);
        nx_p2 = nx_p2 + 1;
      }
      if (planeId == 3)
      {
        fy[padId] = q;
        QpadY.push_back(hit_pair);
        ny = ny + 1 ;
      }
    }

    if ( (nx_p1>0 || nx_p2>0) && ny>0){

    if (nx_p1>0){
    sort(QpadX_p1.begin(),QpadX_p1.end(),sortPairsmwpc1); // el vector se ordena por el primer elemento de cada par
    qmx_p1 = QpadX_p1[0].first;
    padmx_p1 = QpadX_p1[0].second;
    //cout << "qmx1 = " << qmx1 << ", padmx1 = " << padmx1 << endl;
    /*for (Int_t i=0; i<nx_p1; i++){
      cout << "Loop i = " << i << " q = " << QpadX_p1[i].first << endl;
    }*/
    }
    if (nx_p2>0){
    sort(QpadX_p2.begin(),QpadX_p2.end(),sortPairsmwpc1); // el vector se ordena por el primer elemento de cada par
    qmx_p2 = QpadX_p2[0].first;
    padmx_p2 = QpadX_p2[0].second;
    //cout << "qmx1 = " << qmx1 << ", padmx1 = " << padmx1 << endl;
    }
    if (qmx_p1>qmx_p2 && padmx_p1 + 1 < Mw1PadsX && padmx_p2 + 1 < Mw1PadsX){

      qmx1 = qmx_p1;
      padmx1 = padmx_p1;
      qmx2 = qmx_p2;
      padmx2 = padmx_p2;
      qleft1 = fx_p1[padmx1 - 1];
      qright1 = fx_p1[padmx1 + 1];
      qleft2 = fx_p2[padmx2 - 1];
      qright2 = fx_p2[padmx2 + 1];
      planex1 = 1;
      planex2 = 2;
      /*cout << "1>2" << endl;
      cout << "qmx1 = " << qmx1 << ", padmx1 = " << padmx1 << endl;
      cout << "qmx2 = " << qmx2 << ", padmx2 = " << padmx2 << endl;*/
      if ( qleft1 == 0){
	qleft1 = 1;
	}
      if (qright1 == 0){
	(qright1 = 1;
	}
      if (qleft2 == 0){
	qleft2 = 1;
	}
      if (qright2 == 0){
	qright2 = 1;
	}

      if (qmx1 > 10 && qleft1 > 0 && qright1 > 0){
        x1 = GetPositionX(qmx1, padmx1, qleft1, qright1);
      }
      if (qmx2 > 10 && qleft2 > 0 && qright2 > 0){
        x2 = GetPositionX(qmx2, padmx2, qleft2, qright2);
      }
    }

    if (qmx_p2>qmx_p1 && padmx_p1 + 1 < Mw1PadsX && padmx_p2 + 1 < Mw1PadsX){

      qmx1 = qmx_p2;
      padmx1 = padmx_p2;
      qmx2 = qmx_p1;
      padmx2 = padmx_p1;
      qleft1 = fx_p2[padmx1 - 1];
      qright1 = fx_p2[padmx1 + 1];
      qleft2 = fx_p1[padmx2 - 1];
      qright2 = fx_p1[padmx2 + 1];
      planex1 = 2;
      planex2 = 1;
      /*cout << "2>1" << endl;
      cout << "qmx1 = " << qmx1 << ", padmx1 = " << padmx1 << endl;
      cout << "qmx2 = " << qmx2 << ", padmx2 = " << padmx2 << endl;*/
      if ( qleft1 == 0){
	qleft1 = 1;
	}
      if (qright1 == 0){
	(qright1 = 1;
	}
      if (qleft2 == 0){
	qleft2 = 1;
	}
      if (qright2 == 0){
	qright2 = 1;
	}
      if (qmx1 > 10 && qleft1 > 0 && qright1 > 0){
        x1 = GetPositionX(qmx1, padmx1, qleft1, qright1);
      }
      if (qmx2 > 10 && qleft2 > 0 && qright2 > 0){
        x2 = GetPositionX(qmx2, padmx2, qleft2, qright2);
      }
    }

    if (ny>0){

      sort(QpadY.begin(),QpadY.end(),sortPairsmwpc1);
      qmy1 = QpadY[0].first;
      padmy1 = QpadY[0].second;
      //cout << "qmy1 = " << qmy1 << ", padmy1 = " << padmy1 << endl;

      if (padmy1 + 1 < Mw1PadsY)
      {
        // Obtain position Y for 1st charge ----
        qdown1 = fy[padmy1 - 1];
        qup1 = fy[padmy1 + 1];
	if (qdown1 == 0){
		qdown1 = 1;
		}
	if (qup1 == 0){
		qup1 = 1;
		}
        if (qmy1 > 10 && qdown1 > 0 && qup1 > 0){
          y1 = GetPositionY(qmy1, padmy1, qdown1, qup1);
        }
      }

      for (Int_t i =0; i<ny; i++){
        if (QpadY[i].second == padmy1-3 || QpadY[i].second == padmy1-2 || QpadY[i].second == padmy1-1 || QpadY[i].second == padmy1 ||
          QpadY[i].second == padmy1+1 || QpadY[i].second == padmy1+2 || QpadY[i].second == padmy1+3 ){
          QpadY[i].first = 0.;
        }
      }

      sort(QpadY.begin(),QpadY.end(),sortPairsmwpc1);
      qmy2 = QpadY[0].first;
      padmy2 = QpadY[0].second;
      //cout << "qmy2 = " << qmy2 << ", padmy2 = " << padmy2 << endl;

      if (padmy2 + 1 < Mw1PadsY)
      {
        // Obtain position Y for 2nd charge ----
        qdown2 = fy[padmy2 - 1];
        qup2 = fy[padmy2 + 1];
	if (qdown2 == 0){
		qdown2 = 1;
		}
	if (qup2 == 0){
		qup2 = 1;
		}
        if (qmy2 > 10 && qdown2 > 0 && qup2 > 0){
          y2 = GetPositionY(qmy2, padmy2, qdown2, qup2);
        }
      }
    } //if y
    //cout << "x1 = " << x1 << ", y1 = " << y1 << " ,planex1 = " << planex1 << ", x2 = " << x2 << ", y2 = " << y2 << " ,planex2 = " << planex2 << endl;
    AddHitData(x1, y1, planex1);
    AddHitData(x2, y2, planex2);
  } // if all

    if (calData)
        delete[] calData;
    return;
}

// -----   Public method Finish  ------------------------------------------------
void R3BSofMwpc1Cal2Hit::Finish() {}

// -----   Public method Reset   ------------------------------------------------
void R3BSofMwpc1Cal2Hit::Reset()
{
    LOG(DEBUG) << "Clearing Mwpc1HitData Structure";
    if (fMwpcHitDataCA)
        fMwpcHitDataCA->Clear();
}

// -----   Private method AddHitData  --------------------------------------------
R3BSofMwpcHitData* R3BSofMwpc1Cal2Hit::AddHitData(Double_t x, Double_t y, Int_t plane)
{
    // It fills the R3BSofMwpcHitData
    TClonesArray& clref = *fMwpcHitDataCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BSofMwpcHitData(x, y, plane);
}

ClassImp(R3BSofMwpc1Cal2Hit)
