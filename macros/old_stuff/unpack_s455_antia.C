/*
Unipacking macro for s455 offline data analysis
This macro requires an input .lmd file and generates a root file 
 */
typedef struct EXT_STR_h101_t
{
    EXT_STR_h101_unpack_t unpack;
    EXT_STR_h101_TPAT_t unpacktpat;
    EXT_STR_h101_CALIFA_t califa;
    EXT_STR_h101_WRMASTER_t wrm;
    EXT_STR_h101_WRCALIFA_t wrcalifa;
    EXT_STR_h101_SOFTOFW_onion_t tofw;
    EXT_STR_h101_SOFTWIM_onion_t twim;
    EXT_STR_h101_SOFSCI_onion_t sci;
    EXT_STR_h101_WRSOFIA_t wrsofia;

} EXT_STR_h101;

void unpack_s455_antia()
{
    TStopwatch timer;
    timer.Start();

    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y%m%d_%H%M%S");

    const Int_t nev = -1; // number of events to read, -1 - until CTRL+C

    const Int_t expId = 455; // select experiment: 444, 455 or 467
    Int_t NumSofSci = 1;
    Int_t sofiaWR_SE = 0xe00;
    Int_t sofiaWR_ME = 0xf00;

    // Create input -----------------------------------------
//    TString filename = "../../s455_03_273_10_stitched.lmd"; //input lmd file
    TString filename = "../../s455_03_273_10_stitched.lmd"; //input lmd file subrun 1


    // Define output-----------------------------------------
//      TString outputFileName = "../file_src/s455_03_273_10_sq_w_cl_antia.root";
//    TString outputFileName = "../file_src/s455_03_273_1_sq_w_cl_antia.root";
      TString outputFileName = "../file_src/s455_03_273_10_sq_w_cl_antia_newest_parameters.root";

    Bool_t Cal_level = true;          // set true if there exists a file with the calibration parameters
    Bool_t NOTstoremappeddata = false; // if true, don't store mapped data in the root file
    Bool_t NOTstorecaldata = false;    // if true, don't store cal data in the root file
    Bool_t NOTstorehitdata = false;    // if true, don't store hit data in the root file

    //Mapping and Calibration Files --------------------------
    TString califamapfilename = "../parameters/Califa_Mapping_3March2021.par";
    califamapfilename.ReplaceAll("//", "/");

    TString califacalfilename = "../parameters/Califa_CalPar_4March2021.par";
    califacalfilename.ReplaceAll("//", "/");

    //Sofia calibration files, dummy ones ---------------------

    //TString sofiacalfilename = "../parameters/dummy_sofia.par";
    //Sofia calibration files, use the one Antia gave me------
    //TString sofiacalfilename = "../parameters/CalibParam_antia.par";
    //now use params from experiment online
    TString sofiacalfilename = "/u/land/tobias_jenegger/s455_gabri_tobias/parameters/CalibParam_antia.par";

    // UCESB configuration ----------------------------------
    TString ntuple_options = "RAW";
    TString ucesb_dir = getenv("UCESB_DIR");
    TString upexps_dir = "/u/land/fake_cvmfs/9.13/upexps";
    TString ucesb_path;
    if (expId == 455)
    {
        ucesb_path = upexps_dir + "/202103_s455_jentob/202103_s455 --allow-errors --input-buffer=70Mi";

	//only with valid tpat
	//ucesb_path = upexps_dir + "/202103_s455_jentob_tpat/202103_s455 --allow-errors --max-events=10000 --input-buffer=70Mi";
    }
    else
    {
        std::cout << "Experiment was not selected!" << std::endl;
        gApplication->Terminate();
    }
    ucesb_path.ReplaceAll("//", "/");

    // Create source using ucesb for input ------------------
    EXT_STR_h101 ucesb_struct;
    R3BUcesbSource* source =
        new R3BUcesbSource(filename, ntuple_options, ucesb_path, &ucesb_struct, sizeof(ucesb_struct));
    source->SetMaxEvents(nev);

    // Definition of reader ---------------------------------

    //Unpackreader
    R3BUnpackReader* unpackreader =
	           new R3BUnpackReader(&ucesb_struct.unpack, offsetof(EXT_STR_h101, unpack));
    //Tpatreader
    R3BTrloiiTpatReader* unpacktpat =
	           new R3BTrloiiTpatReader(&ucesb_struct.unpacktpat, offsetof(EXT_STR_h101, unpacktpat));
    //CALIFA Febex Reader
    R3BCalifaFebexReader* unpackcalifa = new R3BCalifaFebexReader(&ucesb_struct.califa, offsetof(EXT_STR_h101, califa));

    //Master- Whiterabbitreader
    R3BWhiterabbitMasterReader* unpackWRM =
		   new R3BWhiterabbitMasterReader((EXT_STR_h101_WRMASTER*)&ucesb_struct.wrm, offsetof(EXT_STR_h101, wrm), 0x1000);
     
    //R3BWhiterabbitCalifaReader not used, as wr time information already in CALIFA Febex Reader

    R3BSofTofWReader* unpacktofw = new R3BSofTofWReader((EXT_STR_h101_SOFTOFW_t*)&ucesb_struct.tofw, offsetof(EXT_STR_h101, tofw));

    // R3BSofTwimReader
    R3BSofTwimReader* unpacktwim = new R3BSofTwimReader((EXT_STR_h101_SOFTWIM_t*)&ucesb_struct.twim, offsetof(EXT_STR_h101, twim));

    //R3BSofSciReader
    R3BSofSciReader* unpacksci = new R3BSofSciReader((EXT_STR_h101_SOFSCI_t*)&ucesb_struct.sci, offsetof(EXT_STR_h101, sci),NumSofSci);

    //unpack WR Sofia, this is maybe needed to determine time Sci->TOF
    R3BSofWhiterabbitReader* unpackWRSofia = new R3BSofWhiterabbitReader((EXT_STR_h101_WRSOFIA*)&ucesb_struct.wrsofia, offsetof(EXT_STR_h101, wrsofia),     sofiaWR_SE, sofiaWR_ME);

    // Add readers ------------------------------------------
    source->AddReader(unpackreader);
    source->AddReader(unpacktpat);
    source->AddReader(unpackcalifa);
    unpackcalifa->SetOnline(NOTstoremappeddata);
    source->AddReader(unpackWRM);
    unpackWRM->SetOnline(NOTstoremappeddata);
    source->AddReader(unpacktofw);
    unpacktofw->SetOnline(NOTstoremappeddata);
    unpacktwim->SetOnline(NOTstoremappeddata);
    source->AddReader(unpacktwim);
    unpacksci->SetOnline(NOTstoremappeddata);
    source->AddReader(unpacksci);
    unpackWRSofia->SetOnline(NOTstoremappeddata);
    source->AddReader(unpackWRSofia);


    // Create online run ------------------------------------
    FairRunOnline* run = new FairRunOnline(source);
    run->SetRunId(1);
    run->SetSink(new FairRootFileSink(outputFileName));

    // Runtime data base ------------------------------------
    FairRuntimeDb* rtdb = run->GetRuntimeDb();
    FairParAsciiFileIo* parIo1 = new FairParAsciiFileIo(); // Ascii file

    // Create analysis task ---------------------------------
    if (!Cal_level){
	parIo1->open(califamapfilename, "in"); //CALIFA mapping parameters
    	rtdb->setFirstInput(parIo1);
    	rtdb->print();
    }
    if (Cal_level)
    {
	//parIo1->open(califacalfilename, "in"); //CALIFA mapping and calibration parameters
	TList* parList1 = new TList();
	parList1->Add(new TObjString(sofiacalfilename));  //SOFIA mapping and calibration parameters
	parList1->Add(new TObjString(califacalfilename)); //CALIFA mapping and calibration parameters
	parIo1->open(parList1);
	rtdb->setFirstInput(parIo1);
	rtdb->print();

        // R3BCalifaMapped2CrystalCal ---
        R3BCalifaMapped2CrystalCal* Map2Cal = new R3BCalifaMapped2CrystalCal();
        Map2Cal->SetOnline(NOTstorecaldata);
        run->AddTask(Map2Cal);

        // R3BCalifaCrystalCal2Hit ---
        R3BCalifaCrystalCal2Hit* Cal2Hit = new R3BCalifaCrystalCal2Hit();
	Cal2Hit->SetCrystalThreshold(100.); // 100keV
        Cal2Hit->SelectGeometryVersion(2021);
	Cal2Hit->SetSquareWindowAlg(0.25,0.25);
       	Cal2Hit->SetOnline(NOTstorehitdata);
        run->AddTask(Cal2Hit);

	// R3BSofTwimMapped2Cal ---
	R3BSofTwimMapped2Cal* TwimMap2Cal = new R3BSofTwimMapped2Cal();
	TwimMap2Cal->SetOnline(NOTstorecaldata);
	TwimMap2Cal->SetExpId(expId);
	run->AddTask(TwimMap2Cal);
	
	//R3BSofTwimCal2Hit ---
	R3BSofTwimCal2Hit* TwimCal2Hit = new R3BSofTwimCal2Hit();
	TwimCal2Hit->SetOnline(NOTstorehitdata);
	TwimCal2Hit->SetExpId(expId);
	run->AddTask(TwimCal2Hit);	


        // --- Mapped 2 Tcal for SofSci
        R3BSofSciMapped2Tcal* SofSciMap2Tcal = new R3BSofSciMapped2Tcal();
        SofSciMap2Tcal->SetOnline(NOTstorecaldata);
        run->AddTask(SofSciMap2Tcal);

        // --- Tcal 2 SingleTcal for SofSci
        R3BSofSciTcal2SingleTcal* SofSciTcal2STcal = new R3BSofSciTcal2SingleTcal();
        SofSciTcal2STcal->SetOnline(NOTstorecaldata);
        run->AddTask(SofSciTcal2STcal);
        
	// --- SingleTcal 2 Cal for SofSci
        R3BSofSciSingleTcal2Cal* SofSciSTcal2Cal = new R3BSofSciSingleTcal2Cal();
        SofSciSTcal2Cal->SetOnline(NOTstorecaldata);
        run->AddTask(SofSciSTcal2Cal);
	
	// --- SingleTcal 2 Hit for SofSci
        R3BSofSciSingleTcal2Hit* SofSciSTcal2Hit = new R3BSofSciSingleTcal2Hit();
        SofSciSTcal2Hit->SetOnline(NOTstorehitdata);
        SofSciSTcal2Hit->SetCalParams(675., -1922.);//ToF calibration at Cave-C  tof and offset from spectra
        run->AddTask(SofSciSTcal2Hit);


	// --- Mapped 2 Tcal for SofTofW
	R3BSofTofWMapped2Tcal* SofTofWMap2Tcal = new R3BSofTofWMapped2Tcal();
	SofTofWMap2Tcal->SetOnline(NOTstorecaldata);
	run->AddTask(SofTofWMap2Tcal);


	// --- Tcal 2 SingleTcal for SofTofW
	R3BSofTofWTcal2SingleTcal* SofTofWTcal2STcal = new R3BSofTofWTcal2SingleTcal();
	SofTofWTcal2STcal->SetOnline(NOTstorecaldata);
	run->AddTask(SofTofWTcal2STcal);

	// --- SingleTcal 2 Hit for SofTofW
	R3BSofTofWSingleTCal2Hit* SofTofWSingleTcal2Hit = new R3BSofTofWSingleTCal2Hit();
	SofTofWSingleTcal2Hit->SetOnline(NOTstorehitdata);
	SofTofWSingleTcal2Hit->SetExpId(expId);
	SofTofWSingleTcal2Hit->SetTofLISE(37.48);
	run->AddTask(SofTofWSingleTcal2Hit);
    }

    run->SetSource(source);
    // Initialize -------------------------------------------
    run->Init();
    FairLogger::GetLogger()->SetLogScreenLevel("INFO");

    cout<<"\n\n"<<endl;
    // Run --------------------------------------------------
    run->Run((nev < 0) ? nev : 0, (nev < 0) ? 0 : nev);

    // -----   Finish   -------------------------------------
    cout << endl << endl;
    // Extract the maximal used memory an add is as Dart measurement
    // This line is filtered by CTest and the value send to CDash
    FairSystemInfo sysInfo;
    Float_t maxMemory = sysInfo.GetMaxMemory();
    cout << "MaxMemory: ";
    cout << maxMemory << endl;

    timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();
    Float_t cpuUsage = ctime / rtime;
    cout << "CPU used: " << cpuUsage << endl;

    cout << endl;
    std::cout << "Output file is " << outputFileName << std::endl;
    cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
    cout << "Macro finished successfully." << endl;
    gApplication->Terminate();
}
