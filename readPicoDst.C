
#include <TSystem>

class StMaker;
class StChain;
class StPicoEvent;
class StPicoDstMaker;
class StRefMultCorr;

StChain *chain;
void readPicoDst(const Char_t *inputFile="test.list", Char_t *outputFile="test")
{
        Int_t nEvents = 1000000;
//Load all the System libraries

        gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();

    gSystem->Load("StRefMultCorr");
    gSystem->Load("StPicoEvent");
	gSystem->Load("StPicoDstMaker");
	gSystem->Load("StRecenterPar");
  gSystem->Load("StMyEpdAna");
	chain = new StChain();

	StPicoDstMaker *picoMaker = new StPicoDstMaker(2,inputFile,"picoDst");
	StRecenterPar *recenter = new StRecenterPar("recenter",picoMaker,outputFile);
  PicoAnalyzer *epdpar = new PicoAnalyzer();

	chain->Init();
	cout<<"chain->Init();"<<endl;
	int total = picoMaker->chain()->GetEntries();
	cout << " Total entries = " << total << endl;
	if(nEvents>total) nEvents = total;
	for (Int_t i=0; i<nEvents; i++){

        if(i != 0 && i%50 == 0)
        {
            cout << "." << flush;
        }
        if(i != 0 && i%500 == 0)
        {
            Float_t event_percent = 100.0*i/nEvents;
            cout << " " << i << "(" << event_percent << "%)" << "\n" << " ==> Processing data " << flush;
        }


        chain->Clear();
		int iret = chain->Make(i);

		if (iret) { cout << "Bad return code!" << iret << endl; break;}

		total++;

	}

	cout << "****************************************** " << endl;
	cout << "Work done... now its time to close up shop!"<< endl;
	cout << "****************************************** " << endl;
	chain->Finish();
	cout << "****************************************** " << endl;
	cout << "total number of events  " << nEvents << endl;
	cout << "****************************************** " << endl;

	delete chain;


}
