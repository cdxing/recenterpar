#ifndef StRecenterPar_h
#define StRecenterPar_h

#include "StMaker.h"
#include "TString.h"
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StCut;
class StProManger;
class StCorrection;
class StRefMultCorr;

// run QA
class TProfile;
class TH1F;
class TH2F;

class StRecenterPar : public StMaker {
    public:
        StRecenterPar(const char *name, StPicoDstMaker *picoMaker, char *outName);
        virtual ~StRecenterPar();

        virtual Int_t Init();
        virtual Int_t Make();
        virtual void  Clear(Option_t *opt="");
        virtual Int_t Finish();

        int Centrality(int gRefMult);
        int GetRunIndex(int runId);
        bool isGoodTrigger(StPicoEvent const*) const;

    private:
        static StRefMultCorr *mRefMultCorr;
        StPicoDstMaker *mPicoDstMaker;
        StPicoDst      *mPicoDst;
        StPicoEvent *mPicoEvent;
        StCut *mCut;
        StProManger *mProManger;
        StCorrection *mCorrection;

        Int_t mEnergy;
        TString mOutName;
        TFile *mFile_ReCenterPar;


        ClassDef(StRecenterPar, 1)
};

#endif
