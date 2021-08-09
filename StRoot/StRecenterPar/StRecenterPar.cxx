#include "StRecenterPar.h"
#include "StConstants.h"
#include "StCut.h"
#include "StProManger.h"
#include "StCorrection.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/run/run.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TVector3.h"
#include "TMath.h"
#include "StMessMgr.h"
#include <algorithm>
#include <array>

ClassImp(StRecenterPar)

    //StRefMultCorr* StRecenterPar::mRefMultCorr = NULL;
    //-----------------------------------------------------------------------------
    StRecenterPar::StRecenterPar(const char* name, StPicoDstMaker *picoMaker, char* outName )
: StMaker(name)
{
    mPicoDstMaker = picoMaker;
    mPicoDst = 0;
    mEnergy = 3;
    mOutName = outName;
    mOutName = mOutName + ".root";
}

//----------------------------------------------------------------------------- 
StRecenterPar::~StRecenterPar()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t StRecenterPar::Init() 
{

    //if(!mRefMultCorr)
    //{
    //    mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
    //}

    mCut = new StCut(mEnergy);
    mProManger = new StProManger();
    mCorrection = new StCorrection(mEnergy);

    mFile_ReCenterPar = new TFile(mOutName.Data(),"RECREATE");
    mFile_ReCenterPar->cd();
    mProManger->InitReCenter();

    return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StRecenterPar::Finish() 
{
    {
        mFile_ReCenterPar->cd();
        mProManger->WriteReCenter();
        mFile_ReCenterPar->Close();
    }
    return kStOK;
}

//----------------------------------------------------------------------------- 
void StRecenterPar::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
Int_t StRecenterPar::Make() 
{
    if(!mPicoDstMaker) 
    {
        LOG_WARN << " No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }

    mPicoDst = mPicoDstMaker->picoDst();
    if(!mPicoDst) 
    {
        LOG_WARN << " No PicoDst! Skip! " << endm;
        return kStWarn;
    }

    mPicoEvent = (StPicoEvent*)mPicoDst->event();
    if(!mPicoEvent)
    {
        LOG_WARN << " No PicoEvent! Skip! " << endm;
        return kStWarn;
    }

    // RefMult
    Int_t runId = mPicoEvent->runId();

    //cout << "runID = " << runId << endl;
    Int_t refMult = mPicoEvent->refMult();
    Float_t vz = mPicoEvent->primaryVertex().Z();
    Float_t vx = mPicoEvent->primaryVertex().X();
    Float_t vy = mPicoEvent->primaryVertex().Y();

    Float_t vzvpd = mPicoEvent->vzVpd();
    Int_t TOF_Mul = mPicoEvent->btofTrayMultiplicity();
    Int_t nMatchedToF = mPicoEvent->nBTOFMatch();
    Float_t zdcX = mPicoEvent->ZDCx();

    // vz sign
    Int_t vz_sign;
    if(vz > 0.0)
    {
        vz_sign = 0;
    }
    else
    {
        vz_sign = 1;
    }

    // runIndex
    const int runIndex = GetRunIndex(runId);
    // Event Cut
    if(mCut->passEventCut(mPicoDst)) // event cut
    {
         //mRefMultCorr->init(runId);
         // mRefMultCorr->initEvent(refMult, vz, zdcX);
         //const Int_t cent9 = mRefMultCorr->getCentralityBin9();
         //const Double_t reweight = mRefMultCorr->getWeight();
         const Int_t cent9 = Centrality(refMult);
         const Double_t reweight = 1.0;
         if(cent9 >  8 || cent9 < 0) return 0;

        const Int_t nTracks = mPicoDst->numberOfTracks();
        TVector3 mVertexPos = mPicoDst->event()->primaryVertex();
        float mField = mPicoEvent->bField();
        //cout << "nTracks = " << nTracks << endl;
        for(Int_t i = 0; i < nTracks; i++) // track loop
        {
            StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);
            //StPicoPhysicalHelix helix = track->helix(mField);
            //Float_t dca = helix.geometricSignedDistance(mVertexPos);
            Float_t dca=track->gDCA(mVertexPos).Mag();
            if(mCut->passTrackEP(track,dca)) // track cut   //shaowei
            {
                Float_t pt = track->pMom().Perp();
                Float_t eta = track->pMom().PseudoRapidity();

                if(mCorrection->passTrackFull(track))
                {
                    TVector2 q2Vector_Full = mCorrection->calq2Vector(track);
                    TVector2 q3Vector_Full = mCorrection->calq3Vector(track);
                    mProManger->FillTrackFull(q2Vector_Full,q3Vector_Full,cent9,runIndex,vz_sign,pt);
                }

                for(Int_t j = 0; j < 4; j++) // eta_gap loop
                {
                    if(mCorrection->passTrackEtaEast(track,j,0)) // neg eta sub
                    {
                        TVector2 q2Vector_East = mCorrection->calq2Vector(track);
                        TVector2 q3Vector_East = mCorrection->calq3Vector(track);
                        mProManger->FillTrackEast(q2Vector_East,q3Vector_East,cent9,runIndex,vz_sign,j,pt);
                    }
                    if(mCorrection->passTrackEtaWest(track,j,0)) // pos eta sub
                    {
                        TVector2 q2Vector_West = mCorrection->calq2Vector(track);
                        TVector2 q3Vector_West = mCorrection->calq3Vector(track);
                        mProManger->FillTrackWest(q2Vector_West,q3Vector_West,cent9,runIndex,vz_sign,j,pt);
                    }
                }
            }
        }
    }

    return kStOK;
}

int StRecenterPar::GetRunIndex(int runID)
{
    int runIndex=-999;
    for(int i=0; i<2704; i++)
    {
        if(runID==numbers[i])
        {
            runIndex=i;
        }
    }
    if(runIndex == -999) cout << "Run numbers are not found!!!" << endl;
    return runIndex;
}

int StRecenterPar::Centrality(int gRefMult )
{
    int centrality;
    //int centFull[9]={6,13,25,44,72,113,169,245,295}; // 27gev
    int centFull[9]={4, 9,17,30,50,78, 116,170,205}; // bes2 7.7gev
    if      (gRefMult>=centFull[8]) centrality=8;
    else if (gRefMult>=centFull[7]) centrality=7;
    else if (gRefMult>=centFull[6]) centrality=6;
    else if (gRefMult>=centFull[5]) centrality=5;
    else if (gRefMult>=centFull[4]) centrality=4;
    else if (gRefMult>=centFull[3]) centrality=3;
    else if (gRefMult>=centFull[2]) centrality=2;
    else if (gRefMult>=centFull[1]) centrality=1;
    else if (gRefMult>=centFull[0]) centrality=0;
    else centrality = 9;

    return centrality;
}
