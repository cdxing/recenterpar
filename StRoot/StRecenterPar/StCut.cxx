#include "StCut.h"
#include "StConstants.h"
#include "StRoot/StPicoEvent/StPicoDst.h"  // shaowei  18c
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StMessMgr.h"

ClassImp(StCut)

    //StRefMultCorr* StCut::mRefMultCorr = NULL; // shaowei
    //---------------------------------------------------------------------------------

StCut::StCut(Int_t energy)
{
    mEnergy = energy;
}

//---------------------------------------------------------------------------------

StCut::~StCut()
{
}

//---------------------------------------------------------------------------------

bool StCut::isGoodTrigger(StPicoDst *pico)
{
    array<unsigned int, 6> const triggers = {
        610001,
        610011,
        610021,
        610031,
        610041,
        610051
    };
    StPicoEvent *event = pico->event();
    for(unsigned i=0; i<event->triggerIds().size(); i++)
    {
        //cout << "trigger ID = " << event->triggerIds()[i] << endl;
    }
    for(auto trg: triggers)
    {
        if(event->isTrigger(trg)) return kTRUE;
    }
}
bool StCut::passEventCut(StPicoDst *pico)
{
    //if(!isGoodTrigger(pico)) return kFALSE;
    StPicoEvent *event = pico->event();
    if(!event)
    {
        return kFALSE;
    }

    // initialize StRefMultCorr
    const Float_t vx = event->primaryVertex().X();
    const Float_t vy = event->primaryVertex().Y();
    const Float_t vz = event->primaryVertex().Z();
    // event vertex cut
    // vz cut
    if(fabs(vz) > 70.0)
    {
        return kFALSE;
    }
    if(event->btofTrayMultiplicity()<2)return kFALSE;
    // vr cut
    if(sqrt(vx*vx+vy*vy) > 2.0)
    {
        return kFALSE;
    }
    return kTRUE;
}


bool StCut::passTrackBasic(StPicoTrack *track)
{
    // nHitsFit cut
    if(track->nHitsFit() < TriFlow::mHitsFitTPCMin)
    {
        return kFALSE;
    }

    // nHitsRatio cut
    if(track->nHitsMax() <= TriFlow::mHitsMaxTPCMin)
    {
        return kFALSE;
    }
    if((Float_t)track->nHitsFit()/(Float_t)track->nHitsMax() < TriFlow::mHitsRatioTPCMin)
    {
        return kFALSE;
    }

    // eta cut
    Float_t eta = track->pMom().PseudoRapidity();
    if(fabs(eta) > TriFlow::mEtaMax)
    {
        return kFALSE;
    }

    return kTRUE;
}

//---------------------------------------------------------------------------------

bool StCut::passTrackEP(StPicoTrack *track, float dca)
{
    if(!track) return kFALSE;
    if(!passTrackBasic(track)) return kFALSE;

    // dca cut for event plane reconstruction: 200GeV = 3.0, BES = 1.0
    //if(dca > TriFlow::mDcaEPMax[mEnergy])   // change by shaowei
    if(dca > 1.0)    // old is 1,  change to 3 for test
    {
        return kFALSE;
    }

    // pt cut 0.2 - 2.0 GeV/c
    Float_t pt = track->pMom().Perp();
    Float_t p  = track->pMom().Mag();
    if(!(pt > TriFlow::mPrimPtMin[mEnergy] && pt < TriFlow::mPrimPtMax && p < TriFlow::mPrimMomMax))
    {
        return kFALSE;
    }

    return kTRUE;
}
//---------------------------------------------------------------------------------
