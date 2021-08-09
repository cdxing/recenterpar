#ifndef StCorrection_h
#define StCorrection_h

#include "TObject.h"
#include "TVector2.h"
#include "TString.h"

class StPicoDst;
class StPicoTrack;
class TProfile2D;
class TFile;
class TNtuple;

class StCorrection : public TObject
{
  public:
    StCorrection(Int_t energy);
    ~StCorrection();


    // ReCenter Correction
    bool passTrackEtaEast(StPicoTrack*,Int_t,Int_t);
    bool passTrackEtaWest(StPicoTrack*,Int_t,Int_t);
    bool passTrackFull(StPicoTrack*);

    TVector2 calq2Vector(StPicoTrack*);
    TVector2 calq3Vector(StPicoTrack*);
    Float_t getWeight(StPicoTrack*);

  private:
    Int_t    mEnergy;

    static TString mVStr[2];
    static TString mOrder[2];
    static TString mMethod[2];

  ClassDef(StCorrection,1)
};

#endif
