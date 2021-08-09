#ifndef StCut_h
#define StCut_h

#include "TObject.h"
#include "TString.h"
//#include <array>

class StPicoDst;
class StPicoTrack;
class StPicoBTofPidTraits;
class StRefMultCorr;

using namespace std;


class StCut : public TObject
{
  public:
    StCut(Int_t energy);
    ~StCut();

    bool isGoodTrigger(StPicoDst* );
    bool passEventCut(StPicoDst*);
    bool passTrackBasic(StPicoTrack*);
    bool passTrackEP(StPicoTrack*, float);

  private:
    Int_t mEnergy;

    ClassDef(StCut,1)
};
#endif
