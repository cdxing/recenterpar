#ifndef StProManger_h
#define StProManger_h

#include "TVector2.h"
#include "TString.h"
#include "StMessMgr.h"

class TProfile2D;
class TProfile;
class TH1F;

class StProManger
{
  public:
    StProManger();
    ~StProManger();

    // ReCenter Correction
    void InitReCenter();
    void FillTrackEast(TVector2 q2Vector, TVector2 q3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Float_t pt); // i = vertex pos/neg, j = eta_gap
    void FillTrackWest(TVector2 q2Vector, TVector2 q3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Float_t pt); // i = vertex pos/neg, j = eta_gap
    void FillTrackFull(TVector2 q2Vector, TVector2 q3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Float_t pt); // i = vertex pos/neg

    void WriteReCenter();

  private:

    // ReCenter Correction | x axis is RunIndex, y axis is Centrality
    // Event Plane method
    TProfile2D *p_mq2x_East_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq2y_East_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq2x_West_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq2y_West_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq2x_Full_EP[2];    // 0 = vertex pos/neg
    TProfile2D *p_mq2y_Full_EP[2];    // 0 = vertex pos/neg

    TProfile2D *p_mq3x_East_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq3y_East_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq3x_West_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq3y_West_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq3x_Full_EP[2];    // 0 = vertex pos/neg
    TProfile2D *p_mq3y_Full_EP[2];    // 0 = vertex pos/neg

    static TString mVStr[2];

    ClassDef(StProManger,1)
};

#endif
