#include "StCorrection.h"
#include "StConstants.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StMessMgr.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TF1.h"

Double_t Resolution_Full(Double_t *x_val, Double_t *par)
{
  Double_t y;
  Double_t chi = x_val[0];
  Double_t arg = chi*chi/4.0;
  Double_t norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

  y = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

  return y;
}

ClassImp(StCorrection)

TString StCorrection::mVStr[2] = {"pos","neg"};
TString StCorrection::mOrder[2] = {"2nd","3rd"};
TString StCorrection::mMethod[2] = {"EP","SP"};
//---------------------------------------------------------------------------------

StCorrection::StCorrection(Int_t energy)
{
  mEnergy = energy;
}

//---------------------------------------------------------------------------------

StCorrection::~StCorrection()
{
}

//---------------------------------------------------------------------------------

bool StCorrection::passTrackEtaEast(StPicoTrack *track, Int_t i, Int_t Mode) // neg || i = different eta_gap || Mode = 0 Event Plane Mode, Mode = 1 Flow Mode
{
  Float_t eta = track->pMom().PseudoRapidity();
  
  if(Mode == 0) // Event Plane Mode
  {
    // eta cut
    // eta_gap between two sub event plane is 2*mEta_Gap[i]
    if(!(eta > -1.0*TriFlow::mEtaMax && eta < -1.0*TriFlow::mEta_Gap[i]))
    {
      return kFALSE;
    }

    return kTRUE;
  }
  if(Mode == 1) // Flow Mode
  {
    // eta cut
    // eta_gap between two sub event plane is mEta_Gap[i]
    if(!(eta > -1.0*TriFlow::mEtaMax && eta < 0.0))
    {
      return kFALSE;
    }

    return kTRUE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------

bool StCorrection::passTrackEtaWest(StPicoTrack *track, Int_t i, Int_t Mode) // pos || i = different eta_gap || Mode = 0 Event Plane Mode, Mode = 1 Flow Mode
{
  Float_t eta = track->pMom().PseudoRapidity();

  if(Mode == 0) // Event Plane Mode
  {
    // eta cut
    // eta_gap between two sub event plane is 2*mEta_Gap[i]
    if(!(eta > TriFlow::mEta_Gap[i] && eta < TriFlow::mEtaMax))
    {
      return kFALSE;
    }

    return kTRUE;
  }
  if(Mode == 1) // Flow Mode
  {
    // eta cut
    // eta_gap between two sub event plane is mEta_Gap[i]
    if(!(eta > 0.0 && eta < TriFlow::mEtaMax))
    {
      return kFALSE;
    }

    return kTRUE;
  }

  return kTRUE;
}
//---------------------------------------------------------------------------------

bool StCorrection::passTrackFull(StPicoTrack *track) // Full Event Plane pt Cut
{
  // pt cut 0.2 - 2.0 GeV/c
  Float_t pt = track->pMom().Perp();
  if(!(pt > TriFlow::mPrimPtMin[mEnergy] && pt < TriFlow::mPrimPtMax))
  {
    return kFALSE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------
// Event Plane method and Scalar Product method
TVector2 StCorrection::calq2Vector(StPicoTrack *track)
{
  const Float_t phi = track->pMom().Phi();
  TVector2 q2Vector(0.0,0.0);

  const Float_t q2x = TMath::Cos(2.0*phi);
  const Float_t q2y = TMath::Sin(2.0*phi);
  q2Vector.Set(q2x,q2y);

  return q2Vector;
}

TVector2 StCorrection::calq3Vector(StPicoTrack *track)
{
  const Float_t phi = track->pMom().Phi();
  TVector2 q3Vector(0.0,0.0);

  const Float_t q3x = TMath::Cos(3.0*phi);
  const Float_t q3y = TMath::Sin(3.0*phi);
  q3Vector.Set(q3x,q3y);

  return q3Vector;
}

Float_t StCorrection::getWeight(StPicoTrack *track)
{
  Float_t pt = track->pMom().Perp();
  Float_t w;
  if(pt <= TriFlow::mPrimPtWeight)
  {
    w = pt;
  }
  if(pt > TriFlow::mPrimPtWeight)
  {
    w = TriFlow::mPrimPtWeight;
  }

  return w;
}
//------------------------------------------------------------------------------

