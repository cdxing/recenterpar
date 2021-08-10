#ifndef PicoAnalyzer__
#define PicoAnalyzer__

#include "TObject.h"

class TTree;
class TBranch;
class TLeaf;
class TClonesArray;
class StPicoEvent;
class TH1D;
class TH2D;
class TNtuple;
class TFile;
class StEpdGeom;
class StBbcGeom;
class TRandom3;


class PicoAnalyzer : public TObject {
 public:
  PicoAnalyzer();
  ~PicoAnalyzer();

  void SetPicoDst(TTree*);
  short Init();
  short Make(int iEvent);
  short Finish();

 private:

  void ReadLeaves();   // explictly reads the relevant leaves, from the StPicoEvent object

  double GetBbcPmtPhi(short PmtId);


  int nEvents;

  StEpdGeom* mEpdGeom;
  StBbcGeom* mBbcGeom;

  double mNmipQtB;  // ADC value of MIP peak on rings 6-16 (read out thru QT32Bs)
  double mNmipQtC;  // ADC value of MIP peak on rings 1-5  (read out thru QT32Cs)

  double mnMipThreshold;  // low-signal threshold, to cut out noise basically.

  // the tree, branches and leaves
  TTree* mPicoDst;
  TBranch* mEpdBranch;
  TBranch* mBbcBranch;
  TBranch* mEventBranch;
  // dealing directly with the TLeaf objects is the only way I can get access to the info in the StPicoEvent object
  // a little clumsy but it has its "elegance," I suppose.  And it's pretty fast
  TLeaf* mRefMultNegLeaf;
  TLeaf* mRefMultPosLeaf;
  TLeaf* mVzLeaf;

  // the data objects
  TClonesArray* mEpdHits;
  TClonesArray* mBbcHits;
  int mRefMult;
  double mVz;

  // histograms
  /*
  TH1D* mVzHist;
  TH1D* mRefMultHist;
  TH1D* mDnDetaHist;
  TH1D* mCos1DeltaPsi1;
  TH1D* mCos2DeltaPsi2;
  */
  /*
  TH2D* mPhi1EWhist;
  TH2D* mPhi2EWhist;
  TH2D* mRMvsEpdSumMip;
  TH2D* mBbcAdcSumVsEpdSumMip;
  */

  TH1D* mNmipDists[2][12][31];

  TH1D* mHisto1D[40];
  TH2D* mHisto2D[40];

  TNtuple* mQ1vectorNtuple;
  TNtuple* mQ2vectorNtuple;

  TFile* mHistoFile;
  TFile* mQ1NtupleFile;
  TFile* mQ2NtupleFile;

  TRandom3* mRan;


//  StPicoEvent* mPicoEvent;  // too bad, can't read this due to bloody StThreeVectorF


  ClassDef(PicoAnalyzer, 0)


};


#endif
