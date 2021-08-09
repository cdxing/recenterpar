/*
  Very minimalist example code
*/


#include "tools/PicoAnalyzer.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StPicoEvent/StPicoBbcHit.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StEpdUtil/StEpdGeom.h"
#include "StEpdUtil/StBbcGeom.h"

#include "TChain.h"
//#include "TTree.h"
//#include "TLeaf.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TRandom3.h"

#include "StEpdUtil/StEpdEpFinder.h"

#include <iostream>
#include <fstream>
using namespace std;




//=================================================
PicoAnalyzer::PicoAnalyzer() : mNmipQtB(115), mNmipQtC(160), mnMipThreshold(0.3), mPicoDst(0), mEpdHits(0), mBbcHits(0), mEventClonesArray(0), mQ1NtupleFile(0), mQ2NtupleFile(0), mRunId(0), mRunCollisionSystem(0)
{
  mEpdGeom = new StEpdGeom;
  mBbcGeom = new StBbcGeom;
  mRan = new TRandom3;
  mRan->GetSeed();
}

//=================================================
PicoAnalyzer::~PicoAnalyzer(){
  /* no-op */
}


//=================================================
//void PicoAnalyzer::SetPicoDst(TTree* PicoDst){
void PicoAnalyzer::SetPicoDst(TChain* PicoDst){
  mPicoDst        = PicoDst;

  mEpdHits = new TClonesArray("StPicoEpdHit");
  mBbcHits = new TClonesArray("StPicoBbcHit");
  mTracks  = new TClonesArray("StPicoTrack");
  mEventClonesArray = new TClonesArray("StPicoEvent");

  mPicoDst->SetBranchStatus("*",0);         // turns OFF all branches  (speeds it up :-)
  unsigned int found;
  mPicoDst->SetBranchStatus("EpdHit*",1,&found);   // note you need the asterisk
  cout << "EpdHit Branch returned found= " << found << endl;
  mPicoDst->SetBranchAddress("EpdHit",&mEpdHits);

  mPicoDst->SetBranchStatus("Event*",1,&found);
  cout << "Event Branch returned found= " << found << endl;
  mPicoDst->SetBranchAddress("Event",&mEventClonesArray);


  mPicoDst->SetBranchStatus("Track*",1,&found);
  cout << "Track Branch returned found= " << found << endl;
  mPicoDst->SetBranchAddress("Track",&mTracks);
}

//=================================================
short PicoAnalyzer::Init(){

  // ------------------- for the EP finder ------------------
  mEpFinder = new StEpdEpFinder();
  mEpFinder->SetnMipThreshold(0.3);
  mEpFinder->SetMaxTileWeight(2.0);
  mEpFinder->SetEpdHitFormat(2);     // 2=pico

  // Here, I set eta-based weights
  //    This is only to be used for Au+Au 27 GeV
  // https://drupal.star.bnl.gov/STAR/blog/lisa/optimizing-ep1-resolution-au27au-ring-dependent-weights-eta-dependent-weights
  double lin[9] = {-1.950, -1.900, -1.850, -1.706, -1.438, -1.340, -1.045, -0.717, -0.700};
  double cub[9] = {0.1608, 0.1600, 0.1600, 0.1595, 0.1457, 0.1369, 0.1092, 0.0772, 0.0700};
  TH2D wt("Order1etaWeight","Order1etaWeight",100,1.5,6.5,10,0,10);
  for (int ix=1; ix<101; ix++){
    for (int iy=1; iy<11; iy++){
      double eta = wt.GetXaxis()->GetBinCenter(ix);
      wt.SetBinContent(ix,iy,lin[iy-1]*eta+cub[iy-1]*pow(eta,3));
    }
  }
  mEpFinder->SetEtaWeights(1,wt);

  // --------------------------------------------------------

  mHistoFile = new TFile("EpdHists.root","RECREATE");

  // Miscellaneous One-dimensional histograms
  mHisto1D[0] = new TH1D("Vz","Vz",100,-80,80);
  mHisto1D[1] = new TH1D("RefMult","RefMult",100,-10,600);
  mHisto1D[2] = new TH1D("dNdeta","dNdeta",100,-5.5,5.5);
  mHisto1D[3] = new TH1D("Cos2DeltaPsi2","Cos2DeltaPsi2",1000,-1,1);
  mHisto1D[4] = new TH1D("Cos1DeltaPsi1","Cos1DeltaPsi1",1000,-1,1);
  mHisto1D[5] = new TH1D("BbcCos1DeltaPsi1","BbcCos1DeltaPsi1",1000,-1,1);
  mHisto1D[6] = new TH1D("BbcCos2DeltaPsi2","BbcCos2DeltaPsi2",1000,-1,1);
  //  mHisto1D[7] = new TH1D("TriggerIDs","TriggerIDs",1000000,0,1000000);      // hmm... kind of a big histogram...
  mHisto1D[8]  = new TH1D("FullEventPhi1","Full Event #Psi_{1}",100,0.0,2.0*TMath::Pi());
  mHisto1D[9]  = new TH1D("FullEventPhi2","Full Event #Psi_{2}",100,0.0,2.0*TMath::Pi()/2.0);
  mHisto1D[10] = new TH1D("FullEventPhi3","Full Event #Psi_{3}",100,0.0,2.0*TMath::Pi()/3.0);

  // Miscellaneous Two-dimensional histograms
  mHisto2D[0] = new TH2D("Phi2WestVsPhiEast","Phi2WestVsPhiEast",40,-TMath::Pi()/2.0,TMath::Pi()/2.0,40,-TMath::Pi()/2.0,TMath::Pi()/2.0);
  mHisto2D[0]->GetXaxis()->SetTitle("#phi_{2,East}");     mHisto2D[0]->GetYaxis()->SetTitle("#phi_{2,West}");
  mHisto2D[1] = new TH2D("Phi1WestVsPhiEast","Phi1WestVsPhiEast",40,-TMath::Pi(),TMath::Pi(),40,-TMath::Pi(),TMath::Pi());
  mHisto2D[1]->GetXaxis()->SetTitle("#phi_{1,East}");     mHisto2D[1]->GetYaxis()->SetTitle("#phi_{1,West}");
  mHisto2D[2] = new TH2D("RefMultVsEpdMipSum","RefMultVsEpdMipSum",100,0,5000,100,0,500);
  mHisto2D[2]->GetXaxis()->SetTitle("Epd nMIP Sum");      mHisto2D[2]->GetYaxis()->SetTitle("RefMult");
  mHisto2D[3] = new TH2D("BbcAdcSumVsEpdMipSum","BbcAdcSumVsEpdMipSum",100,0,5000,100,0,50000);
  mHisto2D[3]->GetXaxis()->SetTitle("BBC ADC Sum");      mHisto2D[3]->GetYaxis()->SetTitle("RefMult");
  mHisto2D[4] = new TH2D("BbcPhi1WestVsPhiEast","BbcPhi1WestVsPhiEast",40,-TMath::Pi(),TMath::Pi(),40,-TMath::Pi(),TMath::Pi());
  mHisto2D[4]->GetXaxis()->SetTitle("#phi_{1,East}");     mHisto2D[4]->GetYaxis()->SetTitle("#phi_{1,West}");
  mHisto2D[5] = new TH2D("BbcPhi2WestVsPhiEast","BbcPhi2WestVsPhiEast",40,-TMath::Pi()/2.0,TMath::Pi()/2.0,40,-TMath::Pi()/2.0,TMath::Pi()/2.0);
  mHisto2D[5]->GetXaxis()->SetTitle("#phi_{2,East}");     mHisto2D[5]->GetYaxis()->SetTitle("#phi_{2,West}");
  mHisto2D[6] = new TH2D("VyVx","V_{y} vs V_{x}",400,-1,1,400,1,1);
  mHisto2D[6]->GetXaxis()->SetTitle("V_{x} (cm)");        mHisto2D[6]->GetYaxis()->SetTitle("V_{y} (cm)");
  mHisto2D[7] = new TH2D("nGlobalTracksVsRefMult","Global Mult versus RefMult",100,0,600,100,0,3000);
  mHisto2D[7]->GetXaxis()->SetTitle("RefMult");           mHisto2D[7]->GetYaxis()->SetTitle("nGlobal Tracks");
  mHisto2D[8] = new TH2D("BtofMultVsRefMult","Btof Mult versus RefMult",100,0,600,100,0,3000);
  mHisto2D[8]->GetXaxis()->SetTitle("RefMult");           mHisto2D[8]->GetYaxis()->SetTitle("Btof mult");
  mHisto2D[9] = new TH2D("BtofMultVsGlobalMult","Btof Mult versus Global Mult",100,0,3000,100,0,3000);
  mHisto2D[9]->GetXaxis()->SetTitle("GlobalMult");        mHisto2D[9]->GetYaxis()->SetTitle("Btof mult");
  mHisto2D[10] = new TH2D("EpdWeightEast","Epd East Tile Weight",500,-100.0,100.0,500,-100.0,100.0);
  mHisto2D[10]->GetXaxis()->SetTitle("x (cm)");      mHisto2D[10]->GetYaxis()->SetTitle("y (cm)");
  mHisto2D[11] = new TH2D("EpdWeightWest","Epd West Tile Weight",500,-100.0,100.0,500,-100.0,100.0);
  mHisto2D[11]->GetXaxis()->SetTitle("x (cm)");      mHisto2D[11]->GetYaxis()->SetTitle("y (cm)");
  /*
    For histograms 12,13,14,15:
    Unfortunately, unless you want rotatable lego/surface plots (which probably you don't), the way you
    have to draw those "polar" histograms is as follows
    gPad->DrawFrame(-17,-17,17,17);
    EpdWeightWestPolar->Draw("colz,pol,same");
  */
  mHisto2D[12] = new TH2D("EpdWeightEastPolar","Epd Weight East - Plot as polar",24,0.0,2.0*TMath::Pi(),16,0.5,16.5);   // should plot in polar coordinates.  Radial coordinate is row.
  mHisto2D[13] = new TH2D("EpdWeightWestPolar","Epd Weight West - Plot as polar",24,0.0,2.0*TMath::Pi(),16,0.5,16.5);   // should plot in polar coordinates.  Radial coordinate is row.
  mHisto2D[14] = new TH2D("EpdWeightEastPolarPhiAveraged","Epd Weight East Phi-averaged - Plot as polar",24,0.0,2.0*TMath::Pi(),16,0.5,16.5);   // should plot in polar coordinates.  Radial coordinate is row.
  mHisto2D[15] = new TH2D("EpdWeightWestPolarPhiAveraged","Epd Weight West Phi-averaged - Plot as polar",24,0.0,2.0*TMath::Pi(),16,0.5,16.5);   // should plot in polar coordinates.  Radial coordinate is row.
  //
  mHisto2D[16] = new TH2D("EastEpdRotatedByEastPsi1","Epd East Hit Weight Rotated by -#Psi_{1,E}",500,-100.0,100.0,500,-100.0,100.0);
  mHisto2D[16]->GetXaxis()->SetTitle("x (cm)");      mHisto2D[16]->GetYaxis()->SetTitle("y (cm)");
  mHisto2D[17] = new TH2D("WestEpdRotatedByWestPsi1","Epd West Hit Weight Rotated by -#Psi_{1,W}",500,-100.0,100.0,500,-100.0,100.0);
  mHisto2D[17]->GetXaxis()->SetTitle("x (cm)");      mHisto2D[17]->GetYaxis()->SetTitle("y (cm)");
  mHisto2D[18] = new TH2D("RandomlyRotatedEpdWeightEast","Epd East Hit Weight Randomly Rotated",500,-100.0,100.0,500,-100.0,100.0);  // this is simply nice for normalizing the above histogram offline
  mHisto2D[18]->GetXaxis()->SetTitle("x (cm)");      mHisto2D[18]->GetYaxis()->SetTitle("y (cm)");
  mHisto2D[19] = new TH2D("RandomlyRotatedEpdWeightWest","Epd West Hit Weight Randomly Rotated",500,-100.0,100.0,500,-100.0,100.0);  // this is simply nice for normalizing the above histogram offline
  mHisto2D[19]->GetXaxis()->SetTitle("x (cm)");      mHisto2D[19]->GetYaxis()->SetTitle("y (cm)");
  mHisto2D[20] = new TH2D("EastEpdRotatedByWestPsi1","Epd East Hit Weight Rotated by -#Psi_{1,WEST}",500,-100.0,100.0,500,-100.0,100.0);  // note this is East EPD rotated by _WEST_ angle (more interesting - no "autocorrelation"))
  mHisto2D[20]->GetXaxis()->SetTitle("x (cm)");      mHisto2D[20]->GetYaxis()->SetTitle("y (cm)");
  mHisto2D[21] = new TH2D("WestEpdRotatedByEastPsi1","Epd West Hit Weight Rotated by -#Psi_{1,EAST}",500,-100.0,100.0,500,-100.0,100.0);  // note this is West EPD rotated by _EAST_ angle (more interesting - no "autocorrelation"))
  mHisto2D[21]->GetXaxis()->SetTitle("x (cm)");      mHisto2D[21]->GetYaxis()->SetTitle("y (cm)");

  mHisto2D[22] = new TH2D("EpdWeightEastRotatedFullPsiPolar","Epd Weight East - Rotated by *full event* EP1 - Plot as polar",24,0.0,2.0*TMath::Pi(),16,0.5,16.5);   // should plot in polar coordinates.  Radial coordinate is row.
  mHisto2D[23] = new TH2D("EpdWeightWestRotatedFullPsiPolar","Epd Weight West - Rotated by *full event* EP1 - Plot as polar",24,0.0,2.0*TMath::Pi(),16,0.5,16.5);   // should plot in polar coordinates.  Radial coordinate is row.

  mHisto2D[24] = new TH2D("EpdWeightEastRotatedByWestPolar","Epd Weight East - Rotated by West EP1 - Plot as polar",24,0.0,2.0*TMath::Pi(),16,0.5,16.5);   // should plot in polar coordinates.  Radial coordinate is row.
  mHisto2D[25] = new TH2D("EpdWeightWestRotatedByEastPolar","Epd Weight West - Rotated by East EP1 - Plot as polar",24,0.0,2.0*TMath::Pi(),16,0.5,16.5);   // should plot in polar coordinates.  Radial coordinate is row.


  // One-dimensional profiles
  // raw directed flow (v1)  :-)
  for (int icent=0; icent<9; icent++){
    for (int ewFull=0; ewFull<3; ewFull++){
      for (int order=1; order<4; order++){
	mEpdRawVn[icent][ewFull][order-1]  = new TH1D(Form("EpdRawV%dProfile_EWFull%d_CentID%d",order,ewFull,icent),
						      Form("Raw v%d EWFull=%d - CentId=%d",order,ewFull,icent),100,-5.2,5.2);
	mEpdRawVn[icent][ewFull][order-1]->Sumw2();
      }
      mEpdWeight[icent][ewFull] = new TH1D(Form("EpdWeightProfile_EWFull%d_CentID%d",ewFull,icent),
					   Form("This just needed for normalization at end - EWFull=%d CentId=%d",ewFull,icent),100,-5.2,5.2); // only needed for normalization at the end, then discarded
    }
      //    mTpcRawV1_FullPsiRotate[icent] = new TProfile(Form("TpcRawV1Profile_CentID%d",icent),
      //    Form("Raw v1 from TPC - CentId=%d - This is relative to *full event* EPD plane",icent),100,-5.2,5.2);
  }


  for (int isub=0; isub<mNumberOfEpdSubEvents; isub++){
    for (int order=1; order<5; order++){
      mEastLongDecorProfile[isub][order-1] = new TProfile(Form("EastLongDecorProf_order%d_bin%d",order,isub),Form("East EPD subevent %d - EP correlation of order %d",isub,order),mNumberOfTpcSubEvents,-1.4,1.4);
      mWestLongDecorProfile[isub][order-1] = new TProfile(Form("WestLongDecorProf_order%d_bin%d",order,isub),Form("West EPD subevent %d - EP correlation of order %d",isub,order),mNumberOfTpcSubEvents,-1.4,1.4);
    }
  }

  // Miscellaneous Two-dimensional profiles
  mProfile2D[0] = new TProfile2D("cosDeltaPhi_ij_East","East <cos(#phi_{i}-#phi_{j})>",16,0.5,16.5,16,0.5,16.5);
  mProfile2D[0]->GetXaxis()->SetTitle("Ring i");           mProfile2D[0]->GetYaxis()->SetTitle("Ring j");
  mProfile2D[1] = new TProfile2D("cosDeltaPhi_ij_West","West <cos(#phi_{i}-#phi_{j})>",16,0.5,16.5,16,0.5,16.5);
  mProfile2D[1]->GetXaxis()->SetTitle("Ring i");           mProfile2D[1]->GetYaxis()->SetTitle("Ring j");


  for (int ipmt=0; ipmt<16; ipmt++){
    mBbcXY[0][ipmt] = new TH2D(Form("BBC%dEastNotFires",ipmt+1),Form("When PMT%d does not fire on BBC East",ipmt+1),500,-80,80,500,-80,80);
    mBbcXY[1][ipmt] = new TH2D(Form("BBC%dWestNotFires",ipmt+1),Form("When PMT%d does not fire on BBC West",ipmt+1),500,-80,80,500,-80,80);
  }

  for (int icent=0; icent<7; icent++) mRefMultHist[icent] = new TH1D(Form("RefMult%d",icent),Form("RefMult %d-%d%s",70-10*icent,80-10*icent,"%"),600,-0.5,599.5);
  mRefMultHist[7] = new TH1D("RefMult7","RefMult 5-10%",600,-0.5,599.5);
  mRefMultHist[8] = new TH1D("RefMult8","RefMult 0-5%",600,-0.5,599.5);

  // histograms and profiles of event planes
  // EPD "correction level" (last index):  0=raw; 1=shifted
  // BBC "correction level" (last index):  0=raw; 1=gain; 2=shifted
  for (int n=1; n<4; n++){
    double PhiMax = 2.0*TMath::Pi()/((double)n);
    mEpdEwPsi[n-1][0]  = new TH2D(Form("EpdEwPsi%draw",n),  Form("EPD raw #Psi_{%d,W} vs #Psi_{%d,E}",n,n),  100,0.0,PhiMax,100,0.0,PhiMax);
    mEpdEwPsi[n-1][1]  = new TH2D(Form("EpdEwPsi%dshift",n),Form("EPD shift #Psi_{%d,W} vs #Psi_{%d,E}",n,n),100,0.0,PhiMax,100,0.0,PhiMax);
    //
    mBbcEwPsi[n-1][0]  = new TH2D(Form("BbcEwPsi%draw",n),  Form("BBC raw #Psi_{%d,W} vs #Psi_{%d,E}",n,n),  100,0.0,PhiMax,100,0.0,PhiMax);
    mBbcEwPsi[n-1][1]  = new TH2D(Form("BbcEwPsi%dgain",n), Form("BBC gain #Psi_{%d,W} vs #Psi_{%d,E}",n,n), 100,0.0,PhiMax,100,0.0,PhiMax);
    mBbcEwPsi[n-1][2]  = new TH2D(Form("BbcEwPsi%dshift",n),Form("BBC shift #Psi_{%d,W} vs #Psi_{%d,E}",n,n),100,0.0,PhiMax,100,0.0,PhiMax);
    //
    mEpdEwPsi_midCentral[n-1][0]  = new TH2D(Form("EpdEwPsi%drawMidCent",n),  Form("MidCentral EPD raw #Psi_{%d,W} vs #Psi_{%d,E}",n,n),  100,0.0,PhiMax,100,0.0,PhiMax);
    mEpdEwPsi_midCentral[n-1][1]  = new TH2D(Form("EpdEwPsi%dshiftMidCent",n),Form("MidCentral EPD shift #Psi_{%d,W} vs #Psi_{%d,E}",n,n),100,0.0,PhiMax,100,0.0,PhiMax);

    for (int icor=0; icor<2; icor++){
      mEpdEwPsi[n-1][icor]->GetXaxis()->SetTitle(Form("#Psi_{%d,E}",n));                 mEpdEwPsi[n-1][icor]->GetYaxis()->SetTitle(Form("#Psi_{%d,W}",n));
      mEpdEwPsi_midCentral[n-1][icor]->GetXaxis()->SetTitle(Form("#Psi_{%d,E}",n));      mEpdEwPsi_midCentral[n-1][icor]->GetYaxis()->SetTitle(Form("#Psi_{%d,W}",n));
    }
    for (int icor=0; icor<3; icor++){
      mBbcEwPsi[n-1][icor]->GetXaxis()->SetTitle(Form("#Psi_{%d,E}",n));      mBbcEwPsi[n-1][icor]->GetYaxis()->SetTitle(Form("#Psi_{%d,W}",n));
    }
    //
    mEpdAveCos[n-1][0] = new TProfile(Form("EpdCos%draw",n),  Form("EPD raw #LT cos(%d#Psi_{%d,W}-%d#Psi_{%d,E}) #GT",n,n,n,n)  ,9,-0.5,8.5);
    mEpdAveCos[n-1][1] = new TProfile(Form("EpdCos%dshift",n),Form("EPD shift #LT cos(%d#Psi_{%d,W}-%d#Psi_{%d,E}) #GT",n,n,n,n),9,-0.5,8.5);
    //
    mBbcAveCos[n-1][0] = new TProfile(Form("BbcCos%draw",n),  Form("BBC raw #LT cos(%d#Psi_{%d,W}-%d#Psi_{%d,E}) #GT",n,n,n,n),  9,-0.5,8.5);
    mBbcAveCos[n-1][1] = new TProfile(Form("BbcCos%dgain",n), Form("BBC gain #LT cos(%d#Psi_{%d,W}-%d#Psi_{%d,E}) #GT",n,n,n,n),9,-0.5,8.5);
    mBbcAveCos[n-1][2] = new TProfile(Form("BbcCos%dshift",n),Form("BBC shift #LT cos(%d#Psi_{%d,W}-%d#Psi_{%d,E}) #GT",n,n,n,n),9,-0.5,8.5);
  }

  ReadInSystems();  // read in a text file of what system is being collided in what run.

  cout << "Init done\n";
  return 0;
}

//==============================================
// some things might need resetting when there is a new run
void PicoAnalyzer::NewRun(int runId){
  mRunCollisionSystem = WhichSystem(runId);
  mRunId = runId;
}



//=================================================
short PicoAnalyzer::Make(int iEvent){
  //----------------- get data --------------
  mPicoDst->GetEntry(iEvent);
  StPicoEvent* event = (StPicoEvent*)((*mEventClonesArray)[0]);
  if (event->runId()!=mRunId){
    NewRun(event->runId());        // some things should be reset when there is a new run loaded
    cout << "New run detected: " << mRunId << " and it is collision system #" << mRunCollisionSystem << endl;
  }
  //----- done getting data; have fun! ------

  StThreeVectorF primaryVertex = event->primaryVertex();
  TVector3 PV(primaryVertex.x(),primaryVertex.y(),primaryVertex.z());

  double pi = TMath::Pi();

  int mRefMult = event->refMult();
  int mRunId = event->runId();
  double mVzVpd = event->vzVpd() + 36.4046 - 3.13415;  // numbers from Rosi's email of 3may2018
  int CentId = FindCent(mRefMult);   // returns an integer between 0 (70-80%) and 8 (0-5%)
  int RunYear = 19;
  int RunDay = floor( (mRunId - RunYear*pow(10,6))/pow(10,3) );
  int DayBinId = RunDay-89;
  if (CentId<0) return 0;            // 80-100% - very peripheral


  mHisto1D[0]->Fill(PV.Z());
  mHisto1D[1]->Fill(mRefMult);
  mRefMultHist[CentId]->Fill(mRefMult);
  mHisto2D[6]->Fill(PV.X(),PV.Y());   // Vy versus Vx

  if (fabs(PV.Z())>30.0) return 0;
  if (sqrt(pow(PV.X(),2)+pow(PV.Y(),2))>1.0) return 0;
  mHisto2D[7]->Fill(event->refMult(),              event->numberOfGlobalTracks());
  mHisto2D[8]->Fill(event->refMult(),              event->btofTrayMultiplicity());
  mHisto2D[9]->Fill(event->numberOfGlobalTracks(), event->btofTrayMultiplicity());

  StEpdEpInfo result = mEpFinder->Results(mEpdHits,PV,CentId);  // and now you have all the EP info you could ever want :-)

  for (int order=1; order<4; order++){
    mEpdAveCos[order-1][0]->Fill(double(CentId),cos((double)order*(result.EastRawPsi(order)-result.WestRawPsi(order))));
    mEpdAveCos[order-1][1]->Fill(double(CentId),cos((double)order*(result.EastPhiWeightedAndShiftedPsi(order)-result.WestPhiWeightedAndShiftedPsi(order))));
    mEpdEwPsi[order-1][0]->Fill(result.EastRawPsi(order),result.WestRawPsi(order));
    mEpdEwPsi[order-1][1]->Fill(result.EastPhiWeightedAndShiftedPsi(order),result.WestPhiWeightedAndShiftedPsi(order));
    if ((CentId<=6)&&(CentId>=4)){
      mEpdEwPsi_midCentral[order-1][0]->Fill(result.EastRawPsi(order),result.WestRawPsi(order));
      mEpdEwPsi_midCentral[order-1][1]->Fill(result.EastPhiWeightedAndShiftedPsi(order),result.WestPhiWeightedAndShiftedPsi(order));
    }
  }

  return 0;
}

//=================================================
short PicoAnalyzer::Finish(){

  cout << "In PicoAnalyzer::Finish - calling StEpdEpFinder::Finish()\n";

  mEpFinder->Finish();

  cout << "I have called it\n";

  double tempError[1000];
  for (int i=0; i<1000; i++){
    tempError[i]=0.0;
  }

  for (int icent=0; icent<9; icent++){
    for (int ewFull=0; ewFull<3; ewFull++){
      for (int bin=1; bin<=mEpdWeight[icent][ewFull]->GetXaxis()->GetNbins(); bin++){
	mEpdWeight[icent][ewFull]->SetError(tempError);
      }
      for (int order=1; order<4; order++){
	mEpdRawVn[icent][ewFull][order-1]->Divide(mEpdWeight[icent][ewFull]);
      }
      delete mEpdWeight[icent][ewFull];
    }
  }

  mHistoFile->Write();
  mHistoFile->Close();


  /*
  mCorrectionInputFile->Close();

  for (int ew=0; ew<2; ew++){
    mPhiWeightOutput[ew]->Divide(mPhiAveraged[ew]);
    delete mPhiAveraged[ew];
  }

  mCorrectionOutputFile->Write();
  mCorrectionOutputFile->Close();

  if (mQ1NtupleFile){
    mQ1NtupleFile->Write();
    mQ1NtupleFile->Close();
  }
  if (mQ2NtupleFile){
    mQ2NtupleFile->Write();
    mQ2NtupleFile->Close();
  }

  //  mTpcWeightFile->Write();
  */
  cout << "Finish!!\n\n";

  return 0;
}

//-----------------------------------
// I copied this directly from Isaac
// it's for the 2018 isobars
int PicoAnalyzer::FindCent(int Multiplicity){
  int CentId = -1;
  if(Multiplicity <= 19) CentId = -1; // > 80%
  else if (Multiplicity > 19 && Multiplicity <= 31) CentId = 0; // 70-80%
  else if (Multiplicity > 31 && Multiplicity <= 46) CentId = 1; // 60-70%
  else if (Multiplicity > 46 && Multiplicity <= 67) CentId = 2; // 50-60%
  else if (Multiplicity > 67 && Multiplicity <= 93) CentId = 3; // 40-50%
  else if (Multiplicity > 93 && Multiplicity <= 128) CentId = 4; // 30-40%
  else if (Multiplicity > 128 && Multiplicity <= 172) CentId = 5; // 20-30%
  else if (Multiplicity > 172 && Multiplicity <= 230) CentId = 6; // 10-20%
  else if (Multiplicity > 230 && Multiplicity <= 267) CentId = 7; // 5-10%
  else if (Multiplicity > 267) CentId = 8; // 0-5%
  return CentId;
}



/*
    To make the file that is read in by the method below, there are two ways:

    First, from Dmitry:
    How-To: invoke RunLog, select appropriate filters ("physics" - check, "tpx"
    - check, "filter bad runs": check three times etc), click "select". Then click
    on the number of runs displayed below the select button: "Runs: 3513" (click on
    number). In the middle pane you will get an ascii list like this:

    ...
    19125030 1525539034 Ru 100.00 Ru 100.00
    19125029 1525537184 Ru 100.00 Ru 100.00
    19125028 1525535327 Ru 100.00 Ru 100.00
    19125027 1525534083 Ru 100.00 Ru 100.00
    19125022 1525526621 Zr 100.00 Zr 100.00
    19125021 1525524762 Zr 100.00 Zr 100.00
    19125020 1525522902 Zr 100.00 Zr 100.00
    19125019 1525520656 Zr 100.00 Zr 100.00
    19125018 1525518764 Zr 100.00 Zr 100.00
    ...
    columns: run number, run start time (unixtime), blue beam species, blue energy, yellow beam species, yellow energy


    Next, from Gene, though I didn't try it:
    On RCF:
    mysql -h heston.star.bnl.gov -P 3501 -C RunLog -e "select runNumber,blueSpecies from beamInfo"
    Output that to a file and modify as you please.
*/



void PicoAnalyzer::ReadInSystems(){
  // 2018 run ("run 19") colliding system:  1=Au27Au / 2=Ru200Ru / 3=Zr200Zr / 4-Au27fixedTarget
  for (int day=0; day<365; day++){
    for (int runInDay=0; runInDay<500; runInDay++){
      mCollidingSystem[day][runInDay] = 0;
    }
  }
  int runId,junk,iSys;
  int iReadIn(0);
  ifstream RunSystems("Run2018systems.txt");
  if (!RunSystems){
    cout << "Error opening systems file!\n";
    return;
  }
  do{
    RunSystems >> runId >> junk >> iSys;
    if (RunSystems.eof()) break;
    if (runId!=0){
      iReadIn++;
      mCollidingSystem[(runId%1000000)/1000][runId%1000] = iSys;
    }
    else{
      RunSystems.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
  }while(kTRUE);

  cout << "I have learned the colliding system for " << iReadIn << " runs\n";
  RunSystems.close();
}

//------------------------------------------------------
short PicoAnalyzer::WhichSystem(int runId){
  if (runId/1000000 != 19){
    cout << "Sorry, I only know about 2018 run.  Do some coding!\n";
    return 0;
  }
  short shColSys = mCollidingSystem[(runId%1000000)/1000][runId%1000];
  return shColSys;
}

void PicoAnalyzer::FillPhiWeightHistos(StPicoEpdHit* epdHit, double weight){
  int ew = (epdHit->id()<0)?0:1;

  /*
     mPhiWeightOutput are the histograms that will be used for Phi-Weighting in the next run
     (mPhiAveraged is just used for normalization of mPhiWeightOutput at the end, then it is deleted)
     These histograms are binned in PP and TT, so GetBinContent(10,4) gives the weight for PP10TT04.
  */
  mPhiWeightOutput[ew]->Fill(epdHit->position(),epdHit->tile(),weight);
  if (epdHit->tile()==1){
    for (int pp=1; pp<13; pp++) mPhiAveraged[ew]->Fill(pp,1,weight/12.0);
  }
  else{
    for (int pp=1; pp<13; pp++){
      for (int tt=2*(epdHit->row()-1); tt<2*epdHit->row(); tt++) mPhiAveraged[ew]->Fill(pp,tt,weight/24.0);
    }
  }

  /*
    The histograms below are the same as the ones above, except that they are binned
    in row and phi.  This makes it really nice to make plots in polar mode.

    Unfortunately, unless you want rotatable lego/surface plots (which probably you don't),
    the way you have to draw those "polar" histograms is as follows:
    gPad->DrawFrame(-17,-17,17,17);
    EpdWeightWestPolar->Draw("colz,pol,same");
  */
  double phi = mEpdGeom->TileCenter(epdHit->id()).Phi();
  double phiDraw = (phi>0.0)?phi:phi+2.0*TMath::Pi();
  double dPhi = 2.0*TMath::Pi()/24.0;
  int row = epdHit->row();
  if (row==1){      // special case, tile 1.  Display it as two "tiles"
    mHisto2D[12+ew]->Fill(phiDraw+dPhi/2.0,row,weight);
    mHisto2D[12+ew]->Fill(phiDraw-dPhi/2.0,row,weight);
    for (double phiTemp=dPhi/2.0; phiTemp<2.0*TMath::Pi(); phiTemp+=dPhi)
      mHisto2D[14+ew]->Fill(phiTemp,row,weight/12.0);
  }
  else{
    mHisto2D[12+ew]->Fill(phiDraw,row,weight);
    for (double phiTemp=dPhi/2.0; phiTemp<2.0*TMath::Pi(); phiTemp+=dPhi)
      mHisto2D[14+ew]->Fill(phiTemp,row,weight/24.0);
  }

}

double PicoAnalyzer::v1Weight(int CentId, double eta){
  // https://drupal.star.bnl.gov/STAR/blog/lisa/optimizing-ep1-resolution-au27au-ring-dependent-weights-eta-dependent-weights
  double lin[9] = {-1.950, -1.900, -1.850, -1.706, -1.438, -1.340, -1.045, -0.717, -0.700};
  double cub[9] = {0.1608, 0.1600, 0.1600, 0.1595, 0.1457, 0.1369, 0.1092, 0.0772, 0.0700};

  double etaMagnitude = fabs(eta);   // we have the same weights for both wheels.  Any sign flips (East versus West etc) need to be done in the main code.

  return lin[CentId]*etaMagnitude + cub[CentId]*pow(etaMagnitude,3);
}
