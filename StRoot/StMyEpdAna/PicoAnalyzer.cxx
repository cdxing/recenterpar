#include "PicoAnalyzer.h"
#include "StPicoEpdHit.h"
#include "StPicoBbcHit.h"
#include "../StEpdUtil/StEpdGeom.h"
#include "../StEpdUtil/StBbcGeom.h"
//#include "StPicoEvent.h"

#include "TMath.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TRandom3.h"

#include <iostream>
using namespace std;

ClassImp(PicoAnalyzer)

double PicoAnalyzer::GetBbcPmtPhi(short PmtId){
  unsigned short nBbcTiles;        // how many (1 or 2) tiles are associated with a given BBC phototube
  unsigned short tileNumbers[2];   // what are the tile ids of tiles associated with a given BBC phototube
  int chooseOne=0;
  // here we go with the bloody nightmare of the BBC's "shared phototubes"....
  mBbcGeom->GetTilesOfPmt(abs(PmtId),&nBbcTiles,tileNumbers);
  if (nBbcTiles>1) chooseOne = (mRan->Rndm()<0.5)?0:1;
  double phi = mBbcGeom->TileCenter(tileNumbers[chooseOne]*PmtId/abs(PmtId)).Phi();
  return phi;
}




//=================================================
PicoAnalyzer::PicoAnalyzer() : nEvents(0), mNmipQtB(115), mNmipQtC(160), mnMipThreshold(0.3), mPicoDst(0), mEpdBranch(0), mBbcBranch(0), mEventBranch(0), mEpdHits(0), mBbcHits(0)
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
void PicoAnalyzer::SetPicoDst(TTree* PicoDst){
  mPicoDst        = PicoDst;
  mEpdBranch      = mPicoDst->GetBranch("EpdHit");
  mBbcBranch      = mPicoDst->GetBranch("BbcHit");
  mEventBranch    = mPicoDst->GetBranch("Event");

  mRefMultNegLeaf = mPicoDst->GetLeaf("Event.mRefMultNeg");
  mRefMultPosLeaf = mPicoDst->GetLeaf("Event.mRefMultPos");
  mVzLeaf         = mPicoDst->GetLeaf("Event.mPrimaryVertex.mX3");

  mPicoDst->SetBranchStatus("*",0);         // turns OFF all branches  (speeds it up :-)
  mPicoDst->SetBranchStatus("EpdHit*",1);   // note you need the asterisk
  mPicoDst->SetBranchStatus("BbcHit*",1);
  mPicoDst->SetBranchStatus("Event*",1);

  // this below MIGHT not work to do only once.  might have to do each time (i.e. in make)
  //  Well, it seems to be okay to do it here.
  if (mEpdHits!=0) delete mEpdHits;
  if (mBbcHits!=0) delete mBbcHits;
  mEpdHits = new TClonesArray("StPicoEpdHit");
  mBbcHits = new TClonesArray("StPicoBbcHit");
  mEpdBranch->SetAddress(&mEpdHits);
  mBbcBranch->SetAddress(&mBbcHits);

  cout << "I have set all branches, leaves, etc \n";
}

//=================================================
void PicoAnalyzer::ReadLeaves(){
  mRefMult = mRefMultNegLeaf->GetValue()+mRefMultPosLeaf->GetValue();
  mVz      = mVzLeaf->GetValue();
}

//=================================================
short PicoAnalyzer::Init(){
  mHistoFile = new TFile("EpdHists.root","RECREATE");

  for (int ew=0; ew<2; ew++){
    for (int pp=1; pp<13; pp++){
      for (int tt=1; tt<32; tt++){
	mNmipDists[ew][pp-1][tt-1] = new TH1D(Form("NmipEW%dPP%dTT%d",ew,pp,tt),Form("NmipEW%dPP%dTT%d",ew,pp,tt),200,0,40);
      }
    }
  }

  // One-dimensional histograms
  mHisto1D[0] = new TH1D("Vz","Vz",100,-80,80);
  mHisto1D[1] = new TH1D("RefMult","RefMult",100,-10,600);
  mHisto1D[2] = new TH1D("dNdeta","dNdeta",100,-5.5,5.5);
  mHisto1D[3] = new TH1D("Cos2DeltaPsi2","Cos2DeltaPsi2",1000,-1,1);
  mHisto1D[4] = new TH1D("Cos1DeltaPsi1","Cos1DeltaPsi1",1000,-1,1);
  mHisto1D[5] = new TH1D("BbcCos1DeltaPsi1","BbcCos1DeltaPsi1",1000,-1,1);
  mHisto1D[6] = new TH1D("BbcCos2DeltaPsi2","BbcCos2DeltaPsi2",1000,-1,1);

  // Two-dimensional histograms
  mHisto2D[0] = new TH2D("Phi2WestVsPhiEast","Phi2WestVsPhiEast",40,-TMath::Pi()/2.0,TMath::Pi()/2.0,40,-TMath::Pi()/2.0,TMath::Pi()/2.0);
  mHisto2D[0]->GetXaxis()->SetTitle("#phi_{2,East}");
  mHisto2D[0]->GetYaxis()->SetTitle("#phi_{2,West}");
  mHisto2D[1] = new TH2D("Phi1WestVsPhiEast","Phi1WestVsPhiEast",40,-TMath::Pi(),TMath::Pi(),40,-TMath::Pi(),TMath::Pi());
  mHisto2D[1]->GetXaxis()->SetTitle("#phi_{1,East}");   
  mHisto2D[1]->GetYaxis()->SetTitle("#phi_{1,West}");
  mHisto2D[2] = new TH2D("RefMultVsEpdMipSum","RefMultVsEpdMipSum",100,0,5000,100,0,500);
  mHisto2D[3] = new TH2D("BbcAdcSumVsEpdMipSum","BbcAdcSumVsEpdMipSum",100,0,5000,100,0,50000);
  mHisto2D[4] = new TH2D("BbcPhi1WestVsPhiEast","BbcPhi1WestVsPhiEast",40,-TMath::Pi(),TMath::Pi(),40,-TMath::Pi(),TMath::Pi());
  mHisto2D[4]->GetXaxis()->SetTitle("#phi_{1,East}");
  mHisto2D[4]->GetYaxis()->SetTitle("#phi_{1,West}");
  mHisto2D[5] = new TH2D("BbcPhi2WestVsPhiEast","BbcPhi2WestVsPhiEast",40,-TMath::Pi()/2.0,TMath::Pi()/2.0,40,-TMath::Pi()/2.0,TMath::Pi()/2.0);
  mHisto2D[5]->GetXaxis()->SetTitle("#phi_{2,East}");
  mHisto2D[5]->GetYaxis()->SetTitle("#phi_{2,West}");

  //----------------- Ntuples -----------------
  TString VariableList = "";
  for (int ew=0; ew<2; ew++){
    TString ewstring = (ew==0)?"ER":"WR";
    for (int ring=1; ring<17; ring++){
      VariableList += Form("%s%dQx:%s%dQy",ewstring.Data(),ring,ewstring.Data(),ring);
      if (!((ew==1)&&(ring==16))) VariableList += ":";
    }
  }
  mQ1NtupleFile = new TFile("Q1ntuple.root","RECREATE");
  mQ1vectorNtuple = new TNtuple("RingQ1vectors","RingQ1vectors",VariableList.Data());
  mQ2NtupleFile = new TFile("Q2ntuple.root","RECREATE");
  mQ2vectorNtuple = new TNtuple("RingQ2vectors","RingQ2vectors",VariableList.Data());


  

  cout << "Init done\n";
  return 0;
}

  
//=================================================
short PicoAnalyzer::Make(int iEvent){
  //----------------- get data --------------
  mPicoDst->GetEntry(iEvent);
  ReadLeaves();
  //----- done getting data; have fun! ------

  TVector3 PV(0.0,0.0,mVz);

  mHisto1D[0]->Fill(mVz);
  mHisto1D[1]->Fill(mRefMult);

  StPicoEpdHit* epdHit;
  double Q[2][2][2] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};  // indices are east/west, order (1 or 2), component (x or y)
  double nMipSum(0.0);

  float Q1[64];   // Qx and Qy from all 16 rins on both sides.  Index = 32*EW + 2*(ring-1) + xy.  Where EW=0/1 for East/West.  Ring=1..16.  xy=0/1 for x/y
  float Q2[64];   // Qx and Qy from all 16 rins on both sides.  Index = 32*EW + 2*(ring-1) + xy.  Where EW=0/1 for East/West.  Ring=1..16.  xy=0/1 for x/y
  for (int iii=0; iii<64; iii++){Q1[iii]=0.0;  Q2[iii]=0.0;}

  //--------------- Begin loop over EPD hits --------------
  for (int hit=0; hit<mEpdHits->GetEntries(); hit++){
    epdHit = (StPicoEpdHit*)((*mEpdHits)[hit]);
    double nMip = (epdHit->tile()<10)?epdHit->adc()/mNmipQtC:epdHit->adc()/mNmipQtB;
    if (nMip<mnMipThreshold) continue;

    TVector3 end = mEpdGeom->RandomPointOnTile(epdHit->id());
    TVector3 StraightLine = end-PV;
    double eta = StraightLine.Eta();
    mHisto1D[2]->Fill(eta,nMip);
    nMipSum += nMip;
    int ew = (epdHit->id()<0)?0:1;
    mNmipDists[ew][epdHit->position()-1][epdHit->tile()-1]->Fill(nMip);
    double maxWeight = 3;
    double weight = (nMip>maxWeight)?maxWeight:nMip;
    double phi = mEpdGeom->TileCenter(epdHit->id()).Phi();
    
    int Qindex = 32*ew + 2*(epdHit->row()-1);
    double EWsign = (ew==0)?1.0:-1.0;
    Q1[Qindex]   += EWsign*weight*cos(phi);    // for first-order plane, east and west will point "in opposite directions"
    Q1[Qindex+1] += EWsign*weight*sin(phi);
    Q2[Qindex]   += weight*cos(2.0*phi);
    Q2[Qindex+1] += weight*sin(2.0*phi);



    if ((abs(eta)>0)&&(abs(eta)<4000000.0)){
      for (int order=1; order<3; order++){
	Q[ew][order-1][0] += weight*cos((double)order*phi);
	Q[ew][order-1][1] += weight*sin((double)order*phi);
      }
    }
  }                       // end loop over EpdHits
  //--------------- End loop over EPD hits --------------


  double PsiEP[2][2];     // indices are east/west and order (1 or 2)

//   for (int ew=0; ew<2; ew++){
//     for (int order=1; order<3; order++){
//       PsiEP[ew][order-1] = atan2(Q[ew][order-1][1],Q[ew][order-1][0])/(double)order;
//     }
//   }
//   PsiEP[0][0] += TMath::Pi();  // for first-order plane, East should "point opposite" West...


  double RowWeight[16]={3.8140,0.5399,-0.5231,0.3281,0.8812,1.0733,0.9113,0.9067,1.0962,1.1436,0.9466,0.8508,0.9539,1.0178,0.9952,1.00000000};

  double Q1Xsum[2]={0.0,0.0};  // index = ew
  double Q1Ysum[2]={0.0,0.0};
  double Q2Xsum[2]={0.0,0.0};
  double Q2Ysum[2]={0.0,0.0};
  for (int ew=0; ew<2; ew++){
    for (int jrow=0; jrow<16; jrow++){
      Q1Xsum[ew] += RowWeight[jrow] * Q1[32*ew+2*jrow];
      Q1Ysum[ew] += RowWeight[jrow] * Q1[32*ew+2*jrow+1];
      Q2Xsum[ew] += RowWeight[jrow] * Q2[32*ew+2*jrow];
      Q2Ysum[ew] += RowWeight[jrow] * Q2[32*ew+2*jrow+1];
    }
    PsiEP[ew][0] = atan2(Q1Ysum[ew],Q1Xsum[ew]);
    if (PsiEP[ew][0]<-TMath::Pi()) PsiEP[ew][0] += 2.0*TMath::Pi();
    if (PsiEP[ew][0]>TMath::Pi())  PsiEP[ew][0] -= 2.0*TMath::Pi();
    PsiEP[ew][1] = 0.5*atan2(Q2Ysum[ew],Q2Xsum[ew]);
    if (PsiEP[ew][1]<-TMath::Pi()/2.0) PsiEP[ew][0] += TMath::Pi();
    if (PsiEP[ew][1]>TMath::Pi()/2.0)  PsiEP[ew][0] -= TMath::Pi();
  }




  if ((mRefMult>93)&&(mRefMult<172)){    // according to Isaac, this should be 10%-40%, our "best" centrality.
    mHisto1D[3]->Fill(cos(2.0*(PsiEP[1][1]-PsiEP[0][1])));
    mHisto1D[4]->Fill(cos(1.0*(PsiEP[1][0]-PsiEP[0][0])));
    mHisto2D[0]->Fill(PsiEP[0][1],PsiEP[1][1]);
    mHisto2D[1]->Fill(PsiEP[0][0],PsiEP[1][0]);
    mQ1vectorNtuple->Fill(Q1);
    mQ2vectorNtuple->Fill(Q2);
  }
  mHisto2D[2]->Fill(nMipSum,mRefMult);

  // now look also at BBC stuff
  // 1dhisto[5] = cos2deltapsi2
  // 1dhisto[4] = phi2West vs East

  for (int ew=0; ew<2; ew++){
    for (int order=1; order<3; order++){
      for (int xy=0; xy<2; xy++){
	Q[ew][order-1][xy] = 0.0;
      }
    }
  }

  StPicoBbcHit* bbcHit;
  double BbcAdcSum(0.0);
  
  for (int hit=0; hit<mBbcHits->GetEntries(); hit++){
    bbcHit = (StPicoBbcHit*)((*mBbcHits)[hit]);
    double adc = bbcHit->adc();
    BbcAdcSum += adc;
    int ew = (bbcHit->side()>0)?1:0;

    double phi = GetBbcPmtPhi(bbcHit->id());  // this will take care of the nightmare of BBC's tile-pmt mapping...

    Q[ew][0][0] += adc*cos(phi);
    Q[ew][0][1] += adc*sin(phi);
    Q[ew][1][0] += adc*cos(2*phi);
    Q[ew][1][1] += adc*sin(2*phi);
  }
  mHisto2D[3]->Fill(nMipSum,BbcAdcSum);
  for (int ew=0; ew<2; ew++){
    for (int order=1; order<3; order++){
      PsiEP[ew][order-1] = atan2(Q[ew][order-1][1],Q[ew][order-1][0])/(double)order;
    }
  }
  if ((mRefMult>100)&&(mRefMult<240)){
    mHisto1D[5]->Fill(cos(1.0*(PsiEP[1][0]-PsiEP[0][0])));
    mHisto1D[6]->Fill(cos(2.0*(PsiEP[1][1]-PsiEP[0][1])));

    mHisto2D[4]->Fill(PsiEP[0][0],PsiEP[1][0]);
    mHisto2D[5]->Fill(PsiEP[0][1],PsiEP[1][1]);
  }
  
  return 0;
}

//=================================================
short PicoAnalyzer::Finish(){
  mHistoFile->Write();
  mHistoFile->Close();

  mQ1NtupleFile->Write();
  mQ1NtupleFile->Close();
  mQ2NtupleFile->Write();
  mQ2NtupleFile->Close();

  cout << "Finish!!\n\n";
  /* save histograms to a file if you want */
  return 0;
}


