#include "StProManger.h"
#include "StConstants.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TString.h"
#include "TMath.h"
#include "TH1F.h"

ClassImp(StProManger)

    TString StProManger::mVStr[2] = {"pos","neg"};

//---------------------------------------------------------------------------------

StProManger::StProManger()
{
}

//---------------------------------------------------------------------------------

StProManger::~StProManger()
{
}

//---------------------------------------------------------------------------------

void StProManger::InitReCenter()
{

    for(Int_t i = 0; i < 2; i++) // vertex pos/neg
    {
        for(Int_t j = 0; j < 4; j++) // eta_gap
        {
            TString ProName;
            // Event Plane method
            ProName = Form("qx_2nd_Vertex_%s_EtaGap_%d_East_EP",mVStr[i].Data(),j);
            p_mq2x_East_EP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// neg eta || x axis is RunIndex, y axis is Centrality
            ProName = Form("qy_2nd_Vertex_%s_EtaGap_%d_East_EP",mVStr[i].Data(),j);
            p_mq2y_East_EP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// neg eta || x axis is RunIndex, y axis is Centrality
            ProName = Form("qx_2nd_Vertex_%s_EtaGap_%d_West_EP",mVStr[i].Data(),j);
            p_mq2x_West_EP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// pos eta || x axis is RunIndex, y axis is Centrality
            ProName = Form("qy_2nd_Vertex_%s_EtaGap_%d_West_EP",mVStr[i].Data(),j);
            p_mq2y_West_EP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// pos eta || x axis is RunIndex, y axis is Centrality

            ProName = Form("qx_3rd_Vertex_%s_EtaGap_%d_East_EP",mVStr[i].Data(),j);
            p_mq3x_East_EP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// neg eta || x axis is RunIndex, y axis is Centrality
            ProName = Form("qy_3rd_Vertex_%s_EtaGap_%d_East_EP",mVStr[i].Data(),j);
            p_mq3y_East_EP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// neg eta || x axis is RunIndex, y axis is Centrality
            ProName = Form("qx_3rd_Vertex_%s_EtaGap_%d_West_EP",mVStr[i].Data(),j);
            p_mq3x_West_EP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// pos eta || x axis is RunIndex, y axis is Centrality
            ProName = Form("qy_3rd_Vertex_%s_EtaGap_%d_West_EP",mVStr[i].Data(),j);
            p_mq3y_West_EP[i][j] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// pos eta || x axis is RunIndex, y axis is Centrality

        }

        TString ProName_Full;
        // Event Plane method
        ProName_Full = Form("qx_2nd_Vertex_%s_Full_EP",mVStr[i].Data());
        p_mq2x_Full_EP[i] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// x axis is RunIndex, y axis is Centrality
        ProName_Full = Form("qy_2nd_Vertex_%s_Full_EP",mVStr[i].Data());
        p_mq2y_Full_EP[i] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// x axis is RunIndex, y axis is Centrality
        ProName_Full = Form("qx_3rd_Vertex_%s_Full_EP",mVStr[i].Data());
        p_mq3x_Full_EP[i] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// x axis is RunIndex, y axis is Centrality
        ProName_Full = Form("qy_3rd_Vertex_%s_Full_EP",mVStr[i].Data());
        p_mq3y_Full_EP[i] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5);// x axis is RunIndex, y axis is Centrality
    }
}

//----------------------------------------------------------------------------

//----------------------------------------------------------------------------

void StProManger::FillTrackEast(TVector2 q2Vector, TVector2 q3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Float_t pt) // i = vertex pos/neg, j = eta_gap
{
    const Float_t q2x = q2Vector.X();
    const Float_t q2y = q2Vector.Y();
    const Float_t q3x = q3Vector.X();
    const Float_t q3y = q3Vector.Y();

    Float_t w;
    if(pt <= TriFlow::mPrimPtWeight)
    {
        w = pt;
    }
    if(pt > TriFlow::mPrimPtWeight)
    {
        w = TriFlow::mPrimPtWeight;
    }

    // Event Plane method
    p_mq2x_East_EP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2x,(Double_t)w);
    p_mq2y_East_EP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2y,(Double_t)w);

    p_mq3x_East_EP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q3x,(Double_t)w);
    p_mq3y_East_EP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q3y,(Double_t)w);


}

void StProManger::FillTrackWest(TVector2 q2Vector, TVector2 q3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Float_t pt) // i = vertex pos/neg, j = eta_gap
{
    const Float_t q2x = q2Vector.X();
    const Float_t q2y = q2Vector.Y();
    const Float_t q3x = q3Vector.X();
    const Float_t q3y = q3Vector.Y();

    Float_t w;
    if(pt <= TriFlow::mPrimPtWeight)
    {
        w = pt;
    }
    if(pt > TriFlow::mPrimPtWeight)
    {
        w = TriFlow::mPrimPtWeight;
    }

    // Event Plane method
    p_mq2x_West_EP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2x,(Double_t)w);
    p_mq2y_West_EP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2y,(Double_t)w);

    p_mq3x_West_EP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q3x,(Double_t)w);
    p_mq3y_West_EP[i][j]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q3y,(Double_t)w);


}

void StProManger::FillTrackFull(TVector2 q2Vector, TVector2 q3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Float_t pt) // i = vertex pos/neg
{
    const Float_t q2x = q2Vector.X();
    const Float_t q2y = q2Vector.Y();
    const Float_t q3x = q3Vector.X();
    const Float_t q3y = q3Vector.Y();

    Float_t w;
    if(pt <= TriFlow::mPrimPtWeight)
    {
        w = pt;
    }
    if(pt > TriFlow::mPrimPtWeight)
    {
        w = TriFlow::mPrimPtWeight;
    }
    //cout << "weight = " <<w<< endl;
    //cout << "q2x = " << q2x << endl;
    //cout << "q2y = " << q2y << endl;

    // Event Plane method
    p_mq2x_Full_EP[i]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2x,(Double_t)w);
    p_mq2y_Full_EP[i]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q2y,(Double_t)w);

    p_mq3x_Full_EP[i]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q3x,(Double_t)w);
    p_mq3y_Full_EP[i]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)q3y,(Double_t)w);


}
//----------------------------------------------------------------------------

void StProManger::WriteReCenter()
{
    for(Int_t i = 0; i < 2; i++) // vertex pos/neg
    {
        for(Int_t j = 0; j < 4; j++) // eta_gap
        {
            // Event Plane method
            p_mq2x_East_EP[i][j]->Write();
            p_mq2y_East_EP[i][j]->Write();
            p_mq2x_West_EP[i][j]->Write();
            p_mq2y_West_EP[i][j]->Write();

            p_mq3x_East_EP[i][j]->Write();
            p_mq3y_East_EP[i][j]->Write();
            p_mq3x_West_EP[i][j]->Write();
            p_mq3y_West_EP[i][j]->Write();


        }
        // Event Plane method
        p_mq2x_Full_EP[i]->Write();
        p_mq2y_Full_EP[i]->Write();
        p_mq3x_Full_EP[i]->Write();
        p_mq3y_Full_EP[i]->Write();

    }
}

//----------------------------------------------------------------------------

