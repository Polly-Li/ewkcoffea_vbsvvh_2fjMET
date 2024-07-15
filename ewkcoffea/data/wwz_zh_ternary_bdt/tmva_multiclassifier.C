#include <iostream>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TSystem.h"
#include <vector>

//~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
//
//
// BDT Opposite Flavor V7 (Multi-Classifier)
//
//
//~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*

namespace ewkcoffea
{

    class BDT_OF_v7
    {
        public:
            TMVA::Reader *reader;

            float m_ll;                // Label="m_ll" Title="m_{ll}" Unit="GeV" Internal="m_ll" Type="F" Min="1.20001631e+01" Max="2.46282373e+03"/>
            float m_4l;                // Label="m_4l" Title="m_{4l}" Unit="GeV" Internal="m_4l" Type="F" Min="1.00103073e+02" Max="4.44906885e+03"/>
            float dPhi_4Lep_MET;       // Label="dPhi_4Lep_MET" Title="#Delta#phi(4Lep,p_{T}^{miss})" Unit="" Internal="dPhi_4Lep_MET" Type="F" Min="3.11136246e-05" Max="3.14159226e+00"/>
            float dPhi_WW_MET;         // Label="dPhi_WW_MET" Title="#Delta#phi(Wcands,p_{T}^{miss})" Unit="" Internal="dPhi_WW_MET" Type="F" Min="8.94069672e-06" Max="3.14158940e+00"/>
            float dPhi_W1_MET;         // Label="dPhi_W1_MET" Title="#Delta#phi(l^{W1},p_{T}^{miss})" Unit="" Internal="dPhi_W1_MET" Type="F" Min="0.00000000e+00" Max="3.14158368e+00"/>
            float dPhi_W2_MET;         // Label="dPhi_W2_MET" Title="#Delta#phi(l^{W2},p_{T}^{miss})" Unit="" Internal="dPhi_W2_MET" Type="F" Min="0.00000000e+00" Max="3.14158368e+00"/>
            float dR_Wcands;           // Label="dR_Wcands" Title="#Delta R(l^{W1},l^{W2})" Unit="" Internal="dR_Wcands" Type="F" Min="4.60619330e-02" Max="5.51130724e+00"/>
            float dR_Zcands;           // Label="dR_Zcands" Title="#Delta R(l^{Z1},l^{Z2})" Unit="" Internal="dR_Zcands" Type="F" Min="1.02005213e-01" Max="4.86037207e+00"/>
            float dR_WW_Z;             // Label="dR_WW_Z" Title="#Delta R(Wcands,Zcands)" Unit="" Internal="dR_WW_Z" Type="F" Min="7.74417212e-03" Max="9.68943882e+00"/>
            float MET;                 // Label="MET" Title="p_{T}^{miss}" Unit="GeV" Internal="MET" Type="F" Min="1.14202060e-01" Max="1.62917688e+03"/>
            float MT2;                 // Label="MT2" Title="M_{T2}" Unit="GeV" Internal="MT2" Type="F" Min="1.05712891e-01" Max="3.86917603e+02"/>
            float Pt4l;                // Label="Pt4l" Title="p_{T}^{4l}" Unit="GeV" Internal="Pt4l" Type="F" Min="3.94652113e-02" Max="1.38914185e+03"/>
            float HT;                  // Label="HT" Title="H_{T}" Unit="GeV" Internal="HT" Type="F" Min="0.00000000e+00" Max="2.10278125e+03"/>
            float STLep;               // Label="STLep" Title="#Sigma_{lep,MET} p_{T}" Unit="GeV" Internal="STLep" Type="F" Min="7.89844208e+01" Max="4.27097070e+03"/>
            float leading_Zcand_pt;    // Label="leading_Zcand_pt" Title="p_{T}^{Z1}" Unit="GeV" Internal="leading_Zcand_pt" Type="F" Min="1.20194635e+01" Max="1.63574243e+03"/>
            float subleading_Zcand_pt; // Label="subleading_Zcand_pt" Title="p_{T}^{Z2}" Unit="GeV" Internal="subleading_Zcand_pt" Type="F" Min="1.00015154e+01" Max="7.25295532e+02"/>
            float leading_Wcand_pt;    // Label="leading_Wcand_pt" Title="p_{T}^{W1}" Unit="GeV" Internal="leading_Wcand_pt" Type="F" Min="1.01979656e+01" Max="1.64588867e+03"/>
            float subleading_Wcand_pt; // Label="subleading_Wcand_pt" Title="p_{T}^{W2}" Unit="GeV" Internal="subleading_Wcand_pt" Type="F" Min="1.00002556e+01" Max="1.17922314e+03"/>
            float njets;               // Label="njets" Title="N_{jets}" Unit="" Internal="njets" Type="F" Min="0.00000000e+00" Max="1.10000000e+01"/>
            float cos_helicity_X;      // Label="cos_helicity_X" Title="cos(#theta_{X})" Unit="" Internal="cos_helicity_X" Type="F" Min="6.58120086e-07" Max="9.99999106e-01"/>
            float MT_leading_Wcand;    // Label="MT_leading_Wcand" Title="M_{T}^{W1}" Unit="GeV" Internal="MT_leading_Wcand" Type="F" Min="0.00000000e+00" Max="2.66695630e+03"/>
            float MT_subleading_Wcand; // Label="MT_subleading_Wcand" Title="M_{T}^{W2}" Unit="GeV" Internal="MT_subleading_Wcand" Type="F" Min="0.00000000e+00" Max="1.36826062e+03"/>
            float MT_Wcands;           // Label="MT_Wcands" Title="M_{T}^{Wcands}" Unit="GeV" Internal="MT_Wcands" Type="F" Min="0.00000000e+00" Max="2.27764478e+03"/>
            float MT_4Lep;             // Label="MT_4Lep" Title="M_{T}^{4Lep}" Unit="GeV" Internal="MT_4Lep" Type="F" Min="0.00000000e+00" Max="2.89665723e+03"/>
            float min_dR_W1_jet;       // Label="min_dR_W1_jet" Title="min(#Delta R(l^{W1},j))" Unit="" Internal="min_dR_W1_jet" Type="F" Min="0.00000000e+00" Max="5.68600368e+00"/>
            float min_dR_W2_jet;       // Label="min_dR_W2_jet" Title="min(#Delta R(l^{W2},j))" Unit="" Internal="min_dR_W2_jet" Type="F" Min="0.00000000e+00" Max="5.60851860e+00"/>
            float dPhi_Zcand_MET;
            float leading_jet_pt;
            float subleading_jet_pt;
            float leading_jet_DeepFlav;

            BDT_OF_v7(TString xmlpath);
            ~BDT_OF_v7();
            std::vector<std::vector<float>> Eval(std::vector<float> m_ll_,
                                                 std::vector<float> m_4l_,
                                                 std::vector<float> dPhi_4Lep_MET_,
                                                 std::vector<float> dPhi_WW_MET_,
                                                 std::vector<float> dPhi_W1_MET_,
                                                 std::vector<float> dPhi_W2_MET_,
                                                 std::vector<float> dR_Wcands_,
                                                 std::vector<float> dR_Zcands_,
                                                 std::vector<float> dR_WW_Z_,
                                                 std::vector<float> MET_,
                                                 std::vector<float> MT2_,
                                                 std::vector<float> Pt4l_,
                                                 std::vector<float> HT_,
                                                 std::vector<float> STLep_,
                                                 std::vector<float> leading_Zcand_pt_,
                                                 std::vector<float> subleading_Zcand_pt_,
                                                 std::vector<float> leading_Wcand_pt_,
                                                 std::vector<float> subleading_Wcand_pt_,
                                                 std::vector<float> njets_,
                                                 std::vector<float> cos_helicity_X_,
                                                 std::vector<float> MT_leading_Wcand_,
                                                 std::vector<float> MT_subleading_Wcand_,
                                                 std::vector<float> MT_Wcands_,
                                                 std::vector<float> MT_4Lep_,
                                                 std::vector<float> min_dR_W1_jet_,
                                                 std::vector<float> min_dR_W2_jet_);

    };
}

ewkcoffea::BDT_OF_v7::BDT_OF_v7(TString xmlpath)
{
    reader = new TMVA::Reader("!Color:Silent");
    reader->AddVariable("m_ll", &m_ll);
    reader->AddVariable("m_4l", &m_4l);
    reader->AddVariable("dPhi_4Lep_MET", &dPhi_4Lep_MET);
    reader->AddVariable("dPhi_WW_MET", &dPhi_WW_MET);
    reader->AddVariable("dPhi_W1_MET", &dPhi_W1_MET);
    reader->AddVariable("dPhi_W2_MET", &dPhi_W2_MET);
    reader->AddVariable("dR_Wcands", &dR_Wcands);
    reader->AddVariable("dR_Zcands", &dR_Zcands);
    reader->AddVariable("dR_WW_Z", &dR_WW_Z);
    reader->AddVariable("MET", &MET);
    reader->AddVariable("MT2", &MT2);
    reader->AddVariable("Pt4l", &Pt4l);
    reader->AddVariable("HT", &HT);
    reader->AddVariable("STLep", &STLep);
    reader->AddVariable("leading_Zcand_pt", &leading_Zcand_pt);
    reader->AddVariable("subleading_Zcand_pt", &subleading_Zcand_pt);
    reader->AddVariable("leading_Wcand_pt", &leading_Wcand_pt);
    reader->AddVariable("subleading_Wcand_pt", &subleading_Wcand_pt);
    reader->AddVariable("njets", &njets);
    reader->AddVariable("cos_helicity_X", &cos_helicity_X);
    reader->AddVariable("MT_leading_Wcand", &MT_leading_Wcand);
    reader->AddVariable("MT_subleading_Wcand", &MT_subleading_Wcand);
    reader->AddVariable("MT_Wcands", &MT_Wcands);
    reader->AddVariable("MT_4Lep", &MT_4Lep);
    reader->AddVariable("min_dR_W1_jet", &min_dR_W1_jet);
    reader->AddVariable("min_dR_W2_jet", &min_dR_W2_jet);
    reader->AddSpectator("dPhi_Zcand_MET", &dPhi_Zcand_MET);
    reader->AddSpectator("leading_jet_pt", &leading_jet_pt);
    reader->AddSpectator("subleading_jet_pt", &subleading_jet_pt);
    reader->AddSpectator("leading_jet_DeepFlav", &leading_jet_DeepFlav);
    reader->BookMVA("BDT", xmlpath.Data());
}

ewkcoffea::BDT_OF_v7::~BDT_OF_v7()
{
    delete reader;
}

std::vector<std::vector<float>> ewkcoffea::BDT_OF_v7::Eval(std::vector<float> m_ll_,
                                                           std::vector<float> m_4l_,
                                                           std::vector<float> dPhi_4Lep_MET_,
                                                           std::vector<float> dPhi_WW_MET_,
                                                           std::vector<float> dPhi_W1_MET_,
                                                           std::vector<float> dPhi_W2_MET_,
                                                           std::vector<float> dR_Wcands_,
                                                           std::vector<float> dR_Zcands_,
                                                           std::vector<float> dR_WW_Z_,
                                                           std::vector<float> MET_,
                                                           std::vector<float> MT2_,
                                                           std::vector<float> Pt4l_,
                                                           std::vector<float> HT_,
                                                           std::vector<float> STLep_,
                                                           std::vector<float> leading_Zcand_pt_,
                                                           std::vector<float> subleading_Zcand_pt_,
                                                           std::vector<float> leading_Wcand_pt_,
                                                           std::vector<float> subleading_Wcand_pt_,
                                                           std::vector<float> njets_,
                                                           std::vector<float> cos_helicity_X_,
                                                           std::vector<float> MT_leading_Wcand_,
                                                           std::vector<float> MT_subleading_Wcand_,
                                                           std::vector<float> MT_Wcands_,
                                                           std::vector<float> MT_4Lep_,
                                                           std::vector<float> min_dR_W1_jet_,
                                                           std::vector<float> min_dR_W2_jet_)
{
    std::vector<std::vector<float>> rtn;
    for (unsigned int i = 0; i < m_ll_.size(); ++i)
    {
        m_ll = m_ll_[i];
        m_4l = m_4l_[i];
        dPhi_4Lep_MET = dPhi_4Lep_MET_[i];
        dPhi_WW_MET = dPhi_WW_MET_[i];
        dPhi_W1_MET = dPhi_W1_MET_[i];
        dPhi_W2_MET = dPhi_W2_MET_[i];
        dR_Wcands = dR_Wcands_[i];
        dR_Zcands = dR_Zcands_[i];
        dR_WW_Z = dR_WW_Z_[i];
        MET = MET_[i];
        MT2 = MT2_[i];
        Pt4l = Pt4l_[i];
        HT = HT_[i];
        STLep = STLep_[i];
        leading_Zcand_pt = leading_Zcand_pt_[i];
        subleading_Zcand_pt = subleading_Zcand_pt_[i];
        leading_Wcand_pt = leading_Wcand_pt_[i];
        subleading_Wcand_pt = subleading_Wcand_pt_[i];
        njets = njets_[i];
        cos_helicity_X = cos_helicity_X_[i];
        MT_leading_Wcand = MT_leading_Wcand_[i];
        MT_subleading_Wcand = MT_subleading_Wcand_[i];
        MT_Wcands = MT_Wcands_[i];
        MT_4Lep = MT_4Lep_[i];
        min_dR_W1_jet = min_dR_W1_jet_[i];
        min_dR_W2_jet = min_dR_W2_jet_[i];
        rtn.push_back(reader->EvaluateMulticlass("BDT"));
    }
    return rtn;
}



//~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
//
//
// BDT Same Flavor V7 (Multi-Classifier)
//
//
//~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*

namespace ewkcoffea
{

    class BDT_SF_v7
    {
        public:
            TMVA::Reader *reader;

            float m_ll;                // Label="m_ll" Title="m_{ll}" Unit="GeV" Internal="m_ll" Type="F" Min="1.20001631e+01" Max="2.46282373e+03"/>
            float m_4l;                // Label="m_4l" Title="m_{4l}" Unit="GeV" Internal="m_4l" Type="F" Min="1.00103073e+02" Max="4.44906885e+03"/>
            float dPhi_4Lep_MET;       // Label="dPhi_4Lep_MET" Title="#Delta#phi(4Lep,p_{T}^{miss})" Unit="" Internal="dPhi_4Lep_MET" Type="F" Min="3.11136246e-05" Max="3.14159226e+00"/>
            float dPhi_Zcand_MET;      // Label="dPhi_Zcand_MET" Title="#Delta#phi(4Lep,p_{T}^{miss})" Unit="" Internal="dPhi_Zcand_MET" Type="F" Min="4.30643559e-06" Max="3.14158511e+00"/>
            float dPhi_WW_MET;         // Label="dPhi_WW_MET" Title="#Delta#phi(Wcands,p_{T}^{miss})" Unit="" Internal="dPhi_WW_MET" Type="F" Min="8.94069672e-06" Max="3.14158940e+00"/>
            float dPhi_W1_MET;         // Label="dPhi_W1_MET" Title="#Delta#phi(l^{W1},p_{T}^{miss})" Unit="" Internal="dPhi_W1_MET" Type="F" Min="0.00000000e+00" Max="3.14158368e+00"/>
            float dPhi_W2_MET;         // Label="dPhi_W2_MET" Title="#Delta#phi(l^{W2},p_{T}^{miss})" Unit="" Internal="dPhi_W2_MET" Type="F" Min="0.00000000e+00" Max="3.14158368e+00"/>
            float dR_Wcands;           // Label="dR_Wcands" Title="#Delta R(l^{W1},l^{W2})" Unit="" Internal="dR_Wcands" Type="F" Min="4.60619330e-02" Max="5.51130724e+00"/>
            float dR_Zcands;           // Label="dR_Zcands" Title="#Delta R(l^{Z1},l^{Z2})" Unit="" Internal="dR_Zcands" Type="F" Min="1.02005213e-01" Max="4.86037207e+00"/>
            float dR_WW_Z;             // Label="dR_WW_Z" Title="#Delta R(Wcands,Zcands)" Unit="" Internal="dR_WW_Z" Type="F" Min="7.74417212e-03" Max="9.68943882e+00"/>
            float MET;                 // Label="MET" Title="p_{T}^{miss}" Unit="GeV" Internal="MET" Type="F" Min="1.14202060e-01" Max="1.62917688e+03"/>
            float MT2;                 // Label="MT2" Title="M_{T2}" Unit="GeV" Internal="MT2" Type="F" Min="1.05712891e-01" Max="3.86917603e+02"/>
            float Pt4l;                // Label="Pt4l" Title="p_{T}^{4l}" Unit="GeV" Internal="Pt4l" Type="F" Min="3.94652113e-02" Max="1.38914185e+03"/>
            float HT;                  // Label="HT" Title="H_{T}" Unit="GeV" Internal="HT" Type="F" Min="0.00000000e+00" Max="2.10278125e+03"/>
            float STLep;               // Label="STLep" Title="#Sigma_{lep,MET} p_{T}" Unit="GeV" Internal="STLep" Type="F" Min="7.89844208e+01" Max="4.27097070e+03"/>
            float leading_Zcand_pt;    // Label="leading_Zcand_pt" Title="p_{T}^{Z1}" Unit="GeV" Internal="leading_Zcand_pt" Type="F" Min="1.20194635e+01" Max="1.63574243e+03"/>
            float subleading_Zcand_pt; // Label="subleading_Zcand_pt" Title="p_{T}^{Z2}" Unit="GeV" Internal="subleading_Zcand_pt" Type="F" Min="1.00015154e+01" Max="7.25295532e+02"/>
            float leading_Wcand_pt;    // Label="leading_Wcand_pt" Title="p_{T}^{W1}" Unit="GeV" Internal="leading_Wcand_pt" Type="F" Min="1.01979656e+01" Max="1.64588867e+03"/>
            float subleading_Wcand_pt; // Label="subleading_Wcand_pt" Title="p_{T}^{W2}" Unit="GeV" Internal="subleading_Wcand_pt" Type="F" Min="1.00002556e+01" Max="1.17922314e+03"/>
            float njets;               // Label="njets" Title="N_{jets}" Unit="" Internal="njets" Type="F" Min="0.00000000e+00" Max="1.10000000e+01"/>
            float cos_helicity_X;      // Label="cos_helicity_X" Title="cos(#theta_{X})" Unit="" Internal="cos_helicity_X" Type="F" Min="6.58120086e-07" Max="9.99999106e-01"/>
            float MT_leading_Wcand;    // Label="MT_leading_Wcand" Title="M_{T}^{W1}" Unit="GeV" Internal="MT_leading_Wcand" Type="F" Min="0.00000000e+00" Max="2.66695630e+03"/>
            float MT_subleading_Wcand; // Label="MT_subleading_Wcand" Title="M_{T}^{W2}" Unit="GeV" Internal="MT_subleading_Wcand" Type="F" Min="0.00000000e+00" Max="1.36826062e+03"/>
            float MT_Wcands;           // Label="MT_Wcands" Title="M_{T}^{Wcands}" Unit="GeV" Internal="MT_Wcands" Type="F" Min="0.00000000e+00" Max="2.27764478e+03"/>
            float MT_4Lep;             // Label="MT_4Lep" Title="M_{T}^{4Lep}" Unit="GeV" Internal="MT_4Lep" Type="F" Min="0.00000000e+00" Max="2.89665723e+03"/>
            float min_dR_W1_jet;       // Label="min_dR_W1_jet" Title="min(#Delta R(l^{W1},j))" Unit="" Internal="min_dR_W1_jet" Type="F" Min="0.00000000e+00" Max="5.68600368e+00"/>
            float min_dR_W2_jet;       // Label="min_dR_W2_jet" Title="min(#Delta R(l^{W2},j))" Unit="" Internal="min_dR_W2_jet" Type="F" Min="0.00000000e+00" Max="5.60851860e+00"/>
            float leading_jet_pt;
            float subleading_jet_pt;
            float leading_jet_DeepFlav;

            BDT_SF_v7(TString xmlpath);
            ~BDT_SF_v7();
            std::vector<std::vector<float>> Eval(std::vector<float> m_ll_,
                                                 std::vector<float> m_4l_,
                                                 std::vector<float> dPhi_4Lep_MET_,
                                                 std::vector<float> dPhi_Zcand_MET_,
                                                 std::vector<float> dPhi_WW_MET_,
                                                 std::vector<float> dPhi_W1_MET_,
                                                 std::vector<float> dPhi_W2_MET_,
                                                 std::vector<float> dR_Wcands_,
                                                 std::vector<float> dR_Zcands_,
                                                 std::vector<float> dR_WW_Z_,
                                                 std::vector<float> MET_,
                                                 std::vector<float> MT2_,
                                                 std::vector<float> Pt4l_,
                                                 std::vector<float> HT_,
                                                 std::vector<float> STLep_,
                                                 std::vector<float> leading_Zcand_pt_,
                                                 std::vector<float> subleading_Zcand_pt_,
                                                 std::vector<float> leading_Wcand_pt_,
                                                 std::vector<float> subleading_Wcand_pt_,
                                                 std::vector<float> njets_,
                                                 std::vector<float> cos_helicity_X_,
                                                 std::vector<float> MT_leading_Wcand_,
                                                 std::vector<float> MT_subleading_Wcand_,
                                                 std::vector<float> MT_Wcands_,
                                                 std::vector<float> MT_4Lep_,
                                                 std::vector<float> min_dR_W1_jet_,
                                                 std::vector<float> min_dR_W2_jet_);

    };
}

ewkcoffea::BDT_SF_v7::BDT_SF_v7(TString xmlpath)
{
    reader = new TMVA::Reader("!Color:Silent");
    reader->AddVariable("m_ll", &m_ll);
    reader->AddVariable("m_4l", &m_4l);
    reader->AddVariable("dPhi_4Lep_MET", &dPhi_4Lep_MET);
    reader->AddVariable("dPhi_Zcand_MET", &dPhi_Zcand_MET);
    reader->AddVariable("dPhi_WW_MET", &dPhi_WW_MET);
    reader->AddVariable("dPhi_W1_MET", &dPhi_W1_MET);
    reader->AddVariable("dPhi_W2_MET", &dPhi_W2_MET);
    reader->AddVariable("dR_Wcands", &dR_Wcands);
    reader->AddVariable("dR_Zcands", &dR_Zcands);
    reader->AddVariable("dR_WW_Z", &dR_WW_Z);
    reader->AddVariable("MET", &MET);
    reader->AddVariable("MT2", &MT2);
    reader->AddVariable("Pt4l", &Pt4l);
    reader->AddVariable("HT", &HT);
    reader->AddVariable("STLep", &STLep);
    reader->AddVariable("leading_Zcand_pt", &leading_Zcand_pt);
    reader->AddVariable("subleading_Zcand_pt", &subleading_Zcand_pt);
    reader->AddVariable("leading_Wcand_pt", &leading_Wcand_pt);
    reader->AddVariable("subleading_Wcand_pt", &subleading_Wcand_pt);
    reader->AddVariable("njets", &njets);
    reader->AddVariable("cos_helicity_X", &cos_helicity_X);
    reader->AddVariable("MT_leading_Wcand", &MT_leading_Wcand);
    reader->AddVariable("MT_subleading_Wcand", &MT_subleading_Wcand);
    reader->AddVariable("MT_Wcands", &MT_Wcands);
    reader->AddVariable("MT_4Lep", &MT_4Lep);
    reader->AddVariable("min_dR_W1_jet", &min_dR_W1_jet);
    reader->AddVariable("min_dR_W2_jet", &min_dR_W2_jet);
    reader->AddSpectator("leading_jet_pt", &leading_jet_pt);
    reader->AddSpectator("subleading_jet_pt", &subleading_jet_pt);
    reader->AddSpectator("leading_jet_DeepFlav", &leading_jet_DeepFlav);
    reader->BookMVA("BDT", xmlpath.Data());
}

ewkcoffea::BDT_SF_v7::~BDT_SF_v7()
{
    delete reader;
}

std::vector<std::vector<float>> ewkcoffea::BDT_SF_v7::Eval(std::vector<float> m_ll_,
                                                           std::vector<float> m_4l_,
                                                           std::vector<float> dPhi_4Lep_MET_,
                                                           std::vector<float> dPhi_Zcand_MET_,
                                                           std::vector<float> dPhi_WW_MET_,
                                                           std::vector<float> dPhi_W1_MET_,
                                                           std::vector<float> dPhi_W2_MET_,
                                                           std::vector<float> dR_Wcands_,
                                                           std::vector<float> dR_Zcands_,
                                                           std::vector<float> dR_WW_Z_,
                                                           std::vector<float> MET_,
                                                           std::vector<float> MT2_,
                                                           std::vector<float> Pt4l_,
                                                           std::vector<float> HT_,
                                                           std::vector<float> STLep_,
                                                           std::vector<float> leading_Zcand_pt_,
                                                           std::vector<float> subleading_Zcand_pt_,
                                                           std::vector<float> leading_Wcand_pt_,
                                                           std::vector<float> subleading_Wcand_pt_,
                                                           std::vector<float> njets_,
                                                           std::vector<float> cos_helicity_X_,
                                                           std::vector<float> MT_leading_Wcand_,
                                                           std::vector<float> MT_subleading_Wcand_,
                                                           std::vector<float> MT_Wcands_,
                                                           std::vector<float> MT_4Lep_,
                                                           std::vector<float> min_dR_W1_jet_,
                                                           std::vector<float> min_dR_W2_jet_)
{
    std::vector<std::vector<float>> rtn;
    for (unsigned int i = 0; i < m_ll_.size(); ++i)
    {
        m_ll = m_ll_[i];
        m_4l = m_4l_[i];
        dPhi_4Lep_MET = dPhi_4Lep_MET_[i];
        dPhi_Zcand_MET = dPhi_Zcand_MET_[i];
        dPhi_WW_MET = dPhi_WW_MET_[i];
        dPhi_W1_MET = dPhi_W1_MET_[i];
        dPhi_W2_MET = dPhi_W2_MET_[i];
        dR_Wcands = dR_Wcands_[i];
        dR_Zcands = dR_Zcands_[i];
        dR_WW_Z = dR_WW_Z_[i];
        MET = MET_[i];
        MT2 = MT2_[i];
        Pt4l = Pt4l_[i];
        HT = HT_[i];
        STLep = STLep_[i];
        leading_Zcand_pt = leading_Zcand_pt_[i];
        subleading_Zcand_pt = subleading_Zcand_pt_[i];
        leading_Wcand_pt = leading_Wcand_pt_[i];
        subleading_Wcand_pt = subleading_Wcand_pt_[i];
        njets = njets_[i];
        cos_helicity_X = cos_helicity_X_[i];
        MT_leading_Wcand = MT_leading_Wcand_[i];
        MT_subleading_Wcand = MT_subleading_Wcand_[i];
        MT_Wcands = MT_Wcands_[i];
        MT_4Lep = MT_4Lep_[i];
        min_dR_W1_jet = min_dR_W1_jet_[i];
        min_dR_W2_jet = min_dR_W2_jet_[i];
        rtn.push_back(reader->EvaluateMulticlass("BDT"));
    }
    return rtn;
}
