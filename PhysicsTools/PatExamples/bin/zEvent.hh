#ifndef INC_13TEV_TW_ANALYSIS_ZEVENT_HH
#define INC_13TEV_TW_ANALYSIS_ZEVENT_HH

#include <vector>
#include <string>
#include <math.h>

#include <DataFormats/FWLite/interface/Event.h>
#include <CondFormats/BTauObjects/interface/BTagCalibration.h>
#include <CondTools/BTau/interface/BTagCalibrationReader.h>
#include <TH2F.h>

#include "zLepton.hh"
#include "zJet.hh"
#include "zHLT.hh"
#include "RoccoR.h"

#define LOAD_BRANCH(t, x) {x = 0; t->SetBranchStatus(#x, 1); t->SetBranchAddress(#x, &(x));}
//#define LOAD_BRANCH(t, x) {t->SetBranchAddress(#x, &(x));}

#define EE 0
#define EMu 1
#define MuMu 3

const double lumEraB = 5444711054.147;
const double lumEraC = 2395577279.972;
const double lumEraD = 4255519670.666;
const double lumEraE = 4054218536.816;
const double lumEraF = 3105455590.598;
const double lumEraG = 7544015569.439;
const double lumEraH = 8746001362.368;

const double lumEraBCDEF = (lumEraB + lumEraC + lumEraD + lumEraE + lumEraF);
const double lumEraGH = (lumEraG + lumEraH);

const double lumEraBCDEFGH = (lumEraBCDEF + lumEraGH);

using namespace std;

typedef std::pair<std::string, bool> zFlag;

/*
struct zVertex
{
    float z, dof, RHO;
};
*/
Double_t HTsum(vector<TLorentzVector *> v)
{
    Double_t result = 0;

    for (auto it = v.begin(); it != v.end(); it++)
        result += (*it)->Pt();

    return result;
}

Double_t PTsum(vector<TLorentzVector *> v)
{
    TLorentzVector result;

    for (auto it = v.begin(); it != v.end(); it++)
        result += *(*it);

    return result.Pt();
}

Double_t Et(TLorentzVector v)
{
    return sqrt(v.Mag2() + pow(v.Pt(), 2));
}

class zEvent
{
private:
    vector<zLepton> leptons;
    vector<zJet> jets;

//    vector<zVertex> vertices;
    vector<zHLT> triggers;
    TLorentzVector MET;
    bool isMETok_;
    bool isTrgOk_[4];
    bool BadChargedCandidateFilter_;
    bool BadPFMuonFilter_;

    int puNtrueInteractons;
    double lumiWeight;

public:
    ULong64_t getEvID() const
    {
        return ev_event;
    }

    bool getIsData() const
    {
        return is_data;
    }

    const vector<zHLT> &getTriggers() const
    {
        return triggers;
    }

    const TLorentzVector &getMET() const
    {
        return MET;
    }

    bool isMETok() const
    {
        return isMETok_;
    }

    bool isTrgOk(const int type) const
    {
        return isTrgOk_[type];
    }

    int getPuNtrueInteractons() const
    {
        return puNtrueInteractons;
    }

    const vector<zLepton> &getLeptons() const
    {
        return leptons;
    }

    const vector<zJet> &getJets() const
    {
        return jets;
    }

    void setLumiWeight(double w)
    {
        lumiWeight = w;
    }

private:
    vector<zFlag> event_flags;

public:
    void add_flag(zFlag flag)
    {
        auto it = std::find_if(event_flags.begin(), event_flags.end(),
                               [flag](zFlag item) { return flag.first == item.first; });
        if (it != event_flags.end())
            event_flags.erase(it);

        event_flags.push_back(flag);
    }

    void add_flag(string name, bool value)
    {
        add_flag(make_pair(name, value));
    }

    int get_flag(string flag_name)
    {
        auto it = std::find_if(event_flags.begin(), event_flags.end(),
                               [flag_name](zFlag item) { return flag_name == item.first; });
        if (it != event_flags.end())
            return it->second ? 0 : 1;
        else
            return -1;
    }

    zEvent(TTree *tree, string epoch) : calib("csvv2", "CSVv2_Moriond17_B_H.csv"),
                                        reader(BTagEntry::OP_MEDIUM, "central"),
                                        rc("rcdata.2016.v3")
    {
        reader.load(calib,                // calibration instance
                    BTagEntry::FLAV_B,    // btag flavour
                    "mujets");               // measurement type

        MuIdFile = new TFile *[2];
        MuIdHist = new TH2F *[2];
        MuIsoFile = new TFile *[2];
        MuIsoHist = new TH2F *[2];

        if (epoch == "MC")
        {
            is_data = false;
            MuIdFile[0] = new TFile("MuID_EfficienciesAndSF_BCDEF.root");
            MuIdHist[0] = (TH2F *) ((TDirectoryFile *) MuIdFile[0]->Get(
                    "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta"))->Get(
                    "abseta_pt_ratio");
            MuIsoFile[0] = new TFile("MuIso_EfficienciesAndSF_BCDEF.root");
            MuIsoHist[0] = (TH2F *) ((TDirectoryFile *) MuIsoFile[0]->Get("TightISO_TightID_pt_eta"))->Get(
                    "abseta_pt_ratio");
            MuIdFile[1] = new TFile("MuID_EfficienciesAndSF_GH.root");
            MuIdHist[1] = (TH2F *) ((TDirectoryFile *) MuIdFile[1]->Get(
                    "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta"))->Get(
                    "abseta_pt_ratio");
            MuIsoFile[1] = new TFile("MuIso_EfficienciesAndSF_GH.root");
            MuIsoHist[1] = (TH2F *) ((TDirectoryFile *) MuIsoFile[1]->Get("TightISO_TightID_pt_eta"))->Get(
                    "abseta_pt_ratio");
            EmIdFile = new TFile("egammaEffi.txt_EGM2D.root");
            EmIdHist = (TH2F *) EmIdFile->Get("EGamma_SF2D");
        }
        else
        {
            is_data = true;
//            data_epoch = epoch;
            MuIdHist = NULL;
            MuIsoHist = NULL;
            EmIdHist = NULL;
        }
        tree->SetBranchStatus("*", 0);
//        this->pt_cut_ = 20;
        LOAD_BRANCH(tree, mc_w_sign)
        LOAD_BRANCH(tree, ev_event)
        LOAD_BRANCH(tree, mc_trueNumInteractions)
        LOAD_BRANCH(tree, mc_PU_NumInteractions)
        LOAD_BRANCH(tree, pv_n)
        LOAD_BRANCH(tree, pv_z)
        LOAD_BRANCH(tree, pv_ndof)
        LOAD_BRANCH(tree, pv_normalizedChi2)
        LOAD_BRANCH(tree, pv_isValid)
        LOAD_BRANCH(tree, pv_isFake)
        LOAD_BRANCH(tree, gsf_n)
#ifdef SYNC_EX
        LOAD_BRANCH(tree, gsf_pt)
#else
        LOAD_BRANCH(tree, gsf80_pt)
//        gsf_pt = gsf80_pt;
#endif
        LOAD_BRANCH(tree, gsf_eta)
        LOAD_BRANCH(tree, gsf_phi)
        LOAD_BRANCH(tree, gsf_charge)
        LOAD_BRANCH(tree, gsf_VIDTight)
        LOAD_BRANCH(tree, gsf_dxy_firstPVtx)
        LOAD_BRANCH(tree, gsf_dz_firstPVtx)
        LOAD_BRANCH(tree, gsf_relIso)
        LOAD_BRANCH(tree, gsf_sc_eta)
        LOAD_BRANCH(tree, mu_n)
        LOAD_BRANCH(tree, mu_gt_charge)
        LOAD_BRANCH(tree, mu_gt_pt)
        LOAD_BRANCH(tree, mu_gt_eta)
        LOAD_BRANCH(tree, mu_gt_phi)
        LOAD_BRANCH(tree, mu_gt_d0)
        LOAD_BRANCH(tree, mu_gt_dz)
        LOAD_BRANCH(tree, mu_gt_dz_beamspot)
        LOAD_BRANCH(tree, mu_gt_dz_firstPVtx)
        LOAD_BRANCH(tree, mu_gt_dxy)
        LOAD_BRANCH(tree, mu_gt_dxy_beamspot)
        LOAD_BRANCH(tree, mu_gt_dxy_firstPVtx)
        LOAD_BRANCH(tree, mu_isTightMuon)
        LOAD_BRANCH(tree, mu_isoTrackerBased03)
        LOAD_BRANCH(tree, mu_pfIsoDbCorrected04)
        LOAD_BRANCH(tree, jet_n)
        LOAD_BRANCH(tree, jet_pt)
        LOAD_BRANCH(tree, jet_eta)
        LOAD_BRANCH(tree, jet_phi)
        LOAD_BRANCH(tree, jet_energy)
        LOAD_BRANCH(tree, jet_CSVv2)
        LOAD_BRANCH(tree, jet_isJetIDLoose)
        LOAD_BRANCH(tree, jet_Smeared_pt)
        LOAD_BRANCH(tree, jet_Smeared_energy)
#ifdef SYNC_EX
        LOAD_BRANCH(tree, MET_nominal_Px)
        LOAD_BRANCH(tree, MET_nominal_Py)
#else
        LOAD_BRANCH(tree, MET_T1Txy_Px)
        LOAD_BRANCH(tree, MET_T1Txy_Py)
//        MET_nominal_Px = MET_T1Txy_Px;
//        MET_nominal_Py = MET_T1Txy_Py;
#endif
        LOAD_BRANCH(tree, trig_Flag_HBHENoiseFilter_accept)
        LOAD_BRANCH(tree, trig_Flag_HBHENoiseIsoFilter_accept)
        LOAD_BRANCH(tree, trig_Flag_globalTightHalo2016Filter_accept)
        LOAD_BRANCH(tree, trig_Flag_goodVertices_accept)
        LOAD_BRANCH(tree, trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept)
        LOAD_BRANCH(tree, trig_Flag_BadPFMuonFilter_accept)
        LOAD_BRANCH(tree, trig_Flag_BadChargedCandidateFilter_accept)
        LOAD_BRANCH(tree, trig_Flag_eeBadScFilter_accept)
        LOAD_BRANCH(tree, trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept)
        LOAD_BRANCH(tree, trig_HLT_Ele27_WPTight_Gsf_accept)
        LOAD_BRANCH(tree, trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept)
        LOAD_BRANCH(tree, trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_accept)
        LOAD_BRANCH(tree, trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept)
        LOAD_BRANCH(tree, trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept)
        LOAD_BRANCH(tree, trig_HLT_IsoMu24_accept)
        LOAD_BRANCH(tree, trig_HLT_IsoTkMu24_accept)
        LOAD_BRANCH(tree, trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept)
        LOAD_BRANCH(tree, trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept)
        LOAD_BRANCH(tree, trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept)
        LOAD_BRANCH(tree, trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept)
        LOAD_BRANCH(tree, mu_dB)
        LOAD_BRANCH(tree, mu_isGlobalMuon)
        LOAD_BRANCH(tree, mu_isPFMuon)
        LOAD_BRANCH(tree, mu_gt_normalizedChi2)
        LOAD_BRANCH(tree, mu_numberOfValidMuonHits)
        LOAD_BRANCH(tree, mu_numberOfMatchedStations)
        LOAD_BRANCH(tree, mu_numberOfValidPixelHits)
        LOAD_BRANCH(tree, mu_trackerLayersWithMeasurement)
    }

    Int_t getMc_trueNumInteractions() const
    {
        return mc_trueNumInteractions;
    }

private:
    void reset()
    {
        leptons.clear();
        jets.clear();
        triggers.clear();
//        vertices.clear();
        event_flags.clear();
    }

    void read_electrons()
    {
        for (UInt_t i = 0; i < gsf_eta->size(); i++)
        {
            auto j = static_cast<vector<zLepton>::size_type>(i);
#ifdef SYNC_EX
            if (gsf_pt->at(j) < 20.)
                continue;
#else
            if (gsf80_pt->at(j) < 20.)
                continue;
#endif

            TLorentzVector v;
#ifdef SYNC_EX
            v.SetPtEtaPhiM(gsf_pt->at(j), gsf_eta->at(j), gsf_phi->at(j), 0.000511);
#else
            v.SetPtEtaPhiM(gsf80_pt->at(j), gsf_eta->at(j), gsf_phi->at(j), 0.000511);
#endif
            zLepton thisElectron = zLepton(v, gsf_charge->at(j), gsf_sc_eta->at(j), 0,
                                           gsf_VIDTight->at(j), false, gsf_dxy_firstPVtx->at(j),
                                           gsf_dz_firstPVtx->at(j));
            leptons.push_back(thisElectron);
        }
    }

    void read_muons()
    {
        for (UInt_t i = 0; i < mu_gt_pt->size(); i++)
        {
            auto j = static_cast<vector<zLepton>::size_type>(i);
            if (mu_gt_pt->at(j) < 20.)
                continue;
            TLorentzVector v;
            double pt = mu_gt_pt->at(j);
            if (is_data)
            {
                pt *= rc.kScaleDT(mu_gt_charge->at(j), pt, mu_gt_eta->at(j), mu_gt_phi->at(j), 0, 0);
            }
            else
            {
                pt *= rc.kScaleAndSmearMC(mu_gt_charge->at(j), pt, mu_gt_eta->at(j), mu_gt_phi->at(j),
                                          mu_trackerLayersWithMeasurement->at(j), gRandom->Rndm(), gRandom->Rndm(),
                                          0, 0);
            }
            v.SetPtEtaPhiM(/*mu_gt_pt->at(j)*/pt, mu_gt_eta->at(j), mu_gt_phi->at(j), 0.10566);
            zLepton thisMuon = zLepton(v, mu_gt_charge->at(j), 0,
                                       mu_pfIsoDbCorrected04->at(j),
                                       mu_isTightMuon->at(j), true, 0, 0);
            leptons.push_back(thisMuon);
        }
    }

    void read_jets()
    {
//        cout << "Loading jets" << endl;
        for (UInt_t i = 0; i < jet_pt->size(); i++)
        {
            auto j = static_cast<vector<zJet>::size_type>(i);
            TLorentzVector v;
            v.SetPtEtaPhiE(jet_pt->at(j), jet_eta->at(j), jet_phi->at(j), jet_energy->at(j));
            zJet thisJet = zJet(v, 0, jet_CSVv2->at(j), jet_isJetIDLoose->at(j));
            jets.push_back(thisJet);
        }
//        cout << "Loaded jets" << endl;
        jets.shrink_to_fit();
    }

    void read_vertices()
    {
/*
        for (auto i = 0; i < pv_n; i++)
        {
            auto j = static_cast<vector::size_type>(i);
            zVertex thisVertex;
            thisVertex.z;
        }
*/
    }

    void read_MET()
    {
#ifdef SYNC_EX
        MET = TLorentzVector(MET_nominal_Px, MET_nominal_Py, 0, 0);
#else
        MET = TLorentzVector(MET_T1Txy_Px, MET_T1Txy_Py, 0, 0);
#endif
        isMETok_ = trig_Flag_BadChargedCandidateFilter_accept && trig_Flag_BadPFMuonFilter_accept &&
                   trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept &&
                   trig_Flag_globalTightHalo2016Filter_accept && trig_Flag_goodVertices_accept &&
                   trig_Flag_HBHENoiseFilter_accept && trig_Flag_HBHENoiseIsoFilter_accept;
        if (is_data)
            isMETok_ = isMETok_ && trig_Flag_eeBadScFilter_accept;
    }

    void read_check_trigger()
    {
        // EE
        isTrgOk_[EE] = (trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept + trig_HLT_Ele27_WPTight_Gsf_accept != 0);

        // EMu
        Int_t isTrgOk_EMu = trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept * (is_data ? 1 : 0);
        isTrgOk_EMu += trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept * (is_data ? 1 : 0);
        isTrgOk_EMu += trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept;
        isTrgOk_EMu += trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept;
        isTrgOk_EMu += trig_HLT_Ele27_WPTight_Gsf_accept;
        isTrgOk_EMu += trig_HLT_IsoMu24_accept + trig_HLT_IsoTkMu24_accept;
        isTrgOk_[EMu] = (isTrgOk_EMu != 0);

        Int_t isTrgOk_MuMu = trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept * (is_data ? 1 : 0);
        isTrgOk_MuMu += trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_accept * (is_data ? 1 : 0);
        isTrgOk_MuMu += trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept;
        isTrgOk_MuMu += trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept;
        isTrgOk_MuMu += trig_HLT_IsoMu24_accept + trig_HLT_IsoTkMu24_accept;
        isTrgOk_[MuMu] = (isTrgOk_MuMu != 0);
    }

public:
    /*
    zEvent(edm::EventBase const &ev)
    {
        // Fill event information
#ifndef SYNC_EX
        get_userdata(ev);
        get_triggers(ev);
        get_vertices(ev);
#endif
        edm::Handle<ULong64_t> event_id;
        ev.getByLabel(string("eventInfo:evtInfoEventNumber"), event_id);
        evID = *event_id;
        read_electrons(ev);
        read_muons(ev);
        std::sort(leptons.begin(), leptons.end());
        read_jets(ev);
        clean_jets();
        read_MET(ev);
    }
*/
    bool read_event(TTree *tree, Long64_t id)
    {
        reset();
        tree->GetEntry(id);

        leptons.reserve(gsf_n + mu_n);
        read_electrons();
        read_muons();
        leptons.shrink_to_fit();
        std::sort(leptons.begin(), leptons.end());
        jets.reserve(jet_n);
        read_jets();
        clean_jets();
        read_MET();
        read_check_trigger();
        return true;
    }

private:
    Bool_t trig_Flag_BadPFMuonFilter_accept;
    Bool_t trig_Flag_BadChargedCandidateFilter_accept;
    ULong64_t ev_event;
    Int_t mc_trueNumInteractions;
    Int_t mc_PU_NumInteractions;
    bool is_data;
    Float_t mc_w_sign;
    string data_epoch;

    // Primary vertex
    UInt_t pv_n;
    vector<float> *pv_z;
    vector<float> *pv_ndof;
    vector<float> *pv_normalizedChi2;
    vector<bool> *pv_isValid;
    vector<bool> *pv_isFake;

    // Electron?
    UInt_t gsf_n;
    // vector<float> *gsf_energy;
    // vector<float> *gsf80_energy;
    vector<float> *gsf_pt;
    vector<float> *gsf80_pt;
    vector<float> *gsf_eta;
    vector<float> *gsf_phi;
    vector<int> *gsf_charge;
    vector<bool> *gsf_VIDTight;
    vector<float> *gsf_dxy_firstPVtx;
    vector<float> *gsf_dz_firstPVtx;
    vector<float> *gsf_relIso;
    vector<float> *gsf_sc_eta;

    // Muon
    UInt_t mu_n;
    vector<int> *mu_gt_charge;
    vector<float> *mu_gt_pt;
    vector<float> *mu_gt_eta;
    vector<float> *mu_gt_phi;
    vector<float> *mu_gt_d0;
    vector<float> *mu_gt_dz;
    vector<float> *mu_gt_dz_beamspot;
    vector<float> *mu_gt_dz_firstPVtx;
    vector<float> *mu_gt_dxy;
    vector<float> *mu_gt_dxy_beamspot;
    vector<float> *mu_gt_dxy_firstPVtx;
    vector<bool> *mu_isTightMuon;
    vector<float> *mu_isoTrackerBased03;
    vector<float> *mu_pfIsoDbCorrected04;
    vector<float> *mu_dB;
    vector<bool> *mu_isGlobalMuon;
    vector<bool> *mu_isPFMuon;
    vector<float> *mu_gt_normalizedChi2;
    vector<int> *mu_numberOfValidMuonHits;
    vector<int> *mu_numberOfMatchedStations;
    vector<int> *mu_numberOfValidPixelHits;
    vector<int> *mu_trackerLayersWithMeasurement;

    // Jet
    UInt_t jet_n;
    vector<float> *jet_pt;
    vector<float> *jet_eta;
    vector<float> *jet_phi;
    vector<float> *jet_energy;
    vector<float> *jet_CSVv2;
    vector<bool> *jet_isJetIDLoose;
    vector<float> *jet_Smeared_pt;
    vector<float> *jet_Smeared_energy;

    // MET
    Float_t MET_nominal_Px;
    Float_t MET_nominal_Py;
    Float_t MET_T1Txy_Px;
    Float_t MET_T1Txy_Py;

    // Trigger
    Int_t trig_Flag_HBHENoiseFilter_accept;
    Int_t trig_Flag_HBHENoiseIsoFilter_accept;
    Int_t trig_Flag_globalTightHalo2016Filter_accept;
    Int_t trig_Flag_goodVertices_accept;
    Int_t trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept;
    Int_t trig_Flag_eeBadScFilter_accept;

    //ee triggers
    Int_t trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept;
    Int_t trig_HLT_Ele27_WPTight_Gsf_accept;

    //mumu triggers
    Int_t trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept;
    Int_t trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_accept;
    Int_t trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept;
    Int_t trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept;
    Int_t trig_HLT_IsoMu24_accept;
    Int_t trig_HLT_IsoTkMu24_accept;

    //emu
    Int_t trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept;
    Int_t trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept;
    Int_t trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept;
    Int_t trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept;

    // calibration
    BTagCalibration calib;
    BTagCalibrationReader reader;

    TFile **MuIdFile;
    TH2F **MuIdHist;

    TFile **MuIsoFile;
    TH2F **MuIsoHist;

    TFile *EmIdFile;
    TH2F *EmIdHist;

    RoccoR rc;
/*
    zEvent(edm::EventBase const &ev)
    {
        // Fill event information
#ifndef SYNC_EX
        get_userdata(ev);
        get_triggers(ev);
        get_vertices(ev);
#endif
        edm::Handle<ULong64_t> event_id;
        ev.getByLabel(string("eventInfo:evtInfoEventNumber"), event_id);
        evID = *event_id;
        read_electrons(ev);
        read_muons(ev);
        std::sort(leptons.begin(), leptons.end());
        read_jets(ev);
        clean_jets();
        read_MET(ev);
    }


private:
    void read_electrons(edm::EventBase const &event)
    {
        edm::Handle<std::vector<float> > electronPt;
        event.getByLabel(std::string("electrons:elPt"), electronPt);
        // Handle to the electron eta
        edm::Handle<std::vector<float> > electronEta;
        event.getByLabel(std::string("electrons:elEta"), electronEta);
        // Handle to the electron eta
        edm::Handle<std::vector<float> > electronPhi;
        event.getByLabel(std::string("electrons:elPhi"), electronPhi);
        // Handle to the electron energy
        edm::Handle<std::vector<float> > electronEn;
        event.getByLabel(std::string("electrons:elE"), electronEn);

        // Handle to the electron charge
        edm::Handle<std::vector<float> > electronCharge;
        event.getByLabel(std::string("electrons:elCharge"), electronCharge);
        // Handle to the electron electronSCeta
        edm::Handle<std::vector<float> > electronSCeta;
        event.getByLabel(std::string("electrons:elSCEta"), electronSCeta);

        edm::Handle<std::vector<float> > elvidTight;
        event.getByLabel(std::string("electrons:elvidTight"), elvidTight);

        // Handle to the electron electronD0
        edm::Handle<std::vector<float> > electronD0;
        event.getByLabel(std::string("electrons:elDxy"), electronD0);
        // Handle to the electron electronDz
        edm::Handle<std::vector<float> > electronDz;
        event.getByLabel(std::string("electrons:elDz"), electronDz);

        for (size_t i = 0; i < electronPt->size(); i++)
        {
            TLorentzVector v;
            v.SetPtEtaPhiE(electronPt->at(i), electronEta->at(i), electronPhi->at(i), electronEn->at(i));

            zLepton thisElectron = zLepton(v, electronCharge->at(i), electronSCeta->at(i), 0,
                                           elvidTight->at(i) != 0, false, electronD0->at(i), electronDz->at(i));

            leptons.push_back(thisElectron);
        }
    }

    void read_userdata(edm::EventBase const &event)
    {
        // Handle to the genEventWeight, which I added in eventUserData class
        edm::Handle<double> genEventWeight;
        event.getByLabel(std::string("eventUserData:eventWeight"), genEventWeight);
        weight = *genEventWeight;

        // Handle to the pu, which I added in eventUserData class
        edm::Handle<int> puNtrueIntMC;
        event.getByLabel(std::string("eventUserData:puNtrueInt"), puNtrueIntMC);
        puNtrueInteractons = *puNtrueIntMC;
    }

    void read_triggers(edm::EventBase const &event)
    {
        // Handle to the pass decision of hlt
        edm::Handle<std::vector<float> > HLTdecision;
        event.getByLabel(std::string("TriggerUserData:triggerBitTree"), HLTdecision);
        // Handle to the prescale of hlt
        edm::Handle<std::vector<int> > HLTprescale;
        event.getByLabel(std::string("TriggerUserData:triggerPrescaleTree"), HLTprescale);
        // Handle to the hlt name
        edm::Handle<std::vector<string> > HLTname;
        event.getByLabel(std::string("TriggerUserData:triggerNameTree"), HLTname);

        for (size_t i = 0; i < HLTdecision->size(); i++)
        {
            zHLT thisHLT = zHLT(HLTname->at(i), HLTprescale->at(i), HLTdecision->at(i));
            triggers.push_back(thisHLT);
            // cout << "Trigger " << HLTname->at(i) << " prescale " << HLTprescale->at(i) << " decision "
            //     << HLTdecision->at(i) << endl;
        }
    }

    void read_muons(edm::EventBase const &event)
    {
        // Handle to the muon muonPt
        edm::Handle<std::vector<float> > muonPt;
        event.getByLabel(std::string("muons:muPt"), muonPt);
        // Handle to the muon muonEta
        edm::Handle<std::vector<float> > muonEta;
        event.getByLabel(std::string("muons:muEta"), muonEta);
        // Handle to the muon muonPhi
        edm::Handle<std::vector<float> > muonPhi;
        event.getByLabel(std::string("muons:muPhi"), muonPhi);
        // Handle to the muon energy
        edm::Handle<std::vector<float> > muonEn;
        event.getByLabel(std::string("muons:muE"), muonEn);

        // Handle to the muon muonIso04
        edm::Handle<std::vector<float> > muonIso04;
        event.getByLabel(std::string("muons:muIso04"), muonIso04);
        // Handle to the muon charge
        edm::Handle<std::vector<float> > muonCharge;
        event.getByLabel(std::string("muons:muCharge"), muonCharge);
        // Handle to the muon charge
        edm::Handle<std::vector<float> > muonTight;
        event.getByLabel(std::string("muons:muIsTightMuon"), muonTight);

        for (size_t i = 0; i < muonPt->size(); i++)
        {
            TLorentzVector v;
            v.SetPtEtaPhiE(muonPt->at(i), muonEta->at(i), muonPhi->at(i), muonEn->at(i));

            zLepton thisMuon(v, muonCharge->at(i), 0, muonIso04->at(i), muonTight->at(i) != 0, true, 0, 0);
            leptons.push_back(thisMuon);
        }
    }

    void read_jets(edm::EventBase const &event)
    {
        // Handle to the jet b-tag
        edm::Handle<std::vector<float> > jetBTag;
        event.getByLabel(std::string("jetsAK4CHS:jetAK4CHSCSVv2"), jetBTag);
        // Handle to the jet charge
        edm::Handle<std::vector<float> > jetCharge;
        event.getByLabel(std::string("jetsAK4CHS:jetAK4CHSCharge"), jetCharge);
        // Handle to the jet energy
        edm::Handle<std::vector<float> > jetEn;
        event.getByLabel(std::string("jetsAK4CHS:jetAK4CHSE"), jetEn);
        // Handle to the jet eta
        edm::Handle<std::vector<float> > jetEta;
        event.getByLabel(std::string("jetsAK4CHS:jetAK4CHSEta"), jetEta);
        // Handle to the jet rapidity
        edm::Handle<std::vector<float> > jetRapidity;
        event.getByLabel(std::string("jetsAK4CHS:jetAK4CHSY"), jetRapidity);
        // Handle to the jet mass
        edm::Handle<std::vector<float> > jetMass;
        event.getByLabel(std::string("jetsAK4CHS:jetAK4CHSMass"), jetMass);
        // Handle to the jet phi
        edm::Handle<std::vector<float> > jetPhi;
        event.getByLabel(std::string("jetsAK4CHS:jetAK4CHSPhi"), jetPhi);
        // Handle to the jet pt
        edm::Handle<std::vector<float> > jetPt;
        event.getByLabel(std::string("jetsAK4CHS:jetAK4CHSPt"), jetPt);
        // Handle to the jet Area
        edm::Handle<std::vector<float> > jetArea;
        event.getByLabel(std::string("jetsAK4CHS:jetAK4CHSjetArea"), jetArea);
        ///////////////Jet Loose criteria//////////////////////////////////////////////////////////////////////
        // Handle to the Neutral hadron fraction
        edm::Handle<std::vector<float> > jetneutralHadronEnergyFrac;
        event.getByLabel(std::string("jetsAK4CHS:jetAK4CHSneutralHadronEnergyFrac"), jetneutralHadronEnergyFrac);
        // Handle to the Neutral EM fraction
        edm::Handle<std::vector<float> > jetneutralEmEnergyFrac;
        event.getByLabel(std::string("jetsAK4CHS:jetAK4CHSneutralEmEnergyFrac"), jetneutralEmEnergyFrac);
        // Handle to the Number of constituents
        edm::Handle<std::vector<float> > jetNumConstituents;
        event.getByLabel(std::string("jetsAK4CHS:jetAK4CHSNumConstituents"), jetNumConstituents);
        // Handle to the Charged Hadron fraction
        edm::Handle<std::vector<float> > jetchargedHadronEnergyFrac;
        event.getByLabel(std::string("jetsAK4CHS:jetAK4CHSchargedHadronEnergyFrac"), jetchargedHadronEnergyFrac);
        // Handle to the Charged Multiplicity
        edm::Handle<std::vector<float> > jetchargedMultiplicity;
        event.getByLabel(std::string("jetsAK4CHS:jetAK4CHSchargedMultiplicity"), jetchargedMultiplicity);
        // Handle to the Charged EM fraction
        edm::Handle<std::vector<float> > jetchargedEmEnergyFrac;
        event.getByLabel(std::string("jetsAK4CHS:jetAK4CHSchargedEmEnergyFrac"), jetchargedEmEnergyFrac);
        // Handle to the Neutral particle multiplicity
        edm::Handle<std::vector<float> > jetneutralMultiplicity;
        event.getByLabel(std::string("jetsAK4CHS:jetAK4CHSneutralMultiplicity"), jetneutralMultiplicity);
        /////////////////Jet b-tag/////////////////////////////////////////////////////////////////////////////////


        for (size_t i = 0; i < jetBTag->size(); i++)
        {
            TLorentzVector v;
            v.SetPtEtaPhiE(jetPt->at(i), jetEta->at(i), jetPhi->at(i), jetEn->at(i));

            zJet this_jet = zJet(v, (int) jetCharge->at(i), jetBTag->at(i), jetRapidity->at(i), jetArea->at(i),
                                 jetneutralHadronEnergyFrac->at(i), jetneutralEmEnergyFrac->at(i),
                                 jetchargedHadronEnergyFrac->at(i), jetchargedEmEnergyFrac->at(i),
                                 jetNumConstituents->at(i), jetchargedMultiplicity->at(i),
                                 jetneutralMultiplicity->at(i));
            jets.push_back(this_jet);
        }

        std::sort(jets.begin(), jets.end());
    }

    void read_MET(edm::EventBase const &event)
    {
        edm::Handle<std::vector<float> > MetPx;
        event.getByLabel(std::string("metFull:metFullPx"), MetPx);

        edm::Handle<std::vector<float> > MetPy;
        event.getByLabel(std::string("metFull:metFullPy"), MetPy);

        edm::Handle<bool> BadChargedCandidateFilter;
        event.getByLabel(std::string("BadChargedCandidateFilter"), BadChargedCandidateFilter);
        BadChargedCandidateFilter_ = *BadChargedCandidateFilter;

        edm::Handle<bool> BadPFMuonFilter;
        event.getByLabel(std::string("BadPFMuonFilter"), BadPFMuonFilter);
        BadPFMuonFilter_ = *BadPFMuonFilter;

        edm::Handle<std::vector<string>> METTriggerNameTree;
        event.getByLabel(std::string("METUserData:triggerNameTree"), METTriggerNameTree);

        edm::Handle<std::vector<float>> METTriggerBitTree;
        event.getByLabel(std::string("METUserData:triggerBitTree"), METTriggerBitTree);

        float otherMetFilters = 1;

        for (size_t i = 0; i < METTriggerNameTree->size(); i++)
        {
            if (METTriggerNameTree->at(i) == "Flag_HBHENoiseFilter" ||
                METTriggerNameTree->at(i) == "Flag_HBHENoiseIsoFilter" ||
                METTriggerNameTree->at(i) == "Flag_globalTightHalo2016Filter" ||
                METTriggerNameTree->at(i) == "Flag_goodVertices" ||
                METTriggerNameTree->at(i) == "Flag_EcalDeadCellTriggerPrimitiveFilter")
            {
                otherMetFilters *= METTriggerBitTree->at(i);
            }
        }

        MET = TLorentzVector(MetPx->at(0), MetPy->at(0), 0, 0);
        isMETok_ = ((*BadChargedCandidateFilter) && (*BadPFMuonFilter) && (otherMetFilters == 1));
    }

    void read_vertices(edm::EventBase const &event)
    {
        edm::Handle<std::vector<float> > vtxZ;
        event.getByLabel(std::string("vertexInfo:z"), vtxZ);
        // Handle to the dof of vertex
        edm::Handle<std::vector<int> > dof;
        event.getByLabel(std::string("vertexInfo:ndof"), dof);
        // Handle to the rho
        edm::Handle<std::vector<float> > rhoo;
        event.getByLabel(std::string("vertexInfo:rho"), rhoo);

        for (size_t i = 0; i < vtxZ->size(); i++)
        {
            zVertex this_vtx;
            this_vtx.z = vtxZ->at(i);
            this_vtx.dof = dof->at(i);
            this_vtx.RHO = rhoo->at(i);

            vertices.push_back(this_vtx);
        }
    }
*/

public:
    void clean_jets()
    {
        for (auto thisJet = jets.begin(); thisJet != jets.end(); thisJet++)
        {
            bool jet_flag = std::all_of(leptons.begin(), leptons.end(),
                                        [thisJet](const zLepton &thisLepton) {
                                            return !(thisLepton.is_selected() && thisJet->deltaR(thisLepton) <= 0.4);
                                        });

            thisJet->set_clean_flag(jet_flag);
        }
    }

    void fill_dump_tree(TTree *tree)
    {
        std::vector<TLorentzVector> *JetV, *LepV;
        std::vector<float> *JetCh, *LepCh, *LepDxy, *LepDz, *LepEtaSC, *LepIsoVal;
        std::vector<int> *LepWhere;
        std::vector<bool> *JetB, *JetClean, *JetSel, *JetLoose, *LepSel, *LepIso, *LepTight, *LepMuon;
        std::vector<string> *flags;

        JetV = new std::vector<TLorentzVector>();
        LepV = new std::vector<TLorentzVector>();
        JetCh = new std::vector<float>();
        LepCh = new std::vector<float>();

        JetB = new std::vector<bool>();
        JetClean = new std::vector<bool>();
        JetSel = new std::vector<bool>();
        JetLoose = new std::vector<bool>();
        LepSel = new std::vector<bool>();
        LepWhere = new std::vector<int>();
        LepIso = new std::vector<bool>();
        LepIsoVal = new std::vector<float>();
        LepTight = new std::vector<bool>();
        LepMuon = new std::vector<bool>();
        LepDxy = new std::vector<float>();
        LepDz = new std::vector<float>();
        LepEtaSC = new std::vector<float>();

        flags = new std::vector<string>();

        for (auto thisLepton = leptons.begin(); thisLepton != leptons.end(); thisLepton++)
        {
            LepV->push_back(static_cast<TLorentzVector>(*thisLepton));
            LepCh->push_back(thisLepton->get_charge());
            LepSel->push_back(thisLepton->is_selected());
            LepWhere->push_back(thisLepton->where());
            LepIso->push_back(thisLepton->is_iso());
            LepTight->push_back(thisLepton->is_tight());
            LepMuon->push_back(thisLepton->is_muon());
            LepDxy->push_back(thisLepton->get_d0());
            LepDz->push_back(thisLepton->get_dz());
            LepEtaSC->push_back(thisLepton->get_etaSC());
            LepIsoVal->push_back(thisLepton->get_Iso());
        }

        for (auto thisJet = jets.begin(); thisJet != jets.end(); thisJet++)
        {
            JetV->push_back(static_cast<TLorentzVector>(*thisJet));
            JetCh->push_back(thisJet->get_charge());
            JetSel->push_back(thisJet->is_selected());
            JetClean->push_back(thisJet->is_clean());
            JetB->push_back(thisJet->is_bjet());
            JetLoose->push_back(thisJet->is_loose());
        }

        for (auto thisFlag = event_flags.begin(); thisFlag != event_flags.end(); thisFlag++)
        {
            if (thisFlag->second)
                flags->push_back(thisFlag->first);
        }

        tree->SetBranchAddress("LeptonVec", &LepV);
        tree->SetBranchAddress("LeptonCharge", &LepCh);
        tree->SetBranchAddress("LeptonSelected", &LepSel);
        tree->SetBranchAddress("LeptonWhere", &LepWhere);
        tree->SetBranchAddress("LeptonIsIso", &LepIso);
        tree->SetBranchAddress("LeptonIsTight", &LepTight);
        tree->SetBranchAddress("LeptonIsMuon", &LepMuon);
        tree->SetBranchAddress("LeptonDxy", &LepDxy);
        tree->SetBranchAddress("LeptonDz", &LepDz);
        tree->SetBranchAddress("LeptonEtaSC", &LepEtaSC);
        tree->SetBranchAddress("LeptonIso", &LepIsoVal);

        tree->SetBranchAddress("JetVec", &JetV);
        tree->SetBranchAddress("JetCharge", &JetCh);
        tree->SetBranchAddress("JetSelected", &JetSel);
        tree->SetBranchAddress("JetIsClean", &JetClean);
        tree->SetBranchAddress("JetIsBJet", &JetB);
        tree->SetBranchAddress("JetIsLoose", &JetLoose);

        bool isMetOk = isMETok();
        TLorentzVector *pMet = &MET;

        tree->SetBranchAddress("MetVec", &pMet);
        tree->SetBranchAddress("MetOK", &isMetOk);
        tree->SetBranchAddress("MetBadCCF", &BadChargedCandidateFilter_);
        tree->SetBranchAddress("MetBadPFM", &BadPFMuonFilter_);

        tree->SetBranchAddress("Flags", &flags);
        tree->SetBranchAddress("eventID", &ev_event);

        tree->Fill();
    }

    void fill_BDT_tree(TTree *tree)
    {
        vector<zLepton> selectedLeptons;
        vector<zJet> selectedJets, selectedBJets;
        copy_if(getLeptons().begin(), getLeptons().end(), back_inserter(selectedLeptons),
                [](const zLepton &part) { return part.is_selected(); });

        copy_if(getJets().begin(), getJets().end(), back_inserter(selectedJets),
                [](const zJet &jet) { return jet.is_selected(); });

        copy_if(selectedJets.begin(), selectedJets.end(), back_inserter(selectedBJets),
                [](const zJet &jet) { return jet.is_bjet(); });

        // Apply SFs to MC
        if (!is_data)
        {
            double jet_scalefactor = 1.0;
            for (auto it = selectedBJets.begin(); it != selectedBJets.end(); it++)
            {
                jet_scalefactor *= reader.eval_auto_bounds("central", BTagEntry::FLAV_B, static_cast<float>(it->Eta()),
                                                           static_cast<float>(it->Pt()));
            }

            mc_w_sign *= jet_scalefactor;

            for (auto it = selectedLeptons.begin(); it != selectedLeptons.end(); it++)
            {
                if (it->is_muon())
                {
                    double mu_id_sf[2] = {MuIdHist[0]->GetBinContent(MuIdHist[0]->FindBin(abs(it->Eta()), it->Pt())),
                                          MuIdHist[1]->GetBinContent(MuIdHist[1]->FindBin(abs(it->Eta()), it->Pt()))};
                    double mu_iso_sf[2] = {MuIsoHist[0]->GetBinContent(MuIsoHist[0]->FindBin(abs(it->Eta()), it->Pt())),
                                           MuIsoHist[1]->GetBinContent(
                                                   MuIsoHist[1]->FindBin(abs(it->Eta()), it->Pt()))};


                    mc_w_sign *= mu_id_sf[0] * mu_iso_sf[0] * lumEraBCDEF / lumEraBCDEFGH;
                    mc_w_sign *= mu_id_sf[1] * mu_iso_sf[1] * lumEraGH / lumEraBCDEFGH;
                }
                else
                {
                    double em_id_sf = EmIdHist->GetBinContent(EmIdHist->FindBin(it->get_etaSC(), it->Pt()));
                    mc_w_sign *= em_id_sf;
                }
            }

            mc_w_sign *= lumiWeight;
        }

        TLorentzVector lep1(selectedLeptons.at(0));
        TLorentzVector lep2(selectedLeptons.at(1));
        TLorentzVector bjet1(0, 0, 0, 0);
        TLorentzVector jet1(0, 0, 0, 0);
        TLorentzVector jet2(0, 0, 0, 0);

        TLorentzVector met(getMET());
        if (selectedJets.size() > 0)
            jet1 = selectedJets.at(0);

        if (selectedJets.size() > 1)
            jet2 = selectedJets.at(1);

        if (selectedBJets.size() > 0)
            bjet1 = selectedBJets.at(0);

        Double_t ptllmetj = (lep1 + lep2 + jet1 + met).Pt();
        tree->SetBranchAddress("ptsys", &ptllmetj);

        Double_t dptll_metj = (lep1 + lep2).Pt() - (met + jet1).Pt();
        tree->SetBranchAddress("dpt_ll_metj", &dptll_metj);

        Double_t met_pt = getMET().Pt();
        tree->SetBranchAddress("MET", &met_pt);

        Double_t dptll_met = (lep1 + lep2).Pt() - met_pt;
        tree->SetBranchAddress("dpt_ll_met", &dptll_met);

        Double_t ptlmetj = (lep1 + jet1 + met).Pt();
        tree->SetBranchAddress("pt_lMETj", &ptlmetj);

        Double_t mll = (lep1 + lep2).Mag();
        tree->SetBranchAddress("mll", &mll);

        Double_t cll = Et(lep1 + lep2) / ((lep1 + lep2).Pt());
        tree->SetBranchAddress("cll", &cll);

        Double_t dptl_met = lep1.Pt() - met_pt;
        tree->SetBranchAddress("dpt_l_met", &dptl_met);

        int njet = static_cast<int>(selectedJets.size());
        tree->SetBranchAddress("njets", &njet);

        int nbjet = static_cast<int>(selectedBJets.size());
        tree->SetBranchAddress("nbjets", &nbjet);

        vector<TLorentzVector *> vec;
        vec.push_back(&lep1);
        vec.push_back(&lep2);
        Double_t htll = HTsum(vec);

        vec.push_back(&jet1);
        vec.push_back(&met);
        Double_t htlljmet = HTsum(vec);
        tree->SetBranchAddress("htlljmet", &htlljmet);

        Double_t ptj = jet1.Pt();
        tree->SetBranchAddress("ptj", &ptj);

        Double_t pt_over_ht = ptllmetj / htlljmet;
        tree->SetBranchAddress("pt_over_ht", &pt_over_ht);

        Double_t mlljmet = (lep1 + lep2 + jet1 + getMET()).Mag();
        tree->SetBranchAddress("mlljmet", &mlljmet);

        Double_t ptllj = (lep1 + lep2 + jet1).Pt();
        tree->SetBranchAddress("ptllj", &ptllj);

        Double_t htll_over_ht = htll / htlljmet;
        tree->SetBranchAddress("htll_over_ht", &htll_over_ht);

        Double_t ptll = (lep1 + lep2).Pt();
        tree->SetBranchAddress("ptll", &ptll);

        Double_t drll_metjj = (lep1 + lep2).DeltaR(met + jet1 + jet2);
        tree->SetBranchAddress("dr_llmetjj", &drll_metjj);

        Double_t drll_jj = (lep1 + lep2).DeltaR(jet1 + jet2);
        tree->SetBranchAddress("dr_lljj", &drll_jj);

        Double_t mlj2 = (lep1 + jet2).Mag();
        tree->SetBranchAddress("mlj2", &mlj2);

        Double_t dptl_j = lep1.Pt() - jet1.Pt();
        tree->SetBranchAddress("dpt_l_j", &dptl_j);

        Double_t mlj = (lep1 + jet1).Mag();
        tree->SetBranchAddress("mlj", &mlj);

        Double_t ptl = lep1.Pt();
        tree->SetBranchAddress("ptll", &ptl);

        Double_t ptsl = lep2.Pt();
        tree->SetBranchAddress("ptsl", &ptsl);

        Double_t drl_j = lep1.DeltaR(jet1);
        tree->SetBranchAddress("dr_l_j", &drl_j);

        Double_t ptj2 = jet2.Pt();
        tree->SetBranchAddress("ptj2", &ptj2);

        Double_t ml2jj = (lep2 + jet1 + jet2).Mag();
        tree->SetBranchAddress("ml2jj", &ml2jj);

        Double_t ml2j1 = (lep2 + jet1).Mag();
        tree->SetBranchAddress("ml2j1", &ml2j1);

        Double_t ml2j2 = (lep2 + jet2).Mag();
        tree->SetBranchAddress("ml2j2", &ml2j2);

        float j1tag = -1;
        if (selectedJets.size() > 0)
            j1tag = selectedJets.at(0).get_csv();

        tree->SetBranchAddress("j1csv", &j1tag);

        tree->SetBranchAddress("mc_w_sign", &mc_w_sign);

        int chan = 0;
        if (std::find(event_flags.begin(), event_flags.end(), "EE") != event_flags.end())
            chan = EE;
        else
        {
            if (std::find(event_flags.begin(), event_flags.end(), "EMu") != event_flags.end())
                chan = EMu;
            else if (std::find(event_flags.begin(), event_flags.end(), "MuMu") != event_flags.end())
                chan = MuMu;
        }

        tree->SetBranchAddress("channel", &chan);

        tree->Fill();
    }
};

#endif //INC_13TEV_TW_ANALYSIS_ZEVENT_HH

#pragma clang diagnostic pop