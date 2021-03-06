#ifndef INC_13TEV_TW_ANALYSIS_ZEVENT_HH
#define INC_13TEV_TW_ANALYSIS_ZEVENT_HH

#include <vector>
#include <string>
#include <math.h>

#include <DataFormats/FWLite/interface/Event.h>

#include "zLepton.hh"
#include "zJet.hh"
#include "zHLT.hh"


#define EE 0
#define EMu 1
#define MuMu 3

using namespace std;

typedef std::pair<std::string, bool> zFlag;

struct zVertex
{
    float z, dof, RHO;
};

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

    vector<zVertex> vertices;
    vector<zHLT> triggers;
    TLorentzVector MET;
    bool isMETok_;
    bool BadChargedCandidateFilter_;
    bool BadPFMuonFilter_;

    int puNtrueInteractons;
    double weight;

    ULong64_t evID;
public:
    ULong64_t getEvID() const
    {
        return evID;
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

    const vector<zVertex> &getVertices() const
    {
        return vertices;
    }

    double getWeight() const
    {
        return weight;
    }


private:
    vector<zFlag> event_flags;

public:
    void add_flag(zFlag flag)
    {
        auto it = std::find_if(this->event_flags.begin(), this->event_flags.end(),
                               [flag](zFlag item) { return flag.first == item.first; });
        if (it != this->event_flags.end())
            event_flags.erase(it);

        this->event_flags.push_back(flag);
    }

    void add_flag(string name, bool value)
    {
        add_flag(make_pair(name, value));
    }

    int get_flag(string flag_name)
    {
        auto it = std::find_if(this->event_flags.begin(), this->event_flags.end(),
                               [flag_name](zFlag item) { return flag_name == item.first; });
        if (it != this->event_flags.end())
            return it->second ? 0 : 1;
        else
            return -1;
    }

    bool has_trigger(string name)
    {
        return std::find_if(this->triggers.cbegin(), this->triggers.cend(),
                            [name](zHLT trig) {
                                return trig.HLTName == name && trig.HTLDecision == 1 && trig.HTLPrescale >= 1;
                            }) != this->triggers.cend();
    }

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
        this->evID = *event_id;
        read_electrons(ev);
        read_muons(ev);
        std::sort(this->leptons.begin(), this->leptons.end());
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

            this->leptons.push_back(thisElectron);
        }
    }

    void read_userdata(edm::EventBase const &event)
    {
        // Handle to the genEventWeight, which I added in eventUserData class
        edm::Handle<double> genEventWeight;
        event.getByLabel(std::string("eventUserData:eventWeight"), genEventWeight);
        this->weight = *genEventWeight;

        // Handle to the pu, which I added in eventUserData class
        edm::Handle<int> puNtrueIntMC;
        event.getByLabel(std::string("eventUserData:puNtrueInt"), puNtrueIntMC);
        this->puNtrueInteractons = *puNtrueIntMC;
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
            this->triggers.push_back(thisHLT);
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
            this->leptons.push_back(thisMuon);
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
            this->jets.push_back(this_jet);
        }

        std::sort(this->jets.begin(), this->jets.end());
    }

public:
    void clean_jets()
    {
        for (auto thisJet = this->jets.begin(); thisJet != this->jets.end(); thisJet++)
        {
            bool jet_flag = std::all_of(this->leptons.begin(), this->leptons.end(),
                                        [thisJet](const zLepton &thisLepton) {
                                            return !(thisLepton.is_selected() && thisJet->deltaR(thisLepton) <= 0.4);
                                        });

            thisJet->set_clean_flag(jet_flag);
        }
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

        this->MET = TLorentzVector(MetPx->at(0), MetPy->at(0), 0, 0);
        this->isMETok_ = ((*BadChargedCandidateFilter) && (*BadPFMuonFilter) && (otherMetFilters == 1));
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

            this->vertices.push_back(this_vtx);
        }
    }

    void fill_tree(TTree *tree)
    {
        std::vector<TLorentzVector> *JetV, *LepV;
        std::vector<float> *JetCh, *LepCh, *LepDxy, *LepDz, *LepEtaSC, *LepIso;
        std::vector<int> *LepWhere;
        std::vector<bool> *JetB, *JetClean, *JetSel, *JetLoose, *LepSel, *LepIsIso, *LepTight, *LepMuon;
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
        LepIsIso = new std::vector<bool>();
        LepTight = new std::vector<bool>();
        LepMuon = new std::vector<bool>();
        LepDxy = new std::vector<float>();
        LepDz = new std::vector<float>();
        LepEtaSC = new std::vector<float>();
        LepIso = new std::vector<float>();

        flags = new std::vector<string>();

        for (auto thisLepton = leptons.begin(); thisLepton != leptons.end(); thisLepton++)
        {
            LepV->push_back(static_cast<TLorentzVector>(*thisLepton));
            LepCh->push_back(thisLepton->get_charge());
            LepSel->push_back(thisLepton->is_selected());
            LepWhere->push_back(thisLepton->where());
            LepIsIso->push_back(thisLepton->is_iso());
            LepTight->push_back(thisLepton->is_tight());
            LepMuon->push_back(thisLepton->is_muon());
            LepDxy->push_back(thisLepton->get_d0());
            LepDz->push_back(thisLepton->get_dz());
            LepEtaSC->push_back(thisLepton->get_etaSC());
            LepIso->push_back(thisLepton->get_iso());
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
        tree->SetBranchAddress("LeptonIsIso", &LepIsIso);
        tree->SetBranchAddress("LeptonIsTight", &LepTight);
        tree->SetBranchAddress("LeptonIsMuon", &LepMuon);
        tree->SetBranchAddress("LeptonDxy", &LepDxy);
        tree->SetBranchAddress("LeptonDz", &LepDz);
        tree->SetBranchAddress("LeptonEtaSC", &LepEtaSC);
        tree->SetBranchAddress("LeptonIso", &LepIso);

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
        tree->SetBranchAddress("eventID", &evID);

        tree->Fill();
    }

    void fill_tree_2(TTree *tree)
    {
        vector<zLepton> selectedLeptons;
        vector<zJet> selectedJets;
        vector<zJet> selectedBJets;

        copy_if(this->getLeptons().begin(), this->getLeptons().end(), back_inserter(selectedLeptons),
                [](const zLepton &part) { return part.is_selected(); });

        copy_if(this->getJets().begin(), this->getJets().end(), back_inserter(selectedJets),
                [](const zJet &jet) { return jet.is_selected(); });

        copy_if(selectedJets.begin(), selectedJets.end(), back_inserter(selectedBJets),
                [](const zJet &jet) { return jet.is_bjet(); });

        TLorentzVector lep1(selectedLeptons.at(0));
        TLorentzVector lep2(selectedLeptons.at(0));
        TLorentzVector bjet1(selectedBJets.at(0));
        TLorentzVector jet1(selectedJets.at(0));
        TLorentzVector jet2(0, 0, 0, 0);
        TLorentzVector met(this->getMET());
        if (selectedJets.size() > 1)
        {
            jet2 = selectedJets.at(1);
        }

        Double_t ptllmetj = (lep1 + lep2 + jet1 + met).Pt();
        tree->SetBranchAddress("ptsys", &ptllmetj);

        Double_t dptll_metj = (lep1 + lep2).Pt() - (met + jet1).Pt();
        tree->SetBranchAddress("dpt_ll_metj", &dptll_metj);

        Double_t met_pt = this->getMET().Pt();
        tree->SetBranchAddress("MET", &met_pt);

        Double_t dptll_met = (lep1 + lep2).Pt() - met_pt;
        tree->SetBranchAddress("dpt_ll_met", &dptll_met);

        Double_t ptlmetj = (lep1 + jet1 + met).Pt();
        tree->SetBranchAddress("pt_lMETj", &ptllmetj);

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

        Double_t mlljmet = (lep1 + lep2 + jet1 + this->getMET()).Mag();
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
        tree->SetBranchAddress("ptl", &ptl);

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

        tree->Fill();
    }
};

#endif //INC_13TEV_TW_ANALYSIS_ZEVENT_HH
