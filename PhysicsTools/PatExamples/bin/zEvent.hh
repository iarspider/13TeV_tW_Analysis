#ifndef INC_13TEV_TW_ANALYSIS_ZEVENT_HH
#define INC_13TEV_TW_ANALYSIS_ZEVENT_HH

#include <vector>
#include <string>

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

class zEvent
{
private:
    vector<zLepton> leptons;
    vector<zJet> jets;

    vector<zVertex> vertices;
    vector<zHLT> triggers;
    TLorentzVector MET;
    bool isMETok_;

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
#ifndef TW_SYNC
        get_userdata(ev);
        get_triggers(ev);
        get_vertices(ev);
#endif
        edm::Handle<ULong64_t> event_id;
        ev.getByLabel(string("eventInfo:evtInfoEventNumber"), event_id);
        this->evID = *event_id;
        get_electrons(ev);
        get_muons(ev);
        std::sort(this->leptons.begin(), this->leptons.end());
        get_jets(ev);
        clean_jets();
        get_MET(ev);
    }


private:
    void get_electrons(edm::EventBase const &event)
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

        for (size_t i = 0; i < electronPt->size(); i++)
        {
            TLorentzVector v;
            v.SetPtEtaPhiE(electronPt->at(i), electronEta->at(i), electronPhi->at(i), electronEn->at(i));

            zLepton thisElectron = zLepton(v, electronCharge->at(i), electronSCeta->at(i), 0,
                                           elvidTight->at(i) != 0, false);

            this->leptons.push_back(thisElectron);
        }
    }

    void get_userdata(edm::EventBase const &event)
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

    void get_triggers(edm::EventBase const &event)
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

    void get_muons(edm::EventBase const &event)
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

            zLepton thisMuon(v, muonCharge->at(i), 0, muonIso04->at(i), muonTight->at(i) != 0, true);
            this->leptons.push_back(thisMuon);
        }
    }

    void get_jets(edm::EventBase const &event)
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

    void get_MET(edm::EventBase const &event)
    {
        edm::Handle<std::vector<float> > MetPx;
        event.getByLabel(std::string("metFull:metFullPx"), MetPx);

        edm::Handle<std::vector<float> > MetPy;
        event.getByLabel(std::string("metFull:metFullPy"), MetPy);

        edm::Handle<bool> BadChargedCandidateFilter;
        event.getByLabel(std::string("BadChargedCandidateFilter"), BadChargedCandidateFilter);

        edm::Handle<bool> BadPFMuonFilter;
        event.getByLabel(std::string("BadPFMuonFilter"), BadPFMuonFilter);

        this->MET = TLorentzVector(MetPx->at(0), MetPy->at(0), 0, 0);
        this->isMETok_ = !((*BadChargedCandidateFilter) && (*BadPFMuonFilter));
    }

    void get_vertices(edm::EventBase const &event)
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
        std::vector<float> *JetCh, *LepCh;
        std::vector<bool> *JetB, *JetClean, *JetSel, *JetLoose, *LepSel, *LepGap, *LepIso, *LepTight, *LepMuon;
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
        LepGap = new std::vector<bool>();
        LepIso = new std::vector<bool>();
        LepTight = new std::vector<bool>();
        LepMuon = new std::vector<bool>();

        flags = new std::vector<string>();

        for (auto thisLepton = leptons.begin(); thisLepton != leptons.end(); thisLepton++)
        {
            LepV->push_back(static_cast<TLorentzVector>(*thisLepton));
            LepCh->push_back(thisLepton->get_charge());
            LepSel->push_back(thisLepton->is_selected());
            LepGap->push_back(thisLepton->in_gap());
            LepIso->push_back(thisLepton->is_iso());
            LepTight->push_back(thisLepton->is_tight());
            LepMuon->push_back(thisLepton->is_muon());
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
        tree->SetBranchAddress("LeptonInGap", &LepGap);
        tree->SetBranchAddress("LeptonIsIso", &LepIso);
        tree->SetBranchAddress("LeptonIsTight", &LepTight);
        tree->SetBranchAddress("LeptonIsMuon", &LepMuon);

        tree->SetBranchAddress("JetVec", &JetV);
        tree->SetBranchAddress("JetCharge", &JetCh);
        tree->SetBranchAddress("JetSelected", &JetSel);
        tree->SetBranchAddress("JetIsClean", &JetClean);
        tree->SetBranchAddress("JetIsBJet", &JetB);
        tree->SetBranchAddress("JetIsLoose", &JetLoose);

        bool isMetOk = isMETok();
        TLorentzVector* pMet = &MET;

        tree->SetBranchAddress("MetVec", &pMet);
        tree->SetBranchAddress("MetOK", &isMetOk);

        tree->SetBranchAddress("Flags", &flags);
        tree->SetBranchAddress("eventID", &evID);

        tree->Fill();
    }
};

#endif //INC_13TEV_TW_ANALYSIS_ZEVENT_HH
