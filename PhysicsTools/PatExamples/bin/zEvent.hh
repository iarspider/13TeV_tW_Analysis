#ifndef INC_13TEV_TW_ANALYSIS_ZEVENT_HH
#define INC_13TEV_TW_ANALYSIS_ZEVENT_HH

#include <vector>
#include <string>

#include <DataFormats/FWLite/interface/Event.h>

#include "zElectron.hh"
#include "zMuon.hh"
#include "zJet.hh"
#include "zHLT.hh"

using namespace std;

typedef std::pair<std::string, bool> zFlag;

struct zVertex
{
    float z, dof, RHO;
};

class zEvent
{
private:
    vector<zElectron> electrons;
    vector<zMuon> muons;
    vector<zJet> jets;
    vector<zFlag> selection_pass_flags;
    vector<zVertex> vertices;

public:
    void add_flag(zFlag flag)
    {
        auto it = std::find_if(this->selection_pass_flags.begin(), this->selection_pass_flags.end(),
                               [flag](zFlag item) { return flag.first == item.first; });
        if (it != this->selection_pass_flags.end())
            selection_pass_flags.erase(it);

        this->selection_pass_flags.push_back(flag);
    }

    int get_flag(string flag_name)
    {
        auto it = std::find_if(this->selection_pass_flags.begin(), this->selection_pass_flags.end(),
                               [flag_name](zFlag item) { return flag_name == item.first; });
        if (it != this->selection_pass_flags.end())
            return it->second ? 0 : 1;
        else
            return -1;
    }

    vector<zHLT> triggers;

    TLorentzVector MET;
    bool isMETok;

    int puNtrueInteractons;
    double weight;
public:
    double getWeight() const
    {
        return weight;
    }

private:

    vector<bool> is_clean_jet;

public:
    zEvent(edm::EventBase const &ev)
    {
        // Fill event information
        get_userdata(ev);
        get_triggers(ev);
        get_vertices(ev);
        get_electrons(ev);
        get_muons(ev);
        get_jets(ev);
        clean_jets();
        get_MET(ev);

        // Apply analysis filters
        // check_cuts(ev); // fills selection_pass_flags
    }

    bool has_trigger(string name)
    {
        return std::find_if(this->triggers.cbegin(), this->triggers.cend(),
                            [name](zHLT trig) {
                                return trig.HLTName == name && trig.HTLDecision == 1 && trig.HTLPrescale >= 1;
                            }) != this->triggers.cend();
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

        // Handle to the electron impact parameter
        edm::Handle<std::vector<float> > electronDxy;
        event.getByLabel(std::string("electrons:elDxy"), electronDxy);
        // Handle to the electron veto
        edm::Handle<std::vector<float> > electronisVeto;
        event.getByLabel(std::string("electrons:elvidVeto"), electronisVeto);
        // Handle to the electron missing hits
        edm::Handle<std::vector<float> > electronmissHits;
        event.getByLabel(std::string("electrons:elmissHits"), electronmissHits);
        // Handle to the electron isolation of cone of radius 0.03
        edm::Handle<std::vector<float> > electronIso03;
        event.getByLabel(std::string("electrons:elIso03"), electronIso03);
        // Handle to the electron charge
        edm::Handle<std::vector<float> > electronCharge;
        event.getByLabel(std::string("electrons:elCharge"), electronCharge);

        // Handle to the electron electronSCeta
        edm::Handle<std::vector<float> > electronSCeta;
        event.getByLabel(std::string("electrons:elSCEta"), electronSCeta);
        // Handle to the electron electronfullSigmaEtaEta
        edm::Handle<std::vector<float> > electronfullSigmaEtaEta;
        event.getByLabel(std::string("electrons:elfull5x5siee"), electronfullSigmaEtaEta);
        // Handle to the electron electrondEtaIn
        edm::Handle<std::vector<float> > electrondEtaIn;
        event.getByLabel(std::string("electrons:eldEtaIn"), electrondEtaIn);
        // Handle to the electron electrondPhiIn
        edm::Handle<std::vector<float> > electrondPhiIn;
        event.getByLabel(std::string("electrons:eldPhiIn"), electrondPhiIn);
        // Handle to the electron electronHOverE
        edm::Handle<std::vector<float> > electronHOverE;
        event.getByLabel(std::string("electrons:elHoE"), electronHOverE);
        // Handle to the electron electronooEmooP
        edm::Handle<std::vector<float> > electronooEmooP;
        event.getByLabel(std::string("electrons:elooEmooP"), electronooEmooP);
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

            zElectron thisElectron = zElectron(v, electronCharge->at(i), electronSCeta->at(i), electronDxy->at(i),
                                               electronD0->at(i), electronDz->at(i), electronfullSigmaEtaEta->at(i),
                                               electrondEtaIn->at(i), electrondPhiIn->at(i), electronHOverE->at(i),
                                               electronooEmooP->at(i), electronIso03->at(i),
                                               (int) electronmissHits->at(i), electronisVeto->at(i) == 1);

            this->electrons.push_back(thisElectron);
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

            zMuon thisMuon(v, muonCharge->at(i), muonIso04->at(i), muonTight->at(i) == 1);
            this->muons.push_back(thisMuon);
        }
    }

    void get_jets(edm::EventBase const &event)
    {
        // Handle to the jet CSV
        edm::Handle<std::vector<float> > jetCSV;
        event.getByLabel(std::string("jetsAK4CHS:jetAK4CHSCSVv2"), jetCSV);
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

        for (size_t i = 0; i < jetCSV->size(); i++)
        {
            TLorentzVector v;
            v.SetPtEtaPhiE(jetPt->at(i), jetEta->at(i), jetPhi->at(i), jetEn->at(i));

            zJet this_jet = zJet(v, (int) jetCharge->at(i), jetCSV->at(i), jetRapidity->at(i), jetArea->at(i),
                                 jetneutralHadronEnergyFrac->at(i), jetneutralEmEnergyFrac->at(i),
                                 jetchargedHadronEnergyFrac->at(i), jetchargedEmEnergyFrac->at(i),
                                 jetNumConstituents->at(i), jetchargedMultiplicity->at(i),
                                 jetneutralMultiplicity->at(i));
            this->jets.push_back(this_jet);
        }
    }

    void clean_jets()
    {
        for (size_t i = 0; i < this->jets.size(); i++)
        {
            zJet this_jet = this->jets.at(i);
            bool jet_flag = true;

            for (size_t j = 0; j < this->electrons.size(); j++)
            {
                if (this_jet.deltaR(this->electrons.at(j)) <= 0.4)
                {
                    jet_flag = false;
                    break;
                }
            }

            if (!jet_flag)
            {
                this->is_clean_jet.push_back(false);
                break;
            }

            for (size_t j = 0; j < this->muons.size(); j++)
            {
                if (this_jet.deltaR(this->muons.at(j)) <= 0.4)
                {
                    jet_flag = false;
                    break;
                }
            }

            this->is_clean_jet.push_back(jet_flag);
        }
    }

    void get_MET(edm::EventBase const &event)
    {
        edm::Handle<std::vector<float> > MetPx;
        event.getByLabel(std::string("metFull:metFullPx"), MetPx);

        edm::Handle<std::vector<float> > MetPy;
        event.getByLabel(std::string("metFull:metFullPy"), MetPy);

        this->MET = TLorentzVector(MetPx->at(0), MetPy->at(0), 0, 0);
    }

    void get_vertices(edm::EventBase const &event) {
        edm::Handle<std::vector<float> > vtxZ;
        event.getByLabel(std::string("vertexInfo:z"), vtxZ);
        // Handle to the dof of vertex
        edm::Handle<std::vector<int> > dof;
        event.getByLabel(std::string("vertexInfo:ndof"), dof);
        // Handle to the rho
        edm::Handle<std::vector<float> > rhoo;
        event.getByLabel(std::string("vertexInfo:rho"), rhoo);

        for (size_t i = 0; i < vtxZ->size(); i++) {
            zVertex this_vtx;
            this_vtx.z = vtxZ->at(i);
            this_vtx.dof = dof->at(i);
            this_vtx.RHO = rhoo->at(i);

            this->vertices.push_back(this_vtx);
         }
    }

public:
    size_t NVtx() {
        return this->vertices.size();
    }
    const zVertex& get_vertex(size_t index) {
        return this->vertices.at(index);
    }
};

#endif //INC_13TEV_TW_ANALYSIS_ZEVENT_HH
