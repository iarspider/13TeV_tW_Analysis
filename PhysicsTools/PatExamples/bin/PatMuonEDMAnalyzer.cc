#include <memory>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

#include <TTree.h>
#include <TFile.h>
#include <TSystem.h>
#include <FWCore/FWLite/interface/FWLiteEnabler.h>
#include <PhysicsTools/FWLite/interface/TFileService.h>
#include <PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h>
#include <TROOT.h>

#include "argparse.hpp"

using namespace std;

#include "zEvent.hh"

typedef std::pair<std::vector<zLepton>::iterator, std::vector<zLepton>::iterator> llPair_t;


void MakeBranches(TTree *tree) {
    std::vector<TLorentzVector> JetV, LepV;
    std::vector<float> JetCh, LepCh, LepDxy, LepDz, LepEtaSC;
    std::vector<bool> JetB, JetClean, JetSel, JetLoose, LepSel, LepIso, LepTight, LepMuon;
    std::vector<int> LepWhere;
    std::vector<string> flags;
    TLorentzVector MET;

    tree->Branch("LeptonVec", &LepV);
    tree->Branch("LeptonCharge", &LepCh);
    tree->Branch("LeptonSelected", &LepSel);
    tree->Branch("LeptonWhere", &LepWhere);
    tree->Branch("LeptonIsIso", &LepIso);
    tree->Branch("LeptonIsTight", &LepTight);
    tree->Branch("LeptonIsMuon", &LepMuon);
    tree->Branch("LeptonDxy", &LepDxy);
    tree->Branch("LeptonDz", &LepDz);
    tree->Branch("LeptonEtaSC", &LepEtaSC);
    tree->Branch("LeptonIso", &LepDxy);

    tree->Branch("JetVec", &JetV);
    tree->Branch("JetCharge", &JetCh);
    tree->Branch("JetSelected", &JetSel);
    tree->Branch("JetIsClean", &JetClean);
    tree->Branch("JetIsBJet", &JetB);
    tree->Branch("JetIsLoose", &JetLoose);

    bool isMetOk, BadCCF, BadPFM;
    tree->Branch("MetVec", &MET);
    tree->Branch("MetOK", &isMetOk);
    tree->Branch("MetBadCCF", &BadCCF);
    tree->Branch("MetBadPFM", &BadPFM);

    ULong64_t *evid = NULL;
    tree->Branch("Flags", &flags);
    tree->Branch("eventID", evid);
}

void MakeBDTBranches(TTree *tree) {
    Double_t *temp_double = NULL;
    int *temp_int = NULL;
    Float_t *temp_float = NULL;
    tree->Branch("ptsys", temp_double);
    tree->Branch("dpt_ll_metj", temp_double);
    tree->Branch("MET", temp_double);
    tree->Branch("dpt_ll_met", temp_double);
    tree->Branch("pt_lMETj", temp_double);
    tree->Branch("cll", temp_double);
    tree->Branch("dpt_l_met", temp_double);
    tree->Branch("njets", temp_int);
    tree->Branch("nbjets", temp_int);
    tree->Branch("htlljmet", temp_double);
    tree->Branch("ptj", temp_double);
    tree->Branch("pt_over_ht", temp_double);
    tree->Branch("mlljmet", temp_double);
    tree->Branch("ptllj", temp_double);
    tree->Branch("htll_over_ht", temp_double);
    tree->Branch("ptll", temp_double);
    tree->Branch("dr_llmetjj", temp_double);
    tree->Branch("dr_lljj", temp_double);
    tree->Branch("mlj2", temp_double);
    tree->Branch("dpt_l_j", temp_double);
    tree->Branch("mlj", temp_double);
    tree->Branch("ptll", temp_double);
    tree->Branch("ptsl", temp_double);
    tree->Branch("dr_l_j", temp_double);
    tree->Branch("ptj2", temp_double);
    tree->Branch("ml2jj", temp_double);
    tree->Branch("ml2j1", temp_double);
    tree->Branch("ml2j2", temp_double);
    tree->Branch("mc_w_sign", temp_float);
    tree->Branch("j1csv", temp_float);
    tree->Branch("mll", temp_double);
    tree->Branch("channel", temp_int);
}

int main(int argc, const char *argv[]) {
    // ----------------------------------------------------------------------
    // First Part:
    //
    //  * enable the AutoLibraryLoader
    //  * book the histograms of interest
    //  * open the input file
    // ----------------------------------------------------------------------

    ArgumentParser parser;
    parser.addArgument("--epoch", "-e", 1, true);
    parser.addArgument("--input", "-i", 1, true);
    parser.addArgument("--output", "-o", 1, true);

    parser.parse(static_cast<size_t>(argc), argv);
    // load framework libraries
    gSystem->Load("libFWCoreFWLite");
    FWLiteEnabler::enable();

    char fname[256];

    Float_t counter[4][14] = {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};

    cout << boolalpha;

    string infile = parser.retrieve<string>("input");
    string outfile = parser.retrieve<string>("output");

    fwlite::TFileService fs = fwlite::TFileService(outfile);
    cout << "Input file Name:" << infile << endl;

    TTree *tW_tree = NULL;
    TTree *bdt_tree = fs.make<TTree>("BDT", "BDT");
    MakeBDTBranches(bdt_tree);

    vector<zLepton> selectedLeptons;
    vector<zJet> selectedJets;
    vector<zJet> selectedBJets;

    TFile *inFile = TFile::Open(infile.c_str());

    if (!inFile) {
        cout << "Input Root File not opened for processing " << endl;
        return 1;
    }

    TTree *rTree = (TTree *) inFile->Get("IIHEAnalysis");
    zEvent *event = new zEvent(rTree, parser.retrieve<string>("epoch"));
    vector<float> puMC, puData;

    for (Long64_t evtID = 0; evtID < rTree->GetEntriesFast(); evtID++) {
        counter[EE][0]++;
        counter[EMu][0]++;
        counter[MuMu][0]++;

        event->read_event(rTree, evtID);

        if (!event->getIsData() && (event->get_decay_mode() == -1)) {
            cout << "WARN: Too many leptons in LHE data!" << endl;
            counter[EMu][10]++;
        }
        selectedBJets.clear();
        selectedLeptons.clear();
        selectedJets.clear();

        int event_tag = -1;

        TLorentzVector ll;

        // in_gap will return false for muons; is_iso will return true for electrons
        copy_if(event->getLeptons().begin(), event->getLeptons().end(), back_inserter(selectedLeptons),
                [](const zLepton &part) { return part.is_selected(); });

        copy_if(event->getJets().begin(), event->getJets().end(), back_inserter(selectedJets),
                [](const zJet &jet) { return jet.is_selected(); });

        copy_if(selectedJets.begin(), selectedJets.end(), back_inserter(selectedBJets),
                [](const zJet &jet) { return jet.is_bjet(); });


        if (!event->isMETok()) {
            event->fill_dump_tree(tW_tree);
            continue;
        }
        else
            event->add_flag("pass_met_filters", true);

        if (selectedLeptons.size() >= 2) {
            event->add_flag("is_ll", true);

            auto leading_l = selectedLeptons.at(0);
            auto other_l = selectedLeptons.at(1);
            ll = leading_l + other_l;

            if (!leading_l.is_samesign(other_l) && leading_l.Pt() > 25 && ll.Mag() > 20) {
                if (leading_l.is_muon() && other_l.is_muon())
                    event_tag = MuMu;
                else if (!(leading_l.is_muon() || other_l.is_muon()))
                    event_tag = EE;
                else
                    event_tag = EMu;
            }
        }
        else {
            event->add_flag("no_ll", true);
            event->fill_dump_tree(tW_tree);
            continue;
        }

        counter[event_tag][1]++;

        if (event_tag == EE) {
            event->add_flag("EE", true);
        }
        else if (event_tag == EMu) {
            event->add_flag("EMu", true);
        }
        else if (event_tag == MuMu) {
            event->add_flag("MuMu", true);
        }
        else {
            event->add_flag("not_dilepton", true);
            event->fill_dump_tree(tW_tree);
            continue;
        }

        bool massFlag;

        if (event_tag != EMu)
            massFlag = (ll.Mag() >= 76 && ll.Mag() <= 106);
        else
            massFlag = false;

        if (massFlag) {
            event->fill_dump_tree(tW_tree);
            continue;
        }
        else {
            event->add_flag("pass_dilepton_mass", true);
        }

        counter[event_tag][2]++;
        if (event_tag != EMu) {
            if (event->getMET().Pt() <= 40) {
                event->fill_dump_tree(tW_tree);
                continue;
            }
        }
        event->add_flag("pass_met", true);

        counter[event_tag][3]++;

        if (!event->isTrgOk(event_tag)) {
            event->fill_dump_tree(tW_tree);
            continue;
        }
        event->add_flag("pass_trigger", true);
        counter[event_tag][4] += event->get_mc_w();

        if (selectedJets.size() == 0) {
            counter[event_tag][5] += event->get_mc_w();
            event->fill_BDT_tree(bdt_tree);
            if (event->get_decay_mode() == DECAY_LJ) {
                counter[event_tag][8] += event->get_mc_w();
            }
            if (event->get_decay_mode() == DECAY_JJ) {
                counter[event_tag][12] += event->get_mc_w();
            }
        }

        if ((selectedJets.size() == 1) && (selectedBJets.size() == 0)) {
            counter[event_tag][6] += event->get_mc_w();
            event->fill_BDT_tree(bdt_tree);
            if (event->get_decay_mode() == DECAY_LJ) {
                counter[event_tag][9] += event->get_mc_w();
            }
            if (event->get_decay_mode() == DECAY_JJ) {
                counter[event_tag][13] += event->get_mc_w();
            }
        }

        if (event->get_decay_mode() == DECAY_LJ) {
            counter[event_tag][7] += event->get_mc_w();
        }

        if (event->get_decay_mode() == DECAY_LL) {
            counter[event_tag][11] += event->get_mc_w();
        }
    }

    int i = EMu;
    cout << "=== EMu channel ===" << endl;
    for (int j = 0; j < 14; j++)
        cout << "Counter " << j << " value " << counter[i][j] << endl;

    return 0;
}

                        


