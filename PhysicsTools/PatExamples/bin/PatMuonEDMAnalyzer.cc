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
#include <TROOT.h>

using namespace std;

//#define SYNC_EX
//#define SYNC_TW

#include "zEvent.hh"

typedef std::pair<std::vector<zLepton>::iterator, std::vector<zLepton>::iterator> llPair_t;


void MakeBranches(TTree *tree)
{
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

void MakeBDTBranches(TTree *tree)
{
    Double_t *temp = NULL;
    int *temp_ = NULL;
    tree->Branch("ptsys", temp);
    tree->Branch("dpt_ll_metj", temp);
    tree->Branch("MET", temp);
    tree->Branch("dpt_ll_met", temp);
    tree->Branch("pt_lMETj", temp);
    tree->Branch("cll", temp);
    tree->Branch("dpt_l_met", temp);
    tree->Branch("njets", temp_);
    tree->Branch("nbjets", temp_);
    tree->Branch("htlljmet", temp);
    tree->Branch("ptj", temp);
    tree->Branch("pt_over_ht", temp);
    tree->Branch("mlljmet", temp);
    tree->Branch("ptllj", temp);
    tree->Branch("htll_over_ht", temp);
    tree->Branch("ptll", temp);
    tree->Branch("dr_llmetjj", temp);
    tree->Branch("dr_lljj", temp);
    tree->Branch("mlj2", temp);
    tree->Branch("dpt_l_j", temp);
    tree->Branch("mlj", temp);
    tree->Branch("ptl", temp);
    tree->Branch("dr_l_j", temp);
    tree->Branch("ptj2", temp);
    tree->Branch("ml2jj", temp);
    tree->Branch("ml2j1", temp);
    tree->Branch("ml2j2", temp);
    tree->Branch("mc_w_sign", temp);
}

int main(int argc, char *argv[])
{
    // ----------------------------------------------------------------------
    // First Part:
    //
    //  * enable the AutoLibraryLoader
    //  * book the histograms of interest
    //  * open the input file
    // ----------------------------------------------------------------------
#ifndef SYNC_EX
    if (argc < 2)
    {
        cout << "Usage: "<< argv[0] << "file1 file2 ..."<< endl;
        return 0;
    }
#endif
    // load framework libraries
    gSystem->Load("libFWCoreFWLite");
//    gROOT->ProcessLine("#include <vector>");
    FWLiteEnabler::enable();

    char fname[160];

    //vector<zEvent> events;
//    double cNetEvWt = 0;
    Long64_t counter[4][8] = {{0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0, 0, 0, 0}};

    cout << boolalpha;
#ifdef TW_SYNC
    fstream ee_lep_evid("ee_ihepru_ttbar_lepsel.txt", ee_lep_evid.out | ee_lep_evid.trunc);
    fstream emu_lep_evid("emu_ihepru_ttbar_lepsel.txt", emu_lep_evid.out | emu_lep_evid.trunc);
    fstream mumu_lep_evid("mumu_ihepru_ttbar_lepsel.txt", mumu_lep_evid.out | mumu_lep_evid.trunc);

    fstream ee_jet1_evid("ee_ihepru_ttbar_jetsel1.txt", ee_lep_evid.out | ee_lep_evid.trunc);
    fstream emu_jet1_evid("emu_ihepru_ttbar_jetsel1.txt", emu_lep_evid.out | emu_lep_evid.trunc);
    fstream mumu_jet1_evid("mumu_ihepru_ttbar_jetsel1.txt", mumu_lep_evid.out | mumu_lep_evid.trunc);

    fstream ee_jet2_evid("ee_ihepru_ttbar_jetsel2.txt", ee_lep_evid.out | ee_lep_evid.trunc);
    fstream emu_jet2_evid("emu_ihepru_ttbar_jetsel2.txt", emu_lep_evid.out | emu_lep_evid.trunc);
    fstream mumu_jet2_evid("mumu_ihepru_ttbar_jetsel2.txt", mumu_lep_evid.out | mumu_lep_evid.trunc);

    fstream debug("debug.log", debug.out | debug.trunc);
#endif

#ifdef SYNC_EX
#ifdef SYNC_TT
    fwlite::TFileService fs = fwlite::TFileService("ttbar_sync.root");
#elif defined(SYNC_TW)
    fwlite::TFileService fs = fwlite::TFileService("tW_sync.root");
#endif
#else
    fwlite::TFileService fs = fwlite::TFileService("tW_main.root");
#endif

    TTree *tW_tree = fs.make<TTree>("tW", "tW");
    TTree *bdt_tree = fs.make<TTree>("BDT", "BDT");
    MakeBranches(tW_tree);
    MakeBDTBranches(bdt_tree);

    vector<zLepton> selectedLeptons;
    vector<zJet> selectedJets;
    vector<zJet> selectedBJets;


    vector<llPair_t> llPairs;

    try
    {
#ifdef SYNC_EX
#ifdef SYNC_TT
        sprintf(fname,
                "root://eoscms//eos/cms//store/group/phys_top/Priyanka/ttBar/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttBar/170302_080613/0000/B2GEDMNtuple_7.root");
#elif defined(SYNC_TW)
        sprintf(fname,
                "/afs/cern.ch/work/r/razumov/Reza/CMSSW_8_0_26_patch1/src/tW/PatMuonEDMAnalyzer/sync_tW_in.root");
#endif
#else
        strncpy(fname, argv[1], 159);
        fname[159] = '\0';
#endif
        cout << "File Name:" << fname << endl;

        TFile *inFile = TFile::Open(fname);

        if (!inFile)
        {
            cout << "Input Root File not opened for processing " << endl;
            return 1;
        }

        TTree *rTree = (TTree *) inFile->Get("IIHEAnalysis");
        zEvent *event = new zEvent(rTree, false);

        for (Long64_t evtID = 0; evtID < rTree->GetEntriesFast(); evtID++)
        {
            counter[EE][0]++;
            counter[EMu][0]++;
            counter[MuMu][0]++;

//            cout << "evtID=" << evtID<< endl;
            event->read_event(rTree, evtID);

#ifndef SYNC_EX
/*
            cNetEvWt += event->getWeight();
            if (event->has_trigger("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v9"))
            {
                cTHltEv++;
                event->add_flag(make_pair("trigger", true));
            }
            else
            {
                event->add_flag(make_pair("trigger", false));
                continue;
            }
*/
#endif
            selectedBJets.clear();
            selectedLeptons.clear();
            selectedJets.clear();

            // event.write(fs);

            llPairs.clear();
            int event_tag = -1;

            TLorentzVector ll;

/*
            if (event->getEvID() == 927179)
            {
                cout << "In analyzer" << endl;
                for (auto i = event->getLeptons().begin(); i != event->getLeptons().end(); i++)
                {
                    cout << *i << endl;
                }

                break;
            }
*/
            // in_gap will return false for muons; is_iso will return true for electrons
            copy_if(event->getLeptons().begin(), event->getLeptons().end(), back_inserter(selectedLeptons),
                    [](const zLepton &part) { return part.is_selected(); });

            copy_if(event->getJets().begin(), event->getJets().end(), back_inserter(selectedJets),
                    [](const zJet &jet) { return jet.is_selected(); });

            copy_if(selectedJets.begin(), selectedJets.end(), back_inserter(selectedBJets),
                    [](const zJet &jet) { return jet.is_bjet(); });

            cout << "=== BEGIN EVENT " << evtID << " (" << event->getEvID() << ")" << endl;
            cout << "Raw: # leptons = " << event->getLeptons().size() << ", # jets = " << event->getJets().size()
                 << endl;
            cout << "Sel: # leptons = " << selectedLeptons.size() << ", # jets = " << selectedJets.size();
            cout << ", # bjets = " << selectedBJets.size() << endl;


            if (!event->isMETok())
            {
                event->fill_tree(tW_tree);
                cout << "- Reject step 0.1" << endl;
                continue;
            }
            else
                event->add_flag("pass_met_filters", true);

            if (selectedLeptons.size() >= 2)
            {
                cout << "=== New ll candidate " << evtID << "(" << event->getEvID() << ") ===" << endl;
                event->add_flag("is_ll", true);

                auto leading_l = selectedLeptons.at(0);
                auto other_l = selectedLeptons.at(1);
                ll = leading_l + other_l;

                cout << "Leptons: " << selectedLeptons.size() << endl;
                cout << "Leading: " << leading_l << endl;
                cout << "Second: " << other_l << endl;
                cout << "Dileption: " << "(Pt, Eta, Phi, E) = (" << ll.Pt() << ", " << ll.Eta() << ", " << ll.Phi()
                     << ", "
                     << ll.E() << "), M = " << ll.Mag() << endl;

                if (!leading_l.is_samesign(other_l) && leading_l.Pt() > 25 && ll.Mag() > 20)
                {
                    cout << "+ Accept step 1" << endl;
                    if (leading_l.is_muon() && other_l.is_muon())
                        event_tag = MuMu;
                    else if (!(leading_l.is_muon() || other_l.is_muon()))
                        event_tag = EE;
                    else
                        event_tag = EMu;
                }
            }
            else
            {
                event->add_flag("no_ll", true);
                event->fill_tree(tW_tree);
                cout << "- Reject step 1" << endl;
                continue;
            }

            counter[event_tag][1]++;

            if (event_tag == EE)
            {
#ifdef SYNC_EX
                ee_lep_evid << event->getEvID() << endl;
#endif
                event->add_flag("EE", true);
            }
            else if (event_tag == EMu)
            {
#ifdef SYNC_EX
                emu_lep_evid << event->getEvID() << endl;
#endif
                event->add_flag("EMu", true);
            }
            else if (event_tag == MuMu)
            {
#ifdef SYNC_EX
                mumu_lep_evid << event->getEvID() << endl;
#endif
                event->add_flag("MuMu", true);
            }
            else
            {
                event->add_flag("not_dilepton", true);
                event->fill_tree(tW_tree);
                continue;
            }

            bool massFlag = false;

            if (event_tag != EMu)
                massFlag = (ll.Mag() >= 76 && ll.Mag() <= 106);

            if (massFlag)
            {
                cout << "- Reject step 2" << endl;
                event->fill_tree(tW_tree);
                continue;
            }
            else
            {
                event->add_flag("pass_dilepton_mass", true);
                cout << "+ Accept step 2" << endl;
            }


            counter[event_tag][2]++;

            if (event_tag != EMu)
            {
                if (event->getMET().Pt() <= 40)
                {
                    event->fill_tree(tW_tree);
                    cout << "- Reject step 3" << endl;
                    continue;
                }
            }
            event->add_flag("pass_met", true);
            cout << "+ Accept step 3" << endl;

            counter[event_tag][3]++;

            if (selectedJets.size() >= 2)
            {
                counter[event_tag][4]++;
                cout << "+ Accept step 4.1.1" << endl;
                event->add_flag("2+j", true);
                if (selectedBJets.size() >= 1)
                {
                    counter[event_tag][5]++;
                    cout << "+ Accept step 4.1.2" << endl;
                    event->add_flag("2+j1+t", true);
                    event->fill_tree_2(bdt_tree);
                }

#ifdef SYNC_EX
                if (event_tag == EE)
                {
                    ee_jet2_evid << event->getEvID() << endl;
                }
                else if (event_tag == EMu)
                {
                    emu_jet2_evid << event->getEvID() << endl;
                }
                else
                {
                    mumu_jet2_evid << event->getEvID() << endl;
                }
#endif
            }

            if (selectedJets.size() == 1)
            {
                counter[event_tag][6]++;
                event->add_flag("1j", true);
                cout << "+ Accept step 4.2.1" << endl;
                if (selectedBJets.size() == 1)
                {
                    counter[event_tag][7]++;
                    event->add_flag("1j1t", true);
                    cout << "+ Accept step 4.2.2" << endl;
                    event->fill_tree_2(bdt_tree);
                }
#ifdef SYNC_EX
                if (event_tag == EE)
                {
                    ee_jet1_evid << event->getEvID() << endl;
                }
                else if (event_tag == EMu)
                {
                    emu_jet1_evid << event->getEvID() << endl;
                }
                else
                {
                    mumu_jet1_evid << event->getEvID() << endl;
                }
#endif
            }

            event->fill_tree(tW_tree);
        }
    }
    catch (exception &e)
    {
        cout << "error:" << e.what() << endl;
        cout << "File Name:" << fname << " is not valid!" << endl << endl << endl;
    }

#ifdef SYNC_EX
    ee_lep_evid.close();
    emu_lep_evid.close();
    mumu_lep_evid.close();

    ee_jet1_evid.close();
    emu_jet1_evid.close();
    mumu_jet1_evid.close();

    ee_jet2_evid.close();
    emu_jet2_evid.close();
    mumu_jet2_evid.close();
#endif
    for (int i = 0; i < 4; i++)
    {
        if (i == 2)
            continue;

        cout << "=== ";
        switch (i)
        {
            case EE:
                cout << "EE";
                break;
            case EMu:
                cout << "EMu";
                break;
            case MuMu:
                cout << "MuMu";
                break;
            default:
                cout << "None";
        }
        cout << " channel ===" << endl;
        for (int j = 0; j < 8; j++)
            cout << "Step " << j << " events " << counter[i][j] << endl;
    }

    return 0;
}

                        


