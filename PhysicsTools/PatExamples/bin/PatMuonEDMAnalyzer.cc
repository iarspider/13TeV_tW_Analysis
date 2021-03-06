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

#define SYNC_EX
#define SYNC_TW

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
    tree->Branch("LeptonIso", &LepCh);


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
    Double_t * temp = NULL;
    Int_t* temp_ = NULL;
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

    // load framework libraries
    gSystem->Load("libFWCoreFWLite");
//    gROOT->ProcessLine("#include <vector>");
    FWLiteEnabler::enable();

    char fname[160];

    vector<zEvent> events;
//    double cNetEvWt = 0;
    ulong counter[4][8] = {{0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0}};

    cout << boolalpha;

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

#ifdef SYNC_EX
#ifdef SYNC_TT
    fwlite::TFileService fs = fwlite::TFileService("ttbar_sync.root");
    TTree *tree = fs.make<TTree>("tW", "tW");
    MakeBranches(tree);
#elif defined(SYNC_TW)
    fwlite::TFileService fs = fwlite::TFileService("tW_sync.root");
    TTree *tree = fs.make<TTree>("tW", "tW");
    TTree *bdt_tree = fs.make<TTree>("BDT", "BDT");
    MakeBranches(tree);
    MakeBDTBranches(bdt_tree);
#endif
#else
    fwlite::TFileService fs = fwlite::TFileService("tW_main.root");
    TTree *tree = fs.make<TTree>("tW", "tW");
    MakeBranches(tree);
#endif


#ifndef SYNC_EX
    int n = 2;
    //int s=atoi(argv[2]);
    //int n=atoi(argv[3]);
    //int jobid=atoi(argv[1]);

    // for(int i=s;i<=n;i++)
    for (int i = 1; i <= n; i++)
    {
        //char constName[]="/afs/cern.ch/work/p/ppriyank/private/816Code/CMSSW_8_0_16/src/B2GEDMNtuple_";
        ///////datas
        //char constName[]="root://eoscms//eos/cms/store/group/phys_top/Priyanka/2016B/MuonEG/crab_2016B/170209_115443/0000/B2GEDMNtuple_";
        //char constName[]="root://eoscms//eos/cms/store/group/phys_top/Priyanka/2016C/MuonEG/crab_2016C/170209_120857/0000/B2GEDMNtuple_";
        //char constName[]="root://eoscms//eos/cms/store/group/phys_top/Priyanka/2016D/MuonEG/crab_2016D/170209_121828/0000/B2GEDMNtuple_";
        //char constName[]="root://eoscms//eos/cms/store/group/phys_top/Priyanka/2016E/MuonEG/crab_2016E/170209_130319/0000/B2GEDMNtuple_";
        //  char constName[]="root://eoscms//eos/cms/store/group/phys_top/Priyanka/2016F/MuonEG/crab_2016F/170209_134155/0000/B2GEDMNtuple_";
        //  char constName[]="root://eoscms//eos/cms/store/group/phys_top/Priyanka/2016G/MuonEG/crab_2016G/170209_142904/0000/B2GEDMNtuple_";
        //  char constName[]="root://eoscms//eos/cms/store/group/phys_top/Priyanka/2016H2/MuonEG/crab_2016H2/170209_180745/0000/B2GEDMNtuple_";
        //  char constName[]="root://eoscms//eos/cms/store/group/phys_top/Priyanka/2016H3/MuonEG/crab_2016H3/170209_182434/0000/B2GEDMNtuple_";

        ///////////mcs
        char constName[] = "root://eoscms//eos/cms/store/group/phys_top/Priyanka/tWtop/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_tWtop/170301_124858/0000/B2GEDMNtuple_";
        //char constName[]="root://eoscms//eos/cms/store/group/phys_top/Priyanka/tWtop/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_tWtop/170207_130603/0000/B2GEDMNtuple_";
        //char constName[]="root://eoscms//eos/cms/store/group/phys_top/Priyanka/tWantiTop/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_tWantiTop/170207_130746/0000/B2GEDMNtuple_";
        //char constName[]="root://eoscms//eos/cms/store/group/phys_top/Priyanka/ttBar/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttBar/170207_130836/0000/B2GEDMNtuple_";
        //char constName[]="root://eoscms//eos/cms/store/group/phys_top/Priyanka/DY10to50/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DY10to50/170207_131035/0000/B2GEDMNtuple_";
        //char constName[]="root://eoscms//eos/cms/store/group/phys_top/Priyanka/DYabove50/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYabove50/170207_130929/0000/B2GEDMNtuple_";
        //char constName[]="root://eoscms//eos/cms/store/group/phys_top/Priyanka/WPjets/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WPjets/170207_131406/0000/B2GEDMNtuple_";
        // char constName[]="root://eoscms//eos/cms/store/group/phys_top/Priyanka/WW/WW_TuneCUETP8M1_13TeV-pythia8/crab_WW/170207_131132/0000/B2GEDMNtuple_";
        // char constName[]="root://eoscms//eos/cms/store/group/phys_top/Priyanka/ZZ/ZZ_TuneCUETP8M1_13TeV-pythia8/crab_ZZ/170207_131315/0000/B2GEDMNtuple_";
        // char constName[]="root://eoscms//eos/cms/store/group/phys_top/Priyanka/WZ/WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/170207_131224/0000/B2GEDMNtuple_";
        //char constName[]="root://eoscms//eos/cms/store/group/phys_top/Priyanka/tChnlTop/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/crab_tChnlTop/170207_131459/0000/B2GEDMNtuple_";
        // char constName[]="root://eoscms//eos/cms/store/group/phys_top/Priyanka/tChnlAntiTop/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/crab_tChnlAntiTop/170208_055151/0000/B2GEDMNtuple_";
        // char constName[]="root://eoscms//eos/cms/store/group/phys_top/Priyanka/sChnl/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/crab_sChnl/170207_131653/0000/B2GEDMNtuple_";
        try
        {
            sprintf(fname, "%s%d.root", constName, i);
#else
    try
    {
#endif

#ifdef SYNC_TT
        sprintf(fname,
                "root://eoscms//eos/cms//store/group/phys_top/Priyanka/ttBar/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_ttBar/170302_080613/0000/B2GEDMNtuple_7.root");
#elif defined(SYNC_TW)
        sprintf(fname,
                "root://eoscms//eos/cms//store/group/phys_top/Priyanka/sync/CRAB_UserFiles/crab_sync/170404_071616/0000/B2GEDMNtuple_1.root");
#endif
        cout << "File Name:" << fname << endl;

        TFile *inFile = TFile::Open(fname);

        if (!inFile)
        {
            cout << "Input Root File not opened for processing " << endl;
            return 1;
        }

        fwlite::Event ev(inFile);
        int evtID = 0;

        edm::Handle<std::vector<string>> METTriggerNameTree;

        for (ev.toBegin(); !ev.atEnd(); ++ev, evtID++)
        {
            edm::EventBase const &event = ev;
            events.push_back(zEvent(event));
        }
    }
    catch (exception &e)
    {
        cout << "error:" << e.what() << endl;
        cout << "File Name:" << fname << " is not valid!" << endl << endl << endl;
    }
#ifndef SYNC_EX
    }  // for i
#endif

    counter[EE][0] = events.size();
    counter[EMu][0] = events.size();
    counter[MuMu][0] = events.size();
    vector<zLepton> selectedLeptons;
    vector<zJet> selectedJets;
    vector<zJet> selectedBJets;

    vector<llPair_t> llPairs;

    int iEvent = 0;

    for (auto event = events.begin(); event != events.end(); event++, iEvent++)
    {
#ifndef SYNC_EX
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

        double Vtx_Cut_z = 24.0;
        int Vtx_Cut_ndof = 4;
        double Vtx_Cut_rho = 2.0;
        bool has_pv = false;

        if (event->NVtx() > 0)
        {
            for (size_t iVtx = 0; iVtx < event->NVtx(); iVtx++)
            {
                const zVertex vtx = event->get_vertex(iVtx);

                if ((fabs(vtx.z) < Vtx_Cut_z) && (vtx.dof > Vtx_Cut_ndof) && (vtx.RHO < Vtx_Cut_rho))
                {
                    event->add_flag(make_pair("vertex", true));
                    cVertexEv++;
                    has_pv = true;
                    break;
                }
            }
        }

        if (!has_pv)
        {
            event->add_flag(make_pair("vertex", false));
            continue;
        }
#endif
        selectedBJets.clear();
        selectedLeptons.clear();
        selectedJets.clear();

        // event.write(fs);

        llPairs.clear();
        int event_tag = -1;

        TLorentzVector ll;

        // in_gap will return false for muons; is_iso will return true for electrons
        copy_if(event->getLeptons().begin(), event->getLeptons().end(), back_inserter(selectedLeptons),
                [](const zLepton &part) { return part.is_selected(); });

        copy_if(event->getJets().begin(), event->getJets().end(), back_inserter(selectedJets),
                [](const zJet &jet) { return jet.is_selected(); });

        copy_if(selectedJets.begin(), selectedJets.end(), back_inserter(selectedBJets),
                [](const zJet &jet) { return jet.is_bjet(); });

        cout << "=== BEGIN EVENT " << iEvent << " (" << event->getEvID() << ")" << endl;
        cout << "Raw: # leptons = " << event->getLeptons().size() << ", # jets = " << event->getJets().size() << endl;
        cout << "Sel: # leptons = " << selectedLeptons.size() << ", # jets = " << selectedJets.size();
        cout << ", # bjets = " << selectedBJets.size() << endl;


        if (!event->isMETok())
        {
            event->fill_tree(tree);
            cout << "- Reject step 0.1" << endl;
            continue;
        }
        else
            event->add_flag("pass_met_filters", true);

/*
        if (selectedElectrons.size() >= 2)
        {
            cout << "=== New ee candidate " << iEvent << "(" << event->getEvID() << ") ===" << endl;
            cout << "Selected electrons: " << selectedElectrons.size() << ", muons " << selectedMuons.size()
                 << ", jets "
                 << selectedJets.size() << ", bjets " << selectedBJets.size() << endl;

            event_tagged = true;
        }

        if (selectedMuons.size() >= 2)
        {
            cout << "=== New mumu candidate " << iEvent << "(" << event->getEvID() << ") ===" << endl;
            cout << "Selected electrons: " << selectedElectrons.size() << ", muons " << selectedMuons.size()
                 << ", jets "
                 << selectedJets.size() << ", bjets " << selectedBJets.size() << endl;

            event_tagged = true;
        }

        if (selectedElectrons.size() >= 1 && selectedMuons.size() >= 1)
        {
            cout << "=== New emu candidate " << iEvent << "(" << event->getEvID() << ") ===" << endl;
            cout << "Selected electrons: " << selectedElectrons.size() << ", muons " << selectedMuons.size()
                 << ", jets "
                 << selectedJets.size() << ", bjets " << selectedBJets.size() << endl;
            event_tagged = true;
        }

        if (!event_tagged)
            continue;

        // ee
        size_t ie1 = 0, ie2 = 0, imu1 = 0, imu2 = 0;
        if (selectedElectrons.size() > 0)
        {
            for (auto e1 = selectedElectrons.begin(); e1 != selectedElectrons.end(); e1++, ie1++)
            {
                cout << "Lepton candidate #1 (electron #" << ie1 << "/" << selectedElectrons.size() << "):" << endl;
                cout << "\t" << *e1 << endl;

                if (e1->Pt() > 25.)
                {
                    ie2 = ie1 + 1;
                    for (auto e2 = e1 + 1; e2 != selectedElectrons.end(); e2++, ie2++)
                    {
                        cout << "Lepton candidate #2 (electron #" << ie2 << "/" << selectedElectrons.size() << "):"
                             << endl;
                        cout << "\t" << *e2 << endl;
                        cout << "Dileption object:" << endl;
                        ll = (*e1) + (*e2);
                        cout << "\t(Pt, Eta, Phi, E) = (" << ll.Pt() << ", " << ll.Eta() << ", " << ll.Phi() << ", "
                             << ll.E() << "), M = " << ll.Mag() << endl;


                        if (!e1->is_samesign(*e2) && ll.Mag() > 20.)
                        {
                            cout << "Passed" << endl;
                            eePairs.push_back(make_pair(e1, e2));
                            event_tag = EE;
                        }
                    }
                    if (selectedMuons.size() > 0)
                    {
                        for (auto mu1 = selectedMuons.begin(); mu1 != selectedMuons.end(); mu1++, imu1++)
                        {
                            cout << "Lepton candidate #2 (muon #" << imu1 << "/" << selectedMuons.size() << "):"
                                 << endl;
                            cout << "\t" << *mu1 << endl;
                            cout << "Dileption object:" << endl;
                            ll = (*e1) + (*mu1);
                            cout << "\t(Pt, Eta, Phi, E) = (" << ll.Pt() << ", " << ll.Eta() << ", " << ll.Phi() << ", "
                                 << ll.E() << "), M = " << ll.Mag() << endl;


                            if (!e1->is_samesign(*mu1) && ll.Mag() > 20. && mu1->Pt() < e1->Pt())
                            {
                                emuPairs.push_back(make_pair(e1, mu1));
                                cout << "Passed" << endl;
                                if (event_tag != -1 && event_tag != EMu)
                                {
                                    cout << "! Event " << event->getEvID() << " is already tagged with " << event_tag
                                         << endl;
                                }
                                else
                                {
                                    event_tag = EMu;
                                }
                            }
                        }
                    }
                }
            }
        }

        if (selectedMuons.size() > 0)
        {
            imu1 = 0;

            for (auto mu1 = selectedMuons.begin(); mu1 != selectedMuons.end(); mu1++, imu1++)
            {
                cout << "Lepton candidate #1 (muon #" << ie1 << "/" << selectedMuons.size() << "):" << endl;
                cout << "\t" << *mu1 << endl;

                if (mu1->Pt() > 25.)
                {
                    imu2 = imu1 + 1;
                    for (auto mu2 = mu1 + 1; mu2 != selectedMuons.end(); mu2++)
                    {
                        cout << "Lepton candidate #2 (muon #" << imu2 << "/" << selectedMuons.size() << "):" << endl;
                        cout << "\t" << *mu2 << endl;
                        cout << "Dileption object:" << endl;
                        ll = (*mu1) + (*mu2);
                        cout << "\t(Pt, Eta, Phi, E) = (" << ll.Pt() << ", " << ll.Eta() << ", " << ll.Phi() << ", "
                             << ll.E() << "), M = " << ll.Mag() << endl;


                        if (!mu1->is_samesign(*mu2) && ll.Mag() > 20.)
                        {
                            mumuPairs.push_back(make_pair(mu1, mu2));
                            cout << "Passed" << endl;
                            if (event_tag != -1 && event_tag != MuMu)
                            {
                                cout << "! Event " << event->getEvID() << " is already tagged with " << event_tag
                                     << endl;
                            }
                            else
                            {
                                event_tag = MuMu;
                            }

                        }
                    }

                    ie1 = 0;
                    if (selectedElectrons.size() > 0)
                    {
                        for (auto e1 = selectedElectrons.begin(); e1 != selectedElectrons.end(); e1++, ie1++)
                        {
                            cout << "Lepton candidate #2 (electron #" << ie1 << "/" << selectedElectrons.size() << "):"
                                 << endl;
                            cout << "\t" << *e1 << endl;
                            cout << "Dileption object:" << endl;
                            ll = (*e1) + (*mu1);
                            cout << "\t(Pt, Eta, Phi, E) = (" << ll.Pt() << ", " << ll.Eta() << ", " << ll.Phi() << ", "
                                 << ll.E() << "), M = " << ll.Mag() << endl;


                            if (!e1->is_samesign(*mu1) && ll.Mag() > 20. && mu1->Pt() > e1->Pt())
                            {
                                cout << "Passed" << endl;

                                emuPair_t pair = make_pair(e1, mu1);
                                if (find(emuPairs.begin(), emuPairs.end(), pair) == emuPairs.end())
                                {
                                    if (event_tag != -1 && event_tag != EMu)
                                    {
                                        cout << "! Event " << event->getEvID() << " is already tagged with "
                                             << event_tag << endl;
                                    }
                                    else
                                    {
                                        event_tag = EMu;
                                    }

                                    emuPairs.push_back(pair);
                                }
                            }
                        }
                    }
                }
            }
        }
*/

        if (selectedLeptons.size() >= 2)
        {
            cout << "=== New ll candidate " << iEvent << "(" << event->getEvID() << ") ===" << endl;
            event->add_flag("is_ll", true);

            auto leading_l = selectedLeptons.at(0);
            auto other_l = selectedLeptons.at(1);
            ll = leading_l + other_l;

            cout << "Leptons: " << selectedLeptons.size() << endl;
            cout << "Leading: " << leading_l << endl;
            cout << "Second: " << other_l << endl;
            cout << "Dileption: " << "(Pt, Eta, Phi, E) = (" << ll.Pt() << ", " << ll.Eta() << ", " << ll.Phi() << ", "
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
            event->fill_tree(tree);
            cout << "- Reject step 1" << endl;
            continue;
        }

        counter[event_tag][1]++;

        if (event_tag == EE)
        {
            ee_lep_evid << event->getEvID() << endl;
            event->add_flag("EE", true);
        }
        else if (event_tag == EMu)
        {
            emu_lep_evid << event->getEvID() << endl;
            event->add_flag("EMu", true);
        }
        else if (event_tag == MuMu)
        {
            mumu_lep_evid << event->getEvID() << endl;
            event->add_flag("MuMu", true);
        }
        else
        {
            event->add_flag("not_dilepton", true);
            event->fill_tree(tree);
            continue;
        }

        bool massFlag = false;

        if (event_tag != EMu)
            massFlag = (ll.Mag() >= 76 && ll.Mag() <= 106);
/*
        for (auto ee = eePairs.begin(); ee != eePairs.end(); ee++)
        {
            ll = (*(ee->first)) + (*(ee->second));

            if (ll.Mag() >= 76 && ll.Mag() <= 106)
            {
                massFlag = true;
                break;
            }
        }
*/
        if (massFlag)
        {
            cout << "- Reject step 2" << endl;
            event->add_flag("pass_dilepton_mass", false);
            event->fill_tree(tree);
            continue;
        }
        else
        {
            event->add_flag("pass_dilepton_mass", true);
            cout << "+ Accept step 2" << endl;
        }

        /*
    for (auto mumu = mumuPairs.begin(); mumu != mumuPairs.end(); mumu++)
    {
        ll = (*(mumu->first)) + (*(mumu->second));
        if (ll.Mag() >= 76 && ll.Mag() <= 106)
        {
            massFlag = true;
            break;
        }
    }

    if (massFlag)
        continue;
*/

        counter[event_tag][2]++;

        /*
    if (mumuPairs.size() > 0 || eePairs.size() > 0)
    {
        if (event->getMET().Pt() <= 40)
            continue;
    }

    if (eePairs.size() > 0)
    {
        counter[EE][3]++;
    }

    if (emuPairs.size() > 0)
    {
        counter[EMu][3]++;
    }

    if (mumuPairs.size() > 0)
    {
        counter[MuMu][3]++;
    }
         */

        if (event_tag != EMu)
        {
            if (event->getMET().Pt() <= 40)
            {
                event->add_flag("pass_met_cut", false);
                event->fill_tree(tree);
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
        }

        event->fill_tree(tree);
    }

    ee_lep_evid.close();
    emu_lep_evid.close();
    mumu_lep_evid.close();

    ee_jet1_evid.close();
    emu_jet1_evid.close();
    mumu_jet1_evid.close();

    ee_jet2_evid.close();
    emu_jet2_evid.close();
    mumu_jet2_evid.close();

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
    /*
    cout << "Total number of events processed: " << events.size() << endl;
    cout << "Sum of events Weight : " << cNetEvWt << endl;
    cout << "Total number of events passed HLT cuts: " << cTHltEv << endl;
    cout << "Total number of events passed vertex cuts: " << cVertexEv << endl;

    cout << "Total number of tight electrons: " << cTTE << endl;
    cout << "Total number of events having only one tight electron: " << cEW1e << endl;
    cout << "Total number of events having exactly two tight electrons: " << cEW2e << endl;
    cout << "Total number of tight muons: " << cTTM << endl;
    cout << "Total number of events having only one tight muon: " << cEW1m << endl;
    cout << "Total number of events having exactly two tight muons: " << cEW2m << endl;
    cout << "Total number of Mets passing cuts: " << cMetCut << endl;
    cout << "Total number of Loose jets: " << cTLj << endl;
    cout << "Total number of events having atleast one Loose jet: " << cEwanyLj << endl;
    cout << "Total number of Cleaned jets: " << cToClj << endl;
    cout << "Total number of events having atleast one Cleaned jet: " << cEClj << endl;
    cout << "Total number of events having exactly 1j: " << cE1Lj << endl;
    cout << "Total number of events having exactly 1j1b: " << cE1j1t << endl;
    cout << "Total number of events having exactly 2j: " << cE2Lj << endl;
    cout << "Total number of events having exactly 2j1b: " << cE2j1t << endl;
    cout << "Total number of events having exactly 2j2b: " << cE2j2t << endl;
    cout << "Total number of events having only one tight electron and only one tight muon: " << c1E1M << endl;
    cout << "Total number of events having only two tight electrons: " << c2E << endl;
    cout << "Total number of events having only two tight muons: " << c2M << endl;
    cout << "Total number of events having ooposite charge leptons 1e1mu: " << cOppChrg1e1m << endl;
    cout << "Total number of events having ooposite charge leptons 2e: " << cOppChrg2e << endl;
    cout << "Total number of events having ooposite charge leptons 2mu: " << cOppChrg2m << endl;
    cout << "total number of selected events due to extra loose cuts in 1e 1mu : " << cExLoose1e1m << endl;
    cout << "total number of selected events due to extra loose cuts in 2e : " << cExLoose2e << endl;
    cout << "total number of selected events due to extra loose cuts in 2mu : " << cExLoose2m << endl;
    cout << "total number of selected events after mll cuts in 1e 1mu : " << cmllc1e1m << endl;
    cout << "total number of selected events after mll cuts in 2e : " << cmllc2e << endl;
    cout << "total number of selected events after mll cuts in 2mu : " << cmllc2m << endl;
    cout << "total number of selected events due to met passing cut 1e1mu :" << cEMuMet1e1m << endl;
    cout << "total number of selected events due to met passing cut 2e :" << cEMuMet2e << endl;
    cout << "total number of selected events due to met passing cut 2mu :" << cEMuMet2m << endl;
    cout << "total number of selected events due to atleast 1 Loose jet for emu :" << cEwAtleast1jLemu << endl;
    cout << "total number of selected events due to atleast 1 Cleaned jet for emu :" << cEwAtleast1jClemu << endl;
    cout << "total number of selected events due to 1j jet for emu :" << cSelect1jem << endl;
    cout << "total number of selected events due to 1j1b jet for emu :" << cSelect1j1btem << endl;
    cout << "total number of selected events due to atleast 2 Loose jet for emu :" << cEwAtleast2jLemu << endl;
    cout << "total number of selected events due to atleast 2 Cleaned jet for emu :" << cEwAtleast2jClemu << endl;
    cout << "total number of selected events due to 2j jet for emu :" << cSelect2jem << endl;
    cout << "total number of selected events due to 2j1b jet for emu :" << cSelect2j1btem << endl;
    cout << "total number of selected events due to 2j2b jet for emu :" << cSelect2j2btem << endl;
    cout << "total number of selected events due to 1j jet for 2e :" << cSelect1j2e << endl;
    cout << "total number of selected events due to 1j1b jet for 2e :" << cSelect1j1bt2e << endl;
    cout << "total number of selected events due to 2j jet for 2e :" << cSelect2j2e << endl;
    cout << "total number of selected events due to 2j1b jet for 2e :" << cSelect2j1bt2e << endl;
    cout << "total number of selected events due to 2j2b jet for 2e :" << cSelect2j2bt2e << endl;
    cout << "total number of selected events due to 1j jet for 2mu :" << cSelect1j2m << endl;
    cout << "total number of selected events due to 1j1b jet for 2mu :" << cSelect1j1bt2m << endl;
    cout << "total number of selected events due to 2j jet for 2mu :" << cSelect2j2m << endl;
    cout << "total number of selected events due to 2j1b jet for 2mu :" << cSelect2j1bt2m << endl;
    cout << "total number of selected events due to 2j2b jet for 2mu :" << cSelect2j2bt2m << endl;
     */
    return 0;
}

                        


