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

using namespace std;

#include "zEvent.hh"

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
    FWLiteEnabler::enable();

    int n = 2;
    //int s=atoi(argv[2]);
    //int n=atoi(argv[3]);
    //int jobid=atoi(argv[1]);

    char fname[160];

    std::vector<zEvent> events;
    double cNetEvWt = 0;
    int cTHltEv = 0, cVertexEv = 0;

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
            cout << "File Name:" << fname << endl;

            TFile *inFile = TFile::Open(fname);

            if (!inFile)
            {
                cout << "Input Root File not opened for processing " << endl;
                return 1;
            }

            fwlite::Event ev(inFile);
            int evtID = 0;

            for (ev.toBegin(); !ev.atEnd(); ++ev, evtID++)
            {
                edm::EventBase const &event = ev;
                events.push_back(zEvent(event));
            }
        }
        catch (exception &e)
        {
            cout << "error:" << e.what() << endl;
            cout << "File Name:" << fname << "is not valid!" << endl << endl << endl;
        }

        for (auto event = events.begin(); event != events.end(); event++)
        {
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

        }
    }

    cout << "Total number of events processed: " << events.size() << endl;
    cout << "Sum of events Weight : " << cNetEvWt << endl;
    cout << "Total number of events passed HLT cuts: " << cTHltEv << endl;
    cout << "Total number of events passed vertex cuts: " << cVertexEv << endl;
    /*
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

                        


