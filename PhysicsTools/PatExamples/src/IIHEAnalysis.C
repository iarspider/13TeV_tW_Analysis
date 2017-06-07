#define IIHEAnalysis_cxx

#include "IIHEAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void IIHEAnalysis::Loop() {
//   In a ROOT session, you can do:
//      root> .L IIHEAnalysis.C
//      root> IIHEAnalysis t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
    if (fChain == 0)
        return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0)
            break;

        b_mc_pdgId->GetEntry(jentry);
        b_mc_numberOfDaughters->GetEntry(jentry);
        b_mc_numberOfMothers->GetEntry(jentry);
        b_mc_mother_index->GetEntry(jentry);
        b_mc_mother_pdgId->GetEntry(jentry);
        b_mc_n->GetEntry(jentry);

        // t -> W+ -> q Qbar, tbar -> W- -> l(+-) nu (-+)

        std::vector <UInt_t> lep, nu, W, t, q;
        for(UInt_t i = 0; i < mc_n; i++)
        {
            UInt_t real_i = i;
            for (UInt_t j = 0; j < mc_n; j++)
            {
                if (((mc_numberOfMothers->at(j) == 1) && (mc_mother_index->at(j).at(0) == real_i)) && (mc_mother_pdgId->at(j).at(0) == mc_pdgId->at(real_i)))
                {
                    real_i = j;
                    j = 0;
                }
            }

            switch (abs(mc_pdgId->at(real_i)))
            {
                case 1:
                case 2:
                case 3:
                case 4:
                case 5:
                    q.push_back(real_i);
                    break;
                case 6:
                    t.push_back(real_i);
                    break;
                case 11:
                case 13:
                    lep.push_back(real_i);
                    break;
                case 12:
                case 14:
                    nu.push_back(real_i);
                    break;
                case 24:
                    W.push_back(real_i);
                    break;
                default:
                    break;
            }
        }
    }
}
