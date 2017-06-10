//
// Created by razumov on 6/10/17.
//
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TFileCollection.h>
#include <TH1D.h>

int main(int, char **) {
    TFile f("MyMCTruePileupHistogram", "RECREATE");
    f.cd();
    TChain ch("IIHEAnalysis");
    TFileCollection fc("dum", "", "xrd.list");
    ch.AddFileInfoList((TCollection*)(fc.GetList()));
    Int_t mc_trueNumInteractions;
    ch.SetBranchStatus("*", 0);
    ch.SetBranchStatus("mc_trueNumInteractions", 1);
    ch.SetBranchAddress("mc_trueNumInteractions", &mc_trueNumInteractions);
    TH1D hPileupMC("h", "h", 750, 0, 75);

    for (UInt_t i = 0; i < ch.GetEntries(); i++)
    {
        ch.GetEntry(i);
        hPileupMC.Fill(mc_trueNumInteractions);
    }

    hPileupMC.Write();
    f.Close();
}