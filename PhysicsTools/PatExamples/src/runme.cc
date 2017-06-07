#include "IIHEAnalysis.h"

void runme(TString listname)
{
    TChain *ch = new TChain("IIHEAnalysis");
    TFileCollection fc("dum","",listname);
    ch->AddFileInfoList(fc.GetList());
    TTree* tree = (TTree*) ch->Get("IIHEAnalysis");
    IIHEAnalysis code(tree);
    code.Loop();
}