#ifndef INC_13TEV_TW_ANALYSIS_ZHLT_HH
#define INC_13TEV_TW_ANALYSIS_ZHLT_HH

#include <string>

class zHLT
{
public:
    std::string HLTName;
    int HTLPrescale;
    float HTLDecision;

    zHLT(std::string HLTName, int HTLPrescale, float HTLDecision)
            : HLTName(HLTName), HTLPrescale(HTLPrescale), HTLDecision(HTLDecision)
    {
    }
};
#endif //INC_13TEV_TW_ANALYSIS_ZHLT_HH
