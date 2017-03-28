#ifndef INC_13TEV_TW_ANALYSIS_ZJET_HH
#define INC_13TEV_TW_ANALYSIS_ZJET_HH

#include "zParticle.hh"

class zJet : public zParticle
{
public:
    zJet(const TLorentzVector &v, float charge, float csv, float rapidity, float area, float nhef, float neef,
         float chef,
         float ceef, float n_const, float cmult, float nmult)
            : zParticle(v, charge), csv(csv),
              rapidity(rapidity), area(area),
              NHF(nhef), NEMF(neef), CHF(chef), CEMF(ceef), NumConst(n_const), cmult(cmult), NumNeutralParticle(nmult)
    {
        check_loose();
    }

    bool is_loose()
    {
        return this->loose_flag;
    }

private:
    float csv;
    float rapidity;
    float area;
    float NHF;
    float NEMF;
    float CHF;
    float CEMF;
    float NumConst;
    float cmult;
    float NumNeutralParticle;

    bool loose_flag;

    void check_loose_27()
    {
        std::vector<bool> conditions;
        conditions.push_back(this->NHF < 0.99);
        conditions.push_back(this->NEMF < 0.99);
        conditions.push_back(this->NumConst > 1);

        this->loose_flag = std::all_of(conditions.cbegin(), conditions.cend(), [](bool b) { return b; });
    }

    void check_loose_24()
    {
        std::vector<bool> conditions;
        conditions.push_back(this->CHF > 0);
        conditions.push_back(this->cmult > 0);
        conditions.push_back(this->CEMF < 0.99);

        this->loose_flag =
                this->loose_flag && std::all_of(conditions.cbegin(), conditions.cend(), [](bool b) { return b; });
    }

    void check_loose_30()
    {
        std::vector<bool> conditions;
        conditions.push_back(this->NHF < 0.98);
        conditions.push_back(this->NEMF > 0.01);
        conditions.push_back(this->NumNeutralParticle > 2);

        this->loose_flag = std::all_of(conditions.cbegin(), conditions.cend(), [](bool b) { return b; });
    }

    void check_loose_gt30()
    {
        std::vector<bool> conditions;
        conditions.push_back(this->NEMF < 0.90);
        conditions.push_back(this->NumNeutralParticle > 10);

        this->loose_flag = std::all_of(conditions.cbegin(), conditions.cend(), [](bool b) { return b; });
    }

    void check_loose()
    {
        if (abs(this->Eta()) <= 2.7)
        {
            check_loose_27();
            if (abs(this->Eta() <= 2.4))
            {
                check_loose_24();
            }
        }
        else if (abs(this->Eta()) > 2.7 && abs(this->Eta()) <= 3.0)
        {
            check_loose_30();
        }
        else
        {
            check_loose_gt30();
        }

    }


};

#endif //INC_13TEV_TW_ANALYSIS_ZJET_HH
