#ifndef INC_13TEV_TW_ANALYSIS_ZELECTRON_HH
#define INC_13TEV_TW_ANALYSIS_ZELECTRON_HH

#include "zParticle.hh"

class zElectron : public zParticle
{
private:
    bool is_tight;
    float eta_sc;
    float dxy, d0, dz;
    float full5x5_sigmaIetaIeta, dEtaIn, dPhiIn, HoverE, invE_invP, rel_iso;
    int missing_hits;
    bool is_veto;

    void check_tight()
    {
        if (abs(this->eta_sc) <= 1.479)
        {
            check_tight_eb();
        }
        else
        {
            check_tight_ee();
        }
    }

    void check_tight_eb()
    {
        vector<bool> conditions;
        conditions.push_back(this->full5x5_sigmaIetaIeta < 0.00998);
        conditions.push_back(fabs(this->dEtaIn) < 0.00308);
        conditions.push_back(fabs(this->dPhiIn) < 0.0816);
        conditions.push_back(this->HoverE < 0.0414);
        conditions.push_back(this->rel_iso < 0.0588);
        conditions.push_back(fabs(this->invE_invP) < 0.0129);
        conditions.push_back(this->missing_hits <= 1);
        conditions.push_back(!this->is_veto);
        conditions.push_back(fabs(this->d0) < 0.05);
        conditions.push_back(fabs(this->dz) < 0.10);

        this->is_tight = std::all_of(conditions.cbegin(), conditions.cend(), [](bool b) { return b; });
    }

    void check_tight_ee()
    {
        vector<bool> conditions;
        conditions.push_back(this->full5x5_sigmaIetaIeta < 0.0292);
        conditions.push_back(fabs(this->dEtaIn) < 0.00605);
        conditions.push_back(fabs(this->dPhiIn) < 0.0394);
        conditions.push_back(this->HoverE < 0.0641);
        conditions.push_back(this->rel_iso < 0.0571);
        conditions.push_back(fabs(this->invE_invP) < 0.0129);
        conditions.push_back(this->missing_hits <= 1);
        conditions.push_back(!this->is_veto);
        conditions.push_back(fabs(this->d0) < 0.10);
        conditions.push_back(fabs(this->dz) < 0.20);

        this->is_tight = std::all_of(conditions.cbegin(), conditions.cend(), [](bool b) { return b; });
    }


public:
    zElectron(const TLorentzVector &v, float charge, float eta_sc, float dxy, bool is_veto, float invE_invP,
              float full5x5_sigmaIetaIeta, float dEtaIn, float dPhiIn, float HoverE, int missing_hits,
              float rel_iso, float d0, float dz)
            : zParticle(v, charge),
              eta_sc(eta_sc),
              dxy(dxy), is_veto(is_veto), invE_invP(invE_invP), full5x5_sigmaIetaIeta(full5x5_sigmaIetaIeta),
              dEtaIn(dEtaIn), dPhiIn(dPhiIn), HoverE(HoverE), missing_hits(missing_hits), rel_iso(rel_iso), d0(d0),
              dz(dz)
    {
        this->check_tight();
    }

    zElectron() : zParticle(), eta_sc(0), dxy(0), is_veto(false), invE_invP(0), full5x5_sigmaIetaIeta(0),
                  dEtaIn(0), dPhiIn(0), HoverE(0), missing_hits(0), rel_iso(0), d0(0),
                  dz(0)
    {
    }

    bool get_is_tight()
    {
        return this->is_tight;
    }

    float get_etaSC()
    {
        return this->eta_sc;
    }

};
#endif //INC_13TEV_TW_ANALYSIS_ZELECTRON_HH
