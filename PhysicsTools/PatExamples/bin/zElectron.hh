#ifndef INC_13TEV_TW_ANALYSIS_ZELECTRON_HH
#define INC_13TEV_TW_ANALYSIS_ZELECTRON_HH

#include <ostream>
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

    void check_tight(bool debug = false)
    {
        if (debug)
            cout << "Debug is_tight" << endl;

        if (abs(this->eta_sc) <= 1.479)
            check_tight_eb(debug);
        else
            check_tight_ee(debug);

        if (debug)
            cout << "end" << endl;
    }

    void check_tight_eb(bool debug)
    {
        cout << "EB electron" << endl;
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

        int i = 0;

        for (auto it = conditions.begin(); it != conditions.end(); it++, i++)
            if (!(*it))
                cout << "Condition " << i << " failed" << endl;

        this->is_tight = std::all_of(conditions.cbegin(), conditions.cend(), [](bool b) { return b; });
    }

    void check_tight_ee(bool debug)
    {
        cout << "EE electron" << endl;
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

        int i = 0;

        for (auto it = conditions.begin(); it != conditions.end(); it++, i++)
            if (!(*it))
                cout << "Condition " << i << " failed" << endl;

        this->is_tight = std::all_of(conditions.cbegin(), conditions.cend(), [](bool b) { return b; });
    }


public:
    zElectron(const TLorentzVector &v, float charge, float eta_sc, float dxy, float d0, float dz,
              float full5x5_sigmaIetaIeta, float dEtaIn, float dPhiIn, float HoverE, float invE_invP,
              float rel_iso, int missing_hits, bool is_veto)
            : zParticle(v, charge),
              eta_sc(eta_sc),
              dxy(dxy), d0(d0),
              dz(dz), full5x5_sigmaIetaIeta(full5x5_sigmaIetaIeta),
              dEtaIn(dEtaIn), dPhiIn(dPhiIn), HoverE(HoverE), invE_invP(invE_invP), rel_iso(rel_iso),
              missing_hits(missing_hits), is_veto(is_veto)
    {
        this->check_tight();
    }

    zElectron() = default;

    bool get_istight() const
    {
        return this->is_tight;
    }

    float get_etaSC() const
    {
        return this->eta_sc;
    }

    bool in_gap() const
    {
        return fabs(this->eta_sc) > 1.4442 && fabs(this->eta_sc) < 1.566;
    }

    void debug_tight()
    {
        this->check_tight(true);
    }

    friend ostream &operator<<(ostream &os, const zElectron &electron)
    {
        os << " is_tight: " << electron.is_tight << " eta_sc: "
           << electron.eta_sc << " dxy: " << electron.dxy << " d0: " << electron.d0 << " dz: " << electron.dz
           << " full5x5_sigmaIetaIeta: " << electron.full5x5_sigmaIetaIeta << " dEtaIn: " << electron.dEtaIn
           << " dPhiIn: " << electron.dPhiIn << " HoverE: " << electron.HoverE << " invE_invP: " << electron.invE_invP
           << " rel_iso: " << electron.rel_iso << " missing_hits: " << electron.missing_hits << " is_veto: "
           << electron.is_veto << static_cast<const zParticle &>(electron);
        return os;
    }
};

#endif //INC_13TEV_TW_ANALYSIS_ZELECTRON_HH
