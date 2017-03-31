#ifndef INC_13TEV_TW_ANALYSIS_ZELECTRON_HH
#define INC_13TEV_TW_ANALYSIS_ZELECTRON_HH

#include <ostream>
#include "zParticle.hh"

class zElectron : public zParticle
{
private:
    float eta_sc;
    bool is_tight;

public:
    zElectron(const TLorentzVector &v, float charge, float eta_sc, bool is_tight)
            : zParticle(v, charge),
              eta_sc(eta_sc),
              is_tight(is_tight)
    {
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

    friend ostream &operator<<(ostream &os, const zElectron &electron)
    {
        os << " is_tight_: " << electron.is_tight << ", is_gap " << electron.in_gap() << ", eta_sc: "
           << electron.eta_sc << "; " << static_cast<const zParticle &>(electron);
        return os;
    }
};

#endif //INC_13TEV_TW_ANALYSIS_ZELECTRON_HH
