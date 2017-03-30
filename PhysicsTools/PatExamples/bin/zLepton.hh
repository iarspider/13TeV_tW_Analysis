//
// Created by razumov on 3/30/17.
//

#ifndef INC_13TEV_TW_ANALYSIS_ZLEPTON_HH
#define INC_13TEV_TW_ANALYSIS_ZLEPTON_HH

#include "zParticle.hh"


class zLepton : public zParticle
{
private:
    float eta_sc;
    float iso_;

private:
    bool is_tight;
    bool is_muon_;
public:
    bool is_muon() const
    {
        return is_muon_;
    }

    float is_iso() const
    {
        return (!this->is_muon()) || (iso_ < 0.15);
    }

public:
    zLepton(const TLorentzVector &v, float charge, float eta_sc, float iso, bool is_tight, bool is_muon)
            : zParticle(v, charge),
              eta_sc(eta_sc), iso_(iso),
              is_tight(is_tight), is_muon_(is_muon)
    {
    }

    zLepton() = default;

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
        return (!this->is_muon_) && (fabs(this->eta_sc) > 1.4442 && fabs(this->eta_sc) < 1.566);
    }

    friend ostream &operator<<(ostream &os, const zLepton &lepton)
    {
        if (!lepton.is_muon())
        {
            os << "Electron is_tight: " << lepton.is_tight << ", is_gap " << lepton.in_gap() << ", eta_sc: "
               << lepton.eta_sc << "; " << static_cast<const zParticle &>(lepton);
        }
        else
        {
            os << "Muon is_tight: " << lepton.is_tight << ", iso: " << lepton.iso_ << ", is_iso: " << lepton.is_iso() << "; "
               << static_cast<const zParticle &>(lepton);
        }
        return os;
    }
};

#endif //INC_13TEV_TW_ANALYSIS_ZLEPTON_HH
