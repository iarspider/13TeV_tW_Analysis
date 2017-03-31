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
    bool is_tight_;
    bool is_muon_;
    float pt_cut_, eta_cut_, pt_lead_cut_;

public:
    zLepton(const TLorentzVector &v, float charge, float eta_sc, float iso, bool is_tight, bool is_muon)
            : zParticle(v, charge),
              eta_sc(eta_sc), iso_(iso),
              is_tight_(is_tight), is_muon_(is_muon)
    {
        this->eta_cut_ = 2.4;
        this->pt_cut_ = 20;
        this->pt_lead_cut_ = 25;
    }


public:
    bool is_muon() const
    {
        return is_muon_;
    }

    bool is_tight() const
    {
        return this->is_tight_ && is_iso();
    }

    bool is_iso() const
    {
        return (!this->is_muon()) || (iso_ < 0.15);
    }

    bool in_gap() const
    {
        return (!this->is_muon_) && (fabs(this->eta_sc) > 1.4442 && fabs(this->eta_sc) < 1.566);
    }

    bool is_selected() const
    {
        return this->is_tight() && !this->in_gap() && this->Pt() > this->pt_cut_ && fabs(this->Eta()) < this->eta_cut_;
    }

    bool is_lead() const
    {
        return this->is_selected() && this->Pt() > this->pt_lead_cut_;
    }

    float get_etaSC() const
    {
        return this->eta_sc;
    }

    float getPt_cut_() const
    {
        return pt_cut_;
    }

    void setPt_cut_(float pt_cut_)
    {
        this->pt_cut_ = pt_cut_;
    }

    float getEta_cut_() const
    {
        return eta_cut_;
    }

    void setEta_cut_(float eta_cut_)
    {
        this->eta_cut_ = eta_cut_;
    }

    float getPt_lead_cut_() const
    {
        return pt_lead_cut_;
    }

    void setPt_lead_cut_(float pt_lead_cut_)
    {
        this->pt_lead_cut_ = pt_lead_cut_;
    }

    friend ostream &operator<<(ostream &os, const zLepton &lepton)
    {
        if (!lepton.is_muon())
        {
            os << "Electron is_tight_: " << lepton.is_tight_ << ", is_gap " << lepton.in_gap() << ", eta_sc: "
               << lepton.eta_sc << "; " << static_cast<const zParticle &>(lepton);
        }
        else
        {
            os << "Muon is_tight_: " << lepton.is_tight_ << ", iso: " << lepton.iso_ << ", is_iso: " << lepton.is_iso()
               << "; "
               << static_cast<const zParticle &>(lepton);
        }
        return os;
    }

};

#endif //INC_13TEV_TW_ANALYSIS_ZLEPTON_HH
