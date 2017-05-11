#ifndef INC_13TEV_TW_ANALYSIS_ZJET_HH
#define INC_13TEV_TW_ANALYSIS_ZJET_HH

#include <ostream>
#include "zParticle.hh"

class zJet : public zParticle
{
public:
    zJet(const TLorentzVector &v, float charge, float btag, bool loose_flag)
            : zParticle(v, charge), btag(btag),
              /*rapidity(rapidity), area(area),*/
              loose_flag(loose_flag)
    {
        this->setPt_cut_(30);
        this->setEta_cut_(2.4);
        this->setBtag_cut_(0.8484);
    }

    bool is_loose() const
    {
        return this->loose_flag;
    }

    bool is_bjet() const
    {
        return this->btag > this->btag_cut_;
    }

    bool is_selected() const
    {
        return this->is_loose() && this->is_clean() && (this->Pt() > this->pt_cut_) && (fabs(this->Eta()) < this->eta_cut_);
    }

    float get_csv() const
    {
        return this->btag;
    }

private:
    //friend std::ostream &operator<<(std::ostream &, const zJet &);

    float btag;
//    float rapidity;
//    float area;

    float pt_cut_, eta_cut_, btag_cut_;

    bool loose_flag;
    bool clean_flag;

public:
    bool is_clean() const
    {
        return clean_flag;
    }

    void set_clean_flag(bool clean_flag)
    {
        this->clean_flag = clean_flag;
    }

    friend ostream &operator<<(ostream &os, const zJet &jet)
    {
        return os << "Jet is_clean = "
                  << jet.is_clean() << ", is_loose = " << jet.is_loose() << ", tag =" << jet.is_bjet() << "; "
                  << static_cast<const zParticle & >(jet);
    }

    void setPt_cut_(float pt_cut_)
    {
        this->pt_cut_ = pt_cut_;
    }

    void setEta_cut_(float eta_cut_)
    {
        this->eta_cut_ = eta_cut_;
    }

    void setBtag_cut_(float btag_cut_)
    {
        this->btag_cut_ = btag_cut_;
    }

};
/*
std::ostream &operator<<(std::ostream &strm, const zJet &j)
{
    // return strm << "A(" << a.j << ")";
    return strm << "Jet(" << j.Pt() << ", " << j.Eta() << ", " << j.Phi() << ", " << j.Energy() << ") clean= "
                << j.is_clean() << ", loose = " << j.is_loose() << ", tag =" << j.is_bjet();
}
*/
#endif //INC_13TEV_TW_ANALYSIS_ZJET_HH
