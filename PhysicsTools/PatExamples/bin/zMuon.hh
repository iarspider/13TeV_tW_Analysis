#ifndef INC_13TEV_TW_ANALYSIS_ZMUON_HH
#define INC_13TEV_TW_ANALYSIS_ZMUON_HH

#include "zParticle.hh"

class zMuon: public zParticle {
private:
    float iso;
    bool isTight;

public:
    bool get_istight() const
    {
        return this->isTight;
    }

public:
    zMuon(): zParticle(), iso(0), isTight(false) {}

    zMuon(TLorentzVector v, float charge, float iso, bool isTight): zParticle(v, charge), iso(iso), isTight(isTight) {}

    friend ostream &operator<<(ostream &os, const zMuon &muon)
    {
        os << "Muon tight " << muon.get_istight() << "; " << static_cast<const zParticle &>(muon);
        return os;
    }

};

#endif //INC_13TEV_TW_ANALYSIS_ZMUON_HH
