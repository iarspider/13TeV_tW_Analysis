#ifndef INC_13TEV_TW_ANALYSIS_ZPARTICLE_HH
#define INC_13TEV_TW_ANALYSIS_ZPARTICLE_HH

#include <TLorentzVector.h>

class zParticle : public TLorentzVector
{
private:
    float charge;
public:

    zParticle(TLorentzVector v, float charge) : TLorentzVector(v)
    {
        this->charge = charge;
    }

    zParticle() : TLorentzVector(), charge(0)
    {
    }

    float get_charge()
    {
        return this->charge;
    }

    Double_t deltaR(zParticle &other)
    {
        return this->DeltaR(TLorentzVector(other));
    }

    bool is_samesign(zParticle &other)
    {
        return (this->charge * other.get_charge() > 0);
    }

    bool operator<(const zParticle &rhs) const
    {
        return this->Pt() > rhs.Pt();
    }
};

#endif //INC_13TEV_TW_ANALYSIS_ZPARTICLE_HH
