#ifndef INC_13TEV_TW_ANALYSIS_ZPARTICLE_HH
#define INC_13TEV_TW_ANALYSIS_ZPARTICLE_HH

#include <TLorentzVector.h>

class zParticle : public TLorentzVector
{
private:
    float charge;
public:
    friend std::ostream &operator<<(std::ostream & o, const zParticle & z) {
        o << "Charge" << z.get_charge() << ", (Pt, Eta, Phi, E) = (" << z.Pt() << ", " << z.Eta() << ", " << z.Phi() << ", " << z.E() << ")";
    }

    zParticle(TLorentzVector v, float charge) : TLorentzVector(v)
    {
        this->charge = charge;
    }

    zParticle() : TLorentzVector(), charge(0)
    {
    }

    float get_charge() const
    {
        return this->charge;
    }

    Double_t deltaR(zParticle &other) const
    {
        return this->DeltaR(TLorentzVector(other));
    }

    bool is_samesign(zParticle &other) const
    {
        return (this->charge * other.get_charge() > 0);
    }

    bool operator<(const zParticle &rhs) const
    {
        return this->Pt() > rhs.Pt();
    }
};

#endif //INC_13TEV_TW_ANALYSIS_ZPARTICLE_HH
