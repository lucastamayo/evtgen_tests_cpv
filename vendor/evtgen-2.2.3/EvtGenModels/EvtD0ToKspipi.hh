#ifndef EVTD0TOKSPIPI_HH
#define EVTD0TOKSPIPI_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtDalitzPoint.hh"
#include "EvtGenBase/EvtDalitzReso.hh"
#include "EvtGenBase/EvtDecayAmp.hh"

#include <string>
#include <utility>
#include <vector>

class EvtParticle;

class EvtD0ToKspipi : public EvtDecayAmp {
  public:
    std::string getName() override;
    EvtDecayBase* clone() override;

    void init() override;
    void initProbMax() override;
    void decay( EvtParticle* parent ) override;

  private:
    // Calculate the total amplitude given the Dalitz plot point
    EvtComplex calcTotAmp( const EvtDalitzPoint& point ) const;

    // Set particle IDs and PDG masses
    void setPDGValues();

    // Setup the Dalitz plot resonances and their amplitude coefficients
    void initResonances();

    // Daughter IDs (updated according to decay file ordering)
    int m_d0 = 0;
    int m_d1 = 1;
    int m_d2 = 2;

    // Resonance lineshape and complex amplitude coefficient pair
    typedef std::pair<EvtDalitzReso, EvtComplex> ResAmpPair;

    // Vector of (resonance, coeff) pairs
    std::vector<ResAmpPair> m_resonances;

    // IDs of the relevant particles
    EvtId m_BP;
    EvtId m_BM;
    EvtId m_B0;
    EvtId m_B0B;
    EvtId m_D0;
    EvtId m_D0B;
    EvtId m_KM;
    EvtId m_KP;
    EvtId m_K0;
    EvtId m_K0B;
    EvtId m_KL;
    EvtId m_KS;
    EvtId m_PIM;
    EvtId m_PIP;

    // Flavor of the B parent and D
    EvtId m_bFlavor;
    EvtId m_dFlavor;

    // Masses of the relevant particles
    double m_mD0;
    double m_mKs;
    double m_mPi;
    double m_mK;
};
#endif
