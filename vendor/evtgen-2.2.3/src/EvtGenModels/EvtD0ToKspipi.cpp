#include "EvtGenModels/EvtD0ToKspipi.hh"

#include "EvtGenBase/EvtDecayTable.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <iostream>

std::string EvtD0ToKspipi::getName()
{
    return "D0TOKSPIPI";
}

EvtDecayBase* EvtD0ToKspipi::clone()
{
    return new EvtD0ToKspipi;
}

void EvtD0ToKspipi::init()
{
    // Check that there are 0 arguments
    checkNArg( 0 );

    // Check parent and daughter types
    checkNDaug( 3 );
    checkSpinDaughter( 0, EvtSpinType::SCALAR );
    checkSpinDaughter( 1, EvtSpinType::SCALAR );
    checkSpinDaughter( 2, EvtSpinType::SCALAR );
    checkSpinParent( EvtSpinType::SCALAR );

    // Set the particle IDs and PDG masses
    setPDGValues();

    // Set the EvtId of the three D0 daughters
    const int nDaug = 3;
    std::vector<EvtId> dau;
    for ( int index = 0; index < nDaug; index++ ) {
        dau.push_back( getDaug( index ) );
    }

    // Look for K0bar h+ h-. The order will be K[0SL] h+ h-
    for ( int index = 0; index < nDaug; index++ ) {
        if ( ( dau[index] == m_K0B ) || ( dau[index] == m_KS ) ||
             ( dau[index] == m_KL ) ) {
            m_d0 = index;
        } else if ( dau[index] == m_PIP ) {
            m_d1 = index;
        } else if ( dau[index] == m_PIM ) {
            m_d2 = index;
        } else {
            EvtGenReport( EVTGEN_ERROR, "EvtD0ToKspipi" )
                << "Daughter " << index << " has wrong ID" << std::endl;
        }
    }

    // Setup the Dalitz plot resonances and their amplitude coefficients
    initResonances();
}

void EvtD0ToKspipi::initProbMax()
{
    setProbMax( 6000.0 );
}

void EvtD0ToKspipi::decay( EvtParticle* p )
{
    // Phase space generation and 4-momenta
    p->initializePhaseSpace( getNDaug(), getDaugs() );
    const EvtVector4R p0 = p->getDaug( m_d0 )->getP4();
    const EvtVector4R p1 = p->getDaug( m_d1 )->getP4();
    const EvtVector4R p2 = p->getDaug( m_d2 )->getP4();

    // Squared invariant masses
    const double mSq01 = ( p0 + p1 ).mass2();
    const double mSq02 = ( p0 + p2 ).mass2();
    const double mSq12 = ( p1 + p2 ).mass2();

    // For the decay amplitude
    EvtComplex amp( 0.0, 0.0 );

    // Flavor of the D meson
    m_dFlavor = p->getId();

    // Direct and conjugated Dalitz points
    const EvtDalitzPoint pointD0( m_mKs, m_mPi, m_mPi, mSq02, mSq12, mSq01 );
    const EvtDalitzPoint pointD0b( m_mKs, m_mPi, m_mPi, mSq01, mSq12, mSq02 );

    // Check if the D is from a B+- -> D0 K+- decay with the appropriate model
    EvtParticle* parent = p->getParent();
    if ( parent != nullptr &&
         EvtDecayTable::getInstance()->getDecayFunc( parent )->getName() ==
             "BTODDALITZCPK" ) {
        const EvtId parId = parent->getId();
        if ( ( parId == m_BP ) || ( parId == m_BM ) || ( parId == m_B0 ) ||
             ( parId == m_B0B ) ) {
            m_bFlavor = parId;

            // D0 parent particle is a B meson from the BTODDALITZCPK decay model.
            // D0 decay amplitude combines the interference of D0 and D0bar.
            // Read the D decay parameters from the B decay model.
            // Gamma angle in radians
            const double gamma =
                EvtDecayTable::getInstance()->getDecayFunc( parent )->getArg( 0 );
            // Strong phase in radians
            const double delta =
                EvtDecayTable::getInstance()->getDecayFunc( parent )->getArg( 1 );
            // Ratio between B -> D0 K and B -> D0bar K
            const double rB =
                EvtDecayTable::getInstance()->getDecayFunc( parent )->getArg( 2 );

            // Direct and conjugated amplitudes
            const EvtComplex ampD0 = calcTotAmp( pointD0 );
            const EvtComplex ampD0b = calcTotAmp( pointD0b );

            if ( m_bFlavor == m_BP || m_bFlavor == m_B0 ) {
                // B+ or B0
                const EvtComplex iPhase( 0.0, delta + gamma );
                const EvtComplex expo( exp( iPhase ) );
                amp = ampD0b + rB * expo * ampD0;
            } else {
                // B- or B0bar
                const EvtComplex iPhase( 0.0, delta - gamma );
                const EvtComplex expo( exp( iPhase ) );
                amp = ampD0 + rB * expo * ampD0b;
            }
        }
    } else if ( !parent ) {
        // D0 has no parent particle. Use direct or conjugated amplitude
        if ( m_dFlavor == m_D0 ) {
            amp = calcTotAmp( pointD0 );
        } else {
            amp = calcTotAmp( pointD0b );
        }
    }

    // Set the decay vertex amplitude
    vertex( amp );
}

EvtComplex EvtD0ToKspipi::calcTotAmp( const EvtDalitzPoint& point ) const
{
    // Initialise the total amplitude
    EvtComplex totAmp( 0.0, 0.0 );
    // Check that the Dalitz plot point is OK
    if ( point.isValid() == false ) {
        return totAmp;
    }

    // Add the resonance amplitudes by iterating over the (resonance, coeff) pairs.
    // This includes the BW and LASS lineshapes, as well as the K-matrix contributions
    for ( auto& [res, amp] : m_resonances ) {
        // Evaluate the resonance amplitude and multiply by the coeff
        totAmp += res.evaluate( point ) * amp;
    }
    // Return the total amplitude
    return totAmp;
}

void EvtD0ToKspipi::initResonances()
{
    // Dalitz plot model from combined BaBar and BELLE paper hep-ex/1804.06153

    // Define the Dalitz plot axes
    const EvtCyclic3::Pair AB = EvtCyclic3::AB;
    const EvtCyclic3::Pair AC = EvtCyclic3::AC;
    const EvtCyclic3::Pair BC = EvtCyclic3::BC;

    // Define the particle spin and lineshape types
    const EvtSpinType::spintype vector = EvtSpinType::VECTOR;
    const EvtSpinType::spintype tensor = EvtSpinType::TENSOR;
    const EvtDalitzReso::NumType RBW = EvtDalitzReso::RBW_CLEO_ZEMACH;
    const EvtDalitzReso::NumType KMAT = EvtDalitzReso::K_MATRIX;

    // Define the Dalitz plot
    const EvtDalitzPlot DP( m_mKs, m_mPi, m_mPi, m_mD0, 0, 0 );

    // Clear the internal vector of (resonance, coeff) pairs
    m_resonances.clear();

    // rho BW
    const EvtDalitzReso rhoBW( DP, AB, BC, vector, 0.77155, 0.13469, RBW, 5.0,
                               1.5 );
    const EvtComplex rhoCoeff( 1.0, 0.0 );
    m_resonances.push_back( std::make_pair( rhoBW, rhoCoeff ) );

    // Omega BW
    const EvtDalitzReso omegaBW( DP, AB, BC, vector, 0.78265, 0.00849, RBW, 5.0,
                                 1.5 );
    const EvtComplex omegaCoeff( -0.019829903319132, 0.033339785741436 );
    m_resonances.push_back( std::make_pair( omegaBW, omegaCoeff ) );

    // K*(892)- BW
    const EvtDalitzReso KstarBW( DP, BC, AB, vector, 0.893709298220334,
                                 0.047193287094108, RBW, 5.0, 1.5 );
    const EvtComplex KstarCoeff( -1.255025021860793, 1.176780750003210 );
    m_resonances.push_back( std::make_pair( KstarBW, KstarCoeff ) );

    // K*0(1430)- LASS
    const double LASS_F = 0.955319683174069;
    const double LASS_phi_F = 0.001737032480754;
    const double LASS_R = 1.0;
    const double LASS_phi_R = -1.914503836666840;
    const double LASS_a = 0.112673863011817;
    const double LASS_r = -33.799002116066454;
    const EvtDalitzReso Kstar0_1430LASS = EvtDalitzReso(
        DP, AB, 1.440549945739415, 0.192611512914605, LASS_a, LASS_r, LASS_F,
        LASS_phi_F, LASS_R, LASS_phi_R, -1.0, false );
    const EvtComplex Kstar0_1430Coeff( -0.386469884688245, 2.330315087713914 );
    m_resonances.push_back( std::make_pair( Kstar0_1430LASS, Kstar0_1430Coeff ) );

    // K*2(1430)- BW
    const EvtDalitzReso Kstar2_1430BW( DP, BC, AB, tensor, 1.4256, 0.0985, RBW,
                                       5.0, 1.5 );
    const EvtComplex Kstar2_1430Coeff( 0.914470111251261, -0.885129049790117 );
    m_resonances.push_back( std::make_pair( Kstar2_1430BW, Kstar2_1430Coeff ) );

    // K*(1680)- BW
    const EvtDalitzReso Kstar_1680BW( DP, BC, AB, vector, 1.717, 0.322, RBW,
                                      5.0, 1.5 );
    const EvtComplex Kstar_1680Coeff( -1.560837188791231, -2.916210561577914 );
    m_resonances.push_back( std::make_pair( Kstar_1680BW, Kstar_1680Coeff ) );

    // K*(1410)- BW
    const EvtDalitzReso Kstar_1410BW( DP, BC, AB, vector, 1.414, 0.232, RBW,
                                      5.0, 1.5 );
    const EvtComplex Kstar_1410Coeff( -0.046795079734847, 0.283085379985959 );
    m_resonances.push_back( std::make_pair( Kstar_1410BW, Kstar_1410Coeff ) );

    // K*(892)+ DCS BW
    const EvtDalitzReso Kstar_DCSBW( DP, BC, AC, vector, 0.893709298220334,
                                     0.047193287094108, RBW, 5.0, 1.5 );
    const EvtComplex Kstar_DCSCoeff( 0.121693743404499, -0.110206354657867 );
    m_resonances.push_back( std::make_pair( Kstar_DCSBW, Kstar_DCSCoeff ) );

    // K*0(1430)+ DCS LASS
    const EvtDalitzReso Kstar0_1430_DCSLASS = EvtDalitzReso(
        DP, AC, 1.440549945739415, 0.192611512914605, LASS_a, LASS_r, LASS_F,
        LASS_phi_F, LASS_R, LASS_phi_R, -1.0, false );
    const EvtComplex Kstar0_1430_DCSCoeff( -0.101484805664368, 0.032368302993344 );
    m_resonances.push_back(
        std::make_pair( Kstar0_1430_DCSLASS, Kstar0_1430_DCSCoeff ) );

    // K*2(1430)+ DCS BW
    const EvtDalitzReso Kstar2_1430_DCSBW( DP, AB, AC, tensor, 1.4256, 0.0985,
                                           RBW, 5.0, 1.5 );
    const EvtComplex Kstar2_1430_DCSCoeff( 0.000699701539252, -0.102571188336701 );
    m_resonances.push_back(
        std::make_pair( Kstar2_1430_DCSBW, Kstar2_1430_DCSCoeff ) );

    // K*(1410)+ DCS BW
    const EvtDalitzReso Kstar_1410_DCSBW( DP, BC, AC, vector, 1.414, 0.232, RBW,
                                          5.0, 1.5 );
    const EvtComplex Kstar_1410_DCSCoeff( -0.181330401419455, 0.103990039950039 );
    m_resonances.push_back(
        std::make_pair( Kstar_1410_DCSBW, Kstar_1410_DCSCoeff ) );

    // f2(1270) BW
    const EvtDalitzReso f2_1270BW( DP, AB, BC, tensor, 1.2751, 0.1842, RBW, 5.0,
                                   1.5 );
    const EvtComplex f2_1270Coeff( 1.151785277682948, -0.845612891825272 );
    m_resonances.push_back( std::make_pair( f2_1270BW, f2_1270Coeff ) );

    // rho(1450) BW
    const EvtDalitzReso rho_1450BW( DP, AB, BC, vector, 1.465, 0.400, RBW, 5.0,
                                    1.5 );
    const EvtComplex rho_1450Coeff( -0.597963342540235, 2.787903868470057 );
    m_resonances.push_back( std::make_pair( rho_1450BW, rho_1450Coeff ) );

    // K-matrix pole 1
    const double sProd0 = -0.07;
    const EvtDalitzReso pole1( DP, BC, "Pole1", KMAT, 0, 0, 0, 0, sProd0 );
    const EvtComplex p1Coeff( 3.122415682166643, 7.928823290976309 );
    m_resonances.push_back( std::make_pair( pole1, p1Coeff ) );

    // K-matrix pole 2
    const EvtDalitzReso pole2( DP, BC, "Pole2", KMAT, 0, 0, 0, 0, sProd0 );
    const EvtComplex p2Coeff( 11.139907856904129, 4.948420661321371 );
    m_resonances.push_back( std::make_pair( pole2, p2Coeff ) );

    // K-matrix pole 3
    const EvtDalitzReso pole3( DP, BC, "Pole3", KMAT, 0, 0, 0, 0, sProd0 );
    const EvtComplex p3Coeff( 29.146102368470210, -0.053588781806890 );
    m_resonances.push_back( std::make_pair( pole3, p3Coeff ) );

    // K-matrix pole 4
    const EvtDalitzReso pole4( DP, BC, "Pole4", KMAT, 0, 0, 0, 0, sProd0 );
    const EvtComplex p4Coeff( 6.631556203215280, -8.455370251307063 );
    m_resonances.push_back( std::make_pair( pole4, p4Coeff ) );

    // K-matrix pole 5 is not included since its amplitude coefficient is zero

    // K-matrix slowly varying part
    const EvtComplex fProd11( -4.724094278696236, -6.511009103363590 );
    const EvtComplex fProd12( -23.289333360304212, -12.215597571354197 );
    const EvtComplex fProd13( -1.860311896516422, -32.982507366353126 );
    const EvtComplex fProd14( -13.638752211193912, -22.339804683783186 );
    const EvtComplex fProd15( 0.0, 0.0 );

    const EvtDalitzReso KMSVP( DP, BC, "f11prod", KMAT, fProd12 / fProd11,
                               fProd13 / fProd11, fProd14 / fProd11,
                               fProd15 / fProd11, sProd0 );
    m_resonances.push_back( std::make_pair( KMSVP, fProd11 ) );
}

void EvtD0ToKspipi::setPDGValues()
{
    // Set the EvtIds
    m_BP = EvtPDL::getId( "B+" );
    m_BM = EvtPDL::getId( "B-" );
    m_B0 = EvtPDL::getId( "B0" );
    m_B0B = EvtPDL::getId( "anti-B0" );
    m_D0 = EvtPDL::getId( "D0" );
    m_D0B = EvtPDL::getId( "anti-D0" );
    m_KM = EvtPDL::getId( "K-" );
    m_KP = EvtPDL::getId( "K+" );
    m_K0 = EvtPDL::getId( "K0" );
    m_K0B = EvtPDL::getId( "anti-K0" );
    m_KL = EvtPDL::getId( "K_L0" );
    m_KS = EvtPDL::getId( "K_S0" );
    m_PIM = EvtPDL::getId( "pi-" );
    m_PIP = EvtPDL::getId( "pi+" );

    // Set the particle masses
    m_mD0 = EvtPDL::getMeanMass( m_D0 );
    m_mKs = EvtPDL::getMeanMass( m_KS );
    m_mPi = EvtPDL::getMeanMass( m_PIP );
    m_mK = EvtPDL::getMeanMass( m_KP );
}
