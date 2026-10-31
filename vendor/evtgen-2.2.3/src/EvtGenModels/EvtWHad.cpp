
/***********************************************************************
* Copyright 1998-2020 CERN for the benefit of the EvtGen authors       *
*                                                                      *
* This file is part of EvtGen.                                         *
*                                                                      *
* EvtGen is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation, either version 3 of the License, or    *
* (at your option) any later version.                                  *
*                                                                      *
* EvtGen is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of       *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
* GNU General Public License for more details.                         *
*                                                                      *
* You should have received a copy of the GNU General Public License    *
* along with EvtGen.  If not, see <https://www.gnu.org/licenses/>.     *
***********************************************************************/

#include "EvtGenModels/EvtWHad.hh"

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtTensor4C.hh"

EvtWHad::EvtWHad() :
    mRho_(),
    gamma0_(),
    cK_( 0 ),
    mK_(),
    gammaK_(),
    gKRho_(),
    gKPi_(),
    mPi_( EvtPDL::getMeanMass( EvtPDL::getId( "pi+" ) ) ),
    mPiSq_( mPi_ * mPi_ )
{
    // cK coefficients from Eur. Phys. J. C39, 41 (2005), arXiv:hep-ph/0409080 [hep-ph]

    // rho(770)
    mRho_.push_back( EvtPDL::getMeanMass( EvtPDL::getId( "rho0" ) ) );
    gamma0_.push_back( EvtPDL::getWidth( EvtPDL::getId( "rho0" ) ) );
    cK_.push_back( 1.195 );

    // rho(1450)
    mRho_.push_back( EvtPDL::getMeanMass( EvtPDL::getId( "rho(2S)0" ) ) );
    gamma0_.push_back( EvtPDL::getWidth( EvtPDL::getId( "rho(2S)0" ) ) );
    cK_.push_back( -0.112 );

    // rho(1700)
    mRho_.push_back( EvtPDL::getMeanMass( EvtPDL::getId( "rho(3S)0" ) ) );
    gamma0_.push_back( EvtPDL::getWidth( EvtPDL::getId( "rho(3S)0" ) ) );
    cK_.push_back( -0.083 );

    // rho(2150), PRD 76 092005
    mRho_.push_back( 2.150 );
    gamma0_.push_back( 0.310 );
    cK_.push_back( 0.0 );

    // Storing K resonance information

    // K(892)

    mK_.push_back( EvtPDL::getMeanMass( EvtPDL::getId( "K*0" ) ) );
    gammaK_.push_back( EvtPDL::getWidth( EvtPDL::getId( "K*0" ) ) );
    gKRho_.push_back( 0.0 );
    gKPi_.push_back( 3.26 );

    // K1(1270)

    mK_.push_back( EvtPDL::getMeanMass( EvtPDL::getId( "K_10" ) ) );
    gammaK_.push_back( EvtPDL::getWidth( EvtPDL::getId( "K_10" ) ) );
    gKRho_.push_back( 2.71 );
    gKPi_.push_back( 0.792 );

    // K1(1400)
    mK_.push_back( EvtPDL::getMeanMass( EvtPDL::getId( "K'_10" ) ) );
    gammaK_.push_back( EvtPDL::getWidth( EvtPDL::getId( "K'_10" ) ) );
    gKRho_.push_back( 0.254 );
    gKPi_.push_back( 2.509 );
}

EvtComplex EvtWHad::BWKK( double s, int i ) const
{
    const double m2 = mRho_[i] * mRho_[i];
    const EvtComplex qs = pcm( s );
    const EvtComplex qm = pcm( m2 );
    if ( abs( qm ) < 1e-10 ) {
        return 0;
    }

    const EvtComplex rat = qs / qm;
    const EvtComplex rat3 = rat * rat * rat;
    if ( abs( s ) < 1e-10 ) {
        return 0;
    }

    const EvtComplex gamma = m2 * rat3 * gamma0_[i] / s;
    const EvtComplex I( 0.0, 1.0 );

    const EvtComplex denBW = m2 - s - I * sqrt( s ) * gamma;
    if ( abs( denBW ) < 1e-10 ) {
        return 0;
    }

    return cK_[i] * m2 / denBW;
}

EvtVector4C EvtWHad::WCurrent_KSK( const EvtVector4R& pKS,
                                   const EvtVector4R& pKplus ) const
{
    const double s = ( pKS + pKplus ).mass2();
    const EvtComplex f = BWKK( s, 0 ) + BWKK( s, 1 ) + BWKK( s, 2 );

    return f * ( pKS - pKplus );
}

EvtComplex EvtWHad::pcm( double s ) const
{
    const double mpi2 = 0.0196;    // 0.140*0.140
    if ( abs( s ) < 1e-10 )
        return 0;

    const double pcm2 = 1.0 - 4.0 * mpi2 / s;
    EvtComplex result;

    if ( pcm2 >= 0.0 ) {
        result = EvtComplex( sqrt( pcm2 ), 0.0 );
    } else {
        result = EvtComplex( 0.0, sqrt( -pcm2 ) );
    }

    return result;
}

// =================== W+ -> pi_ current ========================================

EvtVector4C EvtWHad::WCurrent( const EvtVector4R& q1 ) const
{
    return q1;
}

//====================== W+ -> pi+ pi0 current =========================================

EvtVector4C EvtWHad::WCurrent( const EvtVector4R& q1, const EvtVector4R& q2 ) const
{
    return BWr( q1 + q2 ) * ( q1 - q2 );
}

//========================= W+ -> pi+ pi+ pi- current ==============================================

EvtVector4C EvtWHad::WCurrent( const EvtVector4R& q1, const EvtVector4R& q2,
                               const EvtVector4R& q3 ) const
{
    const EvtVector4R Q = q1 + q2 + q3;
    const EvtVector4R q13 = q1 - q3, q23 = q2 - q3;
    const double Q2 = Q.mass2();
    return BWa( Q ) * ( q13 - ( Q * ( Q * q13 ) / Q2 ) * BWr( q2 + q3 ) + q23 -
                        ( Q * ( Q * q23 ) / Q2 ) * BWr( q1 + q3 ) );
}

// ================= W+ -> pi+ pi+ pi- pi- pi+ current with symmetrization ================================

EvtVector4C EvtWHad::WCurrent( const EvtVector4R& q1, const EvtVector4R& q2,
                               const EvtVector4R& q3, const EvtVector4R& q4,
                               const EvtVector4R& q5 ) const
{
    const EvtVector4C term1 = JB( q1, q2, q3, q4, q5 );
    const EvtVector4C term2 = JB( q5, q2, q3, q4, q1 );
    const EvtVector4C term3 = JB( q1, q5, q3, q4, q2 );
    const EvtVector4C term4 = JB( q1, q2, q4, q3, q5 );
    const EvtVector4C term5 = JB( q5, q2, q4, q3, q1 );
    const EvtVector4C term6 = JB( q1, q5, q4, q3, q2 );

    const EvtVector4C V = term1 + term2 + term3 + term4 + term5 + term6;
    return V;
}

// W+ -> pi+ pi+ pi+ pi- pi-
EvtVector4C EvtWHad::WCurrent_5pi( const EvtVector4R& q1, const EvtVector4R& q2,
                                   const EvtVector4R& q3, const EvtVector4R& q4,
                                   const EvtVector4R& q5 ) const
{
    return EvtWHad::WCurrent( q1, q2, q4, q5, q3 );    // WCurrent(++--+)
}

// =========================W+ -> K+ K- pi+ current ==================================================

EvtVector4C EvtWHad::WCurrent_KKP( const EvtVector4R& pKplus,
                                   const EvtVector4R& pKminus,
                                   const EvtVector4R& pPiPlus ) const
{
    const double mA1( 1.239 ), gammaA1( 0.600 );

    const EvtVector4R q = pKplus + pKminus + pPiPlus;
    const double q2 = q.mass2();
    const EvtVector4R pK = pKminus + pPiPlus;
    const double pK2 = pK.mass2();

    const EvtComplex I( 0.0, 1.0 );

    const EvtComplex den1 = 1.0 / ( q2 - mA1 * mA1 + I * mA1 * gammaA1 );
    const EvtComplex den2 = 1.0 / ( pK2 - mK_[0] * mK_[0] +
                                    I * mK_[0] * gammaK_[0] );    //K(892)

    const EvtTensor4C ten = EvtTensor4C::g() -
                            ( 1.0 / q2 ) * EvtGenFunctions::directProd( q, q );

    EvtVector4C vec = den1 * den2 * ( pKminus - pPiPlus );
    vec = ten.cont2( vec );

    return vec;
}

// hadronic hurrent W -> K+ K- pi+ pi+ pi- with identical pi+ symmetry
EvtVector4C EvtWHad::WCurrent_KKPPP( const EvtVector4R& pKplus,
                                     const EvtVector4R& pKminus,
                                     const EvtVector4R& pPi1Plus,
                                     const EvtVector4R& pPi2Plus,
                                     const EvtVector4R& pPiMinus ) const
{
    return EvtWHad::WCurrent_KKPPP_nosym( pKplus, pKminus, pPi1Plus, pPi2Plus,
                                          pPiMinus ) +
           EvtWHad::WCurrent_KKPPP_nosym( pKplus, pKminus, pPi2Plus, pPi1Plus,
                                          pPiMinus );
}

// hadronic hurrent W -> a1(K+ K- pi1+) f0(pi2+ pi-) without identical pi+ symmetry
EvtVector4C EvtWHad::WCurrent_KKPPP_nosym( const EvtVector4R& pKplus,
                                           const EvtVector4R& pKminus,
                                           const EvtVector4R& pPi1Plus,
                                           const EvtVector4R& pPi2Plus,
                                           const EvtVector4R& pPiMinus ) const
{
    const EvtVector4R pf0 = pPi2Plus + pPiMinus;
    const EvtVector4C epsA1 = EvtWHad::WCurrent_KKP( pKplus, pKminus, pPi1Plus );
    const EvtVector4R q = pKplus + pKminus + pPi1Plus + pPi2Plus + pPiMinus;
    return BWa( q ) * epsA1 * BWf( pf0 );
}

// 1=pi+ 2=pi+ 3=pi+ 4=pi+ 5=pi- 6=pi- 7=pi- with symmetrization of the identical particles
EvtVector4C EvtWHad::WCurrent_7pi( const EvtVector4R& p1, const EvtVector4R& p2,
                                   const EvtVector4R& p3, const EvtVector4R& p4,
                                   const EvtVector4R& p5, const EvtVector4R& p6,
                                   const EvtVector4R& p7 ) const
{
    // a1 -> a1(1=pi+ 2=pi+ 3=pi+ 5=pi- 6=pi-) f0(4=pi+ 7=pi-) without symmetrization of the identical particles
    // making p4 symmetric with p1, p2, p3
    //        p7                p5, p6
    EvtVector4C eps;
    eps += EvtWHad::WCurrent_7pi_nosymm( p1, p2, p3, p4, p5, p6, p7 );
    eps += EvtWHad::WCurrent_7pi_nosymm( p1, p2, p4, p3, p5, p6, p7 );
    eps += EvtWHad::WCurrent_7pi_nosymm( p1, p4, p3, p2, p5, p6, p7 );
    eps += EvtWHad::WCurrent_7pi_nosymm( p4, p2, p3, p1, p5, p6, p7 );
    //
    eps += EvtWHad::WCurrent_7pi_nosymm( p1, p2, p3, p4, p5, p7, p6 );
    eps += EvtWHad::WCurrent_7pi_nosymm( p1, p2, p4, p3, p5, p7, p6 );
    eps += EvtWHad::WCurrent_7pi_nosymm( p1, p4, p3, p2, p5, p7, p6 );
    eps += EvtWHad::WCurrent_7pi_nosymm( p4, p2, p3, p1, p5, p7, p6 );
    //
    eps += EvtWHad::WCurrent_7pi_nosymm( p1, p2, p3, p4, p7, p6, p5 );
    eps += EvtWHad::WCurrent_7pi_nosymm( p1, p2, p4, p3, p7, p6, p5 );
    eps += EvtWHad::WCurrent_7pi_nosymm( p1, p4, p3, p2, p7, p6, p5 );
    eps += EvtWHad::WCurrent_7pi_nosymm( p4, p2, p3, p1, p7, p6, p5 );

    return eps;
}

// a1 -> a1(1=pi+ 2=pi+ 3=pi+ 5=pi- 6=pi-) f0(4=pi+ 7=pi-) without symmetrization of the identical particles
EvtVector4C EvtWHad::WCurrent_7pi_nosymm(
    const EvtVector4R& p1, const EvtVector4R& p2, const EvtVector4R& p3,
    const EvtVector4R& p4, const EvtVector4R& p5, const EvtVector4R& p6,
    const EvtVector4R& p7 ) const
{
    const EvtVector4R qTot = p1 + p2 + p3 + p4 + p5 + p6 + p7;
    const EvtVector4C eps1 = EvtWHad::WCurrent_5pi( p1, p2, p3, p5,
                                                    p6 );    // pi+ pi+ pi+ pi- pi-
    const EvtVector4R pf0 = p4 + p7;
    return eps1 * BWa( qTot ) * BWf( pf0 );
};

// hadronic current W+ -> K+ pi+ pi-

EvtVector4C EvtWHad::WCurrent_KPP( const EvtVector4R& pKplus,
                                   const EvtVector4R& pPiPlus,
                                   const EvtVector4R& pPiMinus ) const
{
    const double cK1p = 0.210709, cK1r = -0.0152997, cK2p = 0.0945309,
                 cK2r = 0.504315;
    const double gRho_PiPi = 6.02;

    const EvtVector4R q = pKplus + pPiPlus + pPiMinus;
    const double q2 = q.mass2();

    double pp2( 0.0 );
    EvtVector4C curr( 0, 0, 0, 0 ), curr1;

    // W+ -> K1+(1270) -> K+ rho0 -> K+ pi+ pi-

    pp2 = ( pPiPlus + pPiMinus ).mass2();
    curr1 = ( pPiPlus - pPiMinus ) * Den( q2, mK_[1], gammaK_[1], gKRho_[1] ) *
            Den( pp2, mRho_[0], gamma0_[0], gRho_PiPi );    //K1(1270) and rho(770)
    curr = curr + cK1r * curr1;

    // W+ -> K1+(1270) -> K*(892)0 pi+ -> K+ pi- pi-

    pp2 = ( pKplus + pPiMinus ).mass2();
    curr1 = ( pKplus - pPiMinus ) * Den( q2, mK_[1], gammaK_[1], gKPi_[1] ) *
            Den( pp2, mK_[0], gammaK_[0], gKPi_[0] );    //K1(1270) and K(892)
    curr = curr + cK1p * curr1;

    // W+ -> K1+(1400) -> K+ rho0 -> K+ pi+ pi-

    pp2 = ( pPiMinus + pPiPlus ).mass2();
    curr1 = ( pPiPlus - pPiMinus ) * Den( q2, mK_[2], gammaK_[2], gKRho_[2] ) *
            Den( pp2, mRho_[0], gamma0_[0], gRho_PiPi );    //K1(1400) and rho(770)
    curr = curr + cK2r * curr1;

    // W+ -> K1+(1400) -> K*(892)0 pi+ -> K+ pi- pi+

    pp2 = ( pKplus + pPiMinus ).mass2();
    curr1 = ( pKplus - pPiPlus ) * Den( q2, mK_[2], gammaK_[2], gKPi_[2] ) *
            Den( pp2, mK_[0], gammaK_[0], gKPi_[0] );    //K1(1400) and K(892)
    curr = curr + cK2p * curr1;

    const EvtTensor4C ten = EvtTensor4C::g() -
                            ( 1.0 / q2 ) * EvtGenFunctions::directProd( q, q );
    curr = ten.cont2( curr );

    return curr;
}

EvtComplex EvtWHad::Den( double qSq, const double mR, const double gammaR,
                         const double gR ) const
{
    const EvtComplex I( 0.0, 1.0 );
    const EvtComplex tmp = qSq - mR * mR + I * mR * gammaR;

    if ( abs( tmp ) < 1e-10 )
        return 0;
    return gR / tmp;
}

// a1 -> pi+ pi+ pi- BW
EvtComplex EvtWHad::BWa( const EvtVector4R& q ) const
{
    const double _mA1( 1.26 ), _GA1( 0.4 );
    const double _mA1Sq = _mA1 * _mA1;
    const EvtComplex I( 0.0, 1.0 );
    const double Q2 = q.mass2();
    const double GA1 = _GA1 * pi3G( Q2 ) / pi3G( _mA1Sq );

    const EvtComplex denBA1( _mA1Sq - Q2, -1.0 * _mA1 * GA1 );
    if ( abs( denBA1 ) < 1e-10 )
        return 0;
    return _mA1Sq / denBA1;
}

EvtComplex EvtWHad::BWf( const EvtVector4R& q ) const
{
    const double mf( 0.8 ), Gf( 0.6 );
    const double mfSq = mf * mf;
    const EvtComplex I( 0.0, 1.0 );
    const double Q2 = q.mass2();
    return mfSq / ( mfSq - Q2 - I * mf * Gf );
}

EvtComplex EvtWHad::BWr( const EvtVector4R& q ) const
{
    const double beta( -0.108 );

    const double s = q.mass2();
    const EvtComplex BW_rho = BW( s, mRho_[0], gamma0_[0], mPi_, mPi_ );
    const EvtComplex BW_rhopr = BW( s, mRho_[1], gamma0_[1], mPi_, mPi_ );
    return ( BW_rho + beta * BW_rhopr ) / ( 1.0 + beta );
}

double EvtWHad::pi3G( double Q2 ) const
{
    const double mRhoPi = mRho_[0] + mPi_;
    // Parameterisation of scaling factor for a1 (to 3pi) decay width
    if ( Q2 < mRhoPi * mRhoPi ) {
        const double arg = Q2 - 9. * mPiSq_;
        const double arg2 = arg * arg;
        const double arg3 = arg * arg2;
        return 4.1 * arg3 * ( 1. - 3.3 * arg + 5.8 * arg2 );
    } else {
        const double Q2Sq = Q2 * Q2;
        const double Q2Cu = Q2 * Q2Sq;
        return Q2 * ( 1.623 + 10.38 / Q2 - 9.32 / Q2Sq + 0.65 / Q2Cu );
    }
}

EvtVector4C EvtWHad::JB( const EvtVector4R& p1, const EvtVector4R& p2,
                         const EvtVector4R& p3, const EvtVector4R& p4,
                         const EvtVector4R& p5 ) const
{
    const EvtVector4R Qtot = p1 + p2 + p3 + p4 + p5, Qa = p1 + p2 + p3;
    const EvtTensor4C T = ( 1.0 / Qtot.mass2() ) *
                              EvtGenFunctions::directProd( Qtot, Qtot ) -
                          EvtTensor4C::g();

    const EvtVector4R p13 = p1 - p3, p23 = p2 - p3;
    const EvtVector4R V13 = Qa * ( p2 * p13 ) / Qa.mass2() - p13;
    const EvtVector4R V23 = Qa * ( p1 * p23 ) / Qa.mass2() - p23;

    return BWa( Qtot ) * BWa( Qa ) * BWf( p4 + p5 ) *
           ( T.cont1( V13 ) * BWr( p1 + p3 ) + T.cont1( V23 ) * BWr( p2 + p3 ) );
}

EvtComplex EvtWHad::BW( double s, double m, double gamma, double xm1,
                        double xm2 ) const
{
    const double m2 = m * m;
    const double xmSum = xm1 + xm2;
    const double xmSumSq = xmSum * xmSum;
    if ( s > xmSumSq ) {
        const double xmDiff = xm1 - xm2;
        const double xmDiffSq = xmDiff * xmDiff;
        const double qs = sqrt( fabs( ( s - xmSumSq ) * ( s - xmDiffSq ) ) ) /
                          sqrt( s );
        const double qm = sqrt( fabs( ( m2 - xmSumSq ) * ( m2 - xmDiffSq ) ) ) /
                          m;

        const double qRatio = qm > 0.0 ? qs / qm : 0.0;
        gamma *= m2 / s * ( qRatio * qRatio * qRatio );
    } else
        gamma = 0.;
    const EvtComplex denBW( m2 - s, -1. * sqrt( s ) * gamma );
    return m2 / denBW;
}

EvtVector4C EvtWHad::WCurrent_K4pi( const EvtVector4R& p1, const EvtVector4R& p2,
                                    const EvtVector4R& p3, const EvtVector4R& p4,
                                    const EvtVector4R& p5 ) const
{
    return EvtWHad::WCurrent_K4pi_nosymm( p1, p2, p3, p4, p5 ) +
           EvtWHad::WCurrent_K4pi_nosymm( p1, p2, p3, p5, p4 );
}

// a1 -> K*0 (1=K+ 4=pi-) a1(2=pi+ 3=pi+ 5=pi-)
EvtVector4C EvtWHad::WCurrent_K4pi_nosymm( const EvtVector4R& p1,
                                           const EvtVector4R& p2,
                                           const EvtVector4R& p3,
                                           const EvtVector4R& p4,
                                           const EvtVector4R& p5 ) const
{
    const EvtComplex I( 0, 1 );
    const EvtVector4R pKstar = p1 + p4, pa1 = p2 + p3 + p5;
    EvtComplex denKstar = pKstar * pKstar - mK_[0] * mK_[0] +
                          I * mK_[0] * gammaK_[0];    //K(892)
    if ( abs( denKstar ) < 1e-10 ) {
        denKstar = 1e10;
    }
    const EvtVector4C epsKstar = 1 / denKstar * ( p1 - p4 );
    const EvtVector4C epsA1 = WCurrent( p2, p3, p5 );
    const EvtVector4C eps =
        dual( EvtGenFunctions::directProd( epsKstar, epsA1 ) ).cont2( pKstar - pa1 );
    return eps;
}
