
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

#include "EvtGenBase/EvtDalitzReso.hh"

#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtMatrix.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtdFunction.hh"

#include <assert.h>
#include <cmath>
#include <iostream>
#include <stdlib.h>

#define PRECISION ( 1.e-3 )

using EvtCyclic3::Index;
using EvtCyclic3::Pair;

// single Breit-Wigner
EvtDalitzReso::EvtDalitzReso( const EvtDalitzPlot& dp, Pair pairAng, Pair pairRes,
                              EvtSpinType::spintype spin, double m0, double g0,
                              NumType typeN, double f_b, double f_d ) :
    m_dp( dp ),
    m_pairAng( pairAng ),
    m_pairRes( pairRes ),
    m_spin( spin ),
    m_typeN( typeN ),
    m_m0( m0 ),
    m_g0( g0 ),
    m_massFirst( dp.m( first( pairRes ) ) ),
    m_massSecond( dp.m( second( pairRes ) ) ),
    m_m0_mix( -1. ),
    m_g0_mix( 0. ),
    m_delta_mix( 0. ),
    m_amp_mix( 0., 0. ),
    m_g1( -1. ),
    m_g2( -1. ),
    m_coupling2( Undefined ),
    m_f_b( f_b ),
    m_f_d( f_d ),
    m_kmatrix_index( -1 ),
    m_fr12prod( 0., 0. ),
    m_fr13prod( 0., 0. ),
    m_fr14prod( 0., 0. ),
    m_fr15prod( 0., 0. ),
    m_s0prod( 0. ),
    m_a( 0. ),
    m_r( 0. ),
    m_Blass( 0. ),
    m_phiB( 0. ),
    m_R( 0. ),
    m_phiR( 0. ),
    m_cutoff( -1. ),
    m_scaleByMOverQ( false ),
    m_alpha( 0. )
{
    assert( m_typeN != K_MATRIX && m_typeN != K_MATRIX_I &&
            m_typeN != K_MATRIX_II );    // single BW cannot be K-matrix
}

// Breit-Wigner with electromagnetic mass mixing
EvtDalitzReso::EvtDalitzReso( const EvtDalitzPlot& dp, Pair pairAng, Pair pairRes,
                              EvtSpinType::spintype spin, double m0, double g0,
                              NumType typeN, double m0_mix, double g0_mix,
                              double delta_mix, EvtComplex amp_mix ) :
    m_dp( dp ),
    m_pairAng( pairAng ),
    m_pairRes( pairRes ),
    m_spin( spin ),
    m_typeN( typeN ),
    m_m0( m0 ),
    m_g0( g0 ),
    m_massFirst( dp.m( first( pairRes ) ) ),
    m_massSecond( dp.m( second( pairRes ) ) ),
    m_m0_mix( m0_mix ),
    m_g0_mix( g0_mix ),
    m_delta_mix( delta_mix ),
    m_amp_mix( amp_mix ),
    m_g1( -1. ),
    m_g2( -1. ),
    m_coupling2( Undefined ),
    m_f_b( 0.0 ),
    m_f_d( 1.5 ),
    m_kmatrix_index( -1 ),
    m_fr12prod( 0., 0. ),
    m_fr13prod( 0., 0. ),
    m_fr14prod( 0., 0. ),
    m_fr15prod( 0., 0. ),
    m_s0prod( 0. ),
    m_a( 0. ),
    m_r( 0. ),
    m_Blass( 0. ),
    m_phiB( 0. ),
    m_R( 0. ),
    m_phiR( 0. ),
    m_cutoff( -1. ),
    m_scaleByMOverQ( false ),
    m_alpha( 0. )
{
    // single BW (with electromagnetic mixing) cannot be K-matrix
    assert( m_typeN != K_MATRIX && m_typeN != K_MATRIX_I &&
            m_typeN != K_MATRIX_II );
}

// coupled Breit-Wigner
EvtDalitzReso::EvtDalitzReso( const EvtDalitzPlot& dp, Pair pairAng,
                              Pair pairRes, EvtSpinType::spintype spin,
                              double m0, NumType typeN, double g1, double g2,
                              CouplingType coupling2 ) :
    m_dp( dp ),
    m_pairAng( pairAng ),
    m_pairRes( pairRes ),
    m_spin( spin ),
    m_typeN( typeN ),
    m_m0( m0 ),
    m_g0( -1. ),
    m_massFirst( dp.m( first( pairRes ) ) ),
    m_massSecond( dp.m( second( pairRes ) ) ),
    m_m0_mix( -1. ),
    m_g0_mix( 0. ),
    m_delta_mix( 0. ),
    m_amp_mix( 0., 0. ),
    m_g1( g1 ),
    m_g2( g2 ),
    m_coupling2( coupling2 ),
    m_f_b( 0.0 ),
    m_f_d( 1.5 ),
    m_kmatrix_index( -1 ),
    m_fr12prod( 0., 0. ),
    m_fr13prod( 0., 0. ),
    m_fr14prod( 0., 0. ),
    m_fr15prod( 0., 0. ),
    m_s0prod( 0. ),
    m_a( 0. ),
    m_r( 0. ),
    m_Blass( 0. ),
    m_phiB( 0. ),
    m_R( 0. ),
    m_phiR( 0. ),
    m_cutoff( -1. ),
    m_scaleByMOverQ( false ),
    m_alpha( 0. )
{
    assert( m_coupling2 != Undefined );
    assert( m_typeN != K_MATRIX && m_typeN != K_MATRIX_I &&
            m_typeN != K_MATRIX_II );    // coupled BW cannot be K-matrix
    assert( m_typeN != LASS );           // coupled BW cannot be LASS
    assert( m_typeN != NBW );            // for coupled BW, only relativistic BW
}

// K-Matrix (A&S)
EvtDalitzReso::EvtDalitzReso( const EvtDalitzPlot& dp, Pair pairRes,
                              std::string nameIndex, NumType typeN,
                              EvtComplex fr12prod, EvtComplex fr13prod,
                              EvtComplex fr14prod, EvtComplex fr15prod,
                              double s0prod ) :
    m_dp( dp ),
    m_pairRes( pairRes ),
    m_typeN( typeN ),
    m_m0( 0. ),
    m_g0( 0. ),
    m_massFirst( dp.m( first( pairRes ) ) ),
    m_massSecond( dp.m( second( pairRes ) ) ),
    m_m0_mix( -1. ),
    m_g0_mix( 0. ),
    m_delta_mix( 0. ),
    m_amp_mix( 0., 0. ),
    m_g1( -1. ),
    m_g2( -1. ),
    m_coupling2( Undefined ),
    m_f_b( 0. ),
    m_f_d( 0. ),
    m_kmatrix_index( -1 ),
    m_fr12prod( fr12prod ),
    m_fr13prod( fr13prod ),
    m_fr14prod( fr14prod ),
    m_fr15prod( fr15prod ),
    m_s0prod( s0prod ),
    m_a( 0. ),
    m_r( 0. ),
    m_Blass( 0. ),
    m_phiB( 0. ),
    m_R( 0. ),
    m_phiR( 0. ),
    m_cutoff( -1. ),
    m_scaleByMOverQ( false ),
    m_alpha( 0. )
{
    assert( m_typeN == K_MATRIX || m_typeN == K_MATRIX_I ||
            m_typeN == K_MATRIX_II );
    m_spin = EvtSpinType::SCALAR;
    if ( nameIndex == "Pole1" )
        m_kmatrix_index = 1;
    else if ( nameIndex == "Pole2" )
        m_kmatrix_index = 2;
    else if ( nameIndex == "Pole3" )
        m_kmatrix_index = 3;
    else if ( nameIndex == "Pole4" )
        m_kmatrix_index = 4;
    else if ( nameIndex == "Pole5" )
        m_kmatrix_index = 5;
    else if ( nameIndex == "f11prod" )
        m_kmatrix_index = 6;
    else
        assert( 0 );
}

// LASS parameterization
EvtDalitzReso::EvtDalitzReso( const EvtDalitzPlot& dp, Pair pairRes, double m0,
                              double g0, double a, double r, double B,
                              double phiB, double R, double phiR, double cutoff,
                              bool scaleByMOverQ ) :
    m_dp( dp ),
    m_pairRes( pairRes ),
    m_typeN( LASS ),
    m_m0( m0 ),
    m_g0( g0 ),
    m_massFirst( dp.m( first( pairRes ) ) ),
    m_massSecond( dp.m( second( pairRes ) ) ),
    m_m0_mix( -1. ),
    m_g0_mix( 0. ),
    m_delta_mix( 0. ),
    m_amp_mix( 0., 0. ),
    m_g1( -1. ),
    m_g2( -1. ),
    m_coupling2( Undefined ),
    m_f_b( 0.0 ),
    m_f_d( 1.5 ),
    m_kmatrix_index( -1 ),
    m_fr12prod( 0., 0. ),
    m_fr13prod( 0., 0. ),
    m_fr14prod( 0., 0. ),
    m_fr15prod( 0., 0. ),
    m_s0prod( 0. ),
    m_a( a ),
    m_r( r ),
    m_Blass( B ),
    m_phiB( phiB ),
    m_R( R ),
    m_phiR( phiR ),
    m_cutoff( cutoff ),
    m_scaleByMOverQ( scaleByMOverQ ),
    m_alpha( 0. )
{
    m_spin = EvtSpinType::SCALAR;
}

//Flatte
EvtDalitzReso::EvtDalitzReso( const EvtDalitzPlot& dp, EvtCyclic3::Pair pairRes,
                              double m0 ) :
    m_dp( dp ),
    m_pairRes( pairRes ),
    m_typeN( FLATTE ),
    m_m0( m0 ),
    m_g0( 0. ),
    m_massFirst( dp.m( first( pairRes ) ) ),
    m_massSecond( dp.m( second( pairRes ) ) ),
    m_m0_mix( -1. ),
    m_g0_mix( 0. ),
    m_delta_mix( 0. ),
    m_amp_mix( 0., 0. ),
    m_g1( -1. ),
    m_g2( -1. ),
    m_coupling2( Undefined ),
    m_f_b( 0. ),
    m_f_d( 0. ),
    m_kmatrix_index( -1 ),
    m_fr12prod( 0., 0. ),
    m_fr13prod( 0., 0. ),
    m_fr14prod( 0., 0. ),
    m_fr15prod( 0., 0. ),
    m_s0prod( 0. ),
    m_a( 0. ),
    m_r( 0. ),
    m_Blass( 0. ),
    m_phiB( 0. ),
    m_R( 0. ),
    m_phiR( 0. ),
    m_cutoff( -1. ),
    m_scaleByMOverQ( false ),
    m_alpha( 0. )
{
    m_spin = EvtSpinType::SCALAR;
}

EvtComplex EvtDalitzReso::evaluate( const EvtDalitzPoint& x ) const
{
    if ( m_typeN == NON_RES )
        return EvtComplex( 1.0, 0.0 );

    const double s = x.q( m_pairRes );
    const double m = sqrt( s );

    if ( m_typeN == NON_RES_LIN )
        return s;

    if ( m_typeN == NON_RES_EXP )
        return exp( -m_alpha * s );

    // do use always hash table (speed up fitting)
    if ( m_typeN == K_MATRIX || m_typeN == K_MATRIX_I || m_typeN == K_MATRIX_II )
        return Fvector( s, m_kmatrix_index );

    if ( m_typeN == FLATTE )
        return flatte( s );

    EvtTwoBodyVertex vd( m_massFirst, m_massSecond, m_m0, m_spin );
    vd.set_f( m_f_d );

    const EvtTwoBodyKine kd( m_massFirst, m_massSecond, m );

    if ( m_typeN == LASS )
        return lass( kd, vd );

    EvtComplex amp( 1.0, 0.0 );

    EvtTwoBodyVertex vb( m_m0, m_dp.m( EvtCyclic3::other( m_pairRes ) ),
                         x.bigM(), m_spin );
    vb.set_f( m_f_b );

    const EvtTwoBodyKine kb( m, x.m( EvtCyclic3::other( m_pairRes ) ), x.bigM() );

    EvtComplex prop( 0, 0 );
    if ( m_typeN == NBW ) {
        prop = propBreitWigner( m_m0, m_g0, m );
    } else if ( m_typeN == GAUSS_CLEO || m_typeN == GAUSS_CLEO_ZEMACH ) {
        prop = propGauss( m_m0, m_g0, m );
    } else {
        if ( m_coupling2 == Undefined ) {
            // single BW
            const double g = ( m_g0 <= 0. || vd.pD() <= 0. )
                                 ? -m_g0
                                 : m_g0 * vd.widthFactor( kd );    // running width
            if ( m_typeN == GS_CLEO || m_typeN == GS_CLEO_ZEMACH ) {
                // Gounaris-Sakurai (GS)
                prop = propGounarisSakurai( m_m0, fabs( m_g0 ), vd.pD(), m, g,
                                            kd.p() );
            } else {
                // standard relativistic BW
                prop = propBreitWignerRel( m_m0, g, m );
            }
        } else {
            // coupled width BW
            EvtComplex G1, G2;
            switch ( m_coupling2 ) {
                case PicPic: {
                    G1 = m_g1 * m_g1 * psFactor( m_massFirst, m_massSecond, m );
                    static const double mPic = EvtPDL::getMass(
                        EvtPDL::getId( "pi+" ) );
                    G2 = m_g2 * m_g2 * psFactor( mPic, mPic, m );
                    break;
                }
                case PizPiz: {
                    G1 = m_g1 * m_g1 * psFactor( m_massFirst, m_massSecond, m );
                    static const double mPiz = EvtPDL::getMass(
                        EvtPDL::getId( "pi0" ) );
                    G2 = m_g2 * m_g2 * psFactor( mPiz, mPiz, m );
                    break;
                }
                case PiPi: {
                    G1 = m_g1 * m_g1 * psFactor( m_massFirst, m_massSecond, m );
                    static const double mPic = EvtPDL::getMass(
                        EvtPDL::getId( "pi+" ) );
                    static const double mPiz = EvtPDL::getMass(
                        EvtPDL::getId( "pi0" ) );
                    G2 = m_g2 * m_g2 * psFactor( mPic, mPic, mPiz, mPiz, m );
                    break;
                }
                case KcKc: {
                    G1 = m_g1 * m_g1 * psFactor( m_massFirst, m_massSecond, m );
                    static const double mKc = EvtPDL::getMass(
                        EvtPDL::getId( "K+" ) );
                    G2 = m_g2 * m_g2 * psFactor( mKc, mKc, m );
                    break;
                }
                case KzKz: {
                    G1 = m_g1 * m_g1 * psFactor( m_massFirst, m_massSecond, m );
                    static const double mKz = EvtPDL::getMass(
                        EvtPDL::getId( "K0" ) );
                    G2 = m_g2 * m_g2 * psFactor( mKz, mKz, m );
                    break;
                }
                case KK: {
                    G1 = m_g1 * m_g1 * psFactor( m_massFirst, m_massSecond, m );
                    static const double mKc = EvtPDL::getMass(
                        EvtPDL::getId( "K+" ) );
                    static const double mKz = EvtPDL::getMass(
                        EvtPDL::getId( "K0" ) );
                    G2 = m_g2 * m_g2 * psFactor( mKc, mKc, mKz, mKz, m );
                    break;
                }
                case EtaPic: {
                    G1 = m_g1 * m_g1 * psFactor( m_massFirst, m_massSecond, m );
                    static const double mEta = EvtPDL::getMass(
                        EvtPDL::getId( "eta" ) );
                    static const double mPic = EvtPDL::getMass(
                        EvtPDL::getId( "pi+" ) );
                    G2 = m_g2 * m_g2 * psFactor( mEta, mPic, m );
                    break;
                }
                case EtaPiz: {
                    G1 = m_g1 * m_g1 * psFactor( m_massFirst, m_massSecond, m );
                    static const double mEta = EvtPDL::getMass(
                        EvtPDL::getId( "eta" ) );
                    static const double mPiz = EvtPDL::getMass(
                        EvtPDL::getId( "pi0" ) );
                    G2 = m_g2 * m_g2 * psFactor( mEta, mPiz, m );
                    break;
                }
                case PicPicKK: {
                    static const double mPic = EvtPDL::getMass(
                        EvtPDL::getId( "pi+" ) );
                    G1 = m_g1 * psFactor( mPic, mPic, m );
                    static const double mKc = EvtPDL::getMass(
                        EvtPDL::getId( "K+" ) );
                    static const double mKz = EvtPDL::getMass(
                        EvtPDL::getId( "K0" ) );
                    G2 = m_g2 * psFactor( mKc, mKc, mKz, mKz, m );
                    break;
                }
                default:
                    std::cout
                        << "EvtDalitzReso:evaluate(): PANIC, wrong coupling2 state."
                        << std::endl;
                    assert( 0 );
                    break;
            }
            // calculate standard couple BW propagator
            if ( m_coupling2 != WA76 )
                prop = m_g1 * propBreitWignerRelCoupled( m_m0, G1, G2, m );
        }
    }
    amp *= prop;

    // Compute form-factors (Blatt-Weisskopf penetration factor)
    amp *= vb.formFactor( kb );
    amp *= vd.formFactor( kd );

    // Compute numerator (angular distribution)
    amp *= numerator( x, vb, vd, kb, kd );

    // Compute electromagnetic mass mixing factor
    if ( m_m0_mix > 0. ) {
        EvtComplex prop_mix;
        if ( m_typeN == NBW ) {
            prop_mix = propBreitWigner( m_m0_mix, m_g0_mix, m );
        } else {
            assert( m_g1 < 0. );    // running width only
            const double g_mix = m_g0_mix * vd.widthFactor( kd );
            prop_mix = propBreitWignerRel( m_m0_mix, g_mix, m );
        }
        amp *= mixFactor( prop, prop_mix );
    }

    return amp;
}

EvtComplex EvtDalitzReso::psFactor( const double ma, const double mb,
                                    const double m ) const
{
    if ( m > ( ma + mb ) ) {
        const EvtTwoBodyKine kd( ma, mb, m );
        return EvtComplex( 0, 2 * kd.p() / m );
    } else {
        // analytical continuation
        const double s = m * m;
        const double phaseFactor_analyticalCont =
            -0.5 * ( sqrt( 4 * ma * ma / s - 1 ) + sqrt( 4 * mb * mb / s - 1 ) );
        return EvtComplex( phaseFactor_analyticalCont, 0 );
    }
}

EvtComplex EvtDalitzReso::psFactor( const double ma1, const double mb1,
                                    const double ma2, const double mb2,
                                    const double m ) const
{
    return 0.5 * ( psFactor( ma1, mb1, m ) + psFactor( ma2, mb2, m ) );
}

EvtComplex EvtDalitzReso::propGauss( const double m0, const double s0,
                                     const double m ) const
{
    // Gaussian
    const double gauss = 1. / sqrt( EvtConst::twoPi ) / s0 *
                         exp( -( m - m0 ) * ( m - m0 ) / 2. / ( s0 * s0 ) );
    return EvtComplex( gauss, 0. );
}

EvtComplex EvtDalitzReso::propBreitWigner( const double m0, const double g0,
                                           const double m ) const
{
    // non-relativistic BW
    return sqrt( g0 / EvtConst::twoPi ) / ( m - m0 - EvtComplex( 0.0, g0 / 2. ) );
}

EvtComplex EvtDalitzReso::propBreitWignerRel( const double m0, const double g0,
                                              const double m ) const
{
    // relativistic BW with real width
    return 1. / ( m0 * m0 - m * m - EvtComplex( 0., m0 * g0 ) );
}

EvtComplex EvtDalitzReso::propBreitWignerRel( const double m0,
                                              const EvtComplex& g0,
                                              const double m ) const
{
    // relativistic BW with complex width
    return 1. / ( m0 * m0 - m * m - EvtComplex( 0., m0 ) * g0 );
}

EvtComplex EvtDalitzReso::propBreitWignerRelCoupled( const double m0,
                                                     const EvtComplex& g1,
                                                     const EvtComplex& g2,
                                                     const double m ) const
{
    // relativistic coupled BW
    return 1. / ( m0 * m0 - m * m - ( g1 + g2 ) );
}

EvtComplex EvtDalitzReso::propGounarisSakurai( const double m0, const double g0,
                                               const double k0, const double m,
                                               const double g,
                                               const double k ) const
{
    // Gounaris-Sakurai parameterization of pi+pi- P wave. PRD, Vol61, 112002. PRL, Vol21, 244.
    // Expressions taken from BAD637v4, after fixing the imaginary part of the BW denominator: i M_R Gamma_R(s) --> i sqrt(s) Gamma_R(s)
    return ( 1. + GS_d( m0, k0 ) * g0 / m0 ) /
           ( m0 * m0 - m * m - EvtComplex( 0., m * g ) +
             GS_f( m0, g0, k0, m, k ) );
}

inline double EvtDalitzReso::GS_f( const double m0, const double g0,
                                   const double k0, const double m,
                                   const double k ) const
{
    // m: sqrt(s)
    // m0: nominal resonance mass
    // k: momentum of pion in resonance rest frame (at m)
    // k0: momentum of pion in resonance rest frame (at nominal resonance mass)
    return g0 * m0 * m0 / ( k0 * k0 * k0 ) *
           ( k * k * ( GS_h( m, k ) - GS_h( m0, k0 ) ) +
             ( m0 * m0 - m * m ) * k0 * k0 * GS_dhods( m0, k0 ) );
}

inline double EvtDalitzReso::GS_h( const double m, const double k ) const
{
    return 2. / EvtConst::pi * k / m *
           log( ( m + 2. * k ) / ( 2. * m_massFirst ) );
}

inline double EvtDalitzReso::GS_dhods( const double m0, const double k0 ) const
{
    return GS_h( m0, k0 ) * ( 0.125 / ( k0 * k0 ) - 0.5 / ( m0 * m0 ) ) +
           0.5 / ( EvtConst::pi * m0 * m0 );
}

inline double EvtDalitzReso::GS_d( const double m0, const double k0 ) const
{
    return 3. / EvtConst::pi * m_massFirst * m_massFirst / ( k0 * k0 ) *
               log( ( m0 + 2. * k0 ) / ( 2. * m_massFirst ) ) +
           m0 / ( 2. * EvtConst::pi * k0 ) -
           m_massFirst * m_massFirst * m0 / ( EvtConst::pi * k0 * k0 * k0 );
}

EvtComplex EvtDalitzReso::numerator( const EvtDalitzPoint& x,
                                     const EvtTwoBodyVertex& vb,
                                     const EvtTwoBodyVertex& vd,
                                     const EvtTwoBodyKine& kb,
                                     const EvtTwoBodyKine& kd ) const
{
    EvtComplex ret( 0., 0. );

    // Non-relativistic Breit-Wigner
    if ( NBW == m_typeN ) {
        ret = angDep( x );
    }

    // Standard relativistic Zemach propagator
    else if ( RBW_ZEMACH == m_typeN ) {
        ret = vd.phaseSpaceFactor( kd, EvtTwoBodyKine::AB ) * angDep( x );
    }

    // Standard relativistic Zemach propagator
    else if ( RBW_ZEMACH2 == m_typeN ) {
        ret = vd.phaseSpaceFactor( kd, EvtTwoBodyKine::AB ) *
              vb.phaseSpaceFactor( kb, EvtTwoBodyKine::AB ) * angDep( x );
        if ( m_spin == EvtSpinType::VECTOR ) {
            ret *= -4.;
        } else if ( m_spin == EvtSpinType::TENSOR ) {
            ret *= 16. / 3.;
        } else if ( m_spin != EvtSpinType::SCALAR )
            assert( 0 );
    }

    // Kuehn-Santamaria normalization:
    else if ( RBW_KUEHN == m_typeN ) {
        ret = m_m0 * m_m0 * angDep( x );
    }

    // CLEO amplitude
    else if ( ( RBW_CLEO == m_typeN ) || ( GS_CLEO == m_typeN ) ||
              ( RBW_CLEO_ZEMACH == m_typeN ) || ( GS_CLEO_ZEMACH == m_typeN ) ||
              ( GAUSS_CLEO == m_typeN ) || ( GAUSS_CLEO_ZEMACH == m_typeN ) ) {
        const Index iA = other( m_pairAng );                // A = other(BC)
        const Index iB = common( m_pairRes, m_pairAng );    // B = common(AB,BC)
        const Index iC = other( m_pairRes );                // C = other(AB)

        const double M = x.bigM();
        const double mA = x.m( iA );
        const double mB = x.m( iB );
        const double mC = x.m( iC );
        const double qAB = x.q( combine( iA, iB ) );
        const double qBC = x.q( combine( iB, iC ) );
        const double qCA = x.q( combine( iC, iA ) );

        const double M2 = M * M;
        const double m02 = ( ( RBW_CLEO_ZEMACH == m_typeN ) ||
                             ( GS_CLEO_ZEMACH == m_typeN ) ||
                             ( GAUSS_CLEO_ZEMACH == m_typeN ) )
                               ? qAB
                               : m_m0 * m_m0;
        const double mA2 = mA * mA;
        const double mB2 = mB * mB;
        const double mC2 = mC * mC;

        if ( m_spin == EvtSpinType::SCALAR )
            ret = EvtComplex( 1., 0. );
        else if ( m_spin == EvtSpinType::VECTOR ) {
            ret = qCA - qBC + ( M2 - mC2 ) * ( mB2 - mA2 ) / m02;
        } else if ( m_spin == EvtSpinType::TENSOR ) {
            const double x1 = qBC - qCA + ( M2 - mC2 ) * ( mA2 - mB2 ) / m02;
            const double x2 = M2 - mC2;
            const double x3 = qAB - 2 * M2 - 2 * mC2 + x2 * x2 / m02;
            const double x4 = mA2 - mB2;
            const double x5 = qAB - 2 * mB2 - 2 * mA2 + x4 * x4 / m02;
            ret = x1 * x1 - x3 * x5 / 3.;
        } else
            assert( 0 );
    }

    return ret;
}

double EvtDalitzReso::angDep( const EvtDalitzPoint& x ) const
{
    // Angular dependece for factorizable amplitudes
    // unphysical cosines indicate we are in big trouble
    const double cosTh = x.cosTh(
        m_pairAng, m_pairRes );    // angle between common(reso,ang) and other(reso)
    if ( fabs( cosTh ) > 1. ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "cosTh " << cosTh << std::endl;
        assert( 0 );
    }

    // in units of half-spin
    return EvtdFunction::d( EvtSpinType::getSpin2( m_spin ), 2 * 0, 2 * 0,
                            acos( cosTh ) );
}

EvtComplex EvtDalitzReso::mixFactor( const EvtComplex& prop,
                                     const EvtComplex& prop_mix ) const
{
    const double Delta = m_delta_mix * ( m_m0 + m_m0_mix );
    return 1 / ( 1 - Delta * Delta * prop * prop_mix ) *
           ( 1 + m_amp_mix * Delta * prop_mix );
}

EvtComplex EvtDalitzReso::Fvector( const double s, const int index ) const
{
    assert( index >= 1 && index <= 6 );

    //Define the complex coupling constant
    //The convection is as follow
    //i=0 --> pi+ pi-
    //i=1 --> KK
    //i=2 --> 4pi
    //i=3 --> eta eta
    //i=4 --> eta eta'
    //The first index is the resonace-pole index

    double g[5][5];    // Coupling constants. The first index is the pole index. The second index is the decay channel
    double ma[5];      // Pole masses. The unit is in GeV

    const int solution = ( m_typeN == K_MATRIX )
                             ? 3
                             : ( ( m_typeN == K_MATRIX_I )
                                     ? 1
                                     : ( ( m_typeN == K_MATRIX_II ) ? 2 : 0 ) );
    if ( solution == 0 ) {
        std::cout << "EvtDalitzReso::Fvector() error. Kmatrix solution incorrectly chosen ! "
                  << std::endl;
        abort();
    }

    if ( solution == 3 ) {
        // coupling constants
        //pi+pi- channel
        g[0][0] = 0.22889;
        g[1][0] = 0.94128;
        g[2][0] = 0.36856;
        g[3][0] = 0.33650;
        g[4][0] = 0.18171;
        //K+K- channel
        g[0][1] = -0.55377;
        g[1][1] = 0.55095;
        g[2][1] = 0.23888;
        g[3][1] = 0.40907;
        g[4][1] = -0.17558;
        //4pi channel
        g[0][2] = 0;
        g[1][2] = 0;
        g[2][2] = 0.55639;
        g[3][2] = 0.85679;
        g[4][2] = -0.79658;
        //eta eta channel
        g[0][3] = -0.39899;
        g[1][3] = 0.39065;
        g[2][3] = 0.18340;
        g[3][3] = 0.19906;
        g[4][3] = -0.00355;
        //eta eta' channel
        g[0][4] = -0.34639;
        g[1][4] = 0.31503;
        g[2][4] = 0.18681;
        g[3][4] = -0.00984;
        g[4][4] = 0.22358;

        // Pole masses
        ma[0] = 0.651;
        ma[1] = 1.20360;
        ma[2] = 1.55817;
        ma[3] = 1.21000;
        ma[4] = 1.82206;

    } else if ( solution == 1 ) {    // solnI.txt

        // coupling constants
        //pi+pi- channel
        g[0][0] = 0.31896;
        g[1][0] = 0.85963;
        g[2][0] = 0.47993;
        g[3][0] = 0.45121;
        g[4][0] = 0.39391;
        //K+K- channel
        g[0][1] = -0.49998;
        g[1][1] = 0.52402;
        g[2][1] = 0.40254;
        g[3][1] = 0.42769;
        g[4][1] = -0.30860;
        //4pi channel
        g[0][2] = 0;
        g[1][2] = 0;
        g[2][2] = 1.0;
        g[3][2] = 1.15088;
        g[4][2] = 0.33999;
        //eta eta channel
        g[0][3] = -0.21554;
        g[1][3] = 0.38093;
        g[2][3] = 0.21811;
        g[3][3] = 0.22925;
        g[4][3] = 0.06919;
        //eta eta' channel
        g[0][4] = -0.18294;
        g[1][4] = 0.23788;
        g[2][4] = 0.05454;
        g[3][4] = 0.06444;
        g[4][4] = 0.32620;

        // Pole masses
        ma[0] = 0.7369;
        ma[1] = 1.24347;
        ma[2] = 1.62681;
        ma[3] = 1.21900;
        ma[4] = 1.74932;

    } else if ( solution == 2 ) {    // solnIIa.txt

        // coupling constants
        //pi+pi- channel
        g[0][0] = 0.26014;
        g[1][0] = 0.95289;
        g[2][0] = 0.46244;
        g[3][0] = 0.41848;
        g[4][0] = 0.01804;
        //K+K- channel
        g[0][1] = -0.57849;
        g[1][1] = 0.55887;
        g[2][1] = 0.31712;
        g[3][1] = 0.49910;
        g[4][1] = -0.28430;
        //4pi channel
        g[0][2] = 0;
        g[1][2] = 0;
        g[2][2] = 0.70340;
        g[3][2] = 0.96819;
        g[4][2] = -0.90100;
        //eta eta channel
        g[0][3] = -0.32936;
        g[1][3] = 0.39910;
        g[2][3] = 0.22963;
        g[3][3] = 0.24415;
        g[4][3] = -0.07252;
        //eta eta' channel
        g[0][4] = -0.30906;
        g[1][4] = 0.31143;
        g[2][4] = 0.19802;
        g[3][4] = -0.00522;
        g[4][4] = 0.17097;

        // Pole masses
        ma[0] = 0.67460;
        ma[1] = 1.21094;
        ma[2] = 1.57896;
        ma[3] = 1.21900;
        ma[4] = 1.86602;
    }

    //Now define the K-matrix pole
    double rho1sq, rho2sq, rho4sq, rho5sq;
    EvtComplex rho[5];
    double f[5][5];

    //Initalize the mass of the resonance
    const double mpi = 0.13957;
    const double mK = 0.493677;      //using charged K value
    const double meta = 0.54775;     //using PDG value
    const double metap = 0.95778;    //using PDG value

    //Initialize the matrix to value zero
    EvtComplex K[5][5];
    for ( int i = 0; i < 5; i++ ) {
        for ( int j = 0; j < 5; j++ ) {
            K[i][j] = EvtComplex( 0, 0 );
            f[i][j] = 0;
        }
    }

    //Input the _f[i][j] scattering data
    double s_scatt = 0.0;
    if ( solution == 3 )
        s_scatt = -3.92637;
    else if ( solution == 1 )
        s_scatt = -5.0;
    else if ( solution == 2 )
        s_scatt = -5.0;
    const double sa = 1.0;
    const double sa_0 = -0.15;
    if ( solution == 3 ) {
        f[0][0] = 0.23399;    // f^scatt
        f[0][1] = 0.15044;
        f[0][2] = -0.20545;
        f[0][3] = 0.32825;
        f[0][4] = 0.35412;
    } else if ( solution == 1 ) {
        f[0][0] = 0.04214;    // f^scatt
        f[0][1] = 0.19865;
        f[0][2] = -0.63764;
        f[0][3] = 0.44063;
        f[0][4] = 0.36717;
    } else if ( solution == 2 ) {
        f[0][0] = 0.26447;    // f^scatt
        f[0][1] = 0.10400;
        f[0][2] = -0.35445;
        f[0][3] = 0.31596;
        f[0][4] = 0.42483;
    }
    f[1][0] = f[0][1];
    f[2][0] = f[0][2];
    f[3][0] = f[0][3];
    f[4][0] = f[0][4];

    //Now construct the phase-space factor
    //For eta-eta' there is no difference term
    rho1sq = 1. - pow( mpi + mpi, 2 ) / s;    //pi+ pi- phase factor
    if ( rho1sq >= 0 )
        rho[0] = EvtComplex( sqrt( rho1sq ), 0 );
    else
        rho[0] = EvtComplex( 0, sqrt( -rho1sq ) );

    rho2sq = 1. - pow( mK + mK, 2 ) / s;
    if ( rho2sq >= 0 )
        rho[1] = EvtComplex( sqrt( rho2sq ), 0 );
    else
        rho[1] = EvtComplex( 0, sqrt( -rho2sq ) );

    //using the A&S 4pi phase space Factor:
    //Shit, not continue
    if ( s <= 1 ) {
        const double real = 1.2274 + .00370909 / ( s * s ) - .111203 / s -
                            6.39017 * s + 16.8358 * s * s -
                            21.8845 * s * s * s + 11.3153 * s * s * s * s;
        const double cont32 = sqrt( 1.0 - ( 16.0 * mpi * mpi ) );
        rho[2] = EvtComplex( cont32 * real, 0 );
    } else
        rho[2] = EvtComplex( sqrt( 1. - 16. * mpi * mpi / s ), 0 );

    rho4sq = 1. - pow( meta + meta, 2 ) / s;
    if ( rho4sq >= 0 )
        rho[3] = EvtComplex( sqrt( rho4sq ), 0 );
    else
        rho[3] = EvtComplex( 0, sqrt( -rho4sq ) );

    rho5sq = 1. - pow( meta + metap, 2 ) / s;
    if ( rho5sq >= 0 )
        rho[4] = EvtComplex( sqrt( rho5sq ), 0 );
    else
        rho[4] = EvtComplex( 0, sqrt( -rho5sq ) );

    double smallTerm = 1;    // Factor to prevent divergences.

    // Check if some pole may arise problems.
    for ( int pole = 0; pole < 5; pole++ )
        if ( fabs( pow( ma[pole], 2 ) - s ) < PRECISION )
            smallTerm = pow( ma[pole], 2 ) - s;

    //now sum all the pole
    //equation (3) in the E791 K-matrix paper
    for ( int i = 0; i < 5; i++ ) {
        for ( int j = 0; j < 5; j++ ) {
            for ( int pole_index = 0; pole_index < 5; pole_index++ ) {
                const double A = g[pole_index][i] * g[pole_index][j];
                const double B = ma[pole_index] * ma[pole_index] - s;

                if ( fabs( B ) < PRECISION )
                    K[i][j] += EvtComplex( A, 0 );
                else
                    K[i][j] += EvtComplex( A / B, 0 ) * smallTerm;
            }
        }
    }

    //now add the SVT part
    for ( int i = 0; i < 5; i++ ) {
        for ( int j = 0; j < 5; j++ ) {
            const double C = f[i][j] * ( 1.0 - s_scatt );
            const double D = ( s - s_scatt );
            K[i][j] += EvtComplex( C / D, 0 ) * smallTerm;
        }
    }

    //Fix the bug in the FOCUS paper
    //Include the Alder zero term:
    for ( int i = 0; i < 5; i++ ) {
        for ( int j = 0; j < 5; j++ ) {
            const double E = ( s - ( sa * mpi * mpi * 0.5 ) ) * ( 1.0 - sa_0 );
            const double F = ( s - sa_0 );
            K[i][j] *= EvtComplex( E / F, 0 );
        }
    }

    //This is not correct!
    //(1-ipK) != (1-iKp)
    static thread_local EvtMatrix<EvtComplex> mat;
    mat.setRange(
        5 );    // Try to do in only the first time. DEFINE ALLOCATION IN CONSTRUCTOR.

    for ( int row = 0; row < 5; row++ )
        for ( int col = 0; col < 5; col++ )
            mat( row, col ) = ( row == col ) * smallTerm -
                              EvtComplex( 0., 1. ) * K[row][col] * rho[col];

    EvtMatrix<EvtComplex>* matInverse =
        mat.inverse();    //The 1st row of the inverse matrix. This matrix is {(I-iKp)^-1}_0j
    vector<EvtComplex> U1j;
    for ( int j = 0; j < 5; j++ )
        U1j.push_back( ( *matInverse )[0][j] );

    delete matInverse;

    //this calculates final F0 factor
    EvtComplex value( 0, 0 );
    if ( index <= 5 ) {
        //this calculates the beta_idx Factors
        for ( int j = 0; j < 5; j++ ) {    // sum for 5 channel
            const EvtComplex top = U1j[j] * g[index - 1][j];
            const double bottom = ma[index - 1] * ma[index - 1] - s;

            if ( fabs( bottom ) < PRECISION )
                value += top;
            else
                value += top / bottom * smallTerm;
        }
    } else {
        //this calculates fprod Factors
        value += U1j[0];
        value += U1j[1] * m_fr12prod;
        value += U1j[2] * m_fr13prod;
        value += U1j[3] * m_fr14prod;
        value += U1j[4] * m_fr15prod;

        value *= ( 1 - m_s0prod ) / ( s - m_s0prod ) * smallTerm;
    }

    return value;
}

//replace Breit-Wigner with LASS
EvtComplex EvtDalitzReso::lass( const EvtTwoBodyKine& kd,
                                const EvtTwoBodyVertex& vd ) const
{
    const double m = kd.m();
    const double q = kd.p();
    const double GammaM = m_g0 * vd.widthFactor( kd );    // running width;

    //calculate the background phase motion
    const double cot_deltaB = 1.0 / ( m_a * q ) + 0.5 * m_r * q;
    const double deltaB = atan( 1.0 / cot_deltaB );
    const double totalB = deltaB + m_phiB;

    //calculate the resonant phase motion
    const double deltaR = atan( ( m_m0 * GammaM / ( m_m0 * m_m0 - m * m ) ) );
    const double totalR = deltaR + m_phiR;

    //sum them up
    EvtComplex bkgB, resT;
    bkgB = EvtComplex( m_Blass * sin( totalB ), 0 ) *
           EvtComplex( cos( totalB ), sin( totalB ) );
    resT = EvtComplex( m_R * sin( deltaR ), 0 ) *
           EvtComplex( cos( totalR ), sin( totalR ) ) *
           EvtComplex( cos( 2 * totalB ), sin( 2 * totalB ) );

    EvtComplex T;
    if ( m_cutoff > 0 && m > m_cutoff )
        T = resT;
    else
        T = bkgB + resT;

    if ( m_scaleByMOverQ )
        T *= ( m / q );

    return T;
}

EvtComplex EvtDalitzReso::flatte( const double s ) const
{
    EvtComplex w;

    for ( const auto& param : m_flatteParams ) {
        const double m1 = param.m1();
        const double m2 = param.m2();
        const double g = param.g();
        w += ( g * g *
               sqrtCplx( ( 1 - ( ( m1 - m2 ) * ( m1 - m2 ) ) / ( s ) ) *
                         ( 1 - ( ( m1 + m2 ) * ( m1 + m2 ) ) / ( s ) ) ) );
    }

    const EvtComplex denom = m_m0 * m_m0 - s - EvtComplex( 0, 1 ) * w;

    return EvtComplex( 1.0, 0.0 ) / denom;
}
