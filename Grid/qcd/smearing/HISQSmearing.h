/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/smearing/HISQSmearing.h

Copyright (C) 2023

Author: D. A. Clarke <clarke.davida@gmail.com> 

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*
    @file HISQSmearing.h
    @brief Declares classes related to HISQ smearing 
*/


#pragma once
#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>

NAMESPACE_BEGIN(Grid);



/*!  @brief structure holding the link treatment */
template<class floatT>
struct SmearingParameters{
    SmearingParameters(){}
    floatT c_1;               // 1 link
    floatT c_naik;            // Naik term
    floatT c_3;               // 3 link
    floatT c_5;               // 5 link
    floatT c_7;               // 7 link
    floatT c_lp;              // 5 link Lepage
    SmearingParameters(floatT c1, floatT cnaik, floatT c3, floatT c5, floatT c7, floatT clp) 
        : c_1(c1),
          c_naik(cnaik),
          c_3(c3),
          c_5(c5),
          c_7(c7),
          c_lp(clp){}
};


// I think that "coalesced..." functions are extremely general, which is nice,
// but in the HISQ context it boils down to link reading and writing.
template<class vobj> accelerator_inline
vobj getLink(const vobj & __restrict__ vec,GeneralStencilEntry* SE) {
    return coalescedReadGeneralPermute(vec, SE->_permute, Nd); 
}
#define setLink coalescedWrite 


// figure out the stencil index from mu and nu
accelerator_inline int stencilIndex(int mu, int nu) {
    // Nshifts depends on how you built the stencil
    int Nshifts = 5;
    return Nshifts*nu + Nd*Nshifts*mu;
}



/*!  @brief create fat links from link variables */
template<class Gimpl> 
class Smear_HISQ : public Gimpl {
public:

    GridCartesian* const _grid;

    // Sort out the Gimpl. This handles BCs and part of the precision. 
    INHERIT_GIMPL_TYPES(Gimpl);
    typedef typename Gimpl::GaugeField     GF;
    typedef typename Gimpl::GaugeLinkField LF;
    typedef typename Gimpl::ComplexField   CF;
    typedef typename Gimpl::Scalar ComplexScalar;
    typedef decltype(real(ComplexScalar())) RealScalar;
    typedef iColourMatrix<ComplexScalar> ComplexColourMatrix;

    RealScalar Scut=-1; // Cutoff for U(3) projection eigenvalues, set at initialization
    int HaloDepth=1; 

    SmearingParameters<RealScalar> linkTreatment;

    void initialize() {
        if (sizeof(RealScalar)==4) {
            Scut=1e-5; // Maybe should be higher? e.g. 1e-4
        } else if (sizeof(RealScalar)==8) { 
            Scut=1e-8;
        } else {
            Grid_error("HISQ smearing only implemented for single and double");
        }
        assert(Nc == 3 && "HISQ smearing currently implemented only for Nc==3");
        assert(Nd == 4 && "HISQ smearing only defined for Nd==4");
    }

    Smear_HISQ(GridCartesian* grid, RealScalar c1, RealScalar cnaik, RealScalar c3, RealScalar c5, RealScalar c7, RealScalar clp) 
        : _grid(grid), 
          linkTreatment(c1,cnaik,c3,c5,c7,clp) {
        initialize();
    }

    // Allow to pass a pointer to a C-style array for MILC convenience
    Smear_HISQ(GridCartesian* grid, double* coeff) 
        : _grid(grid), 
          linkTreatment(coeff[0],coeff[1],coeff[2],coeff[3],coeff[4],coeff[5]) {
        initialize();
    }
    Smear_HISQ(GridCartesian* grid, float* coeff) 
        : _grid(grid), 
          linkTreatment(coeff[0],coeff[1],coeff[2],coeff[3],coeff[4],coeff[5]) {
        initialize();
    }

    ~Smear_HISQ() {}

    // Intent: OUT--u_smr (smeared links), 
    //              u_naik (Naik links),
    //          IN--u_thin (thin links)
    void smear(GF& u_smr, GF& u_naik, GF& u_thin) const {

        SmearingParameters lt = this->linkTreatment;
        auto grid = this->_grid;

        // Create a padded cell of extra padding depth=1 and fill the padding.
        PaddedCell Ghost(HaloDepth,grid);
        GF Ughost = Ghost.Exchange(u_thin);

        // This is where auxiliary N-link fields and the final smear will be stored. 
        GF Ughost_fat(Ughost.Grid());
        GF Ughost_3link(Ughost.Grid());
        GF Ughost_5linkA(Ughost.Grid());
        GF Ughost_5linkB(Ughost.Grid());

        // mu-nu plane stencil. We allow mu==nu to make indexing the stencil easier,
        // but these entries will not be used. 
        std::vector<Coordinate> shifts;
        for(int mu=0;mu<Nd;mu++)
        for(int nu=0;nu<Nd;nu++) {
            appendShift<Nd>(shifts,mu);
            appendShift<Nd>(shifts,nu);
            appendShift<Nd>(shifts,shiftSignal::NO_SHIFT);
            appendShift<Nd>(shifts,mu,Back(nu));
            appendShift<Nd>(shifts,Back(nu));
        }

        // A GeneralLocalStencil has two indices: a site and stencil index 
        GeneralLocalStencil gStencil(Ughost.Grid(),shifts);
        typedef decltype(gStencil.GetEntry(0,0)) stencilElement;

        // This is where contributions from the smearing get added together
        Ughost_fat=Zero();

        // This loop handles 3-, 5-, and 7-link constructs, minus Lepage and Naik.
        for(int mu=0;mu<Nd;mu++) {

            // TODO: This approach is slightly memory inefficient. It uses 25% extra memory 
            Ughost_3link =Zero();
            Ughost_5linkA=Zero();
            Ughost_5linkB=Zero();

            // Create the accessors
            autoView(U_v       , Ughost       , AcceleratorRead);
            autoView(U_fat_v   , Ughost_fat   , AcceleratorWrite);
            autoView(U_3link_v , Ughost_3link , AcceleratorWrite);
            autoView(U_5linkA_v, Ughost_5linkA, AcceleratorWrite);
            autoView(U_5linkB_v, Ughost_5linkB, AcceleratorWrite);

            // We infer a type that will be needed in the calculation.
            typedef decltype(getLink(U_v[0](0),gStencil.GetEntry(0,0))) U3matrix;

            int Nsites = U_v.size();
            auto gStencil_v = gStencil.View(AcceleratorRead); 

            accelerator_for(site,Nsites,Simd::Nsimd(),{ // ----------- 3-link constructs
                stencilElement SE0, SE1, SE2, SE3, SE4;
                U3matrix U0, U1, U2, U3, U4, U5, W;
                for(int nu=0;nu<Nd;nu++) {
                    if(nu==mu) continue;
                    int s = stencilIndex(mu,nu);

                    // The stencil gives us support points in the mu-nu plane that we will use to
                    // grab the links we need.
                    SE0 = gStencil_v.GetEntry(s+0,site); int x_p_mu      = SE0->_offset;
                    SE1 = gStencil_v.GetEntry(s+1,site); int x_p_nu      = SE1->_offset;
                    SE2 = gStencil_v.GetEntry(s+2,site); int x           = SE2->_offset;
                    SE3 = gStencil_v.GetEntry(s+3,site); int x_p_mu_m_nu = SE3->_offset;
                    SE4 = gStencil_v.GetEntry(s+4,site); int x_m_nu      = SE4->_offset;

                    // When you're deciding whether to take an adjoint, the question is: how is the
                    // stored link oriented compared to the one you want? If I imagine myself travelling
                    // with the to-be-updated link, I have two possible, alternative 3-link paths I can
                    // take, one starting by going to the left, the other starting by going to the right.
                    U0 = getLink(U_v[x_p_mu     ](nu),SE0);
                    U1 = getLink(U_v[x_p_nu     ](mu),SE1);
                    U2 = getLink(U_v[x          ](nu),SE2);
                    U3 = getLink(U_v[x_p_mu_m_nu](nu),SE3);
                    U4 = getLink(U_v[x_m_nu     ](mu),SE4);
                    U5 = getLink(U_v[x_m_nu     ](nu),SE4);

                    //  "left"          "right"
                    W = U2*U1*adj(U0) + adj(U5)*U4*U3;

                    // Save 3-link construct for later and add to smeared field.
                    setLink(U_3link_v[x](nu), W);

                    // The index operator (x) returns the coalesced read on GPU. The view [] index returns 
                    // a reference to the vector object. The [x](mu) returns a reference to the densely 
                    // packed (contiguous in memory) mu-th element of the vector object. 
                    setLink(U_fat_v[x](mu), U_fat_v(x)(mu) + lt.c_3*W);
                }
            })

            accelerator_for(site,Nsites,Simd::Nsimd(),{ // ----------- 5-link 
                stencilElement SE0, SE1, SE2, SE3, SE4;
                U3matrix U0, U1, U2, U3, U4, U5, W;
                int sigmaIndex = 0;
                for(int nu=0;nu<Nd;nu++) {
                    if(nu==mu) continue;
                    int s = stencilIndex(mu,nu);
                    for(int rho=0;rho<Nd;rho++) {
                        if (rho == mu || rho == nu) continue;

                        SE0 = gStencil_v.GetEntry(s+0,site); int x_p_mu      = SE0->_offset;
                        SE1 = gStencil_v.GetEntry(s+1,site); int x_p_nu      = SE1->_offset;
                        SE2 = gStencil_v.GetEntry(s+2,site); int x           = SE2->_offset;
                        SE3 = gStencil_v.GetEntry(s+3,site); int x_p_mu_m_nu = SE3->_offset;
                        SE4 = gStencil_v.GetEntry(s+4,site); int x_m_nu      = SE4->_offset;

                        U0 = getLink(      U_v[x_p_mu     ](nu ),SE0);
                        U1 = getLink(U_3link_v[x_p_nu     ](rho),SE1);
                        U2 = getLink(      U_v[x          ](nu ),SE2);
                        U3 = getLink(      U_v[x_p_mu_m_nu](nu ),SE3);
                        U4 = getLink(U_3link_v[x_m_nu     ](rho),SE4);
                        U5 = getLink(      U_v[x_m_nu     ](nu ),SE4);

                        W  = U2*U1*adj(U0) + adj(U5)*U4*U3;

                        if(sigmaIndex<3) {
                            setLink(U_5linkA_v[x](rho), W);
                        } else {
                            setLink(U_5linkB_v[x](rho), W);
                        }    

                        setLink(U_fat_v[x](mu), U_fat_v(x)(mu) + lt.c_5*W);
                        sigmaIndex++;
                    }
                }
            })

            accelerator_for(site,Nsites,Simd::Nsimd(),{ // ----------- 7-link
                stencilElement SE0, SE1, SE2, SE3, SE4;
                U3matrix U0, U1, U2, U3, U4, U5, W;
                int sigmaIndex = 0;
                for(int nu=0;nu<Nd;nu++) {
                    if(nu==mu) continue;
                    int s = stencilIndex(mu,nu);
                    for(int rho=0;rho<Nd;rho++) {
                        if (rho == mu || rho == nu) continue;

                        SE0 = gStencil_v.GetEntry(s+0,site); int x_p_mu      = SE0->_offset;
                        SE1 = gStencil_v.GetEntry(s+1,site); int x_p_nu      = SE1->_offset;
                        SE2 = gStencil_v.GetEntry(s+2,site); int x           = SE2->_offset;
                        SE3 = gStencil_v.GetEntry(s+3,site); int x_p_mu_m_nu = SE3->_offset;
                        SE4 = gStencil_v.GetEntry(s+4,site); int x_m_nu      = SE4->_offset;

                        U0 = getLink(U_v[x_p_mu](nu),SE0);
                        if(sigmaIndex<3) {
                            U1 = getLink(U_5linkB_v[x_p_nu](rho),SE1);
                        } else {
                            U1 = getLink(U_5linkA_v[x_p_nu](rho),SE1);
                        }  
                        U2 = getLink(U_v[x](nu),SE2);
                        U3 = getLink(U_v[x_p_mu_m_nu](nu),SE3);
                        if(sigmaIndex<3) {
                            U4 = getLink(U_5linkB_v[x_m_nu](rho),SE4);
                        } else {
                            U4 = getLink(U_5linkA_v[x_m_nu](rho),SE4);
                        }  
                        U5 = getLink(U_v[x_m_nu](nu),SE4);

                        W  = U2*U1*adj(U0) + adj(U5)*U4*U3;

                        setLink(U_fat_v[x](mu), U_fat_v(x)(mu) + lt.c_7*W);
                        sigmaIndex++;
                    }
                }
            })

        } // end mu loop

        // c1, c3, c5, c7 construct contributions
        u_smr = Ghost.Extract(Ughost_fat) + lt.c_1*u_thin;

        // Load up U and V std::vectors to access thin and smeared links.
        std::vector<LF> U(Nd, grid);
        std::vector<LF> V(Nd, grid);
        std::vector<LF> Vnaik(Nd, grid);
        for (int mu = 0; mu < Nd; mu++) {
            U[mu] = PeekIndex<LorentzIndex>(u_thin, mu);
            V[mu] = PeekIndex<LorentzIndex>(u_smr, mu);
        }

        for(int mu=0;mu<Nd;mu++) {

            // Naik
            Vnaik[mu] = lt.c_naik*Gimpl::CovShiftForward(U[mu],mu,
                                    Gimpl::CovShiftForward(U[mu],mu,
                                      Gimpl::CovShiftIdentityForward(U[mu],mu)));

            // LePage
            for (int nu_h=1;nu_h<Nd;nu_h++) {
                int nu=(mu+nu_h)%Nd;
                                // nu, nu, mu, Back(nu), Back(nu)
                V[mu] = V[mu] + lt.c_lp*Gimpl::CovShiftForward(U[nu],nu,
                                          Gimpl::CovShiftForward(U[nu],nu,
                                            Gimpl::CovShiftForward(U[mu],mu,
                                              Gimpl::CovShiftBackward(U[nu],nu,
                                                Gimpl::CovShiftIdentityBackward(U[nu],nu)))))
                                // Back(nu), Back(nu), mu, nu, nu
                              + lt.c_lp*Gimpl::CovShiftBackward(U[nu],nu,
                                          Gimpl::CovShiftBackward(U[nu],nu,
                                            Gimpl::CovShiftForward(U[mu],mu,
                                              Gimpl::CovShiftForward(U[nu],nu,
                                                Gimpl::CovShiftIdentityForward(U[nu],nu)))));
            }
        }

        // Put V back into u_smr.
        for (int mu = 0; mu < Nd; mu++) {
            PokeIndex<LorentzIndex>(u_smr , V[mu]    , mu);
            PokeIndex<LorentzIndex>(u_naik, Vnaik[mu], mu);
        }
    };


    // Intent: OUT--u_proj (U3-projected links)
    //          IN--u_mu (to-be-projected links)
    void projectU3(GF& u_proj, GF& u_mu) const {

        // Open up the views
        autoView(uproj_v, u_proj, AcceleratorWrite);
        autoView(umu_v  , u_mu  , AcceleratorRead);

        // Make sure everyone is using the same Grid
        conformable(u_proj,u_mu);

        // Follow MILC 10.1103/PhysRevD.82.074501, eqs (B2-B3) and (C1-C8)
        accelerator_for(ss,umu_v.size(),Simd::Nsimd(),{
#ifdef GRID_SIMT
            { int blane=acceleratorSIMTlane(Simd::Nsimd());//
#else
            for(int blane=0;blane<Simd::Nsimd();blane++) {
#endif
                RealScalar g1, g2, g0;
                ComplexColourMatrix V;
                auto Vmu = extractLane(blane,umu_v[ss]);
                for (int mu = 0; mu < Nd; mu++) {
                    V()     = Vmu(mu);
                    auto Q  = adj(V)*V;
                    RealScalar c0 =        real(trace(Q))()()();
                    RealScalar c1 = (1/2.)*real(trace(Q*Q))()()();
                    RealScalar c2 = (1/3.)*real(trace(Q*Q*Q))()()();
                    RealScalar S  = (1/3.)*c1-(1/18.)*c0*c0;
                    if (abs(S)<Scut) {
                        g0 = (1/3.)*c0; 
                        g1 = g0; 
                        g2 = g1;
                    } else {
                        RealScalar R     = (1/2.)*c2-(1/3. )*c0*c1+(1/27.)*c0*c0*c0;
                        RealScalar theta = acos(R*pow(S,-1.5));
                        g0 = (1/3.)*c0+2.*sqrt(S)*cos((1/3.)*theta-2*M_PI/3.);
                        g1 = (1/3.)*c0+2.*sqrt(S)*cos((1/3.)*theta          );
                        g2 = (1/3.)*c0+2.*sqrt(S)*cos((1/3.)*theta+2*M_PI/3.);
                    }
                    RealScalar u   = sqrt(g0) + sqrt(g1) + sqrt(g2);
                    RealScalar v   = sqrt(g0*g1) + sqrt(g0*g2) + sqrt(g1*g2);
                    RealScalar w   = sqrt(g0*g1*g2);
                    RealScalar den = w*(u*v-w);
                    RealScalar f0  = (-w*(u*u+v)+u*v*v)/den;
                    RealScalar f1  = (-w-u*u*u+2.*u*v)/den;
                    RealScalar f2  = u/den;

                    auto res = V*(f0 + f1*Q + f2*Q*Q);

                    insertLane(blane,uproj_v[ss](mu),res());
                }
            } 
        });
    };

};



/*!  @brief compute force from link variables */
template<class Gimpl> 
class Force_HISQ : public Gimpl {
public:

    GridCartesian* const _grid;

    // Sort out the Gimpl. This handles BCs and part of the precision. 
    INHERIT_GIMPL_TYPES(Gimpl);
    typedef typename Gimpl::GaugeField     GF;
    typedef typename Gimpl::GaugeLinkField LF;
    typedef typename Gimpl::ComplexField   CF;
    typedef typename Gimpl::Scalar ComplexScalar;
    typedef decltype(real(ComplexScalar())) RealScalar;
    typedef iColourMatrix<ComplexScalar> ComplexColourMatrix;

    RealScalar Scut=-1; // Cutoff for U(3) projection eigenvalues, set at initialization
    int HaloDepth=1; 

    SmearingParameters<RealScalar> linkTreatment;

    void initialize() {
        if (sizeof(RealScalar)==4) {
            Scut=1e-5; // Maybe should be higher? e.g. 1e-4
        } else if (sizeof(RealScalar)==8) { 
            Scut=1e-8;
        } else {
            Grid_error("HISQ force only implemented for single and double");
        }
        assert(Nc == 3 && "HISQ force currently implemented only for Nc==3");
        assert(Nd == 4 && "HISQ force only defined for Nd==4");
    }

    Force_HISQ(GridCartesian* grid, RealScalar c1, RealScalar cnaik, RealScalar c3, RealScalar c5, RealScalar c7, RealScalar clp) 
        : _grid(grid), 
          linkTreatment(c1,cnaik,c3,c5,c7,clp) {
        initialize();
    }

    // Allow to pass a pointer to a C-style array for MILC convenience
    Force_HISQ(GridCartesian* grid, double* coeff) 
        : _grid(grid), 
          linkTreatment(coeff[0],coeff[1],coeff[2],coeff[3],coeff[4],coeff[5]) {
        initialize();
    }
    Force_HISQ(GridCartesian* grid, float* coeff) 
        : _grid(grid), 
          linkTreatment(coeff[0],coeff[1],coeff[2],coeff[3],coeff[4],coeff[5]) {
        initialize();
    }

    ~Force_HISQ() {}


    // Intent: OUT--u_deriv (dW/dV slotted into force)
    //          IN--u_mu (fat links), 
    //              u_force (slot derivative into this force), 
    //              delta (force cutoff)
    // Follow MILC 10.1103/PhysRevD.82.074501
    void ddVprojectU3(GF& u_deriv, GF& u_mu, GF& u_force, RealScalar const delta=5e-5) {

        conformable(u_force,u_mu);
        conformable(u_deriv,u_mu);

        autoView(uderiv_v, u_deriv, AcceleratorWrite);
        autoView(umu_v   , u_mu   , AcceleratorRead);
        autoView(uforce_v, u_force, AcceleratorRead);

        // Follow MILC 10.1103/PhysRevD.82.074501, eqs (B2-B3) and (C1-C8)
        accelerator_for(ss,umu_v.size(),Simd::Nsimd(),{
#ifdef GRID_SIMT
            { int blane=acceleratorSIMTlane(Simd::Nsimd());//
#else
            for(int blane=0;blane<Simd::Nsimd();blane++) {
#endif
                RealScalar g1, g2, g0;
                ComplexColourMatrix V, force;
                auto Vmu     = extractLane(blane,umu_v[ss]);
                auto forcemu = extractLane(blane,uforce_v[ss]);
                for (int mu = 0; mu < Nd; mu++) {
                    V()     = Vmu(mu);
                    auto Q  = adj(V)*V;
                    RealScalar c0 =        real(trace(Q))()()();
                    RealScalar c1 = (1/2.)*real(trace(Q*Q))()()();
                    RealScalar c2 = (1/3.)*real(trace(Q*Q*Q))()()();
                    RealScalar S  = (1/3.)*c1-(1/18.)*c0*c0;
                    if (abs(S)<Scut) {
                        g0 = (1/3.)*c0; 
                        g1 = g0; 
                        g2 = g1;
                    } else {
                        RealScalar R     = (1/2.)*c2-(1/3. )*c0*c1+(1/27.)*c0*c0*c0;
                        RealScalar theta = acos(R*pow(S,-1.5));
                        g0 = (1/3.)*c0+2.*sqrt(S)*cos((1/3.)*theta-2*M_PI/3.);
                        g1 = (1/3.)*c0+2.*sqrt(S)*cos((1/3.)*theta          );
                        g2 = (1/3.)*c0+2.*sqrt(S)*cos((1/3.)*theta+2*M_PI/3.);
                    }

                    if (g0 < delta || g1 < delta || g2 < delta) {
                        // force filter eq (C23)
                        g0 += delta;
                        g1 += delta;
                        g2 += delta;
                        Q   = Q + delta;
                    }
//                if (fabs(Q.determinant()/(g0*g1*g2)-1.0) > 1e-5) { SVD }

                    RealScalar u   = sqrt(g0) + sqrt(g1) + sqrt(g2);
                    RealScalar v   = sqrt(g0*g1) + sqrt(g0*g2) + sqrt(g1*g2);
                    RealScalar w   = sqrt(g0*g1*g2);
                    RealScalar den = w*(u*v-w);
                    RealScalar f0  = (-w*(u*u+v)+u*v*v)/den;
                    RealScalar f1  = (-w-u*u*u+2.*u*v)/den;
                    RealScalar f2  = u/den;

                    auto Qinvsq = f0 + f1*Q + f2*Q*Q;

                    force() = forcemu(mu);
                    auto forcedag = adj(force);

                    RealScalar u2 = u  * u;
                    RealScalar u3 = u2 * u;
                    RealScalar u4 = u3 * u;
                    RealScalar u5 = u4 * u;
                    RealScalar u6 = u5 * u;
                    RealScalar u7 = u6 * u;
                    RealScalar u8 = u7 * u;
                    RealScalar v2 = v  * v;
                    RealScalar v3 = v2 * v;
                    RealScalar v4 = v3 * v;
                    RealScalar v5 = v4 * v;
                    RealScalar v6 = v5 * v;
                    RealScalar w2 = w  * w;
                    RealScalar w3 = w2 * w;
                    RealScalar w4 = w3 * w;
                    RealScalar w5 = w4 * w;
        
                    // eq (C10)
                    auto d = 2*w3*(u*v-w)*(u*v-w)*(u*v-w);
        
                    // eq (C11)
                    auto C00  = ( -w3*u6 + 3*v*w3*u4 + 3*v4*w*u4 - v6*u3 - 4*w4*u3 - 12*v3*w2*u3 + 16*v2*w3*u2 
                                  + 3*v5*w*u2 - 8*v*w4*u - 3*v4*w2*u + w5 + v3*w3 )/d;
                    auto C01  = ( -w2*u7 - v2*w*u6 + v4*u5 + 6*v*w2*u5 - 5*w3*u4 - v3*w*u4 - 2*v5*u3 - 6*v2*w2*u3 
                                  + 10*v*w3*u2 + 6*v4*w*u2 - 3*w4*u - 6*v3*w2*u + 2*v2*w3 )/d;
                    auto C02  = ( w2*u5 + v2*w*u4 - v4*u3 - 4*v*w2*u3 + 4*w3*u2 +3*v3*w*u2 - 3*v2*w2*u + v*w3 )/d;
                    auto C11  = ( -w*u8 - v2*u7 + 7*v*w*u6 + 4*v3*u5 - 5*w2*u5 - 16*v2*w*u4 - 4*v4*u3 + 16*v*w2*u3 
                                  - 3*w3*u2 + 12*v3*w*u2 - 12*v2*w2*u + 3*v*w3 )/d;
                    auto C12  = ( w*u6 + v2*u5 - 5*v*w*u4 - 2*v3*u3 + 4*w2*u3 + 6*v2*w*u2 - 6*v*w2*u + w3 )/d;
                    auto C22  = ( -w*u4 - v2*u3 + 3*v*w*u2 - 3*w2*u )/d;
        
                    // These are all used in the loop over color entries, and we want to avoid recomputing
                    // these products, which should be broadcast to all sites, 3*3*3*3=81 times. 
                    auto Vdag   = adj(V);
                    auto VVdag  = V*Vdag;
                    auto VQ     = V*Q;   
                    auto VQ2    = VQ*Q;
                    auto VQVdag = VQ*Vdag;
                    auto QVdag  = Q*Vdag;
                    auto Q2Vdag = Q*QVdag;
        
                    // eqs (C17-C19)
                    auto PVdag  = ( C00 + C01*Q + C02*Q*Q )*Vdag;
                    auto RVdag  = ( C01 + C11*Q + C12*Q*Q )*Vdag;
                    auto SVdag  = ( C02 + C12*Q + C22*Q*Q )*Vdag;
        
                    // eqs (C20) and (C21)
                    ComplexColourMatrix res = Zero();
                    for (int k = 0; k < 3; k++)
                    for (int l = 0; l < 3; l++)
                    for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++) {
        
                        ComplexScalar deriv = 0.; // dWij/dVkl
                
                        if (k == i) deriv += Qinvsq()()(l,j);
                        if (l == j) deriv += f1*VVdag()()(i,k)+f2*VQVdag()()(i,k);

                        deriv += f2*VVdag()()(i,k)*Q()()(l,j) + V()()(i,j)*PVdag()()(l,k) 
                                 + VQ()()(i,j)*RVdag()()(l,k) + VQ2()()(i,j)*SVdag()()(l,k);
        
                        res()()(l,k) = res()()(l,k) + deriv*force()()(j,i);
                
                        // dWij^+/dVkl
                        deriv = (f1*Vdag()()(i,k)+f2*QVdag()()(i,k))*Vdag()()(l,j) 
                                + f2*Vdag()()(i,k)*QVdag()()(l,j) + Vdag()()(i,j)*PVdag()()(l,k) 
                                + QVdag()()(i,j)*RVdag()()(l,k)+Q2Vdag()()(i,j)*SVdag()()(l,k);
                
                        res()()(l,k) = res()()(l,k) + deriv*forcedag()()(j,i);
                  	}
        
                    insertLane(blane,uderiv_v[ss](mu),res());
                }
            } 
        });
    };


    void ddV_3link(GF& u_deriv, GF& u_mu, GF& u_force) {

        SmearingParameters lt = this->linkTreatment;
        auto grid = this->_grid;

        PaddedCell Ghost(HaloDepth,grid);
        GF Ughost = Ghost.Exchange(u_mu);
        GF Fghost = Ghost.Exchange(u_force);

        GF Ughost_deriv(Ughost.Grid());

        Ughost_deriv = Zero();

        std::vector<Coordinate> shifts;
        for(int mu=0;mu<Nd;mu++)
        for(int nu=0;nu<Nd;nu++) {
            appendShift<Nd>(shifts,mu);
            appendShift<Nd>(shifts,nu);
            appendShift<Nd>(shifts,shiftSignal::NO_SHIFT);
            appendShift<Nd>(shifts,mu,Back(nu));
            appendShift<Nd>(shifts,Back(nu));
            appendShift<Nd>(shifts,Back(mu));
        }

        GeneralLocalStencil gStencil(Ughost.Grid(),shifts);
        typedef decltype(gStencil.GetEntry(0,0)) stencilElement;

        for(int mu=0;mu<Nd;mu++) {

            autoView(U_v      , Ughost      , AcceleratorRead);
            autoView(F_v      , Fghost      , AcceleratorRead);
            autoView(U_deriv_v, Ughost_deriv, AcceleratorWrite);

            typedef decltype(getLink(U_v[0](0),gStencil.GetEntry(0,0))) U3matrix;

            int Nsites = U_v.size();
            auto gStencil_v = gStencil.View(AcceleratorRead);

            accelerator_for(site,Nsites,Simd::Nsimd(),{ 
                stencilElement SE0, SE1, SE2, SE3, SE4;
                U3matrix U0, U1, U2, U3, U4, U5, F0, F1, F2, F3, F4, F5, W;
                for(int nu=0;nu<Nd;nu++) {
                    if(nu==mu) continue;
                    int s = stencilIndex(mu,nu);

                    SE0 = gStencil_v.GetEntry(s+0,site); int x_p_mu      = SE0->_offset;
                    SE1 = gStencil_v.GetEntry(s+1,site); int x_p_nu      = SE1->_offset;
                    SE2 = gStencil_v.GetEntry(s+2,site); int x           = SE2->_offset;
                    SE3 = gStencil_v.GetEntry(s+3,site); int x_p_mu_m_nu = SE3->_offset;
                    SE4 = gStencil_v.GetEntry(s+4,site); int x_m_nu      = SE4->_offset;

                    U0 = getLink(U_v[x_p_mu     ](nu),SE0);
                    U1 = getLink(U_v[x_p_nu     ](mu),SE1);
                    U2 = getLink(U_v[x          ](nu),SE2);
                    U3 = getLink(U_v[x_p_mu_m_nu](nu),SE3);
                    U4 = getLink(U_v[x_m_nu     ](mu),SE4);
                    U5 = getLink(U_v[x_m_nu     ](nu),SE4);

                    F0 = getLink(F_v[x_p_mu     ](nu),SE0);
                    F1 = getLink(F_v[x_p_nu     ](mu),SE1);
                    F2 = getLink(F_v[x          ](nu),SE2);
                    F3 = getLink(F_v[x_p_mu_m_nu](nu),SE3);
                    F4 = getLink(F_v[x_m_nu     ](mu),SE4);
                    F5 = getLink(F_v[x_m_nu     ](nu),SE4);

                    W  =   adj(F2)*U1*adj(U0) +     U2 *adj(F1)*adj(U0) +     U2 *U1*    F0 
                         +     F5 *U4*    U3  + adj(U5)*adj(F4)*    U3  + adj(U5)*U4*adj(F3);                    

                    setLink(U_deriv_v[x](mu), U_deriv_v(x)(mu) + lt.c_3*W);
                }              
            })
        } // end mu loop

        u_deriv = Ghost.Extract(Ughost_deriv);
    }

    void ddV_naik(GF& u_deriv, GF& u_mu, GF& u_force) {

        SmearingParameters lt = this->linkTreatment;
        auto grid = this->_grid;

        PaddedCell Ghost(3,grid);
        GF Ughost = Ghost.Exchange(u_mu);
        GF Fghost = Ghost.Exchange(u_force);

        GF Ughost_deriv(Ughost.Grid());

        Ughost_deriv = Zero();

        std::vector<Coordinate> shifts;
        for(int mu=0;mu<Nd;mu++) {
            appendShift<Nd>(shifts, shiftSignal::NO_SHIFT);
            appendShift<Nd>(shifts, mu);
            appendShift<Nd>(shifts, mu, mu);
            appendShift<Nd>(shifts, Back(mu));
            appendShift<Nd>(shifts, Back(mu), Back(mu));
            appendShift<Nd>(shifts, Back(mu), Back(mu), Back(mu));
        }

        GeneralLocalStencil gStencil(Ughost.Grid(),shifts);
        typedef decltype(gStencil.GetEntry(0,0)) stencilElement;

        autoView(U_v      , Ughost      , AcceleratorRead);
        autoView(F_v      , Fghost      , AcceleratorRead);
        autoView(U_deriv_v, Ughost_deriv, AcceleratorWrite);

        typedef decltype(getLink(U_v[0](0),gStencil.GetEntry(0,0))) U3matrix;

        int Nsites = U_v.size();
        auto gStencil_v = gStencil.View(AcceleratorRead);

        accelerator_for(site,Nsites,Simd::Nsimd(),{ 
            stencilElement SE0, SE1, SE2, SE3, SE4, SE5;
            U3matrix U0, U1, U2, U3, U4, U5, F0, F1, F2, F3, F4, F5, V;
            int s = 0;
            for(int mu=0;mu<Nd;mu++) {

                SE0 = gStencil_v.GetEntry(s+0,site); int x       = SE0->_offset;
                SE1 = gStencil_v.GetEntry(s+1,site); int x_p_mu  = SE1->_offset;
                SE2 = gStencil_v.GetEntry(s+2,site); int x_p_2mu = SE2->_offset;
                SE3 = gStencil_v.GetEntry(s+3,site); int x_m_mu  = SE3->_offset;
                SE4 = gStencil_v.GetEntry(s+4,site); int x_m_2mu = SE4->_offset;
                SE5 = gStencil_v.GetEntry(s+5,site); int x_m_3mu = SE5->_offset;

                U0 = getLink(U_v[x      ](mu),SE0);
                U1 = getLink(U_v[x_p_mu ](mu),SE1);
                U2 = getLink(U_v[x_p_2mu](mu),SE2);
                U3 = getLink(U_v[x_m_mu ](mu),SE3);
                U4 = getLink(U_v[x_m_2mu](mu),SE4);
                U5 = getLink(U_v[x_m_3mu](mu),SE5);

                F0 = getLink(F_v[x      ](mu),SE0);
                F1 = getLink(F_v[x_p_mu ](mu),SE1);
                F2 = getLink(F_v[x_p_2mu](mu),SE2);
                F3 = getLink(F_v[x_m_mu ](mu),SE3);
                F4 = getLink(F_v[x_m_2mu](mu),SE4);
                F5 = getLink(F_v[x_m_3mu](mu),SE5);

                //     ********Forward********   *******Backward********
                V  =  (adj(F2)*    U1 *    U0 )+(adj(U5)*adj(U4)*    F3 ) 
                     +(    U2 *adj(F1)*    U0 )+(adj(U5)*    F4 *adj(U3))
                     +(    U2 *    U1 *adj(F0))+(    F5 *adj(U4)*adj(U3));

                setLink(U_deriv_v[x](mu), U_deriv_v(x)(mu) + lt.c_naik*V);

                s += 6;
            }              
        });
        u_deriv = Ghost.Extract(Ughost_deriv);
    }

};


NAMESPACE_END(Grid);
