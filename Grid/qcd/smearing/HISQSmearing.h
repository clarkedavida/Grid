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
struct SmearingParameters{
    SmearingParameters(){}
    Real c_1;               // 1 link
    Real c_naik;            // Naik term
    Real c_3;               // 3 link
    Real c_5;               // 5 link
    Real c_7;               // 7 link
    Real c_lp;              // 5 link Lepage
    SmearingParameters(Real c1, Real cnaik, Real c3, Real c5, Real c7, Real clp) 
        : c_1(c1),
          c_naik(cnaik),
          c_3(c3),
          c_5(c5),
          c_7(c7),
          c_lp(clp){}
};


/*!  @brief create fat links from link variables */
template<class Gimpl> 
class Smear_HISQ : public Gimpl {

private:
    GridCartesian* const _grid;
    Real const _Scut = 1e-16; // Cutoff for U(3) projection eigenvalues
    SmearingParameters _linkTreatment;

    // figure out the stencil index from mu and nu
    accelerator_inline int stencilIndex(int mu, int nu) const {
        // Nshifts depends on how you built the stencil
        int Nshifts = 6;
        return Nshifts*nu + Nd*Nshifts*mu;
    }

public:

    INHERIT_GIMPL_TYPES(Gimpl);
    typedef typename Gimpl::GaugeField     GF;
    typedef typename Gimpl::GaugeLinkField LF;
    typedef typename Gimpl::ComplexField   CF;

    Smear_HISQ(GridCartesian* grid, Real c1, Real cnaik, Real c3, Real c5, Real c7, Real clp) 
        : _grid(grid), 
          _linkTreatment(c1,cnaik,c3,c5,c7,clp) {
        assert(Nc == 3 && "HISQ smearing currently implemented only for Nc==3");
        assert(Nd == 4 && "HISQ smearing only defined for Nd==4");
    }

    // Allow to pass a pointer to a C-style, double array for MILC convenience
    Smear_HISQ(GridCartesian* grid, double* coeff) 
        : _grid(grid), 
          _linkTreatment(coeff[0],coeff[1],coeff[2],coeff[3],coeff[4],coeff[5]) {
        assert(Nc == 3 && "HISQ smearing currently implemented only for Nc==3");
        assert(Nd == 4 && "HISQ smearing only defined for Nd==4");
    }

    ~Smear_HISQ() {}

    // Intent: OUT--u_smr (smeared links), 
    //              u_naik (Naik links),
    //          IN--u_thin (think links)
    void smear(GF& u_smr, GF& u_naik, GF& u_thin) const {

        SmearingParameters lt = this->_linkTreatment;
        auto grid = this->_grid;

        // Create a padded cell of extra padding depth=1 and fill the padding.
        int depth = 1;
        PaddedCell Ghost(depth,grid);
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
            appendShift<Nd>(shifts,Back(mu));
        }

        // A GeneralLocalStencil has two indices: a site and stencil index 
        GeneralLocalStencil gStencil(Ughost.Grid(),shifts);

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

            // We infer some types that will be needed in the calculation.
            typedef decltype(gStencil.GetEntry(0,0)) stencilElement;
            typedef decltype(coalescedReadGeneralPermute(U_v[0](0),gStencil.GetEntry(0,0)->_permute,Nd)) U3matrix;

            int Nsites = U_v.size();
            auto gStencil_v = gStencil.View(); 

            accelerator_for(site,Nsites,Simd::Nsimd(),{ // ----------- 3-link constructs
                stencilElement SE0, SE1, SE2, SE3, SE4, SE5;
                U3matrix U0, U1, U2, U3, U4, U5, W;
                for(int nu=0;nu<Nd;nu++) {
                    if(nu==mu) continue;
                    int s = this->stencilIndex(mu,nu);

                    // The stencil gives us support points in the mu-nu plane that we will use to
                    // grab the links we need.
                    SE0 = gStencil_v.GetEntry(s+0,site); int x_p_mu      = SE0->_offset;
                    SE1 = gStencil_v.GetEntry(s+1,site); int x_p_nu      = SE1->_offset;
                    SE2 = gStencil_v.GetEntry(s+2,site); int x           = SE2->_offset;
                    SE3 = gStencil_v.GetEntry(s+3,site); int x_p_mu_m_nu = SE3->_offset;
                    SE4 = gStencil_v.GetEntry(s+4,site); int x_m_nu      = SE4->_offset;
                    SE5 = gStencil_v.GetEntry(s+5,site); int x_m_mu      = SE5->_offset;

                    // When you're deciding whether to take an adjoint, the question is: how is the
                    // stored link oriented compared to the one you want? If I imagine myself travelling
                    // with the to-be-updated link, I have two possible, alternative 3-link paths I can
                    // take, one starting by going to the left, the other starting by going to the right.
                    U0 = coalescedReadGeneralPermute(U_v[x_p_mu     ](nu),SE0->_permute,Nd);
                    U1 = coalescedReadGeneralPermute(U_v[x_p_nu     ](mu),SE1->_permute,Nd);
                    U2 = coalescedReadGeneralPermute(U_v[x          ](nu),SE2->_permute,Nd);
                    U3 = coalescedReadGeneralPermute(U_v[x_p_mu_m_nu](nu),SE3->_permute,Nd);
                    U4 = coalescedReadGeneralPermute(U_v[x_m_nu     ](mu),SE4->_permute,Nd);
                    U5 = coalescedReadGeneralPermute(U_v[x_m_nu     ](nu),SE4->_permute,Nd);

                    //  "left"          "right"
                    W = U2*U1*adj(U0) + adj(U5)*U4*U3;

                    // Save 3-link construct for later and add to smeared field.
                    coalescedWrite(U_3link_v[x](nu), W);

                    // The index operator (x) returns the coalesced read on GPU. The view [] index returns 
                    // a reference to the vector object. The [x](mu) returns a reference to the densely 
                    // packed (contiguous in memory) mu-th element of the vector object. On CPU, 
                    // coalescedRead/Write is the identity mapping assigning vector object to vector object.
                    // But on GPU it's non-trivial and maps scalar object to vector object and vice versa.
                    coalescedWrite(U_fat_v[x](mu), U_fat_v(x)(mu) + lt.c_3*W);
                }
            })

            accelerator_for(site,Nsites,Simd::Nsimd(),{ // ----------- 5-link 
                stencilElement SE0, SE1, SE2, SE3, SE4, SE5;
                U3matrix U0, U1, U2, U3, U4, U5, W;
                int sigmaIndex = 0;
                for(int nu=0;nu<Nd;nu++) {
                    if(nu==mu) continue;
                    int s = this->stencilIndex(mu,nu);
                    for(int rho=0;rho<Nd;rho++) {
                        if (rho == mu || rho == nu) continue;

                        SE0 = gStencil_v.GetEntry(s+0,site); int x_p_mu      = SE0->_offset;
                        SE1 = gStencil_v.GetEntry(s+1,site); int x_p_nu      = SE1->_offset;
                        SE2 = gStencil_v.GetEntry(s+2,site); int x           = SE2->_offset;
                        SE3 = gStencil_v.GetEntry(s+3,site); int x_p_mu_m_nu = SE3->_offset;
                        SE4 = gStencil_v.GetEntry(s+4,site); int x_m_nu      = SE4->_offset;

                        U0 = coalescedReadGeneralPermute(      U_v[x_p_mu     ](nu ),SE0->_permute,Nd);
                        U1 = coalescedReadGeneralPermute(U_3link_v[x_p_nu     ](rho),SE1->_permute,Nd);
                        U2 = coalescedReadGeneralPermute(      U_v[x          ](nu ),SE2->_permute,Nd);
                        U3 = coalescedReadGeneralPermute(      U_v[x_p_mu_m_nu](nu ),SE3->_permute,Nd);
                        U4 = coalescedReadGeneralPermute(U_3link_v[x_m_nu     ](rho),SE4->_permute,Nd);
                        U5 = coalescedReadGeneralPermute(      U_v[x_m_nu     ](nu ),SE4->_permute,Nd);

                        W  = U2*U1*adj(U0) + adj(U5)*U4*U3;

                        if(sigmaIndex<3) {
                            coalescedWrite(U_5linkA_v[x](rho), W);
                        } else {
                            coalescedWrite(U_5linkB_v[x](rho), W);
                        }    

                        coalescedWrite(U_fat_v[x](mu), U_fat_v(x)(mu) + lt.c_5*W);
                        sigmaIndex++;
                    }
                }
            })

            accelerator_for(site,Nsites,Simd::Nsimd(),{ // ----------- 7-link
                stencilElement SE0, SE1, SE2, SE3, SE4, SE5;
                U3matrix U0, U1, U2, U3, U4, U5, W;
                int sigmaIndex = 0;
                for(int nu=0;nu<Nd;nu++) {
                    if(nu==mu) continue;
                    int s = this->stencilIndex(mu,nu);
                    for(int rho=0;rho<Nd;rho++) {
                        if (rho == mu || rho == nu) continue;

                        SE0 = gStencil_v.GetEntry(s+0,site); int x_p_mu      = SE0->_offset;
                        SE1 = gStencil_v.GetEntry(s+1,site); int x_p_nu      = SE1->_offset;
                        SE2 = gStencil_v.GetEntry(s+2,site); int x           = SE2->_offset;
                        SE3 = gStencil_v.GetEntry(s+3,site); int x_p_mu_m_nu = SE3->_offset;
                        SE4 = gStencil_v.GetEntry(s+4,site); int x_m_nu      = SE4->_offset;

                        U0 = coalescedReadGeneralPermute(U_v[x_p_mu](nu),SE0->_permute,Nd);
                        if(sigmaIndex<3) {
                            U1 = coalescedReadGeneralPermute(U_5linkB_v[x_p_nu](rho),SE1->_permute,Nd);
                        } else {
                            U1 = coalescedReadGeneralPermute(U_5linkA_v[x_p_nu](rho),SE1->_permute,Nd);
                        }  
                        U2 = coalescedReadGeneralPermute(U_v[x](nu),SE2->_permute,Nd);
                        U3 = coalescedReadGeneralPermute(U_v[x_p_mu_m_nu](nu),SE3->_permute,Nd);
                        if(sigmaIndex<3) {
                            U4 = coalescedReadGeneralPermute(U_5linkB_v[x_m_nu](rho),SE4->_permute,Nd);
                        } else {
                            U4 = coalescedReadGeneralPermute(U_5linkA_v[x_m_nu](rho),SE4->_permute,Nd);
                        }  
                        U5 = coalescedReadGeneralPermute(U_v[x_m_nu](nu),SE4->_permute,Nd);

                        W  = U2*U1*adj(U0) + adj(U5)*U4*U3;

                        coalescedWrite(U_fat_v[x](mu), U_fat_v(x)(mu) + lt.c_7*W);
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

        auto grid = this->_grid;

        LF V(grid), Q(grid), sqrtQinv(grid), id_3(grid);
        CF c0(grid), c1(grid), c2(grid), g0(grid), g1(grid), g2(grid), S(grid), R(grid), theta(grid), 
           u(grid), v(grid), w(grid), den(grid), f0(grid), f1(grid), f2(grid);

        id_3  = 1.;

        // Follow MILC 10.1103/PhysRevD.82.074501, eqs (B2-B3) and (C1-C8)
        for (int mu = 0; mu < Nd; mu++) {
            V  = PeekIndex<LorentzIndex>(u_mu, mu);
            Q  = adj(V)*V;
            c0 =        real(trace(Q));
            c1 = (1/2.)*real(trace(Q*Q));
            c2 = (1/3.)*real(trace(Q*Q*Q));
            S  = (1/3.)*c1-(1/18.)*c0*c0;
            if (norm2(S)<this->_Scut) {
                g0 = (1/3.)*c0; g1 = g0; g2 = g1;
            } else {
                R     = (1/2.)*c2-(1/3. )*c0*c1+(1/27.)*c0*c0*c0;
                theta = acos(R*pow(S,-1.5));
                g0    = (1/3.)*c0+2.*sqrt(S)*cos((1/3.)*theta-2*M_PI/3.);
                g1    = (1/3.)*c0+2.*sqrt(S)*cos((1/3.)*theta          );
                g2    = (1/3.)*c0+2.*sqrt(S)*cos((1/3.)*theta+2*M_PI/3.);
            }
//            if (fabs(Q.determinant()/(g0*g1*g2)-1.0) > 1e-5) { SVD }
            u     = sqrt(g0)    + sqrt(g1)    + sqrt(g2);
            v     = sqrt(g0*g1) + sqrt(g0*g2) + sqrt(g1*g2);
            w     = sqrt(g0*g1*g2);
            den   = w*(u*v-w);
            f0    = (-w*(u*u+v)+u*v*v)/den;
            f1    = (-w-u*u*u+2.*u*v)/den;
            f2    = u/den;

            sqrtQinv = f0*id_3 + f1*Q + f2*Q*Q;

            PokeIndex<LorentzIndex>(u_proj, V*sqrtQinv, mu);
        }
    };

    // Intent: OUT--u_deriv (dW/dV slotted into force)
    //          IN--u_mu (fat links), 
    //              u_force (slot derivative into this force), 
    //              delta (force cutoff)
    // Follow MILC 10.1103/PhysRevD.82.074501
    void ddVprojectU3(GF& u_deriv, GF& u_mu, GF& u_force, Real const delta=5e-5) {

        auto grid = this->_grid;

        LF V(grid), Q(grid), sqrtQinv(grid), id_3(grid), res(grid), force(grid), forcedag(grid), 
           Vdag(grid), VVdag(grid), VQVdag(grid), PVdag(grid), VQ(grid), RVdag(grid), VQ2(grid), 
           SVdag(grid), QVdag(grid), Q2Vdag(grid);

        CF c0(grid), c1(grid), c2(grid), g0(grid), g1(grid), g2(grid), S(grid), R(grid), theta(grid), 
           u(grid), v(grid), w(grid), den(grid), f0(grid), f1(grid), f2(grid), delta(grid),
           u2(grid), u3(grid), u4(grid), u5(grid), u6(grid), u7(grid), u8(grid), 
           v2(grid), v3(grid), v4(grid), v5(grid), v6(grid), 
           w2(grid), w3(grid), w4(grid), w5(grid), d(grid),
           C00(grid), C01(grid), C02(grid), C11(grid), C12(grid), C22(grid), deriv(grid);

        id_3  = 1.;

        // eqs (B2-B3) and (C1-C8)
        for (int mu = 0; mu < Nd; mu++) {
            V  = PeekIndex<LorentzIndex>(u_mu, mu);
            Q  = adj(V)*V;
            c0 =        real(trace(Q));
            c1 = (1/2.)*real(trace(Q*Q));
            c2 = (1/3.)*real(trace(Q*Q*Q));
            S  = (1/3.)*c1-(1/18.)*c0*c0;
            if (norm2(S)<this->_Scut) {
                g0 = (1/3.)*c0; g1 = g0; g2 = g1;
            } else {
                R     = (1/2.)*c2-(1/3. )*c0*c1+(1/27.)*c0*c0*c0;
                theta = acos(R*pow(S,-1.5));
                g0    = (1/3.)*c0+2.*sqrt(S)*cos((1/3.)*theta-2*M_PI/3.);
                g1    = (1/3.)*c0+2.*sqrt(S)*cos((1/3.)*theta          );
                g2    = (1/3.)*c0+2.*sqrt(S)*cos((1/3.)*theta+2*M_PI/3.);
            }
//            if (g0 < delta || g1 < delta || g2 < delta) {
//                // force filter eq (C23)
//                g0 += g0 + delta;
//                g1 += g1 + delta;
//                g2 += g2 + delta;
//                Q  += delta*id_3;
//            }
//            if (fabs(Q.determinant()/(g0*g1*g2)-1.0) > 1e-5) { SVD }
            u     = sqrt(g0)    + sqrt(g1)    + sqrt(g2);
            v     = sqrt(g0*g1) + sqrt(g0*g2) + sqrt(g1*g2);
            w     = sqrt(g0*g1*g2);
            den   = w*(u*v-w);
            f0    = (-w*(u*u+v)+u*v*v)/den;
            f1    = (-w-u*u*u+2.*u*v)/den;
            f2    = u/den;

            sqrtQinv = f0*id_3 + f1*Q + f2*Q*Q;

            force    = PeekIndex<LorentzIndex>(u_force, mu);
            forcedag = adj(force);

            // Ask Peter: is this necessary/helpful?
            u2 = u  * u;
            u3 = u2 * u;
            u4 = u3 * u;
            u5 = u4 * u;
            u6 = u5 * u;
            u7 = u6 * u;
            u8 = u7 * u;
            v2 = v  * v;
            v3 = v2 * v;
            v4 = v3 * v;
            v5 = v4 * v;
            v6 = v5 * v;
            w2 = w  * w;
            w3 = w2 * w;
            w4 = w3 * w;
            w5 = w4 * w;

            // eq (C10)
            d = 2*w3*(u*v-w)*(u*v-w)*(u*v-w);

            // eq (C11)
            C00  = ( -w3*u6 + 3*v*w3*u4 + 3*v4*w*u4 - v6*u3 - 4*w4*u3 - 12*v3*w2*u3 + 16*v2*w3*u2 
                     + 3*v5*w*u2 - 8*v*w4*u - 3*v4*w2*u + w5 + v3*w3 )/d;
            C01  = ( -w2*u7 - v2*w*u6 + v4*u5 + 6*v*w2*u5 - 5*w3*u4 - v3*w*u4 - 2*v5*u3 - 6*v2*w2*u3 
                     + 10*v*w3*u2 + 6*v4*w*u2 - 3*w4*u - 6*v3*w2*u + 2*v2*w3 )/d;
            C02  = ( w2*u5 + v2*w*u4 - v4*u3 - 4*v*w2*u3 + 4*w3*u2 +3*v3*w*u2 - 3*v2*w2*u + v*w3 )/d;
            C11  = ( -w*u8 - v2*u7 + 7*v*w*u6 + 4*v3*u5 - 5*w2*u5 - 16*v2*w*u4 - 4*v4*u3 + 16*v*w2*u3 
                      - 3*w3*u2 + 12*v3*w*u2 - 12*v2*w2*u + 3*v*w3 )/d;
            C12  = ( w*u6 + v2*u5 - 5*v*w*u4 - 2*v3*u3 + 4*w2*u3 + 6*v2*w*u2 - 6*v*w2*u + w3 )/d;
            C22  = ( -w*u4 - v2*u3 + 3*v*w*u2 - 3*w2*u )/d;

            Vd   = adj(V);
            VVd  = V*Vd;
            VQVd = V*Q*Vd;
            VQ   = V*Q;
            VQ2  = V*Q*Q;
            QVd  = Q*Vd;
            Q2Vd = Q*Q*Vd;

            // eqs (C17-C19)
            PVd  = ( C00*id_3 + C01*Q + C02*Q*Q )*Vd;
            RVd  = ( C01*id_3 + C11*Q + C12*Q*Q )*Vd;
            SVd  = ( C02*id_3 + C12*Q + C22*Q*Q )*Vd;

            // eqs (C20) and (C21)
            for (int k = 0; k < 3; k++)
            for (int l = 0; l < 3; l++)
            for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {

                deriv = Zero(); // dWij/dVkl
        
                if (k == i) deriv += Qinvsq()()(l,j);
                if (l == j) deriv += f1*VVd()()(i,k)+f2*VQVd()()(i,k);
        
                deriv += f2*VVd()()(i,k)*Q()()(l,j) + V()()(i,j)*PVd()()(l,k) 
                         + VQ()()(i,j)*RVd()()(l,k) + VQ2()()(i,j)*SVd()()(l,k);

                res()()(l,k) = res()()(l,k) + deriv*force()()(j,i);
        
                // dWij^+/dVkl
                deriv = (f1*Vd()()(i,k)+f2*QVd()()(i,k))*Vd()()(l,j) 
                        + f2*Vd()()(i,k)*QVd()()(l,j) + Vd()()(i,j)*PVd()()(l,k) 
                        + QVd()()(i,j)*RVd()()(l,k)+Q2Vd()()(i,j)*SVd()()(l,k);
        
                res()()(l,k) = res()()(l,k) + deriv*forcedag()()(j,i);
          	}

            PokeIndex<LorentzIndex>(u_deriv, res, mu);
        }
    };

//    void derivative(const GaugeField& Gauge) const {
//    };
};


NAMESPACE_END(Grid);