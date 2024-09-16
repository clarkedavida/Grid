/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/smearing/Test_HISQ_force.cc

Copyright (C) 2024

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
    @file Test_HISQ_force.cc
    @brief test of the HISQ fermion force 
*/


#include <Grid/Grid.h>
#include <Grid/lattice/PaddedCell.h>
#include <Grid/stencil/GeneralLocalStencil.h>
#include <Grid/qcd/smearing/HISQSmearing.h>
using namespace Grid;


#define USE_DOUBLE true 
//#define USE_DOUBLE false 

#if USE_DOUBLE 
    #define PREC double
    typedef LatticeGaugeFieldD LGF;
    typedef StaggeredImplD GIMPL;
    typedef vComplexD COMP;
#else
    #define PREC float
    typedef LatticeGaugeFieldF LGF;
    typedef StaggeredImplF GIMPL;
    typedef vComplexF COMP;
#endif 


// Intent: IN--Umu: thin links
bool testForce(GridCartesian& GRID, LGF Umu, LGF Uforce, LGF Ucontrol) {

    int n_naiks = 1;
    std::array<PREC,GRID_MAX_NAIK> eps_naik = {0,0,0}; 

//    PREC fat7_c1  = 1/8.;                  PREC asqtad_c1  = 1.;
//    PREC fat7_c3  = 1/16.;                 PREC asqtad_c3  = 1/16.;
//    PREC fat7_c5  = 1/64.;                 PREC asqtad_c5  = 1/64.;
//    PREC fat7_c7  = 1/384.;                PREC asqtad_c7  = 1/384.;
//    PREC cnaik    = -1/24.+eps_naik[0]/8;  PREC asqtad_clp = -1/8.;

    PREC fat7_c1  = 1.;  PREC asqtad_c1  = 0.;
    PREC fat7_c3  = 0.;  PREC asqtad_c3  = 1/16.;
    PREC fat7_c5  = 0.;  PREC asqtad_c5  = 0.;
    PREC fat7_c7  = 0.;  PREC asqtad_c7  = 0.;
    PREC cnaik    = 0.;  PREC asqtad_clp = 0.;

    LGF Vmu(&GRID), Wmu(&GRID), Nmu(&GRID);
    Smear_HISQ<GIMPL> fat7(&GRID,fat7_c1,0.,fat7_c3,fat7_c5,fat7_c7,0.);

    fat7.smear(Vmu,Nmu,Umu); // Populate fat7 and Naik links Vmu and Nmu
    fat7.projectU3(Wmu,Vmu); // Populate U(3) projection Wmu

    HISQParameters<PREC> hisq_param(n_naiks  , eps_naik ,
                                    fat7_c1  , fat7_c3  , fat7_c5  , fat7_c7  , 0.,
                                    asqtad_c1, asqtad_c3, asqtad_c5, asqtad_c7, asqtad_clp,
                                    cnaik    , 0.       , 0.);

    HISQReunitSVDParameters<PREC> hisq_reunit_svd(false, false, 1, 1, 1);

    Force_HISQ<GIMPL> hisq_force(&GRID, hisq_param, Wmu, Vmu, Umu, hisq_reunit_svd);

    LGF diff(&GRID);
    hisq_force.ddVprojectU3(Uforce, Umu, Umu, 5e-5);
    bool result;
    diff = Ucontrol-Uforce;
    auto absDiff = norm2(diff)/norm2(Ucontrol);
    if (absDiff < 1e-30) {
        Grid_pass(" |Umu-Usmr|/|Umu| = ",absDiff);
        result = true;
    } else {
        Grid_error(" |Umu-Usmr|/|Umu| = ",absDiff);
        result = false;
    }
//    NerscIO::writeConfiguration(Uforce,"nersc.l8t4b3360.ddVU3");
    return result;
}


int main (int argc, char** argv) {

    // Params for the test.
    int Ns = 8;
    int Nt = 4;
    Coordinate latt_size(Nd,0); latt_size[0]=Ns; latt_size[1]=Ns; latt_size[2]=Ns; latt_size[3]=Nt;
    std::string conf_in  = "nersc.l8t4b3360";
    int threads          = GridThread::GetThreads();

    // Initialize the Grid
    Grid_init(&argc,&argv);
    Coordinate simd_layout = GridDefaultSimd(Nd,COMP::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();
    Grid_log("mpi     = ",mpi_layout);
    Grid_log("simd    = ",simd_layout);
    Grid_log("latt    = ",latt_size);
    Grid_log("threads = ",threads);
    GridCartesian GRID(latt_size,simd_layout,mpi_layout);

    LGF Umu(&GRID), Ucontrol(&GRID);

#if USE_DOUBLE // NerscIO is hard-coded to double.

    // Read the configuration into Umu
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, conf_in);

    bool pass=true;

    NerscIO::readConfiguration(Ucontrol, header, "nersc.l8t4b3360.ddVU3.control");
    pass *= testForce(GRID,Umu,Umu,Ucontrol);

    if(pass){
        Grid_pass("All tests passed.");
    } else {
        Grid_error("At least one test failed.");
    }

#endif

    Grid_finalize();
}