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


//#define USE_DOUBLE true 
#define USE_DOUBLE false 

#if USE_DOUBLE 
    #define PREC double
    typedef LatticeGaugeFieldD LGF;
    typedef PeriodicGimplD GIMPL;
    typedef vComplexD COMP;
#else
    #define PREC float
    typedef LatticeGaugeFieldF LGF;
    typedef PeriodicGimplF GIMPL;
    typedef vComplexF COMP;
#endif 


bool testForce(GridCartesian& GRID, LGF Umu, LGF Uforce, 
               LGF Ucontrol, PREC c1, PREC cnaik, PREC c3, PREC c5, PREC c7, PREC clp) {
    Smear_HISQ<GIMPL> hisq_fat(&GRID,c1,cnaik,c3,c5,c7,clp);
    LGF diff(&GRID);
    hisq_fat.ddVprojectU3(Uforce, Umu, Umu, 5e-5);
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
    pass *= testForce(GRID,Umu,Umu,Ucontrol,0.,0.,0.,0.,0.,0.);

    if(pass){
        Grid_pass("All tests passed.");
    } else {
        Grid_error("At least one test failed.");
    }

#endif

    Grid_finalize();
}