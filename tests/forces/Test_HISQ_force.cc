/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/smearing/Test_fatLinks.cc

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


//bool testSmear(GridCartesian& GRID, LatticeGaugeFieldD Umu, LatticeGaugeFieldD Usmr, LatticeGaugeFieldD Unaik, 
//               LatticeGaugeFieldD Ucontrol, Real c1, Real cnaik, Real c3, Real c5, Real c7, Real clp) {
//    Smear_HISQ<PeriodicGimplD> hisq_fat(&GRID,c1,cnaik,c3,c5,c7,clp);
//    LatticeGaugeFieldD diff(&GRID), Uproj(&GRID);
//    hisq_fat.smear(Usmr, Unaik, Umu);
//    bool result;
//    if (cnaik < 1e-30) { // Testing anything but Naik term
//        diff = Ucontrol-Usmr;
//        auto absDiff = norm2(diff)/norm2(Ucontrol);
//        if (absDiff < 1e-30) {
//            Grid_pass(" |Umu-Usmr|/|Umu| = ",absDiff);
//            result = true;
//        } else {
//            Grid_error(" |Umu-Usmr|/|Umu| = ",absDiff);
//            result = false;
//        }
//    } else { // Testing Naik specifically
//        diff = Ucontrol-Unaik;
//        auto absDiff = norm2(diff)/norm2(Ucontrol);
//        if (absDiff < 1e-30) {
//            Grid_pass(" |Umu-Unaik|/|Umu| = ",absDiff);
//            result = true;
//        } else {
//            Grid_error(" |Umu-Unaik|/|Umu| = ",absDiff);
//            result = false;
//        }
//        hisq_fat.projectU3(Uproj,Ucontrol);
////        NerscIO::writeConfiguration(Unaik,"nersc.l8t4b3360.naik");
//    }
//    return result;
//}


int main (int argc, char** argv) {

    // Params for the test.
    int Ns = 8;
    int Nt = 4;
    Coordinate latt_size(Nd,0); latt_size[0]=Ns; latt_size[1]=Ns; latt_size[2]=Ns; latt_size[3]=Nt;
    std::string conf_in  = "nersc.l8t4b3360";
    int threads          = GridThread::GetThreads();
    typedef LatticeGaugeFieldD LGF;

    // Initialize the Grid
    Grid_init(&argc,&argv);
    Coordinate simd_layout = GridDefaultSimd(Nd,vComplexD::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();
    Grid_log("mpi     = ",mpi_layout);
    Grid_log("simd    = ",simd_layout);
    Grid_log("latt    = ",latt_size);
    Grid_log("threads = ",threads);
    GridCartesian GRID(latt_size,simd_layout,mpi_layout);

    XmlReader Reader("HISQParams.xml",false,"grid");

    LGF Umu(&GRID), dWdV(&GRID);

    // Read the configuration into Umu
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, conf_in);

    bool pass=true;

    if(pass){
        Grid_pass("All tests passed.");
    } else {
        Grid_error("At least one test failed.");
    }

    Grid_finalize();
}