/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/OverlapWilsonPartialFractionZolotarevFermion.h

    Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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

    See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#ifndef OVERLAP_WILSON_PARTFRAC_ZOLOTAREV_FERMION_H
#define OVERLAP_WILSON_PARTFRAC_ZOLOTAREV_FERMION_H

#include <Grid/qcd/action/fermion/FermionCore.h>

NAMESPACE_BEGIN(Grid);

template<class Impl>
class OverlapWilsonPartialFractionZolotarevFermion : public PartialFractionFermion5D<Impl>
{
public:
  INHERIT_IMPL_TYPES(Impl);

  virtual void   Instantiatable(void){};

  void  MomentumSpacePropagator(FermionField &out,const FermionField &in,RealD _m,std::vector<double> twist) {
    this->MomentumSpacePropagatorHw(out,in,_m,twist);
  };

  // Constructors
  OverlapWilsonPartialFractionZolotarevFermion(GaugeField &_Umu,
					       GridCartesian         &FiveDimGrid,
					       GridRedBlackCartesian &FiveDimRedBlackGrid,
					       GridCartesian         &FourDimGrid,
					       GridRedBlackCartesian &FourDimRedBlackGrid,
					       RealD _mass,RealD _M5,
					       RealD lo,RealD hi,const ImplParams &p= ImplParams()):
      
    // b+c=scale, b-c = 0 <=> b =c = scale/2
    PartialFractionFermion5D<Impl>(_Umu,
				   FiveDimGrid,
				   FiveDimRedBlackGrid,
				   FourDimGrid,
				   FourDimRedBlackGrid,_mass,_M5,p)
  {
    assert((this->Ls&0x1)==1); // Odd Ls required

    int nrational=this->Ls;// Odd rational order
    RealD eps = lo/hi;

    Approx::zolotarev_data *zdata = Approx::zolotarev(eps,nrational,0);
    this->SetCoefficientsZolotarev(hi,zdata);
    Approx::zolotarev_free(zdata);

  }
};

NAMESPACE_END(Grid);

#endif
