   /*************************************************************************************
    Grid physics library, www.github.com/paboyle/Grid 
    Source file: ./lib/tensors/Tensor_traits.h
    Copyright (C) 2015
Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Christopher Kelly <ckelly@phys.columbia.edu>
Author: Michael Marshall <michael.marshall@ed.ac.au>
Author: Christoph Lehner <christoph@lhnr.de>
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
#pragma once 

#include <type_traits>

NAMESPACE_BEGIN(Grid);

  // Forward declarations
  template<class T>        class iScalar;
  template<class T, int N> class iVector;
  template<class T, int N> class iMatrix;

  // These are the Grid tensors
  template<typename T>     struct isGridTensor                : public std::false_type { static constexpr bool notvalue = true; };
  template<class T>        struct isGridTensor<iScalar<T> >   : public std::true_type  { static constexpr bool notvalue = false; };
  template<class T, int N> struct isGridTensor<iVector<T, N> >: public std::true_type  { static constexpr bool notvalue = false; };
  template<class T, int N> struct isGridTensor<iMatrix<T, N> >: public std::true_type  { static constexpr bool notvalue = false; };

  template <typename T>  using IfGridTensor    = Invoke<std::enable_if<isGridTensor<T>::value, int> >;
  template <typename T>  using IfNotGridTensor = Invoke<std::enable_if<!isGridTensor<T>::value, int> >;

  // Traits to identify scalars
  template<typename T>     struct isGridScalar                : public std::false_type { static constexpr bool notvalue = true; };
  template<class T>        struct isGridScalar<iScalar<T>>    : public std::true_type  { static constexpr bool notvalue = false; };


  // Traits to identify fundamental data types
  template<typename T>     struct isGridFundamental                : public std::false_type { static constexpr bool notvalue = true; };
  template<>               struct isGridFundamental<vComplexF>     : public std::true_type  { static constexpr bool notvalue = false; };
  template<>               struct isGridFundamental<vComplexD>     : public std::true_type  { static constexpr bool notvalue = false; };
  template<>               struct isGridFundamental<vRealF>        : public std::true_type  { static constexpr bool notvalue = false; };
  template<>               struct isGridFundamental<vRealD>        : public std::true_type  { static constexpr bool notvalue = false; };
  template<>               struct isGridFundamental<ComplexF>      : public std::true_type  { static constexpr bool notvalue = false; };
  template<>               struct isGridFundamental<ComplexD>      : public std::true_type  { static constexpr bool notvalue = false; };
  template<>               struct isGridFundamental<RealF>         : public std::true_type  { static constexpr bool notvalue = false; };
  template<>               struct isGridFundamental<RealD>         : public std::true_type  { static constexpr bool notvalue = false; };
  template<>               struct isGridFundamental<vComplexD2>    : public std::true_type  { static constexpr bool notvalue = false; };
  template<>               struct isGridFundamental<vRealD2>       : public std::true_type  { static constexpr bool notvalue = false; };


//////////////////////////////////////////////////////////////////////////////////
// Want to recurse: GridTypeMapper<Matrix<vComplexD> >::scalar_type == ComplexD.
// Use of a helper class like this allows us to template specialise and "dress"
// other classes such as RealD == double, ComplexD == std::complex<double> with these
// traits.
//
// It is possible that we could do this more elegantly if I introduced a 
// queryable trait in iScalar, iMatrix and iVector and used the query on vtype in 
// place of the type mapper?
//
// Not sure how to do this, but probably could be done with a research effort
// to study C++11's type_traits.h file. (std::enable_if<isGridTensorType<vtype> >)
//
//////////////////////////////////////////////////////////////////////////////////

  // This saves repeating common properties for supported Grid Scalar types
  // TensorLevel    How many nested grid tensors
  // Rank           Rank of the grid tensor
  // count          Total number of elements, i.e. product of dimensions
  // Dimension(dim) Size of dimension dim
  struct GridTypeMapper_Base {
    static constexpr int TensorLevel = 0;
    static constexpr int Rank = 0;
    static constexpr std::size_t count = 1;
    static constexpr int Dimension(int dim) { return 0; }
  };

//////////////////////////////////////////////////////////////////////////////////
// Recursion stops with these template specialisations
//////////////////////////////////////////////////////////////////////////////////

  template<typename T> struct GridTypeMapper {};

  template<> struct GridTypeMapper<RealF> : public GridTypeMapper_Base {
    typedef RealF scalar_type;
    typedef RealD scalar_typeD;
    typedef RealF vector_type;
    typedef RealD vector_typeD;
    typedef RealF tensor_reduced ;
    typedef RealF scalar_object;
    typedef RealD scalar_objectD;
    typedef ComplexF Complexified;
    typedef RealF Realified;
    typedef RealD DoublePrecision;
    typedef RealD DoublePrecision2;
  };
  template<> struct GridTypeMapper<RealD> : public GridTypeMapper_Base {
    typedef RealD scalar_type;
    typedef RealD scalar_typeD;
    typedef RealD vector_type;
    typedef RealD vector_typeD;
    typedef RealD tensor_reduced;
    typedef RealD scalar_object;
    typedef RealD scalar_objectD;
    typedef ComplexD Complexified;
    typedef RealD Realified;
    typedef RealD DoublePrecision;
    typedef RealD DoublePrecision2;
  };
  template<> struct GridTypeMapper<ComplexF> : public GridTypeMapper_Base {
    typedef ComplexF scalar_type;
    typedef ComplexD scalar_typeD;
    typedef ComplexF vector_type;
    typedef ComplexD vector_typeD;
    typedef ComplexF tensor_reduced;
    typedef ComplexF scalar_object;
    typedef ComplexD scalar_objectD;
    typedef ComplexF Complexified;
    typedef RealF Realified;
    typedef ComplexD DoublePrecision;
    typedef ComplexD DoublePrecision2;
  };
  template<> struct GridTypeMapper<ComplexD> : public GridTypeMapper_Base {
    typedef ComplexD scalar_type;
    typedef ComplexD scalar_typeD;
    typedef ComplexD vector_type;
    typedef ComplexD vector_typeD;
    typedef ComplexD tensor_reduced;
    typedef ComplexD scalar_object;
    typedef ComplexD scalar_objectD;
    typedef ComplexD Complexified;
    typedef RealD Realified;
    typedef ComplexD DoublePrecision;
    typedef ComplexD DoublePrecision2;
  };

#if defined(GRID_CUDA) || defined(GRID_HIP)  
  template<> struct GridTypeMapper<std::complex<float> > : public GridTypeMapper_Base {
    typedef std::complex<float> scalar_type;
    typedef std::complex<double> scalar_typeD;
    typedef scalar_type vector_type;
    typedef scalar_typeD vector_typeD;
    typedef scalar_type tensor_reduced;
    typedef scalar_type scalar_object;
    typedef scalar_typeD scalar_objectD;
    typedef scalar_type Complexified;
    typedef RealF Realified;
    typedef scalar_typeD DoublePrecision;
    typedef scalar_typeD DoublePrecision2;
  };
  template<> struct GridTypeMapper<std::complex<double> > : public GridTypeMapper_Base {
    typedef std::complex<double> scalar_type;
    typedef std::complex<double> scalar_typeD;
    typedef scalar_type vector_type;
    typedef scalar_typeD vector_typeD;
    typedef scalar_type tensor_reduced;
    typedef scalar_type scalar_object;
    typedef scalar_typeD scalar_objectD;
    typedef scalar_type Complexified;
    typedef RealD Realified;
    typedef scalar_typeD DoublePrecision;
    typedef scalar_typeD DoublePrecision2;
  };
#endif

  template<> struct GridTypeMapper<Integer> : public GridTypeMapper_Base {
    typedef Integer scalar_type;
    typedef Integer scalar_typeD;
    typedef Integer vector_type;
    typedef Integer vector_typeD;
    typedef Integer tensor_reduced;
    typedef Integer scalar_object;
    typedef Integer scalar_objectD;
    typedef void Complexified;
    typedef void Realified;
    typedef void DoublePrecision;
    typedef void DoublePrecision2;
  };

  template<> struct GridTypeMapper<vRealF> : public GridTypeMapper_Base {
    typedef RealF  scalar_type;
    typedef RealD  scalar_typeD;
    typedef vRealF vector_type;
    typedef vRealD vector_typeD;
    typedef vRealF tensor_reduced;
    typedef RealF  scalar_object;
    typedef RealD  scalar_objectD;
    typedef vComplexF Complexified;
    typedef vRealF Realified;
    typedef vRealD DoublePrecision;
    typedef vRealD2 DoublePrecision2;
  };
  template<> struct GridTypeMapper<vRealD> : public GridTypeMapper_Base {
    typedef RealD  scalar_type;
    typedef RealD  scalar_typeD;
    typedef vRealD vector_type;
    typedef vRealD vector_typeD;
    typedef vRealD tensor_reduced;
    typedef RealD  scalar_object;
    typedef RealD  scalar_objectD;
    typedef vComplexD Complexified;
    typedef vRealD Realified;
    typedef vRealD DoublePrecision;
    typedef vRealD DoublePrecision2;
  };
  template<> struct GridTypeMapper<vRealD2> : public GridTypeMapper_Base {
    typedef RealD  scalar_type;
    typedef RealD  scalar_typeD;
    typedef vRealD2 vector_type;
    typedef vRealD2 vector_typeD;
    typedef vRealD2 tensor_reduced;
    typedef RealD  scalar_object;
    typedef RealD  scalar_objectD;
    typedef vComplexD2 Complexified;
    typedef vRealD2 Realified;
    typedef vRealD2 DoublePrecision;
    typedef vRealD2 DoublePrecision2;
  };
  template<> struct GridTypeMapper<vRealH> : public GridTypeMapper_Base {
    // Fixme this is incomplete until Grid supports fp16 or bfp16 arithmetic types
    typedef RealF  scalar_type;
    typedef RealD  scalar_typeD;
    typedef vRealH vector_type;
    typedef vRealD vector_typeD;
    typedef vRealH tensor_reduced;
    typedef RealF  scalar_object;
    typedef RealD  scalar_objectD;
    typedef vComplexH Complexified;
    typedef vRealH Realified;
    typedef vRealD DoublePrecision;
    typedef vRealD DoublePrecision2;
  };
  template<> struct GridTypeMapper<vComplexH> : public GridTypeMapper_Base {
    // Fixme this is incomplete until Grid supports fp16 or bfp16 arithmetic types
    typedef ComplexF  scalar_type;
    typedef ComplexD  scalar_typeD;
    typedef vComplexH vector_type;
    typedef vComplexD vector_typeD;
    typedef vComplexH tensor_reduced;
    typedef ComplexF  scalar_object;
    typedef ComplexD  scalar_objectD;
    typedef vComplexH Complexified;
    typedef vRealH Realified;
    typedef vComplexD DoublePrecision;
    typedef vComplexD DoublePrecision2;
  };
  template<> struct GridTypeMapper<vComplexF> : public GridTypeMapper_Base {
    typedef ComplexF  scalar_type;
    typedef ComplexD  scalar_typeD;
    typedef vComplexF vector_type;
    typedef vComplexD vector_typeD;
    typedef vComplexF tensor_reduced;
    typedef ComplexF  scalar_object;
    typedef ComplexD  scalar_objectD;
    typedef vComplexF Complexified;
    typedef vRealF Realified;
    typedef vComplexD DoublePrecision;
    typedef vComplexD2 DoublePrecision2;
  };
  template<> struct GridTypeMapper<vComplexD> : public GridTypeMapper_Base {
    typedef ComplexD  scalar_type;
    typedef ComplexD  scalar_typeD;
    typedef vComplexD vector_type;
    typedef vComplexD vector_typeD;
    typedef vComplexD tensor_reduced;
    typedef ComplexD  scalar_object;
    typedef ComplexD  scalar_objectD;
    typedef vComplexD Complexified;
    typedef vRealD Realified;
    typedef vComplexD DoublePrecision;
    typedef vComplexD DoublePrecision2;
  };
  template<> struct GridTypeMapper<vComplexD2> : public GridTypeMapper_Base {
    typedef ComplexD  scalar_type;
    typedef ComplexD  scalar_typeD;
    typedef vComplexD2 vector_type;
    typedef vComplexD2 vector_typeD;
    typedef vComplexD2 tensor_reduced;
    typedef ComplexD  scalar_object;
    typedef ComplexD  scalar_objectD;
    typedef vComplexD2 Complexified;
    typedef vRealD2 Realified;
    typedef vComplexD2 DoublePrecision;
    typedef vComplexD2 DoublePrecision2;
  };
  template<> struct GridTypeMapper<vInteger> : public GridTypeMapper_Base {
    typedef  Integer scalar_type;
    typedef  Integer scalar_typeD;
    typedef vInteger vector_type;
    typedef vInteger vector_typeD;
    typedef vInteger tensor_reduced;
    typedef  Integer scalar_object;
    typedef  Integer scalar_objectD;
    typedef void Complexified;
    typedef void Realified;
    typedef void DoublePrecision;
    typedef void DoublePrecision2;
  };

#define GridTypeMapper_RepeatedTypes \
  using BaseTraits   = GridTypeMapper<T>; \
  using scalar_type  = typename BaseTraits::scalar_type; \
  using vector_type  = typename BaseTraits::vector_type; \
  using scalar_typeD = typename BaseTraits::scalar_typeD; \
  using vector_typeD = typename BaseTraits::vector_typeD; \
  static constexpr int TensorLevel = BaseTraits::TensorLevel + 1

  template<typename T> struct GridTypeMapper<iScalar<T>> {
    GridTypeMapper_RepeatedTypes;
    using tensor_reduced  = iScalar<typename BaseTraits::tensor_reduced>;
    using scalar_object   = iScalar<typename BaseTraits::scalar_object>;
    using scalar_objectD  = iScalar<typename BaseTraits::scalar_objectD>;
    using Complexified    = iScalar<typename BaseTraits::Complexified>;
    using Realified       = iScalar<typename BaseTraits::Realified>;
    using DoublePrecision = iScalar<typename BaseTraits::DoublePrecision>;
    using DoublePrecision2= iScalar<typename BaseTraits::DoublePrecision2>;
    static constexpr int Rank = BaseTraits::Rank + 1;
    static constexpr std::size_t count = BaseTraits::count;
    static constexpr int Dimension(int dim) {
      return ( dim == 0 ) ? 1 : BaseTraits::Dimension(dim - 1); }
  };

  template<typename T, int N> struct GridTypeMapper<iVector<T, N>> {
    GridTypeMapper_RepeatedTypes;
    using tensor_reduced  = iScalar<typename BaseTraits::tensor_reduced>;
    using scalar_object   = iVector<typename BaseTraits::scalar_object,   N>;
    using scalar_objectD  = iVector<typename BaseTraits::scalar_objectD,  N>;
    using Complexified    = iVector<typename BaseTraits::Complexified,    N>;
    using Realified       = iVector<typename BaseTraits::Realified,       N>;
    using DoublePrecision = iVector<typename BaseTraits::DoublePrecision, N>;
    using DoublePrecision2= iVector<typename BaseTraits::DoublePrecision2, N>;
    static constexpr int Rank = BaseTraits::Rank + 1;
    static constexpr std::size_t count = BaseTraits::count * N;
    static constexpr int Dimension(int dim) {
      return ( dim == 0 ) ? N : BaseTraits::Dimension(dim - 1); }
  };

  template<typename T, int N> struct GridTypeMapper<iMatrix<T, N>> {
    GridTypeMapper_RepeatedTypes;
    using tensor_reduced  = iScalar<typename BaseTraits::tensor_reduced>;
    using scalar_object   = iMatrix<typename BaseTraits::scalar_object,   N>;
    using scalar_objectD  = iMatrix<typename BaseTraits::scalar_objectD,  N>;
    using Complexified    = iMatrix<typename BaseTraits::Complexified,    N>;
    using Realified       = iMatrix<typename BaseTraits::Realified,       N>;
    using DoublePrecision = iMatrix<typename BaseTraits::DoublePrecision, N>;
    using DoublePrecision2= iMatrix<typename BaseTraits::DoublePrecision2, N>;
    static constexpr int Rank = BaseTraits::Rank + 2;
    static constexpr std::size_t count = BaseTraits::count * N * N;
    static constexpr int Dimension(int dim) {
      return ( dim == 0 || dim == 1 ) ? N : BaseTraits::Dimension(dim - 2); }
  };

  // Match the index
  template<typename T,int Level> struct matchGridTensorIndex {
    static const bool value = (Level==T::TensorLevel);
    static const bool notvalue = (Level!=T::TensorLevel);
  };
  // What is the vtype
  template<typename T> struct isComplex {
    static const bool value = false;
  };
  template<> struct isComplex<ComplexF> {
    static const bool value = true;
  };
  template<> struct isComplex<ComplexD> {
    static const bool value = true;
  };

  //Get the SIMD vector type from a Grid tensor or Lattice<Tensor>
  template<typename T>
  struct getVectorType{
    typedef T type;
  };
  
  //Query whether a tensor or Lattice<Tensor> is SIMD vector or scalar
  template<typename T, typename V=void> struct isSIMDvectorized : public std::false_type {};
  template<typename U> struct isSIMDvectorized<U, typename std::enable_if< !std::is_same<
    typename GridTypeMapper<typename getVectorType<U>::type>::scalar_type,
    typename GridTypeMapper<typename getVectorType<U>::type>::vector_type>::value, void>::type>
  : public std::true_type {};

  //Get the precision of a Lattice, tensor or scalar type in units of sizeof(float)
  template<typename T>
  class getPrecision{
  public:
    //get the vector_obj (i.e. a grid Tensor) if its a Lattice<vobj>, do nothing otherwise (i.e. if fundamental or grid Tensor)
    typedef typename getVectorType<T>::type vector_obj; 
    typedef typename GridTypeMapper<vector_obj>::scalar_type scalar_type; //get the associated scalar type. Works on fundamental and tensor types
    typedef typename GridTypeMapper<scalar_type>::Realified real_scalar_type; //remove any std::complex wrapper, should get us to the fundamental type

    enum { value = sizeof(real_scalar_type)/sizeof(float) };
  };
NAMESPACE_END(Grid);



