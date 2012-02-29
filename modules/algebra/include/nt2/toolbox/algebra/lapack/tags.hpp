#ifndef __NT2__CORE__ENTITY__COMMON__SETTINGS__INCLUDED__
#define __NT2__CORE__ENTITY__COMMON__SETTINGS__INCLUDED__

////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003-2008 LASMEA UMR 6602 du CNRS.
// Copyright (c) 2007-2008 IEF    UMR 8622 du CNRS.
// All rights reserved.
//
// License information are available in the LICENSE file.
// Additionnal informations are available in the INFOS file.
////////////////////////////////////////////////////////////////////////////////

namespace nt2
{
  //////////////////////////////////////////////////////////////////////////////
  /// @struct dense
  /// @brief  Dense shape tag
  ///
  /// This tag is used to trigger the dense shape settings for the shape class.
  /// A dense shape is supposed to contains no trivial elements and behave as
  /// a standard N-dimensions connex iteration space.
  //////////////////////////////////////////////////////////////////////////////
  struct dense { typedef shape_tag type; };
  NT2_PP_REGISTER_TYPE_ID(dense);

  //////////////////////////////////////////////////////////////////////////////
  /// @struct diagonal
  /// @brief  Diagonal shape tag
  ///
  /// This tag is used to trigger the diagonal shape settings for the shape class.
  /// A diagonal shape is supposed to contains non trivial element on its main
  /// diagonal.
  //////////////////////////////////////////////////////////////////////////////
  struct diagonal { typedef shape_tag type; };
  NT2_PP_REGISTER_TYPE_ID(diagonal);

  //////////////////////////////////////////////////////////////////////////////
  /// @struct packed
  /// @brief  packed shape tag
  ///
  /// This tag is used to trigger the packed shape settings for
  /// the shape class.
  /// A lower_triangular shape is supposed to contains non trivial element
  /// on its lower triangular part.
  //////////////////////////////////////////////////////////////////////////////
  struct packed { typedef shape_tag type; };
  NT2_PP_REGISTER_TYPE_ID(packed);
  
  //////////////////////////////////////////////////////////////////////////////
  /// @struct banded
  /// @brief  banded shape tag
  ///
  /// This tag is used to trigger the banded shape settings for
  /// the shape class.
  /// A banded shape is supposed to contains non trivial element
  /// on its upper triangular part.
  //////////////////////////////////////////////////////////////////////////////
  struct band { typedef shape_tag type; };
  NT2_PP_REGISTER_TYPE_ID(band);

  //////////////////////////////////////////////////////////////////////////////
  /// @struct tridiagonal
  /// @brief  tridiagonal shape tag
  ///
  /// This tag is used to trigger the tridiagonal shape settings for
  /// the shape class.
  /// A tridiagonal shape is supposed to contains non trivial element
  /// on its part above the first subdiagonal.
  //////////////////////////////////////////////////////////////////////////////
  struct tridiagonal { typedef shape_tag type; };
  NT2_PP_REGISTER_TYPE_ID(tridiagonal);

  //////////////////////////////////////////////////////////////////////////////
  /// @struct bidiagonal
  /// @brief  bidiagonal shape tag
  ///
  /// This tag is used to trigger the bidiagonal shape settings for
  /// the shape class.
  /// A bidiagonal shape is supposed to contains non trivial element
  /// on its upper triangular part
  //////////////////////////////////////////////////////////////////////////////
  struct bidiagonal { typedef shape_tag type; };
  NT2_PP_REGISTER_TYPE_ID(bidiagonal);

// *BD 	bidiagonal
// *DI 	diagonal
// *GB 	general band
// GE 	general (i.e., unsymmetric, in some cases rectangular)
// GG 	general matrices, generalized problem (i.e., a pair of general matrices)
// *GT 	general tridiagonal
// HB 	(complex) Hermitian band
// HE 	(complex) Hermitian
// HG 	upper Hessenberg matrix, generalized problem (i.e a Hessenberg and a triangular matrix)
// HP 	(complex) Hermitian, packed storage
// HS 	upper Hessenberg
// OP 	(real) orthogonal, packed storage
// OR 	(real) orthogonal
// PB 	symmetric or Hermitian positive definite band
// PO 	symmetric or Hermitian positive definite
// PP 	symmetric or Hermitian positive definite, packed storage
// PT 	symmetric or Hermitian positive definite tridiagonal
// SB 	(real) symmetric band
// SP 	symmetric, packed storage
// ST 	(real) symmetric tridiagonal
// SY 	symmetric
// TB 	triangular band
// TG 	triangular matrices, generalized problem (i.e., a pair of triangular matrices)
// TP 	triangular, packed storage
// *TR 	triangular (or in some cases quasi-triangular)
// TZ 	trapezoidal
// *UN 	(complex) unitary
// UP 	(complex) unitary, packed storage
  
//   struct bidiagonal { typedef shape_tag type; };
//   NT2_PP_REGISTER_TYPE_ID(bidiagonal);
  
//   struct general_band { typedef shape_tag type; };
//   NT2_PP_REGISTER_TYPE_ID(general_band);
  
//   struct general_tridiagonal { typedef shape_tag type; };
//   NT2_PP_REGISTER_TYPE_ID(general_tridiagonal);
  
//   struct Hermitian { typedef shape_tag type; };
//   NT2_PP_REGISTER_TYPE_ID(Hermitian);
  
//   struct hermitian_band { typedef shape_tag type; };
//   NT2_PP_REGISTER_TYPE_ID(hermitian_band);
  
//   struct unitary { typedef shape_tag type; };
//   NT2_PP_REGISTER_TYPE_ID(unitary);
  
//   struct bidiagonal { typedef shape_tag type; };
//   NT2_PP_REGISTER_TYPE_ID(triangular);
  
//   struct bidiagonal { typedef shape_tag type; };
//   NT2_PP_REGISTER_TYPE_ID(bidiagonal);
//   struct bidiagonal { typedef shape_tag type; };
//   NT2_PP_REGISTER_TYPE_ID(bidiagonal);
//   struct bidiagonal { typedef shape_tag type; };
//   NT2_PP_REGISTER_TYPE_ID(bidiagonal);
//   struct bidiagonal { typedef shape_tag type; };
//   NT2_PP_REGISTER_TYPE_ID(bidiagonal);
//   struct bidiagonal { typedef shape_tag type; };
//   NT2_PP_REGISTER_TYPE_ID(bidiagonal);
//   struct bidiagonal { typedef shape_tag type; };
//   NT2_PP_REGISTER_TYPE_ID(bidiagonal);
//   struct bidiagonal { typedef shape_tag type; };
//   NT2_PP_REGISTER_TYPE_ID(bidiagonal);
//   struct bidiagonal { typedef shape_tag type; };
//   NT2_PP_REGISTER_TYPE_ID(bidiagonal);
//   struct bidiagonal { typedef shape_tag type; };
//   NT2_PP_REGISTER_TYPE_ID(bidiagonal);
//   struct bidiagonal { typedef shape_tag type; };

}

#endif

