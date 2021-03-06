##############################################################################
###   Copyright 2003 and onward LASMEA UMR 6602 CNRS/U.B.P Clermont-Ferrand
###   Copyright 2009 and onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
###
###          Distributed under the Boost Software License, Version 1.0
###                 See accompanying file LICENSE.txt or copy at
###                     http://www.boost.org/LICENSE_1_0.txt
##############################################################################

SET( SOURCES
  colon.cpp
  iscolumn.cpp
  isempty.cpp
  ismatrix.cpp
  isrow.cpp
  isscalar.cpp
  isvector.cpp
  isinside.cpp
  length.cpp
  linspace.cpp
  first_last_index.cpp
  ndims.cpp
  numel.cpp
  ones.cpp
  reshape.cpp
  squeeze.cpp
  ind2sub.cpp
  sub2ind.cpp
  sub2sub.cpp
  size.cpp
  zeros.cpp
  issquare.cpp
  issymetric.cpp
  ishermitian.cpp
  istriangular.cpp
  istril.cpp
  istriu.cpp
  eye.cpp
  triu.cpp
  tril.cpp
  band.cpp
  tri1u.cpp
  tri1l.cpp
  ric.cpp
  rif.cpp
  rows.cpp
  resize.cpp
  expand.cpp
  streaming.cpp
  cic.cpp
  cif.cpp
  cols.cpp
  colvect.cpp
  rowvect.cpp
  fliplr.cpp
  width.cpp
  logspace.cpp
  pispace.cpp
  freqspace1.cpp
  freqspace2.cpp
  sum.cpp
  asum2.cpp
  asum1.cpp
  asump.cpp
  maximum.cpp
  minimum.cpp
  firstnonsingleton.cpp
  tie.cpp
  all.cpp
  any.cpp
  none.cpp
  nbtrue.cpp
  prod.cpp
  mean.cpp
  norm2.cpp
  normp.cpp
#  isequal.cpp
#  isequaln.cpp
  havesamesize.cpp
  isinteger.cpp
  isfloating.cpp
  dot.cpp
# List of test files for toolbox core.settings
   )
set(CMAKE_CXX_FLAGS_NT2TEST "${CMAKE_CXX_FLAGS_NT2TEST} -Wfatal-errors")
nt2_module_add_tests(core.functions.unit ${SOURCES})
