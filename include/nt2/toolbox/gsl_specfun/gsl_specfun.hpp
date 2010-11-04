//////////////////////////////////////////////////////////////////////////////
///   Copyright 2003 and onward LASMEA UMR 6602 CNRS/U.B.P Clermont-Ferrand
///   Copyright 2009 and onward LRI    UMR 8623 CNRS/Univ Paris Sud XI
///
///          Distributed under the Boost Software License, Version 1.0
///                 See accompanying file LICENSE.txt or copy at
///                     http://www.boost.org/LICENSE_1_0.txt
//////////////////////////////////////////////////////////////////////////////
#ifndef NT2_TOOLBOX_GSL_SPECFUN_GSL_SPECFUN_HPP_INCLUDED
#define NT2_TOOLBOX_GSL_SPECFUN_GSL_SPECFUN_HPP_INCLUDED

//<include> please don't modify between these tags
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_sinc.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_beta_inc.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_i1_scaled.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_log_1plusx_mx.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_beta.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_In_scaled.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_log_erfc.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_k1_scaled.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_fermi_dirac_mhalf.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_gamma_inc.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_erfc.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_jl.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_fermi_dirac_3half.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_lambert_Wm1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_hyperg_2F1_conj.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_lnbeta.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_exprel.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_lnchoose.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_atanint.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_conicalP_1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_conicalP_0.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_poch.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_hypot.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_legendre_H3d.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_hyperg_2F1_renorm.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_fermi_dirac_0.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_fermi_dirac_1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_zeta.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_psi_1piy.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_Knu.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_j1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_expm1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_expint_E1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_expint_E2.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_fermi_dirac_half.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_log_1plusx.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_Inu.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_dawson.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_fermi_dirac_2.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_clausen.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_lncosh.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_exprel_n.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_pochrel.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_lndoublefact.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_Klu_scaled_.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_hyperg_0F1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_j0.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_Inu_scaled_.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_legendre_sphPlm.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_fermi_dirac_int.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_hyperg_2F1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_hyperg_2F0.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_yl.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_fermi_dirac_m1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_exp.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_gammainv.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_Shi.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_gamma_inc_P.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_angle_restrict_symm.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_i0_scaled.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_K0_scaled.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_debye_4.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_debye_5.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_debye_6.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_debye_1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_debye_2.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_debye_3.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_legendre_H3d_0.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_legendre_H3d_1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_lnpoch.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_fermi_dirac_inc_0.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_gamma.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_lambert_W0.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_I0.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_I1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_dilog.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_log.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_conicalP_cyl_reg.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_expint_Ei.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_taylorcoeff.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_hzeta.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_fact.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_Yn.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_exprel_2.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_conicalP_half.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_hazard.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_hyperg_U.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_psi_n.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_gammastar.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_angle_restrict_pos.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_lnsinh.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_k0_scaled.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_erf_Z.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_erf.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_erf_Q.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_K1_scaled.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_In.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_doublefact.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_il_scaled.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_hyperg_2F1_conj_renorm.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_Jnu.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_psi_1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_lnfact.hpp>
//#include <nt2/toolbox/gsl_specfun/include/gsl_sf_choose.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_gegenpoly_2.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_legendre_P3.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_Knu_scaled_.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_expint_E1_scaled.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_Y0.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_Y1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_I0_scaled.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_log_abs.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_zetam1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_laguerre_1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_laguerre_3.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_laguerre_2.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_y0.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_psi.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_K0.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_K1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_Ynu.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_gegenpoly_3.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_Jn.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_gegenpoly_1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_cos.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_lngamma.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_y1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_kl_scaled.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_hyperg_1F1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_legendre_Q0.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_legendre_Q1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_I1_scaled.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_expint_Ei_scaled.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_legendre_Pl.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_Kn_scaled.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_expint_E2_scaled.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_gamma_inc_Q.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_gegenpoly_n.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_J1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_J0.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_Kn.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_exp_mult.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_conicalP_sph_reg.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_Ci.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_laguerre_n.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_legendre_P1.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_bessel_Klu.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_legendre_P2.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_sin.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_hyperg_1F1_int.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_legendre_Plm.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_conicalP_mhalf.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_legendre_Ql.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_Chi.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_eta.hpp>
#include <nt2/toolbox/gsl_specfun/include/gsl_sf_expint_3.hpp>
//<\include>

#endif
