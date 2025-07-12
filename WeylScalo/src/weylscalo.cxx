#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <loop_device.hxx>
#include <cx_utils.hxx>

#include <cmath>

#include <cx_derivsGF3D5.hxx>
#include <cx_derivsinline.hxx>
#include <cx_dissinline.hxx>
#include <cx_powerinline.hxx>

namespace WeylScalo {
using namespace Arith;
using namespace Loop;
using namespace CXUtils;

CCTK_DEVICE CCTK_HOST constexpr CCTK_REAL Sqrt(CCTK_REAL x) {
  return std::sqrt(x);
}

CCTK_DEVICE CCTK_HOST constexpr CCTK_REAL Max(CCTK_REAL x, CCTK_REAL y) {
  return std::max(x, y);
}

extern "C" void WeylScalo_calc_psi4(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WeylScalo_calc_psi4;
  DECLARE_CCTK_PARAMETERS;

  for (int d = 0; d < 3; ++d)
    if (cctk_nghostzones[d] < deriv_order / 2)
      CCTK_VERROR("Need at least %d ghost zones", deriv_order / 2);

  const array<CCTK_REAL, 3> invDxyz{1. / CCTK_DELTA_SPACE(0),
                                    1. / CCTK_DELTA_SPACE(1),
                                    1. / CCTK_DELTA_SPACE(2)};

  const GF3D2layout layout2(cctkGH, {0, 0, 0});

  // Input grid functions
  const array<const CCTK_REAL *, 6> gf_gam{gxx, gxy, gxz, gyy, gyz, gzz};
  const array<const CCTK_REAL *, 6> gf_exK{kxx, kxy, kxz, kyy, kyz, kzz};

  // Output grid functions
  CCTK_REAL *gf_Psi4real = Psi4r;
  CCTK_REAL *gf_Psi4imag = Psi4i;

  // Parameters
  // const CCTK_REAL oosqrt2 = 1.0 / sqrt(2);
  constexpr CCTK_REAL ncutoff = 1e-16;

  // Loop
  const Loop::GridDescBaseDevice grid(cctkGH);

  if (calc_derivs_live) {

    // Derivs Lambdas
#include "../wolfram/WeylScalo_derivs1st.hxx"
#include "../wolfram/WeylScalo_derivs2nd.hxx"

    if (elecmag_weyl) {
#include "../wolfram/WeylScalo_set_Psi4.hxx"
    } else {
#include "../wolfram/WeylScalo_set_Psi4_GaussCodazzi.hxx"
    }

  } else {

    // Tile variables for derivatives and so on
    const int ntmps = 72;
    int itmp = 0;
    const GF3D5layout layout5 = CXUtils::get_GF3D5layout<0, 0, 0>(cctkGH);
    CXUtils::GF3D5Factory<CCTK_REAL> fct(layout5, ntmps, itmp);

    const auto tl_dgam = fct.make_smat_vec_gf();
    const auto tl_ddgam = fct.make_smat_smat_gf();
    const auto tl_dexK = fct.make_smat_vec_gf();

    if (itmp != ntmps)
      CCTK_VERROR("Wrong number of temporary variables: ntmps=%d itmp=%d",
                  ntmps, itmp);
    itmp = -1;

    // Define derivs lambdas
    const auto calcderivs1st = [&](const auto &dgf, const auto &gf_) {
      calc_derivs1st<0, 0, 0>(grid, layout5, dgf, layout2, gf_, invDxyz,
                              deriv_order);
    };
    const auto calcderivs2nd = [&](const auto &dgf, const auto &ddgf,
                                   const auto &gf_) {
      calc_derivs2nd<0, 0, 0>(grid, layout5, dgf, ddgf, layout2, gf_, invDxyz,
                              deriv_order);
    };

    calcderivs2nd(tl_dgam, tl_ddgam, gf_gam);
    calcderivs1st(tl_dexK, gf_exK);

    if (elecmag_weyl) {
#include "../wolfram/WeylScalo_set_Psi4_GF3D5.hxx"
    } else {
#include "../wolfram/WeylScalo_set_Psi4_GaussCodazzi_GF3D5.hxx"
    }
  }
}

extern "C" void WeylScalo_Sync(CCTK_ARGUMENTS) {
  // do nothing
}

} // namespace WeylScalo
