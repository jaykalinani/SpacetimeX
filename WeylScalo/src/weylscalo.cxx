#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <loop_device.hxx>
#include <stx_utils.hxx>

#include <cmath>

#include <stx_derivsGF3D5.hxx>
#include <stx_derivsinline.hxx>
#include <stx_dissinline.hxx>
#include <stx_powerinline.hxx>

namespace WeylScalo {
using namespace Arith;
using namespace Loop;
using namespace STXUtils;

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

  // Derivs Lambdas
#include "../wolfram/WeylScalo_derivs1st.hxx"
#include "../wolfram/WeylScalo_derivs2nd.hxx"

  if(elecmag_weyl) {
#include "../wolfram/WeylScalo_set_Psi4.hxx"
  } else {
#include "../wolfram/WeylScalo_set_Psi4_GaussCodazzi.hxx"
  }
}

extern "C" void WeylScalo_Sync(CCTK_ARGUMENTS) {
  // do nothing
}

} // namespace WeylScalo
