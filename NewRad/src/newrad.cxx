#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

#include <cmath>

namespace NewRad {
using namespace Loop;
using namespace std;

namespace {
template <typename T> constexpr T pow2(const T x) { return x * x; }
} // namespace

// Adapted from BSSN_MoL's files NewRad.F and newrad.h. This code was probably
// originally written by Miguel Alcubierre.
void newrad(const cGH *restrict const cctkGH,
            const CCTK_REAL *restrict const var, CCTK_REAL *restrict const rhs,
            const CCTK_REAL var0, //!< value at infinity
            const CCTK_REAL v0    //!< propagation speed
) {
  DECLARE_CCTK_ARGUMENTS;

  constexpr vect<int, dim> DI{1, 0, 0};
  constexpr vect<int, dim> DJ{0, 1, 0};
  constexpr vect<int, dim> DK{0, 0, 1};

  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  const GF3D<const CCTK_REAL, 0, 0, 0> var_(cctkGH, var);
  const GF3D<CCTK_REAL, 0, 0, 0> rhs_(cctkGH, rhs);

  const auto derivx{
      [&](const GF3D<const CCTK_REAL, 0, 0, 0> &u_, const PointDesc &p) {
        const auto I = p.I;
        if (p.NI[0] == 0)
          // interior
          return (u_(I + DI) - u_(I - DI)) / (2 * dx);
        if (p.NI[0] == +1)
          // upper boundary
          return +(3 * u_(I) - 4 * u_(I - DI) + u_(I - 2 * DI)) / (2 * dx);
        if (p.NI[0] == -1)
          // lower boundary
          return -(3 * u_(I) - 4 * u_(I + DI) + u_(I + 2 * DI)) / (2 * dx);
        assert(0);
      }};
  const auto derivy{
      [&](const GF3D<const CCTK_REAL, 0, 0, 0> &u_, const PointDesc &p) {
        const auto I = p.I;
        if (p.NI[1] == 0)
          return (u_(I + DJ) - u_(I - DJ)) / (2 * dy);
        if (p.NI[1] == +1)
          return +(3 * u_(I) - 4 * u_(I - DJ) + u_(I - 2 * DJ)) / (2 * dy);
        if (p.NI[1] == -1)
          return -(3 * u_(I) - 4 * u_(I + DJ) + u_(I + 2 * DJ)) / (2 * dy);
        assert(0);
      }};
  const auto derivz{
      [&](const GF3D<const CCTK_REAL, 0, 0, 0> &u_, const PointDesc &p) {
        const auto I = p.I;
        if (p.NI[2] == 0)
          return (u_(I + DK) - u_(I - DK)) / (2 * dz);
        if (p.NI[2] == +1)
          return +(3 * u_(I) - 4 * u_(I - DK) + u_(I - 2 * DK)) / (2 * dz);
        if (p.NI[2] == -1)
          return -(3 * u_(I) - 4 * u_(I + DK) + u_(I + 2 * DK)) / (2 * dz);
        assert(0);
      }};

  const Loop::GridDescBaseDevice grid(cctkGH);
  grid.loop_bnd_device<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    // The main part of the boundary condition assumes that we have an
    // outgoing radial wave with some speed v0:
    //
    //    var  =  var0 + u(r-v0*t)/r
    //
    // This implies the following differential equation:
    //
    //    d_t var  =  - v^i d_i var  -  v0 (var - var0) / r
    //
    // where  vi = v0 xi/r

    const CCTK_REAL r = sqrt(pow2(p.x) + pow2(p.y) + pow2(p.z));

    // Find local wave speeds
    const CCTK_REAL vx = v0 * p.x / r;
    const CCTK_REAL vy = v0 * p.y / r;
    const CCTK_REAL vz = v0 * p.z / r;
    // CCTK_REAL const vr = sqrt(pow2(vx) + pow2(vy) + pow2(vz));

    // Derivatives
    const CCTK_REAL varx = derivx(var_, p);
    const CCTK_REAL vary = derivy(var_, p);
    const CCTK_REAL varz = derivz(var_, p);

    // Calculate source term
    rhs_(p.I) =
        -vx * varx - vy * vary - vz * varz - v0 * (var_(p.I) - var0) / r;
  });
}

extern "C" void NewRad_Apply(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_NewRad_Apply;

  newrad(cctkGH, chi, chi_rhs, 0, 1);
  newrad(cctkGH, gammatxx, gammatxx_rhs, 0, 1);
  newrad(cctkGH, gammatxy, gammatxy_rhs, 0, 1);
  newrad(cctkGH, gammatxz, gammatxz_rhs, 0, 1);
  newrad(cctkGH, gammatyy, gammatyy_rhs, 0, 1);
  newrad(cctkGH, gammatyz, gammatyz_rhs, 0, 1);
  newrad(cctkGH, gammatzz, gammatzz_rhs, 0, 1);
  newrad(cctkGH, Kh, Kh_rhs, 0, 1);
  newrad(cctkGH, Atxx, Atxx_rhs, 0, 1);
  newrad(cctkGH, Atxy, Atxy_rhs, 0, 1);
  newrad(cctkGH, Atxz, Atxz_rhs, 0, 1);
  newrad(cctkGH, Atyy, Atyy_rhs, 0, 1);
  newrad(cctkGH, Atyz, Atyz_rhs, 0, 1);
  newrad(cctkGH, Atzz, Atzz_rhs, 0, 1);
  newrad(cctkGH, Gamtx, Gamtx_rhs, 0, 1);
  newrad(cctkGH, Gamty, Gamty_rhs, 0, 1);
  newrad(cctkGH, Gamtz, Gamtz_rhs, 0, 1);
  newrad(cctkGH, Theta, Theta_rhs, 0, 1);
  newrad(cctkGH, alphaG, alphaG_rhs, 1, 1);
  newrad(cctkGH, betaG, betaG_rhs, 0, 1);
}

} // namespace NewRad
