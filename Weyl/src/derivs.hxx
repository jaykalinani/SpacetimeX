#ifndef DERIVS_HXX
#define DERIVS_HXX

#include "tensor.hxx"

#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <initializer_list>
#include <ostream>
#include <sstream>
#include <type_traits>

namespace Weyl {
using namespace Loop;
using namespace std;

constexpr int deriv_order = 4;

////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline T deriv1d(const T *restrict const var, const ptrdiff_t di, const T dx) {
  switch (deriv_order) {
  case 2:
    return -1 / T(2) * (var[-di] - var[+di]) / dx;
  case 4:
    return (1 / T(12) * (var[-2 * di] - var[+2 * di]) //
            - 2 / T(3) * (var[-di] - var[+di])) /
           dx;
  default:
    assert(0);
  }
}

template <typename T>
inline T deriv2_1d(const T *restrict const var, const ptrdiff_t di,
                   const T dx) {
  switch (deriv_order) {
  case 2:
    return ((var[-di] + var[+di]) //
            - 2 * var[0]) /
           pow2(dx);
  case 4:
    return (-1 / T(12) * (var[-2 * di] + var[+2 * di]) //
            + 4 / T(3) * (var[-di] + var[+di])         //
            - 5 / T(2) * var[0]) /
           pow2(dx);
  default:
    assert(0);
  }
}

template <typename T>
inline T deriv2_2d(const T *restrict const var, const ptrdiff_t di,
                   const ptrdiff_t dj, const T dx, const T dy) {
  array<T, deriv_order + 1> arrx;
  T *const varx = &arrx[arrx.size() / 2];
  for (int j = -deriv_order / 2; j <= deriv_order / 2; ++j)
    varx[j] = deriv1d(&var[j * dj], di, dx);
  return deriv1d(varx, 1, dy);
}

////////////////////////////////////////////////////////////////////////////////

template <int dir, typename T>
inline T deriv(const GF3D<const T, 0, 0, 0> &gf_, const vect<int, dim> &I,
               const vec3<T, UP> &dx) {
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.offset(DI(dir));
  return deriv1d(&gf_(I), di, dx(dir));
}

template <int dir1, int dir2, typename T>
inline enable_if_t<(dir1 == dir2), T> deriv2(const GF3D<const T, 0, 0, 0> &gf_,
                                             const vect<int, dim> &I,
                                             const vec3<T, UP> &dx) {
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.offset(DI(dir1));
  return deriv2_1d(&gf_(I), di, dx(dir1));
}

template <int dir1, int dir2, typename T>
inline enable_if_t<(dir1 != dir2), T> deriv2(const GF3D<const T, 0, 0, 0> &gf_,
                                             const vect<int, dim> &I,
                                             const vec3<T, UP> &dx) {
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.offset(DI(dir1));
  const ptrdiff_t dj = gf_.offset(DI(dir2));
  return deriv2_2d(&gf_(I), di, dj, dx(dir1), dx(dir2));
}

////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline vec3<T, DN> deriv(const GF3D<const T, 0, 0, 0> &gf_,
                         const vect<int, dim> &I, const vec3<T, UP> &dx) {
  return {
      deriv<0>(gf_, I, dx),
      deriv<1>(gf_, I, dx),
      deriv<2>(gf_, I, dx),
  };
}

template <typename T>
inline mat3<T, DN, DN> deriv2(const GF3D<const T, 0, 0, 0> &gf_,
                              const vect<int, dim> &I, const vec3<T, UP> &dx) {
  return {
      deriv2<0, 0>(gf_, I, dx), deriv2<0, 1>(gf_, I, dx),
      deriv2<0, 2>(gf_, I, dx), deriv2<1, 1>(gf_, I, dx),
      deriv2<1, 2>(gf_, I, dx), deriv2<2, 2>(gf_, I, dx),
  };
}

////////////////////////////////////////////////////////////////////////////////

template <typename T>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs(const cGH *restrict const cctkGH, const GF3D<const T, 0, 0, 0> &gf_,
            const vec3<GF3D<T, 0, 0, 0>, DN> &dgf_) {
  DECLARE_CCTK_ARGUMENTS;

  const vec3<CCTK_REAL, UP> dx([&](int a) { return CCTK_DELTA_SPACE(a); });

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    const auto dval = deriv(gf_, p.I, dx);
    for (int a = 0; a < 3; ++a)
      dgf_(a)(p.I) = dval(a);
  });
}

template <typename T>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs2(const cGH *restrict const cctkGH,
             const GF3D<const T, 0, 0, 0> &gf_,
             const vec3<GF3D<T, 0, 0, 0>, DN> &dgf_,
             const mat3<GF3D<T, 0, 0, 0>, DN, DN> &ddgf_) {
  DECLARE_CCTK_ARGUMENTS;

  const vec3<CCTK_REAL, UP> dx([&](int a) { return CCTK_DELTA_SPACE(a); });

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    const auto dval = deriv(gf_, p.I, dx);
    for (int a = 0; a < 3; ++a)
      dgf_(a)(p.I) = dval(a);
    const auto ddval = deriv2(gf_, p.I, dx);
    for (int a = 0; a < 3; ++a)
      for (int b = a; b < 3; ++b)
        ddgf_(a, b)(p.I) = ddval(a, b);
  });
}

template <typename T, dnup_t dnup>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs(const cGH *restrict const cctkGH,
            const vec3<GF3D<const T, 0, 0, 0>, dnup> &gf_,
            const vec3<vec3<GF3D<T, 0, 0, 0>, DN>, dnup> &dgf_) {
  for (int a = 0; a < 3; ++a)
    calc_derivs(cctkGH, gf_(a), dgf_(a));
}

template <typename T, dnup_t dnup>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs2(const cGH *restrict const cctkGH,
             const vec3<GF3D<const T, 0, 0, 0>, dnup> &gf_,
             const vec3<vec3<GF3D<T, 0, 0, 0>, DN>, dnup> &dgf_,
             const vec3<mat3<GF3D<T, 0, 0, 0>, DN, DN>, dnup> &ddgf_) {
  for (int a = 0; a < 3; ++a)
    calc_derivs2(cctkGH, gf_(a), dgf_(a), ddgf_(a));
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs(const cGH *restrict const cctkGH,
            const mat3<GF3D<const T, 0, 0, 0>, dnup1, dnup2> &gf_,
            const mat3<vec3<GF3D<T, 0, 0, 0>, DN>, dnup1, dnup2> &dgf_) {
  for (int a = 0; a < 3; ++a)
    for (int b = a; b < 3; ++b)
      calc_derivs(cctkGH, gf_(a, b), dgf_(a, b));
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs2(const cGH *restrict const cctkGH,
             const mat3<GF3D<const T, 0, 0, 0>, dnup1, dnup2> &gf_,
             const mat3<vec3<GF3D<T, 0, 0, 0>, DN>, dnup1, dnup2> &dgf_,
             const mat3<mat3<GF3D<T, 0, 0, 0>, DN, DN>, dnup1, dnup2> &ddgf_) {
  for (int a = 0; a < 3; ++a)
    for (int b = a; b < 3; ++b)
      calc_derivs2(cctkGH, gf_(a, b), dgf_(a, b), ddgf_(a, b));
}

} // namespace Weyl

#endif // #ifndef DERIVS_HXX