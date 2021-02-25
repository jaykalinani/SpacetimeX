#ifndef LOOP_HXX
#define LOOP_HXX

#include "mempool.hxx"

#include <vect.hxx>

#include <cctk.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <limits>
#include <ostream>
#include <string>
#include <type_traits>

namespace Loop {
using namespace std;

template <typename F> CCTK_ATTRIBUTE_NOINLINE auto noinline(const F &f) {
  return f();
}

using Arith::vect;

constexpr int dim = 3;

enum class where_t { everywhere, interior, boundary, ghosts_inclusive, ghosts };

struct PointDesc {
  int imin, imax;
  int i, j, k;
  CCTK_REAL x, y, z;
  CCTK_REAL dx, dy, dz;
  int idx;
  static constexpr int di = 1;
  int dj; // stride in the j direction
  int dk; // stride in the k direction
  int np; // number of grid points
  vect<int, dim> I;
  vect<int, dim> NI; // outward boundary normal, or zero
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect<int, dim> DI(int d) const {
    return vect<int, dim>::unit(d);
  }
  vect<CCTK_REAL, dim> X;
  vect<CCTK_REAL, dim> DX;

  friend ostream &operator<<(ostream &os, const PointDesc &p) {
    os << "PointDesc{"
       << "imin,imax:{" << p.imin << "," << p.imax << "}, "
       << "ijk:"
       << "{" << p.i << "," << p.j << "," << p.k << "}, "
       << "xyz:"
       << "{" << p.x << "," << p.y << "," << p.z << "}, "
       << "dxyz:"
       << "{" << p.dx << "," << p.dy << "," << p.dz << "}, "
       << "idx:" << p.idx << ", "
       << "dijk:"
       << "{" << p.di << "," << p.dj << "," << p.dk << "},"
       << "np:" << p.np << ","
       << "nijk:"
       << "{" << p.NI[0] << "," << p.NI[1] << "," << p.NI[2] << "}"
       << "}";
    return os;
  }
};

struct GridDescBase {
  array<int, dim> gsh;
  array<int, dim> lbnd, ubnd;
  array<int, dim> lsh;
  array<int, dim> ash;
  array<int, 2 * dim> bbox;
  array<int, dim> nghostzones;
  array<int, dim> tmin, tmax;

  // for current level
  array<CCTK_REAL, dim> x0;
  array<CCTK_REAL, dim> dx;

  template <typename T, size_t N>
  static void output(ostream &os, const string &str, const array<T, N> &arr) {
    os << str << ":[";
    for (size_t n = 0; n < N; ++n) {
      if (n > 0)
        os << ",";
      os << arr[n];
    }
    os << "]";
  }
  friend ostream &operator<<(ostream &os, const GridDescBase &grid) {
    os << "GridDescBase{";
    output(os, "gsh", grid.gsh);
    output(os, ",lbnd", grid.lbnd);
    output(os, ",ubnd", grid.ubnd);
    output(os, ",lsh", grid.lsh);
    output(os, ",bbox", grid.bbox);
    output(os, ",nghostzones", grid.nghostzones);
    output(os, ",tmin", grid.tmin);
    output(os, ",tmax", grid.tmax);
    os << "}";
    return os;
  }

protected:
  GridDescBase();

public:
  GridDescBase(const cGH *cctkGH);

  template <int CI, int CJ, int CK>
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE PointDesc
  point_desc(const vect<int, dim> &restrict NI, const int imin, const int imax,
             const int i, const int j, const int k) const {
    constexpr int di = 1;
    const int dj = di * (ash[0] - CI);
    const int dk = dj * (ash[1] - CJ);
    const int np = dk * (ash[2] - CK);

    const CCTK_REAL x = x0[0] + (lbnd[0] + i - CCTK_REAL(!CI) / 2) * dx[0];
    const CCTK_REAL y = x0[1] + (lbnd[1] + j - CCTK_REAL(!CJ) / 2) * dx[1];
    const CCTK_REAL z = x0[2] + (lbnd[2] + k - CCTK_REAL(!CK) / 2) * dx[2];
    const int idx = i * di + j * dj + k * dk;
    return PointDesc{imin, imax,      i,     j,         k,   x,  y,
                     z,    dx[0],     dx[1], dx[2],     idx, dj, dk,
                     np,   {i, j, k}, NI,    {x, y, z}, dx};
  }

  template <int CI, int CJ, int CK>
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE PointDesc
  point_desc(const PointDesc &p) const {
    return point_desc<CI, CJ, CK>(p.NI, p.imin, p.imax, p.i, p.j, p.k);
  }

  // Loop over a given box
  template <int CI, int CJ, int CK, typename F>
  CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_box(const F &f, const array<int, dim> &restrict imin,
           const array<int, dim> &restrict imax,
           const array<int, dim> &restrict inormal) const {
    // TODO: convert array to vect (and everywhere nearby)
    static_assert(CI == 0 || CI == 1, "");
    static_assert(CJ == 0 || CJ == 1, "");
    static_assert(CK == 0 || CK == 1, "");

    for (int d = 0; d < dim; ++d)
      if (imin[d] >= imax[d])
        return;

    const auto kernel{[&](const int i, const int j,
                          const int k) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      f(point_desc<CI, CJ, CK>(inormal, imin[0], imax[0], i, j, k));
    }};

    array<bool, dim> bforward;
    for (int d = 0; d < dim; ++d)
      bforward[d] = inormal[d] >= 0;
    bool all_forward = true;
    for (int d = 0; d < dim; ++d)
      all_forward &= bforward[d];

    if (all_forward) {

      for (int k = imin[2]; k < imax[2]; ++k) {
        for (int j = imin[1]; j < imax[1]; ++j) {
#pragma omp simd
          for (int i = imin[0]; i < imax[0]; ++i) {
            kernel(i, j, k);
          }
        }
      }

    } else {
      // At least one direction is reversed; loop from the inside out

      for (int k0 = 0; k0 < imax[2] - imin[2]; ++k0) {
        for (int j0 = 0; j0 < imax[1] - imin[1]; ++j0) {
#pragma omp simd
          for (int i0 = 0; i0 < imax[0] - imin[0]; ++i0) {
            const int i = bforward[0] ? imin[0] + i0 : imax[0] - 1 - i0;
            const int j = bforward[1] ? imin[1] + j0 : imax[1] - 1 - j0;
            const int k = bforward[2] ? imin[2] + k0 : imax[2] - 1 - k0;
            kernel(i, j, k);
          }
        }
      }
    }
  }

  // Box including all points
  template <int CI, int CJ, int CK>
  void box_all(const array<int, dim> &group_nghostzones,
               vect<int, dim> &restrict imin,
               vect<int, dim> &restrict imax) const {
    const array<int, dim> offset{CI, CJ, CK};
    for (int d = 0; d < dim; ++d) {
      int ghost_offset = nghostzones[d] - group_nghostzones[d];
      imin[d] = std::max(tmin[d], ghost_offset);
      imax[d] = std::min(tmax[d], lsh[d] - offset[d] - ghost_offset);
    }
  }

  // Box including all interior points
  template <int CI, int CJ, int CK>
  void box_int(const array<int, dim> &group_nghostzones,
               vect<int, dim> &restrict imin,
               vect<int, dim> &restrict imax) const {
    const array<int, dim> offset{CI, CJ, CK};
    for (int d = 0; d < dim; ++d) {
      imin[d] = std::max(tmin[d], nghostzones[d]);
      imax[d] = std::min(tmax[d], lsh[d] - offset[d] - nghostzones[d]);
    }
  }

  // Loop over all points
  template <int CI, int CJ, int CK, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_all(const array<int, dim> &group_nghostzones, const F &f) const {
    // TODO: call box_all instead
    const array<int, dim> offset{CI, CJ, CK};
    array<int, dim> imin, imax;
    for (int d = 0; d < dim; ++d) {
      int ghost_offset = nghostzones[d] - group_nghostzones[d];
      imin[d] = std::max(tmin[d], ghost_offset);
      imax[d] = std::min(tmax[d], lsh[d] - offset[d] - ghost_offset);
    }
    const array<int, dim> inormal{0, 0, 0};

    loop_box<CI, CJ, CK>(f, imin, imax, inormal);
  }

  // Loop over all interior points
  template <int CI, int CJ, int CK, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_int(const array<int, dim> &group_nghostzones, const F &f) const {
    // TODO: call box_int instead
    const array<int, dim> offset{CI, CJ, CK};
    array<int, dim> imin, imax;
    for (int d = 0; d < dim; ++d) {
      imin[d] = std::max(tmin[d], nghostzones[d]);
      imax[d] = std::min(tmax[d], lsh[d] - offset[d] - nghostzones[d]);
    }
    const array<int, dim> inormal{0, 0, 0};

    loop_box<CI, CJ, CK>(f, imin, imax, inormal);
  }

  // Loop over all outer boundary points. This excludes ghost faces, but
  // includes ghost edges/corners on non-ghost faces. Loop over faces first,
  // then edges, then corners.
  template <int CI, int CJ, int CK, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_bnd(const array<int, dim> &group_nghostzones, const F &f) const {
    const array<int, dim> offset{CI, CJ, CK};

    for (int rank = dim - 1; rank >= 0; --rank) {

      for (int nk = -1; nk <= +1; ++nk) {
        for (int nj = -1; nj <= +1; ++nj) {
          for (int ni = -1; ni <= +1; ++ni) {
            if ((ni == 0) + (nj == 0) + (nk == 0) == rank) {

              if ((ni != 0 && bbox[0 + (ni == -1 ? 0 : 1)]) ||
                  (nj != 0 && bbox[2 + (nj == -1 ? 0 : 1)]) ||
                  (nk != 0 && bbox[4 + (nk == -1 ? 0 : 1)])) {

                const array<int, dim> inormal{ni, nj, nk};

                array<int, dim> imin, imax;
                for (int d = 0; d < dim; ++d) {
                  const int ghost_offset =
                      nghostzones[d] - group_nghostzones[d];
                  const int begin_bnd = ghost_offset;
                  const int begin_int = nghostzones[d];
                  const int end_int = lsh[d] - offset[d] - nghostzones[d];
                  const int end_bnd = lsh[d] - offset[d] - ghost_offset;
                  switch (inormal[d]) {
                  case -1: // lower boundary
                    imin[d] = begin_bnd;
                    imax[d] = begin_int;
                    break;
                  case 0: // interior
                    imin[d] = begin_int;
                    imax[d] = end_int;
                    break;
                  case +1: // upper boundary
                    imin[d] = end_int;
                    imax[d] = end_bnd;
                    break;
                  default:
                    assert(0);
                  }

                  imin[d] = std::max(tmin[d], imin[d]);
                  imax[d] = std::min(tmax[d], imax[d]);
                }

                loop_box<CI, CJ, CK>(f, imin, imax, inormal);
              }
            } // if rank
          }
        }
      }

    } // for rank
  }

  // Loop over all outer ghost points. This includes ghost edges/corners on
  // non-ghost faces. Loop over faces first, then edges, then corners.
  template <int CI, int CJ, int CK, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_ghosts_inclusive(const array<int, dim> &group_nghostzones,
                        const F &f) const {
    const array<int, dim> offset{CI, CJ, CK};

    for (int rank = dim - 1; rank >= 0; --rank) {

      for (int nk = -1; nk <= +1; ++nk) {
        for (int nj = -1; nj <= +1; ++nj) {
          for (int ni = -1; ni <= +1; ++ni) {
            if ((ni == 0) + (nj == 0) + (nk == 0) == rank) {

              if ((ni != 0 && !bbox[0 + (ni == -1 ? 0 : 1)]) ||
                  (nj != 0 && !bbox[2 + (nj == -1 ? 0 : 1)]) ||
                  (nk != 0 && !bbox[4 + (nk == -1 ? 0 : 1)])) {

                const array<int, dim> inormal{ni, nj, nk};

                array<int, dim> imin, imax;
                for (int d = 0; d < dim; ++d) {
                  const int ghost_offset =
                      nghostzones[d] - group_nghostzones[d];
                  const int begin_bnd = ghost_offset;
                  const int begin_int = nghostzones[d];
                  const int end_int = lsh[d] - offset[d] - nghostzones[d];
                  const int end_bnd = lsh[d] - offset[d] - ghost_offset;
                  switch (inormal[d]) {
                  case -1: // lower boundary
                    imin[d] = begin_bnd;
                    imax[d] = begin_int;
                    break;
                  case 0: // interior
                    imin[d] = begin_int;
                    imax[d] = end_int;
                    break;
                  case +1: // upper boundary
                    imin[d] = end_int;
                    imax[d] = end_bnd;
                    break;
                  default:
                    assert(0);
                  }

                  imin[d] = std::max(tmin[d], imin[d]);
                  imax[d] = std::min(tmax[d], imax[d]);
                }

                loop_box<CI, CJ, CK>(f, imin, imax, inormal);
              }
            } // if rank
          }
        }
      }
    } // for rank
  }

  // Loop over all outer ghost points. This excludes ghost edges/corners on
  // non-ghost faces. Loop over faces first, then edges, then corners.
  template <int CI, int CJ, int CK, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_ghosts(const array<int, dim> &group_nghostzones, const F &f) const {
    const array<int, dim> offset{CI, CJ, CK};

    for (int rank = dim - 1; rank >= 0; --rank) {

      for (int nk = -1; nk <= +1; ++nk) {
        for (int nj = -1; nj <= +1; ++nj) {
          for (int ni = -1; ni <= +1; ++ni) {
            if ((ni == 0) + (nj == 0) + (nk == 0) == rank) {

              if ((ni == 0 || !bbox[0 + (ni == -1 ? 0 : 1)]) &&
                  (nj == 0 || !bbox[2 + (nj == -1 ? 0 : 1)]) &&
                  (nk == 0 || !bbox[4 + (nk == -1 ? 0 : 1)])) {

                const array<int, dim> inormal{ni, nj, nk};

                array<int, dim> imin, imax;
                for (int d = 0; d < dim; ++d) {
                  const int ghost_offset =
                      nghostzones[d] - group_nghostzones[d];
                  const int begin_bnd = ghost_offset;
                  const int begin_int = nghostzones[d];
                  const int end_int = lsh[d] - offset[d] - nghostzones[d];
                  const int end_bnd = lsh[d] - offset[d] - ghost_offset;
                  switch (inormal[d]) {
                  case -1: // lower boundary
                    imin[d] = begin_bnd;
                    imax[d] = begin_int;
                    break;
                  case 0: // interior
                    imin[d] = begin_int;
                    imax[d] = end_int;
                    break;
                  case +1: // upper boundary
                    imin[d] = end_int;
                    imax[d] = end_bnd;
                    break;
                  default:
                    assert(0);
                  }

                  imin[d] = std::max(tmin[d], imin[d]);
                  imax[d] = std::min(tmax[d], imax[d]);
                }

                loop_box<CI, CJ, CK>(f, imin, imax, inormal);
              }
            } // if rank
          }
        }
      }
    } // for rank
  }

  template <int CI, int CJ, int CK, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop(where_t where, const array<int, dim> &group_nghostzones,
       const F &f) const {
    switch (where) {
    case where_t::everywhere:
      return noinline(
          [&] { return loop_all<CI, CJ, CK>(group_nghostzones, f); });
    case where_t::interior:
      return noinline(
          [&] { return loop_int<CI, CJ, CK>(group_nghostzones, f); });
    case where_t::boundary:
      return noinline(
          [&] { return loop_bnd<CI, CJ, CK>(group_nghostzones, f); });
    case where_t::ghosts_inclusive:
      return noinline([&] {
        return loop_ghosts_inclusive<CI, CJ, CK>(group_nghostzones, f);
      });
    case where_t::ghosts:
      return noinline(
          [&] { return loop_ghosts<CI, CJ, CK>(group_nghostzones, f); });
    default:
      assert(0);
    }
  }

  template <int CI, int CJ, int CK, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void loop(where_t where,
                                                const F &f) const {
    loop<CI, CJ, CK>(where, nghostzones, f);
  }

  template <typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_idx(where_t where, const array<int, dim> &indextype,
           const array<int, dim> &group_nghostzones, const F &f) const {
    switch (indextype[0] + 2 * indextype[1] + 4 * indextype[2]) {
    case 0b000:
      return noinline(
          [&] { return loop<0, 0, 0>(where, group_nghostzones, f); });
    case 0b001:
      return noinline(
          [&] { return loop<1, 0, 0>(where, group_nghostzones, f); });
    case 0b010:
      return noinline(
          [&] { return loop<0, 1, 0>(where, group_nghostzones, f); });
    case 0b011:
      return noinline(
          [&] { return loop<1, 1, 0>(where, group_nghostzones, f); });
    case 0b100:
      return noinline(
          [&] { return loop<0, 0, 1>(where, group_nghostzones, f); });
    case 0b101:
      return noinline(
          [&] { return loop<1, 0, 1>(where, group_nghostzones, f); });
    case 0b110:
      return noinline(
          [&] { return loop<0, 1, 1>(where, group_nghostzones, f); });
    case 0b111:
      return noinline(
          [&] { return loop<1, 1, 1>(where, group_nghostzones, f); });
    default:
      assert(0);
    }
  }

  template <typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_idx(where_t where, const array<int, dim> &indextype, const F &f) const {
    loop_idx(where, indextype, nghostzones, f);
  }
};

template <typename F>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
loop_idx(const cGH *cctkGH, where_t where, const array<int, dim> &indextype,
         const array<int, dim> &nghostzones, const F &f) {
  GridDescBase(cctkGH).loop_idx(where, indextype, nghostzones, f);
}

template <typename F>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
loop_idx(const cGH *cctkGH, where_t where, const array<int, dim> &indextype,
         const F &f) {
  GridDescBase(cctkGH).loop_idx(where, indextype, f);
}

template <int CI, int CJ, int CK, typename F>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE void loop(const cGH *cctkGH, where_t where,
                                              const F &f) {
  GridDescBase(cctkGH).loop<CI, CJ, CK>(where, f);
}

// Keep these for convenience
template <int CI, int CJ, int CK, typename F>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE void loop_all(const cGH *cctkGH,
                                                  const F &f) {
  loop<CI, CJ, CK>(cctkGH, where_t::everywhere, f);
}

template <int CI, int CJ, int CK, typename F>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE void loop_int(const cGH *cctkGH,
                                                  const F &f) {
  loop<CI, CJ, CK>(cctkGH, where_t::interior, f);
}

template <int CI, int CJ, int CK, typename F>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE void loop_bnd(const cGH *cctkGH,
                                                  const F &f) {
  loop<CI, CJ, CK>(cctkGH, where_t::boundary, f);
}

template <int CI, int CJ, int CK, typename F>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
loop_ghosts_inclusive(const cGH *cctkGH, const F &f) {
  loop<CI, CJ, CK>(cctkGH, where_t::ghosts_inclusive, f);
}

template <int CI, int CJ, int CK, typename F>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE void loop_ghosts(const cGH *cctkGH,
                                                     const F &f) {
  loop<CI, CJ, CK>(cctkGH, where_t::ghosts, f);
}

////////////////////////////////////////////////////////////////////////////////

extern bool CarpetX_poison_undefined_values;

template <typename T, int CI, int CJ, int CK> struct GF3D {
  static_assert(CI == 0 || CI == 1, "");
  static_assert(CJ == 0 || CJ == 1, "");
  static_assert(CK == 0 || CK == 1, "");
  typedef T value_type;
  static constexpr int di = 1;
  const int dj, dk, np;
  const int ni, nj, nk;
  T *restrict const ptr;
  static constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE array<int, dim> indextype() {
    return {CI, CJ, CK};
  }
  GF3D() = delete;
  GF3D(const GF3D &) = default;
  GF3D(GF3D &&) = default;
  GF3D &operator=(const GF3D &) = default;
  GF3D &operator=(GF3D &&) = default;
  GF3D(const cGH *restrict cctkGH, T *restrict ptr)
      : dj(di * (cctkGH->cctk_ash[0] - CI)),
        dk(dj * (cctkGH->cctk_ash[1] - CJ)),
        np(dk * (cctkGH->cctk_ash[2] - CK)), ni(cctkGH->cctk_lsh[0] - CI),
        nj(cctkGH->cctk_lsh[1] - CJ), nk(cctkGH->cctk_lsh[2] - CK), ptr(ptr) {}
  GF3D(const cGH *restrict cctkGH, mempool_t &mempool)
      : dj(di * (cctkGH->cctk_ash[0] - CI)),
        dk(dj * (cctkGH->cctk_ash[1] - CJ)),
        np(dk * (cctkGH->cctk_ash[2] - CK)), ni(cctkGH->cctk_lsh[0] - CI),
        nj(cctkGH->cctk_lsh[1] - CJ), nk(cctkGH->cctk_lsh[2] - CK),
        ptr(mempool.alloc<T>(np)) {}
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE int offset(int i, int j, int k) const {
    // These index checks prevent vectorization. We thus only enable
    // them in debug mode.
#ifdef CCTK_DEBUG
    assert(i >= 0 && i < ni);
    assert(j >= 0 && j < nj);
    assert(k >= 0 && k < nk);
#endif
    return i * di + j * dj + k * dk;
  }
  inline int offset(const vect<int, dim> &I) const {
    return offset(I[0], I[1], I[2]);
  }
  inline T &restrict operator()(int i, int j, int k) const {
    return ptr[offset(i, j, k)];
  }
  inline T &restrict operator()(const vect<int, dim> &I) const {
    return ptr[offset(I)];
  }
};

////////////////////////////////////////////////////////////////////////////////

template <typename T> struct GF3D1 {
  typedef T value_type;
  T *restrict ptr;
#ifdef CCTK_DEBUG
  array<int, dim> imin, imax;
  array<int, dim> ash;
#endif
  static constexpr int di = 1;
  int dj, dk, np;
  int off;
  GF3D1() = delete;
  GF3D1(const GF3D1 &) = default;
  GF3D1(GF3D1 &&) = default;
  GF3D1 &operator=(const GF3D1 &) = default;
  GF3D1 &operator=(GF3D1 &&) = default;
  GF3D1(T *restrict ptr, const array<int, dim> &imin,
        const array<int, dim> &imax, const array<int, dim> &ash)
      : ptr(ptr),
#ifdef CCTK_DEBUG
        imin(imin), imax(imax), ash(ash),
#endif
        dj(di * ash[0]), dk(dj * ash[1]), np(dk * ash[2]),
        off(imin[0] * di + imin[1] * dj + imin[2] * dk) {
  }
  GF3D1(const cGH *restrict cctkGH, const array<int, dim> &indextype,
        const array<int, dim> &nghostzones, T *restrict ptr) {
    for (int d = 0; d < dim; ++d)
      assert(indextype[d] == 0 || indextype[d] == 1);
    for (int d = 0; d < dim; ++d) {
      assert(nghostzones[d] >= 0);
      assert(nghostzones[d] <= cctkGH->cctk_nghostzones[d]);
    }
    array<int, dim> imin, imax;
    for (int d = 0; d < dim; ++d) {
      imin[d] = cctkGH->cctk_nghostzones[d] - nghostzones[d];
      imax[d] = cctkGH->cctk_lsh[d] - indextype[d] -
                (cctkGH->cctk_nghostzones[d] - nghostzones[d]);
    }
    array<int, dim> ash;
    for (int d = 0; d < dim; ++d)
      ash[d] = cctkGH->cctk_ash[d] - indextype[d] -
               2 * (cctkGH->cctk_nghostzones[d] - nghostzones[d]);
    *this = GF3D1(ptr, imin, imax, ash);
  }
  GF3D1(const cGH *restrict cctkGH, const array<int, dim> &indextype,
        const array<int, dim> &nghostzones, mempool_t &mempool)
      : GF3D1(cctkGH, indextype, nghostzones, nullptr) {
    ptr = mempool.alloc<T>(np);
  }
  inline int offset(int i, int j, int k) const {
    // These index checks prevent vectorization. We thus only enable
    // them in debug mode.
#ifdef CCTK_DEBUG
    assert(i >= imin[0] && i < imax[0]);
    assert(j >= imin[1] && j < imax[1]);
    assert(k >= imin[2] && k < imax[2]);
#endif
    return i * di + j * dj + k * dk - off;
  }
  inline int offset(const vect<int, dim> &I) const {
    return offset(I[0], I[1], I[2]);
  }
  inline T &restrict operator()(int i, int j, int k) const {
    return ptr[offset(i, j, k)];
  }
  inline T &restrict operator()(const vect<int, dim> &I) const {
    return ptr[offset(I)];
  }
  inline T &restrict operator()(const PointDesc &p) const {
    return (*this)(p.I);
  }
};

////////////////////////////////////////////////////////////////////////////////

struct GF3D2layout {
#ifdef CCTK_DEBUG
  vect<int, dim> imin, imax;
  vect<int, dim> ash;
#endif
  static constexpr int di = 1;
  int dj, dk, np;
  int off;
  GF3D2layout() = delete;
  GF3D2layout(const GF3D2layout &) = default;
  GF3D2layout(GF3D2layout &&) = default;
  GF3D2layout &operator=(const GF3D2layout &) = default;
  GF3D2layout &operator=(GF3D2layout &&) = default;
  GF3D2layout(const vect<int, dim> &imin, const vect<int, dim> &imax,
              const vect<int, dim> &ash)
      :
#ifdef CCTK_DEBUG
        imin(imin), imax(imax), ash(ash),
#endif
        dj(di * ash[0]), dk(dj * ash[1]), np(dk * ash[2]),
        off(imin[0] * di + imin[1] * dj + imin[2] * dk) {
  }
  GF3D2layout(const vect<int, dim> &imin, const vect<int, dim> &imax)
      : GF3D2layout(imin, imax, imax - imin) {}
  GF3D2layout(const cGH *restrict cctkGH, const vect<int, dim> &indextype,
              const vect<int, dim> &nghostzones) {
    for (int d = 0; d < dim; ++d)
      assert(indextype[d] == 0 || indextype[d] == 1);
    for (int d = 0; d < dim; ++d) {
      assert(nghostzones[d] >= 0);
      assert(nghostzones[d] <= cctkGH->cctk_nghostzones[d]);
    }
    vect<int, dim> imin, imax;
    for (int d = 0; d < dim; ++d) {
      imin[d] = cctkGH->cctk_nghostzones[d] - nghostzones[d];
      imax[d] = cctkGH->cctk_lsh[d] - indextype[d] -
                (cctkGH->cctk_nghostzones[d] - nghostzones[d]);
    }
    vect<int, dim> ash;
    for (int d = 0; d < dim; ++d)
      ash[d] = cctkGH->cctk_ash[d] - indextype[d] -
               2 * (cctkGH->cctk_nghostzones[d] - nghostzones[d]);
    *this = GF3D2layout(imin, imax, ash);
  }
  GF3D2layout(const cGH *restrict cctkGH, const vect<int, dim> &indextype)
      : GF3D2layout(cctkGH, indextype,
                    {cctkGH->cctk_nghostzones[0], cctkGH->cctk_nghostzones[1],
                     cctkGH->cctk_nghostzones[2]}) {}
  bool operator==(const GF3D2layout &other) const {
    bool iseq = true;
#ifdef CCTK_DEBUG
    iseq &= equal_to<vect<int, dim> >()(imin, other.imin) &&
            equal_to<vect<int, dim> >()(imax, other.imax);
    iseq &= equal_to<vect<int, dim> >()(ash, other.ash);
#endif
    iseq &= dj == other.dj && dk == other.dk && np == other.np;
    iseq &= off == other.off;
    return iseq;
  }
  inline int offset(int i, int j, int k) const {
    // These index checks prevent vectorization. We thus only enable
    // them in debug mode.
#ifdef CCTK_DEBUG
    assert(i >= imin[0] && i < imax[0]);
    assert(j >= imin[1] && j < imax[1]);
    assert(k >= imin[2] && k < imax[2]);
#endif
    return i * di + j * dj + k * dk - off;
  }
  inline int offset(const vect<int, dim> &I) const {
    return offset(I[0], I[1], I[2]);
  }
};

template <typename T> struct GF3D2 {
  typedef T value_type;
  // TODO: disallow inf, nan
  // Note: Keep `ptr` as first member, this improves performance a bit (?)
  T *restrict ptr;
  GF3D2layout layout;
  GF3D2() = delete;
  GF3D2(const GF3D2 &) = default;
  GF3D2(GF3D2 &&) = default;
  GF3D2 &operator=(const GF3D2 &) = default;
  GF3D2 &operator=(GF3D2 &&) = default;
  GF3D2(const GF3D2layout &layout, T *restrict ptr)
      : ptr(ptr), layout(layout) {}
  GF3D2(const GF3D2layout &layout, mempool_t &mempool)
      : GF3D2(layout, mempool.alloc<T>(layout.np)) {}
  inline int offset(int i, int j, int k) const {
    return layout.offset(i, j, k);
  }
  inline int offset(const vect<int, dim> &I) const { return layout.offset(I); }
  inline T &restrict operator()(int i, int j, int k) const {
    return ptr[offset(i, j, k)];
  }
  inline T &restrict operator()(const vect<int, dim> &I) const {
    return ptr[offset(I)];
  }
};

////////////////////////////////////////////////////////////////////////////////

template <int NI, int NJ, int NK, int OFF = 0> struct GF3D3layout {
  static_assert(NI >= 0, "");
  static_assert(NJ >= 0, "");
  static_assert(NK >= 0, "");

  static constexpr int ni = NI;
  static constexpr int nj = NJ;
  static constexpr int nk = NK;

  static constexpr int off = OFF;

  static constexpr int di = 1;
  static constexpr int dj = NI * di;
  static constexpr int dk = NJ * dj;
  static constexpr int np = NK * dk;

  constexpr int offset(int i, int j, int k) const {
    return i * di + j * dj + k * dk - OFF;
  }
  constexpr int offset(const vect<int, dim> &I) const {
    return (*this)(I[0], I[1], I[2]);
  }
};

template <int imin0, int imin1, int imin2, int imax0, int imax1, int imax2>
struct makeGF3D3layout {
private:
  static constexpr int ni = imax0 - imin0;
  static constexpr int nj = imax1 - imin1;
  static constexpr int nk = imax2 - imin2;
  static constexpr int di = GF3D3layout<ni, nj, nk>::di;
  static constexpr int dj = GF3D3layout<ni, nj, nk>::dj;
  static constexpr int dk = GF3D3layout<ni, nj, nk>::dk;
  static constexpr int off = imin0 * di + imin1 * dj + imin2 * dk;

public:
  typedef GF3D3layout<ni, nj, nk, off> type;
};
template <int imin0, int imin1, int imin2, int imax0, int imax1, int imax2>
using makeGF3D3layout_t =
    typename makeGF3D3layout<imin0, imin1, imin2, imax0, imax1, imax2>::type;

template <typename T, int NI, int NJ, int NK, int OFF = 0>
struct GF3D3 : GF3D3layout<NI, NJ, NK, OFF> {
  using GF3D3layout<NI, NJ, NK, OFF>::np;
  using GF3D3layout<NI, NJ, NK, OFF>::offset;

  array<T, np> arr;

  constexpr T &restrict operator()(int i, int j, int k) {
    return arr[offset(i, j, k)];
  }
  constexpr const T &restrict operator()(int i, int j, int k) const {
    return arr[offset(i, j, k)];
  }
  constexpr T &restrict operator()(const vect<int, dim> &I) {
    return (*this)(I[0], I[1], I[2]);
  }
  constexpr const T &restrict operator()(const vect<int, dim> &I) const {
    return (*this)(I[0], I[1], I[2]);
  }
};

template <typename T, int NI, int NJ, int NK, int OFF = 0>
struct GF3D3ptr : GF3D3layout<NI, NJ, NK, OFF> {
  using GF3D3layout<NI, NJ, NK, OFF>::np;
  using GF3D3layout<NI, NJ, NK, OFF>::offset;

  T *restrict ptr;

  GF3D3ptr() = delete;
  GF3D3ptr(T *restrict ptr) : ptr(ptr) {}

  constexpr T &restrict operator()(int i, int j, int k) const {
    return ptr[offset(i, j, k)];
  }
  constexpr T &restrict operator()(const vect<int, dim> &I) const {
    return (*this)(I[0], I[1], I[2]);
  }
};

////////////////////////////////////////////////////////////////////////////////

typedef GF3D2layout GF3D5layout;

template <typename T> struct GF3D5 {
#ifdef CCTK_DEBUG
  GF3D5layout layout;
#endif
  T *restrict ptr;
  GF3D5() = delete;
  GF3D5(const GF3D5 &) = default;
  GF3D5(GF3D5 &&) = default;
  GF3D5 &operator=(const GF3D5 &) = default;
  GF3D5 &operator=(GF3D5 &&) = default;
  GF3D5(const GF3D2layout &layout, T *restrict ptr)
      :
#ifdef CCTK_DEBUG
        layout(layout),
#endif
        ptr(ptr) {
  }
  GF3D5(const GF3D2layout &layout, mempool_t &mempool)
      : GF3D5(layout, mempool.alloc<T>(layout.np)) {}
  constexpr T &restrict operator()(const GF3D5layout &layout, int i, int j,
                                   int k) const {
#ifdef CCTK_DEBUG
    assert(layout == this->layout);
#endif
    return ptr[layout.offset(i, j, k)];
  }
  constexpr T &restrict operator()(const GF3D5layout &layout,
                                   const vect<int, dim> &I) const {
    return (*this)(layout, I[0], I[1], I[2]);
  }
};

} // namespace Loop

#endif // #ifndef LOOP_HXX
