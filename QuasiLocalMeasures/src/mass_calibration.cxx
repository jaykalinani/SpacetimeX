#include <cctk.h>
#include <cctk_Functions.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

extern "C" CCTK_INT QuasiLocalMeasures_MeasureInitialMasses(
    CCTK_POINTER_TO_CONST cctkGH_pointer, const CCTK_INT npunctures,
    CCTK_REAL *const measured_mass, CCTK_INT *const valid) {
  DECLARE_CCTK_PARAMETERS;

  if (!cctkGH_pointer || npunctures < 0 || !measured_mass || !valid)
    return -1;
  for (int puncture = 0; puncture < npunctures; ++puncture) {
    measured_mass[puncture] =
        std::numeric_limits<CCTK_REAL>::quiet_NaN();
    valid[puncture] = 0;
  }
  if (!CCTK_IsFunctionAliased("AHFinderDirect_ForceInitialSearch"))
    return -1;

  std::vector<CCTK_INT> found(npunctures);
  std::vector<CCTK_INT> stored_surface(npunctures, -1);
  if (AHFinderDirect_ForceInitialSearch(cctkGH_pointer, npunctures,
                                        found.data(),
                                        stored_surface.data()) < 0)
    return -1;

  auto *const cctkGH = const_cast<cGH *>(
      static_cast<const cGH *>(cctkGH_pointer));
  CallScheduleGroup(cctkGH, "QuasiLocalMeasures_Calculate");
  const auto *const qlm_valid = static_cast<const CCTK_INT *>(
      CCTK_VarDataPtr(cctkGH, 0,
                      "QuasiLocalMeasures::qlm_have_valid_data[0]"));
  const auto *const qlm_mass = static_cast<const CCTK_REAL *>(
      CCTK_VarDataPtr(cctkGH, 0, "QuasiLocalMeasures::qlm_mass[0]"));
  if (!qlm_valid || !qlm_mass)
    return -1;

  for (int puncture = 0; puncture < npunctures; ++puncture) {
    const auto qlm_surface = std::find(
        surface_index, surface_index + num_surfaces,
        stored_surface[puncture]);
    if (found[puncture] && stored_surface[puncture] >= 0 &&
        qlm_surface != surface_index + num_surfaces) {
      const int qlm = int(qlm_surface - surface_index);
      if (qlm_valid[qlm] && std::isfinite(qlm_mass[qlm]) &&
          qlm_mass[qlm] > 0.0) {
        measured_mass[puncture] = qlm_mass[qlm];
        valid[puncture] = 1;
      }
    }
  }
  return 0;
}
