#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Sync.h>

namespace WeylScal4 {

extern "C" void WeylScal4_SyncPsi4(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WeylScal4_SyncPsi4;
  DECLARE_CCTK_PARAMETERS;

  // The Psi4 calculation itself uses this cadence guard.  Keeping the
  // synchronization guard identical avoids an otherwise unconditional
  // collective ghost exchange on iterations where Psi4 was not updated.
  if (cctk_iteration % WeylScal4_psi4_calc_4th_calc_every !=
      WeylScal4_psi4_calc_4th_calc_offset)
    return;

  const int groups[] = {
      CCTK_GroupIndex("WeylScal4::Psi4i_group"),
      CCTK_GroupIndex("WeylScal4::Psi4r_group"),
  };
  for (const int group : groups)
    if (group < 0)
      CCTK_VERROR("Could not find a Psi4 grid-function group");

  const int ierr = CCTK_SyncGroupsI(cctkGH, 2, groups);
  if (ierr < 0)
    CCTK_VERROR("Could not synchronize Psi4 grid-function groups: error %d",
                ierr);
}

} // namespace WeylScal4
