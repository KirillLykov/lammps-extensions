//  (C) Copyright Kirill Lykov 2013.
//
// Distributed under the GNU Software License (See accompanying file LICENSE)

#include "region_complement.h"
#include "stdlib.h"
#include "string.h"
#include "domain.h"
#include "error.h"
#include <assert.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

RegComplement::RegComplement(LAMMPS *lmp, int narg, char **arg) :
Region(lmp, narg, arg), iregion(-1)
{
  if (narg < 3) error->all(FLERR,"Illegal region command");
  options(narg - 3, &arg[3]);

  iregion = domain->find_region(arg[2]);
  if (iregion == -1) error->all(FLERR,"Region complement region ID does not exist");

  bboxflag = 0;
}

/* ---------------------------------------------------------------------- */

RegComplement::~RegComplement()
{
}

/* ----------------------------------------------------------------------
 return 1 if region is dynamic, 0 if static
 dynamic if any sub-region is dynamic, else static
 ------------------------------------------------------------------------- */

int RegComplement::dynamic_check()
{
  return domain->regions[iregion]->dynamic_check();
}

/* ----------------------------------------------------------------------
 inside = 1 if x,y,z is match() with all sub-regions
 else inside = 0
 ------------------------------------------------------------------------- */

int RegComplement::inside(double x, double y, double z)
{
    return !domain->regions[iregion]->match(x,y,z);
}

int RegComplement::surface_interior(double *x, double cutoff)
{
  assert(false);
  return 0;
}

int RegComplement::surface_exterior(double *x, double cutoff)
{
  assert(false);
  return 0;
}
