//  (C) Copyright Kirill Lykov 2013.
//
// Distributed under the GNU Software License (See accompanying file LICENSE)

#include "stdlib.h"
#include "string.h"
#include "region_difference.h"
#include "domain.h"
#include "error.h"
#include <assert.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

RegDifference::RegDifference(LAMMPS *lmp, int narg, char **arg) :
   RegIntersect(lmp, narg, arg)
{
}

/* ---------------------------------------------------------------------- */

RegDifference::~RegDifference()
{
}


/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is match() with all sub-regions
   else inside = 0
------------------------------------------------------------------------- */

int RegDifference::inside(double x, double y, double z)
{
  int ilist;
  Region **regions = domain->regions;
  Region* leftSet = regions[list[0]];
  for (ilist = 1; ilist < nregion; ilist++)
    if ( !(leftSet->match(x, y, z) && !regions[list[ilist]]->match(x, y, z)) )
      break;

  if (ilist == nregion) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   compute contacts with interior of intersection of sub-regions
   (1) compute contacts in each sub-region
   (2) only keep a contact if surface point is match() to all other regions
------------------------------------------------------------------------- */

int RegDifference::surface_interior(double *x, double cutoff)
{
  // not implemented
  assert(false);
  return 0;
}

/* ----------------------------------------------------------------------
   compute contacts with interior of intersection of sub-regions
   (1) flip interior/exterior flag of each sub-region
   (2) compute contacts in each sub-region
   (3) only keep a contact if surface point is not match() to all other regions
   (4) flip interior/exterior flags back to original settings
   this is effectively same algorithm as surface_interior() for RegUnion
------------------------------------------------------------------------- */

int RegDifference::surface_exterior(double *x, double cutoff)
{
  // not implemented
  assert(false);
  return 0;
}
