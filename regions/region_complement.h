//  (C) Copyright Kirill Lykov 2013.
//
// Distributed under the GNU Software License (See accompanying file LICENSE)

#ifdef REGION_CLASS

RegionStyle(complement,RegComplement)

#else

#ifndef LMP_REGION_COMPLEMENT_H
#define LMP_REGION_COMPLEMENT_H

#include "region.h"

namespace LAMMPS_NS {

class RegComplement : public Region {
 public:
  RegComplement(class LAMMPS *, int, char **);
  ~RegComplement();
  int dynamic_check();

  int inside(double, double, double);

  // not implemented:
  int surface_interior(double *, double);
  int surface_exterior(double *, double);
 protected:
  int iregion;
};

}

#endif
#endif
