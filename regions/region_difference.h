//  (C) Copyright Kirill Lykov 2013.
//
// Distributed under the GNU Software License (See accompanying file LICENSE)

#ifdef REGION_CLASS

RegionStyle(difference,RegDifference)

#else

#ifndef LMP_REGION_DIFERENCE_H
#define LMP_REGION_DIFERENCE_H

#include "region_intersect.h"

namespace LAMMPS_NS {

class RegDifference : public RegIntersect {
 public:
  RegDifference(class LAMMPS *, int, char **);
  ~RegDifference();

  int inside(double, double, double);

  // not implemented:
  int surface_interior(double *, double);
  int surface_exterior(double *, double);
};

}

#endif
#endif
