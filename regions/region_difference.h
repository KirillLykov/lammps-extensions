/*
 * region_difference.h
 *
 *  Created on: June 26, 2012
 *      Author: kirill lykov
 */

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
