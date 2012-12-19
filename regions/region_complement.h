/*
 * region_complement.h
 *
 *  Created on: June 26, 2012
 *      Author: kirill lykov
 */

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
