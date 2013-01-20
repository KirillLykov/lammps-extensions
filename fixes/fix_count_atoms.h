//  (C) Copyright Kirill Lykov 2013.
//
// Distributed under the GNU Software License (See accompanying file LICENSE)

#ifdef FIX_CLASS

FixStyle(count/atoms,FixCountAtoms)

#else

#ifndef LMP_FIX_COUNT_ATOMS_H
#define LMP_FIX_COUNT_ATOMS_H

#include "fix.h"
#include <string>

namespace LAMMPS_NS {

class FixCountAtoms : public Fix {
  int m_countOfMesurments;
  class Region* m_region;
  bool m_isActive; //if the subdomain of the current proc doesn't contain region, don't do any computations
  std::string m_fileName;
  int m_atomsCount; // atoms for m_countOfMesurments
  MPI_Comm m_comm;
  int m_root;
  double m_velDir[3];
  double m_avgVel[3];
 public:
  FixCountAtoms(class LAMMPS *, int, char **);
  virtual ~FixCountAtoms();
  int setmask();
  void setup(int);
  void end_of_step();
 private:
  void writeResult();
  MPI_Comm createCommunicator();
};

}

#endif
#endif
