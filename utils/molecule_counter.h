/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Kirill Lykov
------------------------------------------------------------------------- */

#ifndef MOLECULES_COUNTER_H_
#define MOLECULES_COUNTER_H_

#include "pointers.h"

namespace LAMMPS_NS {
/**
 * @class
 *  Counts the number of molecules and creates a mapping between global and local molIndexes.
 *  It was extracted from one of the LAMMPS classes to be used in other places
 */
class MoleculeCounter : protected Pointers
{
  int* molmap; // convert molecule ID to local index
  int m_idlo, m_idhi, m_molnum; // saved from the last run molecules indexes range and number of molecules
public:

  MoleculeCounter(LAMMPS *lmp);
  ~MoleculeCounter();

  void run(int groupbit);
  int getIdlo() const { return m_idlo; }
  int getIdhi() const { return m_idhi; }
  int getMolNum() const { return m_molnum; }

  int getMolIDbyAtom(int atomIndex) const;

  /**
   * @param globalMolIndex
   *  lammps molecule index
   * @return
   *  molecule id in local index, it is in range [0, molnum]
   */
  int getMolID(int globalMolIndex) const;

private:
  int molecules_in_group(int groupbit, int &idlo, int &idhi);
};
}
#endif

