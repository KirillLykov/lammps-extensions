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

#include "MoleculeCounter.h"
#include "memory.h"
#include "error.h"
#include "atom.h"
#include "comm.h"
#include <algorithm>
#include <assert.h>

using namespace LAMMPS_NS;

MoleculeCounter::MoleculeCounter(LAMMPS *lmp)
: Pointers(lmp), molmap(0), m_idlo(0), m_idhi(0), m_molnum(0)
{}

MoleculeCounter::~MoleculeCounter()
{
  memory->destroy(molmap);
}

void MoleculeCounter::run(int groupbit)
{
  m_molnum = molecules_in_group(groupbit, m_idlo, m_idhi);
}

int MoleculeCounter::getMolIDbyAtom(int atomIndex) const
{
  int globalMolIndex = atom->molecule[atomIndex];
  return getMolID(globalMolIndex);
}

int MoleculeCounter::getMolID(int globalMolIndex) const
{
  assert(globalMolIndex >= m_idlo && globalMolIndex <= m_idhi);
  int localMolIndex = molmap[globalMolIndex - m_idlo];
  return localMolIndex;
}

int MoleculeCounter::molecules_in_group(int groupbit, int &idlo, int &idhi)
{
  memory->destroy(molmap);
  molmap = 0;

  // find lo/hi molecule ID for any atom in group
  // warn if atom in group has ID = 0

  int *molecule = atom->molecule;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int lo = INT_MAX;
  int hi = INT_MIN;
  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      // skip molecules with ID=0 because they are free particles
      if (molecule[i] == 0)
        flag = 1;
      else
        lo = std::min(lo, molecule[i]);
      hi = std::max(hi, molecule[i]);
    }

  int flagall;
  MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);
  if (flagall && comm->me == 0)
    error->warning(FLERR,"Atom with molecule ID = 0 included in "
                   "compute molecule group");

  MPI_Allreduce(&lo, &idlo, 1, MPI_INT, MPI_MIN, world);
  MPI_Allreduce(&hi, &idhi, 1, MPI_INT, MPI_MAX, world);
  if (idlo == INT_MAX) return 0;

  // molmap = vector of length nlen
  // set to 1 for IDs that appear in group across all procs, else 0

  int nlen = idhi - idlo + 1;
  memory->create(molmap, nlen, "compute:molmap");
  for (int i = 0; i < nlen; i++) molmap[i] = 0;

  for (int i = 0; i < nlocal; i++)
    if ((mask[i] & groupbit) && molecule[i] != 0) {
      assert((molecule[i]-idlo) >= 0 && (molecule[i]-idlo) < nlen);
      molmap[molecule[i]-idlo] = 1;
    }

  int *molmapall;
  memory->create(molmapall, nlen, "compute:molmapall");
  MPI_Allreduce(molmap, molmapall, nlen, MPI_INT, MPI_MAX, world);

  // nmolecules = # of non-zero IDs in molmap
  // molmap[i] = index of molecule, skipping molecules not in group with -1

  int nmolecules = 0;
  for (int i = 0; i < nlen; i++)
    if (molmapall[i]) molmap[i] = nmolecules++;
    else molmap[i] = -1;
  memory->destroy(molmapall);

  // warn if any molecule has some atoms in group and some not in group

  flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) continue;
    if (molecule[i] < idlo || molecule[i] > idhi) continue;
    if (molmap[molecule[i]-idlo] >= 0) flag = 1;
  }

  MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);
  if (flagall && comm->me == 0)
    error->warning(FLERR,"One or more compute molecules has atoms not in group");

  // if molmap simply stores 1 to Nmolecules, then free it
  //if (idlo == 1 && idhi == nmolecules && nlen == nmolecules) {
  //  memory->destroy(molmap);
  //  molmap = NULL;
  //}
  return nmolecules;
}

