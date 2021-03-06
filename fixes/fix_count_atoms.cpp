//  (C) Copyright Kirill Lykov 2013.
//
// Distributed under the GNU Software License (See accompanying file LICENSE)

#include "fix_count_atoms.h"
#include <fstream>
#include <vector>
#include <iostream>
#include <assert.h>
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "region.h"
#include "error.h"
#include "comm.h"
#include "update.h"
#include "group.h"
#include "math_extra.h"

using namespace LAMMPS_NS;

FixCountAtoms::FixCountAtoms(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), m_countOfMesurments(0), m_region(0),
  m_isActive(false), m_atomsCount(0), m_comm(MPI_COMM_NULL),
  m_root(0), m_firstTimeStep(0)
{
  if (narg < 6) error->all(FLERR,"Illegal fix wall/bb command");

  int iregion = domain->find_region(arg[3]);
  if (iregion == -1) error->all(FLERR,"Region ID does not exist");

  m_region = domain->regions[iregion];

  nevery = force->inumeric(arg[4]);
  m_countOfMesurments = force->inumeric(arg[5]);

  if (narg >= 7)
    m_fileName = std::string(arg[6]);

  memset(m_velDir, 0, 3 * sizeof(m_avgVel[0]));
  memset(m_avgVel, 0, 3 * sizeof(m_avgVel[0]));

  if (narg == 10) {
    m_velDir[0] = force->numeric(arg[7]);
    m_velDir[1] = force->numeric(arg[8]);
    m_velDir[2] = force->numeric(arg[9]);
  }
}

FixCountAtoms::~FixCountAtoms()
{
}

int FixCountAtoms::setmask()
{
  int mask = 0;
  mask |= FixConst::END_OF_STEP;
  return mask;
}

void FixCountAtoms::setup(int)
{
  //if there are several runs then it is not 0
  m_firstTimeStep = update->ntimestep;

  // if subdomain doesn't intersect bounding box for the region
  // there is no need to run end_of_step. I check it using monte-carlo approach
  // due to it seems to be natural in lammps

  igroup = group->find("all");
  if (igroup == -1) error->all(FLERR,"Could not find fix group ID \"all\"");
  int all_groupbit = group->bitmask[igroup];

  double** x = atom->x;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & all_groupbit) {
      if (m_region->match(x[i][0], x[i][1], x[i][2])) {
        m_isActive = true;
        break;
      }
    }
  }
  m_comm = createCommunicator();
  assert(m_isActive ? (m_comm != MPI_COMM_NULL) : (m_comm == MPI_COMM_NULL));
}

void FixCountAtoms::end_of_step()
{
  if (!m_isActive)
    return;

  int countLocal = 0;
  double** x = atom->x;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double localAvgVel[3];
  memset(localAvgVel, 0, 3 * sizeof(localAvgVel[0]));

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (m_region->match(x[i][0], x[i][1], x[i][2])) {
        MathExtra::add3(localAvgVel, atom->v[i], localAvgVel);
        ++countLocal;
      }
    }
  }

  int globalCount = 0;
  assert(m_comm != MPI_COMM_NULL);
  MPI_Reduce(&countLocal, &globalCount, 1, MPI_INT, MPI_SUM, m_root, m_comm);

  double globalAvgVel[3];
  memset(globalAvgVel, 0, 3 * sizeof(globalAvgVel[0]));
  MPI_Reduce(&localAvgVel, &globalAvgVel, 3, MPI_DOUBLE, MPI_SUM, m_root, m_comm);

  MathExtra::scale3(1.0 / globalCount, globalAvgVel);
  MathExtra::add3(m_avgVel, globalAvgVel, m_avgVel);

  int rankInGroup;
  MPI_Comm_rank(m_comm, &rankInGroup);

  if (rankInGroup == m_root) {
    m_atomsCount += globalCount;
    if (update->ntimestep - m_firstTimeStep == nevery * m_countOfMesurments) {
      //time-averaging
      m_atomsCount /= m_countOfMesurments;
      MathExtra::scale3(1.0 / static_cast<double>(m_countOfMesurments), m_avgVel);
      writeResult();
    }
  }
}

void FixCountAtoms::writeResult()
{
  std::fstream file(m_fileName.c_str(), std::fstream::app|std::fstream::out);
  double velInDirection = MathExtra::dot3(m_velDir, m_avgVel);
  if (file.is_open()) {
    file << update->ntimestep << " " << m_atomsCount
      << " " << m_avgVel[0] << " " << m_avgVel[1] << " " << m_avgVel[2] << " " << velInDirection << std::endl;
    file.flush();
  } else {
    std::cout << update->ntimestep << " " << m_atomsCount
        << " " << m_avgVel[0] << " " << m_avgVel[1] << " " << m_avgVel[2] << " " << velInDirection << std::endl;
  }
  file.close();
}

MPI_Comm FixCountAtoms::createCommunicator()
{
  //in sake of communication performance create communicator for active procs
  MPI_Comm communicator = MPI_COMM_NULL;
  MPI_Comm_split(world, m_isActive ? 0 : MPI_UNDEFINED , comm->me, &communicator);
  return communicator;
}
