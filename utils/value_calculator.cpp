//  (C) Copyright Kirill Lykov 2013.
//
// Distributed under the GNU Software License (See accompanying file LICENSE)

#include "ValueCalculator.h"
#include "error.h"
#include "force.h"
#include "math_extra.h"
#include "atom.h"
#include "update.h"
#include "comm.h"
#include <assert.h>
#include "region.h"
#include "group.h"
#include "domain.h"

using namespace LAMMPS_NS;

ValueCalculator::ValueCalculator(class LAMMPS * lmp, int groupbit, int nevery)
: Pointers(lmp), m_groupbit(groupbit), m_countOfMesurments(0), m_region(0),
  m_isActive(false), m_atomsCount(0), m_comm(MPI_COMM_NULL),
  m_root(0), m_firstTimeStep(0), m_nevery(nevery)
{
}

void ValueCalculator::setup()
{
  //if there are several runs then it is not 0
  m_firstTimeStep = update->ntimestep;

  // if subdomain doesn't intersect bounding box for the region
  // there is no need to run end_of_step. I check it using monte-carlo approach
  // due to it seems to be natural in lammps

  int igroup = group->find("all");
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

void ValueCalculator::run()
{
  if (m_isActive) {
    int countLocal = 0;
    double** x = atom->x;

    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double localAvgVel[3];
    memset(localAvgVel, 0, 3 * sizeof(localAvgVel[0]));

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & m_groupbit) {
        if (m_region->match(x[i][0], x[i][1], x[i][2])) {
          calculateLocalValue(i);
          ++countLocal;
        }
      }
    }

    int globalCount = 0;
    assert(m_comm != MPI_COMM_NULL);
    MPI_Allreduce(&countLocal, &globalCount, 1, MPI_INT, MPI_SUM, m_comm);

    calculateGlobalValue(globalCount);

    if (isRoot()) {
      m_atomsCount += globalCount;
      if (update->ntimestep - m_firstTimeStep == m_nevery * m_countOfMesurments) {
        //time-averaging
        m_atomsCount /= m_countOfMesurments;
        writeValue();
      }
    }
  }
  globalAfterRun();
}

bool ValueCalculator::setRegion(const std::string& regionName)
{
  int iregion = domain->find_region(const_cast<char*>(regionName.c_str()));
  m_region = domain->regions[iregion];

  return iregion != -1;
}

bool ValueCalculator::isRoot() const
{
  assert(m_isActive);

  int rankInGroup;
  MPI_Comm_rank(m_comm, &rankInGroup);

  return rankInGroup == m_root;
}

MPI_Comm ValueCalculator::createCommunicator()
{
  //in sake of communication performance create communicator for active procs
  MPI_Comm communicator = MPI_COMM_NULL;
  MPI_Comm_split(world, m_isActive ? 0 : MPI_UNDEFINED , comm->me, &communicator);
  return communicator;
}

// DensityCalculator

void DensityCalculator::globalAfterRun()
{
  double newDensity = 0.0;

  if (isActive() && isRoot() && comm->me == 0) {
    newDensity = m_density;
  } else {
    // Send density from local root to global root
    if (isActive() && isRoot()) {
      MPI_Send(&m_density, 1, MPI_DOUBLE, 0, 0, world);
    }

    if (comm->me == 0) {
       MPI_Status statusProbe;

       MPI_Probe(MPI_ANY_SOURCE, 0, world, &statusProbe);
       MPI_Recv(&newDensity, 1, MPI_DOUBLE, statusProbe.MPI_SOURCE, 0, world, MPI_STATUS_IGNORE);
    }
  }

  // Broadcast to all other processors. It is more efficient that using Send/Recv directly
  MPI_Bcast(&newDensity, 1, MPI_DOUBLE, 0, world);
  if (isActive()) { // to check the correctness of the communication
    assert(m_density == newDensity);  
  }
  m_density = newDensity;
}

// VelocityCalculator

VelocityCalculator::VelocityCalculator(class LAMMPS * lmp, int groupbit, int nevery)
: ValueCalculator(lmp, groupbit, nevery)
{
  memset(m_localAvgVel, 0, 3 * sizeof(m_localAvgVel[0]));
  memset(m_avgVel, 0, 3 * sizeof(m_avgVel[0]));
}

void VelocityCalculator::calculateLocalValue(int i)
{
  MathExtra::add3(m_localAvgVel, atom->v[i], m_localAvgVel);
}

void VelocityCalculator::calculateGlobalValue(int globalCount)
{
  assert(m_countOfMesurments != 0); // TODO it looks like a bug that m_countOfMesurments is always 0

  double globalAvgVel[3];
  memset(globalAvgVel, 0, 3 * sizeof(globalAvgVel[0]));
  MPI_Reduce(&m_localAvgVel, &globalAvgVel, 3, MPI_DOUBLE, MPI_SUM, m_root, m_comm);
  MathExtra::scale3(1.0 / globalCount, globalAvgVel);
  MathExtra::add3(m_avgVel, globalAvgVel, m_avgVel);
}

void VelocityCalculator::writeValue()
{
  MathExtra::scale3(1.0 / static_cast<double>(m_countOfMesurments), m_avgVel);
}
