//  (C) Copyright Kirill Lykov 2013.
//
// Distributed under the GNU Software License (See accompanying file LICENSE)

#include "fix_dump_mesh.h"
#include "update.h"
#include "input.h"
#include "modify.h"
#include "variable.h"
#include "error.h"
#include <string>
#include "../utils/gather_containers.h"
#include "atom.h"
#include "neighbor.h"
#include "comm.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <ios>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixDumpMesh::FixDumpMesh(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), m_nglobalParticles(0)
{
  if (narg < 5) error->all(FLERR,"Illegal fix print command");

  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix dump mesh command: nevery must be positive integer");

  m_fileNameTemplate = std::string(arg[4]);
}

/* ---------------------------------------------------------------------- */

FixDumpMesh::~FixDumpMesh()
{
}

/* ---------------------------------------------------------------------- */

int FixDumpMesh::setmask()
{
  int mask = 0;
  mask |= FixConst::END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDumpMesh::setup(int)
{
  // count n particles of interest
  int nlocalParticles = 0;
  for (int i = 0; i < atom->nlocal; ++i) {
    if (atom->mask[i] & groupbit) {
      ++nlocalParticles;
    }
  }
  MPI_Allreduce(&nlocalParticles, &m_nglobalParticles, 1, MPI_INT, MPI_SUM, world);

  // collect all relevant triangles, assumed that topology is const
  std::vector<int> localTriangulation;
  for (int i = 0; i < atom->nlocal; ++i) {
    if (atom->mask[i] & groupbit) {
      int num_angle = atom->num_angle[i];
      for (int m = 0; m < num_angle; ++m) {
        localTriangulation.push_back(atom->angle_atom1[i][m]);
        localTriangulation.push_back(atom->angle_atom2[i][m]);
        localTriangulation.push_back(atom->angle_atom3[i][m]);
      }
    }
  }
  gatherUnionOfContainers(localTriangulation, world, 0, m_triangulation);
}

/* ---------------------------------------------------------------------- */

void FixDumpMesh::end_of_step()
{
  std::vector<float> localAtoms; // linearized array [tag, x, y ,z]
  localAtoms.reserve(m_nglobalParticles);
  for (int i = 0; i < atom->nlocal; ++i) {
    if (atom->mask[i] & groupbit) {
      localAtoms.push_back( static_cast<float>(atom->tag[i]) );
      localAtoms.push_back( static_cast<float>(atom->x[i][0]) );
      localAtoms.push_back( static_cast<float>(atom->x[i][1]) );
      localAtoms.push_back( static_cast<float>(atom->x[i][2]) );
    }
  }
  std::vector<float> globalUnionVector;
  gatherUnionOfContainers(localAtoms, world, 0, globalUnionVector);

  std::string fileName = m_fileNameTemplate + "." + std::to_string(update->ntimestep/nevery) + ".obj";
  writeObj(fileName, globalUnionVector);
}

/* ---------------------------------------------------------------------- */

struct Point {
  float tag, x, y, z;
};
bool operator< (Point i, Point j) { return (i.tag < j.tag); }

void FixDumpMesh::writeObj(const std::string& fileName, std::vector<float>& positions)
{
  try
  {
    if (comm->me == 0) {
      std::ofstream file(fileName.c_str());
      if (!file.is_open()) {
        error->all(FLERR, "could not open output file");
      } else {
        std::string objName = "Mesh";
        file << "# Generate by FixDumpMesh" << std::endl << "o " + objName << std::endl;
      }

      if (sizeof(Point) != 4 * sizeof(float)) {
        error->one(FLERR, "Internal error: point structure has unexpected padding");
      }

      // write vertices
      Point* pPoints = reinterpret_cast<Point*>( &positions[0] );
      std::sort(pPoints, pPoints + positions.size() / 4);
      for (size_t i = 0; i < positions.size() / 4; ++i) {
        file << "v " << std::fixed << pPoints[i].x << " " << pPoints[i].y << " " << pPoints[i].z << std::endl;
      }

      // write triangles
      for (std::vector<int>::const_iterator it = m_triangulation.begin(); it != m_triangulation.end(); it += 3) {
        file << "f " << *(it) << " " << *(it + 1) << " " << *(it + 2) << std::endl;
      }
    }
  }
  catch(...)
  {
    error->all(FLERR, "Internal error");
  }
}
