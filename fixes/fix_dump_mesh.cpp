//  (C) Copyright Kirill Lykov 2015.
//
// Distributed under the GNU Software License (See accompanying file LICENSE)

#include "fix_dump_mesh.h"
#include "update.h"
#include "input.h"
#include "modify.h"
#include "variable.h"
#include "error.h"
#include "math_extra.h"
#include "rbc_utils.h"
#include "atom.h"
#include "neighbor.h"
#include "comm.h"
#include <fstream>
#include <iostream>
#include <ios>

using namespace LAMMPS_NS;

namespace {
  // aux structure used for easy sorting of linearized array of point
  struct Point {
    float tag, x, y, z;
  };
  bool operator< (Point i, Point j) { return (i.tag < j.tag); }


  // helper for unordered_map object put/get
  void put(int tag, int vertInd, std::unordered_map<int, int>& map)
  {
    map.insert(std::make_pair(tag, vertInd));
  }

  int get(const std::unordered_map<int, int>& map, int tag)
  {
    std::unordered_map<int, int>::const_iterator got = map.find(tag);
    if (got == map.end())
      throw "tags mapping is invalid";
    return got->second;
  }
}

/* ---------------------------------------------------------------------- */

FixDumpMesh::FixDumpMesh(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), m_nglobalParticles(0), m_maxArea(2.0)
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
  std::vector<int> localVertInd2Tag;
  for (int i = 0; i < atom->nlocal; ++i) {
    if (atom->mask[i] & groupbit) {
      localVertInd2Tag.push_back( atom->tag[i] );
    }
  }
  int nlocalParticles = localVertInd2Tag.size();
  MPI_Allreduce(&nlocalParticles, &m_nglobalParticles, 1, MPI_INT, MPI_SUM, world);

  std::vector<int> globalVertInd2Tag;
  allGatherUnionOfContainers(localVertInd2Tag, world, globalVertInd2Tag);
  std::sort(globalVertInd2Tag.begin(), globalVertInd2Tag.end());
  for (int i = 0; i < (int)globalVertInd2Tag.size(); ++i) {
    put(globalVertInd2Tag[i], i, m_tags2VertInd);
  }

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
        // Skip triangles which intersect domain borders
        int vi[] = {get(m_tags2VertInd, *(it)), get(m_tags2VertInd, *(it + 1)), get(m_tags2VertInd, *(it + 2))};
        float bma[3] = {pPoints[vi[0]].x - pPoints[vi[1]].x, pPoints[vi[0]].y - pPoints[vi[1]].y, pPoints[vi[0]].z - pPoints[vi[1]].z};
        float cma[3] = {pPoints[vi[1]].x - pPoints[vi[2]].x, pPoints[vi[1]].y - pPoints[vi[2]].y, pPoints[vi[1]].z - pPoints[vi[2]].z};
        float area = 0.5 * ((bma[1]*cma[2] - bma[2]*cma[1])*(bma[1]*cma[2] - bma[2]*cma[1]) +
                            (bma[2]*cma[0] - bma[0]*cma[2])*(bma[2]*cma[0] - bma[0]*cma[2]) +
                            (cma[0]*bma[1] - bma[0]*cma[1])*(cma[0]*bma[1] - bma[0]*cma[1]));
        if (area < m_maxArea) {
          file << "f " << vi[0] + 1 << " " << vi[1] + 1 << " " << vi[2] + 1 << std::endl;
        }
      }
    }
  }
  catch(...)
  {
    error->all(FLERR, "Internal error");
  }
}
