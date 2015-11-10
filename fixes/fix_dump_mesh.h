//  (C) Copyright Kirill Lykov 2015.
//
// Distributed under the GNU Software License (See accompanying file LICENSE)

#ifdef FIX_CLASS

FixStyle(dump/mesh,FixDumpMesh)

#else

#ifndef LMP_FIX_DUMP_MESH
#define LMP_FIX_DUMP_MESH

#include "fix.h"
#include <string>
#include <vector>
#include <unordered_map>

namespace LAMMPS_NS {

/**
* @class
*   Dumps structure in OBJ geometry format, uses angles for triangulation
*/
class FixDumpMesh : public Fix
{
  int m_nglobalParticles;
  const double m_maxArea;
  std::string m_fileNameTemplate;
  std::unordered_map<int, int> m_tags2VertInd; // mapping from tags to indicies of vertices used for obj
  std::vector<int> m_triangulation; // triangles indices, linearized, tags
public:
  FixDumpMesh(class LAMMPS *, int, char **);
  ~FixDumpMesh();
  int setmask();
  void setup(int);
  void end_of_step();
private:
  void writeObj(const std::string& fileName, std::vector<float>& positions);
};

}

#endif
#endif
