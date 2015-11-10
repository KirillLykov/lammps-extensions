//  (C) Copyright Kirill Lykov 2013.
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
#include <fstream>

namespace LAMMPS_NS {

class FixDumpMesh : public Fix
{
  int m_nglobalParticles;
  std::string m_fileNameTemplate;
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
