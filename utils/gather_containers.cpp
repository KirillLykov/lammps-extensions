//  (C) Copyright Kirill Lykov 2013.
//
// Distributed under the GNU Software License (See accompanying file LICENSE)

#include "gather_containers.h"
#include <algorithm>
#include <assert.h>

void allGatherUnionOfContainers(const SetInt& localSet, MPI_Comm communicator, SetInt& globalUnionSet)
{
  VectorInt localVector;
  std::copy(localSet.begin(), localSet.end(), std::back_inserter(localVector));

  VectorInt globalUnionVector;
  allGatherUnionOfContainers(localVector, communicator, globalUnionVector);

  globalUnionSet.insert(globalUnionVector.begin(), globalUnionVector.end());
}

void allGatherUnionOfContainers(const MapIntInt& localMap, MPI_Comm communicator, MapIntInt& globalUnionMap)
{
  std::vector<int> outLocal;

  size_t sz = 2 * localMap.size();
  outLocal.reserve(sz);
  for (CIterMapInt it = localMap.begin(); it != localMap.end(); ++it) {
    outLocal.push_back(it->first);
    outLocal.push_back(it->second);
  }

  std::vector<int> outGlobal;
  allGatherUnionOfContainers(outLocal, communicator, outGlobal);

  for (size_t i = 0; i < outGlobal.size(); i += 2) {
    std::pair<int, int> p(std::make_pair(outGlobal[i], outGlobal[i + 1]));
    globalUnionMap.insert(p);
  }
}

