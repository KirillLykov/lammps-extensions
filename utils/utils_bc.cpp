//  (C) Copyright Kirill Lykov 2013.
//
// Distributed under the GNU Software License (See accompanying file LICENSE)

#include "gather_containers.h"
#include <algorithm>
#include <assert.h>

// TODO refactoring needed - introduce template method
void gatherUnionOfContainers(const std::vector<double>& localVector,
    MPI_Comm communicator, int root, std::vector<double>& globalUnionVector, size_t resultShift)
{
  int me = -1;
  MPI_Comm_rank(communicator, &me);

  int participants = 0;
  MPI_Comm_size(communicator, &participants); // if it is called on the core which doesn't belong to this comm, application crashes

  std::vector<int> counts(participants, 0);
  int localCount = localVector.size();
  MPI_Gather(&localCount, 1, MPI_INT, &counts[0], 1, MPI_INT, root, communicator);

  int globalCount = 0;
  std::vector<int> displs(participants, 0);
  if (me == root) {
    for (int i = 0; i < participants; ++i) {
      displs[i] += globalCount;
      globalCount += counts[i];
    }

    globalUnionVector.resize(globalCount + resultShift, 0.0);
  }
  std::vector<double>& localWithoutConst = const_cast< std::vector<double>& >(localVector);
  MPI_Gatherv(&localWithoutConst[0], localCount, MPI_DOUBLE, &globalUnionVector[resultShift], &counts[0], &displs[0], MPI_DOUBLE, root, communicator);
}

template<>
void allGatherUnionOfContainers(const VectorInt& localVector, MPI_Comm communicator, VectorInt& globalUnionVector)
{
  int participants = 0;
  MPI_Comm_size(communicator, &participants); // if it is called on the core which doesn't belong to this comm, application crashes

  std::vector<int> counts(participants, 0);
  int localCount = localVector.size();
  MPI_Allgather(&localCount, 1, MPI_INT, &counts[0], 1, MPI_INT, communicator);

  int globalCount = 0;
  std::vector<int> displs(participants, 0);
  for (int i = 0; i < participants; ++i) {
    displs[i] += globalCount;
    globalCount += counts[i];
  }

  globalUnionVector.resize(globalCount, 0);
  VectorInt& localWithoutConst = const_cast<VectorInt&>(localVector);
  MPI_Allgatherv(&localWithoutConst[0], localCount, MPI_INT, &globalUnionVector[0], &counts[0], &displs[0], MPI_INT, communicator);
}

template<>
void allGatherUnionOfContainers(const SetInt& localSet, MPI_Comm communicator, SetInt& globalUnionSet)
{
  VectorInt localVector;
  std::copy(localSet.begin(), localSet.end(), std::back_inserter(localVector));

  VectorInt globalUnionVector;
  allGatherUnionOfContainers(localVector, communicator, globalUnionVector);

  globalUnionSet.insert(globalUnionVector.begin(), globalUnionVector.end());
}

template<>
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

