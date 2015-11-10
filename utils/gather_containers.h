//  (C) Copyright Kirill Lykov 2013.
//
// Distributed under the GNU Software License (See accompanying file LICENSE)

#ifndef GATHER_CONTAINERS_H_
#define GATHER_CONTAINERS_H_

#include "mpi.h"
#include <vector>
#include <set>
#include <map>
#include <cstring>

typedef std::set<int> SetInt;
typedef typename std::set<int>::const_iterator CIterSetInt;

typedef std::map<int, int> MapIntInt;
typedef typename MapIntInt::const_iterator CIterMapInt;
typedef typename MapIntInt::iterator IterMapInt;

/**
 * All gather sets and maps
 */
void allGatherUnionOfContainers(const SetInt& localSet, MPI_Comm communicator, SetInt& globalUnionSet);
void allGatherUnionOfContainers(const MapIntInt& localMap, MPI_Comm communicator, MapIntInt& globalUnionVector);

/**
 * Returns union of vectors from all cores in communicator, do all necessary memory allocation
 * Usage example:
    VectorInt localVector(5 + comm->me*2, comm->me);
    VectorInt globalUnion;

    gatherUnionOfVectors(localVector, world, globalUnion);
    for(int i = 0; i < globalUnion.size(); ++i) {
      std::cout << globalUnion[i] << " ";
    }
    std::cout << std::endl;
    return;
 */
template<class T>
struct MPITrait;

template<>
struct MPITrait<int>
{
  MPI_Datatype dataType;
  MPITrait() : dataType(MPI_INT) {}
};

template<>
struct MPITrait<double>
{
  MPI_Datatype dataType;
  MPITrait() : dataType(MPI_DOUBLE) {}
};

template<>
struct MPITrait<float>
{
  MPI_Datatype dataType;
  MPITrait() : dataType(MPI_FLOAT) {}
};

template<class T>
static void allGatherUnionOfContainers(const std::vector<T>& localVector, MPI_Comm communicator, std::vector<T>& globalUnionVector)
{
  MPITrait<T> trait;
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
  std::vector<T>& localWithoutConst = const_cast<std::vector<T>&>(localVector);
  MPI_Allgatherv(&localWithoutConst[0], localCount, trait.dataType, &globalUnionVector[0], &counts[0], &displs[0], trait.dataType, communicator);
}

template<class T>
static void gatherUnionOfContainers(const std::vector<T>& localVector, MPI_Comm communicator, int root, std::vector<T>& globalUnionVector)
{
  MPITrait<T> trait;
  int participants = 0;
  MPI_Comm_size(communicator, &participants); // if it is called on the core which doesn't belong to this comm, application crashes

  std::vector<int> counts(participants, 0);
  int localCount = localVector.size();
  MPI_Gather(&localCount, 1, MPI_INT, &counts[0], 1, MPI_INT, root, communicator);

  int globalCount = 0;
  std::vector<int> displs(participants, 0);
  for (int i = 0; i < participants; ++i) {
    displs[i] += globalCount;
    globalCount += counts[i];
  }

  globalUnionVector.resize(globalCount, 0);
  std::vector<T>& localWithoutConst = const_cast<std::vector<T>&>(localVector);
  MPI_Gatherv(&localWithoutConst[0], localCount, trait.dataType, &globalUnionVector[0], &counts[0], &displs[0], trait.dataType, root, communicator);
}

#endif /* GATHER_CONTAINERS_H_ */
