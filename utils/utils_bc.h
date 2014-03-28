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

typedef std::vector<int> VectorInt;
typedef typename std::vector<int>::const_iterator CIterVectorInt;

// TODO can be easily generalized for the case Container<T> if needed
template<class Container>
void allGatherUnionOfContainers(const Container& localVector, MPI_Comm communicator, Container& globalUnion);

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
template<>
void allGatherUnionOfContainers(const VectorInt& localVector, MPI_Comm communicator, VectorInt& globalUnionVector);

/**
 * The same but for sets
 */
template<>
void allGatherUnionOfContainers(const SetInt& localSet, MPI_Comm communicator, SetInt& globalUnionSet);

template<>
void allGatherUnionOfContainers(const MapIntInt& localMap, MPI_Comm communicator, MapIntInt& globalUnionVector);

#endif /* GATHER_CONTAINERS_H_ */
