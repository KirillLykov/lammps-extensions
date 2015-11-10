//  (C) Copyright Kirill Lykov 2014.
//
// Distributed under the GNU Software License (See accompanying file LICENSE)

#include "AngleList.h"
#include <assert.h>
#include <memory.h>

namespace LAMMPS_NS {

void AngleList::Iterator::getTriangle(int& i1, int& i2, int& i3) const
{
  if (m_currentTriangle < m_angleList.m_lammpsAnglesCount)
  {
    i1 = m_angleList.m_lammpsAngles[m_currentTriangle][0];
    i2 = m_angleList.m_lammpsAngles[m_currentTriangle][1];
    i3 = m_angleList.m_lammpsAngles[m_currentTriangle][2];
  } else {
    int addInd = m_currentTriangle - m_angleList.m_lammpsAnglesCount;
    assert(addInd < m_angleList.m_addAnglesCount);
    i1 = m_angleList.m_additionalAngles[addInd][0];
    i2 = m_angleList.m_additionalAngles[addInd][1];
    i3 = m_angleList.m_additionalAngles[addInd][2];
  }
}

AngleList::Triangle AngleList::Iterator::getTriangle() const
{
  Triangle tr;
  if (m_currentTriangle < m_angleList.m_lammpsAnglesCount)
  {
    tr.i1 = m_angleList.m_lammpsAngles[m_currentTriangle][0];
    tr.i2 = m_angleList.m_lammpsAngles[m_currentTriangle][1];
    tr.i3 = m_angleList.m_lammpsAngles[m_currentTriangle][2];
  } else {
    int addInd = m_currentTriangle - m_angleList.m_lammpsAnglesCount;
    assert(addInd < m_angleList.m_addAnglesCount);
    tr.i1 = m_angleList.m_additionalAngles[addInd][0];
    tr.i2 = m_angleList.m_additionalAngles[addInd][1];
    tr.i3 = m_angleList.m_additionalAngles[addInd][2];
  }
  return tr;
}

void AngleList::Iterator::getTriangle(int* indexes) const
{
  if (m_currentTriangle < m_angleList.m_lammpsAnglesCount)
  {
    memcpy( indexes, m_angleList.m_lammpsAngles[m_currentTriangle], 3 * sizeof(indexes[0]) );
  } else {
    int addInd = m_currentTriangle - m_angleList.m_lammpsAnglesCount;
    assert(addInd < m_angleList.m_addAnglesCount);
    memcpy( indexes, m_angleList.m_additionalAngles[addInd], 3 * sizeof(indexes[0]) );
  }
}

}
