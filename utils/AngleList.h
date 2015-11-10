//  (C) Copyright Kirill Lykov 2013.
//
// Distributed under the GNU Software License (See accompanying file LICENSE)

#ifndef ANGLELIST_H_
#define ANGLELIST_H_

namespace LAMMPS_NS {
/**
 * @class
 *  Used for iterating though angels in a std-like method.
 *  It might be used with standard lammps angles and with your own angle list, if needed.
 *  Example:
 *    AngleList angleList(neighbor->nanglelist, neighbor->anglelist, myAngles->m_nanglelist, myAngles->m_anglelist);
 *    for (AngleList::Iterator it = angleList.begin(); it != angleList.end(); ++it) { }
 */

class AngleList {

  int m_lammpsAnglesCount;
  int** m_lammpsAngles;
  int m_addAnglesCount;
  int** m_additionalAngles;
public:

  struct Triangle {
    int i1, i2, i3;

    Triangle() : i1(0), i2(0), i3(0)
    {
    }

    Triangle(int inI1, int inI2, int inI3) : i1(inI1), i2(inI2), i3(inI3)
    {
    }

    void set(int inI1, int inI2, int inI3)
    {
      i1 = inI1;
      i2 = inI2;
      i3 = inI3;
    }

    bool operator== (const Triangle& right) const
    {
      return (i1 == right.i1) && (i2 == right.i2) && (i3 == right.i3);
    }

    bool operator!= (const Triangle& right) const
    {
      return !(*this == right);
    }

    bool operator< (const Triangle& right) const
    {
      return (i1 < right.i1) ||
          ((i1 == right.i1) && (i2 < right.i2)) ||
          ((i1 == right.i1) && (i2 == right.i2) && (i3 < right.i3));
    }

    bool operator> (const Triangle& right) const
    {
      return (i1 > right.i1) ||
          ((i1 == right.i1) && (i2 > right.i2)) ||
          ((i1 == right.i1) && (i2 == right.i2) && (i3 > right.i3));
    }
  };

  class Iterator {
    const AngleList& m_angleList;
    int m_currentTriangle;
  public:
    Iterator(const AngleList& angleList, int currentTriangle)
    : m_angleList(angleList), m_currentTriangle(currentTriangle)
    {
    }

    Iterator& operator++()
    {
      ++m_currentTriangle;
      return *this;
    }

    bool operator== (const Iterator& right) const
    {
      return (&m_angleList == &right.m_angleList) &&
             (m_currentTriangle == right.m_currentTriangle);
    }

    bool operator!= (const Iterator& right) const
    {
      return !(*this == right);
    }

    void getTriangle(int& i1, int& i2, int& i3) const;

    Triangle getTriangle() const;

    void getTriangle(int* indexes) const;
  };

  AngleList(int lammpsAnglesCount, int** lammpsAngles, int addAnglesCount = 0, int** additionalAngles = 0)
  : m_lammpsAnglesCount(lammpsAnglesCount), m_lammpsAngles(lammpsAngles),
    m_addAnglesCount(addAnglesCount), m_additionalAngles(additionalAngles)
  {
  }

  Iterator begin() const
  {
    return Iterator(*this, 0);
  }

  Iterator end() const
  {
    return Iterator(*this, m_lammpsAnglesCount + m_addAnglesCount);
  }

private:
  AngleList(const AngleList&);
  AngleList& operator= (const AngleList&);
};

}

#endif /* ANGLELIST_H_ */
