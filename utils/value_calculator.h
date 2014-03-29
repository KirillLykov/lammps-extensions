//  (C) Copyright Kirill Lykov 2013.
//
// Distributed under the GNU Software License (See accompanying file LICENSE)

#ifndef VALUECALCULATOR_H_
#define VALUECALCULATOR_H_

#include <string>
#include <assert.h>
#include "pointers.h"
#include "region.h"

namespace LAMMPS_NS {

/**
 * @class
 *  A base class used for all calculations of statistical properties in a specified region.
 *  Communication is only between active cores, core is active if there are particles of the specified
 *  group which belong to this core and to the region.
 */
class ValueCalculator : protected Pointers
{
protected:
  int m_countOfMesurments;
  class Region* m_region;
  bool m_isActive; //if the subdomain of the current proc doesn't contain region, don't do any computations
  std::string m_fileName;
  int m_atomsCount; // atoms for m_countOfMesurments
  MPI_Comm m_comm;
  int m_root;
  double m_velDir[3];
  bigint m_firstTimeStep;
  int m_groupbit, m_nevery;

public:

  ValueCalculator(class LAMMPS * lmp, int groupbit, int nevery);

  virtual ~ValueCalculator() {}

  virtual void setup();

  virtual void run();

  virtual bool setRegion(const std::string& regionName);

  virtual bool isActive() const { return m_isActive; }

  virtual bool isRoot() const;

protected:
  virtual void calculateLocalValue(int i) = 0;
  virtual void calculateGlobalValue(int globalCount) = 0;
  virtual void writeValue() = 0;
  virtual void globalAfterRun() = 0;

  MPI_Comm createCommunicator();


  ValueCalculator(ValueCalculator&);
  ValueCalculator& operator=(const ValueCalculator&);
};

/**
 * @class
 * Used for the density calculation in the specified region
 * Example:
 *  MyFix::MyFix(LAMMPS *lmp, int narg, char **arg)
 *  : densCalc(lmp, groupbit, 0) {}
 *
 *  void MyFix::setup(int)
 *  {
 *    densCalc.setup();
 *  }
 *
 *   void MyFix::post_integrate()
 *   {
 *     densCalc.run();
 *     currentDensitySum += densCalc.getDensity();
 *   }
 */
class DensityCalculator : public ValueCalculator
{
  double m_volume, m_density;
public:
  DensityCalculator(class LAMMPS * lmp, int groupbit, int nevery)
  : ValueCalculator(lmp, groupbit, nevery), m_volume(0.0), m_density(0.0)
  {}

  void setup()
  {
    ValueCalculator::setup();
    m_volume = m_region->volume_interior();
  }

  void setup(double volume)
  {
    ValueCalculator::setup();
    m_volume = volume;
  }

  double getDensity() const
  {
    return m_density;
  }

protected:
  void calculateLocalValue(int i)
  {
  }

  void calculateGlobalValue(int globalCount)
  {
    assert(m_volume > 0.0);
    m_density = static_cast<double>(globalCount) / m_volume;
  }

  void globalAfterRun();

  void writeValue()
  {
  }
};

/**
 * @class
 *  Computes average velocity for the group of atoms in the specified region
 */
class VelocityCalculator : public ValueCalculator
{
  double m_localAvgVel[3];
  double m_avgVel[3];
protected:

  VelocityCalculator(class LAMMPS * lmp, int groupbit, int nevery);

  virtual void calculateLocalValue(int i);

  void calculateGlobalValue(int globalCount);

  void writeValue();

  void globalAfterRun() {};
};
}

#endif /* VALUECALCULATOR_H_ */
