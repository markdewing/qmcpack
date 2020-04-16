//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_DUMPDERIVOPTIMIZATION_H
#define QMCPLUSPLUS_DUMPDERIVOPTIMIZATION_H

#include "Optimize/OptimizeBase.h"
#include "OhmmsData/ParameterSet.h"

#include <time.h>


template<class T>
class dumpDerivOptimization : public MinimizerBase<T>
{
public:
  typedef T Return_t;
  typedef typename MinimizerBase<T>::ObjectFuncType ObjectFuncType;
  using MinimizerBase<T>::msg_stream;


  ObjectFuncType* TargetFunc;
  int NumParams;
  Return_t Displacement;
  Return_t CostTol;
  Return_t GradTol;
  Return_t deltaG;
  std::vector<Return_t> a_xi, Parms;

  dumpDerivOptimization(ObjectFuncType* atarget = 0)
      : TargetFunc(atarget),
        Displacement(0)
  {
    if (atarget)
      setTarget(atarget);
  }

  ~dumpDerivOptimization() {}

  void setTarget(ObjectFuncType* fn)
  {
    TargetFunc = fn;
    NumParams  = TargetFunc->getNumParams();
    resizeAllArray(NumParams);
    for (int i = 0; i < NumParams; i++)
      Parms[i] = std::real(TargetFunc->Params(i));
  }

  void resizeAllArray(int newSize)
  {
    a_xi.resize(newSize, 0.0);
    Parms.resize(newSize, 0.0);

  }

  bool optimize(ObjectFuncType* fn)
  {
    setTarget(fn);
    return optimize();
  }

  bool optimize()
  {
    std::cout << "Param    Numeric            Analytic       Percent" << std::endl;
#if 0
    double orig = Parms[0];
    for (int i = 0; i < 10; i++) {
      Parms[0] = orig + 0.1*i;
      for (int i = 0; i < NumParams; i++)
        TargetFunc->Params(i) = Parms[i];
      dfunc(Parms, a_xi);
    }
#else
    dfunc(Parms, a_xi);
#endif 
    //make sure wave function has the right parameters for next runs
    for (int i = 0; i < NumParams; i++)
      TargetFunc->Params(i) = Parms[i];
    TargetFunc->Report();
    return true;
  }

  bool get(std::ostream&) const;

  bool put(std::istream&);

  /**  Parse the xml file for parameters
   * @param cur current xmlNode
   * @param a_itmax maximum number of CG iterations
   * @param Displacement used for finite difference cost function
   * @param CG_ortho number of previous search directions we want to maintain orthogonality to.
   * @param a_rich 0=used approximate gradient at new location, 1=calculate new gradient at line minimum
   * @param xybisect Number of times to use bisection before linear interpolation for line minimization minimum finding routine.
   * @param Gfactor max allowed increase in Cost function gradient. Avoids falling in holes.
   * @param xycleanup max steps to tighten the line search
   * @param GradTol gradient to quit line minimization
   * @param a_verbose 0=quiet, 1=normal, 2=chatty, 3+=various degrees of debug (loud!)
   * @param a_lastx_default default step size. Only has a transitive effect.
   * @param a_linmin_maxits Maximum number of steps to try to find a minimum along a_h direction
   * @param deltaG if change in gradient is less than this relax the CG_ortho constraints
   *
   */

  bool put(xmlNodePtr cur)
  {
    return true;
  }

  Return_t func(std::vector<Return_t> _p)
  {
    for (int i = 0; i < NumParams; ++i)
      TargetFunc->Params(i) = _p[i];
    return TargetFunc->Cost();
  }


  void dfunc(std::vector<Return_t> RT, std::vector<Return_t>& FG)
  {
    ///To test we simply output the analytic and numeric gradients of the cost function. Make sure they agree.
    std::vector<Return_t> Dummy(FG);
    TargetFunc->GradCost(Dummy, RT, 1e-5);
    TargetFunc->GradCost(FG, RT, 0.0);

    //std::cout << "Param    Numeric            Analytic       Percent" << std::endl;
    //for (int k = 0; k < NumParams; k++)
    for (int k = 0; k < 1; k++)
    {
      if (Dummy[k] != 0)
        std::cout << RT[k] << "  " << Dummy[k] << "  " << FG[k] << "  " << 100 * (Dummy[k] - FG[k]) / Dummy[k] << std::endl;
      else
        std::cout << RT[k] << "  " << Dummy[k] << "  " << FG[k] << "   inf" << std::endl;
    }
    std::cout << std::endl;
  }

};

#endif
