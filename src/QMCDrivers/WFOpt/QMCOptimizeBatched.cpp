//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCDrivers/WFOpt/QMCOptimizeBatched.h"
#include "Particle/HDFWalkerIO.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#include "Optimize/CGOptimization.h"
#include "Optimize/testDerivOptimization.h"
#include "Optimize/DampedDynamics.h"
#include "QMCDrivers/VMC/VMCBatched.h"
#include "QMCDrivers/WFOpt/QMCCostFunction.h"
#include "QMCHamiltonians/HamiltonianPool.h"

namespace qmcplusplus
{
QMCOptimizeBatched::QMCOptimizeBatched(MCWalkerConfiguration& w,
                                       TrialWaveFunction& psi,
                                       QMCHamiltonian& h,
                                       HamiltonianPool& hpool,
                                       WaveFunctionPool& ppool,
                                       QMCDriverInput&& qmcdriver_input,
                                       VMCDriverInput&& vmcdriver_input,
                                       MCPopulation& population,
                                       SampleStack& samples,
                                       Communicate* comm)
    : QMCDriver(w, psi, h, ppool, comm),
      PartID(0),
      NumParts(1),
      hamPool(hpool),
      optSolver(0),
      vmcEngine(0),
      wfNode(NULL),
      optNode(NULL),
      qmcdriver_input_(qmcdriver_input),
      vmcdriver_input_(vmcdriver_input),
      population_(population),
      samples_(samples)
{
  IsQMCDriver = false;
  //set the optimization flag
  qmc_driver_mode.set(QMC_OPTIMIZE, 1);
  //read to use vmc output (just in case)
  RootName = "pot";
  QMCType  = "QMCOptimizeBatched";
  //default method is cg
  optmethod = "cg";
}

/** Clean up the vector */
QMCOptimizeBatched::~QMCOptimizeBatched()
{
  delete vmcEngine;
  delete optSolver;
}

/** Add configuration files for the optimization
 * @param a root of a hdf5 configuration file
 */
void QMCOptimizeBatched::addConfiguration(const std::string& a)
{
  if (a.size())
    ConfigFile.push_back(a);
}

/** Reimplement QMCDriver::run
 */
bool QMCOptimizeBatched::run()
{
  QMCCostFunction* cost = dynamic_cast<QMCCostFunction*>(optTarget.get());
  cost->makeClones(W, Psi, H);

  //generate samples
  generateSamples();

  NumOfVMCWalkers = W.getActiveWalkers();

  app_log() << "<opt stage=\"setup\">" << std::endl;
  app_log() << "  <log>" << std::endl;
  //reset the rootname
  optTarget->setRootName(RootName);
  optTarget->setWaveFunctionNode(wfNode);
  optTarget->setRng(vmcEngine->getRng());
  app_log() << "   Reading configurations from h5FileRoot " << std::endl;
  //get configuration from the previous run
  Timer t1;
  optTarget->getConfigurations(h5FileRoot);
  optTarget->checkConfigurations();
  app_log() << "  Execution time = " << std::setprecision(4) << t1.elapsed() << std::endl;
  app_log() << "  </log>" << std::endl;
  app_log() << "</opt>" << std::endl;
  app_log() << "<opt stage=\"main\" walkers=\"" << optTarget->getNumSamples() << "\">" << std::endl;
  app_log() << "  <log>" << std::endl;
  optTarget->setTargetEnergy(branchEngine->getEref());
  t1.restart();
  bool success = optSolver->optimize(optTarget.get());
  app_log() << "  Execution time = " << std::setprecision(4) << t1.elapsed() << std::endl;
  ;
  app_log() << "  </log>" << std::endl;
  optTarget->reportParameters();

  int nw_removed = W.getActiveWalkers() - NumOfVMCWalkers;
  app_log() << "   Restore the number of walkers to " << NumOfVMCWalkers << ", removing " << nw_removed << " walkers."
            << std::endl;
  if (nw_removed > 0)
    W.destroyWalkers(nw_removed);
  else
    W.createWalkers(-nw_removed);

  app_log() << "</opt>" << std::endl;
  app_log() << "</optimization-report>" << std::endl;
  MyCounter++;
  return (optTarget->getReportCounter() > 0);
}

void QMCOptimizeBatched::generateSamples()
{
  Timer t1;
  app_log() << "<optimization-report>" << std::endl;

  vmcEngine->qmc_driver_mode_.set(QMC_WARMUP, 1);
  vmcEngine->qmc_driver_mode_.set(QMC_OPTIMIZE, 1);
  vmcEngine->qmc_driver_mode_.set(QMC_WARMUP, 0);

  // TODO - understand what this does
  //vmcEngine->setValue("current", 0); //reset CurrentStep
  app_log() << "<vmc stage=\"main\" blocks=\"" << nBlocks << "\">" << std::endl;
  t1.restart();
  branchEngine->flush(0);
  branchEngine->reset();
  vmcEngine->run();
  app_log() << "  Execution time = " << std::setprecision(4) << t1.elapsed() << std::endl;
  app_log() << "</vmc>" << std::endl;
  //write parameter history and energies to the parameter file in the trial wave function through opttarget
  // TODO - understand what this does - values are not used, but the side effects may be important?
  //FullPrecRealType e, w, var;
  //vmcEngine->Estimators->getEnergyAndWeight(e, w, var);


  h5FileRoot = RootName;
}

void QMCOptimizeBatched::endSection() { vmcEngine->endSection(); }

/** Parses the xml input file for parameter definitions for the wavefunction optimization.
 * @param q current xmlNode
 * @return true if successful
 */
bool QMCOptimizeBatched::put(xmlNodePtr q)
{
  std::string vmcMove("pbyp");
  std::string useGPU("no");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(vmcMove, "move");
  oAttrib.add(useGPU, "gpu");
  oAttrib.put(q);
  xmlNodePtr qsave = q;
  xmlNodePtr cur   = qsave->children;
  int pid          = OHMMS::Controller->rank();
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "mcwalkerset")
    {
      mcwalkerNodePtr.push_back(cur);
    }
    else if (cname.find("optimize") < cname.size())
    {
      const XMLAttrString att(cur, "method");
      if (!att.empty())
        optmethod = att;
      optNode = cur;
    }
    cur = cur->next;
  }
  //no walkers exist, add 10
  if (W.getActiveWalkers() == 0)
    addWalkers(omp_get_max_threads());
  NumOfVMCWalkers = W.getActiveWalkers();
  //create VMC engine
  if (vmcEngine == 0)
  {
    QMCDriverInput qmcdriver_input_copy = qmcdriver_input_;
    VMCDriverInput vmcdriver_input_copy = vmcdriver_input_;
    vmcEngine = new VMCBatched(std::move(qmcdriver_input_copy), std::move(vmcdriver_input_copy), population_, Psi, H,
                               psiPool, samples_, myComm);

    vmcEngine->setUpdateMode(vmcMove[0] == 'p');
  }
  vmcEngine->setStatus(RootName, h5FileRoot, AppendRun);
  vmcEngine->process(qsave);
  if (optSolver == 0)
  {
    if (optmethod == "anneal")
    {
      app_log() << " Annealing optimization using DampedDynamics" << std::endl;
      optSolver = new DampedDynamics<RealType>;
    }
    else if ((optmethod == "flexOpt") | (optmethod == "flexopt") | (optmethod == "macopt"))
    {
      app_log() << "Conjugate-gradient optimization using FlexOptimization" << std::endl;
      app_log() << " This method has been removed. " << std::endl;
      APP_ABORT("QMCOptimizeBatched::put");
    }
    else if (optmethod == "BFGS")
    {
      app_log() << " This method is not implemented correctly yet. " << std::endl;
      APP_ABORT("QMCOptimizeBatched::put");
    }
    else if (optmethod == "test")
    {
      app_log() << "Conjugate-gradient optimization using tester Optimization: " << std::endl;
      optSolver = new testDerivOptimization<RealType>;
    }
    else
    {
      app_log() << " Conjugate-gradient optimization using CGOptimization" << std::endl;
      optSolver = new CGOptimization<RealType>;
    } //set the stream
    optSolver->setOstream(&app_log());
  }
  if (optNode == NULL)
    optSolver->put(qmcNode);
  else
    optSolver->put(optNode);
  bool success = true;
  //allways reset optTarget
  optTarget = std::make_unique<QMCCostFunction>(W, Psi, H, myComm);
  optTarget->setStream(&app_log());
  success = optTarget->put(q);

  return success;
}
} // namespace qmcplusplus
