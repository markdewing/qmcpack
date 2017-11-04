//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"


#include "Utilities/RandomGenerator.h"
#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Utilities/OhmmsInfo.h"
#include "Lattice/ParticleBConds.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "Particle/SymmetricDistanceTableData.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCDrivers/WalkerControlBase.h"


#include <stdio.h>
#include <string>
#include <random>


using std::string;

namespace qmcplusplus
{

// add declaration here so it's accessible for testing
void assignWalkerPopulation(int Cur_pop, int NumContexts, int MyContext, std::vector<int> NumPerNode, std::vector<int> &minus, std::vector<int> &plus);

void output_vector(const std::string &name, std::vector<int> &vec)
{
  std::cout << name;
  for (int i = 0; i < vec.size(); i++) {
    std::cout << vec[i] << " ";
  }
  std::cout << std::endl;
}

// uncomment the std::cout and output_vector lines to see the walker assignments
TEST_CASE("Walker control assign walkers", "[drivers][walker_control]")
{
  int Cur_pop = 8;
  int NumContexts = 4;
  int MyContext = 0;
  std::vector<int> NumPerNode = {4,4,0,0};

  std::vector<int> NewNum = NumPerNode;
  for (int me = 0; me < NumContexts; me++)
  {
    std::vector<int> minus;
    std::vector<int> plus;

    //std::cout << "For processor number " << me << std::endl;
    assignWalkerPopulation(Cur_pop, NumContexts, me, NumPerNode, minus, plus);

    //output_vector("  Minus: ", minus);

    //output_vector("  Plus: ", plus);

    for (int i = 0; i < plus.size(); i++) {
      if (me == plus[i]) NewNum[plus[i]]--;
    }
    for (int i = 0; i < minus.size(); i++) {
      if (me == minus[i]) NewNum[minus[i]]++;
    }
  }
  //output_vector("New num per node: ", NewNum);

  std::vector<int> FairOffset;
  FairDivideLow(Cur_pop,NumContexts,FairOffset);
  for (int i = 0; i < NewNum.size(); i++) {
    int num = FairOffset[i+1] - FairOffset[i];
    REQUIRE(NewNum[i] == num);
  }
}

#if 0
// More of a property-based test than a unit test
// Eventually should make a framework for holding these kinds of tests.
// Constructs random number of processors with a random number of walkers on each node
// Ensures each node ultimately has the number of walkers choosen by FairDivideLow.


TEST_CASE("Walker control assign walkers many", "[drivers][walker_control]")
{
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_int_distribution<int> NumNodes(1,1000);
  std::uniform_int_distribution<int> TotalPop(0,1000);

  for (int nt = 0; nt < 1000; nt++) {
    int NumContexts = NumNodes(mt);
    int Cur_pop = TotalPop(mt);
    std::uniform_int_distribution<int> WalkerPop(0,2*Cur_pop/NumContexts);
    std::vector<int> NumPerNode(NumContexts);
    int current_pop = Cur_pop;
    for (int i = 0; i < NumContexts; i++) {
      int p = WalkerPop(mt);
      p = std::min(current_pop, p);
      current_pop -= p;
      // Make sure all walkers are accounted for on the last node
      if (i == NumContexts-1 && current_pop > 0) {
        p += current_pop;
      }
      NumPerNode[i] = p;
    }
    //std::cout << "NumNodes = " << NumContexts << std::endl;
    //std::cout << "TotalPop = " << Cur_pop << std::endl;
    //output_vector("Start: ",NumPerNode);

    std::vector<int> NewNum = NumPerNode;
    for (int me = 0; me < NumContexts; me++)
    {
      std::vector<int> minus;
      std::vector<int> plus;

      assignWalkerPopulation(Cur_pop, NumContexts, me, NumPerNode, minus, plus);

      for (int i = 0; i < plus.size(); i++) {
        if (me == plus[i]) NewNum[plus[i]]--;
      }
      for (int i = 0; i < minus.size(); i++) {
        if (me == minus[i]) NewNum[minus[i]]++;
      }
    }

    std::vector<int> FairOffset;
    FairDivideLow(Cur_pop,NumContexts,FairOffset);
    for (int i = 0; i < NewNum.size(); i++) {
      int num = FairOffset[i+1] - FairOffset[i];
      REQUIRE(NewNum[i] == num);
    }
  }
}
#endif
}
