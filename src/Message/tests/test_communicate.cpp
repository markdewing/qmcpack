//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "Message/catch_mpi_main.hpp"


#include "Message/Communicate.h"

#include <stdio.h>
#include <string>
#include <sstream>


using std::cout;
using std::endl;

TEST_CASE("split constructor", "[message]")
{

  Communicate *c = OHMMS::Controller;
#if 0
  cout << "initial rank = " << c->rank() << endl;
  cout << "initial size = " << c->size() << endl;
#endif
  int initial_size = c->size();
  int initial_rank = c->rank();
  if (c->size() > 1)
  {
    // Get some segfaults in the destructor - to avoid for now, create object and let it leak
    Communicate *split = new Communicate(*c, 2);
    int new_size = initial_size/2;
    if (initial_rank >= new_size) new_size += initial_size%2;
#if 0
    cout << endl;
    cout << "initial rank = " << c->rank() << endl;
    cout << "rank = " << split.rank() << endl;
    cout << "size = " << split.size() << " new size = " << new_size << endl;
    cout << endl;
#endif
    REQUIRE(split->size() == new_size);

    if (split->rank() == 0) {
      REQUIRE(split->GroupLeaderComm != nullptr);
    } else {
      REQUIRE(split->GroupLeaderComm == nullptr);
    }
  }
}

