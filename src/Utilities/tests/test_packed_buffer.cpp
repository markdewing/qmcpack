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

#include "Utilities/PackedBuffer.h"
#include <stdio.h>
#include <string>
#include <vector>
#include <random>

using std::string;
using std::vector;
namespace qmcplusplus {

TEST_CASE("pack int", "[utilities]")
{
  PackedBuffer p(sizeof(int));
  int i = 1;
  p << i;

  p.reset();
  int j= 0;
  p >> j;
  REQUIRE(i==j);
}

TEST_CASE("pack double", "[utilities]")
{
  PackedBuffer p(sizeof(double));
  double i = 1.1;
  p << i;

  p.reset();
  double j= 0;
  p >> j;
  // No Approx here - expect bitwise equality
  REQUIRE(i==j);
}

TEST_CASE("pack double array", "[utilities]")
{
  PackedBuffer p(3*sizeof(double));
  double tmp[3] = {1.1, 1.2, 1.3};
  p.Pack(tmp, 3);

  p.reset();
  double tmp2[3] = {0.0, 0.0, 0.0};
  p.Unpack(tmp2, 3);
  for (int i = 0; i < 3; i++) {
    REQUIRE(tmp[i] == tmp2[i]);
  }
}


#if 0
// These are more similar to a property-based test than a unit test.
// Need to build a separate way to run these tests - infrequent, usually only
// run when changing the tested code or on new platforms. Much longer-running
// than unit tests.

TEST_CASE("pack int many", "[utilities]")
{
  PackedBuffer p(sizeof(int));
  for (int i = 0; i < 10000; i++) {
    p.reset();
    p << i;

    p.reset();
    int j= 0;
    p >> j;
    REQUIRE(i==j);
  }
}

TEST_CASE("pack vector double many", "[utilities]")
{
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_int_distribution<int> vsize(0,100);
  std::uniform_real_distribution<double> values(-std::numeric_limits<double>::max(),std::numeric_limits<double>::max());
  for (int i = 0; i < 10000; i++) {
    // create vector of random size and populate it will random double values
    int size = vsize(mt);
    vector<double> v(size);
    for (int j = 0; j < size; j++) {
      v[j] = values(mt);
    }

    PackedBuffer p(v.size()*sizeof(v[0]));

    p.Pack(v.data(), v.size());

    p.reset();

    vector<double> w(size);
    p.Unpack(w.data(), size);

    for (int j = 0; j < size; j++) {
      REQUIRE(v[j] == w[j]);
    }
  }
}
#endif



}
