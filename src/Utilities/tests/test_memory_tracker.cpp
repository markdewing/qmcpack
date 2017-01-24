//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, Argonne National Laboratory
//
// File created by: Mark Dewing, markdewing@gmail.com, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Utilities/OhmmsInfo.h"
#include "Utilities/MemoryTracker.h"
#include <stdio.h>
#include <string>
#include <vector>

// Since they aren't public declarations, we put them here in order
//  to test them
bool matches_object_suffix(const std::string &in, std::string &base);
void split_stack(const std::string &in, std::string &top, std::string &rest);

namespace qmcplusplus {


TEST_CASE("test_split_stack", "[utilities]")
{
  std::string t1("stack");
  std::string top;
  std::string rest;
  split_stack(t1, top, rest);
  REQUIRE(top == "stack");
  REQUIRE(rest == "");

  std::string t2("stack/level2");
  top = "";
  rest = "";
  split_stack(t2, top, rest);
  REQUIRE(top == "stack");
  REQUIRE(rest == "level2");

  std::string t3("stack/level2/level3");
  top = "";
  rest = "";
  split_stack(t3, top, rest);
  REQUIRE(top == "stack");
  REQUIRE(rest == "level2/level3");
}

TEST_CASE("test_matches_object_suffix", "[utilities]")
{
  std::string t1("stack");
  std::string base;
  bool matches = matches_object_suffix(t1, base);
  REQUIRE(matches == false);
  REQUIRE(base == "stack");

  std::string t2("stack_");
  base = "";
  matches = matches_object_suffix(t2, base);
  REQUIRE(matches == false);
  REQUIRE(base == "stack_");

  std::string t3("stack_1");
  base = "";
  matches = matches_object_suffix(t3, base);
  REQUIRE(matches == true);
  REQUIRE(base == "stack");

  std::string t4("stack_10");
  base = "";
  matches = matches_object_suffix(t4, base);
  REQUIRE(matches == true);
  REQUIRE(base == "stack");

  std::string t5("stack_10t");
  base = "";
  matches = matches_object_suffix(t5, base);
  REQUIRE(matches == false);
  REQUIRE(base == "stack_10t");
}


}
