//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

/** @file NamedObjectBase.h
 * @brief declaration of NamedObjectBase
 */

#ifndef QMCPLUSPLUS_NAMEDOBJECTBASE_H
#define QMCPLUSPLUS_NAMEDOBJECTBASE_H
#include <string>

namespace qmcplusplus
{

/** Base class for any object which needs to keep an object name
 */
class NamedObjectBase
{
public:

  ///default constructor
  NamedObjectBase() {}

  ///virtual destructor
  virtual ~NamedObjectBase() {}

  ///return the name
  inline const std::string& getName() const
  {
    return myName;
  }

  inline void setName(const std::string& aname)
  {
    myName = aname;
  }

protected:
  /** name of the object */
  std::string myName;
};

}
#endif // QMCPLUSPLUS_NAMEDOBJECTBASE_H
