
#ifndef QMCPLUSPLUS_AFQMC_WAVEFUNCTIONFACTORY_H
#define QMCPLUSPLUS_AFQMC_WAVEFUNCTIONFACTORY_H

#include<iostream>
#include<vector> 
#include<map> 
#include<fstream>
#include "OhmmsData/libxmldefs.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Wavefunctions/WavefunctionBase.hpp"
#include "AFQMC/Wavefunctions/SlaterDetOperations.hpp"
#include "AFQMC/Hamiltonians/HamiltonianBase.hpp"

namespace qmcplusplus
{

namespace afqmc
{

// keep a std::map<*AFQMCInfo,SlaterDetOperations> to construct Wfns, and route all determinant operations through this object in Wfn classes


class WavefunctionFactory
{

  using HamPtr = std::shared_ptr<HamiltonianBase>;

  public:

  WavefunctionFactory(std::map<std::string,AFQMCInfo*>& info):InfoMap(info)
  {
  }

  ~WavefunctionFactory()
  {
  }

  // returns a pointer to the base Wavefunction class associated with a given ID 
  std::shared_ptr<WavefunctionBase> getWavefunction(afqmc::TaskGroup_& TG, const std::string& ID, HamPtr h)
  {
    auto xml = xmlBlocks.find(ID);
    if(xml == xmlBlocks.end())
      return nullptr;
    auto w0 = wavefunctions.find(ID);
    if( w0 == wavefunctions.end() ) {
      auto neww = wavefunctions.insert(std::make_pair(ID,buildWavefunction(TG,xml->second,h)));
      if(!neww.second)
        APP_ABORT(" Error: Problems building new wavefunction in WavefunctionFactory::getWavefunction(string&). \n");
      return (neww.first)->second;  
    } else
      return w0->second;
  }

  // returns the xmlNodePtr associated with ID
  xmlNodePtr getXML(const std::string& ID)
  {
    auto xml = xmlBlocks.find(ID);
    if(xml == xmlBlocks.end())
      return nullptr;
    else
      return xml->second;
  }

  // adds a xml block from which a Wavefunction can be built
  void push(const std::string& ID, xmlNodePtr cur)
  {
    auto xml = xmlBlocks.find(ID);
    if(xml != xmlBlocks.end())
      APP_ABORT("Error: Repeated Wavefunction block in WavefunctionFactory. Wavefunction names must be unique. \n");
    xmlBlocks.insert(std::make_pair(ID,cur));
  }

  protected:

  // reference to container of AFQMCInfo objects 
  std::map<std::string,AFQMCInfo*>& InfoMap;

  // generates a new Wavefunction and returns the pointer to the base class
  std::shared_ptr<WavefunctionBase> buildWavefunction(afqmc::TaskGroup_& TG, xmlNodePtr cur, HamPtr h)
  {
    std::string type;
    ParameterSet m_param;
    m_param.add(type,"filetype","std::string");
    m_param.put(cur);

    if(type == "fcidump" || type == "ascii")
      return fromASCII(TG,cur,h);
    else if(type == "hdf5")
      return fromHDF5(TG,cur,h);
    else {
      app_error()<<"Unknown Wavefunction filetype in WavefunctionFactory::buildWavefunction(): " <<type <<std::endl;
      APP_ABORT(" Error: Unknown Wavefunction filetype in WavefunctionFactory::buildWavefunction(). \n");
    }
  }

  std::shared_ptr<WavefunctionBase> fromASCII(afqmc::TaskGroup_& TG, xmlNodePtr cur, HamPtr h);

  std::shared_ptr<WavefunctionBase> fromHDF5(afqmc::TaskGroup_& TG, xmlNodePtr cur, HamPtr h);

  // MAM: should I store a copy rather than a pointer???
  std::map<std::string,xmlNodePtr> xmlBlocks;

  std::map<std::string,std::shared_ptr<WavefunctionBase>> wavefunctions;

  //std::map<AFQMCInfo,SlaterDetOperations>

};
}
}

#endif
