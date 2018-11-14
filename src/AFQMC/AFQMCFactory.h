
// -*- C++ -*-
/**@file AFQMCFactory.h
 * @brief Top level class for AFQMC. Parses input and performs setup of classes. 
 */

#ifndef QMCPLUSPLUS_AFQMCFACTORY_H
#define QMCPLUSPLUS_AFQMCFACTORY_H

#include<string>
#include<vector>
#include<map>
#include <queue>
#include<algorithm>
#include<Message/MPIObjectBase.h>
#include "OhmmsApp/RandomNumberControl.h"

#include"alf/boost/mpi3/communicator.hpp"

#include "config.h"

#include "AFQMC/Utilities/taskgroup.h"
#include<AFQMC/Walkers/WalkerSetFactory.hpp>
#include<AFQMC/Hamiltonians/HamiltonianFactory.h>
#include<AFQMC/Wavefunctions/WavefunctionFactory.h>
#include<AFQMC/Propagators/PropagatorFactory.h>
#include<AFQMC/Drivers/DriverFactory.h>

#include "OhmmsData/libxmldefs.h"

namespace qmcplusplus
{

namespace afqmc
{

class AFQMCFactory
{

  public:

    ///constructor
    AFQMCFactory(boost::mpi3::communicator& comm_):
        m_series(0),project_title("afqmc"),
        gTG(comm_),
        TGHandler(gTG,-10),
        InfoMap(),
        HamFac(InfoMap),
        WSetFac(InfoMap),
        WfnFac(InfoMap),
        PropFac(InfoMap),
        DriverFac(gTG,TGHandler,InfoMap,WSetFac,PropFac,WfnFac,HamFac) 
    {
      TimerManager.set_timer_threshold(timer_level_coarse);
      setup_timers(AFQMCTimers, AFQMCTimerNames,timer_level_coarse);
    }

    ///destructor
    ~AFQMCFactory() {}

    /* 
     *  Parses xml input and creates all non-executable objects. 
     *  Created objects (pointers actually) are stored in maps based on name in xml block.
     *  Executable sections (drivers) are created with objects already exiting
     *  in the maps. 
     */  
    bool parse(xmlNodePtr cur);

    /*
     *  Parses xml input and creates executable sections, using objects created during parsing. 
     */   
    bool execute(xmlNodePtr cur);

  private:

    int m_series;
    std::string project_title;

    // global TG from which all TGs are constructed
    GlobalTaskGroup gTG;

    // object that manages the TGs. Must be placed here, 
    // since it must be destroyed last
    TaskGroupHandler TGHandler;

    // container of AFQMCInfo objects 
    std::map<std::string,AFQMCInfo> InfoMap; 

    // Hamiltonian factory 
    HamiltonianFactory HamFac; 

    // WalkerHandler factory 
    WalkerSetFactory WSetFac; 

    // Wavefunction factoru
    WavefunctionFactory WfnFac;

    // Propagator factoru
    PropagatorFactory PropFac;

    // driver factory
    DriverFactory DriverFac;

};
}
}

#endif
