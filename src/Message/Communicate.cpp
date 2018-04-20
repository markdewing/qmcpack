//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    




#include <Configuration.h>
#include "Message/Communicate.h"
#include "Message/TagMaker.h"
#include <iostream>
#include <cstdio>
#include <Platforms/sysutil.h>
#include <tau/profiler.h>
#include <Utilities/UtilityFunctions.h>
#include <fstream>

#include "mpi_wrapper/boost/mpi3/communicator.hpp"

#ifdef HAVE_ADIOS
#include <adios.h>
#include <adios_read.h>
#include <adios_error.h>
#include "ADIOS/ADIOS_config.h"
#endif




//static data of TagMaker::CurrentTag is initialized.
int TagMaker::CurrentTag = 1000;

//Global Communicator is created without initialization
Communicate* OHMMS::Controller = new Communicate;
//boost::mpi3::environment* OHMMS::Environment = nullptr;

//default constructor: ready for a serial execution
Communicate::Communicate():
  myMPI(0), d_mycontext(0), d_ncontexts(1), d_groupid(0), d_ngroups(1), GroupLeaderComm(nullptr)
{
}

//Communicate::Communicate(int argc, char **argv):
Communicate::Communicate(boost::mpi3::environment &env):
  GroupLeaderComm(nullptr)
{
  initialize(env);
}

Communicate::~Communicate()
{
  if(GroupLeaderComm!=nullptr) delete GroupLeaderComm;
}

//exclusive:  OOMPI, MPI or Serial

Communicate::Communicate(const mpi_comm_type comm_input):
  myMPI(comm_input), d_groupid(0), d_ngroups(1), GroupLeaderComm(nullptr)
{
  myComm=OOMPI_Intra_comm(myMPI);
  d_mycontext=myComm.Rank();
  d_ncontexts=myComm.Size();
}


Communicate::Communicate(const Communicate& in_comm, int nparts)
{
  std::vector<int> nplist(nparts+1);

  //this is a workaround due to the OOMPI bug with split
  if(nparts>1)
  {
    int p=FairDivideLow(in_comm.rank(), in_comm.size(), nparts, nplist); //group
    int q=in_comm.rank()-nplist[p];//rank within a group
    //int n=comm.size()/nparts;
    //int p=comm.rank()/n;
    //int q=comm.rank()%n;
#if 0
    MPI_Comm row;
    MPI_Comm_split(in_comm.getMPI(),p,q,&row);
    myComm=OOMPI_Intra_comm(row);
#else
    comm = in_comm.comm.split(p,q);
    myComm=OOMPI_Intra_comm(&comm);
#endif
    //comm = row;
    d_groupid=p;
  }
  else
  {
    nplist[0]=0; nplist[1]=in_comm.size();
    myComm=OOMPI_Intra_comm(in_comm.getComm());
    d_groupid=0;
  }
  myMPI = myComm.Get_mpi();
  d_mycontext=myComm.Rank();
  d_ncontexts=myComm.Size();
  d_ngroups=nparts;
  // create a communicator among group leaders.
  MPI_Group parent_group, leader_group;
  MPI_Comm_group(in_comm.getMPI(), &parent_group);
  MPI_Group_incl(parent_group, nparts, nplist.data(), &leader_group);
  MPI_Comm leader_comm;
  MPI_Comm_create(in_comm.getMPI(), leader_group, &leader_comm);
  if(isGroupLeader())
    GroupLeaderComm = new Communicate(leader_comm);
  else
    GroupLeaderComm = nullptr;
  MPI_Group_free(&parent_group);
  MPI_Group_free(&leader_group);
}



//================================================================
// Implements Communicate with OOMPI library
//================================================================

//void Communicate::initialize(int argc, char **argv)
// Empty to keep myriad unit tests happy for now
void Communicate::initialize(int argc, char **argv)
{
}

void Communicate::initialize(boost::mpi3::environment &env)
{
  //myComm = OOMPI_Intra_comm(&OHMMS::Environment->world());
  myComm = OOMPI_Intra_comm(&env.world());
  comm = env.world();
  //OOMPI_COMM_WORLD.Init(argc, argv);
  //myComm = OOMPI_COMM_WORLD;
  myMPI = myComm.Get_mpi();
  //d_mycontext = OOMPI_COMM_WORLD.Rank();
  //d_ncontexts = OOMPI_COMM_WORLD.Size();
  //d_mycontext = OHMMS::Environment->world().rank();
  //d_ncontexts = OHMMS::Environment->world().size();
  d_mycontext = env.world().rank();
  d_ncontexts = env.world().size();
  d_groupid=0;
  d_ngroups=1;
#ifdef __linux__
  for (int proc=0; proc<OHMMS::Controller->size(); proc++)
  {
    if (OHMMS::Controller->rank() == proc)
    {
      fprintf (stderr, "Rank = %4d  Free Memory = %5zu MB\n", proc, freemem());
    }
    barrier();
  }
  barrier();
#endif
  std::string when="qmc."+getDateAndTime("%Y%m%d_%H%M");
}



void Communicate::cleanupMessage(void*)
{
}

void Communicate::abort()
{
  //OOMPI_COMM_WORLD.Abort();
  //OHMMS::Environment->world().abort();
}

void Communicate::barrier()
{
  myComm.Barrier();
}

void Communicate::abort(const char* msg)
{
  std::cerr << msg << std::endl;
  //OOMPI_COMM_WORLD.Abort();
  //OHMMS::Environment->world().abort();
}


