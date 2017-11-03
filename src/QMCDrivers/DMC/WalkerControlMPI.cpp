//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include <QMCDrivers/DMC/WalkerControlMPI.h>
#include <qmc_common.h>
#include <Utilities/IteratorUtility.h>
#include <Utilities/UtilityFunctions.h>
#include <Utilities/NewTimer.h>
#include <Utilities/Timer.h>
#include <Utilities/PackedBuffer.h>

// temporary
#include <QMCDrivers/DMC/PackWrapper.h>

namespace qmcplusplus
{

//#if defined(PRINT_DEBUG)
//#define DMC_BRANCH_START(NOW) NOW
//#define DMC_BRANCH_STOP(TID,TM) TID=TM
//#define DMC_BRANCH_DUMP(IT,T1,T2,T3)  \
//  OhmmsInfo::Debug->getStream()  << "BRANCH " \
//<< std::setw(8) << IT \
//<< " SORT" << std::setw(16) << T1 \
//<< " GSUM" << std::setw(16) << T2 \
//<< " SWAP" << std::setw(16) << T3 << std::endl
//#else
#define DMC_BRANCH_START(NOW)
#define DMC_BRANCH_STOP(TID,TM)
#define DMC_BRANCH_DUMP(IT,T1,T2,T3)
//#endif


//#define MCWALKERSET_MPI_DEBUG

enum DMC_MPI_Timers
{
  DMC_MPI_branch,
  DMC_MPI_prebalance,
  DMC_MPI_loadbalance
};

TimerNameList_t<DMC_MPI_Timers> DMCMPITimerNames =
{
  {DMC_MPI_branch, "WalkerControlMPI::branch"},
  {DMC_MPI_prebalance, "WalkerControlMPI::pre-loadbalance"},
  {DMC_MPI_loadbalance, "WalkerControlMPI::loadbalance"}
};

/** default constructor
 *
 * set SwapMode
 */
WalkerControlMPI::WalkerControlMPI(Communicate* c): WalkerControlBase(c)
{
  SwapMode=1;
  Cur_min=0;
  Cur_max=0;
#ifdef MCWALKERSET_MPI_DEBUG
  char fname[128];
  sprintf(fname,"test.%d",MyContext);
  std::ofstream fout(fname);
#endif
  setup_timers(myTimers, DMCMPITimerNames, timer_level_medium);
}

int
WalkerControlMPI::branch(int iter, MCWalkerConfiguration& W, RealType trigger)
{
  DMC_BRANCH_START(Timer localTimer);
  TinyVector<RealType,3> bTime(0.0);
  myTimers[DMC_MPI_branch]->start();
  myTimers[DMC_MPI_prebalance]->start();
  std::fill(curData.begin(),curData.end(),0);
  //std::fill(NumPerNode.begin(),NumPerNode.end(),0);
  sortWalkers(W);
  //use NumWalkersSent from the previous exchange
  curData[SENTWALKERS_INDEX]=NumWalkersSent;
  //update the number of walkers for this node
  curData[LE_MAX+MyContext]=NumWalkers;
  DMC_BRANCH_STOP(bTime[0],localTimer.elapsed());
  DMC_BRANCH_START(localTimer.restart());
  int nw = copyWalkers(W);
  myComm->allreduce(curData);
  measureProperties(iter);
  W.EnsembleProperty=EnsembleProperty;
  DMC_BRANCH_STOP(bTime[1],localTimer.elapsed());
  DMC_BRANCH_START(localTimer.restart());
  ////update the samples and weights
  //W.EnsembleProperty.NumSamples=curData[WALKERSIZE_INDEX];
  //W.EnsembleProperty.Weight=curData[WEIGHT_INDEX];
  //RealType wgtInv(1.0/curData[WEIGHT_INDEX]);
  //accumData[ENERGY_INDEX]     += curData[ENERGY_INDEX]*wgtInv;
  //accumData[ENERGY_SQ_INDEX]  += curData[ENERGY_SQ_INDEX]*wgtInv;
  //accumData[WALKERSIZE_INDEX] += curData[WALKERSIZE_INDEX];
  //accumData[WEIGHT_INDEX]     += curData[WEIGHT_INDEX];
  Cur_pop=0;
  for(int i=0, j=LE_MAX; i<NumContexts; i++,j++)
  {
    Cur_pop+= NumPerNode[i]=static_cast<int>(curData[j]);
  }
  myTimers[DMC_MPI_prebalance]->stop();
  myTimers[DMC_MPI_loadbalance]->start();
  if(qmc_common.async_swap)
    swapWalkersAsync(W);
  else
    swapWalkersSimple(W);
  myTimers[DMC_MPI_loadbalance]->stop();
  //Do not need to use a trigger.
  //Cur_min=Nmax;
  //Cur_max=0;
  //Cur_pop=0;
  //for(int i=0, j=LE_MAX; i<NumContexts; i++,j++) {
  //  Cur_pop+= NumPerNode[i]=static_cast<int>(curData[j]);
  //  Cur_min = std::min(Cur_min,NumPerNode[i]);
  //  Cur_max = std::max(Cur_max,NumPerNode[i]);
  //}
  //int max_diff = std::max(Cur_max*NumContexts-Cur_pop,Cur_pop-Cur_min*NumContexts);
  //double diff_pop = static_cast<double>(max_diff)/static_cast<double>(Cur_pop);
  //if(diff_pop > trigger) {
  //  swapWalkersSimple(W);
  //  //swapWalkersMap(W);
  //}
  //set Weight and Multiplicity to default values
  MCWalkerConfiguration::iterator it(W.begin()),it_end(W.end());
  while(it != it_end)
  {
    (*it)->Weight= 1.0;
    (*it)->Multiplicity=1.0;
    ++it;
  }
  //update the global number of walkers and offsets
  W.setGlobalNumWalkers(Cur_pop);
  W.setWalkerOffsets(FairOffSet);
  DMC_BRANCH_STOP(bTime[2],localTimer.elapsed());
  DMC_BRANCH_DUMP(iter,bTime[0],bTime[1],bTime[2]);
  myTimers[DMC_MPI_branch]->stop();
  return Cur_pop;
}

/** swap Walkers with Recv/Send
 *
 * The algorithm ensures that the load per node can differ only by one walker.
 * The communication is one-dimensional.
 */
void WalkerControlMPI::swapWalkersSimple(MCWalkerConfiguration& W)
{
  FairDivideLow(Cur_pop,NumContexts,FairOffSet);
  std::vector<int> minus, plus;
  int deltaN;
  for(int ip=0; ip<NumContexts; ip++)
  {
    int dn=NumPerNode[ip]-(FairOffSet[ip+1]-FairOffSet[ip]);
    if(ip == MyContext)
      deltaN=dn;
    if(dn>0)
    {
      plus.insert(plus.end(),dn,ip);
    }
    else
      if(dn<0)
      {
        minus.insert(minus.end(),-dn,ip);
      }
  }
  Walker_t& wRef(*W[0]);
  std::vector<Walker_t*> newW;
  std::vector<Walker_t*> oldW;
#ifdef MCWALKERSET_MPI_DEBUG
  char fname[128];
  sprintf(fname,"test.%d",MyContext);
  std::ofstream fout(fname, std::ios::app);
  //fout << NumSwaps << " " << Cur_pop << " ";
  //for(int ic=0; ic<NumContexts; ic++) fout << NumPerNode[ic] << " ";
  //fout << " | ";
  //for(int ic=0; ic<NumContexts; ic++) fout << FairOffSet[ic+1]-FairOffSet[ic] << " ";
  //fout << " | ";
  for(int ic=0; ic<plus.size(); ic++)
  {
    fout << plus[ic] << " ";
  }
  fout << " | ";
  for(int ic=0; ic<minus.size(); ic++)
  {
    fout << minus[ic] << " ";
  }
  fout << std::endl;
#endif
  int nswap=std::min(plus.size(), minus.size());
  int last=W.getActiveWalkers()-1;
  int nsend=0;
  for(int ic=0; ic<nswap; ic++)
  {
    if(plus[ic]==MyContext)
    {
      PackedBuffer sendBuffer(wRef.byteSize());
      W[last]->putMessage(sendBuffer);
      myComm->getComm()[minus[ic]].Send(toMessage(sendBuffer));
      --last;
      ++nsend;
    }
    if(minus[ic]==MyContext)
    {
      PackedBuffer recvBuffer(wRef.byteSize());
      myComm->getComm()[plus[ic]].Recv(toMessage(recvBuffer));
      Walker_t *awalker= new Walker_t(wRef);
      awalker->getMessage(recvBuffer);
      newW.push_back(awalker);
    }
  }
  //save the number of walkers sent
  NumWalkersSent=nsend;
  if(nsend)
  {
    nsend=NumPerNode[MyContext]-nsend;
    W.destroyWalkers(W.begin()+nsend, W.end());
  }
  //add walkers from other node
  if(newW.size())
    W.insert(W.end(),newW.begin(),newW.end());
}

/** swap Walkers with Irecv/Send
 *
 * The algorithm ensures that the load per node can differ only by one walker.
 * The communication is one-dimensional.
 */
void WalkerControlMPI::swapWalkersAsync(MCWalkerConfiguration& W)
{
  FairDivideLow(Cur_pop,NumContexts,FairOffSet);
  std::vector<int> minus, plus;
  int deltaN;
  for(int ip=0; ip<NumContexts; ip++)
  {
    int dn=NumPerNode[ip]-(FairOffSet[ip+1]-FairOffSet[ip]);
    if(ip == MyContext)
      deltaN=dn;
    if(dn>0)
    {
      plus.insert(plus.end(),dn,ip);
    }
    else
      if(dn<0)
      {
        minus.insert(minus.end(),-dn,ip);
      }
  }
  Walker_t& wRef(*W[0]);
  std::vector<Walker_t*> newW;
  std::vector<Walker_t*> oldW;
#ifdef MCWALKERSET_MPI_DEBUG
  char fname[128];
  sprintf(fname,"test.%d",MyContext);
  std::ofstream fout(fname, std::ios::app);
  //fout << NumSwaps << " " << Cur_pop << " ";
  //for(int ic=0; ic<NumContexts; ic++) fout << NumPerNode[ic] << " ";
  //fout << " | ";
  //for(int ic=0; ic<NumContexts; ic++) fout << FairOffSet[ic+1]-FairOffSet[ic] << " ";
  //fout << " | ";
  for(int ic=0; ic<plus.size(); ic++)
  {
    fout << plus[ic] << " ";
  }
  fout << " | ";
  for(int ic=0; ic<minus.size(); ic++)
  {
    fout << minus[ic] << " ";
  }
  fout << std::endl;
#endif
  int nswap=std::min(plus.size(), minus.size());
  int last=W.getActiveWalkers()-1;
  int nsend=0;
  int countSend = 1;
  PackedBuffer ** sendBuffers = new PackedBuffer*[NumContexts];
  PackedBuffer ** recvBuffers = new PackedBuffer*[NumContexts];
  std::vector<OOMPI_Request> requests(NumContexts);
  std::vector<int> sendCounts(NumContexts,0);
  for(int ip=0; ip < NumContexts; ++ip)
    sendBuffers[ip] = 0;
  for(int ip=0; ip < NumContexts; ++ip)
    recvBuffers[ip] = 0;
  for(int ic=0; ic<nswap; ic++)
  {
    if(plus[ic]==MyContext)
    {
      if((ic < nswap - 1) && (plus[ic] == plus[ic+1]) && (minus[ic] == minus[ic+1]))
      {
        countSend++;
      }
      else
      {
        sendBuffers[minus[ic]] = new PackedBuffer(countSend * wRef.byteSize());
        for(int cs = 0; cs < countSend; ++cs)
        {
          W[last]->putMessage(*(sendBuffers[minus[ic]]));
          --last;
        }
        //OOMPI_COMM_WORLD[minus[ic]].Send(sendBuffer);
        requests[minus[ic]] = myComm->getComm()[minus[ic]].Isend(toMessage(*(sendBuffers[minus[ic]])), plus[ic]);
        nsend += countSend;
        countSend = 1;
      }
    }
    if(minus[ic]==MyContext)
    {
      if((ic < nswap - 1) && (plus[ic] == plus[ic+1]) && (minus[ic] == minus[ic+1]))
      {
        countSend++;
      }
      else
      {
        recvBuffers[plus[ic]] = new PackedBuffer(countSend * wRef.byteSize());
        //OOMPI_COMM_WORLD[plus[ic]].Recv(recvBuffer);
        requests[plus[ic]] = myComm->getComm()[plus[ic]].Irecv(toMessage(*(recvBuffers[plus[ic]])), plus[ic]);
        sendCounts[plus[ic]] = countSend;
        countSend = 1;
      }
    }
  }
  for(int ip = 0; ip < NumContexts; ++ip)
  {
    if(recvBuffers[ip])
    {
      requests[ip].Wait();
      for(int cs = 0; cs < sendCounts[ip]; ++cs)
      {
        Walker_t *awalker= new Walker_t(wRef);
        awalker->getMessage(*(recvBuffers[ip]));
        newW.push_back(awalker);
      }
      delete recvBuffers[ip];
    }
  }
  for(int ip = 0; ip < NumContexts; ++ip)
  {
    if(sendBuffers[ip])
    {
      requests[ip].Wait();
      delete sendBuffers[ip];
    }
  }
  delete[] sendBuffers;
  delete[] recvBuffers;
  //save the number of walkers sent
  NumWalkersSent=nsend;
  if(nsend)
  {
    nsend=NumPerNode[MyContext]-nsend;
    W.destroyWalkers(W.begin()+nsend, W.end());
  }
  //add walkers from other node
  if(newW.size())
    W.insert(W.end(),newW.begin(),newW.end());
}


}

