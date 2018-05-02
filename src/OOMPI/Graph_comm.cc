// -*- c++ -*-
//
// Copyright (c) 2002-2003 Indiana University.  All rights reserved.
// Copyright (c) 1996, 1997, 1998, 2000 University of Notre Dame.
//                         All rights reserved.
// 
// This file is part of the OOMPI software package.  For license
// information, see the LICENSE file in the top level directory of the
// OOMPI source distribution.
//
// $Id$
//
// OOMPI Class library
// Graph communicators
// 

#include <mpi.h>
#include "Comm.h"
#include "Intra_comm.h"
#include "Graph_comm.h"

#if 0

void OOMPI_Graph_comm::do_full_init(MPI_Comm mpi_comm, bool needs_to_be_freed)
{
  if (mpi_comm == MPI_COMM_NULL)
    return;

  int status(0);
  MPI_Topo_test(mpi_comm,&status);

  if (status == MPI_GRAPH)
    MPI_Constructor(mpi_comm, needs_to_be_freed);
  else
    OOMPI_ERROR.Handler(mpi_comm, OOMPI_ERR_COMM);
}


OOMPI_Graph_comm::OOMPI_Graph_comm(const OOMPI_Graph_comm& a)
: OOMPI_Intra_comm(a)
{
}


OOMPI_Graph_comm& 
OOMPI_Graph_comm::operator=(const OOMPI_Graph_comm& a)
{
  (OOMPI_Intra_comm&)* this =(OOMPI_Intra_comm&) a; 
  return *this;
}

OOMPI_Graph_comm::OOMPI_Graph_comm(OOMPI_Intra_comm& intra_comm, 
				   int nnodes, int index[], int edges[], 
				   bool reorder)
: OOMPI_Intra_comm(MPI_COMM_NULL)
{
#if 0
  MPI_Comm mpi_comm;
  if (MPI_Graph_create(intra_comm.comm_wrapper->Get(), nnodes, index, 
		       edges, (int) reorder, &mpi_comm) != MPI_SUCCESS)
    mpi_comm = MPI_COMM_NULL;

  MPI_Constructor(mpi_comm, true); 
#endif
}

OOMPI_Graph_comm::~OOMPI_Graph_comm()
{
}

// ---------- Communicator management

// MPI_Comm_dup
OOMPI_Graph_comm
OOMPI_Graph_comm::Dup(void)
{
  MPI_Comm mpi_comm;
  if (MPI_Comm_dup(Get_mpi(), &mpi_comm) != MPI_SUCCESS)
    return OOMPI_Graph_comm(MPI_COMM_NULL);

  return OOMPI_Graph_comm(mpi_comm, true);
}


// ---------- Process Topologies


// MPI_Graphdims_get
void 
OOMPI_Graph_comm::Dims_get(int& nnodes, int& nedges)
{  
#if 0
  MPI_Graphdims_get(comm_wrapper->Get(), &nnodes, &nedges);
#endif
  return;
}

// MPI_Graphdims_get - returns nnodes
int 
OOMPI_Graph_comm::Num_nodes(void)
{  
  int nnodes, nedges;
#if 0
  if (MPI_Graphdims_get(comm_wrapper->Get(), &nnodes, &nedges) != MPI_SUCCESS)
    return MPI_UNDEFINED;
#endif
  return nnodes;
}

// MPI_Graphdims_get - returns nedges
int 
OOMPI_Graph_comm::Num_edges(void)
{  
  int nnodes, nedges;
#if 0
  if (MPI_Graphdims_get(comm_wrapper->Get(), &nnodes, &nedges) != MPI_SUCCESS)
    return MPI_UNDEFINED;
#endif
  return nedges;
}


// MPI_Graph_get
void
OOMPI_Graph_comm::Get(int maxindex, int maxedges, int index[], int edges[])
{
  if (maxindex == 0)
    maxindex = Num_nodes();
  if (maxedges == 0)
    maxedges = Num_edges();
#if 0
  MPI_Graph_get(comm_wrapper->Get(), maxindex, maxedges, index, edges);
#endif
}

// MPI_Graph_get
void
OOMPI_Graph_comm::Get(int index[], int edges[])
{
  int maxindex, maxedges;
  Dims_get(maxindex, maxedges);
#if 0
  MPI_Graph_get(comm_wrapper->Get(), maxindex, maxedges, index, edges);
#endif
}

// MPI_Graph_get
void 
OOMPI_Graph_comm::Get_index(int index[])
{
  int maxindex, maxedges;
  Dims_get(maxindex, maxedges);
  int* edges = new int [maxedges];
#if 0
  MPI_Graph_get(comm_wrapper->Get(), maxindex, maxedges, index, edges);
#endif
  delete[] edges;
  return;
}

// MPI_Graph_get
void
OOMPI_Graph_comm::Get_edges(int edges[])
{
  int maxindex, maxedges;
  Dims_get(maxindex, maxedges);
  int* index = new int [maxindex];

#if 0
  MPI_Graph_get(comm_wrapper->Get(), maxindex, maxedges, index, edges);
#endif
  delete[] index;
  return;
}

// MPI_Graph_get
int* 
OOMPI_Graph_comm::Get_index(void)
{
  int maxindex, maxedges;
  Dims_get(maxindex, maxedges);
  int* edges = new int [maxedges];
  int* index = new int [maxindex];
  
#if 0
  if (MPI_Graph_get(comm_wrapper->Get(), maxindex, maxedges, index, edges) 
      != MPI_SUCCESS) {
    delete[] edges;
    delete[] index;
    return (int*) 0;
  }
#endif
  delete[] edges;
  return index;
}

// MPI_Graph_get
int* 
OOMPI_Graph_comm::Get_edges(void)
{
  int maxindex(0), maxedges(0);
  Dims_get(maxindex, maxedges);

  int* edges = new int [maxedges];
  int* index = new int [maxindex];
  
#if 0
  if (MPI_Graph_get(comm_wrapper->Get(), maxindex, maxedges, index, edges) 
      != MPI_SUCCESS) {
    delete[] edges;
    delete[] index;
    return (int*) 0;
  }
#endif
  delete[] index;
  return edges;
}


// MPI_Graph_neighbors
int* 
OOMPI_Graph_comm::Neighbors(int maxneighbors, int rank, int neighbors[])
{
#if 0
  if (MPI_Graph_neighbors(comm_wrapper->Get(), rank, maxneighbors, 
			  neighbors) != MPI_SUCCESS)
#endif
    return (int*) 0;
 
  return neighbors;
}

// MPI_Graph_neighbors
int* 
OOMPI_Graph_comm::Neighbors(int rank, int neighbors[])
{
  int maxneighbors = Neighbors_count();
#if 0
  if (MPI_Graph_neighbors(comm_wrapper->Get(), rank, maxneighbors, 
			  neighbors) != MPI_SUCCESS)
    return (int*) 0;
#endif
 
  return neighbors;
}

// MPI_Graph_neighbors
int* 
OOMPI_Graph_comm::Neighbors(int neighbors[])
{
  int rank = OOMPI_Intra_comm::Rank();
  int maxneighbors = Neighbors_count(rank);
  if (MPI_Graph_neighbors(comm_wrapper->Get(), rank, maxneighbors, 
			  neighbors) != MPI_SUCCESS)
    return (int*) 0;
 
  return neighbors;
}

// MPI_Graph_neighbors
int* 
OOMPI_Graph_comm::Neighbors(int maxneighbors, int rank)
{
  int* neighbors = new int [maxneighbors];
  if (MPI_Graph_neighbors(comm_wrapper->Get(), rank, maxneighbors, 
			  neighbors) != MPI_SUCCESS) {
    delete[] neighbors;
    return (int*) 0;
  }
  return neighbors;
}

// MPI_Graph_neighbors
int* 
OOMPI_Graph_comm::Neighbors(int rank)
{
  int maxneighbors = Neighbors_count(rank);
  int* neighbors = new int [maxneighbors];
  if (MPI_Graph_neighbors(comm_wrapper->Get(), rank, maxneighbors, 
			  neighbors) != MPI_SUCCESS) {
    delete[] neighbors;
    return (int*) 0;
  }
  return neighbors;
}

// MPI_Graph_neighbors
int* 
OOMPI_Graph_comm::Neighbors(void)
{
  int rank = OOMPI_Intra_comm::Rank();
  int maxneighbors = Neighbors_count(rank);
  int* neighbors = new int [maxneighbors];
  if (MPI_Graph_neighbors(comm_wrapper->Get(), rank, maxneighbors, 
			  neighbors) != MPI_SUCCESS) {
    delete[] neighbors;
    return (int*) 0;
  }
  return neighbors;
}



// MPI_Graph_neighbors_count
int
OOMPI_Graph_comm::Neighbors_count(int rank)
{
  int nneighbors(0);
  if (MPI_Graph_neighbors_count(comm_wrapper->Get(), rank, &nneighbors) !=
      MPI_SUCCESS)
    return OOMPI_UNDEFINED;

  return nneighbors;
}

// MPI_Graph_neighbors_count
int
OOMPI_Graph_comm::Neighbors_count(void)
{
  int rank = OOMPI_Intra_comm::Rank();
  int nneighbors(0);
  if (MPI_Graph_neighbors_count(comm_wrapper->Get(), rank, &nneighbors) !=
      MPI_SUCCESS)
    return OOMPI_UNDEFINED;

  return nneighbors;
}


//
// Virtual function from OOMPI_Comm
//
bool
OOMPI_Graph_comm::Test_inter()
{
  return OOMPI_Intra_comm::Test_inter();
}
#endif

