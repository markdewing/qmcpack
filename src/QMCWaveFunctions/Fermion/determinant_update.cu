#define DET_BLOCK_SIZE 8

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/time.h>

#include "determinant_update.h"

template<typename T, int BS>
__global__ void
update_inverse_cuda1 (updateJob *updateList,
		      int N, int rowstride)
{
  __shared__ T *A, *Ainv, *u, *Ainv_delta, *Ainv_colk;
  __shared__ int k;
  if (threadIdx.x==0) {
    updateJob job = updateList[blockIdx.y];
    A           = (T*)job.A;
    Ainv        = (T*)job.Ainv;
    u           = (T*)job.newRow;
    Ainv_delta  = (T*)job.AinvDelta;
    Ainv_colk   = (T*)job.AinvColk;
    k           = job.iat;
  }
  __syncthreads();

  // Store the product Ainv * u in shared memory
  T Ainv_delta_tid;
  __shared__ T 
    Ainv_colk_shared[BS], delta[BS];
  Ainv_delta_tid = 0.0f;
 __syncthreads();
  int col = blockIdx.x*BS + threadIdx.x;
  int numblocks = (N+BS-1) / BS;
  int kBlock = k/BS;

  // If the column I need to pull from Ainv is in this thread block
  // domain, do the following
  __syncthreads();
  for (int block=0; block<numblocks; block++) {
    delta[threadIdx.x] = u[block*BS+threadIdx.x] - 
      A[k*rowstride + block*BS + threadIdx.x];
    __syncthreads();
    
    
    int istop = min(BS, N-block*BS);
    for (int i=0; i<istop; i++) {
      int row = block*BS + i;
      T a = Ainv[row*rowstride+col];
      if (col == k)
	Ainv_colk_shared[i] = a;
      Ainv_delta_tid += a*delta[i];
      __syncthreads();
    }
    if (blockIdx.x == kBlock) 
      if (block*BS+threadIdx.x < N)
	Ainv_colk[block*BS+threadIdx.x] = Ainv_colk_shared[threadIdx.x];
    __syncthreads();
  }
  
  // Write the data back to global memory
  if (col < N)
    Ainv_delta[col]    = Ainv_delta_tid;
 __syncthreads();
}


template<typename T, int BS>
__global__ void
update_inverse_cuda2 (updateJob *updateList,
		      int N, int rowstride)
{
  __shared__ T *A, *u, *Ainv, *Ainv_delta, *Ainv_colk;
  int tid = threadIdx.x;
  __shared__ int k;
  if (threadIdx.x==0) {
    updateJob job = updateList[blockIdx.y];
    A          = (T*)job.A;
    u          = (T*)job.newRow;
    Ainv       = (T*)job.Ainv;
    Ainv_delta = (T*)job.AinvDelta;
    Ainv_colk  = (T*)job.AinvColk;
    k          = job.iat;
  }
  __syncthreads();
  
  T Ainv_delta_tid;
  __shared__ T  Ainv_colk_shared[BS];
  int col = blockIdx.x*BS + threadIdx.x;
  // Read the data back from global memory
  Ainv_delta_tid = Ainv_delta[col];
  Ainv_colk_shared[threadIdx.x] = Ainv_colk[col];
  if (col < N)
    A[k*rowstride + col] = u[col];
  __syncthreads();
  
  __shared__ T prefact;
  if (threadIdx.x == 0)
    prefact = -1.0f/(1.0f+Ainv_delta[k]);
		   
  int numblocks = N / BS + ((N % BS) ? 1 : 0);
  __syncthreads();
  for (int block=0; block<numblocks; block++) {
    Ainv_colk_shared[tid] = 
      prefact*Ainv_colk[block*BS+threadIdx.x];
    __syncthreads();
    T *Ainv_row = Ainv+block*BS*rowstride + col;
    int istop = min (BS, N-block*BS);
    if (col < N)
      for (int i=0; i<istop; i++, Ainv_row+=rowstride) 
	*Ainv_row += Ainv_delta_tid*Ainv_colk_shared[i];
    __syncthreads();
  }
}


void
update_inverse_cuda(updateJob updateList[], float dummy,
		    int N, int rowstride, int numWalkers)
{
  const int BS1 = 64;
  const int BS2 = 64;
  int NB1 = (N+BS1-1)/BS1;
  int NB2 = (N+BS2-1)/BS2;

  dim3 dimBlock1(BS1);
  dim3 dimGrid1(NB1, numWalkers);
  dim3 dimBlock2(BS2);
  dim3 dimGrid2(NB2, numWalkers);

  update_inverse_cuda1<float,BS1><<<dimGrid1,dimBlock1>>>
    (updateList, N, rowstride);
  update_inverse_cuda2<float,BS2><<<dimGrid2,dimBlock2>>>
    (updateList, N, rowstride);

  cudaThreadSynchronize();

  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in update_inverse_cuda:\n  %s\n",
  	     cudaGetErrorString(err));
    abort();
  }
}


void
update_inverse_cuda(updateJob updateList[], double dummy,
		    int N, int rowstride, int numWalkers)
{
  const int BS1 = 32;
  const int BS2 = 32;

  int NB1 = N/BS1 + ((N%BS1) ? 1 : 0);
  int NB2 = N/BS2 + ((N%BS2) ? 1 : 0);

  dim3 dimBlock1(BS1);
  dim3 dimGrid1(NB1, numWalkers);
  dim3 dimBlock2(BS2);
  dim3 dimGrid2(NB2, numWalkers);

  update_inverse_cuda1<double,BS1><<<dimGrid1,dimBlock1>>>
    (updateList, N, rowstride);
  update_inverse_cuda2<double,BS2><<<dimGrid2,dimBlock2>>>
    (updateList, N, rowstride);

  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in update_inverse_cuda:\n  %s\n",
  	     cudaGetErrorString(err));
    abort();
  }
}



/////////////////////////////////////////////////
// New version with fewer PCI transfers needed //
/////////////////////////////////////////////////

template<typename T, int BS>
__global__ void
update_inverse_kernel1 (T **data, int *iat, int A_off, int Ainv_off, 
			int newRow_off, int AinvDelta_off, int AinvColk_off,
			int N, int rowstride)
{
  __shared__ T *sdata;
  T *A, *Ainv, *u, *Ainv_delta, *Ainv_colk;
  __shared__ int k;
  if (threadIdx.x==0) {
    sdata = data[blockIdx.y];
    k = iat[blockIdx.y];
  }
  __syncthreads();
  A           = sdata + A_off;
  Ainv        = sdata + Ainv_off;
  u           = sdata + newRow_off;
  Ainv_delta  = sdata + AinvDelta_off;
  Ainv_colk   = sdata + AinvColk_off;

  // Store the product Ainv * u in shared memory
  __shared__ T Ainv_delta_shared[BS], 
    Ainv_colk_shared[BS], delta[BS];
  Ainv_delta_shared[threadIdx.x] = 0.0f;
 __syncthreads();
  int col = blockIdx.x*BS + threadIdx.x;
  int numblocks = N / BS + ((N%BS) ? 1 : 0);
  int kBlock = k/BS;

  // If the column I need to pull from Ainv is in this thread block
  // domain, do the following
  __syncthreads();
  for (int block=0; block<numblocks; block++) {
    delta[threadIdx.x] = u[block*BS+threadIdx.x] - 
      A[k*rowstride + block*BS + threadIdx.x];
    __syncthreads();
    
    
    int istop = min(BS, N-block*BS);
    for (int i=0; i<istop; i++) {
      int row = block*BS + i;
      T a = Ainv[row*rowstride+col];
      if (col == k)
	Ainv_colk_shared[i] = a;
      Ainv_delta_shared[threadIdx.x] += a*delta[i];
      __syncthreads();
    }
    if (blockIdx.x == kBlock) 
      if (block*BS+threadIdx.x < N)
	Ainv_colk[block*BS+threadIdx.x] = Ainv_colk_shared[threadIdx.x];
    __syncthreads();
  }
  
  // Write the data back to global memory
  if (col < N)
    Ainv_delta[col]    = Ainv_delta_shared[threadIdx.x];
 __syncthreads();
}


template<typename T, int BS>
__global__ void
update_inverse_kernel2 (T **data, int *iat, int A_off, int Ainv_off, 
			int newRow_off,	int AinvDelta_off, int AinvColk_off,
			int N, int rowstride)
			
{
  __shared__ T *sdata;
  T *A, *u, *Ainv, *Ainv_delta, *Ainv_colk;
  __shared__ int k;
  int tid = threadIdx.x;
  if (tid == 0) {
    sdata = data[blockIdx.y];
    k = iat[blockIdx.y];
  }
  __syncthreads();
  A           = sdata + A_off;
  Ainv        = sdata + Ainv_off;
  u           = sdata + newRow_off;
  Ainv_delta  = sdata + AinvDelta_off;
  Ainv_colk   = sdata + AinvColk_off;
  
  __shared__ T Ainv_delta_shared[BS];
  __shared__ T  Ainv_colk_shared[BS];
  int col = blockIdx.x*BS + threadIdx.x;
  // Read the data back from global memory
  Ainv_delta_shared[threadIdx.x] = Ainv_delta[col];
  Ainv_colk_shared[threadIdx.x] = Ainv_colk[col];
  if (col < N)
    A[k*rowstride + col] = u[col];
  __syncthreads();
  
  __shared__ T prefact;
  if (threadIdx.x == 0)
    prefact = -1.0f/(1.0f+Ainv_delta[k]);
		   
  int numblocks = N / BS + ((N % BS) ? 1 : 0);
  __syncthreads();
  for (int block=0; block<numblocks; block++) {
    Ainv_colk_shared[tid] = 
      prefact*Ainv_colk[block*BS+threadIdx.x];
    __syncthreads();
    T *Ainv_row = Ainv+block*BS*rowstride + col;
    int istop = min (BS, N-block*BS);
    if (col < N)
      for (int i=0; i<istop; i++, Ainv_row+=rowstride) 
	*Ainv_row += Ainv_delta_shared[tid]*Ainv_colk_shared[i];
    __syncthreads();
  }
}


void
update_inverse_cuda(float **data, int iat[], 
		    int A_off, int Ainv_off, int newRow_off,
		    int AinvDelta_off, int AinvColk_off,
		    int N, int rowstride, int numWalkers)
{
  const int BS1 = 64;
  const int BS2 = 64;
  int NB1 = (N+BS1-1)/BS1;
  int NB2 = (N+BS2-1)/BS2;

  dim3 dimBlock1(BS1);
  dim3 dimGrid1(NB1, numWalkers);
  dim3 dimBlock2(BS2);
  dim3 dimGrid2(NB2, numWalkers);

  update_inverse_kernel1<float,BS1><<<dimGrid1,dimBlock1>>>
    (data, iat, A_off, Ainv_off, newRow_off, AinvDelta_off, AinvColk_off,
     N, rowstride);
  update_inverse_kernel2<float,BS2><<<dimGrid2,dimBlock2>>>
    (data, iat, A_off, Ainv_off, newRow_off, AinvDelta_off, AinvColk_off,
     N, rowstride);

  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in update_inverse_cuda:\n  %s\n",
  	     cudaGetErrorString(err));
    abort();
  }
}

void
update_inverse_cuda(double **data, int iat[], 
		    int A_off, int Ainv_off, int newRow_off,
		    int AinvDelta_off, int AinvColk_off,
		    int N, int rowstride, int numWalkers)
{
  const int BS1 = 32;
  const int BS2 = 32;
  int NB1 = (N+BS1-1)/BS1;
  int NB2 = (N+BS2-1)/BS2;

  dim3 dimBlock1(BS1);
  dim3 dimGrid1(NB1, numWalkers);
  dim3 dimBlock2(BS2);
  dim3 dimGrid2(NB2, numWalkers);

  update_inverse_kernel1<double,BS1><<<dimGrid1,dimBlock1>>>
    (data, iat, A_off, Ainv_off, newRow_off, AinvDelta_off, AinvColk_off,
     N, rowstride);
  update_inverse_kernel2<double,BS2><<<dimGrid2,dimBlock2>>>
    (data, iat, A_off, Ainv_off, newRow_off, AinvDelta_off, AinvColk_off,
     N, rowstride);

  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in update_inverse_cuda:\n  %s\n",
  	     cudaGetErrorString(err));
    abort();
  }
}






template<typename T, int BS>
__global__ void
update_inverse_kernel1 (T **data, int k, int A_off, int Ainv_off, 
			int newRow_off, int AinvDelta_off, int AinvColk_off,
			int N, int rowstride)
{
  __shared__ T *sdata;
  T *A, *Ainv, *u, *Ainv_delta, *Ainv_colk;
  if (threadIdx.x==0) 
    sdata = data[blockIdx.y];
  __syncthreads();
  A           = sdata + A_off;
  Ainv        = sdata + Ainv_off;
  u           = sdata + newRow_off;
  Ainv_delta  = sdata + AinvDelta_off;
  Ainv_colk   = sdata + AinvColk_off;

  // Store the product Ainv * u in shared memory
  __shared__ T Ainv_delta_shared[BS], 
    Ainv_colk_shared[BS], delta[BS];
  Ainv_delta_shared[threadIdx.x] = 0.0f;
 __syncthreads();
  int col = blockIdx.x*BS + threadIdx.x;
  int numblocks = N / BS + ((N%BS) ? 1 : 0);
  int kBlock = k/BS;

  // If the column I need to pull from Ainv is in this thread block
  // domain, do the following
  __syncthreads();
  for (int block=0; block<numblocks; block++) {
    delta[threadIdx.x] = u[block*BS+threadIdx.x] - 
      A[k*rowstride + block*BS + threadIdx.x];
    __syncthreads();
    
    
    int istop = min(BS, N-block*BS);
    for (int i=0; i<istop; i++) {
      int row = block*BS + i;
      T a = Ainv[row*rowstride+col];
      if (col == k)
	Ainv_colk_shared[i] = a;
      Ainv_delta_shared[threadIdx.x] += a*delta[i];
      __syncthreads();
    }
    if (blockIdx.x == kBlock) 
      if (block*BS+threadIdx.x < N)
	Ainv_colk[block*BS+threadIdx.x] = Ainv_colk_shared[threadIdx.x];
    __syncthreads();
  }
  
  // Write the data back to global memory
  if (col < N)
    Ainv_delta[col]    = Ainv_delta_shared[threadIdx.x];
 __syncthreads();
}


template<typename T, int BS>
__global__ void
update_inverse_kernel2 (T **data, int k, int A_off, int Ainv_off, 
			int newRow_off,	int AinvDelta_off, int AinvColk_off,
			int N, int rowstride)
			
{
  __shared__ T *sdata;
  T *A, *u, *Ainv, *Ainv_delta, *Ainv_colk;
  int tid = threadIdx.x;
  if (tid == 0) 
    sdata = data[blockIdx.y];
  __syncthreads();
  A           = sdata + A_off;
  Ainv        = sdata + Ainv_off;
  u           = sdata + newRow_off;
  Ainv_delta  = sdata + AinvDelta_off;
  Ainv_colk   = sdata + AinvColk_off;
  
  __shared__ T Ainv_delta_shared[BS];
  __shared__ T  Ainv_colk_shared[BS];
  int col = blockIdx.x*BS + threadIdx.x;
  // Read the data back from global memory
  Ainv_delta_shared[threadIdx.x] = Ainv_delta[col];
  Ainv_colk_shared[threadIdx.x] = Ainv_colk[col];
  if (col < N)
    A[k*rowstride + col] = u[col];
  __syncthreads();
  
  __shared__ T prefact;
  if (threadIdx.x == 0)
    prefact = -1.0f/(1.0f+Ainv_delta[k]);
		   
  int numblocks = N / BS + ((N % BS) ? 1 : 0);
  __syncthreads();
  for (int block=0; block<numblocks; block++) {
    Ainv_colk_shared[tid] = 
      prefact*Ainv_colk[block*BS+threadIdx.x];
    __syncthreads();
    T *Ainv_row = Ainv+block*BS*rowstride + col;
    int istop = min (BS, N-block*BS);
    if (col < N)
      for (int i=0; i<istop; i++, Ainv_row+=rowstride) 
	*Ainv_row += Ainv_delta_shared[tid]*Ainv_colk_shared[i];
    __syncthreads();
  }
}


void
update_inverse_cuda(float **data, int iat, 
		    int A_off, int Ainv_off, int newRow_off,
		    int AinvDelta_off, int AinvColk_off,
		    int N, int rowstride, int numWalkers)
{
  const int BS1 = 64;
  const int BS2 = 64;
  int NB1 = (N+BS1-1)/BS1;
  int NB2 = (N+BS2-1)/BS2;

  dim3 dimBlock1(BS1);
  dim3 dimGrid1(NB1, numWalkers);
  dim3 dimBlock2(BS2);
  dim3 dimGrid2(NB2, numWalkers);

  update_inverse_kernel1<float,BS1><<<dimGrid1,dimBlock1>>>
    (data, iat, A_off, Ainv_off, newRow_off, AinvDelta_off, AinvColk_off,
     N, rowstride);
  update_inverse_kernel2<float,BS2><<<dimGrid2,dimBlock2>>>
    (data, iat, A_off, Ainv_off, newRow_off, AinvDelta_off, AinvColk_off,
     N, rowstride);

  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in update_inverse_cuda:\n  %s\n",
  	     cudaGetErrorString(err));
    abort();
  }
}

void
update_inverse_cuda(double **data, int iat, 
		    int A_off, int Ainv_off, int newRow_off,
		    int AinvDelta_off, int AinvColk_off,
		    int N, int rowstride, int numWalkers)
{
  const int BS1 = 32;
  const int BS2 = 32;
  int NB1 = (N+BS1-1)/BS1;
  int NB2 = (N+BS2-1)/BS2;

  dim3 dimBlock1(BS1);
  dim3 dimGrid1(NB1, numWalkers);
  dim3 dimBlock2(BS2);
  dim3 dimGrid2(NB2, numWalkers);

  update_inverse_kernel1<double,BS1><<<dimGrid1,dimBlock1>>>
    (data, iat, A_off, Ainv_off, newRow_off, AinvDelta_off, AinvColk_off,
     N, rowstride);
  update_inverse_kernel2<double,BS2><<<dimGrid2,dimBlock2>>>
    (data, iat, A_off, Ainv_off, newRow_off, AinvDelta_off, AinvColk_off,
     N, rowstride);

  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in update_inverse_cuda:\n  %s\n",
  	     cudaGetErrorString(err));
    abort();
  }
}








// The first kernel just computes AinvT * u and also stores the kth
// col of Ainv in global memory
template<typename T, int BS>
__global__ void
update_inverse_cuda1 (T **A_g, T **Ainv_g, T **u_g, 
		      T **Ainv_delta_g, T **Ainv_colk_g, 
		      int N, int rowstride, int k)
{
  __shared__ T *A, *Ainv, *u, *Ainv_delta, *Ainv_colk;
  if (threadIdx.x==0) {
    A           = A_g[blockIdx.y];
    Ainv        = Ainv_g[blockIdx.y];
    u           = u_g[blockIdx.y];
    Ainv_delta  = Ainv_delta_g[blockIdx.y];
    Ainv_colk   = Ainv_colk_g[blockIdx.y];
  }
  __syncthreads();

  // Store the product Ainv * u in shared memory
  T Ainv_delta_tid;
  __shared__ T 
    Ainv_colk_shared[BS], delta[BS];
  Ainv_delta_tid = 0.0f;
 __syncthreads();
  int col = blockIdx.x*BS + threadIdx.x;
  int numblocks = N / BS + ((N%BS) ? 1 : 0);
  int kBlock = k/BS;

  // If the column I need to pull from Ainv is in this thread block
  // domain, do the following
  __syncthreads();
  for (int block=0; block<numblocks; block++) {
    delta[threadIdx.x] = u[block*BS+threadIdx.x] - 
      A[k*rowstride + block*BS + threadIdx.x];
    __syncthreads();
    
    for (int i=0; i<BS; i++) {
      int row = block*BS + i;
      T a = Ainv[row*rowstride+col];
      if (col == k)
	Ainv_colk_shared[i] = a;
      Ainv_delta_tid += a*delta[i];
      __syncthreads();
    }
    if (blockIdx.x == kBlock) 
      Ainv_colk[block*BS+threadIdx.x] = Ainv_colk_shared[threadIdx.x];
    __syncthreads();
  }

  __syncthreads();
  
  // Write the data back to global memory
  Ainv_delta[col]    = Ainv_delta_tid;
 __syncthreads();
}


template<typename T, int BS>
__global__ void
update_inverse_cuda2 (T **A_g, T **Ainv_g, T **u_g,
		      T **Ainv_delta_g, T **Ainv_colk_g, 
		      int N, int rowstride, int k)
{
  __shared__ T *A, *u, *Ainv, *Ainv_delta, *Ainv_colk;
  int tid = threadIdx.x;
  if (threadIdx.x==0) {
    A          = A_g[blockIdx.y];
    u          = u_g[blockIdx.y];
    Ainv       = Ainv_g[blockIdx.y];
    Ainv_delta = Ainv_delta_g[blockIdx.y];
    Ainv_colk  = Ainv_colk_g[blockIdx.y];
  }
  __syncthreads();
  
  T Ainv_delta_tid;
  __shared__ T  Ainv_colk_shared[BS];
  int col = blockIdx.x*BS + threadIdx.x;
  // Read the data back from global memory
  Ainv_delta_tid = Ainv_delta[col];
  Ainv_colk_shared[threadIdx.x] = Ainv_colk[col];
  if (col < N)
    A[k*rowstride + col] = u[col];
  __syncthreads();
  
  __shared__ T prefact;
  if (threadIdx.x == 0)
    prefact = -1.0f/(1.0f+Ainv_delta[k]);
		   
  int numblocks = N / BS + ((N % BS) ? 1 : 0);
  __syncthreads();
  for (int block=0; block<numblocks; block++) {
    Ainv_colk_shared[tid] = 
      prefact*Ainv_colk[block*BS+threadIdx.x];
    __syncthreads();
    T *Ainv_row = Ainv+block*BS*rowstride + col;
    for (int i=0; i<BS; i++, Ainv_row+=rowstride) 
      *Ainv_row += Ainv_delta_tid*Ainv_colk_shared[i];
    __syncthreads();
  }
}


void
update_inverse_cuda(float *A_g[], float *Ainv_g[], float *u_g[], 
		    float *Ainv_delta_g[], float *Ainv_colk_g[], 
		    int N, int rowstride, int iat, int numWalkers)
{
  const int BS1 = 64;
  const int BS2 = 64;

  int NB1 = N/BS1 + ((N%BS1) ? 1 : 0);
  int NB2 = N/BS2 + ((N%BS2) ? 1 : 0);

  dim3 dimBlock1(BS1);
  dim3 dimGrid1(NB1, numWalkers);
  dim3 dimBlock2(BS2);
  dim3 dimGrid2(NB2, numWalkers);

  update_inverse_cuda1<float,BS1><<<dimGrid1,dimBlock1>>>
    (A_g, Ainv_g, u_g, Ainv_delta_g, Ainv_colk_g, N, rowstride, iat);
  update_inverse_cuda2<float,BS2><<<dimGrid2,dimBlock2>>>
    (A_g, Ainv_g, u_g, Ainv_delta_g, Ainv_colk_g, N, rowstride, iat);

  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in update_inverse_cuda:\n  %s\n",
  	     cudaGetErrorString(err));
    abort();
  }
}

void
update_inverse_cuda(double *A_g[], double *Ainv_g[], double *u_g[], 
		    double *Ainv_delta_g[], double *Ainv_colk_g[], 
		    int N, int rowstride, int iat, int numWalkers)
{
  const int BS1 = 32;
  const int BS2 = 32;

  int NB1 = N/BS1 + ((N%BS1) ? 1 : 0);
  int NB2 = N/BS2 + ((N%BS2) ? 1 : 0);

  dim3 dimBlock1(BS1);
  dim3 dimGrid1(NB1, numWalkers);
  dim3 dimBlock2(BS2);
  dim3 dimGrid2(NB2, numWalkers);

  update_inverse_cuda1<double,BS1><<<dimGrid1,dimBlock1>>>
    (A_g, Ainv_g, u_g, Ainv_delta_g, Ainv_colk_g, N, rowstride, iat);
  update_inverse_cuda2<double,BS2><<<dimGrid2,dimBlock2>>>
    (A_g, Ainv_g, u_g, Ainv_delta_g, Ainv_colk_g, N, rowstride, iat);

  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in update_inverse_cuda:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }

}


template<typename T, int BS, int MAXN> 
__global__ void
update_inverse_transpose_cuda(T **A_g, T **AinvT_g, T **u_g, 
			      int N, int row_stride, int elec)
{
  __shared__ float AinvT_row[MAXN], Ainv_colk[MAXN], delta[MAXN];
  int numBlocks = N/blockDim.x + ((N%blockDim.x) ? 1 : 0);
  //int numBlocks = 4;
  __shared__ T *A, *AinvT, *u;
  
  if (threadIdx.x == 0) {
    A     = A_g[blockIdx.x];
    AinvT = AinvT_g[blockIdx.x];
    u     = u_g[blockIdx.x];
  }
  __syncthreads();

  T prefactor;
  __shared__ T sum[BS];
  sum[threadIdx.x] = 0.0f;
  for (int block=0; block<numBlocks; block++) {
    int off = block *BS + threadIdx.x;
    T u1 = u[off];
    delta[off] = u1 - A[elec*row_stride + off];
    Ainv_colk[off] = AinvT[elec*row_stride + off];
    A[elec*row_stride + off] = u1;
    sum[threadIdx.x] += (Ainv_colk[off] * delta[off]);
  }
  __syncthreads();

  for (int s=(BS>>1); s>0; s>>=1) {
    __syncthreads();
    if (threadIdx.x < s)
      sum[threadIdx.x] += sum[threadIdx.x+s];
  }
  __syncthreads();

  prefactor = -1.0f/(1.0f+sum[0]);

  for (int row=0; row<N; row++) {
    // First load row into shared memory
    sum[threadIdx.x] = 0.0; 
    for (int block=0; block<numBlocks; block++) {
      int off = block*BS + threadIdx.x;
      AinvT_row[off] = AinvT[row*row_stride+off];
      sum[threadIdx.x] += (AinvT_row[off] * delta[off]);
    }
    // Now sum across row to get Ainv_delta
    for (int s=BS>>1; s>0; s>>=1) {
      __syncthreads();
      if (threadIdx.x < s)
    	sum[threadIdx.x] += sum[threadIdx.x+s];
    }
    __syncthreads();
    // sum[0] now has the AinvT * delta

    // Add on outer product
    for (int block=0; block<numBlocks; block++) {
      int off = BS*block + threadIdx.x;
      AinvT[row*row_stride + off] = AinvT_row[off] + 
	prefactor*sum[0] *Ainv_colk[off];
    }
    __syncthreads();
  }

}


template<typename T, int BS, int MAXN> 
__global__ void
update_inverse_transpose_cuda_2pass(T **A_g, T **AinvT_g, T **u_g, 
				   int N, int row_stride, int elec)
{
  __shared__ float Ainv_colk[MAXN], delta[MAXN];
  int numBlocks = N/blockDim.x + ((N%blockDim.x) ? 1 : 0);
  //int numBlocks = 4;
  __shared__ T *A, *AinvT, *u;
  
  if (threadIdx.x == 0) {
    A     = A_g[blockIdx.x];
    AinvT = AinvT_g[blockIdx.x];
    u     = u_g[blockIdx.x];
  }
  __syncthreads();

  T prefactor;
  __shared__ T sum[BS];
  sum[threadIdx.x] = 0.0f;
  for (int block=0; block<numBlocks; block++) {
    int off = block *BS + threadIdx.x;
    T u1 = u[off];
    delta[off] = u1 - A[elec*row_stride + off];
    Ainv_colk[off] = AinvT[elec*row_stride + off];
    A[elec*row_stride + off] = u1;
    sum[threadIdx.x] += (Ainv_colk[off] * delta[off]);
  }
  __syncthreads();

  for (int s=(BS>>1); s>0; s>>=1) {
    __syncthreads();
    if (threadIdx.x < s && threadIdx.y == 0)
      sum[threadIdx.x] += sum[threadIdx.x+s];
  }
  __syncthreads();

  prefactor = -1.0f/(1.0f+sum[0]);
  
  __shared__ T sum2[BS][BS+1];

  for (int b1=0; b1 < numBlocks; b1++) {
    sum[threadIdx.x] = 0.0f;
    for (int i=0; i<BS; i++)
      sum2[i][threadIdx.x] = 0.0f;
    // Compute Ainv * delta;
    for (int i=0; i<BS; i++) {
      int row = b1*BS +i;
      for (int b2=0; b2 < numBlocks; b2++) {
	int col = b2*BS + threadIdx.x;
	sum2[i][threadIdx.x] += AinvT[row*row_stride + col] * delta[col];
      }
    }
    __syncthreads();
    for (int i=0; i<BS; i++) 
      sum[threadIdx.x] += prefactor*sum2[threadIdx.x][i];
    
    // Outer product
    for (int i=0; i<BS; i++) {
      int row = b1*BS +i;
      for (int b2=0; b2 < numBlocks; b2++) {
	int col = b2*BS + threadIdx.x;
	AinvT[row*row_stride + col] += sum[i] *Ainv_colk[col];
      }
    }
    
  }
}






template<typename T, int BS>
__global__ void
calc_ratios_transpose (T **AinvT_list, T **new_row_list, 
		       T *ratio_out, int N, int row_stride, int elec,
		       int numMats)
{  
  __shared__ float *AinvT[BS], *new_row[BS];
  int matNum = blockIdx.x*BS + threadIdx.x;

  if (matNum < numMats) {
    AinvT[threadIdx.x]  = AinvT_list[matNum] + row_stride * BS;
    new_row[threadIdx.x] = new_row_list[matNum];
  }
  
  __shared__ float AinvT_phi[BS][BS+1];
  __shared__ float ratio[BS];

  ratio[threadIdx.x] = 0.0;

  int numBlocks = N / BS;
  if (numBlocks*BS < N)
    numBlocks++;
  
  for (int block=0; block<numBlocks; block++) {
    int col = block*BS + threadIdx.x;
    // First, read the data into shared memory
    for (int i=0; i<BS; i++)
      AinvT_phi[i][threadIdx.x] = (AinvT[i])[col] * (new_row[i])[col];
    __syncthreads();
    // Now sum
    for (int i=0; i<BS; i++)
      ratio[threadIdx.x] += AinvT_phi[threadIdx.x][i];
  }
  if (matNum < numMats)
    ratio_out[matNum] = ratio[threadIdx.x];
}



template<typename T, int BS>
__global__ void
calc_ratios (T **Ainv_list, T **new_row_list, 
	     T *ratio, int N, int row_stride, int elec)
{
  int tid = threadIdx.x;

  int col = /*blockIdx.x*BS * */tid;
  __shared__ T *Ainv, *new_row;

  if (tid == 0) {
    Ainv = Ainv_list[blockIdx.x];
    new_row = new_row_list[blockIdx.x];
  }
  __syncthreads();
  __shared__ T new_row_shared[BS];
   
  if (col < N) 
    new_row_shared[tid] = new_row[tid];
    
  __shared__ T Ainv_colk_shared[BS];
  // This is *highly* uncoallesced, but we just have to eat it to allow
  // other kernels to operate quickly.
  if (col < N)
    Ainv_colk_shared[tid] = Ainv[col*row_stride + elec];
  __syncthreads();

  __shared__ T Ainv_new_row[BS];
  if (col < N)
    Ainv_new_row[tid] = Ainv_colk_shared[tid] * new_row_shared[tid];
    
  __syncthreads();
    // Now, we have to dot
  for (unsigned int s=BS/2; s>0; s>>=1) {
    if (tid < s && (tid+s) < N)
      Ainv_new_row[tid] += Ainv_new_row[tid + s];
    __syncthreads();
  }
  if (tid == 0)      ratio[blockIdx.x] = Ainv_new_row[0];
}


void
determinant_ratios_cuda (float *Ainv_list[], float *new_row_list[],
			 float *ratios, int N, int row_stride, int iat,
			 int numWalkers)
{
  dim3 dimBlock(N);
  dim3 dimGrid(numWalkers);

  cudaThreadSynchronize();
  cudaError_t err1 = cudaGetLastError();
  if (err1 != cudaSuccess) {
    fprintf (stderr, "CUDA error before determinant_ratios_cuda:\n  %s\n",
	     cudaGetErrorString(err1));
    abort();
  }

  if (N <= 32) 
    calc_ratios<float,32><<<dimGrid,dimBlock>>>(Ainv_list, new_row_list, ratios, N, row_stride, iat);
  else if (N <= 64)
    calc_ratios<float,64><<<dimGrid,dimBlock>>>(Ainv_list, new_row_list, ratios, N, row_stride, iat);
  else if (N <= 128)
    calc_ratios<float,128><<<dimGrid,dimBlock>>>(Ainv_list, new_row_list, ratios, N, row_stride, iat);
  else if (N <= 256)
    calc_ratios<float,256><<<dimGrid,dimBlock>>>(Ainv_list, new_row_list, ratios, N, row_stride, iat);
  else if (N <= 512)
    calc_ratios<float,512><<<dimGrid,dimBlock>>>(Ainv_list, new_row_list, ratios, N, row_stride, iat);
  else if (N <= 1024)
    calc_ratios<float,1024><<<dimGrid,dimBlock>>>(Ainv_list, new_row_list, ratios, N, row_stride, iat);
  else {
    fprintf (stdout, "Error:  N too large for CUDA evaluation.\n");
    abort();
  }

  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in determinant_ratios_cuda:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }

}

void
determinant_ratios_cuda (double *Ainv_list[], double *new_row_list[],
			 double *ratios, int N, int row_stride, int iat,
			 int numWalkers)
{
  dim3 dimBlock(N);
  dim3 dimGrid(numWalkers);

  if (N <= 32) 
    calc_ratios<double,32><<<dimGrid,dimBlock>>>(Ainv_list, new_row_list, ratios, N, row_stride, iat);
  else if (N <= 64)
    calc_ratios<double,64><<<dimGrid,dimBlock>>>(Ainv_list, new_row_list, ratios, N, row_stride, iat);
  else if (N <= 128)
    calc_ratios<double,128><<<dimGrid,dimBlock>>>(Ainv_list, new_row_list, ratios, N, row_stride, iat);
  else if (N <= 256)
    calc_ratios<double,256><<<dimGrid,dimBlock>>>(Ainv_list, new_row_list, ratios, N, row_stride, iat);
  else if (N <= 512)
    calc_ratios<double,512><<<dimGrid,dimBlock>>>(Ainv_list, new_row_list, ratios, N, row_stride, iat);
  else {
    fprintf (stdout, "Error:  N too large for CUDA evaluation.\n");
    abort();
  }
  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in determinant_ratios_cuda:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }
}


template<typename T, int BS>
__global__ void
calc_ratio_grad_lapl (T **Ainv_list, T **new_row_list, T **grad_lapl_list,
		      T *ratio_grad_lapl, int N, int row_stride, int elec)
{
  int tid = threadIdx.x;
  int NB = N/BS + ((N % BS) ? 1 : 0);

  __shared__ T *Ainv, *new_row, *grad_lapl;

  if (tid == 0) {
    Ainv = Ainv_list[blockIdx.x];
    new_row = new_row_list[blockIdx.x];
    grad_lapl = grad_lapl_list[blockIdx.x];
  }
  __syncthreads();
    
  __shared__ T Ainv_colk_shared[BS];
  __shared__ T ratio_prod[5][BS+1];
  ratio_prod[0][tid] = 0.0f;
  ratio_prod[1][tid] = 0.0f;
  ratio_prod[2][tid] = 0.0f;
  ratio_prod[3][tid] = 0.0f;
  ratio_prod[4][tid] = 0.0f;
  // This is *highly* uncoallesced, but we just have to eat it to allow
  // other kernels to operate quickly.
  __syncthreads();
  for (int block=0; block<NB; block++) {
    int col = block*BS + tid;
    if (col < N)
      Ainv_colk_shared[tid] = Ainv[col*row_stride + elec];
    __syncthreads();
    if (col < N) {
      ratio_prod[0][tid] += Ainv_colk_shared[tid] * new_row[col];
      ratio_prod[1][tid] += Ainv_colk_shared[tid] * grad_lapl[0*row_stride+col];
      ratio_prod[2][tid] += Ainv_colk_shared[tid] * grad_lapl[1*row_stride+col];
      ratio_prod[3][tid] += Ainv_colk_shared[tid] * grad_lapl[2*row_stride+col];
      ratio_prod[4][tid] += Ainv_colk_shared[tid] * grad_lapl[3*row_stride+col];
    }
    __syncthreads();
  }
  // Now, we have to sum
  for (unsigned int s=BS/2; s>0; s>>=1) {
    if (tid < s) {
      ratio_prod[0][tid] += ratio_prod[0][tid + s]; // Value
      ratio_prod[1][tid] += ratio_prod[1][tid + s]; // grad_x
      ratio_prod[2][tid] += ratio_prod[2][tid + s]; // grad_y
      ratio_prod[3][tid] += ratio_prod[3][tid + s]; // grad_z
      ratio_prod[4][tid] += ratio_prod[4][tid + s]; // lapl
    }
    __syncthreads();
  }
  // Subtract off gradient^2 from laplacian
  if (tid == 0) 
    ratio_prod[4][0] -= (ratio_prod[1][0]*ratio_prod[1][0] +
			 ratio_prod[1][0]*ratio_prod[1][0] +
			 ratio_prod[1][0]*ratio_prod[1][0]);
  __syncthreads();
  // Present gradient and laplacian are w.r.t old position.  Divide by 
  // ratio to make it w.r.t. new position
  if (tid < 4)
    ratio_prod[tid+1][0] /= ratio_prod[0][0];

  if (tid < 5) 
    ratio_grad_lapl[5*blockIdx.x+tid] = ratio_prod[tid][0];

}


template<typename T, int BS>
__global__ void
calc_ratio_grad_lapl (T **Ainv_list, T **new_row_list, T **grad_lapl_list,
		      T *ratio_grad_lapl, int N, int row_stride, int *elec_list)
{
  int tid = threadIdx.x;
  int NB = N/BS + ((N % BS) ? 1 : 0);

  __shared__ T *Ainv, *new_row, *grad_lapl;
  __shared__ int elec;

  if (tid == 0) {
    Ainv = Ainv_list[blockIdx.x];
    new_row = new_row_list[blockIdx.x];
    grad_lapl = grad_lapl_list[blockIdx.x];
    elec = elec_list[blockIdx.x];
  }
  __syncthreads();
    
  __shared__ T Ainv_colk_shared[BS];
  __shared__ T ratio_prod[5][BS+1];
  ratio_prod[0][tid] = 0.0f;
  ratio_prod[1][tid] = 0.0f;
  ratio_prod[2][tid] = 0.0f;
  ratio_prod[3][tid] = 0.0f;
  ratio_prod[4][tid] = 0.0f;
  // This is *highly* uncoallesced, but we just have to eat it to allow
  // other kernels to operate quickly.
  __syncthreads();
  for (int block=0; block<NB; block++) {
    int col = block*BS + tid;
    if (col < N)
      Ainv_colk_shared[tid] = Ainv[col*row_stride + elec];
    __syncthreads();
    if (col < N) {
      ratio_prod[0][tid] += Ainv_colk_shared[tid] * new_row[col];
      ratio_prod[1][tid] += Ainv_colk_shared[tid] * grad_lapl[0*row_stride+col];
      ratio_prod[2][tid] += Ainv_colk_shared[tid] * grad_lapl[1*row_stride+col];
      ratio_prod[3][tid] += Ainv_colk_shared[tid] * grad_lapl[2*row_stride+col];
      ratio_prod[4][tid] += Ainv_colk_shared[tid] * grad_lapl[3*row_stride+col];
    }
    __syncthreads();
  }
  // Now, we have to sum
  for (unsigned int s=BS/2; s>0; s>>=1) {
    if (tid < s) {
      ratio_prod[0][tid] += ratio_prod[0][tid + s]; // Value
      ratio_prod[1][tid] += ratio_prod[1][tid + s]; // grad_x
      ratio_prod[2][tid] += ratio_prod[2][tid + s]; // grad_y
      ratio_prod[3][tid] += ratio_prod[3][tid + s]; // grad_z
      ratio_prod[4][tid] += ratio_prod[4][tid + s]; // lapl
    }
    __syncthreads();
  }
  // Subtract off gradient^2 from laplacian
  if (tid == 0) 
    ratio_prod[4][0] -= (ratio_prod[1][0]*ratio_prod[1][0] +
			 ratio_prod[1][0]*ratio_prod[1][0] +
			 ratio_prod[1][0]*ratio_prod[1][0]);
  __syncthreads();
  // Present gradient and laplacian are w.r.t old position.  Divide by 
  // ratio to make it w.r.t. new position
  if (tid < 4)
    ratio_prod[tid+1][0] /= ratio_prod[0][0];

  if (tid < 5) 
    ratio_grad_lapl[5*blockIdx.x+tid] = ratio_prod[tid][0];

}




void
determinant_ratios_grad_lapl_cuda (float *Ainv_list[], float *new_row_list[],
				   float *grad_lapl_list[], float ratios_grad_lapl[], 
				   int N, int row_stride, int iat, int numWalkers)
{
  const int BS = 32;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  
  calc_ratio_grad_lapl<float,BS><<<dimGrid,dimBlock>>>
    (Ainv_list, new_row_list, grad_lapl_list, ratios_grad_lapl, N, row_stride, iat);
}



void
determinant_ratios_grad_lapl_cuda (double *Ainv_list[], double *new_row_list[],
				   double *grad_lapl_list[], double ratios_grad_lapl[], 
				   int N, int row_stride, int iat, int numWalkers)
{
  const int BS = 32;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  
  calc_ratio_grad_lapl<double,BS><<<dimGrid,dimBlock>>>
    (Ainv_list, new_row_list, grad_lapl_list, ratios_grad_lapl, N, row_stride, iat);
}



void
determinant_ratios_grad_lapl_cuda (float *Ainv_list[], float *new_row_list[],
				   float *grad_lapl_list[], float ratios_grad_lapl[], 
				   int N, int row_stride, int iat_list[], int numWalkers)
{
  const int BS = 32;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  
  calc_ratio_grad_lapl<float,BS><<<dimGrid,dimBlock>>>
    (Ainv_list, new_row_list, grad_lapl_list, ratios_grad_lapl, N, row_stride, iat_list);
}



void
determinant_ratios_grad_lapl_cuda (double *Ainv_list[], double *new_row_list[],
				   double *grad_lapl_list[], double ratios_grad_lapl[], 
				   int N, int row_stride, int iat_list[], int numWalkers)
{
  const int BS = 32;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  
  calc_ratio_grad_lapl<double,BS><<<dimGrid,dimBlock>>>
    (Ainv_list, new_row_list, grad_lapl_list, ratios_grad_lapl, N, row_stride, iat_list);
}




template<typename T, int BS>
__global__ void
calc_grad_kernel (T **Ainv_list, T **grad_lapl_list,
		  T *grad, int N, int row_stride, int elec)
{
  int tid = threadIdx.x;
  int NB = N/BS + ((N % BS) ? 1 : 0);

  __shared__ T *Ainv, *grad_lapl;

  if (tid == 0) {
    Ainv = Ainv_list[blockIdx.x];
    grad_lapl = grad_lapl_list[blockIdx.x] + 4*elec*row_stride;
  }
  __syncthreads();
    
  __shared__ T Ainv_colk_shared[BS];
  __shared__ T ratio_prod[3][BS+1];
  ratio_prod[0][tid] = 0.0f;
  ratio_prod[1][tid] = 0.0f;
  ratio_prod[2][tid] = 0.0f;
  // This is *highly* uncoallesced, but we just have to eat it to allow
  // other kernels to operate quickly.
  __syncthreads();
  for (int block=0; block<NB; block++) {
    int col = block*BS + tid;
    if (col < N)
      Ainv_colk_shared[tid] = Ainv[col*row_stride + elec];
    __syncthreads();
    if (col < N) {
      ratio_prod[0][tid] += Ainv_colk_shared[tid] * grad_lapl[0*row_stride+col];
      ratio_prod[1][tid] += Ainv_colk_shared[tid] * grad_lapl[1*row_stride+col];
      ratio_prod[2][tid] += Ainv_colk_shared[tid] * grad_lapl[2*row_stride+col];
    }
    __syncthreads();
  }
  // Now, we have to sum
  for (unsigned int s=BS/2; s>0; s>>=1) {
    if (tid < s) {
      ratio_prod[0][tid] += ratio_prod[0][tid + s]; // grad_x
      ratio_prod[1][tid] += ratio_prod[1][tid + s]; // grad_y
      ratio_prod[2][tid] += ratio_prod[2][tid + s]; // grad_z
    }
    __syncthreads();
  }
  if (tid < 3) 
    grad[3*blockIdx.x+tid] = ratio_prod[tid][0];
}

void
calc_gradient (float *Ainv_list[], float *grad_lapl_list[],
	       float grad[], int N, int row_stride, int elec,
	       int numWalkers)
{
  const int BS = 32;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);

  calc_grad_kernel<float,BS><<<dimGrid,dimBlock>>>
    (Ainv_list, grad_lapl_list, grad, N, row_stride, elec);
}

void
calc_gradient (double *Ainv_list[], double *grad_lapl_list[],
	       double grad[], int N, int row_stride, int elec,
	       int numWalkers)
{
  const int BS = 32;
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);

  calc_grad_kernel<double,BS><<<dimGrid,dimBlock>>>
    (Ainv_list, grad_lapl_list, grad, N, row_stride, elec);
}



#define RATIO_BS 16

template<typename T>
__global__ void
all_ratios_kernel (T **Ainv_list, T **new_mat_list, 
		   T **ratio_list, int N, int row_stride)
{
  __shared__ T *Ainv, *new_mat, *ratio;
  
  if (threadIdx.x == 0 && threadIdx.y == 0) {
    Ainv    = Ainv_list[blockIdx.x];
    new_mat = new_mat_list[blockIdx.x];
    ratio   = ratio_list[blockIdx.x];
  }

  __shared__ float Ainv_block[RATIO_BS][RATIO_BS+1];
  // __shared__ float new_block[RATIO_BS][RATIO_BS+1];
  __shared__ float ratio_block[RATIO_BS][RATIO_BS+1];
  unsigned int numBlocks = N >> 4;
  if (N & 15)
    numBlocks++;

  for (unsigned int yBlock=0; yBlock<numBlocks; yBlock++) {
    ratio_block[threadIdx.y][threadIdx.x] = 0.0f;
    __syncthreads();
    for (unsigned int xBlock=0; xBlock<numBlocks; xBlock++) {
      unsigned int xIndex = yBlock * RATIO_BS + threadIdx.x;
      unsigned int yIndex = xBlock * RATIO_BS + threadIdx.y;
      unsigned int index  = yIndex*row_stride + xIndex;
      if ((xIndex < N) && (yIndex < N))
	Ainv_block[threadIdx.x][threadIdx.y] = Ainv[index];
      __syncthreads();
      xIndex = xBlock * RATIO_BS + threadIdx.x;
      yIndex = yBlock * RATIO_BS + threadIdx.y;
      index  = yIndex*row_stride + xIndex;

      if ((xIndex < N) && (yIndex < N))
	ratio_block[threadIdx.y][threadIdx.x] +=
	  new_mat[index] * Ainv_block[threadIdx.y][threadIdx.x];
      __syncthreads();
    }
    __syncthreads();
    // Now, we have to do the reduction across the ratio_blocks
    
    if (threadIdx.x < 8)
      ratio_block[threadIdx.y][threadIdx.x] +=
	ratio_block[threadIdx.y][threadIdx.x+8];
    if (threadIdx.x < 4)
      ratio_block[threadIdx.y][threadIdx.x] +=
	ratio_block[threadIdx.y][threadIdx.x+4];
    if (threadIdx.x < 2)
      ratio_block[threadIdx.y][threadIdx.x] +=
	ratio_block[threadIdx.y][threadIdx.x+2];
    if (threadIdx.x < 1) 
      ratio_block[threadIdx.y][threadIdx.x] +=
	ratio_block[threadIdx.y][threadIdx.x+1];
    __syncthreads();

    if (threadIdx.y == 0 && (yBlock * RATIO_BS + threadIdx.x) < N)
      ratio[yBlock * RATIO_BS + threadIdx.x] = ratio_block[threadIdx.x][0];
  }      
}




void
calc_all_ratios (float *Ainv_list[], float *new_mat_list[],
		 float *ratio_list[], int N, int row_stride, int num_mats)
{
  dim3 dimBlock(RATIO_BS, RATIO_BS);
  dim3 dimGrid (num_mats);

  all_ratios_kernel<float><<<dimGrid,dimBlock>>>
    (Ainv_list, new_mat_list, ratio_list, N, row_stride);

  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in cal_all_ratios:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }
}


const int MAX_RATIO_ROWS = 20;

template<typename T, int BS>
__global__ void
calc_many_ratios_kernel (T **Ainv_list, T **new_row_list,
			 T **ratio_list, int *num_ratio_list,
			 int N, int row_stride, int *elec_list)
{
  int tid = threadIdx.x;

  __shared__ T *Ainv, *new_rows, *ratios;
  __shared__ int num_ratios, elec;
  if (tid == 0) {
    Ainv = Ainv_list[blockIdx.x];
    new_rows = new_row_list[blockIdx.x];
    num_ratios = num_ratio_list[blockIdx.x];
    ratios   = ratio_list[blockIdx.x];
    elec = elec_list[blockIdx.x];
  }
  __syncthreads();

  int NB = N/BS + ((N%BS) ? 1 : 0);
  __shared__ T Ainv_shared[BS], row[BS];
  // We use BS+1 to avoid bank conflicts in the writing.
  __shared__ T ratio_sum[MAX_RATIO_ROWS][BS+1];
  for (int iratio=0; iratio<num_ratios; iratio++)
    ratio_sum[iratio][tid] = (T)0.0;
  __syncthreads();

  for (int block=0; block<NB; block++) {
    int off = block*BS+tid;
    bool mask = off < N;
    if (mask)
      Ainv_shared[tid] = Ainv[off*row_stride+elec];
    __syncthreads();
    for (int iratio=0; iratio<num_ratios; iratio++) 
      if (mask)
	ratio_sum[iratio][tid] += Ainv_shared[tid] *
	  new_rows[iratio*row_stride + off];
    __syncthreads();
  }
  // now, sum up ratios
  for (int iratio = 0; iratio<num_ratios; iratio++) {
    for (int s=BS>>1; s>0; s>>=1) {
      if (tid < s)
	ratio_sum[iratio][tid] += ratio_sum[iratio][tid+s];
      __syncthreads();
    }
  }
  // Store sums in parallel
  if (tid < num_ratios) 
    ratios[tid] = ratio_sum[tid][0];
}    
  
void
calc_many_ratios (float *Ainv_list[], float *new_row_list[],
		  float* ratio_list[], int num_ratio_list[],
		  int N, int row_stride, int elec_list[],
		  int numWalkers)
{
  const int BS=32;
  
  dim3 dimBlock(BS);
  dim3 dimGrid (numWalkers);
  
  calc_many_ratios_kernel<float,BS><<<dimGrid,dimBlock>>>
    (Ainv_list, new_row_list, ratio_list, num_ratio_list,
     N, row_stride, elec_list);
}

void
calc_many_ratios (double *Ainv_list[], double *new_row_list[],
		  double* ratio_list[], int num_ratio_list[],
		  int N, int row_stride, int elec_list[],
		  int numWalkers)
{
  const int BS=32;
  
  dim3 dimBlock(BS);
  dim3 dimGrid (numWalkers);
  
  calc_many_ratios_kernel<double,BS><<<dimGrid,dimBlock>>>
    (Ainv_list, new_row_list, ratio_list, num_ratio_list,
     N, row_stride, elec_list);

  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in calc_many_ratios:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }
}




#define SCALE_BS 64

__constant__ float GGt[3][3];


template<typename T>
__global__ void
scale_grad_lapl_kernel (T **grad_list, T **hess_list,
			T **grad_lapl_list, T *Linv, int N)
{
  __shared__ float gradBlock[3][SCALE_BS];
  __shared__ float hessBlock[6][SCALE_BS];
  //  __shared__ float outBlock [4][SCALE_BS];
  __shared__ float G[3][3], GGt[3][3];
  __shared__ float *grad, *hess, *out;
  
  if (threadIdx.x == 0) {
    grad = grad_list[blockIdx.y];
    hess = hess_list[blockIdx.y];
    out  = grad_lapl_list[blockIdx.y];
  }

  int i = threadIdx.x/3;
  int j = threadIdx.x%3;

  if (threadIdx.x < 9) 
    G[i][j] = Linv[threadIdx.x];
  __syncthreads();
  
  if (threadIdx.x < 9) {
    GGt[i][j] = (G[i][0] * G[0][j] +
		 G[i][1] * G[1][j] +
		 G[i][2] * G[2][j]);
  }
  
  // Load the gradients into shared memory
  for (int i=0; i<3; i++) {
    unsigned int gIndex = (3 * blockIdx.x+i) * SCALE_BS + threadIdx.x;
    if (gIndex < 3*N)  gradBlock[i][threadIdx.x] = grad[gIndex];
  }
  // Load the hessian into shared memory
  for (int i=0; i<6; i++) {
    unsigned int hIndex = (6 * blockIdx.x+i) * SCALE_BS + threadIdx.x;
    if (hIndex < 6*N)  hessBlock[i][threadIdx.x] = grad[hIndex];
  }

  // Now, loop through the rows that I own and compute the
  // dimensioned gradients and laplacians from the 
  // dimensionless gradients and Hessians.
  int row = blockIdx.x*SCALE_BS;

  float val;
  // x component of gradient
  val = (G[0][0]*gradBlock[0][threadIdx.x] +
	 G[0][1]*gradBlock[1][threadIdx.x] +
	 G[0][2]*gradBlock[2][threadIdx.x]);
  out[row + 0*N + threadIdx.x] = val;

  // y component of gradient
  val = (G[1][0]*gradBlock[0][threadIdx.x] +
	 G[1][1]*gradBlock[1][threadIdx.x] +
	 G[1][2]*gradBlock[2][threadIdx.x]);
  out[row + 1*N + threadIdx.x] = val;

  // z component of gradient
  val = (G[2][0]*gradBlock[0][threadIdx.x] +
	 G[2][1]*gradBlock[1][threadIdx.x] +
	 G[2][2]*gradBlock[2][threadIdx.x]);
  out[row + 2*N + threadIdx.x] = val;

  
  // Hessian = H00 H01 H02 H11 H12 H22
  // Matrix = [0 1 2]
  //          [1 3 4]
  //          [2 4 5]
  // laplacian = Trace(GGt*Hessian)
  val = (GGt[0][0]*hessBlock[0][threadIdx.x] +
	 GGt[0][1]*hessBlock[1][threadIdx.x] +
	 GGt[0][2]*hessBlock[2][threadIdx.x] +

	 GGt[1][0]*hessBlock[1][threadIdx.x] +
	 GGt[1][1]*hessBlock[3][threadIdx.x] +
	 GGt[1][2]*hessBlock[4][threadIdx.x] +

	 GGt[2][0]*hessBlock[2][threadIdx.x] +
	 GGt[2][1]*hessBlock[4][threadIdx.x] +
	 GGt[2][2]*hessBlock[5][threadIdx.x]);

  out[row + 3*N + threadIdx.x] = val;

}


// This function reads the vectors pointed to by grad_list and 
// hess_list.  These are in memory as
// [grad0_x grad0_y grad0_z grad1_x grad1_y ... ] and
// [hess0_xx hess0_xy hess0_xy hess0_yy hess0_yz hess0_zz ...]
// It the writes the data into memory as
// [grad0_x grad1_x ... grad(N-1)_x grad0_y ... grad(N-1)_x lapl0
// lapl1...]

void
scale_grad_lapl(float *grad_list[], float *hess_list[],
		float *grad_lapl_list[], float Linv[], int N,
		int num_walkers)
{
  dim3 dimBlock(SCALE_BS);
  dim3 dimGrid(N/SCALE_BS, num_walkers);
  if (N%SCALE_BS)
    dimGrid.x++;
  
  scale_grad_lapl_kernel<float><<<dimGrid,dimBlock>>>
    (grad_list, hess_list, grad_lapl_list, Linv, N);
}


template<typename T>
__global__ void
all_ratios_grad_lapl_kernel (T **Ainv_list, T **grad_lapl_list, 
			     T **out_list, int N, int row_stride)
{
  __shared__ T *Ainv, *gl_array, *out;
  
  if (threadIdx.x == 0 && threadIdx.y == 0) {
    Ainv     = Ainv_list[blockIdx.x];
    gl_array = grad_lapl_list[blockIdx.x];
    out      = out_list[blockIdx.x];
  }
  __syncthreads();

  __shared__ float Ainv_block[RATIO_BS][RATIO_BS+1];
  __shared__ float grad_lapl_block[4][RATIO_BS][RATIO_BS+1];
  unsigned int numBlocks = N >> 4;
  if (N & 15)
    numBlocks++;

  __syncthreads();
  for (unsigned int yBlock=0; yBlock<numBlocks; yBlock++) {
    __syncthreads();
    grad_lapl_block[0][threadIdx.y][threadIdx.x] = 0.0f;
    grad_lapl_block[1][threadIdx.y][threadIdx.x] = 0.0f;
    grad_lapl_block[2][threadIdx.y][threadIdx.x] = 0.0f;
    grad_lapl_block[3][threadIdx.y][threadIdx.x] = 0.0f;
    __syncthreads();
    for (unsigned int xBlock=0; xBlock<numBlocks; xBlock++) {
      unsigned int xIndex = yBlock * RATIO_BS + threadIdx.x;
      unsigned int yIndex = xBlock * RATIO_BS + threadIdx.y;
      unsigned int index  = yIndex*row_stride + xIndex;
      if ((xIndex < N) && (yIndex < N))
       	Ainv_block[threadIdx.x][threadIdx.y] = Ainv[index];
      __syncthreads();
      xIndex = xBlock * RATIO_BS + threadIdx.x;
      yIndex = yBlock * RATIO_BS + threadIdx.y;
      index  = 4*yIndex*row_stride + xIndex;
      __syncthreads();
      if ((xIndex < N) && (yIndex < N)) {
	grad_lapl_block[0][threadIdx.y][threadIdx.x] +=
	  gl_array[index+0*row_stride] * Ainv_block[threadIdx.y][threadIdx.x];
	grad_lapl_block[1][threadIdx.y][threadIdx.x] +=
	  gl_array[index+1*row_stride] * Ainv_block[threadIdx.y][threadIdx.x];
	grad_lapl_block[2][threadIdx.y][threadIdx.x] +=
	  gl_array[index+2*row_stride] * Ainv_block[threadIdx.y][threadIdx.x];
	grad_lapl_block[3][threadIdx.y][threadIdx.x] +=
	  gl_array[index+3*row_stride] * Ainv_block[threadIdx.y][threadIdx.x];
      }
      __syncthreads();
    }
    // Now, we have to do the reduction across the lapl_blocks
    
    if (threadIdx.x < 8) {
      grad_lapl_block[0][threadIdx.y][threadIdx.x] +=
    	grad_lapl_block[0][threadIdx.y][threadIdx.x+8];
      grad_lapl_block[1][threadIdx.y][threadIdx.x] +=
    	grad_lapl_block[1][threadIdx.y][threadIdx.x+8];
      grad_lapl_block[2][threadIdx.y][threadIdx.x] +=
    	grad_lapl_block[2][threadIdx.y][threadIdx.x+8];
      grad_lapl_block[3][threadIdx.y][threadIdx.x] +=
    	grad_lapl_block[3][threadIdx.y][threadIdx.x+8];
    }
    __syncthreads();
    if (threadIdx.x < 4) {
      grad_lapl_block[0][threadIdx.y][threadIdx.x] +=
    	grad_lapl_block[0][threadIdx.y][threadIdx.x+4];
      grad_lapl_block[1][threadIdx.y][threadIdx.x] +=
    	grad_lapl_block[1][threadIdx.y][threadIdx.x+4];
      grad_lapl_block[2][threadIdx.y][threadIdx.x] +=
    	grad_lapl_block[2][threadIdx.y][threadIdx.x+4];
      grad_lapl_block[3][threadIdx.y][threadIdx.x] +=
    	grad_lapl_block[3][threadIdx.y][threadIdx.x+4];
    }
    __syncthreads();
    if (threadIdx.x < 2) {
      grad_lapl_block[0][threadIdx.y][threadIdx.x] +=
    	grad_lapl_block[0][threadIdx.y][threadIdx.x+2];
      grad_lapl_block[1][threadIdx.y][threadIdx.x] +=
    	grad_lapl_block[1][threadIdx.y][threadIdx.x+2];
      grad_lapl_block[2][threadIdx.y][threadIdx.x] +=
    	grad_lapl_block[2][threadIdx.y][threadIdx.x+2];
      grad_lapl_block[3][threadIdx.y][threadIdx.x] +=
    	grad_lapl_block[3][threadIdx.y][threadIdx.x+2];
    }
    __syncthreads();
    if (threadIdx.x < 1) {
      grad_lapl_block[0][threadIdx.y][threadIdx.x] +=
    	grad_lapl_block[0][threadIdx.y][threadIdx.x+1];
      grad_lapl_block[1][threadIdx.y][threadIdx.x] +=
    	grad_lapl_block[1][threadIdx.y][threadIdx.x+1];
      grad_lapl_block[2][threadIdx.y][threadIdx.x] +=
    	grad_lapl_block[2][threadIdx.y][threadIdx.x+1];
      grad_lapl_block[3][threadIdx.y][threadIdx.x] +=
    	grad_lapl_block[3][threadIdx.y][threadIdx.x+1];
    }
    __syncthreads();

    // unsigned int yIndex = yBlock * RATIO_BS + threadIdx.x;

    // if (threadIdx.y == 0 && yIndex < N) {
    //   out[4*yIndex+0] = grad_lapl_block[0][threadIdx.x][0];
    //   out[4*yIndex+1] = grad_lapl_block[1][threadIdx.x][0];
    //   out[4*yIndex+2] = grad_lapl_block[2][threadIdx.x][0];
    //   out[4*yIndex+3] = grad_lapl_block[3][threadIdx.x][0];
    // }
    //unsigned int yIndex = 4*yBlock*RATIO_BS + 4*threadIdx.y + threadIdx.x;

    unsigned int ix = 16*threadIdx.y + threadIdx.x;
    
    unsigned int yIndex = RATIO_BS * yBlock + (ix >> 2);

    if (ix < 64 && yIndex < N) 
      out[64*yBlock + ix] = grad_lapl_block[ix&3][ix>>2][0];
    // IMPORTANT!!!
    __syncthreads();
  }      
}

void
calc_grad_lapl (float *Ainv_list[], float *grad_lapl_list[],
		float *out_list[], int N, int row_stride, int num_mats)
{
  dim3 dimBlock(RATIO_BS, RATIO_BS);
  dim3 dimGrid (num_mats);

  all_ratios_grad_lapl_kernel<float><<<dimGrid,dimBlock>>>
    (Ainv_list, grad_lapl_list, out_list, N, row_stride);

  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in cal_grad_lapl:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }
}


void
calc_grad_lapl (double *Ainv_list[], double *grad_lapl_list[],
		double *out_list[], int N, int row_stride, int num_mats)
{
  dim3 dimBlock(RATIO_BS, RATIO_BS);
  dim3 dimGrid (num_mats);

  all_ratios_grad_lapl_kernel<double><<<dimGrid,dimBlock>>>
    (Ainv_list, grad_lapl_list, out_list, N, row_stride);

  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in cal_grad_lapl:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }
}

#define COPY_BS 256

template<typename T>
__global__ void
multi_copy (T **dest, T **src, int len)
{
  __shared__ T *mysrc, *mydest;
  if (threadIdx.x ==0) {
    mysrc = src[blockIdx.y];
    mydest = dest[blockIdx.y];
  }
  __syncthreads();
  int i = blockIdx.x * COPY_BS + threadIdx.x;
  if (i < len)
    mydest[i] = mysrc[i];
}


template<typename T>
__global__ void
multi_copy (T **buff, int dest_off, int src_off, int len)
{
  __shared__ T *mysrc, *mydest;
  if (threadIdx.x ==0) {
    T* ptr = buff[blockIdx.y];
    mysrc  = ptr + src_off;
    mydest = ptr + dest_off;
  }
  __syncthreads();
  int i = blockIdx.x * COPY_BS + threadIdx.x;
  if (i < len)
    mydest[i] = mysrc[i];
}



void
multi_copy (float *dest[], float *src[], int len, int num)
{
  dim3 dimBlock(COPY_BS);
  dim3 dimGrid ((len+COPY_BS-1)/COPY_BS, num);
  
  multi_copy<float><<<dimGrid,dimBlock>>>(dest, src, len);
}

void
multi_copy (double *dest[], double *src[], int len, int num)
{
  dim3 dimBlock(COPY_BS);
  dim3 dimGrid (len/COPY_BS, num);
  if (len % COPY_BS)
    dimGrid.x++;
  
  multi_copy<double><<<dimGrid,dimBlock>>>(dest, src, len);
}



void
multi_copy (float *buff[], int dest_off, int src_off, int len, int num)
{
  dim3 dimBlock(COPY_BS);
  dim3 dimGrid ((len+COPY_BS-1)/COPY_BS, num);
  
  multi_copy<float><<<dimGrid,dimBlock>>>(buff, dest_off, src_off, len);
}

void
multi_copy (double *buff[], int dest_off, int src_off, int len, int num)
{
  dim3 dimBlock(COPY_BS);
  dim3 dimGrid (len/COPY_BS, num);
  if (len % COPY_BS)
    dimGrid.x++;
  
  multi_copy<double><<<dimGrid,dimBlock>>>(buff, dest_off, src_off, len);
}



#include <stdlib.h>
#include <time.h>

void
test_all_ratios_kernel()
{
  int N = 128;

  float *A, *A_d, *Ainv, *Ainv_d, *ratio, *ratio_d;

  cudaMalloc ((void**)&A_d,    N*N*sizeof(float));
  cudaMalloc ((void**)&Ainv_d, N*N*sizeof(float));  
  cudaMalloc ((void**)&ratio_d, 1*N*sizeof(float));
  A     = (float *)malloc (N*N*sizeof(float));
  Ainv  = (float *)malloc (N*N*sizeof(float));
  ratio = (float *)malloc (1*N*sizeof(float));

  float ratio2[N];
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++) {
      A[i*N+j] = 1.0f+drand48();
      Ainv[i*N+j] = 1.0f+drand48();
    }
  
  cudaMemcpy (A_d,     A,    N*N*sizeof(float), cudaMemcpyHostToDevice);  
  cudaMemcpy (Ainv_d,  Ainv, N*N*sizeof(float), cudaMemcpyHostToDevice);

  float **A_list, **A_list_d, **Ainv_list, **Ainv_list_d, **ratio_list, **ratio_list_d;
  int numMats = 2000;


  cudaMalloc ((void**)&A_list_d,     numMats*sizeof(float*));
  cudaMalloc ((void**)&Ainv_list_d,  numMats*sizeof(float*));
  cudaMalloc ((void**)&ratio_list_d, numMats*sizeof(float*));
  A_list     = (float **)malloc (numMats*sizeof(float*));
  Ainv_list  = (float **)malloc (numMats*sizeof(float*));
  ratio_list = (float **)malloc (numMats*sizeof(float*));

  for (int i=0; i<numMats; i++) {
    A_list[i] = A_d;
    Ainv_list[i] = Ainv_d;
    ratio_list[i] = ratio_d;
  }

  cudaMemcpy (A_list_d,    A_list,      numMats*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy (Ainv_list_d, Ainv_list,   numMats*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy (ratio_list_d, ratio_list, numMats*sizeof(float*), cudaMemcpyHostToDevice);

  clock_t start = clock();
  for (int i=0; i<1000; i++) 
    calc_all_ratios (Ainv_list_d, A_list_d, ratio_list_d, N, N, numMats);
  clock_t end = clock();
  double time = (double)(end-start)/(double)CLOCKS_PER_SEC;
  fprintf (stderr, "start = %d\n", start);
  fprintf (stderr, "end = %d\n", end);
  double rate = 1000.0/time;
  fprintf (stderr, "Rate = %1.2f generations per second.\n", rate);


  cudaMemcpy (ratio, ratio_d, N*sizeof(float), cudaMemcpyDeviceToHost);

  // for (int i=0; i<N; i++) {
  //   ratio2[i] = 0.0f;
  //   for (int j=0; j<N; j++)
  //     ratio2[i] += A[i*N+j]*Ainv[j*N+i];
  //   fprintf (stderr, "%3d  %10.6f  %10.6f\n", i, ratio2[i], ratio[i]);
  // }
  

}



void
test_all_grad_lapl_kernel()
{
  int N = 128;

  float *A, *A_d, *Ainv, *Ainv_d, *ratio, *ratio_d;

  cudaMalloc ((void**)&A_d,     4*N*N*sizeof(float));
  cudaMalloc ((void**)&Ainv_d,  N*N*sizeof(float));  
  cudaMalloc ((void**)&ratio_d, 4*N*sizeof(float));
  A     = (float *)malloc (4*N*N*sizeof(float));
  Ainv  = (float *)malloc (1*N*N*sizeof(float));
  ratio = (float *)malloc (4*N*sizeof(float));

  float ratio2[4*N];
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++) {
      Ainv[i*N+j] = 1.0f+drand48();
      for (int k=0; k<4; k++)
	A[4*(i*N+j)+k] = 1.0f+drand48();
    }
  
  cudaMemcpy (A_d,     A,    4*N*N*sizeof(float), cudaMemcpyHostToDevice);  
  cudaMemcpy (Ainv_d,  Ainv, 1*N*N*sizeof(float), cudaMemcpyHostToDevice);

  float **A_list, **A_list_d, **Ainv_list, **Ainv_list_d, **ratio_list, **ratio_list_d;
  int numMats = 2000;


  cudaMalloc ((void**)&A_list_d,     numMats*sizeof(float*));
  cudaMalloc ((void**)&Ainv_list_d,  numMats*sizeof(float*));
  cudaMalloc ((void**)&ratio_list_d, numMats*sizeof(float*));
  A_list     = (float **)malloc (numMats*sizeof(float*));
  Ainv_list  = (float **)malloc (numMats*sizeof(float*));
  ratio_list = (float **)malloc (numMats*sizeof(float*));

  for (int i=0; i<numMats; i++) {
    A_list[i] = A_d;
    Ainv_list[i] = Ainv_d;
    ratio_list[i] = ratio_d;
  }

  cudaMemcpy (A_list_d,    A_list,      numMats*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy (Ainv_list_d, Ainv_list,   numMats*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy (ratio_list_d, ratio_list, numMats*sizeof(float*), cudaMemcpyHostToDevice);

  struct timeval tstart, tend;
  gettimeofday(&tstart, NULL);
  for (int i=0; i<100; i++) 
    calc_grad_lapl (Ainv_list_d, A_list_d, ratio_list_d, N, N, numMats);
  cudaMemcpy (ratio, ratio_d, 4*N*sizeof(float), cudaMemcpyDeviceToHost);
  gettimeofday(&tend, NULL);
  double start = (double)tstart.tv_sec + 1.0e-6 * (double)tstart.tv_usec;  
  double end   = (double)tend.tv_sec   + 1.0e-6 * (double)tend.tv_usec;
  fprintf (stderr, "start = %f\n", start);
  fprintf (stderr, "end = %f\n", end);
  double rate = 100.0/(end-start);
  fprintf (stderr, "Rate = %1.2f generations per second.\n", rate);



  for (int i=0; i<N; i++) {
    for (int k=0; k<4; k++)
      ratio2[4*i+k] = 0.0f;
    for (int j=0; j<N; j++)
      for (int k=0; k<4; k++)
	ratio2[4*i+k] += A[(4*i+k)*N+j]*Ainv[j*N+i];
    for (int k=0; k<4; k++)
    fprintf (stderr, "%3d  %10.6f  %10.6f\n", 4*i+k, ratio2[4*i+k], ratio[4*i+k]);
  }
 

}



template<typename T>
__global__ void
woodbury_update_16 (T** Ainv_trans, T** delta,
		    T** Ainv_delta,
		    int N, int rowstride)
{
  T *myAinv, *mydelta, *myAinv_delta;
  int tid = threadIdx.x;

  myAinv       = Ainv_trans[blockIdx.y];
  myAinv_delta = Ainv_delta[blockIdx.y];
  mydelta      =      delta[blockIdx.y];
  int first_row = blockIdx.x*16;
  
  __shared__ T Ainv_s[16][17], delta_s[2][17], Ainv_delta_s[16][17];
  int nb = (N+15)/16;

  for (int row=0; row<16; row++)
    if (tid < 16)
      Ainv_delta_s[row][tid] = 0.0f;
  __syncthreads();

  int col = tid & 15;
  for (int block=0; block<nb; block++) {
    int nend = N - block*16;
    int c = block*16+tid;
    if (tid < 16)
      for (int row=0; row<16; row++) {
      Ainv_s[row][col]  = myAinv[(first_row+row)*rowstride+c];
      }
    __syncthreads();
    for (int irow=0; irow<8; irow++) {
      int odd = tid>15;
      int row = 2*irow + odd;
      delta_s[odd][tid] = mydelta[row*rowstride+c];
      if (row+first_row < N && col < nend)
	for (int k=0; k<16; k++)
	  Ainv_delta_s[row][col] += Ainv_s[col][k] *  delta_s[odd][k];
    }
    __syncthreads();
  }
  int mycol = blockIdx.x*16+tid;
  if (tid < 16 && mycol < N)
    for (int row=0; row<16; row++)
      myAinv_delta[row*rowstride+mycol] = Ainv_delta_s[row][tid];
}


// Require 64 threads per block
template<typename T>
__device__ inline void
block_inverse_16(T A[16][17])
{
  int tid = threadIdx.x;
  __shared__ T Arowk[16], Acolk[16];
  for (int k=0; k<16; k++) {
    T pivotInv = 1.0f/A[k][k];
    if (tid < 16) {
      T tmp = -pivotInv*A[tid][k];
      A[tid][k] = tmp;
      Acolk[tid] = tmp;
      Arowk[tid] = A[k][tid];
    }
    __syncthreads();
    int row= tid >> 4;
    int col = tid & 0x0f;
    for (;row<16;row+=4) 
      A[row][col] += Arowk[col]*Acolk[row];
    __syncthreads();
    A[k][tid] *= pivotInv;
  }
}
    


template<typename T>
__global__ void
woodbury_update_16a (T** Ainv_trans, T** delta,
		     T** Ainv_delta,
		     int N, int rowstride, int kblock)
{
  T *myAinv, *mydelta, *myAinv_delta;
  int tid = threadIdx.x;

  myAinv       = Ainv_trans[blockIdx.y];
  myAinv_delta = Ainv_delta[blockIdx.y];
  mydelta      =      delta[blockIdx.y];

  int first_row = blockIdx.x*16;
  
  __shared__ T Ainv_s[16][17], delta_s[4][17], Ainv_delta_s[16][17];
  int nb = (N+15)/16;

  int row = tid >> 4;
  int col = tid & 0x0000000f;
  for (int irow=0; irow<4; irow++,row+=4) 
    Ainv_delta_s[row][col] = 0.0f;
  __syncthreads();
  
  for (int block=0; block<nb; block++) {
    int c = block*16+ col;
    int row = tid >> 4;
    for (int irow=0; irow<4; irow++) {
      Ainv_s[row][col]  = myAinv[(first_row+row)*rowstride+c];
      row +=4;
    }
    __syncthreads();
    row = tid >> 4;
    for (int irow=0; irow<4; irow++) {
      delta_s[row][col] = mydelta[(row+4*irow)*rowstride+c];
      T mysum = Ainv_delta_s[row+4*irow][col];
      if (row+first_row < N && c < N)
	for (int k=0; k<16; k++)
	  mysum += Ainv_s[col][k] *  delta_s[row][k];
      Ainv_delta_s[row+4*irow][col] = mysum;
    }
    __syncthreads();
  }
  int mycol = blockIdx.x*16+col;
  row = tid >> 4;
  if (mycol < N) {
    for (int irow=0; irow<4; irow++,row+=4)
      myAinv_delta[row*rowstride+mycol] = Ainv_delta_s[row][col];
  }
  __syncthreads();
  if (blockIdx.x == kblock)
    block_inverse_16<T> (Ainv_delta_s);
}




template<typename T>
__global__ void
woodbury_update_32 (T** Ainv_trans, T** delta,
		    T** Ainv_delta,
		    int N, int rowstride)
{
  T *myAinv, *mydelta, *myAinv_delta;
  int tid = threadIdx.x;

  myAinv       = Ainv_trans[blockIdx.y];
  myAinv_delta = Ainv_delta[blockIdx.y];
  mydelta      =      delta[blockIdx.y];
  int first_row = blockIdx.x*32;
  
  __shared__ T Ainv_s[32][33], delta_s[32][33], Ainv_delta_s[32][33];
  int nb = (N+31)/32;

  for (int row=0; row<32; row++)
    if (tid < 32)
      Ainv_delta_s[row][tid] = 0.0f;
  __syncthreads();

  int col = tid;
  for (int block=0; block<nb; block++) {
    int nend = N - block*32;
    int c = block*32+tid;
    for (int row=0; row<32; row++) {
      Ainv_s[row][tid]  = myAinv[(first_row+row)*rowstride+c];
      delta_s[row][tid] = mydelta[row*rowstride+c];
    }
    __syncthreads();
    for (int row=0; row<32; row++) {
      if (row+first_row < N && col < nend)
	for (int k=0; k<32; k++)
	  Ainv_delta_s[row][col] += Ainv_s[row][k] *  delta_s[col][k];
    }
    __syncthreads();
  }
  int mycol = blockIdx.x*32+tid;
  if (mycol < N)
    for (int row=0; row<32; row++)
      myAinv_delta[row*rowstride+mycol] = Ainv_delta_s[row][tid];
}


#ifdef CUDA_TEST_MAIN



// Replaces A with its inverse by gauss-jordan elimination with full pivoting
// Adapted from Numerical Recipes in C
void GJInverse (double *A, int n)
{
  const int maxSize = 2000;

  if (n == 2) { // Special case for 2x2
    double a=A[0]; double b=A[1];
    double c=A[2]; double d=A[3];
    double detInv = 1.0/(a*d-b*c);
    A[0] = d*detInv;
    A[1] = -b*detInv;
    A[2] = -c*detInv;
    A[3] =  a*detInv;
    return;
  }

  int colIndex[maxSize], rowIndex[maxSize], ipiv[maxSize];
  double big, pivInv;
  int icol, irow;
  
  for (int j=0; j<n; j++)
    ipiv[j] = -1;

  for (int i=0; i<n; i++) {
    big = 0.0;
    for (int j=0; j<n; j++) 
      if (ipiv[j] != 0)
	for (int k=0; k<n; k++) {
	  if (ipiv[k] == -1) {
	    if (fabs(A[n*j+k]) >= big) {
	      big = fabs(A[n*j+k]);
	      irow = j; 
	      icol = k;
	    }
	  }
	  else if (ipiv[k] > 0) {
	    fprintf (stderr, "GJInverse: Singular matrix!\n");
	    exit(1);
	  }
	}
    ++(ipiv[icol]); 
    
    if (irow != icol) 
      for (int l=0; l<n; l++) {
	double tmp = A[n*irow+l];
	A[n*irow+l] = A[n*icol+l];
	A[n*icol+l] = tmp;
	// swap (A[n*irow+l], A[n*icol+l]);
      }
			     
    
    rowIndex[i] = irow;
    colIndex[i] = icol;
    if (A[n*icol+icol] == 0.0) { 
      fprintf (stderr, "GJInverse: Singular matrix!\n");
      exit(1);
    }
    pivInv = 1.0/A[n*icol+icol];
    A[n*icol+icol] = 1.0;
    for (int l=0; l<n; l++)
      A[n*icol+l] *= pivInv;
    for (int ll=0; ll<n; ll++)
      if (ll != icol) {
	double dum = A[n*ll+icol];
	A[n*ll+icol] = 0.0;
	for (int l=0; l<n; l++)
	  A[n*ll+l] -= A[n*icol+l]*dum;
      }
  }
  // Now unscramble the permutations
  for (int l=n-1; l>=0; l--) {
    if (rowIndex[l] != colIndex[l])
      for (int k=0; k<n ; k++) {
	double tmp = A[n*k+rowIndex[l]];
	A[n*k+rowIndex[l]] = A[n*k+colIndex[l]];
	A[n*k+colIndex[l]] = tmp;
	// swap (A(k,rowIndex[l]),A(k, colIndex[l]));
      }
  }
}

#include <omp.h>



#define MAT_SIZE 256
#define NUM_MATS 512

void 
test_update()
{
  int const N = MAT_SIZE;
  double *A, *Ainv;
  int numMats = NUM_MATS;
  float *A_h, *Ainv_h, *u_h;
  float *Ainv_d, *Ainv_u_d, *Ainv_colk_d, *u_d;

  A       = (double*)malloc (N*N*sizeof(double));
  Ainv    = (double*)malloc (N*N*sizeof(double));
  Ainv_h  = (float*) malloc (N*N*sizeof(float));
  A_h     = (float*) malloc (N*N*sizeof(float));
  u_h     = (float*) malloc (N*sizeof(float));
  cudaMalloc((void**)&Ainv_d,  N*N*sizeof(float));
  cudaMalloc((void**)&Ainv_d, N*N*sizeof(float));
  cudaMalloc((void**)&u_d, N*sizeof(float));
  cudaMalloc((void**)&Ainv_u_d, N*sizeof(float));
  cudaMalloc((void**)&Ainv_colk_d, N*sizeof(float));
  
  float **AinvList, **Ainv_uList, **AList, 
    **Ainv_colkList, **uList;

  AList        = (float**)malloc(NUM_MATS*sizeof(float*));
  AinvList     = (float**)malloc(NUM_MATS*sizeof(float*));
  Ainv_uList    = (float**)malloc(NUM_MATS*sizeof(float*));
  Ainv_colkList = (float**)malloc(NUM_MATS*sizeof(float*));
  uList         = (float**)malloc(NUM_MATS*sizeof(float*));

  float **AList_d, **AinvList_d, **Ainv_uList_d, **Ainv_colkList_d, **uList_d;
  cudaMalloc((void**)&AList_d,         numMats*sizeof(float*));
  cudaMalloc((void**)&AinvList_d,      numMats*sizeof(float*));
  cudaMalloc((void**)&Ainv_uList_d,    numMats*sizeof(float*));
  cudaMalloc((void**)&Ainv_colkList_d, numMats*sizeof(float*));
  cudaMalloc((void**)&uList_d,         numMats*sizeof(float*));

  fprintf (stderr, "N = %d\n", N);

  
  for (int mat=0; mat<numMats; mat++) {
    cudaMalloc((void**)&(AList[mat])        , N*N*sizeof(float)+1000);
    cudaMalloc((void**)&(AinvList[mat])     , N*N*sizeof(float)+1000);
    cudaMalloc((void**)&(Ainv_uList[mat])   ,   N*sizeof(float)+1000);
    cudaMalloc((void**)&(Ainv_colkList[mat]),   N*sizeof(float)+1000);
    cudaMalloc((void**)&(uList[mat])        ,   N*sizeof(float)+1000);
  }

  fprintf (stderr, "N = %d\n", N);


  cudaMemcpy (AList_d, AList, numMats*sizeof(float*), 
	      cudaMemcpyHostToDevice);
  cudaMemcpy (AinvList_d, AinvList, numMats*sizeof(float*), 
	      cudaMemcpyHostToDevice);
  cudaMemcpy (Ainv_uList_d, Ainv_uList, numMats*sizeof(float*), 
	      cudaMemcpyHostToDevice);
  cudaMemcpy (Ainv_colkList_d, Ainv_colkList, numMats*sizeof(float*), 
	      cudaMemcpyHostToDevice);
  cudaMemcpy (uList_d, uList, numMats*sizeof(float*), 
	      cudaMemcpyHostToDevice);
  
  srand48((long int) 12341313);

  fprintf (stderr, "N = %d\n", N);


  int row = 0;

  for (int mat=0; mat<numMats; mat++) {
    if (mat == 0 ) {
      for (int i=0; i<N; i++) {
	u_h[i] = drand48();
	for (int j=0; j<N; j++) 
	  A[i*N+j] = Ainv[i*N+j] = A_h[i*N+j] = drand48();
      }
      // for (int i=0; i<N; i++)
      // 	u_h[i] = A_h[row*N+i];

      GJInverse(Ainv, N);
      for (int i=0; i<N; i++)
	for (int j=0; j<N; j++) 
	  Ainv_h[i*N+j] = (float)Ainv[i*N+j];
    }
    // for (int i=0; i<N; i++)
    //   u_h[i] = A_h[row*N+i];

    cudaMemcpy (AList[mat], A_h, N*N*sizeof(float), 
		cudaMemcpyHostToDevice);
    cudaMemcpy (AinvList[mat], Ainv_h, N*N*sizeof(float), 
		cudaMemcpyHostToDevice);
    cudaMemcpy (uList[mat], u_h, N*sizeof(float), cudaMemcpyHostToDevice);
  }

  dim3 dimBlock2(64);
  dim3 dimGrid2((N+63)/64, NUM_MATS);


  double start = omp_get_wtime();

  for (int i=0; i<1000; i++) {
    update_inverse_cuda1<float,64><<<dimGrid2,dimBlock2>>>
      (AList_d, AinvList_d, uList_d, Ainv_uList_d, Ainv_colkList_d, N, N, row);
    update_inverse_cuda2<float,64><<<dimGrid2,dimBlock2>>>
      (AList_d, AinvList_d, uList_d, Ainv_uList_d, Ainv_colkList_d, N, N, row);
  }
  cudaThreadSynchronize();
  double end = omp_get_wtime();
  fprintf (stderr, "Rate = %12.8f updates per second.\n",
	   (double)(1000*NUM_MATS)/(end - start));
  cudaMemcpy (Ainv_h, AinvList[0], N*N*sizeof(float),cudaMemcpyDeviceToHost);

  

  /*  for (int j=0; j<16; j++)
    for (int i=0; i<N; i++)
      A[(row+j)*N+i] += delta_h[j*N+i];
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++) {
      double ident = 0.0;
      for (int k=0; k<N; k++)
  	ident += Ainv_h[i*N+k]*A[k*N+j];
      if ((i==j && fabs(ident - 1.0) > 1.0e-4) ||
      	  (i!=j && fabs(ident) > 1.0e-4))
      	fprintf (stderr, "Error in matrix inverse, (%d, %d) = %1.8f\n", i, j, ident);
	}*/
  fprintf (stderr, "Finished.\n");
}

void 
test_update_transpose()
{
  const int N = MAT_SIZE;
  double *A, *Ainv;
  int numMats = NUM_MATS;
  float *A_h, *Ainv_h, *u_h;
  float *Ainv_d, *Ainv_u_d, *Ainv_colk_d, *u_d;

  A       = (double*)malloc (N*N*sizeof(double));
  Ainv    = (double*)malloc (N*N*sizeof(double));
  Ainv_h  = (float*) malloc (N*N*sizeof(float));
  A_h     = (float*) malloc (N*N*sizeof(float));
  u_h     = (float*) malloc (N*sizeof(float));
  cudaMalloc((void**)&Ainv_d,  N*N*sizeof(float));
  cudaMalloc((void**)&Ainv_d, N*N*sizeof(float));
  cudaMalloc((void**)&u_d, N*sizeof(float));
  cudaMalloc((void**)&Ainv_u_d, N*sizeof(float));
  cudaMalloc((void**)&Ainv_colk_d, N*sizeof(float));
  
  float **AinvList, **Ainv_uList, **AList, 
    **Ainv_colkList, **uList;

  AList        = (float**)malloc(NUM_MATS*sizeof(float*));
  AinvList     = (float**)malloc(NUM_MATS*sizeof(float*));
  Ainv_uList    = (float**)malloc(NUM_MATS*sizeof(float*));
  Ainv_colkList = (float**)malloc(NUM_MATS*sizeof(float*));
  uList         = (float**)malloc(NUM_MATS*sizeof(float*));

  float **AList_d, **AinvList_d, **Ainv_uList_d, **Ainv_colkList_d, **uList_d;
  cudaMalloc((void**)&AList_d,         numMats*sizeof(float*));
  cudaMalloc((void**)&AinvList_d,      numMats*sizeof(float*));
  cudaMalloc((void**)&Ainv_uList_d,    numMats*sizeof(float*));
  cudaMalloc((void**)&Ainv_colkList_d, numMats*sizeof(float*));
  cudaMalloc((void**)&uList_d,         numMats*sizeof(float*));

  fprintf (stderr, "N = %d\n", N);

  
  for (int mat=0; mat<numMats; mat++) {
    cudaMalloc((void**)&(AList[mat])        , N*N*sizeof(float));
    cudaMalloc((void**)&(AinvList[mat])     , N*N*sizeof(float));
    cudaMalloc((void**)&(Ainv_uList[mat])   ,   N*sizeof(float));
    cudaMalloc((void**)&(Ainv_colkList[mat]),   N*sizeof(float));
    cudaMalloc((void**)&(uList[mat])        ,   N*sizeof(float));
  }

  fprintf (stderr, "N = %d\n", N);


  cudaMemcpy (AList_d, AList, numMats*sizeof(float*), 
	      cudaMemcpyHostToDevice);
  cudaMemcpy (AinvList_d, AinvList, numMats*sizeof(float*), 
	      cudaMemcpyHostToDevice);
  cudaMemcpy (Ainv_uList_d, Ainv_uList, numMats*sizeof(float*), 
	      cudaMemcpyHostToDevice);
  cudaMemcpy (Ainv_colkList_d, Ainv_colkList, numMats*sizeof(float*), 
	      cudaMemcpyHostToDevice);
  cudaMemcpy (uList_d, uList, numMats*sizeof(float*), 
	      cudaMemcpyHostToDevice);
  
  srand48((long int) 12341313);

  fprintf (stderr, "N = %d\n", N);


  int row = 1;

  for (int mat=0; mat<numMats; mat++) {
    if (mat == 0 ) {
      for (int i=0; i<N; i++) {
	for (int j=0; j<N; j++) 
	  A[i*N+j] = Ainv[i*N+j] = A_h[i*N+j] = drand48();
	//	u_h[i] = drand48();
      }
      for (int j=0; j<N; j++)
	u_h[j] = drand48();//A[N*row+j];
   
      GJInverse(Ainv, N);
      for (int i=0; i<N; i++)
	for (int j=0; j<N; j++) 
	  Ainv_h[j*N+i] = (float)Ainv[i*N+j];
    }

    // for (int i=0; i<N; i++)
    //   u_h[i] = A_h[row*N+i];

    cudaMemcpy (AList[mat], A_h, N*N*sizeof(float), 
		cudaMemcpyHostToDevice);
    cudaMemcpy (AinvList[mat], Ainv_h, N*N*sizeof(float), 
		cudaMemcpyHostToDevice);
    cudaMemcpy (uList[mat], u_h, N*sizeof(float), cudaMemcpyHostToDevice);
  }

  dim3 dimBlock(DET_BLOCK_SIZE);
  dim3 dimGrid(NUM_MATS);

  clock_t start = clock();


  for (int i=0; i<1000; i++) {
    update_inverse_transpose_cuda<float,DET_BLOCK_SIZE,N><<<dimGrid,dimBlock>>>
      (AList_d, AinvList_d, uList_d, N, N, row);
    // update_inverse_transpose_cuda_2pass<float,DET_BLOCK_SIZE,N><<<dimGrid,dimBlock>>>
    //   (AList_d, AinvList_d, uList_d, N, N, row);
  }
  cudaThreadSynchronize();
  clock_t end = clock();
  fprintf (stderr, "Rate = %12.8f updates per second.\n",
	   (double)(1000*NUM_MATS)/((double)(end - start)/(double)CLOCKS_PER_SEC));

  cudaMemcpy (Ainv_h, AinvList[1], N*N*sizeof(float),cudaMemcpyDeviceToHost);

  for (int i=0; i<N; i++)
    A[row*N+i] = u_h[i];
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++) {
      double ident = 0.0;
      for (int k=0; k<N; k++)
  	ident += Ainv_h[k*N+i]*A[k*N+j];
      if ((i==j && fabs(ident - 1.0) > 1.0e-4) ||
      	  (i!=j && fabs(ident) > 1.0e-4))
      	fprintf (stderr, "Error in matrix inverse, (%d, %d) = %1.8f\n", i, j, ident);
    }
  fprintf (stderr, "Finished.\n");
}




void 
test_woodbury()
{
  int const N = MAT_SIZE;
  int M = 16;
  double *A, *Ainv;
  int numMats = NUM_MATS;
  float *A_h, *Ainv_h, *delta_h, *Ainv_delta_h;
  float *Ainv_d, *Ainv_delta_d, *Ainv_colk_d, *delta_d;

  A       = (double*)malloc (N*N*sizeof(double));
  Ainv    = (double*)malloc (N*N*sizeof(double));
  Ainv_h  = (float*) malloc (N*N*sizeof(float));
  A_h     = (float*) malloc (N*N*sizeof(float));
  delta_h     = (float*) malloc (N*M*sizeof(float));
  Ainv_delta_h     = (float*) malloc (N*M*sizeof(float));
  cudaMalloc((void**)&Ainv_d,       N*N*sizeof(float));
  cudaMalloc((void**)&delta_d,      N*M*sizeof(float));
  cudaMalloc((void**)&Ainv_delta_d, N*M*sizeof(float));
  cudaMalloc((void**)&Ainv_colk_d,  N  *sizeof(float));
  
  float **AinvList, **Ainv_deltaList, **AList, 
    **Ainv_colkList, **deltaList;

  AList          = (float**)malloc(NUM_MATS*sizeof(float*));
  AinvList       = (float**)malloc(NUM_MATS*sizeof(float*));
  Ainv_deltaList = (float**)malloc(NUM_MATS*sizeof(float*));
  Ainv_colkList  = (float**)malloc(NUM_MATS*sizeof(float*));
  deltaList      = (float**)malloc(NUM_MATS*sizeof(float*));

  float **AList_d, **AinvList_d, **Ainv_deltaList_d, 
    **Ainv_colkList_d, **deltaList_d;
  cudaMalloc((void**)&AinvList_d,          numMats*sizeof(float*));
  cudaMalloc((void**)&Ainv_deltaList_d,    numMats*sizeof(float*));
  cudaMalloc((void**)&deltaList_d,         numMats*sizeof(float*));

  fprintf (stderr, "N = %d\n", N);
  
  for (int mat=0; mat<numMats; mat++) {
    //    cudaMalloc((void**)&(AList[mat])        , N*N*sizeof(float)+1000);
    cudaMalloc((void**)&(AinvList[mat])      , N*N*sizeof(float)+1000);
    cudaMalloc((void**)&(Ainv_deltaList[mat]), N*M*sizeof(float)+1000);
    cudaMalloc((void**)&(deltaList[mat])     , N*M*sizeof(float)+1000);
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
      fprintf (stderr, "CUDA error in test_woodbury malloc:\n  %s\n",
	       cudaGetErrorString(err));
      abort();
    }

  }

  fprintf (stderr, "N = %d\n", N);


  cudaMemcpy (AinvList_d, AinvList, numMats*sizeof(float*), 
	      cudaMemcpyHostToDevice);
  cudaMemcpy (Ainv_deltaList_d, Ainv_deltaList, numMats*sizeof(float*), 
	      cudaMemcpyHostToDevice);
  cudaMemcpy (deltaList_d, deltaList, numMats*sizeof(float*), 
	      cudaMemcpyHostToDevice);
  
  srand48((long int) 12341313);

  fprintf (stderr, "N = %d    M = %d\n", N, M);
  int row = 0;

  for (int mat=0; mat<numMats; mat++) {
    if (mat == 0 ) {
      for (int i=0; i<N; i++) {
	for (int j=0; j<16; j++)
	  delta_h[j*N+i] = drand48();
	for (int j=0; j<N; j++) 
	  A[i*N+j] = Ainv[i*N+j] = A_h[i*N+j] = drand48();
      }
      // for (int i=0; i<N; i++)
      // 	delta_h[i] = A_h[row*N+i];

      GJInverse(Ainv, N);
      for (int i=0; i<N; i++)
	for (int j=0; j<N; j++) 
	  Ainv_h[i*N+j] = (float)Ainv[i*N+j];
    }
    // for (int i=0; i<N; i++)
    //   delta_h[i] = A_h[row*N+i];

    // cudaMemcpy (AList[mat], A_h, N*N*sizeof(float), 
    // 		cudaMemcpyHostToDevice);
    cudaMemcpy (AinvList[mat], Ainv_h, N*N*sizeof(float), 
		cudaMemcpyHostToDevice);
    cudaMemcpy (deltaList[mat], delta_h, N*M*sizeof(float), cudaMemcpyHostToDevice);
  }

  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in test_woodbury memcopy's:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }


  dim3 dimBlock2(64);
  dim3 dimGrid2((N/16), numMats);
  //dim3 dimGrid2((N/32), numMats);

  double start = omp_get_wtime();
  for (int i=0; i<100; i++) {
    woodbury_update_16a<float><<<dimGrid2,dimBlock2>>>
      (AinvList_d, deltaList_d, Ainv_deltaList_d, N, N, 0);
    // woodbury_update_32<float><<<dimGrid2,dimBlock2>>>
    //   (AinvList_d, deltaList_d, Ainv_deltaList_d, N, N);

    cudaThreadSynchronize();
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
      fprintf (stderr, "CUDA error in woodbury_update_16:\n  %s\n",
	       cudaGetErrorString(err));
      abort();
    }
  }
  double end = omp_get_wtime();
  fprintf (stderr, "Rate = %12.8f updates per second.\n",
	   (double)(100*NUM_MATS)/(end - start));

  fprintf (stderr, "About to copy %ld back\n", N*M*sizeof(float));

  cudaMemcpy (Ainv_delta_h, Ainv_deltaList[0], N*M*sizeof(float),
   	      cudaMemcpyDeviceToHost);
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in test_woodbury memcopy back:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }
  fprintf(stderr, "Copied result back.\n");

  float Ainv_delta[N*M];
  for (int i=0; i<N*M; i++)
    Ainv_delta[i] = 0.0;
  for (int i=0; i<16; i++)
    for (int j=0; j<N; j++)
      for (int k=0; k<N; k++)
	Ainv_delta[i*N+j] += Ainv_h[j*N+k]*delta_h[i*N+k];

  fprintf (stderr, "Ainv_delta_cpu = %1.8e\n", Ainv_delta[51]);
  fprintf (stderr, "Ainv_delta_gpu = %1.8e\n", Ainv_delta_h[51]);

  // cudaMemcpy (Ainv_h, AinvList[0], N*N*sizeof(float),cudaMemcpyDeviceToHost);

  // for (int i=0; i<N; i++)
  //   A[row*N+i] = delta_h[i];
  // for (int i=0; i<N; i++)
  //   for (int j=0; j<N; j++) {
  //     double ident = 0.0;
  //     for (int k=0; k<N; k++)
  // 	ident += Ainv_h[i*N+k]*A[k*N+j];
  //     if ((i==j && fabs(ident - 1.0) > 1.0e-4) ||
  //     	  (i!=j && fabs(ident) > 1.0e-4))
  //     	fprintf (stderr, "Error in matrix inverse, (%d, %d) = %1.8f\n", i, j, ident);
  //   }
  fprintf (stderr, "Finished.\n");
}




// Compile with:
// nvcc -o test_all_ratios -DCUDA_TEST_MAIN ../src/QMCWaveFunctions/Fermion/determinant_update.cu
main()
{
  //test_all_ratios_kernel();
  // test_all_grad_lapl_kernel();
  test_update();
  // test_update_transpose();
  test_woodbury();
}



#endif
