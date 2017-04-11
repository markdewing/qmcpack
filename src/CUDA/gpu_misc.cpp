//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ying Wai Li, yingwaili@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "config.h"
#include "gpu_misc.h"

namespace gpu
{
cudaStream_t kernelStream;
cudaStream_t memoryStream;

cudaEvent_t syncEvent;

cudaEvent_t gradientSyncDiracEvent;
cudaEvent_t gradientSyncOneBodyEvent;
cudaEvent_t gradientSyncTwoBodyEvent;
cudaEvent_t ratioSyncDiracEvent;
cudaEvent_t ratioSyncOneBodyEvent;
cudaEvent_t ratioSyncTwoBodyEvent;
cublasHandle_t cublasHandle;

size_t MaxGPUSpineSizeMB;
int rank;

void
initCUDAStreams()
{
  cudaStreamCreate(&kernelStream);
  cudaStreamCreate(&memoryStream);
}

void
initCUDAEvents()
{
  unsigned int flag = cudaEventDisableTiming;
#ifdef ENABLE_TIMERS
  flag = cudaEventDefault;
#endif
  cudaEventCreateWithFlags(&syncEvent, flag);
  cudaEventCreateWithFlags(&gradientSyncDiracEvent, flag);
  cudaEventCreateWithFlags(&gradientSyncOneBodyEvent, flag);
  cudaEventCreateWithFlags(&gradientSyncTwoBodyEvent, flag);
  cudaEventCreateWithFlags(&ratioSyncDiracEvent, flag);
  cudaEventCreateWithFlags(&ratioSyncOneBodyEvent, flag);
  cudaEventCreateWithFlags(&ratioSyncTwoBodyEvent, flag);
}

void
initCublas()
{
  cublasCreate(&cublasHandle);
}

void
finalizeCUDAStreams()
{
  cudaStreamDestroy(kernelStream);
  cudaStreamDestroy(memoryStream);
}

void
finalizeCUDAEvents()
{
  cudaEventDestroy(syncEvent);
  cudaEventDestroy(gradientSyncDiracEvent);
  cudaEventDestroy(gradientSyncOneBodyEvent);
  cudaEventDestroy(gradientSyncTwoBodyEvent);
  cudaEventDestroy(ratioSyncDiracEvent);
  cudaEventDestroy(ratioSyncOneBodyEvent);
  cudaEventDestroy(ratioSyncTwoBodyEvent);
}

void
finalizeCublas()
{
  cublasDestroy(cublasHandle);
}

void
synchronize()
{
  cudaDeviceSynchronize();
}

void
streamsSynchronize()
{
  cudaEventRecord(syncEvent, 0);
}

#ifdef ENABLE_TIMERS
gpu_timer::gpu_timer()
{
  cudaEventCreate(&startEvent);
  cudaEventCreate(&stopEvent);
}

void gpu_timer::start()
{
  cudaEventRecord(startEvent);
}

void gpu_timer::stop()
{
  cudaEventRecord(stopEvent);
}

double gpu_timer::elapsed()
{
  float ms;
  cudaError_t err = cudaEventElapsedTime(&ms, startEvent, stopEvent);
  if (err == cudaErrorNotReady)
  {
    //printf("GPU section not done\n");
    cudaError_t err2 = cudaEventSynchronize(stopEvent);
    if (err2 != cudaSuccess)
    {
      printf("GPU stop event sync error = %d\n",err);
    }
    err = cudaEventElapsedTime(&ms, startEvent, stopEvent);
    if (err != cudaSuccess)
    {
      printf("GPU event timer try again error = %d\n",err);
    }
  }
  if (err != cudaSuccess)
  {
    printf("GPU event timer error = %d\n",err);
  }
  return ms/1000.0;
}
#endif

}
