//////////////////////////////////////////////////////////////////////////////////////
//// This file is distributed under the University of Illinois/NCSA Open Source License.
//// See LICENSE file in top directory for details.
////
//// Copyright (c) 2019 QMCPACK developers.
////
//// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
////
//// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
////////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_CUDA_ERROR_H
#define QMCPLUSPLUS_CUDA_ERROR_H

#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cuda_runtime_api.h>

#define cudaErrorCheck(ans, cause)                \
  {                                               \
    cudaAssert((ans), __FILE__, __LINE__, cause);	  \
  }

// The cause is largely redundant with the __FILE__ and __LINE__ information
// and it makes CUDA heavy code tedious to write and read
#define cudaCheck(ans)                            \
  {                                               \
    cudaAssert((ans), __FILE__, __LINE__); \
  }                                               \


/// prints CUDA error messages. Always use cudaErrorCheck macro.
inline void cudaAssert(cudaError_t code, const std::string& cause, int line, const char* filename = "")
{
  if (code != cudaSuccess)
  {
    std::ostringstream err;
    err << "cudaAssert: " << cudaGetErrorName(code) << " " << cudaGetErrorString(code) << ", file " << filename
        << ", line " << line << std::endl
        << cause << std::endl;
    std::cerr << err.str();
    throw std::runtime_error(cause);
  }
}

#endif
