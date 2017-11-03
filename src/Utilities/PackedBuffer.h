//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef PACKEDBUFFER_H_
#define PACKEDBUFFER_H_

#include <cstring>
#include <cassert>

// Pack data types into a buffer
// Replacement for OOMPI_Packed (based on MPI_Pack)

#if 1
#define CHECK_BOUNDS(x) assert(x)
#else
#define CHECK_BOUNDS(x)
#endif


class PackedBuffer
{
public:
    char *buf;
    int idx;
    int bufSize;
    PackedBuffer(int byteSize) : idx(0) {
       buf = new char[byteSize]; bufSize=byteSize;
    }

    void reset() {idx = 0; }

    ~PackedBuffer() {
      delete[] buf;
    }

    template<typename T> void Pack(T* start, int count)
    {
      int byteSize = count*sizeof(T);
      CHECK_BOUNDS(idx + byteSize <= bufSize);
      std::memcpy(&buf[idx], start, byteSize);
      idx += byteSize;
    }

    template<typename T> void Unpack(T* start, int count)
    {
      int byteSize = count*sizeof(T);
      CHECK_BOUNDS(idx + byteSize <= bufSize);
      std::memcpy(start, &buf[idx], byteSize);
      idx += byteSize;
    }

    template<typename T> PackedBuffer& operator<<(const T &v)
    {
      CHECK_BOUNDS(idx + sizeof(T) <= bufSize);
      *(T*)(&buf[idx]) = v;
      idx += sizeof(T);
      return *this;
    }

    template<typename T> PackedBuffer& operator>>(T &v)
    {
      CHECK_BOUNDS(idx + sizeof(T) <= bufSize);
      v = *(T*)(&buf[idx]);
      idx += sizeof(T);
      return *this;
    }
};

#endif
