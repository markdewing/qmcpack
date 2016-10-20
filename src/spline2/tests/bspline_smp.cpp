//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp. 
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file bspline_smp.cpp
 * @brief Benchmark einspline. Shared engine among threads.
 */
#include <Configuration.h>
#include <Utilities/OhmmsInfo.h>
#include <Message/Communicate.h>
#include <mpi/collectives.h>
#include <getopt.h>
#include <spline2/tests/einspline_benchmark.h>
using namespace qmcplusplus;

int main(int argc, char** argv)
{
  OHMMS::Controller->initialize(argc,argv);
  Communicate* mycomm=OHMMS::Controller;
  OhmmsInfo Welcome("bspline_smp",mycomm->rank());
  qmcplusplus::Random.init(0,1,11);
  int nx=48,ny=48,nz=48;
  int num_splines=128;
  int nsamples=512;
  int niters=10;
  int opt;
  bool initialize=true;
  while((opt = getopt(argc, argv, "hfg:x:y:z:i:s:p:")) != -1)
  {
    switch(opt)
    {
    case 'h':
      printf("[-g grid| -x grid_x -y grid_y -z grid_z] -s states -p particles -i iterations\n");
      return 1;
    case 'g':
      nx=ny=nz=atoi(optarg);
      break;
    case 'x':
      nx=atoi(optarg);
      break;
    case 'y':
      ny=atoi(optarg);
      break;
    case 'z':
      nz=atoi(optarg);
      break;
    case 's':
      num_splines=atoi(optarg);
      break;
    case 'p':
      nsamples=atoi(optarg);
      break;
    case 'i':
      niters=atoi(optarg);
      break;
    case 'f': //skip spline
      initialize=false;
      break;
    }
  }
  int nthreads=omp_get_max_threads();

  typedef TinyVector<double,3> timer_type;
  timer_type d_timer_t,s_timer_t,z_timer_t,c_timer_t;
  Timer big_clock;
  //einspline3d_benchmark<double> d_bench;
  //d_bench.set(nx,ny,nz,num_splines);
  einspline3d_benchmark<float> s_bench;
  s_bench.set(nx,ny,nz,num_splines);
  double t_init=big_clock.elapsed();
  app_log() << "#einspline benchmark grid = " << nx << " " << ny << " " << nz
            << " num_splines = " << num_splines << " num_samples = " << nsamples
            << " iterations = " << niters
            << " number of operations in millions " << endl;
  app_log() << "#MPI = " << mycomm->size() << "  OMP_NUM_THREADS = " << omp_get_max_threads() << endl;
  app_log() << "#   mpi   openmp    datatype     "
            << "  value_op         vgl_op              vgh_op         value_time       vgl_time         vgh_time" << endl;
  big_clock.restart();
  app_log().flush();

  #pragma omp parallel
  {
    random_position_generator<float> s_pos(nsamples,omp_get_thread_num());
    einspline3d_benchmark<float> s_bench_loc(s_bench);;
    timer_type s_timer;
    for(int i=0; i<niters; ++i)
    {
      s_pos.randomize();
      s_timer+=s_bench_loc.test_all(s_pos.Vpos, s_pos.VGLpos, s_pos.VGHpos);
    }
#pragma omp critical
    {
      s_timer_t += s_timer;
    }
  }

  double t_comp=big_clock.elapsed();
  mpi::reduce(*mycomm,s_timer_t);

  double nops=num_splines*nsamples*1.e-6*mycomm->size()*omp_get_max_threads();;
  double tfac=1.0/static_cast<double>(mycomm->size()*omp_get_max_threads()*niters);

  s_timer_t*=tfac;

  if(mycomm->rank()==0)
  {
    cout.setf(std::ios::scientific, std::ios::floatfield);
    cout.precision(6);
    cout << argv[0]  <<setw(4) << mycomm->size() << setw(4) << omp_get_max_threads() <<  " single Ops  "<< std::fixed << nops/double(s_timer_t[0]) << "   " << nops/double(s_timer_t[1]) << "   " << nops/double(s_timer_t[2])<< endl;
    cout << argv[0] <<setw(4) << mycomm->size() << setw(4) << omp_get_max_threads() <<  " single sec  "<< std::fixed  << s_timer_t << endl;
    cout << " Initialization = " << std::fixed  << t_init << endl;
    cout << " Total time = " << std::fixed << t_comp/omp_get_max_threads() << endl;
  }

  OHMMS::Controller->finalize();
  return 0;
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: OrbitalBase.h 1770 2007-02-17 23:45:38Z jnkim $
 ***************************************************************************/
