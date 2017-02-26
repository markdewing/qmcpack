//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory 
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory

//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file HybridAdoptorBase.h
 *
 * Hybrid adoptor base class
 */
#ifndef QMCPLUSPLUS_HYBRID_ADOPTOR_BASE_H
#define QMCPLUSPLUS_HYBRID_ADOPTOR_BASE_H

#include <Particle/DistanceTableData.h>
#include <QMCWaveFunctions/lcao/SoaSphericalTensor.h>

namespace qmcplusplus
{

template<typename ST>
struct AtomicOrbitalSoA
{
  static const int D=3;
  using AtomicSplineType=typename bspline_traits<ST,1>::SplineType;
  using AtomicBCType=typename bspline_traits<ST,1>::BCType;
  using AtomicSingleSplineType=UBspline_1d_d;
  using PointType=TinyVector<ST,D>;
  using value_type=ST;

  using vContainer_type=aligned_vector<ST>;

  ST cutoff,spline_radius;
  int spline_npoints, BaseN;
  PointType pos;
  const int lmax, lm_tot, NumBands, Npad;
  SoaSphericalTensor<ST> Ylm;
  AtomicSplineType* MultiSpline;

  vContainer_type localV, localG, localL;

  AtomicOrbitalSoA(int Lmax, int Nb):
  Ylm(Lmax), NumBands(Nb), MultiSpline(nullptr), lmax(Lmax),
  lm_tot((Lmax+1)*(Lmax+1)), Npad(getAlignedSize<ST>(Nb))
  {
    localV.resize(Npad*lm_tot);
    localG.resize(Npad*lm_tot);
    localL.resize(Npad*lm_tot);
  }

  //~AtomicOrbitalSoA();

  template<typename PT, typename VT>
  inline void set_info(const PT& R, const VT& cutoff_in, const VT& spline_radius_in, const int& spline_npoints_in)
  {
    pos[0]=R[0];
    pos[1]=R[1];
    pos[2]=R[2];
    cutoff=cutoff_in;
    spline_radius=spline_radius_in;
    spline_npoints=spline_npoints_in;
    BaseN=spline_npoints+2;
  }

  inline void create_spline()
  {
    BCtype_d bc;
    bc.lCode = FLAT;
    bc.rCode = NATURAL;
    Ugrid grid;
    grid.start = 0.0;
    grid.end   = spline_radius;
    grid.num   = spline_npoints;
    MultiSpline = einspline::create(MultiSpline, grid, bc, lm_tot*Npad);
  }

  inline void set_spline(AtomicSingleSplineType* spline, int lm, int ispline)
  {
    einspline::set(MultiSpline, lm*Npad+ispline, spline, 0, BaseN);
  }

  bool read_splines(hdf_archive& h5f)
  {
    //load spline coefficients
    //cutoff, spline radius, spline points, check coeff size
  }

  bool write_splines(hdf_archive& h5f)
  {
    //dump center info, including cutoff, spline radius, spline points, lm
    //name, position, consistency
    //dump spline coefficients    //dump spline coefficients
  }

  template<typename VV>
  inline void evaluate_v(const ST& r, const PointType& dr, VV& myV)
  {
    //evaluate only V
    Ylm.evaluateV(dr[0], dr[1], dr[2]);
    const ST* restrict Ylm_v=Ylm[0];
    einspline::evaluate(MultiSpline,r,localV);
    CONSTEXPR ST czero(0);
    std::fill(myV.begin(),myV.end(),czero);
    for(size_t lm=0; lm<lm_tot; lm++)
    {
      size_t offset=lm*Npad;
      for(size_t ib=0; ib<myV.size(); ib++)
        myV[ib]+=Ylm_v[lm]*localV[offset+ib];
    }
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const ST& r, const PointType& dr, VV& myV, GV& myG, VV& myH)
  {
    //missing
  }

  template<typename VV, typename GV, typename HT>
  void evaluate_vgh(const ST& r, const PointType& dr, VV& myV, GV& myG, HT& myH)
  {
    //Needed to do tensor product here
    //einspline::evaluate_vgh(MultiSpline,r,myV,myG,myH);
    //cheating
    evaluate_v(r, dr, myV);
  }
};

/** adoptor class to match 
 *
 */
template<typename ST>
struct HybridAdoptorBase
{
  using PointType=typename AtomicOrbitalSoA<ST>::PointType;

  std::vector<AtomicOrbitalSoA<ST> > AtomicCenters;
  // I-e distance table
  DistanceTableData* ei_dist;
  //mapping supercell to primitive cell
  std::vector<int> Super2Prim;

  HybridAdoptorBase() { }

  bool read_splines(hdf_archive& h5f)
  {
    // for each center
    // loop load center info, including name, position, check consistency
    // Initial class with index
    // call read_splines
    // push to the vector
  }

  bool write_splines(hdf_archive& h5f)
  {
    //loop for each center
    // write_splines
  }

  template<typename VV>
  inline bool evaluate_v(VV& myV)
  {
    //evaluate only V
    bool inAtom=false;
    const int center_idx=ei_dist->get_first_neighbor_temporal();
    if(center_idx<0) abort();
    auto& myCenter=AtomicCenters[Super2Prim[center_idx]];
    if ( ei_dist->Temp_r[center_idx] < myCenter.cutoff )
    {
      inAtom=true;
      myCenter.evaluate_v(ei_dist->Temp_r[center_idx], ei_dist->Temp_dr[center_idx], myV);
    }
    return inAtom;
  }

  template<typename VV, typename GV>
  inline bool evaluate_vgl(VV& myV, GV& myG, VV& myH)
  {
    //missing
  }

  template<typename VV, typename GV, typename HT>
  inline bool evaluate_vgh(VV& myV, GV& myG, HT& myH)
  {
    bool inAtom=false;
    const int center_idx=ei_dist->get_first_neighbor_temporal();
    if(center_idx<0) abort();
    auto& myCenter=AtomicCenters[Super2Prim[center_idx]];
    if ( ei_dist->Temp_r[center_idx] < myCenter.cutoff )
    {
      inAtom=true;
      myCenter.evaluate_vgh(ei_dist->Temp_r[center_idx], ei_dist->Temp_dr[center_idx], myV, myG, myH);
    }
    return inAtom;
  }
};

}
#endif
