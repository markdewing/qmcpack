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
  vContainer_type l_vals;
  vContainer_type r_power_minus_l;
  AtomicSplineType* MultiSpline;

  vContainer_type localV, localG, localL;

  AtomicOrbitalSoA(int Lmax, int Nb):
  Ylm(Lmax), NumBands(Nb), MultiSpline(nullptr), lmax(Lmax),
  lm_tot((Lmax+1)*(Lmax+1)), Npad(getAlignedSize<ST>(Nb))
  {
    localV.resize(Npad*lm_tot);
    localG.resize(Npad*lm_tot);
    localL.resize(Npad*lm_tot);
    r_power_minus_l.resize(lm_tot);
    l_vals.resize(lm_tot);
    for(int l=0; l<=lmax; l++)
      for(int m=-l; m<=l; m++)
        l_vals[l*(l+1)+m] = l;
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
    einspline_engine<AtomicSplineType> bigtable(MultiSpline);
    return h5f.read(bigtable,"radial_spline");
  }

  bool write_splines(hdf_archive& h5f)
  {
    bool success=true;
    success = success && h5f.write(cutoff, "cutoff_radius");
    success = success && h5f.write(spline_radius, "spline_radius");
    success = success && h5f.write(spline_npoints, "spline_npoints");
    success = success && h5f.write(lmax, "l_max");
    success = success && h5f.write(pos, "position");
    einspline_engine<AtomicSplineType> bigtable(MultiSpline);
    success = success && h5f.write(bigtable,"radial_spline");
    return success;
  }

  //evaluate only V
  template<typename VV>
  inline void evaluate_v(const ST& r, const PointType& dr, VV& myV)
  {
    if (r>0)
      Ylm.evaluateV(-dr[0]/r, -dr[1]/r, -dr[2]/r);
    else
      Ylm.evaluateV(0,0,1);
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

  //evaluate VGL
  template<typename VV, typename GV>
  inline void evaluate_vgl(const ST& r, const PointType& dr, VV& myV, GV& myG, VV& myL)
  {
    const ST rmin=0;
    ST drx, dry, drz, rhatx, rhaty, rhatz, rinv, r_reg;
    if (r>rmin)
    {
      r_reg=r;
      rinv=1.0/r;
      drx=-dr[0];
      dry=-dr[1];
      drz=-dr[2];
      rhatx=drx*rinv;
      rhaty=dry*rinv;
      rhatz=drz*rinv;
    }
    else
    {
      r_reg=0;
      rinv=0;
      drx=-dr[0];
      dry=-dr[1];
      drz=-dr[2];
      rhatx=0;
      rhaty=0;
      rhatz=1;
    }

    Ylm.evaluateVGL(drx, dry, drz);
    const ST* restrict Ylm_v=Ylm[0];
    const ST* restrict Ylm_gx=Ylm[1];
    const ST* restrict Ylm_gy=Ylm[2];
    const ST* restrict Ylm_gz=Ylm[3];

    einspline::evaluate(MultiSpline,r_reg,localV,localG,localL);

    ST* restrict g0=myG.data(0);
    ST* restrict g1=myG.data(1);
    ST* restrict g2=myG.data(2);
    CONSTEXPR ST czero(0), cone(1), chalf(0.5);
    std::fill(myV.begin(),myV.end(),czero);
    std::fill(g0,g0+Npad,czero);
    std::fill(g1,g1+Npad,czero);
    std::fill(g2,g2+Npad,czero);
    std::fill(myL.begin(),myL.end(),czero);

    if(r_reg>0)
    {
      // far from core
      r_power_minus_l[0]=cone;
      ST r_power_temp=cone;
      for(int l=1; l<=lmax; l++)
      {
        r_power_temp*=rinv;
        for(int m=-l, lm=l*l; m<=l; m++,lm++)
          r_power_minus_l[lm]=r_power_temp;
      }

      //std::cout << "debug r " << r << " rx,ry,rz : " << drx << " " << dry << " " << drz << std::endl;
      for(size_t lm=0; lm<lm_tot; lm++)
      {
        size_t offset=lm*Npad;
        const ST& l_val=l_vals[lm];
        const ST& r_power=r_power_minus_l[lm];
        const ST Ylm_rescale=Ylm_v[lm]*r_power;
        //std::cout << "debug lm " << lm << " YlmV : " << Ylm_v[lm] << " YlmG " << Ylm_gx[lm] << " " << Ylm_gy[lm] << " " << Ylm_gz[lm] << std::endl;
        for(size_t ib=0; ib<myV.size(); ib++)
        {
          // value
          const ST Vpart = Ylm_rescale*localV[offset+ib];
          myV[ib] += Vpart;

          // grad
          g0[ib] += localG[offset+ib] * rhatx * Ylm_rescale + localV[offset+ib] * Ylm_gx[lm] * r_power - l_val * rhatx * Vpart * rinv;
          g1[ib] += localG[offset+ib] * rhaty * Ylm_rescale + localV[offset+ib] * Ylm_gy[lm] * r_power - l_val * rhaty * Vpart * rinv;
          g2[ib] += localG[offset+ib] * rhatz * Ylm_rescale + localV[offset+ib] * Ylm_gz[lm] * r_power - l_val * rhatz * Vpart * rinv;

          // laplacian
          ST rhat_dot_G = ( rhatx*Ylm_gx[lm] + rhaty*Ylm_gy[lm] + rhatz*Ylm_gz[lm] ) * r_power;
          //std::cout << "debug " << Ylm_rescale << " rhat_dot_g " << rhat_dot_G << std::endl;
          myL[ib] += (localL[offset+ib] + localG[offset+ib] * 2 * rinv) * Ylm_rescale
                    + localG[offset+ib] * (rhat_dot_G - l_val * Ylm_rescale * rinv )
                    - localV[offset+ib] * l_val * rinv * (Ylm_rescale * rinv + rhat_dot_G );
          //std::cout << "debug ib " << ib << " localVGL : " << localV[offset+ib] << " " << localG[offset+ib] << " " << localL[offset+ib]
          //          << " v " << myV[ib] << " g_xyz " << g0[ib] <<" " << g1[ib] << " " << g2[ib] << " l " << myL[ib] << std::endl;
        }
      }
    }
    else
    {
      // near core, kill divergence
      for(size_t ib=0; ib<myV.size(); ib++)
      {
        // value
        myV[ib] += Ylm_v[0]*localV[ib];

        // laplacian
        myL[ib] += localL[ib] * static_cast<ST>(3) * Ylm_v[0];
      }
      if(lm_tot>0)
      {
        //std::cout << std::endl;
        for(size_t lm=1; lm<4; lm++)
        {
          size_t offset=lm*Npad;
          //std::cout << "debug lm " << lm << " YlmV : " << Ylm_v[lm] << " YlmG " << Ylm_gx[lm] << " " << Ylm_gy[lm] << " " << Ylm_gz[lm] << std::endl;
          for(size_t ib=0; ib<myV.size(); ib++)
          {
            // grad
            g0[ib] += localG[offset+ib] * Ylm_gx[lm];
            g1[ib] += localG[offset+ib] * Ylm_gy[lm];
            g2[ib] += localG[offset+ib] * Ylm_gz[lm];
            //std::cout << "debug ib " << ib << " localVGL : " << localV[offset+ib] << " " << localG[offset+ib] << " " << localL[offset+ib]
            //          << " v " << myV[ib] << " g_xyz " << g0[ib] <<" " << g1[ib] << " " << g2[ib] << " l " << myL[ib] << std::endl;
          }
        }
      }
    }
  }

  template<typename VV, typename GV, typename HT>
  void evaluate_vgh(const ST& r, const PointType& dr, VV& myV, GV& myG, HT& myH)
  {
    //Needed to do tensor product here
    APP_ABORT("AtomicOrbitalSoA::evaluate_vgh");
  }
};

/** adoptor class to match 
 *
 */
template<typename ST>
struct HybridAdoptorBase
{
  using PointType=typename AtomicOrbitalSoA<ST>::PointType;

  // atomic centers
  std::vector<AtomicOrbitalSoA<ST> > AtomicCenters;
  ///table index
  int myTableID;
  //mapping supercell to primitive cell
  std::vector<int> Super2Prim;
  // r, dr for distance table
  DistanceTableData::RealType dist_r;
  DistanceTableData::PosType dist_dr;
  // local r, dr
  ST r;
  PointType dr;

  HybridAdoptorBase() { }

  void set_info(const ParticleSet& ions, ParticleSet& els, const std::vector<int>& mapping)
  {
    myTableID=els.addTable(ions,DT_SOA);
    Super2Prim=mapping;
  }

  bool read_splines(hdf_archive& h5f)
  {
    bool success=true;
    size_t ncenter;
    int Nb;

    success = success && h5f.push("atomic_centers",false);
    success = success && h5f.read(ncenter,"number_of_centers");
    success = success && h5f.read(Nb,"number_of_bands_pad");
    if(!success) return success;
    // read splines of each center
    for(int ic=0; ic<ncenter; ic++)
    {
      std::ostringstream gname;
      gname << "center_" << ic;
      success = success && h5f.push(gname.str().c_str(),false);

      int lmax, spline_npoints;
      PointType pos;
      ST spline_radius, cutoff;
      success = success && h5f.read(cutoff, "cutoff_radius");
      success = success && h5f.read(spline_radius, "spline_radius");
      success = success && h5f.read(spline_npoints, "spline_npoints");
      success = success && h5f.read(pos, "position");
      success = success && h5f.read(lmax, "l_max");
      if(!success) return success;
      AtomicOrbitalSoA<ST> one_center(lmax, Nb);
      if(one_center.Npad!=Nb) success=false;
      one_center.set_info(pos,cutoff,spline_radius,spline_npoints);
      one_center.create_spline();
      success = success && one_center.read_splines(h5f);
      if(!success) return success;
      AtomicCenters.push_back(one_center);

      h5f.pop();
    }
    h5f.pop();
    return success;
  }

  bool write_splines(hdf_archive& h5f)
  {
    bool success=true;
    int ncenter=AtomicCenters.size();
    success = success && h5f.push("atomic_centers",true);
    success = success && h5f.write(ncenter,"number_of_centers");
    success = success && h5f.write(AtomicCenters[0].Npad,"number_of_bands_pad");
    // write splines of each center
    for(int ic=0; ic<AtomicCenters.size(); ic++)
    {
      std::ostringstream gname;
      gname << "center_" << ic;
      success = success && h5f.push(gname.str().c_str(),true);
      success = success && AtomicCenters[ic].write_splines(h5f);
      h5f.pop();
    }
    h5f.pop();
    return success;
  }

  //evaluate only V
  template<typename VV>
  inline bool evaluate_v(const ParticleSet& P, const int iat, VV& myV)
  {
    bool inAtom=false;
    const auto* ei_dist=P.DistTables[myTableID];
    const int center_idx=ei_dist->get_first_neighbor(iat, dist_r, dist_dr);
    if(center_idx<0) abort();
    r=dist_r; dr=dist_dr;
    auto& myCenter=AtomicCenters[Super2Prim[center_idx]];
    if ( r < myCenter.cutoff )
    {
      inAtom=true;
      myCenter.evaluate_v(r, dr, myV);
    }
    return inAtom;
  }

  //evaluate only VGL
  template<typename VV, typename GV>
  inline bool evaluate_vgl(const ParticleSet& P, const int iat, VV& myV, GV& myG, VV& myL)
  {
    bool inAtom=false;
    const auto* ei_dist=P.DistTables[myTableID];
    const int center_idx=ei_dist->get_first_neighbor(iat, dist_r, dist_dr);
    if(center_idx<0) abort();
    r=dist_r; dr=dist_dr;
    auto& myCenter=AtomicCenters[Super2Prim[center_idx]];
    if ( r < myCenter.cutoff )
    {
      inAtom=true;
      myCenter.evaluate_vgl(r, dr, myV, myG, myL);
    }
    return inAtom;
  }

  //evaluate only VGH
  template<typename VV, typename GV, typename HT>
  inline bool evaluate_vgh(const ParticleSet& P, const int iat, VV& myV, GV& myG, HT& myH)
  {
    bool inAtom=false;
    const auto* ei_dist=P.DistTables[myTableID];
    const int center_idx=ei_dist->get_first_neighbor(iat, dist_r, dist_dr);
    if(center_idx<0) abort();
    r=dist_r; dr=dist_dr;
    auto& myCenter=AtomicCenters[Super2Prim[center_idx]];
    if ( r < myCenter.cutoff )
    {
      inAtom=true;
      myCenter.evaluate_vgh(r, dr, myV, myG, myH);
    }
    return inAtom;
  }
};

}
#endif
