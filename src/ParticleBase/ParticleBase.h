//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef OHMMS_PARTICLEBASE_H
#define OHMMS_PARTICLEBASE_H

#include <string>
#include <vector>
#include <map>
#include "ParticleTags.h"

namespace qmcplusplus
{
/**@file ParticleBase.h
 *@brief Declaration of ParticleBase<PT>
 */

/** Container class for particle attributes.
 *
 * The template parameter PT is a trait class which provides typedefs.
 * See, ParticleTraits.h
 */
template<class PT>
class ParticleBase: public ParticleTags
{

public:

  typedef typename PT::ParticleLayout_t      ParticleLayout_t;
  typedef typename PT::Index_t               Index_t;
  typedef typename PT::Scalar_t              Scalar_t;
  typedef typename PT::SingleParticleIndex_t SingleParticleIndex_t;
  typedef typename PT::SingleParticlePos_t   SingleParticlePos_t;
  typedef typename PT::Tensor_t              Tensor_t;
  typedef typename PT::ParticleValue_t       ParticleValue_t;

  //@{containers
  //! Three types of particle attributes are available.
  //! Integer quantity associated with each atom, e.g., ID, GroupID
  typedef typename PT::ParticleIndex_t       ParticleIndex_t;

  //! Scalar quantity associated with each atom, e.g., Local Energy
  typedef typename PT::ParticleScalar_t      ParticleScalar_t;

  //! D-dimensional vector quantity associated with each atom, e.g., R
  typedef typename PT::ParticlePos_t         ParticlePos_t;

  //! D-dimensional tensor quantity associated with each atom, e.g., R
  typedef typename PT::ParticleTensor_t      ParticleTensor_t;

  //! D-dimensional real/complex quantity associated with each atom, e.g., R
  typedef typename PT::ParticleGradient_t     ParticleGradient_t;

  //! D-dimensional real/complex quantity associated with each atom, e.g., R
  typedef typename PT::ParticleLaplacian_t     ParticleLaplacian_t;
  //@}

  //! Mapping between an attribute name and particle attribute
  typedef std::map<std::string,OhmmsObject*>::iterator  PAListIterator;


  /*! \note Lattice and three particle attributes, ID, GroupID, and R are
   * public. Encapsulation of ParticleBase is realized by making
   * a ParticleBase a protected member of PropagatorBase.
   * Modification of particle attributes is only allowed by PropagatorBase
   * and its derived classes.
   */
  //!< ParticleLayout
  ParticleLayout_t  Lattice, PrimitiveLattice;
  //!< unique, persistent ID for each particle
  ParticleIndex_t   ID;
  //!< Species ID
  ParticleIndex_t   GroupID;
  //!< Position
  ParticlePos_t    R;
  //!< Instantaneous Position in unit used by layout/nnlist engines, not added to myAttribList
  ParticlePos_t    curR;
  //!< Default constructor
  ParticleBase();

  ParticleBase(const ParticleBase<PT>& P): Counter(0), LocalNum(0), GlobalNum(0)
  {
    initBase();
    assign(P);
  }

  virtual ~ParticleBase();

  ///return a type id: one of the enum values
  inline int getAttribType(const std::string& tname)
  {
    return AttribTypeMap[tname];
  }

  inline int getNumAttrib() const
  {
    return AttribList.size();
  }
  inline PAListIterator first_attrib()
  {
    return AttribList.begin();
  }
  inline PAListIterator last_attrib()
  {
    return AttribList.end();
  }

  bool hasAttrib(const std::string& attrib_name);

  int addAttribute(const std::string& tname, const std::string& oname);
  int addAttribute(ParticleIndex_t& pa);
  int addAttribute(ParticleScalar_t& pa);
  int addAttribute(ParticlePos_t& pa);
  int addAttribute(ParticleTensor_t& pa);

  ParticleIndex_t*  getIndexAttrib(const std::string& aname);
  ParticleScalar_t* getScalarAttrib(const std::string& aname);
  ParticlePos_t*    getVectorAttrib(const std::string& aname);
  ParticleTensor_t* getTensorAttrib(const std::string& aname);

  inline const ParticleIndex_t*  getIndexAttrib(int i) const
  {
    return INDEX[i];
  }
  inline const ParticleScalar_t* getScalarAttrib(int i) const
  {
    return VAL[i];
  }
  inline const ParticlePos_t*    getVectorAttrib(int i) const
  {
    return POS[i];
  }
  inline const ParticleTensor_t* getTensorAttrib(int i) const
  {
    return TENZOR[i];
  }
  inline ParticleIndex_t*  getIndexAttrib(int i)
  {
    return INDEX[i];
  }
  inline ParticleScalar_t* getScalarAttrib(int i)
  {
    return VAL[i];
  }
  inline ParticlePos_t*    getVectorAttrib(int i)
  {
    return POS[i];
  }
  inline ParticleTensor_t* getTensorAttrib(int i)
  {
    return TENZOR[i];
  }

#if defined(MIXED_PRECISION) || defined(QMC_COMPLEX)
  int addAttribute(ParticleGradient_t& pa);
#endif
#if defined(QMC_COMPLEX)
  int addAttribute(ParticleLaplacian_t& pa);
#endif

  void createBase(unsigned m);
  void createBase(const std::vector<int>& agroup);

  void resize(unsigned m);
  void clear();

  virtual void assign(const ParticleBase<PT>& ptclin)
  {
    resize(ptclin.getLocalNum());
    Lattice = ptclin.Lattice;
    PrimitiveLattice = ptclin.PrimitiveLattice;
    R.InUnit = ptclin.R.InUnit;
    R = ptclin.R;
    ID = ptclin.ID;
    GroupID = ptclin.GroupID;
    if(ptclin.SubPtcl.size())
    {
      SubPtcl.resize(ptclin.SubPtcl.size());
      SubPtcl =ptclin.SubPtcl;
    }
  }

  inline int getLocalNum() const
  {
    return LocalNum;
  }
  inline int getTotalNum() const
  {
    return GlobalNum;
  }

  ///return the number of groups
  inline int groups() const
  {
    return SubPtcl.size()-1;
  }

  ///return the first index of a group i
  inline int first(int igroup) const
  {
    return SubPtcl[igroup];
  }

  ///return the last index of a group i
  inline int last(int igroup) const
  {
    return SubPtcl[igroup+1];
  }

  /// Returns the current counter
  inline int  current() const
  {
    return Counter;
  }

  /// Increments the counter
  inline void advance()
  {
    Counter++;
  }

  /// Sets the counter
  inline void setCounter(int i = 0)
  {
    Counter = i;
  }

protected:

  void initBase();

  //!< Internal counter
  int Counter;
  int LocalNum;
  int GlobalNum;

  ///array to handle a group of distinct particles per species
  ParticleIndex_t             SubPtcl;

  std::map<std::string,int>             AttribTypeMap;
  std::map<std::string,int>             Name2Index;
  std::map<std::string,OhmmsObject*>    AttribList;
  std::vector<ParticleIndex_t*>    INDEX;
  std::vector<ParticleScalar_t*>   VAL;
  std::vector<ParticlePos_t*>      POS;
  std::vector<ParticleTensor_t*>   TENZOR;
#if defined(MIXED_PRECISION) || defined(QMC_COMPLEX)
  std::vector<ParticleGradient_t*>   GRADS;
#endif
#if defined(QMC_COMPLEX)
  std::vector<ParticleLaplacian_t*>  LAPS;
#endif
  std::vector<OhmmsObject*>        myAttribList;
};

template<class T, unsigned D>
inline T* get_first_address(ParticleAttrib<TinyVector<T,D> >& a)
{
  return &(a[0][0]);
}

template<class T, unsigned D>
inline T* get_last_address(ParticleAttrib<TinyVector<T,D> >& a)
{
  return &(a[0][0])+D*a.size();
}

}
#include "ParticleBase/ParticleBase.cpp"
#endif

