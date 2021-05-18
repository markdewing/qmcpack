//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_LCORBITALSETBUILDER_H
#define QMCPLUSPLUS_LCORBITALSETBUILDER_H

#include <vector>
#include "Configuration.h"
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"
#include "QMCWaveFunctions/SPOSetBuilderFactory.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "QMCWaveFunctions/Fermion/MultiSlaterDeterminant.h"
#include "QMCWaveFunctions/Fermion/MultiSlaterDeterminantFast.h"
#include "QMCWaveFunctions/Fermion/ci_configuration.h"
#include "QMCWaveFunctions/Fermion/ci_configuration2.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include "QMCWaveFunctions/Fermion/BackflowBuilder.h"
namespace qmcplusplus
{
/** derived class from WaveFunctionComponentBuilder
 *
 * Builder SlaterDeterminant with LCOrbitalSet
 */
class SlaterDetBuilder : public WaveFunctionComponentBuilder
{
public:
  typedef SlaterDet SlaterDeterminant_t;
  typedef MultiSlaterDeterminant MultiSlaterDeterminant_t;
  /** constructor
   * \param els reference to the electrons
   * \param psi reference to the wavefunction
   * \param ions reference to the ions
   */
  SlaterDetBuilder(Communicate* comm,
                   SPOSetBuilderFactory& factory,
                   ParticleSet& els,
                   TrialWaveFunction& psi,
                   PtclPoolType& psets);

  /** initialize the Antisymmetric wave function for electrons
   *@param cur the current xml node
   *
   */
  WaveFunctionComponent* buildComponent(xmlNodePtr cur) override;

private:
  /// reference to the sposet_builder_factory, should be const once the legacy input style is removed
  SPOSetBuilderFactory& sposet_builder_factory_;
  ///reference to TrialWaveFunction, should go away as the CUDA code.
  TrialWaveFunction& targetPsi;
  ///reference to a PtclPoolType
  PtclPoolType& ptclPool;
  SlaterDeterminant_t* slaterdet_0;
  MultiSlaterDeterminant_t* multislaterdet_0;
  MultiSlaterDeterminantFast* multislaterdetfast_0;

  bool UseBackflow;
  BackflowTransformation* BFTrans;

  /** process a determinant element
   * @param cur xml node
   * @param firstIndex index of the determinant
   * @return firstIndex+number of orbitals
   */
  bool putDeterminant(xmlNodePtr cur, int firstIndex);

  bool createMSD(MultiSlaterDeterminant* multiSD, xmlNodePtr cur);

  bool createMSDFast(std::vector<std::unique_ptr<MultiDiracDeterminant>>& Dets,
                     std::vector<std::vector<size_t>>& C2node,
                     std::vector<ValueType>& C,
                     std::vector<ValueType>& CSFcoeff,
                     std::vector<size_t>& DetsPerCSF,
                     std::vector<RealType>& CSFexpansion,
                     bool& usingCSF,
                     opt_variables_type& myVars,
                     bool& Optimizable,
                     bool& CI_Optimizable,
                     xmlNodePtr cur);


  bool readDetList(xmlNodePtr cur,
                   std::vector<std::vector<ci_configuration>>& uniqueConfgs,
                   std::vector<std::vector<size_t>>& C2nodes,
                   std::vector<std::string>& CItags,
                   std::vector<ValueType>& coeff,
                   bool& optimizeCI,
                   std::vector<int>& nptcls,
                   std::vector<ValueType>& CSFcoeff,
                   std::vector<size_t>& DetsPerCSF,
                   std::vector<RealType>& CSFexpansion,
                   bool& usingCSF);

  bool readDetListH5(xmlNodePtr cur,
                     std::vector<std::vector<ci_configuration>>& uniqueConfgs,
                     std::vector<std::vector<size_t>>& C2nodes,
                     std::vector<std::string>& CItags,
                     std::vector<ValueType>& coeff,
                     bool& optimizeCI,
                     std::vector<int>& nptcls);

  // clang-format off
  template<typename VT,
           std::enable_if_t<(std::is_same<VT, ValueType>::value) &&
                            (std::is_floating_point<VT>::value), int> = 0>
  void readCoeffs(hdf_archive& hin, std::vector<VT>& ci_coeff, size_t n_dets,int ext_level)
  {
    ///Some converters store determinant coeffs in Coeff for the ground state
    ///and some store it in Coeff_0.  The Nth excited state in is in Coeff_N.
    std::string extVar;
    extVar="Coeff_"+std::to_string(ext_level);

    if (!hin.readEntry(ci_coeff,extVar))
      if (ext_level != 0 || !hin.readEntry(ci_coeff,"Coeff"))
          APP_ABORT("Could not read CI coefficients from HDF5");

  }
  template<typename VT,
           std::enable_if_t<(std::is_same<VT, ValueType>::value) &&
                            (std::is_same<VT, std::complex<typename VT::value_type>>::value), int> = 0>
  void readCoeffs(hdf_archive& hin, std::vector<VT>& ci_coeff, size_t n_dets,int ext_level)
  {
    std::string extVar;
    std::vector<double> CIcoeff_real;
    std::vector<double> CIcoeff_imag;
    CIcoeff_imag.resize(n_dets);
    CIcoeff_real.resize(n_dets);
    fill(CIcoeff_imag.begin(), CIcoeff_imag.end(), 0.0);
    ///Determinant coeffs are stored in Coeff_N where N is Nth excited state.
    ///The Ground State is stored in Coeff or Coeff_0.

    std::string ext_var;
    extVar="Coeff_"+std::to_string(ext_level);


    if(!hin.readEntry(CIcoeff_real, extVar))
      if (ext_level != 0 || !hin.readEntry(CIcoeff_real, "Coeff"))
        APP_ABORT("Could not read CI coefficients from HDF5")

    extVar=extVar+"_imag";
    if(!hin.readEntry(CIcoeff_imag, extVar))
      if (ext_level != 0 || !hin.readEntry(CIcoeff_imag, "Coeff_imag"))
        app_log() << "Coeff_imag not found in h5. Set to zero." << std::endl;
 
    for (size_t i = 0; i < n_dets; i++)
      ci_coeff[i] = VT(CIcoeff_real[i], CIcoeff_imag[i]);
  }
  // clang-format on
};

} // namespace qmcplusplus
#endif
