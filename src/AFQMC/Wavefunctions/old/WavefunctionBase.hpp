#ifndef QMCPLUSPLUS_AFQMC_WAVEFUNCTIONBASE_H
#define QMCPLUSPLUS_AFQMC_WAVEFUNCTIONBASE_H

#include "AFQMC/config.h"
#include<Message/MPIObjectBase.h>
#include "AFQMC/Hamiltonians/HamiltonianBase.hpp"
#include "io/hdf_archive.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Walkers/WalkerHandlerBase.hpp"

namespace qmcplusplus
{

class WavefunctionBase: public AFQMCInfo
{

  using WfnPtr = WavefunctionBase*;
  using WSetPtr = std::shared_ptr<WalkerHandlerBase>;

  public:

    WavefunctionBase(AFQMCInfo& info, afqmc::TaskGroup_& tg_, int wlk=1, int wfn=0, int dm=0):
       AFQMCInfo(info),TG(tg_),
       dm_type(dm),walker_type(wlk),wfn_type(wfn),
    {
      parallel = (TG.getNNodesPerTG()>1) || (TG.getNCoresPerTG()>1);
    }

    ~WavefunctionBase() {}

    bool isClosedShell() {return dm_type==0;}

    virtual ComplexType* getOrbMat()=0;
    virtual int getOrbSize()=0;
    virtual std::vector<ComplexType> getCiCoeff()=0;

    virtual int sizeOfInfoForDistributedPropagation() 
    {  
      APP_ABORT("WavefunctionBase::sizeOfInfoForDistributedPropagation() not implemented for this wavefunction type. \n");
      return 0;
    }

    virtual int dm_size(bool full) {
      switch(dm_type) {
        case 0: // closed-shell RHF
          return (full)?(NMO*NMO):(NAEA*NMO);
          break; 
        case 1:
          return (full)?(2*NMO*NMO):((NAEA+NAEB)*NMO);
          break; 
        case 2:
          return (full)?(4*NMO*NMO):((NAEA+NAEB)*2*NMO);
          break;
        default:
          APP_ABORT(" Error: Unknown dm_type in dm_size. \n");
          return -1;     
      }  
    }

  /************************************************************/
  // virtual routines that require implementation
  /************************************************************/
  // A. Single walker routines
    virtual void evaluateLocalEnergy(const ComplexType* SlaterMat, ComplexType& ekin, ComplexType& epot, ComplexType& ovl_alpha, ComplexType& ovl_beta, const int n=-1)=0;
    virtual void evaluateOverlap(const ComplexType* SlaterMat, ComplexType& ovl_alpha, ComplexType& ovl_beta, const int n=-1)=0;
    virtual void evaluateOneBodyMixedDensityMatrix(const ComplexType* SlaterMat, ComplexType* G, ComplexType& oa, ComplexType& ob, bool full=true)=0;
    virtual void evaluateTrialDensityMatrix(ComplexType* G, ComplexType& oa, ComplexType& ob)=0;

    // Two body operators: Not used yet
    virtual void calculateMixedMatrixElementOfTwoBodyOperators(const ComplexType* SlaterMat, const std::vector<s4D<ComplexType> >& vn, const std::vector<IndexType>& vn_indx, ValueSMSpMat&, std::vector<ComplexType>& v, const int n=-1 ) {
      APP_ABORT(" Error: calculateMixedMatrixElementOfTwoBodyOperators() not yet implemented. \n");
    }

  // B. Walker set routines
  protected:
    // reimplement in derived class for good performance
    virtual void dist_evaluateLocalEnergy(WSetPtr wset, std::vector<ComplexType>& oa, std::vector<ComplexType>& ob, std::vector<ComplexType>& eloc, const int n=-1)
    {
      int nw = wset->numWalkers(true);
      if(nw==0) return;
      ComplexType ekin,epot,ovlp_a,ovlp_b;
      std::fill(oa.begin(),oa.end(),ComplexType(0.0,0.0));
      std::fill(ob.begin(),ob.end(),ComplexType(0.0,0.0));
      std::fill(eloc.begin(),eloc.end(),ComplexType(0.0,0.0));
      for(int i=0,cnt=0; i<nw; i++,cnt++) {
        if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue;
        if(cnt%TG.getNCoresPerTG() == TG.getLocalTGRank()) {
          evaluateLocalEnergy(wset->getSM(i),ekin,epot,oa[i],ob[i],n);
          eloc[i] = ekin+epot;
        }
      }
      TG.TG_local().all_reduce_in_place_n(oa.data(),oa.size(),std::plus<>());
      TG.TG_local().all_reduce_in_place_n(ob.data(),ob.size(),std::plus<>());
      TG.TG_local().all_reduce_in_place_n(eloc.data(),eloc.size(),std::plus<>());
    }

  public:
    // important to use multi_array::array_view here for bound safety
    virtual void get_vbias(const ComplexType* G, ComplexType* v, int nW)=0;

    virtual void get_vHS(const ComplexType* G, ComplexType* v, int nW)=0;

    virtual void evaluate_G_for_vbias(WSetPtr wset, ComplexSMVector* buf, int wlksz, int gfoffset, bool transposed, bool full=true)=0; 

    /************************************************************/
    // virtual routines with default implementations that should always work.
    // can be re-implemented for optimal/different behavior
    /************************************************************/
    // Right now it assumes a spin restricted operator which is applied to the same spin sectors
    virtual void calculateMeanFieldMatrixElementOfOneBodyOperators(ValueSMSpMat& vn, std::vector<ComplexType>& v, bool transposed=false) {

      local_buff.resize(dm_size(true));
      ComplexType oa,ob;
      evaluateTrialDensityMatrix(local_buff.data(),oa,ob);
      contractDensityMatrixWithOneBodyOperators(local_buff.data(),vn,v,transposed,true);

    }

    virtual void calculateMixedMatrixElementOfOneBodyOperators(const ComplexType* SlaterMat, ValueSMSpMat& vn, std::vector<ComplexType>& v, bool transposed=false, bool full=true, const int n=-1) {

      local_buff.resize(dm_size(full));
      ComplexType oa,ob;
      evaluateOneBodyMixedDensityMatrix(SlaterMat,local_buff.data(),oa,ob,full);
      contractDensityMatrixWithOneBodyOperators(local_buff.data(),vn,v,transposed,full);

    }

    virtual void contractDensityMatrixWithOneBodyOperators(const ComplexType* G, ValueSMSpMat& vn, std::vector<ComplexType>& v, bool transposed=false, bool full=true)
    {

      const ValueType one = ValueType(1.0);
      const ValueType two = ValueType(2.0);
      const ValueType zero = ValueType(0.0);  
      if(dm_type==0) {
        if(transposed) {
          assert(dm_size(full) == vn.rows());
          SparseMatrixOperators::product_SpMatTV(vn.rows(),vn.cols(),two,vn.values(),vn.column_data(),vn.row_index(),G,zero,v.data());
        } else {
          assert(dm_size(full) == vn.cols());
          SparseMatrixOperators::product_SpMatV(vn.rows(),vn.cols(),two,vn.values(),vn.column_data(),vn.row_index(),G,zero,v.data());
        }
        return;
      } else if(dm_type==1 || dm_type==2) {
        if(transposed) {
          if(full) assert(dm_size(full) == 2*dm_type*vn.rows());
          else assert(dm_size(full) == vn.rows());
          SparseMatrixOperators::product_SpMatTV(vn.rows(),vn.cols(),one,vn.values(),vn.column_data(),vn.row_index(),G,zero,v.data());
          if(full)  
            SparseMatrixOperators::product_SpMatTV(vn.rows(),vn.cols(),one,vn.values(),vn.column_data(),vn.row_index(),G+vn.rows(),one,v.data());
        } else {
          if(full) assert(dm_size(full) == 2*dm_type*vn.cols());
          else assert(dm_size(full) == vn.cols());
          SparseMatrixOperators::product_SpMatV(vn.rows(),vn.cols(),one,vn.values(),vn.column_data(),vn.row_index(),G,zero,v.data());
          if(full)  
            SparseMatrixOperators::product_SpMatV(vn.rows(),vn.cols(),one,vn.values(),vn.column_data(),vn.row_index(),G+vn.cols(),one,v.data());
        }
        return;
      } else {
        APP_ABORT(" Error: Unknown dm_type in contractDensityMatrixWithOneBodyOperators(). \n");
      }              

    }

    // use multi_array later
    virtual void evaluateOneBodyMixedDensityMatrix(WSetPtr wset, ComplexSMVector* buf, int wlksz, int gfoffset, bool transposed, bool full=true) {
      ComplexType oa,ob;
      int nw = wset->numWalkers(true), cnt=0;
      int nw0 = wset->numWalkers(false);
      int sz = 2 + dm_size(full); 
      int wstride = transposed?1:nw0;
      int ax = transposed?wlksz:1;
      // in Buff:  {..., oa, ob, DM, ...} for each walker 
      local_buff.resize(sz);
      for(int i=0; i<nw; i++) {
        if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue;
        if(cnt%ncores_per_TG == core_rank ) {
          evaluateOneBodyMixedDensityMatrix(wset->getSM(i),local_buff.data()+2,local_buff[0],local_buff[1],full);
          //{
          //  boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*(buf->getMutex()));
            BLAS::copy(sz,local_buff.data(),1,buf->values() + ax*cnt + wstride*gfoffset ,wstride);
          //}
        }
        ++cnt;
      }    
    }


    /************************************************************/
    // Evaluation routines that don't need to be implemented 
    /************************************************************/

    void evaluateOneBodyMixedDensityMatrix(const ComplexType* SlaterMat, ComplexMatrix& G, ComplexType& oa, ComplexType& ob, bool full=true) {
        evaluateOneBodyMixedDensityMatrix(SlaterMat,G.data(),oa,ob,full);
    }

    void evaluateLocalEnergy( WSetPtr wset, const int n=-1)
    {
      int nw = wset->numWalkers(true);
      ova_base.resize(nw);
      ovb_base.resize(nw);
      eloc_base.resize(nw);
      if(!parallel) {
        serial_evaluateLocalEnergy(wset,ova_base,ovb_base,eloc_base,n);
      } else {
        dist_evaluateLocalEnergy(wset,ova_base,ovb_base,eloc_base,n);
        TG.local_barrier();
      }
      if(TG.getLocalTGRank() == 0) {   
        for(int i=0; i<nw; i++) 
          wset->setWalker(i,eloc_base[i],ova_base[i],ovb_base[i]);
      }
      TG.local_barrier();
    }

    void evaluateLocalEnergy( WSetPtr wset, std::vector<ComplexType>& oa, std::vector<ComplexType>& ob, std::vector<ComplexType>& eloc, const int n=-1)
    {
      int nw = wset->numWalkers(true);
      oa.resize(nw);
      ob.resize(nw);
      eloc.resize(nw);
      if(parallel) {
        dist_evaluateLocalEnergy(wset,oa,ob,eloc,n);
      } else {
        serial_evaluateLocalEnergy(wset,oa,ob,eloc,n);
      }  
      TG.local_barrier();
    }

    void verifyWalkerData(WSetPtr wset, const int n=-1)
    {
      if(parallel)
        dist_verifyWalkerData(wset,n);
      else
        serial_verifyWalkerData(wset,n);
    }


    void evaluateOverlap(WSetPtr wset, const int n=-1)
    {
      int nw = wset->numWalkers(true);
      ova_base.resize(nw);
      ovb_base.resize(nw);
      if(parallel) { 
        dist_evaluateOverlap(wset,ova_base,ovb_base,n);
        TG.local_barrier();
      } else {
        serial_evaluateOverlap(wset,ova_base,ovb_base,n);
      }
      if(TG.getLocalTGRank() == 0) {  
        for(int i=0; i<nw; i++)  
          wset->setOvlp(i,ova_base[i],ovb_base[i]);
      }
      TG.local_barrier();
    }

    void evaluateOverlap(WSetPtr wset, std::vector<ComplexType>& oa, std::vector<ComplexType>& ob, const int n=-1)
    {
      int nw = wset->numWalkers(true);
      oa.resize(nw);
      ob.resize(nw);
      if(parallel)
        dist_evaluateOverlap(wset,oa,ob,n);
      else
        serial_evaluateOverlap(wset,oa,ob,n);
    }

  protected:

    std::vector<ComplexType> ova_base,ovb_base,eloc_base,local_buff;

    afqmc::TaskGroup_& TG; 

    bool parallel;
    // in both cases below: closed_shell=0, UHF/ROHF=1, GHF=2
    int walker_type;
    int wfn_type;
    // dm_type is the DM type, which should be the largest of walker and wfn types
    int dm_type;

    // used to identify the current step.
    // The main purpose of this is to tellthe different
    // wavefunction objects whether we are in the same step
    // or not. This will allow us to reuse information already
    // calculated in a previous section of the current step.
    // e.g. Not to recalculate density matrices if we are redoing
    // a local energy calculation on the same step. 
    // Specially useful for multideterminant calculations
    int time_stamp;

    virtual void dist_evaluateOverlap(WSetPtr wset, std::vector<ComplexType>& oa, std::vector<ComplexType>& ob, const int n=-1)
    {
      int nw = wset->numWalkers(true);
      if(nw==0) return;
      ComplexType ovlp_a,ovlp_b;
      std::fill(oa.begin(),oa.end(),ComplexType(0.0,0.0));
      std::fill(ob.begin(),ob.end(),ComplexType(0.0,0.0));
      for(int i=0, cnt=0; i<nw; i++, cnt++) {
        if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue;
        if(cnt%TG.getNCoresPerTG() == TG.getLocalTGRank()) 
          evaluateOverlap(wset->getSM(i),oa[i],ob[i],n);
      }
      TG.TG_local().all_reduce_in_place_n(oa.data(),oa.size(),std::plus<>());
      TG.TG_local().all_reduce_in_place_n(ob.data(),ob.size(),std::plus<>());
    } 

    // no need to reimplement this in derived class
    void serial_evaluateOverlap(WSetPtr wset, std::vector<ComplexType>& oa, std::vector<ComplexType>& ob, const int n=-1)
    { 
      int nw = wset->numWalkers(true);
      if(nw==0) return;
      ComplexType ovlp_a,ovlp_b;
      for(int i=0; i<nw; i++) {
        if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue;
        evaluateOverlap(wset->getSM(i),oa[i],ob[i],n);
      }
    }

    // no need to reimplement this in derived class
    void serial_evaluateLocalEnergy(WSetPtr wset, std::vector<ComplexType>& oa, std::vector<ComplexType>& ob, std::vector<ComplexType>& eloc, const int n=-1)
    { 
      int nw = wset->numWalkers(true);
      if(nw==0) return;
      ComplexType ekin,epot,ovlp_a,ovlp_b;   
      for(int i=0; i<nw; i++) {
        if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue;
        evaluateLocalEnergy(wset->getSM(i),ekin,epot,oa[i],ob[i],n);
        eloc[i] = ekin+epot;
      }
    }

    virtual void dist_verifyWalkerData(WSetPtr wset, const int n=-1)
    {  
      serial_verifyWalkerData(wset,n);
    }  


    // right now just dumping messages to app_error
    void serial_verifyWalkerData(WSetPtr wset, const int n=-1)
    {
      int nw = wset->numWalkers(true);
      if(nw==0) return;
      ComplexType ekin,epot,ovlp_a,ovlp_b;
      for(int i=0; i<nw; i++) {
        if(!wset->isAlive(i) || std::abs(wset->getWeight(i)) <= 1e-6) continue;
        evaluateLocalEnergy(wset->getSM(i),ekin,epot,ovlp_a,ovlp_b,n);
        if(std::abs(wset->getEloc(i)-ekin-epot)>1e-8) app_error()<<" diff in verifyWalkerData: eloc (old-new): "<<wset->getEloc(i)<<" "<<(ekin+epot)<<"\n";
        if(std::abs(wset->getOvlpAlpha(i)-ovlp_a)>1e-8) app_error()<<" diff in verifyWalkerData: ovlp_a (old-new): "<<wset->getOvlpAlpha(i)<<" "<<ovlp_a<<"\n";
        if(std::abs(wset->getOvlpBeta(i)-ovlp_b)>1e-8) app_error()<<" diff in verifyWalkerData: ovlp_b (old-new): "<<wset->getOvlpBeta(i)<<" "<<ovlp_b<<"\n";
      }
    }

};
}

#endif
