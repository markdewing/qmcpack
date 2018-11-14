////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Alfredo Correa, correaa@llnl.gov, Lawrence Livermore National Laboratory 
// Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory 
//
// File created by:
// Alfredo Correa, correaa@llnl.gov, Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////

#ifndef SPARSE_CFO_MATRIX_HPP
#define SPARSE_CFO_MATRIX_HPP

#include<array>
#include<cassert>
#include<iostream>
//#include<cstddef>  // ptrdiff_t
#include<vector>
#include<tuple>
#include<numeric>
#include<memory> 
#include<type_traits> // enable_if
#include<algorithm>
#include<utility>

#include "Configuration.h"
#include "AFQMC/Utilities/tuple_iterator.hpp"

#include "mpi.h"

namespace ma{
namespace sparse{

using size_type           = std::size_t;
//using difference_type     = std::ptrdiff_t;
//using index               = std::ptrdiff_t;

template<class Allocator>
struct null_is_root{
	null_is_root(Allocator){}
	bool root(){return true;}
	int size(){return 1;}
	int rank(){return 0;}
	void barrier() {};
};

template<
        class ValType,
        class IndxType = int,
        class IntType = size_type,
	class ValType_alloc = std::allocator<ValType>
>
class csr_matrix_ref {
	public:
	using value_type = ValType;
	using index_type = IndxType;
	using int_type = IntType;
        protected:
        using IndxType_alloc = typename ValType_alloc::template rebind<IndxType>::other;
        using IntType_alloc = typename ValType_alloc::template rebind<IntType>::other;
        using this_t = csr_matrix_ref<ValType,IndxType,IntType,ValType_alloc>;
	using ValTypePtr = typename ValType_alloc::pointer;
	using IndxTypePtr = typename IndxType_alloc::pointer;
	using IntTypePtr = typename IntType_alloc::pointer;
        size_type size1_;
        size_type size2_;
	size_type local_origin1_;
	size_type local_origin2_;
	size_type global_origin1_;
	size_type global_origin2_;
        ValTypePtr data_;
        IndxTypePtr jdata_;
        IntTypePtr pointers_begin_;
        IntTypePtr pointers_end_;
        public:
        using row_type = std::tuple<size_type,ValTypePtr,IndxTypePtr>;
	static const bool sparse = true;
	static const int dimensionality = -2;
        csr_matrix_ref(
                std::tuple<size_type, size_type> const& arr,
                std::tuple<size_type, size_type> const& local,
                std::tuple<size_type, size_type> const& global,
                ValTypePtr __data_,
                IndxTypePtr __jdata_,
                IntTypePtr __pointers_begin_,
                IntTypePtr __pointers_end_
        ) :
                size1_(std::get<0>(arr)), size2_(std::get<1>(arr)),
                local_origin1_(std::get<0>(local)), local_origin2_(std::get<1>(local)),
                global_origin1_(std::get<0>(global)), global_origin2_(std::get<1>(global)),
                data_(__data_),
                jdata_(__jdata_),
                pointers_begin_(__pointers_begin_),
                pointers_end_(__pointers_end_)
        {}
	csr_matrix_ref(this_t const& other) = delete;
	csr_matrix_ref& operator=(this_t const& other) = delete;
	// pointer movement is handled by derived classes
        csr_matrix_ref(this_t&& other) = default; //:csr_matrix_ref() { }
        csr_matrix_ref& operator=(this_t&& other) = default; 
        ~csr_matrix_ref() {}
        auto pointers_begin() const{return pointers_begin_;}
        auto pointers_end() const{return pointers_end_;}
        auto local_origin() const{return std::array<size_type, 2>{{local_origin1_,local_origin2_}};;}
        auto global_origin() const{return std::array<size_type, 2>{{global_origin1_,global_origin2_}};;}
        auto size() const{return size1_;}
	template<typename integer_type>
        auto num_elements(integer_type i) const{
		if(not pointers_begin_)  return size_type(0);
		return static_cast<size_type>(pointers_begin_[i+1]-pointers_begin_[i]);
	}
        auto num_elements() const{
		if(not pointers_begin_)  return size_type(0);
                return static_cast<size_type>(pointers_begin_[size1_]-pointers_begin_[0]);
        }
        auto num_non_zero_elements() const{
                size_type ret = 0;
                for(size_type i = 0; i != size(); ++i)
                        ret += static_cast<size_type>(pointers_end_[i] - pointers_begin_[i]);
                return ret;
        }
        auto shape() const{return std::array<size_type, 2>{{size(),size2_}};}
        auto non_zero_values_data() const{return data_;}
        auto non_zero_indices2_data() const{return jdata_;}
        auto sparse_row(int i) const{
            assert(i >= 0 && i < size1_);
            return std::make_tuple(size_type(pointers_end_[i]-pointers_begin_[i]),
                                   data_+pointers_begin_[i],
                                   jdata_+pointers_begin_[i]);
        }
	auto release_non_zero_values_data() {
		auto t = data_;
		data_ = ValTypePtr(nullptr);
		return t;
	}
	auto release_non_zero_indices2_data() {
		auto t = jdata_; 
		jdata_ = IndxTypePtr(nullptr);
		return t;
	}
	auto release_pointer_begin() {
		auto t = pointers_begin_;
		pointers_begin_ = IntTypePtr(nullptr);
		return t;
	}
	auto release_pointer_end() {
		auto t = pointers_end_;
		pointers_end_ = IntTypePtr(nullptr);
		return t;
	}
        // ??? not sure how to do this better
        auto release_shape() {
                auto t = shape();
                size1_ = size2_ = IntType(0);
                return t;
        }
        friend decltype(auto) size(this_t const& s){return s.size();}
        friend decltype(auto) shape(this_t const& s){return s.shape();}
//        friend auto index_bases(this_t const& s){return s.index_bases();}
	friend auto num_elements(this_t const& s){return s.num_elements();}
	friend auto num_non_zero_elements(this_t const& s){return s.num_non_zero_elements();}
        friend auto non_zero_values_data(this_t const& s){return s.non_zero_values_data();}
        friend auto non_zero_indices2_data(this_t const& s){return s.non_zero_indices2_data();}
        friend auto pointers_begin(this_t const& s){return s.pointers_begin();}
        friend auto pointers_end(this_t const& s){return s.pointers_end();}
};

// View of a sub matrix. Owns pointer_begin and pointer_end 
// Don't know how to make this class derive from the current csr_matrix_ref,
// so I'm just making a new class
template<
        class ValType,
        class IndxType = int,
        class IntType = int
>
class csr_matrix_view {
	public:
	using value_type = ValType;
	using index_type = IndxType;
	using int_type = IntType;
        protected:
        using IndxType_alloc = typename ValType_alloc::template rebind<IndxType>::other;
        using IntType_alloc = typename ValType_alloc::template rebind<IntType>::other;
        using this_t = csr_matrix_ref<ValType,IndxType,IntType,ValType_alloc>;
	using ValTypePtr = typename ValType_alloc::pointer;
	using IndxTypePtr = typename IndxType_alloc::pointer;
	using IntTypePtr = typename IntType_alloc::pointer;
        size_type size1_;
        size_type size2_;
	size_type local_origin1_;
	size_type local_origin2_;
	size_type global_origin1_;
	size_type global_origin2_;
        ValTypePtr data_;
        IndxTypePtr jdata_;
        IntTypePtr pointers_begin_;
        IntTypePtr pointers_end_;
        public:
        using row_type = std::tuple<size_type,ValTypePtr,IndxTypePtr>;
	static const bool sparse = true;
	static const int dimensionality = -2;
        csr_matrix_ref(
                std::tuple<size_type, size_type> const& arr,
                std::tuple<size_type, size_type> const& local,
                std::tuple<size_type, size_type> const& global,
                ValTypePtr __data_,
                IndxTypePtr __jdata_,
                IntTypePtr __pointers_begin_,
                IntTypePtr __pointers_end_
        ) :
                size1_(std::get<0>(arr)), size2_(std::get<1>(arr)),
                local_origin1_(std::get<0>(local)), local_origin2_(std::get<1>(local)),
                global_origin1_(std::get<0>(global)), global_origin2_(std::get<1>(global)),
                data_(__data_),
                jdata_(__jdata_),
                pointers_begin_(__pointers_begin_),
                pointers_end_(__pointers_end_)
        {}
	csr_matrix_ref(this_t const& other) = delete;
	csr_matrix_ref& operator=(this_t const& other) = delete;
	// pointer movement is handled by derived classes
        csr_matrix_ref(this_t&& other) = default; //:csr_matrix_ref() { }
        csr_matrix_ref& operator=(this_t&& other) = default; 
        ~csr_matrix_ref() {}
        auto pointers_begin() const{return pointers_begin_;}
        auto pointers_end() const{return pointers_end_;}
        auto local_origin() const{return std::array<size_type, 2>{{local_origin1_,local_origin2_}};;}
        auto global_origin() const{return std::array<size_type, 2>{{global_origin1_,global_origin2_}};;}
        auto size() const{return size1_;}
	template<typename integer_type>
        auto num_elements(integer_type i) const{
		if(not pointers_begin_)  return size_type(0);
		return static_cast<size_type>(pointers_begin_[i+1]-pointers_begin_[i]);
	}
        auto num_elements() const{
		if(not pointers_begin_)  return size_type(0);
                return static_cast<size_type>(pointers_begin_[size1_]-pointers_begin_[0]);
        }
        auto num_non_zero_elements() const{
                size_type ret = 0;
                for(size_type i = 0; i != size(); ++i)
                        ret += static_cast<size_type>(pointers_end_[i] - pointers_begin_[i]);
                return ret;
        }
        auto shape() const{return std::array<size_type, 2>{{size(),size2_}};}
        auto non_zero_values_data() const{return data_;}
        auto non_zero_indices2_data() const{return jdata_;}
        auto sparse_row(int i) const{
            assert(i >= 0 && i < size1_);
            return std::make_tuple(size_type(pointers_end_[i]-pointers_begin_[i]),
                                   data_+pointers_begin_[i],
                                   jdata_+pointers_begin_[i]);
        }
	auto release_non_zero_values_data() {
		auto t = data_;
		data_ = ValTypePtr(nullptr);
		return t;
	}
	auto release_non_zero_indices2_data() {
		auto t = jdata_; 
		jdata_ = IndxTypePtr(nullptr);
		return t;
	}
	auto release_pointer_begin() {
		auto t = pointers_begin_;
		pointers_begin_ = IntTypePtr(nullptr);
		return t;
	}
	auto release_pointer_end() {
		auto t = pointers_end_;
		pointers_end_ = IntTypePtr(nullptr);
		return t;
	}
        // ??? not sure how to do this better
        auto release_shape() {
                auto t = shape();
                size1_ = size2_ = IntType(0);
                return t;
        }
        friend decltype(auto) size(this_t const& s){return s.size();}
        friend decltype(auto) shape(this_t const& s){return s.shape();}
//        friend auto index_bases(this_t const& s){return s.index_bases();}
	friend auto num_elements(this_t const& s){return s.num_elements();}
	friend auto num_non_zero_elements(this_t const& s){return s.num_non_zero_elements();}
        friend auto non_zero_values_data(this_t const& s){return s.non_zero_values_data();}
        friend auto non_zero_indices2_data(this_t const& s){return s.non_zero_indices2_data();}
        friend auto pointers_begin(this_t const& s){return s.pointers_begin();}
        friend auto pointers_end(this_t const& s){return s.pointers_end();}
};

}


template<
	class ValType,
	class IndxType = int,
	class IntType = size_type,    
	class ValType_alloc = std::allocator<ValType>, 
	class IsRoot = null_is_root<ValType_alloc> 
>
class ucsr_matrix: public csr_matrix_ref<ValType,IndxType,IntType,ValType_alloc>{
	public:
	using value_type = ValType;
	using index_type = IndxType;
	using int_type = IntType;
	protected:
	using this_t = ucsr_matrix<ValType,IndxType,IntType,ValType_alloc,IsRoot>;
	using base = csr_matrix_ref<ValType,IndxType,IntType,ValType_alloc>;
	using IndxType_alloc = typename ValType_alloc::template rebind<IndxType>::other;
	using IntType_alloc = typename ValType_alloc::template rebind<IntType>::other;
	using ValTypePtr = typename ValType_alloc::pointer;
	using IndxTypePtr = typename IndxType_alloc::pointer;
	using IntTypePtr = typename IntType_alloc::pointer;
	using Valloc_ts = std::allocator_traits<ValType_alloc>; 
	using Ialloc_ts = std::allocator_traits<IndxType_alloc>; 
	using Palloc_ts = std::allocator_traits<IntType_alloc>; 
	ValType_alloc Valloc_;
	IndxType_alloc Ialloc_;
	IntType_alloc Palloc_;
        void reset() {
                IsRoot r(Valloc_);
                // calling barrier for safety right now
                r.barrier();
                if(r.root()){
                        size_type tot_sz = base::num_elements();
                        if(base::data_ && base::pointers_begin_ && base::pointers_end_) 
                                for(size_type i = 0; i != base::size1_; ++i)
                                        for(auto p = base::data_ + base::pointers_begin_[i]; 
                                                 p != base::data_ + base::pointers_end_[i]; ++p)
                                                Valloc_.destroy(std::addressof(*p));
                        if(base::jdata_ && base::pointers_begin_ && base::pointers_end_) 
                                for(size_type i = 0; i != base::size1_; ++i)
                                        for(auto p = base::jdata_ + base::pointers_begin_[i]; 
                                                 p != base::jdata_ + base::pointers_end_[i]; ++p)
                                                Ialloc_.destroy(std::addressof(*p));
                        if(base::pointers_begin_ && base::pointers_end_) {
                                for(size_type i = 0; i != base::size1_; ++i){
                                        Palloc_.destroy(std::addressof(base::pointers_begin_[i]));
                                        Palloc_.destroy(std::addressof(base::pointers_end_[i]));
                                }
                                Palloc_.destroy(std::addressof(base::pointers_begin_[base::size1_]));
                        }
                        if(base::data_)   
                                Valloc_.deallocate(base::data_, tot_sz);
                        if(base::jdata_)   
                                Ialloc_.deallocate(base::jdata_, tot_sz);
                        if(base::pointers_begin_)   
                                Palloc_.deallocate(base::pointers_begin_, base::size1_+1);
                        if(base::pointers_end_)   
                                Palloc_.deallocate(base::pointers_end_  , base::size1_);
                }
                base::size1_ = base::size2_ = 0;
                base::data_ = ValTypePtr(nullptr);
                base::jdata_ = IndxTypePtr(nullptr);
                base::pointers_begin_ = IntTypePtr(nullptr);
                base::pointers_end_ = IntTypePtr(nullptr);
                r.barrier();
        }
	public:
        static const bool sparse = true;
        static const int dimensionality = -2;
	template<typename integer_type>
	ucsr_matrix(
		std::tuple<size_type, size_type> const& arr = {0, 0}, 
		integer_type nnzpr_unique = 0,
		ValType_alloc alloc = ValType_alloc{}
	) : 
		csr_matrix_ref<ValType,IndxType,IntType,ValType_alloc>(arr,
			{0,0},
			{0,0},
			ValTypePtr(nullptr),
			IndxTypePtr(nullptr),
			IntTypePtr(nullptr),
			IntTypePtr(nullptr)),
		Valloc_(alloc), 
		Ialloc_(alloc),
		Palloc_(alloc)
	{
		base::data_ = Valloc_.allocate(std::get<0>(arr)*nnzpr_unique);
		base::jdata_ = Ialloc_.allocate(std::get<0>(arr)*nnzpr_unique);
		base::pointers_begin_ = Palloc_.allocate(std::get<0>(arr)+1);
		base::pointers_end_ = Palloc_.allocate(std::get<0>(arr));

		IsRoot r(Valloc_);
		if(r.root()){
			for(size_type i = 0; i != base::size1_; ++i){
				Palloc_ts::construct(Palloc_, std::addressof(base::pointers_begin_[i]), i*nnzpr_unique);
				Palloc_ts::construct(Palloc_, std::addressof(base::pointers_end_[i]), i*nnzpr_unique);
			}
                        Palloc_ts::construct(Palloc_, std::addressof(base::pointers_begin_[base::size1_]), base::size1_*nnzpr_unique);
		}
		r.barrier();
	}
	template<typename integer_type>
        ucsr_matrix(
                std::tuple<size_type, size_type> const& arr = {0, 0},
                std::vector<integer_type> const& nnzpr = std::vector<integer_type>(0),
                ValType_alloc alloc = ValType_alloc{}
        ) :
                csr_matrix_ref<ValType,IndxType,IntType,ValType_alloc>(arr,
			{0,0},
			{0,0},
			ValTypePtr(nullptr),
			IndxTypePtr(nullptr),
			IntTypePtr(nullptr),
			IntTypePtr(nullptr)),
                Valloc_(alloc),
                Ialloc_(alloc),
                Palloc_(alloc)
        {
		size_type sz = size_type(std::accumulate(nnzpr.begin(),nnzpr.end(),integer_type(0)));
                if(std::get<0>(arr)==0) sz=0; // no rows, no capacity
		base::data_ = Valloc_.allocate(sz);
		base::jdata_ = Ialloc_.allocate(sz);
		base::pointers_begin_ = Palloc_.allocate(std::get<0>(arr)+1);
		base::pointers_end_ = Palloc_.allocate(std::get<0>(arr));

		assert(nnzpr.size() >= base::size1_);
                IsRoot r(Valloc_);
                if(r.root()){
			IntType cnter(0);
                        for(size_type i = 0; i != base::size1_; ++i){
                                Palloc_ts::construct(Palloc_, std::addressof(base::pointers_begin_[i]), cnter); 
                                Palloc_ts::construct(Palloc_, std::addressof(base::pointers_end_[i]), cnter);
				cnter += static_cast<IntType>(nnzpr[i]); 
                        }
                        Palloc_ts::construct(Palloc_, std::addressof(base::pointers_begin_[base::size1_]), cnter);
                }
		r.barrier();
        }
	~ucsr_matrix(){
                reset();
        }
	ucsr_matrix(const this_t& other) = delete;  
	ucsr_matrix& operator=(const this_t& other) = delete;  
	ucsr_matrix(this_t&& other):ucsr_matrix({0,0},0,other.Valloc_)
	{ *this = std::move(other); } 
        // Instead of moving allocators, require they are the same right now
	ucsr_matrix& operator=(this_t&& other) {
		if(this != std::addressof(other)) {
                        if(Valloc_ != other.Valloc_ ||
                           Ialloc_ != other.Ialloc_ ||
                           Palloc_ != other.Palloc_ )
                            APP_ABORT(" Error: Can only move assign between csr_matrices with equivalent allocators. \n");
                        reset();
			base::size1_ = std::exchange(other.size1_,IntType(0));
			base::size2_ = std::exchange(other.size2_,IntType(0));
                	base::data_ = std::exchange(other.data_,ValTypePtr(nullptr));
                	base::jdata_ = std::exchange(other.jdata_,IndxTypePtr(nullptr));
                	base::pointers_begin_ = std::exchange(other.pointers_begin_,IntTypePtr(nullptr));
                	base::pointers_end_ = std::exchange(other.pointers_end_,IntTypePtr(nullptr));
		}
		return *this;
	} 
	auto getAlloc() { return Valloc_; }
        template<typename integer_type>
        void reserve(integer_type nnzpr_unique){
                if(base::size1_==0) return;
                IntType minN = IntType(base::pointers_begin_[1]-base::pointers_begin_[0]);
                for(size_type i = 0; i != base::size1_; ++i) 
                        minN = std::min(minN,base::pointers_begin_[i+1]-base::pointers_begin_[i]);
                if( static_cast<IntType>(nnzpr_unique) <= minN) 
                        return;
                this_t other({base::size1_,base::size2_},nnzpr_unique,Valloc_);
                IsRoot r(Valloc_);
                if(r.root()){
                        for(size_type i=0; i<base::size1_; i++)
			{
				size_type disp = static_cast<size_type>(base::pointers_end_[i]-
									base::pointers_begin_[i]); 
				std::copy_n(std::addressof(base::data_[base::pointers_begin_[i]]),
					    disp,
	       				    std::addressof(other.data_[other.pointers_begin_[i]]));
				std::copy_n(std::addressof(base::jdata_[base::pointers_begin_[i]]),
					    disp,
	       				    std::addressof(other.jdata_[other.pointers_begin_[i]]));
				other.pointers_end_[i] = other.pointers_begin_[i] + disp;
			}
                }
                r.barrier();
                *this = std::move(other);
        }
        template<typename integer_type>
        void reserve(std::vector<integer_type>& nnzpr){
                if(base::size1_==0) return;
                bool resz = false;
                assert(nnzpr.size() >= base::size1_);
                for(size_type i=0; i<base::size1_; i++)
                        if(static_cast<IntType>(nnzpr[i]) > 
                                        base::pointers_begin_[i+1]-base::pointers_begin_[i]) {
                                resz=true;
                                break;
                        }
                if(not resz)
                        return;
                this_t other({base::size1_,base::size2_},nnzpr,Valloc_);
                IsRoot r(Valloc_);
                if(r.root()){
                        for(size_type i=0; i<base::size1_; i++)
                        {
                                size_type disp = static_cast<size_type>(base::pointers_end_[i]-
                                                                        base::pointers_begin_[i]);
                                std::copy_n(std::addressof(base::data_[base::pointers_begin_[i]]),
                                            disp,
                                            std::addressof(other.data_[other.pointers_begin_[i]]));
                                std::copy_n(std::addressof(base::jdata_[base::pointers_begin_[i]]),
                                            disp,
                                            std::addressof(other.jdata_[other.pointers_begin_[i]]));
                                other.pointers_end_[i] = other.pointers_begin_[i] + disp;
                        }
                }
                r.barrier();
                *this = std::move(other);
        }
	template<class Pair = std::array<IndxType, 2>, class... Args>
	void emplace(Pair&& indices, Args&&... args){
		using std::get;
                assert(get<0>(indices) >= 0);
                assert(get<0>(indices) < base::size1_);
		if(base::pointers_end_[get<0>(indices)] < base::pointers_begin_[get<0>(indices)+1]) { 
			Valloc_ts::construct(Valloc_,std::addressof(base::data_[base::pointers_end_[get<0>(indices)]]), std::forward<Args>(args)...);
			Ialloc_ts::construct(Ialloc_, std::addressof(base::jdata_[base::pointers_end_[get<0>(indices)]]), get<1>(indices));
			++base::pointers_end_[get<0>(indices)];
		} else   throw std::out_of_range("row size exceeded the maximum");
	}
	protected:
	struct row_reference{
		ucsr_matrix& self_;
		IndxType i_;
		struct element_reference{
			row_reference& self_;
			IndxType j_;
			template<class TT>
			element_reference&&
			operator=(TT&& tt)&&{
				self_.self_.emplace({{self_.i_, j_}}, std::forward<TT>(tt));
				return std::move(*this);
			}
		};
		using reference = element_reference;
		reference operator[](IndxType i)&&{return reference{*this, i};}
	};
	public:
	using reference = row_reference;
	template<typename integer_type>
	reference operator[](integer_type i){return reference{*this, static_cast<IndxType>(i)};}
	ValType_alloc& getValloc() {return Valloc_;}
	IndxType_alloc& getIalloc() {return Ialloc_;}
	IntType_alloc& getPalloc() {return Palloc_;}
};

template<
        class ValType,
        class IndxType = int,
        class IntType = size_type,   
        class ValType_alloc = std::allocator<ValType>,
        class IsRoot = null_is_root<ValType_alloc>
>
class csr_matrix: public ucsr_matrix<ValType,IndxType,IntType,ValType_alloc,IsRoot> 
{
	using base = ucsr_matrix<ValType,IndxType,IntType,ValType_alloc,IsRoot>;
	using this_t = csr_matrix<ValType,IndxType,IntType,ValType_alloc,IsRoot>;
	public:
	using value_type = ValType;
	using index_type = IndxType;
	using int_type = IntType;
	using ValTypePtr = typename ValType_alloc::pointer;
	using IndxTypePtr = typename base::IndxType_alloc::pointer;
	using IntTypePtr = typename base::IntType_alloc::pointer;
        static const bool sparse = true;
        static const int dimensionality = -2;
	template<typename integer_type>
        csr_matrix(
                std::tuple<size_type, size_type> const& arr = {0, 0},
                integer_type nnzpr_unique = 0,
                ValType_alloc alloc = ValType_alloc{}
        ):base(arr,nnzpr_unique,alloc)
	{}
        template<typename integer_type>
        csr_matrix(
                std::tuple<size_type, size_type> const& arr = {0, 0},
                std::vector<integer_type>& nnzpr = std::vector<integer_type>(0),
                ValType_alloc alloc = ValType_alloc{}
        ):base(arr,nnzpr,alloc)
        {}
	csr_matrix(this_t const& ucsr) = delete;
	csr_matrix& operator=(this_t const& ucsr) = delete;
	csr_matrix(this_t&& other):csr_matrix({0,0},0,other.Valloc_) { *this = std::move(other); }
	csr_matrix(ucsr_matrix<ValType,IndxType,IntType,ValType_alloc,IsRoot>&& ucsr):
		csr_matrix({0,0},0,ucsr.getAlloc()) {
		*this = std::move(ucsr);
	}
        csr_matrix& operator=(csr_matrix<ValType,IndxType,IntType,ValType_alloc,IsRoot>&& other) {
                if(this != std::addressof(other)) {
                        if(base::Valloc_ != other.Valloc_ ||
                           base::Ialloc_ != other.Ialloc_ ||
                           base::Palloc_ != other.Palloc_ )
                            APP_ABORT(" Error: Can only move assign between csr_matrices with equivalent allocators. \n");
                        base::reset();
                        base::size1_ = std::exchange(other.size1_,IntType(0));
                        base::size2_ = std::exchange(other.size2_,IntType(0));
                        base::data_ = std::exchange(other.data_,ValTypePtr(nullptr));
                        base::jdata_ = std::exchange(other.jdata_,IndxTypePtr(nullptr));
                        base::pointers_begin_ = std::exchange(other.pointers_begin_,IntTypePtr(nullptr));
                        base::pointers_end_ = std::exchange(other.pointers_end_,IntTypePtr(nullptr));
                }
                return *this;
        }
	csr_matrix& operator=(ucsr_matrix<ValType,IndxType,IntType,ValType_alloc,IsRoot>&& other) {
                if(base::Valloc_ != other.getValloc() ||
                   base::Ialloc_ != other.getIalloc() ||
                   base::Palloc_ != other.getPalloc() )
                        APP_ABORT(" Error: Can only move assign between csr_matrices with equivalent allocators. \n");

                base::reset();
                auto shape_ = other.release_shape();
		base::size1_ = shape_[0]; 
		base::size2_ = shape_[1]; 
                base::data_ = other.release_non_zero_values_data();
                base::jdata_ = other.release_non_zero_indices2_data(); 
                base::pointers_begin_ = other.release_pointer_begin(); 
                base::pointers_end_ = other.release_pointer_end(); 
		using qmcplusplus::make_paired_iterator;
                IsRoot r(base::Valloc_);
		for(size_type p=0; p<base::size1_; p++) {
			if(p%static_cast<size_type>(r.size()) == static_cast<size_type>(r.rank())) {
				auto i1 = base::pointers_begin_[p];
				auto i2 = base::pointers_end_[p];
				std::sort(make_paired_iterator(std::addressof(base::jdata_[i1]),std::addressof(base::data_[i1])),
					  make_paired_iterator(std::addressof(base::jdata_[i2]),std::addressof(base::data_[i2])),
					  [](auto const& a, auto const& b) {
						return std::get<0>(a)<std::get<0>(b);
					  });
			}	
		}	
		r.barrier();
		return *this;
	}
	template<class Pair = std::array<IndxType, 2>, class... Args>
        void emplace(Pair&& indices, Args&&... args){
                using std::get;
                assert(get<0>(indices) >= 0);
                assert(get<0>(indices) < base::size1_);
                if(base::pointers_end_[get<0>(indices)] < base::pointers_begin_[get<0>(indices)+1]) {
			auto loc = std::lower_bound(std::addressof(base::jdata_[base::pointers_begin_[get<0>(indices)]]),
						    std::addressof(base::jdata_[base::pointers_end_[get<0>(indices)]]),
						    get<1>(indices));
			size_type disp = std::distance(std::addressof(base::jdata_[base::pointers_begin_[get<0>(indices)]]),std::addressof(*loc));
			size_type disp_ = std::distance(std::addressof(*loc),std::addressof(base::jdata_[base::pointers_end_[get<0>(indices)]]));
			if( disp_ > 0 && *loc == get<1>(indices)) { 
				// value exists, construct in place 
                        	base::Valloc_ts::construct(base::Valloc_,std::addressof(base::data_[base::pointers_begin_[get<0>(indices)] + disp]), std::forward<Args>(args)...);
			} else {
				// new value, shift back and add in correct place
				if(disp_ > 0) {
					std::move_backward(std::addressof(base::data_[base::pointers_begin_[get<0>(indices)] + disp]), 
					   std::addressof(base::data_[base::pointers_end_[get<0>(indices)]]),
					   std::addressof(base::data_[base::pointers_end_[get<0>(indices)] + 1]));
					std::move_backward(std::addressof(base::jdata_[base::pointers_begin_[get<0>(indices)] + disp]), 
                                           std::addressof(base::jdata_[base::pointers_end_[get<0>(indices)]]),
                                           std::addressof(base::jdata_[base::pointers_end_[get<0>(indices)] + 1]));
				}
                        	++base::pointers_end_[get<0>(indices)];
                        	base::Valloc_ts::construct(base::Valloc_,std::addressof(base::data_[base::pointers_begin_[get<0>(indices)] + disp]), std::forward<Args>(args)...);
                        	base::Ialloc_ts::construct(base::Ialloc_, std::addressof(base::jdata_[base::pointers_begin_[get<0>(indices)] + disp]), get<1>(indices));
			}
                } else throw std::out_of_range("row size exceeded the maximum");
        }	
        // new column index must be larger than all previous column indexes in the row
        // otherwise throws 
	template<class Pair = std::array<IndxType, 2>, class... Args>
        void emplace_back(Pair&& indices, Args&&... args){
                using std::get;
                assert(get<0>(indices) >= 0);
                assert(get<0>(indices) < base::size1_);
                if(base::pointers_end_[get<0>(indices)] < base::pointers_begin_[get<0>(indices)+1]) {
                        // if row is empty or new column index is larger than last column in row
                        if(base::pointers_begin_[get<0>(indices)] == 
                                base::pointers_end_[get<0>(indices)] or 
                                get<1>(indices) > base::jdata_[base::pointers_end_[get<0>(indices)]-1] ) { 
                        	base::Valloc_ts::construct(base::Valloc_,std::addressof(base::data_[base::pointers_end_[get<0>(indices)]]), std::forward<Args>(args)...);
                        	base::Ialloc_ts::construct(base::Ialloc_, std::addressof(base::jdata_[base::pointers_end_[get<0>(indices)]]), get<1>(indices));
                                ++base::pointers_end_[get<0>(indices)];
                        } else // otherwise throw 
                            throw std::runtime_error("inconsistent column index in emplace_back");
                            
                } else throw std::out_of_range("row size exceeded the maximum");
        }
	void remove_empty_spaces() {
		IsRoot r(base::Valloc_);
                if(r.root()){
			for(size_type i=0; i<base::size1_-1; i++) {
                                if(base::pointers_end_[i] == base::pointers_begin_[i+1]) continue;
				auto ni = static_cast<size_type>(base::pointers_end_[i+1]-base::pointers_begin_[i+1]);
				std::move(std::addressof(base::data_[base::pointers_begin_[i+1]]),
					  std::addressof(base::data_[base::pointers_end_[i+1]]),
					  std::addressof(base::data_[base::pointers_end_[i]]));
				std::move(std::addressof(base::jdata_[base::pointers_begin_[i+1]]),
					  std::addressof(base::jdata_[base::pointers_end_[i+1]]),
					  std::addressof(base::jdata_[base::pointers_end_[i]]));
				base::pointers_begin_[i+1] = base::pointers_end_[i];	
				base::pointers_end_[i+1] = base::pointers_begin_[i+1]+ni;	
			}
                }
                r.barrier();
	}
        protected:
        struct row_reference{
                this_t& self_;
                IndxType i_;
                struct element_reference{
                        row_reference& self_;
                        IndxType j_;
                        template<class TT>
                        element_reference&&
                        operator=(TT&& tt)&&{
                                self_.self_.emplace({{self_.i_, j_}}, std::forward<TT>(tt));
                                return std::move(*this);
                        }
                };
                using reference = element_reference;
		template<typename integer_type>
                reference operator[](integer_type i)&&{return reference{*this, static_cast<IndxType>(i)};}
        };

        public:
        using reference = row_reference;
	template<typename integer_type>
        reference operator[](integer_type i){return reference{*this, static_cast<IndxType>(i)};}

	using sub_matrix = csr_matrix_ref<ValType,IndxType,IntType,ValType_alloc>;
	matrix_view operator[](std::array<size_type,4>& arr) {

        template<typename IntT>
	using matrix_view = csr_matrix_view<ValType,IndxType,IntT,ValType_alloc>;
	matrix_view operator[](std::array<size_type,4>& arr) {
		// limited right now
		assert(arr[0]>=0 && arr[1] <= integer_type(base::size1_));
		assert(arr[2]>=0 && arr[3] <= integer_type(base::size2_));
		assert(arr[0] < arr[1]);
		assert(arr[2] < arr[3]);
		// just row partitions right now
		assert(arr[2]==0 && arr[3] == integer_type(base::size2_));
		
	        size_type disp = static_cast<size_type>(base::pointers_begin_[arr[0]]-base::pointers_begin_[0]);
		// Note: This depends on how the view is used/interpreted.
		// e.g. MKL takes as the pointer values ptrb[i]-ptrb[0], so data/jdata
		//      has to be shifted by disp.
		//      If the shift in ptrb is not done, (e.g. some other library),
		//      then data/jdata should be the same (unshifted) value as *this.
		//      cuSparse uses a 3-index format with ptr shift, e.g. ptrb[i]-ptrb[0]
		return matrix_view({base::size1_-(arr[1]-arr[0]),base::size2_},
			{arr[0],arr[2]},
			{0,0},
			base::data_ + disp,
			base::jdata_ + disp,
			base::pointers_begin_ + static_cast<size_type>(arr[0]),
			base::pointers_end_ + static_cast<size_type>(arr[0])
		);
	} 

};

}
}

#endif

