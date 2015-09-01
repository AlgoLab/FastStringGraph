//-----------------------------------------------
// Copyright 2015 Yuri Pirola, Marco Previtali
// Released under the GPL
//-----------------------------------------------
#ifndef SIMPLE_BITVECTOR_HPP
#define SIMPLE_BITVECTOR_HPP

#include <cstdint>
#include <algorithm>
#include <tbb/parallel_for.h>

#include "vector_locker.hpp"

template <class Locker>
class generic_simple_bitvector_t {
public:
  typedef size_t size_type;

private:
  typedef uint64_t data_t;
  const static uint8_t word_size= 8*sizeof(data_t);
  const static uint8_t shift_amount= log_2_t<word_size>::value;
  const static data_t low_mask= (word_size-1);

  const size_type size_;
  data_t* data_;
  Locker locker;

  constexpr static size_type to_word_size(const size_type bit_size) {
    return (bit_size>>shift_amount)+((bit_size & low_mask)? 1 : 0);
  }

public:

  explicit generic_simple_bitvector_t(const size_type size)
    :size_(size),
     data_(new data_t[to_word_size(size)]),
     locker(to_word_size(size))
  {
    reinitialize();
  }

  generic_simple_bitvector_t(const generic_simple_bitvector_t&) = delete;
  generic_simple_bitvector_t(generic_simple_bitvector_t&&) = default;

  ~generic_simple_bitvector_t() {
    //if (data_ != nullptr)
    delete [] data_;
  }

  inline void set(const size_type idx) {
    const size_type arr_pos = idx >> shift_amount;
    const size_type word_mask = 1ull << (idx & low_mask);
    typename Locker::scoped_lock lock(locker.get_mutex_for_el(arr_pos));
    data_[arr_pos] |= word_mask;
  }

  inline size_type size() const {
    return size_;
  }

  inline void reinitialize() {
    //// Serial version:
    // std::fill_n(data_, to_word_size(size_), 0);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, to_word_size(size_)),
                      [=] (const tbb::blocked_range<size_t>& r) {
                        std::fill_n(data_ + r.begin(), r.end() - r.begin(), 0);
                      } );
  }

  inline void swap(generic_simple_bitvector_t<Locker>& other) {
    using std::swap;
    swap(data_, other.data_);
  }

  class ones_enumerator_t {
  private:
    const generic_simple_bitvector_t& bv_;
    const data_t* data;
    data_t cwdata;
    size_type target_i_pos;
    const size_type end_target_i_pos;
    size_type next_aligned_target_i_pos;

    explicit ones_enumerator_t(const generic_simple_bitvector_t& bv,
                               const size_type from= 0)
      : bv_(bv),
        data(bv.data_ + (from >> shift_amount)),
        cwdata((*data) >> (from & low_mask)),
        target_i_pos(from),
        end_target_i_pos(bv.size_),
        next_aligned_target_i_pos((from & ~low_mask) + word_size)
    {
      advance();
    }

    inline void advance();

  public:

    ones_enumerator_t(const ones_enumerator_t&) = default;
    ones_enumerator_t(ones_enumerator_t&&) = default;

    bool has_next() const {
      return target_i_pos <= end_target_i_pos;
    }

    size_type next() {
      const size_type ret= target_i_pos-1;
      advance();
      return ret;
    }

    friend class generic_simple_bitvector_t;

  };

  ones_enumerator_t ones_enumerator(const size_type from= 0) const {
    return ones_enumerator_t(*this, from);
  }


};


template <class Locker>
inline void
swap(generic_simple_bitvector_t<Locker>& bv1,
     generic_simple_bitvector_t<Locker>& bv2) {
  bv1.swap(bv2);
}

// Used for ffsll
#include <cstring>

template <class Locker>
inline void
generic_simple_bitvector_t<Locker>::ones_enumerator_t::advance() {
  if (cwdata == 0) { // no more ones in current word
    ++data;
    target_i_pos = next_aligned_target_i_pos;
    while ((target_i_pos < end_target_i_pos) && (*(data)==0)) {
      ++data;
      target_i_pos += word_size;
    }
    next_aligned_target_i_pos = target_i_pos+word_size;
    if (target_i_pos < end_target_i_pos)
      cwdata= *data;
    else
      target_i_pos= end_target_i_pos+1;
  }
  if (target_i_pos < end_target_i_pos) {
    const size_t first_i_pos= static_cast<size_t>(ffsll(static_cast<long long>(cwdata)))-1;
    //// second version, more portable
    //const size_t first_i_pos= bits::lo(cwdata);
    cwdata >>= first_i_pos;
    cwdata >>= 1;
    target_i_pos += first_i_pos+1;
  }
}

using simple_bitvector_t = generic_simple_bitvector_t<null_vector_locker_t>;
using concurrent_bitvector_t = generic_simple_bitvector_t<spin_vector_locker_t>;

#endif
