//-----------------------------------------------
// Copyright 2015 Yuri Pirola, Marco Previtali
// Released under the GPL
//-----------------------------------------------
#ifndef SUFFMGR_HPP
#define SUFFMGR_HPP

#include <cstdlib>
#include <cassert>

#include "simple_bitvector.hpp"
#include "bwt_interval.hpp"

class SuffMgr
{
public:
  typedef concurrent_bitvector_t bitvect_type;
  typedef bitvect_type::size_type size_type;

private:
  typedef bitvect_type::ones_enumerator_t ones_iterator_type;

public:

  explicit SuffMgr(const size_type bvs_len)
    :bvs_len_(bvs_len), bvs_{{bitvect_type(bvs_len), bitvect_type(bvs_len)}}
  {
  }  

  void swap(SuffMgr& other) {
    using std::swap;
    assert(bvs_len_ == other.bvs_len_);
    swap(bvs_, other.bvs_);
  }

  void add(const TrIntv& intv) {
    assert(begins(intv) < bvs_[0].size());
    assert(ends(intv) < bvs_[1].size());
    assert(endp(intv) < bvs_[1].size());
    bvs_[0].set(begins(intv));
    bvs_[1].set(ends(intv));
    bvs_[1].set(endp(intv));
  }

  void add(const BiIntv& intv) {
    assert(begins(intv) < bvs_[0].size());
    assert(ends(intv) < bvs_[1].size());
    bvs_[0].set(begins(intv));
    bvs_[1].set(ends(intv));
  }

  void reinitialize() {
    bvs_[0].reinitialize();
    bvs_[1].reinitialize();
  }
  
  class enumerator {

  public:
    bool has_next() const {
      return has_next_;
    }

    TrIntv next() {
      TrIntv ret_intv(intv);
      advance();
      return ret_intv;
    }


    enumerator(const enumerator& other)
      :
      last(other.last),
      iters(other.iters),
      intv(other.intv),
      begins_(begins(intv)),
      ends_(ends(intv)),
      endp_(endp(intv)),
      has_next_(other.has_next_)
    {
    }

  private:
    const size_t last;
    std::array<ones_iterator_type, 2> iters;
    TrIntv intv;
    size_t& begins_;
    size_t& ends_;
    size_t& endp_;
    bool has_next_;
    
    
    explicit enumerator(const SuffMgr& smgr)
      :
      last(smgr.bvs_[0].size()),
      iters{ { smgr.bvs_[0].ones_enumerator(), smgr.bvs_[1].ones_enumerator() } },
      intv(),
      begins_(begins(intv)),
      ends_(ends(intv)),
      endp_(endp(intv)),
      has_next_(true)
    {
      advance();
    }

    explicit enumerator(const SuffMgr& smgr, const size_t from, const size_t to)
      :
      last(to),
      iters{ { smgr.bvs_[0].ones_enumerator(from), smgr.bvs_[1].ones_enumerator(from) } },
      intv(),
      begins_(begins(intv)),
      ends_(ends(intv)),
      endp_(endp(intv)),
      has_next_(true)
    {
      advance();
      while (has_next_ && (ends_ < begins_)) { // 'from' is inside a tri-interval
        ends_= endp_;
        endp_= iters[1].next();
      }
    }

    void advance() {
      if (iters[0].has_next()) {
        begins_= iters[0].next();
        ends_= iters[1].next();
        endp_= iters[1].next();
        has_next_= begins_ < last;
      } else {
        has_next_= false;
      }
    }

    friend class SuffMgr;
  };

  enumerator enumerate() const {
    return enumerator(*this);
  }

  enumerator enumerate_in_range(const size_t from, const size_t to) const {
    return enumerator(*this, from, to);
  }



private:

  const size_type bvs_len_;
  std::array<bitvect_type, 2> bvs_;

};


#endif
