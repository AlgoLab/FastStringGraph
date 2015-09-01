//-----------------------------------------------
// Copyright 2015 Yuri Pirola, Marco Previtali
// Released under the GPL
//-----------------------------------------------
#ifndef SEGMENT_HPP
#define SEGMENT_HPP

#include <deque>

#include "bwt_interval.hpp"
#include "hits.hpp"

class Segment
{
public:
  Segment(const size_t& begin, const size_t& bs_len,
          const std::deque<BiIntv>& intvs,
          const std::deque<SmallBiIntv>& int_info)
    : intvs_(intvs), int_info_(int_info), bs_(bs_len, false), begin_(begin)
  {}
    
  Segment(const Segment& other)
    : intvs_(other.intvs_), int_info_(other.int_info_),
      bs_(other.bs_), begin_(other.begin_)
  {}

  Segment(const Segment& os,
          const std::deque<BiIntv>& intvs,
          const std::deque<SmallBiIntv>& int_info)
    : intvs_(intvs), int_info_(int_info), bs_(os.bs_), begin_(os.begin_)
  {}

  inline size_t size() const
  {
    return int_info_.size();
  }

  inline const SmallBiIntv& get_info(const size_t& idx) const
  {
    return int_info_[idx];
  }

  inline const BiIntv& get_interval(const size_t& idx) const
  {
    return intvs_[idx];
  }

  inline const size_t& begin() const
  {
    return begin_;
  }

  inline void reduce(const SmallBiIntv& p)
  {
    std::fill_n(bs_.begin() + begins(p), ends(p)+1-begins(p), true);
  }

  inline bool reduced(const SmallBiIntv& p) const
  {
    for(auto i = begins(p); i <= ends(p); ++i)
      if(!bs_[i])
        return false;
    return true;
  }

  inline void addTrue(Hits& h, const SmallBiIntv& p) const
  {
    for(auto i = begins(p); i <= ends(p); ++i)
      if (!bs_[i])
        h.add_suff(begin_ + i);
  }

private:
  std::deque<BiIntv>       intvs_;
  std::deque<SmallBiIntv>  int_info_;
  std::vector<bool>        bs_;
  const size_t             begin_;
};

#endif
