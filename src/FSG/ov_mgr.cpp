//-----------------------------------------------
// Copyright 2015 Yuri Pirola, Marco Previtali
// Released under the GPL
//-----------------------------------------------
#include "ov_mgr.hpp"

#include <stack>
#include <deque>
#include <cassert>
#include <limits>

// Overlap Manager
OvMgr::OvMgr(const size_type bvs_len, const size_t perm_len)
  : // bvs_len_(bvs_len),
    bvs_ {{ bitvect_type(bvs_len), bitvect_type(bvs_len) }},
    perm_len_(perm_len),
    perm_ {{ perm_type(perm_len), perm_type(perm_len) }}
{
}

// Visitor
segm_start_t* OvMgrVisitor::operator()()
{
  if (startp_ >= om_.perm_len_) {
    return NULL;
  }
  OvMgr::perm_type::size_type perm_beg = startp_;
  OvMgr::bitvect_type::size_type intvs_start = starts_.next();
  while(om_.perm_[0][perm_beg] == 0) ++perm_beg;
  OvMgr::perm_type::size_type i = perm_beg;
  size_t n = om_.perm_[0][i]; // Total number of suffix-intervals
  size_t k = n - om_.perm_[1][i]; // Number of intervals to close
  // push the begins
  while(k > 0)
    {
      ++i;
      n += om_.perm_[0][i];
      k += om_.perm_[0][i];
      k -= om_.perm_[1][i];
    }
  for(size_t i__ = 0; i__ < n-1; ++i__)
    starts_.next();
  startp_ = i+1;
  while(startp_ < om_.perm_len_ && om_.perm_[0][startp_] == 0) ++startp_;
  return new segm_start_t(perm_beg, intvs_start);
}

Segment* SegmConverter::operator()(const segm_start_t* s) const
{
  OvMgr::perm_type::size_type perm_beg = std::get<0>(*s);
  OvMgr::bitvect_type::size_type bpos = std::get<1>(*s);
  delete s;

  OvMgr::bitvect_type::ones_enumerator_t starts_ = om_.bvs_[0].ones_enumerator(bpos);
  OvMgr::bitvect_type::ones_enumerator_t ends_   = om_.bvs_[1].ones_enumerator(bpos);

  std::deque<SmallBiIntv> int_info;
  std::stack<std::deque<SmallBiIntv>::size_type> end_pos;
  std::deque<BiIntv> suffixes;
  OvMgr::perm_type::size_type i = perm_beg;
  size_t n = om_.perm_[0][i]; // Total number of suffix-intervals
  size_t k = n - om_.perm_[1][i]; // Number of intervals to close
  for(OvMgr::perm_el_type i__ = 0; i__ < om_.perm_[0][i]; ++i__)
    {
      end_pos.push(i__);
    }
  int_info.resize(end_pos.size(), SmallBiIntv(0, 0));
  while(k > 0)
    {
      for(OvMgr::perm_el_type i__ = 0; i__ < om_.perm_[1][i]; ++i__)
        {
          assert(perm_beg + std::numeric_limits<BiIntvStorage<SmallBiIntv>::type>::max() >= i);
          ends(int_info[end_pos.top()]) = static_cast<BiIntvStorage<SmallBiIntv>::type>(i - perm_beg);
          end_pos.pop();
        }
      ++i;
      n += om_.perm_[0][i];
      k += om_.perm_[0][i];
      k -= om_.perm_[1][i];
      for(OvMgr::perm_el_type i__ = 0; i__ < om_.perm_[0][i]; ++i__)
        {
          end_pos.push(int_info.size());
          int_info.push_back(SmallBiIntv(i - perm_beg, 0));
        }
    }
  for(OvMgr::perm_el_type i__ = 0; i__ < om_.perm_[1][i]; ++i__)
    {
      assert(perm_beg + std::numeric_limits<BiIntvStorage<SmallBiIntv>::type>::max() >= i);
      ends(int_info[end_pos.top()]) = static_cast<BiIntvStorage<SmallBiIntv>::type>(i - perm_beg);
      end_pos.pop();
    }
  for(size_t i__ = 0; i__ < n; ++i__)
    {
      OvMgr::bitvect_type::size_type intvs_start= starts_.next();
      OvMgr::bitvect_type::size_type intvs_end= ends_.next();
      suffixes.push_back(BiIntv(intvs_start, intvs_end));
    }
  return new Segment(perm_beg, i - perm_beg +1, suffixes, int_info);
}
