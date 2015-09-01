//-----------------------------------------------
// Copyright 2015 Yuri Pirola, Marco Previtali
// Released under the GPL
//-----------------------------------------------
#ifndef OVMGR_HPP
#define OVMGR_HPP

#include <memory>
#include <atomic>

#include "segment.hpp"
#include "simple_bitvector.hpp"

class OvMgrVisitor;
class SegmConverter;

class OvMgr
{
public:
  typedef concurrent_bitvector_t bitvect_type;
  typedef bitvect_type::size_type size_type;
  typedef uint8_t perm_el_type;
  typedef std::vector<std::atomic<perm_el_type> > perm_type;

  friend class OvMgrVisitor;
  friend class SegmConverter;

public:
  explicit OvMgr(const size_type bvs_len, const size_t perm_len);

  inline void add(const BiIntv& suff, const BiIntv& pref)
  {
    bvs_[0].set(begins(suff));
    bvs_[1].set(ends(suff));
    ++perm_[0][begins(pref)];
    ++perm_[1][ends(pref)];
  }

private:
  // const bitvect_type::size_type   bvs_len_;
  std::array<bitvect_type, 2>     bvs_; // suffixes
  const perm_type::size_type      perm_len_;
  std::array<perm_type, 2>        perm_; // overlap intervals on the permutation of the reads
};

typedef std::tuple<OvMgr::perm_type::size_type, OvMgr::bitvect_type::size_type> segm_start_t;

class OvMgrVisitor
{
private:
  const OvMgr& om_;
  OvMgr::perm_type::size_type startp_;
  OvMgr::bitvect_type::ones_enumerator_t starts_;

public:
  OvMgrVisitor(const OvMgr& om)
    : om_(om), startp_(0),
      starts_(om_.bvs_[0].ones_enumerator())
  {}

  OvMgrVisitor(const OvMgrVisitor& other) = default;

  segm_start_t* operator()();
};

class SegmConverter
{
private:
  const OvMgr& om_;

public:
  SegmConverter(const OvMgr& om)
    : om_(om)
  {}

  SegmConverter(const SegmConverter& other) = default;

  Segment* operator()(const segm_start_t*) const;
};

#endif
