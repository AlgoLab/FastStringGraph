//-----------------------------------------------
// Copyright 2015 Yuri Pirola, Marco Previtali
// Released under the GPL
//-----------------------------------------------
#ifndef HITSHPP
#define HITSHPP

#include <vector>

class HitsConverter;

class Hits
{
public:
  const uint64_t        pref_;
  std::vector<uint64_t> suffs_;
  const uint8_t         ext_len_;

  friend class HitsConverter;

  Hits(const uint64_t p, const uint8_t ext_len) noexcept
    : pref_(p), suffs_(), ext_len_(ext_len)
  {}

  Hits(Hits&& other) noexcept
    : pref_(other.pref_), ext_len_(other.ext_len_)
  {
    using std::swap;
    swap(suffs_, other.suffs_);
  }

  inline void add_suff(const uint64_t s);
};

void Hits::add_suff(const uint64_t s)
{
  suffs_.push_back(s);
}


#include <deque>
using HitsVector = std::deque<Hits>;


struct asqg_edge_t {
  std::string pref_id;
  std::string suff_id;
  size_t pref_ov_start;
  size_t pref_ov_end;
  size_t pref_len;
  size_t suff_ov_start;
  size_t suff_ov_end;
  size_t suff_len;
  bool is_reversed;

  asqg_edge_t(std::string&& pref_id_, std::string&& suff_id_) noexcept
  : pref_id(std::forward<std::string>(pref_id_)),
    suff_id(std::forward<std::string>(suff_id_))
  {}
};


template <typename T>
using RawEdgeVectorContainer = std::deque<T>;
using RawEdgeVector = std::tuple<RawEdgeVectorContainer<uint32_t>,
                                 RawEdgeVectorContainer<uint32_t>,
                                 RawEdgeVectorContainer<uint8_t>,
                                 std::vector<bool> >;

inline void
add_rawedge(RawEdgeVector& v,
            const uint32_t id1,
            const uint32_t id2,
            const uint8_t len,
            const bool rev1, const bool rev2) {
  using std::get;
  get<0>(v).push_back(id1);
  get<1>(v).push_back(id2);
  get<2>(v).push_back(len);
  get<3>(v).push_back(rev1);
  get<3>(v).push_back(rev2);
}

#endif
