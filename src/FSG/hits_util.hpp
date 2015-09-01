//-----------------------------------------------
// Copyright 2015 Yuri Pirola, Marco Previtali
// Released under the GPL
//-----------------------------------------------
#ifndef HITSUTIL_HPP
#define HITSUTIL_HPP

#include <iostream>

#include "BWT.h"

#include "segment.hpp"
#include "bwt_interval.hpp"
#include "FSGCommon.h"
#include "ReadInfoTable.h"
#include "ASQG.h"
#include "compressor.hpp"

class HitsProducer
{
private:
  const BWT* pBWT_;
  void extend_segment(Segment* s, HitsVector* hv, uint8_t ext_len) const;

public:
  HitsProducer(const BWT* p) : pBWT_(p) {}
  const HitsVector* operator()(Segment* s) const;
};

size_t convert_hits(const HitsVector* hv,
                    const SuffixArray& SA_,
                    const size_t numreads,
                    RawEdgeVector& raw_edge_vect);

void convert_and_compress_rawedges(cbuffers_t& out_bufs,
                                   cbuffer_t& block_buf,
                                   const RawEdgeVector& rev,
                                   const ReadInfoTable& rInfTab);

#endif
