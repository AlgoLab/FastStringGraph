//-----------------------------------------------
// Copyright 2015 Yuri Pirola, Marco Previtali
// Released under the GPL
//-----------------------------------------------
#include "hits_util.hpp"

//HitsProducer
const HitsVector* HitsProducer::operator()(Segment* s) const
{
  HitsVector* hv = new HitsVector;
  extend_segment(s, hv, 0);
  return hv;
}

void HitsProducer::extend_segment(Segment* s, HitsVector* hv, uint8_t ext_len) const
{
  std::deque<BiIntv>      nIntv[DNA_ALPHABET::size]; // Extended intervals
  std::deque<SmallBiIntv> nInfo[DNA_ALPHABET::size]; // Prefix information
  // Extend all the suffixes in the segment
  AlphaCount64 bOccs, eOccs;
  for(size_t i = s->size(); i>0; --i)
    {
      const BiIntv& _interval = s->get_interval(i-1);
      const SmallBiIntv& _infos = s->get_info(i-1);
      pBWT_->getIntervalFullOcc(begins(_interval)-1, ends(_interval),
                                bOccs, eOccs);
      const auto bOccs_$= bOccs.get('$');
      const auto eOccs_$= eOccs.get('$');
      if(eOccs_$ != bOccs_$)
        {
          for(uint64_t index = bOccs_$; index < eOccs_$; ++index)
            {
              hv->emplace_back(index, ext_len);
              s->addTrue(hv->back(), _infos);
            }
          s->reduce(_infos);
        }
      else if(!s->reduced(_infos))
        for(const char c : ALPHABET_)
          {
            if(bOccs.get(c) != eOccs.get(c))
              {
                const uint8_t c_rank= DNA_ALPHABET::getBaseRank(c);
                const size_t _bp = pBWT_->getPC(c);
                nIntv[c_rank].push_front(BiIntv(_bp + bOccs.get(c), _bp + eOccs.get(c) -1));
                nInfo[c_rank].push_front(_infos);
              }
          }
      // else the arc is reduced, discard it
    }

  // Recursive call on all the new segments
  for(const char c : ALPHABET_) {
    const uint8_t c_rank= DNA_ALPHABET::getBaseRank(c);
    if(!nIntv[c_rank].empty())
      {
        Segment* ns = new Segment(*s, nIntv[c_rank], nInfo[c_rank]);
        extend_segment(ns, hv, ext_len +1);
      }
  }
  delete s;
}


size_t convert_hits(const HitsVector* hv,
                  const SuffixArray& SA_,
                  const size_t numreads,
                  RawEdgeVector& raw_edge_vect)
{
  const size_t numreads_= numreads/2;
  size_t new_edges = 0;
  for(const Hits& h: *hv) {
    uint64_t prefindex = SA_.get(h.pref_).getID();
    const bool prefreverse = (prefindex >= numreads_);
    if (prefreverse) prefindex -= numreads_;
    for(const uint64_t suff: h.suffs_) {
      uint64_t suffindex = SA_.get(suff).getID();
      const bool suffreverse = (suffindex >= numreads_);
      if (suffreverse) suffindex -= numreads_;
      if(suffindex > prefindex) {
        add_rawedge(raw_edge_vect,
                    prefindex, suffindex,
                    h.ext_len_,
                    prefreverse, suffreverse);
        ++new_edges;
      }
    }
  }
  return new_edges;
}


static inline void
format_edge_to_buffer(const asqg_edge_t& edge,
                      cbuffer_t::iterator& buf) {
  *(buf++)= 'E';  *(buf++)= 'D';  *(buf++)= SQG::FIELD_SEP;
  format_to_buffer(edge.pref_id,       buf); *(buf++)= ' ';
  format_to_buffer(edge.suff_id,       buf); *(buf++)= ' ';
  format_to_buffer(edge.pref_ov_start, buf); *(buf++)= ' ';
  format_to_buffer(edge.pref_ov_end,   buf); *(buf++)= ' ';
  format_to_buffer(edge.pref_len,      buf); *(buf++)= ' ';
  format_to_buffer(edge.suff_ov_start, buf); *(buf++)= ' ';
  format_to_buffer(edge.suff_ov_end,   buf); *(buf++)= ' ';
  format_to_buffer(edge.suff_len,      buf); *(buf++)= ' ';
  *(buf++)= edge.is_reversed ? '1' : '0';
  *(buf++)= ' '; *(buf++)= '0';
  *(buf++)= '\n';
}

void convert_and_compress_rawedges(cbuffers_t& out_bufs,
                                   cbuffer_t& block_buf,
                                   const RawEdgeVector& rev,
                                   const ReadInfoTable& rInfTab) {
  using std::get;

  const size_t prev_size= block_buf.size();
  block_buf.resize(GZIP_BLOCK_SIZE);
  cbuffer_t::iterator it= block_buf.begin() + prev_size;

  auto prefindex_it = get<0>(rev).cbegin();
  auto suffindex_it = get<1>(rev).cbegin();
  auto len_it = get<2>(rev).cbegin();
  auto rev_it = get<3>(rev).cbegin();
  while (prefindex_it != get<0>(rev).cend()) {

    if (block_buf.end() - it < 512) {
      out_bufs.emplace_back(GZIP_BLOCK_SIZE);
      block_buf.resize(it - block_buf.begin());
      compress_buffer(block_buf, out_bufs.back());
      block_buf.resize(GZIP_BLOCK_SIZE);
      it = block_buf.begin();
    }

    // Prefix is a suffix (and viceversa)
    const uint32_t prefindex = *suffindex_it; ++suffindex_it;
    const uint32_t suffindex = *prefindex_it; ++prefindex_it;
    const uint8_t len = *len_it; ++len_it;
    const bool prefreverse = *rev_it; ++rev_it;
    const bool suffreverse = *rev_it; ++rev_it;
    asqg_edge_t edge(rInfTab.getReadID(prefindex), rInfTab.getReadID(suffindex));
    edge.pref_len = rInfTab.getReadLength(prefindex);
    edge.suff_len = rInfTab.getReadLength(suffindex);
    const uint8_t ovlen = edge.suff_len - len;
    if(!prefreverse) {
      edge.pref_ov_start = 0;
      edge.pref_ov_end = ovlen - 1;
    } else {
      edge.pref_ov_start = edge.pref_len - ovlen;
      edge.pref_ov_end = edge.pref_len - 1;
    }
    if(!suffreverse) {
      edge.suff_ov_start = edge.suff_len - ovlen;
      edge.suff_ov_end = edge.suff_len - 1;
    } else {
      edge.suff_ov_start = 0;
      edge.suff_ov_end = ovlen - 1;
    }
    edge.is_reversed = prefreverse != suffreverse;

    format_edge_to_buffer(edge, it);
  }
  block_buf.resize(it - block_buf.begin());
}

