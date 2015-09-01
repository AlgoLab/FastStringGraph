//-----------------------------------------------
// Copyright 2015 Yuri Pirola, Marco Previtali
// Released under the GPL
//-----------------------------------------------
#ifndef COMPRESSOR_HPP
#define COMPRESSOR_HPP

#include <cstdint>
#include <string>
#include <vector>
#include <deque>

#include <boost/lexical_cast.hpp>

#define GZIP_BLOCK_SIZE     0x10000  // 64K  good speed/compression tradeoff

using cbuffer_t = std::vector<uint8_t>;
using cbuffers_t = std::deque<cbuffer_t>;


void compress_buffer(cbuffer_t& in_buf, cbuffer_t& out_buf);


template <typename T>
inline void
format_to_buffer(T&& val, cbuffer_t::iterator& buf) {
  format_to_buffer<std::string>(boost::lexical_cast<std::string>(val), buf);
}

template <>
inline void
format_to_buffer(std::string&& str, cbuffer_t::iterator& buf) {
  const size_t str_len= str.length();
  std::copy_n(str.cbegin(), str_len, buf);
  buf += str_len;
}

#endif
