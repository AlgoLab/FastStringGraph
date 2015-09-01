//-----------------------------------------------
// Copyright 2015 Yuri Pirola, Marco Previtali
// Released under the GPL
//-----------------------------------------------
#include "compressor.hpp"

#include <zlib.h>
#include <cassert>

void
compress_buffer(cbuffer_t& in_buf, cbuffer_t& out_buf) {
  out_buf.resize(GZIP_BLOCK_SIZE);
  z_stream zs;
  zs.zalloc = Z_NULL;
  zs.zfree = Z_NULL;
  zs.next_in = &(in_buf[0]);
  zs.avail_in = static_cast<uInt>(in_buf.size());
  zs.next_out = static_cast<uint8_t*>(&out_buf[0]);
  zs.avail_out = static_cast<uInt>(out_buf.size());
#ifndef NDEBUG
  int status =
#endif
    deflateInit2(&zs, Z_BEST_COMPRESSION, Z_DEFLATED, 15+16, 8, Z_DEFAULT_STRATEGY); // +16 to enable gzip header/footer
#ifndef NDEBUG
  assert(status == Z_OK);
  status =
#endif
    deflate(&zs, Z_FINISH);
#ifndef NDEBUG
  assert(status == Z_STREAM_END);
  status =
#endif
    deflateEnd(&zs);
#ifndef NDEBUG
  assert(status == Z_OK);
#endif
  out_buf.resize(zs.total_out);
}

