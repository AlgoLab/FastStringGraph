//-----------------------------------------------
// Copyright 2015 Yuri Pirola, Marco Previtali
// Released under the GPL
//-----------------------------------------------
#include "fsg.h"

#include <deque>
#include <memory>
#include <thread>
#include <getopt.h>
#include <cstdio>
#include <iostream>
#include <iomanip>

#include "suff_mgr.hpp"
#include "ov_mgr.hpp"
#include "bwt_interval.hpp"
#include "hits.hpp"
#include "hits_util.hpp"
#include "segment.hpp"
#include "FSGCommon.h"
#include "gzstream.h"
#include "compressor.hpp"
#include "parallel_util.hpp"
#include "delta_interval.hpp"

#include "config.h"
#include "BWT.h"
#include "SGACommon.h"
#include "ReadInfoTable.h"
#include "SeqReader.h"
#include "Util.h"
#include "ASQG.h"

#include <tbb/tbb.h>
using tbb::parallel_pipeline;
using tbb::make_filter;
using tbb::filter;

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;


using TriDeltaInterval = DeltaInterval<TriIntvStorage<TrIntv>::type, uint8_t>;


//
// Getopt
//
#define SUBPROGRAM "fsg"
static const char *FSG_VERSION_MESSAGE =
  SUBPROGRAM " Version " PACKAGE_VERSION "\n";

static const char *FSG_USAGE_MESSAGE =
  "Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
  "Compute the string graph of the sequences in READSFILE\n"
  "\n"
  "      --help                            display this help and exit\n"
  "\n"
  "      -m, --min-overlap=LEN             minimum overlap required between two reads (default:45)\n"
  "      -p, --prefix=PREFIX               use PREFIX for the names of the index files (default: prefix of the input file)\n"
  "      -d, --sample-rate=N               sample the symbol counts every N symbols in the FM-index. Higher values use significantly\n"
  "                                        less memory at the cost of higher runtime. This value must be a power of 2 (default: 128)\n"
  "\n";

static const char *PROGRAM_IDENT = PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
  static unsigned int numThreads=1;
  static std::string readsFile;
  static std::string outFile;
  static std::string prefix;
  static unsigned int queueSize=512;

  static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;
  static unsigned int minOverlap = DEFAULT_MIN_OVERLAP;
}

static const char *shortopts = "m:p:d:t:q:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
  { "threads",     required_argument, NULL, 't' },
  { "min-overlap", required_argument, NULL, 'm' },
  { "prefix",      required_argument, NULL, 'p' },
  { "sample-rate", required_argument, NULL, 'd' },
  { "queue-size",  required_argument, NULL, 'q' },
  { "help",        no_argument,       NULL, OPT_HELP },
  { "version",     no_argument,       NULL, OPT_VERSION },
  { NULL, 0, NULL, 0}
};


static std::tuple<size_t, size_t>
extend_partial_overlaps(TriDeltaInterval::enumerator_t& src_en,
                        OvMgr& overlaps, SuffMgr& set_sink,
                        const BWT& bwt, const bool are_overlaps) {
  size_t suffnum= 0;
  size_t no_of_overlaps= 0;
  AlphaCount64 bOccs, esOccs, epOccs;
  while (src_en.has_next())
    {
      // Remember that end positions are part of the interval
      TrIntv _interval = src_en.next();
      // AlphaCount64 bOccs = bwt.getFullOcc(begins(_interval)-1);
      // AlphaCount64 esOccs= bwt.getFullOcc(ends(_interval));
      // AlphaCount64 epOccs= bwt.getFullOcc(endp(_interval));
      bwt.getTriIntervalFullOcc(begins(_interval)-1, ends(_interval), endp(_interval),
                                bOccs, esOccs, epOccs);

      // Is it an overlap?
      if(are_overlaps && esOccs.get('$') != epOccs.get('$')) {
        overlaps.add(BiIntv(begins(_interval), ends(_interval)),
                     BiIntv(esOccs.get('$'), epOccs.get('$') -1));
        ++no_of_overlaps;
      }

      // Can we extend the pattern and the suffix with the same nucleotide?
      for(const char c : ALPHABET_)
        if(bOccs.get(c) != esOccs.get(c) && esOccs.get(c) != epOccs.get(c))
          {
            size_t _bp = bwt.getPC(c);
            set_sink.add(TrIntv( _bp + bOccs.get(c),
                                 _bp + esOccs.get(c) -1,
                                 _bp + epOccs.get(c) -1));
            ++suffnum;
          }
    }
  return std::make_tuple(suffnum, no_of_overlaps);
}

static void
step1_compute_overlaps(const BWT& bwt, OvMgr& overlaps,
                       const unsigned int min_overlap,
                       const system_clock::time_point start= timer::now()) {
  using std::swap;

  const size_t numstr = bwt.getNumStrings();
  const size_t bwtlen = bwt.getBWLen();

  size_t suffnum =0; // Number of suffix to extend
  size_t new_suffnum =0; // Number of suffix to extend

  auto stop = timer::now();
  std::cout << "[" << std::setw(12) << duration_cast<milliseconds>(stop-start).count()
            << "] Computing initial suffixes." << std::endl;
  // Generate inital suffix-intervals
  tbb::concurrent_queue<TriDeltaInterval*> di_int_q;
  {
    TriDeltaInterval* _di = new TriDeltaInterval;
    AlphaCount64 occEndSuffixes = bwt.getFullOcc(numstr-1);
    AlphaCount64 occEndBWT = bwt.getFullOcc(bwtlen-1);
    for(const char c : ALPHABET_)
      if(occEndSuffixes.get(c) != 0 && occEndSuffixes.get(c) != occEndBWT.get(c))
        {
          size_t _bp = bwt.getPC(c);
          _di->add(TrIntv(_bp,
                          _bp + occEndSuffixes.get(c) -1,
                          _bp + occEndBWT.get(c) -1));
          ++suffnum;
        }
    di_int_q.push(_di);
  }



  SuffMgr set_sink(bwtlen); // Temporary SetMgr for suffixes of length 1

  size_t no_of_overlaps= 0;
  for(size_t _iter_n =1; suffnum!=0; ++_iter_n)
    {
      stop = timer::now();
      std::cout << "[" << std::setw(12) << duration_cast<milliseconds>(stop-start).count()
                << "] Iteration #" << _iter_n
                << " number of suffix-intervals to extend: " << suffnum << std::endl;

      new_suffnum= 0;

      std::atomic<size_t> conc_new_suffnum; conc_new_suffnum= 0;
      std::atomic<size_t> conc_new_overlaps; conc_new_overlaps= 0;
      parallel_pipeline(
         opt::numThreads,
         make_filter<void, TriDeltaInterval*>(
             filter::serial_in_order,
             [&](tbb::flow_control& fc) -> TriDeltaInterval*
             {
               TriDeltaInterval* _di= nullptr;
               di_int_q.try_pop(_di);
               if (_di == nullptr)  fc.stop();
               return _di;
             }
         )  &
         make_filter<TriDeltaInterval*, void>(
             filter::parallel,
             [&](TriDeltaInterval* di) -> void {
               auto di_en = di->enumerate();
               auto ext_part_res =
                 extend_partial_overlaps(di_en, overlaps, set_sink, bwt,
                                         _iter_n >= min_overlap);
               delete di;
               conc_new_suffnum += std::get<0>(ext_part_res);
               conc_new_overlaps += std::get<1>(ext_part_res);
             }
         )
      );
      new_suffnum = conc_new_suffnum;
      no_of_overlaps += conc_new_overlaps;

      // Initialize source for the next iteration
      size_t start_block= 0;
      const size_t block_len= ((bwtlen/opt::numThreads)+1);
      parallel_pipeline(
           opt::numThreads,
           make_filter<void, SuffMgr::enumerator*>(
                filter::serial_in_order,
                [&](tbb::flow_control& fc) -> SuffMgr::enumerator* {
                  if (start_block == bwtlen) {
                    fc.stop();
                    return NULL;
                  } else {
                    const size_t curr_start= start_block;
                    start_block= std::min(bwtlen, start_block+block_len);
                    return new SuffMgr::enumerator(set_sink.enumerate_in_range(curr_start, start_block));
                  }
                }
           )  &
           make_filter<SuffMgr::enumerator*, void>(
                filter::parallel,
                [&](SuffMgr::enumerator* ptr_src_en) -> void
                {
                  if (ptr_src_en->has_next()) {
                    TriDeltaInterval* _di= new TriDeltaInterval;
                    while(ptr_src_en->has_next())
                      _di->add(ptr_src_en->next());
                    di_int_q.push(_di);
                  }
                  delete ptr_src_en;
                }
           )
      );

      set_sink.reinitialize();
      suffnum = new_suffnum;
    }

  stop = timer::now();
  std::cout << "[" << std::setw(12) << duration_cast<milliseconds>(stop-start).count()
            << "] Overlaps computed: " << no_of_overlaps << std::endl;
}

static inline void
format_seqrecord_to_buffer(const std::string& id,
                           const std::string& dnaseq,
                           cbuffer_t::iterator& buf) {
  *(buf++)= 'V';  *(buf++)= 'T';  *(buf++)= SQG::FIELD_SEP;
  format_to_buffer(id,     buf); *(buf++)= SQG::FIELD_SEP;
  format_to_buffer(dnaseq, buf); *(buf++)= '\n';
}

static inline void
format_seqrecord_to_buffer_or_overflow(const SeqRecord& sr,
                                       const cbuffer_t& buf, cbuffer_t::iterator& it,
                                       cbuffer_t& ov_buf) {
  const size_t id_len= sr.id.length();
  const std::string dnaseq(sr.seq.toString());
  const size_t seq_len= dnaseq.length();
  if ( it + id_len + seq_len + 5 <= buf.end()) {
    format_seqrecord_to_buffer(sr.id, dnaseq, it);
  } else {  // Overflow
    ov_buf.resize(id_len + seq_len + 5);
    cbuffer_t::iterator ov_it= ov_buf.begin();
    format_seqrecord_to_buffer(sr.id, dnaseq, ov_it);
  }
}

static void
prepare_asqg_vertices(std::ostream& asqg_out, const std::string& readsFile,
                      const size_t& numstr, ReadInfoTable& rInfTab)
{
  SeqReader reader(readsFile);
  SeqRecord sr;
  size_t size = 0;

  const size_t N_TOKENS= opt::numThreads+2;

  using block_pool_t= obj_pool_t<cbuffer_t>;
  const auto gzip_block_ctor= [] () -> cbuffer_t* { return new cbuffer_t(GZIP_BLOCK_SIZE); };

  block_pool_t block_pool(gzip_block_ctor);
  cbuffer_t overflow_buffer(0); // The overflow buffer has initial size 0 because no overflow has occurred yet

  parallel_pipeline(
     N_TOKENS,
     make_filter<void, cbuffer_t*>(filter::serial_in_order,
        [&](tbb::flow_control& fc) -> cbuffer_t*
        {
          cbuffer_t* out_buf= block_pool.get();
          out_buf->resize(GZIP_BLOCK_SIZE);
          cbuffer_t::iterator it= out_buf->begin();
          if (!overflow_buffer.empty()) {
            std::copy(overflow_buffer.cbegin(), overflow_buffer.cend(), it);
            it += overflow_buffer.size();
            overflow_buffer.clear();
          }
          while ((out_buf->end() - it >= 256) && // leave 256 bytes to minimize chance of overflow
                 (size++ < numstr) &&
                 reader.get(sr)) {
            rInfTab.add(sr);
            format_seqrecord_to_buffer_or_overflow(sr, *out_buf, it, overflow_buffer);
            if (!overflow_buffer.empty())  break;
          }
          if (it == out_buf->begin()) {
            block_pool.give_back(out_buf);
            fc.stop();
            return nullptr;
          } else {
            out_buf->resize(it - out_buf->begin());
            return out_buf;
          }
        })  &
     make_filter<cbuffer_t*, cbuffer_t*>(filter::parallel,
        [&](cbuffer_t* in_buf) -> cbuffer_t*
        {
          cbuffer_t* out_buf= block_pool.get();
          out_buf->resize(GZIP_BLOCK_SIZE);
          compress_buffer(*in_buf, *out_buf);
          block_pool.give_back(in_buf);
          return out_buf;
        }) &
     make_filter<cbuffer_t*, void>(filter::serial_in_order,
        [&](cbuffer_t* in_buf) -> void
        {
          asqg_out.write(reinterpret_cast<char*>(&((*in_buf)[0])), in_buf->size());
          block_pool.give_back(in_buf);
        })
     );
  assert(overflow_buffer.empty());
}

struct rawedgevect_pool_t {
  using pool_t_= tbb::concurrent_queue<RawEdgeVector*>;

  pool_t_ incomplete_pool_;
  pool_t_ full_pool_;

  rawedgevect_pool_t() = default;
  rawedgevect_pool_t(const rawedgevect_pool_t&) = delete;

  ~rawedgevect_pool_t() {
    RawEdgeVector* p_obj= nullptr;
    while ((p_obj = this->next()) != nullptr)
      delete p_obj;
  }

  RawEdgeVector* get() {
    RawEdgeVector* p_obj= nullptr;
    if (!incomplete_pool_.try_pop(p_obj))
      p_obj= new RawEdgeVector;
    return p_obj;
  }

  void give_back(RawEdgeVector* const p_obj) {
    using std::get;
    if (get<0>(*p_obj).size() > 10000)
      full_pool_.push(p_obj);
    else
      incomplete_pool_.push(p_obj);
  }

  RawEdgeVector* next() {
    RawEdgeVector* p_obj= nullptr;
    if (!full_pool_.try_pop(p_obj))
      if (!incomplete_pool_.try_pop(p_obj))
        p_obj= nullptr;
    return p_obj;
  }

};

static size_t
step2_compute_rawedges(const BWT& bwt, const size_t numstr,
                       const SuffixArray& SA, const OvMgr& overlaps,
                       rawedgevect_pool_t& rawedgevect_pool,
                       const system_clock::time_point start) {
  using std::get;
  std::cout << "[" << std::setw(12) << duration_cast<milliseconds>(timer::now()-start).count()
            << "] Starting the reduction pipeline with " << opt::numThreads
            << " threads." << std::endl;

  size_t n_segments= 0;
  const size_t n_segments_for_message = opt::numThreads * 500000;
  OvMgrVisitor ovmgrvis(overlaps);

  std::atomic<size_t> n_edges; n_edges= 0;

  parallel_pipeline(
     opt::numThreads,
     make_filter<void, segm_start_t*>(filter::serial_in_order,
        [&](tbb::flow_control& fc) -> segm_start_t*
        {
          segm_start_t* s= ovmgrvis();
          if (++n_segments % n_segments_for_message == 0) {
            std::cout << "[" << std::setw(12) << duration_cast<milliseconds>(timer::now()-start).count()
                      << "] Produced " << n_segments << " segments" << std::endl;
          }
          if (s == NULL) {
            std::cout << "[" << std::setw(12) << duration_cast<milliseconds>(timer::now()-start).count()
                      << "] Finished to produce segments (" << n_segments << " produced)" << std::endl;
            fc.stop();
          }
          return s;
        }
     ) &
     make_filter<segm_start_t*, Segment*>(filter::parallel,
                                          SegmConverter(overlaps))
     &
     make_filter<Segment*, const HitsVector*>(filter::parallel,
                                              HitsProducer(&bwt))
     &
     make_filter<const HitsVector*, void>(filter::parallel,
        [&](const HitsVector* hv)
        {
          RawEdgeVector* v= rawedgevect_pool.get();
          const size_t new_edges=
            convert_hits(hv, SA, numstr, *v);
          n_edges += new_edges;
          rawedgevect_pool.give_back(v);
          delete hv;
        }
     ));
  std::cout << "[" << std::setw(12) << duration_cast<milliseconds>(timer::now()-start).count()
            << "] Finished to produce edges" << std::endl;
  std::cout << "[" << std::setw(12) << duration_cast<milliseconds>(timer::now()-start).count()
            << "] Edges computed: " << n_edges << std::endl;
  return n_edges;
}

static void
write_asqg_edges(const ReadInfoTable& rInfTab,
                 rawedgevect_pool_t& rawedgevect_pool,
                 std::ofstream& asqg_out) {
  using block_pool_t= obj_pool_t<cbuffer_t>;
  using blocks_pool_t= obj_pool_t<cbuffers_t>;

  const auto incomplete_block_flush = [&asqg_out] (cbuffer_t* block) -> void {
    cbuffer_t out_buf(GZIP_BLOCK_SIZE);
    compress_buffer(*block, out_buf);
    asqg_out.write(reinterpret_cast<const char*>(&out_buf[0]), out_buf.size());
    delete block;
  };

  const auto empty_block_ctor= [] () -> cbuffer_t* { return new cbuffer_t(0); };

  block_pool_t incomplete_blocks_pool(empty_block_ctor, incomplete_block_flush);
  blocks_pool_t blocks_pool;

  parallel_pipeline(
     opt::numThreads,
     make_filter<void, RawEdgeVector*>(filter::serial_in_order,
        [&](tbb::flow_control& fc) -> RawEdgeVector*
        {
          RawEdgeVector* v= rawedgevect_pool.next();
          if (v == nullptr) fc.stop();
          return v;
        }
     ) &
     make_filter<RawEdgeVector*, cbuffers_t*>(filter::parallel,
        [&](RawEdgeVector* rev) -> cbuffers_t*
        {
          if (rev == nullptr) return nullptr;
          cbuffers_t* out_bufs= blocks_pool.get();
          cbuffer_t* block_buf= incomplete_blocks_pool.get();

          convert_and_compress_rawedges(*out_bufs, *block_buf,
                                        *rev, rInfTab);

          incomplete_blocks_pool.give_back(block_buf);
          delete rev;
          return out_bufs;
        }
     ) &
     make_filter<cbuffers_t*, void>(filter::serial_out_of_order,
        [&](cbuffers_t* in_bufs) -> void
        {
          if (in_bufs == nullptr) return;
          for (const cbuffer_t& buf: *in_bufs) {
            asqg_out.write(reinterpret_cast<const char*>(&buf[0]), buf.size());
          }
          in_bufs->clear();
          blocks_pool.give_back(in_bufs);
        }
      ));
}

//
// Main
//
int fsgMain(int argc, char** argv)
{

  const system_clock::time_point start= timer::now();
  system_clock::time_point stop= start;

  parseFsgOptions(argc, argv);
  std::string indexPrefix;
  if(!opt::prefix.empty())
    indexPrefix = opt::prefix;
  else
    indexPrefix = stripExtension(opt::readsFile);

  tbb::task_scheduler_init init(opt::numThreads);

  size_t numstr= 0;
  size_t bwtlen= 0;

  { // RawEdges block
    rawedgevect_pool_t rawedgevect_pool;

  { // BWT + OvMgr block
    // Initialization: load the BWT
    std::cout << "[" << std::setw(12) << duration_cast<milliseconds>(timer::now()-start).count()
              << "] Reading the BWT..." << std::endl;
    const BWT bwt(std::string(indexPrefix + BWT_EXT), opt::sampleRate);
    numstr = bwt.getNumStrings();
    bwtlen = bwt.getBWLen();
    bwt.printInfo();


    // First step: compute the overlap graph
    OvMgr overlaps(bwtlen, numstr);
    step1_compute_overlaps(bwt, overlaps,
                           opt::minOverlap,
                           start);

  { // SA block

    SuffixArray SA(indexPrefix + ".sai");

    // 2nd step: compute the irreducible edges
    step2_compute_rawedges(bwt, numstr, SA, overlaps,
                           rawedgevect_pool,
                           start);

  } // END -- SA block
  } // END -- BWT + OvMgr block

  { // ReadInfoTable
    // Initialize asgq and readinfotable
    std::cout << "[" << std::setw(12) << duration_cast<milliseconds>(timer::now()-start).count()
              << "] Reading reads and writing the vertices to the ASQG file" << std::endl;
    ReadInfoTable rInfTab(numstr/2);
    std::ofstream asqg_out(std::string(indexPrefix + ASQG_EXT + GZIP_EXT).c_str());
    prepare_asqg_vertices(asqg_out, opt::readsFile, numstr/2, rInfTab);


    std::cout << "[" << std::setw(12) << duration_cast<milliseconds>(timer::now()-start).count()
              << "] Writing the edges to the ASQG file" << std::endl;
    write_asqg_edges(rInfTab, rawedgevect_pool, asqg_out);
  } // END -- ReadInfoTable

  } // END -- RawEdges block

  stop = timer::now();
  std::cout << "[" << std::setw(12) << duration_cast<milliseconds>(stop-start).count()
            << "] Execution completed" << std::endl;

  return 0;
}

//
// Parse options
//
void parseFsgOptions(int argc, char** argv)
{
  bool die = false;
  for(char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
      std::istringstream arg(optarg != NULL ? optarg : "");
      switch(c)
        {
        case 't': arg >> opt::numThreads; break;
        case 'm': arg >> opt::minOverlap; break;
        case 'p': arg >> opt::prefix; break;
        case 'd': arg >> opt::sampleRate; break;
        case 'q': arg >> opt::queueSize; break;
        case '?': die = true; break;
        case OPT_HELP:
          std::cout << FSG_USAGE_MESSAGE;
          exit(EXIT_SUCCESS);
        case OPT_VERSION:
          std::cout << FSG_VERSION_MESSAGE;
          exit(EXIT_SUCCESS);
        }
    }

  if (argc - optind < 1)
    {
      std::cerr << SUBPROGRAM ": missing arguments\n";
      die = true;
    }
  else if (argc - optind > 1)
    {
      std::cerr << SUBPROGRAM ": too many arguments\n";
      die = true;
    }

  if (die)
    {
      std::cout << "\n" << FSG_USAGE_MESSAGE;
      exit(EXIT_FAILURE);
    }

  // Parse the input filenames
  opt::readsFile = argv[optind++];
}
