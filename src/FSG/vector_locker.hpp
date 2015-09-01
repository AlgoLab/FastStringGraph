//-----------------------------------------------
// Copyright 2015 Yuri Pirola, Marco Previtali
// Released under the GPL
//-----------------------------------------------
#ifndef VECTOR_LOCKER_HPP
#define VECTOR_LOCKER_HPP

template <size_t N>
struct log_2_t
{
  enum { value = 1 + log_2_t< (N>>1) >::value };
};

template <>
struct log_2_t<1>
{
  enum { value = 0 };
};


// A lightweight data structure that maintains a vector of mutexes
// corresponding, each one, to N elements of another vector
// Any class modeling the Mutex concept can be used as MutexType
// (see https://software.intel.com/en-us/node/506264 )
template <class MutexType, uint8_t N=4>
class vector_locker_t {
public:
  typedef MutexType mutex_t;
  typedef typename MutexType::scoped_lock scoped_lock;

  explicit vector_locker_t(const size_t no_of_elements)
    :mutexes(new mutex_t[(no_of_elements >> shift_amount)+((no_of_elements & low_mask)? 1 : 0)])
  {
  }

  vector_locker_t(const vector_locker_t<MutexType, N>&) = delete;
  vector_locker_t(vector_locker_t<MutexType, N>&&) = default;

  ~vector_locker_t() {
    delete [] mutexes;
  }

  mutex_t& get_mutex_for_el(const size_t idx) const {
    return mutexes[idx >> shift_amount];
  }

private:
  const static uint8_t shift_amount= log_2_t<N>::value;
  const static size_t low_mask= (N-1);

  mutex_t* mutexes;
};

// A null (=no concurrency) vector locker
#include <tbb/null_mutex.h>
using null_vector_locker_t = vector_locker_t<tbb::null_mutex>;

// A lightweight vector locker
#include <tbb/spin_mutex.h>
using spin_vector_locker_t = vector_locker_t<tbb::spin_mutex>;

#endif
