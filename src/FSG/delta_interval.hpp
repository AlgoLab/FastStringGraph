//-----------------------------------------------
// Copyright 2015 Yuri Pirola, Marco Previtali
// Released under the GPL
//-----------------------------------------------
#ifndef DELTA_INTERVAL_HPP
#define DELTA_INTERVAL_HPP

#include <deque>
#include <cassert>

#include "bwt_interval.hpp"

template <typename VALUES_TYPE, typename CONVERT_TYPE>
class DeltaInterval {
private:

  using container_t = std::deque<CONVERT_TYPE>;

public:

  class enumerator_t {
  private:

    typename container_t::const_iterator it_;
    const typename container_t::const_iterator it_end;
    VALUES_TYPE start_;

    VALUES_TYPE _decode_next() {
      VALUES_TYPE _next=0;
      CONVERT_TYPE _continue = CONTINUE_MASK;
      CONVERT_TYPE _loopn = 0;
      do {
        _next |= (static_cast<VALUES_TYPE>(*it_ & VALUE_MASK) << (_loopn++ * SHIFT_AMOUNT));
        _continue = (*it_ & CONTINUE_MASK);
        ++it_;
      } while(_continue);
      return _next;
    }

    explicit enumerator_t(typename container_t::const_iterator&& cbegin,
                          typename container_t::const_iterator&& cend)
    :it_(cbegin), it_end(cend), start_(0)
    {
    }

  public:

    bool has_next() const {
      return it_ != it_end;
    }

    GeneralTrIntv<VALUES_TYPE> next() {
      GeneralTrIntv<VALUES_TYPE> _ret;
      begins(_ret) = start_ + _decode_next();
      ends(_ret) = begins(_ret) + _decode_next();
      endp(_ret) = ends(_ret) + _decode_next();
      start_ = endp(_ret);
      return _ret;
    }

    friend class DeltaInterval<VALUES_TYPE, CONVERT_TYPE>;

  };

  explicit DeltaInterval()
    :start_(0)
  {
  }

  bool add(const GeneralTrIntv<VALUES_TYPE>& nval)
  {
    _push_val(begins(nval));
    _push_val(ends(nval));
    _push_val(endp(nval));
    return true;
  }

  enumerator_t enumerate() const {
    return enumerator_t(values_.cbegin(), values_.cend());
  }

private:

  static const CONVERT_TYPE   SHIFT_AMOUNT  = sizeof(CONVERT_TYPE)*8 -1;
  static const CONVERT_TYPE   CONTINUE_MASK = (static_cast<CONVERT_TYPE>(0x1) << SHIFT_AMOUNT);
  static const CONVERT_TYPE   VALUE_MASK    = CONTINUE_MASK -1;

  container_t values_;
  VALUES_TYPE start_;

  void _push_val(VALUES_TYPE _x)
  {
    assert( _x >= start_);
    _x -= start_;
    start_ += _x;
    do {
      CONVERT_TYPE _y = _x & VALUE_MASK;
      _x >>= SHIFT_AMOUNT;
      if(_x) _y |= CONTINUE_MASK;
      values_.push_back(_y);
    } while(_x);
  }

};


#endif
