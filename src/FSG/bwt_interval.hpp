//-----------------------------------------------
// Copyright 2015 Yuri Pirola, Marco Previtali
// Released under the GPL
//-----------------------------------------------
#ifndef BWTINTERVAL_HPP
#define BWTINTERVAL_HPP

#include <tuple>

template <typename StorageType>
using GeneralTrIntv = std::tuple<StorageType, StorageType, StorageType>;
// Suffix + Pattern Interval on BWT
using TrIntv = GeneralTrIntv<size_t>;
// Not used now
//using SmallTrIntv = GeneralTrIntv<uint32_t>;

template <typename StorageType>
inline StorageType begins(const GeneralTrIntv<StorageType>& intv) {
  using std::get;
  return get<0>(intv);
}

template <typename StorageType>
inline StorageType ends(const GeneralTrIntv<StorageType>& intv) {
  using std::get;
  return get<1>(intv);
}

template <typename StorageType>
inline StorageType endp(const GeneralTrIntv<StorageType>& intv) {
  using std::get;
  return get<2>(intv);
}

template <typename StorageType>
inline StorageType& begins(GeneralTrIntv<StorageType>& intv) {
  using std::get;
  return get<0>(intv);
}

template <typename StorageType>
inline StorageType& ends(GeneralTrIntv<StorageType>& intv) {
  using std::get;
  return get<1>(intv);
}

template <typename StorageType>
inline StorageType& endp(GeneralTrIntv<StorageType>& intv) {
  using std::get;
  return get<2>(intv);
}


template <typename StorageType>
using GeneralBiIntv = std::tuple<StorageType, StorageType>;
// Suffix Interval on BWT
using BiIntv = GeneralBiIntv<size_t>;

using SmallBiIntv_el_type = uint32_t;
using SmallBiIntv = GeneralBiIntv<SmallBiIntv_el_type>;

template <typename StorageType>
inline StorageType begins(const GeneralBiIntv<StorageType>& intv) {
  using std::get;
  return get<0>(intv);
}

template <typename StorageType>
inline StorageType ends(const GeneralBiIntv<StorageType>& intv) {
  using std::get;
  return get<1>(intv);
}

template <typename StorageType>
inline StorageType& begins(GeneralBiIntv<StorageType>& intv) {
  using std::get;
  return get<0>(intv);
}

template <typename StorageType>
inline StorageType& ends(GeneralBiIntv<StorageType>& intv) {
  using std::get;
  return get<1>(intv);
}


template <typename BiIntvType>
struct BiIntvStorage
{
  using type = typename std::tuple_element<0, BiIntvType>::type;
};

template <typename TriIntvType>
struct TriIntvStorage
{
  using type = typename std::tuple_element<0, TriIntvType>::type;
};

#endif
