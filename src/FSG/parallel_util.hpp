//-----------------------------------------------
// Copyright 2015 Yuri Pirola, Marco Previtali
// Released under the GPL
//-----------------------------------------------
#ifndef PARALLEL_UTIL_HPP
#define PARALLEL_UTIL_HPP

#include <tbb/concurrent_queue.h>

template <class T>
struct default_construct_t {
  template <typename... Args>
  T* operator()(Args... args) const {
    return new T(args...);
  }
};

template <class T>
struct default_destroy_t {
  void operator()(T* p_obj) const {
    if (p_obj != nullptr) delete p_obj;
  }
};

template <class ObjType,
          typename Ctor=std::function<ObjType*()>,
          typename Dtor=std::function<void(ObjType*)> >
class obj_pool_t {
public:
  using object_type= ObjType;
  using constructor_type= Ctor;
  using destructor_type= Dtor;

  obj_pool_t(const constructor_type& ctor=default_construct_t<object_type>(),
             const destructor_type& dtor=default_destroy_t<object_type>())
    :ctor_(ctor), dtor_(dtor)
  {
  }

  ~obj_pool_t() {
    object_type* p_obj= nullptr;
    while (pool_.try_pop(p_obj)) {
      dtor_(p_obj);
    }
  }

  object_type* get() {
    object_type* p_obj= nullptr;
    if (!pool_.try_pop(p_obj))
      p_obj= ctor_();
    return p_obj;
  }

  void give_back(object_type* const p_obj) {
    pool_.push(p_obj);
  }


private:

  const constructor_type ctor_;
  const destructor_type dtor_;
  tbb::concurrent_queue<object_type*> pool_ {};

};


#endif
