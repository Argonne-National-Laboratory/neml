#ifndef HISTORY_H
#define HISTORY_H

#include "math/tensors.h"

#include <string>
#include <map>

namespace neml {

/// Enum type of each allowable object
enum StorageType {
  TYPE_VECTOR   = 0
};

/// Black magic to map a type to the enum
template <class T> constexpr StorageType GetStorageType();
template <> constexpr StorageType GetStorageType<Vector>() {return TYPE_VECTOR;};

/// Black magic to map a type to the right size
template <class T> constexpr size_t GetStorageSize();
template <> constexpr size_t GetStorageSize<Vector>() {return 3;};

class History {
 public:
  History();
  virtual ~History();
  
  void set_data(double * input);
  size_t size() const;
  
  void add_scalar(std::string name);
  double & get_scalar(std::string name);

  void add_array(std::string name, size_t sz);
  size_t array_size(std::string name);
  double * get_array(std::string name);

  template<typename T>
  void add_object(std::string name)
  {
    loc_.insert(std::pair<std::string,size_t>(name, size_));
    size_ += GetStorageSize<T>();
    type_.insert(std::pair<std::string,StorageType>(name, GetStorageType<T>()));
  };

  template<typename T>
  T get_object(std::string name)
  {
    return T(&(storage_[loc_[name]]));
  }

 private:
  size_t size_;
  double * storage_;

  std::map<std::string,size_t> loc_;
  std::map<std::string,size_t> array_size_;
  std::map<std::string,StorageType> type_;

};

} // namespace neml

#endif // HISTORY_H
