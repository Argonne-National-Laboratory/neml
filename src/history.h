#ifndef HISTORY_H
#define HISTORY_H

#include "math/tensors.h"

#include <string>
#include <map>

namespace neml {

/// Enum type of each allowable object
enum StorageType {
  TYPE_VECTOR    = 0,
  TYPE_SCALAR    = 1,
  TYPE_ARRAY     = 2,
  TYPE_RANKTWO   = 3,
  TYPE_SYMMETRIC = 4,
  TYPE_SKEW      = 5
};

/// Black magic to map a type to the enum
template <class T> constexpr StorageType GetStorageType();
template <> constexpr StorageType GetStorageType<Vector>() {return TYPE_VECTOR;};
template <> constexpr StorageType GetStorageType<RankTwo>() {return TYPE_RANKTWO;};
template <> constexpr StorageType GetStorageType<Symmetric>() {return TYPE_SYMMETRIC;};
template <> constexpr StorageType GetStorageType<Skew>() {return TYPE_SKEW;};

/// Black magic to map a type to the right size
template <class T> constexpr size_t GetStorageSize();
template <> constexpr size_t GetStorageSize<Vector>() {return 3;};
template <> constexpr size_t GetStorageSize<RankTwo>() {return 9;};
template <> constexpr size_t GetStorageSize<Symmetric>() {return 6;};
template <> constexpr size_t GetStorageSize<Skew>() {return 3;};

class History {
 public:
  History();
  History(bool store);
  History(const History & other);
  History(History && other);
  virtual ~History();

  /// Copy constructor
  History & operator=(const History & other);
  History & operator=(History && other);

  /// Release ownership of data
  void unown();
  /// Do I own my own data?
  bool store() const {return store_;};
  /// Raw data pointer (const)
  const double * rawptr() const {return storage_;};
  /// Raw data pointer (nonconst)
  double * rawptr() {return storage_;};
  
  /// Set storage to some external pointer
  void set_data(double * input);
  /// Size of storage required
  size_t size() const;
  
  /// Add a scalar parameter
  void add_scalar(std::string name);
  /// Get a scalar parameter
  double & get_scalar(std::string name);
  
  /// Add an arbitrary-sized array
  void add_array(std::string name, size_t sz);
  /// Return the size of an array
  size_t array_size(std::string name);
  /// Return a pointer to the array
  double * get_array(std::string name);

  template<typename T>
  void add_object(std::string name)
  {
    error_if_exists_(name);
    loc_.insert(std::pair<std::string,size_t>(name, size_));
    type_.insert(std::pair<std::string,StorageType>(name, GetStorageType<T>()));
    resize_(GetStorageSize<T>());
  };

  template<typename T>
  T get_object(std::string name)
  {
    error_if_not_exists_(name);
    error_if_wrong_type_(name, GetStorageType<T>());
    return std::move(T(&(storage_[loc_[name]])));
  }

 private:
  void resize_(size_t inc);
  void error_if_exists_(std::string name);
  void error_if_not_exists_(std::string name);
  void error_if_wrong_type_(std::string name, StorageType type);

 private:
  size_t size_;
  bool store_;
  double * storage_;

  std::map<std::string,size_t> loc_;
  std::map<std::string,size_t> array_size_;
  std::map<std::string,StorageType> type_;

};

} // namespace neml

#endif // HISTORY_H