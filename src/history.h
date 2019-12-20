#ifndef HISTORY_H
#define HISTORY_H

#include "math/tensors.h"
#include "math/rotations.h"

#include "windows.h"

#include <string>
#include <map>

namespace neml {

/// Enum type of each allowable object
enum StorageType {
  TYPE_VECTOR    = 0,
  TYPE_SCALAR    = 1,
  TYPE_RANKTWO   = 2,
  TYPE_SYMMETRIC = 3,
  TYPE_SKEW      = 4,
  TYPE_ROT       = 5,
  TYPE_SYMSYM    = 6
};

// Black magic to map a type to the enum
template <class T> constexpr StorageType GetStorageType();
template <> constexpr StorageType GetStorageType<Vector>() {return TYPE_VECTOR;}
template <> constexpr StorageType GetStorageType<RankTwo>() {return TYPE_RANKTWO;}
template <> constexpr StorageType GetStorageType<Symmetric>() {return TYPE_SYMMETRIC;}
template <> constexpr StorageType GetStorageType<Skew>() {return TYPE_SKEW;}
template <> constexpr StorageType GetStorageType<Orientation>() {return TYPE_ROT;}
template <> constexpr StorageType GetStorageType<double>() {return TYPE_SCALAR;}
template <> constexpr StorageType GetStorageType<SymSymR4>() {return TYPE_SYMSYM;}

/// Storage size
const std::map<StorageType,size_t> storage_size =
  { {TYPE_VECTOR,    3},
    {TYPE_SCALAR,    1},
    {TYPE_RANKTWO,   9},
    {TYPE_SYMMETRIC, 6},
    {TYPE_SKEW,      3},
    {TYPE_ROT,       4},
    {TYPE_SYMSYM,    36} };

/// Black magic to map a type to the right size
template <class T> size_t GetStorageSize() {return storage_size.at(GetStorageType<T>());}

/// Map between a type and its derivative
const std::map<StorageType,const std::map<StorageType,StorageType>> derivative_type =
  { {TYPE_SCALAR,
      {{TYPE_SCALAR,    TYPE_SCALAR},
       {TYPE_VECTOR,    TYPE_VECTOR},
       {TYPE_RANKTWO,   TYPE_RANKTWO},
       {TYPE_SYMMETRIC, TYPE_SYMMETRIC},
       {TYPE_SKEW,      TYPE_SKEW},
       {TYPE_ROT,       TYPE_ROT}}}
  };

class NEML_EXPORT History {
 public:
  /// Default constructor (manage own memory)
  History();
  /// Default constructor (option to not manage memory)
  History(bool store);
  /// Copy constructor
  History(const History & other);
  /// Move constructor
  History(const History && other);
  /// Dangerous constructor, only use if you know what you're doing
  History(double * data);
  /// Dangerous constructor, only use if you know what you're doing
  History(const double * data);
  /// Destructor
  virtual ~History();

  /// Copy
  History & operator=(const History & other);
  /// Move
  History & operator=(const History && other);

  /// Explicit deepcopy
  History deepcopy() const;

  /// Do I own my own data?
  bool store() const {return store_;};
  /// Raw data pointer (const)
  const double * rawptr() const {return storage_;};
  /// Raw data pointer (nonconst)
  double * rawptr() {return storage_;};

  /// Set storage to some external pointer
  void set_data(double * input);
  /// Copy data from some external pointer
  void copy_data(const double * const input);
  /// Size of storage required
  size_t size() const;

  /// Convert to store
  void make_store();

  /// Add a generic object
  template<typename T>
  void add(std::string name)
  {
    add(name, GetStorageType<T>(), GetStorageSize<T>());
  }

  /// Add known object
  void add(std::string name, StorageType type, size_t size);

  /// Helper for template magic
  template<class T>
  struct item_return{ typedef T type; };

  /// Get an item (provide with correct class)
  template<class T>
  typename item_return<T>::type get(std::string name) const
  {
    error_if_not_exists_(name);
    error_if_wrong_type_(name, GetStorageType<T>());
    return T(&(storage_[loc_.at(name)]));
  }

  /// Get the location map
  const std::map<std::string,size_t> & get_loc() const {return loc_;};
  /// Get the type map
  const std::map<std::string,StorageType> & get_type() const {return type_;};
  /// Get the name order
  const std::vector<std::string> & get_order() const {return order_;};

  /// Get the location map
  std::map<std::string,size_t> & get_loc() {return loc_;};
  /// Get the type map
  std::map<std::string,StorageType> & get_type() {return type_;};
  /// Get the order
  std::vector<std::string> & get_order() {return order_;};

  const std::vector<std::string> & items() const {return get_order();};

  /// Resize method
  void resize(size_t inc);

  /// Actually increase internal storage
  void increase_store(size_t newsize);

  /// Multiply everything by a scalar
  void scalar_multiply(double scalar);

  /// Add another history to this one
  History & operator+=(const History & other);

  /// Combine another history object through a union
  History & add_union(const History & other);

  /// Make a blank copy
  History copy_blank(std::vector<std::string> exclude = {}) const;

  /// Copy over the order maps
  void copy_maps(const History & other);

  /// Make zero
  void zero();

  /// Make a History appropriate to hold the derivatives of the indicated items
  template<class T>
  History derivative() const
  {
    StorageType dtype = GetStorageType<T>();

    // Precalculate size that will be needed
    size_t nsize = 0;
    for (auto item : order_) {
      StorageType ctype = type_.at(item);
      StorageType ntype = derivative_type.at(ctype).at(dtype);
      nsize += storage_size.at(ntype);
    }
    
    // Cost of moving memory was actually quite high...
    History deriv;
    deriv.increase_store(nsize);

    for (auto item : order_) {
      StorageType ctype = type_.at(item);
      StorageType ntype = derivative_type.at(ctype).at(dtype);
      deriv.add(item, ntype, storage_size.at(ntype));
    }

    deriv.zero();

    return deriv;
  }

  /// Split a history in two
  History split(std::vector<std::string> sep, bool after = true) const;

  /// Quick function to check to see if something is in the vector
  bool contains(std::string name) const;

 private:
  void error_if_exists_(std::string name) const;
  void error_if_not_exists_(std::string name) const;
  void error_if_wrong_type_(std::string name, StorageType type) const;

 private:
  size_t size_;
  size_t storesize_;
  bool store_;
  double * storage_;

  std::map<std::string,size_t> loc_;
  std::map<std::string,StorageType> type_;
  std::vector<std::string> order_;
};

template<>
struct History::item_return<double>{ typedef double & type;};

/// Special case for a double
template<>
inline History::item_return<double>::type History::get<double>(std::string name) const
{
  error_if_not_exists_(name);
  error_if_wrong_type_(name, GetStorageType<double>());
  return storage_[loc_.at(name)];
}

/// Special case for self derivative
template<>
inline History History::derivative<History>() const
{
  History deriv;

  for (auto i1 : order_) {
    StorageType i1_type = type_.at(i1);
    for (auto i2 : order_) {
      StorageType i2_type = type_.at(i2);
      StorageType ntype = derivative_type.at(i1_type).at(i2_type);
      deriv.add(i1+"_"+i2, ntype, storage_size.at(ntype));
    }
  }

  deriv.zero();

  return deriv;
}

} // namespace neml

#endif // HISTORY_H
