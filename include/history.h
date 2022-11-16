#ifndef HISTORY_H
#define HISTORY_H

#include "math/tensors.h"
#include "math/rotations.h"

#include "windows.h"

#include <string>
#include <unordered_map>

namespace neml {

/// Enum type of each allowable object
enum StorageType {
  TYPE_VECTOR    = 0,
  TYPE_SCALAR    = 1,
  TYPE_RANKTWO   = 2,
  TYPE_SYMMETRIC = 3,
  TYPE_SKEW      = 4,
  TYPE_ROT       = 5,
  TYPE_SYMSYM    = 6,
  TYPE_BLANK     = 7
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
const std::unordered_map<StorageType,size_t,std::hash<int>> storage_size =
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
const std::unordered_map<StorageType,const std::unordered_map<StorageType,StorageType,std::hash<int>>,std::hash<int>> derivative_type =
  { {TYPE_SCALAR,
      {{TYPE_SCALAR,    TYPE_SCALAR},
       {TYPE_VECTOR,    TYPE_VECTOR},
       {TYPE_RANKTWO,   TYPE_RANKTWO},
       {TYPE_SYMMETRIC, TYPE_SYMMETRIC},
       {TYPE_SKEW,      TYPE_SKEW},
       {TYPE_ROT,       TYPE_ROT}}},

    {TYPE_SYMMETRIC,
      {{TYPE_SCALAR,    TYPE_SYMMETRIC},
       {TYPE_SYMMETRIC, TYPE_SYMSYM}}}
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

  /// Get a pointer to the raw location of an item
  double * get_data(std::string name)
  {
    error_if_not_exists_(name);
    return &(storage_[loc_.at(name)]);
  }

  /// Get the location map
  const std::unordered_map<std::string,size_t> & get_loc() const {return loc_;};
  /// Get the type map
  const std::unordered_map<std::string,StorageType> & get_type() const {return type_;};
  /// Get the name order
  const std::vector<std::string> & get_order() const {return order_;};

  /// Get the location map
  std::unordered_map<std::string,size_t> & get_loc() {return loc_;};
  /// Get the type map
  std::unordered_map<std::string,StorageType> & get_type() {return type_;};
  /// Get the order
  std::vector<std::string> & get_order() {return order_;};

  /// Return all the items in this object
  const std::vector<std::string> & items() const {return get_order();};

  /// Helper to get the size of a particular object
  size_t size_of_entry(std::string name) const;

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
  History & zero();

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

  /// Derivative with respect to a different history
  History history_derivative(const History & other) const;

  /// Split a history in two
  History split(std::vector<std::string> sep, bool after = true) const;

  /// Extract a subset
  History subset(std::vector<std::string> vars) const;

  /// Reorder based on the provided list of names
  History & reorder(std::vector<std::string> names);

  /// Quick function to check to see if something is in the vector
  inline bool contains(std::string name) const { return loc_.find(name) !=
    loc_.end();};

  /// Postmultiply by various objects
  History postmultiply(const SymSymR4 & T);

  /// This unravels a history derivative into row major storage
  void unravel_hh(const History & base, double * const array);
  
  /// Starting location of an entry
  double * start_loc(std::string name);
  
  /// Nicely formatted names for the flat storage
  std::vector<std::string> formatted_names() const;

 private:
  void error_if_exists_(std::string name) const;
  void error_if_not_exists_(std::string name) const;
  void error_if_wrong_type_(std::string name, StorageType type) const;

 private:
  size_t size_;
  size_t storesize_;
  bool store_;
  double * storage_;

  std::unordered_map<std::string,size_t> loc_;
  std::unordered_map<std::string,StorageType> type_;
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
  return history_derivative(*this);
}

/// NEMLObject that maintains some internal state variables
class NEML_EXPORT HistoryNEMLObject: public NEMLObject {
 public:
  HistoryNEMLObject(ParameterSet & params);
  virtual ~HistoryNEMLObject() {};

  /// Setup the internal state
  virtual void populate_hist(History & h) const = 0;
  /// Initialize the history
  virtual void init_hist(History & h) const = 0;
  /// This should be replaced at some point
  virtual size_t nhist() const;
  /// Dangerous nhist which assumes you cached 
  virtual size_t nh() const;

  /// Setup a flat vector history
  virtual void init_store(double * const h) const;
  /// Now just = nhist
  virtual size_t nstore() const;

  /// Set the object prefix
  void set_variable_prefix(std::string prefix);
  /// Get the prefix
  std::string get_variable_prefix() const;
  /// Assemble a complete variable name
  std::string prefix(std::string basename) const;
  /// Assemble a complete derivative name
  std::string dprefix(std::string a, std::string b) const;

  /// Cache the history to avoid recreating it every time
  void cache_history_();

  /// Quickly setup history
  History gather_history_(double * data) const;
  History gather_history_(const double * data) const;
  History gather_blank_history_() const;

 protected:
  std::string prefix_;
  History stored_hist_;
  
 private:
  bool cached_;
  size_t ncache_;
};

} // namespace neml

#endif // HISTORY_H
