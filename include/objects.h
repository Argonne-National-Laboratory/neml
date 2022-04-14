#ifndef OBJECTS_H
#define OBJECTS_H

#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <memory>
#include <functional>
#include <stdexcept>
#include <algorithm>

#include "windows.h"
#include "nemlerror.h"

namespace neml {

class NEMLObject;

/// Typedef for slip systems
typedef std::vector<std::pair<std::vector<int>,std::vector<int>>> list_systems;
/// Typedef for twin systems
typedef
    std::vector<std::tuple<std::vector<int>,std::vector<int>,std::vector<int>,std::vector<int>>>
    twin_systems;

/// We can avoid this with proper C++14, will need ifdefs
template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

/// This black magic lets us store parameters in a unified map (note that this is a replacement for boost::variant
/// to get rid of a boost dependency which is causing us grief due to spurious memory leaks in boost)
class param_type
{
public:
  template <typename T>
  T get() { return data<T>(); }

  param_type() {}
  param_type(const char * val);
  template <typename T>
  param_type(const T & val) { data<T>() = val; }
  template <typename T>
  param_type(const std::shared_ptr<T> & val) { data<std::shared_ptr<NEMLObject>>() = val; }

  template <typename T>
  const param_type & operator = (const T & val) { data<T>() = val; return *this; }

protected:
  template <typename T>
  T & data();

private:
  // it'd be tempting to use a union here, but then the compiler generated copy constructor would not work anymore
  double double_;
  int int_;
  bool bool_;
  std::vector<double> vec_double_;
  std::shared_ptr<NEMLObject> neml_object_;
  std::vector<std::shared_ptr<NEMLObject>> vec_neml_object_;
  std::string string_;
  list_systems list_systems_;
  twin_systems twin_systems_;
  std::size_t size_t_;
  std::vector<std::size_t> vec_size_t_;
};

template <> inline double & param_type::data<double>() { return double_; }
template <> inline int & param_type::data<int>() { return int_; }
template <> inline bool & param_type::data<bool>() { return bool_; }
template <> inline std::vector<double> & param_type::data<std::vector<double>>() { return vec_double_; }
template <> inline std::shared_ptr<NEMLObject> & param_type::data<std::shared_ptr<NEMLObject>>() { return neml_object_; }
template <> inline std::vector<std::shared_ptr<NEMLObject>> & param_type::data<std::vector<std::shared_ptr<NEMLObject>>>() { return vec_neml_object_; }
template <> inline std::string & param_type::data<std::string>() { return string_; }
template <> inline list_systems & param_type::data<list_systems>() { return list_systems_; }
template <> inline twin_systems & param_type::data<twin_systems>() { return twin_systems_; }
template <> inline std::size_t & param_type::data<std::size_t>() { return size_t_; }
template <> inline std::vector<std::size_t> & param_type::data<std::vector<std::size_t>>() { return vec_size_t_; }

inline param_type::param_type(const char * val) { data<std::string>() = *val; }


/// This is the enum name we assign to each type for the "external" interfaces
/// to use in reconstructing a type from data
enum ParamType {
  TYPE_DOUBLE           = 0,
  TYPE_INT              = 1,
  TYPE_BOOL             = 2,
  TYPE_VEC_DOUBLE       = 3,
  TYPE_NEML_OBJECT      = 4,
  TYPE_VEC_NEML_OBJECT  = 5,
  TYPE_STRING           = 6,
  TYPE_SLIP             = 7,
  TYPE_SIZE_TYPE        = 8,
  TYPE_VEC_SIZE_TYPE    = 9,
  TYPE_TWIN             = 10
};
// This black magic lets us map the actual type of each parameter to the enum
template <class T> constexpr ParamType GetParamType();
template <> constexpr ParamType GetParamType<double>() {return TYPE_DOUBLE;}
template <> constexpr ParamType GetParamType<int>() {return TYPE_INT;}
template <> constexpr ParamType GetParamType<bool>() {return TYPE_BOOL;}
template <> constexpr ParamType GetParamType<std::vector<double>>()
{return TYPE_VEC_DOUBLE;}
template <> constexpr ParamType GetParamType<std::shared_ptr<NEMLObject>>()
{return TYPE_NEML_OBJECT;}
template <> constexpr ParamType GetParamType<NEMLObject>()
{return TYPE_NEML_OBJECT;}
template <> constexpr ParamType GetParamType<std::vector<std::shared_ptr<NEMLObject>>>()
{return TYPE_VEC_NEML_OBJECT;}
template <> constexpr ParamType GetParamType<std::vector<NEMLObject>>()
{return TYPE_VEC_NEML_OBJECT;}
template <> constexpr ParamType GetParamType<std::string>() {return TYPE_STRING;}
template <> constexpr ParamType GetParamType<list_systems>()
{return TYPE_SLIP;}
template <> constexpr ParamType GetParamType<twin_systems>()
{return TYPE_TWIN;}
template <> constexpr ParamType GetParamType<size_t>() {return TYPE_SIZE_TYPE;}
template <> constexpr ParamType GetParamType<std::vector<size_t>>() {return
TYPE_VEC_SIZE_TYPE;}

/// Error if you ask for a parameter that an object doesn't recognize
class NEML_EXPORT UnknownParameter: public NEMLError {
 public:
  UnknownParameter(std::string object, std::string name) :
      NEMLError("Object of type " + object + " has no parameter " + name)
  {
  }
};

/// Error to call if you try a bad cast
class NEML_EXPORT WrongTypeError: public NEMLError {
 public:
  WrongTypeError() :
      NEMLError("Cannot convert object to the correct type!")
  {
  }
};

/// Error to throw if parameters are not completely defined
class NEML_EXPORT UndefinedParameters: public NEMLError {
 public:
  UndefinedParameters(std::string name, std::vector<std::string> unassigned) :
      NEMLError("Parameter set for object " + name + " has undefined parameters"),
      unassigned_(unassigned)
  {
  }

 private:
  std::vector<std::string> unassigned_;
};

/// Error to throw if the class isn't registered
class NEML_EXPORT UnregisteredError: public NEMLError {
 public:
  UnregisteredError(std::string name) :
      NEMLError("Object name " + name + " is not registered with the factory")
  {
  };
};

/// Parameters for objects created through the NEMLObject interface
class NEML_EXPORT ParameterSet {
 public:
  /// Default constructor, needed to push onto stack
  ParameterSet();
  /// Constructor giving object type
  ParameterSet(std::string type);

  virtual ~ParameterSet();

  /// Return the type of object you're supposed to create
  const std::string & type() const;

  /// Add a generic parameter with no default
  template<typename T>
  void add_parameter(std::string name)
  {
    param_names_.push_back(name);
    param_types_[name] = GetParamType<T>();
  }

  /// Immediately assign an input of the right type to a parameter
  void assign_parameter(std::string name, param_type value)
  {
    if (std::find(param_names_.begin(), param_names_.end(), name) == param_names_.end()) {
      throw UnknownParameter(type(), name);
    }
    params_[name] = value;
  }

  /// Add a generic parameter with a default
  template<typename T>
  void add_optional_parameter(std::string name, param_type value)
  {
    add_parameter<T>(name);
    assign_parameter(name, value);
  }

  /// Get a parameter of the given name and type
  template<typename T>
  T get_parameter(std::string name)
  {
    resolve_objects_();
    return params_[name].get<T>();
  }

  /// Assign a parameter set to be used to create an object later
  void assign_defered_parameter(std::string name, ParameterSet value);

  /// Helper method to get a NEMLObject and cast it to subtype in one go
  template<typename T>
  std::shared_ptr<T> get_object_parameter(std::string name)
  {
    auto res = std::dynamic_pointer_cast<T>(get_parameter<std::shared_ptr<NEMLObject>>(name));
    if (res == nullptr) {
      throw WrongTypeError();
    }
    else {
      return res;
    }
  }

  /// Helper to get a vector of NEMLObjects and cast them to subtype in one go
  template<typename T>
  std::vector<std::shared_ptr<T>> get_object_parameter_vector(std::string name)
  {
    std::vector<std::shared_ptr<NEMLObject>> ov =
        get_parameter<std::vector<std::shared_ptr<NEMLObject>>>(name);
    std::vector<std::shared_ptr<T>> nv(ov.size());
    std::transform(std::begin(ov), std::end(ov), std::begin(nv),
                          [](std::shared_ptr<NEMLObject> const & v)
                          {
                            auto res = std::dynamic_pointer_cast<T>(v);

                            if (res == nullptr) {
                              throw WrongTypeError();
                            }
                            else {
                              return res;
                            }
                          });
    return nv;
  }

  /// Get the type of parameter
  ParamType get_object_type(std::string name);

  /// Check if this is an actual parameter
  bool is_parameter(std::string name) const;

  /// Get a list of unassigned parameters
  std::vector<std::string> unassigned_parameters();

  /// Check to make sure this parameter set is ready to go
  bool fully_assigned();
  
  /// Name getter
  const std::vector<std::string> & param_names() const {return param_names_;};

 private:
  /// Run down the chain of deferred objects and actually construct them
  void resolve_objects_();

  std::string type_;

  std::vector<std::string> param_names_;
  std::map<std::string, ParamType> param_types_;
  std::map<std::string, param_type> params_;
  std::map<std::string, ParameterSet> defered_params_;
};

/// Factory that produces NEMLObjects from ParameterSets
class NEML_EXPORT Factory {
 public:
  Factory();

  /// Provide a valid parameter set for the object type
  ParameterSet provide_parameters(std::string type);

  /// Create an object from the parameter set
  std::shared_ptr<NEMLObject> create(ParameterSet & params);

  /// Alternate the makes a unique_ptr
  std::unique_ptr<NEMLObject> create_unique(ParameterSet & params);

  /// Create and cast an object to a type
  template<typename T>
  std::shared_ptr<T> create(ParameterSet & params)
  {
    auto res = std::dynamic_pointer_cast<T>(create(params));
    if (res == nullptr) {
      throw WrongTypeError();
    }
    else {
      return res;
    }
  }

  /// Create and cast an object to a type as a unique_ptr
  template<typename T>
  std::unique_ptr<T> create_unique(ParameterSet & params)
  {
    auto res = std::unique_ptr<T>(dynamic_cast<T*>(create_unique(params).release()));
    if (res == nullptr) {
      throw WrongTypeError();
    }
    else {
      return res;
    }
  }

  /// Register a type with an identifier, create method, and parameter set
  void register_type(std::string type,
                     std::function<std::unique_ptr<NEMLObject>(ParameterSet &)> creator,
                     std::function<ParameterSet()> setup);

  /// Static factor instance
  static Factory * Creator();

 private:
  std::map<std::string, std::function<std::unique_ptr<NEMLObject>(ParameterSet &)>> creators_;
  std::map<std::string, std::function<ParameterSet()>> setups_;
};

/// Little object used for auto registration
template<typename T>
class NEML_EXPORT Register {
 public:
  Register()
  {
    Factory::Creator()->register_type(T::type(), &T::initialize, &T::parameters);
  }
};

/// NEMLObject base calls for serialization
class NEML_EXPORT NEMLObject {
 public:
  NEMLObject(ParameterSet & params);
  virtual ~NEMLObject() {}; // IDK, clang whines
  
  /// Return the current parameter set, including any updates from construction
  // Most of the time the default will be sufficient, but we may need to
  // override for a few objects
  virtual ParameterSet & current_parameters();

  /// Serialize an object to ASCII XML
  std::string serialize(std::string object_name = "object", std::string top_node
                        = "");

 protected:
  ParameterSet current_params_;
};

} //namespace neml

#endif // OBJECTS_H
