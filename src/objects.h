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

#include "boost/variant.hpp"

namespace neml {

/// Typedef for slip systems
typedef std::vector<std::pair<std::vector<int>,std::vector<int>>> list_systems;

/// We can avoid this with proper C++14, will need ifdefs
template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

/// NEMLObjects are current pretty useless.  However, they are a hook
/// for future work on serialization.
class NEMLObject {
 public:
  virtual ~NEMLObject() {};
};

// This version supports the following types of objects as parameters:
//    double
//    int
//    bool
//    vector<double>
//    NEMLObject
//    vector<NEMLObject>
//    string
//    vector<pair<vector<int>,vector<int>>, i.e. groups of slip systems

/// This black magic lets us store parameters in a unified map
typedef boost::variant<double, int, bool, std::vector<double>, 
        std::shared_ptr<NEMLObject>,std::vector<std::shared_ptr<NEMLObject>>,
        std::string, list_systems> param_type;
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
  TYPE_SLIP             = 7
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

/// Error if you ask for a parameter that an object doesn't recognize
class UnknownParameter: public std::exception {
 public:
  UnknownParameter(std::string object, std::string name) :
      object_(object), name_(name)
  {
  
  }

  const char* what() const throw()
  {
    std::stringstream ss;
    
    ss << "Object of type " << object_ << " has no parameter "
        << name_ << "!";

    return ss.str().c_str();
  }

 private:
  std::string object_, name_;
};

/// Error to call if you try a bad cast
class WrongTypeError: public std::exception {
 public:
  WrongTypeError()
  {

  };

  const char * what() const throw ()
  {
    std::stringstream ss;

    ss << "Cannot convert object to the correct type!";

    return ss.str().c_str();
  }
};

/// Parameters for objects created through the NEMLObject interface
class ParameterSet {
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
    return boost::get<T>(params_[name]);
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
class Factory {
 public:
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
class Register {
 public:
  Register()
  {
    Factory::Creator()->register_type(T::type(), &T::initialize, &T::parameters);
  }
};

/// Error to throw if parameters are not completely defined
class UndefinedParameters: public std::exception {
 public:
  UndefinedParameters(ParameterSet & params) :
      params_(params)
  {

  };

  const char* what() const throw()
  {
    std::stringstream ss;

    ss << "Parameter set has undefined parameters:" << std::endl;
    auto ups = params_.unassigned_parameters();
    for (auto it = ups.begin(); it != ups.end(); ++it) {
      ss << *it << " ";
    }

    return ss.str().c_str();
  }

 private:
  ParameterSet & params_;
};

/// Error to throw if the class isn't registered
class UnregisteredError: public std::exception {
 public:
  UnregisteredError(std::string name) :
      name_(name)
  {

  };

  const char * what() const throw ()
  {
    std::stringstream ss;

    ss << "Object named " << name_ << " not registered with factory!";

    return ss.str().c_str();
  };

 private:
  std::string name_;
};

} //namespace neml

#endif // OBJECTS_H
