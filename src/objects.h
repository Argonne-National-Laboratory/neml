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
#include "parameters.h"

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
class NEML_EXPORT NEMLObject {
 public:
  virtual ~NEMLObject() {};
};

/// Error if you ask for a parameter that an object doesn't recognize
class NEML_EXPORT UnknownParameter: public std::exception {
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
class NEML_EXPORT WrongTypeError: public std::exception {
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
class NEML_EXPORT ParameterSet {
 public:
  /// Default constructor, needed to push onto stack
  ParameterSet();

  /// Constructor giving object type
  ParameterSet(std::string type);

  /// Return the type of object you're supposed to create
  const std::string & type() const;

  /// Add a generic parameter with no default
  template<typename T>
  void add_parameter(std::string name);

  /// Immediately assign an input of the right type to a parameter
  template <typename T>
  void assign_parameter(std::string name, const T & value);

  /// Add a generic parameter with a default
  template<typename T>
  void add_optional_parameter(std::string name, const T & value);

  template<typename T>
  void add_optional_parameter(std::string name, const std::shared_ptr<T> & value);

  template<typename T>
  void add_optional_parameter(std::string name, const std::vector<std::shared_ptr<typename T::value_type>> & value);

  /// Add a generic parameter with a string literal value
  template<typename T>
  void add_optional_parameter(std::string name, const char * value)
  {
    add_optional_parameter(name, std::string(value));
  }

  /// Get a parameter of the given name and type
  template<typename T>
  T get_parameter(std::string name)
  {
    resolve_objects_();
    auto res = std::dynamic_pointer_cast<ParamValue<T>>(get_base_parameter(name));
    if (!res)
      throw UnknownParameter(type(), name);
    return res->get();
  }

  const std::shared_ptr<ParamBase> & get_base_parameter(std::string name)
  {
    resolve_objects_();
    auto it = params_.find(name);
    if (it == params_.end())
      throw UnknownParameter(type(), name);
    return it->second;
  }

  /// Assign a parameter set to be used to create an object later
  void assign_defered_parameter(std::string name, const ParameterSet & value);

  /// Helper method to get a NEMLObject and cast it to subtype in one go
  template<typename T>
  std::shared_ptr<T> get_object_parameter(std::string name)
  {
    auto res = std::dynamic_pointer_cast<T>(get_parameter<std::shared_ptr<NEMLObject>>(name));
    if (!res)
      throw WrongTypeError();
    return res;
  }

  /// Helper to get a vector of NEMLObjects and cast them to subtype in one go
  template<typename T>
  std::vector<std::shared_ptr<T>> get_object_parameter_vector(std::string name)
  {
    auto ov = get_parameter<std::vector<std::shared_ptr<NEMLObject>>>(name);
    std::vector<std::shared_ptr<T>> nv(ov.size());
    std::transform(std::begin(ov), std::end(ov), std::begin(nv),
                          [](std::shared_ptr<NEMLObject> const & v)
                          {
                            auto res = std::dynamic_pointer_cast<T>(v);
                            if (res == nullptr)
                              throw WrongTypeError();
                            return res;
                          });
    return nv;
  }

  /// Check if this is an actual parameter
  bool is_parameter(std::string name) const;

  /// Get a list of unassigned parameters
  std::vector<std::string> unassigned_parameters();

  /// Check to make sure this parameter set is ready to go
  bool fully_assigned();

protected:
  template <typename T>
  void assign_parameter_helper(const std::shared_ptr<ParamBase> & p, const T & value);

  template <typename T>
  void assign_parameter_helper(const std::shared_ptr<ParamBase> & p, const std::shared_ptr<T> & value);

  template <typename T>
  void assign_parameter_helper(const std::shared_ptr<ParamBase> & p, const std::vector<std::shared_ptr<T>> & value);

private:
  /// Run down the chain of deferred objects and actually construct them
  void resolve_objects_();

  std::string type_;

  std::map<std::string, std::shared_ptr<ParamBase>> params_;
  std::map<std::string, ParameterSet> defered_params_;
};

template<typename T>
void ParameterSet::add_parameter(std::string name)
{
  params_[name] = std::make_shared<ParamValue<T>>();
}

// farward declare specializations that permit dropping the std::shared_ptr for convenience
template<>
void ParameterSet::add_parameter<NEMLObject>(std::string name);
template<>
void ParameterSet::add_parameter<std::vector<NEMLObject>>(std::string name);

template <typename T>
void ParameterSet::assign_parameter(std::string name, const T & value)
{
  assign_parameter_helper<T>(get_base_parameter(name), value);
}

template <typename T>
void ParameterSet::add_optional_parameter(std::string name, const T & value)
{
  add_parameter<T>(name);
  assign_parameter<T>(name, value);
}

template <typename T>
void ParameterSet::add_optional_parameter(std::string name, const std::shared_ptr<T> & value)
{
  add_optional_parameter<std::shared_ptr<NEMLObject>>(name, value);
}

template <typename T>
void ParameterSet::add_optional_parameter(std::string name, const std::vector<std::shared_ptr<typename T::value_type>> & value)
{
  add_optional_parameter<std::vector<std::shared_ptr<NEMLObject>>>(name, value);
}

// generic template will work for POD types
template <typename T>
void ParameterSet::assign_parameter_helper(const std::shared_ptr<ParamBase> & p, const T & value)
{
  auto pv = std::dynamic_pointer_cast<ParamValue<T>>(p);
  if (!pv)
    throw WrongTypeError();
  pv->assign(value);
}

// neml object smart pointer overload
template <typename T>
void ParameterSet::assign_parameter_helper(const std::shared_ptr<ParamBase> & p, const std::shared_ptr<T> & value)
{
  // check if we can cast to a std::shared_ptr<NEMLObject>
  auto res = std::dynamic_pointer_cast<NEMLObject>(value);
  if (!res)
    throw WrongTypeError();
  // see if we can hold this value type
  auto pv = std::dynamic_pointer_cast<ParamValue<std::shared_ptr<NEMLObject>>>(p);
  if (!pv)
    throw WrongTypeError();
  pv->assign(res);
}

// neml object smart pointer overload
template <typename T>
void ParameterSet::assign_parameter_helper(const std::shared_ptr<ParamBase> & p, const std::vector<std::shared_ptr<T>> & value)
{
  // see if we can hold this value type
  auto pv = std::dynamic_pointer_cast<ParamValue<std::vector<std::shared_ptr<NEMLObject>>>>(p);
  if (!pv)
    throw WrongTypeError();

  // up-cast to NEMLObject pointers
  std::vector<std::shared_ptr<NEMLObject>> nv(value.size());
  std::transform(std::begin(value), std::end(value), std::begin(nv),
                        [](const std::shared_ptr<T> & v)
                        {
                          auto res = std::dynamic_pointer_cast<NEMLObject>(v);
                          if (!res)
                            throw WrongTypeError();
                          return res;
                        });
  pv->assign(nv);
}

/// Factory that produces NEMLObjects from ParameterSets
class NEML_EXPORT Factory {
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
class NEML_EXPORT Register {
 public:
  Register()
  {
     Factory::Creator()->register_type(T::type(), &T::initialize, &T::parameters);
  }
};

/// Error to throw if parameters are not completely defined
class NEML_EXPORT UndefinedParameters: public std::exception {
 public:
  UndefinedParameters(std::string name, std::vector<std::string> unassigned) :
      name_(name), unassigned_(unassigned)
  {

  };

  const char* what() const throw()
  {
    std::stringstream ss;

    ss << "Parameter set for object " << name_ << " has undefined parameters:" << std::endl;

    for (auto it = unassigned_.begin(); it != unassigned_.end(); ++it) {
      ss << "\t" << *it << " ";
    }

    return ss.str().c_str();
  }

 private:
  std::string name_;
  std::vector<std::string> unassigned_;
};

/// Error to throw if the class isn't registered
class NEML_EXPORT UnregisteredError: public std::exception {
 public:
  UnregisteredError(std::string name) :
      name_(name), what_("Object named " + name_ + " not registered with factory!")
  {};

  const char * what() const throw ()
  {
    return what_.c_str();
  };

 private:
   std::string name_;
   std::string what_;
};

} //namespace neml

#endif // OBJECTS_H
