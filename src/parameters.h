#pragma once

#include <memory>
#include <type_traits>

#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"

namespace pybind11
{
class object;
}

namespace neml
{

// Forward declaration
class NEMLObject;

/// base class for all parameter values
class ParamBase
{
public:
  ParamBase() : set_(false) {}
  virtual ~ParamBase() = default;

  template <typename T>
  bool isType() const;

  bool isSet() const { return set_; }

  virtual void setFromXML(rapidxml::xml_node<char> *) = 0;
  virtual void setFromPY(pybind11::object &) = 0;

  virtual std::string typeString() = 0;

protected:
  bool set_;

private:
  template <typename T>
  bool isTypeImpl(T *) const;

  template <typename T>
  bool isTypeImpl(std::shared_ptr<T> *) const;

  template <typename T>
  bool isTypeImpl(std::vector<std::shared_ptr<T>> *) const;
};

/// templated parameter value class that can hold a value of type T
template <typename T>
class ParamValue : public ParamBase
{
public:
  ParamValue() : ParamBase() {}
  virtual ~ParamValue() = default;

  void assign(const T &value)
  {
    set_ = true;
    value_ = value;
  }

  T get() { return value_; }

  virtual void setFromXML(rapidxml::xml_node<char> *);
  virtual void setFromPY(pybind11::object &);

  virtual std::string typeString() { return "<Unknown>"; };

private:
  T value_;
};

template <typename T>
bool ParamBase::isType() const
{
  // turn this into overloading using a dummy parameter
  return isTypeImpl(static_cast<T *>(0));
}

template <typename T>
bool ParamBase::isTypeImpl(T *) const
{
  return dynamic_cast<const ParamValue<T> *>(this) != nullptr;
}

template <typename T>
bool ParamBase::isTypeImpl(std::shared_ptr<T> *) const
{
  return dynamic_cast<const ParamValue<std::shared_ptr<NEMLObject>> *>(this) != nullptr &&
         std::is_base_of<NEMLObject, T>::value;
}

template <typename T>
bool ParamBase::isTypeImpl(std::vector<std::shared_ptr<T>> *) const
{
  return dynamic_cast<const ParamValue<std::vector<std::shared_ptr<NEMLObject>>> *>(this) != nullptr &&
         std::is_base_of<NEMLObject, T>::value;
}

/// the unspecialized version will error out
template <typename T>
void ParamValue<T>::setFromXML(rapidxml::xml_node<char> *)
{
  throw std::runtime_error("Unrecognized object type!");
}

} // namespace neml
