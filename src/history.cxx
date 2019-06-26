#include "history.h"

#include <algorithm>
#include <sstream>
#include <exception>

namespace neml {

History::History() :
    size_(0), store_(true)
{
  storage_ = new double [size_];
}

History::History(bool store) :
    size_(0), store_(store)
{
  if (store) {
    storage_ = new double [size_];
  }
}

History::History(const History & other) :
    size_(other.size()), store_(other.store())
{
  if (store_) {
    storage_ = new double [size_];
    std::copy(other.rawptr(), other.rawptr() + size_, storage_);
  }
  else {
    storage_ = const_cast<double*>(other.rawptr());
  }
  copy_maps_(other);
}

History::History(const History && other) :
    size_(other.size()), store_(other.store())
{
  if (store_) {
    storage_ = new double[size_];
    std::copy(other.rawptr(), other.rawptr() + size_, storage_);
  }
  else {
    storage_ = const_cast<double*>(other.rawptr());
  }
  copy_maps_(other);
}

History::~History()
{
  if (store_) {
    delete [] storage_;
  }
  storage_ = nullptr;
}

History & History::operator=(const History & other)
{
  if (size_ != other.size()) {
    throw std::invalid_argument(
        "History objects in assignment operator do not have the same size");
  }

  if (this != &other) {
    std::copy(other.rawptr(), other.rawptr() + size_, storage_);
  }

  copy_maps_(other);

  return *this;
}

History & History::operator=(const History && other)
{
  if (size_ != other.size()) {
    throw std::invalid_argument(
        "History objects in assignment operator do not have the same size");
  }

  std::copy(other.rawptr(), other.rawptr() + size_, storage_);
  
  copy_maps_(other);

  return *this;
}

History History::deepcopy() const
{
  History nhist;
  nhist.resize(size_);
  std::copy(storage_, storage_+size_, nhist.rawptr());

  nhist.get_loc().insert(loc_.begin(), loc_.end());
  nhist.get_array_size().insert(array_size_.begin(), array_size_.end());
  nhist.get_type().insert(type_.begin(), type_.end());

  return nhist;
}

void History::set_data(double * input)
{
  storage_ = input;
}

void History::copy_data(const double * const input)
{
  std::copy(input, input+size(), storage_);
}

size_t History::size() const 
{
  return size_;
}

void History::add_scalar(std::string name)
{
  error_if_exists_(name);
  loc_.insert(std::pair<std::string,size_t>(name, size_));
  type_.insert(std::pair<std::string,StorageType>(name, TYPE_SCALAR));
  resize(1);
}

double & History::get_scalar(std::string name)
{
  error_if_not_exists_(name);
  error_if_wrong_type_(name, TYPE_SCALAR);
  return storage_[loc_[name]]; 
}

const double & History::get_scalar(std::string name) const
{
  error_if_not_exists_(name);
  error_if_wrong_type_(name, TYPE_SCALAR);
  return storage_[loc_.at(name)]; 
}

void History::add_array(std::string name, size_t sz)
{
  error_if_exists_(name);
  loc_.insert(std::pair<std::string,size_t>(name, size_));
  array_size_.insert(std::pair<std::string,size_t>(name, sz));
  type_.insert(std::pair<std::string,StorageType>(name, TYPE_ARRAY));
  resize(sz);
}

size_t History::array_size(std::string name)
{
  error_if_not_exists_(name);
  return array_size_[name];
}

double * History::get_array(std::string name)
{
  error_if_not_exists_(name);
  error_if_wrong_type_(name, TYPE_ARRAY);
  return &(storage_[loc_[name]]);
}

const double * History::get_array(std::string name) const
{
  error_if_not_exists_(name);
  error_if_wrong_type_(name, TYPE_ARRAY);
  return &(storage_[loc_.at(name)]);
}

void History::resize(size_t inc)
{
  if (store_) {
    double * newstore = new double [size_+inc];
    std::copy(storage_, storage_+size_, newstore);
    delete [] storage_;
    storage_ = newstore;
  }
  size_ += inc;
}

void History::scalar_multiply(double scalar)
{
  for (size_t i = 0; i < size_; i++) {
    storage_[i] *= scalar; 
  }
}

History & History::operator+=(const History & other)
{
  if (size() != other.size()) {
    throw std::runtime_error("Histories to be added do not have the same size!");
  }
  for (size_t i = 0; i < size_; i++) {
    storage_[i] += other.rawptr()[i];
  }
  return *this;
}

void History::copy_maps_(const History & other)
{
  loc_.insert(other.get_loc().begin(), other.get_loc().end());
  array_size_.insert(other.get_array_size().begin(), other.get_array_size().end());
  type_.insert(other.get_type().begin(), other.get_type().end());
}

void History::error_if_exists_(std::string name) const
{
  if (loc_.count(name) != 0) {
    std::stringstream ss;
    ss << "History variable name " << name << " already stored." << std::endl;
    throw std::runtime_error(ss.str());
  }
}

void History::error_if_not_exists_(std::string name) const
{
  if (loc_.count(name) != 1) {
    std::stringstream ss;
    ss << "No history variable named " << name << " is stored." << std::endl;
    throw std::runtime_error(ss.str());
  }
}

void History::error_if_wrong_type_(std::string name, StorageType type) const
{
  if (type != type_.at(name)) {
    std::stringstream ss;
    ss << name << " is not of the type requested." << std::endl;
    throw std::runtime_error(ss.str());
  }
}

} // namespace neml
