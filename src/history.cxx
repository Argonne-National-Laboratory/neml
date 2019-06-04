#include "history.h"

#include <algorithm>
#include <sstream>

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
    size_(other.size()), store_(true)
{
  storage_ = new double [size_];
  std::copy(other.rawptr(), other.rawptr() + size_, storage_);
}

History::History(History && other) :
    size_(other.size()), store_(other.store())
{
  if (store_) {
    storage_ = new double[size_];
    std::copy(other.rawptr(), other.rawptr() + size_, storage_);
  }
  else {
    storage_ = other.rawptr();
    other.unown();
  }
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

  return *this;
}

History & History::operator=(History && other)
{
  if (size_ != other.size()) {
    throw std::invalid_argument(
        "History objects in assignment operator do not have the same size");
  }

  if (other.store()) {
    std::copy(other.rawptr(), other.rawptr() + size_, storage_);
  }
  else {
    storage_ = other.rawptr();
    other.unown();
  }

  return *this;
}

void History::unown()
{
  store_ = false;
}

void History::set_data(double * input)
{
  storage_ = input;
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
  resize_(1);
}

double & History::get_scalar(std::string name)
{
  error_if_not_exists_(name);
  error_if_wrong_type_(name, TYPE_SCALAR);
  return storage_[loc_[name]]; 
}

void History::add_array(std::string name, size_t sz)
{
  error_if_exists_(name);
  loc_.insert(std::pair<std::string,size_t>(name, size_));
  array_size_.insert(std::pair<std::string,size_t>(name, sz));
  type_.insert(std::pair<std::string,StorageType>(name, TYPE_ARRAY));
  resize_(sz);
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

void History::resize_(size_t inc)
{
  if (store_) {
    double * newstore = new double [size_+inc];
    std::copy(storage_, storage_+size_, newstore);
    delete [] storage_;
    storage_ = newstore;
  }
  size_ += inc;
}

void History::error_if_exists_(std::string name)
{
  if (loc_.count(name) != 0) {
    std::stringstream ss;
    ss << "History variable name " << name << " already stored." << std::endl;
    throw std::runtime_error(ss.str());
  }
}

void History::error_if_not_exists_(std::string name)
{
  if (loc_.count(name) != 1) {
    std::stringstream ss;
    ss << "No history variable named " << name << " is stored." << std::endl;
    throw std::runtime_error(ss.str());
  }
}

void History::error_if_wrong_type_(std::string name, StorageType type)
{
  if (type != type_[name]) {
    std::stringstream ss;
    ss << name << " is not of the type requested." << std::endl;
    throw std::runtime_error(ss.str());
  }
}

} // namespace neml
