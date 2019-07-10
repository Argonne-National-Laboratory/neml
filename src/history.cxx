#include "history.h"

#include <algorithm>
#include <sstream>
#include <exception>

namespace neml {

History::History() :
    size_(0), store_(true)
{
  storage_ = new double [size_];
  zero();
}

History::History(bool store) :
    size_(0), store_(store)
{
  if (store) {
    storage_ = new double [size_];
    zero();
  }
}

History::History(const History & other) :
    size_(other.size()), store_(true)
{
  storage_ = new double [size_];
  std::copy(other.rawptr(), other.rawptr() + size_, storage_);
  copy_maps(other);
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
  copy_maps(other);
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

  copy_maps(other);

  return *this;
}

History & History::operator=(const History && other)
{
  if (size_ != other.size()) {
    throw std::invalid_argument(
        "History objects in assignment operator do not have the same size");
  }

  std::copy(other.rawptr(), other.rawptr() + size_, storage_);
  
  copy_maps(other);

  return *this;
}

History History::deepcopy() const
{
  History nhist;
  nhist.resize(size_);
  std::copy(storage_, storage_+size_, nhist.rawptr());
  nhist.copy_maps(*this);

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

void History::make_store()
{
  if (store_) return;

  store_ = true;
  storage_ = new double [size_];
}

void History::add(std::string name, StorageType type, size_t size)
{
  error_if_exists_(name);
  order_.push_back(name);
  loc_.insert(std::pair<std::string,size_t>(name, size_));
  type_.insert(std::pair<std::string,StorageType>(name, type));
  resize(size);
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

History History::copy_blank(std::vector<std::string> exclude) const
{
  History copy;

  for (auto item : order_) {
    if (std::find(exclude.begin(), exclude.end(), item) != exclude.end())
    {
      continue;
    }
    copy.add(item, type_.at(item), storage_size.at(type_.at(item)));
  }

  copy.zero();

  return copy;
}

void History::copy_maps(const History & other)
{
  loc_.insert(other.get_loc().begin(), other.get_loc().end());
  type_.insert(other.get_type().begin(), other.get_type().end());
  order_.assign(other.get_order().begin(), other.get_order().end());
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

void History::zero()
{
  std::fill(storage_, storage_+size_, 0.0);
}

History History::split(std::vector<std::string> sep)
{
  // Check to see if the groups are contiguous and get the offset
  size_t i;
  for (i = 0; i < sep.size(); i++) {
    if (sep[i] != order_[i]) {
      throw std::runtime_error("History items to separate out must be contiguous!");
    }
  }

  History sub(false);

  // Move the maps
  for (size_t j = i; j < order_.size(); j++) {
    sub.add(order_[j], type_.at(order_[j]), storage_size.at(type_.at(order_[j])));
  }

  // Either copy or just split, depending on if we own data
  if (store_) {
    sub.make_store();
    sub.copy_data(&storage_[i]);
  }
  else {
    sub.set_data(&storage_[i]);
  }
  
  return sub;
}

} // namespace neml
