#include "history.h"

#include <algorithm>

namespace neml {

History::History() :
    size_(0), store_(true)
{
  storage_ = new double [0];
}

History::History(bool store) :
    size_(0), store_(store)
{
  if (store) {
    storage_ = new double [0];
  }
}

History::~History()
{
  if (store_) {
    delete [] storage_;
  }
  storage_ = nullptr;
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
  loc_.insert(std::pair<std::string,size_t>(name, size_));
  resize_(1);
}

double & History::get_scalar(std::string name)
{
  return storage_[loc_[name]]; 
}

void History::add_array(std::string name, size_t sz)
{
  loc_.insert(std::pair<std::string,size_t>(name, size_));
  array_size_.insert(std::pair<std::string,size_t>(name, sz));
  resize_(sz);
}

size_t History::array_size(std::string name)
{
  return array_size_[name];
}

double * History::get_array(std::string name)
{
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

} // namespace neml
