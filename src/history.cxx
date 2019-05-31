#include "history.h"

namespace neml {

History::History() :
    size_(0)
{
}

History::~History()
{
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
  size_ += 1;
}

double & History::get_scalar(std::string name)
{
  return storage_[loc_[name]]; 
}

void History::add_array(std::string name, size_t sz)
{
  loc_.insert(std::pair<std::string,size_t>(name, size_));
  array_size_.insert(std::pair<std::string,size_t>(name, sz));
  size_ += sz;
}

size_t History::array_size(std::string name)
{
  return array_size_[name];
}

double * History::get_array(std::string name)
{
  return &(storage_[loc_[name]]);
}

} // namespace neml
