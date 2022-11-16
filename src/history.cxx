#include "history.h"

#include <algorithm>
#include <sstream>

namespace neml {

History::History() :
    size_(0), storesize_(0), store_(true)
{
  storage_ = new double [storesize_];
  zero();
}

History::History(bool store) :
    size_(0), storesize_(0), store_(store)
{
  if (store) {
    storage_ = new double [storesize_];
    zero();
  }
}

History::History(const History & other) :
    size_(other.size()), storesize_(other.size()), store_(other.store())
{
  if (store_) {
    storage_ = new double[storesize_];
    std::copy(other.rawptr(), other.rawptr() + size_, storage_);
  }
  else {
    storage_ = const_cast<double*>(other.rawptr());
  }
  copy_maps(other);
}

History::History(const History && other) :
    size_(other.size()), storesize_(other.size()), store_(other.store())
{
  if (store_) {
    storage_ = new double[storesize_];
    std::copy(other.rawptr(), other.rawptr() + size_, storage_);
  }
  else {
    storage_ = const_cast<double*>(other.rawptr());
  }
  copy_maps(other);
}

History::History(double * data) :
    size_(0), storesize_(0), store_(false)
{
  storage_ = data;
}

History::History(const double * data) :
    size_(0), storesize_(0), store_(false)
{
  storage_ = const_cast<double*>(data);
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

History & History::add_union(const History & other)
{
  // This assumes no overlap
  order_.insert(order_.end(), other.get_order().begin(), other.get_order().end());
  type_.insert(other.get_type().begin(), other.get_type().end());
  loc_.insert(other.get_loc().begin(), other.get_loc().end());
  size_t offset = size();
  resize(other.size());
  for (auto & var : other.get_order()) {
    loc_[var] += offset;
  }

  std::copy(other.rawptr(), other.rawptr()+other.size(), storage_+offset);

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
  if (store_) {
    store_ = false;
    storesize_ = 0;
    delete [] storage_;
    storage_ = nullptr;
  }
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
  storesize_ = size_;
  storage_ = new double [storesize_];
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
    if ((size_ + inc) > storesize_) {
      increase_store((size_ + inc));
    }
  }
  size_ += inc;
}

void History::increase_store(size_t newsize)
{
  storesize_ = newsize;
  double * newstore = new double[newsize];
  std::copy(storage_, storage_+size_, newstore);
  delete [] storage_;
  storage_ = newstore;
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
  if (contains(name)) {
    std::stringstream ss;
    ss << "History variable name " << name << " already stored." << std::endl;
    throw std::runtime_error(ss.str());
  }
}

void History::error_if_not_exists_(std::string name) const
{
  // This is a huge time drain
  if (not contains(name)) {
    std::stringstream ss;
    ss << "No history variable named " << name << " is stored." << std::endl;
    throw std::runtime_error(ss.str());
  }
}

void History::error_if_wrong_type_(std::string name, StorageType type) const
{
  // This is a huge time drain
  if (type != type_.at(name)) {
    std::stringstream ss;
    ss << name << " is not of the type requested." << std::endl;
    throw std::runtime_error(ss.str());
  }
}

History & History::zero()
{
  std::fill(storage_, storage_+size_, 0.0);
  return *this;
}

History History::history_derivative(const History & other) const
{
  History deriv;
  deriv.increase_store(size() * other.size());

  for (auto i1 : order_) {
    StorageType i1_type = type_.at(i1);
    for (auto i2 : other.get_order()) {
      StorageType i2_type = other.get_type().at(i2);
      StorageType ntype = derivative_type.at(i1_type).at(i2_type);
      deriv.add(i1+"_"+i2, ntype, storage_size.at(ntype));
    }
  }

  deriv.zero();

  return deriv;
}

History History::split(std::vector<std::string> sep, bool after) const
{
  // Check to see if the groups are contiguous and get the offset
  size_t i;
  for (i = 0; i < sep.size(); i++) {
    if (sep[i] != order_[i]) {
      throw std::runtime_error("History items to separate out must be contiguous!");
    }
  }
  // Special case where we need to return a zero history
  if (((i == order_.size()) && after) || ((i == 0) && (! after))) {
    if (store_) {
      return History(true);
    }
    else {
      return History(false);
    }
  }

  History res(false);
  
  // Move the maps
  if (after) {
    for (size_t j = i; j < order_.size(); j++) {
      res.add(order_[j], type_.at(order_[j]), storage_size.at(type_.at(order_[j])));
    }
  }
  else {
    for (size_t j = 0; j < i; j++) {
      res.add(order_[j], type_.at(order_[j]), storage_size.at(type_.at(order_[j])));
    }
  }

  // Either copy or just split, depending on if we own data
  if (store_) {
    res.make_store();
    if (after) {
      res.copy_data(&storage_[loc_.at(order_[i])]);
    }
    else {
      res.copy_data(&storage_[loc_.at(order_[0])]);
    }
  }
  else {
    if (after) {
      res.set_data(&storage_[loc_.at(order_[i])]);
    }
    else {
      res.set_data(&storage_[loc_.at(order_[0])]);
    }
  }

  return res;
}

History History::subset(std::vector<std::string> vars) const
{
  History nhist;
  for (auto var : vars) {
    error_if_not_exists_(var);
    size_t sz = storage_size.at(type_.at(var));
    nhist.add(var, type_.at(var), sz);
    std::copy(&storage_[loc_.at(var)], &storage_[loc_.at(var)]+sz, 
              &nhist.rawptr()[nhist.get_loc().at(var)]);

  }
  return nhist;
}

History & History::reorder(std::vector<std::string> vars)
{
  History nhist = subset(vars);

  copy_maps(nhist);

  std::copy(nhist.rawptr(), nhist.rawptr() + size(), storage_);

  return *this;
}

History History::postmultiply(const SymSymR4 & T)
{
  // This is not a guarantee that this object has the right form,
  // but is better than nothing
  if (size_ % 6 != 0)
    throw std::runtime_error("History object does not appear to be suitable "
                             "for postmultiplication by a SymSymR4!");

  History nhist = this->deepcopy();
  double * data = new double [size_];
  std::copy(storage_, storage_+size_, data);
  
  mat_mat(size_ / 6, 6, 6, data, T.data(), nhist.rawptr());

  delete [] data;

  return nhist;
}

size_t History::size_of_entry(std::string name) const
{
  return storage_size.at(type_.at(name));
}

void History::unravel_hh(const History & base, double * const array)
{
  // This assumes storage was correctly allocated

  // Row/col size
  size_t m = base.size();

  for (auto n1 : base.get_order()) {
    size_t s1 = base.size_of_entry(n1);
    size_t o1 = base.get_loc().at(n1);
    for (auto n2: base.get_order()) {
      size_t s2 = base.size_of_entry(n2);
      size_t o2 = base.get_loc().at(n2);
      size_t k = 0;
      size_t start = loc_.at(n1 + "_" + n2);
      for (size_t i = 0; i < s1; i++) {
        for (size_t j = 0; j < s2; j++) {
          array[CINDEX(o1+i,o2+j,m)] = storage_[start + k];
          k++;
        }
      }
    }
  }
}

double * History::start_loc(std::string name)
{
  error_if_not_exists_(name);
  return &(storage_[loc_.at(name)]);
}

std::vector<std::string> History::formatted_names() const
{
  std::vector<std::string> names;
  
  for (auto name : order_)
    if (size_of_entry(name) == 1) {
      names.push_back(name);
    }
    else {
      for (size_t i = 0; i < size_of_entry(name); i++)
        names.push_back(name + "_" + std::to_string(i));
    }

  return names;
}

HistoryNEMLObject::HistoryNEMLObject(ParameterSet & params) :
    NEMLObject(params), prefix_(""), cached_(false)
{
}

size_t HistoryNEMLObject::nhist() const
{
  if (cached_)
    return ncache_;

  History h;
  populate_hist(h);
  return h.size();
}

size_t HistoryNEMLObject::nh() const
{
  return ncache_;
}

void HistoryNEMLObject::init_store(double * const h) const
{
  History hist(h);
  populate_hist(hist);
  init_hist(hist);
}

size_t HistoryNEMLObject::nstore() const
{
  return nhist();
}

void HistoryNEMLObject::set_variable_prefix(std::string prefix)
{
  prefix_ = prefix;
}

std::string HistoryNEMLObject::get_variable_prefix() const
{
  return prefix_;
}

std::string HistoryNEMLObject::prefix(std::string basename) const
{
  return prefix_ + basename;
}

std::string HistoryNEMLObject::dprefix(std::string a, std::string b) const
{
  return prefix(a) + "_" + prefix(b);
}

void HistoryNEMLObject::cache_history_()
{
  populate_hist(stored_hist_);
  ncache_ = nhist();
  cached_ = true;
}

History HistoryNEMLObject::gather_history_(double * data) const
{
  History h = gather_blank_history_();
  h.set_data(data);
  return h;
}

History HistoryNEMLObject::gather_history_(const double * data) const
{
  History h = gather_blank_history_();
  h.set_data(const_cast<double*>(data));
  return h;
}

History HistoryNEMLObject::gather_blank_history_() const
{
  return stored_hist_;
}

} // namespace neml
