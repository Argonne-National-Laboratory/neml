#include "objects.h"

namespace neml {

template <> double & param_type::data<double>() { return double_; }
template <> int & param_type::data<int>() { return int_; }
template <> bool & param_type::data<bool>() { return bool_; }
template <> std::vector<double> & param_type::data<std::vector<double>>() { return vec_double_; }
template <> std::shared_ptr<NEMLObject> & param_type::data<std::shared_ptr<NEMLObject>>() { return neml_object_; }
template <> std::vector<std::shared_ptr<NEMLObject>> & param_type::data<std::vector<std::shared_ptr<NEMLObject>>>() { return vec_neml_object_; }
template <> std::string & param_type::data<std::string>() { return string_; }
template <> list_systems & param_type::data<list_systems>() { return list_systems_; }
template <> std::size_t & param_type::data<std::size_t>() { return size_t_; }
template <> std::vector<std::size_t> & param_type::data<std::vector<std::size_t>>() { return vec_size_t_; }

param_type::param_type(const char * val) { data<std::string>() = *val; }

ParameterSet::ParameterSet() :
    type_("invalid")
{

}

ParameterSet::ParameterSet(std::string type) :
    type_(type)
{

}

ParameterSet::~ParameterSet()
{

}

const std::string & ParameterSet::type() const
{
  return type_;
}

void ParameterSet::assign_defered_parameter(std::string name, ParameterSet value)
{
  defered_params_[name] = value;
}

ParamType ParameterSet::get_object_type(std::string name)
{
  return param_types_[name];
}

bool ParameterSet::is_parameter(std::string name) const
{
  return std::find(param_names_.begin(), param_names_.end(), name) !=
      param_names_.end();
}

std::vector<std::string> ParameterSet::unassigned_parameters()
{
  resolve_objects_();

  std::vector<std::string> uparams;

  for (auto it = param_names_.begin(); it != param_names_.end(); ++it) {
    if (params_.find(*it) == params_.end()) {
      uparams.push_back(*it);
    }
  }

  return uparams;
}

bool ParameterSet::fully_assigned()
{
  resolve_objects_();

  for (auto it = param_names_.begin(); it != param_names_.end(); ++it) {
    if (params_.find(*it) == params_.end()) return false;
  }

  return true;
}

void ParameterSet::resolve_objects_()
{
  for (auto it = defered_params_.begin(); it != defered_params_.end(); ++it) {
    params_[it->first] = Factory::Creator()->create(it->second);
  }
  defered_params_.clear();
}

ParameterSet Factory::provide_parameters(std::string type)
{
  try {
    return setups_[type]();
  }
  catch (std::exception & e) {
    throw UnregisteredError(type);
  }
}

std::shared_ptr<NEMLObject> Factory::create(ParameterSet & params)
{
  if (not params.fully_assigned()) {
    throw UndefinedParameters(params.type(), params.unassigned_parameters());
  }

  try {
    return creators_[params.type()](params);
  }
  catch (std::out_of_range & e) {
    throw UnregisteredError(params.type());
  }
}

std::unique_ptr<NEMLObject> Factory::create_unique(ParameterSet & params)
{
  if (not params.fully_assigned()) {
    throw UndefinedParameters(params.type(), params.unassigned_parameters());
  }

  try {
    return creators_[params.type()](params);
  }
  catch (std::out_of_range & e) {
      throw UnregisteredError(params.type());
  }
}

void Factory::register_type(std::string type,
                            std::function<std::unique_ptr<NEMLObject>(ParameterSet &)> creator,
                            std::function<ParameterSet()> setup)
{
  creators_[type] = creator;
  setups_[type] = setup;
}

Factory * Factory::Creator()
{
  static Factory creator;
  return &creator;
}

} // namespace neml
