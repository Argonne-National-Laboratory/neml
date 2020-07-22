#include "objects.h"

// #include <fenv.h>

namespace neml {

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

Factory::Factory()
{
  // Good spot to place global things
  //feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO);
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
