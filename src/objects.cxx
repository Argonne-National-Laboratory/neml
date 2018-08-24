#include "objects.h"

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
  return setups_[type]();
}

std::shared_ptr<NEMLObject> Factory::create(ParameterSet & params)
{
  if (not params.fully_assigned()) {
    throw std::runtime_error("Parameter set not fully assigned!");
  }

  return creators_[params.type()](params);
}

void Factory::register_type(std::string type,
                            std::function<std::shared_ptr<NEMLObject>(ParameterSet &)> creator,
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
