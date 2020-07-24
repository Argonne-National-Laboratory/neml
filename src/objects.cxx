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

const std::string & ParameterSet::type() const
{
  return type_;
}

void ParameterSet::assign_defered_parameter(std::string name, const ParameterSet & value)
{
  defered_params_[name] = value;
}

bool ParameterSet::is_parameter(std::string name) const
{
  return params_.find(name) != params_.end();
}

std::vector<std::string> ParameterSet::unassigned_parameters()
{
  resolve_objects_();

  std::vector<std::string> uparams;

  for (const auto & p : params_)
    if (!p.second->isSet())
      uparams.push_back(p.first);

  return uparams;
}

bool ParameterSet::fully_assigned()
{
  resolve_objects_();

  for (const auto & p : params_)
    if (!p.second->isSet())
      return false;

  return true;
}

void ParameterSet::resolve_objects_()
{
  for (auto & p : defered_params_)
    add_optional_parameter<NEMLObject>(p.first, Factory::Creator()->create(p.second));
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

// as a shorthand std::shared_ptr<NEMLObject> parameters are declared as NEMLObject parameters
template<>
void ParameterSet::add_parameter<NEMLObject>(std::string name)
{
  add_parameter<std::shared_ptr<NEMLObject>>(name);
}

// as a shorthand std::vector<std::shared_ptr<NEMLObject>> parameters are declared as std::vector<NEMLObject> parameters
template<>
void ParameterSet::add_parameter<std::vector<NEMLObject>>(std::string name)
{
  add_parameter<std::vector<std::shared_ptr<NEMLObject>>>(name);
}

std::shared_ptr<NEMLObject> Factory::create(ParameterSet & params)
{
  if (not params.fully_assigned())
    throw UndefinedParameters(params.type(), params.unassigned_parameters());

  try {
    return creators_[params.type()](params);
  }
  catch (std::out_of_range & e) {
    throw UnregisteredError(params.type());
  }
}

std::unique_ptr<NEMLObject> Factory::create_unique(ParameterSet & params)
{
  if (not params.fully_assigned())
    throw UndefinedParameters(params.type(), params.unassigned_parameters());

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
