#include "parse.h"

namespace neml {

std::shared_ptr<NEMLModel> parse_xml(std::string fname, std::string mname)
{
  // Parse the XML file
  xmlpp::DomParser parser;
  parser.parse_file(fname);
  
  // Grab the root node
  const auto root = parser.get_document()->get_root_node();
  
  // Find the node with the right name
  xmlpp::Node * found = get_child(root, mname);

  // Get the NEMLObject
  std::shared_ptr<NEMLObject> obj = get_object(found);

  // Do a dangerous cast
  auto res = std::dynamic_pointer_cast<NEMLModel>(obj);
  if (res == nullptr) {
    throw std::runtime_error("Top level objects must be of type NEMLModel!");
  }
  else {
    return res;
  }
}

std::unique_ptr<NEMLModel> parse_xml_unique(std::string fname, std::string mname)
{
  // Parse the XML file
  xmlpp::DomParser parser;
  parser.parse_file(fname);
  
  // Grab the root node
  const auto root = parser.get_document()->get_root_node();
  
  // Find the node with the right name
  xmlpp::Node * found = get_child(root, mname);

  // Get the NEMLObject
  std::unique_ptr<NEMLObject> obj = get_object_unique(found);

  // Do a dangerous cast
  auto res = std::unique_ptr<NEMLModel>(dynamic_cast<NEMLModel*>(obj.release()));
  if (res == nullptr) {
    throw std::runtime_error("Top level objects must be of type NEMLModel!");
  }
  else {
    return res;
  }
}

std::unique_ptr<NEMLObject> get_object_unique(const xmlpp::Node * node)
{
  // Special case: could be a ConstantInterpolate
  std::string type = get_type_of_node(node);
  if (type == "none") {
    return make_unique<ConstantInterpolate>(get_double(node));
  }
  else {
    ParameterSet params = get_parameters(node);
    return Factory::Creator()->create_unique(params);
  }
}

std::shared_ptr<NEMLObject> get_object(const xmlpp::Node * node)
{
  // Special case: could be a ConstantInterpolate
  std::string type = get_type_of_node(node);
  if (type == "none") {
    return std::make_shared<ConstantInterpolate>(get_double(node));
  }
  else {
    ParameterSet params = get_parameters(node);
    return Factory::Creator()->create(params);
  }
}

ParameterSet get_parameters(const xmlpp::Node * node)
{
  std::string type = get_type_of_node(node);

  // Needs to have a type at this point
  if (type == "none") {
    throw InvalidType(node->get_name().raw(), type, node->get_line());
  }
  
  ParameterSet pset = Factory::Creator()->provide_parameters(type);

  auto children = node->get_children();
  for (auto it = children.begin(); it != children.end(); ++it) {
    std::string name = (*it)->get_name().raw();
    
    if (name == "text") continue;

    if (not pset.is_parameter(name)) {
      throw UnknownParameterXML(node->get_name().raw(), name, (*it)->get_line());
    }
    // The master switch block
    switch (pset.get_object_type(name)) {
      case TYPE_DOUBLE:
        pset.assign_parameter(name, get_double(*it));
        break;
      case TYPE_INT:
        pset.assign_parameter(name, get_int(*it));
        break;
      case TYPE_BOOL:
        pset.assign_parameter(name, get_bool(*it));
        break;
      case TYPE_VEC_DOUBLE:
        pset.assign_parameter(name, get_vector_double(*it));
        break;
      case TYPE_NEML_OBJECT:
        pset.assign_parameter(name, get_object(*it));
        break;
      case TYPE_VEC_NEML_OBJECT:
        pset.assign_parameter(name, get_vector_object(*it));
        break;
      case TYPE_STRING:
        pset.assign_parameter(name, get_string(*it));
        break;
      default:
        throw std::runtime_error("Unrecognized object type!");
        break;
    }
  }

  return pset;
}

std::vector<std::shared_ptr<NEMLObject>> get_vector_object(
    const xmlpp::Node * node)
{
  std::vector<std::shared_ptr<NEMLObject>> joined;
  
  auto children = node->get_children();
  for (auto it = children.begin(); it != children.end(); ++it) {
    std::string name = (*it)->get_name().raw();
    if (name == "text") continue;
    joined.push_back(get_object(*it));
  }

  return joined;
}

double get_double(const xmlpp::Node* node)
{
  std::string text = get_string(node);
  return std::stod(text);
}

int get_int(const xmlpp::Node* node)
{
  std::string text = get_string(node);
  return std::stoi(text);
}

std::vector<double> get_vector_double(const xmlpp::Node * node)
{
  std::string text = get_string(node);

  return split_string(text);
}

bool get_bool(const xmlpp::Node * node)
{
  std::string text = get_string(node);

  if ((text == "true") || (text == "True") || (text == "T") || (text == "1")) {
    return true;
  }
  else if ((text == "false") || (text == "False") || (text == "F") || (text == "0")) {
    return false;
  }
  else {
    throw InvalidParameter(node->get_name().raw(), node->get_line());
  }  
}

std::string get_string(const xmlpp::Node * node)
{
  auto elem = dynamic_cast<const xmlpp::Element*>(node);

  if (elem->has_child_text()) {
    std::string sval = elem->FIRST_TEXT_FN()->get_content();
    if (sval == "") {
      throw InvalidParameter(node->get_name().raw(), node->get_line());
    }
    else {
      return sval;
    }
  }
  else {
    throw InvalidParameter(node->get_name().raw(), node->get_line());
  }
}

xmlpp::Node * get_child(const xmlpp::Node * node, std::string name)
{
  auto matches = node->get_children(name);

  if (matches.size() > 1) {
    throw DuplicateNode(name, node->get_line());
  }
  else if (matches.size() < 1) {
    throw NodeNotFound(name, node->get_line());
  }
  else {
    return matches.front();
  }
}

std::string get_type_of_node(const xmlpp::Node * node)
{
  auto element = dynamic_cast<const xmlpp::Element*>(node);
  auto attributes = element->get_attributes();
  for (auto it = attributes.begin(); it != attributes.end(); ++it)
  {
    auto attr = *it;
    if (attr->get_name() == "type") {
      return attr->get_value();
    }
  }

  return "none";
}

std::vector<double> split_string(std::string sval)
{
  std::vector<std::string> splits;
  std::stringstream ss(sval);
  std::string temp;
  while (ss >> temp) {
    splits.push_back(temp);
  }
  std::vector<double> value;
  for (auto it = splits.begin(); it != splits.end(); ++it) {
    value.push_back(std::stod(*it));
  }
  return value;
}

} // namespace neml
