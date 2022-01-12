#include "deparse.h"

#include "rapidxml_print.hpp"
#include <iostream>

namespace neml {

std::unique_ptr<xml_document<>> deparse(ParameterSet & input,
                                        std::string top_name,
                                        std::string top_node)
{
  auto doc = neml::make_unique<xml_document<>>();
  
  xml_node<> * top;
  if (top_node == "") 
    top = doc.get();
  else {
    const char * name = doc->allocate_string(top_node.c_str());
    top = doc->allocate_node(node_element, name);
    doc->append_node(top);
  }

  top->append_node(make_object_node(input, top_name, *doc));

  return doc;
}

std::string deparse_to_string(ParameterSet & input, std::string top_name,
                              std::string top_node)
{
  auto doc = deparse(input, top_name, top_node);

  std::string res;
  print(std::back_inserter(res), *doc, 0);

  return res;
}

xml_node<>* make_object_node(ParameterSet & input, std::string name, 
                             xml_document <> & doc)
{
  const char * pname = doc.allocate_string(name.c_str());
  xml_node<>* node = doc.allocate_node(node_element, pname);

  const char* ptype = doc.allocate_string(input.type().c_str());
  xml_attribute<>* attr = doc.allocate_attribute("type", ptype);
  node->append_attribute(attr);

  for (auto name: input.param_names()) {
    switch (input.get_object_type(name)) {
      case TYPE_NEML_OBJECT:
        node->append_node(make_object_node(
                input.get_parameter<std::shared_ptr<NEMLObject>>(name)->current_parameters(),
                name, doc));
        break;
      case TYPE_INT:
        node->append_node(make_simple_node<int>(
                input.get_parameter<int>(name),
                name, doc));
        break;
      case TYPE_SIZE_TYPE:
        node->append_node(make_simple_node<size_t>(
                input.get_parameter<size_t>(name),
                name, doc));
        break;
      case TYPE_DOUBLE:
        node->append_node(make_double_node(input.get_parameter<double>(name),
                                           name, doc));
        break;
      case TYPE_BOOL:
        node->append_node(make_bool_node(input.get_parameter<bool>(name),
                                         name, doc));
        break;
      case TYPE_STRING:
        node->append_node(make_string_node(input.get_parameter<std::string>(name),
                                           name, doc));
        break;
      case TYPE_VEC_DOUBLE:
        node->append_node(make_vec_double_node(input.get_parameter<std::vector<double>>(name),
                                               name, doc));
        break;
      case TYPE_VEC_SIZE_TYPE:
        node->append_node(make_vec_simple_node<size_t>(input.get_parameter<std::vector<size_t>>(name),
                                                       name, doc));
        break;
      case TYPE_SLIP:
        node->append_node(make_slip_node(input.get_parameter<list_systems>(name),
                                         name, doc));
        break;
      case TYPE_TWIN:
        node->append_node(make_twin_node(input.get_parameter<twin_systems>(name),
                                         name, doc));
        break;
      case TYPE_VEC_NEML_OBJECT:
        node->append_node(make_vec_object_node(input.get_parameter<std::vector<std::shared_ptr<NEMLObject>>>(name),
                                               name, doc));
        break;
      default:
        throw std::runtime_error("Internal error: unknown object type in "
                                 "serialization!");
    }
  }

  return node;
}

xml_node<>* make_double_node(const double & input, std::string name,
                             xml_document <> & doc)
{
  const char * pname = doc.allocate_string(name.c_str());

  std::ostringstream ss;
  ss << input;

  const char * data = doc.allocate_string(ss.str().c_str());
  return doc.allocate_node(node_element, pname, data);
}

xml_node<>* make_bool_node(const bool & input, std::string name,
                             xml_document <> & doc)
{
  const char * pname = doc.allocate_string(name.c_str());

  std::string res;
  if (input) 
    res = "true";
  else
    res = "false";

  const char * data = doc.allocate_string(res.c_str());
  return doc.allocate_node(node_element, pname, data);
}

xml_node<>* make_string_node(const std::string & input, std::string name,
                             xml_document <> & doc)
{
  const char * pname = doc.allocate_string(name.c_str());
  const char * data = doc.allocate_string(input.c_str());
  return doc.allocate_node(node_element, pname, data);
}

xml_node<>* make_vec_double_node(const std::vector<double> & input,
                                 std::string name, xml_document <> & doc)
{
  const char * pname = doc.allocate_string(name.c_str());
  std::ostringstream ss;
  for (auto entry : input)
    ss << entry << " ";
  const char * data = doc.allocate_string(ss.str().c_str());
  return doc.allocate_node(node_element, pname, data);
}

xml_node<>* make_slip_node(const list_systems & input,
                           std::string name, xml_document <> & doc)
{
  const char * pname = doc.allocate_string(name.c_str());
  std::ostringstream ss;
  for (auto entry : input) {
    for (auto ni : entry.first)
      ss << ni << " ";
    ss << "; ";
    for (auto ni : entry.second)
      ss << ni << " ";
    ss << ",";
  }

  const char * data = doc.allocate_string(ss.str().c_str());
  return doc.allocate_node(node_element, pname, data);
}

xml_node<>* make_twin_node(const twin_systems & input,
                           std::string name, xml_document <> & doc)
{
  const char * pname = doc.allocate_string(name.c_str());
  std::ostringstream ss;
  for (auto entry : input) {
    for (auto ni : std::get<0>(entry))
      ss << ni << " ";
    ss << ";";
    for (auto ni : std::get<1>(entry))
      ss << ni << " ";
    ss << ";";
    for (auto ni : std::get<2>(entry))
      ss << ni << " ";
    ss << ";";
    for (auto ni : std::get<3>(entry))
      ss << ni << " ";
    ss << ",";
  }

  const char * data = doc.allocate_string(ss.str().c_str());
  return doc.allocate_node(node_element, pname, data);
}

xml_node<>* make_vec_object_node(const std::vector<std::shared_ptr<NEMLObject>> & input,
                                 std::string name, xml_document <> & doc)
{
  const char * pname = doc.allocate_string(name.c_str());
  xml_node<> * vec_node = doc.allocate_node(node_element, pname);
  
  size_t i = 0;
  for (auto obj : input) {
    std::string sub_name = name + std::to_string(i);
    vec_node->append_node(make_object_node(obj->current_parameters(), sub_name, doc));
    i++;
  }

  return vec_node;
}

} // namespace neml
