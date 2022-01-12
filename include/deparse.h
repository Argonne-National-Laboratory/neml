#ifndef DEPARSE_H
#define DEPARSE_H

#include "objects.h"
#include "windows.h"

#include <sstream>

#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"

using namespace rapidxml;

namespace neml {

/// Convert from a parameter set back to an xml object
std::unique_ptr<xml_document<>> deparse(ParameterSet & input,
                                                  std::string top_name =
                                                  "object",
                                                  std::string top_node = "");

/// Convert from a parameter set to an ASCII XML string
std::string deparse_to_string(ParameterSet & input, 
                              std::string top_name = "object",
                              std::string top_node = "");

/// Make an object node from an entry in the parameter tree
xml_node<>* make_object_node(ParameterSet & input, std::string name, 
                             xml_document <> & doc); 

/// Make a node from an entry in the parameter tree using to_string
template <typename T>
xml_node<>* make_simple_node(const T & input, std::string name,
                             xml_document <> & doc)
{
  const char * pname = doc.allocate_string(name.c_str());
  const char * data = doc.allocate_string(std::to_string(input).c_str());
  return doc.allocate_node(node_element, pname, data);
}

/// Make a node for a double
xml_node<>* make_double_node(const double & input, std::string name,
                             xml_document <> & doc);

/// Make a node for a bool
xml_node<>* make_bool_node(const bool & input, std::string name,
                             xml_document <> & doc);

/// Make a node for a string
xml_node<>* make_string_node(const std::string & input, std::string name,
                             xml_document <> & doc);

/// Make a node for a vector of doubles
xml_node<>* make_vec_double_node(const std::vector<double> & input,
                                 std::string name, xml_document <> & doc);

/// Make a node for a vector of some standard type using to_string
template <typename T>
xml_node<>* make_vec_simple_node(const std::vector<T> & input,
                                 std::string name, xml_document <> & doc)
{
  const char * pname = doc.allocate_string(name.c_str());
  std::ostringstream ss;
  for (auto entry : input)
    ss << std::to_string(entry) << " ";
  const char * data = doc.allocate_string(ss.str().c_str());
  return doc.allocate_node(node_element, pname, data);
}

/// Make a node for a list of slip systems
xml_node<>* make_slip_node(const list_systems & input,
                           std::string name, xml_document <> & doc);

/// Make a node for a list of twin systems
xml_node<>* make_twin_node(const twin_systems & input,
                           std::string name, xml_document <> & doc);

/// Make a node for a vector of objects
xml_node<>* make_vec_object_node(const std::vector<std::shared_ptr<NEMLObject>> & input,
                                 std::string name, xml_document <> & doc);

} // namespace neml

#endif
