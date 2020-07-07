#ifndef PARSE_H
#define PARSE_H

#include "objects.h"
#include "models.h"
#include "damage.h"

#include "windows.h"

#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"

#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <exception>

namespace neml {

/// Parse from a string to a shared_ptr
std::shared_ptr<NEMLModel> parse_string(std::string input);

/// Parse from a string to a unique_ptr
std::unique_ptr<NEMLModel> parse_string_unique(std::string input, std::string mname);

/// Parse from file to a shared_ptr
NEML_EXPORT std::shared_ptr<NEMLModel> parse_xml(std::string fname, std::string mname);

/// Parse from file to a unique_ptr
NEML_EXPORT std::unique_ptr<NEMLModel> parse_xml_unique(std::string fname, std::string mname);

/// Extract a NEMLObject from a xml node as a unique_ptr
NEML_EXPORT std::unique_ptr<NEMLObject> get_object_unique(const rapidxml::xml_node<> * node);

/// Extract a NEMLObject from a xml node
std::shared_ptr<NEMLObject> get_object(const rapidxml::xml_node<> * node);

/// Actually get a valid parameter set from a node
 ParameterSet get_parameters(const rapidxml::xml_node<> * node);

/// Extract a vector of NEMLObjects from an xml node
std::vector<std::shared_ptr<NEMLObject>> get_vector_object(const rapidxml::xml_node<> * node);

/// Extract a double from an xml node
double get_double(const rapidxml::xml_node<> * node);

/// Extract an integer parameter
int get_int(const rapidxml::xml_node<> * node);

/// Extract a vector of doubles from an xml node
std::vector<double> get_vector_double(const rapidxml::xml_node<> * node);

/// Extract a bool from an xml node
bool get_bool(const rapidxml::xml_node<> * node);

/// Extract a string from an xml node
std::string get_string(const rapidxml::xml_node<> * node);

/// Extract a slip system from an xml node
list_systems get_slip(const rapidxml::xml_node<> * node);

/// Extract a size_type from an xml node
size_t get_size_type(const rapidxml::xml_node<> * node);

/// Extract a vector of size types from an xml node
std::vector<size_t> get_vector_size_type(const rapidxml::xml_node<> * node);

// Helpers
/// Get a node with a given name
const rapidxml::xml_node<> * get_child(const rapidxml::xml_node<> * node, std::string name);

/// Return the type of a node
std::string get_type_of_node(const rapidxml::xml_node<> * node);

/// Helper to split strings
std::vector<double> split_string(std::string sval);

/// Helper to split lists of size_ts
std::vector<size_t> split_string_size_type(std::string sval);

/// Lame we can't do this with templates
std::vector<int> split_string_int(std::string sval);

/// Helper to strip strings
std::string & strip(std::string & s);

// Exceptions
/// If a node is not found
class NodeNotFound: public std::exception {
 public:
  NodeNotFound(std::string node_name, int line) :
      node_name_(node_name), line_(line)
  {

  };

  const char * what() const throw ()
  {
    std::stringstream ss;

    ss << "Node with name " << node_name_
        << " was not found near line " << line_ << "!";

    return ss.str().c_str();
  };

 private:
  std::string node_name_;
  int line_;
};

/// If a node is not unique (and it should be)
class DuplicateNode: public std::exception {
 public:
  DuplicateNode(std::string node_name, int line) :
      node_name_(node_name), line_(line)
  {

  };

    const char * what() const throw ()
    {
      std::stringstream ss;

      ss << "Multiple nodes with name " << node_name_ << " were found!";

      return ss.str().c_str();
    };

  private:
    std::string node_name_;
    int line_;
};

/// If the object can't be converted
class InvalidType: public std::exception {
 public:
  InvalidType(std::string name, std::string type, std::string ctype) :
      name_(name), type_(type), ctype_(ctype)
  {

  };

  const char * what() const throw ()
  {
    std::stringstream ss;

    ss << "Node with name " << name_ << " and type " << type_
        << "cannot be converted to the correct type " << ctype_ << "!";

    return ss.str().c_str();
  };

 private:
  const std::string name_, type_, ctype_;
};

/// If a parameter doesn't exist
class UnknownParameterXML: public std::exception {
 public:
  UnknownParameterXML(std::string name, std::string param) :
      name_(name), param_(param)
  {

  };

  const char * what() const throw ()
  {
    std::stringstream ss;

    ss << "Object " << name_ << " does not have a parameter called " << param_ << "!";

    return ss.str().c_str();
  };

 private:
  const std::string name_, param_;

};

/// The object isn't in the factory
class UnregisteredXML: public std::exception {
 public:
  UnregisteredXML(std::string name, std::string type) :
      name_(name), type_(type)
  {

  };

  const char * what() const throw ()
  {
    std::stringstream ss;

    ss << "Node named " << name_ << " has an unregistered type of " << type_ << "!";

    return ss.str().c_str();
  };

 private:
  const std::string name_, type_;

};

} // namespace neml

#endif // PARSE_H
