#ifndef PARSE_H
#define PARSE_H

#include "objects.h"
#include "models.h"
#include "damage.h"

#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <exception>

#include <libxml++/libxml++.h>

#ifdef LIBXMLppV3
#define FIRST_TEXT_FN get_first_child_text
#else
#define FIRST_TEXT_FN get_child_text
#endif

namespace neml {

/// Parse from file to a shared_ptr
std::shared_ptr<NEMLModel> parse_xml(std::string fname, std::string mname);

/// Parse from file to a unique_ptr
std::unique_ptr<NEMLModel> parse_xml_unique(std::string fname, std::string mname);

/// Extract a NEMLObject from a xml node as a unique_ptr
std::unique_ptr<NEMLObject> get_object_unique(const xmlpp::Node * node);

/// Extract a NEMLObject from a xml node
std::shared_ptr<NEMLObject> get_object(const xmlpp::Node * node);

/// Actually get a valid parameter set from a node
ParameterSet get_parameters(const xmlpp::Node * node);

/// Extract a vector of NEMLObjects from an xml node
std::vector<std::shared_ptr<NEMLObject>> get_vector_object(const xmlpp::Node * node);

/// Extract a double from an xml node
double get_double(const xmlpp::Node* node);

/// Extract an integer parameter
int get_int(const xmlpp::Node* node);

/// Extract a vector of doubles from an xml node
std::vector<double> get_vector_double(const xmlpp::Node * node);

/// Extract a bool from an xml node
bool get_bool(const xmlpp::Node * node);

/// Extract a string from an xml node
std::string get_string(const xmlpp::Node * node);

// Helpers
/// Get a node with a given name
xmlpp::Node * get_child(const xmlpp::Node * node, std::string name);

/// Return the type of a node
std::string get_type_of_node(const xmlpp::Node *);

/// Helper to split strings
std::vector<double> split_string(std::string sval);

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
  InvalidType(std::string name, std::string type, int line) :
      name_(name), type_(type), line_(line)
  {

  };

  const char * what() const throw ()
  {
    std::stringstream ss;

    ss << "Node with name " << name_ << " and type " << type_ 
        << " near line " << line_ 
        << "cannot be converted to the appropriate type!";

    return ss.str().c_str();
  };

 private:
  std::string name_, type_;
  int line_;
};

/// If a parameter doesn't exist
class UnknownParameterXML: public std::exception {
 public:
  UnknownParameterXML(std::string name, std::string param, int line) :
      name_(name), param_(param), line_(line)
  {

  };

  const char * what() const throw ()
  {
    std::stringstream ss;

    ss << "Object " << name_ << " defined near line " << line_ << 
        " does not have a parameter called " << param_ << "!";

    return ss.str().c_str();
  };

 private:
  std::string name_, param_;
  int line_;
};

/// Nonsensical or invalid parameter
class InvalidParameter: public std::exception {
 public:
  InvalidParameter(std::string name, int line) :
      name_(name), line_(line)
  {

  };

  const char * what() const throw()
  {
    std::stringstream ss;

    ss << "Parameter " << name_ << " near line " << line_ 
        << " is invalid!";

    return ss.str().c_str();
  };

 private:
  std::string name_;
  int line_;
};

} // namespace neml

#endif // PARSE_H
