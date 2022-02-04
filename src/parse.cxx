#include "parse.h"

namespace neml {

std::shared_ptr<NEMLModel> parse_string(std::string input)
{
  // Parse the string to the rapidxml representation
  rapidxml::xml_document<> doc;
  doc.parse<0>(&input[0]);

  // The model is the root node
  const rapidxml::xml_node<> * found = doc.first_node();

  // Get the NEMLObject
  std::shared_ptr<NEMLObject> obj = get_object(found);

  // Do a dangerous cast
  auto res = std::dynamic_pointer_cast<NEMLModel>(obj);
  if (res == nullptr) {
    throw InvalidType(found->name(), get_type_of_node(found), "NEMLModel");
  }
  else {
    return res;
  }
}

std::shared_ptr<NEMLObject> get_object_string(std::string repr)
{
  // Parse the string to the rapidxml representation
  rapidxml::xml_document<> doc;
  doc.parse<0>(&repr[0]);

  // The object, regardless of name
  const rapidxml::xml_node<> * found = doc.first_node();

  return get_object(found);
}

std::unique_ptr<NEMLModel> parse_string_unique(std::string input, std::string mname)
{
  // Parse the string to the rapidxml representation
  rapidxml::xml_document<> doc;
  doc.parse<0>(&input[0]);

  // Grab the root node
  const rapidxml::xml_node<> * root = doc.first_node();

  // Find the node with the right name
  const rapidxml::xml_node<> * found = root->first_node(mname.c_str());
  if (found == 0) throw ModelNotFound(mname);

  // Get the NEMLObject
  std::unique_ptr<NEMLObject> obj = get_object_unique(found);

  // Do a dangerous cast
  auto res = std::unique_ptr<NEMLModel>(dynamic_cast<NEMLModel *>(obj.release()));
  if (res == nullptr) {
    throw InvalidType(found->name(), get_type_of_node(found), "NEMLModel");
  }
  else {
    return res;
  }
}

std::shared_ptr<NEMLModel> parse_xml(std::string fname, std::string mname)
{
  // Parse the XML file
  rapidxml::file <> xmlFile(fname.c_str());
  rapidxml::xml_document<> doc;
  doc.parse<0>(xmlFile.data());

  // Grab the root node
  const rapidxml::xml_node<> * root = doc.first_node();

  // Find the node with the right name
  const rapidxml::xml_node<> * found = root->first_node(mname.c_str());
  if (found == 0) throw ModelNotFound(mname);

  // Get the NEMLObject
  std::shared_ptr<NEMLObject> obj = get_object(found);

  // Do a dangerous cast
  auto res = std::dynamic_pointer_cast<NEMLModel>(obj);
  if (res == nullptr) {
    throw InvalidType(found->name(), get_type_of_node(found), "NEMLModel");
  }
  else {
    return res;
  }
}

std::unique_ptr<NEMLModel> parse_xml_unique(std::string fname, std::string mname)
{
  // Parse the XML file
  rapidxml::file <> xmlFile(fname.c_str());
  rapidxml::xml_document<> doc;
  doc.parse<0>(xmlFile.data());

  // Grab the root node
  const rapidxml::xml_node<> * root = doc.first_node();

  // Find the node with the right name
  const rapidxml::xml_node<> * found = root->first_node(mname.c_str());
  if (found == 0) throw ModelNotFound(mname);

  // Get the NEMLObject
  std::unique_ptr<NEMLObject> obj = get_object_unique(found);

  // Do a dangerous cast
  auto res = std::unique_ptr<NEMLModel>(dynamic_cast<NEMLModel*>(obj.release()));
  if (res == nullptr) {
    throw InvalidType(found->name(), get_type_of_node(found), "NEMLModel");
  }
  else {
    return res;
  }
}

std::unique_ptr<NEMLObject> get_object_unique(const rapidxml::xml_node<> * node) {
  // Special case: could be a ConstantInterpolate
  std::string type = get_type_of_node(node);
  if (type == "none") {
    return make_constant_unique(get_double(node));
  }
  else {
    ParameterSet params = get_parameters(node);
    try {
      return Factory::Creator()->create_unique(params);
    }
    catch (UnregisteredError & e) {
      throw UnregisteredXML(node->name(), type);
    }
  }
}

std::shared_ptr<NEMLObject> get_object(const rapidxml::xml_node<> * node)
{
  // Special case: could be a ConstantInterpolate
  std::string type = get_type_of_node(node);
  if (type == "none") {
    return make_constant(get_double(node));
  }
  else {
    ParameterSet params = get_parameters(node);
    return Factory::Creator()->create(params);
  }
}

ParameterSet get_parameters(const rapidxml::xml_node<> * node)
{
  std::string type = get_type_of_node(node);

  // Needs to have a type at this point
  if (type == "none") {
    throw InvalidType(node->name(), type, "NEMLObject");
  }

  ParameterSet pset;
  try {
    pset = Factory::Creator()->provide_parameters(type);
  }
  catch (UnregisteredError & e) {
    throw UnregisteredXML(node->name(), type);
  }

  for (rapidxml::xml_node<> * child = node->first_node(); child; child = child->next_sibling()) {
    std::string name = (child)->name();

    if (name == "text") continue;

    if (not pset.is_parameter(name)) {
      throw UnknownParameterXML(node->name(), name);
    }

    switch (pset.get_object_type(name)) {
      case TYPE_DOUBLE:
        pset.assign_parameter(name, get_double(child));
        break;
      case TYPE_INT:
        pset.assign_parameter(name, get_int(child));
        break;
      case TYPE_BOOL:
        pset.assign_parameter(name, get_bool(child));
        break;
      case TYPE_VEC_DOUBLE:
        pset.assign_parameter(name, get_vector_double(child));
        break;
      case TYPE_NEML_OBJECT:
        pset.assign_parameter(name, get_object(child));
        break;
      case TYPE_VEC_NEML_OBJECT:
        pset.assign_parameter(name, get_vector_object(child));
        break;
      case TYPE_STRING:
        pset.assign_parameter(name, get_string(child));
        break;
      case TYPE_SLIP:
        pset.assign_parameter(name, get_slip(child));
        break;
      case TYPE_TWIN:
        pset.assign_parameter(name, get_twin(child));
        break;
      case TYPE_SIZE_TYPE:
        pset.assign_parameter(name, get_size_type(child));
        break;
      case TYPE_VEC_SIZE_TYPE:
        pset.assign_parameter(name, get_vector_size_type(child));
        break;
      default:
        throw std::runtime_error("Unrecognized object type!");
        break;
    }
  }

  return pset;
}

std::vector<std::shared_ptr<NEMLObject>> get_vector_object(
    const rapidxml::xml_node<> * node)
{
  if (is_empty(node)) return std::vector<std::shared_ptr<NEMLObject>>();

  std::vector<std::shared_ptr<NEMLObject>> joined;

  // A somewhat dangerous shortcut -- a list of text values should be
  // interpreted as a list of ConstantInterpolates
  if ((rapidxml::count_children(const_cast<rapidxml::xml_node<>*>(node))==1) and 
      (node->first_node()->type() == rapidxml::node_data)) {
    std::vector<double> data = get_vector_double(node);
    for (auto v : data) {
      joined.push_back(make_constant(v));
    }
    return joined;
  }

  for (rapidxml::xml_node<> * child = node->first_node(); child; child = child->next_sibling()) {
    std::string name = (child)->name();
    if (name == "text") continue;
    joined.push_back(get_object(child));
  }

  return joined;
}

// Need to find a way to replace node->get_line() in the throw exception
double get_double(const rapidxml::xml_node<> * node)
{
  try {
    std::string text = get_string(node);
    return std::stod(text);
  }
  catch (std::exception & e) {
    throw InvalidType(node->name(), get_type_of_node(node), "double");
  }
}

int get_int(const rapidxml::xml_node<> * node)
{
  try {
    std::string text = get_string(node);
    return std::stoi(text);
  }
  catch (std::exception & e) {
    throw InvalidType(node->name(), get_type_of_node(node), "int");
  }
}

std::vector<double> get_vector_double(const rapidxml::xml_node<> * node)
{
  if (is_empty(node)) return std::vector<double>();

  try {
    std::string text = get_string(node);
    return split_string(text);
  }
  catch (std::exception & e) {
    throw InvalidType(node->name(), get_type_of_node(node), "vector<double>");
  }
}

bool get_bool(const rapidxml::xml_node<> * node)
{
  std::string text = get_string(node);

  if ((text == "true") || (text == "True") || (text == "T") || (text == "1")) {
    return true;
  }
  else if ((text == "false") || (text == "False") || (text == "F") || (text == "0")) {
    return false;
  }
  else {
    throw InvalidType(node->name(), get_type_of_node(node), "bool");
  }
}

std::string get_string(const rapidxml::xml_node<> * node)
{
  std::string sval = node->first_node()->value();

  if (sval != "")
    return sval;
  else
    throw InvalidType(node->name(), get_type_of_node(node), "string");
}

list_systems get_slip( const rapidxml::xml_node<> * node)
{
  list_systems groups;

  if (node->value_size() == 0)
    return groups;

  std::string text = get_string(node);

  std::stringstream ss(text);
  std::string to;

  // Separate by newlines
  while(std::getline(ss, to,',')) {
    // Delete blank characters at front and back of string
    strip(to);

    // Skip blanks
    if (to == "") continue;

    // Split by semicolon
    std::string dir = to.substr(0, to.find(";"));
    strip(dir);
    std::string nor = to.substr(to.find(";")+1);
    strip(nor);

    // Make into vectors
    auto d = split_string_int(dir);
    auto n = split_string_int(nor);

    groups.push_back(make_pair(d,n));
  }

  return groups;
}

twin_systems get_twin( const rapidxml::xml_node<> * node)
{
  twin_systems groups;

  if (node->value_size() == 0)
    return groups;

  std::string text = get_string(node);

  std::stringstream ss(text);
  std::string to;

  // Separate by newlines
  while(std::getline(ss, to,',')) {
    // Delete blank characters at front and back of string
    strip(to);

    // Skip blanks
    if (to == "") continue;

    // Split by semicolon
    size_t curr = 0;
    size_t next = to.find(";");
    std::string eta1s = to.substr(curr, next-curr);
    strip(eta1s);
    curr = next + 1;

    next = to.find(";", curr);
    std::string K1s = to.substr(curr, next-curr);
    strip(K1s);
    curr = next + 1;
    
    next = to.find(";", curr);
    std::string eta2s = to.substr(curr, next-curr);
    strip(eta2s);
    curr = next + 1;
    
    std::string K2s = to.substr(curr);
    strip(K2s);

    // Make into vectors
    auto eta1 = split_string_int(eta1s);
    auto K1 = split_string_int(K1s);
    auto eta2 = split_string_int(eta2s);
    auto K2 = split_string_int(K2s);

    groups.push_back(make_tuple(eta1,K1,eta2,K2));
  }

  return groups;
}

size_t get_size_type(const rapidxml::xml_node<> * node)
{
  try {
    std::string text = get_string(node);
    return size_t(std::stoul(text));
  }
  catch (std::exception & e) {
    throw InvalidType(node->name(), get_type_of_node(node), "size_t");
  }
}

std::vector<size_t> get_vector_size_type(const rapidxml::xml_node<> * node)
{
  if (is_empty(node)) return std::vector<size_t>();
  try {
    std::string text = get_string(node);
    return split_string_size_type(text);
  }
  catch (std::exception & e) {
    throw InvalidType(node->name(), get_type_of_node(node), "vector<size_t>");
  }
}

std::string get_type_of_node(const rapidxml::xml_node<> * node)
{
  for (auto attributes = node->first_attribute(); attributes; attributes = attributes->next_attribute())
  {

    std::string attr_name = attributes->name();
    if (attr_name == "type")
      return attributes->value();
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

std::vector<size_t> split_string_size_type(std::string sval)
{
  std::vector<std::string> splits;
  std::stringstream ss(sval);
  std::string temp;
  while (ss >> temp) {
    splits.push_back(temp);
  }
  std::vector<size_t> value;
  for (auto it = splits.begin(); it != splits.end(); ++it) {
    value.push_back(size_t(std::stoul(*it)));
  }
  return value;
}

std::vector<int> split_string_int(std::string sval)
{
  std::vector<std::string> splits;
  std::stringstream ss(sval);
  std::string temp;
  while (ss >> temp) {
    splits.push_back(temp);
  }
  std::vector<int> value;
  for (auto it = splits.begin(); it != splits.end(); ++it) {
    value.push_back(std::stoi(*it));
  }
  return value;
}

std::string & strip(std::string & s)
{
  auto noblank = [](char c) { return !std::isspace(c);}; 

  s.erase(s.begin(), std::find_if(s.begin(), s.end(), noblank));
  s.erase(std::find_if(s.rbegin(), s.rend(), noblank).base(), s.end());

  return s;
}

bool is_empty(const rapidxml::xml_node<> * node)
{
  if (node->first_node() == NULL) return true;
  return false;
}

} // namespace neml
