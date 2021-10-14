#ifndef NEML_INTERFACE_H
#define NEML_INTERFACE_H

#include "nemlerror.h"
#include "models.h"

namespace neml {
std::unique_ptr<NEMLModel> parse_xml_unique(std::string fname, std::string mname);
std::unique_ptr<NEMLModel> parse_string_unique(std::string input, std::string mname);
}

#endif // NEML_INTERFACE_H
