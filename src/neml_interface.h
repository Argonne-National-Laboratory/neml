#ifndef NEML_INTERFACE_H
#define NEML_INTERFACE_H

#include "nemlerror.h"
#include "models.h"
#include "math/nemlmath.h"
#include "cp/singlecrystal.h"

namespace neml {
std::unique_ptr<NEMLModel> parse_xml_unique(std::string fname, std::string mname);
std::unique_ptr<NEMLModel> parse_string_unique(std::string input);
}

#endif // NEML_INTERFACE_H
