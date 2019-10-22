#include "nemlerror.h"

#include <iostream>

namespace neml {

void py_error(int ier)
{
  switch(ier) {
    case SUCCESS: return;
    case INCOMPATIBLE_MODELS: throw std::runtime_error("Incompatible submodels");
    case LINALG_FAILURE: throw LinalgError("Generic linear algebra failure");
    case MAX_ITERATIONS: throw std::runtime_error("Maximum iteration count exceeded");
    case KT_VIOLATION: throw std::runtime_error("Integration of rate-independent model resulted in a violation of the Kuhn-Tucker conditions");
    case NODE_NOT_FOUND: throw std::runtime_error("XML node not found");
    case TOO_MANY_NODES: throw std::runtime_error("More than one XML node found");
    case ATTRIBUTE_NOT_FOUND: throw std::runtime_error("XML attribute not found");
    case UNKNOWN_TYPE: throw std::runtime_error("Unknown model type");
    case BAD_TEXT: throw std::runtime_error("Bad text data in XML node");
    case INVALID_TYPE: throw std::runtime_error("Type described by XML node is invalid here");
    case FILE_NOT_FOUND: throw std::runtime_error("File not found");
    case CREEP_PLASTICITY: throw std::runtime_error("Creep models can only be combined with rate independent models");
    case INCOMPATIBLE_KM: throw std::runtime_error("Incompatible lengths in Kocks-Mecking region model: number of models = number of splits + 1");
    case DUMMY_ELASTIC: throw std::runtime_error("Calling for elastic constants from a dummy elastic model.");
    case INCOMPATIBLE_VECTORS: throw std::runtime_error("Inputs do not have the same length.");
    case UNKNOWN_ERROR: throw std::runtime_error("Unknown error");

    default: throw std::runtime_error("Unknown error!");
  }
}

std::string string_error(int ier)
{
  switch(ier) {
    case SUCCESS: return "Success";
    case INCOMPATIBLE_MODELS: return "Incompatible submodels";
    case LINALG_FAILURE: return "Linear algebra call failed";
    case MAX_ITERATIONS: return "Maximum iteration count exceeded";
    case KT_VIOLATION: return "Integration of rate-independent model resulted in a violation of the Kuhn-Tucker conditions";
    case NODE_NOT_FOUND: return "XML node not found";
    case TOO_MANY_NODES: return "More than  one XML node found";
    case ATTRIBUTE_NOT_FOUND: return "XML attribute not found";
    case UNKNOWN_TYPE: return "Unknown model type";
    case BAD_TEXT: return "Bad text data in XML node";
    case INVALID_TYPE: return "Type described by XML node is invalid here";
    case FILE_NOT_FOUND: return "File not found";
    case CREEP_PLASTICITY: return "Creep models can only be combined with rate independent plasticity models";
    case INCOMPATIBLE_KM: return "Incompatible lengths in Kocks-Mecking region model: number of models = number of splits + 1";
    case DUMMY_ELASTIC: return "Calling for elastic constants from a dummy elastic model";
    case INCOMPATIBLE_VECTORS: return "Inputs do not have the same length.";
    case UNKNOWN_ERROR: return "Unknown error";

    default: return "Unknown error";
  }
}

LinalgError::LinalgError(const char * m) :
    std::runtime_error(m)
{

}


} // namespace neml
