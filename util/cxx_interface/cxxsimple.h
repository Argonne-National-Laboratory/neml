#ifndef CXXSIMPLE_H
#define CXXSIMPLE_H

#include "parse.h"

#include <string>

using namespace neml;

int main(int argc, char** argv);
NEMLModel * get_model(std::string xml, std::string name);


#endif // CXXSIMPLE_H
