#include "report.h"

#include <iostream>

#include "parse.h"

int main(int argc, char ** argv) {

  if (argc == 2) {
    std::cout << "Each model name (and its type) in " << argv[1] << " are:" << std::endl;
    neml::print_model_names(argv[1]);
    return 0;
  } else if (argc != 3) {
    std::cout << "Need two command line arguments: An XML filename and a model name from that file." << std::endl;
    return -1;
  }

  auto model = neml::parse_xml(argv[1], argv[2]);

  int n = model->nstore();

  double * ihist = new double[n];
  
  model->init_hist(ihist);

  std::cout << std::endl;
  
  std::cout << "*DEPVAR" << std::endl;
  std::cout << n << std::endl;
  std::cout << "**" << std::endl;

  std::cout << "*INITIAL CONDITIONS, TYPE=SOLUTION" << std::endl;
  for (int i=0; i<n; i++) {
    std::cout << ihist[i] << ", ";
  }
  std::cout << std::endl;

  std::cout << "**" << std::endl;

  std::cout << std::endl;

  delete [] ihist;

  return 0;
}
