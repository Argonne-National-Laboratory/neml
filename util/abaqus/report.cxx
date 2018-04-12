#include "report.h"

#include <iostream>

#include "parse.h"

int main(int argc, char ** argv) {

  if (argc != 3) {
    std::cout << "Need two command line arguments: file and model name." << std::endl;
    return -1;
  }

  int ier;
  auto model = neml::parse_xml(argv[1], argv[2], ier);
  if (ier != 0) {
    std:: cout << "Unable to load model " << argv[2] << " from file " << argv[1] << std::endl;
  }

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
