#include "cxxsimple.h"

#include <algorithm>

NEMLModel * get_model(std::string xml, std::string name)
{
  std::unique_ptr<NEMLModel> model1 = parse_xml_unique(xml, name);
  
  return model1.release(); 
}

int main(int argc, char** argv)
{
  if (argc != 7) {
        printf("Expected 6 arguments:\n");
        printf("\tXML file, model name, max strain, max time, nsteps, temperature.\n");
        return -1;
  }
  
  double e = std::atof(argv[3]);
  double t = std::atof(argv[4]);
  int n = std::atoi(argv[5]);
  double T = std::atof(argv[6]);
  
  NEMLModel * model = get_model(argv[1], argv[2]);

  double * h_n = new double [model->nstore()];
  double * h_np1 = new double [model->nstore()];

  model->init_store(h_n);

  double e_n[6], e_np1[6], s_n[6], s_np1[6], A_np1[36];
  std::fill(e_n, e_n+6, 0.0);
  std::fill(s_n, s_n+6, 0.0);

  double T_np1 = T;
  double T_n = T;

  double t_n = 0.0;
  double t_np1;

  double u_n = 0.0;
  double u_np1;

  double p_n = 0.0;
  double p_np1;

  double estrain[6];


  for (int i=0; i<n; i++) {
    t_np1 = (i+1) * t / ((double) n);
    std::fill(e_np1, e_np1+6, 0.0);
    e_np1[0] = (i+1) * e / ((double) n);

    model->update_sd(e_np1, e_n, T_np1, T_n, t_np1, t_n, s_np1, s_n,
                           h_np1, h_n, A_np1, u_np1, u_n, p_np1, p_n);

    model->elastic_strains(s_np1, T_np1, h_np1, estrain);

    std::copy(e_np1, e_np1+6, e_n);
    std::copy(s_np1, s_np1+6, s_n);
    std::copy(h_np1, h_np1+model->nstore(), h_n);

    t_n = t_np1;
    T_n = T_np1;

    u_n = u_np1;
    p_n = p_np1;

  }
  
  delete model;

  delete [] h_n;
  delete [] h_np1;

  return 0;
}
