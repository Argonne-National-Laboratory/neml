#include "cxxstring.h"

#include <algorithm>
#include <string>

int main(int argc, char** argv)
{
  std::string thing_you_need_to_make = R"(<materials><test_j2iso type="SmallStrainRateIndependentPlasticity"><elastic type="IsotropicLinearElasticModel"><m1>103561.64383561644</m1><m1_type>youngs</m1_type><m2>0.2945205479452055</m2><m2_type>poissons</m2_type></elastic><flow type="RateIndependentAssociativeFlow"><surface type="IsoJ2"/><hardening type="LinearIsotropicHardeningRule"><s0>100.0</s0><K>1000.0</K></hardening></flow></test_j2iso></materials>)";

  if (argc != 5) {
        printf("Expected 4 arguments:\n");
        printf("\tmax strain, max time, nsteps, temperature.\n");
        return -1;
  }
  
  double e = std::atof(argv[1]);
  double t = std::atof(argv[2]);
  int n = std::atoi(argv[3]);
  double T = std::atof(argv[4]);
  std::string model_name("test_j2iso");

  auto model = parse_string_unique(thing_you_need_to_make,model_name);

  double * h_n = new double [model->nstore()];
  double * h_np1 = new double [model->nstore()];

  model->init_hist(h_n);

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

  int ier;

  for (int i=0; i<n; i++) {
    t_np1 = (i+1) * t / ((double) n);
    std::fill(e_np1, e_np1+6, 0.0);
    e_np1[0] = (i+1) * e / ((double) n);

    ier = model->update_sd(e_np1, e_n, T_np1, T_n, t_np1, t_n, s_np1, s_n,
                           h_np1, h_n, A_np1, u_np1, u_n, p_np1, p_n);

    if (ier != 0) {
      throw std::runtime_error("Update went bad");
    }

    ier = model->elastic_strains(s_np1, T_np1, h_np1, estrain);

    if (ier != 0) {
      throw std::runtime_error("Elastic strains went bad");
    }

    std::copy(e_np1, e_np1+6, e_n);
    std::copy(s_np1, s_np1+6, s_n);
    std::copy(h_np1, h_np1+model->nstore(), h_n);

    t_n = t_np1;
    T_n = T_np1;

    u_n = u_np1;
    p_n = p_np1;

  }
  
  delete [] h_n;
  delete [] h_np1;

  return 0;
}
