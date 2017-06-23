#include <memory>
#include <algorithm>

#include "neml_interface.h"

double random(double LO, double HI)
{
  return LO + static_cast<double> (rand()) /( static_cast<double> 
                                              (RAND_MAX/(HI-LO)));
}

int main(int argc, char ** argv) 
{
  int ier; 
  std::unique_ptr<neml::NEMLModel> model = neml::parse_xml(
      "models.xml", "chaboche_600", ier);
  
  int nmodel = 200;
  int nsteps = 150;
  int nhist = model->nhist();

  // Generate history
  double * history_n = new double[nmodel*nhist];
  double * history_np1 = new double[nmodel*nhist];
  for (int i=0; i<nmodel; i++) {
    model->init_hist(&history_n[i*nhist]);
  }

  // Various factors
  double ss_l = -0.0005;
  double ss_h = 0.001;

  // Strains
  double * e_n = new double[nmodel*6];
  double * e_np1 = new double[nmodel*6];
  std::fill(e_n, e_n+nmodel*6, 0.0);
  for (int i=0; i<nmodel; i++) {
    for (int j=0; j<6; j++) {
      e_np1[j+i*6] = random(ss_l, ss_h);
    }
  }

  // Other storage
  double * s_n = new double[nmodel*6];
  std::fill(s_n, s_n+nmodel*6, 0.0);
  double * s_np1 = new double[nmodel*6];
  double * A_np1 = new double[nmodel*36];
  double * u_n = new double[nmodel];
  double * u_np1 = new double[nmodel];
  std::fill(u_n, u_n+nmodel, 0.0);
  double * p_n = new double[nmodel];
  double * p_np1 = new double[nmodel];
  std::fill(p_n, p_n+nmodel, 0.0);

  double T_np1 = 0.0;
  double T_n = 0.0;

  double t_np1 = 0.01;
  double t_n = 0.0;
 
  // Start our hugeo loop
  for (int i = 0; i < nsteps; i++) {
    for (int j = 0; j < nmodel; j++) {
      model->update_sd(&e_np1[j*6], &e_n[j*6], T_np1, T_n, t_np1, t_n,
                    &s_np1[j*6], &s_n[j*6], &history_np1[j*nhist],
                    &history_n[j*nhist], &A_np1[j*36],
                    u_np1[j], u_n[j], p_np1[j], p_n[j]);
    }

    // Switch everything
    std::swap(e_np1, e_n);
    std::swap(s_np1, s_n);
    std::swap(history_np1, history_n);
    std::swap(u_np1, u_n);
    std::swap(p_np1, p_n);

    // Get new strains
    for (int j=0; j<nmodel; j++) {
      for (int k=0; k<6; k++) {
        e_np1[j*6+k] = e_n[j*6+k] + random(ss_l, ss_h);
      }
    }
  }

  delete [] history_n;
  delete [] history_np1;
  delete [] e_n;
  delete [] e_np1;
  delete [] s_n;
  delete [] s_np1;
  delete [] A_np1;
  delete [] p_n;
  delete [] p_np1;
  delete [] u_n;
  delete [] u_np1;

  return 0;
}
