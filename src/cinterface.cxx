#include "cinterface.h"
#include "nemlerror.h"

NEMLMODEL * create_nemlmodel(const char * fname, const char * mname, int * ier)
{
  try {
    // Remember, releasing the unique_ptr means you have to reference count!
    std::unique_ptr<neml::NEMLModel> umodel = neml::parse_xml_unique(fname, mname);
    *ier = 0;

    return umodel.release();
  }
  catch (...) {
    *ier = -1;
    return NULL;
  }
}

void destroy_nemlmodel(NEMLMODEL * model, int * ier)
{
  try {
    delete model;
    *ier = 0;
  }
  catch (...) {
    *ier = -1;
  }
}

double alpha_nemlmodel(NEMLMODEL * model, double T)
{
  try {
    return model->alpha(T);
  }
  catch (...) {
    return -1;
  }
}

void elastic_strains_nemlmodel(NEMLMODEL * model, double * s_np1,
                               double T_np1, double * h_np1,
                               double * e_np1, int * ier)
{
  try {
     model->elastic_strains(s_np1, T_np1, h_np1, e_np1);
  }
  catch (...) {
    *ier = -1;
  }
}

int nstore_nemlmodel(NEMLMODEL * model)
{
  try {
    return model->nstore();
  }
  catch (...) {
    return -1;
  }
}

void init_store_nemlmodel(NEMLMODEL * model, double * store, int * ier)
{
  try {
    model->init_store(store);
  }
  catch (...) {
    *ier = -1;
  }
}

void update_sd_nemlmodel(NEMLMODEL * model, double * e_np1, double * e_n,
                         double T_np1, double T_n,
                         double t_np1, double t_n,
                         double * s_np1, double * s_n,
                         double * h_np1, double * h_n,
                         double * A_np1,
                         double * u_np1, double u_n,
                         double * p_np1, double p_n,
                         int * ier)
{
  try {
    model->update_sd(e_np1, e_n, T_np1, T_n, t_np1, t_n, s_np1, s_n,
                            h_np1, h_n, A_np1, *u_np1, u_n, *p_np1, p_n);
  }
  catch (...) {
    *ier = -1;
  }
}
