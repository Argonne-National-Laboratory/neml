#include "cinterface.h"
#include "nemlerror.h"

NEMLMODEL * create_nemlmodel(const char * fname, const char * mname, int * ier)
{
  try {
    // Remember, releasing the unique_ptr means you have to reference count!
    std::unique_ptr<neml::NEMLModel> umodel = neml::parse_xml(fname, mname, *ier);
    return umodel.release();
  }
  catch (...) {
    *ier = neml::UNKNOWN_ERROR;
    return NULL;
  }
}

void destroy_nemlmodel(NEMLMODEL * model, int * ier)
{
  try {
    delete model;
  }
  catch (...) {
    *ier = neml::UNKNOWN_ERROR;
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
    *ier = model->init_store(store);
  }
  catch (...) {
    *ier = neml::UNKNOWN_ERROR;
  }
}

void update_sd_nemlmodel(NEMLMODEL * model, double * e_np1, double * e_n,
                         double T_np1, double T_n,
                         double t_np1, double t_n,
                         double * s_np1, double * s_n,
                         double * h_np1, double * h_n,
                         double * A_np1, int * ier)
{
  try {
    *ier = model->update_sd(e_np1, e_n, T_np1, T_n, t_np1, t_n, s_np1, s_n,
                            h_np1, h_n, A_np1);
  }
  catch (...) {
    *ier = neml::UNKNOWN_ERROR;
  }
}
