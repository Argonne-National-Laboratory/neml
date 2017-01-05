#include "cinterface.h"

NEMLMODEL * create_nemlmodel(const char * fname, const char * mname, int * ier)
{
  // Gar this is tricky, we need to let C++ know that I don't want it messing
  // with releasing the object.  Really I want a unique_ptr, but I didn't write
  // the interface that way...
  //

}

void destroy_nemlmodel(NEMLMODEL * model, int * ier)
{

}
