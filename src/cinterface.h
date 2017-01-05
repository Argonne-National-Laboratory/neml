#ifndef CINTERFACE_H
#define CINTERFACE_H

// The opaque pointer
typedef struct NEMLMODEL NEMLMODEL;

#ifdef __cplusplus

#include "neml.h"
#include "parse.h"

#include <string>

extern "C" {
#endif

// C interface
NEMLMODEL * create_nemlmodel(const char * fname, const char * mname, int * ier);
void destroy_nemlmodel(NEMLMODEL * model, int * ier);


#ifdef __cplusplus
}
#endif


#endif // CINTERFACE_H
