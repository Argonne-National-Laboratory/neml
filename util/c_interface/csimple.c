#include "csimple.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
      if (argc != 7) {
            printf("Expected 6 arguments:\n");
            printf("\tXML file, model name, max strain, max time, nsteps, temperature.\n");
            return -1;
      }
      
      double e = atof(argv[3]);
      double t = atof(argv[4]);
      int n = atoi(argv[5]);
      double T = atof(argv[6]);

      int ier;
      NEMLMODEL * model = create_nemlmodel(argv[1], argv[2], &ier);
      
      if (ier != 0) {
            printf("Error in creating model.\n");
            return -1;
      }

      double alpha = alpha_nemlmodel(model, 100.0);
      
      // Allocate
      int nstore = nstore_nemlmodel(model);

      double * h_n = malloc(sizeof(double) * nstore);
      init_store_nemlmodel(model, h_n, &ier);
      if (ier != 0) {
            printf("Error in initializing history.\n");
            free(h_n);
            return -1;
      }
      
      double * h_np1 = malloc(sizeof(double) * nstore);

      double e_n[6], e_np1[6], s_n[6], s_np1[6], A_np1[36];
      int i;
      for (i=0; i<6; i++) {
            e_n[i] = 0.0;
            s_n[i] = 0.0;
      }

      double T_np1 = T;
      double T_n = T;

      double t_n = 0.0;
      double t_np1;

      double u_n = 0.0;
      double u_np1;

      double p_n = 0.0;
      double p_np1;

      double ec;
      int j;

      double estrain[6];

      for (i=0; i<n; i++) {
            t_np1 = (i+1) * t / ((double) n);
            for (j=0; j<6; j++) e_np1[j] = 0.0;
            e_np1[0] = (i+1) * e / ((double) n);

            update_sd_nemlmodel(model, e_np1, e_n, T_np1, T_n, t_np1, t_n,
                        s_np1, s_n, h_np1, h_n, A_np1, &u_np1, u_n,
                        &p_np1, p_n, &ier);
            if (ier != 0) {
                  printf("Problem in stress update\n");
                  free(h_n);
                  return -1;
            }

            // printf("Stress: %lf\n", s_np1[0]);

            elastic_strains_nemlmodel(model, s_np1, T_np1, h_np1, estrain, &ier);

            for (j=0; j<6; j++) {
                  s_n[j] = s_np1[j];
                  e_n[j] = e_np1[j];
            }
            for (j=0; j<nstore; j++) {
                  h_n[j] = h_np1[j];
            }
            t_n = t_np1;

            u_n = u_np1;
            p_n = p_np1;
      }


      // Free
      destroy_nemlmodel(model, &ier);
      if (ier != 0) {
            printf("Error in destroying model.\n");
            free(h_n);
            return -1;
      }

      free(h_n);
      free(h_np1);
}
