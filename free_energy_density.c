#include "binary.h"

void Free_Energy_Density(int steps){

	double hphi1;
        double f[nx];

  for(int i=0; i<nx; i++){
    hphi1 = dfdphi[i][Re] * dfdphi[i][Re] * dfdphi[i][Re] * (10.0 - 15.0 * dfdphi[i][Re] + 6.0 * dfdphi[i][Re] *dfdphi[i][Re]);

    f[i] = A * (1 - hphi1) * (dfdc[i][Re] - c_alpha) * (dfdc[i][Re] - c_alpha) 
	+ B * hphi1 * (dfdc[i][Re] - c_beta1) * (dfdc[i][Re] - c_beta1) * (dfdc[i][Re] - c_beta2) * (dfdc[i][Re] - c_beta2) 
	+ (1.0 - chi * dfdc[i][Re]) * P * dfdphi[i][Re] * dfdphi[i][Re] * (1.0 - dfdphi[i][Re]) *  (1.0 - dfdphi[i][Re]);
  }
  
  char fn[100];

  FILE *q;
  sprintf(fn, "With_Constraint_FEDensity.%06d", steps);
 // sprintf(fn, "Without_Constraint_FEDensity.%06d", steps);
  q = fopen(fn, "w");
  for(int i=0; i<nx; i++)
     fprintf(q,"%le\n", f[i]);

  fclose(q);

}
