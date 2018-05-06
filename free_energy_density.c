#include "binary.h"

void Free_Energy_Density(int steps){

	double hphi1;
  double f[nx];
  char fn[100];
  FILE *q;
  //abcd
  for(int i=0; i<nx; i++){
    hphi1 = phi[i][Re] * phi[i][Re] * phi[i][Re] * (10.0 - 15.0 * phi[i][Re] + 6.0 * phi[i][Re] *phi[i][Re]);

    f[i] = A * (1 - hphi1) * (dfdc[i][Re] - c_alpha) * (dfdc[i][Re] - c_alpha) 
	+ B * hphi1 * (dfdc[i][Re] - c_beta1) * (dfdc[i][Re] - c_beta1) * (dfdc[i][Re] - c_beta2) * (dfdc[i][Re] - c_beta2) 
	+ (1.0 - chi * dfdc[i][Re]) * P * phi[i][Re] * phi[i][Re] * (1.0 - phi[i][Re]) *  (1.0 - phi[i][Re]);
  }
  
/*  sprintf(fn, "With_Constraint_FEDensity.%06d", steps);
  q = fopen(fn, "w");
  for(int i=0; i<nx; i++)
     fprintf(q,"%le\n", f[i]);
  fclose(q);
  */
}
