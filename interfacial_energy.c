#include "binary.h"

void Interfacial_energy(double *hphi)
{
	G_bar = 0.0;
	G = 0.0;
  double f[nx];
	unsigned FLAG;

	double cprime[nx], phiprime[nx];

	double dkx, kx[nx];
	double tempreal[nx], tempreal_phi[nx];
	for (int j = 0; j < nx; j++) {
		tempreal[j] = dfdc[j][Re];
		tempreal_phi[j] = phi[j][Re];
	}
/*FILE *Q;
Q = fopen("tempreal","w");
for (int i=0; i<nx; i++)
fprintf(Q,"%lf\n",tempreal_phi[i]);
fclose(Q);*/

	for (int i=0; i<nx; i++){
		G = G + A * (tempreal[i] - c_alpha) * (tempreal[i] - c_alpha) * (1.0 - hphi[i]) 
              + (B * (tempreal[i] - c_beta1) * (tempreal[i] - c_beta1) * (tempreal[i] - c_beta2) * (tempreal[i] - c_beta2 )) * hphi[i] 
              + (1.0 - chi * tempreal[i]) * P * tempreal_phi[i] * tempreal_phi[i] * (1-tempreal_phi[i]) * (1-tempreal_phi[i]);

//		G_bar =  D*c_beta2*0.5;
	}

	InterfacialEnergy = (G - G_bar);

	printf("G is %le\n", G);
//	printf("G_bar is %lf\n",G_bar);
	printf("InterfacialEnergy is %le\n", InterfacialEnergy);


}
