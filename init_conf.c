#include"binary.h"

void Init_Conf()
{

  FILE *fp;
  char fn[100];
  double sum, mean;
 
  sum = 0.0;
/* Phi profile */
  for (int i = 0; i < nx; i++){
     if(fabs((i - nx/2)*dx) <= R*dx){
	// printf("i=%d\n",i);
	phi[i][Re] = 1.0;
     } else {
	//  printf("j=%d\n",i);
        phi[i][Re] = 0.0;
     }
	phi[i][Im] = 0.0;

	sum += phi[i][Re];
  }

  mean = sum * one_by_nx;

  printf("Initial_mean = %le\n", mean);

/* Composition profile  */

  for(int i = 0; i < nx; i++){
    if(fabs(phi[i][Re]) > 1.0e-10){
//if (fabs((i - nx/2)*dx) <= R*dx){

      if (i <= nx/2) 
        comp[i][Re] = c_beta2;
      else 
        comp[i][Re] = c_beta2;
    } else{
      comp[i][Re] = c_alpha;
    }
      comp[i][Im] = 0.0;
   }

  for (int i = 0; i < nx; i++) {
     dfdc[i][Re] = comp[i][Re];
     dfdc[i][Im] = comp[i][Im];
     dfdphi[i][Re] = phi[i][Re];
     dfdphi[i][Im] = phi[i][Im];
  }

  sprintf(fn, "profile.in");  
  if (!(fp = fopen (fn, "w"))) {
     printf ("File:%s could not be opened \n", fn);
     exit (1);
  }

  for (int i = 0; i < nx; i++) 
     fprintf(fp,"%d\t%le\t%le\n",i, comp[i][Re], phi[i][Re]);

  fclose(fp);
}

void Read_Restart()
{
  FILE *fpread;
  char fr[100];

  sprintf (fr,"conf.%06d", initcount);
  fpread = fopen (fr, "r");
  if(fread (&comp[0], sizeof(double), 2 * nx, fpread));
  fclose (fpread);

  for (int i = 0; i < nx; i++) {
     dfdc[i][Re] = comp[i][Re];
     dfdc[i][Im] = comp[i][Im];
     dfdphi[i][Re] = phi[i][Re];
     dfdphi[i][Im] = phi[i][Im];
  }

}
