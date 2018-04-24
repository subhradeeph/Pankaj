#include"binary.h"
void Output_Conf (int steps)
{
 FILE *fpt;
 
 char fn[100];

 sprintf (fn, "conf.%06d", steps);

 fpt = fopen (fn, "w");
 fwrite (&dfdc[0][0], sizeof(double), 2 * nx , fpt);
 fclose (fpt);

 sprintf (fn, "prof_gp.%06d", steps);
 fpt = fopen (fn, "w");

for (int i = 0; i < nx; i++) {
	 fprintf(fpt,"%d\t%le\t%le\n", i, dfdc[i][Re], dfdphi[i][Re]);
	}
  fprintf(fpt,"\n");
 
 fclose(fpt); 

}
