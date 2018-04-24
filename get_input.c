#include"binary.h"
void Get_Input_Parameters (char *fnin, char *fnout)
{
	FILE *fpin, *fpcout;
	char param[100], fn[100];

	if (!(fpcout = fopen (fnout, "w"))) {
		printf ("File:%s could not be opened \n", fnout);
		exit (1);
	}
	fprintf (fpcout, "The name of this file is : %s \n", fnout);
	fprintf (fpcout, "Input is from            : %s \n", fnin);

	if (!(fpin = fopen (fnin, "r"))) {
		printf ("File: %s could not be opened \n", fnin);
		exit (1);
	}

	if(fscanf (fpin, "%s%d", param, &nx));
	if(fscanf (fpin, "%s%lf", param, &dx));
	if(fscanf (fpin, "%s%le",param, &dt));
	if(fscanf (fpin, "%s%d", param, &num_steps));
	if(fscanf (fpin, "%s%d", param, &print_steps));
	if(fscanf (fpin, "%s%lf", param, &R));
	if(fscanf (fpin, "%s%lf%s%lf%s%lf%s%lf", param, &A, param, &B, param, &chi, param, &P));
	if(fscanf (fpin, "%s%lf%s%lf%s%lf", param,  &c_alpha,  param, &c_beta1, param, &c_beta2));
	if(fscanf (fpin, "%s%lf", param, &kappa_c));
	if(fscanf (fpin, "%s%lf", param, &kappa_phi));
	if(fscanf (fpin, "%s%lf", param, &mobility));
	if(fscanf (fpin, "%s%lf", param, &relax_coeff));
	if(fscanf (fpin, "%s%d", param, &initflag));
	if(fscanf (fpin, "%s%d", param, &initcount));
	if(fscanf (fpin, "%s%d", param, &fftw_flag));

	printf("nx=%d\n",nx);
	printf("dx=%lf\n",dx);
	printf("dt=%le\n",dt); 
	printf("Simulation steps = %d\n", num_steps);
	printf("Radius = %lf\n",R);
	printf("Bulk free energy coefficients A=%lf\tB=%lf\tchi=%lf\tP=%lf\n", A, B, chi, P);
	printf("C_alpha = %lf\n", c_alpha);
	printf("C_beta1 = %lf\n", c_beta1);
	printf("C_beta2 = %lf\n", c_beta2);
	printf("Kappa_c = %lf\n", kappa_c);
	printf("Kappa_phi = %lf\n", kappa_phi);
	printf("Mobility = %lf\n", mobility);
	printf("Relax_Coeff = %lf\n", relax_coeff);
	printf("Initflag = %d\n", initflag);
	printf("File read from steps = %d\n", initcount);
	if(kappa_c <= 1.0e-08) {
		printf("Warning: too small or negative values for gradient energy coefficients\n");
		printf("Using default values\n");
		kappa_c = 1.0;
	}
	if(mobility <= 1.0e-08) {
		printf("Warning: too small or negative values for kinetic coefficient\n");
		printf("Using default values\n");
		mobility = 1.0;
	}
	if(initflag == 0) {
		printf("Configuration initialized by me\n");
		initcount = 0;
	} else {
		sprintf(fn,"conf.%06d", initcount);     
		printf("Configuration read from file %s\n",fn);
	}
	printf("FFTW_FLAG = %d\n",fftw_flag);
	fclose (fpin);

	fprintf (fpcout, "nx %d\n", nx);
	fprintf (fpcout, "dx %lf\n",dx);
	fprintf (fpcout, "dt %lf\n",dt);
	fprintf (fpcout, "num_steps %d\n", num_steps);
	fprintf (fpcout, "A   %3.2f   B   %3.2f   P   %3.2f\n", A, B, P);
	fprintf(fpcout,"C_alpha = %lf\n", c_alpha);
	fprintf(fpcout,"C_beta1 = %lf\n", c_beta1);
	fprintf(fpcout,"C_beta2 = %lf\n", c_beta2);
	fprintf (fpcout, "kappa_c %lf\t kappa_phi %lf\n", kappa_c, kappa_phi);
	fprintf (fpcout, "mobility %lf\t relax_coeff %lf\n", mobility, relax_coeff);
	// fprintf (fpcout, "Calc_interface_energy  %d\n", calc_interface_energy);
	fprintf (fpcout, "initflag %d\n", initflag);
	fprintf (fpcout, "initcount %d\n", initcount);
	fprintf (fpcout, "fftw_flag %d\n", fftw_flag);
	fclose(fpcout);
}
