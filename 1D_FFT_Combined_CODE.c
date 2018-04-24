#include"binary.h"
//#include "get_input.c"
//#include "init_conf.c"
//#include "evolve.c"
//#include "out_conf.c"
//#include "interfacial_energy.c"

int main(void)
{

	void Get_Input_Parameters(char *fnin, char *fnout);
	void Init_Conf();
	void Read_Restart();
	void Evolve(double *hphi);//, fftw_complex *dfdc, fftw_complex *dfdphi);
	void Interfacial_energy(double *hphi);

	char finput[15] = "bin1ary";
	char fnin[30], fnout[30];

	FILE *fp;

	unsigned FLAG;

	if (!(fp = fopen(finput, "r"))) {
		printf("File: %s could not opened\n", finput);
		exit(EXIT_FAILURE);
	}

	if (fscanf(fp, "%s", fnin) == 1) {
		printf("Input Parameters Filename: %s\n", fnin);
	}

	if (fscanf(fp, "%s", fnout) == 1) {
		printf("Output Parameters Filename: %s\n", fnout);
	}

	if (!(fpout = fopen (fnout, "w"))) {
		printf ("File:%s could not be opened\n", fnout);
		exit (EXIT_FAILURE);
	}

	fclose(fp);

	Get_Input_Parameters(fnin, fnout);


	comp = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nx);
	dfdc = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nx);
	phi = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nx);
	dfdphi = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nx);

	nx_half = nx / 2;

	one_by_nx = 1 /(double) nx;

	FLAG = FFTW_ESTIMATE;

	if (fftw_flag == 1)
		FLAG = FFTW_MEASURE;

	if (fftw_flag == 2)
		FLAG = FFTW_PATIENT;

	if (fftw_flag == 3)
		FLAG = FFTW_EXHAUSTIVE;

	p_up = fftw_plan_dft_1d(nx, comp, comp, FFTW_FORWARD, FLAG);
	p_dn = fftw_plan_dft_1d(nx, comp, comp, FFTW_BACKWARD, FLAG);

	if (initflag == 0) {
		Init_Conf();
	} else
		Read_Restart();

	sim_time = 0.0;

	double hphi[nx];

	Evolve(hphi);//, dfdc, dfdphi);

	Interfacial_energy(hphi);

	fclose(fpout);

	fftw_destroy_plan(p_up);
	fftw_destroy_plan(p_dn);

	fftw_free(comp);
	fftw_free(dfdc);
	fftw_free(phi);
	fftw_free(dfdphi);
	//fftw_cleaup();
	return 0;
}
