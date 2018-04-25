#include"binary.h"
void Evolve(double *hphi)//, fftw_complex *dfdc, fftw_complex *dfdphi)
{

  void Output_Conf(int steps);
  
  int loop_condition, count;
 
  void Free_Energy_Density(int steps);

  double dkx, kx, kpow2, kpow4;
  
  double rc, fp, rc_new, rphi_new;
  
  double rphi, fpphi;
  
  double total;
  double *tempreal;
  double lhs, rhs, lhse, rhse;
  
  double err, maxerror;

  double sum1, sum2, meanPhi1, meanPhi2;
  
  tempreal = (double *) malloc(sizeof(double) * nx);
  
  dkx = 2.0 * PI / ((double) nx * dx);
  
  for (int i = 0; i < nx; i++) {
      tempreal[i] = comp[i][Re];
  }
  
  loop_condition = 1;
  
  fftw_execute_dft(p_up, comp, comp);
  
  alloycomp = comp[0][Re] * one_by_nx;

//Start the evolution

  for (count = 0; count <= num_steps; count++) {
//printf("time steps %d\n", count);
    if (((count % print_steps) == 0) || (count == num_steps) || (loop_condition == 0)) {
       printf("total_time=%lf\n", sim_time);
       printf("writing configuration to file!\n");
       Output_Conf(count);
       Free_Energy_Density(count);
    }
    if (count > num_steps || loop_condition == 0)
	break;
    
// Evaluate dfdc in real space
   double ctemp;
   double ptemp;
   double gphi, hprime, gprime;

   sum1 = 0.0;
   sum2 = 0.0;

   for (int i = 0; i < nx; i++) {
      ctemp = dfdc[i][Re];	//this dfdc is composition
      ptemp = dfdphi[i][Re];	//this dfdphi is phi
      hphi[i] = ptemp * ptemp * ptemp * (10.0 - 15.0 * ptemp + 6.0 * ptemp * ptemp);
      hprime = 30.0 * (ptemp * ptemp - 2.0 * ptemp * ptemp * ptemp + ptemp * ptemp * ptemp * ptemp);
      gphi = (ptemp * ptemp) * (1.0 - ptemp) * (1.0 - ptemp);
      gprime = 2.0 * ptemp - 6.0 * ptemp * ptemp + 4.0 * ptemp * ptemp * ptemp;

      dfdc[i][Re] = 2.0 * A * (1.0 - hphi[i]) * (ctemp - c_alpha) + 2.0 * B * hphi[i] * (ctemp - c_beta1) * (ctemp - c_beta2) * (2.0 *ctemp - c_beta1 - c_beta2) - chi * P * gphi ;
      dfdc[i][Im] = 0.0;	//this dfdc is derivative of f wrt c

      dfdphi[i][Re] = -1.0 * hprime * A * (ctemp - c_alpha) * (ctemp - c_alpha) +	hprime * B * (ctemp - c_beta1) *
                       (ctemp - c_beta1) * (ctemp - c_beta2) * (ctemp - c_beta2) + (1.0 - chi * ctemp) * P * gprime ;
      dfdphi[i][Im] = 0.0;	// this dfdphi is derivative of f wrt phi
   }
   
   
   fftw_execute_dft(p_up, phi, phi);
   fftw_execute_dft(p_up, dfdc, dfdc);
   fftw_execute_dft(p_up, dfdphi, dfdphi);
   
  for (int i = 0; i < nx; i++) {
    if (i <= nx_half)
       kx = (double) i *dkx;
    else
       kx = (double) (i - nx) * dkx;
    
     kpow2 = kx * kx;
     kpow4 = kpow2 * kpow2;
     
     lhs = 1.0 + 2.0 * mobility * kappa_c * kpow4 * dt;
     
     rc = comp[i][Re];
     fp = dfdc[i][Re];
     rhs = rc - mobility * kpow2 * dt * fp;
     rc_new = rhs / lhs;
     comp[i][Re] = rc_new;
     dfdc[i][Re] = comp[i][Re];

     rc = comp[i][Im];
     fp = dfdc[i][Im];
     rhs = rc - mobility * kpow2 * dt * fp;
     rc_new = rhs / lhs;
     comp[i][Im] = rc_new;
     dfdc[i][Im] = comp[i][Im];

     lhse = 1.0 + 2.0 * relax_coeff * kappa_phi * kpow2 * dt;

     rphi = phi[i][Re];
     fpphi = dfdphi[i][Re];
     rhse = rphi - relax_coeff * dt * fpphi;
     rphi_new = rhse / lhse;
     phi[i][Re] = rphi_new;

     //printf("%lf\n", dfdphi[i][Re]);
     rphi = phi[i][Im];
     fpphi = dfdphi[i][Im];
     rhse = rphi - relax_coeff * dt * fpphi;
     rphi_new = rhse / lhse;
     phi[i][Im] = rphi_new;

  }


//Check for conservation of mass
    total = dfdc[0][Re] * one_by_nx;
    err = fabs(total - alloycomp);
    if (err > COMPERR) {
       printf("ELEMENTS ARE NOT CONSERVED,SORRY!!!!\n");
       printf("error=%lf\n", err);
       exit(0);
    }
    
    fftw_execute_dft(p_dn, dfdc, dfdc);
    fftw_execute_dft(p_dn, dfdphi, dfdphi);

    for (int i = 0; i < nx; i++) {
      dfdc[i][Re] *= one_by_nx;
      dfdc[i][Im] *= one_by_nx;

      dfdphi[i][Re] *= one_by_nx;
      dfdphi[i][Im] *= one_by_nx;
     
      sum1 += dfdphi[i][Re];
      sum2 += dfdphi[i][Im];
    }

  meanPhi1 = sum1 * one_by_nx;
  meanPhi2 = sum2 * one_by_nx;
 // printf("meanPhi1 = %le\n", meanPhi1);

  for (int i=0; i< nx; i++){
      dfdphi[i][Re] -= meanPhi1;
      dfdphi[i][Im] -= meanPhi2;
  }

//Check for bounds 
  for (int i = 0; i < nx; i++) {
     if (dfdc[i][Re] < -0.2 || dfdc[i][Re] > 1.2) {
	printf("Compositions out of bounds. Exiting\n");
	exit(0);
     }
  }

//Check for convergence 
   maxerror = 0.0;
  for (int i = 0; i < nx; i++) {
    err = fabs(tempreal[i] - dfdc[i][Re]);
    if (err > maxerror)
	maxerror = err;
  }

    if (maxerror <= Tolerance) {
      printf("maxerror=%lf\tnumbersteps=%d\n", maxerror, count);
      loop_condition = 0;
    }
     sim_time = sim_time + dt;
  for (int i = 0; i < nx; i++) {
     tempreal[i] = dfdc[i][Re];
  }

  }

}
