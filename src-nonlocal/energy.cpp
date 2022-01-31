// This compute the linear and nonlinear energies, and total wave action


#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>
#include "const.h"
#include "functions.h"
#include <fftw3.h>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace arma;

static double Lin_E,Non_E, Int,k,wave_diss_alpha,wave_diss_nu,energy_diss_alpha,energy_diss_nu, psi_hat2, psi2;
static cx_rowvec psi2_hat(N,fill::zeros), k2(M,fill::zeros);
void define_k2(cx_rowvec & k2_temp);
void embed_N_M( cx_rowvec A, cx_rowvec & B);
void embed_M_N( cx_rowvec B, cx_rowvec & A);
void dealias(cx_rowvec & A);

void energy(double runtime, int out_count, cx_rowvec psi_hat, cx_rowvec psi, double & E_out,cx_rowvec & psi_M,cx_rowvec & psi_hat_M,fftw_plan FFT, fftw_plan IFFT){
    
   /*
    Subroutine that computes the components of the total energy and also energy and waveaction dissipation rates


   */

		k=0.0;
	    wave_diss_alpha = 0.0;
	    wave_diss_nu = 0.0;
	    energy_diss_alpha = 0.0;
	    energy_diss_nu = 0.0;
    
	    Lin_E = 0.0;
	    Non_E = 0.0;
	    Int = 0.0;


  	define_k2(k2);
    embed_N_M(psi_hat,psi_hat_M);
 
    fftw_execute(IFFT);
    psi_M %= conj(psi_M);   
    fftw_execute(FFT);
    psi_hat_M /= double(M);
    psi_hat_M %= pow( 1.0 + (k2/g) ,-0.5);
    embed_M_N(psi_hat_M, psi2_hat);

   
	    for(int i = 0; i < N/2; i++){

	        k = 2.0 * pi * double(i) / Lx;
	        psi_hat2 = pow( abs(psi_hat(i)),2.0);
	        psi2 = pow(abs(psi(i)),2.0);
	        Lin_E += Lx * k*k  * psi_hat2; 
	        Non_E -= Lx* 0.25*pow(abs(psi2_hat(i)),2.0);
	        Int += dx  * psi2;
        
	        wave_diss_alpha += 2.0*alpha*pow(k*k,alphapower) * psi_hat2 ;
	        wave_diss_nu += 2.0*nu*pow(k*k,nupower) * psi_hat2 ;
	        energy_diss_alpha += 2.0*alpha*pow(k*k,alphapower+1.0) * psi_hat2;
	        energy_diss_nu += 2.0*nu*pow(k*k,nupower+1.0) * psi_hat2;

	        k = abs(2.0 * pi * double(-i-1) / Lx);
        
	        psi_hat2 = pow( abs(psi_hat(N-i-1)),2.0);
	        psi2 = pow(abs(psi(N-i-1)),2.0);
        
	        Lin_E += Lx * k*k * psi_hat2;
	        Non_E -= Lx * 0.25*pow(abs(psi2_hat(N-i-1)),2.0);
	        Int += dx * psi2;

	        wave_diss_alpha += 2.0*alpha*pow(k*k,alphapower) * psi_hat2 ;
	        wave_diss_nu += 2.0*nu*pow(k*k,nupower) * psi_hat2 ;
	        energy_diss_alpha += 2.0*alpha*pow(k*k,alphapower+1.0) * psi_hat2;
	        energy_diss_nu += 2.0*nu*pow(k*k,nupower+1.0) * psi_hat2;

	       
	        E_out = Lin_E + Non_E;


	        
	    }
   

	    ostringstream out_D;
	    out_D << "./output/dis." << setw(6) << setfill('0') << out_count << ends;         //creates file name for outputting data at time slice
	    string filenameD = out_D.str();
	    ofstream fout_D(filenameD.c_str());
	    fout_D << scientific;
	    fout_D.precision(12);

	    fout_D << runtime << " " << wave_diss_alpha << " " << wave_diss_nu << " " << energy_diss_alpha << " " << energy_diss_nu << endl;
    
	    fout_D.close();

	    ostringstream out_E;
	    out_E << "./output/energy." << setw(6) << setfill('0') << out_count << ends;			//creates file name for outputting data at time slice
	    string filenameE = out_E.str();
	    ofstream fout_E(filenameE.c_str());
	    fout_E << scientific;
	    fout_E.precision(12);
    
	        fout_E << runtime << " " << Lin_E << " " << Non_E << " " << Lin_E + Non_E << endl;
    
	    fout_E.close();
    
	    ostringstream out_W;
	    out_W << "./output/wave." << setw(6) << setfill('0') << out_count << ends;			//creates file name for outputting data at time slice
	    string filenameW = out_W.str();
	    ofstream fout_W(filenameW.c_str());
	    fout_W << scientific;
	    fout_W.precision(12);

	    fout_W << runtime << " " << Int << " " << pow( abs(psi_hat(0,0)),2.0) << endl;
    
	    fout_W.close();
    

    
	    return;
	}


