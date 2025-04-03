// This compute the linear and nonlinear energies, and total wave action


#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>
#include <fftw3.h>
#include <fstream>
#include <iomanip>
#include "const.h"
#include "nonlinear.h"

using namespace std;
using namespace arma;

static double k,k2, psi_hat2,flux,eflux;

static double k_bin = 2.0 * pi / Lx;

static int ib;
static cx_rowvec psi2_hat(N,fill::zeros), k2_vector(M,fill::zeros);

static cx_rowvec nonlinear_hat(N,fill::zeros),temp_M(M,fill::zeros),E_hat(N,fill::zeros),E_hat_dot(N,fill::zeros), gradx_psi_dot_hat(N,fill::zeros), psi_dot_hat(N,fill::zeros), gradx_psi_hat(N,fill::zeros), rho_hat(N,fill::zeros),rho_dot_hat(N,fill::zeros);

static cx_rowvec lap(N, fill::zeros), gradx(N, fill::zeros);

void define_k2(cx_rowvec & k2_temp);
void embed_N_M( cx_rowvec A, cx_rowvec & B);
void embed_M_N( cx_rowvec B, cx_rowvec & A);
void dealias(cx_rowvec & A);
void generateLaplacian(cx_rowvec & );
void generateGradx(cx_rowvec & );


void computeEnergy(cx_rowvec psi_hat, cx_rowvec psi, double & linear_energy_out, double & potential_energy_out, double & nonlinear_energy_out,cx_rowvec & psi_M,cx_rowvec & psi_hat_M,fftw_plan FFT, fftw_plan IFFT){
    
   /*
    Subroutine that computes the components of the energy
   */

	k=0.0;
	    
	linear_energy_out = 0.0;
	nonlinear_energy_out = 0.0;
	potential_energy_out = 0.0;

	for(int i = 0; i < N/2; i++){

	    k = 2.0 * pi * double(i) / Lx;
	    psi_hat2 = pow( abs(psi_hat(i)),2.0);
	   
	    linear_energy_out -= c * Lx * k * k  * psi_hat2;
	    potential_energy_out +=  Lx * mu  * psi_hat2;
        
		k = abs(2.0 * pi * double(-i-1) / Lx);
	    psi_hat2 = pow( abs(psi_hat(N-i-1)),2.0);
	    
	    linear_energy_out -= c * Lx * k * k * psi_hat2;
	    potential_energy_out += Lx * mu  * psi_hat2;
	        
	}

/* compute the nonlinear term based on model */		


	if( FLAG_MODEL_TYPE == "SHE"){
  		define_k2(k2_vector);
    	embed_N_M(psi_hat,psi_hat_M);
 
    	fftw_execute(IFFT);
    	psi_M %= conj(psi_M);   
    	fftw_execute(FFT);
    	psi_hat_M /= double(M);
    	psi_hat_M %= pow( 1.0 + (beta_param * k2_vector) ,-0.5);
   		embed_M_N(psi_hat_M, psi2_hat);

   
	    for(int i = 0; i < N/2; i++){

	        nonlinear_energy_out += Lx * 0.5 * g * pow(abs(psi2_hat(i)),2.0);
	        nonlinear_energy_out += Lx * 0.5 * g * pow(abs(psi2_hat(N-i-1)),2.0);
	               
	    }
	}
	else if( FLAG_MODEL_TYPE == "NLSE"){
  	
	    for(int i = 0; i < N; i++){

	        nonlinear_energy_out += dx  * 0.5 * g * pow(abs(psi(i)),4.0);
	       
	               
	    }
	}

    
	    return;
	}




void computeWaveaction(cx_rowvec psi_hat, double & waveaction_out){
   /*
    Subroutine that computes the total waveaction
   */
	waveaction_out = 0.0;
        	         
	for(int i = 0; i < N/2; i++){

		waveaction_out +=  Lx * pow( abs(psi_hat(i)),2.0);
    		waveaction_out +=  Lx * pow( abs(psi_hat(N-i-1)),2.0);	        	    
    }
return; 
}


void computeDissipationRate( cx_rowvec psi_hat, double & waveaction_dissipation_rate_alpha, double & waveaction_dissipation_rate_nu, double & energy_dissipation_rate_alpha, double & energy_dissipation_rate_nu ){
    
   /*
    Subroutine that computes the components of the energy and waveaction dissipation rates
   */

		k2=0.0;
	    waveaction_dissipation_rate_alpha = 0.0;
	    waveaction_dissipation_rate_nu = 0.0;
	    energy_dissipation_rate_alpha = 0.0;
	    energy_dissipation_rate_nu = 0.0;
    
	    for(int i = 0; i < N/2; i++){

	        k2 = pow(2.0 * pi * double(i) / Lx, 2.0);
	        psi_hat2 = pow( abs(psi_hat(i)),2.0);
	        
	        if(FLAG_HYPER_DISSIPATION == true){
	        	waveaction_dissipation_rate_nu += 2.0*nu*pow(k2,nupower) * psi_hat2 ;
	        	energy_dissipation_rate_nu -= 2.0 * c * nu*pow(k2,nupower+1.0) * psi_hat2;
	        }
	        if(FLAG_HYPO_DISSIPATION == true && (i != 0 )){
			   		waveaction_dissipation_rate_alpha += 2.0*alpha*pow(k2,alphapower) * psi_hat2;
	        		energy_dissipation_rate_alpha -= 2.0 *  c * alpha*pow(k2,alphapower+1.0) * psi_hat2;		
	        }

	        k2 = pow(2.0 * pi * double(-i-1) / Lx , 2.0);
        
	        psi_hat2 = pow( abs(psi_hat(N-i-1)),2.0);
	                    
	        if(FLAG_HYPER_DISSIPATION == true){
	        	waveaction_dissipation_rate_nu += 2.0*nu*pow(k2,nupower) * psi_hat2 ;
	        	energy_dissipation_rate_nu -= 2.0 * c * nu*pow(k2,nupower+1.0) * psi_hat2;
	        }
	        if(FLAG_HYPO_DISSIPATION == true){
	        	waveaction_dissipation_rate_alpha += 2.0*alpha*pow(k2,alphapower) * psi_hat2 ;
	        	energy_dissipation_rate_alpha -= 2.0 * c * alpha*pow(k2,alphapower+1.0) * psi_hat2;
	        }

	        
	    }
	return;
   
}


void computeSpectrum( cx_rowvec psi_hat, rowvec & wave_spectrum){
   /*
    Subroutine that computes wave action spectrum
   */
    wave_spectrum.zeros();
    
    k = 0.0;

	for(int i =0; i< N/2; i++){
        
		k =   2.0*pi * double(i)/ Lx;
		
		ib = int( k / k_bin + 0.5);

		wave_spectrum(ib) += pow( abs(psi_hat(i)) , 2.0);
	
		k = abs( 2.0*pi * double(-i-1)/ Lx );
		ib = int( k / k_bin + 0.5);
	
		wave_spectrum(ib) += pow( abs(psi_hat(N-i-1)) , 2.0);
	
         
	}

return;

}



void computeFlux(cx_rowvec psi_hat, cx_rowvec L,cx_rowvec & psi_M ,cx_rowvec & psi_hat_M, rowvec & wave_flux ,rowvec & energy_flux,fftw_plan FFT, fftw_plan IFFT ){


	k=0.0;
    wave_flux.zeros();
    energy_flux.zeros();

	generateLaplacian(lap);
    generateGradx(gradx);
	define_k2(k2_vector);
    
	nonlinear_hat = nonlinear(psi_hat,0.0, L , psi_M, psi_hat_M, FFT, IFFT);
	gradx_psi_hat = gradx % psi_hat;
	
	psi_dot_hat = ( - complex<double>(0.0, 1.0) * c * lap % psi_hat  - complex<double>(0.0,1.0) * mu * psi_hat + nonlinear_hat);
	gradx_psi_dot_hat = gradx % psi_dot_hat;



    // compute rho_hat
    
	
	embed_N_M(psi_hat,psi_hat_M);
    fftw_execute(IFFT); 
    temp_M = psi_M; // this is for rho_dot_hat
    psi_M %= conj(psi_M);
    fftw_execute(FFT); 
    psi_hat_M /= double(M);
	if( FLAG_MODEL_TYPE == "SHE"){
		psi_hat_M %= pow(1.0 + beta_param * k2_vector,-0.5);
	}
	embed_M_N(psi_hat_M,rho_hat);
    


    // compute rho_dot_hat
    embed_N_M(psi_dot_hat , psi_hat_M);
    fftw_execute(IFFT); 
    psi_M %= 2.0*conj(temp_M);
    fftw_execute(FFT); 
    psi_hat_M /= double(M);
	if(FLAG_MODEL_TYPE == "SHE"){
		psi_hat_M %= pow(1.0 + beta_param * k2_vector,-0.5);
	}
    embed_M_N(psi_hat_M,rho_dot_hat);


    flux = 0.0;
	eflux=0.0;

    for(int i = 0; i< N/2; i++){


            k =  2.0*pi * double(i)/ Lx;
            ib = int(k/k_bin + 0.5);
         
			flux = -2.0 * real( psi_hat(i) * conj(psi_dot_hat(i)));
			eflux = 2.0 * c * real( gradx_psi_hat(i) * conj(gradx_psi_dot_hat(i)) )  + (mu * flux) - g * real( rho_hat(i) * conj(rho_dot_hat(i)) );

			//integrate from 0 to l
            for(int l=ib; l < N/2 +1; l++){
            	wave_flux(l) += flux;
            	energy_flux(l) += eflux;
			}
    


            k = abs(  2.0*pi * double(-i-1)/ Lx);
            ib = int(k/k_bin + 0.5);
            
		
			flux = -2.0 * real( psi_hat(N-i-1) * conj(psi_dot_hat(N-i-1)));
			eflux = 2.0 * c * real( gradx_psi_hat(N-i-1) * conj(gradx_psi_dot_hat(N-i-1)) )  + (mu * flux) - g * real( rho_hat(N-i-1) * conj(rho_dot_hat(N-i-1)) );

			//integrate from 0 to l
			for(int l=ib; l < N/2 + 1; l++){
            	wave_flux(l) += flux;
            	energy_flux(l) +=  eflux;
			}

        
    }
    return;
}


void generateLaplacian(cx_rowvec & A){

   
	for(int i = 0; i < N/2; i++){

	    A(i) = complex<double>(- pow((2.0 * pi * double(i) / Lx),2.0) ,0.0);
	    A(N-i-1) = complex<double>(- pow((2.0 * pi * double(-i-1) / Lx),2.0)  ,0.0);

	}

    return;
}

void generateGradx(cx_rowvec & A){

	for(int i = 0; i < N/2; i++){

	    A(i) = complex<double>(0.0, 2.0 * pi * double(i) / Lx);
	    A(N-i-1) = complex<double>( 0.0,  2.0 * pi * double(-i-1) / Lx);
		   	   
    }

    return;
}
