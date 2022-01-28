/*

Programme to solve the 1D NLS equation


Author: Jason Laurie
Date: 08/10/2010


*/


#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>
#include "const.h"
#include "functions.h"
#include <time.h>
#include <fstream>
#include <fftw3.h>
#include <chrono>
#include <omp.h>

using namespace std;
using namespace arma;

int main(){
   
	int out_count =0;
	int avg_count = 0;
	double runtime =0.0;
	double start = 0.0;
	double finish = 0.0;
	double E_out = 0.0;
	start = omp_get_wtime();
	int nthreads = 0;
	
	if(num_threads == 0){
		nthreads = omp_get_max_threads();
	}
	else{
		nthreads = num_threads;
	}

	omp_set_num_threads(nthreads);
	
	

	if(nthreads > 1){
		int err = fftw_init_threads();
		if (err==0){
			cout << "thread creation error: " << err << endl;
		}
		else{
			fftw_plan_with_nthreads(nthreads);
			cout << "number of parallel threads used by fftw = " << nthreads << endl;
		}
	}
	else{
		cout << "using a single thread" << endl;
	}
	
    
    cx_rowvec psi_hat(N,fill::zeros), psi_hat_M(M,fill::zeros), psi(N,fill::zeros), psi_M(M,fill::zeros), step1(N,fill::zeros), step2(N,fill::zeros), step3(N,fill::zeros), step4(N,fill::zeros), L(N,fill::zeros), temp(N,fill::zeros),psi_hat_temp(N,fill::zeros);			//Declare arrays
		
		
	cx_rowvec psi_t(N,fill::zeros),NL1(N,fill::zeros), NL2(N,fill::zeros), NL3(N,fill::zeros), F1(N,fill::zeros),F2(N,fill::zeros),F3(N,fill::zeros), Q1(N,fill::zeros),Q2(N,fill::zeros),Q3(N,fill::zeros),Q4(N,fill::zeros),Q5(N,fill::zeros), E1(N,fill::zeros),E2(N,fill::zeros),noise(N,fill::zeros);	
		
		
		
    rowvec famp(N,fill::zeros);
rowvec wave_spec_avg(NBIN,fill::zeros),wave_flux_avg(NBIN,fill::zeros),energy_spec_avg(NBIN,fill::zeros),energy_flux_avg(NBIN,fill::zeros);
    
    fftw_plan FFT;								//delcare plans
	fftw_plan IFFT;
    fftw_plan FFTN;
    fftw_plan IFFTN;
    
    FFTN = fftw_plan_dft_1d(N, (fftw_complex*) psi.memptr(), (fftw_complex*) psi_hat.memptr(), FFTW_FORWARD, FFTW_PATIENT);
    IFFTN = fftw_plan_dft_1d(N, (fftw_complex*) psi_hat.memptr(), (fftw_complex*) psi.memptr(), FFTW_BACKWARD, FFTW_PATIENT);
    FFT = fftw_plan_dft_1d(M, (fftw_complex*) psi_M.memptr(), (fftw_complex*) psi_hat_M.memptr(), FFTW_FORWARD, FFTW_PATIENT);
	IFFT= fftw_plan_dft_1d(M, (fftw_complex*) psi_hat_M.memptr(), (fftw_complex*) psi_M.memptr(), FFTW_BACKWARD, FFTW_PATIENT);

    
	
	//set up random forcing
	unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();  
	mt19937_64 generator(seed1);
	normal_distribution<double> distribution(0.0,1.0);

    
	
    read(psi, psi_hat, runtime, out_count, FFTN);      //reads in initial data
   
    
	system("mkdir ./output");			//Creates the Folders for Outputting

	
	record_parameters();


	linear(L);     //Creates linear operator
   
  	
	if(FLAG_TIMESTEP_ETDRK == true){
		setup_ETDRK(L,Q1,Q2,Q3,Q4,Q5,F1,F2,F3,E1,E2,dt);
	}

	if(FLAG_FORCE_ON == true){
   
		if(FLAG_FORCE_EXP == true){
			initial_forcing2(famp,kf,dk);
		}
		else{
			initial_forcing1(famp,kf,dk); //Initialise Forcing
		} 
	}
	
   // IniForcing(famp,kf); //initialize forcing profile
    //==============Time stepping routine======================================================
    
    for( int time_step = 1; time_step <= nsteps; time_step++ ){
        
        
        
		runtime += dt;                  //increments runtime
			if(FLAG_FORCE_ON == true){
    
		   	/*  This is the section of code that produces the forcing using normal distributed real and imaginary components at a given wave number*/

				noise.zeros();
				
				for(int i = 0; i < N; i++){
		          		if(famp(i) > 0.0){  
		            		noise(i) = famp(i)*pow(0.5,0.5)*complex<double>(distribution(generator),distribution(generator));   	     			
		          		}
		        	}
		      	
			}
     
	
	 
	 
	 
	 
			start = omp_get_wtime();

				if(FLAG_TIMESTEP_ETDRK == true){

					psi_hat += pow(dt,0.5)*noise; // euler noise
  
					if(FLAG_ETDRK_ORDER == 2){
			        	NL1 = nonlinear(psi_hat,0.0, L , psi_M, psi_hat_M, FFT, IFFT);
			        	step1 = E1%psi_hat + Q1%NL1;
			        	NL2 = nonlinear(step1,dt, L , psi_M, psi_hat_M, FFT, IFFT); 
			        	step2 = step1 +  F1 % (NL2 - NL1);
			        	psi_hat = step2;
			        }
					else if(FLAG_ETDRK_ORDER == 3){
						NL1 = nonlinear(psi_hat,0.0, L , psi_M, psi_hat_M, FFT, IFFT);
						step1 = E2%psi_hat + Q1%NL1;
						NL2 = nonlinear(step1,0.5*dt, L , psi_M, psi_hat_M, FFT, IFFT); 
						step2 = E1%psi_hat + Q2%(2.0*NL2-NL1);
						NL3 = nonlinear(step2,dt, L , psi_M, psi_hat_M, FFT, IFFT);
						step3 = E1 % psi_hat + F1 % NL1 + 4.0* F2 % NL2 + F3 % NL3;
						psi_hat = step3;
					}
					else if(FLAG_ETDRK_ORDER == 4){
						NL1 = nonlinear(psi_hat,0.0,L , psi_M, psi_hat_M, FFT, IFFT);
						step1 = E2%psi_hat + Q1%NL1;
						NL2 = nonlinear(step1,0.5*dt,L ,psi_M, psi_hat_M, FFT, IFFT); 
						step2 = E2%psi_hat + Q2%NL1 + Q3%NL2;
						NL3 = nonlinear(step2,0.5*dt,L , psi_M, psi_hat_M, FFT, IFFT);
						step3 = E1%psi_hat + Q4%NL1 + Q5%NL3;
						step4 = E1%psi_hat + F1%NL1 + 2.0*F2%(NL2+NL3) + F3%nonlinear(step3,dt, L , psi_M, psi_hat_M, FFT, IFFT);
						psi_hat = step4;
					}
			    	//==============================================================================================
				}
			    else{   //=====RK2 with noise===================================
					step1 = dt*nonlinear(psi_hat,0.0, L  ,psi_M, psi_hat_M, FFT, IFFT);
					step2 = dt*nonlinear(psi_hat + step1 + pow(dt,0.5)*noise,dt, L , psi_M, psi_hat_M, FFT, IFFT);
					psi_hat += 0.5*(step1 + step2) + pow(dt,0.5)*noise;
					psi_hat %= exp(dt*L);  //transforms data to u=exp(Ldt)v
    	
			    	//=====================================================================
				}
	 
	 
				if( (time_step % outstep) == 0 ){

						out_count++;
						avg_count++; 
						psi_hat_temp = psi_hat; 

						fftw_execute(IFFTN);
         
						spectrum(psi_hat,wave_spec_avg,energy_spec_avg, out_count,avg_count); //computes spectrum
						psi_t = nonlinear(psi_hat,0.0, L , psi_M, psi_hat_M, FFT, IFFT);
				 		flux(psi_hat, psi_t, wave_flux_avg,energy_flux_avg, out_count,avg_count);
						energy(runtime,out_count,psi_hat, psi,E_out, psi_M,psi_hat_M,FFT,IFFT);       //computes energy
						intensity(runtime,out_count,psi); //computes intensity
						
						psi_hat = psi_hat_temp;
						
						cout << "time = " << runtime << " file = " << out_count << " Energy = " << E_out << endl;        //outputs runtime to screen

						ofstream fout_count("./data/curframe.dat");
						fout_count << scientific;
						fout_count.precision(12);
						fout_count << runtime << " " << out_count << endl;    
				    }
				}
	 
	 
	 
	    
   /*============Runge-Kutta order 4  ========================================
        step1 = dt*NL(temp,0.0, L , psi_M, psi_hat_M, FFT, IFFT);
        temp = complex<double>(0.0,0.0);
        temp =0.5*step1 + psi_hat;
        step2 = dt*NL(temp,0.5*dt, L , psi_M, psi_hat_M, FFT, IFFT);
        temp = complex<double>(0.0,0.0);
        temp =0.5*step2 + psi_hat;
        step3 = dt*NL(temp,0.5*dt, L , psi_M, psi_hat_M, FFT, IFFT);
        temp = complex<double>(0.0,0.0);
        temp =step3 + psi_hat;
        step4 = dt*NL(temp,dt, L , psi_M, psi_hat_M, FFT, IFFT);
        
       psi_hat += (step1 + 2.0*step2 + 2.0*step3 + step4)/6.0;
    //=====================================================================
       
        for(int i =0; i < N ; i++){
            psi_hat(i) *=exp(L(i)*dt);          //tranforms data to u=exp(Ldt)v
        }
    
        
        Forcing(psi_hat, famp);  //additive forcing
        
        
        runtime += dt;                  //increments runtime
   
        
        
        if( (time_step % outstep) ==0 ){
        cout << runtime << endl;        //outputs runtime to screen
            out_count++;
            
            Spectrum(psi_hat,out_count);
        Energy(runtime,psi_hat, psi, IFFTN);       //computes energy
        
        Intensity(runtime,out_count,psi); //computes intensity
        }

        fftw_execute(FFTN);
        psi_hat /= (double) N;
        
    }
    
    
  */
    
    
	 
	fftw_destroy_plan(FFT);
	fftw_destroy_plan(IFFT);			//destroy plans
	fftw_destroy_plan(FFTN);
	fftw_destroy_plan(IFFTN);
	
	finish = omp_get_wtime();                       //end clock
	cout << "time taken for code is = " << finish-start << endl;    //output time
    

	return 0;
	
}
	

