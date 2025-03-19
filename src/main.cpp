/*
Programme to solve the 1D Schrödinger-Helmholtz / nonlinear Schrödinger equation

"SHE"

i PSI_t = c * PSI_xx    +  g * PSI (1 - beta_param * d^2_xx)^-1 |PSI|^2  + mu * PSI    + dissipation + forcing

or

"NLS"

i PSI_t = c * PSI_xx    +  g * PSI |PSI|^2  +   mu * PSI    + dissipation + forcing

Author: Jason Laurie
Date: 27/09/2023

*/


#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>
#include <time.h>
#include <fstream>
#include <fftw3.h>
#include <chrono>
#include <omp.h>
#include "const.h"
#include "functions.h"
#include "linear.h"
#include "forcing.h"
#include "compute.h"
#include "print.h"

using namespace std;
using namespace arma;

int main(){
   
	//define variables counting files and run time 

	int file_number_start = 0;
	int file_number = 0;
	int avg_number = 0;
	double run_time_start = 0.0;
	double run_time = 0.0;
	
	//clock variables
	double start = 0.0;
	double finish = 0.0;
	
	
	//to record diagonostics (conservation laws)
	double waveaction = 0.0;
	double linear_energy = 0.0;
	double potential_energy = 0.0;
	double nonlinear_energy = 0.0;

	//to record dissipation rates
	double waveaction_dissipation_rate_alpha = 0.0;
	double waveaction_dissipation_rate_nu = 0.0;
	double energy_dissipation_rate_alpha = 0.0;
	double energy_dissipation_rate_nu = 0.0;
	
	//arrays for recording spectra and flux
	rowvec wave_spectrum(N/2 +1,fill::zeros), wave_spectrum_avg(N/2 +1 ,fill::zeros), wave_flux(N/2 +1,fill::zeros),energy_flux(N/2 +1,fill::zeros),wave_flux_avg(N/2 +1,fill::zeros),energy_flux_avg(N/2 +1,fill::zeros);
	
	//start clock
	start = omp_get_wtime();


	/*================= openmp (not in use)================
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
	=============================*/
    
	//declare  arrays
    cx_rowvec psi_hat(N,fill::zeros), psi_hat_M(M,fill::zeros), psi(N,fill::zeros), psi_M(M,fill::zeros), L(N,fill::zeros), temp(N,fill::zeros),psi_hat_temp(N,fill::zeros);	
		
		
	cx_rowvec psi_t(N,fill::zeros),NL1(N,fill::zeros), NL2(N,fill::zeros), NL3(N,fill::zeros), F1(N,fill::zeros),F2(N,fill::zeros),F3(N,fill::zeros), Q1(N,fill::zeros),Q2(N,fill::zeros),Q3(N,fill::zeros),Q4(N,fill::zeros),Q5(N,fill::zeros), E1(N,fill::zeros),E2(N,fill::zeros),noise(N,fill::zeros);	
	
    rowvec famp(N,fill::zeros);

    //declare fftw plans
    fftw_plan FFT;							
	fftw_plan IFFT;
    fftw_plan FFTN;
    fftw_plan IFFTN;
    
    FFTN = fftw_plan_dft_1d(N, (fftw_complex*) psi.memptr(), (fftw_complex*) psi_hat.memptr(), FFTW_FORWARD, FFTW_PATIENT);
    IFFTN = fftw_plan_dft_1d(N, (fftw_complex*) psi_hat.memptr(), (fftw_complex*) psi.memptr(), FFTW_BACKWARD, FFTW_PATIENT);
    FFT = fftw_plan_dft_1d(M, (fftw_complex*) psi_M.memptr(), (fftw_complex*) psi_hat_M.memptr(), FFTW_FORWARD, FFTW_PATIENT);
	IFFT= fftw_plan_dft_1d(M, (fftw_complex*) psi_hat_M.memptr(), (fftw_complex*) psi_M.memptr(), FFTW_BACKWARD, FFTW_PATIENT);

    
	
	//set up random forcing (white noise)
	unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();  
	mt19937_64 generator(seed1);
	normal_distribution<double> distribution(0.0,1.0);

    
	//read in initial data
    readData(psi, psi_hat, run_time_start, file_number_start, FFTN);      
   
    
	system("mkdir ./output");			//Creates the Folders for Outputting

	//records model parameters in file and outputs them to terminal
	recordParameters();

	//create the linear operator	
	linearOperator(L);     
   
  	//set up ETDRK-4 routine
	if(FLAG_TIMESTEP_METHOD == "ETDRK4"){
		setupETDRK(L,Q1,Q2,Q3,Q4,Q5,F1,F2,F3,E1,E2,dt);
	}

	//initialise forcing profile
	if(FLAG_FORCING_TYPE != "off"){
		cout << "initialising forcing" << endl;
		initialiseForcing(famp);
	}
	else{
		cout << "no forcing routine" << endl;
	}
	
	
	

    //============== time stepping loop ======================================================
    
    for( int time_step = 1; time_step <= total_steps; time_step++ ){
        
		//increments runtime
		run_time += dt;                  
		
		
		if(FLAG_FORCING_TYPE != "off" ){
		/*  This is the section of code that produces the forcing using normal distributed real and imaginary components at a given wave number*/

			noise.zeros();
				
			for(int i = 0; i < N; i++){
		    	if(famp(i) > 0.0){  
		       		noise(i) = famp(i)*pow(0.5,0.5)*complex<double>(distribution(generator),distribution(generator));   	     			
		    	}
		    }
		      	
		}
     

		if(FLAG_TIMESTEP_METHOD == "ETDRK4"){

			psi_hat += pow(dt,0.5)*noise; // add addative Euler method noise
	
			NL1 = nonlinear(psi_hat,0.0,L , psi_M, psi_hat_M, FFT, IFFT);	
			NL2 = nonlinear(E2%psi_hat + Q1%NL1,0.5*dt,L ,psi_M, psi_hat_M, FFT, IFFT); 	
			NL3 = nonlinear(E2%psi_hat + Q2%NL1 + Q3%NL2,0.5*dt,L , psi_M, psi_hat_M, FFT, IFFT);
			psi_hat = E1%psi_hat + F1%NL1 + 2.0*F2%(NL2+NL3) + F3%nonlinear(E1%psi_hat + Q4%NL1 + Q5%NL3,dt, L , psi_M, psi_hat_M, FFT, IFFT);

			if(FLAG_HYPO_DISSIPATION == true && alphapower < 0){
				psi_hat(0) = 0.0;
			}
				
		}
		else if( FLAG_TIMESTEP_METHOD == "RK2"){   
			//Runge-Kutta order 2 method  with noise
			NL1 = dt*nonlinear(psi_hat,0.0, L  ,psi_M, psi_hat_M, FFT, IFFT);
			NL2 = dt*nonlinear(psi_hat + NL1 + pow(dt,0.5)*noise,dt, L , psi_M, psi_hat_M, FFT, IFFT);
			psi_hat += 0.5*(NL1 + NL2) + pow(dt,0.5)*noise;
			
			//transforms data to u=exp(Ldt)v (integrating factor)
			psi_hat %= exp(dt*L);  

			if(FLAG_HYPO_DISSIPATION == true && alphapower < 0){
				psi_hat(0) = 0.0;
			}
    	
			    	
		}
		else{ //checks to make sure no other method has been defined
			cout << "FLAG_TIMETEP_METHOD invalid" << endl;
			exit(1);
		}
	 
	 	//printing of diagonistics and state
		if( int(run_time / output_time) == file_number + 1){

			//increments file number
			file_number++;
			//increments average number for averaged spectra
			avg_number++; 

			//saves state			
			psi_hat_temp = psi_hat; 

			//get physical state
			fftw_execute(IFFTN);
         
			//prints wave function			
			printWavefunction(file_number_start + file_number,psi); //computes intensity

			if(FLAG_OUTPUT_DIAGONOSTICS == true){
				
				//compute and prints waveaction
				computeWaveaction( psi_hat, waveaction);
				printWaveaction( run_time_start + run_time, file_number_start + file_number,  waveaction);

				//compute and prints energy
				computeEnergy(psi_hat,  psi,linear_energy,  potential_energy,  nonlinear_energy, psi_M, psi_hat_M,FFT, IFFT);
				printEnergy(run_time_start + run_time, file_number_start + file_number,  linear_energy, potential_energy, nonlinear_energy);

						
				//compute and print dissipation rates		
				computeDissipationRate(psi_hat, waveaction_dissipation_rate_alpha, waveaction_dissipation_rate_nu,energy_dissipation_rate_alpha,  energy_dissipation_rate_nu);
				printDissipationRate(run_time_start + run_time, file_number_start + file_number, waveaction_dissipation_rate_alpha, waveaction_dissipation_rate_nu, energy_dissipation_rate_alpha, energy_dissipation_rate_nu);

				//compute and print spectrum
				computeSpectrum(psi_hat, wave_spectrum);
				printSpectrum(file_number_start + file_number, avg_number, wave_spectrum, wave_spectrum_avg);

				//compute and print flux
				computeFlux(psi_hat,  L, psi_M , psi_hat_M,  wave_flux , energy_flux, FFT, IFFT );
				printFlux(file_number_start + file_number, avg_number, wave_flux,  energy_flux, wave_flux_avg, energy_flux_avg );
			}		
			
			//restore state (just in case)
			psi_hat = psi_hat_temp;
						
			//records some basic info to terminal
			cout << "time = " << run_time_start + run_time << " file = " << file_number_start + file_number << " Energy = " << linear_energy + potential_energy + nonlinear_energy << endl;        

			//overwrite curframe.dat
			ofstream fout_count("./data/curframe.dat");
			fout_count << scientific; fout_count.precision(12);
			fout_count << run_time_start + run_time  << " " << file_number_start + file_number << endl;    
		}
	}
	 
	
	//destroy fftw plans 
	fftw_destroy_plan(FFT);
	fftw_destroy_plan(IFFT);			
	fftw_destroy_plan(FFTN);
	fftw_destroy_plan(IFFTN);
	
	//end clock and display time
	finish = omp_get_wtime();                      
	cout << "Code has taken " << finish-start << " seconds" << endl;    
    

	return 0;
	
}
	

