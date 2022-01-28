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

static int ib = 0;
static rowvec wave_flux(NBIN),energy_flux(NBIN);
static double k_bin,k,fwave;


void flux(cx_rowvec A, cx_rowvec B, rowvec & wave_flux_avg, rowvec & energy_flux_avg, int out_count, int avg_count){

    
    k=0.0;
    wave_flux.zeros();
    energy_flux.zeros();
    k_bin = 2.0*pi / Lx;
    fwave=0.0;
   
        for(int i =0; i< N/2; i++){
        
        
        
            k =  2.0*pi * double(i)/ Lx;
            ib = int(k/k_bin);
            
            fwave = -2.0 * ( real(A(i))*real(B(i)) + imag(A(i))*imag(B(i)) );
  
            wave_flux(ib) += fwave;
            energy_flux(ib) += k*k*fwave;
    

            k = abs(  2.0*pi * double(-i-1)/ Lx);
            ib = int(k/k_bin);
            
            fwave = -2.0 * ( real(A(N-i-1))*real(B(N-i-1)) + imag(A(N-i-1))*imag(B(N-i-1)));

            wave_flux(ib) += fwave;
            energy_flux(ib) += k*k*fwave;

        
    }

    wave_flux_avg += wave_flux;
    energy_flux_avg += energy_flux;

    
    ostringstream out_flux;
    out_flux << "./output/flux." << setw(6) << setfill('0') << out_count << ends;			//creates file name for outputting data at time slice
    string filename = out_flux.str();
    ofstream fout_flux(filename.c_str());
    fout_flux << scientific;
    fout_flux.precision(12);
    
    for(int i=0; i < NBIN; i++){
            
            fout_flux << (i*k_bin) << " " << wave_flux(i) << " " << energy_flux(i) << " " << wave_flux_avg(i) / double(avg_count) << " " << energy_flux_avg(i) / double(avg_count) <<  endl;
    }
    
    fout_flux.close();
   
    
    return;
}
