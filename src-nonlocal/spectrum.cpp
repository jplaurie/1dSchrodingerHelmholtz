

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
static rowvec wave_spec(NBIN),energy_spec(NBIN);
static double k_bin,k,psi2;


void spectrum( cx_rowvec A, rowvec & wave_spec_avg, rowvec & energy_spec_avg,int out_count, int avg_count){

    
k=0.0;
psi2=0.0;
    
wave_spec.zeros();
energy_spec.zeros();

k_bin = 2.0*pi / Lx;

	for(int i =0; i< N/2; i++){
        
		k =   (2.0*pi * double(i)/ Lx);
		
		ib = int(k/k_bin);
	
		psi2=pow( abs(A(i)) , 2.0);
		wave_spec(ib) += psi2;
		energy_spec(ib) += k*k*psi2;
      
		k = abs(2.0*pi * double(-i-1)/ Lx );
		ib = int(k/k_bin);
		psi2 = pow( abs(A(N-i-1)) , 2.0);
		wave_spec(ib) += psi2;
		energy_spec(ib) += k*k*psi2;
         
	}



wave_spec_avg += wave_spec;
energy_spec_avg += energy_spec;

    
ostringstream out_spec;
out_spec << "./output/spec." << setw(6) << setfill('0') << out_count << ends;			//creates file name for outputting data at time slice
string filename = out_spec.str();
ofstream fout_spec(filename.c_str());
fout_spec << scientific;
fout_spec.precision(12);
    
for(int i=0; i < NBIN; i++){
            
	fout_spec << (i*k_bin) << " " << wave_spec(i) << " " << energy_spec(i) << " " << wave_spec_avg(i) / double(avg_count) << " " << energy_spec_avg(i) / double(avg_count) <<  endl;
}
    
fout_spec.close();
        
return;

}
