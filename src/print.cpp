
#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>
#include "const.h"
#include <fftw3.h>
#include <fstream>
#include <iomanip>
using namespace std;
using namespace arma;


static double k_bin = (2.0*pi/Lx) ;


void printWavefunction(int file_number, cx_rowvec psi){
/* prints wave function in physical space */    
ostringstream out_I;
out_I << "./data/psi." << setw(6) << setfill('0') << file_number << ends;			//creates file name for outputting data at time slice
string filename = out_I.str();
ofstream fout_I(filename.c_str());
fout_I << scientific;
fout_I.precision(12);
    
for(int i=0; i < N; i++){
		fout_I << i*dx << " " << real(psi(i)) << " " << imag(psi(i)) << endl;
	}

fout_I.close();

return;
    
}


void printDissipationRate(double run_time, int file_number,  double waveaction_dissipation_rate_alpha, double waveaction_dissipation_rate_nu, double energy_dissipation_rate_alpha, double energy_dissipation_rate_nu){
 /* prints dissipation rates from dissipation terms */   

	ostringstream out_D;
	out_D << "./output/dissipation." << setw(6) << setfill('0') << file_number << ends;         //creates file name for outputting data at time slice
	string filenameD = out_D.str();
	ofstream fout_D(filenameD.c_str());
	fout_D << scientific;
	fout_D.precision(12);

	fout_D << run_time << " " << waveaction_dissipation_rate_alpha << " " << waveaction_dissipation_rate_nu << " " << energy_dissipation_rate_alpha << " " << energy_dissipation_rate_nu << endl;
    
	fout_D.close();

return;
    
}

void printEnergy(double run_time, int file_number,  double linear_energy, double potential_energy, double nonlinear_energy){
    /* prints energy */

	ostringstream out_E;
	out_E << "./output/energy." << setw(6) << setfill('0') << file_number << ends;			//creates file name for outputting data at time slice
	string filenameE = out_E.str();
	ofstream fout_E(filenameE.c_str());
	fout_E << scientific;
	fout_E.precision(12);
    
    fout_E << run_time << " " << linear_energy << " " << potential_energy << " " << nonlinear_energy << endl;
    
	fout_E.close();

return;
    
}


void printWaveaction(double run_time, int file_number,  double waveaction){
    /* prints total wave action */

	ostringstream out_WA;
	out_WA << "./output/waveaction." << setw(6) << setfill('0') << file_number << ends;			//creates file name for outputting data at time slice
	string filenameWA = out_WA.str();
	ofstream fout_WA(filenameWA.c_str());
	fout_WA << scientific;
	fout_WA.precision(12);
    
    fout_WA << run_time << " " << waveaction << endl;
    
	fout_WA.close();

return;
    
}



void printSpectrum(int file_number, int avg_number, rowvec wave_spectrum, rowvec & wave_spectrum_avg){
/* prints spectrum and averaged spectrum */
    
wave_spectrum_avg += wave_spectrum;
    
ostringstream out_spec;
out_spec << "./output/spectrum." << setw(6) << setfill('0') << file_number << ends;			//creates file name for outputting data at time slice
string filename = out_spec.str();
ofstream fout_spec(filename.c_str());
fout_spec << scientific;
fout_spec.precision(12);
    
for(int i=0; i < N/2+1; i++){
            
	fout_spec << (i*k_bin) << " " << wave_spectrum(i) << " " << wave_spectrum_avg(i) / double(avg_number) <<  endl;
}
    
fout_spec.close();
        
return;

}



void printFlux(int file_number, int avg_number, rowvec wave_flux, rowvec energy_flux, rowvec & wave_flux_avg,rowvec & energy_flux_avg ){
/* prints flux an averaged flux */
    
    wave_flux_avg += wave_flux;
    energy_flux_avg += energy_flux;

    
    ostringstream out_flux;
    out_flux << "./output/flux." << setw(6) << setfill('0') << file_number << ends;           //creates file name for outputting data at time slice
    string filename = out_flux.str();
    ofstream fout_flux(filename.c_str());
    fout_flux << scientific;
    fout_flux.precision(12);
    
    for(int i=0; i < N/2 +1; i++){
            
            fout_flux << (i*k_bin) << " " << wave_flux(i) << " " << energy_flux(i) << " " << wave_flux_avg(i) / double(avg_number) << " " << energy_flux_avg(i) / double(avg_number) <<  endl;
    }
    
    fout_flux.close();
        
return;

}
