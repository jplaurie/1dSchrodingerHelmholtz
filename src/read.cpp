//this code reads in the initial data from file ./Initial.dat

#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>
#include <fftw3.h>
#include <fstream>
#include <iomanip> 
#include "const.h"
#include "print.h"

using namespace std;
using namespace arma;

static double ignore_in, real_part, imag_part;

void readData(cx_rowvec & psi,cx_rowvec & psi_hat, double & run_time_start, int & file_number_start, fftw_plan FFTN){
    
/*
    reads in initial data located in file with number given by curframe.dat
*/

ignore_in = 0.0;
real_part = 0.0;
imag_part = 0.0;
    
fstream filein("./data/curframe.dat");
filein.precision(12);
filein >> run_time_start >> file_number_start;
    

if(file_number_start >= 0 && FLAG_INITIAL_CONDITION == "FILE"){

	ostringstream in_data;
	in_data << "./data/psi." << setw(6) << setfill('0') << file_number_start << ends;
	string filename = in_data.str();
	ifstream filein2(filename.c_str());
	filein2.precision(12);
    
	for(int i = 0; i < N; i++){
			filein2 >> ignore_in >> real_part >> imag_part;
			psi(i)= complex<double>(real_part,imag_part);
	}    

	cout << "Loading data from file psi." << setw(6) << setfill('0') << file_number_start << endl;
	cout << "Runtime of simulation = " << run_time_start << endl;
	cout << "Loading successful" << endl;
}
else if(file_number_start < 0 || FLAG_INITIAL_CONDITION == "ZERO"){


	file_number_start = 0;
	run_time_start = 0.0;
	psi.zeros();
	/*ofstream fout_read("./data/psi.000000");
	fout_read.precision(12);

	for(int i=0; i < N; i++){
		
			fout_read << i*dx << " " << real(psi(i)) << " " << imag(psi(i)) << endl;
			
	}
	fout_read.close();
*/
	printWavefunction(0 , psi );

	cout << "Simulation starting...from zero state" << endl;
	cout << "time = " << run_time_start << endl;
	cout << "file_number = 0" << endl; 

}
else{
	cout << "Problem with setting initial condition" << endl;
	exit(1);
}
    
fftw_execute(FFTN);
psi_hat /= double(N);

return;
}


void recordParameters(){


ofstream fout_data("./parameters.txt");

fout_data << "N = " << N <<  endl;
fout_data << "Lx = " << Lx << endl; 
fout_data << "c = " << c << " g = " << g << " mu = " << mu << " beta = " << beta_param << endl;;
fout_data <<"================ Dissipation =====================" << endl;
fout_data << "nupower = " << nupower << " nu = " << nu << endl;
fout_data << "alphapower = " << alphapower << " alpha = " << alpha << endl;
fout_data <<"================ Timestepping =========================" << endl;
fout_data << "FLAG_TIMESTEP_METHOD = " << FLAG_TIMESTEP_METHOD << endl;
fout_data << "dt = " << dt << endl;
fout_data << "total_steps = " << total_steps << endl;
fout_data << "output_time = " << output_time << endl;
fout_data << "FLAG_INITIAL_CONDITION = " << FLAG_INITIAL_CONDITION << endl;
fout_data <<"================ Forcing =========================" << endl;
fout_data << "FLAG_FORCING_TYPE = " << FLAG_FORCING_TYPE << endl;
if(FLAG_FORCING_TYPE == "annulus"){
	fout_data << "FLAG_FORCING_RESCALE = " << FLAG_FORCING_RESCALE << endl;
	fout_data << "force_amplitude = "  << force_amplitude << " kf = " << kf << " dk = " << dk << endl;
}
else if(FLAG_FORCING_TYPE == "exponential"){
	fout_data << "FLAG_FORCING_RESCALE = " << FLAG_FORCING_RESCALE << endl;
	fout_data << "force_amplitude = "  << force_amplitude << " force_power = " << force_power << endl;
}
else if(FLAG_FORCING_TYPE == "gaussian"){
	fout_data << "FLAG_FORCING_RESCALE = " << FLAG_FORCING_RESCALE << endl;
	fout_data << "force_amplitude = "  << force_amplitude << " kf = " << kf << " sigmaf = "<<  sigmaf << endl;
}


fout_data.close();

cout << "N = " << N <<  endl;
cout << "Lx = " << Lx << endl; 
cout << "c = " << c << " g = " << g << " mu = " << mu << " beta = " << beta_param << endl;;
cout <<"================ Dissipation =====================" << endl;
cout << "nupower = " << nupower << " nu = " << nu << endl;
cout << "alphapower = " << alphapower << " alpha = " << alpha << endl;
cout <<"================ Timestepping =========================" << endl;
cout << "FLAG_TIMESTEP_METHOD = " << FLAG_TIMESTEP_METHOD << endl;
cout << "dt = " << dt << endl;
cout << "total_steps = " << total_steps << endl;
cout << "output_time = " << output_time << endl;
cout << "FLAG_INITIAL_CONDITION = " << FLAG_INITIAL_CONDITION << endl;
cout <<"================ Forcing =========================" << endl;
cout << "FLAG_FORCING_TYPE = " << FLAG_FORCING_TYPE << endl;
if(FLAG_FORCING_TYPE == "annulus"){
	cout << "FLAG_FORCING_RESCALE = " << FLAG_FORCING_RESCALE << endl;
	cout << "force_amplitude = "  << force_amplitude << " kf = " << kf << " dk = " << dk << endl;
}
else if(FLAG_FORCING_TYPE == "exponential"){
	cout << "FLAG_FORCING_RESCALE = " << FLAG_FORCING_RESCALE << endl;
	cout << "force_amplitude = "  << force_amplitude << " force_power = " << force_power << endl;
}
else if(FLAG_FORCING_TYPE == "gaussian"){
	cout << "FLAG_FORCING_RESCALE = " << FLAG_FORCING_RESCALE << endl;
	cout << "force_amplitude = "  << force_amplitude << " kf = " << kf << " sigmaf = "<<  sigmaf << endl;
}
return;
}
